#!/usr/bin/env python3
"""Reproduce a traveling-wave solution from Nizovtseva et al. (2016)."""
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OOFEM = ROOT / "build/debug_magneto/oofem"
TAU = 0.2                         # dimensionless relaxation time t_R M
T = 0.2
XMIN, XMAX = -12.0, 12.0
LAMBDA = math.sqrt(2.0 + 9.0 * TAU) / 2.0


def exact(x: np.ndarray, t: float) -> np.ndarray:
    # Nizovtseva et al., Eq. (29), Q=A=B=+1, c_0=0.
    z = np.clip(LAMBDA * x - 1.5 * t, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(z))


def make_input(nx: int, dt: float, outname: str = "case.out") -> str:
    ny = 1
    nn, ne = (nx + 1) * (ny + 1), nx * ny
    nsteps = round(T / dt)
    node = lambda i, j: j * (nx + 1) + i + 1
    lines = [
        outname, "Nizovtseva et al. hyperbolic Allen-Cahn traveling wave",
        f"nltransienttransportproblem nsteps {nsteps} deltat {dt:.16g} alpha 1.0 rtol 1.e-8 nsmax 100 nmodules 1 wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0 wavedamping {1.0/TAU:.16g}",
        "vtkxml tstep_all domain_all primvars 1 6", "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}",
    ]
    xnodes = []
    vnodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = XMIN + (XMAX - XMIN) * i / nx
            xnodes.append(x)
            u0 = float(exact(np.array([x]), 0.0)[0])
            vnodes.append(1.5 * u0 * (1.0 - u0))
            lines.append(f"node {node(i,j)} coords 3 {x:.16g} {float(j):.1f} 0.0")
    for i in range(nx):
        lines.append(f"quad1ht {i+1} nodes 4 {node(i,0)} {node(i+1,0)} {node(i+1,1)} {node(i,1)} nip 9")
    lines += [
        "SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
        f"allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 hyperbolic 1 rho {TAU:.16g} surfacecapacity 1.0",
    ]
    for j in range(ny + 1):
        for i in range(nx + 1):
            n = node(i, j)
            x = XMIN + (XMAX - XMIN) * i / nx
            u0 = float(exact(np.array([x]), 0.0)[0])
            lines.append(f"InitialCondition {n} dofs 1 10 Conditions 1 u {u0:.17g} set {n+1}")
    # The scalar transport model has one equation at every node, numbered in node order.
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {ne})}}")
    lines += [f"Set {n+1} nodes 1 {n}" for n in range(1, nn + 1)]
    lines[2] += f" initialvelocityvalues {nn} " + " ".join(f"{v:.17g}" for v in vnodes)
    return "\n".join(lines) + "\n"


def read_final(folder: Path) -> tuple[np.ndarray, np.ndarray]:
    pvd = folder / "case.out.m0.pvd"
    root = ET.parse(pvd).getroot()
    vtu = folder / root.findall(".//DataSet")[-1].attrib["file"]
    tree = ET.parse(vtu).getroot()
    coords = np.fromstring(tree.find(".//Points/DataArray").text, sep=" ").reshape(-1, 3)
    values = next(a for a in tree.findall(".//PointData/DataArray") if a.attrib.get("Name") == "Temperature")
    u = np.fromstring(values.text or "", sep=" ")
    return coords[:, 0], u


def run(nx: int, dt: float, tag: str) -> tuple[np.ndarray, np.ndarray, dict]:
    folder = HERE / f"nizovtseva_{tag}"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "case.in").write_text(make_input(nx, dt))
    result = subprocess.run([str(OOFEM), "-f", str(folder / "case.in")], cwd=folder,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (folder / "run.log").write_text(result.stdout)
    if result.returncode or "ANALYSIS FINISHED" not in result.stdout:
        raise RuntimeError(f"OOFEM failed for {tag}; see {folder/'run.log'}")
    x, u = read_final(folder)
    ref = exact(x, T)
    row = {"tag": tag, "nx": nx, "h": (XMAX-XMIN)/nx, "dt": dt, "T": T,
           "rms_error_vs_exact": float(np.sqrt(np.mean((u-ref)**2))),
           "max_error_vs_exact": float(np.max(np.abs(u-ref)))}
    print(row, flush=True)
    return x, u, row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="run one coarse case to check the input")
    args = parser.parse_args()
    if args.quick:
        _, _, row = run(80, 0.001, "quick")
        print(row)
        return
    if not OOFEM.exists():
        raise FileNotFoundError(OOFEM)

    rows = []
    temporal = []
    for dt in (0.004, 0.002, 0.001):
        x, u, row = run(160, dt, f"time_n160_dt{dt:g}")
        rows.append(row)
        temporal.append((dt, u))
    for coarse, fine in zip(temporal, temporal[1:]):
        delta = float(np.sqrt(np.mean((coarse[1]-fine[1])**2)))
        next(row for row in rows if row["dt"] == coarse[0])["rms_to_next_finer"] = delta
    r0, r1 = rows[0]["rms_to_next_finer"], rows[1]["rms_to_next_finer"]
    rows[0]["observed_temporal_order"] = math.log(r0/r1, 2.0)

    for nx in (80, 160, 320):
        _, _, row = run(nx, 0.001, f"space_n{nx}_dt0p001")
        rows.append(row)
    space_rows = [r for r in rows if r["tag"].startswith("space_")]
    space_rows.sort(key=lambda r: r["nx"])
    for a, b in zip(space_rows, space_rows[1:]):
        b["observed_spatial_order"] = math.log(a["rms_error_vs_exact"]/b["rms_error_vs_exact"], 2.0)

    fields = sorted({k for r in rows for k in r})
    with (HERE / "nizovtseva_traveling_wave_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
