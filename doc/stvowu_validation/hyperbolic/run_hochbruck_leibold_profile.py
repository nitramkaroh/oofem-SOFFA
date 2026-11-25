#!/usr/bin/env python3
"""Use the exact profile from Hochbruck--Leibold in the OOFEM wave model.

The source terms are the residual of OOFEM's current bulk-inertial,
first-order dynamic-surface operator. This is a code-verification adaptation,
not a reproduction of the paper's second-order acoustic boundary equation.
"""
from __future__ import annotations

import csv
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OOFEM = ROOT / "build/debug_magneto/oofem"
T = 0.2


def make_input(n: int, dt: float) -> str:
    nn = (n + 1) ** 2
    quads = n * n
    edges = 4 * n
    ne = quads + edges
    steps = round(T / dt)
    node = lambda i, j: j * (n + 1) + i + 1
    lines = [
        "case.out", "Hochbruck-Leibold exact profile adapted to OOFEM operator",
        f"nltransienttransportproblem nsteps {steps} deltat {dt:.16g} alpha 1.0 rtol 1.e-9 nsmax 60 nmodules 1 wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0 wavedamping 1.0",
        "vtkxml tstep_all domain_all primvars 1 6", "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}",
    ]
    velocities = []
    for j in range(n + 1):
        y = j / n
        for i in range(n + 1):
            x = i / n
            nd = node(i, j)
            lines.append(f"node {nd} coords 3 {x:.16g} {y:.16g} 0.0")
            velocities.append(2.0 * math.pi * x * y)
    eid = 1
    for j in range(n):
        for i in range(n):
            lines.append(f"quad1ht {eid} nodes 4 {node(i,j)} {node(i+1,j)} {node(i+1,j+1)} {node(i,j+1)} nip 9")
            eid += 1
    # Boundary lines in a consistent orientation, as in the rectangular dynamic-boundary cases.
    for i in range(n):
        lines.append(f"line1ht {eid} nodes 2 {node(i,0)} {node(i+1,0)} nip 3"); eid += 1
    for j in range(n):
        lines.append(f"line1ht {eid} nodes 2 {node(n,j)} {node(n,j+1)} nip 3"); eid += 1
    for i in range(n, 0, -1):
        lines.append(f"line1ht {eid} nodes 2 {node(i,n)} {node(i-1,n)} nip 3"); eid += 1
    for j in range(n, 0, -1):
        lines.append(f"line1ht {eid} nodes 2 {node(0,j)} {node(0,j-1)} nip 3"); eid += 1
    lines += [
        "SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
        "allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 hyperbolic 1 rho 1.0 surfacecapacity 1.0 manufactured 1 mmsbase 0.0 mmsamplitude 1.0 mmsrate 6.283185307179586 mmsbeta 0.0 mmslx 1.0 mmsly 1.0 mmstemporal 1 mmsspatial 2",
    ]
    for j in range(n + 1):
        for i in range(n + 1):
            nd = node(i, j)
            lines.append(f"InitialCondition {nd} dofs 1 10 Conditions 1 u 0.0 set {nd+1}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {ne})}}")
    lines += [f"Set {nd+1} nodes 1 {nd}" for nd in range(1, nn+1)]
    lines[2] += f" initialvelocityvalues {nn} " + " ".join(f"{v:.17g}" for v in velocities)
    return "\n".join(lines) + "\n"


def final_field(folder: Path) -> tuple[np.ndarray, np.ndarray]:
    pvd = folder / "case.out.m0.pvd"
    root = ET.parse(pvd).getroot()
    vtu = folder / root.findall(".//DataSet")[-1].attrib["file"]
    tree = ET.parse(vtu).getroot()
    coords = np.fromstring(tree.find(".//Points/DataArray").text, sep=" ").reshape(-1, 3)
    arr = next(a for a in tree.findall(".//PointData/DataArray") if a.attrib.get("Name") == "Temperature")
    return coords, np.fromstring(arr.text or "", sep=" ")


def run(n: int, dt: float, tag: str) -> tuple[np.ndarray, np.ndarray, dict]:
    folder = HERE / f"hl_profile_{tag}"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "case.in").write_text(make_input(n, dt))
    result = subprocess.run([str(OOFEM), "-f", str(folder / "case.in")], cwd=folder,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (folder / "run.log").write_text(result.stdout)
    if result.returncode or "ANALYSIS FINISHED" not in result.stdout:
        raise RuntimeError(f"OOFEM failed for {tag}; see {folder/'run.log'}")
    coords, u = final_field(folder)
    t = T
    exact = np.sin(2.0 * math.pi * t) * coords[:, 0] * coords[:, 1]
    row = {"tag": tag, "n": n, "h": 1.0/n, "dt": dt, "T": T,
           "rms_error_vs_exact_profile": float(np.sqrt(np.mean((u-exact)**2))),
           "max_error_vs_exact_profile": float(np.max(np.abs(u-exact)))}
    print(row, flush=True)
    return coords, u, row


def main() -> None:
    rows = []
    for dt in (0.01, 0.005, 0.0025):
        _, _, row = run(20, dt, f"time_n20_dt{dt:g}")
        rows.append(row)
    # Re-read the final fields to calculate successive-step differences on the same mesh.
    fields = [(dt, final_field(HERE / f"hl_profile_time_n20_dt{dt:g}")[1]) for dt in (0.01, 0.005, 0.0025)]
    for (dt0, u0), (dt1, u1) in zip(fields, fields[1:]):
        next(r for r in rows if r["dt"] == dt0)["rms_to_next_finer"] = float(np.sqrt(np.mean((u0-u1)**2)))
    d0, d1 = rows[0]["rms_to_next_finer"], rows[1]["rms_to_next_finer"]
    rows[0]["observed_temporal_order"] = math.log(d0/d1, 2.0)

    for n in (10, 20, 40):
        _, _, row = run(n, 0.0025, f"space_n{n}_dt0p0025")
        rows.append(row)
    space = [r for r in rows if r["tag"].startswith("space_")]
    space.sort(key=lambda r: r["n"])
    for a, b in zip(space, space[1:]):
        b["observed_spatial_order"] = math.log(a["rms_error_vs_exact_profile"]/b["rms_error_vs_exact_profile"], 2.0)

    columns = sorted({key for row in rows for key in row})
    with (HERE / "hochbruck_leibold_profile_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
