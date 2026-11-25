#!/usr/bin/env python3
"""No-MMS verification of the hyperbolic bulk/Newmark path in OOFEM.

Two bulk-only Neumann problems are run on a one-cell-thick Q1 strip:
1. a linear cosine eigenmode (mobility=0), checked against its exact
   consistent-mass semi-discrete solution;
2. a spatially uniform nonlinear Allen-Cahn oscillator, checked against an
   independent high-accuracy integration of rho*u''=M*(u-u**3).

The script writes inputs, logs, VTU files and a CSV below this directory.
"""
from __future__ import annotations

import csv
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OOFEM = ROOT / "build/debug_magneto/oofem"
NX, NY, T = 40, 1, 0.2
BETA, GAMMA = 0.25, 0.5


def make_case(case: str, dt: float) -> tuple[str, np.ndarray, np.ndarray]:
    steps = round(T / dt)
    nn = (NX + 1) * (NY + 1)
    ne = NX * NY
    lines = [
        f"{case}_dt{dt:g}.out", f"Hyperbolic bulk check: {case}",
        f"nltransienttransportproblem nsteps {steps} deltat {dt:.16g} alpha 1.0 rtol 1.e-11 nsmax 50 nmodules 1 wave 1 newmarkbeta {BETA} newmarkgamma {GAMMA} initialvelocity {0.3 if case == 'oscillator' else 0.0}",
        "vtkxml tstep_all domain_all primvars 1 6", "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic 1 nltf 1 nset 2",
    ]
    node_id = lambda i, j: j * (NX + 1) + i + 1
    for j in range(NY + 1):
        for i in range(NX + 1):
            x, y = i / NX, j / NY
            lines.append(f"node {node_id(i,j)} coords 3 {x:.16g} {y:.16g} 0.0")
    eid = 1
    for j in range(NY):
        for i in range(NX):
            n1, n2 = node_id(i,j), node_id(i+1,j)
            n3, n4 = node_id(i+1,j+1), node_id(i,j+1)
            lines.append(f"quad1ht {eid} nodes 4 {n1} {n2} {n3} {n4} nip 9")
            eid += 1
    lines.extend(["SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1"])
    mobility = 0.0 if case == "linear" else 1.0
    lines.append(f"allencahn 1 d 1.0 k 1.0 c 1.0 mobility {mobility} hyperbolic 1 rho 1.0 surfacecapacity 1.0")
    lines.append("InitialCondition 1 dofs 1 10 Conditions 1 u 0.0 set 2" if case == "linear" else
                 "InitialCondition 1 dofs 1 10 Conditions 1 u 0.2 set 2")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {ne})}}")
    lines.append(f"Set 2 nodes {nn} " + " ".join(str(i) for i in range(1, nn+1)))
    if case == "linear":
        # Override the common initial condition with one value per node.
        lines[6] = f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}"
        lines = [s for s in lines if not s.startswith("InitialCondition 1 ")]
        lines = [s for s in lines if not s.startswith("Set 1 elementranges") and not s.startswith("Set 2 nodes")]
        # Insert nodal initial conditions and one-node sets after material.
        matpos = next(i for i, s in enumerate(lines) if s.startswith("allencahn ")) + 1
        icsets=[]
        for j in range(NY+1):
            for i in range(NX+1):
                n=node_id(i,j); val=math.cos(math.pi*i/NX)
                set_id=n+1
                icsets.append(f"InitialCondition {n} dofs 1 10 Conditions 1 u {val:.17g} set {set_id}")
        lines[matpos:matpos] = icsets
        lines += [f"Set 1 elementranges {{(1 {ne})}}"]
        lines += [f"Set {n+1} nodes 1 {n}" for n in range(1,nn+1)]
        initial = np.array([math.cos(math.pi*i/NX) for j in range(NY+1) for i in range(NX+1)])
    else:
        initial = np.full(nn, 0.2)
    text = "\n".join(lines) + "\n"
    return text, initial, np.array([i/NX for j in range(NY+1) for i in range(NX+1)])


def final_vtu(path: Path) -> np.ndarray:
    pvd = path.with_suffix(".out.m0.pvd")
    # Input's declared output basename is the full path's .out basename.
    if not pvd.exists():
        pvd = path.parent / (path.stem + ".out.m0.pvd")
    root = ET.parse(pvd).getroot()
    datasets = root.findall(".//DataSet")
    vtu = pvd.parent / datasets[-1].attrib["file"]
    tree = ET.parse(vtu)
    arrays = tree.findall(".//PointData/DataArray")
    if not arrays:
        raise RuntimeError(f"No point data in {vtu}")
    arr = next((a for a in arrays if a.attrib.get("Name", "").lower() in ("u", "temperature", "unknown")), arrays[0])
    return np.fromstring(arr.text or "", sep=" ")


def main() -> None:
    if not OOFEM.exists():
        raise FileNotFoundError(OOFEM)
    OUT.mkdir(parents=True, exist_ok=True)
    rows=[]
    dts=(0.04,0.02,0.01)
    for case in ("linear","oscillator"):
        for dt in dts:
            folder=OUT / f"{case}_dt{dt:g}"
            folder.mkdir(exist_ok=True)
            input_text, initial, x=make_case(case,dt)
            inp=folder / "case.in"
            # Use a local basename so OOFEM places all outputs with this case.
            input_text=input_text.replace(f"{case}_dt{dt:g}.out", "case.out", 1)
            inp.write_text(input_text)
            run=subprocess.run([str(OOFEM), "-f", str(inp)], cwd=folder,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (folder/"run.log").write_text(run.stdout)
            if run.returncode or "ANALYSIS FINISHED" not in run.stdout:
                raise RuntimeError(f"OOFEM failed for {case} dt={dt}; see {folder/'run.log'}")
            u=final_vtu(folder/"case.out")
            if u.size != initial.size:
                raise RuntimeError(f"Unexpected field size {u.size}, expected {initial.size}")
            if case == "linear":
                # Exact eigenvalue of the 1D consistent-mass Q1 FE operator.
                h=1/NX
                lam=6*(1-math.cos(math.pi/NX))/(h*h*(2+math.cos(math.pi/NX)))
                ref=initial*math.cos(math.sqrt(lam)*T)
                err=u-ref
                label="error_vs_exact_semidiscrete"
            else:
                def ode(t,y): return (y[1], y[0]-y[0]**3)
                ref=solve_ivp(ode,(0,T),(0.2,0.3),rtol=2e-13,atol=2e-14,
                              method="DOP853",dense_output=True).sol(T)[0]
                err=u-ref
                label="error_vs_independent_ode"
            rows.append({"case":case,"dt":dt,"T":T,"max_abs_error":float(np.max(np.abs(err))),
                         "rms_error":float(np.sqrt(np.mean(err**2))),"final_min":float(u.min()),
                         "final_max":float(u.max()),"reference":float(ref if np.isscalar(ref) else np.max(np.abs(ref))),
                         "reference_type":label,"observed_order":""})
            print(rows[-1])
    for case in ("linear", "oscillator"):
        ordered=sorted((r for r in rows if r["case"] == case), key=lambda r:r["dt"], reverse=True)
        for coarse, fine in zip(ordered, ordered[1:]):
            fine["observed_order"] = math.log(coarse["rms_error"]/fine["rms_error"], 2)
        ordered[-1]["observed_order"] = ""
    with (OUT/"bulk_reference_metrics.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
