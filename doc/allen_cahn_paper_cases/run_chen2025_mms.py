#!/usr/bin/env python3
"""Generate/run OOFEM cases for Chen et al. (2025), Section 6.2 MMS test."""
from pathlib import Path
import math
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "doc/allen_cahn_paper_cases/chen2025_mms"
OOFEM = ROOT / "build/opt_magneto/oofem"
L = 2.0 * math.pi


def make_case(name: str, n: int, dt: float, final_time: float) -> Path:
    case = OUT / name
    case.mkdir(parents=True, exist_ok=True)
    steps = round(final_time / dt)
    dt = final_time / steps
    node_id = lambda i, j: j * (n + 1) + i + 1
    nodes = []
    for j in range(n + 1):
        for i in range(n + 1):
            x, y = L * i / n, L * j / n
            nodes.append((x, y, 0.05 * math.sin(x) * math.sin(y)))
    elements = []
    eid = 1
    for j in range(n):
        for i in range(n):
            a, b = node_id(i, j), node_id(i + 1, j)
            c, d = node_id(i + 1, j + 1), node_id(i, j + 1)
            elements.extend([(eid, a, b, c), (eid + 1, a, c, d)])
            eid += 2
    boundary = [node_id(i, 0) for i in range(n + 1)]
    boundary += [node_id(i, n) for i in range(n + 1)]
    boundary += [node_id(0, j) for j in range(1, n)]
    boundary += [node_id(n, j) for j in range(1, n)]
    lines = [
        f"{name}.out",
        "Chen et al. (2025) Sec. 6.2 exact-solution Allen-Cahn MMS; P1 triangles, homogeneous Dirichlet data",
        f"nltransienttransportproblem nsteps {steps} deltat {dt:.16g} alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1",
        "vtkxml tstep_all domain_all primvars 1 6",
        "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {len(nodes)} nelem {len(elements)} ncrosssect 1 nmat 1 nbc 1 nic {len(nodes)} nltf 1 nset {len(nodes)+2}",
    ]
    for idx, (x, y, _) in enumerate(nodes, 1):
        lines.append(f"node {idx} coords 3 {x:.16g} {y:.16g} 0")
    for eid, a, b, c in elements:
        lines.append(f"tr1ht {eid} nodes 3 {a} {b} {c}")
    lines += [
        f"SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
        # Paper equation: u_t - eps^2 Delta(u) + (u^3-u) = g; eps=0.01.
        # Thus k=eps^2, mobility=1 in OOFEM's c*u_t-k*Delta(u)=mobility*(u-u^3)+source form.
        "allencahn 1 d 1.0 k 0.0001 c 1.0 mobility 1.0 manufactured 1 mmsbase 0.0 mmsamplitude 0.05 mmsrate 0.1 mmsbeta 0.0 mmslx 3.141592653589793 mmsly 3.141592653589793 mmstemporal 0 mmsspatial 1",
        f"BoundaryCondition 1 loadTimeFunction 1 dofs 1 10 values 1 0.0 set 2",
    ]
    for idx, (_, _, val) in enumerate(nodes, 1):
        sid = idx + 2
        lines.append(f"InitialCondition {idx} dofs 1 10 Conditions 1 u {val:.16g} set {sid}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append("Set 1 elementranges {(1 " + str(len(elements)) + ")}")
    lines.append(f"Set 2 nodes {len(boundary)} " + " ".join(map(str, boundary)))
    for idx in range(1, len(nodes) + 1):
        lines.append(f"Set {idx+2} nodes 1 {idx}")
    input_path = case / "case.in"
    input_path.write_text("\n".join(lines) + "\n")
    return input_path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "generate":
        # Temporal convergence: hold spatial error small; spatial convergence: hold dt small.
        cases = [(f"time_dt_{dt:g}", 32, dt, 1.0)
                 for dt in (0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.0009765625)]
        cases += [(f"space_n_{n}", n, 0.0025, 1.0) for n in (8, 16, 32, 64)]
        for args in cases:
            print(make_case(*args))
        return
    if len(sys.argv) > 1 and sys.argv[1] == "run_all":
        for case in sorted(OUT.iterdir()):
            if (case / "case.in").exists():
                with (case / "run.log").open("w") as log:
                    subprocess.run([str(OOFEM), "-f", "case.in"], cwd=case,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
                print(f"finished {case.name}")
        return
    if len(sys.argv) < 2:
        raise SystemExit("usage: run_chen2025_mms.py generate | <case-directory>")
    case = Path(sys.argv[1]).resolve()
    subprocess.run([str(OOFEM), "-f", str(case / "case.in")], cwd=case, check=True)


if __name__ == "__main__":
    main()
