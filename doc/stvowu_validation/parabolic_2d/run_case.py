#!/usr/bin/env python3
"""Generate and run a no-source Q1 Allen-Cahn bulk/dynamic-boundary case.

This is an exploratory Xiao--Chen Case 4 parameterization. Their paper omits
epsilon for Case 4, so epsilon is explicit here and is never inferred from it.
"""
from __future__ import annotations
import argparse, math, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OOFEM = ROOT / "build/debug_magneto/oofem"
HERE = Path(__file__).resolve().parent

def make_input(n: int, dt: float, end_time: float, out_dir: Path,
               output_every: int, eps: float, reverse_boundary: bool = False,
               permute_boundary: bool = False) -> Path:
    steps = round(end_time / dt)
    if abs(steps * dt - end_time) > 1.e-12:
        raise ValueError("end_time must be an integer multiple of dt")
    nodes = (n + 1) ** 2
    elems = n * n + 4 * n
    prefix = "case"
    lines = [
        f"{prefix}.out", f"Xiao-Chen Case 4, unpenalized exploratory epsilon={eps:g}",
        f"nltransienttransportproblem nsteps {steps} deltat {dt:.14g} alpha 1.0 rtol 1.e-10 nsmax 40 nmodules 1",
        f"vtkxml tstep_step {output_every} domain_all primvars 1 6", "domain heattransfer",
        f"OutputManager tstep_step {output_every} dofman_all element_all",
        f"ndofman {nodes} nelem {elems} ncrosssect 1 nmat 1 nbc 0 nic {nodes} nltf 1 nset {nodes + 1}"
    ]
    for j in range(n + 1):
        for i in range(n + 1):
            node = j * (n + 1) + i + 1
            lines.append(f"node {node} coords 3 {i/n:.14g} {j/n:.14g} 0")
    eid = 1
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i + 1; b = a + 1
            d = a + n + 1; c = d + 1
            lines.append(f"quad1ht {eid} nodes 4 {a} {b} {c} {d} nip 9"); eid += 1
    boundary = []
    boundary += [(i + 1, i + 2) for i in range(n)]
    boundary += [(j * (n + 1) + n + 1, (j + 1) * (n + 1) + n + 1) for j in range(n)]
    boundary += [(n * (n + 1) + i + 1, n * (n + 1) + i) for i in range(n, 0, -1)]
    boundary += [(j * (n + 1) + 1, (j - 1) * (n + 1) + 1) for j in range(n, 0, -1)]
    expected={(i+1,i+2) for i in range(n)}
    expected|={(j*(n+1)+n+1,(j+1)*(n+1)+n+1) for j in range(n)}
    expected|={(n*(n+1)+i+1,n*(n+1)+i) for i in range(n,0,-1)}
    expected|={(j*(n+1)+1,(j-1)*(n+1)+1) for j in range(n,0,-1)}
    canonical={tuple(sorted(e)) for e in boundary}
    if len(boundary)!=4*n or len(canonical)!=4*n or canonical!={tuple(sorted(e)) for e in expected}:
        raise ValueError('Boundary line mesh must cover every square edge exactly once')
    if reverse_boundary:
        boundary = [(b, a) for a, b in boundary]
    if permute_boundary:
        boundary = boundary[::2] + boundary[1::2]
    for a, b in boundary:
        lines.append(f"line1ht {eid} nodes 2 {a} {b} nip 3"); eid += 1
    lines += [
        "SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
        f"allencahn 1 d 1.0 k 0.01 c 1.0 mobility {0.01 / eps**2:.14g} surfacebeta 1.0"
    ]
    # Piecewise constant +/-1 square droplet, matching the figure's initial shape.
    for j in range(n + 1):
        y = j / n
        for i in range(n + 1):
            x = i / n
            value = 1.0 if 0.25 <= x <= 0.75 and 0.0 <= y <= 0.5 else -1.0
            node = j * (n + 1) + i + 1
            lines.append(f"InitialCondition {node} dofs 1 10 Conditions 1 u {value:g} set {node + 1}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {elems})}}")
    for node in range(1, nodes + 1):
        lines.append(f"Set {node + 1} nodes 1 {node}")
    out_dir.mkdir(parents=True, exist_ok=True)
    inp = out_dir / "case.in"
    inp.write_text("\n".join(lines) + "\n")
    return inp

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--dt", type=float, required=True)
    p.add_argument("--end-time", type=float, default=.1)
    p.add_argument("--eps", type=float, default=.1)
    p.add_argument("--output-every", type=int, default=1)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reverse-boundary", action="store_true")
    p.add_argument("--permute-boundary", action="store_true")
    a = p.parse_args()
    out_dir=a.out.resolve()
    inp = make_input(a.n, a.dt, a.end_time, out_dir, a.output_every, a.eps,
                     a.reverse_boundary, a.permute_boundary)
    run = subprocess.run([str(OOFEM), "-f", str(inp)], cwd=out_dir,
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (out_dir / "run.log").write_text(run.stdout)
    print(run.stdout[-2500:])
    if run.returncode or "0 error(s) and 0 warning(s)" not in run.stdout:
        raise SystemExit(f"OOFEM failed; inspect {out_dir / 'run.log'}")
    print(f"Completed n={a.n}, dt={a.dt:g}, T={a.end_time:g}; output={out_dir}")

if __name__ == "__main__":
    main()
