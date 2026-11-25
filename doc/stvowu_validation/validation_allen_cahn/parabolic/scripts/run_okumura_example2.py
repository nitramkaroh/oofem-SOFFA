#!/usr/bin/env python3
"""Reproduce Okumura dissertation (2021), Ch. 4, Computation example 2, in OOFEM.

Problem (thesis eqs. (4.1), (4.91), section 6.2, printed p. 70-71):

    u_t = u_xx - F'(u)          in (0, L) x (0, T],      L = 1,
    eps_ex * u_t(0,t)  =  u_x(0,t)  - F'(u(0,t)),
    eps_ex * u_t(L,t) = -u_x(L,t)  - F'(u(L,t)),

    F(s) = (gamma/4) (s^2 - 1)^2,  gamma = 100,  eps_ex = 10,
    u(x,0) = exp(-500 (x - 0.5)^2),
    K = 100 (dx = 0.01),  N = 700 (dt = 1e-3),  T = 0.7.

Implementation in OOFEM (see ../../README.md for the full derivation):

* the 1D interval is extruded to a unit-depth 2D strip of 100x1 Q1 cells
  (``quad1ht``, 3x3 Gauss); a y-independent field then reduces exactly to
  the 1D problem (verified a posteriori: max variation across the depth is
  below export precision);
* the two 0D endpoints of the 1D interval are represented by two
  ``line1ht`` edges of unit height carrying the ``allencahn`` material with
  capacity ``c=10`` (= eps_ex) on the endpoint nodes, which assembles the
  first-order dynamic law 10 psi_t = -d_n u - F'(psi) on each end (the
  tangential Laplacian of a 0D endpoint is zero);
* top and bottom edges of the strip are natural Neumann edges, matching the
  y-independence of the solution.

Bulk material:  allencahn 1 d 1.0 k 1.0 c 1.0 mobility 100.0
End material:   allencahn 2 d 1.0 k 1.0 c 10.0 mobility 100.0 surfacebeta 1.0

Run:  python3 run_okumura_example2.py
Outputs are written to ../results/okumura_ch4_case2/ (input, log, .out,
.pvd and every VTU frame) and to ../logs/.
"""
from __future__ import annotations
import math, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]          # repo root
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'okumura_ch4_case2'
LOG = HERE.parent / 'logs'
OOFEM = ROOT / 'build' / 'debug_magneto' / 'oofem'

NX, NY = 100, 1
DT, T = 0.001, 0.7
NSTEPS = round(T / DT)
NN = (NX + 1) * (NY + 1)
NE = NX * NY + 2


def make_input() -> str:
    lines = [
        'okumura_ch4_case2.out',
        'Okumura thesis (2021), Chapter 4, Computation Example 2',
        f'nltransienttransportproblem nsteps {NSTEPS} deltat {DT} alpha 1.0 rtol 1.e-10 nsmax 30 nmodules 1',
        'vtkxml tstep_all domain_all primvars 1 6',
        'domain heattransfer',
        'OutputManager tstep_all dofman_all element_all',
        f'ndofman {NN} nelem {NE} ncrosssect 2 nmat 2 nbc 0 nic {NX + 1} nltf 1 nset {NX + 3}',
    ]
    for j in range(NY + 1):
        for i in range(NX + 1):
            node = j * (NX + 1) + i + 1
            lines.append(f'node {node} coords 3 {i / NX:.12g} {j / NY:.12g} 0.0')
    eid = 1
    for j in range(NY):
        for i in range(NX):
            a = j * (NX + 1) + i + 1
            b = a + 1
            d = a + NX + 1
            c = d + 1
            lines.append(f'quad1ht {eid} nodes 4 {a} {b} {c} {d} nip 9')
            eid += 1
    # Dynamic endpoints of the 1D interval, represented as unit-height edges.
    lines.append(f'line1ht {eid} nodes 2 1 {NX + 2} nip 3')
    eid += 1
    lines.append(f'line1ht {eid} nodes 2 {NX + 1} {NN} nip 3')
    lines += [
        'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
        'SimpleTransportCS 2 thickness 1.0 area 1.0 mat 2 set 2',
        'allencahn 1 d 1.0 k 1.0 c 1.0 mobility 100.0',
        'allencahn 2 d 1.0 k 1.0 c 10.0 mobility 100.0 surfacebeta 1.0',
    ]
    for i in range(NX + 1):
        value = math.exp(-500.0 * (i / NX - .5) ** 2)
        n0 = i + 1
        n1 = NX + 1 + i + 1
        lines.append(f'InitialCondition {i + 1} dofs 1 10 Conditions 1 u {value:.16g} set {i + 3}')
    lines += [
        'ConstantFunction 1 f(t) 1.0',
        f'Set 1 elementranges {{(1 {NX})}}',
        f'Set 2 elements 2 {NX + 1} {NX + 2}',
    ]
    for i in range(NX + 1):
        n0 = i + 1
        n1 = NX + 1 + i + 1
        lines.append(f'Set {i + 3} nodes 2 {n0} {n1}')
    return '\n'.join(lines) + '\n'


def main() -> None:
    if not OOFEM.exists():
        sys.exit(f'OOFEM executable not found: {OOFEM}')
    OUT.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    inp = make_input()
    (OUT / 'case.in').write_text(inp)
    (HERE.parent / 'input' / 'okumura_ch4_case2.in').write_text(inp)
    p = subprocess.run([str(OOFEM), '-f', str(OUT / 'case.in')], cwd=OUT,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT / 'run.log').write_text(p.stdout)
    (LOG / 'okumura_ch4_case2.log').write_text(p.stdout)
    tail = p.stdout.splitlines()[-12:]
    print('\n'.join(tail))
    if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout:
        sys.exit(f'OOFEM failed; inspect {OUT / "run.log"}')
    print(f'OK: Nx={NX}, dt={DT}, T={T}; raw outputs in {OUT}')


if __name__ == '__main__':
    main()
