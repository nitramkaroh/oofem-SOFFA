#!/usr/bin/env python3
"""Spatial refinement study for the Okumura Ch.4 Example 2 reproduction.

Same problem and time discretization as run_okumura_example2.py
(dt = 1e-3, T = 0.7, backward Euler), on Q1 meshes with NX = 50, 100, 200,
400 cells along the interval (unit-depth strip). The nodal grids nest
(x_i^50 = x_{2i}^100 = x_{4i}^200 = x_{8i}^400), so postprocessing can
compare the coarser solutions against the finest (NX = 400) solution
nodally, without any interpolation.

Run:  python3 run_okumura_refinement.py [--only 50,400]
Raw outputs: ../results/okumura_refinement_n{50,100,200}/
"""
from __future__ import annotations
import argparse, math, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]          # repo root
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results'
LOG = HERE.parent / 'logs'
OOFEM = ROOT / 'build' / 'debug_magneto' / 'oofem'

NY = 1
DT, T = 0.001, 0.7
NSTEPS = round(T / DT)


def make_input(nx: int) -> str:
    ny = NY
    nn = (nx + 1) * (ny + 1)
    ne = nx * ny + 2
    lines = [
        'case.out',
        'Okumura Ch.4 Example 2, spatial refinement (Q1 strip)',
        (f'nltransienttransportproblem nsteps {NSTEPS} deltat {DT} alpha 1.0 '
         f'rtol 1.e-10 nsmax 30 nmodules 1'),
        'vtkxml tstep_all domain_all primvars 1 6',
        'domain heattransfer',
        'OutputManager tstep_all dofman_all element_all',
        f'ndofman {nn} nelem {ne} ncrosssect 2 nmat 2 nbc 0 nic {nx + 1} '
        f'nltf 1 nset {nx + 3}',
    ]
    for j in range(ny + 1):
        for i in range(nx + 1):
            node = j * (nx + 1) + i + 1
            lines.append(f'node {node} coords 3 {i / nx:.12g} {j / ny:.12g} 0.0')
    eid = 1
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i + 1
            b = a + 1
            d = a + nx + 1
            c = d + 1
            lines.append(f'quad1ht {eid} nodes 4 {a} {b} {c} {d} nip 9')
            eid += 1
    lines.append(f'line1ht {eid} nodes 2 1 {nx + 2} nip 3')
    eid += 1
    lines.append(f'line1ht {eid} nodes 2 {nx + 1} {nn} nip 3')
    lines += [
        'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
        'SimpleTransportCS 2 thickness 1.0 area 1.0 mat 2 set 2',
        'allencahn 1 d 1.0 k 1.0 c 1.0 mobility 100.0',
        'allencahn 2 d 1.0 k 1.0 c 10.0 mobility 100.0 surfacebeta 1.0',
    ]
    for i in range(nx + 1):
        value = math.exp(-500.0 * (i / nx - .5) ** 2)
        n0 = i + 1
        n1 = nx + 1 + i + 1
        lines.append(f'InitialCondition {i + 1} dofs 1 10 Conditions 1 u '
                     f'{value:.16g} set {i + 3}')
    lines += [
        'ConstantFunction 1 f(t) 1.0',
        f'Set 1 elementranges {{(1 {nx})}}',
        f'Set 2 elements 2 {nx + 1} {nx + 2}',
    ]
    for i in range(nx + 1):
        n0 = i + 1
        n1 = nx + 1 + i + 1
        lines.append(f'Set {i + 3} nodes 2 {n0} {n1}')
    return '\n'.join(lines) + '\n'


def run(nx: int) -> None:
    folder = RESULTS / f'okumura_refinement_n{nx}'
    folder.mkdir(parents=True, exist_ok=True)
    inp = make_input(nx)
    (folder / 'case.in').write_text(inp)
    p = subprocess.run([str(OOFEM), '-f', str(folder / 'case.in')], cwd=folder,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (folder / 'run.log').write_text(p.stdout)
    (LOG / f'okumura_refinement_n{nx}.log').write_text(p.stdout)
    if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout:
        sys.exit(f'OOFEM failed for n={nx}; inspect {folder / "run.log"}')
    print(f'OK n={nx}: {folder}', flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', type=str, default=None,
                    help='comma-separated subset of 50,100,200,400')
    args = ap.parse_args()
    if not OOFEM.exists():
        sys.exit(f'OOFEM executable not found: {OOFEM}')
    only = None if args.only is None else {int(s) for s in args.only.split(',')}
    for nx in (50, 100, 200, 400):
        if only is not None and nx not in only:
            continue
        run(nx)


if __name__ == '__main__':
    main()
