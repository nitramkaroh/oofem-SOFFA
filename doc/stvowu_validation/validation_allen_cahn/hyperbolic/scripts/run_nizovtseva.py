#!/usr/bin/env python3
"""Reproduce the Nizovtseva et al. (2016) hyperbolic Allen-Cahn traveling wave.

Reference: I. G. Nizovtseva, P. K. Galenko, D. V. Alexandrov, "The
hyperbolic Allen-Cahn equation: exact solutions", J. Phys. A: Math. Theor.
49 (2016) 435201, eqs. (1)-(4) and the bounded traveling wave (29),
Table 1 (sign choice k = 1: Q = A = B = +1, c > 0), printed p. 6.

Dimensional PDE (checked against the PDF):

    tau_R u_tt + u_t = D Delta u + M (u - u^3),   F(u) = 1/4 (u^2 - 1)^2.

Dimensionless variables t' = M t, x' = x sqrt(M/D), tau = tau_R M.
With D = M = 1 (so x' = x, t' = t) and tau = tau_R = 0.2 the bounded
traveling wave (29) with c_0 = 0 reads

    u(x, t) = 1 / (1 + exp[ lambda (x - c t) ]),
    lambda = sqrt(2 + 9 tau) / 2,     c = 3 / (2 lambda) = 3/sqrt(2+9 tau)
           = 1.5389675281...,

which is verified here to satisfy the PDE to machine precision
(see check_exact_residual.py). The level u = 1/2 (the interface) moves as
x_0(t) = c t.

Translation to OOFEM (see the master README for the derivation): the
assembled bulk equation is  rho u_tt + rho*delta u_t - k Delta u - m(u-u^3) = 0,
so the input uses rho = 0.2, wavedamping = 5 (rho*delta = 1), k = D = 1,
mobility = m = M = 1, Newmark beta = 1/4, gamma = 1/2, and the exact
nodal initial velocity u_t(x,0) = 1.5 u(x,0) (1 - u(x,0)).

Geometry: strip [-12, 12] x [0, 1], one layer of Q1 cells along y, natural
Neumann on all edges (the exact normal derivative is exponentially small at
x = +/-12; the truncation effect is quantified in the report). No boundary
elements: the traveling wave does not carry a dynamic boundary condition,
and none of the available literature sources provides an exact solution for
the coupled inertial-bulk / dynamic-surface model implemented in OOFEM.

Run:  python3 run_nizovtseva.py            (runs n = 80, 160, 320)
      python3 run_nizovtseva.py --only 320
Outputs in ../results/nizovtseva_n{80,160,320}/ (case.in, run.log, .out,
.pvd, every VTU frame); the n320 input is also copied to ../input/.
"""
from __future__ import annotations
import argparse, math, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]          # repo root
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results'
OOFEM = ROOT / 'build' / 'debug_magneto' / 'oofem'

TAU = 0.2                 # = tau_R with M = 1
T = 0.2
D = MOBILITY = DAMPING = 1.0
RHO = TAU * MOBILITY      # 0.2
BETA, GAMMA = 0.25, 0.5
LAMBDA = math.sqrt(2.0 + 9.0 * TAU) / 2.0
SPEED = 1.5 / LAMBDA      # = 3/sqrt(2+9 tau)
DT = 0.001
XMIN, XMAX = -12.0, 12.0


def exact(x, t):
    z = LAMBDA * (x - SPEED * t)
    return 1.0 / (1.0 + math.exp(max(min(z, 700.0), -700.0)))


def exact_velocity(x, t):
    u = exact(x, t)
    return 1.5 * MOBILITY * u * (1.0 - u)


def make_input(nx: int) -> str:
    ny = 1
    nn, ne = (nx + 1) * 2, nx
    nsteps = round(T / DT)
    assert abs(nsteps * DT - T) < 1e-12
    node = lambda i, j: j * (nx + 1) + i + 1
    lines = [
        'case.out',
        'Nizovtseva et al. (2016), Eq. (29), bounded traveling wave, bulk only',
        (f'nltransienttransportproblem nsteps {nsteps} deltat {DT:.16g} alpha 1.0 '
         f'rtol 1.e-10 nsmax 100 nmodules 1 wave 1 newmarkbeta {BETA} newmarkgamma {GAMMA} '
         f'initialvelocity 0.0 wavedamping {DAMPING / RHO:.16g}'),
        'vtkxml tstep_all domain_all primvars 1 6',
        'domain heattransfer',
        'OutputManager tstep_all dofman_all element_all',
        f'ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn + 1}',
    ]
    values_u, values_v = [], []
    for j in range(2):
        for i in range(nx + 1):
            x = XMIN + (XMAX - XMIN) * i / nx
            lines.append(f'node {node(i, j)} coords 3 {x:.17g} {float(j):.1f} 0.0')
            values_u.append(exact(x, 0.0))
            values_v.append(exact_velocity(x, 0.0))
    for i in range(nx):
        lines.append(f'quad1ht {i + 1} nodes 4 {node(i, 0)} {node(i + 1, 0)} {node(i + 1, 1)} {node(i, 1)} nip 9')
    lines += [
        'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
        f'allencahn 1 d 1.0 k {D} c 1.0 mobility {MOBILITY} hyperbolic 1 rho {RHO} surfacecapacity 1.0',
    ]
    for n, u in enumerate(values_u, 1):
        lines.append(f'InitialCondition {n} dofs 1 10 Conditions 1 u {u:.17g} set {n + 1}')
    lines.append('ConstantFunction 1 f(t) 1.0')
    lines.append(f'Set 1 elementranges {{(1 {ne})}}')
    lines += [f'Set {n + 1} nodes 1 {n}' for n in range(1, nn + 1)]
    lines[2] += f' initialvelocityvalues {nn} ' + ' '.join(f'{v:.17g}' for v in values_v)
    return '\n'.join(lines) + '\n'


def run(nx: int) -> None:
    folder = RESULTS / f'nizovtseva_n{nx}'
    folder.mkdir(parents=True, exist_ok=True)
    inp = make_input(nx)
    (folder / 'case.in').write_text(inp)
    if nx == 320:
        (HERE.parent / 'input' / 'nizovtseva_n320.in').write_text(inp)
    p = subprocess.run([str(OOFEM), '-f', str(folder / 'case.in')], cwd=folder,
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (folder / 'run.log').write_text(p.stdout)
    (HERE.parent / 'logs' / f'nizovtseva_n{nx}.log').write_text(p.stdout)
    if p.returncode or 'ANALYSIS FINISHED' not in p.stdout:
        sys.exit(f'OOFEM failed for n={nx}; inspect {folder / "run.log"}')
    print(f'OK n={nx}: {folder}')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', type=int, default=None)
    args = ap.parse_args()
    if not OOFEM.exists():
        sys.exit(f'OOFEM executable not found: {OOFEM}')
    for nx in (80, 160, 320):
        if args.only is not None and nx != args.only:
            continue
        run(nx)


if __name__ == '__main__':
    main()
