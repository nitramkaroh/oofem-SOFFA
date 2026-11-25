#!/usr/bin/env python3
"""PDE-residual verification of the Nizovtseva et al. (2016) Eq. (29)
traveling-wave profile used by validate_nizovtseva_bulk.py.

Verifies, without any OOFEM simulation, that

    f(x,t) = 1 / (1 + exp( lam*x - s*t )),   s = 3/2,   lam = sqrt(2 + 9*tau)/2

satisfies the 1-D hyperbolic Allen-Cahn equation

    tau * f_tt + f_t = f_xx + (f - f^3)        (D = M = 1, damping coeff 1)

i.e. rho = tau*M = 0.2 and wavedamping = 5 in the OOFEM input
(rho*zeta = 1), and that the initial velocity v0 = +s*f*(1-f) matches f_t.

Symbolic result (sympy, optional): the residual factors as
    g(1-g) [ A (1-g) + B g ],   g = f
with A = -lam^2 + tau s^2 + s - 1, B = lam^2 - tau s^2 + s - 2,
which vanish simultaneously iff s = 3/2 and lam^2 = (9 tau + 2)/4.
"""
import csv
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

TAU = 0.2
LAM = np.sqrt(2 + 9 * TAU) / 2
S = 1.5
XMIN, XMAX, TEND, NX, NT = -12.0, 12.0, 0.2, 2401, 41


def pde_residual(tau, lam, s, x, t):
    phi = lam * x - s * t
    f = 1.0 / (1.0 + np.exp(np.clip(phi, -60, 60)))
    gp = -f * (1 - f)          # d/dphi
    gpp = f * (1 - f) * (1 - 2 * f)
    f_t = -s * gp
    f_tt = s * s * gpp
    f_xx = lam * lam * gpp
    return tau * f_tt + f_t - f_xx - (f - f ** 3)


def symbolic_factor():
    try:
        import sympy as sp
    except ImportError:
        return None
    x, t, tau, s, lam = sp.symbols('x t tau s lam', positive=True)
    E = sp.exp(lam * x - s * t)      # E = e^phi, f = 1/(1+E)
    f = 1 / (1 + E)
    res = tau * sp.diff(f, t, 2) + sp.diff(f, t) - sp.diff(f, x, 2) - (f - f ** 3)
    A = -lam ** 2 + tau * s ** 2 + s - 1
    B = lam ** 2 - tau * s ** 2 + s - 2
    # with g = f: f_t = -s g' , f_tt = s^2 g'', g' = -g(1-g), g'' = g(1-g)(1-2g)
    # => residual = g(1-g)[A(1-g) + B g] = E (A E + B) / (1+E)^3
    expected = E * (A * E + B) / (1 + E) ** 3
    ok = sp.simplify(sp.together(res - expected)) == 0
    sol = sp.solve([A, B], [s, lam ** 2], dict=True)
    return {'residual_equals_factor_form': bool(ok), 'solutions': sol}


def main():
    x = np.linspace(XMIN, XMAX, NX)
    t = np.linspace(0, TEND, NT)
    X, T = np.meshgrid(x, t)

    rows = []
    rows.append({'case': 'repo params (tau=0.2, lam=sqrt(3.8)/2, s=1.5)',
                 'tau': TAU, 'lam': LAM, 's': S,
                 'max_abs_residual': float(np.max(np.abs(pde_residual(TAU, LAM, S, X, T))))})
    for name, tau, lam, s in [
        ('control: s=2.0 (wrong speed)', TAU, LAM, 2.0),
        ('control: lam*1.1 (wrong width)', TAU, 1.1 * LAM, S),
        ('control: tau=0.3, lam unchanged', 0.3, LAM, S),
    ]:
        rows.append({'case': name, 'tau': tau, 'lam': lam, 's': s,
                     'max_abs_residual': float(np.max(np.abs(pde_residual(tau, lam, s, X, T))))})

    # initial-velocity consistency: v0 = s*f*(1-f) vs f_t(x,0)
    f0 = 1.0 / (1.0 + np.exp(LAM * x))
    ft0 = S * f0 * (1 - f0)  # f_t = -s*g' = +s*g(1-g)
    rows.append({'case': 'v0 = s*f0*(1-f0) vs f_t(x,0)', 'tau': '', 'lam': '', 's': '',
                 'max_abs_residual': float(np.max(np.abs(S * f0 * (1 - f0) - ft0)))})

    out = HERE / 'nizovtseva_pde_residual_check.csv'
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['case', 'tau', 'lam', 's', 'max_abs_residual'])
        w.writeheader()
        w.writerows(rows)

    sym = symbolic_factor()
    print(f'grid: x in [{XMIN},{XMAX}] ({NX} pts), t in [0,{TEND}] ({NT} pts)')
    for r in rows:
        print(f"  {r['case']}: max|res| = {r['max_abs_residual']:.3e}")
    if sym is not None:
        print('symbolic: residual == g(1-g)[A(1-g)+Bg] ?', sym['residual_equals_factor_form'])
        print('symbolic solutions (s, lam^2):', sym['solutions'])
    print(f'wrote {out.name}')
    print(f'wave speed c = s/lam = {S / LAM!r}')


if __name__ == '__main__':
    main()
