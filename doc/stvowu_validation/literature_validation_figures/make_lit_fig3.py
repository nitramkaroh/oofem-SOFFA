#!/usr/bin/env python3
"""
make_lit_fig3.py — OOFEM code verification against closed-form solutions of
the hyperbolic Allen-Cahn bulk dynamics (2x2 panel).

(a) Nizovtseva--Galenko--Alexandrov (2016) exact traveling wave,
    Eq. (29) (Q=A=B=+1, c0=0), tau u_tt + u_t - Laplacian u = u - u^3,
    tau = 0.2, strip [-12,12]x[0,1], T = 0.2, dt = 0.001, h refinement.
    Measured local orders 1.995 / 1.976.
(b) Stationary hyperbolic Allen-Cahn front tanh(x/sqrt(2)),
    [-8,8], T = 0.2, dt = 0.001, h refinement. Measured 2.006 / 1.887.
(c) Linear cosine mode, Neumann strip: error vs the exact SEMIDISCRETE Q1
    (consistent-mass) solution — a time-stepping test, dt refinement.
    Measured 1.997 (archived).
(d) Spatially uniform Allen-Cahn oscillator u_tt = u - u^3: error vs an
    independent ODE integration, dt refinement. Measured 2.016 (archived);
    the third step is close to the VTU output-precision floor.

Panels (a,b) are re-derived from the raw VTU files and asserted against the
archived metrics; (c,d) use the archived bulk_reference_metrics.csv as the
source of record (see README).

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_lit_fig3.py
"""
import csv
from pathlib import Path

import figcommon_lit as fc
import matplotlib.pyplot as plt

print('[1] recomputing (a,b) from the archived VTU files; (c,d) from CSV')
niz = fc.niz_space_study()
kink = fc.kink_space_study()
bulk = fc.bulk_reference_rows()

(h1, e1, _), (h2, e2, _), (h3, e3, _) = niz
(k1, g1, _), (k2, g2, _), (k3, g3, _) = kink
lin = bulk['linear']      # [(dt, rms, order_or_None)] dt descending
osc = bulk['oscillator']
o_niz = (fc.log2(e1, e2), fc.log2(e2, e3))
o_kink = (fc.log2(g1, g2), fc.log2(g2, g3))
o_lin = lin[1][2]
o_osc = osc[1][2]
print(f'  orders  niz: {o_niz[0]:.4f}/{o_niz[1]:.4f}  kink: {o_kink[0]:.4f}/{o_kink[1]:.4f}'
      f'  linear: {o_lin:.4f}  oscillator: {o_osc:.4f}')

print('[2] writing data CSV')
out = Path(__file__).parent
with open(out / 'data_fig3_bulk_exact_solutions.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['panel', 'test', 'resolution', 'rms_error'])
    for h, e, _ in niz:
        w.writerow(['a', 'nizovtseva_traveling_wave_space', f'{h:g}', f'{e:.10e}'])
    for h, e, _ in kink:
        w.writerow(['b', 'stationary_kink_space', f'{h:g}', f'{e:.10e}'])
    for dt, e, _ in lin:
        w.writerow(['c', 'linear_cosine_mode_time', f'{dt:g}', f'{e:.10e}'])
    for dt, e, _ in osc:
        w.writerow(['d', 'allencahn_oscillator_time', f'{dt:g}', f'{e:.10e}'])

print('[3] rendering')
plt = fc.style_fig()
fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2))
ax = axes.ravel()

def panel(a, pts, xlabel, title, pair_orders):
    """pts: [(x, y)]; pair_orders: list of order or None per consecutive pair."""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    a.loglog(xs, ys, 'o-', color=fc.C_OOF, ms=5, lw=1.2, zorder=3)
    for (x1, y1), (x2, y2), o in zip(pts[:-1], pts[1:], pair_orders):
        if o is not None:
            fc.annotate_order(a, x1, x2, y1, y2, f'$\\approx$ {o:.2f}')
    a.set_xlabel(xlabel)
    a.set_title(title, fontsize=8.5)
    a.invert_xaxis()
    a.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
    a.set_ylabel(r'error (RMS) at $T$')

panel(ax[0], [(h1, e1), (h2, e2), (h3, e3)],
      r'mesh size $h$', '(a) Nizovtseva et al. (2016) traveling wave — space',
      list(o_niz))
panel(ax[1], [(k1, g1), (k2, g2), (k3, g3)],
      r'mesh size $h$', '(b) stationary front tanh($x/\\sqrt{2}$) — space',
      list(o_kink))
# (c),(d): only the first pair is annotated — the second pair of the
# oscillator reaches the VTU output-precision floor (see caption)
panel(ax[2], [(r[0], r[1]) for r in lin],
      r'$\Delta t$ (s)', '(c) linear cosine mode — time (vs exact semidiscrete)',
      [o_lin, None])
panel(ax[3], [(r[0], r[1]) for r in osc],
      r'$\Delta t$ (s)', '(d) Allen--Cahn oscillator — time (vs independent ODE)',
      [o_osc, None])

fig.suptitle('OOFEM code verification — closed-form solutions of the '
             'hyperbolic Allen--Cahn bulk', fontsize=9.5, y=1.005)
fig.tight_layout()
fig.savefig(out / 'lit_fig3_bulk_exact_solutions.pdf', bbox_inches='tight')
fig.savefig(out / 'lit_fig3_bulk_exact_solutions.png', bbox_inches='tight')
print('[4] written: lit_fig3_bulk_exact_solutions.pdf / .png')
print('DONE')
