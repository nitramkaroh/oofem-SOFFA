#!/usr/bin/env python3
"""
make_lit_fig1.py — OOFEM on the Hochbruck--Leibold exact profile: temporal refinement.

Published profile u(t,x) = sin(2 pi t) x y with initial velocity 2 pi x y
(Hochbruck & Leibold 2021, Sec. 4.5, their Fig. 1 test).  The forcing in the
archived OOFEM runs was manufactured for OOFEM's own operator (bulk-inertial
hyperbolic Allen-Cahn + first-order dynamic surface law), i.e. this is an
ADAPTED exact-profile test of OOFEM's time integration, not a reproduction
of the paper's second-order acoustic boundary equation (see README).

Panel (a): RMS error vs Delta t (fixed mesh, h ~ 0.16, 179 nodes) — the
error is pinned by the spatial floor of the fixed mesh.
Panel (b): difference between two successive refinements vs Delta t —
local order 2.02 measured from the two archived differences.

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_lit_fig1.py
"""
import csv
from pathlib import Path

import figcommon_lit as fc
import matplotlib.pyplot as plt

print('[1] recomputing final-field errors from the archived VTU files')
study = fc.hl_time_study()
# (dt, rms, max, diff) with diff of the last entry None
(d1, e1, m1, f1), (d2, e2, m2, f2), (d3, e3, m3, f3) = study
floor = e3
order = fc.log2(f1, f2)
print(f'  temporal local order (from successive differences): {order:.4f}')
assert 1.9 < order < 2.1, order

print('[2] writing data CSV')
out = Path(__file__).parent
with open(out / 'data_fig1_hl_profile_temporal.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['dt', 'rms_error_vs_exact_profile', 'max_error_vs_exact_profile',
                'diff_to_next_finer', 'T', 'nodes'])
    for (dt, e, m, fd) in study:
        w.writerow([f'{dt:g}', f'{e:.10e}', f'{m:.10e}',
                    f'{fd:.10e}' if fd is not None else '', fc.T, 179])

print('[3] rendering')
plt = fc.style_fig()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.6))

dt = [d1, d2, d3]
err = [e1, e2, e3]
ax1.loglog(dt, err, 'o-', color=fc.C_OOF, ms=5, lw=1.2, zorder=3)
ax1.axhline(floor, color=fc.C_FLOOR, ls=':', lw=1.0, zorder=2)
ax1.annotate('spatial floor (fixed $h \\approx$ 0.16)', xy=(d1, floor),
             xytext=(-8, 8), textcoords='offset points', ha='right',
             va='bottom', fontsize=7.5, color='0.3')
ax1.set_xlabel(r'$\Delta t$ (s)')
ax1.set_ylabel(r'RMS error vs exact profile at $T$')
ax1.set_title('(a) error pinned by the spatial floor')
ax1.invert_xaxis()
ax1.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
ax1.set_xlim(d3 * 0.85, d1 * 1.15)

df = [d1, d2]
ax2.loglog(df, [f1, f2], 's-', color=fc.C_OOF, ms=5, lw=1.2, zorder=3)
fc.annotate_order(ax2, d1, d2, f1, f2, f'local order $\\approx$ {order:.2f}')
ax2.set_xlabel(r'$\Delta t$ (s)')
ax2.set_ylabel(r'$\|u_{\Delta t}-u_{\Delta t/2}\|_{\mathrm{RMS}}$')
ax2.set_title('(b) successive-refinement differences')
ax2.invert_xaxis()
ax2.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
ax2.set_xlim(d2 * 0.85, d1 * 1.15)

fig.suptitle('OOFEM vs the Hochbruck--Leibold exact profile — time refinement '
             '(adapted operator, unit disk)', fontsize=9.5, y=1.02)
fig.tight_layout()
fig.savefig(out / 'lit_fig1_hl_profile_temporal.pdf', bbox_inches='tight')
fig.savefig(out / 'lit_fig1_hl_profile_temporal.png', bbox_inches='tight')
print('[4] written: lit_fig1_hl_profile_temporal.pdf / .png')
print('DONE')
