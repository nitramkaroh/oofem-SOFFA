#!/usr/bin/env python3
"""
make_lit_fig2.py — OOFEM on the Hochbruck--Leibold exact profile: spatial
refinement on the polygonally approximated unit disk.

Same adapted exact-profile test as Fig. 1, now with fixed dt = 0.0025 and
gmsh-triangulated polygonal unit disks (60 / 179 / 502 nodes,
h_target 0.3 / 0.16 / 0.09).

Panel (a): RMS error vs mesh size with the measured local log-2 ratios
(0.33, 0.24) annotated — no order is claimed.
Panel (b): maximum pointwise error, essentially flat (approx 4.9e-2),
indicating a geometric (boundary-approximation) error floor.

The reduced spatial rate is the boundary-approximation effect that the
bulk-surface literature emphasizes for non-polygonal domains
(Kovacs & Lubich 2017, Introduction, 3rd remark; Hipp & Kovacs 2021,
Sec. 2, nonconforming bulk-surface FEM).  See README.

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_lit_fig2.py
"""
import csv
from pathlib import Path

import figcommon_lit as fc
import matplotlib.pyplot as plt

print('[1] recomputing final-field errors from the archived VTU files')
study = fc.hl_space_study()   # [(h, nodes, rms, max)]
(h1, n1, e1, x1), (h2, n2, e2, x2), (h3, n3, e3, x3) = study
o12 = fc.log2(e1, e2)
o23 = fc.log2(e2, e3)
print(f'  local spatial ratios: {o12:.4f}, {o23:.4f}')

print('[2] writing data CSV')
out = Path(__file__).parent
with open(out / 'data_fig2_hl_profile_spatial.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['h_target', 'nodes', 'rms_error_vs_exact_profile',
                'max_error_vs_exact_profile', 'dt', 'T'])
    for h, n, e, x in study:
        w.writerow([f'{h:g}', n, f'{e:.10e}', f'{x:.10e}', 0.0025, fc.T])

print('[3] rendering')
plt = fc.style_fig()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.6))

h = [h1, h2, h3]
ax1.loglog(h, [e1, e2, e3], 'o-', color=fc.C_OOF, ms=5, lw=1.2, zorder=3)
fc.annotate_order(ax1, h1, h2, e1, e2, f'local $\\approx$ {o12:.2f}')
fc.annotate_order(ax1, h2, h3, e2, e3, f'local $\\approx$ {o23:.2f}')
ax1.set_xlabel(r'mesh size $h$ (refinement to the right)')
ax1.set_ylabel(r'RMS error vs exact profile at $T$')
ax1.set_title('(a) spatial refinement, polygonal disk')
ax1.invert_xaxis()
ax1.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
ax1.set_xlim(h1 * 1.25, h3 * 0.8)

ax2.plot(h, [x1, x2, x3], 'o-', color=fc.C_OOF, ms=5, lw=1.2, zorder=3)
ax2.annotate(f'$\\approx$ 4.9e-2 (variation $<$ 1%)', xy=(h2, x2),
             xytext=(0, 8), textcoords='offset points', ha='center',
             fontsize=7.5, color='0.2', bbox=fc.ABOX)
ax2.set_xlabel(r'mesh size $h$ (refinement to the right)')
ax2.set_ylabel('maximum pointwise error')
ax2.set_title('(b) geometric floor')
ax2.invert_xaxis()
ax2.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
ax2.set_xlim(h1 * 1.25, h3 * 0.8)
ax2.set_ylim(0.045, 0.050)

fig.suptitle('OOFEM vs the Hochbruck--Leibold exact profile — space refinement '
             '(adapted operator, polygonal unit disk)', fontsize=9.5, y=1.02)
fig.tight_layout()
fig.savefig(out / 'lit_fig2_hl_profile_spatial.pdf', bbox_inches='tight')
fig.savefig(out / 'lit_fig2_hl_profile_spatial.png', bbox_inches='tight')
print('[4] written: lit_fig2_hl_profile_spatial.pdf / .png')
print('DONE')
