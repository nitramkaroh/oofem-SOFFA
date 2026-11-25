#!/usr/bin/env python3
"""
make_cmp_fig1.py - OOFEM solution vs the literature (Hochbruck--Leibold 2021)
exact profile, side by side on the polygonal unit disk.

Panel (a): OOFEM solution field at T (finest archived run, 502 nodes,
h ~ 0.09, dt = 0.0025) as a filled triangulated plot, with the published
profile u = sin(2 pi T) x y drawn as black contour lines.
Panel (b): pointwise difference u_h - u_exact on the same triangulation
(symmetric colour scale) - where the two solutions agree, the panel is
uniform; the deviation is concentrated at the polygonal boundary.

The forcing of the archived run was manufactured for OOFEM's own operator
(adapted exact-profile test; see README).

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_cmp_fig1.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

import figcommon_lit as fc

print('[1] reading the finest archived HL-profile run (space_h0.09_dt0p0025)')
d = 'hl_disk_space_h0.09_dt0p0025'
coords, u, tri = fc.final_tri(fc.HYPER / d)
rms, mx = fc.rms_max(coords, u, fc.hl_exact)
rows = {r['tag']: r for r in fc._csv_rows('hochbruck_leibold_disk_metrics.csv')}
fc._check(f'{d} rms', rms,
          float(rows['space_h0.09_dt0p0025']['rms_error_vs_exact_profile']))
fc._check(f'{d} max', mx,
          float(rows['space_h0.09_dt0p0025']['max_error_vs_exact_profile']))

exact = fc.hl_exact(coords)
diff = u - exact
print(f'  max |diff| = {mx:.4e}   RMS = {rms:.4e}')

out = Path(__file__).parent
plt = fc.style_fig()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.9))

T = Triangulation(coords[:, 0], coords[:, 1], tri)
levels = np.arange(-0.2, 0.2001, 0.05)

pc = ax1.tripcolor(T, u, cmap='viridis', shading='gouraud')
fig.colorbar(pc, ax=ax1, fraction=0.046, pad=0.04)
ax1.tricontour(T, exact, levels=levels, colors='k', linewidths=0.8)
ax1.plot(coords[:, 0], coords[:, 1], '.', ms=1, color='0.4')
ax1.set_title('(a) OOFEM field (shaded) vs literature profile (contours)')
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
ax1.set_aspect('equal')

pd = ax2.tripcolor(T, diff, cmap='RdBu_r', shading='gouraud',
                   vmin=-mx, vmax=mx)
fig.colorbar(pd, ax=ax2, fraction=0.046, pad=0.04)
ax2.tricontour(T, exact, levels=levels, colors='0.5', linewidths=0.6, linestyles=':')
ax2.set_title(f'(b) difference $u_{{h}}-u_{{\\mathrm{{exact}}}}$, '
              f'max error = {mx:.2e}')
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_aspect('equal')

fig.suptitle('OOFEM vs the Hochbruck--Leibold (2021) exact profile at '
             r'$T=0.2$ (adapted operator, unit disk)', fontsize=9.5, y=1.02)
fig.tight_layout()
fig.savefig(out / 'cmp_fig1_hl_profile_field.pdf', bbox_inches='tight')
fig.savefig(out / 'cmp_fig1_hl_profile_field.png', bbox_inches='tight')
print('[2] written: cmp_fig1_hl_profile_field.pdf / .png')
print('DONE')
