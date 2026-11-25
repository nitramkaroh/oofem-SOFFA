#!/usr/bin/env python3
"""
make_cmp_fig2.py — OOFEM solution vs the literature exact traveling wave
(Nizovtseva, Galenko & Alexandrov 2016, Eq. (29)), side by side.

Panel (a): the exact profile u(x, T) = 1/(1+exp(lambda x - 1.5 T))
(lambda = sqrt(2+9 tau)/2, tau = 0.2) as a black dashed curve; the three
archived OOFEM solutions (nx = 80/160/320, dt = 0.001, T = 0.2) as blue
curves — they coincide with the literature solution within the plot
resolution.
Panel (b): pointwise differences u_h - u_exact for the three meshes.

All final fields are re-read from the archived VTU files and asserted
against the archived metrics before plotting.

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_cmp_fig2.py
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import figcommon_lit as fc

print('[1] reading the three archived Nizovtseva runs (final fields)')
runs = [('nizovtseva_space_n80_dt0p001', 'space_n80_dt0p001', 80),
        ('nizovtseva_space_n160_dt0p001', 'space_n160_dt0p001', 160),
        ('nizovtseva_space_n320_dt0p001', 'space_n320_dt0p001', 320)]
curves = []
for d, tag, nx in runs:
    coords, u = fc.check_field(d, fc.niz_exact, 'nizovtseva_traveling_wave_metrics.csv',
                               'tag', tag, 'rms_error_vs_exact', 'max_error_vs_exact')
    sel = coords[:, 1] < 0.5          # lower row of the strip (field is 1-D)
    x = np.sort(coords[sel, 0])
    uu = u[sel][np.argsort(coords[sel, 0])]
    curves.append((x, uu, nx))

xg = np.linspace(-12.0, 12.0, 2001)
exact = fc.niz_exact(np.column_stack([xg, np.zeros_like(xg)]))
maxdiff = {nx: float(np.max(np.abs(uu - fc.niz_exact(np.column_stack([x, np.zeros_like(x)])))))
           for x, uu, nx in curves}
print('  max |diff| per mesh:', {k: f'{v:.3e}' for k, v in maxdiff.items()})

out = Path(__file__).parent
with open(out / 'data_cmp_fig2_nizovtseva_wave.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['x', 'exact', 'oofem_n80', 'oofem_n160', 'oofem_n320'])
    for i, xi in enumerate(xg[::10]):
        row = [f'{xi:.8e}']
        for (x, uu, nx) in curves:
            j = np.searchsorted(x, xi)
            row.append(f'{uu[min(j, len(uu) - 1)]:.10e}')
        row.insert(1, f'{exact[i * 10]:.10e}')
        w.writerow(row)

plt = fc.style_fig()
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 6.0), sharex=True)

alphas = (0.35, 0.60, 1.00)
widths = (1.0, 1.2, 1.5)
for (x, uu, nx), a, lw in zip(curves, alphas, widths):
    lab = f'OOFEM, $n_x$ = {nx}' if nx == 320 else None
    ax1.plot(x, uu, color=fc.C_OOF, alpha=a, lw=lw, label=lab)
    ax2.plot(x, uu - fc.niz_exact(np.column_stack([x, np.zeros_like(x)])),
             color=fc.C_OOF, alpha=a, lw=lw)
ax1.plot(xg, exact, 'k--', lw=1.3, label=r'exact literature solution $u(x,T)$')
ax1.set_ylabel(r'$u$ at $T = 0.2$')
ax1.set_title('(a) OOFEM vs Nizovtseva et al. (2016) exact traveling wave')
ax1.legend(loc='upper left', frameon=False, fontsize=7.5)
ax1.grid(True, lw=0.4, color='0.8')

ax2.axhline(0.0, color='0.6', lw=0.6)
ax2.set_ylabel(r'$u_h - u_{\mathrm{exact}}$')
ax2.set_xlabel(r'$x$  (strip $[-12,12]$, $\Delta t = 10^{-3}$ s)')
mc, mf = maxdiff[80], maxdiff[320]
ax2.set_title(f'(b) pointwise difference, max error '
               rf'{mc:.1e} (coarsest) $\rightarrow$ {mf:.1e} (finest)')
ax2.grid(True, lw=0.4, color='0.8')
ax2.set_xlim(-12.5, 12.5)

fig.suptitle('OOFEM vs the literature exact traveling wave — hyperbolic '
             r'Allen--Cahn bulk, $\tau u_{tt}+u_t-\Delta u = u-u^3$, $\tau=0.2$',
             fontsize=9.5, y=1.00)
fig.tight_layout()
fig.savefig(out / 'cmp_fig2_nizovtseva_wave.pdf', bbox_inches='tight')
fig.savefig(out / 'cmp_fig2_nizovtseva_wave.png', bbox_inches='tight')
print('[2] written: cmp_fig2_nizovtseva_wave.pdf / .png')
print('DONE')
