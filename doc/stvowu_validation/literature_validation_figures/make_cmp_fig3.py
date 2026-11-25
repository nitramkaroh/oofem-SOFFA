#!/usr/bin/env python3
"""
make_cmp_fig3.py — OOFEM solution vs the literature stationary front
tanh(x/sqrt(2)) of the hyperbolic Allen-Cahn equation, side by side.

Panel (a): exact profile u(x) = tanh(x/sqrt(2)) as a black dashed curve;
the three archived OOFEM solutions (Nx = 40/80/160, dt = 0.001, T = 0.2,
front on [-8,8]) as blue curves — coinciding with the literature profile
within the plot resolution.
Panel (b): pointwise differences u_h - u_exact for the three meshes.

All final fields are re-read from the archived VTU files and asserted
against the archived metrics before plotting.

Read-only w.r.t. pre-existing files; writes only to this directory.
Run:  python3 make_cmp_fig3.py
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import figcommon_lit as fc

print('[1] reading the three archived kink runs (final fields)')
runs = [('kink_n40', '40', 40), ('kink_n80', '80', 80), ('kink_n160', '160', 160)]
curves = []
for d, key, nx in runs:
    coords, u = fc.check_field(d, fc.kink_exact, 'stationary_kink_space_refinement.csv',
                               'nx', key, 'rms_error', 'max_error')
    sel = coords[:, 1] < 0.5          # lower row (field is 1-D in x)
    order = np.argsort(coords[sel, 0])
    curves.append((coords[sel, 0][order], u[sel][order], nx))

xg = np.linspace(-8.0, 8.0, 2001)
exact = np.tanh(xg / np.sqrt(2.0))
maxdiff = {nx: float(np.max(np.abs(uu - fc.kink_exact(np.column_stack([x, np.zeros_like(x)])))))
           for x, uu, nx in curves}
print('  max |diff| per mesh:', {k: f'{v:.3e}' for k, v in maxdiff.items()})

out = Path(__file__).parent
with open(out / 'data_cmp_fig3_kink_profile.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['x', 'exact', 'oofem_nx40', 'oofem_nx80', 'oofem_nx160'])
    for i, xi in enumerate(xg[::10]):
        row = [f'{xi:.8e}', f'{exact[i * 10]:.10e}']
        for (x, uu, nx) in curves:
            j = np.searchsorted(x, xi)
            row.append(f'{uu[min(j, len(uu) - 1)]:.10e}')
        w.writerow(row)

plt = fc.style_fig()
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 6.0), sharex=True)

alphas = (0.35, 0.60, 1.00)
widths = (1.0, 1.2, 1.5)
for (x, uu, nx), a, lw in zip(curves, alphas, widths):
    lab = f'OOFEM, $N_x$ = {nx}' if nx == 160 else None
    ax1.plot(x, uu, color=fc.C_OOF, alpha=a, lw=lw, label=lab)
    ax2.plot(x, uu - fc.kink_exact(np.column_stack([x, np.zeros_like(x)])),
             color=fc.C_OOF, alpha=a, lw=lw)
ax1.plot(xg, exact, 'k--', lw=1.3, label=r'exact literature profile $\tanh(x/\sqrt{2})$')
ax1.set_ylabel(r'$u$ at $T = 0.2$')
ax1.set_title('(a) OOFEM vs the stationary hyperbolic Allen--Cahn front')
ax1.legend(loc='lower center', frameon=False, fontsize=7.5)
ax1.grid(True, lw=0.4, color='0.8')

ax2.axhline(0.0, color='0.6', lw=0.6)
ax2.set_ylabel(r'$u_h - u_{\mathrm{exact}}$')
ax2.set_xlabel(r'$x$  (domain $[-8,8]$, $\Delta t = 10^{-3}$ s)')
mc, mf = maxdiff[40], maxdiff[160]
ax2.set_title(f'(b) pointwise difference, max error '
               rf'{mc:.1e} (coarsest) $\rightarrow$ {mf:.1e} (finest)')
ax2.grid(True, lw=0.4, color='0.8')
ax2.set_xlim(-8.5, 8.5)

fig.suptitle('OOFEM vs the literature stationary front — hyperbolic '
             r'Allen--Cahn bulk, $u_{tt}-\Delta u = u-u^3$', fontsize=9.5, y=1.00)
fig.tight_layout()
fig.savefig(out / 'cmp_fig3_kink_profile.pdf', bbox_inches='tight')
fig.savefig(out / 'cmp_fig3_kink_profile.png', bbox_inches='tight')
print('[2] written: cmp_fig3_kink_profile.pdf / .png')
print('DONE')
