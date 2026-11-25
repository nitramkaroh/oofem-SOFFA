#!/usr/bin/env python3
"""
make_fig1.py — Main validation figure (3 panels) at t* = T = 0.02.

  (a) OOFEM solution u_h on the 32x16 mesh (T1-A rerun, step 80)
  (b) Independent semidiscrete DOP853 reference (n128x64, rtol 2e-12),
      step 80, restricted to the n32 nodes by EXACT index map (nested
      meshes coincide nodally; zero interpolation error)
  (c) |u_h - u_ref| on the n32 nodes, separate colour scale

Panels (a) and (b) share one colour scale; panel (c) has its own.
The reference is a semidiscrete (space-discretised) DOP853 solution of
the same finite-element system — NOT an exact PDE solution. The
same-mesh OOFEM-vs-DOP853 floor at T is 2.916141e-08 (bulk FE RMS),
7-8 orders of magnitude below the displayed error, so the displayed
error is the spatial discretisation error of OOFEM relative to the
semidiscrete family.

All plotted values are cross-checked against the audited archive
(t1b_metrics_T.csv, floors_recomputed.csv) before writing any output.
Read-only with respect to every pre-existing file; writes only to this
directory (fig1_*.pdf/png and data/fig1_fields_T0p02.csv).

Run:  python3 make_fig1.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import figcommon as fc

# ------------------------------------------------------------- assertions --
def check(name, val, ref, rtol=2e-6):
    ok = abs(val - ref) <= rtol * max(abs(ref), 1e-30)
    print(f'  check {name:32s} {val:.6e}  (archive {ref:.6e})  {"OK" if ok else "MISMATCH"}')
    assert ok, name

print('[1] loading fields (t* = T = 0.02, step 80)')
u_oof = fc.oofem_field('n32', 80)              # OOFEM, 32x16, T1-A rerun
u_ref_fine = fc.ref_at('n128', 80)             # DOP853 n128x64 semidiscrete
m128_32 = fc.restrict_map(128, 64, 32, 16)
u_ref = u_ref_fine[m128_32]                    # exact restriction to n32 nodes
diff = u_oof - u_ref                           # signed error (for norms)
err = np.abs(diff)                             # absolute error (for panel c)

arch = fc.archive_T_metrics()
print('[2] cross-checks vs audited archive')
check('bulk FE RMS',  fc.fe_rms_bulk(diff, 32, 16),
      arch[('n32x16 bulk_fe_rms @T', 'vs_proxy128')])
check('surf FE RMS',  fc.fe_rms_surf(diff, 32, 16),
      arch[('n32x16 surf_fe_rms @T', 'vs_proxy128')])
bm = fc.bnd_mask(32, 16)
check('max (interior nodes)', float(err[~bm].max()),
      arch[('n32x16 bulk_max @T', 'vs_proxy128')])
check('same-mesh floor n32@T', fc.fe_rms_bulk(u_oof - fc.ref_at('n32', 80), 32, 16),
      fc.archive_floor('n32x16', 0.02, 'bulk_fe_rms'))

# ------------------------------------------------------------------- data --
import csv
from pathlib import Path
X, Y = fc.grid(32, 16)
out = Path(__file__).parent
(out / 'data').mkdir(exist_ok=True)
with open(out / 'data' / 'fig1_fields_T0p02.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['x', 'y', 'u_oofem_n32', 'u_dop853_n128_restricted', 'abs_diff'])
    for j in range(17):
        for i in range(33):
            k = j * 33 + i
            w.writerow([f'{X[j, i]:.6f}', f'{Y[j, i]:.6f}',
                        f'{u_oof[k]:.10e}', f'{u_ref[k]:.10e}', f'{err[k]:.10e}'])
print('[3] data written: data/fig1_fields_T0p02.csv')

# ------------------------------------------------------------------ plot --
print('[4] rendering')
plt.rcParams.update({
    'font.size': 8.5, 'axes.linewidth': 0.8, 'axes.titlesize': 9,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 150, 'savefig.dpi': 300,
})
CMAP = plt.get_cmap('viridis')
vmin = min(u_oof.min(), u_ref.min())
vmax = max(u_oof.max(), u_ref.max())

fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7),
                         gridspec_kw={'width_ratios': [1, 1, 1.06]})

def field_panel(ax, C, cmap, vmin_, vmax_, cbar):
    p = ax.pcolormesh(X, Y, C.reshape(17, 33), cmap=cmap,
                      vmin=vmin_, vmax=vmax_, shading='gouraud')
    ax.set_aspect('equal')
    ax.set_title(cbar)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    return p

U = u_oof.reshape(17, 33)
R = u_ref.reshape(17, 33)
E = err.reshape(17, 33)

p1 = field_panel(axes[0], U, CMAP, vmin, vmax,
                 r'(a)  $u_h$ — OOFEM (32$\times$16)')
p2 = field_panel(axes[1], R, CMAP, vmin, vmax,
                 r'(b)  $u_{\mathrm{ref}}$ — DOP853 n128$\times$64, restricted')
p3 = field_panel(axes[2], E, CMAP, 0.0, E.max(),
                 r'(c)  $|u_h - u_{\mathrm{ref}}|$')

for ax in axes:
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

for p, ax in ((p1, axes[0]), (p2, axes[1]), (p3, axes[2])):
    cb = fig.colorbar(p, ax=ax, fraction=0.042, pad=0.05)
    cb.ax.tick_params(labelsize=7.5)

fig.text(0.005, 0.012,
         r'$t^{*}=T=0.02$  (step 80);  panels (a), (b) share the scale;  '
         r'(c) is the pointwise absolute difference on the 32$\times$16 nodes.',
         fontsize=7.5, color='0.3')

fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(out / 'fig1_oofem_vs_dop853_fields_T.pdf')
fig.savefig(out / 'fig1_oofem_vs_dop853_fields_T.png')
print('[5] written: fig1_oofem_vs_dop853_fields_T.pdf / .png')
print('DONE')
