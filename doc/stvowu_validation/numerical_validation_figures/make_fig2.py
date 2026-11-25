#!/usr/bin/env python3
"""
make_fig2.py — Spatial error vs mesh size, FE RMS norms, at t* = T = 0.02.

Two panels (bulk domain / boundary Gamma), log-log.

Series (both at the SAME time T, same physical problem):
  * OOFEM chain (blue, circles, solid):
        E_h = ||u_h - u_{ref,128}|_h||_FE-RMS   for h = 0.25, 0.125, 0.0625
    where u_{ref,128} is the DOP853 semidiscrete n128x64 solution restricted
    to the coarser nodes (exact index map).
  * Semidiscrete family level-differences (vermilion, squares, dashed):
        h = 0.0625 : ||U_32 - U_64|_32||  (evaluated on the 32x16 mesh)
        h = 0.03125: ||U_64 - U_128|_64|| (evaluated on the 64x32 mesh)
    h is the mesh size of the evaluation mesh (convention verified by the
    2026-10-27 audit).

Local log2 ratios of CONSECUTIVE points are annotated as measured values
(0.72, 1.05 for the OOFEM chain; 3.09 for the family at T). No theoretical
slope is drawn: the audit found 0/80 time steps with both consecutive
family ratios in [1.5, 2.5], so no asymptotic order is claimed.

Read-only w.r.t. pre-existing files; writes only to this directory
(fig2_*.pdf/png and data/fig2_convergence_data.csv).

Run:  python3 make_fig2.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import figcommon as fc

# ------------------------------------------------------------- assertions --
print('[1] recomputing E(T) and cross-checking vs audited archive')
u128 = fc.ref_at('n128', 80)
u64 = fc.ref_at('n64', 80)
u32 = fc.ref_at('n32', 80)
m8, m16, m32 = (fc.restrict_map(128, 64, nx, ny) for nx, ny in ((8, 4), (16, 8), (32, 16)))
m64 = fc.restrict_map(128, 64, 64, 32)
i32_64 = fc.restrict_map(64, 32, 32, 16)

arch = fc.archive_T_metrics()
def check(name, val, key, rtol=2e-6):
    a = arch[key]
    ok = abs(val - a) <= rtol * max(abs(a), 1e-30)
    print(f'  check {name:30s} {val:.6e}  (archive {a:.6e})  {"OK" if ok else "MISMATCH"}')
    assert ok, name

o8, o16, o32 = fc.oofem_field('n8', 80), fc.oofem_field('n16', 80), fc.oofem_field('n32', 80)
eo8 = o8 - u128[m8]; eo16 = o16 - u128[m16]; eo32 = o32 - u128[m32]
ef3264 = u32 - u64[i32_64]; ef64128 = u64 - u128[m64]

check('oofem n8  bulk', fc.fe_rms_bulk(eo8, 8, 4), ('n8x4 bulk_fe_rms @T', 'vs_proxy128'))
check('oofem n16 bulk', fc.fe_rms_bulk(eo16, 16, 8), ('n16x8 bulk_fe_rms @T', 'vs_proxy128'))
check('oofem n32 bulk', fc.fe_rms_bulk(eo32, 32, 16), ('n32x16 bulk_fe_rms @T', 'vs_proxy128'))
check('oofem n8  surf', fc.fe_rms_surf(eo8, 8, 4), ('n8x4 surf_fe_rms @T', 'vs_proxy128'))
check('oofem n16 surf', fc.fe_rms_surf(eo16, 16, 8), ('n16x8 surf_fe_rms @T', 'vs_proxy128'))
check('oofem n32 surf', fc.fe_rms_surf(eo32, 32, 16), ('n32x16 surf_fe_rms @T', 'vs_proxy128'))
check('fam 32->64 bulk', fc.fe_rms_bulk(ef3264, 32, 16), ('family 32->64 bulk_fe_rms @T', 'family_level_diff'))
check('fam 64->128 bulk', fc.fe_rms_bulk(ef64128, 64, 32), ('family 64->128 bulk_fe_rms @T', 'family_level_diff'))
check('fam 32->64 surf', fc.fe_rms_surf(ef3264, 32, 16), ('family 32->64 surf_fe_rms @T', 'family_level_diff'))
check('fam 64->128 surf', fc.fe_rms_surf(ef64128, 64, 32), ('family 64->128 surf_fe_rms @T', 'family_level_diff'))

# ------------------------------------------------------------- data table --
import csv
from pathlib import Path
out = Path(__file__).parent
(out / 'data').mkdir(exist_ok=True)
h_oof = [0.25, 0.125, 0.0625]
E_bulk_oof = [fc.fe_rms_bulk(e, nx, ny) for e, (nx, ny) in
              zip((eo8, eo16, eo32), ((8, 4), (16, 8), (32, 16)))]
E_surf_oof = [fc.fe_rms_surf(e, nx, ny) for e, (nx, ny) in
              zip((eo8, eo16, eo32), ((8, 4), (16, 8), (32, 16)))]
h_fam = [0.0625, 0.03125]
E_bulk_fam = [fc.fe_rms_bulk(ef3264, 32, 16), fc.fe_rms_bulk(ef64128, 64, 32)]
E_surf_fam = [fc.fe_rms_surf(ef3264, 32, 16), fc.fe_rms_surf(ef64128, 64, 32)]

with open(out / 'data' / 'fig2_convergence_data.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['series', 'norm', 'h', 'E_FE_RMS'])
    for h, eb, es in zip(h_oof, E_bulk_oof, E_surf_oof):
        w.writerow(['oofem_vs_dop853_n128', 'bulk', f'{h:g}', f'{eb:.10e}'])
        w.writerow(['oofem_vs_dop853_n128', 'surf', f'{h:g}', f'{es:.10e}'])
    for h, eb, es in zip(h_fam, E_bulk_fam, E_surf_fam):
        w.writerow(['semidiscrete_family_level_diff', 'bulk', f'{h:g}', f'{eb:.10e}'])
        w.writerow(['semidiscrete_family_level_diff', 'surf', f'{h:g}', f'{es:.10e}'])
print('[2] data written: data/fig2_convergence_data.csv')

# ------------------------------------------------------------------ plot --
print('[3] rendering')
plt.rcParams.update({
    'font.size': 8.5, 'axes.linewidth': 0.8, 'axes.titlesize': 9,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 150, 'savefig.dpi': 300,
})
C_OOF = '#0072B2'   # Okabe-Ito blue      — OOFEM
C_FAM = '#D55E00'   # Okabe-Ito vermilion — semidiscrete family
# CVD separation of this pair: Oklab Delta-E x100 = 32.7 (>= 8 required)

def log2(a, b):
    return float(np.log2(a / b))

ABOX = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85)
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.9))
for ax, (Eb, Ef) in zip(axes, ((E_bulk_oof, E_bulk_fam), (E_surf_oof, E_surf_fam))):
    ax.loglog(h_oof, Eb, 'o-', color=C_OOF, ms=5, lw=1.2, zorder=3,
              label=r'$u_h$: OOFEM vs DOP853 n128 (restricted)')
    ax.loglog(h_fam, Ef, 's--', color=C_FAM, ms=5, lw=1.2, zorder=3,
              label=r'family level-diff $U_{32}-U_{64}$, $U_{64}-U_{128}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.invert_xaxis()
    ax.set_xticks(h_oof + [0.03125])
    ax.set_xticklabels([f'{h:g}' for h in h_oof + [0.03125]])
    ax.minorticks_on()
    ax.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
    ax.set_xlabel('mesh size $h$ (refinement to the right)')
    # local log2 ratios of consecutive points, as measured at T
    r1, r2 = log2(Eb[0], Eb[1]), log2(Eb[1], Eb[2])
    ax.annotate(f'local order $\\approx$ {r1:.2f}',
                xy=((h_oof[0] + h_oof[1]) / 2, np.sqrt(Eb[0] * Eb[1])),
                ha='center', va='bottom', fontsize=7.5, color='0.2', bbox=ABOX)
    ax.annotate(f'local order $\\approx$ {r2:.2f}',
                xy=((h_oof[1] + h_oof[2]) / 2, np.sqrt(Eb[1] * Eb[2])),
                ha='center', va='bottom', fontsize=7.5, color='0.2', bbox=ABOX)
    rf = log2(Ef[0], Ef[1])
    ax.annotate(f'local ratio $\\approx$ {rf:.2f} (single, at $T$)',
                xy=((h_fam[0] + h_fam[1]) / 2, np.sqrt(Ef[0] * Ef[1])),
                xytext=(0, -18), textcoords='offset points',
                ha='center', va='top', fontsize=7.5, color='0.2', bbox=ABOX)
    ax.set_xlim(0.29, 0.027)
    ax.set_title('bulk $\\Omega$' if ax is axes[0] else 'boundary $\\Gamma$')
    ax.set_ylabel(r'$E(h)=\|u_h-u_{\mathrm{ref}}\|_{\mathrm{FE\,RMS}}$  at $t^{*}=T$')

# single shared legend (both panels show the same two series, same styles);
# series carry label= so the auto legend gets real handles
axes[0].legend(loc='lower left', frameon=False, fontsize=7.5)
fig.tight_layout()
fig.savefig(out / 'fig2_convergence_FE_RMS_T.pdf', bbox_inches='tight')
fig.savefig(out / 'fig2_convergence_FE_RMS_T.png', bbox_inches='tight')
print('[4] written: fig2_convergence_FE_RMS_T.pdf / .png')
print('DONE')
