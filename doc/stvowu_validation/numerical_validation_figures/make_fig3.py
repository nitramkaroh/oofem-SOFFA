#!/usr/bin/env python3
"""
make_fig3.py — Boundary-time-history validation figure.

E_Gamma(t) = ||u_h(t) - u_ref(t)||_{L2(Gamma)}  as FE-RMS norm over the
boundary nodes (surface variable of the model is the trace u|_Gamma, so
this is the model-consistent boundary error norm).

Two series, all 80 time steps (t = k*2.5e-4, k = 1..80):
  * OOFEM n32x16 (T1-A rerun) vs DOP853 n128x64 restricted to n32 nodes
    (blue, solid)
  * semidiscrete family level-difference U_64 - U_128|_64 on the 64x32
    boundary (vermilion, dashed) — reference-solver context only

The endpoint values at T are cross-checked against the audited archive
before writing any output. Read-only w.r.t. pre-existing files; writes
only to this directory (fig3_*.pdf/png and
data/fig3_boundary_time_series.csv).

Run:  python3 make_fig3.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import figcommon as fc

# -------------------------------------------------------------- time loops --
print('[1] computing boundary time series (80 steps)')
u128_all = fc.ref('n128')['U']
u64_all = fc.ref('n64')['U']
m32 = fc.restrict_map(128, 64, 32, 16)
m64 = fc.restrict_map(128, 64, 64, 32)

t = fc.ref('n128')['t_eval'][1:]          # stored evaluation times
assert np.abs(t - np.arange(1, 81) * fc.DT).max() < 1e-15
E_oof = np.empty(80)
E_fam = np.empty(80)
for idx, k in enumerate(range(1, 81)):
    o = fc.oofem_field('n32', k)
    e_oof = o - u128_all[k][m32]
    E_oof[idx] = fc.fe_rms_surf(e_oof, 32, 16)
    e_fam = u64_all[k] - u128_all[k][m64]
    E_fam[idx] = fc.fe_rms_surf(e_fam, 64, 32)

print('[2] cross-check endpoints vs audited archive (t = T)')
arch = fc.archive_T_metrics()
def check(name, val, key, rtol=2e-6):
    a = arch[key]
    ok = abs(val - a) <= rtol * max(abs(a), 1e-30)
    print(f'  check {name:34s} {val:.6e}  (archive {a:.6e})  {"OK" if ok else "MISMATCH"}')
    assert ok, name
check('E_Gamma(T) OOFEM n32 vs proxy128', E_oof[-1], ('n32x16 surf_fe_rms @T', 'vs_proxy128'))
check('E_Gamma(T) family 64->128',        E_fam[-1], ('family 64->128 surf_fe_rms @T', 'family_level_diff'))

# ------------------------------------------------------------------- data --
import csv
from pathlib import Path
out = Path(__file__).parent
(out / 'data').mkdir(exist_ok=True)
with open(out / 'data' / 'fig3_boundary_time_series.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['t', 'E_gamma_oofem_n32_vs_proxy128', 'E_gamma_family_64_128_on_n64'])
    for k, (eo, ef) in enumerate(zip(E_oof, E_fam)):
        w.writerow([f'{t[k]:.6f}', f'{eo:.10e}', f'{ef:.10e}'])
print('[3] data written: data/fig3_boundary_time_series.csv')

# ------------------------------------------------------------------ plot --
print('[4] rendering')
plt.rcParams.update({
    'font.size': 8.5, 'axes.linewidth': 0.8, 'axes.titlesize': 9,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'figure.dpi': 150, 'savefig.dpi': 300,
})
C_OOF = '#0072B2'
C_FAM = '#D55E00'

fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.semilogy(t, E_oof, '-', color=C_OOF, lw=1.3,
            label=r'OOFEM n32$\times$16 vs DOP853 n128 (restricted)')
ax.semilogy(t, E_fam, '--', color=C_FAM, lw=1.2,
            label=r'family $U_{64}-U_{128}|_{64}$ (n64 boundary)')
ax.set_yscale('log')
ax.grid(True, which='both', lw=0.4, color='0.8', zorder=0)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.set_xlabel('$t$')
ax.set_ylabel(r'$E_{\Gamma}(t)=\|u_h-u_{\mathrm{ref}}\|_{\mathrm{FE\,RMS},\,\Gamma}$')
ax.legend(loc='upper right', frameon=True, framealpha=0.92,
          edgecolor='0.8', facecolor='white', fontsize=7.8)
ax.set_xlim(0.0, 0.021)
ax.set_ylim(8e-6, 4e-3)
fig.tight_layout()
fig.savefig(out / 'fig3_boundary_error_time_series.pdf')
fig.savefig(out / 'fig3_boundary_error_time_series.png')
print('[5] written: fig3_boundary_error_time_series.pdf / .png')
print('DONE')
