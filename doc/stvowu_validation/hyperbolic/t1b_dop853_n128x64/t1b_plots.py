#!/usr/bin/env python3
"""T1-B figures: regenerate the 3 PNGs in results/ from the result CSVs
(deterministic; metrics themselves come from t1b_validate.py)."""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

C1, C2 = '#2a78d6', '#eb6834'          # dataviz default slots 1-2
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'text.color': INK, 'axes.edgecolor': '#c9c8c4',
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': '#fcfcfb', 'figure.facecolor': '#fcfcfb',
                     'grid.color': '#c9c8c4', 'grid.alpha': 0.55, 'font.size': 9})


def style_ax(ax):
    ax.grid(True, lw=0.6)


def metrics():
    m = {}
    for row in csv.DictReader(open(RES / 't1b_metrics_T.csv')):
        m[row['quantity']] = row
    return m


def orders():
    return {r['order']: float(r['value'])
            for r in csv.DictReader(open(RES / 't1b_family_orders_T.csv'))}


def fam_series():
    return list(csv.DictReader(open(RES / 't1b_orders_time_64_128.csv')))


def moving_max():
    rows = [r for r in csv.DictReader(open(RES / 't1b_moving_max_T.csv'))
            if not r['mesh'].startswith('fit')]
    fit = next(r for r in csv.DictReader(open(RES / 't1b_moving_max_T.csv'))
               if r['mesh'].startswith('fit'))
    return rows, float(fit['err'].split('^')[1].rstrip(')'))


# ------------------------------------------------------------- fig 1
m = metrics()
e = {k: m[k]['vs_proxy128'] for k in
     ('n8x4 bulk_fe_rms @T', 'n16x8 bulk_fe_rms @T', 'n32x16 bulk_fe_rms @T')}
fam = {k: m[k]['family_level_diff'] for k in
       ('family 32->64 bulk_fe_rms @T', 'family 64->128 bulk_fe_rms @T')}
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
h_oof, y_oof = [0.25, 0.125, 0.0625], [float(e[k]) for k in
                                       ('n8x4 bulk_fe_rms @T', 'n16x8 bulk_fe_rms @T',
                                        'n32x16 bulk_fe_rms @T')]
h_fam, y_fam = [0.0625, 0.03125], [float(fam[k]) for k in
                                   ('family 32->64 bulk_fe_rms @T',
                                    'family 64->128 bulk_fe_rms @T')]
ax.plot(h_oof, y_oof, 'o-', ms=7, lw=1.8, color=C1,
        label='OOFEM vs proxy n128x64 (bulk FE RMS)')
ax.plot(h_fam, y_fam, 's--', ms=7, lw=1.8, color=C2,
        label='DOP853 family level diff (bulk FE RMS)')
# short slope-2 reference wedge in the empty top-left region
y0 = 4.0e-4
ax.loglog([0.025, 0.050], [y0, y0 * (0.050 / 0.025) ** 2], '--',
          color='#8a8984', lw=1.2)
ax.text(0.031, y0 * (0.031 / 0.025) ** 2 * 1.15, 'slope 2 (h$^2$)',
        color=INK2, fontsize=8)
ax.set_xlim(0.025, 0.30)
ax.set_ylim(3e-5, 3e-3)
for x, y, lab in zip(h_oof, y_oof, ('n8x4', 'n16x8', 'n32x16')):
    ax.annotate(lab, (x, y), xytext=(5, 8), textcoords='offset points',
                fontsize=8, color=INK)
for x, y, lab in zip(h_fam, y_fam, ('32$\\to$64', '64$\\to$128')):
    ax.annotate(lab, (x, y), xytext=(10, -2), textcoords='offset points',
                fontsize=8, color=INK)
style_ax(ax)
ax.set_xlabel('h (element size)')
ax.set_ylabel('bulk FE RMS error at T = 0.02')
ax.set_title('T1-B: spatial error vs mesh size at T')
ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_error_vs_h_T.png')
plt.close(fig)

# ------------------------------------------------------------- fig 2
fam_rows = fam_series()
o = orders()
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
steps = [int(r['step']) for r in fam_rows]
ax.plot(steps, [float(r['order_fam_bulk_fe_rms']) for r in fam_rows], '-',
        color=C1, lw=1.8, label='bulk FE RMS')
ax.plot(steps, [float(r['order_fam_surf_fe_rms']) for r in fam_rows], '-',
        color=C2, lw=1.8, label='surface FE RMS')
ax.axhline(2.0, color='#8a8984', ls='--', lw=1.2)
ax.text(12, 2.07, 'expected 2 (h$^2$)', color=INK2, fontsize=8)
style_ax(ax)
ax.set_xlabel('time step s  (t = s$\\cdot$DT)')
ax.set_ylabel('order $\\log_2(e_{32\\to64}/e_{64\\to128})$')
ax.set_title('T1-B: family convergence order vs time (same-time, DOP853-only)')
ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_family_order_64_128_time.png')
plt.close(fig)

# ------------------------------------------------------------- fig 3
mv, p_fit = moving_max()
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
xs = [float(r['r_c']) for r in mv]
ys = [float(r['abs_err']) for r in mv]
ax.loglog(xs, ys, 'o-', ms=7, lw=1.8, color=C1,
          label='|err| at first interior corner node (T)')
xr = np.logspace(np.log10(min(xs) * 0.8), np.log10(max(xs) * 1.3), 50)
ax.loglog(xr, np.log10(ys[0]) + p_fit * (np.log10(xr) - np.log10(xs[0])),
          '--', color='#8a8984', lw=1.2,
          label=f'fit: err $\\propto$ r$_c^{p_fit:.2f}$ (descriptive, 4 pts)')
offsets = {(0, 8), (5, -12), (-95, 6), (5, 8)}
for (x, y, lab), off in zip(zip(xs, ys, (r['mesh'] for r in mv)), offsets):
    ax.annotate(lab, (x, y), xytext=off, textcoords='offset points',
                fontsize=7.5, color=INK)
style_ax(ax)
ax.set_xlabel('r$_c$ = distance of node to nearest corner of the active pair')
ax.set_ylabel('|err| at T = 0.02')
ax.set_title('T1-B: moving maximum vs distance from corner (NOT proof of singularity)')
ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_moving_max_rc.png')
plt.close(fig)
print('figures regenerated from CSVs')
