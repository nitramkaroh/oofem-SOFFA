#!/usr/bin/env python3
"""Figure 1 - OOFEM validation: parabolic Allen-Cahn with dynamic boundary
conditions (Okumura 2021, Ch. 4, Example 2).

Panels (identical scales in comparable panels, explicit units):
  (a) literature: energy E(t) digitized from Okumura (2021) Fig. 4.4
  (b) OOFEM: field snapshots u(x, t) at six times
  (c) energy overlay: OOFEM (exact 3x3 Gauss postprocessing of the Q1 field)
      vs the digitized literature curve, over 0.01 <= t <= 0.4
  (d) difference E_OOFEM - E_paper(t) over the same interval

Colors: the dataviz skill's validated default palette, light mode
(slot 1 blue #2a78d6, slot 2 orange #eb6834); literature/analytic
references in primary ink #0b0b0b; surfaces/ink from the same reference.
Digitized data are labeled as digitized.

Run:  python3 make_figure1.py
Writes ../figures/figure_1_parabolic.{png,pdf} (figures/ is
validation_allen_cahn/figures, created if missing).
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BASE = HERE.parent                      # .../validation_allen_cahn/parabolic
FIG = HERE.parents[1] / 'figures'

# --- dataviz skill reference (light mode) -------------------------------------
SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'
BLUE = '#2a78d6'
ORANGE = '#eb6834'
# blue sequential ramp, steps 100..700 (one hue, light -> dark)
BLUE_RAMP = ['#cde2fb', '#9ec5f4', '#6da7ec', '#2a78d6', '#1c5cab', '#0d366b']

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': SURFACE,
    'axes.edgecolor': AXIS, 'axes.labelcolor': INK,
    'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7,
    'font.family': 'DejaVu Sans', 'font.size': 9.5,
    'axes.titlesize': 10.5, 'axes.titleweight': 'bold',
    'legend.frameon': False, 'lines.linewidth': 1.6,
    'savefig.facecolor': 'white',
})


def style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3.5, width=0.9)


def load_csv(name):
    return np.genfromtxt(BASE / name, delimiter=',', names=True)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    energy = load_csv('okumura_example2_energy.csv')
    paper = load_csv('literature_okumura_fig44_digitized.csv')
    diff = load_csv('okumura_example2_energy_difference.csv')
    # read with the csv module: genfromtxt sanitizes the dots out of the
    # column names (u_t_0.001 -> u_t_0001)
    with (BASE / 'okumura_example2_profiles.csv').open() as f:
        prof_rows = list(csv.reader(f))
    prof_header = prof_rows[0]
    prof_x = np.array([float(r[0]) for r in prof_rows[1:]])
    prof_cols = {c: np.array([float(r[i]) for r in prof_rows[1:]])
                 for i, c in enumerate(prof_header[1:], 1)}
    SNAP_T = [0.001, 0.1, 0.2, 0.3, 0.4, 0.7]
    with (BASE / 'okumura_example2_literature_comparison.csv').open() as f:
        summary = {k: v for k, v in csv.reader(f)}

    t_o, e_o = energy['time'], energy['OOFEM_Q1_total_energy_3x3gauss']
    t_p, e_p = paper['time'], paper['paper_energy_digitized_from_Fig4_4']
    rmse = float(summary['digitized_curve_RMSE'])
    mae = float(summary['digitized_curve_MAE'])
    emax = float(summary['digitized_curve_max_abs_error'])

    fig, axs = plt.subplots(2, 2, figsize=(11.2, 8.6))
    fig.suptitle('OOFEM validation: parabolic Allen-Cahn equation with\n'
                 'dynamic boundary conditions  (Okumura 2021, Ch. 4, Example 2)',
                 fontsize=13, fontweight='bold', color=INK, y=0.995)

    # (a) literature (digitized)
    ax = axs[0, 0]
    ax.plot(t_p, e_p, color=INK, lw=1.8)
    ax.set_title('(a) Literature energy  [digitized from Okumura (2021), Fig. 4.4]')
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('E = J(u) + F(u(0)) + F(u(L))')
    ax.set_xlim(0.0, 0.45)
    ax.set_ylim(0.0, 80.0)

    # (b) OOFEM snapshots
    ax = axs[0, 1]
    for tval, c in zip(SNAP_T, BLUE_RAMP):
        col = f'u_t_{tval:g}'
        ax.plot(prof_x, prof_cols[col], color=c, lw=1.5, label=f't = {tval:g}')
    ax.set_title('(b) OOFEM field  u(x, t)   (Q1, N = 100, $\\Delta t$ = 10$^{-3}$)')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.03, 1.05)
    ax.legend(loc='center right', fontsize=8, handlelength=1.6)

    # (c) energy overlay (comparable range with panel a)
    ax = axs[1, 0]
    mask = (t_o >= 0.01) & (t_o <= 0.4)
    ax.plot(t_p[(t_p >= 0.01) & (t_p <= 0.4)],
            e_p[(t_p >= 0.01) & (t_p <= 0.4)],
            color=INK, lw=1.8, label='literature (digitized)')
    ax.plot(t_o[mask], e_o[mask], color=ORANGE, lw=1.6, label='OOFEM (this work)')
    ax.set_title('(c) Energy: OOFEM vs literature')
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('E')
    ax.set_xlim(0.0, 0.45)
    ax.set_ylim(0.0, 80.0)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(0.012, 74.5,
            f'0.01 $\\leq$ t $\\leq$ 0.4:  RMSE = {rmse:.3f}\n'
            f'MAE = {mae:.3f},  max |diff| = {emax:.3f}\n'
            '(scale 0-90 in the source figure)',
            fontsize=8.5, color=INK2, va='top',
            bbox=dict(boxstyle='round,pad=0.35', fc='white',
                      ec=GRID, lw=0.8))

    # (d) difference
    ax = axs[1, 1]
    d = diff['difference_OOFEM_minus_paper']
    ax.axhline(0.0, color=MUTED, lw=0.9, ls=':')
    ax.plot(diff['time'], d, color=BLUE, lw=1.3)
    ax.set_title('(d) Difference  E$_{\\mathrm{OOFEM}}$ - E$_{\\mathrm{paper}}$'
                 '  (digitized)')
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('$\\Delta$E')
    ax.set_xlim(0.01, 0.4)
    lo, hi = min(-2.6, d.min() - 0.3), max(2.6, d.max() + 0.3)
    ax.set_ylim(lo, hi)
    ax.text(0.395, hi - 0.15,
            f'max $|\\Delta$E$|$ = {emax:.3f} over 0.01 $\\leq$ t $\\leq$ 0.4\n'
            'digitization uncertainty of the paper curve is not quantified;\n'
            'see comparison report for scheme differences (BE-Q1 vs FD)',
            fontsize=8, color=INK2, va='top', ha='right')

    for ax in axs.ravel():
        style(ax)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(FIG / 'figure_1_parabolic.pdf')
    fig.savefig(FIG / 'figure_1_parabolic.png', dpi=300)
    print('wrote', FIG / 'figure_1_parabolic.png')
    print(f'RMSE={rmse:.4f} MAE={mae:.4f} max={emax:.4f}')


if __name__ == '__main__':
    main()
