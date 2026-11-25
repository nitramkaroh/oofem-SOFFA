#!/usr/bin/env python3
"""Figure 2 - OOFEM validation: hyperbolic Allen-Cahn traveling wave
(Nizovtseva, Galenko, Alexandrov 2016, J. Phys. A 49, 435201, Eq. (29),
bounded wave, tau = 0.2, D = M = 1).

Panels (identical scales in comparable panels, explicit units):
  (a) exact profiles u(x, t) = 1/(1 + exp[lambda (x - c t)])
  (b) OOFEM profiles, finest mesh (n = 320)
  (c) overlay exact vs OOFEM at five times
  (d) interface position x_0(t) (u = 1/2): exact c t vs OOFEM n = 320,
      with an inset of the spatial RMS error vs h and the O(h^2) reference

The traveling-wave test exercises the inertial bulk model only: no dynamic
boundary condition is present in this benchmark (declared in the report).

Colors: dataviz skill validated default palette, light mode: slot 1 blue
#2a78d6, slot 2 orange #eb6834 (one-hue ramps derived by blending with the
chart surface #fcfcfb); exact references in the documented ink grays
#0b0b0b -> #898781.

Run:  python3 make_figure2.py
Writes ../figures/figure_2_hyperbolic.{png,pdf}.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BASE = HERE.parent                      # .../validation_allen_cahn/hyperbolic
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
SURFACE_RGB = np.array([252, 252, 251])
ORANGE_RGB = np.array([235, 104, 52])
MUTED_RGB = np.array([137, 135, 129])
INK_RGB = np.array([11, 11, 11])


def ramp(base, steps):
    out = []
    for t in steps:
        v = (t * base + (1 - t) * SURFACE_RGB).round().astype(int)
        out.append('#%02x%02x%02x' % tuple(v))
    return out


# one-hue sequential ramps, light -> dark for increasing time
GRAY_FULL = ramp(MUTED_RGB, [0.25, 0.48, 0.68, 0.85, 1.00])  # t = 0 ... 0.2
GRAY_FULL[-1] = INK                                          # last: full ink
ORANGE_FULL = ramp(ORANGE_RGB, [0.30, 0.52, 0.72, 0.87, 1.00])

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

TAU = 0.2
LAMBDA = np.sqrt(2.0 + 9.0 * TAU) / 2.0
SPEED = 1.5 / LAMBDA
SNAP_T = [0.0, 0.05, 0.10, 0.15, 0.20]
XMIN, XMAX = -8.0, 8.0


def exact(x, t):
    return 1.0 / (1.0 + np.exp(np.clip(LAMBDA * (x - SPEED * t), -700, 700)))


def style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3.5, width=0.9)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    # OOFEM profiles (n = 320)
    prof = {}
    with (BASE / 'nizovtseva_profiles.csv').open() as f:
        for row in csv.DictReader(f):
            key = (float(row['time']), int(row['n']))
            prof.setdefault(key, ([], []))
            prof[key][0].append(float(row['x']))
            prof[key][1].append(float(row['u']))
    # metrics
    metrics = []
    with (BASE / 'nizovtseva_metrics.csv').open() as f:
        for row in csv.DictReader(f):
            metrics.append(row)
    m320 = [m for m in metrics if m['n'] == '320'][0]
    rms320, mx320 = float(m320['rms_full']), float(m320['max_full'])
    # interface
    it = []
    with (BASE / 'nizovtseva_interface_n320.csv').open() as f:
        for row in csv.DictReader(f):
            it.append((float(row['time']), float(row['x0_OOFEM_u1_2'])))
    it_t = np.array([a for a, b in it])
    it_x = np.array([b for a, b in it])

    fig, axs = plt.subplots(2, 2, figsize=(11.2, 8.6))
    fig.suptitle('OOFEM validation: hyperbolic Allen-Cahn traveling wave\n'
                 '($\\tau_R u_{tt} + u_t = \\Delta u + u - u^3$, $\\tau$ = 0.2; '
                 'Nizovtseva et al. 2016, Eq. (29))\n'
                 'bulk inertial model - no dynamic boundary condition in this test',
                 fontsize=12.5, fontweight='bold', color=INK, y=0.995)

    xg = np.linspace(XMIN, XMAX, 2401)

    # (a) exact profiles
    ax = axs[0, 0]
    for tval, c in zip(SNAP_T, GRAY_FULL):
        ax.plot(xg, exact(xg, tval), color=c, lw=1.5, label=f't = {tval:g}')
    ax.set_title('(a) Exact traveling wave  u(x, t)')
    ax.set_xlabel('x (dimensionless)')
    ax.set_ylabel('u')
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(-0.03, 1.05)
    ax.legend(loc='lower right', fontsize=8, handlelength=1.6)

    # (b) OOFEM profiles, n = 320
    ax = axs[0, 1]
    for tval, c in zip(SNAP_T, ORANGE_FULL):
        xs, us = prof[(tval, 320)]
        ax.plot(xs, us, color=c, lw=1.5, label=f't = {tval:g}')
    ax.set_title('(b) OOFEM  u(x, t)   (Q1, n = 320, $\\Delta t$ = 10$^{-3}$)')
    ax.set_xlabel('x (dimensionless)')
    ax.set_ylabel('u')
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(-0.03, 1.05)
    ax.legend(loc='lower right', fontsize=8, handlelength=1.6)

    # (c) overlay exact vs OOFEM
    ax = axs[1, 0]
    for tval, cg, co in zip(SNAP_T, GRAY_FULL, ORANGE_FULL):
        xs, us = prof[(tval, 320)]
        ax.plot(xg, exact(xg, tval), color=cg, lw=1.5)
        ax.plot(xs, us, color=co, lw=1.5, ls=(0, (4, 2)))
    ax.set_title('(c) Overlay: exact (solid) vs OOFEM n = 320 (dashed)')
    ax.set_xlabel('x (dimensionless)')
    ax.set_ylabel('u')
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(-0.03, 1.05)
    ax.legend([plt.Line2D([], [], color=INK, lw=1.5),
               plt.Line2D([], [], color=ORANGE, lw=1.5, ls=(0, (4, 2)))],
              ['exact', 'OOFEM (n = 320)'], loc='lower right', fontsize=8)
    ax.text(0.015, 1.01,
            f'at T = 0.2:  RMS = {rms320:.2e},  max = {mx320:.2e}\n'
            'deviation below the visible thickness of the lines',
            fontsize=8, color=INK2, va='top')

    # (d) interface position + convergence inset
    ax = axs[1, 1]
    tgrid = np.linspace(0.0, 0.2, 41)
    ax.plot(tgrid, SPEED * tgrid, color=INK, lw=1.8,
            label='exact  $x_0(t) = c t$')
    ax.plot(it_t, it_x, color=ORANGE, lw=1.5, label='OOFEM  (n = 320)')
    ax.set_title('(d) Interface  $x_0(t)$  (level u = 1/2)')
    ax.set_xlabel('t (dimensionless)')
    ax.set_ylabel('$x_0$')
    ax.set_xlim(0.0, 0.21)
    ax.set_ylim(0.0, 0.33)
    ax.set_xticks([0.0, 0.05, 0.10, 0.15, 0.20])
    ax.legend(loc='lower right', fontsize=9)
    x0err = float(it_x[-1]) - SPEED * float(it_t[-1])
    sperr = abs(float(it_x[-1]) / float(it_t[-1]) / SPEED - 1.0)
    ax.text(0.008, 0.150,
            f'$c$ = {SPEED:.5f}   (exact speed)\n'
            f'$|x_0(T) - cT|$ = {abs(x0err):.2e},  '
            f'rel. speed error = {sperr:.1e}',
            fontsize=8.5, color=INK2, va='top')
    # inset: spatial convergence
    # inset in the top-left region, away from the x_0(t) line
    axins = ax.inset_axes([0.05, 0.52, 0.46, 0.42])
    hs = np.array([float(m['h']) for m in metrics])
    rms = np.array([float(m['rms_full']) for m in metrics])
    axins.loglog(hs, rms, 'o-', color=BLUE, lw=1.4, ms=5)
    ref = rms[0] * (hs[0] / hs) ** 2
    axins.loglog(hs, ref, '--', color=MUTED, lw=1.1)
    o1 = np.log2(rms[0] / rms[1])
    o2 = np.log2(rms[1] / rms[2])
    axins.set_xticks([0.3, 0.15, 0.075])
    axins.set_xticklabels(['0.3', '0.15', '0.075'], fontsize=7.5)
    axins.set_yticks([1e-5, 1e-6])
    axins.set_yticklabels(['1e-5', '1e-6'], fontsize=7.5)
    axins.tick_params(length=2.5, width=0.7)
    for s in axins.spines.values():
        s.set_linewidth(0.7)
    axins.set_ylabel('RMS', fontsize=8, color=INK2)
    axins.set_title(f'RMS vs h  (Q1, $\\Delta t$ = 10$^{-3}$): order '
                    f'$\\approx$ {o1:.2f} / {o2:.2f}', fontsize=8, color=INK2)

    for ax in axs.ravel():
        style(ax)
    fig.tight_layout(rect=(0, 0, 1, 0.935))
    fig.savefig(FIG / 'figure_2_hyperbolic.pdf')
    fig.savefig(FIG / 'figure_2_hyperbolic.png', dpi=300)
    print('wrote', FIG / 'figure_2_hyperbolic.png')
    print(f'rms_full: {rms}')
    print(f'orders: {o1:.3f}, {o2:.3f}')
    print(f'x0(T) err = {x0err:.3e}, rel speed err = {sperr:.3e}')


if __name__ == '__main__':
    main()
