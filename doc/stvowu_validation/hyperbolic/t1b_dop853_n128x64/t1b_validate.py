#!/usr/bin/env python3
"""T1-B FASE 3: validation of the n128x64 proxy + comparisons (READ-ONLY on
all archives; outputs only into results/). No new integration/simulation.

Conventions identical to T1: bulk FE RMS = sqrt(e'Mb e / AREA), surf FE RMS =
sqrt(e'Ms e / PERIM) (corrected Ms), bulk max on interior nodes, surf max on
boundary; nested-grid restriction by EXACT nodal coincidence (assert d2 < 1e-26);
same-time comparisons only (OOFEM step s <-> U[s], t = s*DT).
"""
from __future__ import annotations
import csv, sys, time, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
HYPERBOLIC = HERE.parent
RES = HERE / 'results'
RES.mkdir(exist_ok=True)
sys.path.insert(0, str(HYPERBOLIC))
sys.path.insert(0, str(HYPERBOLIC / 'coupled_independent_jacobian_fix'))
import run_coupled_independent as orig
import check_jacobian_fix as fj

DT = 2.5e-4
NSTEP = 80
TIMES = np.linspace(0, orig.T, NSTEP + 1)
AREA, PERIM = 2.0, 6.0
FIX = HYPERBOLIC / 'coupled_independent_spatial_jacobian_fix'
F8 = HYPERBOLIC / 'coupled_independent_space_n8x4_dt0p00025'
F16 = HYPERBOLIC / 'coupled_independent_space_n16x8_dt0p00025'
F32NEW = HYPERBOLIC / 't1a_n32x16_rerun' / 'oofem'
T0 = time.perf_counter()

# ------------------------------------------------------------ load references
U128 = np.load(HERE / 'reference_n128x64_dt0p00025.npz')['U']
U64 = np.load(FIX / 'reference_n64x32_dt0p00025.npz')['U']
U32 = np.load(FIX / 'reference_n32x16_dt0p00025.npz')['U']
c128, c64, c32, c16, c8 = (orig.mesh(*m)[0] for m in ((128, 64), (64, 32), (32, 16), (16, 8), (8, 4)))


def cmap(nc, nf):
    idx = np.empty(len(nc), dtype=int)
    for i, (x, y) in enumerate(nc):
        d2 = (nf[:, 0] - x) ** 2 + (nf[:, 1] - y) ** 2
        j = int(np.argmin(d2))
        assert d2[j] < 1e-26, 'nested-grid coincidence broken'
        idx[i] = j
    return idx


i64 = cmap(c64, c128)
i32 = cmap(c32, c128)
i16 = cmap(c16, c128)
i8 = cmap(c8, c128)
i8_64 = cmap(c8, c64)
i16_64 = cmap(c16, c64)
i32_64 = cmap(c32, c64)
U128_64 = U128[:, i64]
print(f'loaded: U128 {U128.shape}, U64 {U64.shape}, U32 {U32.shape}; '
      f'nested maps exact (d2<1e-26)', flush=True)

# 1) finiteness / completeness / final time
assert U128.shape == (NSTEP + 1, len(c128))
assert np.isfinite(U128).all()
assert abs(TIMES[-1] - 0.02) < 1e-15
print('finiteness / completeness / final time: OK', flush=True)

# 2) quadrature invariants on n128 (from the verified cache)
d = np.load(HERE / 'matrices_n128x64.npz')
N128 = len(c128)
from scipy.sparse import coo_matrix
M128 = {}
for pfx in ('Mb', 'Kb', 'Ms', 'Ks'):
    M128[pfx] = coo_matrix((d[pfx + '_data'], (d[pfx + '_row'], d[pfx + '_col'])),
                           shape=(N128, N128)).tocsc()
ones = np.ones(N128)
inv_rows = [
    ('total Mb1 = AREA (2.0)', abs(float(ones @ (M128['Mb'] @ ones)) - AREA) < 1e-12,
     f'{float(ones @ (M128["Mb"] @ ones)):.15g}'),
    ('total Ms1 = PERIM (6.0)', abs(float(ones @ (M128['Ms'] @ ones)) - PERIM) < 1e-12,
     f'{float(ones @ (M128["Ms"] @ ones)):.15g}'),
    ('row sums Kb = 0', float(np.abs(M128['Kb'] @ ones).max()) < 1e-12,
     f'{float(np.abs(M128["Kb"] @ ones).max()):.3e}'),
    ('row sums Ks = 0', float(np.abs(M128['Ks'] @ ones).max()) < 1e-12,
     f'{float(np.abs(M128["Ks"] @ ones).max()):.3e}'),
]
for name, ok, val in inv_rows:
    print(f'quadrature n128: {name}: {val} -> {ok}', flush=True)
    assert ok

# ------------------------------------------------------------ per-mesh ctx
def ctx_of(nx, ny):
    f = fj.assemble_fixed(nx, ny)
    return f[3], f[4], f[6]   # boundary, Mb, Ms

B8, Mb8, Ms8 = ctx_of(8, 4)
B16, Mb16, Ms16 = ctx_of(16, 8)
B32, Mb32, Ms32 = ctx_of(32, 16)
B64, Mb64, Ms64 = ctx_of(64, 32)
print('mesh contexts assembled (n8..n64)', flush=True)


def metrics(Mb, Ms, bnd, e):
    m = np.zeros(len(e), dtype=bool); m[bnd] = True
    return (orig.rms_fe(Mb, e, AREA), float(np.abs(e[~m]).max()),
            orig.rms_fe(Ms, e, PERIM), float(np.abs(e[m]).max()))


def read_oofem_step(folder, step, ctarget):
    """Read OOFEM VTU step `step`, reordering nodal values into `ctarget`'s
    mesh node order by EXACT coordinate match (T1 methodology; no assumption
    on OOFEM's internal node numbering)."""
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    ds = next(x for x in pvd.findall('.//DataSet')
              if abs(float(x.attrib['timestep']) - step * DT) < 1e-12)
    vtu = ET.parse(folder / ds.attrib['file']).getroot()
    pts = np.fromstring(vtu.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)[:, :2]
    u = np.fromstring(next(a for a in vtu.findall('.//PointData/DataArray')
                           if a.attrib.get('Name') == 'Temperature').text or '', sep=' ')
    assert len(pts) == len(u) == len(ctarget), f'VTU node count mismatch (step {step})'
    perm = np.empty(len(pts), dtype=int)
    for i, (x, y) in enumerate(pts):
        d2 = (ctarget[:, 0] - x) ** 2 + (ctarget[:, 1] - y) ** 2
        j = int(np.argmin(d2))
        assert d2[j] < 1e-26, f'VTU point off mesh at step {step}'
        perm[i] = j
    out = np.empty(len(ctarget))
    out[perm] = u
    assert np.isfinite(out).all()
    return out


# ------------------------------------------------- 3) common-time T metrics
def zones(c):
    x, y = c[:, 0], c[:, 1]
    rA = np.minimum(np.hypot(x, y - 1.0), np.hypot(x - 2.0, y))
    rB = np.minimum(np.hypot(x, y), np.hypot(x - 2.0, y - 1.0))
    zA = rA < 0.5
    return zA, (rB < 0.5) & ~zA, ~(zA | (rB < 0.5))


# OOFEM values (read once, reused by sections 3-5)
OOF = {}


def oof(folder, ctarget, s):
    if (folder, s) not in OOF:
        OOF[(folder, s)] = read_oofem_step(folder, s, ctarget)
    return OOF[(folder, s)]


S = 80
u8, u16, u32 = oof(F8, c8, S), oof(F16, c16, S), oof(F32NEW, c32, S)
e8_128 = metrics(Mb8, Ms8, B8, u8 - U128[S, i8])
e16_128 = metrics(Mb16, Ms16, B16, u16 - U128[S, i16])
e32_128 = metrics(Mb32, Ms32, B32, u32 - U128[S, i32])
e8_64 = metrics(Mb8, Ms8, B8, u8 - U64[S, i8_64])
e16_64 = metrics(Mb16, Ms16, B16, u16 - U64[S, i16_64])
e32_64 = metrics(Mb32, Ms32, B32, u32 - U64[S, i32_64])
# family level differences at T (DOP853-only, implementation-independent)
e32_64fam = metrics(Mb32, Ms32, B32, U32[S] - U64[S, i32_64])
e64_128fam = metrics(Mb64, Ms64, B64, U64[S] - U128[S, i64])

MT = ['bulk_fe_rms', 'bulk_max', 'surf_fe_rms', 'surf_max']
with (RES / 't1b_metrics_T.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['quantity', 'vs_proxy128', 'vs_proxy64', 'family_level_diff'])
    for k, name in enumerate(MT):
        w.writerow([f'n8x4 {name} @T', f'{e8_128[k]:.6e}', f'{e8_64[k]:.6e}', ''])
    for k, name in enumerate(MT):
        w.writerow([f'n16x8 {name} @T', f'{e16_128[k]:.6e}', f'{e16_64[k]:.6e}', ''])
    for k, name in enumerate(MT):
        w.writerow([f'n32x16 {name} @T', f'{e32_128[k]:.6e}', f'{e32_64[k]:.6e}', ''])
    for k, name in enumerate(MT):
        w.writerow([f'family 32->64 {name} @T', '', '', f'{e32_64fam[k]:.6e}'])
    for k, name in enumerate(MT):
        w.writerow([f'family 64->128 {name} @T', '', '', f'{e64_128fam[k]:.6e}'])
orders = {}
for k, name in enumerate(MT):
    orders[f'8->16 {name} @T (vs proxy128)'] = float(np.log2(e8_128[k] / e16_128[k]))
    orders[f'16->32 {name} @T (vs proxy128)'] = float(np.log2(e16_128[k] / e32_128[k]))
    orders[f'family 32->64 {name} @T'] = float(np.log2(e32_64fam[k] / e64_128fam[k]))
with (RES / 't1b_family_orders_T.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['order', 'value'])
    for k, v in orders.items():
        w.writerow([k, f'{v:.4f}'])
        print(f'order {k}: {v:.4f}', flush=True)
# proxy-level stability
with (RES / 't1b_proxy_stability_T.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mesh', 'metric', 'err_vs_proxy64', 'err_vs_proxy128',
                'abs_diff', 'family_64_128_at_this_level'])
    for tag, e64v, e128v in (('n8x4', e8_64, e8_128), ('n16x8', e16_64, e16_128),
                             ('n32x16', e32_64, e32_128)):
        for k, name in enumerate(MT):
            w.writerow([tag, name, f'{e64v[k]:.6e}', f'{e128v[k]:.6e}',
                        f'{abs(e128v[k] - e64v[k]):.6e}', ''])

# ------------------------------------------------- 4) zones at T
zrows = []
for tag, c, e in (('family 64->128 @T (n64 nodes)', c64, U64[S] - U128[S, i64]),
                  ('OOFEM n32x16 vs proxy128 @T', c32, u32 - U128[S, i32])):
    zA, zB, zR = zones(c)
    tot = float(np.sum(e ** 2))
    for zname, z in (('cornerA (0,1)/(2,0)', zA), ('cornerB (0,0)/(2,1)', zB),
                     ('remainder', zR)):
        zrows.append({'dataset': tag, 'zone': zname,
                      'nodal_rms': float(np.sqrt(np.mean(e[z] ** 2))),
                      'max': float(np.abs(e[z]).max()),
                      'energy_frac': float(np.sum(e[z] ** 2) / tot)})
with (RES / 't1b_zones_T.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['dataset', 'zone', 'nodal_rms', 'max', 'energy_frac'])
    w.writeheader()
    for r in zrows:
        w.writerow({k: (f'{r[k]:.6e}' if isinstance(r[k], float) else r[k]) for k in r})
for r in zrows:
    print(f'zone {r["dataset"]} {r["zone"]}: energy={r["energy_frac"]:.3f} '
          f'rms={r["nodal_rms"]:.3e}', flush=True)

# ------------------------------------------------- 5) moving maximum (T)
def first_interior_corner_node(c):
    x, y = c[:, 0], c[:, 1]
    bnd = np.zeros(len(c), dtype=bool)
    bnd[(x == 0) | (x == 2.0) | (y == 0) | (y == 1.0)] = True
    inter = ~bnd
    rA = np.hypot(x, y - 1.0)
    rB = np.hypot(x - 2.0, y)
    r = np.minimum(rA, rB)
    r[~inter] = np.inf
    j = int(np.argmin(r))
    return j, float(r[j])


mv = []
for tag, c, bnd, err in (('n8x4 OOFEM', c8, B8, u8 - U128[S, i8]),
                         ('n16x8 OOFEM', c16, B16, u16 - U128[S, i16]),
                         ('n32x16 OOFEM-new', c32, B32, u32 - U128[S, i32]),
                         ('n64 family diff', c64, B64, U64[S] - U128[S, i64])):
    j, rj = first_interior_corner_node(c)
    xm, ym = c[j, 0], c[j, 1]
    assert not (xm == 0 or xm == 2.0 or ym == 0 or ym == 1.0), \
        'first-corner node is on the boundary (unexpected)'
    mv.append({'mesh': tag, 'x': c[j, 0], 'y': c[j, 1], 'r_c': rj,
               'err': float(err[j]), 'abs_err': float(abs(err[j]))})
rs = np.array([m['r_c'] for m in mv]); es = np.array([m['abs_err'] for m in mv])
p_fit = float(np.polyfit(np.log(rs), np.log(es), 1)[0])
with (RES / 't1b_moving_max_T.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(mv[0].keys()))
    w.writeheader()
    for m in mv:
        w.writerow({k: (f'{v:.10e}' if isinstance(v, float) else v) for k, v in m.items()})
    w.writerow({'mesh': 'fit (4 pts, descriptive)', 'x': '', 'y': '', 'r_c': '',
                'err': f'err ~ r_c^{p_fit:.3f}', 'abs_err': ''})
print(f'moving max: ' + '  '.join(f'{m["mesh"]} r_c={m["r_c"]:.4f} err={m["abs_err"]:.3e}'
                                  for m in mv), flush=True)
print(f'  4-point log-log slope (descriptive, NOT proof): p = {p_fit:.3f}', flush=True)

# ------------------------------------------------- 6) time series (same-time)
# family order(s) = log2( ||U32-U64||_{n32} / ||U64-U128||_{n64} ), both at t=s*DT
fam = []
oofrows = []
for s in range(1, NSTEP + 1):
    e64f = metrics(Mb64, Ms64, B64, U64[s] - U128[s, i64])
    e32f = metrics(Mb32, Ms32, B32, U32[s] - U64[s, i32_64])
    fam.append({'step': s, 't': f'{s * DT:.5f}',
                **{f'fam32_64_{n}': f'{e32f[k]:.6e}' for k, n in enumerate(MT)},
                **{f'fam64_128_{n}': f'{e64f[k]:.6e}' for k, n in enumerate(MT)},
                **{f'order_fam_{n}': f'{np.log2(e32f[k] / e64f[k]):.4f}'
                   for k, n in enumerate(MT)}})
    if s % 10 == 0 or s in (1, 39):
        e8o = metrics(Mb8, Ms8, B8, oof(F8, c8, s) - U128[s, i8])
        e16o = metrics(Mb16, Ms16, B16, oof(F16, c16, s) - U128[s, i16])
        e32o = metrics(Mb32, Ms32, B32, oof(F32NEW, c32, s) - U128[s, i32])
        oofrows.append({'step': s, 't': f'{s * DT:.5f}',
                        **{f'e8_{n}': f'{e8o[k]:.6e}' for k, n in enumerate(MT)},
                        **{f'e16_{n}': f'{e16o[k]:.6e}' for k, n in enumerate(MT)},
                        **{f'e32_{n}': f'{e32o[k]:.6e}' for k, n in enumerate(MT)},
                        **{f'ord8_16_{n}': f'{np.log2(e8o[k] / e16o[k]):.4f}'
                           for k, n in enumerate(MT)},
                        **{f'ord16_32_{n}': f'{np.log2(e16o[k] / e32o[k]):.4f}'
                           for k, n in enumerate(MT)}})
# full family series (all 80 steps)
with (RES / 't1b_orders_time_64_128.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(fam[0].keys()))
    w.writeheader(); w.writerows(fam)
# OOFEM series vs proxy128 at selected steps (1,10,20,30,39,40,50,60,70,80)
with (RES / 't1b_orders_time_oofem_proxy128.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(oofrows[0].keys()))
    w.writeheader(); w.writerows(oofrows)
print('time series written (family all 80 steps; OOFEM selected steps)', flush=True)
print(f'  family order 64->128 bulk FE RMS: step1={fam[0]["order_fam_bulk_fe_rms"]}, '
      f'step20={fam[19]["order_fam_bulk_fe_rms"]}, step40={fam[39]["order_fam_bulk_fe_rms"]}, '
      f'step80={fam[79]["order_fam_bulk_fe_rms"]}', flush=True)

# ------------------------------------------------- 7) energy of n128 trajectory
V128 = np.load(HERE / 'reference_n128x64_dt0p00025.npz')['V']
quads128 = d['quads']
edges128 = d['edges']
ens = []
for s in range(NSTEP + 1):
    e = fj.energy_fixed(U128[s], V128[s], c128, quads128, edges128)
    ens.append((TIMES[s], *e))
    if s % 20 == 0:
        print(f'  energy t={TIMES[s]:.5f}: total={e[0]:.8f}', flush=True)
with (RES / 't1b_energy_n128.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['time', 'total', 'bulk_kinetic',
                'bulk_potential (gradient + double well)', 'surface'])
    for r in ens:
        w.writerow([f'{r[0]:.5f}'] + [f'{x:.10e}' for x in r[1:]])

# ------------------------------------------------- 8) figures (static PNG)
# Colors: skill default categorical slots 1-2 (validated by the dataviz skill;
# no Node runtime available to re-run the validator here). Secondary encoding:
# marker shape + line style, so identity is never color-alone.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C1, C2 = '#2a78d6', '#eb6834'          # blue, orange
INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'text.color': INK, 'axes.edgecolor': '#c9c8c4',
                     'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.facecolor': '#fcfcfb', 'figure.facecolor': '#fcfcfb',
                     'grid.color': '#c9c8c4', 'grid.alpha': 0.55, 'font.size': 9})


def style_ax(ax):
    ax.grid(True, lw=0.6)


# fig 1: bulk FE RMS error vs h at T — OOFEM chain + DOP853 family diffs
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
hmin, hmax = 0.025, 0.30
ax.loglog([hmin, hmax], [(h / 0.25) ** 2 for h in (hmin, hmax)],
          '--', color='#8a8984', lw=1.2)
ax.text(hmax * 0.55, (hmax / 0.25) ** 2 * 1.15, 'h$^2$', color=INK2)
h_oof, y_oof = [0.25, 0.125, 0.0625], [e8_128[0], e16_128[0], e32_128[0]]
h_fam, y_fam = [0.0625, 0.03125], [e32_64fam[0], e64_128fam[0]]
ax.plot(h_oof, y_oof, 'o-', ms=7, lw=1.8, color=C1,
        label='OOFEM vs proxy n128x64 (bulk FE RMS)')
ax.plot(h_fam, y_fam, 's--', ms=7, lw=1.8, color=C2,
        label='DOP853 family level diff (bulk FE RMS)')
for x, y, lab in zip(h_oof, y_oof, ('n8x4', 'n16x8', 'n32x16')):
    ax.annotate(lab, (x, y), xytext=(5, 8), textcoords='offset points',
                fontsize=8, color=INK)
for x, y, lab in zip(h_fam, y_fam, ('32$\\to$64', '64$\\to$128')):
    ax.annotate(lab, (x, y), xytext=(5, -12), textcoords='offset points',
                fontsize=8, color=INK)
style_ax(ax)
ax.set_xlabel('h (element size)')
ax.set_ylabel('bulk FE RMS error at T = 0.02')
ax.set_title('T1-B: spatial error vs mesh size at T')
ax.legend(loc='lower left', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_error_vs_h_T.png')
plt.close(fig)

# fig 2: family order (32->64 over 64->128) vs time, FE RMS bulk/surface
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
steps = [r['step'] for r in fam]
ax.plot(steps, [float(r['order_fam_bulk_fe_rms']) for r in fam], '-',
        color=C1, lw=1.8, label='bulk FE RMS')
ax.plot(steps, [float(r['order_fam_surf_fe_rms']) for r in fam], '-',
        color=C2, lw=1.8, label='surface FE RMS')
ax.axhline(2.0, color='#8a8984', ls='--', lw=1.2)
ax.text(12, 2.07, 'expected 2 (h$^2$)', color=INK2, fontsize=8)
style_ax(ax)
ax.set_xlabel('time step s  (t = s$\\cdot$DT)')
ax.set_ylabel('order $\\log_2(e_{32\\to64}/e_{64\\to128})$')
ax.set_title('T1-B: family convergence order vs time (same-time, DOP853-only)')
ax.legend(loc='best', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_family_order_64_128_time.png')
plt.close(fig)

# fig 3: moving maximum — |err| at first interior corner node vs r_c (4 points)
fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=150)
xs, ys = [m['r_c'] for m in mv], [m['abs_err'] for m in mv]
ax.loglog(xs, ys, 'o-', ms=7, lw=1.8, color=C1,
          label='|err| at first interior corner node (T)')
xr = np.logspace(np.log10(min(xs) * 0.8), np.log10(max(xs) * 1.3), 50)
ax.loglog(xr, np.log10(ys[0]) + p_fit * (np.log10(xr) - np.log10(xs[0])),
          '--', color='#8a8984', lw=1.2,
          label=f'fit: err $\\propto$ r$_c^{p_fit:.2f}$ (descriptive, 4 pts)')
for x, y, lab in zip(xs, ys, (f'${m_["mesh"]}$' for m_ in mv)):
    ax.annotate(lab, (x, y), xytext=(5, 8), textcoords='offset points',
                fontsize=7.5, color=INK)
style_ax(ax)
ax.set_xlabel('r$_c$ = distance of node to nearest corner of the active pair')
ax.set_ylabel('|err| at T = 0.02')
ax.set_title('T1-B: moving maximum vs distance from corner (NOT proof of singularity)')
ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
fig.tight_layout()
fig.savefig(RES / 't1b_moving_max_rc.png')
plt.close(fig)
print('figures written', flush=True)

print(f'\nDONE wall={time.perf_counter() - T0:.1f}s', flush=True)
