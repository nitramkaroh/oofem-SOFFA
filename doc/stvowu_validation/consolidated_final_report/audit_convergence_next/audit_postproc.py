#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_postproc.py — Audit indipendente in sola lettura (audit_convergence_next)

Scopo
-----
Ricalcolo indipendente (post-processing) di tutte le grandezze numeriche
riferite al caso accoppiato iperbolico, a partire dai DATI GRI già archiviati:

  * riferimenti DOP853 (corretti, post-Jacobian-fix):
      - n16x8 : coupled_independent_spatial_jacobian_fix/tolcheck_n16x8_rtol2e-12.npz
      - n32x16: coupled_independent_spatial_jacobian_fix/reference_n32x16_dt0p00025.npz
      - n64x32: coupled_independent_spatial_jacobian_fix/reference_n64x32_dt0p00025.npz
      - n128x64: t1b_dop853_n128x64/reference_n128x64_dt0p00025.npz  (U e V)
  * run OOFEM (VTU, Float64, campo "Temperature"):
      - n8x4, n16x8 : coupled_independent_space_n{8x4,16x8}_dt0p00025/
      - n32x16      : t1a_n32x16_rerun/oofem/   (rerun bit-identico steps 1-39)
  * matrici/quadrature n128: t1b_dop853_n128x64/matrices_n128x64.npz

NON modifica, NON sovrascrive e NON avvia nulla: nessuna simulazione OOFEM,
nessuna integrazione DOP853, nessun pacchetto nuovo, niente matrici dense di
grandi dimensioni (solo assembly Q1 sparse di n<=64, già usate in T1-B).

Scrive SOLO in: audit_convergence_next/results/*.csv
Stampa un digest finale.

Livelli di indipendenza (documentati nell'audit_report.md):
  * mesh, quadrature, norme, zone, fit, energia, parsing VTU: codice NUOVO qui;
  * assembly Mb/Ms: doppia implementazione (qui + check_jacobian_fix.assemble_fixed,
    verificato nella campagna Jacobian-fix); le differenze massime sono riportate.
"""
import csv
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path('/home/marco/Documents/Work/oofem-SOFFA/doc/stvowu_validation')
HYPER = ROOT / 'hyperbolic'
FIX = HYPER / 'coupled_independent_spatial_jacobian_fix'
JACFIX = HYPER / 'coupled_independent_jacobian_fix'  # check_jacobian_fix.py (T1)
T1B = HYPER / 't1b_dop853_n128x64'
F8 = HYPER / 'coupled_independent_space_n8x4_dt0p00025'
F16 = HYPER / 'coupled_independent_space_n16x8_dt0p00025'
F32A = HYPER / 't1a_n32x16_rerun' / 'oofem'
OUT = Path(__file__).resolve().parent / 'results'
OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------- costanti caso --
LX, LY = 2.0, 1.0
DT, NSTEP, TEND = 2.5e-4, 80, 0.02
AREA, PERIM = 2.0, 6.0
S = NSTEP  # step finale (t = T)

sys.path.insert(0, str(HYPER))
sys.path.insert(0, str(JACFIX))
import run_coupled_independent as orig      # noqa: E402  (solo costanti/funzioni)
import check_jacobian_fix as fj             # noqa: E402  (assembly verificato)

CHECKS = []  # righe per integrity_checks.csv


def log(check, detail, value, ok):
    CHECKS.append((check, detail, value, 'OK' if ok else 'FAIL'))
    print(f'  [{"OK " if ok else "FAIL"}] {check}: {detail} = {value}', flush=True)
    return ok


def write_csv(name, header, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f'  -> results/{name}: {len(rows)} righe', flush=True)


# ------------------------------------------------------------------- mesh ---
def mesh(nx, ny):
    """Griglia Q1 canonica, node(i,j) = j*(nx+1)+i (x-fasta), come T1/T1-B."""
    x = np.linspace(0.0, LX, nx + 1)
    y = np.linspace(0.0, LY, ny + 1)
    X, Y = np.meshgrid(x, y, indexing='ij')
    coords = np.column_stack([X.ravel(order='F'), Y.ravel(order='F')])
    quads = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            quads.append((a, a + 1, a + nx + 2, a + nx + 1))
    edges = []
    for i in range(nx):  # basso (y=0) e alto (y=LY)
        edges.append((i, i + 1))
        a = ny * (nx + 1) + i
        edges.append((a, a + 1))
    for j in range(ny):  # sinistro (x=0) e destro (x=LX)
        edges.append((j * (nx + 1), (j + 1) * (nx + 1)))
        a = j * (nx + 1) + nx
        edges.append((a, a + nx + 1))
    bnd = (coords[:, 0] == 0) | (coords[:, 0] == LX) | (coords[:, 1] == 0) | \
          (coords[:, 1] == LY)
    return coords, np.array(quads, dtype=int), np.array(edges, dtype=int), bnd


def restrict(c_coarse, c_fine, tag):
    """Mappa di restrizione per coincidenza esatta di coordinate (d2 < 1e-26)."""
    out = np.empty(len(c_coarse), dtype=int)
    for i, (x, y) in enumerate(c_coarse):
        d2 = (c_fine[:, 0] - x) ** 2 + (c_fine[:, 1] - y) ** 2
        j = int(np.argmin(d2))
        n_match = int((d2 < 1e-26).sum())
        if d2[j] >= 1e-26 or n_match != 1:
            log('restrizione', f'{tag} nodo {i} ({x},{y})',
                f'd2min={d2[j]:.3e} n_match={n_match}', False)
        out[i] = j
    return out


# ------------------------------------------------------------ quadrature mie --
# Gauss 3 punti su [-1,1] (esatte per polinomi <= grado 5)
GX = np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)])
GW = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


def qshape(x, y):
    N = np.array([.25 * (1 - x) * (1 - y), .25 * (1 + x) * (1 - y),
                  .25 * (1 + x) * (1 + y), .25 * (1 - x) * (1 + y)])
    dN = np.array([[-.25 * (1 - y), -.25 * (1 - x)],
                   [.25 * (1 - y), -.25 * (1 + x)],
                   [.25 * (1 + y), .25 * (1 + x)],
                   [-.25 * (1 + y), .25 * (1 - x)]])
    return N, dN  # dN = [dN/dxi, dN/deta] riga per nodo


def assemble_mine(coords, quads, edges, bnd):
    """Assemblaggio Q1 indipendente (bulk 3x3, linee con J = L/2 corretto)."""
    nd = len(coords)
    Mb = np.zeros((nd, nd))
    Ms = np.zeros((nd, nd))
    for q in quads:
        xy = coords[q]
        for a in range(3):
            for b in range(3):
                N, dN = qshape(GX[a], GX[b])
                J = xy.T @ dN
                det = np.linalg.det(J)
                w = det * GW[a] * GW[b]
                Mb[np.ix_(q, q)] += np.outer(N, N) * w
    for e in edges:
        xy = coords[e]
        L = float(np.linalg.norm(xy[1] - xy[0]))
        for a in range(3):
            s = GX[a]
            N = np.array([(1 - s) / 2, (1 + s) / 2])
            dN = np.array([-.5, .5]) * (2.0 / L)
            w = 0.5 * L * GW[a]
            Ms[np.ix_(e, e)] += np.outer(N, N) * w
    return Mb, Ms


def rms_fe(M, e, measure):
    v = float(e @ (M @ e))
    return float(np.sqrt(max(v / measure, 0.0)))


def metrics4(Mb, Ms, bnd, e):
    m = np.zeros(len(e), dtype=bool)
    m[bnd] = True
    return (rms_fe(Mb, e, AREA), float(np.abs(e[~m]).max()),
            rms_fe(Ms, e, PERIM), float(np.abs(e[m]).max()))


# ------------------------------------------------------------------ energia --
def my_energy(u, v, coords, quads, edges):
    """E = ek + ep + es + esp con Gauss 3x3 (Q1) e 3 punti (P1 linea).
    Esatta per quartici Q1/P1. Parametri tutti = 1."""
    ek = ep = es = esp = 0.0
    for a in range(3):
        for b in range(3):
            x, y = GX[a], GX[b]
            wa = GW[a] * GW[b]
            N, dN = qshape(x, y)
            e = quads
            xy = coords[e]
            x0, y0 = xy[:, :, 0], xy[:, :, 1]
            J11 = dN[:, 0] @ x0.T
            J12 = dN[:, 0] @ y0.T
            J21 = dN[:, 1] @ x0.T
            J22 = dN[:, 1] @ y0.T
            det = J11 * J22 - J12 * J21
            uel, vel = u[e], v[e]
            uq = (N[0] * uel[:, 0] + N[1] * uel[:, 1]
                  + N[2] * uel[:, 2] + N[3] * uel[:, 3])
            vq = (N[0] * vel[:, 0] + N[1] * vel[:, 1]
                  + N[2] * vel[:, 2] + N[3] * vel[:, 3])
            coef_x = (dN[:, 0][:, None] * J22[None, :]
                      - dN[:, 1][:, None] * J21[None, :]) / det[None, :]
            coef_y = (-dN[:, 0][:, None] * J12[None, :]
                      + dN[:, 1][:, None] * J11[None, :]) / det[None, :]
            gux = np.einsum('ie,ei->e', coef_x, uel)
            guy = np.einsum('ie,ei->e', coef_y, uel)
            w = det * wa
            ek += 0.5 * float(np.sum(vq * vq * w))
            ep += 0.5 * float(np.sum((gux * gux + guy * guy) * w))
            es += 0.25 * float(np.sum((uq * uq - 1.0) ** 2 * w))
    for a in range(3):
        s = GX[a]
        e = edges
        xy = coords[e]
        L = np.linalg.norm(xy[:, 1] - xy[:, 0], axis=1)
        ul = u[e]
        uq = 0.5 * (1 - s) * ul[:, 0] + 0.5 * (1 + s) * ul[:, 1]
        du = (ul[:, 1] - ul[:, 0]) / L
        w = 0.5 * L * GW[a]
        esp += 0.5 * float(np.sum(du * du * w)) \
            + 0.25 * float(np.sum((uq * uq - 1.0) ** 2 * w))
    return ek + ep + es + esp, ek, ep + es, esp


# ------------------------------------------------------------------- VTU -----
VTOF = {}


def read_oofem_step(folder, step, ctarget, tag):
    key = (str(folder), int(step))
    if key in VTOF:
        return VTOF[key]
    pvd = ET.parse(Path(folder) / 'case.out.m0.pvd').getroot()
    ds = next(x for x in pvd.findall('.//DataSet')
              if abs(float(x.attrib['timestep']) - step * DT) < 1e-12)
    vtu = ET.parse(Path(folder) / ds.attrib['file']).getroot()
    pts = np.fromstring(vtu.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)
    u = np.fromstring(next(a for a in vtu.findall('.//PointData/DataArray')
                           if a.attrib.get('Name') == 'Temperature').text or '',
                      sep=' ')
    assert len(pts) == len(u) == len(ctarget), f'VTU mismatch {tag} step {step}'
    perm = np.empty(len(pts), dtype=int)
    worst = 0.0
    for i, (x, y) in enumerate(pts[:, :2]):
        d2 = (ctarget[:, 0] - x) ** 2 + (ctarget[:, 1] - y) ** 2
        j = int(np.argmin(d2))
        assert d2[j] < 1e-26 and (d2 < 1e-26).sum() == 1, f'VTU off-mesh {tag}'
        worst = max(worst, float(d2[j]))
        perm[i] = j
    out = np.empty(len(ctarget))
    out[perm] = u
    assert np.isfinite(out).all()
    VTOF[key] = out
    log('VTU', f'{tag} step {step}', f'node={len(pts)} max_d2={worst:.1e}', True)
    return out


# ================================================================ setup =======
t_start = time.time()
print('== setup: mesh + assembly', flush=True)

MESHS = {}
for nx, ny, tag in ((8, 4, 'n8'), (16, 8, 'n16'), (32, 16, 'n32'),
                    (64, 32, 'n64'), (128, 64, 'n128')):
    MESHS[tag] = mesh(nx, ny)

NXNY = {'n8': (8, 4), 'n16': (16, 8), 'n32': (32, 16), 'n64': (64, 32)}
for tag in ('n8', 'n16', 'n32', 'n64'):
    c, q, e, b = MESHS[tag]
    Mb_m, Ms_m = assemble_mine(c, q, e, b)
    f = fj.assemble_fixed(*NXNY[tag])
    f_mb, f_ms = f[4], f[6]
    d_mb = float(np.abs(Mb_m - f_mb.toarray()).max() / np.abs(f_mb.toarray()).max())
    d_ms = float(np.abs(Ms_m - f_ms.toarray()).max() / np.abs(f_ms.toarray()).max())
    MESHS[tag] = (c, q, e, b, f_mb, f_ms)  # uso l'assembly verificato (fj)
    log('assembly_xcheck', f'{tag} Mb vs fj.assemble_fixed', f'max_rel_diff={d_mb:.3e}',
        d_mb < 1e-12)
    log('assembly_xcheck', f'{tag} Ms vs fj.assemble_fixed', f'max_rel_diff={d_ms:.3e}',
        d_ms < 1e-12)

# matrici n128 dalla cache verificata (sha256 documentato in T1-B)
d128 = np.load(T1B / 'matrices_n128x64.npz')
from scipy.sparse import coo_matrix  # noqa: E402
N128 = MESHS['n128'][0].shape[0]
M128 = {}
for pfx in ('Mb', 'Kb', 'Ms', 'Ks'):
    M128[pfx] = coo_matrix((d128[pfx + '_data'],
                            (d128[pfx + '_row'], d128[pfx + '_col'])),
                           shape=(N128, N128)).tocsc()
if 'coords' in d128.files:
    d_c128 = float(np.abs(d128['coords'] - MESHS['n128'][0]).max())
    log('quadrature_n128', 'coords npz == mesh canonica (mio builder)',
        f'max_diff={d_c128:.2e}', d_c128 < 1e-14)
ones128 = np.ones(N128)
inv = [
    ('Mb1=AREA', abs(float(ones128 @ (M128['Mb'] @ ones128)) - AREA) < 1e-12,
     f'{float(ones128 @ (M128["Mb"] @ ones128)):.15g}'),
    ('Ms1=PERIM', abs(float(ones128 @ (M128['Ms'] @ ones128)) - PERIM) < 1e-12,
     f'{float(ones128 @ (M128["Ms"] @ ones128)):.15g}'),
    ('rowsum Kb=0', float(np.abs(M128['Kb'] @ ones128).max()) < 1e-12,
     f'{float(np.abs(M128["Kb"] @ ones128).max()):.3e}'),
    ('rowsum Ks=0', float(np.abs(M128['Ks'] @ ones128).max()) < 1e-12,
     f'{float(np.abs(M128["Ks"] @ ones128).max()):.3e}'),
]
for name, ok, val in inv:
    log('quadrature_n128', name, val, ok)

# ---------------------------------------------------------------- references --
print('== references: caricamento + integrità', flush=True)


refs = {}
REF_FILES = {
    'n16': FIX / 'tolcheck_n16x8_rtol2e-12.npz',
    'n32': FIX / 'reference_n32x16_dt0p00025.npz',
    'n64': FIX / 'reference_n64x32_dt0p00025.npz',
    'n128': T1B / 'reference_n128x64_dt0p00025.npz',
}
for tag, path in REF_FILES.items():
    d = np.load(path)
    U = d['U']
    V = d['V'] if 'V' in d.files else None
    t = d['times']
    n = U.shape[1]
    c = MESHS[tag][0]
    assert U.shape[0] == NSTEP + 1 and len(t) == NSTEP + 1, f'{tag} shape'
    d_t = float(np.abs(t - np.arange(NSTEP + 1) * DT).max())
    log('ref_integrita', f'{tag} tempi', f'N={n} max|t-canonical|={d_t:.2e}',
        d_t < 1e-15)
    # IC (convenzione x-fasta)
    X, Y = c.T
    u0 = 0.2 + 0.1 * np.cos(np.pi * X / LX) * np.cos(np.pi * Y / LY)
    d_ic = float(np.abs(U[0] - u0).max())
    log('ref_integrita', f'{tag} IC esatta', f'max|U0-IC|={d_ic:.2e}', d_ic < 1e-14)
    log('ref_integrita', f'{tag} finitezza U', 'all_finite='
        f'{bool(np.isfinite(U).all())}', bool(np.isfinite(U).all()))
    if V is not None:
        log('ref_integrita', f'{tag} V0=0 e finitezza V',
            f'max|V0|={float(np.abs(V[0]).max()):.2e} finite={bool(np.isfinite(V).all())}',
            float(np.abs(V[0]).max()) < 1e-14 and bool(np.isfinite(V).all()))
    refs[tag] = (U, V)
    for k in ('nfev', 'rtol'):
        if k in d.files:
            print(f'    {tag} {k} = {d[k]}', flush=True)

U16, _ = refs['n16']
U32, _ = refs['n32']
U64, _ = refs['n64']
U128, V128 = refs['n128']

# mappe di restrizione (coincidenza esatta)
I8_16 = restrict(MESHS['n8'][0], MESHS['n16'][0], 'n8->n16')
I16_32 = restrict(MESHS['n16'][0], MESHS['n32'][0], 'n16->n32')
I32_64 = restrict(MESHS['n32'][0], MESHS['n64'][0], 'n32->n64')
I64_128 = restrict(MESHS['n64'][0], MESHS['n128'][0], 'n64->n128')
I8_128 = I64_128[I32_64[I16_32[I8_16]]]
I16_128 = I64_128[I32_64[I16_32]]
I32_128 = I64_128[I32_64]
log('restrizione', 'mappe concatenate n8->...->n128', 'costruite (d2<1e-26)', True)

# ------------------------------------------------------- OOFEM: passi da ----
# leggere: 39 e 80 (T) per n8/n16/n32 + la griglia di 10 istanti di T1-B
TEN_STEPS = [1, 10, 20, 30, 39, 40, 50, 60, 70, 80]
print('== OOFEM: lettura VTU (passi 39, 80 e griglia 10 istanti)', flush=True)
UOOF = {}  # (tag, step) -> u ordinato sulla mesh canonica
for tag, folder in (('n8', F8), ('n16', F16), ('n32', F32A)):
    c = MESHS[tag][0]
    for step in sorted(set(TEN_STEPS) | {80}):
        UOOF[(tag, step)] = read_oofem_step(folder, step, c, f'OOFEM {tag}')

print(f'setup completato in {time.time() - t_start:.1f} s', flush=True)

# ================================================================= S2: @T =====
print('== S2: metriche a T ricalcolate vs t1b_metrics_T.csv', flush=True)
MB = {t: MESHS[t][4] for t in ('n8', 'n16', 'n32', 'n64')}
MS = {t: MESHS[t][5] for t in ('n8', 'n16', 'n32', 'n64')}
BN = {t: MESHS[t][3] for t in ('n8', 'n16', 'n32', 'n64')}
C = {t: MESHS[t][0] for t in ('n8', 'n16', 'n32', 'n64')}

I8_64 = I32_64[I16_32[I8_16]]
I16_64 = I32_64[I16_32]
u8, u16, u32 = UOOF[('n8', S)], UOOF[('n16', S)], UOOF[('n32', S)]
E = {}
for tag, u in (('n8', u8), ('n16', u16), ('n32', u32)):
    i128 = {'n8': I8_128, 'n16': I16_128, 'n32': I32_128}[tag]
    i64 = {'n8': I8_64, 'n16': I16_64, 'n32': I32_64}[tag]
    E[f'{tag}_128'] = metrics4(MB[tag], MS[tag], BN[tag], u - U128[S, i128])
    E[f'{tag}_64'] = metrics4(MB[tag], MS[tag], BN[tag], u - U64[S, i64])
E['fam32_64'] = metrics4(MB['n32'], MS['n32'], BN['n32'], U32[S] - U64[S, I32_64])
E['fam64_128'] = metrics4(MB['n64'], MS['n64'], BN['n64'], U64[S] - U128[S, I64_128])

# valori archiviati
arch = {}
with (T1B / 'results' / 't1b_metrics_T.csv').open() as f:
    r = csv.reader(f)
    next(r)
    for row in r:
        label = row[0]
        for col, key in ((1, 'vs_proxy128'), (2, 'vs_proxy64'),
                         (3, 'family_level_diff')):
            if row[col].strip():
                arch[(label, key)] = float(row[col])

MT = ['bulk_fe_rms', 'bulk_max', 'surf_fe_rms', 'surf_max']
rows = []


def cmp_row(group, label, key, mine):
    a = arch.get((label, key))
    if a is None:
        rows.append([group, label, f'{mine:.6e}', '', '', 'ARCH-MISSING'])
        return
    rd = abs(mine - a) / a
    rows.append([group, label, f'{mine:.6e}', f'{a:.6e}', f'{rd:.3e}',
                 'OK' if rd < 5e-5 else 'CHECK'])
    log('metriche_T', f'{label} ({key})', f'mine={mine:.6e} arch={a:.6e} '
        f'rel_diff={rd:.2e}', rd < 5e-5)


for tag in ('n8', 'n16', 'n32'):
    for k, m in enumerate(MT):
        cmp_row('OOFEM_vs_proxy128', f'{tag} {m} @T', 'vs_proxy128', E[f'{tag}_128'][k])
        cmp_row('OOFEM_vs_proxy64', f'{tag} {m} @T', 'vs_proxy64', E[f'{tag}_64'][k])
for fam in ('fam32_64', 'fam64_128'):
    label = 'family 32->64' if fam == 'fam32_64' else 'family 64->128'
    for k, m in enumerate(MT):
        cmp_row('family_diff', f'{label} {m} @T', 'family_level_diff', E[fam][k])

# ordini a T
orders = {}
for k, m in enumerate(MT):
    orders[f'8->16 {m} @T (vs proxy128)'] = float(np.log2(E['n8_128'][k] / E['n16_128'][k]))
    orders[f'16->32 {m} @T (vs proxy128)'] = float(np.log2(E['n16_128'][k] / E['n32_128'][k]))
    orders[f'family 32->64 {m} @T'] = float(np.log2(E['fam32_64'][k] / E['fam64_128'][k]))
arch_ord = {}
with (T1B / 'results' / 't1b_family_orders_T.csv').open() as f:
    r = csv.reader(f)
    next(r)
    for row in r:
        arch_ord[row[0]] = float(row[1])
ord_rows = []
for k, v in orders.items():
    a = arch_ord.get(k)
    rd = abs(v - a) if a is not None else float('nan')
    ord_rows.append([k, f'{v:.4f}', f'{a:.4f}' if a else '', f'{rd:.3e}',
                     'OK' if rd < 5e-4 else 'CHECK'])
    log('ordini_T', k, f'mine={v:.4f} arch={a:.4f}', rd < 5e-4)
write_csv('metrics_T_recomputed.csv',
          ['group', 'quantity', 'mine', 'archived', 'rel_diff', 'status'], rows)
write_csv('orders_T_recomputed.csv',
          ['order', 'mine', 'archived', 'abs_diff', 'status'], ord_rows)

# ================================================================= S3: bound ==
print('== S3: stabilità proxy + bound di triangolo', flush=True)
ps_rows = []
for tag, i128, i64 in (('n8x4', I8_128, I8_64),
                       ('n16x8', I16_128, I16_64),
                       ('n32x16', I32_128, I32_64)):
    tt = 'n8' if tag == 'n8x4' else 'n16' if tag == 'n16x8' else 'n32'
    e128v = E[f'{tt}_128']
    e64v = E[f'{tt}_64']
    bound = U128[S, i128] - U64[S, i64]
    b_bulk = rms_fe(MB[tt], bound, AREA)
    b_surf = rms_fe(MS[tt], bound, PERIM)
    # bound per metric con la NORMA GIUSTA: per i max serve la norma max
    # della differenza di famiglia SULLO STESSO insieme di nodi
    # (interiori per bulk_max, bordo per surf_max): e128-e64 = U64|_k - U128|_k
    bnd_tt = BN[tt]
    b_bulk_max = float(np.abs(bound[~bnd_tt]).max())
    b_surf_max = float(np.abs(bound[bnd_tt]).max())
    bmap = {'bulk_fe_rms': b_bulk, 'bulk_max': b_bulk_max,
            'surf_fe_rms': b_surf, 'surf_max': b_surf_max}
    for k, m in enumerate(MT):
        ad = abs(e128v[k] - e64v[k])
        b = bmap[m]
        ps_rows.append([tag, m, f'{b:.6e}', f'{ad:.6e}',
                        'SI' if b >= ad else 'NO(!)', f'{b / ad:.2f}' if ad else ''])
    log('proxy_stability', f'{tag} bound bulk FE', f'bound={b_bulk:.6e} '
        f'abs_diff={abs(e128v[0] - e64v[0]):.6e}', b_bulk >= abs(e128v[0] - e64v[0]))
    log('proxy_stability', f'{tag} bound surf FE', f'bound={b_surf:.6e} '
        f'abs_diff={abs(e128v[2] - e64v[2]):.6e}', b_surf >= abs(e128v[2] - e64v[2]))
    log('proxy_stability', f'{tag} bound bulk max', f'bound={b_bulk_max:.6e} '
        f'abs_diff={abs(e128v[1] - e64v[1]):.6e}',
        b_bulk_max >= abs(e128v[1] - e64v[1]))
    log('proxy_stability', f'{tag} bound surf max', f'bound={b_surf_max:.6e} '
        f'abs_diff={abs(e128v[3] - e64v[3]):.6e}',
        b_surf_max >= abs(e128v[3] - e64v[3]))
write_csv('proxy_stability_bound.csv',
          ['mesh', 'metric', 'triangle_bound_mine', 'abs_diff_mine',
           'bound>=abs_diff', 'ratio'], ps_rows)

# ========================================================== S4: zone/max mov. ==
print('== S4: zone + massimo mobile + fit', flush=True)


def zones(c):
    x, y = c[:, 0], c[:, 1]
    rA = np.minimum(np.hypot(x, y - 1.0), np.hypot(x - 2.0, y))
    rB = np.minimum(np.hypot(x, y), np.hypot(x - 2.0, y - 1.0))
    zA = rA < 0.5
    return zA, (rB < 0.5) & ~zA, ~(zA | (rB < 0.5))


zrows = []
arch_z = {}
with (T1B / 'results' / 't1b_zones_T.csv').open() as f:
    for rrow in csv.DictReader(f):
        arch_z[(rrow['dataset'], rrow['zone'])] = rrow
for tag, c, e in (('family 64->128 @T (n64 nodes)', C['n64'], U64[S] - U128[S, I64_128]),
                  ('OOFEM n32x16 vs proxy128 @T', C['n32'], u32 - U128[S, I32_128])):
    zA, zB, zR = zones(c)
    tot = float(np.sum(e ** 2))
    for zname, z in (('cornerA (0,1)/(2,0)', zA), ('cornerB (0,0)/(2,1)', zB),
                     ('remainder', zR)):
        mine = float(np.sum(e[z] ** 2) / tot)
        a = arch_z.get((tag, zname), {}).get('energy_frac')
        af = float(a) if a else None
        rd = abs(mine - af) / af if af else float('nan')
        zrows.append([tag, zname, f'{mine:.6e}', f'{af:.6e}' if af else '',
                      f'{rd:.3e}', 'OK' if rd < 5e-5 else 'CHECK'])
        log('zone', f'{tag} {zname}', f'mine={mine:.6e} arch={af}', rd < 5e-5)
write_csv('zones_recomputed.csv',
          ['dataset', 'zone', 'energy_frac_mine', 'archived', 'rel_diff',
           'status'], zrows)

# massimo mobile
def first_interior_corner_node(c):
    x, y = c[:, 0], c[:, 1]
    inter = ~((x == 0) | (x == 2.0) | (y == 0) | (y == 1.0))
    r = np.minimum(np.hypot(x, y - 1.0), np.hypot(x - 2.0, y))
    r[~inter] = np.inf
    j = int(np.argmin(r))
    return j, float(r[j])


mv = []
for tag, c, bnd, err in (('n8x4 OOFEM', C['n8'], BN['n8'], u8 - U128[S, I8_128]),
                         ('n16x8 OOFEM', C['n16'], BN['n16'], u16 - U128[S, I16_128]),
                         ('n32x16 OOFEM-new', C['n32'], BN['n32'],
                          u32 - U128[S, I32_128]),
                         ('n64 family diff', C['n64'], BN['n64'],
                          U64[S] - U128[S, I64_128])):
    j, rj = first_interior_corner_node(c)
    mv.append({'mesh': tag, 'x': float(c[j, 0]), 'y': float(c[j, 1]),
               'r_c': rj, 'err': float(err[j]), 'abs_err': float(abs(err[j]))})
rs4 = np.array([m['r_c'] for m in mv])
es4 = np.array([m['abs_err'] for m in mv])
p_fit4 = float(np.polyfit(np.log(rs4), np.log(es4), 1)[0])
# fit SU 3 PUNTI SOLO OOFEM (punti 1-3) — decomposizione metodologica
p_fit3 = float(np.polyfit(np.log(rs4[:3]), np.log(es4[:3]), 1)[0])
mv_rows = [[m['mesh'], f'{m["x"]:.6f}', f'{m["y"]:.6f}', f'{m["r_c"]:.6e}',
            f'{m["err"]:.10e}', f'{m["abs_err"]:.10e}'] for m in mv]
mv_rows.append(['fit 4pt (misto: 3 OOFEM + 1 family diff)', '', '', '',
                f'p={p_fit4:.4f}', ''])
mv_rows.append(['fit 3pt (solo OOFEM n8/n16/n32)', '', '', '',
                f'p={p_fit3:.4f}', ''])
write_csv('moving_max_recomputed.csv',
          ['mesh', 'x', 'y', 'r_c', 'err', 'note'], mv_rows)
log('moving_max', 'fit 4pt (serie mista, come T1-B)', f'p={p_fit4:.4f}', True)
log('moving_max', 'fit 3pt (solo OOFEM, stesso punto fisico per mesh)',
    f'p={p_fit3:.4f}', True)

# nodo fisico FISSO (1.75, 0.25) — presente in n8, n16, n32 (e n64, n128)
fn = np.array([1.75, 0.25])
fix_rows = []
for tag, u, i128 in (('n8', u8, I8_128), ('n16', u16, I16_128),
                     ('n32', u32, I32_128)):
    j = np.argmin(np.sum((C[tag] - fn) ** 2, axis=1))
    err = float(u[j] - U128[S, i128[j]])
    fix_rows.append([f'OOFEM {tag}', f'{err:.6e}'])
# differenze di famiglia allo stesso punto
i32f = np.argmin(np.sum((C['n32'] - fn) ** 2, axis=1))
i64f = I32_64[i32f]
i128f = I64_128[i64f]
d32_64_fn = float(abs(U32[S, i32f] - U64[S, i64f]))
d64_128_fn = float(abs(U64[S, i64f] - U128[S, i128f]))
i16f = int(np.argmin(np.sum((C['n16'] - fn) ** 2, axis=1)))
assert np.allclose(C['n16'][i16f], fn, atol=1e-14)
d16_32_fn = float(abs(U16[S, i16f] - U32[S, i32f]))
fix_rows.append(['fam |U32-U64| @fn', f'{d32_64_fn:.6e}'])
fix_rows.append(['fam |U64-U128| @fn', f'{d64_128_fn:.6e}'])
fix_rows.append(['fam |U16-U32| @fn', f'{d16_32_fn:.6e}'])
fix_rows.append(['Richardson 3-livello @fn (T)',
                 f'{np.log2(d16_32_fn / d64_128_fn) / 2:.4f}'])
write_csv('fixed_node_175_025_T.csv', ['object', 'value'], fix_rows)
for r in fix_rows:
    log('fissocorner', f'nodo (1.75,0.25) {r[0]}', r[1], True)

# ================================================================= S5: local ==
print('== S5: Richardson puntuale (nodi fissi, 3-4 soluzioni)', flush=True)
# nodi fissi: P* sulla griglia n16 (presenti in n16..n128), Q* solo n32 (n32..n128)
PTS = [
    ('P1 (1.875,0.125)', 1.875, 0.125, 'n16'),
    ('P2 (1.75,0.25)', 1.75, 0.25, 'n16'),
    ('P3 (1.625,0.375)', 1.625, 0.375, 'n16'),
    ('P4 (1.5,0.25)', 1.5, 0.25, 'n16'),
    ('P5 (1.25,0.25)', 1.25, 0.25, 'n16'),
    ('P6 (1.0,0.25)', 1.0, 0.25, 'n16'),
    ('P7 (1.0,0.5)', 1.0, 0.5, 'n16'),
    ('P8 (0.5,0.5)', 0.5, 0.5, 'n16'),
    ('Q1 (1.9375,0.0625)', 1.9375, 0.0625, 'n32'),
    ('Q2 (1.8125,0.1875)', 1.8125, 0.1875, 'n32'),
]
STEP5 = [10, 20, 40, 60, 80]
ric = []
for pname, px, py, level in PTS:
    j32 = int(np.argmin(np.sum((C['n32'] - [px, py]) ** 2, axis=1)))
    assert np.allclose(C['n32'][j32], [px, py], atol=1e-14)
    j64 = int(I32_64[j32])
    j128 = int(I64_128[j64])
    j16 = None
    if level == 'n16':
        j16 = int(np.argmin(np.sum((C['n16'] - [px, py]) ** 2, axis=1)))
        assert np.allclose(C['n16'][j16], [px, py], atol=1e-14)
    rc = float(min(np.hypot(px - 2.0, py), np.hypot(px, py - 1.0)))
    for s in STEP5:
        d3264 = float(abs(U32[s, j32] - U64[s, j64]))
        d64128 = float(abs(U64[s, j64] - U128[s, j128]))
        if j16 is not None:
            d1632 = float(abs(U16[s, j16] - U32[s, j32]))
            p_low = float(np.log2(d1632 / d3264))
            p2 = float(np.log2(d1632 / d64128) / 2)
        else:
            d1632, p_low, p2 = float('nan'), float('nan'), float('nan')
        p_high = float(np.log2(d3264 / d64128))
        ric.append([pname, f'{rc:.6f}', level, s, f'{s * DT:.6f}',
                    f'{d1632:.6e}' if j16 is not None else '',
                    f'{d3264:.6e}', f'{d64128:.6e}',
                    f'{p_low:.4f}' if j16 is not None else '',
                    f'{p_high:.4f}', f'{p2:.4f}' if j16 is not None else ''])
write_csv('fixed_node_richardson.csv',
          ['point', 'r_c', 'level', 'step', 't', 'D16_32', 'D32_64', 'D64_128',
           'p_low=log2(D16_32/D32_64)', 'p_high=log2(D32_64/D64_128)',
           'p_2step=log2(D16_32/D64_128)/2'], ric)
# digest per punto a T
for pname, px, py, level in PTS:
    j32 = int(np.argmin(np.sum((C['n32'] - [px, py]) ** 2, axis=1)))
    j64 = int(I32_64[j32]); j128 = int(I64_128[j64])
    d3264 = abs(U32[S, j32] - U64[S, j64]); d64128 = abs(U64[S, j64] - U128[S, j128])
    line = (f'  {pname}: D32_64={d3264:.2e} D64_128={d64128:.2e} '
            f'p_high={np.log2(d3264 / d64128):.2f}')
    if level == 'n16':
        j16 = int(np.argmin(np.sum((C['n16'] - [px, py]) ** 2, axis=1)))
        d1632 = abs(U16[S, j16] - U32[S, j32])
        line += (f' D16_32={d1632:.2e} p_low={np.log2(d1632 / d3264):.2f} '
                 f'p_2step={np.log2(d1632 / d64128) / 2:.2f}')
    print(line, flush=True)

# ================================================================= S6: floor ==
print('== S6: floor solutore (OOFEM vs DOP853 stessa mesh)', flush=True)
floors = [
    ('OOFEM n16x8 vs DOP853 n16x8', 'n16', 39, '2.0995e-8 (static audit §1.5)'),
    ('OOFEM n16x8 vs DOP853 n16x8', 'n16', 80, '2.4441e-8 (static audit §1.5)'),
    ('OOFEM n32x16 (T1-A) vs DOP853 n32x16', 'n32', 39, '2.5972e-8 (static audit §1.5)'),
    ('OOFEM n32x16 (T1-A) vs DOP853 n32x16', 'n32', 80, 'NON ARCHIVIATO (nuovo)'),
]
fl_rows = []
for label, tag, step, archnote in floors:
    u = UOOF[(tag, step)]
    U = (U16 if tag == 'n16' else U32)
    e = u - U[step]
    m = metrics4(MB[tag], MS[tag], BN[tag], e)
    for k, mm in enumerate(MT):
        fl_rows.append([label, f't={step * DT:.5f}', mm, f'{m[k]:.6e}',
                        archnote if k == 0 else ''])
    log('floor', f'{label} step {step} bulk_fe', f'{m[0]:.6e}', True)
    print(f'    {label} step {step}: bulk_fe={m[0]:.3e} bulk_max={m[1]:.3e} '
          f'surf_fe={m[2]:.3e} surf_max={m[3]:.3e}', flush=True)
# provenienza: U16 tolcheck DEVE essere il riferimento CORRETTO — il floor
# lo dimostra (se fosse il riferimento buggy la discrepanza sarebbe O(1e-3)).
u16_80 = UOOF[('n16', 80)]
dbug = float(np.abs(u16_80 - np.load(F16 / 'independent_reference.npz')['U'][-1]).max())
log('provenienza', 'OOFEM n16 vs riferimento BUGGY originale (max nodale)',
    f'{dbug:.3e} (atteso O(1e-3), NON usato dai calcoli T1-B)', dbug > 1e-4)
write_csv('floors_recomputed.csv', ['case', 't', 'metric', 'mine', 'archived'], fl_rows)

# ======================================================= S7: serie temporali ==
print('== S7: differenze di famiglia su 80 step (+ n16 nuovo)', flush=True)
ts_rows = []
arch_ts = {}
with (T1B / 'results' / 't1b_orders_time_64_128.csv').open() as f:
    for rrow in csv.DictReader(f):
        arch_ts[int(rrow['step'])] = rrow
D1632 = np.empty(NSTEP)
D3264 = np.empty(NSTEP)
D64128 = np.empty(NSTEP)
for s in range(1, NSTEP + 1):
    d16 = rms_fe(MB['n16'], U16[s] - U32[s, I16_32], AREA)
    d32 = rms_fe(MB['n32'], U32[s] - U64[s, I32_64], AREA)
    d64 = rms_fe(MB['n64'], U64[s] - U128[s, I64_128], AREA)
    D1632[s - 1], D3264[s - 1], D64128[s - 1] = d16, d32, d64
    a = arch_ts[s]
    ok32 = abs(d32 - float(a['fam32_64_bulk_fe_rms'])) / float(
        a['fam32_64_bulk_fe_rms']) < 5e-5
    ok64 = abs(d64 - float(a['fam64_128_bulk_fe_rms'])) / float(
        a['fam64_128_bulk_fe_rms']) < 5e-5
    p_hi = float(np.log2(d32 / d64))
    okp = abs(p_hi - float(a['order_fam_bulk_fe_rms'])) < 5e-4
    ts_rows.append([s, f'{s * DT:.6f}', f'{d16:.6e}', f'{d32:.6e}', f'{d64:.6e}',
                    f'{np.log2(d16 / d32):.4f}', f'{p_hi:.4f}',
                    a['fam32_64_bulk_fe_rms'], a['fam64_128_bulk_fe_rms'],
                    a['order_fam_bulk_fe_rms'],
                    'OK' if (ok32 and ok64 and okp) else 'CHECK'])
    if not (ok32 and ok64 and okp):
        log('serie_tempo', f'step {s}', 'mismatch vs archivio', False)
write_csv('family_time_series_recomputed.csv',
          ['step', 't', 'D16_32_mine', 'D32_64_mine', 'D64_128_mine',
           'p_low_mine', 'p_high_mine', 'D32_64_arch', 'D64_128_arch',
           'p_high_arch', 'status'], ts_rows)
log('serie_tempo', '80/80 step coerenti con t1b_orders_time_64_128.csv',
    'D32_64, D64_128, p_high', all(r[-1] == 'OK' for r in ts_rows))

i_pk = int(np.argmax(D64128))
p_lo = np.log2(D1632 / D3264)
p_hi = np.log2(D3264 / D64128)
both_in = int(np.sum((p_lo >= 1.5) & (p_lo <= 2.5) & (p_hi >= 1.5) & (p_hi <= 2.5)))
close = int(np.sum(np.abs(p_lo - p_hi) < 0.25))
summ = [
    ['D32_64(s): monotono crescente?', f'{bool(np.all(np.diff(D3264) > 0))}',
     f'{D3264[0]:.3e} -> {D3264[-1]:.3e}'],
    ['D64_128(s): monotono?', f'{bool(np.all(np.diff(D64128) > 0))}',
     f'picco step {i_pk + 1} (t={(i_pk + 1) * DT:.5f}) = {D64128[i_pk]:.3e}; '
     f'valore a T = {D64128[-1]:.3e}; D(T)/D(picco) = {D64128[-1] / D64128[i_pk]:.3f}'],
    ['p_low(t) range', f'{p_lo.min():.3f} .. {p_lo.max():.3f}',
     f'a T: {p_lo[-1]:.3f}'],
    ['p_high(t) range', f'{p_hi.min():.3f} .. {p_hi.max():.3f}',
     f'a T: {p_hi[-1]:.3f}'],
    ['step con p_low e p_high entrambi in [1.5,2.5]', str(both_in),
     'di 80'],
    ['step con |p_low - p_high| < 0.25', str(close), 'di 80'],
]
write_csv('family_time_summary.csv', ['quantity', 'value', 'note'], summ)
for r in summ:
    log('serie_tempo', r[0], f'{r[1]} {r[2]}', True)

# ================================================================= S8: energia
print('== S8: energia n128 (81 step, quadrature indipendenti)', flush=True)
c128, q128, e128, b128 = MESHS['n128']
EN = []
for s in range(NSTEP + 1):
    EN.append(my_energy(U128[s], V128[s], c128, q128, e128))
EN = np.array(EN)
arch_en = np.loadtxt(T1B / 'results' / 't1b_energy_n128.csv', delimiter=',',
                     skiprows=1)
de = np.diff(EN[:, 0])
mono = bool(np.all(de < 0))
xchk = []
for s in (0, 40, 80):
    ef = fj.energy_fixed(U128[s], V128[s], c128, q128, e128)
    xchk.append(max(abs(EN[s, k] - ef[k]) for k in range(4)))
# arch_en colonne: 0=time, 1=total, 2=bulk_kinetic, 3=bulk_potential, 4=surface
en_rows = [[f'{s * DT:.6f}', f'{EN[s, 0]:.10e}', f'{arch_en[s, 1]:.10e}',
            f'{abs(EN[s, 0] - arch_en[s, 1]):.3e}', f'{EN[s, 1]:.6e}',
            f'{arch_en[s, 2]:.6e}'] for s in range(NSTEP + 1)]
write_csv('energy_recomputed.csv',
          ['t', 'E_tot_mine', 'E_tot_arch', 'abs_diff_total',
           'kinetic_mine', 'kinetic_arch'], en_rows)
es_rows = [
    ['E0', f'{EN[0, 0]:.10e}'], ['ET', f'{EN[-1, 0]:.10e}'],
    ['variazione relativa totale', f'{(EN[-1, 0] / EN[0, 0] - 1) * 100:.5f}%'],
    ['monotona decrescente (80 passi)', str(mono),
     f'max dE_consecutivo = {de.max():.3e}'],
    ['cinetica t=0 -> T', f'{EN[0, 1]:.3e} -> {EN[-1, 1]:.10e}'],
    ['max|E_mine - E_arch| su 81 step',
     f'{float(np.abs(EN[:, 0] - arch_en[:, 1]).max()):.3e}'],
    ['max|E_mine - fj.energy_fixed| (s=0,40,80)', f'{max(xchk):.3e}'],
]
write_csv('energy_summary.csv', ['quantity', 'value'], es_rows)
for r in es_rows:
    log('energia', r[0], r[1], mono if r[0].startswith('monotona') else True)

# ================================================================= S9: nodale =
print('== S9: decomposizione nodale n32 a T (floor + famiglia)', flush=True)
floor_n = u32 - U32[S]
fam_n = U32[S] - U128[S, I32_128]
e_n = u32 - U128[S, I32_128]
assert float(np.abs(e_n - (floor_n + fam_n)).max()) == 0.0
inter = ~BN['n32']
corr = float(np.corrcoef(np.abs(e_n[inter]), np.abs(fam_n[inter]))[0, 1])
nd_rows = [
    ['max|floor_nodale| = max|u32 - U32|', f'{float(np.abs(floor_n).max()):.3e}'],
    ['frazione nodi interni con |floor| < 1e-7',
     f'{float(np.mean(np.abs(floor_n[inter]) < 1e-7)):.6f}'],
    ['correlazione Pearson |e| vs |famiglia| (nodi interni)', f'{corr:.6f}'],
    ['nodi interni con |e - fam| > 1e-7',
     str(int(np.sum(np.abs(e_n[inter] - fam_n[inter]) > 1e-7)))],
    ['argmax|e| nodo',
     f'{C["n32"][int(np.argmax(np.abs(e_n)))].tolist()} '
     f'valore={float(np.abs(e_n).max()):.6e}'],
    ['argmax|fam| nodo',
     f'{C["n32"][int(np.argmax(np.abs(fam_n)))].tolist()} '
     f'valore={float(np.abs(fam_n).max()):.6e}'],
]
write_csv('node_decomp_n32_T.csv', ['quantity', 'value'], nd_rows)
for r in nd_rows:
    log('nodale_n32', r[0], r[1], True)

# ============================================================== S10: 10 tempi ==
print('== S10: ordini OOFEM a 10 istanti vs t1b_orders_time_oofem_proxy128.csv',
      flush=True)
with (T1B / 'results' / 't1b_orders_time_oofem_proxy128.csv').open() as f:
    reader = csv.DictReader(f)
    cols10 = reader.fieldnames
    arch10 = {int(r['step']): r for r in reader}
s10_rows = []
for step in TEN_STEPS:
    m8 = metrics4(MB['n8'], MS['n8'], BN['n8'], UOOF[('n8', step)] - U128[step, I8_128])
    m16 = metrics4(MB['n16'], MS['n16'], BN['n16'], UOOF[('n16', step)] - U128[step, I16_128])
    m32 = metrics4(MB['n32'], MS['n32'], BN['n32'], UOOF[('n32', step)] - U128[step, I32_128])
    mine = {
        'ord8_16_bulk_fe': np.log2(m8[0] / m16[0]),
        'ord16_32_bulk_fe': np.log2(m16[0] / m32[0]),
        'ord8_16_surf_fe': np.log2(m8[2] / m16[2]),
        'ord16_32_surf_fe': np.log2(m16[2] / m32[2]),
    }
    a = arch10[step]
    for k, v in mine.items():
        ak = [c for c in cols10 if c.startswith(k)]
        if ak:
            av = float(a[ak[0]])
            rd = abs(v - av)
            s10_rows.append([step, f'{step * DT:.5f}', ak[0], f'{v:.4f}',
                             f'{av:.4f}', f'{rd:.3e}',
                             'OK' if rd < 5e-3 else 'CHECK'])
            log('ordini_10t', f'step {step} {ak[0]}', f'mine={v:.4f} arch={av:.4f}',
                rd < 5e-3)
write_csv('oofem_orders_10times_recomputed.csv',
          ['step', 't', 'column', 'mine', 'archived', 'abs_diff', 'status'], s10_rows)

# ============================================================== integrità fin =
write_csv('integrity_checks.csv', ['check', 'detail', 'value', 'status'], CHECKS)
n_fail = sum(1 for c in CHECKS if c[3] == 'FAIL')
print(f'\n=== FINE: {len(CHECKS)} check, FAIL = {n_fail}, '
      f'tempo totale {time.time() - t_start:.1f} s ===', flush=True)
