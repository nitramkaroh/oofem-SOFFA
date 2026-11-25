"""
figcommon.py — shared, strictly READ-ONLY data access layer for the
publication validation figures in this directory.

Data sources (all pre-existing; nothing here is written except data/ and
figure outputs of the make_fig*.py scripts):

  OOFEM VTU exports (field "Temperature" = the shared bulk/surface field u;
  the surface variable of this implementation is the trace u|_Gamma):
    n8x4 : hyperbolic/coupled_independent_space_n8x4_dt0p00025/case.out.m0.NN.vtu
    n16x8: hyperbolic/coupled_independent_space_n16x8_dt0p00025/case.out.m0.NN.vtu
    n32x16: hyperbolic/t1a_n32x16_rerun/oofem/case.out.m0.NN.vtu
            (T1-A rerun — the ONLY valid n32 source: in the original T1
            directory the files for step >= 40 are 0-byte stubs; steps
             1-39 of the rerun are bit-identical to the original, input
             file byte-identical. See hyperbolic/t1a_n32x16_rerun/results/t1a_report.md)

  Semidiscrete DOP853 references (independent integrator, rtol = 2e-12,
  atol = 2e-14; 81 steps, t = k*2.5e-4, k = 0..80; U and V arrays):
    n16  : hyperbolic/coupled_independent_spatial_jacobian_fix/tolcheck_n16x8_rtol2e-12.npz
    n32  : hyperbolic/coupled_independent_spatial_jacobian_fix/reference_n32x16_dt0p00025.npz
    n64  : hyperbolic/coupled_independent_spatial_jacobian_fix/reference_n64x32_dt0p00025.npz
    n128 : hyperbolic/t1b_dop853_n128x64/reference_n128x64_dt0p00025.npz
           (n128 = the T1-B "proxy"; nfev 2117; 21.0 min; RSS ~109 MB)

Conventions (all independently re-verified by the audit of 2026-10-27,
consolidated_final_report/audit_convergence_next/audit_report.md — 172/172
checks OK):
  * tensor mesh on Omega = [0,2] x [0,1], x-fast row-major: node (i,j) = j*(nx+1)+i
  * nested meshes share nodes EXACTLY (max distance 0.0) -> "interpolation"
    between levels is pure index restriction (map below), zero interpolation
    error
  * FE RMS norms (the E(h) used throughout this dossier):
        bulk : sqrt( e^T M_b e / |Omega| ),  |Omega| = 2.0
        surf : sqrt( e^T M_s e / |Gamma| ),  |Gamma| = 6.0
    M_b by 3x3 Gauss quadrature per Q1 quad (exact for the quartic products
    involved), M_s as 1D P1 line elements on the 4 boundary edges
    (local matrix (L/12)[[2,1],[1,2]], exact).
  * The quadratic forms are accumulated element-by-element (no dense
    matrices are built).
"""
from pathlib import Path
import numpy as np

DOSSIER = Path('/home/marco/Documents/Work/oofem-SOFFA/doc/stvowu_validation')
HYPER = DOSSIER / 'hyperbolic'
AUD = DOSSIER / 'consolidated_final_report' / 'audit_convergence_next'
FIX = HYPER / 'coupled_independent_spatial_jacobian_fix'
T1B = HYPER / 't1b_dop853_n128x64'

RUN = {
    'n8':  HYPER / 'coupled_independent_space_n8x4_dt0p00025',
    'n16': HYPER / 'coupled_independent_space_n16x8_dt0p00025',
    'n32': HYPER / 't1a_n32x16_rerun' / 'oofem',
}
REF = {
    'n16':  FIX / 'tolcheck_n16x8_rtol2e-12.npz',
    'n32':  FIX / 'reference_n32x16_dt0p00025.npz',
    'n64':  FIX / 'reference_n64x32_dt0p00025.npz',
    'n128': T1B / 'reference_n128x64_dt0p00025.npz',
}
NX = {'n8': 8, 'n16': 16, 'n32': 32, 'n64': 64, 'n128': 128}
NY = {'n8': 4, 'n16': 8, 'n32': 16, 'n64': 32, 'n128': 64}
LX, LY = 2.0, 1.0
DT, NSTEP, T = 2.5e-4, 80, 0.02
AREA, PERIM = 2.0, 6.0

# ---------------------------------------------------------------- geometry --
def grid(nx, ny):
    """Canonical x-fast tensor grid: X, Y of shape (ny+1, nx+1)."""
    xs = LX * np.arange(nx + 1) / nx
    ys = LY * np.arange(ny + 1) / ny
    return np.meshgrid(xs, ys)  # (ny+1, nx+1) each


def quads(nx, ny):
    """Q1 quads (a, a+1, a+nx+2, a+nx+1) — same pattern as the OOFEM input
    generator run_coupled_independent.py; reproduced exactly by the audit."""
    npr = nx + 1
    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    a = jj * npr + ii  # first node of quad (i, j) = j*(nx+1)+i, row-major
    return np.stack([a, a + 1, a + npr + 1, a + npr], axis=-1).reshape(-1, 4).astype(int)


def boundary_edges(nx, ny):
    """Boundary edges: bottom, right, top (reversed), left (reversed) —
    identical to the audited edge list (audit bug #2 fixed form)."""
    npr = nx + 1
    e = [(i, i + 1) for i in range(nx)]
    e += [(nx + j * npr, nx + (j + 1) * npr) for j in range(ny)]
    e += [(ny * npr + i + 1, ny * npr + i) for i in range(nx - 1, -1, -1)]
    e += [((j + 1) * npr, j * npr) for j in range(ny - 1, -1, -1)]
    return e


def bnd_mask(nx, ny):
    """Boolean mask of boundary nodes, flat (row-major) order."""
    X, Y = grid(nx, ny)
    return ((X == 0) | (X == LX) | (Y == 0) | (Y == LY)).ravel()


def restrict_map(fine_nx, fine_ny, coarse_nx, coarse_ny):
    """Exact nodal restriction map coarse -> fine (uniform tensor grids).
    Verified at call sites: coords_fine[map] == coords_coarse (max diff 0.0)."""
    r = fine_nx // coarse_nx
    assert r * coarse_nx == fine_nx and r * coarse_ny == fine_ny
    ffpr = fine_nx + 1
    II, JJ = np.meshgrid(np.arange(coarse_nx + 1), np.arange(coarse_ny + 1))
    return ((r * JJ) * ffpr + r * II).ravel()


# -------------------------------------------------------------- VTU reader --
def load_vtu(path):
    """Read an OOFEM .vtu (ascii Float64, field 'Temperature').
    Returns (u, pts) in canonical x-fast order; asserts exact coordinate
    coincidence with the canonical tensor grid (max distance 0.0)."""
    raw = Path(path).read_bytes()
    npts = int(raw.split(b'NumberOfPoints="', 1)[1].split(b'"', 1)[0])
    # npts = (nx+1)(ny+1) with nx = 2*ny  =>  npts = (2k+1)(k+1) = 2k^2+3k+1
    #  =>  k = (sqrt(8*npts+1) - 3)/4, exact for all grid sizes used here
    k = int(round((np.sqrt(8 * npts + 1) - 3) / 4))
    ny, nx = k, 2 * k
    assert (nx + 1) * (ny + 1) == npts, f'{path}: not a 2:1 tensor grid'

    i0 = raw.find(b'<Points>')
    i1 = raw.find(b'>', raw.find(b'<DataArray', i0)) + 1
    i2 = raw.find(b'</DataArray>', i1)
    pts = np.fromstring(raw[i1:i2].decode(), sep=' ').reshape(-1, 3)

    j0 = raw.find(b'Name="Temperature"')
    j1 = raw.find(b'>', raw.rfind(b'<DataArray', 0, j0)) + 1
    j2 = raw.find(b'</DataArray>', j1)
    u = np.fromstring(raw[j1:j2].decode(), sep=' ')
    assert len(u) == npts, f'{path}: field length mismatch'

    X, Y = grid(nx, ny)
    key2idx = {(round(float(p[0]), 12), round(float(p[1]), 12)): k
               for k, p in enumerate(pts)}
    pos = np.empty(npts, dtype=int)
    for j in range(ny + 1):
        for i in range(nx + 1):
            pos[j * (nx + 1) + i] = key2idx[(round(float(X[j, i]), 12),
                                             round(float(Y[j, i]), 12))]
    pts_c, u_c = pts[pos], u[pos]
    dmax = float(np.abs(pts_c[:, :2] - np.column_stack([X.ravel(), Y.ravel()])).max())
    assert dmax == 0.0, f'{path}: coordinate mismatch {dmax}'
    return u_c, pts_c


# ------------------------------------------------------------- references --
_REF_CACHE = {}


def ref(name):
    if name not in _REF_CACHE:
        _REF_CACHE[name] = dict(np.load(REF[name]))
    return _REF_CACHE[name]


def ref_at(name, step):
    """U[step] of the DOP853 reference on mesh `name`."""
    return ref(name)['U'][step]


# ------------------------------------------------------------------ norms --
def _qshape(x, y):
    """Q1 shape functions and derivatives on the reference quad (factor .25,
    same as run_coupled_independent.py — audit bug #1 fixed form)."""
    N = np.array([(1 - x) * (1 - y), (1 + x) * (1 - y),
                  (1 + x) * (1 + y), (1 - x) * (1 + y)]) * 0.25
    # dN[0,:] = dN/dxi, dN[1,:] = dN/deta (same values as run_coupled_independent.py,
    # dN = .25*([[-(1-eta),-(1-xi)],[(1-eta),-(1+xi)],[(1+eta),(1+xi)],[-(1+eta),(1-xi)]],
    # transposed here to coordinate-direction rows)
    dN = np.array([[-(1 - y), (1 - y), (1 + y), -(1 + y)],
                   [-(1 - x), -(1 + x), (1 + x), (1 - x)]]) * 0.25
    return N, dN


_GX = np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)])
_GW = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


def bulk_qf(e, nx, ny):
    """e^T M_b e accumulated element-by-element, 3x3 Gauss (exact here)."""
    X, Y = grid(nx, ny)
    xy = np.column_stack([X.ravel(), Y.ravel()])
    q = quads(nx, ny)
    el = e[q]
    x0, y0 = xy[q][:, :, 0], xy[q][:, :, 1]
    s = 0.0
    for a in range(3):
        for b in range(3):
            N, dN = _qshape(_GX[a], _GX[b])  # N (4,); dN[0,:] = dN/dxi, dN[1,:] = dN/deta
            J11 = dN[0, :] @ x0.T
            J12 = dN[0, :] @ y0.T
            J21 = dN[1, :] @ x0.T
            J22 = dN[1, :] @ y0.T
            det = J11 * J22 - J12 * J21
            # MASS matrix: M_ij = w * det * N_i * N_j  (exact with 3x3 Gauss
            # for the quartic integrands of interest; audited assembler used
            # the same form, verified to 0.0 vs check_jacobian_fix)
            Ml = ((N[:, None] * N[None, :])[:, :, None]
                  * det[None, None, :] * _GW[a] * _GW[b])
            s += float(np.einsum('ei,ije,ej->', el, Ml, el))
    return s


def surf_qf(e, nx, ny):
    """e^T M_s e accumulated boundary-edge-by-edge (1D P1, exact)."""
    X, Y = grid(nx, ny)
    xy = np.column_stack([X.ravel(), Y.ravel()])
    s = 0.0
    for a, b in boundary_edges(nx, ny):
        dx, dy = xy[b, 0] - xy[a, 0], xy[b, 1] - xy[a, 1]
        L = np.hypot(dx, dy)
        v = np.array([e[a], e[b]])
        s += float((L / 6.0) * (2.0 * v[0] * v[0] + 2.0 * v[1] * v[1]
                                + 2.0 * v[0] * v[1]))
    return s


def fe_rms_bulk(e, nx, ny):
    return float(np.sqrt(bulk_qf(e, nx, ny) / AREA))


def fe_rms_surf(e, nx, ny):
    return float(np.sqrt(surf_qf(e, nx, ny) / PERIM))


def selfcheck_quadrature():
    """Constant-field invariants: bulk_qf(1) = AREA, surf_qf(1) = PERIM."""
    for nx, ny in ((8, 4), (32, 16), (64, 32)):
        n = (nx + 1) * (ny + 1)
        assert abs(bulk_qf(np.ones(n), nx, ny) - AREA) < 1e-12
        assert abs(surf_qf(np.ones(n), nx, ny) - PERIM) < 1e-12


def oofem_field(tag, step):
    """OOFEM field on mesh `tag` at step `step` (1..80), canonical order."""
    return load_vtu(RUN[tag] / f'case.out.m0.{step}.vtu')[0]


# ----------------------------------------------------------------- archive --
def archive_T_metrics():
    """All stored values of t1b_metrics_T.csv keyed by (label, column).
    Labels are the full CSV row labels, e.g. 'n32x16 bulk_fe_rms @T'."""
    import csv
    arch = {}
    with (T1B / 'results' / 't1b_metrics_T.csv').open() as f:
        r = csv.reader(f)
        next(r)  # header
        for row in r:
            label = row[0]
            for col, key in ((1, 'vs_proxy128'), (2, 'vs_proxy64'),
                             (3, 'family_level_diff')):
                if len(row) > col and row[col].strip():
                    arch[(label, key)] = float(row[col])
    return arch


def archive_floor(case_sub, t, metric):
    """Value from the audit's floors_recomputed.csv ('mine' column).
    case_sub: substring of the case column, e.g. 'n32x16'."""
    import csv
    with (AUD / 'results' / 'floors_recomputed.csv').open() as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            case, tcol, met, mine = row[0], row[1], row[2], row[3]
            if case_sub in case and tcol == f't={t:.5f}' and met == metric:
                return float(mine)
    raise KeyError(f'floor not found: {case_sub} t={t} {metric}')
