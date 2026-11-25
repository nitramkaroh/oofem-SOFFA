#!/usr/bin/env python3
"""T1 — Spatial convergence of the coupled hyperbolic Allen-Cahn system.

Uses ONLY archived OOFEM outputs (coupled_independent_space_n{8x4,16x8,32x16}_dt0p00025)
and the corrected independent DOP853 reference (Jacobian-fixed line assembly).
No OOFEM re-run, no overwrite of archived results: all outputs are written to
this directory only.

Stages (run individually, each cached; `all` runs them in order):
  check     mesh comparability + assembly/Jacobian checks on the fine mesh + IC/params vs case.in
  tolcheck  DOP853 n16x8 at rtol 2e-12 vs 2e-14 (reference temporal-error bound)
  control32 DOP853 corrected n32x16 (direct Newmark temporal-error floor at dt=2.5e-4)
  proxy64   DOP853 corrected n64x32 trajectory (spatial proxy reference)
  tolcheck64 targeted fine-mesh tolerance sensitivity (9 pts, t<=0.002, ~1 min)
  compare   per-mesh errors vs proxy, orders, CSV + PNG
  all       the above

Conventions (identical to the Jacobian-fix verification):
  bulk FE RMS   = sqrt(e' Mb e / (LX*LY))            over all nodes
  surface FE RMS= sqrt(e' Ms e / (2(LX+LY)))         over boundary nodes, corrected Ms
  nodal max bulk / max surface (interior / boundary masks)
  proxy: DOP853 (rtol 2e-12, atol 2e-14) on n64x32, nodal restriction (nested grids)
"""
from __future__ import annotations
import csv, math, re, resource, sys, time, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

HERE = Path(__file__).resolve().parent
HYPERBOLIC = HERE.parent
sys.path.insert(0, str(HYPERBOLIC))
sys.path.insert(0, str(HERE.parent / 'coupled_independent_jacobian_fix'))
import run_coupled_independent as orig            # read-only
import check_jacobian_fix as fj                   # read-only (corrected assembly)

LX, LY, T = orig.LX, orig.LY, orig.T              # 2.0, 1.0, 0.02
DT = 2.5e-4
NSTEP = round(T / DT)                              # 80
TIMES = np.linspace(0, T, NSTEP + 1)
AREA = LX * LY
PERIM = 2 * (LX + LY)
RTOL, ATOL = 2e-12, 2e-14

LOG_ROWS = []


def log(stage, wall, note=''):
    LOG_ROWS.append({'stage': stage, 'wall_s': round(wall, 2), 'note': note})
    print(f'[{stage}] {wall:.2f}s {note}', flush=True)


def write_log():
    with (HERE / 'run_log.csv').open('a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['stage', 'wall_s', 'note'])
        if f.tell() == 0:
            w.writeheader()
        w.writerows(LOG_ROWS)


# ------------------------------------------------------------ stage: check --
def parse_case_in(folder):
    """Material constants, problem line, per-node IC, node coords, line nodes."""
    text = (folder / 'case.in').read_text()
    mat = dict(re.findall(r'(\w+)\s+([-\d.eE]+)',
                          re.search(r'^allencahn 1 (.+)$', text, re.M).group(1)))
    prob = re.search(r'^nltransienttransportproblem (.+)$', text, re.M).group(1)
    ic = {}
    sets = {}
    for m in re.finditer(r'^Set (\d+) nodes 1 (\d+)$', text, re.M):
        sets[int(m.group(1))] = int(m.group(2))          # set -> node (1-based)
    for m in re.finditer(r'^InitialCondition \d+ dofs 1 10 Conditions 1 u ([-\d.eE]+) set (\d+)$',
                         text, re.M):
        node = sets[int(m.group(2))]
        ic[node] = float(m.group(1))
    coords = np.zeros((0, 2))
    for m in re.finditer(r'^node (\d+) coords 3 ([-\d.eE]+) ([-\d.eE]+) ([-\d.eE]+)$',
                         text, re.M):
        coords = np.vstack([coords, [float(m.group(2)), float(m.group(3))]])
    linenodes = set()
    for m in re.finditer(r'^line1ht \d+ nodes 2 (\d+) (\d+) nip 3$', text, re.M):
        linenodes.update([int(m.group(1)), int(m.group(2))])
    return mat, prob, ic, coords, linenodes


def validated_mesh(nx, ny, rows):
    """Regenerate the mesh and cross-validate it against the archived case.in
    (node coordinates and line-element node set). mesh.npz is used when valid;
    the n32x16 archive is a 0-byte file (corrupted write) and is left untouched.
    """
    coords, quads, edges, boundary = orig.mesh(nx, ny)
    folder = HYPERBOLIC / f'coupled_independent_space_n{nx}x{ny}_dt0p00025'
    mat, prob, ic, in_coords, linenodes = parse_case_in(folder)
    ok_c = len(in_coords) == len(coords) and len(coords) > 0 \
        and bool(np.abs(in_coords - coords).max() < 1e-13)
    check_row(rows, f'n{nx}x{ny} regenerated node coords == case.in node coords',
              ok_c, f'max|diff|={float(np.abs(in_coords - coords).max()):.2e}'
              if len(in_coords) == len(coords) else 'no/mismatched nodes')
    in_bnd = np.array(sorted(linenodes)) - 1            # 1-based -> 0-based
    check_row(rows, f'n{nx}x{ny} regenerated boundary == case.in line1ht node set',
              np.array_equal(np.sort(boundary), in_bnd))
    mpath = folder / 'mesh.npz'
    if mpath.stat().st_size == 0:
        check_row(rows, f'n{nx}x{ny} mesh.npz is a 0-byte file (corrupted archive, left untouched; '
                        f'regenerated mesh used after validation vs case.in)',
                  True, '0 bytes')
        return coords, quads, edges, boundary
    mesh = np.load(mpath)
    check_row(rows, f'n{nx}x{ny} mesh.npz == regenerated mesh',
              np.array_equal(mesh['coords'], coords)
              and np.array_equal(mesh['quads'], quads)
              and np.array_equal(mesh['edges'], edges)
              and np.array_equal(mesh['boundary'], boundary))
    return coords, quads, edges, boundary


def check_row(rows, name, ok, value=''):
    rows.append({'check': name, 'pass': bool(ok), 'value': value})
    return ok


def stage_check():
    t0 = time.perf_counter()
    rows = []
    ctx = {}
    # 1) meshes validated against archived case.in (and mesh.npz where valid)
    for nx, ny in [(8, 4), (16, 8), (32, 16)]:
        ctx[(nx, ny)] = validated_mesh(nx, ny, rows)
    ctx[(64, 32)] = orig.mesh(64, 32)
    # 2) nested grids: every coarse node is a fine node (exact)
    cf = ctx[(64, 32)][0]
    for nx, ny in [(8, 4), (16, 8), (32, 16)]:
        c, _, _, _ = ctx[(nx, ny)]
        d = np.linalg.norm(c[:, np.newaxis, :] - cf[np.newaxis, :, :], axis=2)
        check_row(rows, f'n{nx}x{ny} nodes subset of n64x32 nodes (atol 1e-13)',
                  d.min(axis=1).max() < 1e-13, f'max min-dist={d.min(axis=1).max():.2e}')
    # 3) boundary fully covered by line elements (closed dynamic boundary)
    for nx, ny in [(8, 4), (16, 8), (32, 16), (64, 32)]:
        _, _, edges, boundary = ctx[(nx, ny)]
        onb = set(np.unique(edges.ravel()).tolist()) <= set(boundary.tolist())
        check_row(rows, f'n{nx}x{ny} line elements cover the whole boundary '
                        f'({len(edges)} edges = 2({nx}+{ny}), all nodes on boundary)',
                  onb and len(edges) == 2 * (nx + ny))
    # 4) corrected-assembly Jacobian/quadrature checks on the FINE mesh
    c, quads, edges, boundary = ctx[(64, 32)]
    nd = len(c)
    _, _, Ms, Ks = fj.assemble_fixed(64, 32)[4:]
    Mb = lil_matrix((nd, nd)); Kb = lil_matrix((nd, nd))
    for ids in quads:
        xy = c[ids]
        for a, xi in enumerate(orig.GAUSS_X):
            for b, eta in enumerate(orig.GAUSS_X):
                N, dN = orig.qshape(xi, eta)
                J = xy.T @ dN; det = np.linalg.det(J)
                grad = dN @ np.linalg.inv(J)
                w = det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
                Mb[np.ix_(ids, ids)] += orig.RHO * np.outer(N, N) * w
                Kb[np.ix_(ids, ids)] += orig.KB * (grad @ grad.T) * w
    ones = np.ones(nd)
    Mb1, Ms1, Kb1, Ks1 = Mb @ ones, Ms @ ones, Kb @ ones, Ks @ ones
    check_row(rows, 'total Mb 1 = area (bulk mass, 3x3 Gauss)',
              abs(float(ones @ Mb1) - AREA) < 1e-12,
              f'{float(ones @ Mb1):.15g} vs {AREA}')
    check_row(rows, 'total Ms 1 = perimeter (surface mass, corrected J=L/2)',
              abs(float(ones @ Ms1) - PERIM) < 1e-12,
              f'{float(ones @ Ms1):.15g} vs {PERIM}')
    check_row(rows, 'row sums Kb = 0 (bulk diffusion constant invariance)',
              np.abs(Kb1).max() < 1e-12, f'{np.abs(Kb1).max():.3e}')
    check_row(rows, 'row sums Ks = 0 (surface diffusion constant invariance)',
              np.abs(Ks1).max() < 1e-12, f'{np.abs(Ks1).max():.3e}')
    # single horizontal edge: exact P1 matrices (the Jacobian factors themselves)
    e0 = edges[0]
    L = float(np.linalg.norm(c[e0[1]] - c[e0[0]]))
    Ms_e = np.zeros((2, 2)); Ks_e = np.zeros((2, 2))
    for a, s in enumerate(orig.GAUSS_X):
        N = np.array([(1 - s) / 2, (1 + s) / 2])
        dN = np.array([-.5, .5]) * (2.0 / L)
        w = 0.5 * L * orig.GAUSS_W[a]
        Ms_e += np.outer(N, N) * w
        Ks_e += np.outer(dN, dN) * w
    exp_M = (L / 6) * np.array([[2, 1], [1, 2]])
    exp_K = (1 / L) * np.array([[1, -1], [-1, 1]])
    check_row(rows, f'single-edge Ms_e = (L/6)[[2,1],[1,2]] (L={L})',
              np.abs(Ms_e - exp_M).max() < 1e-14, f'max|diff|={np.abs(Ms_e - exp_M).max():.2e}')
    check_row(rows, f'single-edge Ks_e = (1/L)[[1,-1],[-1,1]] (L={L})',
              np.abs(Ks_e - exp_K).max() < 1e-13, f'max|diff|={np.abs(Ks_e - exp_K).max():.2e}')
    # reaction partition of unity (u = 0.2 constant)
    u02 = np.full(nd, 0.2)
    rb, rs = fj.residual_reaction_fixed(u02, c, quads, edges)
    f02 = orig.MB * (0.2 - 0.2 ** 3)
    check_row(rows, 'bulk reaction: sum rb = m(0.2)*area (partition of unity)',
              abs(rb.sum() - f02 * AREA) < 1e-12, f'{float(rb.sum()):.15g} vs {f02 * AREA:.15g}')
    check_row(rows, 'surface reaction: sum rs = m(0.2)*perimeter (corrected weight L/2)',
              abs(rs.sum() - f02 * PERIM) < 1e-12, f'{float(rs.sum()):.15g} vs {f02 * PERIM:.15g}')
    # surface energy at constant u: gradient 0, potential = F(u)*perimeter
    tot, ek, epes, esp = fj.energy_fixed(u02, np.zeros(nd), c, quads, edges)
    exp_s = orig.MS * 0.25 * (0.2 ** 2 - 1) ** 2 * PERIM
    check_row(rows, 'surface energy at u=0.2: esp = F(0.2)*perimeter (grad=0)',
              abs(esp - exp_s) < 1e-12, f'{esp:.15g} vs {exp_s:.15g}')
    # 5) archived inputs: params identical to reference + IC = cosine formula at every node
    expect_mat = {'d': '1.0', 'k': '1.0', 'c': '1.0', 'mobility': '1.0',
                  'surfacebeta': '1.0', 'manufactured': '0',
                  'hyperbolic': '1', 'rho': '1.0', 'surfacecapacity': '1.0'}
    for nx, ny in [(8, 4), (16, 8), (32, 16)]:
        folder = HYPERBOLIC / f'coupled_independent_space_n{nx}x{ny}_dt0p00025'
        mat, prob, ic, _, _ = parse_case_in(folder)
        check_row(rows, f'n{nx}x{ny} case.in material params == reference constants',
                  all(mat.get(k) == v for k, v in expect_mat.items()),
                  str(mat))
        ok_prob = ('nsteps 80 deltat 0.00025' in prob and 'rtol 1.e-8' in prob
                   and 'wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0' in prob)
        check_row(rows, f'n{nx}x{ny} case.in problem line (80 steps, dt 2.5e-4, rtol 1e-8, '
                        f'wave 1, Newmark 1/4 1/2, v0 0)', ok_prob)
        coords, _, _, _ = ctx[(nx, ny)]
        err_ic = max(abs(ic[n + 1] - (0.2 + 0.1 * math.cos(math.pi * coords[n, 0] / 2)
                                     * math.cos(math.pi * coords[n, 1])))
                     for n in range(len(coords)))
        check_row(rows, f'n{nx}x{ny} case.in IC = 0.2+0.1 cos(pi x/2) cos(pi y) at every node',
                  err_ic < 1e-14, f'max|diff|={err_ic:.2e}')
    with (HERE / 'assembly_and_comparability_checks.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['check', 'pass', 'value'])
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if r['pass'] is not True]
    assert not bad, f'checks failed: {bad}'
    log('check', time.perf_counter() - t0, f'{len(rows)} checks, all passed')


# ------------------------------------------------------ stage: DOP853 runs --
def dop853(nx, ny, times, rtol=RTOL, atol=ATOL):
    """Corrected-reference DOP853 with configurable tolerances (same assembly
    and ODE structure as check_jacobian_fix.reference_fixed)."""
    coords, quads, edges, boundary, Mb, Kb, Ms, Ks = fj.assemble_fixed(nx, ny)
    mb_lu = orig.splu(Mb)
    C = orig.ETA * Mb + Ms
    K = Kb + Ks

    def rhs(t, z):
        u = z[:len(coords)]; v = z[len(coords):]
        rb, rs = fj.residual_reaction_fixed(u, coords, quads, edges)
        a = mb_lu.solve(-(C @ v) - (K @ u) + rb + rs)
        return np.r_[v, a]

    y0 = np.r_[orig.initial_u(coords), np.zeros(len(coords))]
    sol = solve_ivp(rhs, (0, float(times[-1])), y0, t_eval=times,
                    method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    U = sol.y[:len(coords), :].T
    assert np.isfinite(U).all(), 'non-finite reference state (NaN/Inf): stopping'
    return (coords, quads, edges, boundary, Mb, Kb, Ms, Ks), U, \
        sol.y[len(coords):, :].T, int(sol.nfev)


def cached_run(name, nx, ny, times, rtol=RTOL, atol=ATOL):
    out = HERE / name
    if out.exists():
        d = np.load(out)
        nfev = int(d['nfev']) if 'nfev' in d.files else 'n/a (pre-diagnostic run)'
        log(name, 0.0, f'cached (nfev={nfev})')
        return d['U'], d['times']
    rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()
    _, U, V, nfev = dop853(nx, ny, times, rtol, atol)
    wall = time.perf_counter() - t0
    rss1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # hard internal wall cap: 5-min budget for a single full attempt
    assert wall < 300, f'full run took {wall:.0f}s (> 300 s budget): stopping, no retry'
    np.savez(out, U=U, V=V, times=times, nx=nx, ny=ny, rtol=rtol, atol=atol, nfev=nfev)
    log(name, wall,
        f'nx={nx} ny={ny} rtol={rtol:g} atol={atol:g} tpoints={len(times)} '
        f'nfev={nfev} (~{nfev // 12} steps of 12 stages) '
        f'rss_peak={rss1 / 1024:.1f}MB (+{rss1 - rss0}kB)')
    return U, times


def stage_tolcheck():
    t0 = time.perf_counter()
    U12, _ = cached_run('tolcheck_n16x8_rtol2e-12.npz', 16, 8, TIMES, 2e-12, 2e-14)
    U14, _ = cached_run('tolcheck_n16x8_rtol2e-14.npz', 16, 8, TIMES, 2e-14, 2e-16)
    d = float(np.max(np.abs(U12 - U14)))
    with (HERE / 'tolerance_check_n16x8.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rtol', 'atol', 'max_abs_diff_vs_tighter'])
        w.writerow(['2e-12', '2e-14', f'{d:.6e}'])
    log('tolcheck', time.perf_counter() - t0, f'max|U(rtol 2e-12) - U(rtol 2e-14)| = {d:.3e}')
    return d


def stage_control32():
    t0 = time.perf_counter()
    cached_run('reference_n32x16_dt0p00025.npz', 32, 16, TIMES)
    log('control32', time.perf_counter() - t0, 'reference trajectory n32x16 stored')


def stage_proxy64():
    # cost probe: 9 time points (t <= 0.002) to bound the full run before committing
    probe_n = 9
    t0 = time.perf_counter()
    _, U, V, nfev_p = dop853(64, 32, TIMES[:probe_n])
    w = time.perf_counter() - t0
    est = w * (NSTEP + 1) / probe_n
    LOG_ROWS.append({'stage': 'proxy64_probe', 'wall_s': round(w, 2),
                     'note': f'est full {est:.0f}s (nfev {nfev_p} for t<=0.002)'})
    print(f'[proxy64] probe {probe_n} pts in {w:.1f}s -> full estimate {est:.0f}s', flush=True)
    assert est < 270, (f'full proxy run estimated {est:.0f}s (> 270 s margin under the '
                        f'300 s budget): stop and decide')
    cached_run('reference_n64x32_dt0p00025.npz', 64, 32, TIMES)


def stage_tolcheck64():
    """Targeted tolerance-sensitivity check on the fine mesh itself (limited
    budget by design: 9 points, t <= 0.002, ~1 min on 1 core) — compares the
    stored rtol 2e-12 reference against a tighter rtol 2e-14 integration.
    A second FULL run would be ~3 min and is not launched automatically.
    """
    t0 = time.perf_counter()
    U14, t14 = cached_run('tolcheck64_partial_rtol2e-14.npz', 64, 32,
                          TIMES[:9], 2e-14, 2e-16)
    ref = np.load(HERE / 'reference_n64x32_dt0p00025.npz')
    U12 = ref['U'][:len(t14)]
    assert np.allclose(t14, ref['times'][:len(t14)], atol=1e-15), 'time grids differ'
    d = float(np.max(np.abs(U12 - U14)))
    with (HERE / 'tolerance_check_n64x32_partial.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['mesh', 'interval', 'rtol_base', 'rtol_tight', 'max_abs_diff', 'note'])
        w.writerow(['n64x32', 't in [0, 0.002] (9 pts)', '2e-12 (atol 2e-14)',
                    '2e-14 (atol 2e-16)', f'{d:.6e}',
                    'scipy clamps rtol 2e-14 to 2.2e-14 (machine epsilon)'])
    log('tolcheck64', time.perf_counter() - t0,
        f'max|U(2e-12) - U(2e-14)| on t<=0.002 = {d:.3e}')
    return d


# --------------------------------------------------------- stage: compare --
def read_step_oofem(folder, step, coords_expect):
    """Archived OOFEM field at step (1-based). Raises a documented error if the
    VTU is a 0-byte file (corrupted archive, left untouched)."""
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    want = step * DT
    ds = next(d for d in pvd.findall('.//DataSet')
              if abs(float(d.attrib['timestep']) - want) < 1e-12)
    path = folder / ds.attrib['file']
    if path.stat().st_size == 0:
        raise FileNotFoundError(f'{path.name} is a 0-byte file (corrupted archive, '
                                f'left untouched)')
    vtu = ET.parse(path).getroot()
    pts = np.fromstring(vtu.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)[:, :2]
    assert np.allclose(pts, coords_expect, rtol=0, atol=1e-13), 'VTU coords != mesh coords'
    arr = next(a for a in vtu.findall('.//PointData/DataArray')
               if a.attrib.get('Name') == 'Temperature')
    u = np.fromstring(arr.text or '', sep=' ')
    assert np.isfinite(u).all(), 'non-finite archived field'
    return u


def fine_to_coarse_map(nc, nf):
    """Index of each coarse node in the fine node array (exact coordinate match)."""
    idx = np.empty(len(nc), dtype=int)
    for i, (x, y) in enumerate(nc):
        d2 = (nf[:, 0] - x) ** 2 + (nf[:, 1] - y) ** 2
        j = int(np.argmin(d2))
        assert d2[j] < 1e-26, f'coarse node {i} ({x}, {y}) has no matching fine node'
        idx[i] = j
    return idx


def stage_compare():
    """Dual-time comparison. The n32x16 archive is corrupted from step 40
    (0-byte VTUs, left untouched): the t=0.02 field is NOT available there.
      interval 8x4 -> 16x8   at t = 0.02    (step 80; both archives complete)
      interval 16x8 -> 32x16 at t = 0.00975 (step 39; last complete n32x16 step)
    Proxy: DOP853 n64x32 trajectory, nodal restriction (nested grids, exact
    coordinate match). Norms: FE RMS with Mb (area) / corrected Ms (perimeter),
    nodal max with interior/boundary masks (same conventions as the
    Jacobian-fix verification).
    """
    t0 = time.perf_counter()
    # cached_run returns (U, times): U first
    U64, _ = cached_run('reference_n64x32_dt0p00025.npz', 64, 32, TIMES)
    U32, _ = cached_run('reference_n32x16_dt0p00025.npz', 32, 16, TIMES)
    fine_coords = orig.mesh(64, 32)[0]
    n32_coords = orig.mesh(32, 16)[0]

    plan = [((8, 4), 80), ((16, 8), 80), ((16, 8), 39), ((32, 16), 39)]
    rows = []
    ctx = {}
    per_node = {}
    for (nx, ny), step in plan:
        key = (nx, ny)
        if key not in ctx:
            folder = HYPERBOLIC / f'coupled_independent_space_n{nx}x{ny}_dt0p00025'
            # mesh.npz is a 0-byte (corrupted) file for n32x16: use the
            # regenerated mesh, validated against case.in
            vrows = []
            coords, _, _, boundary = validated_mesh(nx, ny, vrows)
            bad = [r for r in vrows if r['pass'] is not True]
            assert not bad, f'mesh validation failed for n{nx}x{ny}: {bad}'
            full = fj.assemble_fixed(nx, ny)      # (..., Mb, Kb, Ms, Ks)
            Mb, Ms = full[4], full[6]
            ctx[key] = dict(folder=folder, coords=coords, boundary=boundary,
                            Mb=Mb, Ms=Ms,
                            idx64=fine_to_coarse_map(coords, fine_coords),
                            idx32=fine_to_coarse_map(coords, n32_coords))
        c = ctx[key]
        # TIME CONVENTION (verified): OOFEM pvd DataSet m0.k has timestep k*DT
        # (m0.1 -> 2.5e-4 ... m0.80 -> 0.02), and the reference U[k] is at
        # t = k*DT (TIMES = linspace(0, T, NSTEP+1), U[0] = IC). So OOFEM step
        # s is compared with U[s] at the SAME time t = s*DT, the same
        # convention used by the accepted Jacobian-fix verification.
        # (An earlier draft used U[step-1], a one-dt offset; that made the
        # 'temporal floor' equal to the proxy's own one-step drift and is
        # corrected here. Documented in spatial_verification_report.md.)
        u_oofem = read_step_oofem(c['folder'], step, c['coords'])
        u64c = U64[step][c['idx64']]
        e = u_oofem - u64c
        mask = np.zeros(len(e), dtype=bool); mask[c['boundary']] = True
        temporal = None
        if key == (32, 16):
            e_t = u_oofem - U32[step]      # same mesh AND same time: solver floor
            temporal = (orig.rms_fe(c['Mb'], e_t, AREA), float(np.abs(e_t).max()))
        rows.append({
            'mesh': f'n{nx}x{ny}', 'step': step, 't_compare': f'{step * DT:g}',
            'bulk_fe_rms': orig.rms_fe(c['Mb'], e, AREA),
            'bulk_max': float(np.abs(e[~mask]).max()),
            'surf_fe_rms': orig.rms_fe(c['Ms'], e, PERIM),
            'surf_max': float(np.abs(e[mask]).max()),
            # proxy sanity: level difference between the two proxy meshes at
            # the same comparison time
            'proxy_level_diff_fe': orig.rms_fe(
                c['Mb'], U32[step][c['idx32']] - u64c, AREA),
            'temporal_floor_fe': temporal[0] if temporal else '',
            'temporal_floor_max': temporal[1] if temporal else '',
            'note': ''})
        # per-node detail at each mesh's primary comparison step
        if step == 80 or key == (32, 16):
            per_node[f'n{nx}x{ny}'] = (e, u64c, c['coords'], c['boundary'], step)
    # documented missing field: n32x16 at step 80 (corrupted archive)
    rows.append({'mesh': 'n32x16', 'step': 80, 't_compare': f'{80 * DT:g}',
                 'bulk_fe_rms': '', 'bulk_max': '', 'surf_fe_rms': '', 'surf_max': '',
                 'proxy_level_diff_fe': '', 'temporal_floor_fe': '',
                 'temporal_floor_max': '',
                 'note': 'NOT AVAILABLE: case.out.m0.80.vtu is a 0-byte file '
                         '(archive corrupted from step 40; left untouched)'})

    # orders over both intervals (each at its own comparison time)
    orders = []
    for metric in ('bulk_fe_rms', 'bulk_max', 'surf_fe_rms', 'surf_max'):
        for a, b, step in [((8, 4), (16, 8), 80), ((16, 8), (32, 16), 39)]:
            ra = next(r for r in rows if r['mesh'] == f'n{a[0]}x{a[1]}'
                      and r['step'] == step)
            rb = next(r for r in rows if r['mesh'] == f'n{b[0]}x{b[1]}'
                      and r['step'] == step)
            orders.append({'interval': f"{a[0]}x{a[1]}->{b[0]}x{b[1]}",
                           't_compare': f'{step * DT:g}', 'metric': metric,
                           'err_coarse': ra[metric], 'err_fine': rb[metric],
                           'order': math.log(ra[metric] / rb[metric], 2.0)
                           if rb[metric] > 0 else float('nan')})
    with (HERE / 'convergence_orders.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['interval', 't_compare', 'metric',
                                          'err_coarse', 'err_fine', 'order'])
        w.writeheader()
        w.writerows(orders)

    fields = ['mesh', 'step', 't_compare', 'bulk_fe_rms', 'bulk_max', 'surf_fe_rms',
              'surf_max', 'proxy_level_diff_fe', 'temporal_floor_fe',
              'temporal_floor_max', 'note']
    with (HERE / 'spatial_convergence_final.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})
    # per-node errors (reproducible detail)
    for tag, (e, up, coords, boundary, step) in per_node.items():
        bset = set(boundary.tolist())
        with (HERE / f'per_node_{tag}_step{step}_vs_proxy_n64x32.csv').open('w',
                                                                           newline='') as f:
            w = csv.writer(f)
            w.writerow(['node', 'x', 'y', f'u_oofem_t{step * DT:g}',
                        f'u_proxy_t{step * DT:g}', 'error', 'on_boundary'])
            for i in range(len(e)):
                w.writerow([i, coords[i, 0], coords[i, 1],
                            f'{up[i] + e[i]:.15g}', f'{up[i]:.15g}',
                            f'{e[i]:.15g}', i in bset])

    # plot: two panels, one per comparison time
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        styles = {'bulk_fe_rms': ('o-', 'tab:blue'), 'surf_fe_rms': ('s-', 'tab:red'),
                  'bulk_max': ('o--', 'tab:cyan'), 'surf_max': ('s--', 'tab:orange')}

        def panel(ax, pts, step, title):
            for key, (mk, col) in styles.items():
                hx = [LX / nx for nx, _ in pts]
                vals = [next(r for r in rows if r['mesh'] == f'n{nx}x{ny}'
                             and r['step'] == step)[key] for nx, ny in pts]
                ax.semilogy(hx, vals, mk, color=col, label=key)
            ax.set_xlabel('h_x = 2/nx')
            ax.set_title(title, fontsize=9)
            ax.grid(True, which='both', alpha=0.3)
            ax.legend(fontsize=7)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), dpi=150)
        panel(ax1, [(8, 4), (16, 8)], 80, 't = 0.02 (steps 1-80 intact)')
        panel(ax2, [(16, 8), (32, 16)], 39,
              't = 0.00975 (n32x16 archive corrupted from step 40)')
        ax1.set_ylabel('error at t_compare (vs DOP853 proxy n64x32)')
        fig.suptitle('Coupled hyperbolic Allen-Cahn: spatial convergence '
                     '(OOFEM archived, dt=2.5e-4)', fontsize=10)
        fig.tight_layout()
        fig.savefig(HERE / 'spatial_convergence.png')
        plt.close(fig)
    except Exception as ex:  # plot must never kill the stage
        print(f'[compare] plot failed: {ex}', flush=True)
    log('compare', time.perf_counter() - t0,
        'dual-time comparison (n32x16 step 80 missing): '
        'spatial_convergence_final.csv + convergence_orders.csv + per-node CSVs + PNG')
    return rows


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    t_all = time.perf_counter()
    if stage in ('check', 'all'):
        stage_check()
    if stage in ('tolcheck', 'all'):
        stage_tolcheck()
    if stage in ('control32', 'all'):
        stage_control32()
    if stage in ('proxy64', 'all'):
        stage_proxy64()
    if stage in ('tolcheck64', 'all'):
        stage_tolcheck64()
    if stage in ('compare', 'all'):
        stage_compare()
    write_log()
    print(f'total wall: {time.perf_counter() - t_all:.1f}s', flush=True)


if __name__ == '__main__':
    main()
