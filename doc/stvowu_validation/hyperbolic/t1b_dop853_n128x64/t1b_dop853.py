#!/usr/bin/env python3
"""T1-B: independent DOP853 reference on n128x64 (proxy refinement) — driver.

Constraint set (user-authorized 2026-10-27):
  * FASE 0 static audit -> see t1b_fase0_cost.py + fase0_time_report.txt.
  * FASE 1 `pilot`:  single process/core, hard 60 s cap; measure wall, RSS,
    nfev, integrator behaviour. NOT a proof the full run fits the budget.
  * FASE 2 `full`:   single full integration [0,0.02] (81 pts). In-process
    guards: NaN/Inf, |u|<5 (anomalous growth), wall<40 min. Outer `timeout
    2700` + `ulimit -v` enforced by the shell wrapper. No auto-retry.
  * FASE 3 `validate`: finiteness/completeness/final time + comparisons vs
    available levels (norms/conventions identical to T1).

Assembly = check_jacobian_fix.assemble_fixed / residual_reaction_fixed
(corrected line Jacobian w=(L/2)*Gauss, dN=[-.5,.5]*(2/L)); independent of
OOFEM (pure-Python, no OOFEM code/binary). Tolerances rtol=2e-12 atol=2e-14
(identical to the accepted n64x32 proxy). Only t_eval (81 pts) is stored.

Assembly cache (why): the pure-Python assembly of n128x64 takes ~28 s. The
bounded 60 s pilot cannot afford assembly + the first DOP853 step (~30 s,
initial step-size exploration, ~50 RHS) in one process. The `cache` stage
(assembled fresh once) serializes the four sparse matrices + mesh; the pilot
loads it (~2 s) and measures a real integration step; the `full` stage
RE-ASSEMBLES FRESH and asserts EXACT equality with the cache before
integrating, so no unverified shortcut enters the reference.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, resource, sys, time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

HERE = Path(__file__).resolve().parent
HYPERBOLIC = HERE.parent
sys.path.insert(0, str(HYPERBOLIC))
sys.path.insert(0, str(HYPERBOLIC / 'coupled_independent_jacobian_fix'))
import run_coupled_independent as orig   # read-only (mesh, constants, splu)
import check_jacobian_fix as fj          # read-only (corrected assembly)

LX, LY, T = orig.LX, orig.LY, orig.T          # 2.0, 1.0, 0.02
DT = 2.5e-4
NSTEP = round(T / DT)                          # 80
TIMES = np.linspace(0, T, NSTEP + 1)
AREA = LX * LY
PERIM = 2 * (LX + LY)
RTOL, ATOL = 2e-12, 2e-14
U_BOUND = 5.0                                  # anomalous-growth tripwire (phys |u|~<=1.5)

NX, NY = 128, 64
NPZ = HERE / 'reference_n128x64_dt0p00025.npz'
CACHE_NPZ = HERE / 'matrices_n128x64.npz'
CACHE_META = HERE / 'matrices_n128x64.json'
COST_CSV = HERE / 't1b_cost_log.csv'


def rss_kb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


class Guard(RuntimeError):
    pass


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def save_cache(ctx):
    coords, quads, edges, boundary, Mb, Kb, Ms, Ks = ctx
    arrays = {}
    for pfx, M in (('Mb', Mb), ('Kb', Kb), ('Ms', Ms), ('Ks', Ks)):
        Mc = M.tocoo()
        arrays[pfx + '_data'] = Mc.data
        arrays[pfx + '_row'] = Mc.row
        arrays[pfx + '_col'] = Mc.col
    np.savez(CACHE_NPZ, **arrays, coords=coords, quads=quads, edges=edges,
             boundary=boundary)
    meta = {'sha256': sha256_of(CACHE_NPZ), 'nx': NX, 'ny': NY, 'n': len(coords),
            'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
            'source': 'fj.assemble_fixed(128, 64) fresh (corrected line Jacobian)'}
    CACHE_META.write_text(json.dumps(meta, indent=2))
    return meta


def load_cache_ctx():
    if not (CACHE_NPZ.exists() and CACHE_META.exists()):
        raise Guard('matrix cache missing: run the `cache` stage first')
    meta = json.loads(CACHE_META.read_text())
    if meta['sha256'] != sha256_of(CACHE_NPZ):
        raise Guard('matrix cache sha256 mismatch: refusing to use it (no retry)')
    d = np.load(CACHE_NPZ)
    N = meta['n']
    M = lambda pfx: coo_matrix((d[pfx + '_data'], (d[pfx + '_row'], d[pfx + '_col'])),
                               shape=(N, N)).tocsc()
    coords, quads, edges = d['coords'], d['quads'], d['edges']
    boundary = d['boundary']
    assert len(coords) == N
    return (coords, quads, edges, boundary, M('Mb'), M('Kb'), M('Ms'), M('Ks')), meta


def build_rhs(ctx):
    coords, quads, edges, boundary, Mb, Kb, Ms, Ks = ctx
    mb_lu = orig.splu(Mb)
    C = orig.ETA * Mb + Ms
    K = Kb + Ks
    N = len(coords)
    state = {'nfev': 0, 't0': time.perf_counter(), 'rss_peak': rss_kb()}

    def rhs(t, z):
        state['nfev'] += 1
        # in-process safety guards (cheap; every RHS call)
        if not np.isfinite(z).all():
            raise Guard(f'NaN/Inf in state at nfev={state["nfev"]}: aborting (no retry)')
        if np.max(np.abs(z[:N])) > U_BOUND:
            raise Guard(f'|u|_max={np.max(np.abs(z[:N])):.3e} > {U_BOUND}: anomalous growth, '
                        f'aborting (no retry)')
        state['rss_peak'] = max(state['rss_peak'], rss_kb())
        u = z[:N]
        v = z[N:]
        rb, rs = fj.residual_reaction_fixed(u, coords, quads, edges)
        a = mb_lu.solve(-(C @ v) - (K @ u) + rb + rs)
        return np.r_[v, a]

    return rhs, state, N


def append_cost(row):
    new = not COST_CSV.exists()
    with COST_CSV.open('a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if new:
            w.writeheader()
        w.writerow(row)
    print(f'[cost] {row}', flush=True)


# -------------------------------------------------------------------- cache --
def cmd_cache(_):
    t0 = time.perf_counter()
    print(f'[cache] fresh assembly n{NX}x{NY} (corrected line Jacobian) ...', flush=True)
    ctx = fj.assemble_fixed(NX, NY)
    meta = save_cache(ctx)
    print(f'[cache] done in {time.perf_counter()-t0:.1f}s; sha256={meta["sha256"][:16]}... '
          f'N={meta["n"]}', flush=True)


# -------------------------------------------------------------------- pilot --
def cmd_pilot(_):
    t0 = time.perf_counter()
    ctx, meta = load_cache_ctx()
    t_load = time.perf_counter() - t0
    N = len(ctx[0])
    print(f'[pilot] cache loaded in {t_load:.2f}s (sha256 {meta["sha256"][:16]}...), '
          f'N={N}', flush=True)
    if t_load > 10:
        raise Guard(f'cache load took {t_load:.1f}s (> 10 s): degraded environment, '
                    f'stopping (no retry)')
    rhs, state, _ = build_rhs(ctx)
    y0 = np.r_[orig.initial_u(ctx[0]), np.zeros(N)]
    # Bounded pilot: 1 output step (t <= 2.5e-4), including DOP853 initial
    # step-size exploration. Hard 60 s cap enforced by the shell wrapper
    # (timeout 60); soft in-process cap 55 s.
    ptimes = TIMES[:2]
    if time.perf_counter() - t0 > 55:
        raise Guard('soft 55 s pilot cap hit before integration: stopping (no retry)')
    t1 = time.perf_counter()
    sol = solve_ivp(rhs, (0, float(ptimes[-1])), y0, t_eval=ptimes,
                    method='DOP853', rtol=RTOL, atol=ATOL)
    t_int = time.perf_counter() - t1
    if not sol.success:
        raise RuntimeError(f'pilot integration failed: {sol.message}')
    assert np.isfinite(sol.y).all(), 'pilot produced non-finite state'
    nfev = int(sol.nfev)
    ms_per_rhs = t_int / max(nfev, 1) * 1000.0
    est_full_nfev = nfev + (NSTEP - 1) * 13.6   # first step + n64 steady-state rate
    est_full_s = 28.0 + 0.02 + est_full_nfev * t_int / max(nfev, 1)
    row = {'stage': 'pilot', 'wall_s': f'{time.perf_counter()-t0:.2f}',
           'load_s': f'{t_load:.2f}', 'int_s': f'{t_int:.2f}', 'n': N, 'nfev': nfev,
           'ms_per_rhs': f'{ms_per_rhs:.1f}', 'rss_peak_mb': f'{state["rss_peak"]/1024:.1f}',
           'n_accepted': len(sol.t), 'first_step_h': f'{float(sol.t[1]) if len(sol.t) > 1 else float(sol.t[-1]):.6e}',
           'note': (f'est full: nfev~{est_full_nfev:.0f} -> ~{est_full_s/60:.1f} min '
                    f'(cost model, NOT a budget proof)')}
    append_cost(row)
    print(f'[pilot] int={t_int:.1f}s: nfev={nfev}, {ms_per_rhs:.0f} ms/RHS, '
          f'n_accepted_steps={len(sol.t)}, rss_peak={state["rss_peak"]/1024:.1f} MB; '
          f'est full ~{est_full_s/60:.1f} min', flush=True)
    print('[pilot] NOTE: pilot success does NOT prove the full run fits the budget.',
          flush=True)


# --------------------------------------------------------------------- full --
def cmd_full(_):
    t_start = time.perf_counter()
    SOFT_WALL = 40 * 60   # in-process cap under the 45-min outer timeout
    # 1) cache must exist and be intact
    cache_ctx, meta = load_cache_ctx()
    # 2) FRESH assembly + exact-equality verification against the cache
    print(f'[full] fresh assembly n{NX}x{NY} (cache verification) ...', flush=True)
    t0 = time.perf_counter()
    fresh = fj.assemble_fixed(NX, NY)
    t_asm = time.perf_counter() - t0
    for i, name in enumerate(('Mb', 'Kb', 'Ms', 'Ks')):
        diff = fresh[4 + i] - cache_ctx[4 + i]
        d = float(abs(diff).max()) if diff.nnz else 0.0
        print(f'[full] cache check {name}: max|fresh - cache| = {d:.3e}', flush=True)
        if d > 1e-18:
            raise Guard(f'{name} fresh assembly differs from cache (max {d:.3e}): '
                        f'non-deterministic assembly suspected, stopping (no retry)')
    print(f'[full] fresh assembly verified exact in {t_asm:.1f}s; integrating ...',
          flush=True)
    ctx = fresh
    N = len(ctx[0])
    rhs, state, _ = build_rhs(ctx)
    y0 = np.r_[orig.initial_u(ctx[0]), np.zeros(N)]

    def guarded(t, z):
        if time.perf_counter() - t_start > SOFT_WALL:
            raise Guard(f'wall {SOFT_WALL}s exceeded: aborting (no retry)')
        return rhs(t, z)

    t0 = time.perf_counter()
    sol = solve_ivp(guarded, (0, float(TIMES[-1])), y0, t_eval=TIMES,
                    method='DOP853', rtol=RTOL, atol=ATOL)
    wall = time.perf_counter() - t0
    if not sol.success:
        raise RuntimeError(f'full integration failed: {sol.message}')
    U = sol.y[:N, :].T
    V = sol.y[N:, :].T
    assert np.isfinite(U).all() and np.isfinite(V).all(), 'non-finite trajectory: stopping'
    nfev = int(sol.nfev)
    np.savez(NPZ, U=U, V=V, times=sol.t, t_eval=TIMES, nx=NX, ny=NY,
             rtol=RTOL, atol=ATOL, nfev=nfev, n=N,
             wall_s=wall, rss_peak_kb=state['rss_peak'],
             final_t=float(sol.t[-1]), cache_sha256=meta['sha256'])
    row = {'stage': 'full', 'asm_s': f'{t_asm:.2f}', 'wall_s': f'{wall:.2f}', 'n': N,
           'nfev': nfev, 'nfev_per_step': f'{nfev/NSTEP:.2f}',
           'rss_peak_mb': f'{state["rss_peak"]/1024:.1f}',
           'final_t': f'{float(sol.t[-1]):.17g}',
           'note': f'saved {NPZ.name} ({U.nbytes/1e6:.1f} MB U + {V.nbytes/1e6:.1f} MB V)'}
    append_cost(row)
    print(f'[full] done in {wall:.1f}s ({wall/60:.1f} min): nfev={nfev}, '
          f'rss_peak={state["rss_peak"]/1024:.1f} MB, final_t={float(sol.t[-1]):.17g}',
          flush=True)


# ----------------------------------------------------------------- validate --
def cmd_validate(_):
    d = np.load(NPZ)
    U = d['U']
    print(f'[validate] {NPZ.name}: U={U.shape}, finite={np.isfinite(U).all()}, '
          f'final_t={float(d["final_t"]):.17g} (target {T}), nfev={int(d["nfev"])}')
    assert U.shape == (NSTEP + 1, int(d['n'])), 'trajectory shape mismatch'
    assert np.isfinite(U).all(), 'non-finite trajectory'
    assert abs(float(d['final_t']) - T) < 1e-12, 'did not reach final time'
    assert np.array_equal(d['t_eval'], TIMES), 't_eval mismatch vs canonical grid'
    print('[validate] finiteness / completeness / final-time / grid: OK')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['cache', 'pilot', 'full', 'validate'])
    args = ap.parse_args()
    {'cache': cmd_cache, 'pilot': cmd_pilot, 'full': cmd_full,
     'validate': cmd_validate}[args.stage](args)


if __name__ == '__main__':
    main()
