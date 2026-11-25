#!/usr/bin/env python3
"""T1-B FASE 0: static cost measurement for the DOP853 n128x64 reference.

NO integration is performed. Measures: assembly time/nnz, splu factor time and
size, per-RHS-evaluation wall time (30 frozen-state calls), baseline/peak
memory (RSS + VmSize), IC sanity, times sanity. Read-only on all modules.
"""
import resource
import sys
import time
from pathlib import Path

import numpy as np

HYPER = Path('/home/marco/Documents/Work/oofem-SOFFA/doc/stvowu_validation/hyperbolic')
sys.path.insert(0, str(HYPER))
sys.path.insert(0, str(HYPER / 'coupled_independent_jacobian_fix'))
sys.path.insert(0, str(HYPER / 'coupled_independent_spatial_jacobian_fix'))
import run_coupled_independent as orig  # read-only
import check_jacobian_fix as fj         # read-only (corrected assembly)

DT = 2.5e-4
NSTEP = 80
TIMES = np.linspace(0, orig.T, NSTEP + 1)


def vm_kb():
    with open('/proc/self/status') as f:
        for line in f:
            if line.startswith('VmSize:'):
                return int(line.split()[1])
    return -1


def rss_kb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


NX, NY = 128, 64
print(f'=== FASE 0 static cost: n{NX}x{NY} ===')
print(f'baseline after imports: VmSize={vm_kb()/1024:.1f} MB  RSS_peak={rss_kb()/1024:.1f} MB')

t0 = time.perf_counter()
coords, quads, edges, boundary, Mb, Kb, Ms, Ks = fj.assemble_fixed(NX, NY)
t_asm = time.perf_counter() - t0
N = len(coords)
print(f'assembly: {t_asm:.2f}s  N={N}  quads={len(quads)}  edges={len(edges)}  '
      f'boundary_nodes={len(boundary)}')
for name, M in (('Mb', Mb), ('Kb', Kb), ('Ms', Ms), ('Ks', Ks)):
    assert M.format == 'csc'
    assert M.nnz < N * N, 'matrix denser than expected: STOP'
    print(f'  {name}: format={M.format} shape={M.shape} nnz={M.nnz} '
          f'(dense would be {N*N:,} entries -> ratio {M.nnz/(N*N)*100:.4f}%)')

# symmetry sanity (all four must be symmetric: verified assembly)
for name, M in (('Mb', Mb), ('Kb', Kb), ('Ms', Ms), ('Ks', Ks)):
    d = abs(M - M.T).max()
    print(f'  {name} max|A - A^T| = {d:.3e}')

t0 = time.perf_counter()
mb_lu = orig.splu(Mb)
t_splu = time.perf_counter() - t0
L, Uu = mb_lu.L, mb_lu.U
print(f'splu(Mb): {t_splu:.2f}s  L.nnz={L.nnz:,}  U.nnz={Uu.nnz:,}  '
      f'factor mem ~ {8*(L.nnz+Uu.nnz)/1e6:.1f} MB')

C = orig.ETA * Mb + Ms
K = Kb + Ks
print(f'C=ETA*Mb+Ms nnz={C.nnz:,}  K=Kb+Ks nnz={K.nnz:,} (both csc={C.format},{K.format})')


def rhs(t, z):
    u = z[:N]
    v = z[N:]
    rb, rs = fj.residual_reaction_fixed(u, coords, quads, edges)
    a = mb_lu.solve(-(C @ v) - (K @ u) + rb + rs)
    return np.r_[v, a]


y0 = np.r_[orig.initial_u(coords), np.zeros(N)]
print(f'IC: u0 in [{y0[:N].min():.6f}, {y0[:N].max():.6f}]  finite={np.isfinite(y0).all()}  '
      f'v0 all zero={not y0[N:].any()}')
print(f'times: {len(TIMES)} pts, t[1]={TIMES[1]:.17g} (dt={DT:g}), '
      f't[-1]={TIMES[-1]:.17g} (T={orig.T}), max|t_k - k*dt|={max(abs(t - k * DT) for k, t in enumerate(TIMES)):.2e}')

# time 30 frozen-state RHS evaluations (first one includes lazy setup; report both)
z = y0.copy()
t0 = time.perf_counter()
for _ in range(30):
    out = rhs(0.0, z)
t30 = time.perf_counter() - t0
assert np.isfinite(out).all()
ms = t30 / 30 * 1000
print(f'RHS: 30 evals in {t30:.3f}s -> {ms:.2f} ms/eval  '
      f'(state finite: {np.isfinite(out).all()})')
# one isolated eval for variance
t1 = time.perf_counter()
rhs(0.0, z)
t1 = time.perf_counter() - t1
print(f'RHS single repeat: {t1*1000:.2f} ms')

for nfev in (1088, 1300, 1600):
    tot = t_asm + t_splu + nfev * t30 / 30
    print(f'  est total if nfev={nfev}: {tot:.0f}s = {tot/60:.1f} min')
print(f'memory now: VmSize={vm_kb()/1024:.1f} MB  RSS_peak={rss_kb()/1024:.1f} MB')
print('trajectory storage needed: U,V at 81 t_eval pts = '
      f'{2*81*N*8/1e6:.1f} MB (no dense_output, no internal steps)')
print('DONE (no integration performed)')
