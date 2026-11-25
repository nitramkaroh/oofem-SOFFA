#!/usr/bin/env python3
"""Verification of the line-element Jacobian fix for the coupled hyperbolic
independent reference (see ../hyperbolic_coupled_discrepancy_audit.md).

Constraints honoured:
  * OOFEM is NOT re-run: the archived VTU histories of the three temporal runs
    (coupled_independent_time_n16x8_dt{0.002,0.001,0.0005}) are read as-is.
  * The original reference script run_coupled_independent.py is NOT modified;
    it is imported read-only and its archived trajectory
    (independent_reference.npz) is reused for the "before" numbers.
  * Only the DOP853 integration of the corrected reference is executed,
    sequentially, for the three archived time steps.
  * All new outputs are written under this directory only.

Correction (line elements only, xi in [-1,1], element length L):
    weight   w = (L/2) * GAUSS_W[a]        (was: L * GAUSS_W[a])
    dN/ds    [-.5,.5] * (2/L)              (was: [-.5,.5] / L)
No other model term, coefficient or initial condition is changed.
"""
from __future__ import annotations
import csv, math, sys, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

HERE = Path(__file__).resolve().parent
HYPERBOLIC = HERE.parent
sys.path.insert(0, str(HYPERBOLIC))
import run_coupled_independent as orig  # read-only; original reference

LX, LY, T = orig.LX, orig.LY, orig.T


# ---------------------------------------------------------------- assembly --
def assemble_fixed(nx, ny):
    """Same as orig.assemble, with the line-element isoparametric Jacobian."""
    coords, quads, edges, boundary = orig.mesh(nx, ny)
    nd = len(coords)
    Mb = lil_matrix((nd, nd)); Kb = lil_matrix((nd, nd))
    Ms = lil_matrix((nd, nd)); Ks = lil_matrix((nd, nd))
    for ids in quads:  # bulk: unchanged
        xy = coords[ids]
        for a, xi in enumerate(orig.GAUSS_X):
            for b, eta in enumerate(orig.GAUSS_X):
                N, dN = orig.qshape(xi, eta)
                J = xy.T @ dN; det = np.linalg.det(J)
                grad = dN @ np.linalg.inv(J)
                w = det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
                Mb[np.ix_(ids, ids)] += orig.RHO * np.outer(N, N) * w
                Kb[np.ix_(ids, ids)] += orig.KB * (grad @ grad.T) * w
    for ids in edges:  # FIXED: J = L/2, physical dN = dN/dxi * 2/L
        xy = coords[ids]
        length = np.linalg.norm(xy[1] - xy[0])
        for a, s in enumerate(orig.GAUSS_X):
            N = np.array([(1 - s) / 2, (1 + s) / 2])
            dN = np.array([-.5, .5]) * (2.0 / length)
            w = 0.5 * length * orig.GAUSS_W[a]
            Ms[np.ix_(ids, ids)] += orig.CS * np.outer(N, N) * w
            Ks[np.ix_(ids, ids)] += orig.KS * np.outer(dN, dN) * w
    return (coords, quads, edges, boundary,
            Mb.tocsc(), Kb.tocsc(), Ms.tocsc(), Ks.tocsc())


def residual_reaction_fixed(u, coords, quads, edges):
    """Same as orig.residual_reaction, with the line weight factor 1/2."""
    rb = np.zeros(len(u)); rs = np.zeros(len(u))
    for ids in quads:  # bulk: unchanged
        xy = coords[ids]; ul = u[ids]
        for a, xi in enumerate(orig.GAUSS_X):
            for b, eta in enumerate(orig.GAUSS_X):
                N, dN = orig.qshape(xi, eta)
                det = np.linalg.det(xy.T @ dN)
                uq = N @ ul
                rb[ids] += orig.MB * N * (uq - uq ** 3) * det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
    for ids in edges:  # FIXED
        xy = coords[ids]; length = np.linalg.norm(xy[1] - xy[0]); ul = u[ids]
        for a, s in enumerate(orig.GAUSS_X):
            N = np.array([(1 - s) / 2, (1 + s) / 2])
            uq = N @ ul
            rs[ids] += orig.MS * N * (uq - uq ** 3) * 0.5 * length * orig.GAUSS_W[a]
    return rb, rs


def energy_fixed(u, v, coords, quads, edges):
    """Same as orig.energy, with the line terms corrected for consistency."""
    ek = ep = es = esp = 0.
    for ids in quads:
        xy = coords[ids]; ul = u[ids]; vl = v[ids]
        for a, xi in enumerate(orig.GAUSS_X):
            for b, eta in enumerate(orig.GAUSS_X):
                N, dN = orig.qshape(xi, eta)
                J = xy.T @ dN; det = np.linalg.det(J)
                grad = dN @ np.linalg.inv(J)
                uq = N @ ul; vq = N @ vl
                ek += .5 * orig.RHO * vq * vq * det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
                ep += .5 * orig.KB * (grad.T @ ul) @ (grad.T @ ul) * det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
                es += .25 * orig.MB * (uq * uq - 1.) ** 2 * det * orig.GAUSS_W[a] * orig.GAUSS_W[b]
    for ids in edges:
        length = np.linalg.norm(coords[ids[1]] - coords[ids[0]])
        ul = u[ids]
        for a, s in enumerate(orig.GAUSS_X):
            N = np.array([(1 - s) / 2, (1 + s) / 2])
            dN = np.array([-.5, .5]) * (2.0 / length)
            w = 0.5 * length * orig.GAUSS_W[a]
            uq = N @ ul
            esp += .5 * orig.KS * (dN @ ul) ** 2 * w + .25 * orig.MS * (uq * uq - 1.) ** 2 * w
    return ek + ep + es + esp, ek, ep + es, esp


# ------------------------------------------------------------- references --
def reference_fixed(nx, ny, times):
    """DOP853 of the corrected semidiscrete system; same structure/tolerances
    as orig.reference (C = ETA*Mb + Ms, K = Kb + Ks, Mb only in the inverse).
    """
    ctx = assemble_fixed(nx, ny)
    coords, quads, edges, boundary, Mb, Kb, Ms, Ks = ctx
    mb_lu = orig.splu(Mb)
    C = orig.ETA * Mb + Ms
    K = Kb + Ks

    def rhs(t, z):
        u = z[:len(coords)]; v = z[len(coords):]
        rb, rs = residual_reaction_fixed(u, coords, quads, edges)
        a = mb_lu.solve(-(C @ v) - (K @ u) + rb + rs)
        return np.r_[v, a]

    y0 = np.r_[orig.initial_u(coords), np.zeros(len(coords))]
    sol = solve_ivp(rhs, (0, float(times[-1])), y0, t_eval=times,
                    method='DOP853', rtol=2e-12, atol=2e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    return ctx, sol.y[:len(coords), :].T, sol.y[len(coords):, :].T


# ----------------------------------------------------------------- reading --
def read_oofem_steps(folder):
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    ts = []; us = []; pts_first = None
    for ds in pvd.findall('.//DataSet'):
        vtu = ET.parse(folder / ds.attrib['file']).getroot()
        pts = np.fromstring(vtu.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)[:, :2]
        arr = next(a for a in vtu.findall('.//PointData/DataArray')
                   if a.attrib.get('Name') == 'Temperature')
        if pts_first is None:
            pts_first = pts
        ts.append(float(ds.attrib['timestep']))
        us.append(np.fromstring(arr.text or '', sep=' '))
    return np.array(ts), np.array(us), pts_first


def error_row(u_oofem, u_ref, Mb, Ms, boundary):
    e = u_oofem - u_ref
    mask = np.zeros(len(u_oofem), dtype=bool); mask[boundary] = True
    return {
        'bulk_rms': orig.rms_fe(Mb, e, LX * LY),
        'surf_rms': orig.rms_fe(Ms, e, 2 * (LX + LY)),
        'bulk_max': float(np.max(np.abs(e[~mask]))),
        'surf_max': float(np.max(np.abs(e[mask]))),
        'all_max': float(np.max(np.abs(e))),
    }


def main():
    summary = []
    for dt in (0.002, 0.001, 0.0005):
        tag = f'time_n16x8_dt{dt:g}'
        folder = HYPERBOLIC / f'coupled_independent_{tag}'
        mesh = np.load(folder / 'mesh.npz')
        coords, quads, edges, boundary = mesh['coords'], mesh['quads'], mesh['edges'], mesh['boundary']
        nsteps = round(T / dt)
        times = np.linspace(0, T, nsteps + 1)

        # OOFEM archived history (untouched); vtkxml starts at t=dt (no t=0 file)
        ts, Uoofem, pts = read_oofem_steps(folder)
        assert len(ts) == nsteps and np.allclose(ts, times[1:], atol=1e-12), 'vtu time mismatch'
        assert np.allclose(ts[-1], T, atol=1e-10)
        assert np.allclose(pts, coords, rtol=0, atol=1e-13), 'vtu coordinates mismatch'

        # "before": original reference trajectory as archived
        zref = np.load(folder / 'independent_reference.npz')
        Uo = zref['U']
        assert np.allclose(zref['times'], times, atol=1e-12)
        assert np.allclose(zref['coords'], coords, rtol=0, atol=1e-13)

        ctx_o = orig.assemble(16, 8)
        Mb_o, Kb_o, Ms_o, Ks_o = ctx_o[4:]  # assemble returns (coords,quads,edges,boundary,Mb,Kb,Ms,Ks)

        # "after": corrected reference, DOP853 only
        print(f'[{tag}] running corrected DOP853 reference ...', flush=True)
        ctx_f, Uf, Vf = reference_fixed(16, 8, times)
        _, _, Ms_f, Ks_f = ctx_f[4:]

        # per-step comparison (OOFEM file k corresponds to reference index k+1)
        rows = []
        for i in range(1, nsteps + 1):
            u = Uoofem[i - 1]
            b = error_row(u, Uo[i], Mb_o, Ms_o, boundary)
            a = error_row(u, Uf[i], Mb_o, Ms_f, boundary)
            before_surf_correctnorm = orig.rms_fe(Ms_f, u - Uo[i], 2 * (LX + LY))
            after_surf_on_orig_norm = orig.rms_fe(Ms_o, u - Uf[i], 2 * (LX + LY))
            rows.append({'step': i, 't': ts[i - 1],
                         'bulk_rms_before': b['bulk_rms'],
                         'surf_rms_before': b['surf_rms'],          # original norm (2x Ms)
                         'surf_rms_before_correctnorm': before_surf_correctnorm,
                         'bulk_rms_after': a['bulk_rms'],
                         'surf_rms_after': a['surf_rms'],           # correct norm (1x Ms)
                         'surf_rms_after_on_orig_norm': after_surf_on_orig_norm,
                         'bulk_max_before': b['bulk_max'], 'surf_max_before': b['surf_max'],
                         'bulk_max_after': a['bulk_max'], 'surf_max_after': a['surf_max'],
                         'all_max_before': b['all_max'], 'all_max_after': a['all_max']})

        out_csv = HERE / f'per_step_{tag}.csv'
        with out_csv.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

        # energy of the corrected reference (internally consistent)
        quads_l = orig.mesh(16, 8)[1]
        ens = np.array([energy_fixed(Uf[i], Vf[i], coords, quads_l, edges) for i in range(nsteps + 1)])
        with (HERE / f'independent_energy_fixed_{tag}.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['time', 'total', 'bulk_kinetic', 'bulk_potential_gradient', 'surface'])
            for t, e in zip(times, ens):
                w.writerow([t, *e])
        np.savez(HERE / f'independent_reference_fixed_{tag}.npz',
                 times=times, U=Uf, V=Vf, boundary=boundary,
                 Ms=Ms_f.tocoo(), Ks=Ks_f.tocoo())

        fin = rows[-1]
        summary.append({'dt': dt, **{k: fin[k] for k in
            ('bulk_rms_before', 'surf_rms_before', 'surf_rms_before_correctnorm',
             'bulk_rms_after', 'surf_rms_after', 'bulk_max_before', 'surf_max_before',
             'bulk_max_after', 'surf_max_after')}})
        print(f'[{tag}] done: after bulk={fin["bulk_rms_after"]:.6e} '
              f'surf={fin["surf_rms_after"]:.6e} '
              f'(before bulk={fin["bulk_rms_before"]:.6e} surf={fin["surf_rms_before"]:.6e})', flush=True)

    # orders of convergence (final time, corrected reference)
    srt = sorted(summary, key=lambda r: -r['dt'])
    for key in ('bulk_rms_after', 'surf_rms_after'):
        for i in range(len(srt) - 1):
            srt[i][key + '_order_to_next'] = math.log(srt[i][key] / srt[i + 1][key], 2.0)
    srt[0]['bulk_rms_before_order_to_next'] = math.log(
        srt[0]['bulk_rms_before'] / srt[1]['bulk_rms_before'], 2.0)
    srt[0]['surf_rms_before_order_to_next'] = math.log(
        srt[0]['surf_rms_before'] / srt[1]['surf_rms_before'], 2.0)

    with (HERE / 'summary_final.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in summary for k in r)))
        w.writeheader(); w.writerows(summary)
    print('\nfinal summary:')
    for r in srt:
        print(' ', r, flush=True)


if __name__ == '__main__':
    main()
