#!/usr/bin/env python3
"""Postprocess the Nizovtseva traveling-wave OOFEM runs.

Reads the raw VTU histories in ../results/nizovtseva_n{80,160,320}/ and
produces (in the parent directory):

  nizovtseva_metrics.csv        per-mesh error vs the exact formula at T,
                                observed spatial orders, interface data
  nizovtseva_profiles.csv       u(x,t) at t = 0, .05, .10, .15, .20 for each mesh
  nizovtseva_interface_n320.csv x_0(t) of the u = 1/2 level (n=320) vs exact c t
  nizovtseva_error_n320_at_T.csv nodal error u_OOFEM - u_exact at T (n=320)

Metrics: RMS nodal and max nodal error of the collapsed (y-averaged) field
against the exact formula at the unique uniform x-grid; an "inner" metric
restricted to |x| <= 6 away from the truncated Neumann ends. The interface
is the linearly interpolated u = 1/2 crossing, defined identically for the
exact formula and the FE field.
"""
from __future__ import annotations
import csv, math
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RESULTS = BASE / 'results'

TAU = 0.2
T = 0.2
LAMBDA = math.sqrt(2.0 + 9.0 * TAU) / 2.0
SPEED = 1.5 / LAMBDA
SNAPSHOTS = (0.0, 0.05, 0.10, 0.15, 0.20)


def exact(x, t):
    z = np.clip(LAMBDA * (np.asarray(x) - SPEED * t), -700., 700.)
    return 1.0 / (1.0 + np.exp(z))


def read_vtu(path):
    root = ET.parse(path).getroot()
    points = np.fromstring(root.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)
    pd = root.find('.//PointData')
    arr = next(a for a in pd.findall('DataArray') if a.attrib.get('Name') == 'Temperature')
    u = np.fromstring(arr.text or '', sep=' ')
    return points, u


def load_history(folder: Path):
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    frames = []
    for ds in pvd.findall('.//DataSet'):
        pts, u = read_vtu(folder / ds.attrib['file'])
        frames.append((float(ds.attrib['timestep']), pts, u))
    return frames


def collapse_y(pts, u):
    x = pts[:, 0]
    xu = np.unique(x)
    uu = np.array([np.mean(u[np.isclose(x, z)]) for z in xu])
    return xu, uu


def zero_crossing(x, u):
    for i in range(len(x) - 1):
        if (u[i] - .5) * (u[i + 1] - .5) <= 0 and u[i] != u[i + 1]:
            return float(x[i] + (.5 - u[i]) * (x[i + 1] - x[i]) / (u[i + 1] - u[i]))
    return float('nan')


def main() -> None:
    rows, profiles, prof_headers = [], [], ['x']
    for nx in (80, 160, 320):
        folder = RESULTS / f'nizovtseva_n{nx}'
        frames = load_history(folder)
        assert abs(frames[-1][0] - T) < 1e-10
        # final field and error vs exact
        xf, uf = collapse_y(frames[-1][1], frames[-1][2])
        ue = exact(xf, T)
        err = uf - ue
        inner = np.abs(xf) <= 6.0 + 1e-12
        row = {
            'n': nx, 'h': 24.0 / nx, 'dt': 0.001, 'T': T,
            'rms_full': float(np.sqrt(np.mean(err ** 2))),
            'max_full': float(np.max(np.abs(err))),
            'rms_inner_absx_le6': float(np.sqrt(np.mean(err[inner] ** 2))),
            'max_inner_absx_le6': float(np.max(np.abs(err[inner]))),
            'u_min': float(np.min(uf)), 'u_max': float(np.max(uf)),
            'x0_exact': SPEED * T,
            'x0_numerical': zero_crossing(xf, uf),
            'speed_numerical': zero_crossing(xf, uf) / T,
        }
        rows.append(row)
        # error vs exact at T, for the finest mesh (figure panel D)
        if nx == 320:
            with (BASE / 'nizovtseva_error_n320_at_T.csv').open('w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['x', 'u_OOFEM', 'u_exact', 'error_OOFEM_minus_exact'])
                for a, b, c in zip(xf, uf, ue):
                    w.writerow([a, b, c, b - c])
            # interface tracking at every saved frame
            recs = []
            for t, pts, u in frames:
                xi, ui = collapse_y(pts, u)
                recs.append((t, zero_crossing(xi, ui)))
            with (BASE / 'nizovtseva_interface_n320.csv').open('w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['time', 'x0_OOFEM_u1_2', 'x0_exact_c_t'])
                for t, x0 in recs:
                    w.writerow([t, x0, SPEED * t])
        # profiles at selected times (t = 0 from the exact initial condition,
        # because the PVD starts at t = dt)
        snap = {}
        for t, pts, u in frames:
            for s in SNAPSHOTS:
                if abs(t - s) < 1e-10:
                    snap[s] = collapse_y(pts, u)
        if 0.0 not in snap:
            x0g = np.linspace(-12, 12, nx + 1)
            snap[0.0] = (x0g, exact(x0g, 0.0))   # exact initial condition
        for s in SNAPSHOTS:
            xs, us = snap[s]
            profiles.append((s, nx, xs, us))
        print(row, flush=True)
    # observed spatial orders (RMS full and inner)
    for a, b in zip(rows, rows[1:]):
        b['order_rms_full'] = math.log(a['rms_full'] / b['rms_full'], 2)
        b['order_rms_inner'] = math.log(a['rms_inner_absx_le6'] / b['rms_inner_absx_le6'], 2)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (BASE / 'nizovtseva_metrics.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with (BASE / 'nizovtseva_profiles.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'n', 'x', 'u'])
        for (t, nx, xs, us) in profiles:
            for xx, uu in zip(xs, us):
                w.writerow([f'{t:g}', nx, xx, uu])
    print('metrics written to', BASE / 'nizovtseva_metrics.csv')


if __name__ == '__main__':
    main()
