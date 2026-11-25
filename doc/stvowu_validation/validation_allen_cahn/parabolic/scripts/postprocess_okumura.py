#!/usr/bin/env python3
"""Postprocess the Okumura Ch.4 example 2 OOFEM run.

Computes, from the raw VTU frames in ../results/okumura_ch4_case2/:

* total energy  E(t) = int_Omega [1/2 |grad u|^2 + F(u)] dA + sum_ends F(psi),
  F(s) = (gamma/4)(s^2-1)^2, gamma = 100; the bulk integral uses a 3x3 Gauss
  rule, which is exact for the Q1 quartic potential and for the quadratic
  gradient term; each endpoint edge contributes F at its two nodes
  (unit-length 2-point Gauss);
* spatial profiles at selected times;
* 1D-consistency checks (variation across the extruded depth);
* energy monotonicity;
* comparison against the digitized literature curve
  (okumura_ch4_example2_paper_fig44_digitized.csv, retained copy of the
  digitization of Fig. 4.4 of the dissertation) over 0.01 <= t <= 0.4.

Outputs (written next to this script's parent, i.e. ../):
  okumura_example2_energy.csv
  okumura_example2_profiles.csv
  okumura_example2_literature_comparison.csv
"""
from __future__ import annotations
import csv, re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT = BASE / 'results' / 'okumura_ch4_case2'
DIGITIZED = BASE / 'literature_okumura_fig44_digitized.csv'

NX, DT, T = 100, 0.001, 0.7
GAMMA = 100.0
SNAPSHOTS = (0.001, 0.1, 0.2, 0.3, 0.4, 0.7)

# 3x3 Gauss on [-1,1], exact for polynomials of degree <= 5 per direction.
GP = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])
GW = np.array([5 / 9, 8 / 9, 5 / 9])
# 2-point Gauss for the unit endpoint edges.
G2 = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])
W2 = np.array([1.0, 1.0])


def read_frame(path: Path):
    piece = ET.parse(path).getroot().find('.//Piece')
    xyz = np.fromstring(piece.find('./Points/DataArray').text, sep=' ').reshape(-1, 3)
    vals = np.fromstring(piece.find('./PointData/DataArray').text, sep=' ')
    u = np.empty((2, NX + 1))
    for (x, y), v in zip(xyz[:, :2], vals):
        u[round(y), round(x * NX)] = v
    return u


def bulk_energy(u: np.ndarray) -> float:
    e = 0.0
    for i in range(NX):
        q = np.array([u[0, i], u[0, i + 1], u[1, i + 1], u[1, i]])
        # Cell-center derivatives: du/dx = dN/dxi * (2/NX) * (1/2)(dU_b + dU_t)
        du_dx = 0.25 * (u[0, i + 1] - u[0, i] + u[1, i + 1] - u[1, i]) * 2 * NX
        du_dy = 0.25 * (u[1, i + 1] - u[0, i + 1] + u[1, i] - u[0, i]) * 2
        for ixi, a in enumerate(GP):
            for jeta, b in enumerate(GP):
                z = (0.25 * (1 - a) * (1 - b) * q[0] + 0.25 * (1 + a) * (1 - b) * q[1]
                     + 0.25 * (1 + a) * (1 + b) * q[2] + 0.25 * (1 - a) * (1 + b) * q[3])
                e += 0.25 / NX * (0.5 * (du_dx ** 2 + du_dy ** 2)
                                  + GAMMA / 4 * (z * z - 1) ** 2) * GW[ixi] * GW[jeta]
    return e


def end_energy(u: np.ndarray) -> float:
    e = 0.0
    for end in (u[:, 0], u[:, -1]):
        for b, w in zip(G2, W2):
            z = 0.5 * (1 - b) * end[0] + 0.5 * (1 + b) * end[1]
            e += 0.5 * (GAMMA / 4 * (z * z - 1) ** 2) * w
    return e


def main() -> None:
    files = sorted(OUT.glob('okumura_ch4_case2.out.m0.*.vtu'),
                   key=lambda p: int(re.search(r'm0\.(\d+)', p.name).group(1)))
    assert len(files) == round(T / DT), f'expected {T / DT} frames, found {len(files)}'
    times, energies, fields, y_spread = [], [], [], []
    for path in files:
        idx = int(re.search(r'm0\.(\d+)', path.name).group(1))
        u = read_frame(path)
        y_spread.append(float(np.max(np.abs(u[0] - u[1]))))
        times.append(idx * DT)
        energies.append(bulk_energy(u) + end_energy(u))
        fields.append(u[0].copy())
    times = np.array(times)
    energies = np.array(energies)

    with (BASE / 'okumura_example2_energy.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'OOFEM_Q1_total_energy_3x3gauss'])
        w.writerows(zip(times, energies))

    ids = [int(np.argmin(np.abs(times - t))) for t in SNAPSHOTS]
    snap_t = times[ids]
    with (BASE / 'okumura_example2_profiles.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['x'] + [f'u_t_{t:g}' for t in snap_t])
        for i in range(NX + 1):
            w.writerow([i / NX] + [fields[k][i] for k in ids])

    # Comparison with the digitized literature energy curve.
    paper = np.genfromtxt(DIGITIZED, delimiter=',', names=True)
    pt, pe = paper['time'], paper['paper_energy_digitized_from_Fig4_4']
    mask = (times >= 0.01) & (times <= 0.4)
    ref = np.interp(times[mask], pt, pe)
    err = energies[mask] - ref
    summary = {
        'frames': len(times),
        'time_range': f'{times[0]} .. {times[-1]}',
        'energy_initial': energies[0],
        'energy_final': energies[-1],
        'max_energy_increment': float(np.max(np.diff(energies))),
        'energy_monotone_nonincreasing': bool(np.all(np.diff(energies) <= 1e-12)),
        'max_y_variation_across_depth': max(y_spread),
        'u_final_min': float(np.min(fields[-1])),
        'u_final_max': float(np.max(fields[-1])),
        'digitized_comparison_interval': '0.01 <= t <= 0.4',
        'digitized_curve_RMSE': float(np.sqrt(np.mean(err ** 2))),
        'digitized_curve_MAE': float(np.mean(np.abs(err))),
        'digitized_curve_max_abs_error': float(np.max(np.abs(err))),
        'energy_scale_of_paper_figure': '0 .. 90',
    }
    with (BASE / 'okumura_example2_literature_comparison.csv').open('w', newline='') as f:
        w = csv.writer(f)
        for k, v in summary.items():
            w.writerow([k, v])
    # Per-time difference series for the error panel of Figure 1.
    with (BASE / 'okumura_example2_energy_difference.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time', 'OOFEM_energy', 'paper_energy_digitized', 'difference_OOFEM_minus_paper'])
        for t, e in zip(times[mask], energies[mask]):
            w.writerow([t, e, float(ref[np.argmin(np.abs(times[mask] - t))]), e - float(ref[np.argmin(np.abs(times[mask] - t))])])
    print(summary)


if __name__ == '__main__':
    main()
