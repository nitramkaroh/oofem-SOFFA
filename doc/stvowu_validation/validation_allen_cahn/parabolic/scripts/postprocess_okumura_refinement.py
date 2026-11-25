#!/usr/bin/env python3
"""Postprocess the Okumura spatial refinement runs (n = 50, 100, 200, 400).

The nodal grids nest (x_i^50 = x_{2i}^100 = x_{4i}^200 = x_{8i}^400), so
each coarser solution is compared against the finest (n = 400) solution
nodally, without any interpolation.

Export precision. The OOFEM VTU writer prints point data with the C++
default stream precision (6 significant digits; vtkxmlexportmodule.C,
writeVTKPointData: fileStream << scientific << value). For the O(1) field
u the export quantum is therefore ~1e-6: mesh-to-mesh differences below
~5e-6 (a few quantization steps) are export rounding, not solution
differences. Observed convergence orders are consequently reported only
where the differences clear that floor (t = 0.3, 0.5, interface in
transit); at t = 0.7 (near the final state) the meshes agree within the
export precision and the row is reported as an agreement check without
an order. No data are filtered or interpolated; the guard only withholds
an order estimate from numbers that are below the resolution of the
exported values.

Output: ../okumura_refinement_convergence.csv
"""
from __future__ import annotations
import csv, math
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RESULTS = BASE / 'results'
DT, T = 0.001, 0.7
T_EVAL = (0.3, 0.5)      # rows where an observed order may be assigned
T_CHECK = (0.7,)         # agreement-check rows only
MESHES = (50, 100, 200, 400)
REF = 400
QUANTUM = 1e-6           # 6 significant digits on O(1) values
ORDER_FLOOR = 5 * QUANTUM  # below this, a row cannot support an order


def read_bottom_row(path: Path, nx: int) -> np.ndarray:
    piece = ET.parse(path).getroot().find('.//Piece')
    xyz = np.fromstring(piece.find('./Points/DataArray').text, sep=' ').reshape(-1, 3)
    vals = np.fromstring(piece.find('./PointData/DataArray').text, sep=' ')
    u = np.empty(nx + 1)
    for (x, y), v in zip(xyz[:, :2], vals):
        if abs(y) < 0.5:
            u[round(x * nx)] = v
    return u


def frame_path(nx: int, t: float) -> Path:
    k = int(round(t / DT))
    return RESULTS / f'okumura_refinement_n{nx}' / f'case.out.m0.{k}.vtu'


def main() -> None:
    fields = {}
    for nx in MESHES:
        for t in T_EVAL + T_CHECK:
            fields[(nx, t)] = read_bottom_row(frame_path(nx, t), nx)

    rows = []
    for t in T_EVAL + T_CHECK:
        ref = fields[(REF, t)]
        prev_rms = None
        prev_ok = False
        for nx in MESHES:
            u = fields[(nx, t)]
            if nx == REF:
                rms = mx = 0.0
                note = 'reference (finest) mesh'
            else:
                err = u - ref[::(REF // nx)]
                rms = float(np.sqrt(np.mean(err ** 2)))
                mx = float(np.max(np.abs(err)))
                if mx < 5 * QUANTUM:
                    note = ('differences within export precision '
                            '(6 significant digits); no order assigned')
                elif mx < ORDER_FLOOR:
                    note = ''
                else:
                    note = ''
            order = ''
            if note == '' and prev_ok and prev_rms >= ORDER_FLOOR:
                order = f'{math.log(prev_rms / rms, 2):.3f}'
            rows.append({'nx': nx, 'h': 1.0 / nx, 't_eval': t,
                         'rms_vs_n400_ref': rms, 'max_vs_n400_ref': mx,
                         'observed_order': order, 'note': note})
            prev_rms = rms
            prev_ok = (note == '')
    out = BASE / 'okumura_refinement_convergence.csv'
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print('wrote', out)


if __name__ == '__main__':
    main()
