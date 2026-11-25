#!/usr/bin/env python3
"""Shared read-only data layer for the literature-comparison figures.

All data come from pre-existing OOFEM run directories in
`doc/stvowu_validation/hyperbolic/`.  For every run used by figures 1, 2
and 3(a,b) the final VTU field is re-read and the error against the
published exact profile is recomputed, then asserted to agree with the
archived metrics CSV (relative 2e-6).  The two bulk-reference panels
3(c,d) (linear cosine mode, Allen-Cahn oscillator) use the archived CSV
`bulk_reference_metrics.csv` directly: their references (the exact
semidiscrete Q1 solution and an independent ODE integration of
u_tt = u - u^3, respectively) are not re-derived here.

Read-only w.r.t. every pre-existing file; writes only to this directory.
"""
import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

HYPER = Path(__file__).resolve().parents[1] / 'hyperbolic'
T = 0.2                          # final time of every archived run below
TAU = 0.2                        # Nizovtseva relaxation time
LAMBDA = math.sqrt(2.0 + 9.0 * TAU) / 2.0   # Nizovtseva et al. (2016), Eq. (29)


# ------------------------------------------------------------- VTU layer --
def final_field(folder):
    """(coords (N,2), u) from the last VTU of folder/case.out.m0.pvd."""
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    vtu = folder / pvd.findall('.//DataSet')[-1].attrib['file']
    root = ET.parse(vtu).getroot()
    coords = np.fromstring(root.find('.//Points/DataArray').text,
                           sep=' ').reshape(-1, 3)[:, :2]
    arr = next(a for a in root.findall('.//PointData/DataArray')
               if a.attrib.get('Name') == 'Temperature')
    u = np.fromstring(arr.text or '', sep=' ')
    return coords, u


def rms_max(coords, u, exact_fn):
    e = u - exact_fn(coords)
    return float(np.sqrt(np.mean(e ** 2))), float(np.max(np.abs(e)))


def final_tri(folder):
    """(coords (N,2), u, triangles (M,3) int) from the last VTU of a triangular mesh."""
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    vtu = folder / pvd.findall('.//DataSet')[-1].attrib['file']
    root = ET.parse(vtu).getroot()
    coords = np.fromstring(root.find('.//Points/DataArray').text,
                           sep=' ').reshape(-1, 3)[:, :2]
    arr = next(a for a in root.findall('.//PointData/DataArray')
               if a.attrib.get('Name') == 'Temperature')
    u = np.fromstring(arr.text or '', sep=' ')
    # OOFEM VTU layout: UnstructuredGrid/Piece/Cells with named DataArrays.
    # OOFEM writes M "cumulative" offset entries (offsets[-1] == len(conn)):
    # cell j (0-based) is conn[start[j]:start[j+1]] with start = [0]+offsets.
    # The mesh also contains boundary line elements (surface model) — keep
    # triangles only (VTK_TRIANGLE = 5).
    cells_el = root.find('.//UnstructuredGrid/Piece/Cells')
    arrays = {a.attrib.get('Name'): a for a in cells_el.findall('DataArray')}
    types = np.fromstring(arrays['types'].text, sep=' ', dtype=np.int64)
    offs = np.fromstring(arrays['offsets'].text, sep=' ', dtype=np.int64)
    cells = np.fromstring(arrays['connectivity'].text, sep=' ', dtype=np.int64)
    start = np.concatenate([[0], offs])
    assert start[-1] == len(cells), (start[-1], len(cells))
    tris = [tuple(cells[start[j]:start[j + 1]]) for j in range(len(offs))
            if types[j] == 5]
    return coords, u, np.array(tris, dtype=int)


def check_field(d, exact_fn, csv_name, key_col, key_val,
                rms_col='rms_error_vs_exact_profile', max_col=None):
    """Read one archived run, assert the error against the metrics CSV, return (coords, u)."""
    coords, u = final_field(HYPER / d)
    rms, mx = rms_max(coords, u, exact_fn)
    rows = {r[key_col]: r for r in _csv_rows(csv_name)}
    _check(f'{d} rms', rms, float(rows[key_val][rms_col]))
    if max_col:
        _check(f'{d} max', mx, float(rows[key_val][max_col]))
    return coords, u


# --------------------------------------------------- exact profiles (T) --
def hl_exact(coords):
    """Hochbruck & Leibold (2021), Sec. 4.5 test profile, at t = T."""
    return math.sin(2.0 * math.pi * T) * coords[:, 0] * coords[:, 1]


def niz_exact(coords):
    """Nizovtseva, Galenko & Alexandrov (2016), Eq. (29); Q=A=B=+1, c0=0, at t = T."""
    z = np.clip(LAMBDA * coords[:, 0] - 1.5 * T, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(z))


def kink_exact(coords):
    """Stationary hyperbolic Allen-Cahn front tanh(x/sqrt(2))."""
    return np.tanh(coords[:, 0] / math.sqrt(2.0))


# -------------------------------------------------- archived metric rows --
def _csv_rows(path):
    with open(HYPER / path, newline='') as f:
        return list(csv.DictReader(f))


def _check(name, recomputed, archived, rtol=2e-6):
    ok = abs(recomputed - archived) <= rtol * max(abs(archived), 1e-300)
    print(f'  check {name:34s} {recomputed:.9e}  (csv {archived:.9e})  '
          f'{"OK" if ok else "MISMATCH"}')
    assert ok, name


# Figure 1: HL profile, time study (fixed mesh h ~ 0.16, 179 nodes).
def hl_time_study():
    """[(dt, rms, max, diff_to_next)] recomputed from VTU, checked vs CSV."""
    rows = {r['tag']: r for r in _csv_rows('hochbruck_leibold_disk_metrics.csv')}
    out = []
    # the CSV tags are: time_h016_dt0.01 / dt0.005 / dt0.0025
    tag_of = {'0.01': 'time_h016_dt0.01', '0.005': 'time_h016_dt0.005',
              '0.0025': 'time_h016_dt0.0025'}
    for dt, d in ((0.01, 'hl_disk_time_h016_dt0.01'),
                  (0.005, 'hl_disk_time_h016_dt0.005'),
                  (0.0025, 'hl_disk_time_h016_dt0.0025')):
        tag = tag_of[f'{dt:g}']
        coords, u = final_field(HYPER / d)
        rms, mx = rms_max(coords, u, hl_exact)
        r = rows[tag]
        _check(f'hl time dt={dt:g} rms', rms, float(r['rms_error_vs_exact_profile']))
        _check(f'hl time dt={dt:g} max', mx, float(r['max_error_vs_exact_profile']))
        out.append((dt, rms, mx,
                    float(r['rms_to_next_finer']) if r['rms_to_next_finer'] else None))
    return out


# Figure 2: HL profile, space study (dt = 0.0025, polygonal unit disk).
def hl_space_study():
    rows = {r['tag']: r for r in _csv_rows('hochbruck_leibold_disk_metrics.csv')}
    out = []
    for d, tag, h, nodes in (('hl_disk_space_h0.3_dt0p0025', 'space_h0.3_dt0p0025', 0.3, 60),
                             ('hl_disk_space_h0.16_dt0p0025', 'space_h0.16_dt0p0025', 0.16, 179),
                             ('hl_disk_space_h0.09_dt0p0025', 'space_h0.09_dt0p0025', 0.09, 502)):
        coords, u = final_field(HYPER / d)
        rms, mx = rms_max(coords, u, hl_exact)
        assert len(coords) == nodes, (d, len(coords), nodes)
        r = rows[tag]
        _check(f'hl space h={h:g} rms', rms, float(r['rms_error_vs_exact_profile']))
        _check(f'hl space h={h:g} max', mx, float(r['max_error_vs_exact_profile']))
        out.append((h, nodes, rms, mx))
    return out


# Figure 3a: Nizovtseva traveling wave, space study (dt = 0.001).
def niz_space_study():
    rows = {r['tag']: r for r in _csv_rows('nizovtseva_traveling_wave_metrics.csv')}
    out = []
    for d, tag, h, nx in (('nizovtseva_space_n80_dt0p001', 'space_n80_dt0p001', 0.3, 80),
                          ('nizovtseva_space_n160_dt0p001', 'space_n160_dt0p001', 0.15, 160),
                          ('nizovtseva_space_n320_dt0p001', 'space_n320_dt0p001', 0.075, 320)):
        coords, u = final_field(HYPER / d)
        rms, mx = rms_max(coords, u, niz_exact)
        assert len(coords) == (nx + 1) * 2, (d, len(coords))
        r = rows[tag]
        _check(f'niz space h={h:g} rms', rms, float(r['rms_error_vs_exact']))
        _check(f'niz space h={h:g} max', mx, float(r['max_error_vs_exact']))
        out.append((h, rms, mx))
    return out


# Figure 3b: stationary kink, space study (T = 0.2, dt = 0.001).
def kink_space_study():
    rows = {r['nx']: r for r in _csv_rows('stationary_kink_space_refinement.csv')}
    out = []
    for d, nx, h in (('kink_n40', 40, 0.4), ('kink_n80', 80, 0.2), ('kink_n160', 160, 0.1)):
        coords, u = final_field(HYPER / d)
        rms, mx = rms_max(coords, u, kink_exact)
        r = rows[str(nx)]
        _check(f'kink nx={nx} rms', rms, float(r['rms_error']))
        _check(f'kink nx={nx} max', mx, float(r['max_error']))
        out.append((h, rms, mx))
    return out


# Figure 3c/d: bulk references (archived CSV is the source of record).
def bulk_reference_rows():
    """{case: [(dt, rms, observed_order_or_None)]} from bulk_reference_metrics.csv."""
    rows = _csv_rows('bulk_reference_metrics.csv')
    out = {}
    for case in ('linear', 'oscillator'):
        sel = sorted((r for r in rows if r['case'] == case),
                     key=lambda r: float(r['dt']), reverse=True)
        assert [float(r['dt']) for r in sel] == [0.04, 0.02, 0.01], sel
        out[case] = [(float(r['dt']), float(r['rms_error']),
                      float(r['observed_order']) if r['observed_order'] else None)
                     for r in sel]
    return out


# ------------------------------------------------------------ plot style --
C_OOF = '#0072B2'      # Okabe-Ito blue   — OOFEM data
C_FLOOR = '#999999'    # gray — reference lines (spatial floors)
ABOX = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85)


def log2(a, b):
    return float(math.log2(a / b))


def style_fig():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.size': 8.5, 'axes.linewidth': 0.8, 'axes.titlesize': 9,
        'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
        'figure.dpi': 150, 'savefig.dpi': 300,
    })
    return plt


def annotate_order(ax, x1, x2, y1, y2, label):
    ax.annotate(label, xy=((x1 + x2) / 2, math.sqrt(y1 * y2)),
                ha='center', va='bottom', fontsize=7.5, color='0.2', bbox=ABOX)
