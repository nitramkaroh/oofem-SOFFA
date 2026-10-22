#!/usr/bin/env python3
"""Create side-by-side source-paper/OOFEM energy-figure comparisons.

The OOFEM data are read from VTK outputs. Energy is evaluated on each P1
triangle with degree-four quadrature, using the Allen-Cahn energy from the
papers. Run from the repository root; results are saved in this script's
figures/ directory. Requires numpy and matplotlib.
"""
from pathlib import Path
import re
import textwrap
import argparse
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
DATA = Path("/tmp/ac_paperplots")

# Six-point, degree-four Dunavant rule. Weights sum to one; multiply by area.
a, b = 0.445948490915965, 0.108103018168070
c, d = 0.091576213509771, 0.816847572980459
BARY = np.array([[a,a,b],[a,b,a],[b,a,a],[c,c,d],[c,d,c],[d,c,c]])
WEIGHTS = np.array([0.223381589678011]*3 + [0.109951743655322]*3)


def vtk_data(path):
    root = ET.parse(path).getroot()
    arrays = {el.get("Name"): np.fromstring(el.text or "", sep=" ")
              for el in root.iter("DataArray") if el.get("Name")}
    # The coordinate array has no name.
    pts = next(el for el in root.iter("Points") if True)
    xyz = np.fromstring(next(iter(pts)).text or "", sep=" ").reshape(-1, 3)[:, :2]
    u = arrays["Temperature"]
    conn = arrays["connectivity"].astype(int)
    offsets = arrays["offsets"].astype(int)
    cells = np.vstack([conn[(0 if i == 0 else offsets[i-1]):off]
                       for i, off in enumerate(offsets)])
    if cells.shape[1] != 3:
        raise ValueError(f"Expected P1 triangles, got cells {cells.shape}")
    return xyz, cells, u


def energy(xy, cells, u, eps):
    tri = xy[cells]
    x0, x1, x2 = tri[:, 0], tri[:, 1], tri[:, 2]
    det = (x1[:,0]-x0[:,0])*(x2[:,1]-x0[:,1]) - (x2[:,0]-x0[:,0])*(x1[:,1]-x0[:,1])
    area = np.abs(det) * 0.5
    grads = np.empty((len(cells), 3, 2))
    grads[:,0,:] = np.column_stack((x1[:,1]-x2[:,1], x2[:,0]-x1[:,0])) / det[:,None]
    grads[:,1,:] = np.column_stack((x2[:,1]-x0[:,1], x0[:,0]-x2[:,0])) / det[:,None]
    grads[:,2,:] = np.column_stack((x0[:,1]-x1[:,1], x1[:,0]-x0[:,0])) / det[:,None]
    ue = u[cells]
    gu = np.einsum("ei,eij->ej", ue, grads)
    e_grad = 0.5 * np.sum(area * np.einsum("ei,ei->e", gu, gu))
    uq = ue @ BARY.T
    f = 0.25 * (uq*uq - 1.0)**2
    e_pot = (1/eps**2) * np.sum(area * (f @ WEIGHTS))
    return float(e_grad + e_pot)


def load_curve(case_dir, eps):
    case_dir = Path(case_dir)
    case_text = (case_dir / "case.in").read_text()
    initial = np.array([float(v) for v in re.findall(
        r"^InitialCondition\s+\d+.*?\bu\s+([-+0-9.eE]+)\s+set\s+\d+", case_text, re.M)])
    pvd = next(case_dir.glob("*.out.m0.pvd"))
    collection = ET.parse(pvd).getroot()
    sets = collection.findall(".//DataSet")
    times, files = [], []
    for ds in sets:
        times.append(float(ds.get("timestep")))
        files.append(case_dir / ds.get("file"))
    xy, cells, first_u = vtk_data(files[0])
    if len(initial) != len(xy):
        raise ValueError(f"Need {len(xy)} initial values, found {len(initial)} in {case_dir}")
    ts = [0.0] + times
    es = [energy(xy, cells, initial, eps)]
    kept_times = []
    for t, f in zip(times, files):
        x, tri, u = vtk_data(f)
        value = energy(x, tri, u, eps)
        if not np.isfinite(value):
            print(f"WARNING {case_dir}: non-finite energy at t={t:g}; keeping only the finite prefix")
            break
        kept_times.append(t)
        es.append(value)
    return np.array([0.0] + kept_times), np.array(es)


def draw_comparison(original_name, series, title, note, output_png):
    original = Image.open(FIG / original_name).convert("RGB")
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.2, 7.3), gridspec_kw={"width_ratios":[1,1.15]})
    fig.suptitle(title, fontsize=13, fontweight="bold", y=.985)
    left.imshow(original)
    left.set_title("Figura originale dell’articolo", fontsize=10)
    left.axis("off")
    for eps, t, e, label, style, color in series:
        if len(t) == 1:
            right.plot(t, e, marker="x", linestyle="none", color=color, markersize=7, label=label)
        else:
            right.plot(t, e, linestyle=style, color=color, linewidth=1.8, label=label)
    right.set_title("Energia calcolata da OOFEM", fontsize=10)
    right.set_xlabel("t")
    right.set_ylabel("E(u)")
    right.set_xlim(left=0)
    right.grid(True, alpha=.25)
    right.legend(frameon=True, fontsize=8)
    fig.text(.5, .02, textwrap.fill(note, width=145), ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0,.16,1,.94))
    fig.savefig(output_png, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA, help="directory containing wang/ and akrivis/ run folders")
    data_root = parser.parse_args().data_root
    wang = []
    for eps, dt, color in ((.1,.005,"#1f77b4"),(.05,.00125,"#ff7f0e"),(.01,.0001,"#e6ad00")):
        t,e = load_curve(data_root/"wang"/f"eps_{eps:g}", eps)
        suffix = (", solo t=0; NaN al primo passo" if len(t) == 1 else (f", si ferma a t={t[-1]:.3g}" if t[-1] < 0.999 else ""))
        wang.append((eps,t,e,fr"$\varepsilon={eps:g}$ (OOFEM $h=0.025$, $\Delta t={dt:g}$){suffix}","-",color))
        print(f"Wang eps={eps:g}: E0={e[0]:.8g}, Eend={e[-1]:.8g} at t={t[-1]:g}, max increase={(np.maximum(0,np.diff(e)).max() if len(e)>1 else float('nan')):.3g}")
    draw_comparison("wang_figure1_original.png", wang,
        "Wang, Wei & Chen (2025), Fig. 1 — stabilità dell’energia",
        "Paper: Ω=[0,1]², cos(2πx)cos(2πy), h=1/20, Δt=0.01; il testo indica T=20 ma l’asse della Fig. 1 arriva a t=1. "
        "OOFEM: mesh uniforme raffinata h=0.025 (P1), backward Euler. ε=0.1 e 0.05 completano t=1 con energia decrescente; ε=0.01 genera valori non finiti al primo passo. E(0)≈1616 è sopra il limite verticale 800 della figura del paper, quindi quel valore è tagliato nel grafico originale. Non è una replica dello schema del paper.",
        FIG/"wang_figure1_comparison.png")

    akrivis = []
    for eps, dt, color in ((.04,.0004,"black"),(.02,.0001,"red"),(.01,.0001,"blue")):
        t,e = load_curve(data_root/"akrivis"/f"eps_{eps:g}", eps)
        suffix = f", si ferma a t={t[-1]:.3g}" if t[-1] < .0999 else ""
        akrivis.append((eps,t,e,fr"$\varepsilon={eps:g}$ (OOFEM BE, $h=0.05$, $\Delta t={dt:g}$){suffix}","-",color))
        print(f"Akrivis eps={eps:g}: E0={e[0]:.8g}, Eend={e[-1]:.8g} at t={t[-1]:g}, max increase={(np.maximum(0,np.diff(e)).max() if len(e)>1 else float('nan')):.3g}")
    draw_comparison("akrivis_figure3_original.png", akrivis,
        "Akrivis & Li (2022), Fig. 3 — evoluzione dell’energia",
        "Paper: τ=ε², h=ε/2, BDF5 con avvio Radau IIA. OOFEM: backward Euler, mesh uniforme n=40 (h=0.05), doppio raffinamento rispetto ai run precedenti ma più grossolana della mesh del paper (h=ε/2). "
        "ε=0.04 e 0.02 completano t=0.1 con energia decrescente; ε=0.01 ha energia non monotona (incrementi locali fino a 4.01) e valori non finiti a t≈0.0221. Il confronto resta qualitativo: schema e mesh non coincidono.",
        FIG/"akrivis_figure3_comparison.png")

    with PdfPages(FIG/"confronto_grafici_originali_oofem.pdf") as pdf:
        for name in ("wang_figure1_comparison.png", "akrivis_figure3_comparison.png"):
            im = plt.imread(FIG/name)
            fig = plt.figure(figsize=(12.2,7.3))
            ax = fig.add_axes([0,0,1,1]); ax.imshow(im); ax.axis("off")
            pdf.savefig(fig); plt.close(fig)

if __name__ == "__main__":
    main()
