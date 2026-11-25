#!/usr/bin/env python3
"""Regenerate the Calatroni--Colli refinement figure from archived metrics."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "comparison_report_figures" / "calatroni_colli_refinement.png"

with (HERE / "spatial_convergence.csv").open() as f:
    spatial = list(csv.DictReader(f))
with (HERE / "calatroni_colli_temporal_refinement.csv").open() as f:
    temporal = list(csv.DictReader(f))

fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0))
axs[0].loglog(
    [float(r["h"]) for r in spatial[:2]],
    [float(r["nodal_rms_vs_40"]) for r in spatial[:2]],
    "-o", linewidth=2,
)
axs[0].set_xlabel("Passo spaziale h")
axs[0].set_ylabel("RMS nodale vs mesh 40×40")
axs[0].set_title("Raffinamento spaziale")
axs[0].annotate("riferimento: mesh 40×40\n(non soluzione esatta)", xy=(.06, 3e-5), fontsize=8)

axs[1].loglog(
    [float(r["dt_coarse"]) for r in temporal],
    [float(r["nodal_rms_vs_next_finer"]) for r in temporal],
    "-s", linewidth=2, color="#c44e52",
)
axs[1].set_xlabel("Passo temporale Δt")
axs[1].set_ylabel("RMS tra campi finali consecutivi")
axs[1].set_title("Raffinamento temporale")
axs[1].annotate("rapporto errori ≈ 1.97\nordine osservato ≈ 0.98", xy=(.0029, 1.38e-4), fontsize=8)
for ax in axs:
    ax.grid(True, which="both", alpha=.25)
fig.suptitle("Confronto numerico con la soluzione debole unica (test 2D parabolico)")
fig.tight_layout()
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=180, bbox_inches="tight")
