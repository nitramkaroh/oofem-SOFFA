# Publication validation figures — OOFEM Allen–Cahn (coupled hyperbolic model, dynamic boundary)

Date: 2026-10-27. All figures in this directory are generated **exclusively from
pre-existing numerical data** in the `doc/stvowu_validation` dossier. No OOFEM
simulation was run, no reference integration was performed, no OOFEM code or
original script was modified. Every script in this directory is read-only with
respect to all pre-existing files; new outputs go only to this directory.

## Contents

| File | Role |
|---|---|
| `figcommon.py` | shared read-only data layer (VTU/npz readers, exact restriction maps, FE-RMS norms, archive cross-checks) |
| `make_fig1.py` | source of `fig1_oofem_vs_dop853_fields_T.pdf` / `.png` + `data/fig1_fields_T0p02.csv` |
| `make_fig2.py` | source of `fig2_convergence_FE_RMS_T.pdf` / `.png` + `data/fig2_convergence_data.csv` |
| `make_fig3.py` | source of `fig3_boundary_error_time_series.pdf` / `.png` + `data/fig3_boundary_time_series.csv` |
| `captions.md` | publication-ready English captions for the three figures |
| `final_report.md` | what the figures do and do not support |

## 1. Physical problem (identical for every curve in every figure)

Coupled hyperbolic Allen–Cahn on Ω = [0,2] × [0,1], Γ = ∂Ω:

- bulk:  ρ u_tt − k Δu = m (u − u³)  (Newmark β = 1/4, γ = 1/2)
- surface (first order in time) on Γ, with the surface variable equal to the
  **trace** u|_Γ (shared DOFs)
- parameters RHO = KB = MB = CS = KS = MS = 1
- T = 0.02 s, Δt = 2.5·10⁻⁴ s, 80 steps; initial condition
  u₀ = 0.2 + 0.1·cos(πx/2)cos(πy); homogeneous (no imposed) boundary data
- meshes: Q1 tensor grids with 2:1 aspect (nx = 2·ny), x-fast row-major
  numbering node(i,j) = j(nx+1)+i; nested levels coincide **nodally and
  exactly** (verified max coordinate distance 0.0), so restriction between
  levels is a pure index map with zero interpolation error

## 2. Data sources (all pre-existing)

OOFEM solutions (VTU export of the shared field `Temperature` = u, one file per
time step):

| mesh | h | source directory |
|---|---|---|
| 8×4   | 0.25    | `hyperbolic/coupled_independent_space_n8x4_dt0p00025/` |
| 16×8  | 0.125   | `hyperbolic/coupled_independent_space_n16x8_dt0p00025/` |
| 32×16 | 0.0625  | `hyperbolic/t1a_n32x16_rerun/oofem/` (T1-A rerun — see note) |

*Provenance note (n32):* in the original T1 directory
(`hyperbolic/coupled_independent_space_n32x16_dt0p00025/`) the VTU files for
steps ≥ 40 are **0-byte stubs** (step 39, 43 074 B, is the last valid one). The
T1-A rerun directory is the only valid n32 source; its steps 1–39 are
bit-identical to the original run and its input file is byte-identical
(`hyperbolic/t1a_n32x16_rerun/results/t1a_report.md`).

Semidiscrete reference solutions (independent DOP853 integrator of the
*space-discretised* systems, rtol = 2·10⁻¹², atol = 2·10⁻¹⁴, 81 steps,
arrays U, V at t = k·2.5·10⁻⁴):

| mesh | source file |
|---|---|
| 16×8  | `hyperbolic/coupled_independent_spatial_jacobian_fix/tolcheck_n16x8_rtol2e-12.npz` |
| 32×16 | `hyperbolic/coupled_independent_spatial_jacobian_fix/reference_n32x16_dt0p00025.npz` |
| 64×32 | `hyperbolic/coupled_independent_spatial_jacobian_fix/reference_n64x32_dt0p00025.npz` |
| 128×64 | `hyperbolic/t1b_dop853_n128x64/reference_n128x64_dt0p00025.npz` (T1-B "proxy") |

**These references are semidiscrete finite-element solutions, NOT exact PDE
solutions.** They are used as an independent-solver benchmark family; the
128×64 member serves as the error proxy. The OOFEM-vs-DOP853 difference on an
identical mesh (the "same-mesh floor") at T is 2.916141·10⁻⁸ in bulk FE RMS
(audited value, `consolidated_final_report/audit_convergence_next/results/floors_recomputed.csv`),
7–8 orders of magnitude below the displayed errors — so the displayed
error is the spatial discretisation error relative to the semidiscrete family.

## 3. Error definitions (all norms measured, not interpolated)

FE-RMS norms, accumulated element-by-element (no dense matrices):

- bulk:  ‖e‖ = sqrt( eᵀ M_b e / |Ω| ),  |Ω| = 2.0
- surface: ‖e‖ = sqrt( eᵀ M_s e / |Γ| ), |Γ| = 6.0

M_b by 3×3 Gauss quadrature per Q1 quad (exact for the quartic integrands),
M_s as 1-D P1 line elements on the four boundary edges. The assembler and the
quadrature were independently re-verified by the 2026-10-27 audit
(`consolidated_final_report/audit_convergence_next/`, 172/172 checks, matching
the original `check_jacobian_fix.assemble_fixed` to machine precision).

Quantities plotted:

- **OOFEM chain** (blue): E_h = ‖u_h − U₁₂₈|_h‖, h = 0.25 / 0.125 / 0.0625,
  with U₁₂₈|_h the exact index restriction of the DOP853 128×64 solution to the
  coarser nodes.
- **Family level-differences** (vermilion): ‖U₃₂ − U₆₄|₃₂‖ evaluated on the
  32×16 mesh (h = 0.0625) and ‖U₆₄ − U₁₂₈|₆₄‖ evaluated on the 64×32 mesh
  (h = 0.03125) — h is the mesh size of the *evaluation* mesh (audited
  convention).
- Pointwise max convention (Fig. 1 only, stated in the caption): interior max
  over non-boundary nodes, boundary max over boundary nodes.

## 4. Cross-checks against the audited archive

Every `make_fig*.py` re-derives its plotted values from the raw VTU/npz files
and asserts agreement (relative 2·10⁻⁶, in practice exact to all printed
digits) with the audited archive **before writing any output**:

- `hyperbolic/t1b_dop853_n128x64/results/t1b_metrics_T.csv` — all OOFEM-vs-proxy
  and family level-difference norms at T (16 values used across the figures)
- `consolidated_final_report/audit_convergence_next/results/floors_recomputed.csv`
  — same-mesh floor n32@T

`figcommon.py` additionally self-checks the quadrature (constant-field
invariants ‖1‖²_bulk = |Ω|, ‖1‖²_surf = |Γ| on 8×4, 32×16, 64×32) and the exact
nodal coincidence of every restriction map used.

## 5. Reproducibility

```
cd doc/stvowu_validation/numerical_validation_figures
python3 make_fig1.py
python3 make_fig2.py
python3 make_fig3.py
```

Environment used: Python 3.12.3, NumPy 1.26.4, matplotlib 3.6.3 (Agg backend).
Each run re-reads the raw sources, re-runs all archive assertions, rewrites the
PDF/PNG pair and the `data/` CSV for that figure. Runtime is a few seconds
(fig3 re-reads 80 small VTU files). No network, no packages, no writes outside
this directory.

## 6. Style

- vector PDF + 300 dpi PNG for each figure
- Okabe–Ito blue `#0072B2` = OOFEM; Okabe–Ito vermilion `#D55E00` =
  reference-solver family (CVD-safe pair, Oklab ΔE×100 = 32.7 ≥ 8)
- sequential fields on viridis (perceptually uniform); panels (a)/(b) of Fig. 1
  share one colour scale, panel (c) has its own 0–max scale
- no theoretical slopes drawn anywhere; only measured local log₂ ratios are
  annotated, and only as local ratios (see `final_report.md`)
