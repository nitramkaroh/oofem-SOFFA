# Final report — validation figures

Date: 2026-10-27. Directory: `doc/stvowu_validation/numerical_validation_figures/`.

## 1. Figures generated (3, all from existing data)

| Figure | Files | Content |
|---|---|---|
| **Fig. 1** — field comparison at t\* = T | `fig1_oofem_vs_dop853_fields_T.pdf` / `.png`, `make_fig1.py`, `data/fig1_fields_T0p02.csv` | 3 panels on the 32×16 nodes: (a) OOFEM solution, (b) DOP853 n128 semidiscrete reference restricted exactly to the same nodes, (c) \|u_h − u_ref\| on its own scale |
| **Fig. 2** — error vs mesh size at T | `fig2_convergence_FE_RMS_T.pdf` / `.png`, `make_fig2.py`, `data/fig2_convergence_data.csv` | FE RMS norms, log-log, bulk + boundary panels; OOFEM chain (h = 0.25/0.125/0.0625) and semidiscrete family level-differences (h = 0.0625/0.03125); measured local log₂ ratios annotated, **no slopes** |
| **Fig. 3** — boundary error in time | `fig3_boundary_error_time_series.pdf` / `.png`, `make_fig3.py`, `data/fig3_boundary_time_series.csv` | E_Γ(t) = ‖u_h − u_ref‖_{FE RMS, Γ} over all 80 steps, OOFEM n32 vs proxy128, plus the family 64→128 boundary curve as reference-solver context |

Shared read-only data layer: `figcommon.py`. Captions: `captions.md`.
Data sources and error definitions: `README.md`.

## 2. Verification performed (task §7)

- **Every plotted value comes from an identifiable file.** Raw sources: 3
  OOFEM VTU series (80 steps each), 4 DOP853 npz references; nothing else.
- **OOFEM and reference at the same time and physical setup.** All figures use
  step 80 / t = 0.02 (Fig. 1, 2) or the common 80-step time grid (Fig. 3);
  identical Ω, Γ, parameters, initial condition, Δt.
- **Consistent coordinates and mappings.** VTU files reordered by exact
  coordinate coincidence (asserted max distance 0.0); nested meshes coincide
  nodally, so all restrictions are exact index maps (asserted max coordinate
  distance 0.0 for every map used).
- **Error norms correct.** FE RMS norms with element-by-element mass matrices
  (3×3 Gauss in the bulk, 1-D P1 on the boundary); quadrature validated by
  constant-field invariants; the whole assembler reproduces the audited
  archive values (16/16 cross-checks in `t1b_metrics_T.csv`, floor
  2.916141·10⁻⁸ from `floors_recomputed.csv`) **before** any figure is written.
- **Regenerable.** `python3 make_fig1.py && python3 make_fig2.py &&
  python3 make_fig3.py` re-derives everything from the raw files, re-runs the
  assertions, and rewrites all PDF/PNG/CSV outputs.
- **No pre-existing file modified.** Only this directory was created; `git
  status` shows the new directory as untracked, nothing else touched.

## 3. Claims the figures DO support

1. **Agreement with an independent solver family (Figs. 1–3).** OOFEM and an
   independently integrated DOP853 solution of the same space-discretised
   system agree to FE RMS 3.8·10⁻⁴ (bulk) / 8.4·10⁻⁴ (boundary) at T on the
   32×16 mesh, with the difference concentrated at the two diagonal corners.
   The same-mesh OOFEM–DOP853 floor (2.9·10⁻⁸) is 7–8 orders below this, so
   the displayed error is spatial, not a solver/implementation discrepancy.
2. **Observed spatial convergence of both families (Fig. 2).** Error decreases
   monotonically with refinement for the OOFEM chain (h = 0.25 → 0.0625) and
   for the semidiscrete family (h = 0.0625 → 0.015625, factor ≈ 8.5 in the
   bulk at T), with no visible plateau in the available range.
3. **Small, stable boundary error in time (Fig. 3).** E_Γ(t) stays ≤ 1.2·10⁻³
   over the whole window and tracks the family level-difference shape, with
   no late-time drift or blow-up.

## 4. Claims the figures do NOT support (stated in captions where relevant)

1. **No asymptotic convergence order.** Annotated numbers are *local* log₂
   ratios of consecutive points at a single time level. The audited full-time
   analysis found 0/80 steps with both consecutive family ratios in [1.5, 2.5]
   (e.g. p_low(T) = 0.50 vs p_high(T) = 3.09); no order is claimed and no
   theoretical slope is drawn.
2. **No exact-PDE accuracy.** The DOP853 references are semidiscrete
   (space-discretised) solutions; the "proxy" is itself a discretisation.
   Accuracy against the exact PDE solution is NOT VERIFIABLE with the data
   available (no known exact solution).
3. **No angular-singularity proof.** The corner concentration visible in
   Fig. 1 is an *observation* (supported by the audited zone analysis, 65–75%
   of the error energy in the corner zones), not a demonstrated singularity.
4. **No universality.** Figures cover one parameter set, one final time
   T = 0.02 (plus the 80-step history), and meshes h ∈ [0.015625, 0.25]. They
   do not validate OOFEM in general.
5. **Measured vs post-processed, clearly separated.** All curves are computed
   quantities (mass-matrix norms, pointwise differences); the only
   post-processing step is the exact nodal restriction (zero interpolation
   error), which is documented in the captions and README.

## 5. Known limitations carried into the figures

- The OOFEM chain stops at h = 0.0625 (no OOFEM 64×32 run exists; the
  recommended next step is specified but **not** executed in
  `consolidated_final_report/audit_convergence_next/next_step.md`).
- The 3.09 family ratio at T is a single local ratio in a transient regime
  (the 64→128 difference is 3.2× smaller at T than at its step-44 peak); it
  is annotated as "single, at T" and not extrapolated.
- The n32 OOFEM field comes from the T1-A rerun because the original T1
  directory contains 0-byte VTU stubs for steps ≥ 40 (documented in
  `README.md` §2; rerun steps 1–39 are bit-identical to the original).
