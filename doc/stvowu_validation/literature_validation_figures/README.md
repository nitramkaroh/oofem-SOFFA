# Literature-comparison validation figures — OOFEM Allen–Cahn (coupled hyperbolic model)

Date: 2026-10-27. Directory: `doc/stvowu_validation/literature_validation_figures/`.

Three publication-ready figures for the thesis numerical-validation section
("validate the code with results available in the literature for evolution
equations with dynamic boundary conditions — Kovács–Lubich, Hipp–Kovács,
Hochbruck–Leibold"). Each figure **superimposes the OOFEM solution on the
literature solution** so that their coincidence is visible directly, plus a
pointwise-difference panel that quantifies it. **Every plotted value comes
from pre-existing archived OOFEM data** in `doc/stvowu_validation/hyperbolic/`.
No OOFEM simulation was run, no reference integration was performed, no
pre-existing file was modified. This directory is deliberately kept
**separate** from `numerical_validation_figures/` (OOFEM-vs-DOP853
semidiscrete comparison) — see
`../literature_analysis_2026/06_scope_separation_dop853.md`.

## Figures (primary deliverable — OOFEM vs literature, superposed)

| Figure | Files | Content |
|---|---|---|
| **Fig. 1** — OOFEM field vs the Hochbruck–Leibold exact profile | `cmp_fig1_hl_profile_field.pdf` / `.png`, `make_cmp_fig1.py` | (a) OOFEM solution at T on the finest polygonal unit disk (502 nodes, h ≈ 0.09, Δt = 2.5·10⁻³) as a shaded field, with the published profile u = sin(2πT)·xy (Hochbruck & Leibold 2021, §4.5) drawn as black contours; the two coincide to within the shading resolution. (b) pointwise difference u_h − u_exact: the interior is essentially zero, the deviation is concentrated at the polygonal boundary; max error 4.88·10⁻², RMS 1.42·10⁻². |
| **Fig. 2** — OOFEM curves vs the Nizovtseva exact traveling wave | `cmp_fig2_nizovtseva_wave.pdf` / `.png`, `make_cmp_fig2.py`, `data_cmp_fig2_nizovtseva_wave.csv` | (a) the exact solution u(x,T) = 1/(1+exp(λx − 1.5T)), λ = √(2+9τ)/2, τ = 0.2 (Nizovtseva, Galenko & Alexandrov 2016, Eq. (29)) as a black dashed curve; the three archived OOFEM solutions (n_x = 80/160/320, Δt = 10⁻³, T = 0.2, strip [−12,12]×[0,1]) lie on it within the plot resolution. (b) pointwise differences: max error 1.16·10⁻⁴ (coarsest) → 7.19·10⁻⁶ (finest). |
| **Fig. 3** — OOFEM curves vs the stationary hyperbolic Allen–Cahn front | `cmp_fig3_kink_profile.pdf` / `.png`, `make_cmp_fig3.py`, `data_cmp_fig3_kink_profile.csv` | (a) the exact front u(x) = tanh(x/√2) of u_tt − Δu = u − u³ as a black dashed curve; the three archived OOFEM solutions (N_x = 40/80/160, Δt = 10⁻³, T = 0.2, domain [−8,8]) lie on it. (b) pointwise differences: max error 7.27·10⁻⁵ (coarsest) → 6.59·10⁻⁶ (finest). |

Shared read-only data layer: `figcommon_lit.py`. Captions: `captions.md`.

## Supplementary (convergence-view of the same archived runs)

| Figure | Files | Content |
|---|---|---|
| S1 — HL exact profile, time refinement | `lit_fig1_hl_profile_temporal.pdf` / `.png`, `make_lit_fig1.py`, `data_fig1_hl_profile_temporal.csv` | OOFEM on the published profile (fixed 179-node mesh): error pinned by the spatial floor; successive-refinement differences show measured local temporal order **2.02** |
| S2 — HL exact profile, space refinement | `lit_fig2_hl_profile_spatial.pdf` / `.png`, `make_lit_fig2.py`, `data_fig2_hl_profile_spatial.csv` | Polygonal unit disks (60/179/502 nodes): measured local spatial ratios **0.33 / 0.24** (no order claimed) and a flat maximum error ≈ 4.9·10⁻² — the boundary-approximation (geometric) floor |
| S3 — bulk code verification, 2×2 | `lit_fig3_bulk_exact_solutions.pdf` / `.png`, `make_lit_fig3.py`, `data_fig3_bulk_exact_solutions.csv` | (a) Nizovtseva wave, space, orders 1.995/1.976; (b) tanh front, space, 2.006/1.887; (c) linear cosine mode vs exact semidiscrete, time, 1.997; (d) Allen–Cahn oscillator vs independent ODE, time, 2.016 |

## Data provenance

| Figure | Raw sources (all pre-existing) |
|---|---|
| 1 | Final VTU field: `hyperbolic/hl_disk_space_h0.09_dt0p0025/case.out.m0.*.vtu` (unit disk, 502 nodes, T = 0.2) |
| 2 | 3 final VTU fields: `hyperbolic/nizovtseva_space_n{80,160,320}_dt0p001/case.out.m0.*.vtu` |
| 3 | 3 final VTU fields: `hyperbolic/kink_n{40,80,160}/case.out.m0.*.vtu` |
| S1 | 3 final VTU fields: `hyperbolic/hl_disk_time_h016_dt{0.01,0.005,0.0025}/case.out.m0.*.vtu` (179 nodes) |
| S2 | 3 final VTU fields: `hyperbolic/hl_disk_space_h{0.3,0.16,0.09}_dt0p0025/case.out.m0.*.vtu` |
| S3(a,b) | same VTU fields as Figs. 2–3 |
| S3(c,d) | `hyperbolic/bulk_reference_metrics.csv` (source of record; the linear-mode reference is the exact **semidiscrete** Q1 solution, the oscillator reference an independent ODE integration — neither is re-derived here) |

**Re-derivation and cross-check.** Every script re-reads the final VTU of
each run, recomputes the nodal RMS and maximum error against the published
exact profile, and **asserts agreement (relative 2·10⁻⁶) with the archived
metrics CSV** before writing any output: `hochbruck_leibold_disk_metrics.csv`,
`nizovtseva_traveling_wave_metrics.csv`,
`stationary_kink_space_refinement.csv`. Regenerate the primary set with:

```
cd doc/stvowu_validation/literature_validation_figures
python3 make_cmp_fig1.py && python3 make_cmp_fig2.py && python3 make_cmp_fig3.py
```

(…and `make_lit_fig1.py`, `make_lit_fig2.py`, `make_lit_fig3.py` for the
supplementary set.)

Environment: Python 3.12.3, NumPy 1.26.4, matplotlib 3.6.3 (Agg). No network,
no packages, no writes outside this directory.

## What these figures DO show

1. **Figs. 2 and 3 are direct comparisons of OOFEM's output with exact
   solutions published in the literature.** The OOFEM curves coincide with
   the literature solution within the plot resolution; the difference panels
   quantify the coincidence (≤ 1.2·10⁻⁴ on the coarsest mesh, ≤ 7.2·10⁻⁶ on
   the finest) and show that the discrepancy shrinks with mesh refinement
   (measured local spatial orders 1.995/1.976 for the wave and 2.006/1.887
   for the front — see S3 for the convergence view).
2. **Fig. 1 shows OOFEM's field reproducing the Hochbruck–Leibold test
   profile** u(t,x) = sin(2πt)·xy at T = 0.2: the shaded OOFEM field and the
   literature contours are indistinguishable in the interior; the difference
   is a boundary layer of the polygonal disk (max 4.88·10⁻²), which is the
   boundary-approximation effect the literature flags for non-polygonal
   domains (Kovács & Lubich 2017, Introduction, 3rd remark; Hipp & Kovács
   2021, Sec. 2).
3. **The behavior is consistent with the cited works**: second-order spatial
   behavior on bulk 1-D domains (same O(h²)-class as Hipp–Kovács 2021,
   Table 1, pure second-order column, and Hochbruck–Leibold 2020, Fig. 5.1),
   and second-order time behavior on the fixed mesh (S1, local order 2.02 —
   the same order class as the second-order IMEX of Hochbruck–Leibold 2021
   and the Gauss–RK order-2s schemes of Hipp–Kovács 2021, Thms. 7.1–7.2).

## What these figures do NOT show (stated in the captions)

1. **No reproduction of a P1–P4 benchmark of OOFEM's exact coupled model.**
   None of the four literature papers contains a numerical benchmark of the
   OOFEM model (hyperbolic nonlinear bulk + first-order-in-time dynamic
   surface with Wentzell flux coupling) — see the compatibility analysis in
   `../literature_analysis_2026/03_compatibility_matrix.md`.
2. **Fig. 1 is an ADAPTED exact-profile test.** The forcing in the `hl_disk`
   runs was manufactured for OOFEM's own operator (bulk-inertial hyperbolic +
   first-order surface, `wavedamping 1.0`, `manufactured 1`); the OOFEM
   surface law is **first order** in time, while the Hochbruck–Leibold
   boundary equation is **second order** (kinetic/acoustic). The figure
   verifies that OOFEM reproduces the paper's test profile under the adapted
   operator; it is not a reproduction of the paper's boundary benchmark.
   The P4 profile (e^{−t}x₁²x₂²) has not been implemented in OOFEM.
3. **Nizovtseva et al. (2016) is an additional reference** (Fig. 2), not part
   of the P1–P4 corpus of the literature analysis; the thesis must add the
   citation. The traveling wave is an exact solution of the linearly-damped
   hyperbolic Allen–Cahn bulk equation (no dynamic boundary condition); it
   verifies the bulk machinery, not the surface law.
4. **No asymptotic order is claimed in the primary figures.** The
   convergence-set figures (S1–S3) annotate measured local log₂ ratios only;
   S2 explicitly shows a reduced (non-converged-in-rate) spatial behavior on
   the polygonal disk.
5. **P1 (Kovács–Lubich 2017) contributes theory only** (parabolic problems,
   no numerical experiments at all — verified page by page); it is used here
   only for the boundary-approximation remark and as parabolic context.
6. **The OOFEM-vs-DOP853 semidiscrete comparison is a separate comparison**
   (Rule 8) and is not part of these literature figures: see
   `../numerical_validation_figures/`.

## Style

- vector PDF + 300 dpi PNG; Okabe–Ito blue `#0072B2` for OOFEM data, black
  dashed for the literature solution, symmetric diverging scale for
  difference panels; measured quantities only in the titles.
