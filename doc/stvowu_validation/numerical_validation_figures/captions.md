# Captions (publication-ready English)

Conventions: u_h denotes the OOFEM solution; U_N denotes the DOP853
semidiscrete solution on the N×N/2 tensor mesh; "FE RMS" denotes the
finite-element L²-type RMS norm with the mass matrix of the evaluation mesh
(‖e‖ = sqrt(eᵀMe/|Ω|) in the bulk, sqrt(eᵀMe/|Γ|) on the boundary).

---

**Figure 1 — Field comparison at t\* = T = 0.02 (step 80).**
*(a) OOFEM solution u_h on the 32×16 mesh. (b) Semidiscrete DOP853
reference U₁₂₈ (128×64 mesh, independent integrator, rtol = 2·10⁻¹²)
restricted to the 32×16 nodes by exact nodal coincidence (no interpolation).
(c) Pointwise absolute difference |u_h − U₁₂₈| on the 32×16 nodes, on a
separate scale. Panels (a) and (b) share one colour scale. Problem: coupled
hyperbolic Allen–Cahn on Ω = [0,2]×[0,1] with trace-coupled boundary,
homogeneous initial-boundary data except the prescribed initial condition;
identical parameters, mesh topology and time level in both solvers. The
reference is a semidiscrete (space-discretised) solution, not an exact PDE
solution; the same-mesh OOFEM–DOP853 difference at T is 2.9·10⁻⁸ in bulk FE
RMS, so the displayed error is the spatial error of OOFEM relative to the
semidiscrete family. The difference is concentrated at the two diagonal
corners (0,1) and (2,0). Bulk FE RMS of the difference: 3.82·10⁻⁴; boundary
FE RMS: 8.37·10⁻⁴.*

---

**Figure 2 — Error versus mesh size at t\* = T = 0.02 (FE RMS norms).**
*Blue: OOFEM error E_h = ‖u_h − U₁₂₈|_h‖ for h = 0.25, 0.125, 0.0625
(8×4, 16×8, 32×16 meshes), with U₁₂₈ restricted exactly to the coarser nodes.
Vermilion: level-differences of the semidiscrete DOP853 family,
‖U₃₂ − U₆₄|₃₂‖ (evaluated on the 32×16 mesh, h = 0.0625) and
‖U₆₄ − U₁₂₈|₆₄‖ (evaluated on the 64×32 mesh, h = 0.03125); h denotes the
mesh size of the evaluation mesh. Both panels: same time level and physical
problem; left panel — bulk Ω, right panel — boundary Γ. Annotated numbers are
measured local log₂ ratios of consecutive points at T only; no asymptotic
convergence order is claimed or indicated by any slope, because the
consecutive local ratios (0.72 and 1.05 for the OOFEM chain; 3.09 for the
family in the bulk) are not mutually consistent with a stable order.*

---

**Figure 3 — Boundary error in time, E_Γ(t) = ‖u_h(t) − U(t)‖_{FE RMS, Γ}.**
*Blue: OOFEM 32×16 solution against the DOP853 128×64 semidiscrete reference
restricted to the 32×16 boundary nodes, all 80 time steps
(t = k·2.5·10⁻⁴, k = 1…80). Vermilion (dashed): semidiscrete family
level-difference ‖U₆₄ − U₁₂₈|₆₄‖ on the 64×32 boundary, shown as
reference-solver context. The surface variable of the model is the trace
u|_Γ, so the boundary FE RMS norm is the model-consistent boundary error
norm. Both curves grow from ~2·10⁻⁶ at the first step (t = 2.5·10⁻⁴), peak
between t ≈ 0.010 and 0.013, and remain ≤ 1.2·10⁻³ over the whole window.
The OOFEM curve and the family curve have qualitatively the same shape; the
OOFEM error is at most about 8× the family level-difference at any time
after the initial growth (2.9× at the OOFEM peak, t ≈ 0.013; 5.1× at T). The
family curve is non-monotone near the end of the window (peak at t ≈ 0.010,
minimum near t ≈ 0.0185), reflecting transient cancellation in the
level-difference rather than an error floor.*
