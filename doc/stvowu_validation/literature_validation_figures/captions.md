# Captions (publication-ready English)

Conventions: u_h denotes the OOFEM solution; "RMS" denotes the nodal
root-mean-square error over all mesh nodes at the final time T = 0.2 s;
"local order" denotes a measured log₂ ratio of two consecutive error values
— it is not a demonstrated asymptotic convergence rate.

---

## Primary figures (OOFEM vs literature, superposed)

**Figure 1 — OOFEM solution and the Hochbruck–Leibold exact profile at T = 0.2 (polygonal unit disk).**
*The published test profile u(t,x) = sin(2πt)·xy of Hochbruck & Leibold
(2021, Sec. 4.5) is imposed on the finest archived mesh (502 nodes,
h ≈ 0.09, Δt = 2.5·10⁻³ s) by manufacturing the source terms of OOFEM's own
operator (bulk-inertial hyperbolic Allen–Cahn with a first-order-in-time
dynamic surface law); the profile is therefore an exact solution of the
*adapted* problem — OOFEM's surface law is first order in time, whereas the
paper's boundary equation is second order — so the figure verifies that
OOFEM reproduces the paper's test profile under the adapted operator, not
the paper's boundary benchmark. (a) The OOFEM field (shaded) and the
literature profile (black contours) coincide to within the shading
resolution in the interior of the domain. (b) Pointwise difference
u_h − u_exact: the deviation is confined to the polygonal boundary
(max 4.88·10⁻², RMS 1.42·10⁻²), the boundary-approximation effect emphasized
for non-polygonal domains in the bulk–surface literature (Kovács & Lubich
2017, Introduction, 3rd remark; Hipp & Kovács 2021, Sec. 2).*

---

**Figure 2 — OOFEM solution and the exact traveling wave of Nizovtseva, Galenko & Alexandrov (2016) at T = 0.2.**
*Exact solution u(x,t) = 1/(1+exp(λx − 1.5t)) of the damped hyperbolic
Allen–Cahn equation τ u_tt + u_t − Δu = u − u³ (τ = 0.2), Eq. (29) of
Nizovtseva et al. (2016) (Q = A = B = +1, c₀ = 0), at T = 0.2 s on the strip
[−12,12]×[0,1]. (a) Black dashed curve: the literature solution; blue
curves: the three archived OOFEM solutions (n_x = 80/160/320, Δt = 10⁻³ s),
which coincide with the literature solution within the plot resolution.
(b) Pointwise differences u_h − u_exact: the discrepancy is localized at
the wave crest and decreases with refinement — maximum error
1.16·10⁻⁴ (coarsest) → 7.19·10⁻⁶ (finest); measured local spatial orders
1.995 and 1.976 (Fig. S3a).*

---

**Figure 3 — OOFEM solution and the stationary hyperbolic Allen–Cahn front at T = 0.2.**
*Exact stationary front u(x) = tanh(x/√2) of the hyperbolic Allen–Cahn
equation u_tt − Δu = u − u³ on [−8,8]. (a) Black dashed curve: the exact
front; blue curves: the three archived OOFEM solutions (N_x = 40/80/160,
Δt = 10⁻³ s, T = 0.2 s), which coincide with it within the plot resolution.
(b) Pointwise differences u_h − u_exact: maximum error 7.27·10⁻⁵ (coarsest)
→ 6.59·10⁻⁶ (finest); measured local spatial orders 2.006 and 1.887 (the
second ratio is affected by the truncation of the front and by output
precision; Fig. S3b).*

---

## Supplementary figures (convergence view of the same archived runs)

**Figure S1 — OOFEM on the Hochbruck–Leibold exact profile: temporal refinement.**
*Same adapted exact-profile test as Fig. 1, on a fixed polygonal unit disk
(179 nodes). (a) RMS error versus Δt: the error is pinned by the spatial
floor of the fixed mesh (dotted line). (b) Difference between two
successive Δt refinements: local order ≈ 2.02, consistent with a
second-order time discretization, the same order class as the schemes
analyzed in the cited works (second-order IMEX in Hochbruck & Leibold 2021;
Gauss–Runge–Kutta order 2s in Hipp & Kovács 2021, Theorems 7.1–7.2).
T = 0.2 s; Δt = 10⁻², 5·10⁻³, 2.5·10⁻³ s.*

**Figure S2 — OOFEM on the Hochbruck–Leibold exact profile: spatial refinement on the polygonal unit disk.**
*Same adapted test with fixed Δt = 2.5·10⁻³ s and gmsh-triangulated
polygonal unit disks of 60, 179 and 502 nodes (h ≈ 0.30, 0.16, 0.09).
(a) RMS error versus mesh size; the measured local log₂ ratios are ≈ 0.33
and ≈ 0.24 — no convergence order is claimed. (b) Maximum pointwise error,
essentially flat (≈ 4.9·10⁻²): a geometric error floor. The reduced spatial
rate is the boundary-approximation effect emphasized in the bulk–surface
literature for non-polygonal domains (Kovács & Lubich 2017, Introduction;
Hipp & Kovács 2021, Sec. 2), here dominated by the polygonal approximation
of the circular boundary and the tangential derivatives at the vertices.
The case is exploratory and is not evidence for a spatial convergence rate
of the coupled model.*

**Figure S3 — OOFEM code verification against closed-form solutions of the hyperbolic Allen–Cahn bulk dynamics.**
*(a) Exact traveling wave of Nizovtseva et al. (2016), Eq. (29)
(τ = 0.2): measured local spatial orders 1.995 and 1.976. (b) Stationary
hyperbolic Allen–Cahn front tanh(x/√2): measured local orders 2.006 and
1.887. (c) Linear cosine mode with natural Neumann conditions: RMS error
versus Δt against the exact solution of the semidiscrete Q1
(consistent-mass) system — a time-stepping test; local order 1.997.
(d) Spatially uniform Allen–Cahn oscillator u_tt = u − u³: RMS error versus
Δt against an independent ODE integration; local order 2.016 (the third
step is close to the output-precision floor). Annotated numbers are
measured local log₂ ratios; no asymptotic order is claimed.*
