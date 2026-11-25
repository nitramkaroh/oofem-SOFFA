# Benchmark selection and literature compatibility matrix

Date: 2026-10-27. All statements below are anchored to the PDFs on disk;
see `bibliography.md` for full citations and paths.

## Compatibility classes

| Class | Meaning |
|---|---|
| **A** | Directly quantitatively comparable: same PDE class, same type of boundary condition, and a published exact/reference solution or published quantitative data for a specified problem. |
| **B** | Partially compatible: same or close PDE/BC class, but no reference solution for the *nonlinear Allen–Cahn* problem (e.g. linear experiments, different data, theory-only). Usable for context or for parts of the model. |
| **C** | Different problem class (different PDE or boundary-condition type); context only. |
| **D** | No benchmark content relevant here (software description, general methodology); context only. |

## Matrix

| Ref. | PDE | Boundary condition | Reference data published? | Class | Use in this validation |
|---|---|---|---|---|---|
| [O1] Okumura 2021, Ch. 4 | parabolic Allen–Cahn, u_t = u_xx − F'(u), regular quartic F | **dynamic** (1st order) at both ends | yes — full parameters (K=100, N=700, T=0.7, γ=100, ε_ex=10, IC), Fig. 4.3 field, **Fig. 4.4 energy** | **A (parabolic)** | **Benchmark 1**: reproduced in OOFEM; energy and field compared |
| [N1] Nizovtseva et al. 2016 | hyperbolic Allen–Cahn, τ_R u_tt + u_t = DΔu + M(u−u³) | none (Cauchy problem; traveling wave) | yes — **exact bounded traveling wave**, eq. (29) + Table 1 | **A (hyperbolic, bulk)** | **Benchmark 2**: reproduced in OOFEM (bulk only) |
| [C1] Calatroni–Colli 2012 | parabolic Allen–Cahn, **singular** potentials | dynamic | no (existence/uniqueness only) | B | context: well-posedness of the same dynamic-boundary structure |
| [X1] Xiao–Chen 2025 | parabolic Allen–Cahn | dynamic (Allen–Cahn type) | yes, but 2D droplet contour/energy data, FD scheme | B | exploratory only in the pre-existing corpus; not a 1D-comparable benchmark |
| [K1] Kovács–Lubich 2017 | parabolic, dynamic BC (linear + semilinear theory) | dynamic | experiments for linear/semilinear heat, not Allen–Cahn | B | context: boundary regularity and error theory behind the OOFEM boundary element |
| [M1] Makki–Miranville–Petcu 2022 | coupled Cahn–Hilliard/Allen–Cahn | dynamic | yes (their own scheme, coupled system) | C | out of scope (single equation here) |
| [H1] Hipp–Kovács 2021 | wave (linear) with dynamic BC | dynamic | L² error theory only | C | context: inertial + dynamic BC theory |
| [L1] Hochbruck–Leibold 2021 | semilinear wave, IMEX | dynamic | exact *linear* wave profiles | C | pre-existing OOFEM replicas documented in `../hyperbolic/` |
| [L2] Hochbruck–Leibold 2020 | semilinear acoustic wave | kinetic/dynamic | exact *linear* profiles | C | as [L1] |
| [P1] Patzák 2012 | — | — | — | D | software description |
| [H2] Hairer et al. 1993 | — | — | — | D | DOP853 methodology (pre-existing coupled-model audit) |
| [R1] Roache 2002 | — | — | — | D | exact-solution verification methodology |

## Rationale

**Benchmark 1 (parabolic, dynamic BC).** [O1] is the only Class-A source for
the parabolic Allen–Cahn problem *with a dynamic boundary condition*: the
dissertation publishes the complete numerical recipe (mesh, time step,
horizon, γ, ε_ex, initial condition) **and** a quantitative reference curve
(Fig. 4.4, energy E(t), axis scale 0–90) plus the field evolution (Fig. 4.3).
This is the strongest available anchor for the OOFEM parabolic model with the
first-order dynamic law ε_ex ψ_t = −∂_n u − F'(ψ). [C1] and [K1] confirm the
well-posedness/error theory of the same structure but publish nothing
comparable quantitatively. [X1] publishes 2D droplet data for a different
discretization; a prior exploratory comparison exists in the pre-existing
corpus but is not a clean 1D benchmark, so it is not used for the headline
verification.

**Benchmark 2 (hyperbolic, traveling wave).** [N1] is the only Class-A source
for the hyperbolic Allen–Cahn equation: it provides an **exact bounded
traveling wave** (eq. 29, sign set of Table 1) for the dimensional PDE with
τ = 0.2, D = M = 1. The OOFEM hyperbolic model contains the inertial bulk term
ρ u_tt + ρδ u_t = k Δu + m(u − u³); with ρ = τ_R M, ρδ = M, k = D, m = M it
is exactly the [N1] PDE. The wave carries **no boundary condition**, so the
test is run bulk-only (natural Neumann on the truncated strip, whose
truncation effect is quantified in the report).

**Explicit declaration (integrity requirement).** Among all documents
available, **none** provides an exact or published reference solution for the
hyperbolic Allen–Cahn problem *with a dynamic boundary condition* (the
dynamic/kinetic-boundary wave papers [H1], [L1], [L2] are theory or linear
profile results). Therefore, in Benchmark 2 the dynamic boundary condition is
**not exercised and not validated**, and this is stated in the figure, in the
executive summary and in the conclusions. The OOFEM coupled
inertial-bulk/dynamic-surface model was separately verified against a
corrected independent high-order (DOP853) semidiscrete reference in the
pre-existing audit `../hyperbolic/coupled_independent_jacobian_fix/`
(orders 1.95–2.03 vs that reference, residuals 5e-8/1.4e-7); that result is
reported as complementary evidence of the coupled *code path*, **not** as a
validation against literature.

## What was deliberately not done

- No benchmark was constructed by digitizing a curve and treating it as an
  "exact solution"; the only digitized data used are the [O1] Fig. 4.4 energy
  curve, clearly labeled as digitized everywhere it appears, and compared with
  the stated RMSE/MAE/max metrics.
- No parameter set was tuned to make the curves match: the OOFEM inputs use
  the published values verbatim (parabolic: N=100 cells, Δt=10⁻³, T=0.7,
  γ=100, ε_ex=10; hyperbolic: τ=0.2, D=M=1, Δt=10⁻³, T=0.2).
- The hyperbolic domain truncation [−12, 12] is a numerical necessity (the
  Cauchy problem has no boundary); it is disclosed and quantified, not hidden.
