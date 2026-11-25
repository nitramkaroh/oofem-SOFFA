# Bibliography — OOFEM Allen–Cahn validation

All references below were read from the PDFs on disk (paths given); page and
equation anchors refer to the printed numbering of each document. No
bibliographic datum in this file was taken from memory or invented.

Relative paths are taken from the campaign root
`doc/stvowu_validation/validation_allen_cahn/`; a `../` prefix refers to the
pre-existing corpus in `doc/stvowu_validation/`.

## Primary benchmark sources

**[O1] Makoto Okumura**, *Design and analysis of structure-preserving schemes
for parabolic partial differential equations with dynamic boundary
conditions*, Dottorato di Ricerca (PhD dissertation), Graduate School of
Information Science and Technology, Osaka University, January 2021.
DOI 10.18910/82279.
PDF on disk: `/home/marco/Documents/Work/Riccardo/refs/31954_Dissertation.pdf`.

Used in the parabolic benchmark (all anchors re-verified against the PDF on
2026-10-27; PDF page = print page + 10):
- Ch. 4 (§1, print p. 40; PDF p. 49–50), governing parabolic Allen–Cahn
  problem: u_t = u_xx − F'(u) in (0,L)×(0,∞), eq. (4.1); dynamic boundary
  conditions u_t(0,t) = u_x(0,t) − F'(u(0,t)) and u_t(L,t) = −u_x(L,t) −
  F'(u(L,t)), eqs. (4.2)–(4.3);
- energy dissipation, print p. 41 (PDF p. 50): J(u) = ∫₀ᴸ (|u_x|²/2 + F(u)) dx
  and d/dt {J(u(t)) + F(u(0,t)) + F(u(L,t))} ≤ 0, eq. (4.5);
- computational section of Ch. 4 (§6, print p. 70; PDF p. 79): same dynamic
  law multiplied by ε_ex, eq. (4.91); K = 100 (Δx = 10⁻²), L = 1,
  ε_ex = 10, F(s) = (γ/4)(s²−1)² with γ = 100;
- **Computation example 2** (§6.2, print p. 71; PDF p. 80): N = 700
  (Δt = 10⁻³), T = 0.7, u(x,0) = exp{−500(x − 0.5)²};
- Fig. 4.3 (field evolution) and **Fig. 4.4 (energy decay, axis scale 0–90)**,
  print p. 71 (PDF p. 80) — Fig. 4.4 digitized, see
  `parabolic/literature_okumura_fig44_digitized.csv`.

**[N1] I. G. Nizovtseva, P. K. Galenko, D. V. Alexandrov**, "The hyperbolic
Allen–Cahn equation: exact solutions", *J. Phys. A: Math. Theor.* **49**
(2016) 435201.
PDF on disk:
`/home/marco/Documents/Work/Riccardo/refs/Nizovtseva_2016_J._Phys._A__Math._Theor._49_435201_new.pdf`.

Used in the hyperbolic benchmark (anchors re-verified against the PDF on
2026-10-27; print = PDF page):
- §2, p. 4: hyperbolic Allen–Cahn equation (1), τ_R ∂²_t f + ∂_t f = D∇²f − M
  f′(f), with free-energy density f(f) = (1/4)(f²−1)², eq. (2), giving
  ∂²_t f + ∂_t f = D∇²f + M(f−f³), eq. (3); dimensionless traveling-wave ODE
  (1−τc²)f″ + cf′ + f − f³ = 0, eq. (4), τ = τ_R M;
- p. 5: velocity c = ±3/√(2+9τ), eq. (26); **Table 1** (sign counters Q, A, B,
  k = 1..8, solution k = 1: Q = A = B = +1);
- p. 6: **bounded traveling wave, eq. (29)** — with the sign choice of
  solution k = 1 (Q = A = B = +1), c₀ = 0 and the positive velocity:
  u(x,t) = 1/(1+exp[λ(x−ct)]), λ = √(2+9τ)/2, c = 3/√(2+9τ);
  boundedness class discussed in §5.2.1 (p. 10).

## Secondary sources (contextual / not directly benchmarked)

**[C1] L. Calatroni, P. Colli**, "Global solution to the Allen–Cahn equation
with singular potentials and dynamic boundary conditions",
arXiv:1206.6738v1 [math.AP], 28 Jun 2012.
PDF on disk: `/home/marco/Documents/Work/Riccardo/refs/1206.6738v1_new.pdf`.
Well-posedness theory for the same parabolic problem with dynamic boundary
condition and **singular** potentials; contains no numerical benchmark with an
exact or published reference solution. Used only to contextualize the
regular-quartic case of [O1].

**[X1] M. Xiao, R. Chen**, "A first-order energy stable scheme for the
Allen–Cahn equation with the Allen–Cahn type dynamic boundary condition",
*J. Comput. Appl. Math.* **460** (2025) 116409.
PDF on disk:
`/home/marco/Documents/Work/Riccardo/refs/1-s2.0-S0377042724006575-main.pdf`.
Energy-stable finite-difference scheme with numerical experiments; its test
cases (e.g. the "case 4" droplet relaxation) were explored in the pre-existing
corpus (`../literature_comparison/xiao_chen_case4_*`) but are **not** part of
the two main benchmarks of this validation, because the published data are
contours/energies of a 2D droplet, not a 1D profile comparable to a Q1 strip
solution.

**[K1] G. Kovács, C. Lubich**, "Numerical analysis of parabolic problems with
dynamic boundary conditions", *IMA J. Numer. Anal.* **37**(1), 1–39, 2017,
DOI 10.1093/imanum/drw015 (arXiv:1501.01882).
PDF on disk: `/home/marco/Documents/Work/Riccardo/refs/imanum_2Fdrw015.pdf`
(= `../literature_analysis_2026/pdfs/P1_kovacs_lubich_2017_imajna.pdf`).
A priori error analysis for parabolic problems with dynamic boundary
conditions (linear and semilinear); the numerical experiments are for the
linear heat equation and a semilinear test, not for the nonlinear
Allen–Cahn dynamic-boundary problem of [O1]. Contextualizes the boundary
regularity theory behind the OOFEM boundary element.

**[M1] A. Makki, A. Miranville, M. Petcu**, "A numerical analysis of the
coupled Cahn–Hilliard/Allen–Cahn system with dynamic boundary conditions",
*Int. J. Numer. Anal. Model.* **19**(5), 2022 (HAL hal-05015316).
PDF on disk: `/home/marco/Documents/Work/Riccardo/refs/MMP1-rewritten2.pdf`.
Coupled system; out of scope of the present single-equation benchmarks.

**[H1] N. Hipp, G. Kovács**, "Finite element error analysis of wave equations
with dynamic boundary conditions: L² estimates", *IMA J. Numer. Anal.*
**41**(1), 638–728, 2021, DOI 10.1093/imanum/drz073.
PDF on disk: `../literature_analysis_2026/pdfs/P2_hipp_kovacs_2021_imajna.pdf`.
Wave-equation (inertial) theory with dynamic boundary conditions; theory
paper for the linear problem — it motivates the hyperbolic model but contains
no Allen–Cahn traveling-wave benchmark.

**[L1] M. Hochbruck, M. Leibold**, "An implicit–explicit time discretization
scheme for second-order semilinear wave equations with application to dynamic
boundary conditions", *Numer. Math.* **147**(4), 869–899, 2021,
DOI 10.1007/s00211-021-01184-w.
PDF on disk: `../literature_analysis_2026/pdfs/P3_hochbruck_leibold_2021_nummath.pdf`.
**[L2] M. Hochbruck, M. Leibold**, "Finite element discretization of
semilinear acoustic wave equations with kinetic boundary conditions",
*ETNA* **53**, 522–540, 2020, DOI 10.1553/etna_vol53s522.
PDF on disk: `../literature_analysis_2026/pdfs/P4_hochbruck_leibold_2020_etna.pdf`.
IMEX schemes for semilinear wave equations with dynamic/kinetic boundary
conditions; their exact profile solutions are for the *linear* (or weakly
nonlinear) wave problem, not for the hyperbolic Allen–Cahn equation. Pre-
existing OOFEM replicas of the L1/L2 profiles and their honest limitations are
documented in `../hyperbolic/` (see the pre-existing corpus README).

**None of [K1], [H1], [L1], [L2] provides an exact or published reference
solution for the coupled inertial-bulk / dynamic-surface Allen–Cahn model
implemented in OOFEM.** This is the basis of the explicit declaration in
`validation/comparison_report.md`: in the hyperbolic traveling-wave test the
dynamic boundary condition was **not** validated.

## Software and methodology references (context only)

**[P1] T. Patzák**, "OOFEM — an object-oriented simulation framework for
finite element computations", *Acta Polytechnica* **52**(6), 59–66, 2012.
The OOFEM framework used here is the local 2.6.devel source build
(see master README for the exact build).

**[H2] E. Hairer, S. N. Nørsett, G. Wanner**, *Solving Ordinary Differential
Equations I: Nonstiff Problems*, 2nd ed., Springer, 1993.
Methodology behind the DOP853 independent references used in the pre-existing
coupled-model audit (`../hyperbolic/coupled_independent_jacobian_fix/`).

**[R1] C. T. Roache**, "Code verification by the method of manufactured
solutions", *J. Fluids Eng.* **124**(1), 4–10, 2002.
Methodology reference for the exact-solution verification used here.
