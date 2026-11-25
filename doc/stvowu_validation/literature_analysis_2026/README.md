# Reference literature analysis — OOFEM Allen–Cahn validation (2026)

Date: 2026-10-27. Directory: `doc/stvowu_validation/literature_analysis_2026/`.

This dossier analyzes the seven collected publications (P1–P7) against the OOFEM
coupled hyperbolic Allen–Cahn model with dynamic boundary, and determines — for
each paper — what it actually contains, what can and cannot be compared with the
OOFEM results, and where a genuinely comparable reference problem would have to
come from.

## Files

| File | Content |
|---|---|
| `README.md` | this file: scope, method, rule compliance, PDF provenance |
| `01_oofem_model_statement.md` | the exact OOFEM model as assembled by the code (Rule 3): discrete system, strong form, surface law, provenance in source files |
| `02_paper_by_paper.md` | P1–P7: what each paper actually contains (theory vs. numerical experiments), with page / equation / table / figure anchors (Rule 4), and per-paper verdict |
| `03_compatibility_matrix.md` | compatibility matrix (PDE class, boundary-condition type, domain, norms, time scheme) and the quantitative-comparability verdicts (Rules 2, 5, 6) |
| `04_existing_oofem_replica_attempts.md` | the pre-existing OOFEM implementations of the Hochbruck–Leibold exact profiles (`hl_disk`, `hl_profile`), their metrics and honest limitations (Rule 7) |
| `05_limitations_and_search_directions.md` | the Rule-7 gap report: what none of P1–P7 provides, what a truly comparable benchmark must have, and concrete search directions |
| `06_scope_separation_dop853.md` | separation of the OOFEM-versus-literature comparison from the OOFEM-versus-DOP853 comparison (Rule 8); DOP853 provenance (SciPy 1.11.4); MMS statement (P7) |

## Corpus and PDF provenance

| ID | Paper | PDF on disk | PDF version read |
|---|---|---|---|
| P1 | Kovács & Lubich, *Numerical analysis of parabolic problems with dynamic boundary conditions*, IMA J. Numer. Anal. 37(1), 1–39, 2017, DOI 10.1093/imanum/drw015 | `pdfs/P1_kovacs_lubich_2017_imajna.pdf` | arXiv:1501.01882v1 ("Version of 9 April 2018"), 30 pp, **read in full** |
| P2 | Hipp & Kovács, *Finite element error analysis of wave equations with dynamic boundary conditions: L² estimates*, IMA J. Numer. Anal. 41(1), 638–728, 2021, DOI 10.1093/imanum/drz073 | `pdfs/P2_hipp_kovacs_2021_imajna.pdf` | arXiv:1901.01792v2, 44 pp, read pp. 1–14 and 23–44 (all model, theorem and experiment sections; only interior proof pages 15–22 not read page-by-page) |
| P3 | Hochbruck & Leibold, *An implicit-explicit time discretization scheme for second-order semilinear wave equations with application to dynamic boundary conditions*, Numer. Math. 147(4), 869–899, 2021, DOI 10.1007/s00211-021-01184-w | `pdfs/P3_hochbruck_leibold_2021_nummath.pdf` | open-access journal version, 31 pp, **read in full** |
| P4 | Hochbruck & Leibold, *Finite element discretization of semilinear acoustic wave equations with kinetic boundary conditions*, ETNA 53, 522–540, 2020, DOI 10.1553/etna_vol53s522 | `pdfs/P4_hochbruck_leibold_2020_etna.pdf` | open-access journal version, 19 pp, **read in full** |
| P5 | Patzák, *OOFEM — an object-oriented simulation framework for finite element computations*, Acta Polytechnica 52(6), 59–66, 2012, DOI 10.14311/1678 (as provided) | not available | not readable: download from ojs.cvut.cz timed out in this environment (recorded in the session task log); **citation only**, per the task instruction this software-description paper is not treated as a numerical benchmark |
| P6 | Hairer, Nørsett, Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed., Springer, 1993, DOI 10.1007/978-3-540-78862-1 | not on disk | book; cited for the adaptive Runge–Kutta methodology only (see `06_scope_separation_dop853.md`); not a source of published benchmark results |
| P7 | Roache, *Code verification by the method of manufactured solutions*, J. Fluids Eng. 124(1), 4–10, 2002, DOI 10.1115/1.1436090 | not on disk | background on the MMS methodology; the exact statement of where (and where not) MMS was used in the present dossier is in `06_scope_separation_dop853.md` |

P1–P4 are the CRC 1173 preprints referenced in `hyperbolic/README.md` in their
published/open-access form.

## Method (Rule 1)

The PDFs were read page by page and the **actual numerical experiments** were
extracted (governing equations, boundary conditions, initial data, exact or
reference solutions, meshes, time steps, norms, figures), not the abstracts.
Where a statement in this dossier quotes a paper, it gives the equation / table /
figure number and, where known from the reading, the page.

## Compliance with the eight rules

| Rule | Where satisfied |
|---|---|
| 1. Read the PDFs, inspect the actual numerical experiments | §"Method" above; per-paper evidence in `02_paper_by_paper.md` |
| 2. Preserve the distinction parabolic / hyperbolic / acoustic, and of the different boundary-condition types | `03_compatibility_matrix.md` (rows 1–3); per-paper verdicts |
| 3. Identify the precise equations relevant to the OOFEM model | `01_oofem_model_statement.md` (from the OOFEM source and the independent reference implementation) |
| 4. Source page / equation / table / figure for each candidate comparison | `02_paper_by_paper.md`, per paper |
| 5. Quantitative comparisons only when problem and quantities are compatible | `03_compatibility_matrix.md` §3 — result: no paper passes this test for the full coupled model |
| 6. Theoretical-only papers → contextualize, do not fabricate comparison curves | `03_compatibility_matrix.md` §4 (contextualization-only usage per paper) |
| 7. If no sufficiently comparable benchmark, report the limitation and search for a more appropriate reference problem | `05_limitations_and_search_directions.md`; the pre-existing OOFEM replica attempts are documented in `04_existing_oofem_replica_attempts.md` |
| 8. Keep OOFEM-vs-literature separate from OOFEM-vs-DOP853 | `06_scope_separation_dop853.md` |

## Standing constraints honored

No OOFEM code or original script was modified; no new OOFEM simulation or new
DOP853 integration was run; all statements about OOFEM results in this dossier
reference pre-existing archived data (paths given). Where a paper's printed
content is internally inconsistent (one instance found in P4 §5), it is flagged
explicitly rather than silently corrected.
