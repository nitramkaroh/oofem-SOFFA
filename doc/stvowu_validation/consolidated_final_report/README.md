# Rapporto consolidato di validazione OOFEM–Allen–Cahn (post T1-B)

Data: 2026-10-27 — Versione 1.0. Directory dedicata al solo consolidamento
documentale: nessun risultato numerico è stato eseguito o modificato per produrre
questo rapporto; tutti i valori sono riferimenti a artefatti già archiviati.

## File di questa directory

| File | Contenuto |
|---|---|
| `validation_report_consolidated.md` | Report principale (sezioni A–J: sintesi esecutiva, modello e implementazione, matrice di validazione, correttezza/riproducibilità con storia degli errori, convergenza spaziale del caso accoppiato, angoli e regolarità, bilancio energetico, inventario artefatti, limiti residui, conclusioni e raccomandazioni). |
| `validation_matrix.csv` | Matrice di validazione: 18 righe (C1–C18) × 11 colonne, una riga per test, con riferimento, metrica, risultato, criterio, verdetto, limitazioni e artefatti fonte. |
| `coupled_spatial_convergence.csv` | Dati della convergenza spaziale del caso accoppiato (sezione E del report): errori vs proxy, differenze di famiglia, ordini, stabilità del proxy, massimo mobile, zone angolari, serie temporale degli ordini, floor del solutore, energia T1-B. Formato lungo: un valore per riga; le serie complete a 80 step restano nei CSV fonte (non duplicati). |
| `README.md` | Questo file. |

## Come tracciare ogni numero alle fonti

Tutti i percorsi sono relativi alla radice del dossier
`doc/stvowu_validation/`. Ogni valore citato nel report ha una colonna/riga
`artefatti_fonte` (nel CSV) o una citazione inline (nel report).

### Caso accoppiato iperbolico (T1, T1-A, T1-D, T1-B)

- **Risultati T1-B (proxy DOP853 n128×64, 21.0 min, RSS ≈ 109 MB, nfev 2117):**
  `hyperbolic/t1b_dop853_n128x64/results/` — `t1b_metrics_T.csv`,
  `t1b_family_orders_T.csv`, `t1b_proxy_stability_T.csv`,
  `t1b_moving_max_T.csv`, `t1b_zones_T.csv`,
  `t1b_orders_time_64_128.csv`, `t1b_orders_time_oofem_proxy128.csv`,
  `t1b_energy_n128.csv`; figure `t1b_error_vs_h_T.png`,
  `t1b_family_order_64_128_time.png`, `t1b_moving_max_rc.png`.
  Riferimento (10.9 MB): `hyperbolic/t1b_dop853_n128x64/reference_n128x64_dt0p00025.npz`
  (U, V: 81×8385); matrici `matrices_n128x64.npz` (sha256 68c85bfd…1330);
  log `run.log`, `validate.log`, `resource_report.txt`;
  report `hyperbolic/t1b_dop853_n128x64/t1b_report.md`;
  script `t1b_dop853.py`, `t1b_validate.py`, `t1b_plots.py`, `run_full.sh`.
- **T1-A (rerun OOFEM n32×16, steps 1–39 bit-identici all'originale):**
  `hyperbolic/t1a_n32x16_rerun/results/t1a_report.md` (+ CSV di confronto).
- **T1-D (audit analitico pre-T1-B):** `hyperbolic/t1d_analytic_audit/t1d_report.md`.
- **T1 originale + audit statico + erratum sul fit h²+h:**
  `hyperbolic/coupled_independent_spatial_jacobian_fix/static_audit_report.md`
  (§1.5 floor: 2.0995e-8 / 2.4441e-8 / 2.5972e-8; §4 erratum),
  `hyperbolic/coupled_independent_spatial_jacobian_fix/jacobian_fix_verification_report.md`
  (verifica post-fix: ordini 1.972/1.948/1.975/2.027).

### Caso parabolico e letterario

- **Caso parabolico Calatroni–Colli (C6):**
  `literature_comparison/calatroni_colli_temporal_refinement.csv`,
  `literature_comparison/spatial_convergence.csv`,
  `literature_comparison/plot_calatroni_colli_refinement.py`;
  audit delle metriche: `validation_metrics_audit.md` (row Calatroni:
  NOT VERIFIABLE sui raw, solo metriche aggregate conservate).
- **Parabolico 2D (C7, C8, energia secante):**
  `parabolic_2d/parabolic_2d_validation.md`,
  `parabolic_2d/spatial_refinement.csv`, `parabolic_2d/spatial_refinement_T0p0625.csv`,
  `parabolic_2d/temporal_refinement.csv`, `parabolic_2d/energy_balance_secant.csv`.
- **Okumura (C4) e Xiao–Chen (C5):** `literature_comparison/README.md`,
  `literature_comparison/okumura_ch4_example2_literature_comparison.txt`,
  `literature_comparison/xiao_chen_case4_*.csv`.
- **Nizovtseva e kink (C1, C2, C13):**
  `hyperbolic/hyperbolic_bulk_validation.md`,
  `hyperbolic/nizovtseva_bulk_v2_spatial.csv`,
  `hyperbolic/nizovtseva_bulk_v2_temporal.csv`,
  `hyperbolic/nizovtseva_bulk_v2_energy.csv`,
  `hyperbolic/stationary_kink_space_refinement.csv`.
- **Newmark semidiscreto (C3):** `hyperbolic/bulk_reference_metrics.csv`,
  `hyperbolic/README.md`.
- **Lee 2024 (C14):** `validation_metrics_audit.md` (L2 6.0675e-3 vs 2.61689e-3,
  fattore 2.32× irrisolto).

### Energia, riproducibilità e audit

- **Bilanci energetici:** `energy_validation_final.md` (fasi 1–5, CSV legacy
  non autorevoli, limiti di secanti/ricostruzione), `energy_theory.md`,
  `energy_balance.csv`, `energy_summary.csv`.
- **Riproducibilità:** `reproducibility_runs.md` (3 rerun OOFEM bit-identici),
  `reproducibility_audit.md` (per-test), `hyperbolic/t1b_dop853_n128x64/`
  (assembly fresca == cache: 0.000e+00).
- **Audit iniziale con criteri pre-registrati e plan T1–T5:**
  `audit_iniziale_2026-10-27.md` (§7–8 criteri [1.5, 2.5] / [0.75, 1.25]).
- **Storia degli errori (13 voci):** sezione D.5 di
  `validation_report_consolidated.md`; voci originali in
  `static_audit_report.md` §4 (erratum), `t1a_report.md`, `t1d_report.md`,
  `t1b_report.md`, `validation_metrics_audit.md`, `energy_validation_final.md`.

### Note di tracciabilità

- **Nessun file pesante è stato copiato** in questa directory: i `.npz`, i raw
  VTK/VTU e i log restano dove sono; qui ci sono solo riferimenti.
- I **raw VTK di diversi test sono in `/tmp`** (non versionati, non durevoli):
  dove il report lo dice, i dati originari sono *NOT VERIFIABLE* e il fatto è
  dichiarato esplicitamente invece di ricostruire i valori.
- **CSV energetici legacy** (Gauss 2×2, `energy_balance.csv` ecc.) sono
  marcati *non autorevoli* in `energy_validation_final.md`; i valori citati nel
  report consolidato provengono dai ricalcoli moderni (3×3/5×5).
- I file marcati **SOPPIANTATI** (es. `validation_matrix_final.md` con la riga
  "Sistema iperbolico accoppiato FAIL") non sono usati per i verdetti correnti;
  sono conservati per la storia e segnalati in D.5.
- Build di riferimento: `rel-2.1-5430-g2cc1a7706-dirty` (debug). Python
  3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3 (nessun lockfile).
