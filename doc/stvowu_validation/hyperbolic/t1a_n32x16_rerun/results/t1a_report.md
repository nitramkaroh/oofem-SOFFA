# T1-A — Riesecuzione OOFEM n32×16 (recupero T=0.02): rapporto

**Data:** 2026-10-27
**Stato:** COMPLETATO entro tutti i limiti dichiarati. Nessuna seconda esecuzione.
**Directory nuove:** `t1a_n32x16_rerun/oofem/` (esecuzione) e `t1a_n32x16_rerun/results/`
(output di post-processing, log risorse, questo rapporto). Nessun archivio preesistente
modificato o sovrascritto.

---

## 1. Esecuzione — tempo e memoria effettivi

| voce | valore |
|---|---|
| binario | `build/debug_magneto/oofem` (stesso del run originale, `run_coupled_independent.py`) |
| comando | `oofem -f <nuova_dir>/case.in`, cwd = nuova directory (identico all'originale) |
| input | `case.in` **byte-identico** all'archivio originale (verificato con `cmp`; la rigenerazione da `make_input(32,16,2.5e-4,'space_n32x16_dt0p00025')` è anch'essa byte-identica) |
| inviluppo | 1 processo, 1 core (`taskset -c 0`), cap memoria virtuale 8 GB (`ulimit -v 8388608`), timeout hard 60 s (`timeout -s TERM --kill-after=5`), memoria monitorata con GNU `/usr/bin/time -v` |
| **real time (OOFEM)** | **10 s** (wall clock 10.25 s; budget 60 s) |
| **CPU** | 99% di un core (user 10.18 s, sys 0.02 s) |
| **picco RSS** | **32768 kB = 32 MB** |
| esito | exit 0, "ANALYSIS FINISHED", 0 errori, 0 warning |
| swap | 0 |

Confronto dei tempi reali OOFEM (dai rispettivi run.log): n8×4 **0 s** (<1 s),
n16×8 **2 s**, n32×16 (nuovo) **10 s** — scala come atteso (~×3.7 DOF).

## 2. Integrità del nuovo output

- 80/80 VTU: parsabili, 561 punti, coordinate identiche alla mesh (atol 1e-13),
  `Temperature` finito; pvd con 80 DataSets, tempi t = k·2.5×10⁻⁴ esatti,
  **ultimo t = 0.02 raggiunto**. (Dettaglio: `t1a_new_run_integrity.csv`.)
- **Conferma di impostazioni identiche**: per gli step 1–39 (intatti nell'archivio
  originale) il campo del nuovo run è **bit-identico** a quello originale:
  max\|u_new − u_old\| = **0.000×10⁰** su 39 step. Esecuzione deterministica a 1 core.

## 3. Metriche sincronizzate a tempo comune T=0.02 (proxy n64×32 DOP853, norme verificate)

Stesse norme e convenzioni della verifica principale (FE RMS con M_b/M_s corretti,
max interni/confine, restrizione per coincidenza nodale esatta, U[80] a t=0.02):

| mesh | bulk FE RMS | bulk max | surf FE RMS | surf max |
|---|---|---|---|---|
| n8×4 | 1.3286×10⁻³ | 8.345×10⁻⁴ | 2.5727×10⁻³ | 6.069×10⁻³ |
| n16×8 | 8.188×10⁻⁴ | 1.880×10⁻³ | 1.9948×10⁻³ | 4.744×10⁻³ |
| n32×16 (nuovo) | 4.141×10⁻⁴ | 3.123×10⁻³ | 9.265×10⁻⁴ | 2.230×10⁻³ |

**Floor del solvere a T** (nuovo n32 OOFEM vs DOP853 n32 in cache, stessa mesh e
stesso tempo; nessuna nuova integrazione): bulk FE RMS **2.9×10⁻⁸**, bulk max
4.4×10⁻⁷, surf FE RMS 9.3×10⁻⁸ → 6–7 ordini sotto l'errore spaziale.
(`t1a_solver_floor.csv`.)

## 4. Ordini spaziali osservati — entrambi gli intervalli a t=0.02

`order = log2(err_grosso/err_fine)`; criterio preliminare [1.5, 2.5].

| intervallo | bulk FE RMS | bulk max | surf FE RMS | surf max |
|---|---|---|---|---|
| 8×4 → 16×8 | **0.698** | **−1.172** (non monotono) | **0.367** | 0.355 |
| 16×8 → 32×16 | **0.983** | **−0.732** (non monotono) | **1.106** | 1.089 |

(`t1a_convergence_orders.csv`.)

Osservazioni:
- Gli ordini 8→16 sono invariati (stessi dati di prima); gli ordini 16→32 a T sono
  superiori a quelli misurati a t=0.00975 (0.761/0.718 → 0.983/1.106 FE RMS):
  l'ordine dipende dal tempo (l'errore si accumula, D3 della verifica principale).
- **bulk max resta non monotono** (8.3×10⁻⁴ → 1.88×10⁻³ → 3.12×10⁻³) su entrambi
  gli intervalli a T: segnale flaggato, localizzato agli angoli (CSV per-nodo).
- FE RMS: monotono decrescente su entrambi gli intervalli.

## 5. Fit e(h) = a·h² + b·h a tempo comune (3 livelli, h = 0.25/0.125/0.0625)

| metrica | a (h²) | b (h) | residuo max |
|---|---|---|---|
| bulk FE RMS | **−8.77×10⁻³** | +7.52×10⁻³ | 5.18% |
| surf FE RMS | **−3.72×10⁻²** | +1.97×10⁻² | 16.96% |

(`t1a_fit_3point.csv`.) **Limiti dichiarati:** tre punti con due parametri non
consentono alcuna verifica indipendente della bontà del modello (1 solo grado di
libertà residuo); il coefficiente di h² risulta **negativo** in entrambe le variabili,
cioè il modello h²+h a coefficienti positivi **non è supportato** nemmeno dai dati a
tempo unico. **Il fit non è presentato, e non deve essere usato, come prova di
singolarità angolare** (né a suo favore né a suo sfavore): è solo una descrizione
locale dei tre punti. L'ipotesi di singolarità angolare resta un'interpretazione
coerente con le evidenze osservate (concentrazione energetica angolare ~74%,
simmetria 180° esatta, decadimento di potenza verso l'angolo) ma non è dimostrata
dal dataset attuale.

## 6. Analisi di zona a T (solo dati esistenti; metriche nodali non pesate)

| zona (r<0.5) | n8: RMS / energia | n16: RMS / energia | n32 nuovo: RMS / energia | ordine n16→n32 (nodal) |
|---|---|---|---|---|
| A: (0,1)/(2,0) | 3.97×10⁻³ / 0.708 | 2.26×10⁻³ / 0.725 | 9.92×10⁻⁴ / **0.739** | 1.185 |
| B: (0,0)/(2,1) | 1.76×10⁻³ / 0.139 | 1.03×10⁻³ / 0.151 | 4.57×10⁻⁴ / 0.157 | 1.170 |
| resto | 9.71×10⁻⁴ / 0.153 | 5.30×10⁻⁴ / 0.124 | 2.14×10⁻⁴ / 0.104 | 1.308 |

L'energia di errore si concentra progressivamente nella zona A (70.8% → 72.5% →
73.9%) col raffinamento; gli ordini di zona a T (~1.2) non differiscono
sostanzialmente dal resto del dominio. (`t1a_region_analysis_T.csv`.)

## 7. Verdetto aggiornato

| voce | verdetto | motivazione |
|---|---|---|
| Convergenza spaziale **bulk** | **NON VALIDATA** (criterio preliminare) | Ordini a T: 0.698 / 0.983 ∉ [1.5, 2.5]; bulk max non monotono su entrambi gli intervalli. |
| Convergenza spaziale **superficie** | **NON VALIDATA** (criterio preliminare) | Ordini a T: 0.367 / 1.106 (FE RMS) e 0.355 / 1.089 (max) ∉ [1.5, 2.5]. |
| **Validazione complessiva dell'accoppiamento** | **SUPPORTATA PARZIALMENTE** | Il solvere OOFEM riproduce il riferimento indipendente a ~3×10⁻⁸ FE RMS a T (6–7 ordini sotto l'errore spaziale); il nuovo run è bit-identico agli step intatti dell'archivio (riproducibilità dell'implementazione confermata). L'ordine quadratico non è raggiunto: il comportamento osservato (FE RMS ~1.0–1.1 sull'intervallo fine a T, concentrazione angolare stabile ~74%) resta compatibile con una limitazione di regolarità della discretizzazione accoppiata, ma ciò resta **ipotesi**, non evidenza (fit h²+h non supportato, §5). |

## 8. Limiti interpretativi (dichiarati)

1. Due soli intervalli di raffinamento a un unico tempo (T=0.02): l'ordine dipende
   dal tempo (0.76→0.98 bulk tra t=0.00975 e 0.02); un singolo tempo non
   caratterizza la dipendenza temporale dell'ordine.
2. Il fit 3 punti / 2 parametri (§5) non ha validazione indipendente; coefficiente
   di h² negativo → modello h²+h a coefficienti positivi non supportato.
3. Le metriche di zona sono nodali (non pesate FE) con raggio fisso 0.5; gli ordini
   di zona sono quindi indicative.
4. Nessun run OOFEM su n64×32 (non autorizzato): il proxy resta DOP853-only;
   l'errore OOFEM a n64 è sconosciuto.
5. L'unico campo esportato è `u` (`Temperature`); nessun controllo su `v`.
6. La corruzione dell'archivio originale n32×16 (step 40–80) resta tale: l'archivio
   è rimasto intatto; il nuovo run è in directory separata.

## 9. File prodotti (tutti in directory nuove)

`oofem/`: case.in, run.log, resource_report.txt, case.out.m0.pvd, 80× VTU.
`results/`: `t1a_new_run_integrity.csv`, `t1a_metrics_common_time.csv`,
`t1a_convergence_orders.csv`, `t1a_solver_floor.csv`, `t1a_fit_3point.csv`,
`t1a_region_analysis_T.csv`, `per_node_n32x16_new_step80_T0p02_vs_proxy_n64x32.csv`,
questo rapporto.
