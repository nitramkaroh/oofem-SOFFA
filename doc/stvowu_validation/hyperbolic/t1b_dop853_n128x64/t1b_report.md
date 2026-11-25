# T1-B — Proxy indipendente DOP853 su n128×64: rapporto finale

**Data:** 2026-10-27 · **Directory:** `hyperbolic/t1b_dop853_n128x64/`
**Sistema:** Allen–Cahn iperbolico bulk–surface accoppiato, Ω=[0,2]×[0,1], T=0.02, DT=2.5e-4 (80 passi),
u₀=0.2+0.1·cos(πx/2)cos(πy), v₀=0; RHO=KB=MB=CS=KS=MS=1, ETA=0.
**Scopo:** fornire un quarto livello spaziale indipendente da OOFEM per discriminare
convergenza spaziale vs limitazione di regolarità localizzata (verdetto T1-D: *NON VALIDATA*).

---

## 1. Script e parametri effettivamente utilizzati

| Componente | File | Parametri |
|---|---|---|
| Assemblaggio (corretto) | `../coupled_independent_jacobian_fix/check_jacobian_fix.py` (read-only) | `assemble_fixed`, `residual_reaction_fixed`, `energy_fixed`; quadratura Gauss 3 punti; Jacobiano di linea corretto (dN=[−.5,.5]·(2/L), w=0.5·L·w_a) |
| Costi FASE 0 | `t1b_fase0_cost.py` | misura statica (nessuna integrazione) |
| Driver | `t1b_dop853.py` | DOP853 (scipy `solve_ivp`), **rtol=2e-12, atol=2e-14** (identici al proxy n64 accettato), t_eval=81 punti canonici, `splu(Mb)` riusato, C=ETA·Mb+Ms=Ms, K=Kb+Ks |
| Validazione FASE 3 | `t1b_validate.py` | post-processing read-only; norme FE RMS (Mb/area, Ms/perimetro), max su nodi interni/bordo; restrizione per coincidenza nodale esatta (d²<1e-26); confronti a istante identico (step s ↔ U[s]) |
| Grafici | `t1b_plots.py` | rigenera i 3 PNG dai CSV (deterministico) |
| Wrapper bounded | `pilot_run.sh`, `run_full.sh`, `run_validate.sh` | `taskset -c 0`, `ulimit -v 2457600` (2.4 GB), `OMP/OPENBLAS/MKL_NUM_THREADS=1`, `timeout -s TERM --kill-after=5`, `/usr/bin/time -v` |

Indipendenza dal codice OOFEM: confermata (assemblaggio puro Python, nessun binario/codice/
output OOFEM nel percorso di integrazione; OOFEM entra solo come *dato letto in sola lettura*
nella FASE 3). Il Jacobiano di linea errato (versione originale) non è presente nel percorso
del riferimento: esiste solo in `orig.*`, usato esclusivamente per i confronti "prima".
Mesh n128×64: **N=8385 nodi** (129×65), 8192 quad, 384 spigoli di bordo, h=0.015625,
tutti gli elementi quadrati.

## 2. FASE 0 — Audit statico (esito: PASS)

- Stima costo RHS: **591 ms/valutazione** (misura su stato congelato, 30 ripetizioni).
- Stima totale: nfev 1088–1600 → **11.2–16.2 min**; memoria stimata ~100–300 MB
  (baseline VmSize 125.7 MB / RSS 61.4 MB; dopo assemblaggio+RHS 315.9/88.1 MB).
- Sparse: Mb/Kb nnz=74305 (0.11% di dense), Ms/Ks nnz=1152, tutte CSC, simmetria
  max|A−Aᵀ|=0. Nessuna conversione in dense; nessuna traiettoria ad alta risoluzione
  (solo 81 punti t_eval → 10.9 MB).
- Timeout e monitoraggio RAM: applicabili realmente (dimostrati nei run successivi).
- Verifiche quadratura: invariati su n128 (vedi §5.2) + 34 check sullo stesso percorso
  di assemblaggio (T1, fix dir): matrici di linea P1 esatte (1e-14), partizione
  dell'unità per la reazione, Mb·1=area, Ms·1=perimetro, somme di riga Kb/Ks=0.

## 3. FASE 1 — Pilot limitato (60 s, 1 processo/1 core)

| Misurata | Valore |
|---|---|
| Tempo totale | **33.2 s / 60 s** |
| nfev per il primo passo di output | 56 (esplorazione del passo iniziale DOP853) |
| Costo RHS | 590 ms (conferma la stima FASE 0) |
| RSS picco | 80.2 MB |
| n_accepted / first_step_h | 2 / 2.5e-4 |
| Stima run completo | ~11.6 min (cost model, **non** prova di budget) |

Il pilot non avrebbe contenuto assemblaggio (~28 s) + primo passo (~33 s) entro 60 s:
è stata usata una **cache delle matrici** (assemblata una volta, `matrices_n128x64.npz`,
sha256=68c85bfd…1330) caricata dal pilot; la `full` **riassembla fresco e verifica
uguaglianza esatta** con la cache (max|fresh−cache| = 0.000e+00 su Mb, Kb, Ms, Ks)
prima di integrare — nessuna scorciatoia non verificata entra nel riferimento.

## 4. FASE 2 — Esecuzione completa (1 processo/1 core)

**Esito: SUCCESSO, nessun arresto, nessun retry.**

| Voce | Limite | Valore effettivo |
|---|---|---|
| Tempo reale | 45 min (timeout 2700 s) | **21.0 min** (1259.3 s; soft wall 40 min non raggiunto) |
| Memoria | 2 GB aggiuntivi (ulimit 2.4 GB virt) | **RSS picco 109 MB** (`/usr/bin/time`; 98.8 MB in-process) |
| Parallelismo | 1 core | `taskset -c 0`, 99% CPU, user 1282.3 s |
| nfev | — | **2117** (26.46/passi) |
| Costo RHS | — | 595 ms (stima: 591) |
| final_t | 0.02 | **0.02 (esatto)** |
| Guard | — | nessuna violazione (NaN/Inf ad ogni RHS; tripwire \|u\|_max>5; crescita anomala) |
| Output | — | `reference_n128x64_dt0p00025.npz` (10.9 MB: U,V 81×8385, times, nfev, wall, rss_peak, final_t, cache_sha256) |

**Scostamento dalla stima (dichiarato):** il modello di costo usava il tasso stazionario
del run n64 (13.6 eval/step); a n128 DOP853 richiede **26.5 eval/step** (nfev 2117 vs ~1130
stimate) → 21.0 min vs 11.6 min stimati. Il costo per RHS era previsto correttamente.
Il run rientrava comunque nel budget (45 min / 2 GB) con ampio margine.
Nota di esecuzione: il run è girato come processo di fondo sotto il wrapper a limiti duri
(2700 s + ulimit + 1 core) perché la stima 11–16 min supera il limite dei 10 min dei
comandi in primo piano; l'autorizzazione FASE 2 prevedeva esplicitamente fino a 45 min.

## 5. FASE 3 — Validazione del proxy

### 5.1 Finitezza, completezza, tempo finale
U=(81, 8385), tutti finiti; `final_t` = 0.02 (scarto < 1e-12); `t_eval` identico alla griglia
canonica (81 punti a step DT). **OK.**

### 5.2 Quadratura, Jacobiani, termini di bordo, accoppiamento (su n128)
- Mb·1 = 2.0 (area) **esatto**; Ms·1 = 6.0 (perimetro) **esatto** (Ms corretto);
- somme di riga Kb = 1.8e-15, Ks = 0.0; simmetria esatta (FASE 0);
- identità esatte P1 per un elemento di linea + partizione dell'unità (T1, stesso percorso);
- accoppiamento: C = Ms (ETA=0), K = Kb+Ks, inversa su Mb — struttura identica al proxy n64
  accettato; coerenza cross-livello: la traiettoria n128 riproduce n64 con differenza
  FE RMS 4.87e-5 (il livello più fine della famiglia).
**OK** — il proxy è internamente consistente e costruito con la stessa quadratura verificata.

### 5.3 Errori a T=0.02 (convenzioni T1; dettaglio in `results/t1b_metrics_T.csv`)

| Quantità (bulk) | n8×4 | n16×8 | n32×16 (new) | fam. 32→64 | fam. 64→128 |
|---|---|---|---|---|---|
| FE RMS (vs proxy128) | 1.299e-3 | 7.903e-4 | 3.822e-4 | 4.141e-4 | **4.874e-5** |
| max interno (vs proxy128) | 8.35e-4 | 1.910e-3 | 2.957e-3 | 3.123e-3 | 7.485e-4 |

Surface FE RMS (vs proxy128): 2.509e-3 / 1.900e-3 / 8.372e-4; fam. 32→64: 9.265e-4; fam. 64→128: 1.642e-4.

**Ordini a T** (`results/t1b_family_orders_T.csv`):
bulk FE RMS: 8→16 **0.717**, 16→32 **1.048**, fam. 32→64 **3.087**;
bulk max: −1.195, −0.631, fam. 2.061; surf FE RMS: 0.401, 1.183, fam. 2.496; surf max: 0.462, 1.218, fam. 2.208.

### 5.4 Osservazioni chiave (dati, non interpretazioni)
1. **La famiglia continua a convergere fino a h=0.03125**: la differenza 64→128
   (4.87e-5 FE RMS) è 8.5× più piccola della 32→64 (4.14e-4). Nessun plateau osservato.
2. **L'errore OOFEM n32-new (3.82e-4) coincide con la differenza di famiglia 32→64 (4.14e-4)**
   entro l'8%: a n32 l'errore misurato è dominato dall'errore FE di famiglia; il residuo
   specifico-OOFEM resta al livello del floor 1e-8 di T1.
3. **Gli ordini osservati crescono col raffinamento** (0.72 → 1.05 → ~3.1) invece di
   stabilizzarsi in una finestra costante; la finestra accettazione [1.5, 2.5] non è
   coperta in modo costante su h ∈ [0.031, 0.25].
4. **Concentrazione angolare persistente al livello più fine**: 65.8% dell'energia della
   differenza 64→128 in zona A ((0,1)/(2,0)), 28.0% in B, 6.1% resto; OOFEM n32-new vs proxy128:
   74.6% / 14.2% / 11.2% (`results/t1b_zones_T.csv`).
5. **Massimo mobile (primo nodo interno verso l'angolo, T)**: 8.35e-4 → 1.910e-3 → 2.957e-3 →
   7.485e-4 per r_c = 0.354 → 0.177 → 0.088 → 0.044 (`results/t1b_moving_max_T.csv`).
   Non monotona (cresce poi cade); pendenza log-log su 4 punti **p = −0.016 ≈ 0**;
   segno negativo a tutti i livelli. Non è una legge di potenza: **nessuna singolarità
   può essere dichiarata dimostrata** su questa base.
6. **Serie temporale dell'ordine di famiglia 64→128 (dato nuovo, tutti i 80 passi**,
   `results/t1b_orders_time_64_128.csv`): l'ordine cresce monotono da **−0.45 (s1)** a
   **3.09 (T)** (bulk FE RMS; surface: −0.89 → 2.50, con picco 2.97 a s75). La differenza
   e64_128 **non è monotona nel tempo**: cresce da 3.4e-7 (s1) a un massimo ~1.55e-4
   (s45, t≈0.0113) e poi **decresce** a 4.87e-5 (T); e32_64 cresce invece monotono
   (2.5e-7 → 4.14e-4).
7. **Stabilità al livello di riferimento** (|err_vs64 − err_vs128|,
   `results/t1b_proxy_stability_T.csv`): FE RMS scambia di ≤ 3.2e-5 (2.2–8.3% in relativo);
   max bulk ≤ 1.7e-4; surf max n32: 4.0e-4 (~22%). Gli ordini si spostano di ≤ 0.07.
   I risultati non dipendono significativamente dalla scelta del proxy di riferimento.
8. **Energia della traiettoria n128** (corretta, `results/t1b_energy_n128.csv`): totale
   monotono decrescente 1.932710 → 1.927464 (−0.27% su T); superficie dominante (1.434),
   potenziale bulk 0.491, cinetica 2.4e-3. Nessuna crescita anomala.

### 5.5 Interpretazioni (esplicite, NON dimostrate)
- (I1) La sequenza è in **regime pre-asintotico dominato dal layer angolare**: coerente con
  ordini crescenti col raffinamento, concentrazione >65% in zona A a ogni livello,
  non monotonia del massimo puntuale. Non è una dimostrazione.
- (I2) La non monotonia del massimo puntuale è la **firma classica di limitazione di
  regolarità localizzata** ai vertici; resta compatibile anche con artefatti di campionamento
  del layer in posizione relativa fissa (h√2). **Singolarità: non dimostrata.**
- (I3) Il decadimento di e64_128 nel tempo finale (s45→T) è compatibile con lisciamiento/
  dissipazione dei tratti ad alta curvatura che distinguono i livelli (energia decrescente);
  non provato.
- (I4) Il valore ~3.1 dell'ordine sul pair più fine è un artefatto dell'estimatore in regime
  pre-asintotico (e64_128 in decadenza tarda), **non** evidenza di superconvergenza.

### 5.6 Limiti dell'analisi
- Non esistono riferimenti DOP853 n8/n16 a dt=2.5e-4: i primi 3 punti della catena del
  massimo mobile usano i valori OOFEM (≡ DOP853 a floor ~1e-8 a parità di mesh, T1);
  il 4° punto usa la differenza di famiglia 64→128 (mix dichiarato).
- Non esiste un livello oltre n128: l'errore assoluto di n128 non è misurabile;
  e64_128 è un surrogato (bound) per l'errore di n64 → l'ordine più fine (3.09) è un
  estimatore, non una misura.
- Gli ordini "OOFEM" misurano OOFEM vs proxy128 (incl. la catena di famiglia), non
  solo l'errore OOFEM-puro; la decomposizione è possibile solo entro l'8% (§5.4-2).
- Un solo istante per gli ordini OOFEM a 3 livelli (T); le serie temporali coprono la
  famiglia (T1-D aveva già campionato gli istanti intermedi per i livelli OOFEM).

## 6. Verdetti (separati, vocabolario del progetto)

| Domanda | Verdetto | Base |
|---|---|---|
| Riproducibilità (OOFEM) | **SUPPORTATA** (invariato, rafforzato) | T1 + §5.4-2: n32-new vs proxy128 = errore di famiglia, residuo al floor 1e-8 |
| Accordo col riferimento semidiscreto | **SUPPORTATA** (rafforzato) | Floor 1e-8 (T1) + coerenza cross-livello n64↔n128 (4.87e-5 FE RMS) |
| Convergenza spaziale — esistenza | **SOSTENUTA** | Errore decrescente a ogni livello fino a h=0.031; nessun plateau; 64→128 = minimo |
| Convergenza spaziale — tasso (criterio [1.5, 2.5]) | **NON ASSESTATO nel caso tested** | Ordini bulk FE RMS a T: 0.72, 1.05, ~3.1 (crescenti, non costanti); max non monotono |
| Singolarità ai vertici | **NON DEMOSTRATA** | Concentrazione angolare + non monotonia coerenti ma pende su 4 punti p≈0, non legge di potenza |

**Sintesi:** rispetto al verdetto T1-D (*convergenza NON VALIDATA*), T1-B dimostra che la
sequenza **converge** (nessuna saturazione fino a h=0.031) e che l'errore OOFEM a n32 è
dominato dall'errore FE di famiglia; ciò che **non** emerge è un regime asintotico di
secondo ordine costante sull'intervallo h ∈ [0.031, 0.25]: il comportamento è quello
pre-asintotico di un problema con layer angolari. Il verdetto complessivo sulla convergenza
spaziale è **SUPPORTATO PARZIALMENTE**.

## 7. Stima iniziale vs consumo effettivo (riassunto)

| Voce | Stima (FASE 0/1) | Effettivo |
|---|---|---|
| Tempo | 11.2–16.2 min (pilot: 11.6) | **21.0 min** |
| CPU | 1 core | user 1282.3 s, 99% CPU |
| RAM | ~100–300 MB | **109 MB** (RSS picco) |
| RHS | 591 ms/eval, ~1130–1600 eval | 595 ms/eval, **2117 eval** |
| Disco | ~11 MB | 10.9 MB npz + 2.8 MB cache |

## 8. Note di esecuzione
- Due esecuzioni fallite del **post-processing FASE 3** (bug dello script di analisi:
  mappa di indicazione n128 usata su U64; array-indici trattato come boolean mask).
  Corretti; l'integrazione FASE 2 **non è mai stata riavviata**.
- `t1b_cost_log.csv`: la riga `full` originariamente scritta con colonne disallineate
  (field mismatch del DictWriter) è stata riscritta con schema coerente; i valori sono
  identici a `run.log`/`resource_report.txt`.
- Nessun archivio originale modificato; tutti gli output in questa directory.

## 9. Inventario file
```
t1b_fase0_cost.py, fase0_time_report.txt          FASE 0 (costi statici)
t1b_dop853.py                                     driver (cache/pilot/full/validate)
matrices_n128x64.npz/.json, cache_time_report.txt cache matrici (sha256 verificato)
pilot_run.sh, pilot.log, pilot_resource_report.txt FASE 1
run_full.sh, run.log, resource_report.txt          FASE 2 (limiti duri)
reference_n128x64_dt0p00025.npz                    il riferimento (10.9 MB)
t1b_cost_log.csv                                   stima vs effettivo
run_validate.sh, t1b_validate.py, validate.log,    FASE 3
validate_resource_report.txt
t1b_plots.py                                       rigenerazione grafici
t1b_report.md                                      questo rapporto
results/
  t1b_metrics_T.csv, t1b_family_orders_T.csv,
  t1b_proxy_stability_T.csv, t1b_zones_T.csv,
  t1b_moving_max_T.csv,
  t1b_orders_time_64_128.csv   (80 passi, tutti i campi)
  t1b_orders_time_oofem_proxy128.csv (passi 1,10,20,30,39,40,50,60,70,80)
  t1b_energy_n128.csv
  t1b_error_vs_h_T.png, t1b_family_order_64_128_time.png, t1b_moving_max_rc.png
```

**STOP** — T1-B completato. Nessun T1-C, nessun ulteriore raffinamento, nessuna
simulazione OOFEM avviata; in attesa di nuova autorizzazione.
