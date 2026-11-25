# Audit indipendente della convergenza spaziale — caso accoppiato iperbolico Allen–Cahn

Data: 2026-10-27 — Directory: `consolidated_final_report/audit_convergence_next/`
Strumento: `audit_postproc.py` (post-processing di sola lettura, 11.9 s, 172 check, 0 FAIL — log `rerun_full.log`)
Vincoli rispettati: nessuna simulazione OOFEM, nessuna integrazione DOP853, nessun lavoro in parallelo, nessun pacchetto installato, nessuna matrice dense di grandi dimensioni, nessun MMS/WED, nessun dato eliminato, nessun file preesistente modificato o sovrascritto (tutti i nuovi artefatti sono in questa directory; verificato con `git status`: l'intero `consolidated_final_report/` risulta untracked, i soli file modificati nell'intero dossier restano i due già modificati a inizio sessione, estranei a questo lavoro).

Etichette usate (obbligatorie): **VERIFICATO**, **SUPPORTATO DALLE EVIDENZE**, **NON DIMOSTRATO**, **NON VERIFICABILE CON I DATI DISPONIBILI**.

---

## 1. Obiettivo e perimetro

Audit scientifico e riproducibile dei risultati **già disponibili** per il caso accoppiato
iperbolico (bulk `ρu_tt − kΔu = m(u−u³)`, Newmark β=1/4, γ=1/2; superficie a 1° ordine;
DOF traccia condivisi; RHO=KB=MB=CS=KS=MS=1; Ω=[0,2]×[0,1]; T=0.02; DT=2.5e-4; 80 passi;
IC `u₀=0.2+0.1·cos(πx/2)cos(πy)`), per stabilire:

- (A) se la dimostrazione della convergenza spaziale sia solida a livello di protocollo;
- (B) perché l'ordine asintotico non sia dimostrato e cosa mancherebbe esattamente;
- (C) se i massimi locali vicino agli angoli ammettano una lettura a singolarità;
- (D) se dati e post-processing siano indipendentemente verificabili;
- (E) quale sia il **singolo** passo minimo per rafforzare la dimostrazione.

Nessuna nuova simulazione o integrazione è stata eseguita: tutto è post-processing di
artefatti archiviati.

## 2. Fonti lette (niente di nuovo eseguito)

| Fonte | Contenuto usato |
|---|---|
| `hyperbolic/t1b_dop853_n128x64/reference_n128x64_dt0p00025.npz` | proxy n128 (81×8385 U,V; nfev 2117) |
| `hyperbolic/t1b_dop853_n128x64/results/*.npz` (n32, n64) e `hyperbolic/coupled_independent_spatial_jacobian_fix/tolcheck_n16x8_rtol2e-12.npz` | riferimenti DOP853 n32 (81×561), n64 (81×2145, nfev 1088), n16 **corretto** (81×153) |
| `hyperbolic/coupled_independent_space_n{16x8,32x16}_dt0p00025/case.out.m0.NN.vtu` (+PVD) | campi OOFEM n8/n16/n32 (n8 letto da `coupled_independent_space_n8x4…`; 30 file: 3 mesh × 10 passi {1,10,20,30,39,40,50,60,70,80}) |
| `hyperbolic/t1b_dop853_n128x64/results/t1b_{metrics_T,orders_time_64_128,orders_time_oofem_proxy128,proxy_stability_T,moving_max_T,zones_T,energy_n128}.csv` | archivio di riferimento per i confronti |
| `hyperbolic/t1b_dop853_n128x64/t1b_validate.py` (l. 60–256) | convenzioni originali: `zones()`, `metrics()`, `read_oofem_step()`, moving-max, energia |
| `hyperbolic/coupled_independent_jacobian_fix/check_jacobian_fix.py` | `assemble_fixed`, `energy_fixed` (cross-check indipendente) |
| `hyperbolic/run_coupled_independent.py` | `qshape` (derivazioni `.25`), builder mesh (pattern quads/edges) |
| `hyperbolic/coupled_independent_spatial_jacobian_fix/static_audit_report.md` §1.5, §4 | valori floor archiviati; erratum fit h²+h |
| `consolidated_final_report/{validation_report_consolidated.md, validation_matrix.csv, coupled_spatial_convergence.csv, README.md}` | punto di partenza consolidato |
| `hyperbolic/coupled_independent_spatial_jacobian_fix/independent_reference.npz` | riferimento n16 **buggy** originale (solo per provenienza) |

## 3. Metodo

`audit_postproc.py` ricalcola **da zero** (codice indipendente, nessuna importazione delle
funzioni di validazione originali) e poi incrocia con l'archivio:

1. builder mesh x-fast (n8→n128), forma Q1 con quadratura 3×3 Gauss, assemblaggio M_b, K_b, M_s, K_s;
2. lettura dei 4 riferimenti npz (integrità: forma, tempi, IC, V0, finitezza, rtol, nfev);
3. lettura dei 30 VTU OOFEM (riordino per coincidenza esatta di coordinate);
4. mappe di restrizione per **coincidenza nodale esatta** (d²<1e-26; in pratica 0.0);
5. metriche 4-tuplo (FE RMS massa-pesata normalizzata per AREA=2.0 / PERIM=6.0; max interno; max bordo);
6. energia con quadratura propria (3×3 Gauss nel bulk, 3 punti sulla superficie; esatte per quartici Q1/P1);
7. confronti contro ogni CSV archiviato (20 metriche a T, 12 ordini a T, 40 ordini a 10 istanti,
   80 passi di famiglia, 81 passi di energia, 6 frazioni di zona, 4 punti del massimo mobile,
   16 floor) e contro `fj.assemble_fixed`/`fj.energy_fixed` come seconda implementazione.

Ogni check è tracciato in `results/integrity_checks.csv` (172 righe: 166 del run base + 6
assert aggiunti dopo la correzione del bound dei max, §3.1).

### 3.1 Bug dello script di audit stesso — rilevati e corretti (documentati per trasparenza)

Quattro bug della **prima stesura** dello script di audit sono stati rilevati proprio dai
cross-check (nessuno ha contaminato i risultati finali: i CSV scritti provengono dal run
corretto, 172/172 OK). Nessuna modifica è stata fatta a file preesistenti.

| # | Bug | Rilevamento | Correzione |
|---|---|---|---|
| 1 | `qshape` con derivazioni `.5` invece di `.25` → J 2×, M_b 4× | confronto con la massa elementare esatta (M00=1/144, M01=1/288, M0,2=1/576 su quadrato 0.25) | `.25` come in `run_coupled_independent.py` |
| 2 | edges di bordo destro `(a, a+1)` (orizzontali) invece di `(a, a+nx+1)` (verticali) → ΣM_s=13.06 ≠ PERIM=6.0 | integrazione esatta 1D (∫N²=L/3) e diff di insiemi di edges | `(a, a+nx+1)` |
| 3 | composizione mappe `I64_128[I32_64[I16_32[I8_16]]]` (=I8_128) usata come mappa n8→n64 → IndexError | fallimento del cross-check (indice 2160 su 2145) | `I8_64 = I32_64[I16_32[I8_16]]` |
| 4 | S3: per le righe `max` del bound di triangolo usata la norma FE-RMS della differenza di famiglia (norma sbagliata) → 4 falsi `NO(!)` | revisione del CSV dopo il run (5 righe `NO` matematicamente impossibili per una disuguaglianza di triangolo) | norma max della differenza di famiglia sullo stesso insieme di nodi; dopo la correzione **12/12 SI**, di cui 5 uguaglianze esatte (ratio 1.000) |

## 4. Verifica A — protocollo della convergenza spaziale

**VERDICTO: VERIFICATO.**

| Proprietà del protocollo | Evidenza | Esito |
|---|---|---|
| Confronti a **tempo uguale** | tutti i riferimenti: `max|t − griglia_canonica| = 0.0` (16 check `ref_integrita` in `integrity_checks.csv`); i VTU OOFEM letti ai passi esatti 1/10/20/30/39/40/50/60/70/80 | OK |
| **Definizione di metrica coerente** | assemblaggio indipendente vs `fj.assemble_fixed`: M_b `max_rel_diff = 0.000e+00` e M_s `1.665e-16` su tutte e 4 le mesh (8 check `assembly_xcheck`); invarianti di quadratura su n128: M_b·1=2.0=AREA, M_s·1=6.0=PERIM, rowsum K_b=1.832e-15, rowsum K_s=0.0 (5 check `quadrature_n128`) | OK |
| **Stessa regione fisica** (mesh nidificate, stesso punto) | mappe di restrizione per coincidenza esatta: `d² < 1e-26` (in pratica 0.0); 30/30 file VTU riordinati con `max_d2 = 0.0` (30 check `VTU`); coordinate npz n128 == builder canonico: `max_diff = 0.0` | OK |
| **Stabilità del riferimento fine (proxy)** | non più solo "piccola differenza": ora **bound rigoroso per disuguaglianza di triangolo** `|‖e_128‖ − ‖e_64‖| ≤ ‖U_64|_k − U_128|_k‖` con la norma corretta per ogni metrica — 12/12 SI in `results/proxy_stability_bound.csv` (ratio 1.26/1.72/1.43 bulk FE, 1.23/1.55/1.80 surf FE; ratio 1.00 su 5 righe max = bound teso, 1.21 su n32 surf max). Cambiamento relativo dell'errore al cambio di proxy 64→128: 2.2% (n8), 3.5% (n16), 7.7% (n32) — tutti ≤ 8% (criterio T1) | OK |
| Metriche a T riprodotte | 24/24 valori OOFEM-vs-proxy **identici a tutte le cifre archiviate** in `t1b_metrics_T.csv` (colonne `vs_proxy128`/`vs_proxy64`); 8/8 differenze di famiglia con diff relativa ≤ 1.754e-7; 12/12 ordini a T con diff ass ≤ 4.017e-05 (`results/metrics_T_recomputed.csv`, `results/orders_T_recomputed.csv`) | OK |

### A.1 Esistenza della convergenza (distinta dall'ordine)

- **Catena OOFEM (FE RMS bulk a T):** 1.299267e-3 (h=0.25) → 7.902543e-4 (h=0.125) → 3.822430e-4 (h=0.0625): monotono decrescente; ordini a T 0.7173 e 1.0478 — **SUPPORTATO DALLE EVIDENZE** per h=0.25→0.0625, ma gli ordini sono **sotto** l'intervallo pre-registrato [1.5,2.5] (`audit_iniziale_2026-10-27.md` §7).
- **Famiglia DOP853 (FE RMS bulk a T):** 4.141362e-4 (32→64) → 4.874007e-5 (64→128): riduzione 8.5× senza plateau; la serie 80-passi (nuova, con D16_32 ricomparsa dal riferimento tolcheck n16) mostra D32_64 monotono crescente e D64_128 non monotona (§4-B). — **SUPPORTATO DALLE EVIDENZE** per h=0.0625→0.015625.
- **Accordo OOFEM–semidiscreto a mesh uguale (floor):** 2.099453e-8 (n16, step 39), 2.444096e-08 (n16, step 80), 2.597176e-08 (n32, step 39) — tutti riproducibili e coerenti con l'archivio (`static_audit_report.md` §1.5); **nuovo**: n32 step 80 = **2.916141e-08** (primo rilievo a h=0.0625 a T; non archiviato prima). 6–7 ordini sotto l'errore spaziale. — **VERIFICATO** (a h=0.125 e 0.0625); a h=0.03125 **NON VERIFICABILE CON I DATI DISPONIBILI** (manca il run OOFEM n64).

## 5. Verifica B — ordine asintotico: perché NON è dimostrato

**VERDICTO: NON DIMOSTRATO** (e ora il motivo è quantificato, non solo dichiarato).

1. **Il valore 3.0869 è un singolo rapporto locale, non un ordine.**
   `results/family_time_series_recomputed.csv` (colonne `D32_64_mine`, `D64_128_mine`, `p_high_mine`): a T il rapporto `log2(4.141362e-4/4.874007e-5) = 3.0869` è il **rapporto di due differenze calcolate in un solo istante**. Il rapporto *precedente* della stessa catena, p_low(T) = log2(D16_32/D32_64) = **0.5006**, non gli è affine: i due rapporti consecutivi del **stesso** istante differiscono di 2.6 unità.
2. **Contaminazione temporale del denominatore.**
   D64_128(s) **non è monotona** in tempo: picco allo step 44 (t=0.011) = 1.550e-4, valore a T = 4.874e-5, quindi **D(T) = 0.315 × picco** (`results/family_time_summary.csv`, riga 2). Il denominatore del rapporto 3.0869 è già 3× più piccolo del suo massimo: il rapporto è una stima locale in un regime transiente, non il limite h→0.
3. **Nessun passo soddisfa il criterio pre-registrato.**
   Su **80/80** passi almeno uno dei due rapporti consecutivi (p_low, p_high) è fuori da [1.5, 2.5] (`family_time_summary.csv`, riga 5: "step con p_low e p_high entrambi in [1.5,2.5] = 0 di 80"). Gli unici 22 passi con |p_low−p_high|<0.25 sono i primi (entrambi ≈ −0.4, regime di avvio).
4. **Ordine sulla catena OOFEM a 10 istanti: mai in [1.5,2.5].**
   `results/oofem_orders_10times_recomputed.csv`: ord16_32 bulk FE varia 0.6956→1.0478 e ord8_16 0.4713→0.7173 tra step 1 e 80; nessun valore entra in [1.5,2.5] a **nessun** istante. La catena OOFEM cresce verso l'ordine 2 ma non lo raggiunge entro h=0.0625.
5. **Richardson puntuale a nodo fisso (3–4 livelli, 10 punti × 5 istanti = 50 righe nuove):**
   `results/fixed_node_richardson.csv`. Gli ordini locali p_high salgono da 2.7 (step 1) a 4–14 (step 80) a seconda del punto: questo **non è** evidenza di ordine elevato, ma il segno che D64_128 → 0 sul nodo considerato (n64 già convergente lì), cioè il terzo livello è sotto il livello di risoluzione della differenza — il rapporto 3-livello (es. 5.3007 al nodo (1.75,0.25) a T) è pre-asintotico/transiente.
6. **Localizzazione angolare del residuo di famiglia.**
   La differenza 64→128 a T è concentrata: 65.84% dell'energia in zona A (coppia angoli attivi), 28.04% in zona B, 6.12% resto (`results/zones_recomputed.csv`, righe 1–3; replica `t1b_zones_T.csv` a ≤7.9e-8 relativa). L'errore OOFEM n32 vs proxy128: 74.61% in zona A (righe 4–6). Il "plateau angolare" (§C) e la non monotonia locale sono la stessa fisica: finché il residuo sta negli angoli, i rapporti globali non stabilizzano.

**Cosa mancherebbe per dimostrare l'ordine:** un 4° punto della catena OOFEM (n64×32) **e**
che i due rapporti consecutivi entrino in [1.5,2.5] a T (e nella maggioranza dei 10 istanti),
con floor a h=0.03125 confermato ≤ O(1e-7). Tutto e solo questo è specificato in
`next_step.md` (non eseguito).

## 6. Verifica C — massimi locali e possibile singolarità angolare

**VERDICTO: ipotesi di singolarità NON DIMOSTRATA; l'osservazione di concentrazione angolare è SUPPORTATA DALLE EVIDENZE; la non monotonia del massimo mobile è spiegata (VERIFICATO) e non richiede una singolarità.**

1. **Il "massimo mobile" non è una sequenza a punto fisso.**
   I 4 punti di `t1b_moving_max_T.csv` (riprodotto in `results/moving_max_recomputed.csv`)
   sono su **4 nodi diversi** (il primo nodo interno adiacente a ciascun angolo attivo,
   r_c = 0.35355/0.17678/0.08839/0.04419), di cui il 4° è una **differenza di famiglia**
   DOP853 64→128, non un errore OOFEM. Il fit log-log a 4 punti misti p = −0.0159 è quindi
   un fit di **punti non comparabili**: rifacendo il fit solo sui 3 punti OOFEM si ottiene
   **p = −0.9125** — il segno della "legge di potenza" cambia con la scelta del 4° punto.
   Un'osservazione che cambia segno con la selezione del campione non è evidenza di alcuna
   legge di scala: **NON DIMOSTRATO** (e dichiarato come tale, come già in T1-B).
2. **A nodo FISSO l'errore OOFEM decade monotono.**
   Al nodo (1.75, 0.25) (presente su tutte e 3 le mesh, r_c=0.35355):
   |−8.346605e-4| → |3.342395e-4| → |2.723954e-5| (`results/fixed_node_175_025_T.csv`),
   con cambio di segno n8→n16 (comportamento già documentato in T1-D). La crescita apparente
   8.35e-4 → 1.91e-3 → 2.96e-3 del "massimo mobile" è interamente dovuta al cambio di
   **punto di valutazione** (si avvicina all'angolo), non a un degrado locale: **VERIFICATO**.
3. **Concentrazione angolare (l'unica parte supportata dell'ipotesi).**
   Zone A/B (§5.6): 65.84%/28.04% della differenza di famiglia 64→128 e 74.61%/14.19%
   dell'errore OOFEM n32 stanno nei due coni r<0.5 attorno agli angoli. Questo è un fatto
   misurato e riprodotto — **SUPPORTATO DALLE EVIDENZE** come *concentrazione spaziale*
   dell'errore/residuo. **NON** è, di per sé, una dimostrazione di singolarità (es. tipo
   r^α con α<1): non c'è un profilo r^α misurato a r costante su una famiglia di mesh (che
   richiederebbe il 4° punto OOFEM e un'analisi ancora non fatta) e la soluzione del modello
   in questo test non ha motivo noto di essere singolare (IC liscia, dominio poligonale,
   condizioni omogenee — ma il bulk è iperbolico con Newmark: la regolarità non è teorema
   invocato qui).
4. **Distinzione richiesta:** *osservazione* (concentrazione angolare, non monotonia del
   massimo mobile) → documentata; *ipotesi* (singolarità angolare che degrada l'ordine
   locale) → **NON DIMOSTRATA**; *tesi dimostrata* → nessuna.

## 7. Verifica D — verifica indipendente di dati e post-processing

**VERDICTO: VERIFICATO. Nessun problema di dati o di post-processing trovato nell'archivio; 4 bug dello script di audit (primo draft) rilevati e corretti (§3.1).**

- **Integrità riferimenti (16 check):** forme 81×{153,561,2145,8385}; tempi == griglia
  canonica (0.0); IC esatta (0.0); V0 = 0; tutti finiti; rtol = 2e-12; nfev n64 = 1088,
  n128 = 2117. Il conteggio nodi **8385** (non 8449) è confermato direttamente sulle shape.
- **Provenienza del riferimento n16 (nuovo, §5.4 del report consolidato):** il tolcheck n16
  usato da T1-B è il riferimento **corretto**, non quello buggy: OOFEM n16 vs riferimento
  buggy originale = **3.452e-3** (massimo nodale, ordine atteso O(1e-3) del bug di
  Giacobiano — `integrity_checks.csv` riga `provenienza`), mentre vs il tolcheck il floor è
  **2.444e-8** (6 ordini sotto): se il tolcheck fosse il buggy, il "floor" sarebbe O(1e-3).
- **Confronti completi con l'archivio:** 20 metriche a T + 12 ordini a T + 40 ordini a 10
  istanti (40/40, abs_diff ≤ 4.822e-05) + serie 80 passi (D32_64 e D64_128 identici a
  `t1b_orders_time_64_128.csv`; D16_32 **nuovo** — non esiste in archivio) + 81 passi di
  energia + 6 frazioni di zona + 4 punti moving-max + 16 floor. Max diff energia
  `|E_mine − E_arch| = 4.992e-11` su 81 step (limitato dalle 7–10 cifre archiviate) e
  `|E_mine − fj.energy_fixed| = 7.216e-15` a s=0,40,80.
- **Coerenza con la storia degli errori** (nessuna contraddizione nuova):
  (i) bug di Giacobiano → la catena dei riferimenti corretti è confermata dalla provenienza
  di cui sopra; (ii) off-by-one T1 (step 39 vs 40) → l'audit legge **entrambi** i passi e
  riporta floor a entrambi; (iii) erratum fit h²+h → il 3.0869 attuale non è un fit ma un
  singolo rapporto e §5-B lo qualifica come non dimostrativo (la retrazione resta);
  (iv) corruzione archivio n32 → l'audit usa i VTU T1-A (steps 1–39 bit-identici
  all'originale) e il floor n32 a T (2.916e-8) è coerente con i floor n16 (2.1–2.4e-8):
  nessuna firma di corruzione residua; (v) 8449→8385 → shape 81×8385 confermata;
  (vi) inflazione √2 → le distanze r_c dei punti angolari sono quelle corrette
  (es. (1.75,0.25)→(2,0): 0.353553 = √(0.25²+0.25²)).

## 8. Le quattro proprietà non equivalenti (richiesta dal task)

| Proprietà | Stato in questo audit |
|---|---|
| 1. Convergenza **osservata** della famiglia FE (l'errore diminuisce al raffinare) | **SUPPORTATA DALLE EVIDENZE** — FE RMS monotono: OOFEM fino a h=0.0625, famiglia DOP853 fino a h=0.015625; non è una prova di ordine |
| 2. Accordo **OOFEM–semidiscreto** (stessa mesh) | **VERIFICATO** — floor 2.1e-8…2.9e-8, 6–7 ordini sotto l'errore spaziale |
| 3. **Ordine asintotico** (p=2 su un intervallo di h, rapporto stabile) | **NON DIMOSTRATO** — 0/80 passi nel criterio; rapporti consecutivi discordi a T (0.5006 vs 3.0869); denominatore non monotono in tempo |
| 4. **Accuratezza vs soluzione esatta della PDE** | **NON VERIFICABILE CON I DATI DISPONIBILI** — nessuna soluzione esatta nota per questo caso; il proxy n128 è una discretizzazione (la sua "verità" è limitata dal suo errore, bound di triangolo §4) |

## 9. Energia (controllata nell'audit, non oggetto principale)

**VERIFICATO** (conferma, non estensione): 81 step ricalcolati con quadratura indipendente;
monotona decrescente (max ΔE consecutivo = −2.737e-7); E0 = 1.9327101889 →
ET = 1.9274638141 (−0.27145%); cinetica 0 → 2.4766801329e-3; max|Δ vs archivio| = 4.992e-11,
max|Δ vs `fj.energy_fixed`| = 7.216e-15 (`results/energy_recomputed.csv`,
`results/energy_summary.csv`). Nessuna identità discreta è dimostrata (Newmark nonlineare):
la monotonia resta una proprietà osservata.

## 10. Criterio di decisione (task §5)

**Opzione 2 — la dimostrazione è rafforzata con post-processing di dati già disponibili**
(eseguito in questo audit: bound di triangolo rigoroso, convergenza a nodo fisso,
Richardson 3-livello, serie D16_32 a 80 passi, floor n32 a T, decomposizione nodale,
provenienza del riferimento n16, verifica indipendente completa). Il **residuo** — un
ordine asintotico dimostrato — non è raggiungibile con i soli dati esistenti perché la
catena OOFEM si ferma a h=0.0625: serve **un solo** nuovo dato, specificato senza eseguirlo
in `next_step.md` (run OOFEM n64×32, criterio di accettazione pre-registrato).

## 11. Limiti dell'audit

- Tutto è riferito a **un singolo istante finale T=0.02** (più 10 istanti per gli ordini
  OOFEM e 80 per la famiglia) e a un'unica combinazione di parametri.
- Il proxy n128 non è la soluzione esatta: i verdetti su "convergenza" sono relativi alla
  famiglia di discretizzazioni (proprietà 1 e 2 della tabella §8), mai alla PDE (proprietà 4).
- I 4 campi OOFEM provengono da 4 run separati (n8/n16/n32 in T1, n32 T1-A): la
  riproducibilità bit-esatta è documentata ma non ri-verificata qui per i campi n8/n16
  (solo i VTU T1-A di n32 hanno il rerun formale; per n8/n16 l'audit verifica solo la
  coerenza geometrica e i valori a tutte le cifre archiviate).
- L'energia del bulk iperbolico usa la ricostruzione di v già discussa in T1-B
  (pavimento ~7e-7): l'audit conferma i valori, non rimuove quel limite.
- Nessun calcolo costoso: 172 check in 11.9 s; picco di memoria modesto (campi 81×8385).

## 12. Verdetto sintetico

| # | Affermazione | Etichetta |
|---|---|---|
| V1 | Protocollo (tempi uguali, metriche coerenti, regione fisica identica, proxy stabile con bound rigoroso) | **VERIFICATO** |
| V2 | Convergenza osservata famiglia OOFEM FE RMS, h=0.25→0.0625 | **SUPPORTATO DALLE EVIDENZE** |
| V3 | Convergenza osservata famiglia DOP853 FE RMS, h=0.0625→0.015625 (8.5×, no plateau) | **SUPPORTATO DALLE EVIDENZE** |
| V4 | Accordo OOFEM–semidiscreto a mesh uguale (floor ≤ 2.9e-8) | **VERIFICATO** |
| V5 | Ordine asintotico p=2 | **NON DIMOSTRATO** (motivo quantificato §5) |
| V6 | Il rapporto 3.0869 come stima di ordine | **NON DIMOSTRATO** (rapporto locale singolo, denominatore non monotono) |
| V7 | Non monotonia temporale di D64_128 (picco step 44, D(T)=0.315×picco) | **VERIFICATO** |
| V8 | Concentrazione angolare di errore/residuo (65.8%/74.6% in zona A) | **SUPPORTATO DALLE EVIDENZE** |
| V9 | Singolarità angolare | **NON DIMOSTRATO** (ipotesi; fit fragile p=−0.016/−0.912) |
| V10 | La non monotonia del massimo mobile è scelta del punto di valutazione, non degrado locale | **VERIFICATO** |
| V11 | Dati e post-processing indipendentemente verificabili; nessun problema in archivio | **VERIFICATO** |
| V12 | Bilancio energetico T1-B (monotono, −0.27145%) | **VERIFICATO** |
| V13 | Accuratezza vs soluzione esatta della PDE | **NON VERIFICABILE CON I DATI DISPONIBILI** |
| V14 | Floor del solutore a h=0.03125 | **NON VERIFICABILE CON I DATI DISPONIBILI** (gap → `next_step.md`) |

## 13. Inventario di questa directory

| File | Contenuto |
|---|---|
| `audit_postproc.py` | script di audit (sola lettura, ~770 righe) |
| `rerun_full.log` | log completo dell'ultimo run (172 check, 0 FAIL, 11.9 s) |
| `audit_report.md` | questo report |
| `findings.csv` | tabella sintetica verifica → evidenza → conclusione → etichetta |
| `next_step.md` | unico passo successivo raccomandato (non eseguito, richiede autorizzazione) |
| `results/integrity_checks.csv` | 172 check con esito |
| `results/metrics_T_recomputed.csv` (24), `orders_T_recomputed.csv` (12) | metriche e ordini a T vs archivio |
| `results/proxy_stability_bound.csv` (12) | bound di triangolo per mesh×metrica (12/12 SI) |
| `results/zones_recomputed.csv` (6), `moving_max_recomputed.csv` (4+2 fit) | zone angolari e massimo mobile |
| `results/fixed_node_175_025_T.csv`, `fixed_node_richardson.csv` (50) | convergenza a nodo fisso + Richardson 3–4 livelli su 10 punti × 5 istanti |
| `results/floors_recomputed.csv` (16) | floor OOFEM vs DOP853 stessa mesh (incl. n32@T nuovo) |
| `results/family_time_series_recomputed.csv` (80), `family_time_summary.csv` | serie di famiglia a 80 passi (D16_32 nuovo) + statistiche |
| `results/energy_recomputed.csv` (81), `energy_summary.csv` | energia ricalcolata indipendente |
| `results/node_decomp_n32_T.csv` | decomposizione nodale |e| = floor + famiglia a n32@T |
| `results/oofem_orders_10times_recomputed.csv` (40) | ordini catena OOFEM a 10 istanti vs archivio |
