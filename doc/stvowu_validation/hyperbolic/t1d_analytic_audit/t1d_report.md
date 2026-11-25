# T1-D — Audit analitico: conclusioni sostenute dai dati esistenti e prossimo test più informativo

**Data:** 2026-10-27
**Ambito:** esclusivamente analisi su dati già archiviati/calcolati. Nessuna nuova
integrazione o simulazione, nessuna modifica o sovrascrittura di archivi originali,
nessun processo parallelo, nessuna conversione sparse→dense (1 core, ~30 s di CPU
totali). Directory nuova: `t1d_analytic_audit/`.

---

## 1. Dimensione caratteristica h e calcolo degli ordini — **VERIFICATA**

| mesh | h_x = 2/nx | h_y = 1/ny | lati quadrilateri | lunghezze linee | quadrati? | linee = h? |
|---|---|---|---|---|---|---|
| n8×4 | 0.25 | 0.25 | [0.25, 0.25] | [0.25, 0.25] | sì (a 1e-14) | sì |
| n16×8 | 0.125 | 0.125 | [0.125, 0.125] | [0.125, 0.125] | sì | sì |
| n32×16 | 0.0625 | 0.0625 | [0.0625, 0.0625] | [0.0625, 0.0625] | sì | sì |

- `h` è **non ambigua**: h_x = h_y in entrambe le direzioni, e coincide con la lunghezza
  di ogni lato (quadrilateri e linee di superficie), quindi è la dimensione
  caratteristica sia per bulk che per superficie.
- Rapporto di raffinamento esatto: h8/h16 = h16/h32 = **2** → la base **log2** usata in
  `order = log2(e_grosso/e_fine)` è la base corretta.
- Ogni ordine è un **rapporto a tempo uguale** (le due mesh della coppia sono valutate
  a t = s·dt, §2). Nessun ordine usa istanti diversi.
- Conclusione: nessuna ambiguità geometrica nella definizione di h né nella base degli
  ordini. (`t1d_h_geometry.csv`.)

## 2. Dipendenza temporale degli ordini — istanti archiviati soltanto

Metodo: a ogni step s ∈ 1…80 (1…39 per n32×16 originale; step 80 da T1-A) le due mesh
della coppia sono prese **allo stesso tempo** t = s·dt (convenzione verificata in T1 e
nell'audit statico §1.1: OOFEM step s ↔ U[s], Δt ≤ 3.5×10⁻¹⁸). Dove si confronta un
ordine a un istante con un ordine a un altro istante (es. 0.76 → 0.98 tra t=0.00975 e
t=0.02) si tratta di **confronto temporale, esplicitamente segnato** — non è un ordine di
raffinamento.

Range su tutti gli istanti archiviati (serie complete:
`t1d_order_time_series_816.csv`, 80 righe; `t1d_order_time_series_1632.csv`, 39 righe):

| intervallo (t pari) | metrica | min (step) | max (step) | step 1 | ultimo |
|---|---|---|---|---|---|
| 8→16 | bulk FE | 0.570 (19) | 0.698 (80) | 0.686 | 0.698 |
| 8→16 | bulk max | −1.308 (3) | −0.810 (1) | −0.810 | −1.172 |
| 8→16 | surf FE | 0.185 (19) | 0.367 (80) | 0.312 | 0.367 |
| 8→16 | surf max | −0.078 (13) | 0.355 (80) | 0.223 | 0.355 |
| 16→32 | bulk FE | 0.730 (24) | 0.781 (2) | 0.772 | 0.761 (39) |
| 16→32 | bulk max | −1.704 (1) | −1.402 (39) | −1.704 | −1.402 (39) |
| 16→32 | surf FE | 0.625 (3) | 0.718 (39) | 0.631 | 0.718 (39) |
| 16→32 | surf max | 0.506 (8) | 0.753 (39) | 0.586 | 0.753 (39) |

Step 80 a T=0.02 (da T1-A, `t1a_convergence_orders.csv`): 16→32 bulk FE **0.983**,
surf FE **1.106**, bulk max **−0.732**, surf max **1.089**.

**Risultati:**

1. **Nessun recupero verso l'ordine 2 a nessun istante**: l'ordine bulk FE 8→16 non
   supera mai 0.70 nell'intero intervallo di tempo; 16→32 sta in 0.73–0.78 per
   t ≤ 0.00975. Il mancato raggiungimento del criterio [1.5, 2.5] **non è un artefatto
   della scelta dell'istante**.
2. **La concentrazione angolare è presente già al primo step** (verifica nuova,
   `t1d_early_steps_zones.csv`): la quota di energia di errore nella zona A (r<0.5 da
   (0,1)/(2,0)) è 0.735–0.766 a step 1 per tutte e tre le mesh, e 0.71–0.75 a step 5 e
   20 — quasi identica a quella a T (0.708–0.739). Al step 1 l'errore totale è
   ~10⁻⁷ (u(t≈dt) ≈ u₀, ancora regolare): la struttura angolare dell'errore **non è una
   accumolazione tarda nel tempo**, è presente dal primo passo.
3. **Debole dipendenza temporale**: 8→16 ha un minimo depresso a t≈0.00475–0.005
   (0.570 bulk FE, 0.185 surf FE) e risale verso T; 16→32 **cresce col tempo** verso T
   (bulk FE 0.76→0.98, surf FE 0.72→1.11 tra t=0.00975 e t=0.02) [confronto temporale,
   segnato].
4. **Il bulk max è non monotono a ogni istante archiviato** (ordini negativi in tutta la
   serie, entrambi gli intervalli): la non monotonia non è sporadica. Il meccanismo è
   quantificato in §4.

Limite dichiarato: la floor temporale (OOFEM vs DOP853 a mesh uguale) è misurata a
step 39/80 (2.1–2.9×10⁻⁸), non a step 1; essa vale il 3–6% dell'errore di step 1, quindi
le quantità a step 1 hanno una contaminazione temporale di quella entità possibile.
Tutte le affermazioni robuste di questo § usano step ≥ 5 (errori ≥ 1.5×10⁻⁵,
100× sopra la floor).

## 3. Concentrazione agli angoli vs. prova di singolarità

**Fatti osservati (dai dati, con file di riferimento):**
- ~71–74% dell'energia di errore nella zona angolare A su tutte e tre le mesh
  (0.708/0.725/0.739 a T), stabile col raffinamento (`t1a_region_analysis_T.csv`).
- **La stessa concentrazione appare nell'auto-differenza di livello del proxy**
  (DOP853 n32×16 vs DOP853 n64×32, nessuna OOFEM coinvolta): zona A 0.722
  (t=0.00975) / 0.739 (t=0.02); B 0.184/0.157; resto 0.094/0.104
  (`t1d_proxy_selfdiff_zones.csv`). Poiché la floor OOFEM (~10⁻⁸) è 4 ordini sotto
  l'errore spaziale, OOFEM-vs-proxy ≈ famiglia-vs-famiglia: **la struttura dell'errore
  è proprietà della famiglia FE per questo problema, non dell'implementazione OOFEM.**
- Simmetria 180° esatta (max|Δ| = 2.2×10⁻¹⁶, precisione macchina).
- Decadimento di potenza |err| ~ r_c^p verso l'angolo, p ≈ −2 negli interni
  (regione-dipendente; audit statico §1.3).
- **Convergenza puntuale a nodo fisso** (verifica nuova, `t1d_fixed_nodes.csv`), a T:
  nodo (0.25, 0.75): 8.34×10⁻⁴ (n8) → 3.34×10⁻⁴ (n16) → 2.74×10⁻⁵ (n32);
  nodo (0.125, 0.875): 1.88×10⁻³ (n16) → 6.08×10⁻⁴ (n32); con cambio di segno al primo
  raffinamento. I gemelli 180° sono uguali a ~10⁻¹⁷.
- L'errore del **primo nodo interno** (più vicino all'angolo) **cresce** al
  raffinamento: 8.34×10⁻⁴ → 1.88×10⁻³ → 3.12×10⁻³ a T, scaling err ~ r_c^−0.952
  (§4).

**Non dimostrato (limiti interpretativi, dichiarati):**
- La "singolarità" (soluzione esatta non regolare nell'angolo) **non può essere
  stabilita da dati solo-FE**: il proxy è esso stesso una discretizzazione; la soluzione
  esatta è sconosciuta; ogni "errore" è una differenza rispetto al proxy. In particolare
  l'"errore" di OOFEM-n32 rispetto al proxy n64 è in sostanza la **differenza di livello
  di famiglia 32→64**, non l'errore rispetto alla soluzione esatta.
- La concentrazione in A è una condizione necessaria (ma non sufficiente) dell'ipotesi
  di limite non regolare: **non è presentata qui come dimostrazione di singolarità.**
- Il fit e(h) = a·h² + b·h **non è supportato** (a < 0 a tempo uguale, T1-A §5):
  nessuna decomposizione Ch + Ch² è identificabile con i dati disponibili.

**Cosa i dati sostengono:** (i) la soluzione FE converge puntualmente a punti fissi
(compresi punti a distanza finta dall'angolo, r = 0.177); (ii) l'errore di
discretizzazione non svanisce al tendere dei nodi all'angolo (cresce ~ r^−0.95 a 3
livelli); (iii) la struttura è di livello famiglia, non dell'implementazione.
L'ipotesi di **limite non regolare (singolare) agli angoli (0,1)/(2,0)** resta
l'interpretazione principale, coerente con tutti i fatti sopra — ma resta
**interpretazione** finché non esiste una risoluzione indipendente (ordini di famiglia
a h più fini, §5.1).

## 4. Massimi non monotoni: nodi dominanti

Argmax interno per ogni (mesh, istante) disponibile (`t1d_max_nodes.csv`):

| mesh (step, t) | nodo argmax interno | err | r_c | == nodo più vicino all'angolo? |
|---|---|---|---|---|
| n8×4 (80, 0.02) | (0.25, 0.75) | −8.345×10⁻⁴ | 0.3536 | sì (twin 180° a pari errore) |
| n16×8 (80, 0.02) | (1.875, 0.125) | −1.880×10⁻³ | 0.1768 | sì (twin a pari errore) |
| n32×16 (39, 0.00975) | (0.0625, 0.9375) | −1.281×10⁻³ | 0.0884 | sì (twin a pari errore) |
| n32×16 nuovo (80, 0.02) | (1.9375, 0.0625) | −3.123×10⁻³ | 0.0884 | sì (twin a pari errore) |

- **Il nodo dominante è in tutti i casi il nodo interno più vicino a un angolo**
  (con il twin 180° a errore identico): i dati lo consentono e lo identificano.
- r_c = h√2 si dimezza a ogni raffinamento; |err| su quel nodo a T scala come
  **err ~ r_c^−0.952 ≈ 1/r_c** (3 punti: r = 0.3536/0.1768/0.0884).
- **Meccanismo (fatti + costruzione): massimo mobile.** L'argmax "insegue" l'angolo al
  raffinamento: il bulk max osservato è misurato su nodi diversi a mesh diverse. La non
  monotonia del bulk max (−1.17, −0.73 ordini) **non è una violazione di convergenza**:
  coesiste con la convergenza puntuale a nodo fisso (§3: 8.3×10⁻⁴ → 3.3×10⁻⁴ →
  2.7×10⁻⁵ su (0.25, 0.75)). È la firma classica FE di una limitazione di regolarità
  localizzata — ma, come in §3, resta interpretazione finché il comportamento r→0 non è
  risolto (richiede r_c < 0.088, §5.1).
- Il **surf max è monotono decrescente** (6.07×10⁻³ → 4.74×10⁻³ → 2.23×10⁻³ a T):
  la non monotonia riguarda solo il bulk max.
- Limite: il trend r_c^−0.95 è su 3 livelli (r_c ≥ 0.0884); se l'errore continui a
  crescere o saturi per r→0 **non è identificabile coi dati esistenti** (dichiarato).

## 5. Test candidati — valore informativo e costo (nessuno eseguito)

### 5.1 Proxy DOP853 n128×64 (T1-B)
- **Valore informativo: il più alto — è l'unico test che discrimina le due spiegazioni
  in competizione dell'ordine ridotto:**
  (i) *effetto finite-h / pre-asintotico*: se l'ordine di famiglia 64→128 recupera verso
  2, gli ordini 0.6–1.1 misurati sono un effetto di h finito e la discretizzazione è
  sana;
  (ii) *limite di regolarità del problema accoppiato*: se 64→128 resta ≈ 1, il limite
  è del problema (non regularità agli angoli).
  Aggiunge inoltre: raffinamento del proxy per tutti i confronti (oggi il proxy è solo
  2× più fine del più fine OOFEM); dati per-nodo a r_c = 0.0442 e 0.0221 → risponde
  alla domanda §4 (crescita continua verso r→0 o saturazione); serie di ordini di
  famiglia a istanti precoci su un livello nuovo.
- **Costo stimato** (stima, da confermare con pilot): 8449 DOF, 81 step, matrici sparse
  (stesso codice via di n64×32, già eseguito senza problemi): **~30–60 min su 1 core,
  RAM dell'ordine di qualche centinaio di MB**. Directory nuova, npz nuovo, nessun file
  preesistente toccato.
- **Precondizioni**: autorizzazione esplicita del costo; **pilot vincolato** (primi 10
  step, monitoraggio RSS, arresto se fuori budget — protocollo FASE 2); 1 processo.

### 5.2 Test bulk-only (senza superficie)
- **Valore informativo: medio.** Isola l'equazione bulk (ρu_tt − kΔu = m(u−u³)) dalla
  coupling: se bulk-only raggiunge ordine ~2 e il caso accoppiato no, la riduzione è
  attribubile all'accoppiamento (comportamento angolare del sistema accoppiato); se
  anche bulk-only riduce, il problema è nel setup/discretizzazione bulk.
- **Costo stimato**: 3 run OOFEM (n8/n16/n32: ~0/2/10 s, misurati) + 3 run DOP853
  bulk-only (45–561 DOF: secondi) + variante di riferimento + directory input nuove:
  **~15–30 min** (CPU + setup).
- **Rischio dichiarato**: la fattibilità dell'input **non è verificata
  staticamente** — va prima controllato se il modello OOFEM accetta un caso senza
  elementi di superficie (nbc=0, termini di superficie assenti). Se non lo accetta, il
  test richiederebbe modifiche a modello/formato di input, vietate dalla vincolistica
  (punto 3) senza aver prima provato un errore specifico. Non immediatamente esegubile
  senza quel controllo.

### 5.3 Ulteriore raffinamento OOFEM (n64×32)
- **Valore informativo: limitato e parziale.** Aggiunge: (i) la **floor OOFEM a n64**
  (OOFEM-n64 vs DOP853-n64, stessa mesh e tempo): verifica che l'accordo ~10⁻⁸ persista
  a h più fine — assicurazione utile sul verdetto di accordo; (ii) un quarto punto
  per-nodo del trend del massimo mobile (r_c = 0.0442). **Non aggiunge un nuovo ordine
  di convergenza**: OOFEM-n64 vs proxy-n64 è la floor del solvere, non un errore
  (servirebbe la differenza di livello del proxy, che è esattamente 5.1).
- **Costo stimato**: il più basso — run OOFEM ~40–90 s (misurato n32: 10 s, RSS 32 MB;
  n64 ~4× DOF, RAM attesa 100–200 MB) + ~5 s di post-processing; 1 core, directory nuova.

### 5.4 Ranking
1. **5.1 (DOP853 n128×64)**: unico test decisivo sulla domanda interpretativa; nessuna
   riesecuzione OOFEM né modifica di codice; migliora anche tutti i confronti esistenti
   (proxy più fine). **Raccomandato come prossimo test.**
2. 5.3 (OOFEM n64): assicurazione economica sul verdetto di accordo a h più fine e un
   punto per-nodo in più; secondaria; esegubile se 5.1 non viene autorizzata (costo
   ~1–2 min), ma non risponde alla domanda di convergenza.
3. 5.2 (bulk-only): ha senso **dopo** l'esito di 5.1 (per localizzare la causa se il
   limite di regolarità fosse confermato); bloccata dalla fattibilità input da verificare.

## 6. Verdetti (mantenuti separati, come richiesto)

| verdetto | stato | evidenza |
|---|---|---|
| **Riproducibilità** | **SUPPORTATA** | Rierogazione T1-A bit-identica agli step intatti 1–39 (max\|Δ\| = 0), deterministica a 1 core; 16/16 metriche ricalcolate end-to-end con diff 0 (audit statico §1.1). |
| **Accordo col riferimento semidiscreto** | **SUPPORTATA** | Floor 2.1/2.4/2.6/2.9×10⁻⁸ (bulk FE RMS) su 4 coppie (mesh, step), a mesh, tempo e convenzioni identiche e verificate; 6–7 ordini sotto l'errore spaziale → la discretizzazione temporale non contamina le metriche spaziali. |
| **Convergenza spaziale** | **NON VALIDATA** (criterio preliminare [1.5, 2.5]) | A T: FE 0.698/0.983 (bulk), 0.367/1.106 (surf); l'ordine non raggiunge 1.5 a **nessun** istante archiviato (§2); bulk max non monotono a ogni istante (meccanismo quantificato, §4); FE RMS monotono decrescente. Il fallimento è stabile nel tempo → non artefatto della scelta dell'istante. |

## 7. Gap di dati dichiarati esplicitamente

1. **Nessuna soluzione esatta** del problema accoppiato (vincolo 1: no MMS; non
   disponibile in letteratura per questo sistema): ogni "errore" è rispetto al proxy FE
   DOP853 n64×32; l'"errore" a n32 è in sostanza la differenza di livello di famiglia
   32→64.
2. Dati per-nodo OOFEM solo per n8@80, n16@80, n32@39, n32_new@80 (no per-nodo n16/n8 a
   step 39 → il controllo a nodo fisso a t=0.00975 è possibile solo su n32).
3. Un solo campo esportato: `u` (`Temperature`); nessun `v`.
4. Archivio originale n32×16: step 40–80 corrotti (intatto, non modificato); T=0.02
   recuperato da T1-A in directory separata.
5. Nessun run OOFEM su n64×32: l'errore OOFEM a h = 0.03125 è sconosciuto.
6. Il trend err ~ r_c^−0.95 del primo nodo interno è su 3 livelli (r_c ≥ 0.0884): il
   comportamento r→0 non è identificabile coi dati esistenti.
7. La floor temporale non è misurata a step 1 (valore a step 39/80 usato come
   superiore, §2).

## 8. File

**Usati (sola lettura, nessun archivio modificato):**
- `hyperbolic/coupled_independent_space_n{8x4,16x8,32x16}_dt0p00025/` (pvd + VTU);
- `hyperbolic/t1a_n32x16_rerun/oofem/` (VTU/pvd del rerun) e `…/results/` (CSV T1-A);
- `hyperbolic/coupled_independent_spatial_jacobian_fix/`:
  `reference_n64x32_dt0p00025.npz`, `reference_n32x16_dt0p00025.npz` (cache DOP853),
  3 CSV per-nodo, `spatial_convergence_final.csv`, `convergence_orders.csv`,
  `static_audit_report.md`, `spatial_verification_report.md`;
- moduli `run_coupled_independent.py` (`mesh`, `rms_fe`) e
  `check_jacobian_fix.py` (`assemble_fixed`) — senza modifiche.

**Prodotti in `t1d_analytic_audit/` (directory nuova):**
`t1d_h_geometry.csv`, `t1d_order_time_series_816.csv`,
`t1d_order_time_series_1632.csv`, `t1d_early_steps_zones.csv`,
`t1d_proxy_selfdiff_zones.csv`, `t1d_max_nodes.csv`, `t1d_fixed_nodes.csv`,
questo rapporto.

---

## Raccomandazione

**Prossimo test con il massimo valore informativo: T1-B — proxy DOP853 n128×64**
(~30–60 min, 1 core, RAM ~centinaia di MB da confermare con pilot vincolato; nessuna
riesecuzione OOFEM, nessuna modifica di codice, directory nuova). È l'unico test
disponibile che discrimina *effetto finite-h* da *limite di regolarità angolare*,
risponde alla domanda del massimo mobile (r_c fino a 0.0221) e raffina il proxy per
tutti i confronti. Se non autorizzato: alternativa economica T-5.3 (OOFEM n64×32,
~1–2 min) che rafforza solo il verdetto di accordo; il verdetto di convergenza spaziale
resterebbe invariato (NON VALIDATA) con i limiti dichiarati.

**Fine di T1-D. In attesa di istruzioni; nessun test avviato.**
