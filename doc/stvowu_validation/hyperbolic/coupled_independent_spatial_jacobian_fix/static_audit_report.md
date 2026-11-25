# Audit statico e di post-processing di T1 — rapporto

**Data:** 2026-10-27
**Ambito:** solo ispezione statica e post-processing su dati già archiviati/calcolati.
Nessuna nuova integrazione, nessuna modifica/sovrascrittura di dati originali,
nessuna matrice densa, nessun processo parallelo (1 core, ~10 s di CPU totali).
File nuovi prodotti: `archive_integrity_audit.csv`, `reproducibility_check.csv`,
`fit_sensitivity_audit.csv`, `region_analysis_audit.csv`, questo rapporto.

---

## 1. Verifiche eseguite ed evidenze

### 1.1 Norme, interpolazione, tempi PVD/VTU, sincronizzazione — **CONFERMATI**

| verifica | esito |
|---|---|
| Definizione norme | `rms_fe(M,e,measure) = sqrt(eᵀMe/measure)` (FE L² RMS); bulk max sui nodi interni, surf max sui nodi di confine (masche); AREA=2.0, PERIM=6.0 |
| Interpolazione | coincidenza nodale **esatta** sulle griglie annidate: `d² < 1e-26` per nodo (assert in `fine_to_coarse_map`); distanza minima massima = 0.00×10⁰; nessuna interpolazione con funzioni forma |
| Griglie temporali | `TIMES[s] = s·dt` **esatto** (max diff 0.00e+00); tempi PVD vs `TIMES[1..80]`: max diff **3.469×10⁻¹⁸** su tutti e tre gli archivi; `read_step_oofem` risolve il DataSet con `|t_pvd − step·dt| < 1e-12` |
| Ordinamento nodi | `array_equal(orig.mesh, assemble_fixed) = True`, max diff 0.00e+00 su n8/n16/n32/n64: campo OOFEM (ordine VTU) e riferimento (indice `U[s]`) usano lo stesso ordinamento |
| Semantica del campo | i VTU contengono un unico campo, `Temperature` (= u, variabile primaria); l'identificazione è sostenuta dalla floor 10⁻⁸ (un campo diverso darebbe differenze O(0.1)) e dal controllo IC nodo-per-nodo di `case.in` |
| **Riproducibilità end-to-end** | tutte e 16 le metriche di `spatial_convergence_final.csv` ricalcolate da zero (riparsing VTU, nuova restrizione, indici a tempo uguale): **max \|ri_calcolato − salvato\| = 0.000×10⁰** |

Dettagli: `reproducibility_check.csv`.

### 1.2 Integrità archivio n32×16 — **MAPPATA ESAATTAMENTE** (sola lettura)

| archivio | VTU integri | corrotti (0 byte) | pvd | mesh.npz | run.log | case.in |
|---|---|---|---|---|---|---|
| n8×4 | 80/80 | — | 80 DataSets, t=k·dt esatti | 3318 B ok | 28154 B ok | ok |
| n16×8 | 80/80 | — | 80 DataSets, t=k·dt esatti | 8694 B ok | 28154 B ok | ok |
| n32×16 | **39/80 (step 1–39, t ≤ 0.00975)** | **41 file consecutivi: step 40–80 (t = 0.01000…0.02000)** | 80 DataSets, t=k·dt esatti (intatto, riferisce tutti i file) | **0 B (corrotto)** | **0 B (corrotto)** | 95890 B ok |

Tutti i VTU integri: parsabili, numero punti corretto, coordinate identiche alla mesh
rigenerata (atol 1e-13), `Temperature` finito. Unico campo in ogni VTU: `Temperature`.
Dati per file (240 righe): `archive_integrity_audit.csv`. **Nessun dato modificato.**

### 1.3 Sensibilità dei fit angolari — **IL FIT h²+h NON È ROBUSTO (recesso)**

Fit `e = a·h² + b·h` sui 3 punti FE RMS (n8, n16 a t=0.02; n32 a t=0.00975 — **tempi misti**,
l'unico set a 3 mesh disponibile):

| metrica | a (h²) | b (h) | residuo max |
|---|---|---|---|
| bulk FE RMS | **−2.042×10⁻³** | 5.900×10⁻³ | **71.9%** |
| surf FE RMS | **−2.524×10⁻²** | 1.680×10⁻² | **68.5%** |

Il coefficiente di h² è **negativo** e il fit è scarno: sul dato corretto a tempo uguale il
modello h²+h non è supportato. Il fit "residuo 9%" citato in §6 di
`spatial_verification_report.md` (a=+9.48×10⁻³, b=+2.83×10⁻³) era calcolato sui valori
**pre-correzione** (con l'offset di un passo) e **non è riproducibile sui dati corretti:
si recede da esso come evidenza** (erratum §4 di questo audit).

**Quantità osservate robuste** (stesso tempo, due mesh per intervallo):

| intervallo | bulk FE RMS | surf FE RMS |
|---|---|---|
| 8→16 @ t=0.02 | ordine 0.698 | ordine 0.367 |
| 16→32 @ t=0.00975 | ordine 0.761 | ordine 0.718 |

Contrasto log-log a un solo esponente (3 punti, tempi misti): p = 1.331 (bulk), 1.094 (surf).

**Esponente \|err\| ~ r_c^p** (nodi interni, distanza dal paio di angoli (0,1)/(2,0),
soglia \|err\|>10⁻⁶), in funzione del taglio di regione:

| mesh | r_c<0.25 | r_c<0.33 | r_c<0.5 | r_c<0.75 |
|---|---|---|---|---|
| n32×16 (n=16/34/80/146) | −2.73 | −1.95 | −2.65 | −2.10 |
| n16×8 (n=6/16/44) | −5.15* | −2.21 | — | −2.85 |
| n8×4 (n=8) | — | — | — | −3.78* |

\* pochi punti, instabile. Bordo (n32): −0.14/−0.19/−0.30/−0.62 (molto più piano).
Distanza dai **4** angoli (n32 interno): −2.05 / −2.77. L'esponente dipende da taglio e
soglia: si può affermare solo che **\|err\| decaye verso l'angolo con legge di potenza di
indice ≈ 2 (range −1.9…−2.8 interno)**; l'indice esatto non è identificabile coi dati
disponibili. (Il valore −1.75 precedentemente citato era su dati pre-correzione e con
soglia diversa.)

**Simmetria 180° (fatto osservato, non ipotesi):** max\|err(p) − err(2−p)\| = 2.2×10⁻¹⁶
(n32), 6.6×10⁻¹⁶ (n16), 6.7×10⁻¹⁶ (n8) — precisione macchina.

Dettagli: `fit_sensitivity_audit.csv`.

**Osservato vs ipotesi interpretative (distinguere esplicitamente):**
- *Osservato:* ordini a tempo uguale 0.698/0.761 (bulk) e 0.367/0.718 (surf) FE RMS;
  FE RMS monotono decrescente; simmetria 180° esatta; ~72% dell'energia di errore nodale
  concentrata nella zona angolare A (§1.4); decadimento di potenza verso l'angolo
  (indice interno ≈ −2, regione-dipendente); floor solvere 10⁻⁸ (§1.5).
- *Ipotesi interpretativa (non dimostrata):* singolarità angolare della soluzione
  accoppiata che limita l'ordine L². Coerente con tutte le osservazioni sopra, ma il fit
  h²/h che avrebbe potuto confermarla **non è identificabile** coi dati disponibili
  (servirebbero 3 mesh a un unico tempo — gap che T1-A colmerebbe, §3). Il fit h²+h
  *non supporta attualmente* la decomposizione.

### 1.4 Errori per regione (solo dati archiviati; metriche nodali non pesate — dichiarato)

Zone: A = r<0.5 dagli angoli (0,1)/(2,0); B = r<0.5 dagli angoli (0,0)/(2,1); resto.

| mesh (step, t) | zona A: RMS / energia | zona B: RMS / energia | resto: RMS / energia |
|---|---|---|---|
| n8×4 (80, 0.02) | 3.97×10⁻³ / **0.708** | 1.76×10⁻³ / 0.139 | 9.71×10⁻⁴ / 0.153 |
| n16×8 (80, 0.02) | 2.26×10⁻³ / **0.725** | 1.03×10⁻³ / 0.151 | 5.30×10⁻⁴ / 0.124 |
| n32×16 (39, 0.00975) | 5.21×10⁻⁴ / **0.722** | 2.63×10⁻⁴ / 0.184 | 1.08×10⁻⁴ / 0.094 |

- L'energia di errore in A è **stabile (~72%)** su tutte e tre le mesh: fatto osservato.
- A tempo uguale (n8→n16 @0.02, nodal RMS): ordini 0.81 (A), 0.77 (B), 0.87 (resto) —
  l'ordine ridotto **non è esclusivo della zona angolare** nemmeno nell'interno; gli angoli
  contribuiscono per la quota di energia, non per un ordine locale diverso.
- Anche la zona B (paio "inattivo") porta 15–18% di energia: coerente con la simmetria
  180°, che accoppia B con se stessa.
- n16→n32 non confrontabile a tempo uguale (0.02 vs 0.00975): **dato mancante**.

Dettagli: `region_analysis_audit.csv`.

### 1.5 Floor OOFEM–DOP853 — **CONFERMATO a mesh, tempo e convenzioni identiche**

| mesh/step (t) | bulk FE RMS (ricalcolato) |
|---|---|
| n16×8 step 39 (0.00975) | 2.0995×10⁻⁸ |
| n16×8 step 80 (0.02) | 2.4441×10⁻⁸ |
| n32×16 step 39 (0.00975) | 2.5972×10⁻⁸ |

Identico ai valori salvati (diff 0.0). Verificato: stessa mesh (array_equal, diff 0,
stesso ordinamento), stesso tempo (Δt ≤ 3.5×10⁻¹⁸), stessa convenzione di campo
(`Temperature`), stessa norma FE (stesso `Mb`/`Ms` dalla stessa via di assemblaggio),
stesse AREA/PERIM. Il riferimento della floor n16×8 è `tolcheck_n16x8_rtol2e-12.npz`
(DOP853 rtol 2×10⁻¹², stessa mesh). La floor 10⁻⁸ è quindi 6–7 ordini sotto l'errore
spaziale (≥2×10⁻⁴): l'errore temporale non contribuisce alle metriche spaziali.

---

## 2. Limiti — dati dichiaratamente mancanti

1. **n32×16 a t=0.02 (step 40–80) manca** (corruzione archivio): nessun ordine a 3 mesh
   su tempo unico; nessun fit h²+h a tempo unico (non identificabile); nessuna analisi di
   zona su n32 a T; nessuna floor n32 a step 80.
2. **Nessun run OOFEM su n64×32** (per design, nessuna nuova simulazione autorizzata):
   il proxy è solo DOP853; l'errore OOFEM a n64 è sconosciuto.
3. Dati per-nodo solo per n8@step80, n16@step80, n32@step39.
4. I VTU esportano solo `u` (`Temperature`), non `v`: nessun controllo sul campo velocità.
5. Gli esponenti r_c e gli ordini di zona dipendono da taglio/soglia (riportati come range);
   le metriche di zona sono nodali (non pesate FE).
6. I due intervalli di raffinamento restano a tempi diversi (0.02 / 0.00975): gli ordini
   non sono confrontabili tra intervalli.

---

## 3. Raccomandazione su T1-A (motivata) + costo residuo

**Raccomando T1-A** (riesecuzione OOFEM n32×16 in nuova directory, parametri identici,
archivio corrotto intatto) — è l'unico test proposto che colma un **gap di dato reale**
identificato da questo audit, a costo minimo:

- (i) ripristina la tabella a 3 mesh su tempo unico (T=0.02) → ordini puliti per tutte le
  4 metriche su un solo tempo;
- (ii) rende **identificabile** il fit h²+h a tempo unico — l'audit ha mostrato che la
  versione a tempi misti è instabile (a<0, residuo ~70%): è l'unico modo residuo di
  testare l'ipotesi di termine O(h)/singolarità angolare sul dataset corrente;
- (iii) abilita l'analisi di zona n16→n32 a tempo uguale (ordini per zona a T);
- (iv) aggiunge la floor n32 a step 80.

Non cambia i verdetti correnti (gli ordini già misurati a t=0.00975: 0.761/0.718 FE RMS)
ma rimuove il caveat dual-time e affina la diagnosi.

**Stima del costo residuo (documentata):** il run OOFEM originale n16×8 (306 DOF) ha
consumato **2 s di real time** (suo `run.log`: "Real time consumed: 000h:00m:02s").
n32×16 (1122 DOF, ×3.7) → stimati **8–30 s di real time**, 1 processo, RAM di ordine di
decine di MB (binario OOFEM già compilato, nessuna modifica al codice, nuova directory
`coupled_independent_space_n32x16_dt0p00025_rerun_t1`, nome unico, parametri documentati).
Post-processing: ~3 s (stage `compare` con riferimenti in cache). **Totale ≈ 1 min di CPU,
RAM trascurabile** — ben dentro il budget di sicurezza.

T1-B (proxy n128×64, ~30–60 min) e T1-C (bulk-only) restano non autorizzati e non avviati.
Nota di pertinenza: dopo T1-A, T1-B conserverebbe utilità per misurare l'ordine 32→64
della famiglia FE (test di recupero verso O(h²) a h più fini).

Se T1-A non viene autorizzato: i verdetti correnti restano validi con il limite dual-time
documentato; il fit h²+h va trattato come recesso (§4).

---

## 4. Erratum a `spatial_verification_report.md`

§6, punto D5(iii) e paragrafo "Interpretazione": il fit
`err(h) ≈ 9.48×10⁻³·h² + 2.83×10⁻³·h (residuo 9%)` era calcolato sui valori
pre-correzione (convenzione temporale con offset di un passo) e **non è riproducibile sui
dati corretti a tempo uguale** (a<0, residuo ~72%, §1.3). Si recede da esso come evidenza.
Le osservazioni robuste restano gli ordini a tempo uguale, la monotonia FE RMS, la
simmetria 180° esatta e la concentrazione energetica angolare (~72%, §1.4). L'ipotesi di
singolarità angolare resta un'interpretazione coerente con l'evidenza, ma non è più
sostenuta dal fit h²/h.
