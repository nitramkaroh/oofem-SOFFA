# T1 — Verifica della convergenza spaziale, sistema Allen–Cahn iperbolico bulk–superficie

**Data:** 2026-10-27
**Stato:** COMPLETATO — **criterio preliminare di ordine NON soddisfatto**; diagnosi completata;
prossimi test minimi proposti (richiedono autorizzazione esplicita).
**Verdetto complessivo:** vedi §7.

Questo rapporto è il deliverable finale di T1: "completare la verifica della convergenza
spaziale del sistema Allen–Cahn iperbolico bulk-superficie, utilizzando esclusivamente i
risultati OOFEM già archiviati e il riferimento indipendente DOP853 corretto".

---

## 1. Setup del test

**Sistema** (invariato rispetto ai test precedenti, nessuna modifica a matematica, IC, BC, parametri):

- Bulk: `ρ u_tt − k Δu = m(u − u³)` in Ω = [0,2]×[0,1], iperbolico, Newmark β=1/4, γ=1/2.
- Superficie: `c_Γ ψ_t − k Δ_Γ ψ + k ∂ₙu = β m(ψ − ψ³)` su Γ = ∂Ω, primo ordine, senza inerzia superficiale.
- Accoppiamento: DOF traccia condivisi (nbc = 0, Neumann naturale), `ψ = u|_Γ`.
- Costanti: ρ = k = m = c_Γ = k_s = 1, η = 0.
- IC: `u₀ = 0.2 + 0.1 cos(πx/2) cos(πy)`, `v₀ = 0` (simmetria rotazionale 180°: u(2−x, 1−y) = u(x,y)).
- `T = 0.02`, `dt = 2.5×10⁻⁴`, 80 step. OOFEM: rtol 1e-8, 80 step (archiviati, non rieseguiti).

**Mesh (annidate 8 | 16 | 32 | 64):**

| mesh | nodi | DOF (u,v) | h |
|---|---|---|---|
| n8×4 | 45 | 90 | 0.25 |
| n16×8 | 153 | 306 | 0.125 |
| n32×16 | 561 | 1122 | 0.0625 |
| n64×32 (proxy) | 2145 | 4290 | 0.03125 |

**Riferimento indipendente (DOP853, assemblaggio con fattore di quadratura corretto):**
la traiettoria fine n64×32 è usata come *proxy* (quasi) indipendente dalla mesh: l'errore di
OOFEM(h) è misurato sui nodi di h (sottoinsiuto esatto dei nodi n64, vedi §3.2).

**Convenzione di interpolazione:** le griglie sono annidate — ogni nodo della mesh grossolana
è un nodo *esatto* della mesh fine (distanza minima verificata = 0.00×10⁰, atol 1e-13,
`assembly_and_comparability_checks.csv`). La restrizione è quindi per **coincidenza nodale
esatta**, senza interpolazione con funzioni forma. Nessuna interpolazione di forma è
necessaria né usata; gli array nodali di mesh diverse non vengono mai confrontati senza
questa restrizione.

**Norme** (stesse della verifica di correzione del Jacobiano, con Ms corretto):

- bulk FE RMS = `sqrt(eᵀ M_b e / AREA)`, AREA = 2.0;
- superficie FE RMS = `sqrt(eᵀ M_s e / PERIM)`, PERIM = 6.0;
- bulk max = max |e| sui nodi interni; surf max = max |e| sui nodi di confine;
- `M_b`, `M_s` dall'assemblaggio indipendente corretto (verificati in §2.1).

---

## 2. Deliverable 1 e 2 — Script, modifiche, parametri, costo

### 2.1 Script utilizzati

| file | ruolo | stato |
|---|---|---|
| `check_spatial_convergence.py` (qui) | script principale T1 (stages: `check`, `tolcheck`, `control32`, `proxy64`, `tolcheck64`, `compare`) | **nuovo** |
| `../run_coupled_independent.py` | mesh, costanti, `rms_fe`, IC | **solo lettura** (invariato) |
| `../coupled_independent_jacobian_fix/check_jacobian_fix.py` | `assemble_fixed` (fattori di quadratura corretti), `residual_reaction_fixed`, `energy_fixed` | **solo lettura** (invariato) |

Nessuna modifica al codice OOFEM, ai benchmark o al riferimento matematico: costanti, IC, BC,
assemblaggio e fattori di quadratura sono identici alla verifica del Jacobiano già accettata.

**Modifiche apportate nello script T1 (solo diagnostiche/affidabilità, nessuna modifica matematica):**
1. `dop853`/`cached_run`: restituzione di `nfev`, log del picco RSS (`ru_maxrss`),
   `assert` di finitezza (NaN/Inf → stop) e cap hard `wall < 300 s` (auto-stop, no retry).
2. `stage_proxy64`: pre-probe su 9 punti temporali con stima del costo completo e
   `assert est < 270 s` prima del lancio (compliance FASE 2).
3. `stage_tolcheck64`: verifica di accuratezza **parziale** (9 punti, t ≤ 0.002) per evitare
   un secondo calcolo completo costoso (compliance FASE 4).
4. `stage_compare`: ridisegno *dual-time* (vedi §5.2) per aggirare la corruzione
   dell'archivio n32×16; mappature nodali separate `idx64`/`idx32`; riga `NOT AVAILABLE`
   documentata per il campo n32×16 step 80; floor del solvere su mesh e tempo uguali;
   CSV per-nodo e PNG a due pannelli.
5. `read_step_oofem`: rilevamento e errore documentato per VTU a 0 byte.
6. `validated_mesh`: rigenerazione mesh + cross-validazione contro `case.in`
   (l'`mesh.npz` di n32×16 è 0 byte — archivio corrotto, lasciato intatto).
7. **Correzione off-by-one nella convenzione temporale** (documentata in §5.4): la prima
   stesura del `compare` usava `U[step−1]` (proxy a t = (step−1)·dt) confrontando con il
   campo OOFEM a t = step·dt; corretto a `U[step]` (stesso tempo, convenzione della
   verifica Jacobiano). Solo lo script di confronto è stato toccato; nessun dato,
   assemblaggio o parametro è stato alterato.

### 2.2 Parametri del riferimento DOP853

- Metodo: DOP853 (scipy 1.11.4), `rtol = 2×10⁻¹²`, `atol = 2×10⁻¹⁴` (nessuna tolleranza
  allentata per terminare prima).
- Intervallo: t ∈ [0, 0.02], 81 punti (`t_eval`), nessun salvataggio di storia temporale
  oltre i punti richiesti; stato ODE = 2×(nodi) float64.
- Rilevato: scipy clamp di `rtol = 2×10⁻¹⁴` a 2.2×10⁻¹⁴ (epsilon macchina) nella run
  di controllo più stretta — documentato, non modificato.

### 2.3 Costo effettivo (da `run_log.csv`; 1 core, `taskset -c 0`, cap `ulimit -v` 3 GB)

| stage | wall (s) | note |
|---|---|---|
| `check` | 13.77 | 34 controlli di comparabilità/assemble, tutti passati |
| `tolcheck` n16×8 | 2.72 + 4.47 | due run DOP853 (rtol 2e-12 e 2e-14), 81 punti |
| `tolcheck` (confronto) | 7.18 | max\|ΔU\| = 2.658×10⁻¹⁴ |
| `control32` (n32×16) | 30.41 | traiettoria di controllo completa |
| `proxy64` probe | 28.98 | stima run completo: 261 s < budget 300 s |
| `proxy64` **completo** | **167.9** | n64×32, nfev = 1088 (~90 step × 12 stadi), RSS picco 77.1 MB |
| `tolcheck64` parziale | 43.48 | 9 punti (t ≤ 0.002), nfev = 245, max\|ΔU\| = 2.109×10⁻¹⁵ |
| `compare` | 2.79 + 2.67 | post-processing (prima stesura + riesecuzione dopo correzione off-by-one) |

**Totale ≈ 5.2 min di CPU su 1 core; picco di RAM 77.1 MB** (budget: +2 GB; margine ~25×).
Pilot (FASE 2): RHS n64×32 = 151.52 ms (mediana, 10 valutazioni a stato fisso e limitato);
nessuna matrice densa; nessun assemblaggio all'interno del RHS (operatori sparsi assemblati
una volta, riutilizzati; `splu(M_b)` una volta).

---

## 3. Deliverable 3 — Tabella degli errori (bulk e superficie separati)

### 3.1 Errore di OOFEM rispetto al proxy n64×32, a tempo uguale (t = step·dt)

| mesh | step | t | bulk FE RMS | bulk max | surf FE RMS | surf max |
|---|---|---|---|---|---|---|
| n8×4 | 80 | 0.02 | 1.3286×10⁻³ | 8.345×10⁻⁴ | 2.5727×10⁻³ | 6.069×10⁻³ |
| n16×8 | 80 | 0.02 | 8.188×10⁻⁴ | 1.880×10⁻³ | 1.9948×10⁻³ | 4.744×10⁻³ |
| n16×8 | 39 | 0.00975 | 3.557×10⁻⁴ | 4.845×10⁻⁴ | 9.287×10⁻⁴ | 2.466×10⁻³ |
| n32×16 | 39 | 0.00975 | 2.099×10⁻⁴ | 1.281×10⁻³ | 5.645×10⁻⁴ | 1.463×10⁻³ |
| n32×16 | 80 | 0.02 | **NON DISPONIBILE** (archivio corrotto, §5.2) | | | |

Convergenza FE RMS (monotona): bulk 1.33×10⁻³ → 8.2×10⁻⁴ → 2.1×10⁻⁴;
superficie 2.57×10⁻³ → 2.0×10⁻³ → 5.6×10⁻⁴.

### 3.2 Floor del solvere: OOFEM vs DOP853 su *stessa mesh e stesso tempo*

| mesh | step | t | bulk FE RMS | bulk max | surf FE RMS | surf max |
|---|---|---|---|---|---|---|
| n16×8 | 39 | 0.00975 | 2.099×10⁻⁸ | 4.8×10⁻⁸ | 3.567×10⁻⁸ | 1.24×10⁻⁷ |
| n16×8 | 80 | 0.02 | 2.444×10⁻⁸ | 6.7×10⁻⁸ | 5.151×10⁻⁸ | 2.09×10⁻⁷ |
| n32×16 | 39 | 0.00975 | 2.597×10⁻⁸ | 3.7×10⁻⁷ | — | — |

**OOFEM riproduce il riferimento indipendente a ~10⁻⁸–10⁻⁷**, cioè 6–7 ordini di
grandezza sotto gli errori spaziali misurati: l'implementazione (assemblaggio, accoppiamento,
Newmark, stepping temporale) è corretta entro il rumore numerico dei due codici.

### 3.3 Controlli di accuratezza del riferimento (FASE 4)

| mesh | intervallo | rtol base → stretta | max\|ΔU\| |
|---|---|---|---|
| n16×8 | [0, 0.02] completo | 2×10⁻¹² → 2×10⁻¹⁴ | 2.658×10⁻¹⁴ |
| n64×32 | [0, 0.002] (9 punti) | 2×10⁻¹² → 2×10⁻¹⁴ | 2.109×10⁻¹⁵ |

Errore temporale del riferimento ≤ ~3×10⁻¹⁴: 8–9 ordini sotto il segnale spaziale
(≥ 2×10⁻⁴) → trascurabile, dimostrato senza secondi run completi costosi (check parziale
sulla mesh fine + intervallo completo sulla mesh n16×8).

### 3.4 Livelli della famiglia FE (sanità del proxy, FE RMS bulk)

| quantità @ t | valore |
|---|---|
| OOFEM n16\|n8 − OOFEM n8, t=0.02 | 5.045×10⁻⁴ (bulk), 7.411×10⁻⁴ (surf) |
| PROXY n32\|n8 − PROXY n64\|n8, t=0.02 | 4.191×10⁻⁴ (bulk), 8.746×10⁻⁴ (surf) |
| OOFEM n32\|n16 − OOFEM n16, t=0.00975 | 1.861×10⁻⁴ (bulk), 3.741×10⁻⁴ (surf) |
| PROXY n32\|n16 − PROXY n64\|n16, t=0.00975 | 1.799×10⁻⁴ (bulk), 5.358×10⁻⁴ (surf) |

Le differenze di livello sono uguali (a fattore 1.1–1.5) nelle due implementazioni
indipendenti: il comportamento sub-quadratico è una proprietà della **famiglia FE**,
non uno scostamento OOFEM-vs-riferimento (diagnosi D1).

---

## 4. Deliverable 4 — Ordini di convergenza

`order = log(err_coarse/err_fine, 2)` (raffinamento ×2 su entrambe le direzioni).
**Criterio preliminare dichiarato: ordine ∈ [1.5, 2.5].**

| intervallo | t | metrica | err grosso | err fine | **ordine** |
|---|---|---|---|---|---|
| 8×4 → 16×8 | 0.02 | bulk FE RMS | 1.329×10⁻³ | 8.188×10⁻⁴ | **0.698** |
| 8×4 → 16×8 | 0.02 | bulk max | 8.345×10⁻⁴ | 1.880×10⁻³ | **−1.172** (non monotono) |
| 8×4 → 16×8 | 0.02 | surf FE RMS | 2.573×10⁻³ | 1.995×10⁻³ | **0.367** |
| 8×4 → 16×8 | 0.02 | surf max | 6.069×10⁻³ | 4.744×10⁻³ | **0.355** |
| 16×8 → 32×16 | 0.00975 | bulk FE RMS | 3.557×10⁻⁴ | 2.099×10⁻⁴ | **0.761** |
| 16×8 → 32×16 | 0.00975 | bulk max | 4.845×10⁻⁴ | 1.281×10⁻³ | **−1.402** (non monotono) |
| 16×8 → 32×16 | 0.00975 | surf FE RMS | 9.287×10⁻⁴ | 5.645×10⁻⁴ | **0.718** |
| 16×8 → 32×16 | 0.00975 | surf max | 2.466×10⁻³ | 1.463×10⁻³ | **0.753** |

**Riscontro al criterio: NON soddisfatto.** Ordini FE RMS 0.37–0.76 ( fuori [1.5, 2.5]);
bulk max non monotono (cresce col raffinamento). Gli errori FE RMS sono monotoni
decrescenti su entrambi gli intervalli. I due intervalli sono a tempi diversi (0.02 e
0.00975) — conseguenza della corruzione dell'archivio n32×16 (§5.2) — e non sono
direttamente confrontabili tra loro.

---

## 5. Note operative

### 5.1 Controlli preliminari (tutti passati, 34/34 in `assembly_and_comparability_checks.csv`)
- Mesh rigenerata == coordinate di `case.in` (max|diff| = 0) per n8/n16/n32;
  nodi di confine == insieme `line1ht`.
- Annidamento: nodi n8, n16, n32 ⊂ nodi n64 (distanza minima 0.00×10⁰).
- Elementi di linea coprono l'intero confine (24/48/96/192 archi).
- `Σ M_b 1 = AREA = 2.0`, `Σ M_s 1 = PERIM = 6.0` (fattore di quadratura corretto L/2),
  righe di K_b, K_s a somma nulla, `M_s` di singolo arco = (L/6)[[2,1],[1,2]],
  `K_s` = (1/L)[[1,−1],[−1,1]], partition-of-unity delle reazioni, energia superficiale
  F(0.2)·PERIM, parametri `case.in` == costanti del riferimento, IC identica a ogni nodo.

### 5.2 Corruzione dell'archivio n32×16 (dati, non codice)
Nell'archivio `../coupled_independent_space_n32x16_dt0p00025/` i file
`case.out.m0.40.vtu` … `case.out.m0.80.vtu` (41 file), `mesh.npz` e `run.log` sono **0 byte**;
gli step 1–39 (t ≤ 0.00975) sono integri e finiti. L'archivio è stato **lasciato intatto**
(vincolo 4). Conseguenza: il confronto n16→32 a T=0.02 non è eseguibile; è stato adottato il
disegno **dual-time**: intervallo 8→16 a t=0.02 (entrambi integri) e 16→32 a t=0.00975
(step 39, tutti interi, finiti). La riga mancante è riportata esplicitamente nella tabella.

### 5.3 Adesione alla protocollo FASE 0–5
- **FASE 0** (ispezione statica, nessuna integrazione ODE): fattori di quadratura corretti
  verificati (§5.1); inventario risorse: nessuna matrice densa, nessun assemblaggio nel RHS,
  stato ODE 2×4290 float64, costo singolo RHS n64 misurato 151.52 ms → stima run completo
  121–261 s ≪ budget; stima memoria dal vettore e dalle strutture sparse (non solo dal
  vettore soluzione): picco atteso < 400 MB. Nessun errore concreto trovato → nessun codice
  modificato in questa fase.
- **FASE 1** (strategia a memoria controllata): operatori sparsi, nessuna moltiplicazione
  sparso→denso, assemblaggio una volta, nessun ri-build di invarianti nel RHS, nessuna
  storia temporale, nessun salvataggio di stato per RHS, 1 processo / 1 core (`taskset -c 0`),
  cap virtual memory `ulimit -v 3145728` (3 GB, margine ~25× sul picco reale 77.1 MB).
  I limiti sono imposti e verificati (picco misurato con `ru_maxrss`); il timeout (60 s
  pilot / 300 s completo) non è l'unica protezione: esistono assert di finitezza, cap di
  wall nel codice e il cap di memoria del sistema operativo.
- **FASE 2** (pilot obbligatorio): eseguito prima del run completo (10 valutazioni RHS a
  stato fisso e fisicamente limitato, nessuna integrazione); stima 261 s < 300 s → green.
- **FASE 3** (riferimento fine): eseguito dentro i limiti (167.9 s; nfev 1088; RSS picco
  77.1 MB; nessun NaN/Inf; nessuna crescita anomala). Tolleranze non allentate.
- **FASE 4** (accuratezza del riferimento): dimostrata con check mirati a budget limitato
  (§3.3); nessun secondo run completo lanciato.
- **FASE 5** (confronti senza riesecuzione OOFEM): solo archivi in sola lettura;
  compatibilità verificata prima delle metriche (§5.1); interpolazione per coincidenza
  nodale esatta (§1); norme e convenzioni documentate (§1); bulk e superficie separati.

### 5.4 Correzione off-by-one nella convenzione temporale del confronto (trasparenza)
La prima stesura di `compare` confrontava il campo OOFEM allo step s (t = s·dt, verificato
sul pvd: m0.80 → t = 0.02) con il proxy all'indice `s−1` (t = (s−1)·dt). Il riscontro: la
"floor temporale" misurata su n32 (1.006×10⁻⁵ FE RMS) coincideva a 4 cifre con il drift
del proxy su un singolo step (1.007×10⁻⁵) — cioè la floor era il drift, non la
differenza OOFEM–riferimento. Correzione: indice `U[s]` (stesso tempo, convenzione della
verifica Jacobiano già accettata); `compare` rieseguito (2.67 s, post-processing).
Effetti: (a) tabelle aggiornate con confronto a tempo uguale; (b) la floor del solvere è
stata raffinata da ~10⁻⁵ (che era il drift) a **~10⁻⁸ (valore vero)**; (c) gli ordini
cambiano per < 0.01 e le conclusioni restano invariate. Il bug era nello script di
confronto, non nei dati né nel riferimento.

---

## 6. Diagnosi (il criterio non è soddisfatto → §3 del mandato: diagnosticare)

Tutte le diagnosi sono post-processing su dati già archiviati/calcolati (nessuna nuova
simulazione).

- **D1 — Proprietà della famiglia FE.** Le differenze di livello OOFEM−OOFEM e
  proxy−proxy coincidono a fattore 1.1–1.5 su entrambi gli intervalli e in entrambe le
  variabili (§3.4). L'ordine sub-quadratico è quindi una proprietà della discretizzazione
  FE del sistema accoppiato, comune a due implementazioni indipendenti — non un difetto
  dell'implementazione OOFEM.
- **D2 — L'errore temporale non è dominante.** Floor del solvere 2.1–5.2×10⁻⁸ (n16×8) e
  2.6×10⁻⁸ (n32×16); errore temporale del riferimento ≤ 3×10⁻¹⁴. La drift su un step
  (8.6×10⁻⁶…1.5×10⁻⁵ FE RMS) è ≥ 15× sotto gli errori spaziali misurati.
- **D3 — Accumulo, non instabilità.** L'errore n32×16 rispetto al proxy cresce
  regolarmente col tempo: 2.8×10⁻⁶ (t=0.00125) → 1.5×10⁻⁵ → 6.2×10⁻⁵ → 2.0×10⁻⁴
  (t=0.00975) (bulk FE RMS); nessun blow-up.
- **D4 — Accenno di perdita di regolarità agli angoli.** Sul proxy n64×32 a t=0.00975,
  lungo il bordo sinistro vicino all'angolo (0,1): max |d²u/dy²| = 0.854 contro 0.365 a
  metà bordo (2.3×), con profilo regolare (nessun kink visibile a n64): la curvatura è
  amplificata in prossimità dell'angolo.
- **D5 — Localizzazione dell'errore agli angoli.** (i) i massimi errori interni su n32×16
  sono ai nodi (0.0625, 0.9375) e (1.9375, 0.0625) — adiacenti agli angoli (0,1) e (2,0),
  con simmetria 180° esatta e segno negativo (massima discrepanza 180°: 2.2×10⁻¹⁶,
  precisione macchina); (ii) per i nodi interni con 0 < r_c < 0.5 (distanza dall'angolo
  attivo) |err| decade con legge di potenza (indice ≈ −2, range −1.9…−2.8 a seconda del
  taglio di regione; vedi audit): l'errore è fortemente concentrato in prossimità degli
  angoli, e ~72% dell'energia di errore nodale sta nella zona angolare (stabile su tutte
  e tre le mesh); (iii) **ERRATUM**: il fit `err(h) ≈ 9.48×10⁻³·h² + 2.83×10⁻³·h
  (residuo 9%)` citato in questa stesura era calcolato sui valori pre-correzione (offset
  temporale di un passo) e **non è riproducibile sui dati corretti a tempo uguale**
  (coefficiente di h² negativo, residuo ~72%): si recede da esso come evidenza. Il fit
  h²+h a 3 punti non è identificabile a tempo unico coi dati disponibili (manca n32 a
  T=0.02). Vedi `static_audit_report.md` §1.3 e §4.

**Interpretazione:** il sistema accoppiato bulk–superficie ammette, in questo caso test,
una singolarità agli angoli del dominio (il dato di bordo effettivo di tipo Robin porta,
al cambio di bordo, un kink nel derivato normale; il bulk eredita una componente
singolare). Con elementi Q1/P1 ciò limiterebbe l'ordine di convergenza L² globale a < 2
(osservato 0.37–0.76) e produrrebbe la non-monotonicità della norma max vicino agli
angoli. La convergenza esiste (FE RMS monotono) ma è più lenta del criterio
preliminare. *Nota (audit 2026-10-27):* questa interpretazione resta un'ipotesi
coerente con l'evidenza osservata (ordini a tempo uguale, simmetria 180° esatta,
concentrazione energetica angolare ~72%), ma il fit h²/h che avrebbe potuto
confermarla non è identificabile coi dati disponibili (erratum D5(iii) e
`static_audit_report.md` §1.3); il test minimo per confermarla/infamarla resta T1-A
(3 mesh a tempo unico). La validazione
dell'implementazione OOFEM non è in discussione: il solvere riproduce il riferimento
indipendente a ~10⁻⁸ (D2, §3.2).

---

## 7. Deliverable 6 — Verdetti separati

| voce | verdetto | motivazione |
|---|---|---|
| Convergenza spaziale **bulk** | **NON VALIDATO** (rispetto al criterio preliminare) | Ordini 0.698 / 0.761 ∉ [1.5, 2.5]; bulk max non monotono (−1.172 / −1.402). FE RMS monotono decrescente. |
| Convergenza spaziale **superficie** | **NON VALIDATO** (rispetto al criterio preliminare) | Ordini 0.367 / 0.718 (FE RMS), 0.355 / 0.753 (max) ∉ [1.5, 2.5]; FE RMS monotono decrescente. |
| **Validazione complessiva dell'accoppiamento** | **SUPPORTATA PARZIALMENTE** | L'implementazione OOFEM (assemblaggio, accoppiamento bulk–superficie, Newmark, stepping) è corretta: riproduce il riferimento indipendente DOP853 a ~10⁻⁸–10⁻⁷ a mesh e tempo uguali, 6–7 ordini sotto l'errore spaziale. La deviazione dall'ordine atteso è diagnosticata come proprietà matematica della discretizzazione FE accoppiata (regolarità angolare, D1/D4/D5), non come difetto OOFEM. L'ordine di convergenza spaziale quadratico, comunque, **non** è dimostrato su questo caso test. |

Nessun "perfettamente validato": i risultati sono quelli della tabella, con le
incertezze e i limiti dichiarati (proxy come surrogato di riferimento; due intervalli a
tempi diversi per corruzione d'archivio; ordine max non monotono).

---

## 8. File prodotti (tutti in questa directory, riproducibili)

| file | contenuto |
|---|---|
| `check_spatial_convergence.py` | script T1 (stages documentati) |
| `spatial_verification_report.md` | questo rapporto |
| `spatial_convergence_final.csv` | tabella errori + proxy level diff + floor (riga n32 step 80 documentata come non disponibile) |
| `convergence_orders.csv` | 8 ordini (4 metriche × 2 intervalli) |
| `per_node_n8x4_step80_vs_proxy_n64x32.csv`, `per_node_n16x8_step80_…`, `per_node_n32x16_step39_…` | errori per nodo (localizzazione) |
| `spatial_convergence.png` | grafico a due pannelli (t=0.02: n8,n16; t=0.00975: n16,n32; 4 metriche) |
| `assembly_and_comparability_checks.csv` | 34 controlli preliminari (tutti True) |
| `tolerance_check_n16x8.csv`, `tolerance_check_n64x32_partial.csv` | accuratezza del riferimento |
| `reference_n64x32_dt0p00025.npz` (nfev 1088), `reference_n32x16_dt0p00025.npz`, `tolcheck_*.npz` | traiettorie di riferimento (cache) |
| `run_log.csv` | cronologia costi di ogni stage |

Archivi OOFEM utilizzati (sola lettura, non modificati):
`../coupled_independent_space_n{8x4,16x8,32x16}_dt0p00025/` (n32×16 corrotto da step 40,
lasciato intatto).

---

## 9. Test minimi successivi proposti (nessuno eseguito — richiedono autorizzazione esplicita)

- **T1-A — Recuperare T=0.02 su n32×16 (richiede riesecuzione OOFEM).** Rieseguire
  *solo* il caso OOFEM n32×16 in una **nuova directory**
  (`coupled_independent_space_n32x16_dt0p00025_rerun_t1`), parametri identici; l'archivio
  corrotto resta intatto. Costo stimato: ~10–30 s di CPU (dalla durata del run originale)
  + pochi secondi di post-processing; RAM trascurabile. Consegna: i due intervalli di
  raffinamento a un unico tempo (T=0.02) e una tabella di convergenza pulita. Non cambia
  l'ordine atteso (già misurabile a t=0.00975) ma rimuove il disegno dual-time.
- **T1-B — Secondo intervallo della famiglia FE (nessun OOFEM; richiede autorizzazione di costo).**
  Proxy DOP853 n128×64 (8385 nodi, 16770 DOF): ordine 32→64 della famiglia FE. Se
  l'ordine recupera verso ~2, l'effetto è di h-finito; se resta ~0.8, la limitazione di
  regolarità angolare è confermata. Costo stimato: ~30–60 min su 1 core, ~300 MB RAM
  (scala da 151.5 ms/RHS su n64 a ~0.5–0.7 s/RHS su n128, ~150–250 step).
- **T1-C — Controllo bulk-only (richiede riesecuzione OOFEM + variante di riferimento).**
  Problema bulk senza superficie né accoppiamento (Neumann naturale) su n8/n16/n32,
  OOFEM + riferimento DOP853: separa la discretizzazione bulk pura (attesa O(h²))
  dall'effetto dell'accoppiamento angolare. Costo stimato: 3 run OOFEM ~10–60 s ciascuno,
  3 run DOP853 ~5–30 s ciascuno, + implementazione della variante di riferimento
  (nessuna modifica al codice OOFEM esistente: nuovo caso .in + nuovo riferimento).

**Nessuna di queste azioni è stata intrapresa.** In attesa di decisione.
