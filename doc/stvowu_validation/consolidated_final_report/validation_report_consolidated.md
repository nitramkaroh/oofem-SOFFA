# Rapporto consolidato di validazione — modelli Allen–Cahn bulk–superficie in OOFEM

**Data:** 2026-10-27
**Versione:** 1.0 (stato delle evidenze al completamento di T1-B)
**Perimetro:** esclusivamente documentale e analitico. Nessuna nuova simulazione, integrazione,
raffinamento o modifica a codice/risultati archiviati è stata eseguita per questo rapporto.
Tutti i numeri provengono da artefatti archiviati, qui referenziati (sezione H e lista finale).

**Vincolistica osservata:** nessuna esecuzione OOFEM; nessuna integrazione DOP853; nessun
raffinamento di mesh; nessun processo in background/parallelismo/benchmark; nessuna modifica
ai risultati numerici archiviati; nessun sovrascrittura di rapporti, script, CSV, VTU, PVD,
log o cache; nessuna correzione retroattiva dei dati originali; nessuna conversione
sparse→dense; nessuna installazione di dipendenze. I dati non disponibili sono registrati
come mancanti, non ricostruiti per supposizione.

**Lingua e convenzioni:** rapporti di campagna in italiano; norme e definizioni come nei
rapporti di riferimento (sezione E). Il riferimento DOP853 è in ogni caso un **riferimento
numerico per il sistema semidiscreto**, non la soluzione esatta della PDE.

---

## A. Sintesi esecutiva

La campagna copre tre famiglie di modelli implementati nel materiale `allencahnmat` di OOFEM:

1. **Parabolica bulk–superficie** `c u_t − kΔu = m(u−u³)` + `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)`
   (bordo di primo ordine, capacità superficiale, diffusione tangenziale, flusso normale accoppiato);
2. **Iperbolica `wave 1`**: bulk inerziale smorzato `ρu_tt + ρζu_t − kΔu = m(u−u³)`
   (Newmark β=1/4, γ=1/2, full-Newton) + superficie di **primo ordine** (stessa legge del parabolico,
   `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)`);
3. **Caso accoppiato dinamico bulk–superficie** (2D, Ω=[0,2]×[0,1]) — oggetto di T1, T1-A, T1-D, T1-B.

**Stato sintetico delle evidenze (dettaglio in sezione J):**

| Oggetto | Verdetto |
|---|---|
| Correttezza dell'implementazione/assemblaggio | **SUPPORTATA** |
| Riproducibilità | **SUPPORTATA** |
| Accordo col riferimento semidiscreto indipendente (temporale) | **SUPPORTATA** (nel caso testato) |
| Convergenza spaziale del caso accoppiato (esistenza, fino a h=0.03125) | **SUPPORTATA** |
| Ordine asintotico spaziale del caso accoppiato | **NON DIMOSTRATA** |
| Singolarità angolare | **NON DEMOSTRATA** |
| **Validazione complessiva** | **SUPPORTATA PARZIALMENTE** |

Fatti principali:

- Il **bulk iperbolico** è confermato contro la soluzione analitica di Nizovtseva (onda
  viaggiante, ordini ≈2 spaziali e temporali, velocità d'onda a 9.8·10⁻⁵ relativa) e contro il
  kink stazionario tanh (ordini 2.006/1.887): l'equazione bulk, l'avvio Newmark e la
  discretizzazione Q1 sono coerenti in quei casi.
- La **discrepanza originaria** del cross-check accoppiato (errori OOFEM-vs-ODE quasi costanti,
  ordini ~0) è stata **identificata, dimostrata e corretta sul lato del riferimento Python
  indipendente** (fattori di Jacobiano degli elementi linea: Ms ×2, Ks ×½, rs ×2); il codice
  OOFEM non era in causa. Dopo la correzione, la convergenza temporale accoppiata mostra
  ordini 1.972/1.948 (bulk) e 1.975/2.027 (superficie). I documenti che riportano ancora
  "FAIL" per il sistema accoppiato (matrice "finale" e report parabolico del 26-10) sono
  **superati** e lo si dichiara qui esplicitamente.
- La **convergenza spaziale del caso accoppiato** esiste e prosegue fino alla fine della
  famiglia proxy: la differenza di livello famiglia 64→128 (h=0.03125→0.015625) è
  4.87·10⁻⁵ (bulk FE RMS), 8.5× sotto la 32→64 (4.14·10⁻⁴): nessun plateau. L'errore OOFEM
  n32×16 rispetto al proxy n128 (3.82·10⁻⁴) è entro l'8% della differenza di famiglia 32→64:
  a quel livello la componente OOFEM (floor ~10⁻⁸) è trascurabile rispetto alla differenza
  di famiglia. **Il tasso asintotico non è però stabilito**: gli ordini 0.72 / 1.05 / 3.09
  (bulk FE RMS) non sono costanti, dipendono dal tempo, e il valore ~3.1 è un singolo
  rapporto fra due differenze di livello, **non** una dimostrazione di convergenza di
  terzo ordine.
- La **concentrazione dell'errore agli angoli (0,1)/(2,0)** (66–75% dell'energia di errore)
  è una proprietà della famiglia FE per questo problema (compare anche nell'auto-differenza
  del proxy, senza OOFEM), presente dal primo passo; l'ipotesi di limite non regolare
  (singolarità) resta un'interpretazione coerente ma **non dimostrata**: il massimo mobile
  a T scala come r_c^(−0.016) su 4 punti (essenzialmente piatto), non come legge di potenza.
- I **bilanci energetici** sono coerenti con le leggi di dissipazione continue in entrambi i
  rami (monotonia, residui che calano col raffinamento, additività a livello di roundoff),
  ma la dissipazione è ricostruita a posteriori (differenze secanti, quadratura dichiarata,
  velocità non esportata dal solver): la monotonia da sola non costituisce validazione
  completa del solver.
- La **riproducibilità** è bit-identica dove testata (3 run OOFEM ri-eseguiti; T1-A steps 1–39
  identici; assembly fresca == cache in T1-B con diff 0.000e+00; 16/16 metriche ricalcolate
  con diff 0). Restano limiti ambientali (build dirty, raw in `/tmp`, nessun lockfile).

**Attività future raccomandate (non eseguite, in ordine di valore informativo/costo, §J.3):**
1. run OOFEM n64×32 del caso accoppiato (floor a h=0.03125 + quarto punto per-nodo);
2. export della velocità `waveVelocity` nei VTU (richiede approvazione esplicita, modifica
   codice mirata e non numerica) per rendere il bilancio iperbolico una verifica indipendente;
3. test bulk-only per isolare l'effetto dell'accoppiamento sull'ordine ridotto
   (fattibilità input da verificare prima).

---

## B. Modello matematico e implementazione

### B.1 Equazioni implementate (confermate dal codice)

Verifica termine per termine del sorgente (ispezione statica, nessun MMS):

| Componente | Equazione continua | Implementazione OOFEM (verificata) |
|---|---|---|
| Bulk parabolico | `c u_t − kΔu = m(u−u³)` | capacità `Σ N Nᵀ c dV`, rigidezza `Σ Bᵀ d B dV`, sorgente `Σ N m(u−u³) dV`, tangente `Σ N Nᵀ m(3u²−1) dV` (full-Newton) |
| Superficie parabolica | `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)` | `line1ht` isoparametrico (J = L/2), 3 punti di Gauss (esatti per la cubica), β = `surfacebeta`; **in modalità parabolica la linea condivide la `c` del bulk** (nessun `c_Γ` indipendente in input — semantica documentata, non un difetto) |
| Bulk iperbolico | `ρu_tt + ρζu_t − kΔu = m(u−u³)` | Newmark media-accelerazione: LHS `K + (1/(βdt²) + ζγ/(βdt))·ρM_bulk` (solo bulk), RHS con `a + ζv`; Newton completo con tangente aggiornata; accelerazione iniziale coerente `M_bulk a₀ = F − K u₀ − M_Γ v₀` (correzione del 26-10: la zero-inizializzazione dava errore di avvio di 1° ordine) |
| Superficie iperbolica (1° ordine) | `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)` | riga LHS `K_Γ + γ/(βdt)·M_Γ`, RHS `M_Γ v`; **nessuna inerzia superficiale** (`BulkWaveMassAssembler` esclude le linee) |
| Capacità | — | parabolico: `c` su bulk e superficie; iperbolico: `rho` (inerzia bulk), `surfacecapacity` (superficie), `c` non usata |
| Accoppiamento `k∂ₙu` | flusso normale implicito | DOF di tracciato **condivisi** tra quad e line (la `line1ht` usa gli stessi ID nodo del bordo dei quadrilateri): la continuità della traccia è cinematica; il termine entra come flusso naturale nelle righe di bordo del bulk (nbc 0, Neumann) — costruzione identica al riferimento indipendente `K = Kb + Ks` |

Riferimenti di codice: `allencahnmat.C` (materiali), `transportelement.C` (residuo, matrici,
tangente, quadratura), `nltransienttransportproblem.C` (modello transiente, ramo `wave 1`,
Newmark, `applyIC`), `line1_ht.C` / `fei3dlinelin.C` (elementi linea), `quad1_ht.C`
(quadratura), `initialcondition.C` (condizioni iniziali).

### B.2 Integrazione temporale, condizioni iniziali e al contorno

- **Parabolico:** backward Euler (BDF1) totalmente implicito, Newton con tangente aggiornata
  alla prima iterazione del passo. Ordine atteso: 1.
- **Iperbolico:** Newmark media-accelerazione (β=1/4, γ=1/2) su stato (u, v), Newton completo.
  Ordine atteso: 2.
- **IC:** `u₀` da campo, `v₀` da `initialvelocity` (valori per DOF in ordine nodale), `a₀`
  consistente risolto dal residuo semidiscreto (non imposta a zero).
- **BC:** nessun Dirichlet nei benchmark; Neumann omogenea naturale su tutto il contorno
  (tratto come condizione di troncamento, dove rilevante, e quantificato: Nizovtseva §5).
- **Quadratura:** `quad1ht` default 2×2 — **sotto-integra la reazione cubica Q1** (grado 4);
  i benchmark usano `nip 9` (3×3, esatta) ed è verificato che RHS e assemblaggio usino la
  stessa regola. `line1ht` 3 punti (esatto per la cubica). **Caveat: gli input devono
  dichiarare `nip 9`** (vero in tutti i benchmark di validazione).

### B.3 Solo documentato, non confermato dai test

- I due **bachi noti nel codice** (documentati in `validation_audit.md`, non toccati dalla
  campagna): (i) velocità/accelerazione prescritte sui DOF non onorate nel ramo Newmark
  (Dirichlet tempo-dipendente + `wave 1`); (ii) sorgente MMS iperbolica usa `c·u_t` invece
  di `ρζ·u_t`. Nessun benchmark di validazione li esercita (nessun test usa Dirichlet
  tempo-dipendente in wave mode; i casi MMS sono esclusi per vincolo).
- **Inerzia superficiale di 2° ordine non implementata**: il modello iperbolico è
  bulk-inerziale + superficie-parabolica. I preprint che descrivono superficie inerziale
  (Hipp–Kovács 2018/53, Hochbruck–Leibold 2019/26, 2020/20) **non** descrivono il modello
  implementato: nessun loro benchmark è usato come validazione del modello accoppiato.
- Coefficienti indipendenti `c_Γ`, `k_Γ`, `m_Γ` non isolati uno per volta in nessun test
  (tutti i benchmark usano valori unitari o Okumura `c=10` agli estremi).

---

## C. Matrice di validazione

Legenda classi di riferimento: **A** soluzione analitica esatta non fabbricata; **A-sem**
soluzione esatta del sistema FEM semidiscreto; **B** dati pubblicati (digitalizzati);
**C** solutore indipendente; **D** raffinamento contro riferimento numerico più fine;
**E** confronto qualitativo; **F** bilancio energetico.
Il **criterio** è "pre-registrato" solo se definito prima dei risultati in
`audit_iniziale_2026-10-27.md` §7–8 (ordine in [1.5, 2.5] per schemi 2° ordine,
[0.75, 1.25] per 1° ordine; per T1: monotonia, ordini in [1.5, 2.5] su entrambi gli
intervalli, sanity del proxy). Gli altri test non hanno criteri pre-registrati e sono
marcati **descrittivo/esplorativo**. Versione tabellare: `validation_matrix.csv`.

| # | Test | Equazione / componente | Riferimento (classe) | Metrica | Risultato | Criterio | Verdetto | Limiti principali |
|---|---|---|---|---|---|---|---|---|
| C1 | Onda viaggiante di Nizovtseva | bulk iperbolico | formula Eq. (29), verificata come soluzione esatta (residuo 1.1e-16) (A) | RMS/max nodale, velocità d'onda | N=80/160/320: RMS 2.9183e-5 / 7.3231e-6 / 1.8621e-6; ordini 1.9946 / 1.9755; temporali (N=320) 2.0086 / 1.9165; velocità 1.539119 vs 1.538968 (9.8e-5 rel.); IC mismatch 3.61e-16 | pre-registrato [1.5, 2.5] | **SUPPORTATA** (nel caso testato) | solo bulk (niente superficie); troncamento Neumann quantificato ([−8,8] vs [−12,12]: max 1.7e-4 vs 1.3e-5); RMS nodale, non norma FE |
| C2 | Kink stazionario tanh(x/√2) | bulk iperbolico, equilibrio | profilo esatto stazionario (A) | RMS nodale su [−8,8], T=0.2 | N=40/80/160: 3.0955e-5 / 7.7074e-6 / 2.0844e-6; ordini 2.0059 / 1.8866 | pre-registrato [1.5, 2.5] | **SUPPORTATA** (nel caso testato) | Neumann ai bordi solo approssimata (derivata esatta ≠ 0); pavimento VTU sull'ultimo intervallo |
| C3 | Newmark bulk vs semidiscreto (modo coseno; oscillatore uniforme) | bulk iperbolico, integrazione temporale | eigenmodo Q1 semidiscreto (A-sem); ODE DOP853 indipendente (C) | RMS finale a mesh fissa | coseno: 3.4715e-4 / 8.6981e-5 / 2.1759e-5, p=1.9968; oscillatore: 7.3064e-6 / 1.8064e-6 / 5.0638e-7, p=2.0161 | pre-registrato [1.5, 2.5] | **SUPPORTATA** (nel caso testato) | sistema semidiscreto, non PDE continua; nessun bordo; ultimo intervallo oscillatore vicino al pavimento |
| C4 | Okumura, diss. 2021 Ch.4 Ex. 2 (1D) | parabolico bulk + estremi dinamici | energia pubblicata, **curva raster digitalizzata** (B) | RMSE/MAE/max su 391 campioni, 0.01≤t≤0.4, asse energia 0–90 | RMSE 0.9881, MAE 0.8804, max 1.8323; energia OOFEM 77.42 → 2.22e-6, monotona; parametri e IC replicati | **descrittivo** (nessun criterio pre-registrato; accordo entro incertezza raster dichiarata) | **SUPPORTATA PARZIALMENTE** | incertezza raster ~0.1–1 unità; energia FEM con Gauss 2×2 (non esatta per il quartico Q1); **mismatch 391 vs 372 campioni non riconciliato**; FD (paper) vs Q1 (OOFEM) |
| C5 | Xiao–Chen 2025, Caso 4 | parabolico 2D + bordo dinamico | figure raster; **ε omesso nel paper** (E) | morfologia u=0.5 a T=0.1; energia | forma coerente a vista (intersezioni 0.2578/0.7422, sommità 0.4922 vs ~0.25/0.75, 0.5); energia 144.16 → 49.67 monotona | **descrittivo/esplorativo** (ε=0.1 è assunzione esplicita, non dato degli autori) | **SUPPORTATA PARZIALMENTE** (solo qualitativo) | nessuna validazione quantitativa possibile (ε assente, dati solo raster, h=0.05 vs 1/300); n20 non spazialmente convergente a T=0.1 |
| C6 | Caso parabolico Calatroni–Colli | parabolico 2D bulk–superficie | teoria 2012 (esistenza/unicità) + raffinamenti (D) | RMS nodali spaziali vs mesh 40×40; differenze temporali consecutive | spaziale 10/20/40 a T=0.1, dt=0.0025: RMS 1.27871e-4 / 2.59027e-5 (rapporto 4.94, ordine apparente 2.30); temporale dt=0.005/0.0025: 2.28476e-4 / 1.15949e-4 (rapporto 1.97, ordine osservato ≈0.98, coerente con BE) | **descrittivo** (il paper non fornisce valori esatti: nessun criterio quantitativo esistente) | **SUPPORTATA PARZIALMENTE** (descrittivo) | **NON è una verifica A**: i teoremi non danno valori della soluzione; i CSV conservano solo metriche aggregate e i raw non sono archiviati → errori originari NOT VERIFIABLE; il raffinamento da solo non certifica l'equazione |
| C7 | Raffinamento spaziale interno 2D parabolico | parabolico 2D | mesh N=40 / N=60 (D) | L2 bulk/traccia su griglia comune | T=0.0625 (dt=6.25e-4): bulk .11731→.05461→.02994 (N=20/30/40 vs 60), convergente; T=0.1 (dt=5e-4): bulk .07998→.06211 ma **traccia .08963→.08861 quasi stazionaria** | **descrittivo** (riferimento numerico, non esatto) | **SUPPORTATA PARZIALMENTE** | N=40/60 non esatte; la sequenza a T=0.1 non dimostra robustezza spaziale; nessun ordine asintotico certificato |
| C8 | Raffinamento temporale parabolico (BDF1) | parabolico 2D, mesh 20×20 | dt=1e-4 riferimento numerico (D) | L2 mass-weighted bulk/traccia; differenze adiacenti | p ≈ 0.999 (bulk: 8.6004e-5 / 4.3031e-5), coerente con 1° ordine | pre-registrato [0.75, 1.25] | **SUPPORTATA** (nel caso testato) | verso il sistema semidiscreto, non la PDE; attenzione: nel report 2D i valori "vs dt fine" erano etichettati come "differenze adiacenti" (errore descrittivo, ordine non intaccato) |
| C9 | Caso accoppiato iperbolico — temporale | bulk inerziale + superficie 1° ordine | ODE Q1/P1 semidiscreta DOP853 indipendente, **corretta** (C) | FE RMS bulk/superficie a T=0.02, mesh 16×8, dt=0.002/0.001/0.0005 | dopo la correzione del riferimento: bulk 7.50e-7 / 1.91e-7 / 4.96e-8 (ordini 1.972 / 1.948); superficie 2.27e-6 / 5.78e-7 / 1.42e-7 (ordini 1.975 / 2.027); ΔE(T)−ΔE(OOFEM) = 6.5e-6 ∝ O(dt²) | pre-registrato [1.5, 2.5] | **SUPPORTATA** (nel caso testato) | singola combinazione di parametri (ρ=1, ζ=0, c_Γ=1, β=1, IC data, T=0.02); il riferimento è numerico, non la soluzione esatta della PDE |
| C10 | Caso accoppiato iperbolico — spaziale (T1, T1-A, T1-D, T1-B) | bulk inerziale + superficie 1° ordine, 2D | proxy DOP853 n64×32, poi n128×64 (D) | FE RMS massa-pesata (normalizzata) + max nodale, bulk e superficie, a tempi uguali | errori a T: n8 1.299e-3, n16 7.903e-4, n32 3.822e-4 (bulk FE); differenze di famiglia 4.141e-4 (32→64), 4.874e-5 (64→128); ordini 0.7173 / 1.0478 (OOFEM) e rapporto famiglia 3.0869; convergenza monotona FE RMS fino a h=0.015625 | pre-registrato per T1: monotonia + ordini in [1.5, 2.5] su entrambi gli intervalli | **SUPPORTATA PARZIALMENTE**: esistenza **SOSTENUTA**; ordine **NON DIMOSTRATA**; singolarità **NON DEMOSTRATA** | nessuna soluzione esatta; il proxy è una discretizzazione; bulk max non monotono (massimo mobile); vedi sezione E |
| C11 | Energia — caso accoppiato (T1-B n128; 16×8) | iperbolico accoppiato | legge continua (F) | E_tot(t), monotonia, R(T) | T1-B: 1.93271019 → 1.92746381 (−0.27%), monotono su 80 passi; 16×8: 1.93160616 → 1.92995451, R(T) ≈ 1.8e-7 (velocità ricostruita) | **descrittivo** (nessuna identità discreta dimostrata per Newmark nonlineare) | **SUPPORTATA PARZIALMENTE** | dissipazione a secanti; velocità non esportata; pavimento di ricostruzione ~7e-7 |
| C12 | Energia — parabolica (20×20 tre dt; 16×8; 2D secante) | parabolico | legge continua (F) | E(t), residuo normalizzato, additività | E0=1.28445597 → 1.20912/1.20902/1.20897; residuo 4.231e-3 → 2.289e-3 → 1.201e-3; 16×8: R(T) −2.947e-4 / −1.525e-4 / −7.766e-5 (ordine 0.987); 2D secante: identità entro 5.6e-17, max residuo 1.73e-3 | **descrittivo** | **SUPPORTATA PARZIALMENTE** | secanti ≠ dissipazione esatta; quadratura legacy 2×2 (energy.csv); pannello superiore della figura secante dimensionalmente errato (non usare) |
| C13 | Energia — Nizovtseva bulk | iperbolico bulk | legge continua (F) | E(t), R(T) | E0=3.16242072, ET=3.08549010, R=3.465e-6 | **descrittivo** | **SUPPORTATA PARZIALMENTE** | u_t stimata con np.gradient: la cinetica interna è **NOT VERIFIABLE** |
| C14 | Lee et al. 2024, equilibrio t=5 | parabolico 1D | profilo esatto tanh(x/(√2 ε)) (A) | L2 continuo, max nodale | max 1.30468e-2 (consistente); **L2 6.0675e-3 (audit) vs 2.61689e-3 (ricalcolo): 2.32× irrisolto**, nessuna derivazione persistente del valore vecchio | **descrittivo** | **SUPPORTATA PARZIALMENTE** | discrepanza L2 irrisolta: il valore L2 non deve essere citato senza chiarimento |
| C15 | Superficie con inerzia di 2° ordine | — | — | — | non implementata nel modello `wave 1` | — | **NON APPLICABILE** | limite di modello, non di validazione |
| C16 | Wang / Akrivis (figure raster energetiche) | parabolico 2D | figure raster (B) | — | raw VTK solo in /tmp, non versionati; alcuni run NaN | **descrittivo** | **NON APPLICABILE** (dati non archiviati → non riproducibili) | — |
| C17 | Suite MMS (Hochbruck–Leibold adattato, Chen 2025 §6.2, beta_mms interno) | — | sorgenti costruite | — | **escluse per vincolo "no MMS"**; le metriche relative (es. ordini stagnanti 0.33/0.24 sul disco H-L) non sono evidenza di validazione | — | **NON APPLICABILE** (esclusione da protocollo) | — |

Nota di coerenza: la matrice "finale" archiviata (`validation_matrix_final.md`) riporta
"Sistema iperbolico accoppiato: **FAIL**" e "Convergenza spaziale: PARTIAL" con la nota
"la causa del mismatch non è determinata". Entrambi gli stati sono **superati**: la causa
era nel riferimento Python (sezione D.2) e la convergenza spaziale è stata completata in
T1/T1-A/T1-B (sezione E). Questo rapporto è lo stato più recente e autorevole.

---

## D. Correttezza, riproducibilità e storia degli errori

### D.1 Correttezza dell'implementazione

- **Verifica della formula esatta di Nizovtseva** (nuova, in `audit_iniziale_2026-10-27.md`
  §1, `nizovtseva_pde_residual_check.py`): il residuo della Eq. (29) si annulla solo per
  s=3/2 e λ²=(9τ+2)/4; max|residuo| = 1.1e-16 con τ=0.2; `v₀ = 1.5 f₀(1−f₀)` esatto. La
  formula è usata come riferimento con questo status verificato.
- **Ispezione statica termine per termine** (sezione B.1): equazioni, segni, coefficienti,
  accoppiamento, quadratura coerenti con il modello.
- **Coerenza OOFEM dopo la correzione del riferimento** (`jacobian_fix_verification_report.md`
  §4, 9 righe termine per termine): residuo, matrici, capacità, reazione, flusso, termini di
  bordo — tutti coerenti; l'unica spiegazione della discrepanza originaria è il riferimento
  Python (sotto).
- **Comportamento coerente**: ordini ≈2 dove il problema è regolare e il riferimento è
  adeguato (C1–C4, C9); invarianza esatta a inversione/permutazione delle linee (test 8×8,
  diff 0); simmetria 180° a precisione macchina (2.2e-16) nel caso accoppiato.

### D.2 Il difetto del riferimento Python indipendente (identificato, corretto, documentato)

**Problema:** nello script di riferimento originale (`run_coupled_independent.py`, 232 righe,
conservato intatto e non modificato) l'assemblaggio degli elementi linea usava
`dN = [−0.5, 0.5]/length` e `w = length·GAUSS_W`. Correttamente: `dN = [−0.5, 0.5]·(2/length)`,
`w = 0.5·length·GAUSS_W` (e analogamente per `rs`). Effetto: **Ms ×2, Ks ×½, rs ×2** — il
riferimento risolveva la diversa equazione di bordo
`2 c_Γ ψ_t − (k/2)Δ_Γψ + k∂ₙu = 2βm(ψ−ψ³)`.

**Conseguenza osservata (pre-correzione):** errori OOFEM-vs-ODE quasi costanti al
raffinamento (bulk 4.6937e-4 / 4.6967e-4 / 4.6975e-4; superficie 2.1058e-3 / 2.1070e-3 /
2.1073e-3, dt 0.002/0.001/0.0005), ordini ~0/negativi → la matrice "finale" aveva assegnato
FAIL al sistema accoppiato.

**Correzione** (`check_jacobian_fix.py`, verificata contro le matrici esatte di un'edge P1:
`Ms_e = (L/6)[[2,1],[1,2]]`, `Ks_e = (1/L)[[1,−1],[−1,1]]`, ΣMs = 6.0 = perimetro, ΣMb = 2.0 =
area, nnz 144) e **risultati post-correzione**: bulk 7.4998e-7 / 1.9119e-7 / 4.9563e-8
(ordini 1.972 / 1.948); superficie 2.27e-6 / 5.78e-7 / 1.42e-7 (ordini 1.975 / 2.027);
energia del riferimento corretta monotona, ΔE(T) − ΔE(OOFEM) = 6.5e-6 coerente con O(dt²).
Le metriche superficiali "prima" erano gonfiate di √2 dalla Ms errata (plateau corrette
1.4891/1.4899/1.4901e-3).

**Conseguenza sulle conclusioni:** il verdetto "accoppiato FAIL" è **recepito come errato**
(causa lato riferimento, non OOFEM); i documenti che lo riportano sono marcati superati
(§H). Il codice OOFEM non è stato modificato in nessun punto della campagna.

### D.3 Riproducibilità verificata (bit-identica dove testata)

| Prova | Esito | Fonte |
|---|---|---|
| Rerun OOFEM parabolico 20×20 (19 s), Nizovtseva N=160 (12 s), accoppiato n16×8 dt=0.0005 (1 s) | campi finali identici agli archiviati (RMS/max diff = 0), 0 errori, 0 warning, input byte-identici | `reproducibility_runs.md` |
| T1-A: rerun OOFEM n32×16 (case.in byte-identico, `cmp`) | steps 1–39 **bit-identici** all'originale (max\|Δ\| = 0.000e+00); 80/80 VTU completi, 10 s, RSS 32 MB, exit 0 | `t1a_report.md` §3 |
| T1-B: assembly matrici fresca vs cache n128×64 | diff massima 0.000e+00 su Kb, Ks, Mb, Ms (sha256 matrici 68c85bfd…1330) | `t1b_report.md`, `run.log` |
| Riferimento DOP853 corretto rieseguito in memoria | bit-identico (max\|ΔU\| = 0, max\|ΔV\| = 0, max\|ΔE\| = 0) | `audit_iniziale_2026-10-27.md` §3.4 |
| 16/16 metriche T1 ricalcolate end-to-end (norme, nesting, PVD, ordinamento nodi) | diff 0.000e+00; nesting esatto (d² < 1e-26); tempi PVD Δ ≤ 3.5e-18 | `static_audit_report.md` §1.1 |
| Determinismo T1-B (1 core, taskset) | RSS 109 MB, 21.0 min, nfev 2117, exit 0, tutti i guard intatti | `t1b_report.md`, `resource_report.txt` |

**Limiti ambientali:** build `rel-2.1-5430-g2cc1a7706-dirty` (debug, non release immutabile);
raw di diversi test solo in `/tmp` (non archivio durevole); nessun lockfile/requirements;
percorsi assoluti in alcuni runner (Dipendenze: Python 3.12.3, NumPy 1.26.4, SciPy 1.11.4,
Matplotlib 3.6.3). La ripetibilità bitwise vale per questo eseguibile e questi input.

### D.4 Dipendenze non corrette, archivi incompleti, problemi di post-processing

Rimasti aperti (nessuna correzione eseguita in questa campagna, per vincolo):

1. **CSV energetici legacy non riproducibili** dai loro stessi output:
   `energy_validation_finer.csv` (16×8, dt=0.00025): dichiara E(T)=1.9299885203,
   R=3.47e-5, **17 incrementi positivi**; ricalcolo dagli stessi `.out`:
   E(T)=1.9299547624, R=6.87e-7, **0 incrementi**. A dt=0.0005, cinetica legacy vs
   high-precision differisce di 1.8135e-5 (potenziale: −1.4e-8). I due CSV
   (`energy_time_series.csv`, `energy_validation_finer.csv`) e `energy_validation.png`
   sono marcati **non autorevoli** (non evidenza sul solver).
2. `energy.csv` (legacy): quadratura Gauss 2×2 sul potenziale quartico Q1 → sottointegrata,
   solo diagnostica qualitativa.
3. `parabolic_2d/energy_balance_secant.png` pannello superiore: confronta −ΔE/Δt (tasso)
   con dissipazioni integrate non divise per Δt (4967.6 vs 0.4959 al primo passo) —
   **dimensionalmente errato, non va usato**; il CSV e il pannello del residuo sono coerenti.
4. Link rotto nel report 2D: `energy_interface.csv` indicato nel padre, reale in
   `parabolic_2d/runs/energy_n20_dt0p0001_every_step/`.
5. **Discrepanze metriche irrisolte**: Okumura 391 vs 372 campioni (RMSE 0.98810 vs
   0.99288) mai riconciliate; Lee L2 2.32× (C14).
6. **Archivi incompleti**: n32×16 originale step 40–80 a 0 byte (recuperato da T1-A;
   originale lasciato intatto); `energy_validation_runs_highprecision/parabolic_n32x16_dt0p0005`
   4/40 VTU; `final_reproduction/energy_metrics_verified/` vuoto;
   `energy_parabolic_single_test_run/` solo input; suite legacy `nizovtseva_*` v1 e
   `hl_disk_*` solo stato finale (non re-derivabili).
7. **Velocità non esportata** nei VTU (solo `u`/`Temperature`): ogni cinetica iperbolica è
   ricostruita (np.gradient o ricorrenza Newmark su spostamenti testuali a ~8 cifre);
   pavimento di ricostruzione ~7e-7 a dt=0.00025.
8. Runner che sovrascrivono output (nessun manifest/hash input↔result) — rischio di riuso
   silenzioso; `postprocess_energy.py` è l'unica eccezione positiva (modalità esclusiva).

### D.5 Storia completa dei problemi, correzioni e conseguenze (sezione dedicata)

La storia degli errori **non è eliminata**: ogni voce include stato attuale.

| # | Problema | Dove | Correzione / stato | Conseguenza sulle conclusioni |
|---|---|---|---|---|
| 1 | Jacobiano degli elementi linea nel riferimento Python (Ms ×2, Ks ×½, rs ×2) | `run_coupled_independent.py` (originale, intatto) | corretto in `check_jacobian_fix.py`; riferimento buggy conservato | il verdetto "accoppiato FAIL" (matrice "finale", report 26-10) è **recepito come errato**; l'accordo temporale accoppiato è SUPPORTATA (C9). OOFEM non modificato |
| 2 | Convenzione temporale con offset di un passo nei confronti T1 (U[step−1] invece di U[step]) | post-processing T1 | corretta; Δt verificato ≤ 3.5e-18 | floor raffinato ~1e-5 → ~1e-8; ordini cambiati < 0.01; le tabelle T1 usano la convenzione corretta |
| 3 | Fit `err(h) = a·h² + b·h` (a=+9.48e-3, b=+2.83e-3, residuo 9%) in `spatial_verification_report.md` §6 D5(iii) | T1 | **ritirato** (erratum `static_audit_report.md` §4): calcolato sui dati pre-correzione; sui dati corretti a tempo uguale a<0, residuo ~72%; rifatto in T1-A §5 a tempo unico: bulk a=−8.77e-3 (residuo 5.18%), surf a=−3.72e-2 (16.96%) | il fit h²+h **non è supportato** dai dati: nessuna decomposizione "termine O(h) + O(h²)" identificabile; non può essere usato come prova di singolarità |
| 4 | Corruzione archivio n32×16 originale: step 40–80 (41 VTU + mesh.npz + run.log) a 0 byte | `hyperbolic/coupled_independent_space_n32x16_dt0p00025/` | archivio lasciato **intatto**; T=0.02 recuperato dal rerun T1-A in directory separata (`t1a_n32x16_rerun/`, steps 1–39 bit-identici) | nessun dato perduto; il rerun è la fonte primaria per n32 a T (sezione E) |
| 5 | Stima DOF n128×64 errata in T1-D (8449) | `t1d_report.md` §5.1 | misurata in T1-B FASE 0: **8385 nodi** (129×65) | nessuna conseguenza sui risultati; le stime di costo restano valide (confermate in FASE 1) |
| 6 | `t1b_cost_log.csv`: riga "full" scritta con chiavi colonna diverse da "pilot" (valori disallineati, es. wall_s=27.58 che era asm_s) | T1-B | riga riscritta con schema a 14 colonne coerente, valori da `run.log`/`resource_report.txt`; dichiarato in `t1b_report.md` §8 | nessuna conseguenza scientifica (solo tracciamento costi) |
| 7 | Due bug di script in `t1b_validate.py` (FASE 3): (a) indicizzazione di U64 (larghezza 2145) con mappe n128 (larghezza 8385) → IndexError; (b) `bnd` (array di indici) trattato come maschera booleana nel moving-max → AssertionError | T1-B | entrambi corretti prima della run finale di validazione (mappe `i8_64/i16_64/i32_64`; assert a coordinate) | i dati di FASE 2 (integrazione) non sono mai stati ritoccati; i CSV di risultato finali derivano dalla run corretta |
| 8 | Stima nfev T1-B: ~1130 previsti, 2117 effettivi | T1-B FASE 2 | documentata: DOP853 a n128 usa 26.5 RHS/step (vs 13.6 a n64); il costo per RHS era accurato (595 vs 591 ms) | la run è rimasta entro tutti i limiti (21.0 min < 45, 109 MB < 2 GB, 1 core) |
| 9 | Report 2D: "differenze consecutive in norma di massa" erano in realtà errori rispetto a dt=1e-4 | `parabolic_2d_validation.md` §3.3 | **non corretto nel file originale** (per vincolo): marcatto come errore descrittivo in `validation_metrics_audit.md`; l'ordine ~1 resta corretto ricalcolato dalle vere coppie adiacenti | nessun impatto sui verdetti; la citazione deve usare la formulazione corretta |
| 10 | `OOFEM_AllenCahn_Validation_Report.md` (26-10) e `validation_matrix_final.md/.csv` anteriori alla risoluzione: "La causa resta da isolare" / "accoppiato FAIL" | root | **non modificati** (vincolo): superati da `hyperbolic_coupled_discrepancy_audit.md` + `jacobian_fix_verification_report.md` e dichiarati obsoleti qui | questo rapporto è lo stato autorevole; i documenti obsoleti restano in archivio per tracciabilità |
| 11 | `literature_comparison/README.md`: "il solver parabolico è validato… nel caso 1D" senza qualifiche | letteratura | segnalato in `audit_iniziale_2026-10-27.md` §4.1, non riformulato (vincolo) | lo stato corretto è C4: SUPPORTATA PARZIALMENTE con riserve raster/quadratura |
| 12 | Metriche superficiali "prima" gonfiate di √2 (2.1058e-3 vs plateau corretta 1.49e-3) | T1 / report | spiegato dal Ms errato del riferimento (voce 1); annotate in `audit_iniziale_2026-10-27.md` §4.1 | le tabelle "prima" vanno lette con questa qualificazione |
| 13 | Figura `energy_validation.png` con provenienza non verificata (residuo legacy rumoroso) | root | marcata **non autorevole** negli audit; non rigenerata | non usare per conclusioni sul residuo iperbolico |

---

## E. Convergenza spaziale del caso accoppiato (sezione centrale)

**Caso:** Ω = [0,2]×[0,1], Γ = ∂Ω, T = 0.02, dt = 2.5e-4 (80 step), RHO=KB=MB=CS=KS=MS=1,
ETA=ζ=0, IC u₀ = 0.2 + 0.1·cos(πx/2)·cos(πy), v₀ = 0, nessun BC (Neumann naturale).
Simmetria 180° u(2−x,1−y)=u(x,y) a 2.2e-16. I tre oggetti distinti — **(1)** differenza di
livello di famiglia (proxy vs proxy), **(2)** errore OOFEM vs semidiscreto a mesh uguale
("floor" del solutore), **(3)** errore vs soluzione esatta (non disponibile) — sono
separati in ogni affermazione qui sotto.

### E.1 Mesh, h e tempi

| Livello | mesh | nodi | h (verificato: h_x = h_y, lati quadrati, linee = h) |
|---|---|---|---|
| 1 | n8×4 | 45 | 0.25 |
| 2 | n16×8 | 153 | 0.125 |
| 3 | n32×16 | 561 | 0.0625 |
| 4 | n64×32 (proxy) | 2145 | 0.03125 |
| 5 | n128×64 (proxy, T1-B) | 8385 = 129×65 | 0.015625 |

Raffinamento esatto ×2 (rapporti verificati a 1e-14, `t1d_h_geometry.csv`) → base log2
corretta per gli ordini. Mesh annidate per coincidenza esatta di nodi (d² < 1e-26), senza
interpolazione. **Tempi di confronto:** OOFEM step s ↔ U[s] del proxy a t = s·dt (convenzione
verificata: Δt ≤ 3.5e-18; correzione dell'off-by-one, voce D.5-2). **Tutti gli ordini sono
rapporti a tempo uguale.**

### E.2 Definizioni delle norme

- **bulk FE RMS** = √(eᵀ M_b e / AREA), AREA = 2.0 (massa Q1 del livello coerente,
  normalizzata per area);
- **surf FE RMS** = √(eᵀ M_s e / PERIM), PERIM = 6.0, con **Ms corretta** (fattori
  0.5·L·w, dN = (2/L)·[−.5,.5]);
- **bulk max** su nodi interni (bordo escluso), **surf max** su nodi di bordo;
- restrizione per coincidenza esatta di coordinate (nessuna interpolazione di forma).

### E.3 Errori a T = 0.02 (vs proxy n128×64) e differenze di livello di famiglia

Fonti: `t1b_metrics_T.csv` (T1-B); n32 = rerun T1-A (bit-identico a step 1–39 dell'originale).

| Quantità | n8×4 (h=.25) | n16×8 (.125) | n32×16 (.0625) | famiglia 32→64 (.03125) | famiglia 64→128 (.015625) |
|---|---:|---:|---:|---:|---:|
| bulk FE RMS | 1.299267e-3 | 7.902543e-4 | 3.822430e-4 | 4.141362e-4 | **4.874007e-5** |
| bulk max | 8.346605e-4 | 1.910358e-3 | 2.957438e-3 | 3.123389e-3 | 7.484989e-4 |
| surf FE RMS | 2.508695e-3 | 1.900346e-3 | 8.372465e-4 | 9.264513e-4 | 1.642138e-4 |
| surf max | 5.873672e-3 | 4.263484e-3 | 1.832372e-3 | 2.229713e-3 | 4.826234e-4 |

Le tre colonne "n8/n16/n32" sono OOFEM-vs-proxy128 (oggetto (2)+(1) sommati); le due colonne
"famiglia" sono DOP853-vs-DOP853 (oggetto (1) puro).

### E.4 Ordini a T

| Quantità | 8→16 (OOFEM) | 16→32 (OOFEM) | rapporto famiglia 32→64 / 64→128 |
|---|---:|---:|---:|
| bulk FE RMS | 0.7173 | 1.0478 | 3.0869 |
| bulk max | −1.1946 | −0.6305 | 2.0610 |
| surf FE RMS | 0.4007 | 1.1825 | 2.4961 |
| surf max | 0.4622 | 1.2183 | 2.2079 |

Controlli di coerenza (eseguiti per questo rapporto): log₂(1.299267e-3/7.902543e-4) = 0.7173 ✓;
log₂(7.902543e-4/3.822430e-4) = 1.0478 ✓; log₂(4.141362e-4/4.874007e-5) = 3.0869 ✓;
idem per le altre 8 righe ✓.

**Interpretazione obbligatoria:** le prime due colonne misurano (1)+(2); la terza colonna è
un **singolo rapporto fra due differenze di livello di famiglia** (diagnosi del regime della
famiglia FE). Il criterio pre-registrato [1.5, 2.5] non è soddisfatto da nessuna delle
prime due colonne; il valore 3.0869 **non è un ordine di convergenza misurato su tre
soluzioni** e **non dimostra convergenza di terzo ordine**: è non costante nel tempo
(§E.5), è costruito su due differenze e usa un solo rapporto.

### E.5 Dipendenza dal tempo (T1-D + T1-B)

- **Nessun recupero verso l'ordine 2 a nessun istante archiviato** (T1-D §2): l'ordine
  8→16 bulk FE sta in 0.570–0.698 su tutti gli 80 step; 16→32 in 0.73–0.78 per t ≤ 0.00975
  e arriva a 0.983 solo a T (confronto temporale, segnato).
- **Serie di famiglia 64→128 (T1-B, solo DOP853)**: ordine bulk FE −0.4462 (s1) → −0.1140
  (s20) → 0.7156 (s45) → 3.0869 (s80); la differenza e64→128 è **non monotona** nel tempo
  (picco ~1.55e-4 a s45, poi decresce a 4.87e-5 a T). Il valore a T è quindi il risultato
  di una storia temporale non monotona, non di un regime asintotico spaziale stabile.
- L'errore FE RMS **è monotono decrescente** al raffinamento a ogni istante; il **bulk max
  non è monotono** (argmax mobile, §F).

### E.6 Stabilità del livello di riferimento (robustezza rispetto alla scelta del proxy)

Errori degli stessi campi OOFEM calcolati vs proxy n64 e vs proxy n128 (T1-B
`t1b_proxy_stability_T.csv`): differenze assolute 2.93e-5 / 2.85e-5 / 3.19e-5 (bulk FE RMS,
n8/n16/n32) = **2.2% / 3.5% / 7.7%** (≤8%); i cambi d'ordine sono ≤0.07. Il raffinare il
proxy da n64 a n128 **non cambia le conclusioni**.

### E.7 Errori regionali e massimo mobile

- **Zone** (r < 0.5 dagli angoli; A = (0,1)/(2,0) "coppia attiva", B = (0,0)/(2,1)): a T,
  energia di errore in A: **65.8%** per la famiglia 64→128, **74.6%** per OOFEM n32-vs-proxy128
  (B: 28.0% / 14.2%; resto: 6.1% / 11.2%).
- **Massimo mobile a T** (primo nodo interno vicino all'angolo attivo; T1-B):
  r_c = 0.3536 → 8.347e-4 (n8); 0.1768 → 1.910e-3 (n16); 0.0884 → 2.957e-3 (n32);
  0.0442 → 7.485e-4 (famiglia n64). Pendenza log-log a 4 punti: **p = −0.016 ≈ 0** —
  non una legge di potenza.

### E.8 I tre oggetti, separatamente

1. **Differenze di famiglia (1):** 32→64 = 4.14e-4, 64→128 = 4.87e-5 (bulk FE RMS a T):
   la famiglia **sta ancora convergendo** a h = 0.015625 (rapporto 8.5, nessun plateau
   visibile). La struttura è della famiglia FE per questo problema, non dell'implementazione:
   l'auto-differenza del proxy (n32 vs n64, senza OOFEM) mostra la stessa concentrazione
   angolare (T1-D: zona A 0.722/0.739) e la stessa non monotonia del max.
2. **OOFEM vs semidiscreto, stessa mesh (2):** floor 2.10e-8 / 2.44e-8 / 2.60e-8 (bulk FE
   RMS, n16@39, n16@80, n32@39) — 6–7 ordini sotto l'errore spaziale. A n32 a T, l'errore
   OOFEM-vs-proxy128 (3.82e-4) è **entro l'8%** della differenza di famiglia 32→64
   (4.14e-4): la componente OOFEM misurabile è al livello della floor e l'errore visibile a
   n32 è in sostanza la differenza di famiglia 32→64. **Non esiste nessun run OOFEM n64×32**:
   la floor a h = 0.03125 resta sconosciuta (attività futura 1, §J.3).
3. **Errore vs soluzione esatta (3):** non disponibile per il caso accoppiato (vincolo no-MMS;
   non in letteratura per questo sistema). **Nessuna affermazione di accuratezza assoluta
   è fatta in questo rapporto.**

### E.9 Perché la convergenza (esistenza) è sostenuta ma il tasso asintotico no

**Sostenuto (convergenza esistenza, fino a h = 0.03125 per OOFEM e a h = 0.015625 per la
famiglia):**
- FE RMS monotono decrescente a ogni raffinamento e a ogni istante archiviato;
- la famiglia prosegue a raffinare sotto h = 0.03125 (64→128 = 4.87e-5, 8.5× sotto 32→64):
  **nessun plateau** → l'errore non satura;
- convergenza puntuale a nodi fissi (T1-D: 8.34e-4 → 3.34e-4 → 2.74e-5 su (0.25, 0.75));
- l'errore OOFEM a n32 è dominato dalla famiglia, non da un difetto del solutore (§E.8-2);
- la scelta del proxy non cambia le conclusioni (§E.6).

**Non stabilito (tasso asintotico):**
- gli ordini (1)+(2) a T sono 0.72 e 1.05 (bulk FE) — fuori dal criterio [1.5, 2.5],
  crescenti col raffinamento, non costanti;
- il rapporto di famiglia 3.09 (bulk FE) è un singolo rapporto di differenze, non costante
  nel tempo (−0.45 → 3.09) e **non costituisce prova di terzo ordine** (né, per converso,
  di secondo: il valore non è stabile né in h né in t);
- il bulk max è non monotono a ogni istante archiviato (argmax mobile, §F);
- mancano 2–3 livelli aggiuntivi di h per identificare un regime (richiederebbe
  OOFEM n64×32 e/o proxy n256×128, non eseguiti);
- l'eventuale limite di regolarità del problema (singolarità angolare) resta
  indistinguibile da un effetto pre-asintotico di h finito con i dati esistenti.

**Formulazione ammissibile:** "La famiglia FE Q1 converge puntualmente e in FE RMS per il
caso accoppiato testato, senza plateau fino a h = 0.015625; il tasso di convergenza non è
identificato: gli ordini misurati non sono costanti e nessun regime asintotico è
stabilito." Non è ammissibile: "il caso accoppiato converge con ordine 2" né "con ordine 3".

---

## F. Analisi angolare e regolarità

**Osservazioni (fatti, con fonti):**
1. ~71–75% dell'energia di errore nella zona A (angoli (0,1)/(2,0)) su tutte le mesh a T
   (0.708/0.725/0.739, T1-A; 0.746 OOFEM n32 e 0.658 famiglia n64, T1-B): **stabile col
   raffinamento**.
2. La stessa concentrazione è presente **al passo 1** (0.735–0.766, T1-D): non è
   accumulación tarda nel tempo.
3. Compare anche nell'**auto-differenza del proxy** (n32 vs n64 DOP853, senza OOFEM):
   zona A 0.722 (t=0.00975) / 0.739 (t=0.02) → proprietà della famiglia FE per questo
   problema.
4. Simmetria 180° esatta (2.2e-16); i gemelli di ogni nodo hanno errore identico a ~1e-17.
5. Verso l'angolo, |err| ~ r_c^p con **p ≈ −2 negli interni** (range −1.9…−2.8,
   regione/slice-dipendente; molto più piatto al bordo, −0.14…−0.62) — range, non legge
   unica (audit statico §1.3).
6. **Massimo mobile**: l'argmax interno è sempre il primo nodo vicino all'angolo attivo, e
   "insegue" l'angolo al raffinamento (r_c = h√2); la non monotonia del bulk max (ordini
   −1.19/−0.63) è la conseguenza di misurare il max su nodi diversi a mesh diverse, e
   coesiste con la convergenza puntuale a nodo fisso.
7. A 4 punti (T1-B, r_c = 0.3536…0.0442) il massimo a T è **essenzialmente piatto**
   (p = −0.016); il comportamento r → 0 (crescita continua o saturazione) **non è
   identificabile** con i dati esistenti (servirebbe r_c < 0.044).

**Interpretazioni plausibili (non dimostrate):**
- limite non regolare (singolare) della soluzione esatta negli angoli (0,1)/(2,0) del
  sistema accoppiato: coerente con tutti i fatti 1–7 (concentrazione stabile dal passo 1,
  p ≈ −2 interno, max mobile);
- in alternativa o in aggiunta: effetto pre-asintotico di h finito (la famiglia non ha
  ancora raggiunto il regime).

**Dichiarazioni non ammesse con i dati esistenti:**
- "è dimostrata una singolarità angolare" — **no**: il proxy è esso stesso una
  discretizzazione, la soluzione esatta è sconosciuta, ogni "errore" è una differenza di
  livello; la concentrazione è condizione necessaria ma non sufficiente;
- "gli ordini ridotti sono causati dalla singolarità" — **no**: causa non attribuita
  automaticamente; il fit h²+h (meccanismo atteso per un termine O(h) singolare) è
  **ritirato** (a<0, voce D.5-3);
- alcuna affermazione sul comportamento r → 0 (punto 7).

---

## G. Bilancio energetico

**Leggi continue (coerenti con il codice, `energy_theory.md`):**
- parabolico: Ẽ_pot = −∫_Ω c u_t² dx − ∫_Γ c_Γ ψ_t² dS ≤ 0, con
  E_pot = ∫_Ω(k/2|∇u|² + mF(u))dx + ∫_Γ(k/2|∇_Γψ|² + βmF(ψ))dS, F(u) = (u²−1)²/4;
- iperbolico `wave 1`: E_tot = E_pot + ∫_Ω(ρ/2)u_t²dx,
  Ẽ_tot = −∫_Ω ρζ u_t² dx − ∫_Γ c_Γ ψ_t² dS ≤ 0.
Il segno del termine di bordo +k∂ₙu si cancella nel bilancio accoppiato; le formule sono
confermate termine per termine sul sorgente (sezione B). **Nessuna identità discreta è
dimostrata nel codice per Newmark nonlineare** (α=1/BE parabolico incluso): la legge
continua non stabilisce da sola una legge di dissipazione a passo finito.

**Osservato (tendenze, tutti run moderni):**
- T1-B (n128, 80 passi): E_tot 1.93271019 → 1.92746381 (−0.27%), **monotono decrescente**,
  cinetica 0 → 2.48e-3, potenziale bulk 4.894e-1 → 4.911e-1 (cresce: conversione
  cinetica, coerente con la legge del totale);
- parabolico 20×20: E 1.28445597 → 1.20912/1.20902/1.20897 (tre dt), nessun incremento;
- 16×8: R(T) = E(T)−E(0)+D(T) = −2.947e-4 / −1.525e-4 / −7.766e-5 (parabolico, dt
  0.002/0.001/0.0005, ordine curva 0.987); 3.692e-6 / 9.156e-7 / 3.473e-7 (iperbolico,
  poi plateau ~7e-7 a dt=0.00025 = pavimento della ricostruzione);
- 2D secante: identità `ΔE + D_secante = R` a 5.6e-17 su 1000 righe, max |R|/|ΔE| = 1.73e-3;
- Nizovtseva: E0=3.16242072 → ET=3.08549010, R=3.465e-6;
- **nessun incremento positivo** dell'energia totale in nessun run moderno.

**Limiti delle misure (dichiarati, non superati dalla monotonia):**
1. **Dissipazione per differenze secanti** (Δu/Δt): è un'approssimazione a passo finito,
   non la dissipazione esatta integrata;
2. **Quadratura**: i CSV recenti usano 3×3/5×5 (esatti per il quartico Q1 su questi elementi);
   il legacy `energy.csv` usa 2×2 (sottointegrato) — valori solo qualitativi;
3. **Velocità ricostruita, non esportata**: ricorrenza Newmark v_{n+1} = 2(u_{n+1}−u_n)/Δt − v_n
   su spostamenti testuali a ~8 cifre (divisione per Δt amplifica l'arrotondamento) o
   `np.gradient` (Nizovtseva): la cinetica interna è **NOT VERIFIABLE**; pavimento ~7e-7;
4. **CSV legacy non riproducibili** (D.4-1): i grandi "aumenti" e residui archiviati in
   `energy_validation_finer.csv` / `energy_time_series.csv` sono artefatti di
   post-processing, **non evidenza di violazione del solver**;
5. **Figura secante**: pannello superiore dimensionalmente incoerente (D.4-3) — non usare.

**Conclusione ammissibile:** le implementazioni energetiche sono coerenti con le leggi
continue in entrambi i rami e i residui moderni scendono col raffinamento coerentemente
collo schema; la monotonia osservata è un controllo interno, **non una validazione
completa né indipendente del solver** (manca: residuo variazionale indipendente, identità
discreta, velocità exportata).

---

## H. Riproducibilità e inventario degli artefatti

**Directory dedicate (risultati):** `hyperbolic/coupled_independent_space_n{8x4,16x8,32x16}_dt0p00025/`
(n32 originale: step 40–80 corrotti, intatti per tracciabilità), `hyperbolic/t1a_n32x16_rerun/`,
`hyperbolic/t1b_dop853_n128x64/` (+ `results/`), `hyperbolic/t1d_analytic_audit/`,
`hyperbolic/coupled_independent_spatial_jacobian_fix/` (reference npz n64/n32, report T1/fix,
audit statico), `parabolic_2d/runs/`, `literature_comparison/`, `energy_validation_runs*`.

**Artefatti integri e verificati** (riproducibilità D.3): run OOFEM iperbolici (Niz v2
n80/160/320, kink n40/80/160, modi/oscillatore, accoppiato n8/n16, dt triple), run
parabolici (20×20, suite 2D, 16×8 energia), CSV energetici moderni, riferimenti DOP853
corretti (npz bit-identici), output T1-A (80 VTU), output T1-B (npz 10.9 MB, matrici
sha256, run.log, resource_report, 7 CSV, 3 PNG), report e audit della campagna.
I file pesanti (VTU/npz) sono **referenziati, non copiati**.

**Legacy non riproducibili / non autorevoli** (marcati, non modificati):
`energy_time_series.csv`, `energy_validation_finer.csv` (D.4-1); `energy.csv` (quadratura
2×2); `energy_validation.png` (provenienza non verificata); raw Wang/Akrivis (solo /tmp);
VTU intermedi Okumura (solo stato finale); serie N40/dt Xiao–Chen (solo /tmp);
derivazione L2 di Lee (non persistente); raw Calatroni–Colli (solo CSV aggregati);
suite `nizovtseva_*` v1, `hl_disk_*` (solo stato finale);
`energy_parabolic_single_test_run/` (solo input); `final_reproduction/energy_metrics_verified/`
(vuota); `energy_validation_runs_highprecision/parabolic_n32x16_dt0p0005` (4/40 VTU).

**Corretti/superati** (marcati in D.5): riferimento Python buggy (intatto, sostituito
funzionalmente da `check_jacobian_fix.py`); `validation_matrix_final.md/.csv` e
`OOFEM_AllenCahn_Validation_Report.md` §4 (stato "FAIL"/"causa da isolare" superato);
T1 §6 D5(iii) (fit ritirato); convenzione off-by-one (corretta nei post-processing).

**Script di riferimento e post-processing** (percorsi completi nel README della directory
di questo rapporto): `run_coupled_independent.py`, `check_jacobian_fix.py`,
`run_coupled_refinement.py`, `run_bulk_references.py`, `run_stationary_kink.py`,
`validate_nizovtseva_bulk.py`, `t1b_dop853.py`, `t1b_validate.py`, `t1b_plots.py`,
`run_case.py`/`run_suite.py`/`analyze.py` (parabolico 2D), `analyze_energy*.py`,
`run_energy_balance.py`, `run_okumura_1d_strip.py`, `run_xiao_chen_droplet.py`,
`digitize_okumura_fig44.py`, `postprocess_energy.py`.

---

## I. Limiti residui

1. **Nessuna soluzione esatta** del problema accoppiato 2D (vincolo no-MMS; non disponibile
   in letteratura): ogni errore spaziale accoppiato è rispetto al proxy FE; l'accuratezza
   assoluta non è misurata.
2. **Ordine asintotico spaziale non stabilito** (E.9): ordini non costanti, regime non
   identificato; il criterio [1.5, 2.5] non soddisfatto su nessuna coppia OOFEM-vs-proxy.
3. **Singolarità angolare non dimostrata** (F): interpretazione coerente ma non risolta;
   comportamento r → 0 non identificabile con i dati esistenti.
4. **Limiti del riferimento numerico**: il proxy DOP853 (rtol 2e-12, atol 2e-14) è una
   discretizzazione, non la soluzione della PDE; la sua floor temporale a mesh uguale è
   ~1e-8 (n16/n32); il proxy n128 non è a sua volta verificato contro un livello più fine.
5. **Riproducibilità ambientale**: build dirty (rel-2.1-5430-g2cc1a7706-dirty, debug), raw
   in /tmp per diversi test, nessun lockfile, runner con sovrascritture e percorsi fissi:
   la ripetizione bitwise è garantita solo per questo eseguibile e questi input.
6. **Benchmark da grafici digitalizzati** (Okumura): incertezza raster ~0.1–1 unità di
   energia; mismatch 391/372 campioni non riconciliato; quadratura legacy 2×2 sull'energia
   FEM; FD vs Q1.
7. **Metriche energetiche**: dissipazione a secanti, quadratura dichiarata, velocità
   ricostruita (pavimento ~7e-7), nessuna identità discreta dimostrata per Newmark
   nonlineare, nessun residuo variazionale indipendente; CSV legacy non riproducibili.
8. **Nessuna capacità inerziale superficiale** nei modelli iperbolici implementati
   (limite di modello, non di validazione): i risultati non si estendono a superfici di
   2° ordine.
9. (aggiuntivo) **Singola combinazione di parametri** per l'accoppiato (ρ=1, ζ=0, c_Γ=1,
   β=1, IC e T dati): nessuna generalizzazione ad altri parametri/IC; VTU a ~8 cifre
   significative; nessun run OOFEM a n64×32 (floor a h=0.03125 sconosciuta); i due bachi
   di codice noti (§B.3) restano aperti (non esercitati dai benchmark).

---

## J. Conclusione e raccomandazioni

### J.1 Verdetti separati

| # | Oggetto | Verdetto | Base |
|---|---|---|---|
| 1 | Correttezza dell'implementazione/assemblaggio | **SUPPORTATA** | ispezione termine per termine (B.1); formula esatta verificata (D.1); coerenza OOFEM post-fix (D.2); ordini ≈2 nei casi regolari (C1–C3, C9); invarianza linee (diff 0); simmetria 180° (2.2e-16). Limite: nessun residuo variazionale indipendente sul sistema nonlineare accoppiato |
| 2 | Riproducibilità | **SUPPORTATA** | D.3 (bit-identica su 3 rerun OOFEM + T1-A + T1-B + riferimento DOP853; 16/16 metriche diff 0). Limite: ambiente D.5 / I-5 |
| 3 | Accordo col riferimento semidiscreto | **SUPPORTATA** (nel caso testato) | C9: ordini 1.972/1.948/1.975/2.027; floor ~1e-8; ΔE 6.5e-6 ∝ O(dt²). Limite: combinazione di parametri singola; riferimento numerico, non soluzione esatta |
| 4 | Convergenza spaziale (esistenza) del caso accoppiato | **SUPPORTATA** (fino a h=0.03125; famiglia fino a h=0.015625) | E.9: FE RMS monotono, nessun plateau, convergenza puntuale, componente OOFEM alla floor a n32, stabilità del proxy ≤8% |
| 5 | Ordine asintotico spaziale | **NON DIMOSTRATA** | E.4/E.5: ordini 0.72/1.05 (OOFEM) non costanti; rapporto famiglia 3.09 = singolo rapporto di differenze, non costante nel tempo; bulk max non monotono; **non è affermato un ordine 2 né un ordine 3** |
| 6 | Singolarità angolare | **NON DEMOSTRATA** | F: fatti coerenti, interpretazione non risolta; p = −0.016 su 4 punti; fit h²+h ritirato |
| — | **Validazione complessiva** | **SUPPORTATA PARZIALMENTE** | bulk iperbolico forte (C1–C3); parabolico 1D favorevole con riserve (C4); 2D parabolico qualitativo/partial (C5–C7); accoppiato: tempo SUPPORTATA, spazio SUPPORTATA (esistenza) con ordine e singolarità NON DIMOSTRATI; energia coerente ma non indipendente (G) |

### J.2 Affermazioni non ammesse con le evidenze attuali

- "il solver è validato integralmente" / "il modello Allen–Cahn è validato in generale";
- "il caso accoppiato converge con ordine 2" (o con qualsiasi ordine asintotico);
- "è dimostrata una singolarità angolare" o "gli ordini ridotti sono causati da essa";
- "l'energia valida il solver" (monotonia ≠ validazione completa);
- "l'accordo col DOP853 è accordo con la soluzione esatta della PDE";
- estensioni ad altri parametri, IC, geometrie o a superficie inerziale.

### J.3 Attività future raccomandate (massimo 3, in ordine di valore informativo/costo;
**nessuna eseguita**)

1. **Run OOFEM n64×32 del caso accoppiato** (stesso input, dt, T; ~40–90 s, 1 core,
   ~100–200 MB, directory nuova): fornisce (i) la floor OOFEM a h = 0.03125
   (OOFEM-n64 vs DOP853-n64 a mesh/tempo uguali) — verifica che l'accordo ~1e-8 persista
   al livello più fine, rafforzando il verdetto 3; (ii) un quarto punto per-nodo del
   massimo mobile (r_c = 0.0442 sul lato OOFEM). Non aggiunge un ordine di convergenza
   (per quello serve la differenza di famiglia, già in T1-B).
2. **Export della velocità `waveVelocity` nei VTU** (modifica di codice mirata e non
   numerica, **richiede approvazione esplicita**): ricalcolo del bilancio iperbolico 16×8
   con la velocità del solver; accettazione proposta: R(T) sotto il pavimento di
   ricostruzione (~7e-7) e 0 incrementi positivi. Trasforma l'energia da coerenza
   interna a verifica indipendente del ramo Newmark.
3. **Test bulk-only** (ρu_tt − kΔu = m(u−u³) senza elementi di superficie, n8/n16/n32 +
   riferimento DOP853 bulk-only): se il bulk-only raggiunge ~2 e l'accoppiato no, la
   riduzione dell'ordine è attribubile all'accoppiamento/comportamento angolare; se anche
   il bulk-only riduce, il limite è nel setup/discretizzazione bulk. **Rischio:** la
   fattibilità dell'input (modello senza linee) va prima verificata staticamente; in caso
   contrario servirebbero modifiche a modello/formato, fuori vincolistica senza errore
   specifico documentato.

---

## Versione, data e controlli di coerenza

- **Data:** 2026-10-27; **Versione:** 1.0; **Autore:** campagna di validazione
  Allen–Cahn OOFEM (post T1-B).
- **Controlli eseguiti per questo rapporto:** (i) ogni ordine della tabella E.4
  ricalcolato come log2 del rapporto degli errori adiacenti — tutti coerenti entro
  l'arrotondamento a 4 cifre; (ii) i tempi di confronto sono tutti a tempo uguale
  (convenzione verificata, Δt ≤ 3.5e-18); (iii) le definizioni di norma (E.2) coincidono
  con gli script `t1b_validate.py` / `check_jacobian_fix.py`; (iv) i valori della sezione
  E trascritti dai CSV di T1-B/T1-A citati, senza modifiche; (v) le discrepanze tra
  documenti precedenti e dati più recenti sono segnalate in D.5 (in particolare: stato
  "FAIL" accoppiato → superato; fit h²+h → ritirato; 8449 → 8385 DOF; √2 sulle metriche
  "prima").

### Lista degli artefatti sorgente principali

| Argomento | Percorso (relativo a `doc/stvowu_validation/`) |
|---|---|
| Report parabolico principale (26-10, parti superate marcate) | `OOFEM_AllenCahn_Validation_Report.md` |
| Audit iniziale 27-10 (inventario, lacune, piano T1–T5, criteri pre-registrati) | `audit_iniziale_2026-10-27.md` |
| Audit metriche quantitative (Lee, Okumura, Calatroni, energia) | `validation_metrics_audit.md` |
| Audit di riproducibilità (stati per test) | `reproducibility_audit.md` |
| Rerun operative (3 run bit-identici) | `reproducibility_runs.md` |
| Matrice "finale" (stato superato, conservata) | `validation_matrix_final.md` / `.csv` |
| Matrice lavoro (criteri proposti) | `validation_matrix.md` |
| Inventario artefatti | `validation_inventory.md` |
| Audit energetico finale | `energy_validation_final.md`; teoria: `energy_theory.md` |
| Parabolico 2D (Xiao–Chen, raffinamenti) | `parabolic_2d/parabolic_2d_validation.md` |
| Confronti di letteratura (Okumura, Xiao–Chen, Calatroni–Colli CSV) | `literature_comparison/README.md` |
| Bulk iperbolico (Nizovtseva, kink, modi) | `hyperbolic/hyperbolic_bulk_validation.md`, `hyperbolic/README.md`, `hyperbolic/nizovtseva_pde_residual_check.{py,csv}` |
| Discrepanza accoppiata (audit termini) | `hyperbolic/hyperbolic_coupled_discrepancy_audit.md` |
| T1: verifica spaziale + audit statico + erratum | `hyperbolic/coupled_independent_spatial_jacobian_fix/{spatial_verification_report.md, static_audit_report.md}` |
| Correzione Jacobiano del riferimento | `hyperbolic/coupled_independent_jacobian_fix/jacobian_fix_verification_report.md` (directory di lavoro `coupled_independent_spatial_jacobian_fix/check_jacobian_fix.py`) |
| T1-A: rerun n32×16 | `hyperbolic/t1a_n32x16_rerun/results/t1a_report.md` (+ 7 CSV) |
| T1-D: audit analitico | `hyperbolic/t1d_analytic_audit/t1d_report.md` (+ 7 CSV) |
| T1-B: proxy n128×64, FASE 0–3 | `hyperbolic/t1b_dop853_n128x64/t1b_report.md`, `results/t1b_*.csv`, `run.log`, `resource_report.txt`, `reference_n128x64_dt0p00025.npz`, `matrices_n128x64.npz` |
| Questo rapporto e file tabellari | `consolidated_final_report/{validation_report_consolidated.md, validation_matrix.csv, coupled_spatial_convergence.csv, README.md}` |

**Fine del rapporto consolidato.**
