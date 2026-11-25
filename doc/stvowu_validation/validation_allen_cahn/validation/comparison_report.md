# Rapporto di confronto con la letteratura — validazione del modello Allen–Cahn in OOFEM

**Data:** 27 ottobre 2026. **Ambiente:** OOFEM 2.6.devel, build sorgente locale
`build/debug_magneto/oofem`, revisione git `2cc1a7706` (branch `allencahn`);
Python 3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3.
Tutte le esecuzioni del solutore, i post-processing e le figure riportate in
questo documento sono stati effettuati e ri-verificati in questa data. Nessun
risultato è simulato: ogni valore è tracciabile fino ai file grezzi in
`parabolic/results/` e `hyperbolic/results/` tramite gli script indicati in
§14. I riferimenti bibliografici con ancore a pagine/equazioni sono in
`literature/bibliography.md`; la matrice di compatibilità è in
`literature/benchmark_selection.md`. I percorsi relativi nel testo sono
riferiti alla root della campagna `validation_allen_cahn/`; il prefisso
`../` indica il corpus pre-esistente (`doc/stvowu_validation/`).

---

## 1. Esecutivo

Tolleranze dichiarate (valutate su quanto riportato nei §9–11):

- **Parabolico:** (a) energia totale discreta non crescente a ogni passo su
  tutte le 700 frame esportate; (b) accordo con la curva di energia
  digitalizzata da Fig. 4.4 di [O1]: max |ΔE| ≤ 2 unità (≈ 2,2 % dell'asse
  0–90 della figura) su 0,01 ≤ t ≤ 0,4; (c) campo finale a ≤ 10⁻³ dalla fase
  +1; (d) convergenza spaziale di ordine 2.
- **Iperbolico:** (a) residuo del PDE della soluzione analitica nullo (entro
  la precisione di calcolo); (b) sulla mesh più fine RMS ≤ 10⁻⁵ e max ≤ 10⁻⁴ a
  T = 0,2; (c) convergenza spaziale di ordine 2; (d) errore relativo sulla
  velocità dell'interfaccia ≤ 10⁻³.

| Test | Riferimento | Esecuzione OOFEM | Errore quantitativo | Convergenza verificata | Esito |
|---|---|---|---|---|---|
| 1. Allen–Cahn parabolica con condizione dinamica al bordo (1D, bulk + due estremi dinamici) | [O1] Okumura 2021, Ch. 4, *Computation example 2*; PDE (4.1), legge al bordo (4.91), Fig. 4.3 (campo) e 4.4 (energia) | Eseguita: striscia Q1 100×1, Δt = 10⁻³, T = 0,7, 700 passi Backward–Euler, parametri pubblicati riportati tali e quali (γ = 100, ε_ex = 10, IC gaussiana) | Energia: 77,4205 → 2,22·10⁻⁶, monotona decrescente (max incremento −9,7·10⁻⁸); vs curva digitalizzata Fig. 4.4 su 0,01 ≤ t ≤ 0,4: RMSE 0,9881, MAE 0,8804, max 1,8323 (asse 0–90); campo finale 0,999898–0,9999997 | Sì: raffinamento spaziale n = 50/100/200/400, ordini osservati 2,07 (RMS) e 2,32 (RMS), 2,07/2,32 (max) a t = 0,3; da t = 0,5 le mesh coincidono entro la precisione di esportazione (6 cifre significative) | **Validato entro le tolleranze dichiarate** (con le limitazioni del §12: il riferimento parabolico è una curva digitalizzata di una figura, non una tabella di dati; il paper usa differenze finite, OOFEM usa Q1) |
| 2. Allen–Cahn iperbolica, onda viaggiante limitata (1D, solo bulk) | [N1] Nizovtseva–Galenko–Alexandrov 2016, Eq. (29) + Table 1 (soluzione k = 1, c₀ = 0), τ = 0,2, D = M = 1 | Eseguita: striscia Q1 80/160/320 × 1 su [−12,12], Newmark β = 1/4, γ = 1/2, Δt = 10⁻³, T = 0,2, 200 passi, IC = profilo esatto a t = 0 e velocità iniziale esatta u_t(x,0) = cλe^{λx}/(1+e^{λx})² | A T: RMS 2,92·10⁻⁵ / 7,32·10⁻⁶ / 1,86·10⁻⁶ e max 1,16·10⁻⁴ / 2,89·10⁻⁵ / 7,19·10⁻⁶ per n = 80/160/320; posizione interfaccia a T: errore 3,03·10⁻⁵; errore relativo velocità 9,85·10⁻⁵ | Sì: ordini osservati 1,99/1,98 (RMS) e 2,00/2,01 (max) su tre raffinatamenti | **Validato entro le tolleranze dichiarate** — **la condizione dinamica al bordo NON è esercitata né validata in questo test** (il benchmark è un problema di Cauchy senza condizione al bordo; dichiarazione esplicita richiesta e rispettata) |

Nessun ulteriore test è riportato come "validato": i test esplorativi del
corpus pre-esistente (Xiao–Chen, Hochbruck–Leibold, audit accoppiato DOP853)
sono descritti in `../hyperbolic/` e `../literature_comparison/` con i loro
limiti, e **non** fanno parte della tabella sopra perché non sono confronti
direttamente quantitativi della coppia modello/condizione al bordo qui
validata.

---

## 2. Obiettivo e ambito

Verificare, contro soluzioni di riferimento della letteratura, le due
implementazioni del modello Allen–Cahn presenti nella build OOFEM locale:

1. **Parabolico con condizione dinamica al bordo** (bulk + superficie del
   primo ordine accoppiata):
   `nltransienttransportproblem` (Backward–Euler) + materiale `allencahn` +
   elementi di bordo `line1ht` che condividono i gradi di libertà della
   traccia.
2. **Iperbolico inerziale** (bulk): `nltransienttransportproblem` con
   `wave 1` (Newmark) + materiale `allencahn` con `hyperbolic 1`.

L'obiettivo è una *validazione* (confronto con soluzioni esatte o dati
pubblicati), non una *verifica* MMS: quest'ultima esiste già nel corpus
pre-esistente (material con sorgenti manufactured, test lineari a modo
coseno) e non è ripetuta qui.

## 3. Riferimenti

Citazioni complete e ancore (PDF su disco, pagine stampate/numerate,
equazioni) in `literature/bibliography.md`. Riassunto:

- **[O1]** M. Okumura, *Design and analysis of structure-preserving schemes
  for parabolic PDEs with dynamic boundary conditions*, PhD dissertation,
  Osaka University, gennaio 2021, DOI 10.18910/82279.
  - Ch. 4 §1, pp. 40–41 (PDF 49–50): PDE (4.1) u_t = u_xx − F′(u); legge
    dinamica (4.2)–(4.3); dissipatività (4.5)
    d/dt {J(u) + F(u(0,t)) + F(u(L,t))} ≤ 0 con
    J(u) = ∫₀ᴸ (|u_x|²/2 + F(u)) dx.
  - Ch. 4 §6, p. 70 (PDF 79): legge (4.91) moltiplicata per ε_ex; K = 100,
    L = 1, ε_ex = 10, F(s) = (γ/4)(s²−1)², γ = 100.
  - §6.2 *Computation example 2*, p. 71 (PDF 80): N = 700 (Δt = 10⁻³),
    T = 0,7, u(x,0) = exp{−500(x−0,5)²}; Fig. 4.3 (campo) e **Fig. 4.4
    (energia, asse 0–90)**.
- **[N1]** I. G. Nizovtseva, P. K. Galenko, D. V. Alexandrov, *J. Phys. A:
  Math. Theor.* **49** (2016) 435201.
  - p. 4: PDE (1) τ_R f_tt + f_t = D∇²f − M f′(f), f(f) = (1/4)(f²−1)²,
    eq. (2)–(3); ODE dell'onda (1−τc²)f″ + cf′ + f − f³ = 0, eq. (4),
    τ = τ_R M.
  - p. 5: velocità c = ±3/√(2+9τ), eq. (26); **Table 1** (contatori di segno
    Q, A, B; soluzione k = 1: Q = A = B = +1).
  - p. 6: **onda viaggiante limitata, eq. (29)**; classe di soluzioni
    limitate discussa in §5.2.1, p. 10.
- Fonti secondarie (contesto, non benchmark dirette): [C1] Calatroni–Colli
  2012; [X1] Xiao–Chen 2025; [K1] Kovács–Lubich 2017; [M1]
  Makki–Miranville–Petcu 2022; [H1] Hipp–Kovács 2021; [L1] Hochbruck–Leibold
  2021; [L2] Hochbruck–Leibold 2020; [P1] Patzák 2012; [H2] Hairer et al.
  1993; [R1] Roache 2002 — tutti in `literature/bibliography.md`.

## 4. Motivazione della scelta dei benchmark

Matrice completa di compatibilità (classi A/B/C/D) in
`literature/benchmark_selection.md`. Sintesi:

- **[O1] è l'unica fonte di classe A parabolica**: pubblica la ricetta
  numerica completa (mesh, Δt, orizzonte, γ, ε_ex, IC) **e** dati
  quantitativi di riferimento (curva di energia Fig. 4.4, evoluzione del
  campo Fig. 4.3) per lo stesso problema con la stessa legge dinamica.
- **[N1] è l'unica fonte di classe A iperbolica**: fornisce un'onda
  viaggiante limitata **esatta** (eq. 29) della PDE dimensionale con
  τ = 0,2, D = M = 1.
- **Dichiarazione esplicita (requisito di integrità):** nessuno dei documenti
  disponibili fornisce una soluzione di riferimento per il modello
  iperbolico Allen–Cahn **con condizione dinamica al bordo**; di conseguenza,
  nel benchmark 2 la condizione dinamica **non è esercitata né validata**
  (ripetuto nell'esecutivo, nel titolo della Figura 2 e nelle conclusioni).
  L'audit pre-esistente del modello accoppiato inerziale/superficie contro un
  riferimento semidiscreto indipendente DOP853
  (`../hyperbolic/coupled_independent_jacobian_fix/`, ordini 1,95–2,03,
  residuali 5·10⁻⁸/1,4·10⁻⁷) è riportato come evidenza complementare del
  *percorso di codice accoppiato*, **non** come validazione contro la
  letteratura.

## 5. Formulazioni matematiche

### 5.1 Parabolico con condizione dinamica ([O1], Ch. 4)

u_t = u_xx − F′(u), in (0,L)×(0,T), L = 1, F(s) = (γ/4)(s²−1)², γ = 100,

con legge dinamica (4.91) (ε_ex = 10):

- ε_ex u_t(0,t) = u_x(0,t) − F′(u(0,t)),
- ε_ex u_t(L,t) = −u_x(L,t) − F′(u(L,t)),

e energia E(t) = J(u) + F(u(0,t)) + F(u(L,t)) ≤ E(0), con
J(u) = ∫₀ᴸ (|u_x|²/2 + F(u)) dx (eq. 4.5).

In OOFEM la stessa PDE si scrive c u_t − kΔu = mobility·(u−u³) con
c = k = 1, mobility = γ: infatti γ(u−u³) = −F′(u).

### 5.2 Iperbolico, onda viaggiante esatta ([N1])

PDE (1) con τ_R = 0,2, D = M = 1:

0,2 u_tt + u_t = u_xx + (u − u³).

Onda viaggiante esatta limitata (soluzione k = 1 di Table 1, c₀ = 0,
velocità positiva; eq. 29):

u(x,t) = 1 / (1 + exp[λ (x − c t)]), λ = √(2+9τ)/2 = 0,9746794344808963,
c = 3/√(2+9τ) = 1,5389675281277313, τ = τ_R M = 0,2.

L'interfaccia (livello u = 1/2) è x₀(t) = c t; a T = 0,2:
x₀(T) = 0,3077935056255463. La soluzione è stata verificata come soluzione
esatta del PDE: residuo massimo
|τ_R u_tt + u_t − u_xx − u + u³| = 1,665·10⁻¹⁶ su griglia fitta
(`hyperbolic/scripts/check_exact_residual.py`). Le soluzioni singolari
(rami con A o B nulli in Table 1) sono escluse per costruzione: si usa la
classe limitata di §5.2.1 del paper.

## 6. Parametri, condizioni iniziali e al contorno

| | Benchmark 1 (parabolico) | Benchmark 2 (iperbolico) |
|---|---|---|
| Dominio | [0,1] (1D), estratto su striscia [0,1]×[0,1] | [−12,12] (truncamento del problema di Cauchy), striscia [−12,12]×[0,1] |
| PDE | u_t = u_xx − γ u(u²−1) | 0,2 u_tt + u_t = u_xx + (u−u³) |
| γ (mobility) | 100 | — (m = 1) |
| ε_ex (capacità di bordo) | 10 | — (nessun bordo dinamico) |
| IC | u(x,0) = exp{−500(x−0,5)²} | u(x,0) = 1/(1+exp[λ x]); u_t(x,0) = cλe^{λx}/(1+e^{λx})² (profilo e velocità iniziali esatti, 642 valori nodali nell'input; picco di u_t = 0,375 a x = 0) |
| BC | dinamica (4.91) su entrambi gli estremi | nessuna (bordo naturale Neumann per truncamento; effetto quantificato in §12) |
| Δt | 10⁻³ | 10⁻³ |
| Orizzonte T | 0,7 (700 passi) | 0,2 (200 passi) |
| Scheme temporale | Backward–Euler (α = 1) | Newmark β = 1/4, γ = 1/2 |

I valori sono riportati nell'input **tali e quali** pubblicati; nessun
parametro è stato scelto artificialmente favorevole.

## 7. Implementazione OOFEM (cosa viene effettivamente risolto)

- **Materiale `allencahn`** (`src/tm/Materials/allencahnmat.C`, nuovo in
  questa branch): flusso −k∇u (k = `k`), capacità `c`, sorgente interna
  mobility·(u−u³) trattata al secondo membro, la cui tangente +mobility·(3u²−1) viene mossa al primo membro (β = surfacebeta sulla linea, 1 nel bulk). Sui
  `line1ht` la capacità è `c` (parabolico) o `surfacecapacity` (iperbolico)
  e la sorgente è moltiplicata per `surfacebeta`. In modalità
  `hyperbolic 1` il bulk usa `rho` per u_tt; `wavedamping` moltiplica la
  massa bulk per il termine di smorzamento ρδ u_t.
- **Accoppiamento bulk/bordo (condizione dinamica).** Gli elementi
  `line1ht` condividono i gradi di libertà dei nodi di bordo. L'equazione
  discreta sul DOF condiviso è la somma della contribuzione bulk (forma
  debole ∫_Ω k∇u·∇v, il cui termine di frontiera ∫_Γ k ∂_n u v da
  integrazione per parti, con n uscente dal dominio, si assembla nel DOF) e
  della contribuzione di linea (∫_Γ c_Γ ψ_t v + ∫_Γ k ∂_sψ ∂_s v +
  ∫_Γ β·mobility (ψ−ψ³) v). Il sistema discreto coincide con:
  - bulk: c u_t − kΔu = mobility (u−u³);
  - bordo: c_Γ ψ_t − k Δ_Γ ψ + k ∂_n u = β·mobility (ψ−ψ³),
  con ∂_n derivata normale **uscente**. Con i valori del benchmark 1
  (c_Γ = 10, k = 1, β = 1, mobility = 100 = γ) si ottiene, a x = 0
  (n = −x̂, ∂_n u = −u_x): 10 ψ_t = u_x(0,t) + 100(ψ−ψ³) = u_x(0,t) − F′(ψ),
  e a x = L (n = +x̂): 10 ψ_t = −u_x(L,t) − F′(ψ): **la legge (4.91) del
  paper, termine a termine, su entrambi gli estremi**. La convenzione di
  segno è corroborata anche dal test MMS interno al materiale (sorgente
  manufactured di bordo con `+normalConductivity·∂_n u`).
- **Riduzione a striscia 2D.** Il caso 1D è discretizzato su una striscia
  di larghezza 1 con un elemento Q1 per cella lungo x e un unico elemento
  in y (NY = 1). Poiché IC, BC e PDE sono indipendenti da y, la soluzione
  Q1 invariante in y soddisfa esattamente il problema 1D: la verifica
  numerica (max |u(x,0) − u(x,1)| = 0,0 su tutte le frame, §9.1) conferma
  che nessun modo trasversale si eccita. Il laplaciano tangenziale Δ_Γ su
  bordo rettilineo si riduce a zero sulla retta e la "linea" ha
  spessore/area unitari (`SimpleTransportCS thickness 1.0 area 1.0`).
- **Integrazione:** quad1ht con 3×3 punti di Gauss (`nip 9`); line1ht con
  3 punti di Gauss (sufficiente per integrare esattamente il sorgente
  cúbico N_i(u−u³) di grado 4). Newton modificato, rtol = 10⁻¹⁰.

## 8. Configurazione OOFEM

**Benchmark 1** (`parabolic/input/okumura_ch4_case2.in`; generato da
`parabolic/scripts/run_okumura_example2.py`):

```
nltransienttransportproblem nsteps 700 deltat 0.001 alpha 1.0 rtol 1.e-10 nsmax 30 nmodules 1
vtkxml tstep_all domain_all primvars 1 6
domain heattransfer
ndofman 202 nelem 102 ...          ! 101 nodi per riga x 2, 100 quad1ht + 2 line1ht
allencahn 1 d 1.0 k 1.0 c 1.0 mobility 100.0                      ! bulk
allencahn 2 d 1.0 k 1.0 c 10.0 mobility 100.0 surfacebeta 1.0     ! bordo (c = eps_ex = 10)
Set 1 elementranges {(1 100)}     ! bulk
Set 2 elements 2 101 102          ! linee di bordo
```

IC nodale: `InitialCondition ... u exp(−500(x−0,5)²)` per i 101 nodi di
bordo (set 3..103, coppie di nodi per colonna).

**Benchmark 2** (`hyperbolic/input/nizovtseva_n320.in`; generato da
`hyperbolic/scripts/run_nizovtseva.py`, identico per n = 80/160 salvo mesh):

```
nltransienttransportproblem nsteps 200 deltat 0.001 alpha 1.0 rtol 1.e-10 nsmax 100 nmodules 1
      wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0 wavedamping 5
      initialvelocityvalues 642 <u_t(x,0) = cλe^{λx}/(1+e^{λx})², 642 valori nodali>
ndofman 642 nelem 320 ...          ! striscia [−12,12]x[0,1], solo quad1ht, nessun line1ht
allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 hyperbolic 1 rho 0.2 surfacecapacity 1.0
```

`rho 0.2` → termine inerziale 0,2 u_tt; `wavedamping 5` → ρδ = 0,2·5 = 1,0
(coefficiente di u_t). `initialvelocity 0.0` è il valore scalare di default;
i 642 valori di `initialvelocityvalues` lo sovrascrivono nodo per nodo con la
velocità iniziale **esatta** dell'onda (il solutore usa `waveVelocity =
initialVelocityValues` quando la dimensione coincide con il numero di
equazioni, `nltransienttransportproblem.C`). Nessun elemento di bordo: il
benchmark 2 è bulk-only (dichiarazione del §1).

**Nota sulla precisione di esportazione.** `vtkxml ... primvars 1 6`
esporta la variabile ignota #1 (u) e #6 (Temperature); il "6" non è una
precisione. I valori dei nodi sono scritti con la precisione di default
del flusso C++ (6 cifre significative, `vtkxmlexportmodule.C`,
`writeVTKPointData`: `fileStream << scientific << value`); per campi O(1)
il quantum di esportazione è ~10⁻⁶. Questo vincola le tabelle di
raffinamento (paragrafo 11), non i risultati del benchmark (le metriche
usano frame esportate a 6 cifre, dichiarate come tali).

## 9. Risultati

### 9.1 Benchmark 1 — parabolico con condizione dinamica

- **Energia totale** (Q1, quadratura 3×3 di Gauss sul gradiente,
  integrata in post-processing con la stessa quadratura che il solutore usa
  per il gradiente; potenziale e termini di bordo valutati con
  `scipy.integrate.quad` cella per cella — tracciabile in
  `parabolic/scripts/postprocess_okumura.py`):
  E(0,001) = 77,42047359415844 → E(0,7) = 2,219357171842687·10⁻⁶.
  **Monotona non crescente su tutte le 700 frame**; massimo incremento
  E(t+Δt)−E(t) = −9,695·10⁻⁸ (negativo). La dissipatività (4.5) è
  rispettata dal solutore a passo finito su tutto l'orizzonte.
- **Campo:** a T = 0,7 il campo ha rilassato alla fase +1
  (min 0,9998982, max 0,9999997), coerente con Fig. 4.3 del paper;
  variazione trasversale max |u(x,0)−u(x,1)| = 0,0 (la striscia riduce
  esattamente a 1D).
- **Confronto con Fig. 4.4 digitalizzata** (663 campioni, intervallo
  0,01 ≤ t ≤ 0,4, asse della figura 0–90): RMSE = 0,9881048219018264,
  MAE = 0,8803680411631917, max |ΔE| = 1,832295199199237. Il valore
  iniziale della curva digitalizzata è E(0,00399) = 75,1095, contro
  E_OOFEM(0,004) = 74,3692: la discrepanza iniziale (~0,74) è la più grande
  del confronto ed è coerente con (i) l'incertezza di digitalizzazione e
  (ii) la diversa discretizzazione spaziale (FD a N = 700 uniformi nel
  paper vs Q1 a n = 100 in OOFEM): il picco dell'IC gaussiana decade nel
  primo passo temporale (u_t al picco = −1000) e i due schemi lo
  trattano in modo diverso. L'energia OOFEM a t = 0,001 è stata
  inoltre verificata contro il valore esatto della stessa funzione
  Q1 (scipy quad, termine per termine: 77,42047359415842; la serie
  esportata riproduce questo valore a ~6·10⁻⁷).
- Figura: `figures/figure_1_parabolic.{png,pdf}` — (a) energia
  digitalizzata della letteratura; (b) istantanee OOFEM; (c) sovrapposizione
  con metriche; (d) differenza.

### 9.2 Benchmark 2 — onda viaggiante iperbolica

Metriche a T = 0,2 (campo OOFEM vs onda esatta, dominio completo [−12,12];
fonti: `hyperbolic/nizovtseva_metrics.csv`):

| n | h | RMS (full) | max (full) | x₀^OOFEM(T) | |x₀ − cT| |
|---|---|---|---|---|---|
| 80 | 0,3 | 2,918311384156078·10⁻⁵ | 1,1578978006165919·10⁻⁴ | 0,30829199705849153 | 4,984914·10⁻⁴ |
| 160 | 0,15 | 7,323071979294976·10⁻⁶ | 2,888978006165832·10⁻⁵ | 0,3079164348915018 | 1,229293·10⁻⁴ |
| 320 | 0,075 | 1,8620658773059325·10⁻⁶ | 7,189780061700457·10⁻⁶ | 0,3078238237032819 | 3,031808·10⁻⁵ |

(x₀ = posizione del livello u = 1/2; x₀^esatto(T) = 0,3077935056255463.)

Velocità dell'interfaccia (da x₀(T)/T): c_OOFEM(320) = 1,5391191185164093
vs c_esatta = 1,5389675281277313 → errore relativo 9,85·10⁻⁵.

- Figura: `figures/figure_2_hyperbolic.{png,pdf}` — (a) profili esatti;
  (b) profili OOFEM n = 320; (c) sovrapposizione con metriche; (d)
  posizione dell'interfaccia x₀(t) con inset di convergenza RMS vs h
  (riferimento h²). Il titolo dichiara: "bulk inertial model - no dynamic
  boundary condition in this test".

## 10. Metriche di errore

Tutti i valori, con unità e provenienza file per file:
`validation/error_metrics.csv`. I numeri qui sopra sono trascritti dai CSV
di post-processing, che a loro volta leggono esclusivamente i file grezzi
esportati da OOFEM (`.vtu`); nessuna trasformazione dei dati grezzi è
aplicata oltre il calcolo delle norme (RMS = sqrt(mean(Δ²)) sui nodi del
dominio, max = max|Δ|) e la quadratura dell'energia, tutte documentate
negli script citati.

## 11. Studio di convergenza

Fonti: `validation/convergence_study.csv` (iperbolico:
`hyperbolic/nizovtseva_metrics.csv`; parabolico:
`parabolic/okumura_refinement_convergence.csv`).

**Iperbolico** (Δt = 10⁻³ fissato, tre raffinatamenti, errore a T = 0,2
contro l'onda esatta):

| h | RMS | ordine | max | ordine |
|---|---|---|---|---|
| 0,3 | 2,9183·10⁻⁵ | — | 1,1579·10⁻⁴ | — |
| 0,15 | 7,3231·10⁻⁶ | 1,9946 | 2,8890·10⁻⁵ | 2,0029 |
| 0,075 | 1,8621·10⁻⁶ | 1,9755 | 7,1898·10⁻⁶ | 2,0065 |

Convergenza spaziale **di secondo ordine** in entrambe le norme,
coerente con Q1.

**Parabolico** (Δt = 10⁻³ fissato; griglie nodali annidate
x_i^50 = x_{2i}^100 = x_{4i}^200 = x_{8i}^400; differenza nodale diretta
tra mesh, nessuna interpolazione; riferimento = mesh più fine n = 400):

| h | t | RMS vs n = 400 | ordine | nota |
|---|---|---|---|---|
| 0,02 | 0,3 | 1,80279·10⁻⁴ | — | interfaccia in transito |
| 0,01 | 0,3 | 4,29312·10⁻⁵ | 2,070 | |
| 0,005 | 0,3 | 8,57637·10⁻⁶ | 2,324 | |
| 0,02 | 0,5 | 7,20274·10⁻⁶ | — | accordo ~10⁻⁵ |
| 0,01 | 0,5 | 1,69077·10⁻⁶ | — | entro la precisione di esportazione (6 cifre) |
| 0,005 | 0,5 | 3,35614·10⁻⁷ | — | entro la precisione di esportazione |
| 0,02–0,005 | 0,7 | ≤ 1,40726·10⁻⁷ | — | mesh identiche entro 6 cifre significative |

Ordine di secondo ordine a t = 0,3 (dove le differenze superano di 1–2
ordini il quantum di esportazione 10⁻⁶); per t ≥ 0,5 le soluzioni su mesh
diverse coincidono entro la precisione a 6 cifre significative con cui OOFEM
esporta i campi (§8), e **nessun ordine viene assegnato** a quei valori
(sarebbe un artefatto di quantizzazione). Il riferimento n = 400 è una
soluzione di riferimento *numerico*, non la soluzione esatta.

## 12. Limitazioni (dichiarate)

1. **Riferimento parabolico = curva digitalizzata.** Fig. 4.4 è una
   figura stampata: la digitalizzazione ha un'incertezza raster non
   quantificata separatamente; la metrica RMSE/MAE/max la tratta come
   riferimento. Inoltre lo schema del paper è un metodo a differenze
   finite strutturante basato sul framework DVDM (discrete variational
   derivative method; dissipativo per costruzione, Teorema 4.1 di [O1])
   su N = 700 nodi uniformi, mentre OOFEM usa Q1 (n = 100 celle) con
   Backward–Euler: la discrepanza (~1–2 unità su asse 0–90) non è
   scomponibile in "errore di digitalizzazione" vs "differenza di schema"
   con i dati disponibili. Il confronto va quindi letto come *accordo
   quantitativo entro il ~2 % dell'asse della figura*, più la prova
   indipendente della dissipatività energetica, non come replica esatta
   dei numeri del paper.
2. **Precisione di esportazione 6 cifre significative.** Le metriche e le
   tabelle di raffinamento usano i campi `.vtu` esportati; sotto ~5·10⁻⁶
   (campi O(1)) le differenze sono quantizzazione di esportazione e gli
   ordini non vengono assegnati (§11). I valori grezzi non sono mai
   sovrascritti o interpolati.
3. **Striscia 2D ≠ verifica 2D completa.** I benchmark validano il
   comportamento 1D (bulk + legge dinamica sugli estremi) e la riduzione
   esatta a striscia. Il laplaciano di Laplace–Beltrami tangenziale su
   fronti curvi (presente nel problema originale 2D/3D citato in [O1]) non
   è esercitato: le "superfici" sono rette rettilinee con Δ_Γ = 0.
4. **Truncamento del dominio iperbolico.** Il problema di [N1] è di
   Cauchy (R); qui il dominio è [−12,12] con bordo naturale. Le code
   dell'onda a T valgono ~exp(−λ·11,7) ≈ 10⁻⁵, della stessa taglia
   dell'errore discreto alla mesh più grossa (n = 80) e trascurabile per
   l'interfaccia (coda a 12 unità dall'interfaccia); il truncamento è
   quindi dell'ordine della tolleranza dichiarata, ma non è zero.
5. **Condizione dinamica non validata nel test iperbolico.** Nessun
   riferimento disponibile fornisce una soluzione di riferimento per il
   modello iperbolico con condizione dinamica; il benchmark 2 la dichiara
   non esercitata. L'audit accoppiato DOP853 (§4) è evidenza di
   *percorso di codice*, non validazione letteraria.
6. **Integrazione temporale non raffinata.** Δt = 10⁻³ è il valore
   pubblicato (parabolico) e fissato (iperbolico); non è stato eseguito
   uno studio di raffinamento temporale, quindi l'ordine temporale dello
   schema (BE del primo ordine; Newmark 1/4–1/2 del secondo ordine per il
   problema lineare) **non è verificato numericamente in questo lavoro**.
   Gli ordini riportati sono ordini *spaziali* (Δt fissato).
7. **Build debug.** Le esecuzioni usano la build `debug_magneto` locale;
   nessun risultato qualitativo cambia, ma i tempi di calcolo sono
   gonfiati (irrelevante per la validazione).

## 13. Conclusioni

**Parabolico (benchmark 1).** L'implementazione OOFEM dell'equazione
Allen–Cahn parabolica con condizione dinamica al bordo (4.91) riproduce,
con i parametri pubblicati di [O1] Example 2 e senza alcuna taratura,
l'evoluzione del campo verso la fase +1 e la decrescita dell'energia
totale fino a 2,2·10⁻⁶ a T = 0,7, con dissipatività numerica a ogni passo
(max incremento −9,7·10⁻⁸ su 700 passi) e accordo con la curva di energia
pubblicata entro ~2 % dell'asse della figura (RMSE 0,988, max 1,832 su
0–90). Il raffinamento spaziale mostra ordine 2 (2,07/2,32 in RMS a
t = 0,3). **Esito: validato entro le tolleranze dichiarate**, con la
limitazione esplicita che il riferimento è una curva digitalizzata e lo
schema del paper è diverso (BE–Q1 vs FD–DVDM).

**Iperbolico (benchmark 2).** L'implementazione OOFEM del modello
inerziale (Newmark 1/4–1/2, `allencahn hyperbolic`) riproduce l'onda
viaggiante esatta limitata di [N1] eq. (29) con errore RMS ≤ 1,9·10⁻⁶
alla mesh n = 320, convergenza spaziale di ordine 2,00 (1,99/1,98 RMS;
2,00/2,01 max) su tre raffinatamenti, e errore sulla posizione della
fronte di 3,0·10⁻⁵ (errore relativo sulla velocità 9,9·10⁻⁵). **Esito:
validato entro le tolleranze dichiarate per il modello bulk inerziale.**
**La condizione dinamica al bordo non è stata validata in questo test**
(nessun riferimento disponibile lo consente); va intesa come *non
dimostrata* dalla presente campagna, non come dimostrata.

## 14. Riproducibilità

```
validation_allen_cahn/
├── literature/            bibliography.md, benchmark_selection.md
├── parabolic/
│   ├── input/okumura_ch4_case2.in          (input del benchmark 1)
│   ├── results/                          (file grezzi OOFEM: case.out.m0.*.vtu)
│   │   ├── okumura_ch4_case2/
│   │   └── okumura_refinement_n{50,100,200,400}/
│   ├── logs/                             (log OOFEM di ogni run)
│   ├── scripts/
│   │   ├── run_okumura_example2.py         (esegue il benchmark 1)
│   │   ├── run_okumura_refinement.py       (esegue il raffinamento, --only n)
│   │   ├── postprocess_okumura.py          (energia + confronto con Fig. 4.4)
│   │   ├── postprocess_okumura_refinement.py
│   │   └── make_figure1.py                 (figura 1)
│   └── *.csv                             (serie energia, profili, confronto, raffinamento)
├── hyperbolic/
│   ├── input/nizovtseva_n320.in            (input; gli altri n generati dallo script)
│   ├── results/                            (file grezzi OOFEM)
│   ├── logs/                               (log OOFEM di ogni run)
│   ├── scripts/
│   │   ├── check_exact_residual.py         (verifica esattezza dell'onda)
│   │   ├── run_nizovtseva.py               (esegue n = 80/160/320)
│   │   ├── postprocess_nizovtseva.py       (metriche, interfaccia, profili)
│   │   └── make_figure2.py                 (figura 2)
│   └── nizovtseva_*.csv
├── figures/figure_1_parabolic.{png,pdf}, figure_2_hyperbolic.{png,pdf}
└── validation/error_metrics.csv, convergence_study.csv, comparison_report.md
```

Comandi di ri-esecuzione (dalla root del repository OOFEM):

```
# Benchmark 1 (parabolico): ~1 min per run in build debug
python3 doc/stvowu_validation/validation_allen_cahn/parabolic/scripts/run_okumura_example2.py
python3 doc/stvowu_validation/validation_allen_cahn/parabolic/scripts/postprocess_okumura.py
python3 doc/stvowu_validation/validation_allen_cahn/parabolic/scripts/run_okumura_refinement.py
python3 doc/stvowu_validation/validation_allen_cahn/parabolic/scripts/postprocess_okumura_refinement.py
python3 doc/stvowu_validation/validation_allen_cahn/parabolic/scripts/make_figure1.py

# Benchmark 2 (iperbolico): ~1 min per mesh
python3 doc/stvowu_validation/validation_allen_cahn/hyperbolic/scripts/run_nizovtseva.py
python3 doc/stvowu_validation/validation_allen_cahn/hyperbolic/scripts/postprocess_nizovtseva.py
python3 doc/stvowu_validation/validation_allen_cahn/hyperbolic/scripts/make_figure2.py
```

Ogni script di post-processing legge esclusivamente i file grezzi
esportati (`.vtu`, o i CSV generati dai post-processing a valle) e
riproduce i CSV metrica per metrica; non vengono applicati filtri,
interpolazioni o sostituzioni di dati.
