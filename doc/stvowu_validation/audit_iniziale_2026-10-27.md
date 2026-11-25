# Audit iniziale della campagna di validazione Allen–Cahn (bulk + dinamica superficiale)

**Data:** 27 ottobre 2026
**Scopo:** inventario, verifica delle evidenze pregresse, conferma della correzione del riferimento Python, identificazione precisa delle lacune e piano minimo ordinato di test — **senza eseguire nuove simulazioni** (vincolo rispettato: in questa sessione sono state eseguite solo (a) ispezione statica del codice, (b) rieseguzione **in memoria, senza scrittura** del riferimento DOP853 corretto per la verifica di riproducibilità, (c) verifica simbolica/numerica della soluzione analitica di Nizovtseva, senza alcuna simulazione OOFEM).
**Stato a fine audit:** FERMATO, in attesa di decisione prima di qualsiasi nuova simulazione.

---

## 1. Fase I — Verifica del modello implementato (ispezione statica del codice)

| Componente | Equazione attesa | Implementazione OOFEM | Verifica effettuata | Esito | Evidenza |
|---|---|---|---|---|---|
| Bulk parabolico | `c u_t − kΔu = m(u−u³)` | matrice capacità `Σ NNᵀc dV`, rigidezza `Σ BᵀdB dV`, sorgente `Σ N m(u−u³) dV`, tangente `Σ NNᵀ m(3u²−1) dV` (IntSourceLHS, full-Newton) | ispezione termine per termine | **OK** | [transportelement.C:327-387,444-465](../src/tm/Elements/transportelement.C#L327-L387), [allencahnmat.C:162-190](../src/tm/Materials/allencahnmat.C#L162-L190) |
| Superficie parabolica | `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)` | `line1ht` isoparametrico (J=L/2), 3-pt Gauss (esatta per la cubica), β=`surfacebeta`; **attenzione**: in modalità parabolica la superficie condivide la `c` del bulk (non esiste `c_Γ` indipendente in input) | ispezione + verifica aritmetica dei fattori | **OK** (semantica documentata) | [allencahnmat.C:166-170](../src/tm/Materials/allencahnmat.C#L166-L170), [line1_ht.C:54-106](../src/tm/Elements/line1_ht.C#L54-L106), [fei3dlinelin.C:114](../src/oofemlib/fei3dlinelin.C#L114) |
| Bulk iperbolico | `ρu_tt + ρζu_t − kΔu = m(u−u³)` | Newmark: LHS `K + (1/(βdt²)+ζγ/(βdt))·ρM` (bulk-only), RHS `a+ζv`, accelerazione iniziale consistente `M_bulk a₀ = F − Ku₀ − M_Γv₀`; `wave 1` forza full-Newton | ispezione termine per termine | **OK** | [nltransienttransportproblem.C:53-87,473-531,667-710](../src/tm/EngineeringModels/nltransienttransportproblem.C#L53-L87) |
| Superficie iperbolica (1ª ordine) | `c_Γ ψ_t − kΔ_Γψ + k∂ₙu = βm(ψ−ψ³)` | riga LHS `K_Γ + γ/(βdt)·M_Γ`, RHS `M_Γ v`; **nessuna inerzia superficiale** (la massa d'onda esclude gli elementi riga) | ispezione | **OK** (limite di modello: la superficie inerziale di 2ª ordine non è implementata — non è una lacuna di validazione) | [nltransienttransportproblem.C:53-87](../src/tm/EngineeringModels/nltransienttransportproblem.C#L53-L87) |
| Capacità | — | parabolico: `c` su bulk **e** superficie; iperbolico: `rho` (inerzia bulk), `surfaceCapacity` (superficie), `c` non usata | ispezione | **OK** | [allencahnmat.C:162-190](../src/tm/Materials/allencahnmat.C#L162-L190) |
| Accoppiamento `k∂ₙu` | flusso normale implicito | DOF di tracciato condivisi tra quad e line; il termine entra come flusso naturale nelle righe di bordo del bulk (nbc 0, Neumann) — identico alla costruzione del riferimento `K = Kb + Ks` | ispezione + coerenza col riferimento | **OK** | §4 di [hyperbolic_coupled_discrepancy_audit.md](hyperbolic/hyperbolic_coupled_discrepancy_audit.md) |
| Quadratura | quadratura esatta per la reazione cubica | `quad1ht` default 2×2 (sotto-integra `N(u−u³)`, grado 4/direzione; la tangente grado 3 è esatta anche 2×2); `nip` letto da `Element::initializeFrom`; i casi di validazione usano `nip 9` (3×3, esatta); `line1ht` default 3-pt (esatta); RHS e assemblaggio usano la **stessa** regola (`giveDefaultIntegrationRulePtr() = integrationRulesArray[0]`) | ispezione + analisi del grado polinomiale | **OK con caveat**: gli input devono dichiarare `nip 9` (veri fatto nei benchmark) | [quad1_ht.C:61,84](../src/tm/Elements/quad1_ht.C#L61), [line1_ht.C:63](../src/tm/Elements/line1_ht.C#L63), [element.C:697](../src/oofemlib/element.C#L697), [element.h:836](../src/oofemlib/element.h#L836) |
| Condizioni iniziali | `u₀`, `v₀`, `a₀` consistente | `InitialCondition u` + `initialvelocity` + solve di `a₀` | ispezione | **OK** | [initialcondition.C](../src/oofemlib/initialcondition.C), [nltransienttransportproblem.C:473-531](../src/tm/EngineeringModels/nltransienttransportproblem.C#L473-L531) |
| Soluzione analitica di Nizovtseva (Eq. 29) | `τf_tt + f_t = Δf + (f−f³)` soddisfatta da `f = 1/(1+e^{λx−1.5t})`, `λ=√(2+9τ)/2` | — (dato di riferimento) | **verifica simbolica + numerica nuova** in questa sessione: residuo identico a `g(1−g)[A(1−g)+Bg]`, che si annulla per ogni g solo per `s=3/2, λ²=(9τ+2)/4`; max\|residuo\| = 1.1e-16 con τ=0.2, controlli non-vuoti (0.02–0.14), `v₀ = 1.5 f₀(1−f₀) = f_t(x,0)` esatto, velocità d'onda `c = 1.5/λ = 1.5389675281` | **OK — verificato** | [nizovtseva_pde_residual_check.py](hyperbolic/nizovtseva_pde_residual_check.py), [nizovtseva_pde_residual_check.csv](hyperbolic/nizovtseva_pde_residual_check.csv) |
| Leggi energetiche continue | `dE/dt = −∫c u_t² − ∫c_Γψ_t²` (parabolico); `dE_tot/dt = −∫ρζu_t² − ∫c_Γψ_t²` (iperbolico) | formule di post-processing coerenti con le leggi; **velocità non esportata dal solver** (ricostruita) | audit dei 4 script energetici + codice C++ | **OK** le formule; **LIMITATO** la misura (vedi §6, G5) | [energy_theory.md](energy_theory.md), audit energetico §1-3 |

**Bachi noti nel codice** (documentati in `validation_audit.md`, non toccati in questa sessione): (1) velocità/accelerazione prescritte sui DOF non onorate nel ramo Newmark (Dirichlet tempo-dipendente + `wave 1`); (2) sorgente MMS iperbolica usa `c·u_t` invece di `ρζ·u_t`. **Nessun benchmark di validazione li riguarda** (nessun test usa Dirichlet tempo-dipendente in wave mode; i casi MMS sono esclusi), ma vanno riportati nel report finale come limitazioni note.

## 2. Fase II — Matrice bibliografica

| Riferimento | Equazione | Caso numerico | Dati disponibili | Parametri completi | Confronto OOFEM possibile | Lacune |
|---|---|---|---|---|---|---|
| Okumura (diss. 2021, Ch.4 Ex.2) / Okumura–Furihata DCDS 2020 | Allen–Cahn 1D parabolico + BC dinamica | L=1, c_Γ=10, γ=K=100, u₀=e^{−500(x−0.5)²}, Δt=1e-3, T=0.7 | **curva di energia digitalizzata dal raster** (Fig. 4.4; 663 pt; crop+calibrazione manuali; asse t∈[0,0.7], E∈[0,90]); nessuna tabella | sì (dal testo della dissertazione) | sì, solo curva energia: 391 pt in 0.01≤t≤0.4, RMSE 0.9881 / MAE 0.8804 / max 1.8323 | incertezza raster (~0.1–1 unità); energia legacy calcolata con Gauss 2×2 (non esatta per il potenziale quartico); FD (paper) vs Q1 (OOFEM); **mismatch 391 vs 372 campioni non riconciliato** (vedi §4) |
| Nizovtseva 2016, Eq. (29) | AC iperbolico 1D `τf_tt+f_t=Δf+(f−f³)` | onda viaggiante τ=0.2, [−12,12]×[0,1], T=0.2 | **formula analitica** (verificata come soluzione esatta in questa sessione) | sì | **sì — il confronto più forte**: RMS 2.92e-5/7.32e-6/1.86e-6 (N80/160/320), ordini 1.995/1.976 spaziali, 2.009/1.917 temporali, velocità 9.8e-5 rel., troncamento e sensibilità IC quantificati | solo bulk (niente superficie dinamica); troncamento Neumann (quantificato); energia con `np.gradient` (NOT VERIFIABLE per la cinetica interna) |
| Xiao–Chen 2025, Case 4 | AC parabolico 2D + BC dinamica | Ω=[0,1]², goccia, γ=M=0.01, τ=1e-4, h=1/300 | **solo figure raster**; **ε non specificato nel paper** | **no** (ε assente) | **no** — solo qualitativo (contorni, monotonia energia); **ε=0.1 è nostra assunzione esplicita**, documentata in ogni script e report | ε assunto; h=0.05 vs 1/300; non convergenza spaziale a T=0.1 |
| Lee et al. 2024, Fig. 7 | AC parabolico 1D, equilibrio | kink tanh, [−1,1], P1, NX=200 | formula analitica (soluzione stazionaria) | sì (ε derivato 0.009606) | parziale — max nodale riconciliato; **L2: 6.0675e-3 (audit) vs 2.61689e-3 (ricalcolo), discrepanza 2.32× irrisolta** | nessun script/CSV metriche persistente |
| Wang, Wei & Chen 2025, Fig. 1 | AC parabolico 2D, energia | cos(2πx)cos(2πy), h=1/20, Δt=0.01 | solo figura energia raster | parzialmente | no — solo side-by-side; raw VTK solo in /tmp (non archiviati) | non riproducibile; NaN a ε=0.01; schema/mesh diversi |
| Akrivis & Li 2022, Fig. 3 | AC parabolico 2D, energia | τ=ε², h=ε/2, BDF5+Radau | solo figura energia raster | parzialmente | no — solo side-by-side; Radau non implementato | non riproducibile; ε=0.01 non valido |
| Hochbruck–Leibold 2019/26, 2020/20 | onda semilineare + BC cinetica | disco, profilo sin(2πt)x₁x₂ | profilo esatto + **sorgente ricostruita = MMS** | — | **escluso** ("Non ammesso (MMS)"); ordini spaziali stagnanti (0.33/0.24) | non è il benchmark del paper (inerzia superficiale, BC diverse) |
| Chen et al. 2025 §6.2 | AC MMS | u=0.05e^{−0.1t}sin x sin y, ε=0.01 | MMS del paper | — | **escluso come MMS** | — |
| beta_mms / suite MMS interne | AC accoppiato MMS | u=0.2+0.1e^{−t}[…] | MMS interno (Stefanelli–Voso–Wu) | — | **escluso** (code verification) | — |
| Calatroni–Colli 2012; Kovács–Lubich 2017; Makki–Miranville–Petcu 2022; Hipp–Kovács 2018/53 | teoria / modelli diversi | — | nessuna formula numerica utilizzabile | — | no | solo riferimenti teorici |

## 3. Evidenze già sufficienti (verificate in questa sessione o già documentate e coerenti)

1. **Bulk iperbolico — onda viaggiante di Nizovtseva.** Ordini spaziali 1.995/1.976, temporali 2.009/1.917, velocità d'onda 9.8e-5 rel., troncamento di dominio quantificato ([−8,8] vs [−12,12]: 1.7e-4 vs 1.3e-5), sensibilità IC (+1% velocità → 1.65e-4), **riproducibilità bit-identica** (3.61e-16 mismatch IC; ri-esecuzione `final_reproduction/nizovtseva` identica). In questa sessione: **la formula è verificata come soluzione esatta della PDE** (simbolico + numerico, §1). → **VALIDATO NEL CASO TESTATO** (bulk, quel profilo, quel dominio).
2. **Integrazione temporale bulk (Newmark)** su sistemi semidiscreti: modo coseno (p≈1.997, mobility 0), oscillatore uniforme (p≈2.02 vs DOP853), Nizovtseva temporale (2.01/1.92). → coerente con l'ordine 2 atteso.
3. **Sistema accoppiato iperbolico — convergenza temporale.** Dopo la correzione del riferimento (causa identificata e dimostrata lato Python, §4.2): bulk RMS 7.50e-7/1.91e-7/4.96e-8 (ordini 1.972/1.948), superficie 2.27e-6/5.78e-7/1.42e-7 (ordini 1.975/2.027), dt 0.002/0.001/0.0005, mesh 16×8, T=0.02; energia del riferimento corretto monotona, ΔE(T)−ΔE(OOFEM) = 6.5e-6 coerente con O(dt²).
4. **Riferimento Python corretto: conservato e riproducibile** — verificato in questa sessione senza sovrascrivere nulla: rieseguito DOP853 in memoria, **bit-identico** (max|ΔU|=0, max|ΔV|=0, energia max|ΔE|=0) rispetto a `independent_reference_fixed_time_n16x8_dt0.001.npz`/CSV; lo script originale buggy (`run_coupled_independent.py`, 232 righe) è intatto e non modificato; tutti gli artefatti di `coupled_independent_jacobian_fix/` sono presenti.
5. **Bilanci energetici — leggi e formule.** Le quattro implementazioni di post-processing sono coerenti con le leggi di dissipazione continue in entrambi i rami; quadrature 3×3/5×5 (esatte per il quartico su Q1); i residui moderni ad alta precisione sono coerenti con l'ordine dello schema: parabolico R(T) = −2.95e-4/−1.52e-4/−7.77e-5 (1ª ordine, backward Euler), iperbolico 3.69e-6/9.16e-7/3.47e-7 (2ª ordine, poi plateau ~7e-7 a dt=0.00025 = pavimento della ricostruzione velocità). Nessun incremento positivo dell'energia totale nei run moderni.
6. **Riproducibilità dei run OOFEM archiviati:** 3 ri-esecuzioni bit-identiche (parabolico n20, Niz n160, accoppiato dt=0.0005), 0 errori/0 warning.
7. **Esclusioni coerenti:** tutti i casi MMS (H-L, Chen 2025, beta_mms) sono correttamente esclusi dalle conclusioni di validazione.

## 4. Evidenze da verificare / documenti obsoleti / discrepanze irrisolte

### 4.1 Documenti "finali" obsoleti (aggiornare, nessuna simulazione)
- **`validation_matrix_final.md/.csv` e `OOFEM_AllenCahn_Validation_Report.md` predata la risoluzione del 26-27/10**: la matrice riporta ancora *Sistema iperbolico accoppiato: FAIL* e il report "La causa resta da isolare". Entrambi sono superati da [hyperbolic_coupled_discrepancy_audit.md](hyperbolic/hyperbolic_coupled_discrepancy_audit.md) (causa: Jacobiano delle linee nel riferimento Python, Ms 2×, Ks ½×, rs 2×) e da [jacobian_fix_verification_report.md](hyperbolic/coupled_independent_jacobian_fix/jacobian_fix_verification_report.md).
- **"Mancano il run e il riferimento n32×16" è metà falso**: `hyperbolic/coupled_independent_space_n32x16_dt0p00025/` **contiene il run OOFEM completo** (case.in, .out, 80 VTU, mesh.npz, log); manca solo l'`independent_reference.npz` (e in ogni caso il riferimento archiviato lì sarebbe quello buggy — da rigenerare con la versione corretta).
- `literature_comparison/README.md` afferma "il solver parabolico è validato… nel caso 1D" senza qualificazioni: incoerente con la matrice (PARTIAL) e con le riserve su raster/quadrature. Da riformulare.
- Le metriche superficiali "prima" (2.1058e-3 ecc.) citate nel report sono gonfiate di √2 dalla Ms errata del riferimento (plateforme corrette 1.4891/1.4899/1.4901e-3): da annotare.
- I due bachi codice noti (§1) non hanno riga di stato nella matrice.

### 4.2 Discrepanze numeriche irrisolte (risolvibili a costo zero, solo post-processing)
- **Okumura: 391 campioni (RMSE 0.9881) vs 372 (RMSE 0.99288)** tra `literature_comparison/` e `final_reproduction/okumura/literature_metrics.csv` — stessa casistica, metriche diverse, mai riconciliate.
- **Lee 2024: L2 6.0675e-3 vs 2.61689e-3** (2.32×) tra `validation_audit.md` e il ricalcolo di `validation_metrics_audit.md`; provenienza del valore vecchio non recuperabile; nessun script metriche persistente.
- **CSV energetici legacy** (`energy_time_series.csv`, `energy_validation_finer.csv`): **non riproducibili** dai loro stessi output (ricalcolo: R(T)=6.87e-7, 0 incrementi positivi, contro 3.47e-5 e 17 dichiarati) → da marcare esplicitamente "non riproducibili, non sono evidenza" (già fatto negli audit; da propagare nel report finale).
- Figura `parabolic_2d/energy_balance_secant.png`: pannello superiore dimensionalmente incoerente (≈4967.6 vs 0.4959); il report stesso dice "non va usato", ma il PNG resta archiviato accanto al CSV. Da marcare.
- Link rotto nel report: `parabolic_2d/energy_interface.csv` → i file reali sono in `parabolic_2d/runs/*/`.

### 4.3 Stub e record incompleti (da marcare, non da "completare" a meno di necessità)
- `energy_parabolic_single_test_run/` — solo input, nessun output mai prodotto.
- `final_reproduction/energy_metrics_verified/` — directory vuota.
- `energy_validation_runs_highprecision/parabolic_n32x16_dt0p0005/` — 4 di 40 VTU.
- `final_reproduction/energy/parabolic_n16x8_dt0p001/` — 0 VTU (il CSV energia è però completo).
- Suite legacy `nizovtseva_*` v1 e `hl_disk_*` — solo stato finale, senza `.out`/storia: non re-derivabili (v1 superata da v2; H-L esclusa come MMS).
- Wang/Akrivis — raw VTK solo in `/tmp`, non versionati: **non riproducibili**.

## 5. Benchmark incompleti (in sintesi)

| Benchmark | Stato | Cosa manca |
|---|---|---|
| Accoppiato iperbolico — **spaziale** vs riferimento corretto | **incompleto (lacuna principale)** | riferimento DOP853 corretto su n8×4/n16×8/n32×16 mai eseguito; per la convergenza spaziale serve anche un livello-proxy fine (n64×32) |
| Accoppiato iperbolico — temporale | completo (dopo fix) | — |
| Bulk iperbolico (Niz, kink, coseno, oscillatore) | completo | — |
| Okumura 1D | completo ma metriche non riconciliate | riconciliazione 391/372 |
| Xiao–Chen 2D | qualitativo per costruzione (ε assente nel paper) | — (limite di dato, non di esecuzione) |
| Lee 2024 | metriche non persistenti | script/CSV metriche + riconciliazione L2 |
| Bilanci energetici | moderno completo; legacy non riproducibile | marcatura; (opzionale) export velocità |
| Wang / Akrivis | non riproducibili (raw in /tmp) | da dichiarare NON APPLICABILE/IN ATTESA nel report finale |

## 6. Lacune precise

- **G1 (decisiva):** convergenza **spaziale** del sistema accoppiato iperbolico non dimostrata. Attenzione al design: il confronto "OOFEM(h) vs DOP853 sullo stesso mesh" misura solo l'errore temporale (stesso sistema semidiscreto); per l'errore spaziale serve un riferimento (quasi) indipendente dal mesh. Design minimo: traiettoria DOP853 su **n64×32** come proxy di riferimento (storia completa, nodi dei mesh fini a sottoinsieme dei nodi del proxy), confronto con i run OOFEM archiviati n8×4/n16×8/n32×16 (dt=2.5e-4, T=0.02, già completi) — **nessuna riesecuzione OOFEM necessaria**.
- **G2:** documenti finali obsoleti (§4.1).
- **G3:** discrepanze metriche Okumura e Lee (§4.2).
- **G4:** energia cinetica iperbolica basata su velocità ricostruita (pavimento ~7e-7); `postprocess_energy.py` rifiuta la stima e lascia `kinetic` vuoto nei VTU attuali.
- **G5:** generalità: l'accoppiato è validato (temporalmente) su un'unica combinazione di parametri (ρ=1, ζ=0, c_Γ=1, β=1, T=0.02, quella IC).
- **G6 (limite di modello, non di validazione):** superficie inerziale di 2ª ordine non implementata.

## 7. Piano minimo ordinato di test (priorità, costo, criteri di accettazione)

Regola rispettata: prima ogni gruppo — test proposti, lacuna risolta, costo, check ridondanza — e un test alla volta. **Nessun sweep esteso.**

### T1 (priorità 1) — Convergenza spaziale dell'accoppiato con riferimento corretto [risolve G1]
- **Cosa:** (a) generare la traiettoria DOP853 corretta su n64×32 (dt non rilevante: rtol 2e-12/atol 2e-14, T=0.02); (b) opzionale ma economico: riferimenti corretti su n8×4/n16×8/n32×16 come controllo di consistenza interna; (c) confronto nodi-per-nodi (i mesh sono nidificati) dei 3 run OOFEM archiviati vs proxy n64×32, RMS pesati FE (massa del mesh coerente, normalizzata per area/perimetro) + max nodale, bulk e superficie, a T=0.02.
- **Costo stimato:** DOP853 n64×32 ≈ 1–10 min CPU (stime dai run 16×8 ≈ 3-4 s; lavoro per RHS ×≈14.7); post-processing: secondi. **Nessuna simulazione OOFEM.**
- **Criteri di accettazione:** (i) errore proxy-monotono decrescente a ogni raffinamento; (ii) ordini in [1.5, 2.5] su entrambi gli intervalli (n8→n16, n16→n32) per bulk e superficie; (iii) sanity del proxy: l'errore temporale interno a n32×16 (~1e-8, dal T temporale) è trascurabile rispetto a ‖u₃₂−u₆₄‖; (iv) dichiarati: 3 livelli OOFEM + 1 proxy, norme FE e nodali, dt=2.5e-4.
- **Ridondanza:** nessuna — è l'unico dato spaziale mancante sull'accoppiato.
- **T1b (condizionata):** se l'ordine risulta < 1.5 o ambiguo → aggiungere il livello OOFEM n64×32 (≈ 1–2 min) + confronto: totale comunque < 15 min.

### T2 (priorità 2) — Coerenza documentale [risolve G2, parte di G3] — nessuna simulazione
Aggiornare `validation_matrix_final.md/.csv` e `OOFEM_AllenCahn_Validation_Report.md`: stato accoppiato FAIL→(esito di T1); "run mancante"→"riferimento mancante (run esistente)"; √2 sulle metriche "prima"; status Okumura qualificato; CSV legacy non riproducibili; figura secant marcata; link rotto; righe per i due bachi codice; stub marcati; Wang/Akrivis → NON APPLICABILE (dati non archiviati). Costo: tempo di scrittura. Accettazione: ogni riga della matrice conciliabile con un artefatto archiviato; nessuna affermazione "validato" senza qualifiche.

### T3 (priorità 3) — Riconciliazione metriche Okumura e Lee [risolve G3] — nessuna simulazione
Post-processing sui CSV/VTU archiviati: definire una sola finestra di confronto e una sola norma per benchmark, documentarle, o dichiarare esplicitamente le due finestre. Per Lee: script metriche persistente + scelta motivata della norma (1D trapezoidale). Costo: minuti.

### T4 (priorità 4) — Consolidamento energetico [risolve G4]
- (a) dichiarare formalmente il pavimento di ricostruzione (~7e-7, dt=0.00025) come incertezza di R(T) iperbolico; marcare i CSV legacy (già fatto negli audit → propagare).
- (b) **opzionale, richiede approvazione esplicita (modifica codice OOFEM, mirata e non numerica):** export di `waveVelocity` nei VTU → ricalcolo del bilancio 16×8 con velocità del solver. Accettazione: R(T) sotto il pavimento di ricostruzione, 0 incrementi positivi. Senza (b), l'energia resta evidenza di coerenza interna, non di validazione indipendente.

### T5 (priorità 5, facoltativa, basso valore) — Rafforzare il 2D parabolico
Convergenza spaziale interna a T=0.0625 con N=80 (≈ 5–10 min). Criteri solo interni (nessun target esterno esiste per Xiao–Chen): decrescita dell'errore di livello. Fare **solo** se richiesto: non cambia lo status (resterà qualitativo vs letteratura a causa di ε).

**Costo totale del piano minimo (T1+T2+T3+T4a):** < 15 min di CPU + lavoro documentale. Timeout suggeriti: 10 min per T1 (margine 4× sull'estimato).

## 8. Criteri di accettazione generali (Fase V)

- Ordine di convergenza: dichiarare sempre numero di livelli e norma; accettare [1.5, 2.5] per schemi 2ª ordine (Newmark) e [0.75, 1.25] per 1ª ordine (backward Euler); ordini stimati su almeno 2 intervalli quando possibile.
- Distinguere e riportare norma FE (massa-pesata, normalizzata) e max nodale; dichiarare la misura (area/perimetro).
- Nei test spaziali verificare che l'errore temporale non domini (a dt=2.5e-4 il pavimento temporale è ~1e-8, trascurabile).
- Nessuna soglia universale di "accordo": i confronti con dati digitalizzati (Okumura) si esprimono come accordo favorevole entro l'incertezza di digitalizzazione dichiarata.
- Energia: R(T) = E(T)−E(0)+D(T) deve scalare con l'ordine dello schema e restare sotto il pavimento documentato della ricostruzione; la monotonia da sola non è validazione.
- Stati di accettazione del report finale: VALIDATO NEL CASO TESTATO / SUPPORTATO PARZIALMENTE / NON VALIDATO / NON APPLICABILE / IN ATTESA. Mai "perfettamente validato"; mai "solver validato integralmente".

## 9. Conclusioni che è già lecito formulare

1. Il modello parabolico (bulk + superficie dinamica) è assemblato coerentemente con le equazioni continue; l'integrazione temporale (backward Euler) converge con ordine 1 sul sistema semidiscreto; le leggi energetiche sono implementate correttamente e i residui moderni sono coerenti con l'ordine dello schema; il confronto 1D con i dati digitalizzati di Okumura è favorevole entro l'incertezza di digitalizzazione (RMSE 0.99 su scala di energia 90). → **parabolico: SUPPORTATO (1D quantizzato, 2D qualitativo)**.
2. Il **bulk iperbolico è VALIDATO NEL CASO TESTATO** contro la soluzione analitica di Nizovtseva (ordini ≈2 spaziali e temporali, velocità d'onda a 1e-4, soluzione verificata come esatta in questa sessione, troncamento e sensibilità quantificati, riproducibilità bit-identica).
3. Il **sistema accoppiato iperbolico (bulk inerziale + superficie 1ª ordine) è SUPPORTATO PARZIALMENTE**: la convergenza temporale è verificata contro un riferimento semidiscreto indipendente corretto (ordini 1.95–2.03, residui finali 5e-8/1.4e-7); la discrepanza originaria è identificata, dimostrata e documentata come difetto del riferimento Python, non di OOFEM.
4. I bilanci energetici dimostrano coerenza delle implementazioni con le leggi di dissipazione continue in entrambi i rami, con incertezza nota legata alla ricostruzione della velocità.
5. I casi MMS sono correttamente esclusi; la riproducibilità dei run archiviati è bit-identica (3 casi) e quella del riferimento corretto è bit-identica (verificata oggi).

## 10. Conclusioni che NON è ancora lecito formulare

1. Validazione **spaziale** del sistema accoppiato iperbolico (pende da T1).
2. "Il solver è validato integralmente" / "il modello Allen–Cahn è validato in generale" — non supportato da alcuna evidenza presente.
3. Validazione del 2D parabolico contro Xiao–Chen (resterà qualitativa: ε non nel paper, dati solo raster, mesh 6× più grossolana).
4. Generalizzazione del risultato accoppiato ad altre combinazioni di parametri/IC (singola combinazione testata).
5. Bilancio energetico iperbolico come verifica **indipendente** del solver (velocità ricostruita; da esportare per diventare tale — richiede modifica codice con approvazione).
6. "Parabolico 1D validato" senza qualifiche (incertezza raster + quadratura 2×2 legacy nell'energia).
7. Riproduzione dei casi Wang/Akrivis (dati non archiviati → NON APPLICABILE finché non ri-archiviati).
8. Dinamica superficiale inerziale di 2ª ordine: non implementata, non validabile (limite di modello).

---

## 11. Fermo e richiesta di decisione

Come richiesto, mi fermo qui. **Nessuna nuova simulazione è stata avviata.** Prossimo passo proposto, in attesa di decisione: **T1** (riferimento DOP853 corretto su n64×32 + confronto con i 3 run spaziali accoppiati già archiviati; nessuna riesecuzione OOFEM; costo stimato < 15 min di CPU). In alternativa o successivamente: T2 (documentazione) e T3 (metriche), che non richiedono simulazioni.
