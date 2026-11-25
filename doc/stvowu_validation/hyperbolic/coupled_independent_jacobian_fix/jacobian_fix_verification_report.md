# Verifica della correzione del Jacobiano negli elementi di bordo del riferimento indipendente

**Data:** 27 ottobre 2026
**Scopo:** verificare se la discrepanza OOFEM ↔ riferimento ODE indipendente (errori costanti in dt, 4.69e-4 bulk / 2.11e-3 superficie a T=0.02, v. [../hyperbolic_coupled_discrepancy_audit.md](../hyperbolic_coupled_discrepancy_audit.md)) era causata dagli errori di quadratura degli elementi di bordo nel riferimento Python.

**Vincoli rispettati:**
- OOFEM non è stato rieseguito: letti solo i campi VTU (Float64) archiviati in `../coupled_independent_time_n16x8_dt{0.002,0.001,0.0005}`.
- `run_coupled_independent.py` e tutti i risultati precedenti non sono modificati (importato in sola lettura; la traiettoria originale `independent_reference.npz` archiviata è riutilizzata per i valori "prima").
- Eseguita solo l'integrazione DOP853 del riferimento corretto (3 run, sequenziale, ~11 s totali), nessun MMS/WED, nessun sweep.
- Output nuovi solo in questa directory.

---

## 1. Modifiche al riferimento (file e righe)

Il riferimento originale non è toccato. La versione corretta vive in `check_jacobian_fix.py` (stessa directory) e differisce dall'originale **solo** nei tre fattori isoparametrici degli elementi riga:

| termine | originale (`run_coupled_independent.py`) | corretto (`check_jacobian_fix.py`) |
|---|---|---|
| peso quadratura, linea (Ms, Ks) | riga 62: `w = length*GAUSS_W[a]` | `w = 0.5*length*GAUSS_W[a]` |
| derivata fisica, linea (Ks) | riga 62: `dN = np.array([-.5,.5])/length` | `dN = np.array([-.5,.5])*(2.0/length)` |
| peso quadratura, reazione (rs) | riga 80: `… * length*GAUSS_W[a]` | `… * 0.5*length*GAUSS_W[a]` |
| termini energetici di superficie | righe 169-172 (stesso schema) | `energy_fixed`: `w = 0.5*length*GAUSS_W[a]`, `dN = [-.5,.5]*(2.0/length)` |

Nessun altra modifica: mesh, costanti (RHO=KB=MB=CS=KS=MS=1, ETA=0), punti di Gauss, IC, tolleranze DOP853 (rtol 2e-12, atol 2e-14), struttura dell'ODE (`Mb⁻¹(−(ETA·Mb+Ms)v − (Kb+Ks)u + rb + rs)`) e definizioni di norma invariate.

## 2. Dimostrazione dei fattori di quadratura corretti

Mappatura isoparametrica P1 di un elemento riga di lunghezza L, ξ ∈ [−1,1]:
x(ξ) = (x₁+x₂)/2 + (L/2)·ξ·t̂ ⇒ **J = L/2**;
N = [(1−ξ)/2, (1+ξ)/2], ∂N/∂x = (∂N/∂ξ)·(1/J) = [−½,½]·(2/L) = **[−1,1]/L**;
∫_e f ds = (L/2)·Σₐ wₐ f(ξₐ).

Il codice originale usava peso `L·wₐ` (2×) e `dN=[−½,½]/L` (½×): di conseguenza Ms e rs uscivano 2× e Ks ½× il valore dell'elemento P1 (verificato aritmeticamente nell'audit). Verifiche numeriche dell'assemblaggio corretto (nessuna simulazione):

| check | risultato | atteso |
|---|---|---|
| quadratura manuale di un elemento, Ms_e | [[0.041667, 0.020833],[0.020833, 0.041667]] | (L/6)[[2,1],[1,2]], L=0.125 ✓ |
| quadratura manuale di un elemento, Ks_e | [[8, −8],[−8, 8]] | (1/L)[[1,−1],[−1,1]] ✓ |
| reazione di linea con u costante 0.2 (nodo non d'angolo) | L/2·m(0.2) per elemento | partizione dell'unità ✓ |
| Σ elementi Ms | 6.0 = perimetro (orig: 12.0 = 2×) | ∫_Γ 1 dS ✓ |
| righe di Ks | somma 0 | proprietà di rigidità ✓ |
| Σ elementi Mb | 2.0 = area | ✓ |
| nnz di Ms | 144 = 4·48 − 48 (ciclo chiuso, 48 nodi angolari condivisi da 2 elementi) | ✓ |

(NB: il blocco assemblato su un nodo d'angolo contiene i contributi di **entrambi** gli elementi adiacenti; il confronto col blocco mono-elemento va fatto sulla quadratura manuale, non sul blocco assemblato.)

## 3. Risultati: RMS "prima" (riferimento originale) vs "dopo" (riferimento corretto)

Confronto per-passo con i campi OOFEM archiviati; i valori "prima" riproducono **bit-per-bit** il CSV archiviato `../coupled_independent_temporal.csv` (es. 4.6937105848466384e-4), quindi il collaudo è equivalente alla pipeline originale. Tabella finale (t = T = 0.02):

| dt | bulk RMS prima | bulk RMS dopo | rapporto | ordine | superficie RMS prima (norma orig, 2×Ms) | superficie RMS dopo (norma corretta) | ordine (dopo) |
|---|---|---|---|---|---|---|---|
| 0.002 | 4.6937e-4 | 7.4998e-7 | — | — | 2.1058e-3 | 2.2734e-6 | — |
| 0.001 | 4.6967e-4 | 1.9119e-7 | 3.923 | 1.972 | 2.1070e-3 | 5.7830e-7 | 1.975 (bulk/surf) |
| 0.0005 | 4.6975e-4 | 4.9563e-8 | 3.858 | 1.948 | 2.1073e-3 | 1.4190e-7 | 2.027 (surf) |

Ordine "prima" (da CSV archiviato): bulk −0.0009, superficie −0.0008 (piattaforma, leggermente crescente).
Max nodali finali: bulk 4.37e-6 → 2.67e-7; superficie 1.01e-5 → 6.08e-7 (dopo).

Note sulla metrica superficiale "prima": la norma originale usava la Ms sbagliata (2×) con misura 2(LX+LY), gonfiando il valore di √2; con la norma corretta la piattaforma "prima" è 1.4891/1.4899/1.4901e-3 (rapporto esatto √2 = 1.4141 contro i valori con norma originale), quindi anche la piattaforma di superficie si spiega interamente con i fattori errati, non con la norma.

### Andamento temporale (dt = 0.001, estratto; file completo `per_step_time_n16x8_dt0.001.csv`)

| passo (t) | bulk RMS prima | bulk RMS dopo | superficie RMS dopo |
|---|---|---|---|
| 1 (0.001) | 1.5124e-6 | 2.1753e-8 | 4.886e-8 |
| 5 (0.005) | 3.5858e-5 | 7.1309e-8 | 2.624e-7 |
| 10 (0.01) | 1.3390e-4 | 1.2012e-7 | 4.318e-7 |
| 15 (0.015) | 2.8177e-4 | 1.4555e-7 | 4.592e-7 |
| 20 (0.02) | 4.6967e-4 | 1.9119e-7 | 5.7830e-7 |

L'errore "dopo" cresce da 2.2e-8 a 1.9e-7 su T (accumulo del truncation error locale O(dt²) — coerente con Newmark), mentre quello "prima" cresce di 3 ordini di grandezza: la piattaforma era l'accumulo in tempo del bias di sistema del riferimento.

### Energia del riferimento corretto (dt=0.001)

`independent_energy_fixed_time_n16x8_dt0.001.csv`: totale 1.93160616 → 1.92996107, monotono decrescente, cinetica 0 → 3.56e-3, superficie 1.44273 → 1.43779. Confronto con il run OOFEM archiviato dt=0.0005 citato nel rapporto (1.93160616 → 1.92995451): ΔE(T) = 6.5e-6, coerente con O(dt²).

## 4. Coerenza termine per termine del riferimento corretto con OOFEM

| termine | riferimento corretto | OOFEM |
|---|---|---|
| inerzia bulk | Mb (RHO=1), Q1, det·w | MassMatrix con density (iperbolico), quad1ht isoparametrico |
| smorzamento bulk | ETA=0 | wavedamping default 0.0 |
| diffusione bulk | Kb (K=Kb+Ks nel residuo) | TangentStiffnessMatrix |
| reazione bulk | rb, 3×3 Gauss | computeInternalSourceRhsVectorAt, nip 9 |
| massa superficiale | Ms = c_Γ(L/6)[[2,1],[1,2]] per elemento (Σ=perimetro) | line1ht, J=L/2, surfacecapacity |
| diffusione tangenziale | Ks = k/L[[1,−1],[−1,1]] per elemento | line1ht, dNdx=2/L |
| reazione superficiale | rs = ∫N βm(ψ−ψ³) dS, 3-pt Gauss | line1ht nip 3, surfacebeta |
| accoppiamento k∂ₙu | implicito nelle righe di bordo di Kb | implicito (nbc 0, Neumann) |
| IC | u₀ coseno, v₀=0, a₀ da rhs(0) | InitialCondition u, initialvelocity 0.0, a₀ risolta (M_bulk a₀ = F−Ku₀−Ms v₀) |
| ODE | Mb⁻¹(−Ms v − (Kb+Ks)u + rb + rs), DOP853 | Newmark β=1/4, γ=1/2 (bulk 2a, riga 1a) |

## 5. Verdetto (distinto nei tre livelli richiesti)

1. **Causa statica individuata: SÌ.** Il plateau dt-indipendente (4.69e-4/2.11e-3, ordine ≈0) scompare interamente con la correzione dei soli fattori di quadratura delle linee: il bias era il riscalaggio errato (2×, ½×, 2×) di massa, rigidità e reazione di superficie nel riferimento Python, che integrava un problema al bordo diverso dal modello implementato.

2. **Riferimento corretto verificato: SÌ.** Con la correzione, OOFEM converge verso il riferimento con ordine **1.95–2.03** (atteso: 2 di Newmark) a ogni dimezzamento di dt, e il residuo finale scende a 5.0e-8 (bulk) / 1.4e-7 (superficie) per dt=0.0005 — livello coerente con il truncation error O(dt²) e la tolleranza non lineare OOFEM (rtol 1e-8). Nessun residuo sistemico in dt: il sistema semidiscreto del riferimento corretto e quello di OOFEM sono lo stesso sistema (stessa mesh, stesse matrici, stessi punti di quadratura, stesse IC).

3. **Accoppiamento iperbolico validato: parzialmente (solo la parte temporale di questo test).** Il cross-check semidiscreto indipendente valida ora la dinamica **accoppiata bulk-inerziale/superficie-di-primo-ordine** per integrazione temporale, su mesh 16×8, T=0.02, con questa IC e questi coefficienti. Non è una validazione completa: (i) resta aperto il buco di raffinamento **spaziale** del sistema accoppiato (mancano il run e il riferimento n32×16, già segnalati nel rapporto) — si consiglia di rieseguire la suite spaziale con il riferimento corretto, ora che il bias è rimosso; (ii) il test copre un'unica combinazione di parametri/IC; (iii) la superficie iperbolica inerziale non è un'equazione implementata (limite di modello, non di validazione). **Nessuna dichiarazione di validazione integrale del solver** è supportata da questi risultati.

## 6. Evidenze conservate (questa directory)

- `check_jacobian_fix.py` — riferimento corretto + driver di confronto (importa l'originale in sola lettura)
- `per_step_time_n16x8_dt{0.002,0.001,0.0005}.csv` — RMS/max per passo, prima/dopo
- `independent_reference_fixed_time_n16x8_dt{...}.npz` — traiettorie U, V del riferimento corretto (+ Ms, Ks)
- `independent_energy_fixed_time_n16x8_dt{...}.csv` — energia del riferimento corretto
- `summary_final.csv` — tabella finale e ordini
