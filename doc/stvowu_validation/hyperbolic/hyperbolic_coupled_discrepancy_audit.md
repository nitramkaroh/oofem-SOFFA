# Audit della discrepanza iperbolica accoppiata OOFEM vs riferimento ODE semidiscreto

**Data:** 26 ottobre 2026
**Problema esaminato:** confronto tra il sistema iperbolico accoppiato bulk-superficie di OOFEM (`wave 1`) e il riferimento ODE semidiscreto indipendente (`run_coupled_independent.py`, DOP853). Errori RMS quasi costanti col raffinamento temporale:

| dt | bulk L2 RMS | superficie L2 RMS |
|---|---|---|
| 0.002 | 4.6937e-4 | 2.1058e-3 |
| 0.001 | 4.6967e-4 | 2.1070e-3 |
| 0.0005 | 4.6975e-4 | 2.1073e-3 |

(fonte: `coupled_independent_temporal.csv`; errori *crescenti* dello 0.08% mentre dt scende 4× — l'ordine stimato è ≈0/negativo).

**Metodo:** solo ispezione statica del codice e verifica aritmetica delle formule di assemblaggio. Nessuna nuova simulazione, nessuna modifica a solver o dati di riferimento.

---

## 1. Sistemi confrontati (equazioni semidiscrete attese)

Modello dichiarato nel rapporto, con `wavedamping` non specificato (default 0.0, [nltransienttransportproblem.h:84](../../../src/tm/EngineeringModels/nltransienttransportproblem.h#L84)) e η=0 nel riferimento:

- bulk: ρ a + Kb u = rb
- bordo: ρ a + c_Γ v + (Kb + Ks) u = rb + rs

con a = ü, v = u̇, rb = ∫Ω N m(u−u³) dx, rs = ∫Γ N βm(ψ−ψ³) dS, e l'accoppiamento k∂ₙu che entra automaticamente nelle righe di bordo di Kb u. Condizioni iniziali u₀ = 0.2 + 0.1 cos(πx/2) cos(πy), v₀ = 0 in entrambi.

## 2. Verifica termine per termine

| Termine | OOFEM | Riferimento Python | Esito |
|---|---|---|---|
| Inerzia bulk ρ a | `MassMatrix` con `Capacity → density` in modo iperbolico ([allencahnmat.C:166-169](../../../src/tm/Materials/allencahnmat.C#L166-L169)); escluso per gli elementi riga ([nltransienttransportproblem.C:75-86](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L75-L86)) | `Mb` con RHO=1, assembly Q1 standard con `det·w` ([run_coupled_independent.py:56-58](run_coupled_independent.py#L56-L58)) | **coincidenti** |
| Smorzamento bulk ρζ v | assente: `wavedamping` default 0.0, non dato in input | assente: `ETA=0.`, `C = ETA*Mb + Ms` (riga 92) | **coincidenti** |
| Diffusione bulk k | `TangentStiffnessMatrix` → `computeConductivityMatrix`, BC matrix nulla (nbc 0) ([transportelement.C:107-113](../../../src/tm/Elements/transportelement.C#L107-L113)) | `Kb` con `grad = dN @ inv(J)`, peso `det·w·w` — isoparametrico corretto | **coincidenti** |
| Reazione bulk m(u−u³) | ∫N·f via `computeInternalSourceRhsVectorAt`, GP temperature dallo status ([allencahnmat.C:183-190](../../../src/tm/Materials/allencahnmat.C#L183-L190), [transportelement.C:380-387](../../../src/tm/Elements/transportelement.C#L380-L387)); full-Newton (wave → `nrsolverFullNRM`, [nltransienttransportproblem.C:145-149](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L145-L149)) | `rb += MB·N·(uq−uq³)·det·w·w` (righe 72-75) | **coincidenti** |
| Dinamica superficiale c_Γ ṽ | riga di bordo: `γ/(βdt)·M_line` a sinistra, `M_line·v` a destra ([nltransienttransportproblem.C:67-72, 684-701](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L667-L710)) | `Ms·v` in `C` | **diverso di fattore 2** (vedi §3) |
| Diffusione tangenziale | `line1ht` standard isoparametrico: J = L/2 ([fei3dlinelin.C:113-114](../../../src/oofemlib/fei3dlinelin.C#L113-L114)), volume = detJ·w·area, area=1 ([line1_ht.C:100-106](../../../src/tm/Elements/line1_ht.C#L100-L106)) | `Ks` con `dN=[-.5,.5]/length`, `w=length*GAUSS_W[a]` (righe 60-64) | **diverso di fattore 1/2** (vedi §3) |
| Reazione superficiale βm(ψ−ψ³) | come sopra, β=surfacebeta=1 | `rs += MS·N·(uq−uq³)·length·w_a` (riga 80) | **diverso di fattore 2** |
| Accoppiamento k∂ₙu | implicito nelle righe di bordo di Kb u (Neumann naturale, nbc 0) | idem (stessa Kb) | **coincidenti** |
| CI posizione | `InitialCondition … u <val>` su dof 10, stessa formula coseno (case.in) | stessa formula (riga 84-86) | **coincidenti** |
| CI velocità | `initialvelocity 0.0` → `waveVelocity=0` su tutte le equazioni | `v0=0` (riga 98) | **coincidenti** |
| Accelerazione iniziale | risolta con M_bulk a₀ = F − K u₀ − Ms v₀ ([nltransienttransportproblem.C:473-531](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L473-L531)) | implicita da `rhs(0, z0)` | **coerenti**¹ |
| Quadratura | quad `nip 9` = 3×3; line `nip 3` = 3 punti Gauss | stessi punti/pesi 3-pt (riga 23-24) | **stessi punti** |
| Norme di confronto | — | RMS FE con Mb/Ms su stessi nodi, stesso tempo T (con assert su coordinate 1e-13 e tempo 1e-10) | procedura **corretta**² |

¹ Il riferimento usa qui le sue Ms/Ks (sbagliate): a₀ è coerente col proprio sistema, ma diverso da quello di OOFEM.
² Dettaglio minore: la Ms del riferimento (2×) usata in `rms_fe` con misura 2(LX+LY) gonfia di √2 la sola metrica superficiale ruotando il valore assoluto, non il plateau.

## 3. Differenza individuata: quadratura degli elementi riga nel riferimento

[run_coupled_independent.py:60-64](run_coupled_independent.py#L60-L64) (e analoghi alle righe 76-80 per `rs` e 168-172 per l'energia superficiale):

```python
for a,s in enumerate(GAUSS_X):                       # s ∈ [−1, 1], pesi che sommano a 2
    N=np.array([(1-s)/2,(1+s)/2]); dN=np.array([-.5,.5])/length; w=length*GAUSS_W[a]
    Ms[...] += CS*np.outer(N,N)*w
    Ks[...] += KS*np.outer(dN,dN)*w
```

Su un elemento riga isoparametrico con ξ ∈ [−1,1] il Jacobiano è **L/2**, non L: il peso corretto è `(L/2)·w_a` e la derivata fisica ∂N/∂s = `[−.5,.5]·(2/L)`. Verifica aritmetica diretta (nessuna simulazione: solo valutazione delle somme di quadratura su un elemento):

| quantità | riferimento / valore corretto |
|---|---|
| massa superficiale Ms | **2.0×** |
| rigidezza tangenziale Ks | **0.5×** |
| reazione superficiale rs | **2.0×** |

(per Ms: CS·L·∫₋₁¹ NNᵀ dξ = 2 × [CS·(L/6)[[2,1],[1,2]]]; per Ks: (1/4 dN_corr)·2·w_corr = ½ della rigidezza corretta k/L [[1,−1],[−1,1]]).

Il riferimento integra quindi **un'equazione di bordo diversa**:

> 2 c_Γ ψ_t − (k/2) Δ_Γ ψ + k ∂ₙu = 2 βm (ψ − ψ³)   (riferimento)
>
> c_Γ ψ_t − k Δ_Γ ψ + k ∂ₙu = βm (ψ − ψ³)   (modello dichiarato, OOFEM)

La riscalatura non è omogenea (2, 1/2, 2) quindi non è una ridefinizione di parametri: è un problema al contorno differente, con stato stazionario diverso. I termini bulk (Mb, Kb, rb, energia cinetica) sono invece assemblati correttamente nello script — coerentemente col fatto che i test bulk (Nizovtseva, kink, modo coseno, oscillatore) passano.

### Coerenza con la firma osservata degli errori

- Errore costante in dt ⇒ differenza **spaziale** (sistemi semidiscreti diversi), non di integrazione temporale: il riferimento è esatto in tempo (DOP853, rtol 2e-12) e Newmark ha errore O(dt²) che scomparirebbe ×4 per dimezzamento.
- Errori *leggermente crescenti* riducendo dt: OOFEM converge alla propria traiettoria semidiscreta allontanandosi dal riferimento "sbagliato" per un bias O(dt²) che sparisce.
- Il raffinamento temporale **interno** OOFEM (confronto run vs run più fine, `coupled_time_refinement.csv`) ha ordine 1.95: lo schema OOFEM è auto-consistente e converge.
- L'energia superficiale calcolata dal riferimento (`energy()`, righe 168-172) usa gli stessi pesi sbagliati: i suoi diagnostici energetici sono internamente coerenti e non rivelano il bug.

## 4. File e righe interessati

**Prima causa (riferimento Python):**
- `doc/stvowu_validation/hyperbolic/run_coupled_independent.py:60-64` — assembly Ms, Ks: `w=length*GAUSS_W[a]`, `dN=[-.5,.5]/length`
- `... :78-80` — rs con lo stesso peso doppio
- `... :168-172` — termini energetici superficiali con gli stessi fattori

**Codice OOFEM verificato come corretto per questi termini (non modificato):**
- [allencahnmat.C:162-190](../../../src/tm/Materials/allencahnmat.C#L162-L190) — Capacity iperbolico (rho/surfacecapacity), sorgente e tangente
- [nltransienttransportproblem.C:53-87](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L53-L87) — LHS Newmark (bulk 1/(βdt²), riga γ/(βdt), massa bulk-only)
- [nltransienttransportproblem.C:667-710](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L667-L710) — RHS dinamico (bulk a+ζv, riga v)
- [nltransienttransportproblem.C:473-531](../../../src/tm/EngineeringModels/nltransienttransportproblem.C#L473-L531) — accelerazione iniziale consistente
- [line1_ht.C:58-63, 100-106](../../../src/tm/Elements/line1_ht.C#L100-L106) — 3-pt Gauss, volume con Jacobiano
- [fei3dlinelin.C:49-70, 113-114](../../../src/oofemlib/fei3dlinelin.C#L113-L114) — N=(1∓ξ)/2, dNdx con 2/L, J=L/2

## 5. Diagnostici minimi (da eseguire solo su decisione; non eseguiti qui)

1. **Correzione a costo nullo del riferimento** (righe 62 e 80, e 170-172 per coerenza):
   `dN = np.array([-.5, .5]) * (2.0/length)`, `w = 0.5*length*GAUSS_W[a]`.
   Rieseguire solo l'integrazione DOP853 della suite temporale (i run OOFEM archiviati restano validi): predizione — gli RMS vs indipendente cadono di 1-2 ordini di grandezza e dimezzano in ×4 (ordine ≈2 di Newmark) a ogni dimezzamento di dt.
2. **Controllo matriciale senza alcuna simulazione**: confronto numerico elemento-per-elemento delle Ms/Ks/rs dello script contro i valori analitici dell'elemento P1 corretto (Ms = c_Γ L/6 [[2,1],[1,2]], Ks = k/L [[1,−1],[−1,1]]) — già confermato aritmeticamente in questa audit; basta renderlo test persistente.
3. **Triangolo di conferma opzionale**: con il riferimento corretto, rieseguire anche la suite spaziale (`--suite space`, nuovo run n32x16 mancante): ci si aspetta accordo che migliora con h, chiudendo contemporaneamente il buco di raffinamento spaziale segnalato nel rapporto.

## 6. Conclusione provvisoria

**La discrepanza nasce dal riferimento Python, non da OOFEM né dalla procedura di confronto.** L'assemblaggio manuale degli elementi riga nel riferimento dimentica il fattore di Jacobiano ½ della trasformazione isoparametrica (ξ ∈ [−1,1]), producendo massa e reazione superficiale doppie e rigidezza tangenziale dimezzata — un problema al contorno differente da quello implementato (correttamente, secondo lo standard isoparametrico verificato in `line1ht`/`FEI3dLineLin`) da OOFEM. Questo spiega esattamente la firma osservata: errori costanti in dt, ordine ≈0/negativo, leggera crescita col raffinamento, e contemporaneo ordine ~1.95 nella convergenza interna OOFEM.

Limiti di questa conclusione: la prova riguarda i termini superficiali; gli errori residui dopo la correzione (attesi O(dt²) + tolleranze) vanno riesaminati — in particolare la norma superficiale del riferimento usa la sua Ms erronea, da correggere insieme. Non ho verificato con una simulazione che la correzione chiuda numericamente i residui: il test 5.1 è il controllo minimo decisivo.
