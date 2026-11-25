# Unico passo successivo raccomandato — run OOFEM n64×32 (NON eseguito; richiede autorizzazione)

Data: 2026-10-27 — Derivato da: `audit_report.md` (verdetto §5: opzione 2, gap residuo §5.5).
**Questo documento è una proposta, non un esito: nessun run è stato avviato.**

## 1. Cosa fare (una sola cosa)

Eseguire **un solo** run OOFEM del caso accoppiato iperbolico identico a T1, con mesh
**n64×32** (2145 nodi, h=0.03125), DT=2.5e-4, T=0.02 (80 passi), stesso input di
`hyperbolic/coupled_independent_space_n32x16_dt0p00025/` (stessi parametri, IC, BC,
exporter VTU a ogni passo). Nessuna altra simulazione, nessuna integrazione DOP853,
nessun raffinamento temporale, nessun sweep.

## 2. Perché questo e non altro (criterio di scelta, task §5)

Il gap residuo dopo l'esaustivo post-processing (audit §10) è esattamente uno: la
catena OOFEM si ferma a h=0.0625. Il 4° punto (h=0.03125) è l'unico dato in grado di:

1. produrre **due rapporti consecutivi** sulla catena OOFEM a T e sui 10 istanti
   (oggi ce n'è uno solo per intervallo: 0.7173, 1.0478 — entrambi sotto [1.5,2.5]);
2. misurare il **floor a h=0.03125** (oggi NON VERIFICABILE: 2.1e-8/2.4e-8 a n16,
   2.6e-8/2.9e-8 a n32 — la sequenza dice ~3e-8, ma va misurata);
3. permettere un Richardson 3-livello **OOFEM** (8/16/32/64) invece che solo
   semi-discreto.

Alternative scartate (e perché):

- **Nuova integrazione DOP853 n256**: non indispensabile — i riferimenti n64 e n128 sono
  già su disco (`reference_n64x32_dt0p00025.npz`, `reference_n128x64_dt0p00025.npz`,
  rtol 2e-12) e la stabilità del proxy è ora coperta da bound rigorosi
  (`results/proxy_stability_bound.csv`, 12/12 SI). Il gap è OOFEM, non riferimento.
  (Rispetta il vincolo: nessuna DOP853 prima di aver dimostrato che è indispensabile.)
- **Nuovi tempi terminali / parametri**: non risponde alla domanda su h (l'ordine
  asintotico è una proprietà in h a T fissato); la dipendenza temporale è già
  caratterizzata su 80 passi (`results/family_time_series_recomputed.csv`).
- **Altro post-processing**: esausto — tutto ciò che i dati esistenti ammettono è già
  in `results/` (audit §10).

## 3. Stima del costo (verificabile prima dell'esecuzione)

| Voce | Stima | Base |
|---|---|---|
| CPU | **≈ 40–90 s** (1 core) | T1-A: n32×16 = 10.25 s wall / 10.18 s CPU (`t1a_n32x16_rerun/results/t1a_report.md` §1); n64 = 4× DOF (2145 vs 561) |
| RAM | < 1 GB | T1-B DOP853 su sistema 8385-DOF: RSS ≈ 109 MB; OOFEM n64 (2145 DOF) è minore |
| Disco | VTU 81 file × ~0.5 MB ≈ < 100 MB | dai run n32 esistenti |
| Budget proposto | **5 min wall, 1 core, nessuna parall** | 5× il peggio stimato |

Nessun rischio di blocco macchina: stesso run della serie T1, stessa build
(`rel-2.1-5430-g2cc1a7706-dirty`), stesso runner.

## 4. Dove (nessuna sovrascrittura)

- Esecuzione: `hyperbolic/t1c_n64x32_oofem/oofem/` (nuova)
- Risultati: `hyperbolic/t1c_n64x32_oofem/results/` (nuova)
- Post-processing: `hyperbolic/t1c_n64x32_oofem/validate_t1c.py` — riusa
  `audit_postproc.py` (funzioni `read_oofem_step`, `metrics4`, restrizione n64→n128
  già implementate e verificate) senza toccare script originali.

## 5. Criteri di accettazione pre-registrati (da applicare DOPO il run)

1. **Integrità run:** 80 step, exit code 0, tempi esportati == griglia canonica
   (|t−k·DT|<1e-12), VTU riordinabili per coincidenza esatta di coordinate con
   max_d2=0.0 (come i 30 check VTU dell'audit).
2. **Floor (criterio chiave):** OOFEM n64 vs riferimento DOP853 n64 (stessa mesh,
   `reference_n64x32_dt0p00025.npz`): bulk FE RMS a step 39 e 80.
   - atteso: ≤ O(1e-7) (coerenza con 2.1e-8…2.9e-8 di n16/n32);
   - se > 1e-6: problema del solutore/exporter — si ferma la catena di conclusioni
     e si diagnostica (nessun verdetto di convergenza).
3. **Catena a 4 punti a T (bulk e surf, FE RMS + max):**
   - monotonia FE RMS 8→16→32→64;
   - i **due** rapporti consecutivi (16→32, 32→64) entrambi in [1.5, 2.5] a T
     **e** in ≥ 6 dei 10 istanti {0.00025…0.02};
   - stabilità: errore vs proxy128 e vs proxy64 del 4° punto differiscono ≤ 8%
     (criterio T1, esteso al nuovo punto).
4. **Richardson 3-livello OOFEM a T:** p_2step = log2(e16/e64)/2 coerente con i due
   rapporti (|p_2step − media| < 0.5).

**Se 2+3 passano** → l'ordine 2 è **dimostrato sull'intervallo h∈[0.03125, 0.25]**
per questo caso (T=0.02, parametri dati; i limiti di generalizzazione vanno dichiarati
come già in C10). Il verdetto C10 della matrice passa da "ordine NON DIMOSTRATA" a
"ordine 2 SOSTENUTO/VERIFICATO sull'intervallo" con la catena a 4 punti.

**Se 2 passa ma 3 fallisce** → NON DIMOSTRATO confermato sulla catena completa a 4
punti: la mancanza di ordine entro h=0.03125 non è più un'ipotesi ma un fatto
misurato su tutta la catena disponibile, e il sospetto si sposta sulla concentrazione
angolare (§C) o sul regime pre-asintotico a T=0.02 — con l'evidenza esatta in mano.

**Se 2 fallisce** → stop della catena di conclusioni; si isola floor vs errore
spaziale (la decomposizione nodale E.8-2, già implementata in `audit_postproc.py` §S9,
dice dove sta la discrepanza).

## 6. Cosa NON fa questo passo (limiti dichiarati)

- Non dimostra l'ordine per altri T, parametri o mesh (resta la generalità da C10).
- Non fornisce la soluzione esatta della PDE (proprietà 4 della tabella §8: resta
  NON VERIFICABILE CON I DATI DISPONIBILI).
- Non dimostra la regolarità/singolarità angolare (C4): per quella servirebbe un
  profilo r^α a r costante sulla famiglia completa — proposta separata, non inclusa.

## 7. Stato

**NON ESEGUITO.** In attesa di autorizzazione esplicita. Stima: ~1–2 min di CPU;
nessuna modifica a codice, input o risultati preesistenti; nuovi artefatti solo nelle
due directory `t1c_n64x32_oofem/` di cui sopra.
