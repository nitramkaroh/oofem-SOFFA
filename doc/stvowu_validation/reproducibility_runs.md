# Ripetizioni operative Allen–Cahn OOFEM

## Ambito

Sono state effettuate tre esecuzioni sequenziali, tutte in directory nuove sotto `/tmp/oofem_repro_audit_20261026/`. Nessun run archiviato è stato sovrascritto; gli input copiati sono byte-per-byte identici agli originali. Binario usato: `build/debug_magneto/oofem` del working tree corrente. Ogni processo ha terminato con `ANALYSIS FINISHED`, 0 errori e 0 warning. Non sono stati eseguiti raffinamenti, scansioni o altre simulazioni.

Le directory temporanee seguenti conservano input, log e output completi per la sessione corrente; per archiviazione durevole vanno copiate in un archivio versionato prima della pulizia automatica del sistema.

## Test 1 — Parabolico 2D accoppiato

**Scelta e costo stimato.** Caso già eseguito `space_n20_dt0p0005`: Allen–Cahn parabolica bulk-superficie Q1, 20×20 quad, 80 elementi linea, 441 nodi/DOF, 200 passi, \(dt=0.0005\), \(T=0.1\), 0 BC imposte e bordo dinamico completo. Log storico: circa 19 s. Stima prima del run: decine di secondi; memoria modesta, mesh da 441 DOF.

**Preflight.** Destinazione `/tmp/oofem_repro_audit_20261026/parabolic_n20_dt0p0005` verificata assente. L’output cadence dell’input è ogni 200 passi, quindi un VTU finale oltre a PVD e testo.

**Comando eseguito.**

```sh
mkdir -p /tmp/oofem_repro_audit_20261026/parabolic_n20_dt0p0005
cp doc/stvowu_validation/parabolic_2d/runs/space_n20_dt0p0005/case.in \
  /tmp/oofem_repro_audit_20261026/parabolic_n20_dt0p0005/case.in
cd /tmp/oofem_repro_audit_20261026/parabolic_n20_dt0p0005
/home/marco/Documents/Work/oofem-SOFFA/build/debug_magneto/oofem -f case.in > run.log 2>&1
```

**Esito e metriche.** Log conservato in `/tmp/oofem_repro_audit_20261026/parabolic_n20_dt0p0005/run.log`; `ANALYSIS FINISHED`, 0 errori/0 warning, tempo reale 19 s. Creati `case.in`, `run.log`, `case.out`, `case.out.m0.pvd`, `case.out.m0.200.vtu`. Il campo/coordinate al passo finale coincidono esattamente col VTU archiviato (RMS e massimo della differenza nodale: 0). L’analisi energetica sul nuovo VTU dà:

| Quantità | Archiviata | Ripetizione |
|---|---:|---:|
| \(E_{bulk}(T)\) | 21.297949516661088 | 21.297949516661088 |
| \(E_{surface}(T)\) | 28.416249050730585 | 28.416249050730585 |
| \(E_{total}(T)\) | 49.714198567391676 | 49.714198567391676 |

Gli scarti spaziali rispetto al riferimento numerico N=40 nel CSV, RMS bulk 0.0799777541, RMS traccia 0.0896266004, massimo bulk 0.25815093, restano identici perché il nuovo campo N=20 è identico all’archiviato. N=40 è un riferimento di griglia, non soluzione esatta. \(\varepsilon=0.1\) è un valore esplorativo assunto per Xiao–Chen, non specificato nel paper; il run è quindi una riproduzione del caso OOFEM, non una replica parametrica esatta dell’articolo.

**Esito di riproducibilità: riuscito con riserva sul riferimento.** Input, comando e metriche sono rigenerabili; il valore di \(\varepsilon\) limita la comparabilità bibliografica.

## Test 2 — Bulk iperbolico Nizovtseva

**Scelta e costo stimato.** Input esistente `nizovtseva_bulk_v2_space_n160`: 322 nodi/DOF, 160 quad Q1, \(dt=0.001\), 200 passi fino a \(T=0.2\), `rho=0.2`, `wavedamping=5`, mobilità/diffusione 1. Log storico: circa 12 s. Stima: circa 10–15 s e output contenuto (200 VTU).

**Preflight.** Destinazione `/tmp/oofem_repro_audit_20261026/hyperbolic_nizovtseva_n160` assente; input e riferimento analitico disponibili. Non serve generare una suite o variare parametri.

**Comando eseguito.**

```sh
mkdir -p /tmp/oofem_repro_audit_20261026/hyperbolic_nizovtseva_n160
cp doc/stvowu_validation/hyperbolic/nizovtseva_bulk_v2_space_n160/case.in \
  /tmp/oofem_repro_audit_20261026/hyperbolic_nizovtseva_n160/case.in
cd /tmp/oofem_repro_audit_20261026/hyperbolic_nizovtseva_n160
/home/marco/Documents/Work/oofem-SOFFA/build/debug_magneto/oofem -f case.in > run.log 2>&1
```

**Esito e metriche.** Log conservato in `/tmp/oofem_repro_audit_20261026/hyperbolic_nizovtseva_n160/run.log`; `ANALYSIS FINISHED`, 0 errori/0 warning, tempo reale 12 s. Creati `case.in`, `run.log`, `case.out`, PVD e 200 VTU (tempi 0.001–0.2). L’input è identico a quello archivato e il campo finale coincide bit-per-bit (RMS/max differenza vs raw originale: 0). Ricalcolo indipendente sul campo finale, collassando le due righe y identiche e usando l’Eq. (29) ai nodi x unici:

| Metrica a \(T=0.2\) | CSV archiviato | Nuova esecuzione |
|---|---:|---:|
| RMS nodale | 7.323071979296714e-6 | 7.323071979296714e-6 |
| Massimo nodale | 2.888978006165832e-5 | 2.888978006165832e-5 |

**Esito di riproducibilità: riuscito.** Il controllo riguarda il bulk iperbolico e non la dinamica di bordo.

## Test 3 — Sistema iperbolico bulk-inerziale / bordo dinamico

**Scelta e costo stimato.** Caso esistente `coupled_independent_time_n16x8_dt0.0005`: 153 nodi/DOF, 128 quad + 48 linee, 40 passi a \(dt=0.0005\), \(T=0.02\), Newmark \((1/4,1/2)\), nessun BC. Log storico: circa 1 s; stima inferiore a 5 s, memoria modesta. Il riferimento `independent_reference.npz` è presente.

**Preflight.** Destinazione `/tmp/oofem_repro_audit_20261026/hyperbolic_coupled_n16x8_dt0p0005` assente. Il run produce 40 VTU e PVD; non viene rilanciata né ricostruita l’integrazione DOP853 del riferimento, che è già archiviato.

**Comando eseguito.**

```sh
mkdir -p /tmp/oofem_repro_audit_20261026/hyperbolic_coupled_n16x8_dt0p0005
cp doc/stvowu_validation/hyperbolic/coupled_independent_time_n16x8_dt0.0005/case.in \
  /tmp/oofem_repro_audit_20261026/hyperbolic_coupled_n16x8_dt0p0005/case.in
cd /tmp/oofem_repro_audit_20261026/hyperbolic_coupled_n16x8_dt0p0005
/home/marco/Documents/Work/oofem-SOFFA/build/debug_magneto/oofem -f case.in > run.log 2>&1
```

**Esito e metriche.** Log conservato in `/tmp/oofem_repro_audit_20261026/hyperbolic_coupled_n16x8_dt0p0005/run.log`; `ANALYSIS FINISHED`, 0 errori/0 warning, tempo reale 1 s. Creati input, log, `.out`, PVD e 40 VTU. Il campo finale coincide esattamente col raw archiviato (RMS/max differenza nodale: 0), ma il confronto contro il riferimento ODE indipendente riproduce la discrepanza già segnalata:

| Metrica a \(dt=0.0005\), \(T=0.02\) | CSV archiviato | Nuova esecuzione |
|---|---:|---:|
| RMS bulk FEM vs ODE | 4.697547032139687e-4 | 4.697547032139687e-4 |
| RMS superficie FEM vs ODE | 2.107343741036224e-3 | 2.107343741036224e-3 |
| Massimo DOF bulk/superficie | 3.45184784626007e-3 | 3.45184784626007e-3 |

Il valore non migliora rispetto agli altri passi già archiviati; il run conferma la riproducibilità del dato, non la correttezza della discrepanza né la convergenza al riferimento. Non è attribuito automaticamente al solver.

**Esito di riproducibilità: riuscito operativamente, anomalia scientifica non risolta.** Questo era il terzo e ultimo run; non sono state avviate altre verifiche.

## Test non eseguiti

- Nessun altro benchmark parabolico: la ripetizione selezionata copre già il sistema accoppiato 2D e dispone di input, runner e metriche archiviate.
- Nessun raffinamento temporale/spaziale: non autorizzato dalla selezione e non necessario per verificare l’identità del run.
- Nessun test isolato di superficie iperbolica: nel modello `wave 1` la superficie è del primo ordine; non esiste un benchmark separato adatto.
- Nessuna ulteriore suite accoppiata: completato un caso a basso costo; il risultato ha riprodotto l’anomalia ODE/OOFEM, quindi mi fermo come richiesto.

## Versione e limiti

I tre run usano lo stesso eseguibile `build/debug_magneto/oofem` presente nel workspace; gli output identificano un build Debug del working tree, non una release pulita immutabile. I log e i raw delle ripetizioni sono stati conservati in `/tmp/oofem_repro_audit_20261026/`; l’input di ogni run è stato copiato senza modifiche e verificato byte-per-byte contro la sorgente.

Questa procedura dimostra che, con questo eseguibile e questi input, i risultati archiviati si ripetono. Non sostituisce la fissazione di un commit/build pulito né risolve la discrepanza del riferimento indipendente per il sistema accoppiato.