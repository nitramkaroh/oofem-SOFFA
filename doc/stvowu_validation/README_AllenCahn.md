# Guida di riproduzione — Allen–Cahn OOFEM

Questa guida riguarda le prove Allen–Cahn OOFEM descritte in [`OOFEM_AllenCahn_Validation_Report.md`](OOFEM_AllenCahn_Validation_Report.md). La [matrice](validation_matrix_final.md) riassume la copertura e i limiti per componente. Le istruzioni sotto separano test bulk, sistemi con bordo dinamico, analisi energetiche e confronti con riferimenti.

## Ambiente

- Repository alla radice `oofem-SOFFA`.
- Eseguibile atteso dai runner: `build/debug_magneto/oofem`; alcuni test storici usano `build/opt_magneto/oofem`. Gli script non selezionano un binario comune né lo ricompilano.
- Gli header `.out` registrano la revisione incorporata. I risultati esaminati provengono da un working tree modificato (`rel-2.1-5430-g2cc1a7706-dirty`); per una riproduzione scientifica pubblicabile annotare commit pulito, opzioni di compilazione e versione esatta dell’eseguibile.
- Dipendenze Python usate dagli analizzatori: Python 3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3. Gli script energetici usano anche `xml.etree.ElementTree` della libreria standard. Il benchmark ODE usa `scipy.integrate.solve_ivp` con DOP853.
- Il confronto Okumura digitalizzato richiede `pdftoppm`, Pillow e il PDF della tesi citata nel report. Il file PDF è esterno al repository; occorre ottenerlo legittimamente e passare il suo percorso al digitalizzatore.

Non è presente un lockfile Python che fissi le dipendenze. Verificare la presenza del binario e registrare le versioni prima di confrontare numeri fra macchine.

## Test bulk iperbolici

**Modo coseno e oscillatore uniforme**

```sh
python3 doc/stvowu_validation/hyperbolic/run_bulk_references.py
```

Il runner genera input, log e VTU nei folder `hyperbolic/linear_dt*` e `hyperbolic/oscillator_dt*`, poi riscrive `bulk_reference_metrics.csv`. Usa la soluzione Q1 semidiscreta esatta per il modo lineare e un’integrazione DOP853 per l’ODE uniforme. Per conservare gli artefatti esistenti, eseguire in una copia di lavoro pulita o salvarne una copia prima: il runner riusa le directory e sovrascrive file.

**Kink stazionario**

```sh
python3 doc/stvowu_validation/hyperbolic/run_stationary_kink.py
```

Genera i tre input e risultati in `kink_n40`, `kink_n80`, `kink_n160` e riscrive `stationary_kink_space_refinement.csv`. Il riferimento è il profilo analitico tanh; il dominio finito applica Neumann naturale approssimata agli estremi. Eseguire solo con una copia dei risultati se si desidera conservarli.

**Onda viaggiante Nizovtseva**

```sh
python3 doc/stvowu_validation/hyperbolic/validate_nizovtseva_bulk.py --suite all
python3 doc/stvowu_validation/hyperbolic/analyze_nizovtseva_bulk.py
```

La prima riga esegue la suite spaziale/temporale e i controlli di dominio/dato iniziale; la seconda ricalcola CSV e figure dai PVD/VTU. Il riferimento è la soluzione analitica dell’articolo. Il runner scrive in `hyperbolic/nizovtseva_bulk_v2_*` e aggiorna i CSV/PNG nella cartella: usare una copia del repository o archiviare prima i file correnti. Per rigenerare solo un singolo caso, usare gli argomenti CLI descritti da `--help` (senza lanciare `--suite all`).

## Bulk e bordo dinamico parabolico

**Caso 2D singolo** — il runner genera ed esegue un input, con cartella di output esplicita:

```sh
python3 doc/stvowu_validation/parabolic_2d/run_case.py \
  --n 20 --dt 0.0005 --end-time 0.1 --eps 0.1 --output-every 200 \
  --out /percorso/nuovo/parabolic_n20_dt0p0005
```

Per un confronto col caso archivato, usare gli stessi parametri descritti in `parabolic_2d_validation.md`. La scelta \(\varepsilon=0.1\) del caso Xiao–Chen è esplorativa e non è un parametro dichiarato nel paper.

Per metriche e figure della suite, il comando previsto è:

```sh
python3 doc/stvowu_validation/parabolic_2d/analyze.py --root /percorso/cartella/runs
```

L’analizzatore aggiorna CSV e PNG nella directory dello script e può generare/riscrivere la serie energetica secante. Controllare che la cartella passata contenga i run attesi e salvare prima una copia dei risultati correnti. Il report usa riferimenti N=40/N=60 e passi più fini, non soluzioni esatte.

## Bilanci energetici

**Bilancio 20×20 a tre passi**

```sh
python3 doc/stvowu_validation/run_energy_balance.py
python3 doc/stvowu_validation/update_report_assets.py
```

Il primo script esegue tre simulazioni e scrive raw in `/tmp/stvowu_energy_balance`; riscrive inoltre `energy_balance.csv`. Il secondo dipende da quelle directory temporanee e rigenera traiettoria/figure. La riproduzione è quindi specifica per macchina se non si conserva `/tmp`; non lanciare il secondo comando senza avere gli output corrispondenti. Per mantenere risultati esistenti, usare una copia di lavoro isolata.

**Raffinamento energetico 16×8**

```sh
python3 doc/stvowu_validation/run_energy_validation.py
python3 doc/stvowu_validation/analyze_energy_validation.py
```

Il runner principale esegue dieci configurazioni (raffinamenti temporali e mesh) e sovrascrive input/raw e CSV/PNG di riepilogo. Il suo default scrive in `energy_validation_runs/`; l'analizzatore usa invece prevalentemente `energy_validation_runs_highprecision/` (con un caso in `energy_validation_runs/`). **Non sono quindi una catena run→analisi univoca**: concatenarli non ricalcola necessariamente i CSV dal run appena eseguito. Scegliere e annotare esplicitamente la cartella dati per ogni caso ed eseguire soltanto in una copia/checkout dedicato. Le variabili `OOFEM_ENERGY_OUT` e `OOFEM_ENERGY_PREFIX` cambiano il runner, ma non vengono applicate automaticamente all'analizzatore. I CSV storici iperbolici non sono tutti coerenti con il ricalcolo dai raw; consultare `validation_metrics_audit.md` prima di usare i dati.

Il post-processor [`postprocess_energy.py`](postprocess_energy.py) è non distruttivo rispetto al CSV di destinazione (rifiuta overwrite) e può essere applicato a input compatibili. Per il ramo iperbolico i VTU correnti non esportano `Velocity`, quindi fornisce i contributi potenziali ma non cinetica o energia totale.

## Ripetizioni operative archiviate

Le ripetizioni già documentate e verificate sono descritte in [`reproducibility_runs.md`](reproducibility_runs.md). Per un nuovo run manuale, copiare l’input esistente in una directory nuova, sostituire il prefisso output solo se necessario e invocare:

```sh
/path/to/oofem -f case.in
```

Mantenere nella stessa directory input, log, `.out`, PVD e VTU. Non riutilizzare un folder di output se non si intende sovrascrivere i risultati.

## Test regressivi e limiti

Gli input `tests/tm/allencahn*.in` sono regressioni di esecuzione; da soli non calcolano errori quantitativi o bilanci energetici. Per i confronti numerici usare i runner e i riferimenti elencati nel rapporto finale.

Alcuni benchmark letterari dipendono da PDF/figure esterni o da raw output precedentemente collocati sotto `/tmp`; i CSV aggregati non bastano sempre a ricostruire errori e norme. Registrare percorso e checksum dei riferimenti esterni, input effettivo, commit/build e script di post-processing con ogni nuova esecuzione.
