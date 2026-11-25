# Audit finale dei test energetici Allen–Cahn

## Perimetro e verdetto

Audit dei documenti, input, CSV e figure energetiche presenti nel working tree al 26 ottobre 2026. Non sono state eseguite simulazioni né modificate sorgenti o input. Il repository non associa esplicitamente gli artefatti alle etichette “fasi 5.1–5.6”; la tabella seguente le mappa alle sei famiglie di evidenza energetica riconoscibili dai file.

**Conclusione:** le equazioni continue e le formule energetiche sono coerenti con il ramo parabolico e con il ramo iperbolico bulk-inerziale/superficie-parabolica implementati. Le serie ricalcolate più recenti rispettano la somma degli addendi a precisione macchina e mostrano residui che calano nel raffinamento temporale iniziale. Questo **non valida il solver**: la velocità iperbolica non è esportata, alcuni CSV legacy non si riproducono dai loro stessi output, e un grafico confronta quantità con unità diverse. Non è stata trovata una violazione del bilancio imputabile con certezza al solver.

Significato degli stati: **PASS** soddisfa il criterio quantitativo esplicito indicato; **PARTIAL** fornisce evidenza utile ma incompleta; **FAIL** viola un criterio verificabile (specificando se riguarda un artefatto o il solver); **NOT TESTED** non ha dati sufficienti/interpretabili per l'osservabile richiesto.

Classificazione applicata ai singoli controlli, non al solver nel suo complesso:

| Fase / prova | Stato | Criterio verificato e conclusione |
|---|---|---|
| 5.1, legge energetica e mappatura materiale/assemblaggio | **PARTIAL** | Segni e coefficienti di flusso, reazione, capacità e inerzia concordano con le formule. Manca un confronto indipendente del residuo FEM assemblato contro la variazione della funzionale. |
| 5.2, decadimento parabolico e CSV componenti | **PARTIAL** | Le serie affidabili esaminate decrescono; [energy_components_parabolic_example.csv](energy_components_parabolic_example.csv) è però un solo frame e [energy.csv](energy.csv) usa quadratura insufficiente per il potenziale quartico. Nessuno dei due, da solo, verifica il bilancio. |
| 5.3, bilancio parabolico e raffinamento temporale | **PASS** limitato al criterio diagnostico | Nei tre passi 0.005, 0.0025, 0.00125 il residuo relativo cala 0.004231, 0.002289, 0.001201; nel test 16×8 il residuo e la differenza di curva calano quasi di un fattore 2 a ogni dimezzamento. È un PASS del trend di post-processing, non una validazione del solver. |
| 5.4, parabolico 2D, energia e bilancio secante | **PARTIAL** | La somma bulk+superficie è coerente e l'energia cala in 1000 passi. Il residuo usa una velocità secante, non la dissipazione esatta; il pannello superiore di `energy_balance_secant.png` è dimensionalmente errato (**FAIL della figura**, non del run). |
| 5.5, post-processing iperbolico senza velocità | **NOT TESTED** per energia totale | [energy_components_hyperbolic_example.csv](energy_components_hyperbolic_example.csv) ha `kinetic` e `total_energy` vuoti (`missing_velocity`): consente un controllo dei quattro addendi potenziali, non del bilancio iperbolico. |
| 5.5, caso iperbolico singolo 16×8 | **PARTIAL** | Il totale è ricostruito e il residuo relativo finale è circa \(1.80\cdot10^{-7}\), ma la velocità deriva da spostamenti testuali arrotondati, non da un export diretto dello stato interno. |
| 5.6, raffinamento parabolico 16×8 | **PASS** limitato al criterio diagnostico | Gli errori \(L_\infty\) tra curve successive sono 6.624e-5 e 3.343e-5 (ordine indicativo 0.987); il modulo del residuo finale cala a ogni livello. |
| 5.6, raffinamento iperbolico 16×8 | **PARTIAL** | I tre passi 0.002, 0.001, 0.0005 migliorano, ma la serie più fine disponibile e l'incertezza sulla velocità non stabiliscono un regime asintotico robusto del residuo. |
| CSV iperbolico legacy [energy_validation_finer.csv](energy_validation_finer.csv) | **FAIL** della riproducibilità del post-processing | Ricalcolando gli stessi `.out` con la routine energetica corrente, il run 16×8, \(dt=0.00025\), dà \(E(T)=1.9299547624\), \(R(T)=6.87\cdot10^{-7}\), senza incrementi positivi; il CSV riporta 1.9299885203, 3.47e-5 e 17 incrementi positivi. È una discrepanza verificata del dato derivato, non una prova di errore del solver. |

## Formule e implementazione

Ponendo \(F(u)=\tfrac14(u^2-1)^2\), \(F'(u)=u^3-u\), la legge coerente con il codice è:

\[
E_{pot}=\int_\Omega\left(\frac{k}{2}|\nabla u|^2+mF(u)\right)dx
+\int_\Gamma\left(\frac{k}{2}|\nabla_\Gamma\psi|^2+\beta mF(\psi)\right)dS.
\]

Nel parabolico:

\[
\dot E_{pot}=-\int_\Omega c\,u_t^2dx-\int_\Gamma c\,\psi_t^2dS\le0.
\]

Nel ramo `wave 1`, con inerzia solo nel bulk:

\[
E_{tot}=E_{pot}+\int_\Omega\frac{\rho}{2}u_t^2dx,
\qquad
\dot E_{tot}=-\int_\Omega\rho\zeta u_t^2dx-\int_\Gamma c_\Gamma\psi_t^2dS\le0.
\]

Il sorgente `AllenCahnMaterial::computeInternalSourceVector` restituisce \(m(u-u^3)\), moltiplicato per `surfacebeta` sulle linee; `computeFlux3D` restituisce \(-k\nabla u\); `giveCharacteristicValue(Capacity)` usa `c` nel parabolico e `rho` bulk / `surfacecapacity` linea nel ramo iperbolico. `BulkWaveMassAssembler` esclude le linee dall'inerzia. Queste scelte danno i segni della variazione di \(E_{pot}\) e il termine di bordo \(+k\partial_nu\) che si cancella nel bilancio accoppiato. Le formule sono documentate in [energy_theory.md](energy_theory.md) e corrispondono alla sorgente C++ per i casi non-MMS esaminati.

Nei casi di raffinamento il materiale e la sezione hanno coefficienti unitari, spessore/area unitari e la stessa discretizzazione sul contorno; le omissioni di fattori geometrici non alterano quei CSV. I contributi energetici usano Gauss 5×5 sui quadrilateri e 5 punti sulle linee: adeguato al potenziale quartico e ai gradienti sui rettangoli affini in esame. Le dissipazioni usano matrici di massa con quadratura sufficiente per questi elementi.

La legge continua riguarda il modello regolare non forzato e non stabilisce da sola una legge di dissipazione a passo finito. `alpha=1` è il backward Euler dei casi parabolici; Newmark \((\beta_N,\gamma_N)=(1/4,1/2)\) non ha nel codice un'identità discreta dimostrata per il sistema nonlineare accoppiato. Un aumento numerico va confrontato con residuo, passo e accuratezza del post-processing, non dichiarato automaticamente bug o accettato perché l'energia media scende.

## Ramo parabolico

- [energy_balance.csv](energy_balance.csv) documenta tre passi \(0.005,0.0025,0.00125\): il residuo normalizzato diminuisce da 4.231e-3 a 2.289e-3 e 1.201e-3. Il calcolo è un diagnostico a passi finiti basato su differenze secanti e quadratura dichiarata; non è l'identità continua esatta.
- `energy_trajectory.csv` ha 81 righe, additività \(E=E_{bulk}+E_\Gamma\) entro \(2.3\cdot10^{-16}\) e nessun incremento positivo. `energy_trajectory.png` rappresenta queste componenti. Prova di decadimento per quel caso, non verifica del residuo né prova indipendente.
- La serie 16×8 di `energy_refinement_study.csv` ha \(R(T)=-2.947\cdot10^{-4},-1.525\cdot10^{-4},-7.766\cdot10^{-5}\); le differenze successive della curva totale hanno rapporto 1.982 e ordine indicativo 0.987. La scala è coerente con il primo ordine temporale del caso.
- `parabolic_2d/runs/energy_n20_dt0p0001_every_step/energy_interface.csv` contiene 1001 campioni: \(E_{total}=E_{bulk}+E_{surface}\) entro \(1.5\cdot10^{-14}\) e nessun incremento nei dati. `energy_time.png` è coerente con questa serie, con l'ipotesi dichiarata \(\varepsilon=0.1\); non è una verifica indipendente contro il paper.
- `energy_components_parabolic_example.csv` contiene un solo istante: verifica la somma dei componenti per quel frame, non la legge temporale.
- **Errore del post-processing legacy:** `analyze_energy.py` integra il potenziale Q1 bulk con Gauss 2×2. \((u_h^2-1)^2\) è quartico e 2×2 non è esatto in generale; il CSV `energy.csv` è dunque una diagnostica qualitativa, non un valore energetico accurato. Il post-processing più recente usa una quadratura più adeguata.

Nel [CSV secante 2D](parabolic_2d/energy_balance_secant.csv) l'identità discreta `balance_residual_secant = delta_E + dissipation_secant_sum` vale entro \(5.6\cdot10^{-17}\); il massimo \(|R|/|\Delta E|\) è 1.73e-3. Ma le colonne di dissipazione sono integrali per passo, mentre [energy_balance_secant.png](parabolic_2d/energy_balance_secant.png) traccia \(-\Delta E/\Delta t\) come tasso e non divide la dissipazione per \(\Delta t\). Al primo passo le curve comparabili sarebbero circa 4967.6 e 4959.0; la figura attuale mostra invece 4967.6 contro 0.496. Il pannello superiore è quindi fuorviante; il pannello inferiore usa il residuo CSV e resta interpretabile. Inoltre il report 2D collega `parabolic_2d/energy_interface.csv`, percorso inesistente: la serie effettiva è nella cartella `runs/` indicata sopra.

## Ramo iperbolico

La sola energia potenziale può scendere o salire per conversione cinetica. Nel caso singolo la cinetica sale da 0 a 0.00355156 e il potenziale scende di 0.00520321; il totale scende da 1.93160616 a 1.92995451. Questo è coerente con la legge del totale, ma non basta a convalidare il solver.

La relazione Newmark del caso \((\beta_N,\gamma_N)=(1/4,1/2)\) è

\[
v_{n+1}=\frac{2(u_{n+1}-u_n)}{\Delta t}-v_n.
\]

La ricorrenza è algebricamente coerente con lo stato interno per velocità iniziale zero. Tuttavia `Velocity` non è nei VTU; il file `.out` stampa gli spostamenti a circa 8 cifre significative, non a piena precisione come afferma il docstring di `read_text_history`. La divisione per \(\Delta t\) amplifica l'arrotondamento. L'uguaglianza cinematica non equivale quindi a un confronto dei valori ricostruiti con la velocità effettivamente memorizzata dal solver.

| Fonte | \(dt\) | \(E(T)\) | \(R(T)\) | Incrementi positivi del totale |
|---|---:|---:|---:|---:|
| Serie ricalcolata high-precision 16×8 | 0.002 | 1.9299595493 | 3.692e-6 | 0 |
| Serie ricalcolata high-precision 16×8 | 0.001 | 1.9299554140 | 9.156e-7 | 0 |
| Serie ricalcolata high-precision 16×8 | 0.0005 | 1.9299545074 | 3.473e-7 | 0 |
| Ricalcolo dagli output già esistenti 16×8 | 0.00025 | 1.9299547624 | 6.867e-7 | 0 |

Nel triplo 16×8 iniziale il residuo diminuisce, ma al passo 0.00025 non migliora e ha già raggiunto un pavimento confrontabile con l'errore di ricostruzione. Le differenze delle curve tra \(dt\) e \(dt/2\) hanno ordine indicativo 2.19 per il totale e 2.00 per il potenziale; non sono sufficienti per affermare un ordine asintotico del bilancio totale.

**Anomalia di post-processing identificata:** a \(dt=0.0005\), [energy_time_series.csv](energy_time_series.csv) e [energy_validation_hp_hyperbolic_time_series.csv](energy_validation_hp_hyperbolic_time_series.csv) hanno potenziali finali che differiscono solo di circa \(1.4\cdot10^{-8}\), ma energie cinetiche che differiscono di \(1.8135\cdot10^{-5}\). Gli input e gli spostamenti testuali OOFEM a questo passo coincidono; la discrepanza è nel calcolo della cinetica, non una differenza osservata nella soluzione del solver. In `energy_validation_finer.csv` le serie 8×4 riportano 28–81 aumenti per run, massimi fino a \(1.46\cdot10^{-4}\); ricalcolando gli stessi `.out` tali aumenti scompaiono, salvo un incremento residuo di \(5.2\cdot10^{-9}\) al passo 0.000125. Quelli archiviati non sono evidenza affidabile di violazione energetica del solutore.

`energy_validation.png` visualizza il residuo rumoroso della serie legacy, non quello del ricalcolo high-precision più recente; non è stata trovata nel working tree una pipeline attiva che ne registri la provenienza. [energy_refinement_study.png](energy_refinement_study.png) e [energy_hyperbolic_single_test.png](energy_hyperbolic_single_test.png) sono invece coerenti con i rispettivi CSV recenti. Il primo grafico va tenuto come figura legacy/non autorevole finché non viene rigenerato da una sorgente esplicita.

## Coerenza figure e CSV

- **Coerenti nel perimetro controllato:** `energy_trajectory.png` ↔ `energy_trajectory.csv`; `energy_time.png` ↔ il CSV 2D nella cartella `runs/`; `energy_refinement_study.png` ↔ `energy_refinement_study.csv`; `energy_hyperbolic_single_test.png` ↔ il CSV del caso singolo. Le somme energia totale/componenti verificate sui CSV recenti hanno errore massimo dell'ordine di \(2.3\cdot10^{-16}\) (circa \(1.5\cdot10^{-14}\) nel CSV 2D).
- **Non autorevole:** `energy_validation.png` rispecchia l'andamento legacy ma non la serie iperbolica ricalcolata; non usarlo per concludere sul residuo del solver.
- **Figura non coerente dimensionalmente:** pannello superiore di `energy_balance_secant.png`, come sopra. Il CSV secante e il pannello del residuo sono coerenti.
- **Link rotto:** il report `parabolic_2d_validation.md` punta a un `energy_interface.csv` nella directory sbagliata.

## Anomalie residue e verifiche mirate

1. **Esportazione/precisione di `Velocity`.** Aggiungere, in una verifica futura autorizzata, un confronto tra `giveUnknownComponent(VM_Velocity)` e la velocità usata dal post-processing; mantenere alta precisione e verificare l'inizializzazione. Non è stato fatto qui.
2. **Provenienza dei CSV legacy.** Versionare lo script che ha prodotto `energy_time_series.csv` e `energy_validation_finer.csv`; rigenerare i dati solo dopo aver chiarito la formula cinetica e conservare separati raw output e post-processing. Non rigenerato in questo audit.
3. **Quadratura di `analyze_energy.py`.** Sostituire il controllo 2×2 con una regola esatta/di ordine adeguato, o etichettare `energy.csv` come stima sottintegrata; il file non è stato modificato.
4. **Figura secante e link 2D.** Confrontare rate con rate oppure incrementi integrati con dissipazioni integrate; correggere il percorso CSV del report. Nessuna correzione applicata.
5. **Geometrie non rettangolari e sezione.** Il generatore di raffinamento usa quadrilateri rettangolari allineati; `energy_components` richiede un trattamento attento di thickness/area e non valida automaticamente copertura/duplicati delle linee. `run_energy_validation.py` trasforma il gradiente con `inv(J.T)`, formula indistinguibile da `inv(J)` sulle mesh rettangolari testate ma da verificare su quad ruotati/distorti. Il postprocessor assegna area superficiale zero se il parametro manca: verificare che il default non nasconda contributi di bordo negli input futuri.
6. **Energia discreta del solver.** Anche con velocità esportata, servono un residuo variazionale indipendente, prova su quadratura e sensibilità temporale prima di concludere su stabilità/dissipazione discreta Newmark. La verifica di un'identità continua e la monotonia osservata non sostituiscono questi controlli.

La conclusione più forte ammessa è: la formulazione energetica continua è compatibile con i coefficienti assemblati; alcune curve paraboliche e iperboliche sono coerenti con il bilancio al livello di risoluzione studiato; una parte dei grandi aumenti/residui iperbolici archiviati è spiegata dal post-processing della velocità. Restano non testate una verifica indipendente del residuo FEM e una valutazione conclusiva della dissipazione discreta del solver.