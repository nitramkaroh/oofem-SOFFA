# Verifica del ramo iperbolico OOFEM

## Cosa risolve il ramo attuale

Il flag `wave 1` applica Newmark medio ai termini inerziali del bulk. Il materiale Allen--Cahn fornisce la capacità bulk `rho` e la reazione `mobility*(u-u^3)`. Sulle linee di bordo, invece, `surfacecapacity` moltiplica una derivata prima nel tempo: il modello OOFEM corrente è quindi **inerziale nel bulk e parabolico/dinamico al bordo**, non il problema d'onda con inerzia anche sulla superficie studiato dai preprint CRC 1173.

La prima prova quantitativa ha rivelato che inizializzare l'accelerazione a zero produceva errore di avvio del primo ordine. `applyIC` ora calcola l'accelerazione iniziale coerente risolvendo il residuo iniziale sulla matrice d'inerzia del bulk e includendo il contributo di velocità delle linee dinamiche. Il build Debug è stato rigenerato e i test sotto sono stati ripetuti con quel binario.

## Verifiche eseguite

`run_bulk_references.py` crea due problemi bulk-only su una striscia Q1 con Neumann naturale, senza termini sorgente fabbricati:

1. **Modo coseno lineare.** Si pone `mobility=0`, `u(x,0)=cos(pi*x)`, `u_t(x,0)=0`, `rho=k=1`. Il confronto è con la soluzione esatta del sistema FEM semidiscreto Q1 a massa consistente; così l'errore di tempo è separato dall'errore spaziale. Per `dt = 0.04, 0.02, 0.01`, l'errore RMS scende da `3.47e-4` a `8.70e-5` e `2.18e-5`; l'ordine osservato al dimezzamento è 2.00.
2. **Oscillatore Allen--Cahn non lineare.** Si usa un campo spazialmente uniforme `u(0)=0.2`, `u_t(0)=0.3`, `rho=mobility=1`. Il sistema FEM si riduce a `u_tt=u-u^3`, confrontato con un'integrazione indipendente DOP853 dell'ODE. Per gli stessi passi l'errore assoluto finale è `7.31e-6`, `1.81e-6`, `5.06e-7`; il primo ordine osservato è 2.02. I valori del terzo passo iniziano a risentire dell'arrotondamento dell'output VTU.

`run_stationary_kink.py` verifica inoltre il profilo stazionario `u(x)=tanh(x/sqrt(2))`, la forma kink riportata nella letteratura sull'Allen--Cahn iperbolico. Su `[-8,8]` con Neumann naturale approssimato agli estremi, `T=0.2` e `dt=0.001`, l'errore RMS rispetto al profilo esatto è `3.10e-5`, `7.71e-6`, `2.08e-6` sulle mesh `Nx=40,80,160`; gli ordini spaziali osservati sono 2.01 e 1.89. Il secondo è influenzato dalla troncatura ai bordi e dalla precisione VTU.

`run_coupled_refinement.py` ripete il caso bulk inerziale + bordo Allen--Cahn dinamico 2D già archiviato in `../allencahn_wave_dynamic.in`, a `T=0.02`. Il confronto dei campi finali per `dt = 0.002, 0.001, 0.0005` dà ordine temporale 1.95, vicino all'ordine due atteso; il riferimento è la differenza fra due passi, non una soluzione esatta. Tutti i run terminano senza errori OOFEM.

Metriche, input, log e VTU sono conservati nelle sottocartelle `linear_dt*`, `oscillator_dt*`, `kink_n*` e `coupled_dt*`; CSV riassuntivi: `bulk_reference_metrics.csv`, `stationary_kink_space_refinement.csv` e `coupled_time_refinement.csv`.

## Che cosa dicono le nuove fonti

- Hipp–Kovács, *Finite element error analysis of wave equations with dynamic boundary conditions: L2 estimates* (CRC 1173 Preprint 2018/53), fornisce modelli d'onda con inerzia anche al bordo (`mu*u_tt`) e test numerici di convergenza. È un riferimento per FEM d'onda bulk-superficie, ma la legge temporale al bordo non coincide con quella implementata qui. [PDF](</home/marco/Documents/Work/Riccardo/refs/CRC1173_Preprint_2018-53.pdf>).
- Hochbruck–Leibold, *Finite element discretization of semilinear acoustic wave equations with kinetic boundary conditions* (CRC 1173 Preprint 2019/26), studia `u_tt` nel bulk e sulla frontiera, con termini di smorzamento e reazioni semilineari specifiche; il test su disco unitario usa una soluzione esatta fabbricata. Il modello e la geometria non sono una replica diretta di Allen--Cahn OOFEM. [PDF](</home/marco/Documents/Work/Riccardo/refs/CRC1173_Preprint_2019-26.pdf>).
- Hochbruck–Leibold, *An implicit-explicit time discretization scheme for second-order semilinear wave equations with application to dynamic boundary conditions* (CRC 1173 Preprint 2020/20), tratta analogamente un'onda acustica con condizione cinetica e test temporali contro una soluzione esatta. [PDF](</home/marco/Documents/Work/Riccardo/refs/CRC1173_Preprint_2020-20 (00000003).pdf>).
- Nizovtseva, Galenko e Alexandrov, *The hyperbolic Allen--Cahn equation: exact solutions*, J. Phys. A 49 (2016), 435201, riporta profili kink/tanh viaggianti ([DOI](https://doi.org/10.1088/1751-8113/49/43/435201)). Il controllo kink qui usa il limite statico, con velocità zero, della forma tanh; verifica l'Allen--Cahn iperbolico nel bulk, ma non l'accoppiamento specifico con bordo dinamico OOFEM.

**Esito circoscritto:** il percorso Newmark e la reazione iperbolica bulk sono verificati quantitativamente e il caso accoppiato mostra convergenza temporale di ordine circa due dopo la correzione dell'accelerazione iniziale. Le fonti CRC non validano direttamente l'intero modello misto bulk-inerziale/bordo-parabolico. Per poter dire che l'intero problema iperbolico OOFEM è validato contro letteratura serve un riferimento con la stessa legge di bordo o un benchmark numerico pubblicato comparabile; questi tre preprint, da soli, non lo forniscono.

## Replica dell'onda viaggiante di Nizovtseva et al. (2016)

L'audit aggiornato del riferimento, della velocità iniziale, del dominio troncato e della convergenza temporale/spaziale è in [`hyperbolic_bulk_validation.md`](hyperbolic_bulk_validation.md). Per questa valutazione dettagliata fa fede quel report e il nuovo runner `validate_nizovtseva_bulk.py`; la sequenza temporale storica qui sotto non va usata da sola per stimare l'ordine, perché il precedente confronto era vicino al pavimento spaziale.

`run_nizovtseva_traveling_wave.py` confronta direttamente il bulk Allen--Cahn iperbolico con la soluzione chiusa dell'Eq. (29), scegliendo `tau=0.2`, `D=M=1`, dominio `[-12,12]x[0,1]`, e velocità iniziale esatta `u_t=1.5*u*(1-u)`. Il damping input `wavedamping=5` realizza il termine fisico `u_t` dopo aver diviso l'equazione per `rho=tau=0.2`. Non è MMS: si confronta la traiettoria con la formula analitica e la condizione di Neumann naturale incide solo sulle code esponenzialmente piccole al bordo troncato.

Gli errori RMS a `T=0.2` sulle griglie `h=0.30, 0.15, 0.075` sono rispettivamente `2.9183e-5`, `7.3231e-6`, `1.8621e-6`, con ordini spaziali `1.995` e `1.976`. Questo è un confronto quantitativo positivo per l'equazione bulk e per l'inizializzazione della velocità; non verifica l'equazione dinamica di superficie.

## Profilo esatto di Hochbruck--Leibold: prova adattata

`run_hochbruck_leibold_disk.py` usa il profilo pubblicato `u=sin(2*pi*t)*x*y` e i suoi dati iniziali su un disco unitario triangolato. Perché il campo sia soluzione del problema OOFEM, la forzante è ricalcolata per l'operatore OOFEM (inerzia bulk e legge Allen--Cahn di primo ordine al bordo). È quindi una verifica analitica adattata e non la replica del benchmark dell'articolo: Hochbruck--Leibold usano inerzia al bordo e termini bulk ulteriori.

A `T=0.2`, sulle mesh da 60, 179 e 502 nodi gli errori RMS sono `2.121e-2`, `1.687e-2` e `1.424e-2`; il massimo resta circa `4.9e-2`. Il raffinamento spaziale è insufficiente e dominato dall'approssimazione poligonale del bordo e dalle derivate tangenziali ai vertici. Il caso resta esplorativo e non è evidenza di validazione del benchmark H--L. Una replica diretta richiede l'implementazione dell'equazione acustica di secondo ordine sulla superficie e dei suoi termini bulk.

Metriche in `nizovtseva_traveling_wave_metrics.csv` e `hochbruck_leibold_disk_metrics.csv`; figure riassuntive in `../comparison_report_figures/`. Le due simulazioni si rigenerano con i rispettivi script `run_*.py`.
