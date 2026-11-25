# Validazione parabolica 2D Allen–Cahn con bordo dinamico

**Esito:** il test conferma la coerenza dell'assemblaggio bulk-superficie per i coefficienti unitari del caso, la convergenza temporale di primo ordine e l'invarianza rispetto a inversione/riordinamento delle linee di bordo. La convergenza spaziale è supportata dalla sequenza fino a 60×60 a (T=0{,}0625), ma il caso a (T=0{,}1) resta sensibile alla mesh. Il confronto con Xiao–Chen è soltanto qualitativo: il paper non dichiara (ε) per il Caso 4 e non rende disponibili dati numerici tabulati.

Nessuna sorgente artificiale o MMS è stata usata; nessuna modifica al codice C++ è stata fatta; nessun test WED è incluso.

## 1. Audit del riferimento Xiao–Chen

Riferimento: M. Xiao, R. Chen, “A first-order energy stable scheme for the Allen–Cahn equation with the Allen–Cahn type dynamic boundary condition”, *Journal of Computational and Applied Mathematics* 460 (2025), 116409, [DOI 10.1016/j.cam.2024.116409](https://doi.org/10.1016/j.cam.2024.116409). È stato esaminato il PDF originale fornito, in particolare le equazioni (2.7)–(2.10), (2.15)–(2.20), i parametri (5.1) e il Caso 4/Fig. 5.18–5.21.

| Voce | Quanto specifica il riferimento |
|---|---|
| Equazione | φ_t = −γμ nel bulk; μ = −Δφ + F'(φ); φ|_Γ = ψ; ψ_t = −Mμ_Γ; μ_Γ = −Δ_Γψ + G'(ψ) + ∂_nφ. Il Caso 4 mostrato in Fig. 5.19 è il ramo **senza** i due termini globali di penalizzazione della massa. |
| Potenziale | Il doppio pozzo classico è (F(z)=\frac{1}{4\varepsilon^2}(z^2-1)^2), con derivata ((z^3-z)/\varepsilon^2). Il modello ammette anche una scelta distinta di (G); il Caso 4 non ristampa una formula separata per (G). L'input esplorativo qui usa il doppio pozzo anche sulla superficie, coerente con la scelta classica descritta nell'articolo. |
| Coefficienti pubblicati | I valori predefiniti sono γ = M = 0,01, τ = 10⁻⁴, h = 1/300; i coefficienti di penalizzazione valgono solo per i rami penalizzati. |
| Parametro d'interfaccia | **ε non è specificato per il Caso 4** né nell'elenco (5.1), né nel testo del Caso 4, né nella didascalia di Fig. 5.19. Non gli attribuisco il valore usato in OOFEM. |
| Dominio e dato iniziale | Ω = [0,1]². Il testo colloca una goccia quadrata di lato 0,5 e centro (0,5; 0,25); Fig. 5.18 mostra valore +1 nella goccia a contatto col lato inferiore e −1 all'esterno. L'immagine non fornisce una formula nodale più dettagliata per i punti sugli spigoli. |
| Bordo e tempi | Condizione Allen–Cahn dinamica sull'intera frontiera; Fig. 5.19 riporta T = 0,1; 0,5; 1; 2; 4; 6. |
| Metodo e dati | Differenze finite centrate di ordine spaziale due e schema stabilizzato di primo ordine nel tempo. Per il Caso 4 non sono pubblicati campi nodali, tabelle di energia o coordinate dell'interfaccia: sono disponibili figure raster, non dati quantitativi grezzi. |

**Confronto consentito:** morfologia e ordine di grandezza della forma a T = 0,1, dichiarando che OOFEM assume ε = 0,1. **Non consentito:** errore quantitativo sul campo, energia o posizione dell'interfaccia attribuito al paper, finché ε e i dati numerici di riferimento non sono noti. La figura digitalizzata non è stata convertita in misure: una lettura pixel per pixel non risolverebbe l'incertezza sul parametro mancante.

## 2. Equazione implementata e mapping

Con il doppio pozzo (F=G=(z^2-1)^2/(4\varepsilon^2)), il ramo senza penalizzazione del paper equivale a

\[
\begin{aligned}
\phi_t &= \gamma\Delta\phi + \frac{\gamma}{\varepsilon^2}(\phi-\phi^3) &&\text{in }\Omega,\\
\psi_t &= M\Delta_\Gamma\psi + \frac{M}{\varepsilon^2}(\psi-\psi^3)-M\partial_n\phi &&\text{su }\Gamma,\\
\psi&=\phi|_\Gamma.
\end{aligned}
\]

Per γ = M = 0,01, ε = 0,1 (ipotesi esplorativa), i coefficienti OOFEM sono `c=1`, `k=0.01`, `mobility=1`, `surfacebeta=1` e capacità superficiale 1. Il materiale usa lo stesso `k` sui quad e sulle linee: coincide con i coefficienti bulk e tangenziale perché in questo riferimento γ = M. Il termine normale non è assemblato come carico separato. Nella forma debole combinata il flusso normale bulk e il termine (M\partial_n\phi) superficiale si cancellano proprio quando i coefficienti sono uguali, come in questo caso.

`allencahnmat` restituisce flusso (-k\nabla u) e sorgente (m(u-u^3)); `TransportElement` assembla diffusione con segno positivo e sottrae la sorgente al residuo. Le linee aggiungono massa superficiale, diffusione tangenziale e reazione al sistema globale. La variabile superficiale non ha gradi di libertà distinti: usa gli stessi DOF nodali del bulk e coincide con la traccia. Di conseguenza, l'errore di traccia coincide con quello del campo al bordo; non esiste un errore indipendente su una seconda variabile ψ.

Il campo bulk è discretizzato con quadrilateri bilineari Q1, 9 punti di Gauss; ogni lato del quadrato è discretizzato una volta con segmenti lineari e 3 punti di Gauss. Questa quadratura integra esattamente, su elementi rettangolari regolari, i polinomi della reazione cubica Q1 nel residuo/tangente. Il generatore controlla che le 4N linee coprano ogni lato una volta sola. Non ci sono Dirichlet o sorgenti: tutte le linee rappresentano la frontiera dinamica del modello. Il passo temporale OOFEM è backward Euler (BDF1) totalmente implicito; la soluzione non lineare è ottenuta con Newton modificato, aggiornando la tangente nella prima iterazione del passo.

La verifica del codice ha riguardato [allencahnmat.C](../../../src/tm/Materials/allencahnmat.C) (`computeFlux3D`, `giveCharacteristicValue`, `computeInternalSourceVector`), [transportelement.C](../../../src/tm/Elements/transportelement.C) (residuo, matrice di capacità e quadratura), [nltransienttransportproblem.C](../../../src/tm/EngineeringModels/nltransienttransportproblem.C) e [nonstationarytransportproblem.C](../../../src/tm/EngineeringModels/nonstationarytransportproblem.C) (BDF/Newton e integrazione temporale), oltre a [quad1_ht.C](../../../src/tm/Elements/quad1_ht.C) e [line1_ht.C](../../../src/tm/Elements/line1_ht.C). Non è stato identificato in questo caso un errore di assemblaggio che giustifichi una patch; non ho modificato questi sorgenti.

## 3. Test eseguiti e metriche

Gli input contengono `manufactured` disattivato per default, non hanno carichi o sorgenti imposti e non chiamano alcun percorso MMS. Le esecuzioni completate hanno raggiunto il tempo finale e terminano con `ANALYSIS FINISHED`; i risultati grezzi sono conservati sotto `runs/`.

### 3.1 Confronto con Fig. 5.19

Il run OOFEM su 20×20, Δt = 10⁻⁴, T = 0,1 usa ε = 0,1 **assunto**. Il confronto [figura a lato](xiao_chen_qualitative_T0p1.png) mostra la stessa goccia appoggiata al bordo inferiore, con appiattimento iniziale del vertice e interfaccia diffusa. È **evidenza qualitativa (E)**, non una replica parametrica o una validazione quantitativa. La goccia iniziale è rappresentata assegnando ±1 ai nodi secondo la figura; il suo interpolante Q1 smussa il salto su circa un elemento, mentre il paper usa differenze finite con h = 1/300.

### 3.2 Raffinamento spaziale a T = 0,1

Tre mesh Q1, N×N, con identici coefficienti e Δt = 5·10⁻⁴. Si confronta il campo ricostruito su una griglia comune con N = 40 (riferimento numerico più fine, classe D). L'interfaccia è definita internamente come livello u = 0; (y_{top}) è l'intersezione con x = 0,5, mentre (x_L,x_R) sono le intersezioni col bordo inferiore.

| N | L2 bulk vs N=40 | L2 traccia vs N=40 | max bulk | E_bulk(T) | E_surface(T) | E_total(T) | (y_{top}) | (x_L,x_R) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 0,07998 | 0,08963 | 0,2582 | 21,29795 | 28,41625 | 49,71420 | 0,5250 | 0,2250; 0,7750 |
| 30 | 0,06211 | 0,08861 | 0,2246 | 20,75819 | 28,91298 | 49,67117 | 0,5167 | 0,2500; 0,7500 |
| 40 | 0 | 0 | 0 | 21,18762 | 29,14523 | 50,33284 | 0,5125 | 0,2375; 0,7625 |

Il bulk migliora passando 20→30, ma la norma L2 della traccia quasi non cambia (0,08963→0,08861) e il massimo resta elevato. Anche le energie variano di circa 1,2–1,3% nel totale e di fino a 2,5% nel contributo superficiale. **Questa sequenza a T=0,1 non dimostra robustezza spaziale.**

Il campo supera leggermente l'intervallo iniziale [−1,1] sulla mesh più grossolana: a T=0,1, min/max sono −1,01099/1,02186 per N=20, −1,00144/1,00287 per N=30 e −1,00021/1,00042 per N=40. Nella sequenza a T=0,0625 il sovraelongamento cala da circa 2,1% (N=20) a circa 0,01% (N=60). Lo registro come effetto di discretizzazione/risoluzione; non lo trasformo in un test del principio del massimo, che non è stato dimostrato per questo schema accoppiato.

Per distinguere questo esito dalla diversa risoluzione iniziale, è stata eseguita una seconda sequenza coerente, con identico Δt = 6,25·10⁻⁴ e T = 0,0625 fino a 60×60. Rispetto a N=60, la L2 bulk/traccia è rispettivamente 0,1173/0,1320 (N=20), 0,05461/0,06707 (N=30), 0,02994/0,03387 (N=40), 0 (N=60). Gli scarti diminuiscono regolarmente. E_total a T = 0,0625 è 58,4632; 58,9530; 60,0641; 60,4805 per N = 20, 30, 40, 60. Questo **supporta convergenza numerica spaziale (D)**, ma non è un ordine asintotico certificato e non sostituisce un riferimento indipendente.

Le metriche e i campi sono in [`spatial_refinement.csv`](spatial_refinement.csv), [`spatial_refinement_T0p0625.csv`](spatial_refinement_T0p0625.csv), [`spatial_fields_T0p1.png`](spatial_fields_T0p1.png) e [`spatial_fields_T0p0625.png`](spatial_fields_T0p0625.png). Le differenze L2 spaziali sono valutate su una griglia cartesiana comune con quadratura trapezoidale; sono norme della soluzione FE ricostruita e della traccia, non errori rispetto alla soluzione esatta.

### 3.3 Raffinamento temporale

Con N = 20 fissato e T = 0,1, i campi finali a Δt = 4·10⁻⁴, 2·10⁻⁴, 10⁻⁴ danno differenze L2 bulk rispetto al run più fine di 1,290·10⁻⁴ e 4,303·10⁻⁵. Le differenze tra coppie adiacenti hanno rapporto 1,999 e ordine apparente 0,9994, coerente con backward Euler. Le differenze L2 sulla traccia rispetto a Δt = 10⁻⁴ sono 1,447·10⁻⁴ e 4,831·10⁻⁵; la norma di traccia usa la matrice di massa lineare Q1. Il run Δt = 5·10⁻⁴ usato nel raffinamento spaziale differisce dal caso Δt = 10⁻⁴ di L2 bulk 1,720·10⁻⁴ e L2 traccia 1,929·10⁻⁴.

È **convergenza temporale verso la soluzione semidiscreta (D)**, non errore rispetto alla soluzione esatta della PDE. Dati: [`temporal_refinement.csv`](temporal_refinement.csv).

### 3.4 Geometria e accoppiamento

Nel test 8×8, T = 0,01, ogni segmento di bordo è invertito e l'ordine degli elementi di bordo è permutato. Il campo nodale finale è identico alla numerazione standard: RMS e massimo scarto sono entrambi zero nei dati esportati. La prova conferma l'invarianza rispetto all'orientamento delle linee per questa geometria rettangolare e quadratura simmetrica. Non verifica normali su bordi curvi o mesh non conformi; inoltre la normale non compare esplicitamente nell'assemblaggio parabolico no-MMS.

Non c'è ambiguità tra DOF superficiali e bulk: la variabile al bordo è la traccia nodale condivisa. Le linee di bordo sono generate sui quattro lati e validate senza doppioni. Il test controlla una geometria semplice, ma non dimostra l'accuratezza del flusso normale su una superficie curva o con coefficienti bulk/superficie diversi.

### 3.5 Energia e dissipazione

Per il ramo senza penalità il funzionale corretto è

\[
E(\phi,\psi)=\int_\Omega\left[\tfrac12|\nabla\phi|^2+\frac{(\phi^2-1)^2}{4\varepsilon^2}\right]dx+
\int_\Gamma\left[\tfrac12|\nabla_\Gamma\psi|^2+\frac{(\psi^2-1)^2}{4\varepsilon^2}\right]dS.
\]

La legge continua è dE/dt = −γ‖μ‖²Ω − M‖μΓ‖²Γ ≤ 0. Per il mapping usato è anche dE/dt = −‖φ_t‖²Ω/γ − ‖ψ_t‖²Γ/M. Non ho esteso questa proprietà a passi temporali arbitrari: backward Euler completamente implicito su un potenziale non convesso non implica qui, senza ulteriore analisi, stabilità energetica incondizionata.

Sul caso N = 20, Δt = 10⁻⁴, ε = 0,1, l'energia è stata ricalcolata a ogni passo con 3 punti di Gauss (esatti per il potenziale quartico Q1) e le componenti bulk/superficie sono entrambe non crescenti in tutti i 1000 passi. E_bulk passa da 63,71889 a 21,28574; E_surface da 81,33333 a 28,40098; E_total da 145,05222 a 49,68672.

Come controllo del bilancio ho confrontato ΔE con la dissipazione secante ∫(Δu/Δt)²/γ nel bulk e ∫Γ(Δu/Δt)²/M sul bordo, usando la matrice di massa Q1. Il residuo relativo medio è 8,54·10⁻⁴ e il massimo 1,73·10⁻³; l'approssimazione secante non è la dissipazione esatta integrata nel tempo, quindi questa è una verifica numerica di consistenza, non una prova della legge continua. Figure e dati: [`energy_time.png`](energy_time.png), [`energy_balance_secant.png`](energy_balance_secant.png), [`energy_interface.csv`](runs/energy_n20_dt0p0001_every_step/energy_interface.csv), [`energy_balance_secant.csv`](energy_balance_secant.csv). È **controllo energetico (F)**; non è una validazione indipendente del campo.

## 4. Giudizio per componente

| Componente | Esito |
|---|---|
| Mapping dell'equazione, coefficienti e forma debole bulk-superficie | Coerente per γ=M=0,01; ε=0,1 assunto; `surfacebeta=1`. Non è una verifica quantitativa indipendente del termine normale. |
| Campo bulk 2D e dinamica superficiale | Supportati da convergenza numerica a T=0,0625; non robusti alla mesh dimostrati a T=0,1. La superficie è la traccia condivisa, non un campo indipendente. |
| Accuratezza temporale | Verificata quantitativamente come primo ordine verso il riferimento numerico più fine sulla mesh fissata; ordine apparente 0,9994. |
| Energia | Monotonia osservata per questo caso e Δt in ogni step; bilancio secante con residuo relativo max 0,173%. Non si dichiara stabilità energetica per ogni Δt. |
| Geometria / orientamento | Invarianza numerica esatta nel test rettangolare 8×8 invertito e permutato; geometrie curve e accoppiamento con coefficienti diversi non verificati. |
| Xiao–Chen, Caso 4 | Solo confronto qualitativo. ε mancante e figure senza dati grezzi impediscono confronto quantitativo affidabile. |
| Soluzione analitica 2D o solver indipendente | Non disponibile in questi test; nessuna validazione A o C per il problema accoppiato 2D. |

**Conclusione:** il solver mostra una dinamica 2D coerente e una buona convergenza temporale; la prova spaziale più fine supporta convergenza ma lascia errori di riferimento non nulli, e il test a T=0,1 resta sensibile alla mesh. Il confronto con la letteratura non convalida quantitativamente il campo perché il Caso 4 non dichiara ε. Il risultato complessivo è quindi una **verifica numerica parziale**, non una validazione indipendente completa del modello parabolico 2D.

## 5. Riproducibilità e artefatti

- [`run_case.py`](run_case.py): genera l'input per N, Δt, T, ε, cadenza di output e variante d'orientamento; verifica che le linee di bordo coprano ciascun lato una volta.
- [`run_suite.py`](run_suite.py): sequenza riproducibile di tutti i run del report (può richiedere diversi minuti).
- [`analyze.py`](analyze.py): ricalcola energie, norme, metriche d'interfaccia, CSV e figure dagli output VTU.
- [`make_qualitative_comparison.py`](make_qualitative_comparison.py): genera il pannello qualitativo a T=0,1.
- Gli input generati, i log e gli stati finali VTU sono in [`runs/`](runs/). La prova energetica conserva la serie CSV a tutti i 1000 passi e il campo finale; i VTU intermedi sono rigenerabili dal runner e non sono archiviati perché voluminosi. `run_suite.py` salta i casi già completi; per rigenerare una storia completa di VTU, lanciare `run_case.py` in una nuova cartella di output.

Per rieseguire un caso singolo, ad esempio la mesh 40×40 a T=0,1:

```sh
python3 doc/stvowu_validation/parabolic_2d/run_case.py \\
  --n 40 --dt 0.0005 --end-time 0.1 --eps 0.1 --output-every 200 \\
  --out doc/stvowu_validation/parabolic_2d/runs/space_n40_dt0p0005
```

Per rigenerare metriche e figure:

```sh
python3 doc/stvowu_validation/parabolic_2d/analyze.py
python3 doc/stvowu_validation/parabolic_2d/make_qualitative_comparison.py
```

