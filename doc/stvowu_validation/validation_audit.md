# Audit di validazione del solutore Allen–Cahn OOFEM

Data dell'audit: 26 ottobre 2026. Il perimetro è il codice e gli artefatti presenti nel working tree corrente; questo working tree aveva già modifiche non commesse prima dell'audit. In questa fase non ho modificato sorgenti C++ né input/CSV archiviati. Le repliche qui descritte sono state scritte in `/tmp/oofem_ac_validation_audit/`, per non sostituire i risultati esistenti. Nessun test WED è stato lanciato.

## Esito in breve

Il codice implementa un Allen–Cahn scalare nel bulk con una dinamica di Allen–Cahn sulla frontiera espressa nella formulazione debole attraverso elementi `line1ht` che condividono i nodi con gli elementi bulk. In modalità `wave 1`, l'inerzia è nel bulk; il bordo rimane del primo ordine nel tempo. Non è il modello acustico con inerzia superficiale dei preprint CRC.

Le evidenze più forti disponibili sono circoscritte:

- **Benchmark naturale Okumura, parabolico 1D:** il run ripetuto è riproducibile e la curva energetica concorda con la digitalizzazione della tesi (RMSE 0.988, MAE 0.880, massimo 1.832 su scala circa 0–90). È un confronto con dati pubblicati, non con una soluzione tabulata esatta. La dinamica dei due estremi è rappresentata da linee verticali di altezza unitaria in una striscia Q1.
- **Onda Nizovtseva, iperbolico bulk:** il run ripetuto sulla mesh intermedia dà RMS 7.3231 × 10⁻⁶ contro la formula esatta; il CSV archiviato riporta ordini spaziali 1.995 e 1.976. È un controllo diretto molto buono dell'equazione bulk e dell'avvio della velocità, con coda troncata e Neumann naturale.
- **Newmark bulk:** il modo lineare e l'oscillatore uniforme sono stati confrontati, rispettivamente, con la soluzione esatta Q1 semidiscreta e con un'integrazione DOP853 indipendente. Le repliche a Δt=0.01 danno RMS 2.176 × 10⁻⁵ e 5.064 × 10⁻⁷.
- **Accoppiamento bulk-superficie iperbolico:** esiste solo un raffinamento temporale D, non un confronto indipendente completo. L'ordine 1.95 sul caso corrente dimostra una tendenza di convergenza ma non dimostra che il termine di bordo e il suo accoppiamento abbiano i coefficienti/segni corretti.
- **Energia parabolica:** i tre run ripetuti hanno energia monotona e residuo del bilancio discreto normalizzato che scende da 4.23 × 10⁻³ a 1.20 × 10⁻³. Il valore è una diagnostica numerica a passi finiti, non una prova di identità energetica esatta.
- **Calatroni–Colli e Xiao–Chen:** il primo è un confronto teorico con conseguenze di esistenza/unicità più raffinamenti D; non fornisce valori esatti per sovrapposizione. Il secondo è E, con ipotesi esplicita ε=0.1 non pubblicata e mesh molto più grossolana.

**Verdetto circoscritto:** la componente bulk parabolica ha un benchmark pubblicato favorevole nel caso 1D e un controllo analitico del kink di equilibrio; il bulk iperbolico ha una soluzione esatta non fabbricata favorevole. Il problema parabolico 2D e il sistema iperbolico accoppiato non sono ancora validati contro una soluzione indipendente del medesimo modello. La convergenza temporale e spaziale già osservata non colma da sola questa lacuna.

## 1. Codice esaminato e modello ricostruito

### Classi e flusso di assemblaggio

- `src/tm/Materials/allencahnmat.{h,C}` — coefficienti, flusso diffusivo, capacità, reazione Allen–Cahn e sorgente MMS opzionale. La classe restituisce il flusso convenzionale, pari a meno k per gradiente di u, la capacità e il termine sorgente; vedi `computeFlux3D`, `giveCharacteristicValue`, `computeInternalSourceVector`.
- `src/tm/EngineeringModels/nltransienttransportproblem.{h,C}` — modello transiente nonlineare. Il ramo ordinario assembla il passo implicito e la tangente della reazione; il ramo opt-in `wave 1` usa Newmark medio-accelerazione, con velocità/accelerazione memorizzate dal modello.
- `src/tm/EngineeringModels/nonstationarytransportproblem.C` — `TransportExternalForceAssembler` somma carichi e sorgente interna; gli assemblatori standard combinano capacità e tangente.
- `src/tm/Elements/transportelement.C` — matrici di diffusione e massa, vettori sorgente, tangente della sorgente e aggiornamento del valore ai punti di Gauss.
- `src/tm/Elements/quad1_ht.C` e `src/tm/Elements/line1_ht.C` — discretizzazioni Q1 bulk e lineari della frontiera. `line1_ht` ha tre punti di Gauss nel codice; i benchmark Q1 documentati impostano esplicitamente `nip 9` sui quadrilateri, adeguato all'integrazione dei termini cubici.

Il campo bulk e quello di frontiera **non sono due incognite separate**: una linea `line1ht` usa gli stessi ID di nodo del bordo dei quadrilateri. La continuità della traccia è quindi imposta cinematicamente dalla condivisione dei DOF. La classe non verifica però che le linee selezionate formino una partizione non duplicata del bordo fisico: questa è una responsabilità del modello/input.

### Forma debole effettivamente assemblata

Siano (k) la conducibilità, (c) la capacità bulk, (m) la mobilità, (c_\Gamma) la capacità superficiale, (m_\Gamma=m\,\texttt{surfacebeta}), ρ la densità bulk e ζ=`wavedamping`. Per il caso scalare con doppio pozzo, gli assemblaggi corrispondono a:

**Parabolico, bulk**
\[
 c u_t-k\Delta u=m(u-u^3)\quad\text{in }\Omega.
\]

**Parabolico, bordo dinamico**
\[
 c_\Gamma u_t-k_\Gamma\Delta_\Gamma u+k\partial_nu
     =m_\Gamma(u-u^3)\quad\text{su }\Gamma.
\]

Il flusso normale non appare come carico separato nella forma debole combinata: il termine negativo k per flusso normale del bulk si cancella con il termine positivo k per flusso normale dell’equazione superficiale. Il risultato assemblato somma massa/reazione/diffusione bulk e massa/reazione/diffusione superficiale sulle stesse incognite di traccia. Per riprodurre questa forma servono linee di bordo coerenti e condivise; una linea staccata crea una variabile indipendente e un problema diverso.

**Iperbolico OOFEM (`wave 1`)**
\[
 \rho u_{tt}+\rho\zeta u_t-k\Delta u=m(u-u^3)\quad\text{in }\Omega,
\]
\[
 c_\Gamma u_t-k_\Gamma\Delta_\Gamma u+k\partial_nu
     =m_\Gamma(u-u^3)\quad\text{su }\Gamma.
\]
L'inerzia è esclusivamente bulk; `surfacecapacity` moltiplica la velocità u_t sulle linee. `wavedamping` è un **rapporto** che moltiplica la matrice di massa bulk, non il coefficiente fisico autonomo: il coefficiente di smorzamento bulk è ρζ. L'equazione dell'onda non è attivata automaticamente dal materiale: occorrono sia `hyperbolic 1` nel materiale per la capacità bulk ρ, sia `wave 1` nel modello per il residuo Newmark.

### Segni e tangente nonlineare

Il flusso pari a meno k per gradiente u è assemblato da `computeInternalForcesVector` con segno che produce il termine positivo di diffusione nel residuo. La sorgente m(u-u³) è assemblata sul lato destro e quindi sottratta nel residuo. La derivata del residuo di reazione è m(3u²−1), restituita per `IntSource`; il Newton completo è attivato in modalità onda. Per la modalità parabolica la tangente completa viene formata dall'assemblaggio di conducibilità, derivata della sorgente e capacità. L'aggiornamento dello stato ai punti di Gauss precede la lettura di (u) usata nella sorgente/tangente.

Con `alpha=1` e `bdforder=1`, l'equazione ordinaria è discretizzata con backward Euler. Il ramo `wave 1` costruisce il predittore/correttore Newmark; `beta=1/4`, `gamma=1/2` è la scelta della media accelerazione usata dai test. Il ramo Newmark aggiorna la tangente completa a ogni iterazione. I test qui citati non verificano le altre coppie β,γ.

### Parametri configurabili e limiti

- `k` controlla diffusione bulk e, se assegnato un materiale a una linea, anche la diffusione tangenziale. Un materiale distinto sulle linee permette una diffusione tangenziale kΓ diversa da k bulk; non esiste però una keyword specifica `surfacek`.
- `capacity` è usata nel problema parabolico sul tipo di elemento cui è assegnato il materiale. Nel ramo iperbolico, la capacità bulk diventa `rho` e quella delle linee `surfacecapacity`.
- `mobility` scala la reazione bulk e di linea; `surfacebeta` dà un fattore moltiplicativo aggiuntivo alla reazione superficiale, non alla diffusione o capacità.
- I benchmark naturali usano condizioni naturali di Neumann dove non sono assegnate linee dinamiche o vincoli: il codice non applica un flusso normale imposto implicito al posto della condizione naturale.
- Le normali non entrano nel residuo debole ordinario. Una scelta d'orientamento errata delle linee non dovrebbe cambiare il problema di diffusione 1D; l'orientamento viene usato soltanto dal codice di sorgente MMS per ricostruire una normale.

## 2. Audit dei test, riproducibilità e classi A–F

| Prova | Classe | Modello / osservabile | Riscontro audit | Limite e conclusione ammissibile |
|---|---|---|---|---|
| Okumura, Example 2 | **B** | Parabolico 1D, dati e parametri pubblicati; energia digitale Fig. 4.4 | Ripetuti 700 passi in cartella temporanea. Energia da 77.42047 a 2.2194 × 10⁻⁶, massimo incremento negativo. Stessa analisi della curva digitalizzata: 391 campioni, RMSE 0.98810, MAE 0.88037, max 1.83230. | Benchmark di letteratura con curva raster, non soluzione analitica; valido per quel caso 1D estruso, non per variazioni tangenziali 2D. |
| Energia naturale 2D | **F** | Funzionale Allen–Cahn: energia bulk più energia tangenziale e potenziale sulla frontiera | Tre run ripetuti; tutti monotoni. E0=1.28445597; ET=1.20911974, 1.20901760 e 1.20896702. Residui normalizzati 4.2314 × 10⁻³, 2.2885 × 10⁻³ e 1.2005 × 10⁻³. | Controllo energetico diagnosticato, non confronto puntuale. Il residuo discreto include errori di temporizzazione/quadratura e tende a zero con Δt. Non implica conservazione della massa, che infatti cresce da 0.6 a circa 0.662. |
| Parabolico, mesh 10/20/40 | **D** | Dato naturale liscio, differenze nodali con mesh 40 | Run ripetuto in `/tmp`; RMS 1.27871 × 10⁻⁴ e 2.59027 × 10⁻⁵, identici ai CSV archiviati. | La mesh più fine è un riferimento numerico, non esatto; confronto di griglia solo con lo stesso schema/famiglia. |
| Calatroni–Colli | teoria + **D** | Compatibilità del potenziale polinomiale e raffinamenti spazio/tempo | Equazione polinomiale compatibile con il quadro teorico dichiarato; i CSV mostrano tendenza temporale/spaziale. | Teoremi di esistenza/unicità non forniscono valori della soluzione. Non è una verifica A e i raffinamenti da soli non certificano l'equazione. Regolarità del dominio e condizioni del teorema vanno ricontrollate sulla specifica mesh/dato. |
| Xiao–Chen, Caso 4 | **E** | Goccia quadrata, contorno u=0.5, forma a T=0.1 | La morfologia è coerente a vista con la figura. | Il riferimento omette ε per quel caso; è stata assunta ε=0.1. Mesh 20×20 contro h circa 1/300; i CSV mostrano sensitività spaziale. Nessuna validazione quantitativa contro l'articolo. |
| Modo coseno lineare bulk | **A** (semidiscreto) | Newmark contro soluzione esatta del sistema Q1 semidiscreto consistente | Run Δt=.01 ripetuto: RMS 2.17590 × 10⁻⁵, coincide col CSV. | Verifica temporale Newmark dopo discretizzazione spaziale, non accuratezza spaziale rispetto al PDE continuo né bordo dinamico. |
| Oscillatore Allen–Cahn uniforme | **C** | Bulk u_tt=u-u^3, confronto DOP853 indipendente | Run Δt=.01 ripetuto: RMS 5.06384 × 10⁻⁷, riferimento u(0.2)=0.2641955936. | Controlla reazione e integrazione temporale su soluzione spazialmente uniforme; non diffusione né frontiera. |
| Kink stazionario tanh | **A** | u=tanh(x/sqrt(2)), OOFEM a T=0.2 | Run (N=80) ripetuto RMS 7.70736 × 10⁻⁶, uguale al CSV. CSV (N=40,80,160) dà ordini 2.006 e 1.887. | Soluzione analitica bulk, ma il dominio [-8,8] usa Neumann naturale mentre la derivata esatta ai bordi non è esattamente zero; non c'è bordo dinamico. |
| Onda viaggiante Nizovtseva et al. | **A** | Eq. esatta (29), τ=.2, T=0.2, velocità iniziale analitica | Run (N_x=160) ripetuto: RMS 7.32307 × 10⁻⁶, max 2.88898 × 10⁻⁵. CSV raffinamenti h=.3,.15,.075: ordini 1.995,1.976. | Evidenza diretta forte per il bulk; Neumann ai lati troncati approssima le code esponenzialmente piccole; non verifica l'equazione di superficie. Il controllo temporale a h=.15 è sul plateau spaziale e non stima robustamente l'ordine temporale. |
| Sistema iperbolico accoppiato | **D** | Soluzioni finali Δt=.002,.001,.0005, T=0.02 | È archiviato ordine apparente 1.950; nel campione ripetuto il run Δt=.0005 termina e produce campo nel range [0.10354, 0.29819]. | La tripla è convergenza interna contro la soluzione più fine, senza soluzione indipendente o analitica; non dimostra segni/scala della capacità superficiale. |
| Profilo Hochbruck–Leibold adattato | **Non ammesso (MMS)** | Profilo sin(2 pi t) x y con `manufactured 1` e forzante ricalcolata per OOFEM | Il file/script usa sorgenti artificiali; non è prova valida secondo il vincolo di questo audit. | Da escludere da conteggi e conclusioni. Le metriche precedenti RMS 1.4 × 10⁻² non sono considerate evidenza di validazione. |

### Altri input Allen–Cahn nel repository (fuori da stvowu_validation)

Ho esteso l'inventario anche a `tests/tm/` e `doc/allen_cahn_paper_cases/`, così l'audit non si limita alla cartella indicata nel contesto.

| Input / studio | Categoria e modello | Audit / limite |
|---|---|---|
| `tests/tm/allencahn.in` | Smoke test parabolico bulk, Q1 20×20, 10 passi | Ripetuto in directory temporanea; OOFEM termina regolarmente. Non c'è soluzione di riferimento o metrica d'errore: è solo regressione/esecuzione. |
| `tests/tm/allencahn_dynamic_boundary.in` | Smoke test parabolico bulk + linee dinamiche, Q1 20×20, 40 passi | Ripetuto in directory temporanea; termina. Non confronta valori indipendenti e non prova il segno dell'accoppiamento. |
| `tests/tm/allencahn_mixed_dirichlet_neumann.in`, `allencahn_triangular_dynamic.in`, `allencahn_paper_mms_sine.in` | MMS / sorgente costruita | Esclusi dall'uso e non lanciati in questo audit, per il vincolo no-MMS. Il precedente resoconto di `ctest -R allencahn` (6/6) includeva test MMS e non vale come evidenza ammessa in questa valutazione. |
| `doc/allen_cahn_paper_cases/lee2024_equilibrium_t5/` | **A**, equilibrio analitico bulk: u*(x)=tanh(x/(sqrt(2) epsilon)); rilassamento da 0.5 tanh senza sorgente, epsilon=0.009606, h=0.01, dt=1e-5, T=5 | Ho ricalcolato metriche dal VTU finale archiviato (500000 passi): errore L2 continuo ricostruito 6.0675 × 10⁻³, massimo nodale 1.30468 × 10⁻². È un confronto A del bulk con l’equilibrio analitico, non del bordo dinamico; il paper usa uno schema compatto FD di ordine 4, OOFEM P1 + backward Euler, quindi non è confronto a pari discretizzazione. Non ho ripetuto 500000 passi per il costo; ho verificato gli output archiviati. |
| `doc/allen_cahn_paper_cases/akrivis_ellipse_be_*.in` e generazione paper-parameter | **B/D**, frontiera ellittica naturale e confronto temporale, senza sorgente MMS | La prova pubblicata confronta schemi BDF5/Radau; i casi `be` sono backward Euler e la mesh delle sequenze originarie è grossolana. Il generatore con h=ε/2 e BDF5 non riproduce i primi quattro passi Radau; il run ε=0.01 è incompleto per costo. È supporto per convergenza temporale, non una replica completa o una validazione del bordo dinamico. |
| `doc/allen_cahn_paper_cases/reproduction_inputs/` e `refined_inputs_n40/` | **B/E**, energia di Wang Fig. 1 e Akrivis–Li Fig. 3, senza sorgente artificiale | Il report archiviato dice 4 casi finiti/dissipativi e 2 casi ε=0.01 non finiti/non monotoni. Schemi e mesh differiscono dagli articoli, i riferimenti sono figure e non curve numeriche tabulate. I VTK in `/tmp/ac_refined` non sono presenti nel filesystem corrente: ho controllato il report/input ma non ho ricalcolato queste sei curve. Non contano come validazione diretta. |
| `doc/allen_cahn_paper_cases/chen2025_mms/` | MMS per Wang et al., soluzione sinusoidale e sorgente apposita | Escluso. I tassi riportati possono servire per la suite di verifica solo se il vincolo MMS viene revocato; qui non li considero nella conclusione. |

La formula di Lee, u*=tanh(x/(sqrt(2) epsilon)), è la soluzione stazionaria del bulk u_xx = epsilon^(-2)(u^3-u). Il dato iniziale è a metà ampiezza, mentre il confronto finale è col kink completo. Il profilo è esatto sul dominio infinito; sui limiti finiti [-1,1] il Neumann naturale è un’approssimazione, ma la derivata alle estremità è esponenzialmente piccola. L’errore residuo è coerente con h≈epsilon e P1, che non risolvono finemente l’interfaccia.

### Integrità e ripetibilità

I casi ripetuti qui sono stati lanciati dal binario `build/debug_magneto/oofem` in directory nuove sotto `/tmp/oofem_ac_validation_audit/`; le metriche sono state ricalcolate dai VTU o dalle stesse routine di postprocessing copiate in `/tmp`. I CSV originali e le cartelle di run nel repository non sono stati sovrascritti. La replica indipendente di Okumura ha prodotto 700 VTU e ha rigenerato l'energia; la replica energia 2D e il raffinamento spaziale hanno dato gli stessi numeri archiviati. La mesh intermedia Nizovtseva e i campioni Newmark/kink/accoppiato sono stati ripetuti e corrispondono ai CSV.

**Nota di processo sull'MMS:** durante il campionamento ho lanciato una volta in directory temporanea l'input del test adattato Hochbruck–Leibold prima di rilevare `manufactured 1`. Questo run è escluso integralmente dai riscontri e non è usato come prova o risultato. L'output temporaneo è in `/tmp/oofem_ac_validation_audit/hochbruck_disk_h016`; i risultati archiviate originali non sono stati toccati. Il test non verrà ripetuto nell'ambito di questo audit.

## 3. Possibili bug e rischi

### Dimostrati dal codice (condizionali ai casi d'uso)

1. **Velocità/accelerazione dei DOF prescritti non entra nel residuo Newmark con il valore imposto.** In `giveUnknownComponent`, DOF senza equation number possono consultare il BC per velocity/acceleration, ma il ramo di assemblaggio Newmark costruisce direttamente `currentVelocity/currentAcceleration` da `initialVelocity` e zero accelerazione per `eq<=0`; lo stesso fallback compare nel calcolo iniziale. Se si usano condizioni di Dirichlet dipendenti dal tempo in `wave 1`, i contributi di massa/capacità sui DOF vincolati sono quindi incoerenti con la derivata del valore prescritto. I test corrente sono Neumann naturali e non attivano il caso. È una lacuna funzionale dimostrabile per quel caso d'uso; va isolata in un test prima di una correzione.
2. **Sorgente MMS iperbolica parametrizzata con `capacity*u_t`, mentre il ramo onda assembla `rho*wavedamping*u_t`.** `computeManufacturedSource` restituisce ρu_tt + c u_t nel bulk; il residuo `wave 1` è ρu_tt+ρζu_t. Sono equivalenti solo se (c=\rho\zeta). È dimostrato per il percorso MMS iperbolico con coefficienti non corrispondenti; non influenza i benchmark non-MMS ammessi. Da tenere fuori dalla validazione qui, ma da documentare/controllare prima di riutilizzare quel flag.

### Sospetti / non dimostrati

- **Topologia delle linee dinamiche:** non c'è un controllo applicativo che ogni tratto fisico di Γ sia presente una sola volta come `line1ht`, con nodi condivisi ai bulk. Input incompleto o duplicato cambia massa e diffusione superficiali senza errore. È una precondizione di mesh, non un bug dimostrato nei casi esaminati.
- **Quadratura insufficiente se si omette `nip`:** `quad1_ht` ha una regola predefinita a 4 punti. I vettori (N_i(u-u^3)) per Q1 e la tangente (N_iN_j(3u^2-1)) non sono integrati esattamente da 2x2; i benchmark attuali impostano `nip 9`, ma input futuri potrebbero non farlo. Va quantificato come errore di quadratura, non descritto come errore algebraico certo.
- **Sensibilità all'orientamento ai vertici:** il debole lineare è invariante a inversione delle linee, ma il generatore MMS disk-specifico ricostruisce normali edgewise e il profilo (xy) ha salti di derivata tangenziale ai vertici del disco poligonale. L'elevato errore spaziale del test H–L non separa gli effetti di geometria, sorgente singolare ai vertici e implementazione. Inoltre quel test è MMS e resta escluso.
- **Matrice di massa e vincoli non standard:** i test Newmark usano massa consistente e nessun Dirichlet. Non sono ancora verificati l'eliminazione dei DOF prescritti, condizioni non omogenee, massa concentrata (`lumpedcapa`) o mesh con equation numbering non banale.
- **Regime di reazione nonlineare e Newton:** non esiste ancora un test diretto della tangente tramite differenze finite del residuo; la convergenza delle simulazioni non certifica da sola la derivata `IntSource`.
- **Dissipazione iperbolica:** non è archiviato un bilancio comprendente energia cinetica bulk, energia Allen–Cahn e dissipazioni bulk/superficie. Non è quindi verificata la scala/segno dello smorzamento né la dissipazione del bordo nel caso onda accoppiato.

## 4. Limiti delle conclusioni attuali

- Le classi **A** bulk (kink e onda viaggiante) verificano fortemente la reazione/diffusione bulk e l'integrazione temporale Newmark nei casi provati, ma non la dinamica superficiale.
- Il test **B** di Okumura verifica un particolare problema parabolico 1D con capacità agli estremi e accordo con grafico digitale; non basta a convalidare la diffusione tangenziale 2D.
- Le classi **C**, **D**, **E** e **F** sono informazioni complementari, non equivalenti a un confronto esatto/publicato punto-a-punto.
- La prova Hochbruck–Leibold con sorgente costruita non è ammessa dal vincolo no-MMS e viene esclusa. I loro benchmark inoltre descrivono un operatore al bordo diverso.
- La teoria di Calatroni–Colli sostiene proprietà del problema continuo sotto ipotesi, ma il raffinamento del codice non è una prova del teorema di convergenza per questo schema.

Il piano quantitativo per colmare le lacune è in [`validation_matrix.md`](validation_matrix.md).
