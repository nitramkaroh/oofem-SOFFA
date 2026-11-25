# Validazione del ramo Allen–Cahn iperbolico nel bulk

**Verdetto: PASS per il benchmark bulk di onda viaggiante di Nizovtseva et al., con il perimetro e le riserve descritti qui.** Il profilo esatto, la propagazione, la convergenza spaziale Q1, l’ordine temporale Newmark e il bilancio energetico sono stati controllati quantitativamente. Questo esito non si estende al sistema bulk-superficie accoppiato né prova ogni opzione del solutore.

Nessun MMS, carico o termine sorgente artificiale è stato usato. Non è stato modificato codice C++ durante questa verifica e non è stato eseguito lavoro WED.

## 1. Riferimento e soluzione esatta

Fonte originale: I. G. Nizovtseva, P. K. Galenko e D. V. Alexandrov, “The hyperbolic Allen–Cahn equation: exact solutions”, *Journal of Physics A: Mathematical and Theoretical* 49 (2016), 435201, [DOI 10.1088/1751-8113/49/43/435201](https://doi.org/10.1088/1751-8113/49/43/435201). Ho controllato il PDF fornito, in particolare le equazioni (1)–(4) e la soluzione (29), pagina stampata 6. Il paper dà una soluzione analitica esatta e illustra i profili, ma non specifica una discretizzazione numerica o un benchmark a dominio finito con dati tabulati: questo è un confronto OOFEM contro soluzione analitica, non una replica di un calcolo numerico pubblicato. La formula è esatta sulla retta infinita; la condizione Neumann su un dominio finito è una scelta del test e il suo effetto è misurato sotto.

L’equazione dimensionale del paper è

\[
\tau_R f_{tt}+f_t=D\Delta f+M(f-f^3),
\]

con potenziale bistabile \(F(f)=\tfrac14(f^2-1)^2\), per cui \(-M F'(f)=M(f-f^3)\).

con tempo e spazio adimensionali \(t'=Mt\), \(x'=x\sqrt{M/D}\), e \(\tau=\tau_R M\). Per la soluzione (29), scegliendo i segni \(A=B=Q=+1\) e \(c_0=0\),

\[
f(x,t)=\frac{1}{1+\exp\!\left[\frac{\sqrt{2+9\tau}}{2}\,x'-\frac32t'\right]}.
\]

La velocità del profilo è positiva; il livello \(f=1/2\) si muove verso x crescente. La velocità iniziale da passare al codice è la derivata temporale della formula,

\[
f_t(x,0)=\frac32 M f(x,0)\,[1-f(x,0)],
\qquad
f_{tt}(x,0)=\frac94 M^2f(x,0)[1-f(x,0)][1-2f(x,0)].
\]

Il fattore spaziale è \(\sqrt{2+9\tau}/2\), non \(\sqrt{(2+9\tau)/2}\). Lo script usa la forma del paper.

### Parametri e traduzione in OOFEM

Il benchmark usa \(D=M=1\), \(\tau_R=0.2\), dunque \(\tau=0.2\), \(x'=x\), \(t'=t\), \(\lambda=\sqrt{3.8}/2\), e velocità analitica \(c=3/(2\lambda)=1.5389675\). In forma dimensionale la PDE è \(0.2 f_{tt}+f_t=\Delta f+(f-f^3)\).

La lettura del codice dà l’equazione bulk assemblata

\[
\rho u_{tt}+\rho\,\delta u_t-k\Delta u-m(u-u^3)=0,
\]

dove `rho` è la massa inerziale, `waveDamping` è \(\delta\), `k` è la diffusività e `mobility` è \(m\). Pertanto l’input usa `rho=0.2`, `wavedamping=5`, `k=1`, `mobility=1`: il prodotto \(\rho\,\texttt{wavedamping}=1\) riproduce il coefficiente di \(f_t\) del paper. `c=1` è la capacità parabolica e non è il coefficiente inerziale quando `hyperbolic=1`.

In [allencahnmat.C](../../../src/tm/Materials/allencahnmat.C), il ramo iperbolico restituisce capacità \(\rho\), flusso \(-k\nabla u\) e sorgente \(m(u-u^3)\). [TransportElement](../../../src/tm/Elements/transportelement.C) assembla la diffusione e la reazione. [NLTransientTransportProblem](../../../src/tm/EngineeringModels/nltransienttransportproblem.C) applica Newmark con \(\beta=1/4\), \(\gamma=1/2\), alpha 1 e Newton completo; il tangente è aggiornato alle iterazioni non lineari. Il termine di smorzamento si forma come \(\rho\delta M u_t\). Il ramo calcola l’accelerazione iniziale coerente dal residuo semidiscreto; non parte da accelerazione zero.

## 2. Geometria, condizioni e controlli iniziali

Il riferimento OOFEM è una striscia rettangolare \([-12,12]\times[0,1]\), con 80, 160 o 320 elementi Q1 lungo x e un solo strato lungo y. La soluzione esatta è indipendente da y; le due file di nodi hanno quindi gli stessi valori. Non ci sono elementi di bordo né carichi: la forma debole applica Neumann naturale su tutti i lati. Sul lato y la condizione esatta è omogenea; sui due estremi x la derivata esatta non è esattamente zero, ma a \(x=\pm12\) è molto piccola. La sua influenza è misurata con una prova di dominio troncato sotto.

Tutti i run terminano a \(T=0.2\). Il profilo iniziale è campionato direttamente dalla formula esatta in ogni nodo. `initialvelocityvalues` contiene un valore per equazione/DOF nell’ordine nodale dell’input. Il confronto indipendente tra questi valori e \(f_t(x,0)\) dà mismatch massimo \(3.61\cdot10^{-16}\) in tutti i run: la velocità iniziale nodale è corretta entro l’arrotondamento macchina.

La formula dà accelerazione iniziale analitica riportata sopra. Il solver invece calcola l’accelerazione iniziale coerente con il residuo Q1 semidiscreto, includendo reazione, diffusione e smorzamento. L’accelerazione interna non è esportata nel VTU, perciò non è stata confrontata nodalmente con \(f_{tt}(x,0)\). Il test di ordine temporale e il modo lineare/oscillatore indipendenti verificano l’avvio Newmark senza sostituire tale misura diretta.

## 3. Convergenza spaziale

Le tre mesh usano \(\Delta t=0.001\), lo stesso dominio, lo stesso tempo finale e le coordinate fisiche della soluzione analitica. L’errore è calcolato sui nodi x unici (le due righe y duplicate sono collassate), uniformemente distribuiti: RMS nodale e massimo nodale. La metrica interna include i nodi \(|x|\le6\), lontani dagli estremi. Non è un errore rispetto a un’altra simulazione.

| N in x | h | RMS full | max full | RMS \(|x|\le6\) | max \(|x|\le6\) | ordine RMS full / interno |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 0,300 | 2,9183·10⁻⁵ | 1,1579·10⁻⁴ | 4,1011·10⁻⁵ | 1,1579·10⁻⁴ | — |
| 160 | 0,150 | 7,3231·10⁻⁶ | 2,8890·10⁻⁵ | 1,0309·10⁻⁵ | 2,8890·10⁻⁵ | 1,995 / 1,992 |
| 320 | 0,075 | 1,8621·10⁻⁶ | 7,1898·10⁻⁶ | 2,5811·10⁻⁶ | 7,1898·10⁻⁶ | 1,976 / 1,998 |

Gli ordini sono compatibili con Q1 di secondo ordine per questo profilo regolare. RMS interno maggiore dell’RMS sull’intero dominio è atteso: l’errore si concentra nell’interfaccia, mentre le code quasi costanti abbassano la media full-domain. I risultati sono in [`nizovtseva_bulk_v2_spatial.csv`](nizovtseva_bulk_v2_spatial.csv).

## 4. Convergenza temporale e propagazione

La mesh è fissata a N=320, \(h=0.075\); si confrontano i campi a T=0.2 sugli stessi nodi. Newmark average acceleration è di ordine due. La tabella riporta anche la differenza RMS fra un passo e quello dimezzato; gli ordini da queste differenze sono le stime temporali rilevanti, poiché il contributo spaziale è comune alle coppie.

| Δt | RMS vs esatta | max vs esatta | RMS vs passo dimezzato | ordine dalle differenze |
|---:|---:|---:|---:|---:|
| 0,040 | 1,9364·10⁻⁶ | 4,7250·10⁻⁶ | 2,1093·10⁻⁶ | 2,009 |
| 0,020 | 1,4383·10⁻⁶ | 4,9345·10⁻⁶ | 5,2419·10⁻⁷ | 1,917 |
| 0,010 | 1,7402·10⁻⁶ | 6,5898·10⁻⁶ | 1,3885·10⁻⁷ | — |
| 0,005 | 1,8346·10⁻⁶ | 7,0898·10⁻⁶ | — | — |

L’errore rispetto alla formula non diminuisce monotonicamente: a passi piccoli il residuo spaziale domina e si somma/cancella con l’errore temporale. Per questo non uso tale colonna per stimare l’ordine. Le differenze tra campi a mesh fissa sono compatibili con secondo ordine nelle prime due coppie; la coppia 0,01/0,005 è vicina al pavimento dovuto alla precisione numerica dell’output VTU e non viene usata per un’affermazione d’ordine. Non emerge instabilità: tutti i run raggiungono T e mantengono il profilo limitato tra circa 1,5·10⁻⁵ e 0,999992.

Il livello \(u=1/2\) definisce l’interfaccia in modo identico nella formula e nei campi FE, con interpolazione lineare tra nodi. Per N=320, \(x_{0,\mathrm{exact}}=0.30779351\), \(x_{0,\mathrm{OOFEM}}=0.30782382\); la velocità misurata è 1,539119 contro 1,538968 esatta, errore relativo circa \(9.8\cdot10^{-5}\). Sulle mesh N=80 e 160 la velocità misurata è rispettivamente 1,541460 e 1,539582, mostrando convergenza nella direzione e nella velocità di propagazione.

Metriche temporali: [`nizovtseva_bulk_v2_temporal.csv`](nizovtseva_bulk_v2_temporal.csv).

## 5. Sensibilità al dominio e ai dati iniziali

La prova di troncamento tiene \(h=0.1\), \(\Delta t=0.001\) e confronta \([-8,8]\) con \([-12,12]\) sugli stessi nodi fisici in \([-8,8]\). Il max errore sul dominio corto è \(1.716\cdot10^{-4}\), contro \(1.279\cdot10^{-5}\) per la simulazione lunga; l’RMS sul comune dominio è \(2.260\cdot10^{-5}\) contro \(3.981\cdot10^{-6}\). Limitando il confronto a \(|x|\le6\), i due RMS coincidono a \(4.592\cdot10^{-6}\) e i massimi a \(1.279\cdot10^{-5}\). L’effetto della condizione Neumann troncata è quindi localizzato vicino a \(x=\pm8\); il dominio \([-12,12]\) riduce il mismatch alle estremità abbastanza da non contaminare l’interfaccia.

Come controllo della mappatura della velocità, un run N=160 con la velocità iniziale aumentata del 1% dà a T=0.2 RMS \(1.653\cdot10^{-4}\) rispetto all’onda non perturbata e posizione dell’interfaccia 0,309816, invece di 0,307916 nel run con velocità esatta. Il cambiamento è una risposta fisica ai dati iniziali e mostra che l’array di velocità influenza la propagazione; non è un criterio di stabilità o un limite d’accettazione.

Risultati: [`nizovtseva_bulk_v2_domain.csv`](nizovtseva_bulk_v2_domain.csv), [`nizovtseva_bulk_v2_initial_sensitivity.csv`](nizovtseva_bulk_v2_initial_sensitivity.csv).

## 6. Energia del sistema iperbolico smorzato

Per l’equazione \(\rho u_{tt}+\eta u_t-k\Delta u-m(u-u^3)=0\), con Neumann omogenea, l’energia continua corretta è

\[
E(t)=\int_\Omega\left[\frac{\rho}{2}u_t^2+\frac{k}{2}|\nabla u|^2+\frac{m}{4}(u^2-1)^2\right]dx,
\qquad
\frac{dE}{dt}=-\eta\int_\Omega u_t^2\,dx\le0.
\]

Nel test \(\rho=0.2\), \(\eta=1\), \(k=m=1\). L’energia è integrata sui Q1 con quadratura 3×3, esatta per il potenziale quartico; la velocità è ricostruita dalle storie di spostamento VTU con differenze finite centrate nel tempo. Il PVD comincia dal primo passo \(t=0.001\), quindi il frame iniziale è ricostruito dai dati analitici iniziali dell’input e usa la velocità analitica prescritta. L’energia osservata decresce in tutti i passi: 3,162421 a 3,085490. La dissipazione integrata numericamente è 0,07693409 e il residuo finale \(E(T)-E(0)+\int_0^T\eta\|u_t\|^2dt=3.47\cdot10^{-6}\), pari a \(1.10\cdot10^{-6}\) dell’energia iniziale.

È un controllo energetico post-processato: la derivata temporale è approssimata dai campi stampati e non coincide con la velocità interna Newmark. Non lo interpreto come identità discreta esatta né prova autonoma di stabilità incondizionata. Dati per tutti i tempi: [`nizovtseva_bulk_v2_energy.csv`](nizovtseva_bulk_v2_energy.csv).

## 7. Evidenza indipendente già presente nel repository

La prova [`run_bulk_references.py`](run_bulk_references.py) e [`bulk_reference_metrics.csv`](bulk_reference_metrics.csv) è stata riusata come controllo complementare, non duplicata. Il modo coseno lineare con `mobility=0` è confrontato con la soluzione esatta del sistema FEM Q1 a massa consistente: al dimezzamento di Δt l’RMS va da 3,4715·10⁻⁴ a 8,6981·10⁻⁵ e 2,1759·10⁻⁵, ordine 1,997. L’oscillatore spazialmente uniforme con \(u(0)=0.2\), \(u_t(0)=0.3\), \(\rho=m=1\) si riduce a \(u_{tt}=u-u^3\), confrontato con integrazione DOP853 indipendente: errori 7,3064·10⁻⁶, 1,8064·10⁻⁶ e 5,0638·10⁻⁷ per Δt 0,04, 0,02 e 0,01; ordine 2,016 sul primo dimezzamento. Questo isola la reazione non lineare e sostiene il comportamento Newmark senza diffusione spaziale.

## 8. Artefatti e riproducibilità

- [`validate_nizovtseva_bulk.py`](validate_nizovtseva_bulk.py) genera input, esegue OOFEM e salva log, PVD e tutti i VTU temporali in cartelle nuove `nizovtseva_bulk_v2_*`. Le prove precedenti non sono state sovrascritte.
- Un input spaziale rappresentativo è [`nizovtseva_bulk_v2_space_n320/case.in`](nizovtseva_bulk_v2_space_n320/case.in); ogni cartella `nizovtseva_bulk_v2_*` contiene il proprio input e `run.log` oltre agli output.
- [`analyze_nizovtseva_bulk.py`](analyze_nizovtseva_bulk.py) ricostruisce CSV e grafici dagli output grezzi.
- [Confronto profilo, raffinamenti ed energia](nizovtseva_bulk_validation.png); [dominio troncato e dato iniziale perturbato](nizovtseva_bulk_boundary_initial.png).
- Per rieseguire l’intera suite e poi rigenerare analisi/figure:

```sh
python3 doc/stvowu_validation/hyperbolic/validate_nizovtseva_bulk.py --suite all
python3 doc/stvowu_validation/hyperbolic/analyze_nizovtseva_bulk.py
```

Le metriche spaziali sono RMS nodali su griglie uniformi, non norme integrate FE; il massimo è il massimo nodale. La verifica di energia è una quadratura FE del campo ricostruito con velocità post-processata. Gli output grezzi completi sono conservati per ciascun run, perciò entrambe le analisi sono ripetibili.

## Conclusione

**PASS — ramo iperbolico bulk, benchmark Nizovtseva.** Formula e normalizzazione sono allineate; la velocità iniziale è derivata dalla formula e immessa correttamente; il campo converge spazialmente vicino al secondo ordine, Newmark mostra ordine vicino a due dove risolto, l’interfaccia propaga nel verso e alla velocità attesi, il troncamento è quantificato e il bilancio di energia è coerente entro l’errore di post-processing. L’oscillatore indipendente conferma la reazione cubica. Restano fuori da questo PASS il bordo dinamico accoppiato, geometrie multidimensionali generiche, un confronto diretto dell’accelerazione interna iniziale e una dimostrazione di stabilità energetica discreta per passi arbitrari.
