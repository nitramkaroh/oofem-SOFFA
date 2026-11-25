# Bilanci energetici dei modelli Allen–Cahn in OOFEM

## Campo e coefficienti letti dal codice

Sia \(\Omega\subset\mathbb R^d\) un dominio regolare e \(\Gamma=\partial\Omega\). Il campo di linea usato dagli elementi `line1ht` è la traccia del campo bulk, \(\psi=u|_\Gamma\): i gradi di libertà sono condivisi, non esiste un secondo campo superficiale indipendente.

Dal materiale `AllenCahnMaterial` (`src/tm/Materials/allencahnmat.C`) risultano:

- flusso \(-k\nabla u\), quindi diffusione isotropa con coefficiente `k`;
- sorgente bulk \(m(u-u^3)\);
- sorgente sulle linee \(\beta m(\psi-\psi^3)\), con `surfacebeta` \(\beta\);
- nel ramo parabolico capacità `c` sia nel bulk sia sulle linee;
- nel ramo iperbolico capacità `rho` nel bulk e `surfacecapacity` \(c_\Gamma\) sulle linee.

La diffusività tangenziale sulle linee è ancora `k`: il codice non offre un coefficiente superficiale indipendente. Se il modello temporale è `nltransienttransportproblem` con `wave 1`, Newmark tratta il bulk come secondo ordine e le linee come primo ordine. `wavedamping` \(\zeta\), predefinito a zero, moltiplica la massa bulk e genera smorzamento \(\rho\zeta u_t\). Con `wave 0` il ramo parabolico usa lo schema a un passo controllato da `alpha` (o BDF se configurato).

Le equazioni seguenti descrivono la parte non forzata (`manufactured 0`), senza carichi esterni e con le linee che ricoprono una volta ciascuna la frontiera dinamica prevista.

## Equazioni effettive e forma debole

Poniamo \(f(s)=s-s^3\) e \(F(s)=\tfrac14(s^2-1)^2\), cosicché \(F'(s)=s^3-s=-f(s)\).

### Ramo parabolico

La forma forte equivalente assemblata è

\[
\begin{aligned}
c u_t-k\Delta u &=m f(u) &&\text{in }\Omega,\\
c\psi_t-k\Delta_\Gamma\psi+k\partial_nu&=\beta m f(\psi)&&\text{su }\Gamma,\\
\psi&=u|_\Gamma.
\end{aligned}
\]

Qui \(\partial_nu=\nabla u\cdot n\), con \(n\) normale uscente. La condizione superficiale è dinamica: non si impone separatamente una condizione di Neumann omogenea. Eventuali condizioni di Dirichlet imposte su parti del bordo cambiano il bilancio e vanno trattate come vincoli/carichi; la derivazione sotto assume l'intera frontiera dinamica e nessun lavoro esterno.

Moltiplicando la prima equazione per una funzione test \(v\), integrando per parti e sommando l'equazione superficiale si ottiene, per \(v\in H^1(\Omega)\) e la sua traccia:

\[
\int_\Omega c u_t v+k\nabla u\cdot\nabla v-mf(u)v\,dx
+\int_\Gamma c\psi_t v+k\nabla_\Gamma\psi\cdot\nabla_\Gamma v-\beta m f(\psi)v\,dS=0.
\]

Infatti l'integrazione per parti del termine bulk produce \(+k\partial_nu\,v\) sul bordo; questo è esattamente il termine di accoppiamento nella seconda equazione. Nella forma assemblata non compare un termine normale esplicito sugli elementi di linea: esso è implicito nel termine al bordo della forma debole bulk. Questa formulazione equivale al coefficiente normale \(k\), non a un coefficiente autonomo.

### Ramo iperbolico bulk / parabolico superficiale

Con \(\zeta=\texttt{wavedamping}\), \(\rho=\texttt{rho}\) e \(c_\Gamma=\texttt{surfacecapacity}\):

\[
\begin{aligned}
\rho u_{tt}+\rho\zeta u_t-k\Delta u&=m f(u)&&\text{in }\Omega,\\
c_\Gamma\psi_t-k\Delta_\Gamma\psi+k\partial_nu&=\beta m f(\psi)&&\text{su }\Gamma,\\
\psi&=u|_\Gamma.
\end{aligned}
\]

Servono dati iniziali \(u(0)=u_0\), \(u_t(0)=v_0\); la superficie iniziale è la traccia \(u_0|_\Gamma\). Non c'è inerzia superficiale: `BulkWaveMassAssembler` esclude le linee dalla matrice di massa inerziale e `NewmarkWaveLhsAssembler` assegna loro il termine di velocità di primo ordine.

La forma debole è

\[
\begin{aligned}
&\int_\Omega \rho(u_{tt}+\zeta u_t)v+k\nabla u\cdot\nabla v-mf(u)v\,dx\\
&\quad+\int_\Gamma c_\Gamma\psi_t v+k\nabla_\Gamma\psi\cdot\nabla_\Gamma v-\beta m f(\psi)v\,dS=0.
\end{aligned}
\]

Lo schema Newmark implementato usa i parametri `newmarkbeta` e `newmarkgamma` (valori predefiniti \(1/4,1/2\)); il ramo richiede `alpha 1`. La positività di questi parametri è controllata, ma il codice non impone la relazione specifica necessaria per ogni proprietà di stabilità/accuratezza dello schema.

## Energia potenziale comune

Per \(k,m,\beta>0\), l'energia potenziale coerente con la forma debole è

\[
\boxed{\quad
\mathcal E_{pot}(u)=\int_\Omega\left[\frac{k}{2}|\nabla u|^2+mF(u)\right]dx
+\int_\Gamma\left[\frac{k}{2}|\nabla_\Gamma\psi|^2+\beta mF(\psi)\right]dS.
\quad}
\]

Il fattore \(\beta\) moltiplica solo il potenziale superficiale. La stessa \(k\) moltiplica energia di gradiente bulk e tangenziale: ciò segue dall'uso della medesima matrice costitutiva nel materiale di dominio e di linea.

Per una variazione \(w\), la derivata prima è

\[
\begin{aligned}
D\mathcal E_{pot}(u)[w]
={}&\int_\Omega[-k\Delta u+mF'(u)]w\,dx
+\int_\Gamma k\partial_nu\,w\,dS\\
&+\int_\Gamma[-k\Delta_\Gamma\psi+\beta mF'(\psi)]w\,dS,
\end{aligned}
\]

assumendo bordo regolare e nessun bordo della varietà \(\Gamma\) (per un dominio chiuso, \(\Gamma\) non ha bordo). Poiché \(F'=-f\), le equazioni sopra sono rispettivamente gradient flow in bulk/superficie oppure dinamica inerziale in bulk accoppiata a gradient flow superficiale. Il termine \(+k\partial_nu\) ha segno positivo: deriva dalla variazione dell'energia bulk. Non va cambiato segno isolandolo dal termine di bordo prodotto dall'integrazione per parti.

## Derivazione del bilancio parabolico

Moltiplichiamo la PDE bulk per \(u_t\), quella superficiale per \(\psi_t\), integriamo e sommiamo. Per il bulk:

\[
\int_\Omega c u_t^2 dx-\int_\Omega k\Delta u\,u_t dx
=\int_\Omega mf(u)u_t dx.
\]

L'integrazione per parti dà

\[
-\int_\Omega k\Delta u\,u_t dx
=\frac{d}{dt}\int_\Omega\frac{k}{2}|\nabla u|^2dx
-\int_\Gamma k\partial_nu\,\psi_t dS.
\]

Pertanto

\[
\frac{d}{dt}\int_\Omega\left(\frac{k}{2}|\nabla u|^2+mF(u)\right)dx
+\int_\Omega c u_t^2dx-\int_\Gamma k\partial_nu\,\psi_t dS=0,
\]

perché \(mF'(u)u_t=-mf(u)u_t\). Moltiplicando l'equazione superficiale per \(\psi_t\), usando integrazione per parti tangenziale:

\[
\frac{d}{dt}\int_\Gamma\left(\frac{k}{2}|\nabla_\Gamma\psi|^2+\beta mF(\psi)\right)dS
+\int_\Gamma c\psi_t^2dS+\int_\Gamma k\partial_nu\,\psi_t dS=0.
\]

I due termini normali si cancellano con segno opposto. Quindi

\[
\boxed{\frac{d}{dt}\mathcal E_{pot}(u)
=-\int_\Omega c u_t^2dx-\int_\Gamma c\psi_t^2dS\le0.}
\]

La quantità monotona è l'energia potenziale totale bulk più superficie; la massa \(\int_\Omega u\) non è in generale conservata.

## Derivazione del bilancio iperbolico

Moltiplichiamo la PDE bulk per \(u_t\). Il termine inerziale soddisfa

\[
\int_\Omega \rho u_{tt}u_tdx=\frac{d}{dt}\int_\Omega\frac{\rho}{2}u_t^2dx.
\]

Ripetendo l'integrazione per parti precedente:

\[
\frac{d}{dt}\left[\int_\Omega\frac{\rho}{2}u_t^2dx+\mathcal E_{pot}(u)\right]
+\int_\Omega\rho\zeta u_t^2dx-\int_\Gamma k\partial_nu\,\psi_t dS=0.
\]

La moltiplicazione dell'equazione superficiale per \(\psi_t\) dà il bilancio superficiale già derivato, con dissipazione \(\int_\Gamma c_\Gamma\psi_t^2\) e il termine normale positivo. Sommando:

\[
\boxed{\frac{d}{dt}\mathcal E_{tot}(u,u_t)
=-\int_\Omega\rho\zeta u_t^2dx-\int_\Gamma c_\Gamma\psi_t^2dS\le0,}
\]

con

\[
\mathcal E_{tot}(u,u_t)=\underbrace{\int_\Omega\frac{\rho}{2}u_t^2dx}_{\text{energia cinetica bulk}}
+\underbrace{\mathcal E_{pot}(u)}_{\text{energia potenziale bulk e superficie}}.
\]

La sola energia potenziale non è necessariamente monotona nel modello iperbolico: può convertirsi in energia cinetica. È l'energia totale a essere non crescente. Se \(\zeta=0\), la dissipazione proviene comunque dal flusso Allen–Cahn superficiale; l'energia totale è costante solo se anche \(\psi_t=0\) (o in assenza di dissipazione superficiale). Se si elimina la dinamica superficiale oppure si pone \(c_\Gamma=0\), occorre rivalutare l'ipotesi di positività e il bilancio nel limite corrispondente.

## Ipotesi necessarie e condizioni al contorno

Le identità continue richiedono soluzioni sufficientemente regolari per le integrazioni per parti e la moltiplicazione per \(u_t\) (ad esempio \(u\) con derivate spaziali seconde e temporali necessarie, traccia compatibile e \(\psi\) regolare lungo \(\Gamma\)); coefficienti costanti positivi \(c,k,m,\beta\) per il parabolico e \(\rho,c_\Gamma,k,m,\beta>0\), \(\zeta\ge0\) per l'iperbolico; nessuna forzante né lavoro al bordo. La frontiera considerata dinamica deve essere la frontiera su cui si integra la superficie.

Con bordi misti, condizioni di Dirichlet o carichi esterni compaiono termini di potenza/lavoro aggiuntivi; la monotonia non segue senza analizzarli. Spigoli poligonali sono ammessi nella forma debole con regolarità adatta, ma non necessariamente nella derivazione classica punto per punto.

## Discretizzazione ed energia numerica

`TransportElement` assembla la diffusione dalla matrice tangente del materiale, la capacità dalla risposta `Capacity` e la reazione da `computeInternalSourceVector`; `Line1_ht::computeVolumeAround` integra le linee con il fattore area della sezione. Gli elementi di linea e bulk usano gli stessi nodi di bordo: la formulazione FEM è conforme rispetto alla traccia. La forma semidiscreta, se assemblata coerentemente e senza forzanti, conserva la struttura del bilancio: nel caso parabolico \(\dot E_h=-\dot U^T C_h\dot U\); nel caso iperbolico \(\dot E_h=-\dot U^T C_{d,h}\dot U\), dove \(C_h\) include capacità bulk e di bordo nel parabolico, mentre \(C_{d,h}\) comprende smorzamento bulk e capacità superficiale nel modello misto.

Questa proprietà semidiscreta non implica automaticamente monotonia dell'energia calcolata a ogni passo. Lo schema a un passo del ramo parabolico dipende da `alpha`; in generale valori diversi dal caso totalmente implicito possono non preservare monotonia discreta per ogni passo. Newmark (predefinito \(\beta_N=1/4,\gamma_N=1/2\)) integra l'inerzia bulk, ma il codice non dimostra né impone un'identità energetica discreta per il sistema nonlineare accoppiato. La convergenza Newton del passo non prova una legge di dissipazione discreta. Il confronto numerico deve quindi distinguere: identità energetica continua, bilancio semidiscreto e andamento dell'energia ricostruita dalla soluzione a passi finiti.

## Quantità da misurare

Per il parabolico salvare separatamente energia gradiente bulk, potenziale bulk, energia gradiente superficiale e potenziale superficiale, oltre a \(\mathcal E_{pot}\), e verificare il bilancio integrato

\[
\mathcal E_{pot}(t)-\mathcal E_{pot}(0)
+\int_0^t\!\left(\int_\Omega c u_t^2dx+\int_\Gamma c\psi_t^2dS\right)ds=0.
\]

Per l'iperbolico aggiungere \(\tfrac12\int_\Omega\rho u_t^2\) e controllare

\[
\mathcal E_{tot}(t)-\mathcal E_{tot}(0)
+\int_0^t\!\left(\int_\Omega\rho\zeta u_t^2dx+\int_\Gamma c_\Gamma\psi_t^2dS\right)ds=0.
\]

Riportare residuo assoluto e relativo del bilancio, contributi di dissipazione, minimo incremento passo-passo dell'energia e sensibilità a mesh, quadratura e passo temporale. Derivate temporali discrete e quadratura devono essere dichiarate: ricostruzioni incoerenti possono generare un residuo anche quando l'identità continua è corretta.

## Incongruenze e limiti rilevati nel codice

1. **`knormal` non è un coefficiente fisico indipendente nella soluzione non MMS.** Viene inizializzato e usato nella sorgente manufactured, ma la forma debole effettivamente assemblata produce il termine normale con coefficiente `k` bulk. Impostare `knormal != k` non modifica il problema non forzato, mentre cambia la sorgente MMS: è una discrepanza fra interfaccia dei parametri e problema risolto.
2. **La capacità parabolica superficiale non è separata.** Il ramo parabolico usa `c` anche sulle linee. Solo il ramo iperbolico espone `surfacecapacity`.
3. **La diffusione superficiale non è separata.** Il coefficiente tangenziale delle linee è `k`, uguale a quello bulk; non si può riprodurre con i soli parametri attuali un modello con \(k_\Gamma\ne k\).
4. **Il fattore `surfacebeta` modifica soltanto il potenziale/reazione superficiale**, non la diffusione o la capacità. Il funzionale energetico deve quindi pesare di \(\beta\) il potenziale superficiale, non l'intera energia di superficie.
5. **Copertura della frontiera.** Ogni bordo dinamico deve essere rappresentato da elementi `line1ht` una sola volta. Bordi mancanti, duplicati o segmenti interni aggiungono/tolgono termini di energia e dissipazione. L'audit qui identifica la regola d'assemblaggio, ma non certifica automaticamente ogni input geometrico.
6. **Monotonia discreta non garantita dal solo assemblaggio.** Gli schemi temporali possono introdurre variazioni energetiche finite; non si deve interpretare ogni piccolo aumento numerico come violazione della legge continua senza studio di raffinamento, ma nemmeno attribuire al metodo una proprietà di dissipazione non dimostrata.

Di conseguenza, i bilanci in riquadro sono quelli del modello continuo coerente con la forma debole assemblata, sotto le ipotesi dichiarate. Non costituiscono da soli una verifica numerica dell'implementazione o della monotonia a passo finito.
