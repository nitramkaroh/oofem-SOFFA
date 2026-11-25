# Piano di implementazione del calcolo energetico FEM

## Scopo e modello a cui si riferisce

Questo piano descrive un post-processing esterno della soluzione già calcolata. Non cambia le equazioni né la procedura del solver. Le formule seguono il bilancio in [`energy_theory.md`](energy_theory.md): i nodi degli elementi `line1ht` sono gli stessi nodi bulk, quindi la variabile superficiale è la traccia di `Temperature`.

Per gli input non forzati, i coefficienti da leggere dal record `allencahn` sono `k`, `mobility` \(m\), `surfacebeta` \(\beta\), `c`, `rho`, `surfacecapacity` e il flag `hyperbolic`. Il parametro temporale `wavedamping` \(\zeta\) appartiene al record del problema, non al materiale. I valori omessi hanno i default dichiarati in `allencahnmat.h` e `nltransienttransportproblem.h`. Il campo `d` presente in alcuni vecchi input non è letto da `AllenCahnMaterial`: per la diffusione il codice legge `k`. Per riproducibilità lo script dovrà estrarre i valori effettivi dal record e registrare anche i default assunti.

Il parametro `knormal` non entra nel bilancio della soluzione non MMS: il termine normale si genera dalla forma debole bulk con coefficiente `k`. `knormal` compare nella sorgente manufactured. Non si deve aggiungere all’energia un termine separato con derivata normale.

## Contributi energetici da integrare

Per ogni istante, con \(F(u)=\tfrac14(u^2-1)^2\):

| Contributo | Formula | Dominio e dati FEM |
|---|---|---|
| Gradiente bulk | \(E_{\nabla,b}=\int_\Omega \tfrac{k}{2}|\nabla u_h|^2\,dx\) | Elementi di volume; campo nodale `Temperature`; gradienti \(B u_e\); Jacobiano e thickness/volume della sezione. |
| Potenziale bulk | \(E_{F,b}=\int_\Omega \tfrac{m}{4}(u_h^2-1)^2\,dx\) | Elementi di volume; \(u_h=N u_e\); stesso peso geometrico del bulk. |
| Gradiente superficiale | \(E_{\nabla,\Gamma}=\int_\Gamma \tfrac{k}{2}|\nabla_\Gamma\psi_h|^2\,dS\) | Solo elementi `line1ht` che rappresentano la frontiera dinamica; valori nodali condivisi; gradiente lungo la linea. |
| Potenziale superficiale | \(E_{F,\Gamma}=\int_\Gamma \tfrac{\beta m}{4}(\psi_h^2-1)^2\,dS\) | Gli stessi elementi di frontiera una sola volta; \(\psi_h=N_\Gamma u_e\). |
| Cinetica bulk iperbolica | \(E_{kin}=\int_\Omega\tfrac{\rho}{2}v_h^2\,dx\), \(v_h=u_{h,t}\) | Solo elementi di volume e velocità nodale Newmark. Non esiste energia cinetica superficiale. |

L’energia potenziale totale è la somma dei primi quattro contributi. Nel ramo parabolico coincide con l’energia totale. Nel ramo iperbolico l’energia totale è \(E_{tot}=E_{kin}+E_{pot}\). Capacità \(c\) e \(c_\Gamma\) non pesano l’energia: moltiplicano la dissipazione. La dissipazione continua attesa è \(\int_\Omega c u_t^2+\int_\Gamma c\psi_t^2\) nel parabolico e \(\int_\Omega\rho\zeta u_t^2+\int_\Gamma c_\Gamma\psi_t^2\) nell’iperbolico.

## Assemblaggio OOFEM da cui ricavare i coefficienti e i controlli

| Parte | Routine/file | Cosa conferma per il post-processing |
|---|---|---|
| Coefficienti e reazione | `AllenCahnMaterial::initializeFrom`, `giveCharacteristicValue`, `computeInternalSourceVector` in `src/tm/Materials/allencahnmat.C` | `k` è la diffusività; `mobility` moltiplica \(u-u^3\); sulle linee il fattore aggiuntivo è `surfacebeta`; capacità e ramo iperbolico distinguono bulk e linea. |
| Flusso e tangente diffusiva | `computeFlux3D`, `computeTangent3D` in `allencahnmat.C` | \(q=-k\nabla u\); la matrice costitutiva restituisce \(kI\), anche sugli elementi di linea. |
| Massa/capacità | `TransportElement::computeCapacitySubMatrix` in `src/tm/Elements/transportelement.C`; `giveCharacteristicValue` | La capacità FEM è integrata con \(N^TN\). Per energia cinetica usare `rho` sulla regione bulk, senza assegnare massa inerziale alle linee. |
| Diffusione | `TransportElement::computeConductivitySubMatrix` | Integra \(B^TkB\). È un controllo utile per il termine di gradiente, ma non fornisce direttamente il valore scalare dell’energia a un dato istante. |
| Reazione assemblata | `computeInternalForcesVector`, `computeInternalSourceRhsSubVectorAt`, `computeIntSourceLHSSubMatrix` | Verificano il segno e la quadratura della reazione. Il residuo usa la sorgente come termine di carico positivo. L’energia potenziale va comunque calcolata dalla formula di \(F\), non integrando il residuo Newton. |
| Misura bulk | `Quad1_ht::computeVolumeAround` / `Tr1_ht::computeVolumeAround` (e rispettivi elementi) | Peso \(|J|w\) moltiplicato per `CS_Thickness`; per altri elementi 3D/assialsimmetrici usare la loro propria misura geometrica. |
| Misura superficiale | `Line1_ht::computeVolumeAround` in `src/tm/Elements/line1_ht.C` | Peso \(|J_\Gamma|w\) moltiplicato per `CS_Area`. Questo fattore va incluso in entrambi i contributi di linea. |
| Velocità/accelerazione | `NLTransientTransportProblem::giveUnknownComponent`, `updateWaveKinematics` in `src/tm/EngineeringModels/nltransienttransportproblem.C` | Nel ramo Newmark la velocità è uno stato interno del problema, aggiornato con `newmarkgamma`; non coincide in generale con una derivata centrata approssimata del solo campo esportato. |

Le routine di assemblaggio di massa, diffusione e reazione sono riportate per associare i coefficienti e i segni; non vanno usate come shortcut del post-processing energetico. In particolare, la matrice di capacità non è l’energia cinetica e la matrice di diffusione non è l’energia di gradiente finché non si valuta \(\tfrac12 U^TKU\) con lo stesso coefficiente e la stessa quadratura.

## Recupero dei campi

### Dati attualmente disponibili

I file VTU dei casi esistenti ispezionati riportano `PointData/Temperature` e le coordinate dei nodi; non riportano `Velocity` o `Acceleration`. Gli input (`.in`) contengono mesh, tipologia di elemento, cross section e materiale. Il file PVD associa ogni frame al tempo. Perciò:

- i quattro termini potenziali/di gradiente si possono ricostruire dai VTU e dall’input;
- la cinetica iperbolica **non si può ricostruire esattamente dal solo VTU corrente**;
- `initialvelocity` / `initialvelocityvalues` danno solo la condizione iniziale, non la storia di \(v_h\).

Il file `.out` ordinario non va assunto come archivio di velocità: la configurazione VTK esaminata esporta soltanto la variabile primaria. Per energia cinetica affidabile serve un’esportazione della velocità a ogni frame o un postprocessor eseguito dentro OOFEM che interroghi `giveUnknownComponent(VM_Velocity, ...)`. Un’alternativa provvisoria è stimare \(v_h\) dai campi a tempi consecutivi, con formule coerenti e accuratezza dichiarata; questa derivata può degradare l’ordine, specialmente ai bordi temporali, e non è la scelta preferita per un bilancio fine. L’accelerazione non serve per l’energia cinetica.

### Ricostruzione sui punti di Gauss

Per ogni elemento bulk, ricostruire \(u_h=N(\xi)u_e\) e \(\nabla u_h=B(\xi)u_e\) dai valori `Temperature` nell’ordine locale dei nodi definito dalla connettività OOFEM. Nel ramo iperbolico interpolare analogamente le velocità nodali \(v_h=Nv_e\). Per ogni elemento linea, parametrizzare la linea e calcolare \(\psi_h=N_\Gamma u_e\), \(\nabla_\Gamma\psi_h=B_\Gamma u_e\); su una linea retta quest’ultimo è la derivata rispetto all’arco. L’orientamento della linea può cambiare il segno del gradiente, ma non i contributi energetici quadratici.

Per ogni punto di quadratura accumulare:

\[
\begin{array}{ll}
E_{\nabla,b} \mathrel{+}= \frac{k}{2}(\nabla u_h\cdot\nabla u_h)|J|w\,t, &
E_{F,b} \mathrel{+}= \frac{m}{4}(u_h^2-1)^2|J|w\,t,\\[2mm]
E_{\nabla,\Gamma} \mathrel{+}= \frac{k}{2}|\nabla_\Gamma\psi_h|^2|J_\Gamma|w\,A, &
E_{F,\Gamma} \mathrel{+}= \frac{\beta m}{4}(\psi_h^2-1)^2|J_\Gamma|w\,A,\\[2mm]
E_{kin} \mathrel{+}= \frac{\rho}{2}v_h^2|J|w\,t. &
\end{array}
\]

Qui \(t\) è lo spessore/cross-section thickness richiesto dall’elemento bulk 2D e \(A\) è `CS_Area` della linea. In 3D la misura bulk è il volume 3D; non moltiplicare nuovamente per thickness. Non confondere area di sezione della linea e spessore del bulk.

## Quadratura

La quadratura del postprocessor deve seguire la geometria e l’ordine FEM e deve essere almeno sufficientemente accurata per il potenziale quartico. Per un quadrilatero Q1 affine rettangolare, \(u_h\) è bilineare e \((u_h^2-1)^2\) può avere grado quattro in ciascuna coordinata: una regola Gauss 3×3 è esatta per questo integrando; 2×2 non lo è in generale. Per elementi Q1 non affini, il gradiente contiene l’inversa del Jacobiano e non è in generale polinomiale: usare più punti e documentare uno studio di convergenza di quadratura. Per una linea lineare, il potenziale quartico ha grado al più quattro in coordinata locale: Gauss a 3 punti è esatta; una quadratura inferiore può alterare l’energia superficiale.

Per la riproduzione dell’energia discreta, registrare anche il numero effettivo di punti impostato con `nip` e la regola usata dall’elemento. In `TransportElement`, le matrici di capacità/conducibilità e le forze interne usano `integrationRulesArray[0]`, mentre il vettore della sorgente assemblato da `computeInternalSourceRhsSubVectorAt` itera `giveDefaultIntegrationRulePtr()`. Queste vie possono non avere identica quadratura: l’implementazione dell’energia deve riportare sia una quadratura di post-processing accurata sia, se si confronta la derivata dell’energia con il residuo FEM, la differenza rispetto alla quadratura effettivamente usata nell’assemblaggio della reazione.

## Regole per non contare due volte

1. Integrare gradienti/potenziale bulk solo sugli elementi di area/volume.
2. Integrare contributi superficiali solo sugli elementi `line1ht` designati come frontiera dinamica; verificare che ogni segmento sia presente una volta e che non ci siano linee interne o bordi duplicati.
3. Non reintegrare le facce del bordo degli elementi bulk: il bulk contribuisce già una volta attraverso i suoi elementi volumetrici.
4. Non aggiungere \(k\partial_nu\) come energia o come ulteriore carico di bordo. È l’accoppiamento nella variazione dell’energia bulk e si cancella nella derivazione del bilancio, non è un quinto addendo energetico.
5. Le coordinate/nodi condivisi tra elementi bulk e linee non significano che i valori superficiali vadano sommati ai valori bulk nei nodi. Sono campioni dello stesso campo; i contributi si sommano solo dopo integrazione con le rispettive misure.
6. Non assegnare energia cinetica alle linee nel modello iperbolico misto.

## Ambiguità da risolvere o documentare nei dati di input

- Verificare per ogni mesh se `line1ht` ricopre l’intero bordo o soltanto porzioni dinamiche; l’integrale energetico superficiale segue esattamente quelle linee.
- Registrare `thickness`, `area`, unità e convenzione geometrica della sezione. Questi fattori cambiano l’energia assoluta.
- Registrare il mapping fra `Temperature` e il campo matematico \(u\); evitare normalizzazioni o conversioni non lineari dei valori esportati.
- Accertare che gli snapshot siano stati scritti per ogni passo o annotare eventuali frame mancanti, necessari per stime temporali.
- Per il ramo iperbolico ottenere velocità nodali Newmark. Nei dof vincolati, velocità/condizioni cinematiche devono essere quelle effettivamente applicate dal modello.
- Se esistono coefficienti diversi per materiale/regione, l’energia va integrata elemento per elemento col coefficiente locale, non con un valore globale estratto dal primo materiale.
- `mobility` è il coefficiente della reazione così come implementata. Non sostituirlo con un parametro di interfaccia \(\varepsilon\) senza una mappatura fisica derivata dall’input.

## Implementazione minima proposta (fase successiva)

Creare uno script standalone di post-processing, senza modificare il solver, che:

1. legga l’input OOFEM per nodi, connettività, classi di elemento, cross section, materiali e parametri temporali;
2. legga il PVD/VTU per associare ogni frame al tempo e ottenere coordinate e `Temperature`;
3. ricostruisca shape functions e gradienti per ogni classe di elemento supportata (inizialmente `quad1ht` + `line1ht`, poi gli altri tipi effettivamente usati);
4. integri i cinque contributi con regole dichiarate, coefficienti per materiale e fattori di misura coerenti con `computeVolumeAround`;
5. legga la velocità da un’esportazione dedicata, se disponibile; in assenza, scriva `kinetic_energy=NA` invece di presentare una stima come dato esatto;
6. esporti CSV con tempo, cinque contributi, energia potenziale, energia totale e metadati di mesh, quadratura e coefficienti;
7. faccia controlli automatici: copertura/duplicati delle linee, coerenza delle coordinate PVD-VTU, energia non negativa per coefficienti positivi, test d’invarianza rispetto all’orientamento linea e confronto della quadratura a ordini crescenti.

Il primo rilascio minimo può calcolare i quattro contributi potenziali dal VTU e CSV; il calcolo dell’energia cinetica iperbolica va abilitato soltanto dopo avere una velocità nodale verificabile. Nessuno dei passaggi richiede di cambiare la fisica o il comportamento di OOFEM.

## Prima implementazione disponibile

È stato aggiunto [`postprocess_energy.py`](postprocess_energy.py). Legge il nome base dei risultati dalla prima riga dell’input OOFEM, quindi carica il PVD e i relativi VTU ASCII. Integra `quad1ht` e `line1ht` con Gauss 3×3 nel bulk e Gauss 3 punti sulle linee. La prima versione richiede un solo `allencahn`, un solo `SimpleTransportCS`, `thickness` esplicito e mesh bidimensionali con celle VTU quad/line. Non avvia OOFEM e rifiuta di sovrascrivere un CSV esistente.

Comando:

```bash
python3 doc/stvowu_validation/postprocess_energy.py \
  doc/stvowu_validation/hyperbolic/coupled_dt0.001/case.in \
  --output doc/stvowu_validation/energy_components_hyperbolic_example.csv
```

Le formule integrate sono quelle nella tabella sopra. I valori sono espressi nelle unità coerenti usate nell’input OOFEM, senza conversione o fattore di normalizzazione aggiuntivo; in 2D le energie bulk sono moltiplicate per `thickness` e quelle superficiali per `area` della sezione. Il CSV contiene `bulk_gradient`, `bulk_potential`, `surface_gradient`, `surface_potential`, `kinetic`, `potential_energy` e `total_energy`. Nel parabolico `kinetic` è vuoto/non applicabile e `total_energy=potential_energy`. Nel modello iperbolico `kinetic` e `total_energy` sono vuoti se non è disponibile un campo di velocità scalare, con `kinetic_status=missing_velocity`.

Esecuzioni minime su risultati già esistenti, senza simulazioni:

- `parabolic_2d/runs/orientation_base_n8_dt0p001/case.in` → [`energy_components_parabolic_example.csv`](energy_components_parabolic_example.csv), un frame; energia totale potenziale risultante \(6.1992360950\times10^{-1}\).
- `hyperbolic/coupled_dt0.001/case.in` → [`energy_components_hyperbolic_example.csv`](energy_components_hyperbolic_example.csv), 20 frame; i quattro contributi potenziali sono disponibili, ma `Velocity` non è presente nei VTU e cinetica/totale restano vuote.

Limiti di questa versione: un solo materiale/sezione, elementi bidimensionali quad/line ASCII, campo primario scalare chiamato `Temperature`, quadratura fissa 3×3/3 punti e nessuna lettura della velocità Newmark da `.out`. Per chiudere l’energia totale iperbolica occorre esportare la velocità nodale effettiva a ogni frame (oppure aggiungere in seguito un postprocessor OOFEM non invasivo). Il CSV attuale non stima tale campo da differenze temporali.
