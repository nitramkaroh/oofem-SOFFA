# Confronto OOFEM con i riferimenti Allen–Cahn a bordo dinamico

## Equazione e corrispondenza del solver

Per il potenziale convenzionale \(F(u)=\tfrac14(u^2-1)^2\), la specializzazione parabolica implementata in `allencahnmat` è

\[
 u_t-\Delta u=u-u^3\quad\text{in }\Omega,\qquad
 u_t=\Delta_\Gamma u-(u^3-u)-\partial_nu\quad\text{su }\Gamma,
\]

con traccia continua tra dominio e bordo. Questa coincide con la forma Allen–Cahn bulk/surface discussa da Kovács–Lubich per capacità superficiale, diffusione e potenziali unitari. Nel materiale OOFEM sono stati usati `c=k=mobility=surfacebeta=1`; `surfacebeta` scala la reazione di bordo. Il coefficiente indipendente della diffusione tangenziale non è attualmente configurabile: la corrispondenza parametrica diretta copre il caso unitario.

La forma debole FEM include i contributi di massa e diffusione delle linee di bordo e il flusso normale accoppiato alla soluzione bulk. Il riferimento teorico analizza elementi finiti P1, mentre i test qui usano quadrilateri Q1 di OOFEM: la struttura variazionale è la stessa, ma non la discretizzazione spaziale identica. Per la forma con parametro d'interfaccia \(\varepsilon\), il mapping dei coefficienti è `c=1`, `k=0.01`, `mobility=0.01/epsilon^2`, `surfacebeta=1`; ciò corrisponde a \(\gamma=M=0.01\) e allo stesso doppio pozzo nel bulk e sul bordo.

## Test eseguiti

1. **Dissipazione ed errore di bilancio temporale.** `run_energy_balance.py` ha eseguito OOFEM su mesh Q1 20×20 con dato iniziale liscio \(u_0(x,y)=0.1+0.05\cos(2\pi x)\), fino a \(T=0.1\), con \(\Delta t=0.005, 0.0025, 0.00125\). Nessun dato fabbricato o MMS è usato. La funzionale monitorata è \(E=\int_\Omega[\tfrac12|\nabla u|^2+F(u)]+\int_\Gamma[\tfrac12|\nabla_\Gamma u|^2+F(u)]\).
2. **Convergenza temporale.** Si confrontano i campi finali sulle stesse coordinate e con la matrice di massa bulk+bordo; \(\Delta t=0.00125\) è il riferimento numerico più fine, non una soluzione esatta.
3. **Convergenza spaziale.** `run_spatial_convergence.py` ha ripetuto la stessa evoluzione naturale sulle mesh 10×10, 20×20 e 40×40 con passo comune \(\Delta t=0.0025\), fino a \(T=0.1\). I valori nodali delle mesh annidate sono confrontati con la mesh 40×40, che è un riferimento numerico.

I risultati sono in `metrics.csv`, `energy_trajectory.csv`, `energy_decay.png` e `spatial_convergence.csv`. Le esecuzioni grezze sono mantenute in `/tmp/stvowu_energy_balance` e `/tmp/stvowu_literature_spatial`; si possono rigenerare con:

```sh
python3 doc/stvowu_validation/run_energy_balance.py
python3 doc/stvowu_validation/literature_comparison/postprocess.py
python3 doc/stvowu_validation/literature_comparison/run_spatial_convergence.py
```

## Risultati sintetici

- Energia iniziale: 1.28445597. A \(T=0.1\), l'energia vale 1.20911974, 1.20901760 e 1.20896702 per i tre passi temporali. In nessuno dei passi temporali l'energia discreta misurata aumenta.
- Il residuo del bilancio energetico normalizzato scende da 4.23e-3 a 1.20e-3 quando \(\Delta t\) è ridotto di un fattore quattro; il residuo è una misura diagnostica del bilancio completamente discreto, non una violazione osservata della monotonia dell'energia.
- Il rapporto delle differenze tra campi finali successivi in norma di massa è 1.978, che corrisponde all'ordine temporale empirico 0.984, compatibile con il backward Euler di OOFEM.
- Rispetto alla mesh 40×40, l'RMS nodale è 1.279e-4 sulla 10×10 e 2.590e-5 sulla 20×20; i massimi scarti sono rispettivamente 2.075e-4 e 4.230e-5. Il rapporto tra i due RMS è 4.94 (ordine apparente 2.30), indicativo ma non una stima rigorosa perché la mesh 40×40 è anch'essa approssimata.
- La massa Allen–Cahn non è conservata: in questo run passa da 0.6000 a circa 0.6620. Non si tratta del modello conservativo/penalizzato di alcuni lavori recenti.

## Confronto con gli articoli forniti

- **Okumura, dissertation (2021), Chapter 4, Computation Example 2:** questa nuova fonte fornisce un benchmark naturale completamente parametrizzato per Allen–Cahn 1D con bordo dinamico. Usa \(L=1\), \(K=100\) (\(\Delta x=0.01\)), capacità dinamica agli estremi \(\epsilon_{ex}=10\), potenziale \(F(s)=\frac{\gamma}{4}(s^2-1)^2\) con \(\gamma=100\), dato iniziale \(u_0(x)=\exp[-500(x-0.5)^2]\), \(\Delta t=10^{-3}\), \(T=0.7\). Non è MMS: il dato iniziale è naturale e il lavoro pubblica sia il campo evoluto (Fig. 4.3) sia l'energia totale (Fig. 4.4).

  Abbiamo riprodotto questo caso con il solver OOFEM 2D su una striscia di profondità unitaria, 100 Q1 lungo \(x\) e una sola fila lungo \(y\). Il campo iniziale non dipende da \(y\); i bordi dinamici sono assegnati solo ai due lati \(x=0,1\), mentre i lati superiore e inferiore hanno Neumann naturale. La riga di materiale degli estremi usa capacità `c=10`, il bulk `c=1`; entrambi usano `k=1`, `mobility=100`. Così la soluzione FEM si riduce esattamente al problema 1D del capitolo, senza cambiare le equazioni né aggiungere termini sorgente.

  Il campo evolve dal picco gaussiano fino alla fase uniforme \(+1\), come nella Fig. 4.3; la variazione massima lungo la profondità estrusa è zero alla precisione dei dati esportati. L'energia totale bulk più i due estremi scende da 77.42 al primo output fino a \(2.22\times10^{-6}\) a \(T=0.7\), senza alcun incremento tra passi. Abbiamo digitalizzato la curva pubblicata in Fig. 4.4 e confrontato il tratto \(0.01\le t\le0.4\): RMSE 0.988, MAE 0.880 e massimo scarto 1.832, su un asse energetico 0–90. La digitalizzazione raster ha un'incertezza dell'ordine di alcune decime di unità; l'accordo delle curve è molto stretto e le piccole differenze sono compatibili con la differenza tra differenze finite e Q1 FEM e con la lettura della figura.

  Artefatti: `run_okumura_1d_strip.py`, `analyze_okumura_1d_strip.py`, `digitize_okumura_fig44.py`, `okumura_ch4_case2.in`, `okumura_ch4_case2_final.vtu`, `okumura_ch4_example2_oofem.png`, `okumura_ch4_example2_energy_overlay.png`, `okumura_ch4_example2_energy.csv`, `okumura_ch4_example2_paper_fig44_digitized.csv`, `okumura_ch4_example2_literature_comparison.txt`. Il run si ripete con `python3 doc/stvowu_validation/literature_comparison/run_okumura_1d_strip.py`; poi lanciare i due script di analisi e digitalizzazione.

- **Okumura–Furihata (2020):** è il confronto più vicino per il principio strutturale: Allen–Cahn con bordo dinamico, schema a differenze finite che preserva la dissipazione e figure di evoluzione/energia. I test OOFEM confermano qualitativamente la dissipazione e la coerenza temporale. Non è stata ottenuta una replica punto-a-punto delle figure: la pagina pubblica consultabile espone abstract e figure, ma il testo integrale (che contiene dati iniziali e parametri) reindirizza a una pagina di accesso. Di conseguenza non si attribuiscono ai nostri test gli stessi dati iniziali dell'articolo.
- **Kovács–Lubich (2017):** fornisce il quadro FEM bulk-superficie e tratta l'equazione Allen–Cahn con condizione dinamica. Contiene anche esperimenti numerici, ma nel caso Allen–Cahn scelto usa la soluzione fabbricata \(u=e^{-4t}x_1x_2\) sul disco unitario con sorgenti costruite appositamente. Come richiesto, non usiamo quel caso come validazione fisica del dato iniziale naturale; resta un riferimento di discretizzazione e analisi. La forma debole è pertinente, ma i coefficienti dei potenziali e della diffusione vanno mappati con cura (l'articolo consente potenziali generali e non fissa un unico benchmark nonlineare naturale).
- **Xiao–Chen (2025):** il PDF conferma che il lavoro contiene sia modelli penalizzati sia risultati *senza* i termini di penalità. In particolare, il Caso 4 usa \(\Omega=[0,1]^2\), una goccia quadrata centrata in \((0.5,0.25)\) di lato 0.5 e mostra snapshot a \(T=0.1,0.5,1,2,4,6\) per il ramo senza penalità (Fig. 5.19), oltre all'energia (Fig. 5.21). Questo ramo è strutturalmente compatibile con il nostro Allen–Cahn non conservativo; il ramo con penalità non lo è, perché OOFEM non implementa i termini globali che impongono la conservazione separata delle masse.

  **Il PDF omette \(\varepsilon\) per il Caso 4.** L'elenco dei parametri numerici (5.1) indica \(M=\gamma=0.01\), \(\tau=10^{-4}\), \(h=1/300\) e i parametri di penalità, ma non \(\varepsilon\). Il potenziale classico (2.7) dipende esplicitamente da \(\varepsilon\), e le figure non forniscono dati tabulati. Perciò non è possibile una replica a parametri identici né una validazione quantitativa contro il paper.

## Replica esplorativa del Caso 4 non penalizzato

Non ci siamo fermati all'assenza di \(\varepsilon\): abbiamo eseguito il Caso 4 con l'ipotesi dichiarata \(\varepsilon=0.1\), usando il mapping sopra e il dato nodale letto da Fig. 5.18 (\(+1\) nella goccia quadrata che tocca il bordo inferiore, \(-1\) altrove). La discretizzazione Q1 OOFEM è su \(20\times20\) celle (\(h=0.05\), contro \(h=1/300\) nel paper); il run principale usa \(\Delta t=10^{-4}\), \(T=0.1\). La scelta di \(\varepsilon\) è una nostra ipotesi esplicita, non un dato attribuito agli autori.

- La forma a \(T=0.1\) è coerente con Fig. 5.19: il contorno interpolato \(u=0.5\) taglia il bordo inferiore a \(x=0.2578,0.7422\) e ha sommità circa \(y=0.4922\); nella figura pubblicata si leggono approssimativamente \(x=0.25,0.75\) e \(y=0.5\). Questo è un accordo **qualitativo e visivo**, non un errore quantitativo contro dati esatti.
- Il controllo temporale specifico del caso ha trovato differenza RMS nodale \(1.35\times10^{-4}\) fra \(\Delta t=5\times10^{-4}\) e \(2.5\times10^{-4}\), e \(8.12\times10^{-5}\) fra \(2.5\times10^{-4}\) e \(10^{-4}\). Invece \(\Delta t=10^{-3}\) produce una soluzione fortemente asimmetrica (errore RMS 0.182 rispetto a \(5\times10^{-4}\), errore massimo 0.758): quel passo non è accettabile, benché Newton converga. Il passo usato dal paper, \(10^{-4}\), evita questo difetto temporale.
- La verifica spaziale a \(\Delta t=5\times10^{-4}\) fra \(20\times20\) e \(40\times40\) dà RMS nodale 0.0788 e massimo 0.258; il contorno \(u=0.5\) si sposta di circa 0.010 nelle coordinate d'intersezione. La mesh \(20\times20\) non è quindi spazialmente convergente. La mesh più fine si avvicina alla forma di Fig. 5.19, ma servono ulteriori raffinamenti per stabilire il limite spaziale.
- Con la stessa ipotesi \(\varepsilon=0.1\), l'energia Q1 calcolata con quadratura di Gauss diminuisce da 144.16 a 49.67 tra i primi output e \(T=0.1\), senza incrementi per il run a \(\Delta t=10^{-4}\). Non confrontiamo i valori assoluti con Fig. 5.21 perché \(\varepsilon\) non è noto e la mesh è molto più grossolana.

Gli input di riproduzione sono in `run_xiao_chen_droplet.py` e `xiao_chen_case4_n20_dt0p0001_eps0p1.in`; il campo finale è archiviato in `xiao_chen_case4_n20_dt0p0001_eps0p1_final.vtu`. Il grafico e i dati estratti sono in `xiao_chen_case4_exploratory_n20_dt0.0001.png`, `xiao_chen_case4_contour_n20_dt0.0001.csv`, `xiao_chen_case4_time_refinement.csv`, `xiao_chen_case4_space_refinement.csv` e `xiao_chen_case4_energy_n20_dt0p0001_eps0p1.csv`. Per ripetere il run OOFEM con il valore esplorativo dichiarato: `python3 doc/stvowu_validation/literature_comparison/run_xiao_chen_droplet.py --n 20 --dt 0.0001 --eps 0.1`. I risultati temporanei completi VTU erano sotto `/tmp/stvowu_xiao_chen_case4_*`.

## Esito della validazione

**Esito aggiornato: sì, il solver parabolico Allen–Cahn di OOFEM è validato contro un benchmark naturale e parametrizzato di letteratura nel caso 1D (replicato tramite la striscia 2D uniforme nello spessore).** La prova più forte è la curva di energia di Okumura: sovrapposizione molto stretta con Fig. 4.4, energia monotona, parametri e dati iniziali replicati.

Questa validazione copre il bulk Allen–Cahn, la capacità dinamica e la reazione ai due estremi, e il loro accoppiamento nella forma debole. Non certifica ancora tutti i comportamenti del problema 2D con variazioni lungo il bordo e Laplaciano tangenziale non nullo. Il test della goccia Xiao–Chen resta una verifica 2D qualitativa perché il paper non specifica \(\varepsilon\) per il Caso 4 e la raffinazione spaziale è ancora sensibile. Il PDF di Makki–Miranville–Petcu descrive invece un sistema accoppiato Cahn–Hilliard/Allen–Cahn con più variabili e leggi di conservazione; non è il modello scalare implementato da `allencahnmat`, quindi non lo presentiamo come confronto diretto.

## Riferimenti

- M. Okumura, D. Furihata, *A structure-preserving scheme for the Allen–Cahn equation with a dynamic boundary condition*, DCDS 40(8), 4927–4960 (2020), DOI: 10.3934/dcds.2020206.
- B. Kovács, C. Lubich, *Numerical analysis of parabolic problems with dynamic boundary conditions*, IMA J. Numer. Anal. 37(1), 1–39 (2017), DOI: 10.1093/imanum/drw015.
- M. Xiao, R. Chen, *A first-order energy stable scheme for the Allen–Cahn equation with the Allen–Cahn type dynamic boundary condition*, JCAM 460, 116409 (2025), DOI: 10.1016/j.cam.2024.116409.
- M. Okumura, *Design and analysis of structure-preserving schemes for parabolic partial differential equations with dynamic boundary conditions*, Osaka University doctoral dissertation (2021), [institutional record](https://doi.org/10.18910/82279).
- A. Makki, A. Miranville, M. Petcu, *A numerical analysis of the coupled Cahn–Hilliard/Allen–Cahn system with dynamic boundary conditions*, International Journal of Numerical Analysis and Modeling 19(5) (2022).
