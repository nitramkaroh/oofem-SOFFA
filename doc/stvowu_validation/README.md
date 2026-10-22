# Validazione Allen–Cahn con bordo dinamico

`report.pdf` contiene la verifica parabolica OOFEM diretta e il confronto con il prototipo WED Python. I backend sono distinti: il WED non è eseguito nativamente dall'eseguibile C++.

## Riproduzione

Prerequisiti Python: NumPy, SciPy e Matplotlib. I driver OOFEM usano `build/debug_magneto/oofem`.

Da root del repository:

```sh
python3 doc/stvowu_validation/run_wed_homogeneous.py
python3 doc/stvowu_validation/run_wed_nonuniform.py
python3 doc/stvowu_validation/run_wed_mesh_refinement.py
python3 doc/stvowu_validation/run_wed_horizon_sensitivity.py
python3 doc/stvowu_validation/run_wed_robustness.py
python3 doc/stvowu_validation/run_wed_small_epsilon_time_check.py
python3 doc/stvowu_validation/run_surface_beta_study.py
python3 doc/stvowu_validation/run_energy_balance.py
python3 doc/stvowu_validation/run_mms_refinement.py
python3 doc/stvowu_validation/update_report_assets.py
python3 doc/stvowu_validation/plot_recommended.py
python3 doc/stvowu_validation/plot_wed_refinements.py
python3 doc/stvowu_validation/plot_checklist_robustness.py
python3 doc/stvowu_validation/plot_paper_figures.py
cd doc/stvowu_validation && pdflatex -interaction=nonstopmode -halt-on-error report.tex
pdflatex -interaction=nonstopmode -halt-on-error numerical_section.tex
pdflatex -interaction=nonstopmode -halt-on-error numerical_section_en.tex
```

I driver producono i dati in `*.csv`, i grafici in `*.png` e i log/input dei casi archiviati in `archive/`. Le figure e le tabelle sono associate qui:

- errore vs epsilon e curva media di bordo: `plot_recommended.py` da `wed_nonuniform.csv` e `surface_beta.csv`;
- residuo energetico vs passo: `plot_recommended.py` da `energy_balance.csv`; energia nel tempo: `update_report_assets.py` da `energy_decay.in` e dai VTU del caso;
- convergenza mesh WED e passo WED: `plot_wed_refinements.py` da `wed_mesh_refinement.csv` e `wed_nonuniform_time_refinement.csv`;
- snapshot bulk/traccia: `update_report_assets.py` dagli output del run diretto non uniforme;
- tabella MMS mesh: `run_mms_refinement.py` da soluzione fabbricata analitica.

Parametri principali: WED non uniforme usa mesh $8\times4$, riferimento OOFEM $\Delta t=10^{-3}$, $S_{\max}=12$, $\Delta s=0.1$, tolleranza di ottimizzazione $2\,10^{-5}$ e massimo 100 iterazioni. La mesh WED usa $\varepsilon=0.05$ e gli stessi parametri temporali sulle griglie $8\times4$, $16\times8$, $32\times16$; la terza è il riferimento numerico fine. Quadratura: 3 punti di Gauss su `line1ht`, `nip 9` sui quadrilateri non lineari Q1. Il MMS diretto usa $T=0.02$ e $\Delta t=5\,10^{-4}$.

## Stato e limiti

- Il trend WED in epsilon è stato esteso fino a 0.0125; il passo $\Delta s$ è stato dimezzato per $\varepsilon=0.05$.
- `run_wed_robustness.py` ripete l'intera serie epsilon con passo fisico costante $\Delta t=\varepsilon\Delta s=0.005$ e confronta le mesh con $\Delta s=0.05$, tolleranza $2\,10^{-6}$.
- `run_wed_small_epsilon_time_check.py`, dopo il run di robustezza, dimezza ancora $\Delta s$ a $\varepsilon=0.0125$ mantenendo $S_{\max}$, mesh e tolleranza, e riporta sia gli errori rispetto a OOFEM sia la differenza diretta fra le due traiettorie WED.
- `numerical_section_en.pdf` è la sezione autonoma in inglese pronta per l'inserimento nel paper, con Tabelle A--F e Figure 1--6; `numerical_section.pdf` è la versione italiana. I PNG usati sono in `paper_figures/`.
- L'orizzonte è confrontato a $S_{\max}=6,12,18$ su una finestra fisica comune.
- Il raffinamento WED confronta due mesh coarse con la soluzione sulla mesh $32\times16$; quest'ultima è un riferimento numerico, non esatto.
- Il MMS verifica separatamente il solver diretto OOFEM.
- Il bilancio energetico usa differenze all'indietro per $u_t$ e matrici FEM di massa per la dissipazione; l'energia include termini bulk e superficiali.
- Gli studi beta presentano separatamente OOFEM e WED. Non è implementato un confronto con un problema privo di dinamica superficiale; porre $\beta=0$ non definirebbe da solo tale problema.
- Il caso ondoso è fuori ambito. Gli ordini osservati sono empirici sui casi riportati, non risultati teorici universali.

## Riproducibilità archiviata

Una ripetizione in directory temporanea pulita del caso WED non uniforme con `epsilon=0.05` ha riprodotto gli errori integrati bulk/bordo (`1.1511e-3`, `2.0198e-3`). Input e log sono conservati in `archive/inputs/clean_repro/` e `archive/logs/clean_repro/`; output grezzi compatti sono in `archive/results/clean_repro/`. Gli archivi per gli altri casi contengono gli input e i log OOFEM/WED disponibili. Ambiente usato: Python 3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3; revisione sorgente `0680f2da49050e2ef72407897d7b292ce66cd054`, build Debug `build/debug_magneto/oofem`. Per rigenerare il PDF serve anche `pdflatex`.

La replica pulita riguarda un caso rappresentativo e usa la build esistente: non è una build isolata da zero. I file VTU completi non sono duplicati nell'archivio; i risultati tabellari sono conservati nei CSV della cartella.
