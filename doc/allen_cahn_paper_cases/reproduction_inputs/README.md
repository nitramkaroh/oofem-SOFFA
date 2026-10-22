# Input per i confronti delle figure

I sei input riproducono i casi di energia delle Figure 1 (Wang et al.) e 3 (Akrivis–Li) con il materiale Allen–Cahn attuale. Le mesh sono triangolari strutturate con P1; ogni input contiene nodi, elementi, condizioni iniziali e bordo.

Per eseguire un caso, copiare l'input in una cartella di lavoro vuota ed eseguire, per esempio:

```sh
/path/to/oofem -f wang_fig1_eps0p1.in
```

Le scelte di passo sono quelle usate nei grafici comparativi. Per Wang, `eps0p01` è incluso intenzionalmente: il run restituisce valori non finiti già al primo passo. Per Akrivis–Li, le curve sono esplorative: il paper usa BDF5/Radau IIA e `h=epsilon/sqrt(2)`, mentre questi input usano backward Euler e lato di griglia 0.1. Il caso `eps0p01` mostra crescita dell'energia e overshoot, quindi non è una riproduzione valida.

Per rigenerare i grafici, preparare una directory dati con la struttura `wang/eps_0.1`, `wang/eps_0.05`, `wang/eps_0.01`, `akrivis/eps_0.04`, `akrivis/eps_0.02`, `akrivis/eps_0.01`; in ciascuna cartella mettere l’input come `case.in`, eseguire OOFEM lì e poi lanciare `python3 doc/allen_cahn_paper_cases/compare_paper_figures.py --data-root /percorso/dati`. Senza opzioni, lo script legge `/tmp/ac_paperplots/`. I file VTK non sono inclusi nel repository.

## Run raffinati del confronto energetico

Gli input a 40 divisioni per lato sono generati da `python3 doc/allen_cahn_paper_cases/refine_case_inputs.py --n 40` e salvati in `doc/allen_cahn_paper_cases/refined_inputs_n40/`. I sei run richiesti per aggiornare il confronto sono stati eseguiti con il binario `build/opt_magneto/oofem`; input, log e output VTK sono in `/tmp/ac_refined/`. Per rigenerare il PDF comparativo:

```sh
python3 doc/allen_cahn_paper_cases/compare_paper_figures.py --data-root /tmp/ac_refined
```

I run con valori non finiti sono riportati solo fino all’ultimo stato finito; l’analisi OOFEM può comunque dichiarare `ANALYSIS FINISHED` dopo aver propagato NaN.
## Akrivis--Li Fig. 3 paper-parameter inputs

Generate with `python3 doc/allen_cahn_paper_cases/refine_case_inputs.py --akrivis-paper --output /tmp/ac_akrivis_paper`.
These inputs use a uniform structured P1 triangle mesh with `h=epsilon/2`, `deltat=epsilon^2`, BDF5 and full Newton updates. BDF order is ramped from 1 to 5; the paper instead uses three-stage Radau IIA for its first four steps, which is not implemented yet. Mesh nodes/connectivity are not supplied by the article, so only the mesh scale and problem data are matched. VTK fields are saved about every 1% of the run to control output size.
