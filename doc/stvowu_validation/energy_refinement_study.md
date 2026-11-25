# Studio temporale di convergenza del bilancio energetico

## Scopo e criterio

La legge continua non forzata è quella in [energy_theory.md](energy_theory.md). Nel ramo parabolico l'energia potenziale soddisfa \(dE_{pot}/dt=-D\le0\); nel ramo iperbolico è il totale \(E_{tot}=E_{pot}+E_{kin}\) a soddisfare \(dE_{tot}/dt=-D_{bulk}-D_\Gamma\le0\). La teoria continua non garantisce da sola monotonia discreta per ogni schema e passo: qui il confronto si basa sulle curve e sul residuo integrato

\[
R(t)=E(t)-E(0)+D_{bulk}(t)+D_\Gamma(t),
\]

non sulla sola monotonia.

## Casi e configurazioni

Sono stati riutilizzati risultati di simulazioni già presenti; per questo studio non è stato lanciato alcun run. Per ciascun ramo sono considerate esattamente tre configurazioni, tutte sulla stessa mesh 16×8, sul dominio 2×1 e fino a \(T=0.02\). I relativi input sono:

| Ramo | \(\Delta t\) | Passi | Input esistente |
|---|---:|---:|---|
| Parabolico | 0.002 | 10 | [case.in](energy_validation_runs_highprecision/parabolic_n16x8_dt0p002/case.in) |
| Parabolico | 0.001 | 20 | [case.in](energy_validation_runs_highprecision/parabolic_n16x8_dt0p001/case.in) |
| Parabolico | 0.0005 | 40 | [case.in](energy_validation_runs_highprecision/parabolic_n16x8_dt0p0005/case.in) |
| Iperbolico | 0.002 | 10 | [case.in](energy_validation_runs_highprecision/hyperbolic_n16x8_dt0p002/case.in) |
| Iperbolico | 0.001 | 20 | [case.in](energy_validation_runs_highprecision/hyperbolic_n16x8_dt0p001/case.in) |
| Iperbolico | 0.0005 | 40 | [case.in](energy_validation_runs_highprecision/hyperbolic_n16x8_dt0p0005/case.in) |

All'interno di ciascun ramo sono invariati mesh, dati iniziali, coefficienti, condizioni al contorno e tolleranze; cambiano solo `nsteps` e `deltat`. La condizione iniziale è \(u_0=0.2+0.1\cos(\pi x/2)\cos(\pi y)\), senza sorgenti esterne. Il ramo parabolico usa `alpha=1`; quello iperbolico usa `rho=1`, `surfacecapacity=1`, `wavedamping=0.2`, `newmarkbeta=0.25` e `newmarkgamma=0.5`. Non sono stati raffinati spazio o parametri.

## Metodo e dati

Le storie nodali salvate sono state rielaborate con le routine energetiche FEM già esistenti in [run_energy_validation.py](run_energy_validation.py): integrazione Gauss a 5 punti per gli addendi energetici, dissipazione parabolica da differenze all'indietro, e dissipazione iperbolica con velocità Newmark ricostruita e interpolazione lineare per passo. Sono confrontati 11 istanti comuni \(t=0,0.002,\ldots,0.02\). Il file dati contiene tutti i contributi separati, le dissipazioni cumulative, il residuo e gli scarti rispetto alla serie \(\Delta t=0.0005\): [CSV](energy_refinement_study.csv). Il grafico mostra le curve totali, i residui normalizzati e gli errori delle curve rispetto al passo più fine: [grafico](energy_refinement_study.png).

## Bilancio energetico

| Ramo | \(\Delta t\) | \(E(0)\) | \(E(T)\) | \(D_{bulk}(T)\) | \(D_\Gamma(T)\) | \(R(T)\) | \(\max_t|R(t)|\) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Parabolico | 0.002 | 1.9316061622 | 1.9061825476 | 0.0066730036 | 0.0184559281 | -2.94683e-4 | 2.94683e-4 |
| Parabolico | 0.001 | 1.9316061622 | 1.9061163082 | 0.0067700730 | 0.0185673011 | -1.52480e-4 | 1.52480e-4 |
| Parabolico | 0.0005 | 1.9316061622 | 1.9060828813 | 0.0068209484 | 0.0186246690 | -7.76636e-5 | 7.76636e-5 |
| Iperbolico | 0.002 | 1.9316061622 | 1.9299595493 | 1.13363e-5 | 0.0016389684 | 3.69173e-6 | 3.69173e-6 |
| Iperbolico | 0.001 | 1.9316061622 | 1.9299554140 | 1.13479e-5 | 0.0016403160 | 9.15565e-7 | 9.15565e-7 |
| Iperbolico | 0.0005 | 1.9316061622 | 1.9299545074 | 1.13508e-5 | 0.0016406514 | 3.47269e-7 | 3.67927e-7 |

Il residuo parabolico finale, normalizzato per \(E(0)\), scende da \(1.526\cdot10^{-4}\) a \(7.894\cdot10^{-5}\) e \(4.021\cdot10^{-5}\). Nel ramo iperbolico scende da \(1.911\cdot10^{-6}\) a \(4.740\cdot10^{-7}\) e \(1.798\cdot10^{-7}\). Il residuo iperbolico più fine ha un massimo leggermente precedente a \(T\), perciò \(\max|R|\) è riportato separatamente dal valore finale.

## Differenza tra curve e ordine indicativo

Le norme seguenti sono massimi assoluti sulle 11 istantanee comuni. L'ordine è ricavato dal rapporto tra le differenze successive, \(p=\log_2(\|E_{dt}-E_{dt/2}\|/\|E_{dt/2}-E_{dt/4}\|)\); è un indicatore limitato a questi tre livelli, non un ordine certificato.

| Ramo / quantità | \(\|dt-dt/2\|_\infty\) | \(\|dt/2-dt/4\|_\infty\) | Rapporto | \(p\) indicativo |
|---|---:|---:|---:|---:|
| Parabolico, energia totale | 6.624e-5 | 3.343e-5 | 1.982 | 0.987 |
| Iperbolico, energia totale | 4.135e-6 | 9.066e-7 | 4.561 | 2.189 |
| Iperbolico, energia potenziale | 3.250e-6 | 8.122e-7 | 4.001 | 2.000 |

Per il residuo finale si confrontano invece i moduli \(|R(T)|\) direttamente:

| Ramo | \(|R(T,dt)|\) | \(|R(T,dt/2)|\) | \(|R(T,dt/4)|\) | Ordine dai primi due / dagli ultimi due |
|---|---:|---:|---:|---:|
| Parabolico | 2.947e-4 | 1.525e-4 | 7.766e-5 | 0.950 / 0.974 |
| Iperbolico | 3.692e-6 | 9.156e-7 | 3.473e-7 | 2.012 / 1.398 |

**Interpretazione.** Nel parabolico le differenze delle curve e il residuo calano quasi linearmente con \(\Delta t\), coerentemente con un comportamento di primo ordine per questo test. Nell'iperbolico la parte potenziale mostra un rapporto prossimo a 4, compatibile con secondo ordine; per la curva totale la stima è \(p\approx2.19\), ma la componente cinetica ricostruita influenza il rapporto. Inoltre la riduzione del residuo non mantiene lo stesso rapporto tra gli ultimi due livelli. Il dato suggerisce una riduzione dell'errore temporale, ma non basta a stabilire un regime asintotico del bilancio iperbolico né a provare una proprietà dissipativa discreta generale dello schema.

Nel caso iperbolico l'energia potenziale e quella cinetica vanno lette insieme al totale: la legge teorica riguarda quest'ultimo. La monotonia osservata delle curve non è stata usata da sola come criterio di accettazione.

## Discrepanze e limiti

- La velocità non è esportata direttamente nei risultati OOFEM. Per Newmark con \(\beta_N=1/4\), \(\gamma_N=1/2\) è stata ricostruita tramite \(v_{n+1}=2(u_{n+1}-u_n)/\Delta t-v_n\), coerente algebricamente con l'aggiornamento interno. Gli spostamenti testuali hanno precisione limitata, che viene amplificata dalla ricorrenza; ciò può contaminare soprattutto energia cinetica e residuo iperbolico.
- I tre livelli temporali non consentono di separare completamente errore di integrazione, arrotondamento dei campi e quadratura/post-processing. Le dissipazioni sono integrate con la ricostruzione temporale dichiarata e non costituiscono una misura esatta del residuo non lineare a ogni iterazione.
- Il calcolo usa una sola mesh, quindi non identifica errori spaziali, copertura errata del bordo o dipendenza dalla quadratura geometrica. Non sono state fatte prove spaziali in questa fase.
- I numeri provengono dai sei output preesistenti nella cartella `energy_validation_runs_highprecision`; lo studio non ha avviato nuove simulazioni. Il precedente test singolo iperbolico è documentato in [energy_hyperbolic_single_test.md](energy_hyperbolic_single_test.md).

Un eventuale seguito dovrebbe prima esportare la velocità Newmark nodale a precisione adeguata e ripetere il bilancio temporale; solo dopo, se il residuo residuo lo giustifica, pianificare separatamente una verifica spaziale a passo temporale sufficientemente piccolo.