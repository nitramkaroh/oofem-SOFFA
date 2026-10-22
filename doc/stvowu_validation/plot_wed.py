#!/usr/bin/env python3
from pathlib import Path
import csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).parent
rows=list(csv.DictReader((p/'wed_homogeneous.csv').open()))
fig,ax=plt.subplots(figsize=(5.6,3.8))
for t in (0.05,0.1,0.2):
 r=[x for x in rows if abs(float(x['time'])-t)<1e-12]
 ax.loglog([float(x['epsilon']) for x in r],[float(x['absolute_error']) for x in r],'o-',label=f't={t:g}')
ax.set_xlabel('epsilon'); ax.set_ylabel('errore assoluto vs soluzione analitica'); ax.grid(True,which='both',alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(p/'wed_epsilon_convergence.png',dpi=180)
