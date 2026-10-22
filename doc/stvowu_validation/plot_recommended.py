#!/usr/bin/env python3
"""Plot epsilon, surface-beta and energy-balance studies for the report."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent

def rows(name):
    with (P/name).open() as f: return list(csv.DictReader(f))
fig,ax=plt.subplots(2,2,figsize=(10,7.2))
wed=rows('wed_nonuniform.csv')
eps=[float(r['epsilon']) for r in wed]
ax[0,0].loglog(eps,[float(r['bulk_L2_space_time']) for r in wed],'o-',label='bulk')
ax[0,0].loglog(eps,[float(r['boundary_L2_space_time']) for r in wed],'s-',label='bordo')
ax[0,0].set_xlabel(r'$\varepsilon$'); ax[0,0].set_ylabel(r'errore $L^2$ spazio-tempo'); ax[0,0].set_title('WED non uniforme vs OOFEM, t <= 0.1'); ax[0,0].legend(); ax[0,0].grid(True,which='both',alpha=.3)
beta=rows('surface_beta.csv')
for b in sorted({float(r['beta']) for r in beta}):
 q=[r for r in beta if float(r['beta'])==b]; ax[0,1].plot([float(r['time']) for r in q],[float(r['boundary_mean']) for r in q],label=fr'$\beta={b:g}$')
ax[0,1].set_xlabel('tempo'); ax[0,1].set_ylabel(r'$M_\Gamma(t)$'); ax[0,1].set_title('OOFEM diretto: effetto di beta'); ax[0,1].legend(); ax[0,1].grid(True,alpha=.3)
wb=rows('wed_surface_beta.csv')
for b in sorted({float(r['beta']) for r in wb}):
 q=[r for r in wb if float(r['beta'])==b]; ax[1,0].plot([float(r['physical_time']) for r in q],[float(r['boundary_mean']) for r in q],label=fr'$\beta={b:g}$')
ax[1,0].set_xlabel('tempo fisico'); ax[1,0].set_ylabel(r'$M_\Gamma(t)$'); ax[1,0].set_title('WED (epsilon=0.05)'); ax[1,0].legend(); ax[1,0].grid(True,alpha=.3)
eb=rows('energy_balance.csv')
ax[1,1].loglog([float(r['dt']) for r in eb],[float(r['normalized_abs_residual']) for r in eb],'o-')
ax[1,1].set_xlabel(r'$\Delta t$'); ax[1,1].set_ylabel(r'$|R_E|/E(0)$'); ax[1,1].set_title('Residuo bilancio energia non uniforme'); ax[1,1].grid(True,which='both',alpha=.3)
fig.tight_layout(); fig.savefig(P/'recommended_tests.png',dpi=190)
