import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).parent
with open(p/'spatial_convergence.csv') as f: spatial=list(csv.DictReader(f))
with open(p/'energy.csv') as f: energy=list(csv.DictReader(f))
fig,ax=plt.subplots(1,2,figsize=(9.2,3.6))
h=[2/int(r['nx']) for r in spatial]; e=[float(r['nodal_rms']) for r in spatial]
ax[0].loglog(h,e,'o-',label='RMS nodale')
ax[0].loglog(h,[e[0]*(x/h[0])**2 for x in h],'--',label='riferimento $h^2$')
ax[0].set_xlabel('passo spaziale h'); ax[0].set_ylabel('errore RMS'); ax[0].grid(True,which='both',alpha=.3); ax[0].legend()
t=[float(r['time']) for r in energy]; E=[float(r['discrete_total_energy']) for r in energy]
ax[1].plot(t,E,'-'); ax[1].set_xlabel('tempo'); ax[1].set_ylabel('energia discreta'); ax[1].grid(True,alpha=.3)
fig.tight_layout(); fig.savefig(p/'validation_results.png',dpi=180)
