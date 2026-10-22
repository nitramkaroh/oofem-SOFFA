#!/usr/bin/env python3
"""Publication-style plots for the targeted checklist robustness checks."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent
def read(name):
    with (P/name).open() as f:return list(csv.DictReader(f))

fig,ax=plt.subplots(1,3,figsize=(13,3.7))
base=read('wed_nonuniform.csv'); fixed=read('wed_epsilon_fixed_dt.csv')
ax[0].loglog([float(r['epsilon']) for r in base],[float(r['bulk_L2_space_time']) for r in base],'o-',label='bulk, ds=0.1')
ax[0].loglog([float(r['epsilon']) for r in base],[float(r['boundary_L2_space_time']) for r in base],'s--',label='bordo, ds=0.1')
ax[0].loglog([float(r['epsilon']) for r in fixed],[float(r['bulk_L2_space_time_vs_OOFEM']) for r in fixed],'o-',label='bulk, dt=0.005')
ax[0].loglog([float(r['epsilon']) for r in fixed],[float(r['boundary_L2_space_time_vs_OOFEM']) for r in fixed],'s--',label='bordo, dt=0.005')
ax[0].set_xlabel('epsilon');ax[0].set_ylabel('errore spazio-tempo');ax[0].set_title('Trend in epsilon');ax[0].grid(True,which='both',alpha=.3);ax[0].legend(fontsize=7)

h=read('wed_horizon_sensitivity.csv')
ax[1].plot([float(r['Smax']) for r in h],[float(r['space_time_RMS_vs_S18']) for r in h],'o-',label='RMS vs Smax=18')
ax[1].plot([float(r['Smax']) for r in h],[float(r['max_nodal_difference_vs_S18']) for r in h],'s--',label='massima differenza nodale')
ax[1].set_xlabel('Smax');ax[1].set_ylabel('differenza su t in [0, 0.1]');ax[1].set_title('Sensibilita all’orizzonte');ax[1].grid(True,alpha=.3);ax[1].legend(fontsize=7)

sp=read('wed_mesh_refinement.csv'); tight=read('wed_mesh_robustness.csv')
for rows,label,marker in ((sp,'baseline: ds=0.1, tol=2e-5','o'),(tight,'raffinato: ds=0.05, tol=2e-6','s')):
    fine_key='bulk_L2_space_time_vs_fine_WED_same_epsilon' if 'bulk_L2_space_time_vs_fine_WED_same_epsilon' in rows[0] else 'bulk_L2_space_time_vs_fine_WED'
    vals=[r for r in rows if float(r[fine_key])>0]
    ax[2].loglog([float(r['h']) for r in vals],[float(r[fine_key]) for r in vals],marker+'-',label=label)
ax[2].set_xlabel('h');ax[2].set_ylabel('errore bulk vs WED 32x16');ax[2].set_title('Robustezza mesh WED');ax[2].grid(True,which='both',alpha=.3);ax[2].legend(fontsize=7)
fig.tight_layout();fig.savefig(P/'wed_checklist_robustness.png',dpi=200)
