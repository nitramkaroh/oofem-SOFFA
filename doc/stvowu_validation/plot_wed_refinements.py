#!/usr/bin/env python3
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent
def read(name):
 with (P/name).open() as f:return list(csv.DictReader(f))
fig,ax=plt.subplots(1,2,figsize=(9,3.7))
t=read('wed_nonuniform_time_refinement.csv')
ax[0].loglog([float(r['physical_dt']) for r in t],[float(r['bulk_L2_space_time']) for r in t],'o-',label='bulk')
ax[0].loglog([float(r['physical_dt']) for r in t],[float(r['boundary_L2_space_time']) for r in t],'s-',label='bordo')
ax[0].set_xlabel('passo fisico WED');ax[0].set_ylabel('errore L2 spazio-tempo');ax[0].set_title('Raffinamento temporale, epsilon=0.05');ax[0].legend();ax[0].grid(True,which='both',alpha=.3)
ax[0].set_xticks([0.00125,0.0025],['0.00125','0.0025'])
h=read('wed_mesh_refinement.csv')
h=[r for r in h if float(r['bulk_L2_space_time_vs_fine_WED_same_epsilon'])>0]
ax[1].loglog([float(r['h']) for r in h],[float(r['bulk_L2_space_time_vs_fine_WED_same_epsilon']) for r in h],'o-',label='bulk')
ax[1].loglog([float(r['h']) for r in h],[float(r['boundary_L2_space_time_vs_fine_WED_same_epsilon']) for r in h],'s-',label='bordo')
ax[1].set_xlabel('h');ax[1].set_ylabel('errore L2 spazio-tempo vs WED fine');ax[1].set_title('Convergenza WED (rif. 32x16)');ax[1].legend();ax[1].grid(True,which='both',alpha=.3)
ax[1].set_xticks([0.125,0.25],['0.125','0.25'])
fig.tight_layout();fig.savefig(P/'wed_refinements.png',dpi=190)
