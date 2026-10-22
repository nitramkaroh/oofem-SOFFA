#!/usr/bin/env python3
"""Build the six self-contained figures requested for the paper section."""
from pathlib import Path
import csv, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent
OUT=P/'paper_figures';OUT.mkdir(exist_ok=True)
def read(name):
    with (P/name).open() as f:return list(csv.DictReader(f))

# Figure 1: bulk snapshots and boundary trace.
shutil.copy2(P/'field_and_boundary_snapshots.png',OUT/'figure1_field_boundary.png')

# Figure 2: epsilon convergence with separate fixed-ds and fixed-physical-dt studies.
a=read('wed_nonuniform.csv');b=read('wed_epsilon_fixed_dt.csv')
fig,axs=plt.subplots(1,2,figsize=(9.2,3.7))
for ax,key,fixedkey,title in ((axs[0],'bulk_L2_space_time','bulk_L2_space_time_vs_OOFEM','Bulk'),(axs[1],'boundary_L2_space_time','boundary_L2_space_time_vs_OOFEM','Bordo')):
    ax.loglog([float(r['epsilon']) for r in a],[float(r[key]) for r in a],'o-',label='passo riscalato fisso: ds=0.1')
    ax.loglog([float(r['epsilon']) for r in b],[float(r[fixedkey]) for r in b],'s--',label='passo fisico fisso: dt=0.005')
    ax.set_xlabel('epsilon');ax.set_ylabel('errore L2 spazio-tempo');ax.set_title(title);ax.grid(True,which='both',alpha=.3);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'figure2_epsilon.png',dpi=210);plt.close(fig)

# Figure 3: trunction sensitivity; errors are relative to Smax=18.
h=read('wed_horizon_sensitivity.csv'); h=h[:2]
fig,axs=plt.subplots(1,2,figsize=(8.5,3.5))
for ax,key,title in ((axs[0],'space_time_RMS_vs_S18','RMS spazio-tempo'),(axs[1],'max_nodal_difference_vs_S18','Massimo scarto nodale')):
    ax.bar([str(int(float(r['Smax']))) for r in h],[float(r[key]) for r in h],color=['#4477AA','#66CCEE'])
    ax.set_xlabel('Smax (riferimento: 18)');ax.set_ylabel('scarto su t in [0, 0.1]');ax.set_title(title);ax.grid(axis='y',alpha=.3)
    for i,r in enumerate(h):ax.text(i,float(r[key]),f"{float(r[key]):.2e}",ha='center',va='bottom',fontsize=8)
fig.suptitle('Sensibilità al troncamento WED: epsilon=0.05, ds=0.1, tol=2e-5')
fig.tight_layout();fig.savefig(OUT/'figure3_horizon.png',dpi=210);plt.close(fig)

# Figure 4: both bulk and boundary mesh trends, comparing baseline and tighter settings.
base=read('wed_mesh_refinement.csv');tight=read('wed_mesh_robustness.csv')
fig,axs=plt.subplots(1,2,figsize=(9.2,3.7))
for ax,bkey,tkey,title in ((axs[0],'bulk_L2_space_time_vs_fine_WED_same_epsilon','bulk_L2_space_time_vs_fine_WED','Bulk'),(axs[1],'boundary_L2_space_time_vs_fine_WED_same_epsilon','boundary_L2_space_time_vs_fine_WED','Bordo')):
    for rows,key,label,mark in ((base,bkey,'ds=0.1, tol=2e-5','o'),(tight,tkey,'ds=0.05, tol=2e-6','s')):
        vals=[r for r in rows if float(r[key])>0]
        ax.loglog([float(r['h']) for r in vals],[float(r[key]) for r in vals],mark+'-',label=label)
    ax.set_xlabel('h');ax.set_ylabel('errore spazio-tempo vs WED 32x16');ax.set_title(title);ax.grid(True,which='both',alpha=.3);ax.legend(fontsize=8)
fig.suptitle('Convergenza spaziale WED, epsilon=0.05; riferimento numerico fine')
fig.tight_layout();fig.savefig(OUT/'figure4_mesh.png',dpi=210);plt.close(fig)

# Figure 5: total and component energies from corrected postprocessing.
shutil.copy2(P/'energy_trajectory.png',OUT/'figure5_energy.png')

# Figure 6: direct and WED beta series shown separately with their distinct final times.
direct=read('surface_beta.csv');wed=read('wed_surface_beta.csv')
fig,axs=plt.subplots(1,2,figsize=(9.2,3.7))
for beta in sorted({float(r['beta']) for r in direct}):
    q=[r for r in direct if float(r['beta'])==beta]
    axs[0].plot([float(r['time']) for r in q],[float(r['boundary_mean']) for r in q],label=f'beta={beta:g}')
for beta in sorted({float(r['beta']) for r in wed}):
    q=[r for r in wed if float(r['beta'])==beta]
    axs[1].plot([float(r['physical_time']) for r in q],[float(r['boundary_mean']) for r in q],label=f'beta={beta:g}')
for ax,title,xlabel in ((axs[0],'OOFEM diretto, T=0.02','tempo'),(axs[1],'minimizzatore WED, T=0.2','tempo fisico')):
    ax.set_title(title);ax.set_xlabel(xlabel);ax.set_ylabel('media di u sul bordo');ax.grid(True,alpha=.3);ax.legend(fontsize=8)
fig.suptitle('Risposta al parametro beta; serie di backend e durate non confrontabili punto a punto')
fig.tight_layout();fig.savefig(OUT/'figure6_beta.png',dpi=210);plt.close(fig)
