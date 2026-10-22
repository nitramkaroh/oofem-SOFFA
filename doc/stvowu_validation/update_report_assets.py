#!/usr/bin/env python3
"""Regenerate energy history, beta summary, and requested report figures."""
from pathlib import Path
import csv, importlib.util, math, xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('energy',HERE/'run_energy_balance.py')
energy=importlib.util.module_from_spec(spec); spec.loader.exec_module(energy)

def values(path):
    root=ET.parse(path).getroot(); return np.fromstring(root.find('.//PointData/DataArray').text,sep=' ')

def save_energy_history():
    case=Path('/tmp/stvowu_energy_balance/dt_0p00125/case.in')
    nodes,elems,u0=energy.geometry(case); folder=case.parent
    rows=[]; eb,eg=energy.energy(u0,nodes,elems); rows.append((0.,eb,eg,eb+eg))
    for step in range(1,81):
        u=values(folder/f'energy_balance_0p00125.out.m0.{step}.vtu')
        eb,eg=energy.energy(u,nodes,elems); rows.append((step*.00125,eb,eg,eb+eg))
    with (HERE/'energy_trajectory.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['time','bulk_energy','boundary_energy','total_energy']); w.writerows(rows)
    fig,ax=plt.subplots(figsize=(6.7,4.1)); t=np.array(rows)
    ax.plot(t[:,0],t[:,3],lw=2,label='totale'); ax.plot(t[:,0],t[:,1],label='bulk'); ax.plot(t[:,0],t[:,2],label='bordo')
    ax.set(xlabel='tempo',ylabel='energia',title='Dissipazione energetica, $\\Delta t=0.00125$'); ax.grid(True,alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(HERE/'energy_trajectory.png',dpi=190); plt.close(fig)

def save_beta_table():
    from run_surface_beta_study import read_mesh, surface_stats
    nodes,edges=read_mesh(); base=Path('/tmp/stvowu_beta_study'); dt=.0005; rows=[]
    for beta in (.5,1.,2.,5.):
        label=str(beta).replace('.','p'); folder=base/f'beta_{label}'; case=folder/'case.in'
        ns,elems,u0=energy.geometry(case); mb,mg=energy.assembled(ns,elems)
        prev=u0; db=dg=0.
        for step in range(1,41):
            u=values(folder/f'beta_{label}.out.m0.{step}.vtu'); du=u-prev
            db+=float(du@mb@du)/dt; dg+=float(du@mg@du)/dt; prev=u
        mean,e=surface_stats(prev,nodes,edges,beta)
        rows.append((beta,mean,e,db,dg,db+dg))
    with (HERE/'surface_beta_summary.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['beta','boundary_mean_T','boundary_energy_T','bulk_dissipation','boundary_dissipation','total_dissipation']); w.writerows(rows)

def save_snapshots():
    case=Path('/tmp/stvowu_wed_nonuniform/direct.in'); nodes,elems,u0=energy.geometry(case)
    ids=sorted(nodes); xy=np.array([nodes[i] for i in ids]); idmap={n:i for i,n in enumerate(ids)}
    quads=[conn for typ,conn in elems if typ=='quad1ht']; tris=[]
    for conn in quads:
        a,b,c,d=[idmap[n] for n in conn]; tris.extend(((a,b,c),(a,c,d)))
    tri=mtri.Triangulation(xy[:,0],xy[:,1],np.array(tris))
    out=Path('/tmp/stvowu_wed_nonuniform'); times=(0.,.02,.1,.2); fields=[u0]
    for t in times[1:]: fields.append(values(out/f'direct_nonuniform.out.m0.{round(t/.001)}.vtu'))
    levels=np.linspace(min(float(u.min()) for u in fields),max(float(u.max()) for u in fields),21)
    trace_data=[boundary_trace(xy,u) for u in fields]
    all_trace=np.concatenate([v for _,v in trace_data]); lo=float(all_trace.min()); hi=float(all_trace.max()); pad=.04*(hi-lo)
    fig,axs=plt.subplots(2,4,figsize=(11,5.2))
    for j,(t,u) in enumerate(zip(times,fields)):
        cf=axs[0,j].tricontourf(tri,u,levels=levels,cmap='coolwarm'); axs[0,j].set_title(f'$t={t:g}$'); axs[0,j].set_aspect('equal'); axs[0,j].set_xlabel('x')
        if j==0: axs[0,j].set_ylabel('y')
        axs[1,j].plot(*trace_data[j],color='black'); axs[1,j].set_ylim(lo-pad,hi+pad); axs[1,j].set_xlim(trace_data[0][0][0],trace_data[0][0][-1]); axs[1,j].set_title(f'$u_\\Gamma$, $t={t:g}$'); axs[1,j].grid(True,alpha=.25); axs[1,j].set_xlabel('ascissa curvilinea')
        if j==0: axs[1,j].set_ylabel('$u_\\Gamma$')
    fig.subplots_adjust(left=.07,right=.90,bottom=.10,top=.86,wspace=.28,hspace=.42)
    cax=fig.add_axes([.925,.56,.018,.27]); fig.colorbar(cf,cax=cax,label='$u(x,t)$')
    fig.suptitle('Campo nel bulk e traccia sul bordo'); fig.savefig(HERE/'field_and_boundary_snapshots.png',dpi=190); plt.close(fig)

def save_mesh_plot():
    with (HERE/'spatial_convergence.csv').open() as f: rows=list(csv.DictReader(f))
    h=np.array([2./int(r['nx']) for r in rows]); err=np.array([float(r['nodal_rms']) for r in rows])
    fig,ax=plt.subplots(figsize=(5.4,4.0)); ax.loglog(h,err,'o-',label='RMS nodale MMS')
    ax.set(xlabel='$h$',ylabel='errore RMS',title='Convergenza spaziale FEM'); ax.grid(True,which='both',alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(HERE/'spatial_convergence.png',dpi=190); plt.close(fig)

def boundary_trace(xy,u):
    xmax=float(xy[:,0].max()); ymax=float(xy[:,1].max()); tol=1.e-10; pts=[]
    for i,(x,y) in enumerate(xy):
        if abs(y)<tol: s=x
        elif abs(x-xmax)<tol: s=xmax+y
        elif abs(y-ymax)<tol: s=xmax+ymax+(xmax-x)
        elif abs(x)<tol: s=2*xmax+ymax+(ymax-y)
        else: continue
        pts.append((s,float(u[i])))
    pts.sort(); out=[]
    for p in pts:
        if not out or abs(p[0]-out[-1][0])>tol: out.append(p)
    a=np.array(out); return a[:,0],a[:,1]

if __name__=='__main__':
    save_energy_history(); save_beta_table(); save_snapshots(); save_mesh_plot()
