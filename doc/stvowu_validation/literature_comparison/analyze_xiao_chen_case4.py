#!/usr/bin/env python3
"""Extract the exploratory OOFEM Case 4 field and its half-level contour."""
from pathlib import Path
import argparse, csv, xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,default=20);parser.add_argument('--dt',type=float,default=.0001);args=parser.parse_args();n=args.n;dt=args.dt
tag=f'{dt:g}'
src=Path(f'/tmp/stvowu_xiao_chen_case4_n{n}_dt{tag}_eps0.1/xiao_chen_case4.out.m0.{round(.1/dt)}.vtu')
if n==20 and dt==.001 and not src.exists():src=Path('/tmp/stvowu_xiao_chen_case4/xiao_chen_case4.out.m0.100.vtu')
if not src.exists() and n==20 and dt==.0001:
    src=Path(__file__).parent/'xiao_chen_case4_n20_dt0p0001_eps0p1_final.vtu'
root=ET.parse(src).getroot(); piece=root.find('.//Piece')
npts=int(piece.get('NumberOfPoints')); assert npts==(n+1)**2
xyz=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)
values=np.fromstring(piece.find('./PointData/DataArray').text,sep=' ')
x=np.unique(np.round(xyz[:,0],12)); y=np.unique(np.round(xyz[:,1],12))
u=np.empty((n+1,n+1))
for (xx,yy),value in zip(xyz[:,:2],values):
    u[int(round(yy*n)),int(round(xx*n))]=value
fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
cf=ax.contourf(x,y,u,levels=np.linspace(-1,1,21),cmap='viridis',vmin=-1,vmax=1,extend='both')
ax.contour(x,y,u,levels=[0],colors='white',linewidths=1.5)
ax.set(xlim=(0,1),ylim=(0,1),xlabel='x',ylabel='y',title=f'OOFEM Xiao–Chen Case 4, T=0.1\nexploratory ε=0.1, dt={dt:g}; Q1 mesh {n}×{n}')
ax.set_aspect('equal'); fig.colorbar(cf,ax=ax,label='u')
out=Path(__file__).parent/f'xiao_chen_case4_exploratory_n{n}_dt{tag}.png'; fig.savefig(out,dpi=180)
# Contour at u=0.5 approximated by linear interpolation on vertical/horizontal grid lines.
rows=[]
for j in range(n+1):
    ix=np.where((u[j,:-1]-.5)*(u[j,1:]-.5)<=0)[0]
    for i in ix:
        a=(.5-u[j,i])/(u[j,i+1]-u[j,i]) if u[j,i+1]!=u[j,i] else 0.
        rows.append({'direction':'horizontal','fixed_coord':y[j],'crossing':x[i]+a*(x[i+1]-x[i])})
for i in range(n+1):
    iy=np.where((u[:-1,i]-.5)*(u[1:,i]-.5)<=0)[0]
    if iy.size:
        j=iy[-1]; a=(.5-u[j,i])/(u[j+1,i]-u[j,i]) if u[j+1,i]!=u[j,i] else 0.
        rows.append({'direction':'vertical','fixed_coord':x[i],'crossing':y[j]+a*(y[j+1]-y[j])})
with (Path(__file__).parent/f'xiao_chen_case4_contour_n{n}_dt{tag}.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['direction','fixed_coord','crossing']);w.writeheader();w.writerows(rows)
print('range',float(u.min()),float(u.max()))
print('saved',out)
print('horizontal u=0.5 crossings:',[(r['fixed_coord'],r['crossing']) for r in rows if r['direction']=='horizontal'])
print('vertical u=0.5 crossings:',[(r['fixed_coord'],r['crossing']) for r in rows if r['direction']=='vertical'])
