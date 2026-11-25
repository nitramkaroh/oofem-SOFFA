#!/usr/bin/env python3
"""Postprocess OOFEM reproduction of Okumura thesis Ch.4, Example 2."""
from pathlib import Path
import csv, re, xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

N=100; DT=.001; GAMMA=100.; OUT=Path('/tmp/stvowu_okumura_ch4_case2')
files=sorted(OUT.glob('okumura_ch4_case2.out.m0.*.vtu'),key=lambda p:int(re.search(r'm0\.(\d+)',p.name).group(1)))
gauss=(-1/np.sqrt(3),1/np.sqrt(3)); times=[]; energies=[]; fields=[]; y_spread=[]
for path in files:
    idx=int(re.search(r'm0\.(\d+)',path.name).group(1)); t=idx*DT
    piece=ET.parse(path).getroot().find('.//Piece')
    xyz=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)
    vals=np.fromstring(piece.find('./PointData/DataArray').text,sep=' ')
    # Coordinate-based assembly avoids assuming VTK point ordering.
    u=np.empty((2,N+1))
    for (x,y),v in zip(xyz[:,:2],vals):u[round(y)*(1),round(x*N)]=v
    y_spread.append(float(np.max(np.abs(u[0]-u[1]))))
    E=0.
    for i in range(N):
        q=np.array([u[0,i],u[0,i+1],u[1,i+1],u[1,i]])
        for a in gauss:
            for b in gauss:
                Ns=np.array([.25*(1-a)*(1-b),.25*(1+a)*(1-b),.25*(1+a)*(1+b),.25*(1-a)*(1+b)])
                dxi=np.array([-.25*(1-b),.25*(1-b),.25*(1+b),-.25*(1+b)])*2*N
                deta=np.array([-.25*(1-a),-.25*(1+a),.25*(1+a),.25*(1-a)])*2
                z=Ns@q
                E+=.25/N*(.5*((dxi@q)**2+(deta@q)**2)+GAMMA/4*(z*z-1)**2)
    # Two vertical dynamic ends; each has unit length and zero tangential gradient.
    for end in (u[:,0],u[:,-1]):
        for b in gauss:
            z=.5*(1-b)*end[0]+.5*(1+b)*end[1]
            E+=.5*(GAMMA/4*(z*z-1)**2)
    times.append(t); energies.append(E); fields.append(u[0].copy())

base=Path(__file__).parent
with (base/'okumura_ch4_example2_energy.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['time','OOFEM_Q1_total_energy']);w.writerows(zip(times,energies))
with (base/'okumura_ch4_example2_profiles.csv').open('w',newline='') as f:
    snapshots=(.001,.1,.2,.3,.4,.7)
    ids=[int(np.argmin(np.abs(np.asarray(times)-t))) for t in snapshots]
    w=csv.writer(f);w.writerow(['x']+[f'u_t_{times[k]:g}' for k in ids])
    for i in range(N+1):w.writerow([i/N]+[fields[k][i] for k in ids])

fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
im=axs[0].imshow(np.array(fields),origin='lower',aspect='auto',extent=[0,1,times[0],times[-1]],cmap='coolwarm',vmin=-1,vmax=1)
axs[0].set(xlabel='x',ylabel='t',title='OOFEM: u(x,t), Okumura Example 2');fig.colorbar(im,ax=axs[0],label='u')
axs[1].plot(times,energies,color='navy');axs[1].set(xlabel='t',ylabel='Total energy',title='Bulk + two dynamic endpoints');axs[1].grid(alpha=.25)
fig.savefig(base/'okumura_ch4_example2_oofem.png',dpi=180)
print('steps',len(times),'time range',times[0],times[-1])
print('energy first/final',energies[0],energies[-1],'max increment',max(np.diff(energies)))
print('max variation across extruded y',max(y_spread))
print('u(T=.7) range',float(fields[-1].min()),float(fields[-1].max()))
