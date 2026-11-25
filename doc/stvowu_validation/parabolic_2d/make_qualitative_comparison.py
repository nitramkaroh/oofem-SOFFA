#!/usr/bin/env python3
"""Make a non-quantitative side-by-side plot for Xiao--Chen Fig. 5.19 at T=.1."""
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

HERE=Path(__file__).resolve().parent
REF=HERE/'reference_case4_fig519_T0p1.png'  # crop from the supplied original PDF, p.18
VTU=HERE/'runs/energy_n20_dt0p0001_every_step/case.out.m0.1000.vtu'
n=20
piece=ET.parse(VTU).getroot().find('.//Piece')
xyz=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)
values=np.fromstring(piece.find('./PointData/DataArray').text,sep=' ')
u=np.empty((n+1,n+1))
for (x,y,_),v in zip(xyz,values): u[int(round(y*n)),int(round(x*n))]=v
xy=np.linspace(0,1,n+1)
fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
axs[0].imshow(Image.open(REF));axs[0].axis('off');axs[0].set_title('Xiao–Chen (2025), Fig. 5.19, T=0.1')
im=axs[1].contourf(xy,xy,u,levels=np.linspace(-1.05,1.05,22),cmap='viridis',vmin=-1,vmax=1,extend='both')
axs[1].contour(xy,xy,u,levels=[0],colors='white',linewidths=1)
axs[1].set(xlim=(0,1),ylim=(0,1),xlabel='x',ylabel='y',title='OOFEM Q1, 20×20, Δt=10⁻⁴\nε=0.1 assumed (paper does not report ε)')
axs[1].set_aspect('equal');fig.colorbar(im,ax=axs[1],label='u')
fig.savefig(HERE/'xiao_chen_qualitative_T0p1.png',dpi=180)
