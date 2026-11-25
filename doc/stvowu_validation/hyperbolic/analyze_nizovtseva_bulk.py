#!/usr/bin/env python3
"""Recompute Nizovtseva benchmark metrics and figures from retained VTU histories."""
from pathlib import Path
import importlib.util, math
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('nizovalid',HERE/'validate_nizovtseva_bulk.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def last(tag):
    frames=m.load_history(HERE/f'nizovtseva_bulk_v2_{tag}')
    return frames[-1]

def field(frame):
    _,pts,u,_=frame; x=pts[:,0]; xu=np.unique(x)
    return xu,np.array([np.mean(u[np.isclose(x,z)]) for z in xu])

sp=[]
for n in (80,160,320):
    tag=f'space_n{n}'; fr=last(tag); x,u=field(fr)
    row=m.metrics(x,u,float(x.min()),float(x.max()),tag,n,.001)
    sp.append(row)
for a,b in zip(sp,sp[1:]):
    b['order_full']=math.log(a['rms_full']/b['rms_full'],2)
    b['order_inner_absx_le6']=math.log(a['rms_inner_absx_le6']/b['rms_inner_absx_le6'],2)
m.write_csv(HERE/'nizovtseva_bulk_v2_spatial.csv',sp)

ti=[]; fields=[]
for dt in (.04,.02,.01,.005):
    tag=f'time_n320_dt{dt:g}'; fr=last(tag); x,u=field(fr)
    ti.append(m.metrics(x,u,-12.,12.,tag,320,dt)); fields.append((x,u))
ti.sort(key=lambda r:r['dt'],reverse=True)
for i in range(len(ti)-1):
    diff=fields[i][1]-fields[i+1][1]; x=fields[i][0]; inner=np.abs(x)<=6
    ti[i]['rms_to_next_finer']=float(np.sqrt(np.mean(diff**2)))
    ti[i]['rms_inner_to_next_finer']=float(np.sqrt(np.mean(diff[inner]**2)))
for i in range(len(ti)-2):
    a,b=ti[i]['rms_to_next_finer'],ti[i+1]['rms_to_next_finer']
    ti[i]['order_from_adjacent_differences']=math.log(a/b,2)
m.write_csv(HERE/'nizovtseva_bulk_v2_temporal.csv',ti)

short=last('domain_m8_p8'); long=last('domain_m12_p12')
xs,us=field(short); xl,ul=field(long)
common=(xl>=-8-1e-10)&(xl<=8+1e-10)
assert np.allclose(xs,xl[common])
err_s=us-m.exact(xs,m.T); err_l=ul[common]-m.exact(xl[common],m.T)
delta=us-ul[common]; inner=np.abs(xs)<=6
boundary=[
 {'tag':'truncated_domain_minus8_plus8','xmin':-8,'xmax':8,'h':.1,
  'rms_exact_common':float(np.sqrt(np.mean(err_s**2))),'max_exact_common':float(np.max(np.abs(err_s))),
  'rms_exact_inner_absx_le6':float(np.sqrt(np.mean(err_s[inner]**2))),'max_exact_inner_absx_le6':float(np.max(np.abs(err_s[inner]))),
  'rms_vs_long_common':float(np.sqrt(np.mean(delta**2))),'max_vs_long_common':float(np.max(np.abs(delta)))},
 {'tag':'reference_domain_minus12_plus12','xmin':-12,'xmax':12,'h':.1,
  'rms_exact_common':float(np.sqrt(np.mean(err_l**2))),'max_exact_common':float(np.max(np.abs(err_l))),
  'rms_exact_inner_absx_le6':float(np.sqrt(np.mean(err_l[inner]**2))),'max_exact_inner_absx_le6':float(np.max(np.abs(err_l[inner]))),
  'rms_vs_long_common':0.,'max_vs_long_common':0.}]
m.write_csv(HERE/'nizovtseva_bulk_v2_domain.csv',boundary)

energy=np.genfromtxt(HERE/'nizovtseva_bulk_v2_energy.csv',delimiter=',',names=True)
dt=np.diff(energy['time']); incr=np.diff(energy['energy'])
monotone=bool(np.all(incr<=1e-10))

fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
# Profile
x,u=field(last('space_n320')); ax[0,0].plot(x,u,label='OOFEM Q1, N=320, dt=.001')
ax[0,0].plot(x,m.exact(x,m.T),'--',label='Nizovtseva Eq. (29)')
ax[0,0].set(xlim=(-1,1),xlabel='x',ylabel='u(x,T)',title='Profilo a T=0.2'); ax[0,0].legend()
# Spatial
hs=np.array([r['h'] for r in sp]); er=np.array([r['rms_full'] for r in sp]); ei=np.array([r['rms_inner_absx_le6'] for r in sp])
ax[0,1].loglog(hs,er,'o-',label='RMS intero dominio'); ax[0,1].loglog(hs,ei,'s-',label='RMS |x|≤6')
ax[0,1].loglog(hs,er[-1]*(hs/hs[-1])**2,'k:',label='pendenza h²')
ax[0,1].set(xlabel='h',ylabel='RMS verso esatta',title='Raffinamento spaziale'); ax[0,1].legend()
# Temporal pair differences
pairs=[(ti[i]['dt'],ti[i]['rms_to_next_finer']) for i in range(len(ti)-1)]
ax[1,0].loglog([p[0] for p in pairs],[p[1] for p in pairs],'o-',label='RMS tra passi adiacenti')
ax[1,0].loglog([p[0] for p in pairs[:2]],[pairs[0][1]*(p[0]/pairs[0][0])**2 for p in pairs[:2]],'k:',label='pendenza dt²')
ax[1,0].set(xlabel='dt del passo più grosso',ylabel='RMS differenza campi',title='Raffinamento temporale (N=320)'); ax[1,0].legend()
# Energy and residual
cum=energy['energy_balance_residual']-energy['energy']+energy['energy'][0]
ax[1,1].plot(energy['time'],energy['energy']-energy['energy'][0],label='E(t)−E(0)')
ax[1,1].plot(energy['time'],-cum,'--',label='−∫₀ᵗ η‖uₜ‖²')
ax[1,1].set(xlabel='t',ylabel='variazione di energia',title=f'Bilancio, E monotona: {monotone}'); ax[1,1].legend(loc='best')
fig.savefig(HERE/'nizovtseva_bulk_validation.png',dpi=180)
# Boundary truncation and initial-velocity sensitivity, in physical coordinates.
fig2,ax2=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
for tag,label in [('domain_m8_p8','dominio [-8,8]'),('domain_m12_p12','dominio [-12,12]')]:
    xx,uu=field(last(tag)); ax2[0,0].plot(xx,uu,label=label)
    ax2[0,1].semilogy(xx,np.maximum(np.abs(uu-m.exact(xx,m.T)),1.e-10),label=label)
xx=np.linspace(-8,8,801); ax2[0,0].plot(xx,m.exact(xx,m.T),'k--',label='esatta')
ax2[0,0].set(xlim=(-8,8),xlabel='x',ylabel='u(x,T)',title='Profilo con bordo troncato'); ax2[0,0].legend()
ax2[0,1].set(xlim=(-8,8),xlabel='x',ylabel='|u−u esatta|',title='Errore in funzione di x'); ax2[0,1].legend()
base_x,base_u=field(last('space_n160'))
pert_x,pert_u=field(last('initial_velocity_plus1pct'))
ax2[1,0].plot(base_x,base_u,label='velocità iniziale esatta')
ax2[1,0].plot(pert_x,pert_u,'--',label='velocità iniziale +1%')
ax2[1,0].set(xlim=(-1,1),xlabel='x',ylabel='u(x,T)',title='Profilo dopo perturbazione iniziale'); ax2[1,0].legend()
ax2[1,1].plot(base_x,pert_u-base_u)
ax2[1,1].set(xlim=(-1,1),xlabel='x',ylabel='u perturbato−u base',title='Variazione del campo per +1% in velocità')
fig2.savefig(HERE/'nizovtseva_bulk_boundary_initial.png',dpi=180)
print('spatial=',sp)
print('temporal=',ti)
print('domain=',boundary)
print('energy monotone=',monotone,'max energy increase=',float(np.max(incr)),'max relative residual=',float(np.max(np.abs(energy['energy_balance_residual']))/energy['energy'][0]))
