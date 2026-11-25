#!/usr/bin/env python3
"""Run non-forced parabolic and hyperbolic Allen--Cahn energy checks.

Q1 bulk / P1 boundary FEM matrices and energies are independently assembled
from the structured mesh description. The OOFEM histories are retained in
energy_validation_runs/ for reproducibility.
"""
from __future__ import annotations
import csv, math, os, re, subprocess, time, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/oofem_energy_mpl')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OOFEM=ROOT/'build/debug_magneto/oofem'
OUT=Path(os.environ.get('OOFEM_ENERGY_OUT', HERE/'energy_validation_runs'))
OUTPUT_PREFIX=os.environ.get('OOFEM_ENERGY_PREFIX', 'energy_validation')
LX,LY,T=2.,1.,.02
K,M,C,RHO,CS,BETA=1.,1.,1.,1.,1.,1.
DELTA=.2                 # physical bulk damping coefficient in rho*delta*u_t
GLX,GLW=np.polynomial.legendre.leggauss(5)

def mesh(nx,ny):
    xy=np.array([(i*LX/nx,j*LY/ny) for j in range(ny+1) for i in range(nx+1)])
    nd=lambda i,j:j*(nx+1)+i
    quads=np.array([[nd(i,j),nd(i+1,j),nd(i+1,j+1),nd(i,j+1)] for j in range(ny) for i in range(nx)])
    edges=[]
    edges += [[nd(i,0),nd(i+1,0)] for i in range(nx)]
    edges += [[nd(nx,j),nd(nx,j+1)] for j in range(ny)]
    edges += [[nd(i,ny),nd(i-1,ny)] for i in range(nx,0,-1)]
    edges += [[nd(0,j),nd(0,j-1)] for j in range(ny,0,-1)]
    return xy,quads,np.asarray(edges,dtype=int),np.unique(np.asarray(edges).ravel())

def qshape(rr,ss):
    N=.25*np.asarray([(1-rr)*(1-ss),(1+rr)*(1-ss),(1+rr)*(1+ss),(1-rr)*(1+ss)],dtype=float)
    dN=.25*np.asarray([[-(1-ss),-(1-rr)],[(1-ss),-(1+rr)],[(1+ss),(1+rr)],[-(1+ss),(1-rr)]],dtype=float)
    return N,dN

def fields(xy):
    x,y=xy.T
    return .2+.1*np.cos(np.pi*x/LX)*np.cos(np.pi*y/LY)

def matrices(xy,quads,edges):
    n=len(xy); Mb=np.zeros((n,n)); Ms=np.zeros((n,n))
    # Three point tensor Gauss integrates all entries and the cubic source.
    z,w=np.polynomial.legendre.leggauss(3)
    for ids in quads:
        X=xy[ids]
        for i,r in enumerate(z):
            for j,s in enumerate(z):
                N,dN=qshape(r,s); J=X.T@dN; det=np.linalg.det(J)
                Mb[np.ix_(ids,ids)] += np.outer(N,N)*det*w[i]*w[j]
    for ids in edges:
        length=np.linalg.norm(xy[ids[1]]-xy[ids[0]])
        for r,ww in zip(z,w):
            N=np.array([(1-r)/2,(1+r)/2])
            Ms[np.ix_(ids,ids)] += np.outer(N,N)*length*ww/2
    return Mb,Ms

def energy_parts(u,v,xy,quads,edges):
    # Five-point Gauss integrates Q1 gradient and quartic potential exactly.
    ebg=ebp=ek=esg=esp=0.
    for ids in quads:
        X=xy[ids]; U=u[ids]; V=v[ids]
        for i,r in enumerate(GLX):
            for j,s in enumerate(GLX):
                N,dN=qshape(r,s); J=X.T@dN; det=np.linalg.det(J)
                grad=dN@np.linalg.inv(J.T); uh=N@U; vh=N@V; gu=U@grad
                w=det*GLW[i]*GLW[j]
                ebg += .5*K*(gu@gu)*w
                ebp += M*.25*(uh*uh-1.)**2*w
                ek += .5*RHO*vh*vh*w
    z,w=np.polynomial.legendre.leggauss(5)
    for ids in edges:
        X=xy[ids]; length=np.linalg.norm(X[1]-X[0]); U=u[ids]
        dN=np.array([-1.,1.])/length
        for r,ww in zip(z,w):
            N=np.array([(1-r)/2,(1+r)/2]); uh=N@U; tangent=dN@U
            esg += .5*K*tangent*tangent*length*ww/2
            esp += BETA*M*.25*(uh*uh-1.)**2*length*ww/2
    return ebg,ebp,ek,esg,esp

def input_text(nx,ny,dt,kind,tag):
    xy,quads,edges,boundary=mesh(nx,ny); nn=len(xy); nsteps=round(T/dt); ne=len(quads)+len(edges)
    wave=kind=='hyperbolic'
    header=(f'nltransienttransportproblem nsteps {nsteps} deltat {dt:.16g} alpha 1.0 rtol 1.e-8 nsmax 80 nmodules 1'
            + (f' wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0 wavedamping {DELTA}' if wave else ''))
    lines=[f'{tag}.out',f'Non-forced {kind} Allen-Cahn energy check',header,
      'vtkxml tstep_all domain_all primvars 1 6','domain heattransfer','OutputManager tstep_all dofman_all element_all',
      f'ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}']
    for i,(x,y) in enumerate(xy,1): lines.append(f'node {i} coords 3 {x:.16g} {y:.16g} 0.0')
    for e,ids in enumerate(quads,1): lines.append(f'quad1ht {e} nodes 4 '+' '.join(str(i+1) for i in ids)+' nip 9')
    for e,ids in enumerate(edges,len(quads)+1): lines.append(f'line1ht {e} nodes 2 '+' '.join(str(i+1) for i in ids)+' nip 3')
    lines += ['SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
       f'allencahn 1 d 1.0 k {K} c {C} mobility {M} surfacebeta {BETA} manufactured 0'+
       (f' hyperbolic 1 rho {RHO} surfacecapacity {CS}' if wave else '')]
    u0=fields(xy)
    for i,u in enumerate(u0,1): lines.append(f'InitialCondition {i} dofs 1 10 Conditions 1 u {u:.17g} set {i+1}')
    lines += ['ConstantFunction 1 f(t) 1.0',f'Set 1 elementranges {{(1 {ne})}}']
    lines += [f'Set {i+1} nodes 1 {i}' for i in range(1,nn+1)]
    return '\n'.join(lines)+'\n',xy,quads,edges,boundary

def read_vtu(file):
    root=ET.parse(file).getroot(); piece=root.find('.//Piece')
    xyz=np.fromstring(piece.find('./Points/DataArray').text,sep=' ').reshape(-1,3)[:,:2]
    arrays=piece.findall('./PointData/DataArray')
    data={a.attrib.get('Name',''):np.fromstring(a.text or '',sep=' ') for a in arrays}
    if 'Temperature' not in data: raise RuntimeError(f'No Temperature array in {file}')
    return xyz,data['Temperature']

def read_text_history(file, expected_steps, expected_nodes, dt):
    """Read full precision nodal values from OOFEM's text output.

    VTU ASCII output is rounded more aggressively. Differentiating its fields
    to recover Newmark velocity can magnify that rounding when dt is small.
    """
    content=Path(file).read_text()
    blocks=re.split(r'Output for time\s+([\deE.+-]+)',content)
    history=[]
    for i in range(1,len(blocks),2):
        time_value=float(blocks[i]); body=blocks[i+1]
        values=re.findall(r'dof\s+10\s+d\s+([\deE.+-]+)',body)
        if not values:
            continue
        if len(values)!=expected_nodes:
            raise RuntimeError(f'Expected {expected_nodes} nodal values at t={time_value:g}; got {len(values)}')
        history.append((time_value,np.asarray(values,dtype=float)))
    if len(history)!=expected_steps:
        raise RuntimeError(f'Expected {expected_steps} text-output steps; got {len(history)}')
    for i,(time_value,_) in enumerate(history,1):
        if abs(time_value-i*dt)>1.e-10:
            raise RuntimeError(f'Text-output time mismatch at step {i}: {time_value:g}')
    return history

def run_case(kind,nx,ny,dt):
    tag=f'{kind}_n{nx}x{ny}_dt{dt:g}'.replace('.','p')
    folder=OUT/tag; folder.mkdir(parents=True,exist_ok=True)
    inp,xy,quads,edges,boundary=input_text(nx,ny,dt,kind,tag)
    (folder/'case.in').write_text(inp)
    tic=time.perf_counter(); p=subprocess.run([str(OOFEM),'-f',str(folder/'case.in')],cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (folder/'oofem.log').write_text(p.stdout)
    if p.returncode or not (folder/f'{tag}.out.m0.pvd').exists():
        raise RuntimeError(f'{kind} {nx}x{ny} dt={dt} failed: {folder}/oofem.log')
    elapsed=time.perf_counter()-tic; nsteps=round(T/dt); u0=fields(xy); U=[u0]
    text_history=read_text_history(folder/f'{tag}.out',nsteps,len(xy),dt)
    for s in range(1,nsteps+1):
        x,u=read_vtu(folder/f'{tag}.out.m0.{s}.vtu')
        if not np.allclose(x,xy,rtol=0,atol=1e-12): raise RuntimeError('VTU node ordering differs from input')
        U.append(text_history[s-1][1])
    U=np.asarray(U); V=np.zeros_like(U)
    if kind=='hyperbolic':
        for j in range(nsteps): V[j+1]=2*(U[j+1]-U[j])/dt-V[j]
    else:
        V[1:]=(U[1:]-U[:-1])/dt
    Mb,Ms=matrices(xy,quads,edges); series=[]
    Db=Ds=0.
    for j in range(nsteps+1):
        parts=energy_parts(U[j],V[j],xy,quads,edges)
        ebg,ebp,ek,esg,esp=parts
        if kind=='parabolic':
            ek=0.0
        Et=ebg+ebp+ek+esg+esp
        if j:
            if kind=='parabolic':
                du=U[j]-U[j-1]
                db=C*(du@Mb@du)/dt
                ds=CS*(du@Ms@du)/dt
            else:
                # Newmark average acceleration makes velocity linear over
                # each step; integrate its square exactly for that segment.
                v0,v1=V[j-1],V[j]
                integral_b=dt/3*(v0@Mb@v0+v0@Mb@v1+v1@Mb@v1)
                integral_s=dt/3*(v0@Ms@v0+v0@Ms@v1+v1@Ms@v1)
                db=RHO*DELTA*integral_b
                ds=CS*integral_s
            Db+=db; Ds+=ds
        else: db=ds=0.
        if j==0: E0=Et
        residual=Et-E0+Db+Ds
        series.append({'kind':kind,'mesh':f'{nx}x{ny}','dt':dt,'time':j*dt,'bulk_gradient':ebg,'bulk_potential':ebp,
          'bulk_kinetic':ek,'surface_gradient':esg,'surface_potential':esp,'energy_total':Et,
          'dissipation_bulk_cumulative':Db,'dissipation_surface_cumulative':Ds,'balance_residual':residual,
          'normalized_balance_residual':residual/max(abs(E0),1e-14)})
    return series,elapsed

def main():
    if not OOFEM.exists(): raise FileNotFoundError(OOFEM)
    OUT.mkdir(exist_ok=True); allrows=[]; summary=[]
    # Time refinement at fixed mesh, and spatial refinement at fixed dt.
    cases=[]
    for kind in ('parabolic','hyperbolic'):
        cases += [(kind,16,8,dt) for dt in (.002,.001,.0005)]
        cases += [(kind,nx,ny,.0005) for nx,ny in ((8,4),(32,16))]
    for kind,nx,ny,dt in cases:
        series,elapsed=run_case(kind,nx,ny,dt); allrows.extend(series)
        end=series[-1]; summary.append({'kind':kind,'mesh':f'{nx}x{ny}','dt':dt,'T':T,
          'E0':series[0]['energy_total'],'ET':end['energy_total'],
          'bulk_dissipation':end['dissipation_bulk_cumulative'],'surface_dissipation':end['dissipation_surface_cumulative'],
          'balance_residual':end['balance_residual'],'relative_balance_residual':end['normalized_balance_residual'],
          'max_abs_residual':max(abs(r['balance_residual']) for r in series),'runtime_seconds':elapsed})
        print(summary[-1],flush=True)
    for name,rows in [(f'{OUTPUT_PREFIX}_time_series.csv',allrows),(f'{OUTPUT_PREFIX}_summary.csv',summary)]:
        with (HERE/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for row,kind in enumerate(('parabolic','hyperbolic')):
        for dt in (.002,.001,.0005):
            sub=[r for r in allrows if r['kind']==kind and r['mesh']=='16x8' and r['dt']==dt]
            axs[row,0].plot([r['time'] for r in sub],[r['energy_total']-sub[0]['energy_total'] for r in sub],label=f'dt={dt:g}')
            axs[row,1].plot([r['time'] for r in sub],[r['normalized_balance_residual'] for r in sub],label=f'dt={dt:g}')
        axs[row,0].set_title(f'{kind}: E(t)-E(0)'); axs[row,1].set_title(f'{kind}: normalized balance residual')
        for ax in axs[row]: ax.set_xlabel('time'); ax.grid(True,alpha=.3); ax.legend()
    fig.savefig(HERE/f'{OUTPUT_PREFIX}_temporal.png',dpi=180); plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for row,kind in enumerate(('parabolic','hyperbolic')):
        for meshname in ('8x4','16x8','32x16'):
            sub=[r for r in allrows if r['kind']==kind and r['mesh']==meshname and r['dt']==.0005]
            if not sub: continue
            axs[row,0].plot([r['time'] for r in sub],[r['energy_total']-sub[0]['energy_total'] for r in sub],label=meshname)
            axs[row,1].plot([r['time'] for r in sub],[r['normalized_balance_residual'] for r in sub],label=meshname)
        axs[row,0].set_title(f'{kind}: spatial refinement, E(t)-E(0)')
        axs[row,1].set_title(f'{kind}: spatial refinement, balance residual')
        for ax in axs[row]: ax.set_xlabel('time'); ax.grid(True,alpha=.3); ax.legend(title='mesh')
    fig.savefig(HERE/f'{OUTPUT_PREFIX}_mesh.png',dpi=180); plt.close(fig)

if __name__=='__main__': main()
