#!/usr/bin/env python3
"""Independent semidiscrete reference for OOFEM's inertial-bulk/AC-surface model.

The reference assembles the Q1 bulk and P1 boundary energy forms in Python,
then integrates the resulting finite-dimensional ODE with SciPy DOP853. It does
not reuse OOFEM element/assembly code. OOFEM inputs and full VTU histories are
written to new coupled_independent_* folders.
"""
from __future__ import annotations
import argparse, csv, math, subprocess, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OOFEM=ROOT/'build/debug_magneto/oofem'
LX,LY,T=2.,1.,.02
RHO=1.; ETA=0.; KB=1.; MB=1.; CS=1.; KS=1.; MS=1.
BETA,GAMMA=.25,.5
GAUSS_X=np.array([-math.sqrt(3/5),0.,math.sqrt(3/5)])
GAUSS_W=np.array([5/9,8/9,5/9])


def qshape(xi,eta):
    N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
    dN=.25*np.array([[-(1-eta),-(1-xi)],[(1-eta),-(1+xi)],[(1+eta),(1+xi)],[-(1+eta),(1-xi)]])
    return N,dN


def mesh(nx,ny):
    xs=np.linspace(0,LX,nx+1); ys=np.linspace(0,LY,ny+1)
    coords=np.array([(x,y) for y in ys for x in xs])
    node=lambda i,j:j*(nx+1)+i
    quads=[]
    for j in range(ny):
        for i in range(nx): quads.append([node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)])
    edges=[]
    for i in range(nx): edges.append([node(i,0),node(i+1,0)])
    for j in range(ny): edges.append([node(nx,j),node(nx,j+1)])
    for i in range(nx,0,-1): edges.append([node(i,ny),node(i-1,ny)])
    for j in range(ny,0,-1): edges.append([node(0,j),node(0,j-1)])
    boundary=np.unique(np.asarray(edges).ravel())
    return coords,np.asarray(quads,dtype=int),np.asarray(edges,dtype=int),boundary


def assemble(nx,ny):
    coords,quads,edges,boundary=mesh(nx,ny); nd=len(coords)
    Mb=lil_matrix((nd,nd)); Kb=lil_matrix((nd,nd)); Ms=lil_matrix((nd,nd)); Ks=lil_matrix((nd,nd))
    for ids in quads:
        xy=coords[ids]
        for a,xi in enumerate(GAUSS_X):
            for b,eta in enumerate(GAUSS_X):
                N,dN=qshape(xi,eta); J=xy.T@dN; det=np.linalg.det(J); grad=dN@np.linalg.inv(J); w=det*GAUSS_W[a]*GAUSS_W[b]
                Mb[np.ix_(ids,ids)] += RHO*np.outer(N,N)*w
                Kb[np.ix_(ids,ids)] += KB*(grad@grad.T)*w
    for ids in edges:
        xy=coords[ids]; tangent=xy[1]-xy[0]; length=np.linalg.norm(tangent)
        for a,s in enumerate(GAUSS_X):
            N=np.array([(1-s)/2,(1+s)/2]); dN=np.array([-.5,.5])/length; w=length*GAUSS_W[a]
            Ms[np.ix_(ids,ids)] += CS*np.outer(N,N)*w
            Ks[np.ix_(ids,ids)] += KS*np.outer(dN,dN)*w
    return coords,quads,edges,boundary,Mb.tocsc(),Kb.tocsc(),Ms.tocsc(),Ks.tocsc()


def residual_reaction(u,coords,quads,edges):
    rb=np.zeros(len(u)); rs=np.zeros(len(u))
    for ids in quads:
        xy=coords[ids]; ul=u[ids]
        for a,xi in enumerate(GAUSS_X):
            for b,eta in enumerate(GAUSS_X):
                N,dN=qshape(xi,eta); det=np.linalg.det(xy.T@dN); uq=N@ul
                rb[ids] += MB*N*(uq-uq**3)*det*GAUSS_W[a]*GAUSS_W[b]
    for ids in edges:
        xy=coords[ids]; length=np.linalg.norm(xy[1]-xy[0]); ul=u[ids]
        for a,s in enumerate(GAUSS_X):
            N=np.array([(1-s)/2,(1+s)/2]); uq=N@ul
            rs[ids] += MS*N*(uq-uq**3)*length*GAUSS_W[a]
    return rb,rs


def initial_u(coords):
    x,y=coords[:,0],coords[:,1]
    return .2+.1*np.cos(math.pi*x/LX)*np.cos(math.pi*y/LY)


def reference(nx,ny,times):
    coords,quads,edges,boundary,Mb,Kb,Ms,Ks=assemble(nx,ny)
    mb_lu=splu(Mb)
    C=ETA*Mb+Ms; K=Kb+Ks
    def rhs(t,z):
        u=z[:len(coords)]; v=z[len(coords):]
        rb,rs=residual_reaction(u,coords,quads,edges)
        a=mb_lu.solve(-(C@v)-(K@u)+rb+rs)
        return np.r_[v,a]
    y0=np.r_[initial_u(coords),np.zeros(len(coords))]
    sol=solve_ivp(rhs,(0,float(times[-1])),y0,t_eval=times,method='DOP853',rtol=2e-12,atol=2e-14)
    if not sol.success: raise RuntimeError(sol.message)
    return coords,edges,boundary,sol.y[:len(coords),:].T,sol.y[len(coords):,:].T,(Mb,Kb,Ms,Ks)


def make_input(nx,ny,dt,tag):
    coords,quads,edges,boundary=mesh(nx,ny); nn=len(coords); ne=len(quads)+len(edges); nsteps=round(T/dt)
    lines=['case.out',f'Coupled hyperbolic bulk / Allen-Cahn surface independent-reference test {tag}',
      # The assembled residual reaches a ~1e-10 roundoff plateau on the
      # refined meshes; a looser nonlinear tolerance avoids false failure.
      f'nltransienttransportproblem nsteps {nsteps} deltat {dt:.16g} alpha 1.0 rtol 1.e-8 nsmax 80 nmodules 1 wave 1 newmarkbeta {BETA} newmarkgamma {GAMMA} initialvelocity 0.0',
      'vtkxml tstep_all domain_all primvars 1 6','domain heattransfer','OutputManager tstep_all dofman_all element_all',
      f'ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}']
    for n,(x,y) in enumerate(coords,1): lines.append(f'node {n} coords 3 {x:.16g} {y:.16g} 0.0')
    for e,ids in enumerate(quads,1): lines.append(f'quad1ht {e} nodes 4 '+' '.join(str(i+1) for i in ids)+' nip 9')
    for e,ids in enumerate(edges,len(quads)+1): lines.append(f'line1ht {e} nodes 2 '+' '.join(str(i+1) for i in ids)+' nip 3')
    lines += ['SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
      'allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 surfacebeta 1.0 manufactured 0 hyperbolic 1 rho 1.0 surfacecapacity 1.0']
    u0=initial_u(coords)
    for n,u in enumerate(u0,1): lines.append(f'InitialCondition {n} dofs 1 10 Conditions 1 u {u:.17g} set {n+1}')
    lines += ['ConstantFunction 1 f(t) 1.0',f'Set 1 elementranges {{(1 {ne})}}']
    lines += [f'Set {n+1} nodes 1 {n}' for n in range(1,nn+1)]
    return '\n'.join(lines)+'\n',coords,quads,edges,boundary


def read_final(folder):
    pvd=ET.parse(folder/'case.out.m0.pvd').getroot(); ds=pvd.findall('.//DataSet')
    vtu=folder/ds[-1].attrib['file']; root=ET.parse(vtu).getroot()
    coords=np.fromstring(root.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)[:,:2]
    arr=next(a for a in root.findall('.//PointData/DataArray') if a.attrib.get('Name')=='Temperature')
    u=np.fromstring(arr.text or '',sep=' ')
    return float(ds[-1].attrib['timestep']),coords,u


def run_oofem(nx,ny,dt,tag):
    folder=HERE/f'coupled_independent_{tag}'; folder.mkdir(exist_ok=True)
    text,coords,quads,edges,boundary=make_input(nx,ny,dt,tag); (folder/'case.in').write_text(text)
    p=subprocess.run([str(OOFEM),'-f',str(folder/'case.in')],cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (folder/'run.log').write_text(p.stdout)
    if p.returncode or 'ANALYSIS FINISHED' not in p.stdout: raise RuntimeError(f'OOFEM failed: {folder}/run.log')
    t,x,u=read_final(folder)
    if abs(t-T)>1e-10 or not np.allclose(x,coords,rtol=0,atol=1e-13): raise RuntimeError('OOFEM time/coordinates mismatch')
    np.savez(folder/'mesh.npz',coords=coords,quads=quads,edges=edges,boundary=boundary)
    return folder,coords,quads,edges,boundary,u


def rms_fe(M,e,measure): return float(math.sqrt(max(float(e@(M@e))/measure,0.)))
def maxerr(e): return float(np.max(np.abs(e)))


def compare_row(nx,ny,dt,tag,coords,boundary,u,uref,Mb,Ms):
    eb=u-uref; mask=np.zeros(len(u),dtype=bool); mask[boundary]=True
    return {'tag':tag,'nx':nx,'ny':ny,'dt':dt,'T':T,
      'bulk_L2_RMS_vs_independent':rms_fe(Mb,eb,LX*LY),
      'surface_L2_RMS_vs_independent':rms_fe(Ms,eb,2*(LX+LY)),
      'bulk_max_nodal':maxerr(eb[np.flatnonzero(~mask)]),
      'surface_max_nodal':maxerr(eb[boundary]),'all_dof_max':maxerr(eb)}


def energy(u,v,coords,quads,edges):
    ek=ep=es=esp=0.
    for ids in quads:
        xy=coords[ids]; ul=u[ids]; vl=v[ids]
        for a,xi in enumerate(GAUSS_X):
            for b,eta in enumerate(GAUSS_X):
                N,dN=qshape(xi,eta); J=xy.T@dN; det=np.linalg.det(J); grad=dN@np.linalg.inv(J); uq=N@ul; vq=N@vl
                ek+=.5*RHO*vq*vq*det*GAUSS_W[a]*GAUSS_W[b]
                ep+=.5*KB*(grad.T@ul)@(grad.T@ul)*det*GAUSS_W[a]*GAUSS_W[b]
                es+=.25*MB*(uq*uq-1.)**2*det*GAUSS_W[a]*GAUSS_W[b]
    for ids in edges:
        length=np.linalg.norm(coords[ids[1]]-coords[ids[0]]); ul=u[ids]
        for a,s in enumerate(GAUSS_X):
            N=np.array([(1-s)/2,(1+s)/2]); dN=np.array([-.5,.5])/length; uq=N@ul
            esp+=.5*KS*(dN@ul)**2*length*GAUSS_W[a]+.25*MS*(uq*uq-1.)**2*length*GAUSS_W[a]
    return ek+ep+es+esp,ek,ep+es,esp


def save_reference(tag,nx,ny,dt):
    nsteps=round(T/dt); times=np.linspace(0,T,nsteps+1)
    coords,edges,boundary,U,V,mats=reference(nx,ny,times)
    folder=HERE/f'coupled_independent_{tag}'
    np.savez(folder/'independent_reference.npz',times=times,coords=coords,U=U,V=V,boundary=boundary)
    quads=mesh(nx,ny)[1]
    energies=np.array([energy(U[i],V[i],coords,quads,edges) for i in range(len(times))])
    with (folder/'independent_energy.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['time','total','bulk_kinetic','bulk_potential_gradient','surface']);
        for t,e in zip(times,energies): w.writerow([t,*e])
    return U[-1],(coords,quads,edges,boundary,*mats),times,U,V


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--suite',choices=['all','space','time'],default='all'); args=ap.parse_args()
    if not OOFEM.exists(): raise FileNotFoundError(OOFEM)
    if args.suite in ('all','space'):
        rows=[]; snapshots=[]
        for nx,ny in ((8,4),(16,8),(32,16)):
            dt=.00025; tag=f'space_n{nx}x{ny}_dt0p00025'
            folder,coords,quads,edges,boundary,u=run_oofem(nx,ny,dt,tag)
            uref,ctx,times,U,V=save_reference(tag,nx,ny,dt)
            Mb,Kb,Ms,Ks=ctx[-4:]
            rows.append(compare_row(nx,ny,dt,tag,coords,boundary,u,uref,Mb,Ms))
            snapshots.append((nx,ny,coords,u,boundary,Mb,Ms))
        # Internal spatial convergence: common-node comparison to next mesh.
        for i in range(2):
            n0,ny0,x0,u0,b0,Mb0,Ms0=snapshots[i]; n1,ny1,x1,u1,b1,Mb1,Ms1=snapshots[i+1]
            lookup={tuple(np.round(x,12)):j for j,x in enumerate(x1)}
            idx=np.array([lookup[tuple(np.round(x,12))] for x in x0]); d=u0-u1[idx]
            mask=np.zeros(len(x0),bool); mask[b0]=True
            rows[i]['bulk_RMS_to_next_finer']=rms_fe(Mb0,d,LX*LY)
            rows[i]['surface_RMS_to_next_finer']=rms_fe(Ms0,d,2*(LX+LY))
            rows[i]['bulk_max_to_next_finer']=maxerr(d[~mask]); rows[i]['surface_max_to_next_finer']=maxerr(d[b0])
        for key in ('bulk_RMS_to_next_finer','surface_RMS_to_next_finer'):
            a,b=rows[0][key],rows[1][key]; rows[0][key+'_order']=math.log(a/b,2.)
        with (HERE/'coupled_independent_spatial.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r))); w.writeheader(); w.writerows(rows)
        print('spatial',rows,flush=True)
    if args.suite in ('all','time'):
        rows=[]
        nx,ny=16,8
        for dt in (.002,.001,.0005):
            tag=f'time_n{nx}x{ny}_dt{dt:g}'
            folder,coords,quads,edges,boundary,u=run_oofem(nx,ny,dt,tag)
            uref,ctx,times,U,V=save_reference(tag,nx,ny,dt)
            Mb,Kb,Ms,Ks=ctx[-4:]
            rows.append(compare_row(nx,ny,dt,tag,coords,boundary,u,uref,Mb,Ms))
        rows.sort(key=lambda r:r['dt'],reverse=True)
        for key in ('bulk_L2_RMS_vs_independent','surface_L2_RMS_vs_independent'):
            for i in range(len(rows)-1):
                rows[i][key+'_order_to_next']=math.log(rows[i][key]/rows[i+1][key],2.)
        with (HERE/'coupled_independent_temporal.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r))); w.writeheader(); w.writerows(rows)
        print('temporal',rows,flush=True)

if __name__=='__main__': main()
