#!/usr/bin/env python3
"""Finite-element-in-space WED minimizer for OOFEM Allen-Cahn input meshes.

Supports 2D quad1ht/line1ht and tri3ht/line1ht scalar meshes with one
Allen-Cahn material. It minimizes the paper's weighted functional on
0 <= s=t/epsilon <= Smax using P1 time elements and a free terminal value.
This is an isolated research solver; it does not change OOFEM's transient engine.
"""
from __future__ import annotations
import argparse, math, re
import warnings
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags, bmat, eye
from scipy.sparse.linalg import spsolve, MatrixRankWarning


def read_oofem_mesh(path):
    lines=Path(path).read_text().splitlines()
    nodes={}; elems=[]; initials={}; mat=None
    for line in lines:
        s=line.split()
        if not s or s[0].startswith('#'): continue
        if s[0]=='node': nodes[int(s[1])]=np.array([float(s[4]),float(s[5])])
        elif s[0] in ('quad1ht','quad1hmt','quad1mt','tri3ht','tri3mt','line1ht','line1mt'):
            count=int(s[s.index('nodes')+1]); ids=list(map(int,s[s.index('nodes')+2:s.index('nodes')+2+count])); elems.append((s[0],ids))
        elif s[0]=='allencahn':
            kv=dict(zip(s[2::2],s[3::2])); mat={k:float(v) for k,v in kv.items() if k in {'k','c','mobility','surfacebeta'}}
        elif s[0]=='InitialCondition':
            ids=s.index('u'); val=float(s[ids+1]); setid=int(s[s.index('set')+1]); initials[setid-1]=val
    if not nodes or not elems or mat is None: raise ValueError('input needs nodes, transport elements and one allencahn material')
    node_ids=sorted(nodes); index={n:i for i,n in enumerate(node_ids)}
    xy=np.array([nodes[n] for n in node_ids]); u0=np.empty(len(node_ids))
    for n in node_ids:
        if n not in initials: raise ValueError(f'missing scalar initial condition for node {n}')
        u0[index[n]]=initials[n]
    return xy,index,elems,u0,mat


def quadrature(kind, X):
    if kind.startswith('quad1'):
        # 3x3 Gauss exactly integrates the quartic potential and its cubic
        # derivative for bilinear Q1 fields. The 2x2 rule is insufficient.
        a=math.sqrt(3/5); points=((-a,5/9),(0.,8/9),(a,5/9)); out=[]
        for r,wr in points:
            for s,ws in points:
                N=.25*np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])
                dN=.25*np.array([[-(1-s),-(1-r)],[(1-s),-(1+r)],[(1+s),(1+r)],[-(1+s),(1-r)]])
                J=X.T@dN; det=np.linalg.det(J)
                if det<=0: raise ValueError('quad has non-positive Jacobian')
                grad=dN@np.linalg.inv(J)
                out.append((N,grad,det*wr*ws))
        return out
    if kind.startswith('tri3'):
        # Three-point degree-two rule on the reference triangle.
        bary=((2/3,1/6,1/6),(1/6,2/3,1/6),(1/6,1/6,2/3)); dN=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
        J=X.T@dN; det=np.linalg.det(J)
        if det<=0: raise ValueError('triangle has non-positive Jacobian')
        grad=dN@np.linalg.inv(J)
        return [(np.array(N),grad,det/6) for N in bary]
    if kind.startswith('line1'):
        # Three-point Gauss integrates quartic line potentials exactly.
        length=np.linalg.norm(X[1]-X[0]); a=math.sqrt(3/5); t=(X[1]-X[0])/length
        points=((-a,5/9),(0.,8/9),(a,5/9))
        return [(np.array([(1-r)/2,(1+r)/2]),np.outer(np.array([-1.,1.])/length,t),length*w/2) for r,w in points]
    raise ValueError(f'unsupported element {kind}')


def assemble(path):
    xy,index,elems,u0,mat=read_oofem_mesh(path)
    n=len(xy); rows=[]; cols=[]; mv=[]; kv=[]; qrows=[]; qcols=[]; qvals=[]; pweights=[]; nq=0
    k=mat.get('k',1.); cap=mat.get('c',1.); mob=mat.get('mobility',1.); beta=mat.get('surfacebeta',1.)
    for kind, ids0 in elems:
        ids=np.array([index[i] for i in ids0]); X=xy[ids]
        for N,grad,w in quadrature(kind,X):
            mm=cap*np.outer(N,N)*w
            kk=k*(grad@grad.T)*w
            for a,ia in enumerate(ids):
                for b,ib in enumerate(ids):
                    rows.append(ia); cols.append(ib); mv.append(mm[a,b]); kv.append(kk[a,b])
            for a,ia in enumerate(ids): qrows.append(nq); qcols.append(ia); qvals.append(N[a])
            pweights.append(mob*w*(beta if kind.startswith('line1') else 1.)); nq+=1
    M=coo_matrix((mv,(rows,cols)),shape=(n,n)).tocsr(); K=coo_matrix((kv,(rows,cols)),shape=(n,n)).tocsr()
    Nq=coo_matrix((qvals,(qrows,qcols)),shape=(nq,n)).tocsr()
    return M,K,Nq,np.asarray(pweights),u0


def solve_wed(input_path, eps, smax, ds, maxiter=3000, tol=2e-5):
    M,K,Nq,pweights,u0=assemble(input_path); n=len(u0); nt=round(smax/ds)
    if nt<2 or not math.isclose(nt*ds,smax,rel_tol=1e-10,abs_tol=1e-12): raise ValueError('Smax must be an integer multiple of ds')
    h=smax/nt; tau=np.array([.5*(1-1/math.sqrt(3)),.5*(1+1/math.sqrt(3))]); weights=np.array([.5,.5])
    def objective(z, want_hessian=False):
        U=np.empty((nt+1,n)); U[0]=u0; U[1:]=z.reshape(nt,n)
        G=np.zeros_like(U); J=0.
        blocks={}
        def add_block(a,b,value):
            if a<0 or b<0: return
            key=(a,b); blocks[key]=value if key not in blocks else blocks[key]+value
        for i in range(nt):
            dU=(U[i+1]-U[i])/h
            for tq,wq in zip(tau,weights):
                s=(i+tq)*h; wt=h*wq*math.exp(-s); u=(1-tq)*U[i]+tq*U[i+1]
                Mu=M@dU; Ku=K@u
                energy=.5*float(dU@Mu)+.5*eps*float(u@Ku)
                uq=Nq@u
                pot=.25*float(pweights@((uq*uq-1)**2))
                gp=Nq.T@(pweights*(uq**3-uq))
                J+=wt*(energy+eps*pot)
                gu=eps*(Ku+gp)
                G[i] += wt*((-1/h)*Mu+(1-tq)*gu)
                G[i+1] += wt*((1/h)*Mu+tq*gu)
                if want_hessian:
                    Kpot=Nq.T@diags(pweights*(3*uq*uq-1))@Nq
                    Hsp=K+Kpot
                    Htt=M/(eps*h*h)
                    ii=i-1; jj=i
                    add_block(ii,ii,wt*(Htt+(1-tq)**2*Hsp))
                    add_block(jj,jj,wt*(Htt+tq*tq*Hsp))
                    add_block(ii,jj,wt*(-Htt+tq*(1-tq)*Hsp))
                    add_block(jj,ii,wt*(-Htt+tq*(1-tq)*Hsp))
        # Scaling the objective does not change its minimizer and prevents
        # small epsilon values from meeting the optimizer's absolute gradient
        # tolerance before the path has converged.
        if not want_hessian:
            return J/eps,G[1:].ravel()/eps
        matrix=bmat([[blocks.get((i,j),None) for j in range(nt)] for i in range(nt)],format='csc')
        return J/eps,G[1:].ravel()/eps,matrix
    z=np.tile(u0,(nt,1)).ravel()
    direction=np.random.default_rng(17).normal(size=z.size); direction/=np.linalg.norm(direction)
    fdstep=1.e-6
    f0,g0,H0=objective(z,True)
    fp,gp=objective(z+fdstep*direction); fm,gm=objective(z-fdstep*direction)
    gradcheck=abs((fp-fm)/(2*fdstep)-float(g0@direction))/(1+abs(float(g0@direction)))
    if gradcheck>2.e-5: raise RuntimeError(f'WED gradient check failed: relative directional error {gradcheck:g}')
    hesscheck=np.linalg.norm((gp-gm)/(2*fdstep)-H0@direction)/(1+np.linalg.norm(H0@direction))
    if hesscheck>2.e-4: raise RuntimeError(f'WED Hessian check failed: relative directional error {hesscheck:g}')
    success=False; message='maximum Newton iterations reached'; Jval=np.inf; gradnorm=np.inf
    for iteration in range(maxiter):
        Jval,g,H=objective(z,True); gradnorm=np.linalg.norm(g,np.inf)
        if gradnorm<tol:
            success=True; message='Newton gradient tolerance satisfied'; break
        direction_new=None; alpha=1.0
        for shift in (0.,1.e-10,1.e-8,1.e-6,1.e-4,1.e-2,1.):
            try:
                A=H if shift==0 else H+shift*eye(H.shape[0],format='csc')
                with warnings.catch_warnings():
                    warnings.simplefilter('error',MatrixRankWarning)
                    candidate=spsolve(A,-g)
                if np.all(np.isfinite(candidate)) and float(g@candidate)<0:
                    direction_new=candidate; break
            except Exception:
                continue
        if direction_new is None:
            direction_new=-g
        slope=float(g@direction_new); accepted=False
        for _ in range(30):
            trial=z+alpha*direction_new; Jtrial,_=objective(trial)
            if np.isfinite(Jtrial) and Jtrial<=Jval+1.e-4*alpha*slope:
                z=trial; accepted=True; break
            alpha*=.5
        if not accepted:
            message='Newton line search failed'; break
    Jval,g,H=objective(z,True); gradnorm=np.linalg.norm(g,np.inf)
    if gradnorm<tol: success=True; message='Newton gradient tolerance satisfied'
    class Result: pass
    result=Result(); result.success=success; result.nit=iteration+1; result.fun=Jval; result.message=message
    U=np.vstack((u0,z.reshape(nt,n)))
    return result,U,gradnorm,gradcheck,hesscheck


def write_vtk_series(input_path, prefix, times, solution):
    xy,index,elems,u0,mat=read_oofem_mesh(input_path)
    vtk_type={'quad1ht':9,'quad1hmt':9,'quad1mt':9,'tri3ht':5,'tri3mt':5,'line1ht':3,'line1mt':3}
    cells=[]; types=[]; offsets=[]; nconn=0
    for kind,ids in elems:
        cells.extend(index[i] for i in ids); nconn+=len(ids); offsets.append(nconn); types.append(vtk_type[kind])
    outprefix=Path(prefix); outprefix.parent.mkdir(parents=True,exist_ok=True); datasets=[]
    for step,(time,u) in enumerate(zip(times,solution)):
        name=f'{outprefix.name}.m0.{step}.vtu'; dest=outprefix.parent/name
        with dest.open('w') as f:
            f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian"><UnstructuredGrid>\n')
            f.write(f'<Piece NumberOfPoints="{len(xy)}" NumberOfCells="{len(types)}"><Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">')
            f.write(' '.join(f'{x:.12g} {y:.12g} 0' for x,y in xy)); f.write('</DataArray></Points><Cells>')
            f.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+' '.join(map(str,cells))+'</DataArray>')
            f.write('<DataArray type="Int32" Name="offsets" format="ascii">'+' '.join(map(str,offsets))+'</DataArray>')
            f.write('<DataArray type="UInt8" Name="types" format="ascii">'+' '.join(map(str,types))+'</DataArray></Cells>')
            f.write('<PointData Scalars="Temperature"><DataArray type="Float64" Name="Temperature" format="ascii">')
            f.write(' '.join(f'{v:.12g}' for v in u)); f.write('</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>\n')
        datasets.append(f'<DataSet timestep="{time:.12g}" group="" part="0" file="{name}"/>')
    with Path(str(outprefix)+'.pvd').open('w') as f:
        f.write('<?xml version="1.0"?><VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>')
        f.write(''.join(datasets)); f.write('</Collection></VTKFile>\n')


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('input'); p.add_argument('--eps',type=float,required=True); p.add_argument('--Smax',type=float,default=8.0,help='terminal time in rescaled variable s=t/epsilon'); p.add_argument('--ds',type=float,default=0.1); p.add_argument('--output',default='wed_path.npz'); p.add_argument('--maxiter',type=int,default=3000); p.add_argument('--tol',type=float,default=2.e-5,help='infinity-norm stopping tolerance for the scaled WED gradient')
    a=p.parse_args(); result,U,gn,gc,hc=solve_wed(a.input,a.eps,a.Smax,a.ds,a.maxiter,a.tol)
    times=a.eps*np.linspace(0,a.Smax,U.shape[0])
    np.savez(a.output,time=times,rescaled_time=np.linspace(0,a.Smax,U.shape[0]),solution=U,epsilon=a.eps,objective=result.fun,success=result.success,message=result.message,iterations=result.nit,gradient_inf=gn,gradient_check_relative=gc,hessian_check_relative=hc)
    write_vtk_series(a.input,Path(a.output).with_suffix(''),times,U)
    print(f'WED success={result.success}; iterations={result.nit}; objective={result.fun:.12g}; |grad|inf={gn:.4e}; grad check={gc:.3e}; Hessian check={hc:.3e}; message={result.message}')
    if not result.success: raise SystemExit(2)

if __name__=='__main__': main()
