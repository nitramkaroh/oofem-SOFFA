#!/usr/bin/env python3
"""Measure the Allen-Cahn energy balance and its time-step sensitivity."""
from pathlib import Path
import csv, math, re, subprocess, time
import numpy as np
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]; OOFEM=ROOT/'build/debug_magneto/oofem'
BASE=Path(__file__).with_name('energy_decay.in'); OUT=Path('/tmp/stvowu_energy_balance'); T=.1

def geometry(inp):
    nodes={}; elems=[]; ic={}; sets={}
    for line in inp.read_text().splitlines():
        s=line.split()
        if s[:1]==['node']: nodes[int(s[1])]=np.array([float(s[4]),float(s[5])])
        elif s[:1] in (['quad1ht'],['line1ht']):
            q=s.index('nodes'); elems.append((s[0],[int(x) for x in s[q+2:q+2+int(s[q+1])]]))
        elif s[:1]==['InitialCondition']:
            ic[int(s[s.index('set')+1])]=float(s[s.index('u')+1])
        elif s[:1]==['Set'] and 'nodes' in s:
            q=s.index('nodes'); sets[int(s[1])]=[int(x) for x in s[q+2:q+2+int(s[q+1])]]
    u0=np.zeros(len(nodes));
    for sid,val in ic.items():
        for nid in sets[sid]: u0[nid-1]=val
    return nodes,elems,u0

def qpts(kind,X):
    a=1/math.sqrt(3)
    if kind=='quad1ht':
        for r in (-a,a):
            for s in (-a,a):
                N=.25*np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])
                dN=.25*np.array([[-(1-s),-(1-r)],[(1-s),-(1+r)],[(1+s),(1+r)],[-(1+s),(1-r)]])
                J=X.T@dN; yield N,dN@np.linalg.inv(J),np.linalg.det(J)
    else:
        length=np.linalg.norm(X[1]-X[0]);
        for r in (-a,a): yield np.array([(1-r)/2,(1+r)/2]),np.array([-1.,1.])/length,length/2

def assembled(nodes,elems):
    n=len(nodes); Mb=np.zeros((n,n)); Mg=np.zeros((n,n))
    for kind,ids in elems:
        ii=np.array(ids)-1; X=np.array([nodes[k] for k in ids]); M=Mb if kind=='quad1ht' else Mg
        for N,_,w in qpts(kind,X): M[np.ix_(ii,ii)]+=np.outer(N,N)*w
    return Mb,Mg

def energy(U,nodes,elems):
    eb=eg=0.
    for kind,ids in elems:
        ii=np.array(ids)-1; X=np.array([nodes[k] for k in ids]); uv=U[ii]
        # The standard 2x2 rule used for the assembled Allen-Cahn reaction
        # does not integrate the quartic energy F(u_h) exactly for bilinear Q1
        # fields. Use five-point Gauss in each direction (three points on
        # line elements) for the diagnostic energy functional.
        if kind=='quad1ht':
            z,wg=np.polynomial.legendre.leggauss(5)
            def energy_qpts():
                for r,wr in zip(z,wg):
                    for s,ws in zip(z,wg):
                        N=.25*np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])
                        dN=.25*np.array([[-(1-s),-(1-r)],[(1-s),-(1+r)],[(1+s),(1+r)],[-(1+s),(1-r)]])
                        J=X.T@dN
                        yield N,dN@np.linalg.inv(J),np.linalg.det(J)*wr*ws
        else:
            z,wg=np.polynomial.legendre.leggauss(3)
            length=np.linalg.norm(X[1]-X[0])
            def energy_qpts():
                for r,wr in zip(z,wg):
                    yield np.array([(1-r)/2,(1+r)/2]),np.array([-1.,1.])/length,length*wr/2
        for N,grad,w in energy_qpts():
            u=N@uv; gu=uv@grad if kind=='quad1ht' else (uv@grad)
            e=.5*float(gu@gu) if np.ndim(gu) else .5*float(gu*gu)
            e+=.25*(u*u-1)**2
            if kind=='quad1ht': eb+=e*w
            else: eg+=e*w
    return eb,eg

def load_values(file):
    root=ET.parse(file).getroot(); piece=root.find('.//Piece')
    return np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))

def main():
    OUT.mkdir(exist_ok=True); rows=[]
    for dt in (.005,.0025,.00125):
        nsteps=round(T/dt); label=str(dt).replace('.','p'); folder=OUT/f'dt_{label}'; folder.mkdir(exist_ok=True)
        lines=BASE.read_text().splitlines(); lines[0]=f'energy_balance_{label}.out'; lines[2]=f'nltransienttransportproblem nsteps {nsteps} deltat {dt} alpha 1.0 rtol 1.e-10 nsmax 30 nmodules 1'
        # Use exact quadrature for the quartic potential/cubic reaction.
        lines=[line+' nip 9' if line.startswith('quad1ht ') and 'nip' not in line.split() else line for line in lines]
        inp=folder/'case.in'; inp.write_text('\n'.join(lines)+'\n'); start=time.perf_counter()
        p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); runtime=time.perf_counter()-start
        (folder/'oofem.log').write_text(p.stdout)
        if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout: raise RuntimeError(f'energy run dt={dt} failed')
        nodes,elems,u0=geometry(inp); Mb,Mg=assembled(nodes,elems); e0=sum(energy(u0,nodes,elems))
        previous=u0; diss_b=diss_g=0.
        for step in range(1,nsteps+1):
            u=load_values(folder/f'energy_balance_{label}.out.m0.{step}.vtu'); du=u-previous
            diss_b+=float(du@Mb@du)/dt; diss_g+=float(du@Mg@du)/dt; previous=u
        eb,eg=energy(previous,nodes,elems); ef=eb+eg; residual=ef-e0+diss_b+diss_g
        rows.append((dt,nsteps,e0,ef,diss_b,diss_g,residual,abs(residual)/max(abs(e0),1.e-14),runtime))
        print(f'dt={dt:g}: E0={e0:.9g}, ET={ef:.9g}, Db={diss_b:.9g}, Dg={diss_g:.9g}, R={residual:.4e}, |R|/E0={abs(residual)/e0:.4e}, runtime={runtime:.2f}s')
    with (Path(__file__).with_name('energy_balance.csv')).open('w') as f:
        w=csv.writer(f); w.writerow(['dt','steps','E0','ET','dissipation_bulk','dissipation_boundary','residual','normalized_abs_residual','runtime_seconds']); w.writerows(rows)
if __name__=='__main__': main()
