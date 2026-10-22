#!/usr/bin/env python3
"""Run the direct Allen-Cahn surface-potential parameter study."""
from pathlib import Path
import csv, math, re, subprocess, time
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OOFEM=ROOT/'build/debug_magneto/oofem'
BASE=Path(__file__).with_name('allencahn_rectangle_dynamic.in')
OUT=Path('/tmp/stvowu_beta_study')

def read_mesh():
    lines=BASE.read_text().splitlines(); nodes={}; edges=[]
    for line in lines:
        s=line.split()
        if s[:1]==['node']: nodes[int(s[1])]=np.array([float(s[4]),float(s[5])])
        if s[:1]==['line1ht']:
            j=s.index('nodes'); edges.append([int(x) for x in s[j+2:j+4]])
    return nodes,edges

def values(path, step):
    root=ET.parse(path/f'beta_{path.name.split("_")[-1]}.out.m0.{step}.vtu').getroot(); piece=root.find('.//Piece')
    return np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))

def surface_stats(u, nodes, edges, beta=1.0):
    idmap={n:i for i,n in enumerate(sorted(nodes))}; total_length=energy=mean_num=0.0; a=1/math.sqrt(3)
    for edge in edges:
        ii=[idmap[n] for n in edge]; X=np.array([nodes[n] for n in edge]); length=np.linalg.norm(X[1]-X[0]); ue=u[ii]
        for r in (-a,a):
            N=np.array([(1-r)/2,(1+r)/2]); dN=np.array([-1.,1.]); val=float(N@ue); du=float(dN@ue)/length
            w=length/2; total_length+=w; mean_num+=w*val; energy+=w*(.5*du*du+beta*.25*(val*val-1)**2)
    return mean_num/total_length,energy

def main():
    OUT.mkdir(exist_ok=True); nodes,edges=read_mesh(); rows=[]
    original=BASE.read_text().splitlines()
    for beta in (.5,1.,2.,5.):
        label=str(beta).replace('.','p'); folder=OUT/f'beta_{label}'; folder.mkdir(exist_ok=True)
        lines=original.copy(); lines[0]=f'beta_{label}.out'; lines[2]='nltransienttransportproblem nsteps 40 deltat 0.0005 alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1'
        lines=[re.sub(r'manufactured 1', 'manufactured 0', x) for x in lines]
        lines=[re.sub(r'(allencahn 1 .*?mobility 1\.0)',rf'\1 surfacebeta {beta:g}',x) if x.startswith('allencahn 1 ') else x for x in lines]
        inp=folder/'case.in'; inp.write_text('\n'.join(lines)+'\n')
        start=time.perf_counter(); p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); runtime=time.perf_counter()-start
        (folder/'oofem.log').write_text(p.stdout)
        if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout: raise RuntimeError(f'OOFEM beta={beta} failed; inspect {folder}/oofem.log')
        U0=np.array([float(next(x.split('u ')[1].split()[0] for x in lines if x.startswith(f'InitialCondition {i} '))) for i in range(1,len(nodes)+1)])
        mu0,e0=surface_stats(U0,nodes,edges)
        final=None
        for step in range(1,41):
            u=values(folder,step); mean,e=surface_stats(u,nodes,edges,beta)
            rows.append((beta,step*.0005,mean,e,runtime))
            if step==40: final=(mean,e)
        print(f'beta={beta:g}: M_gamma(T)={final[0]:.9g}; E_gamma(T)={final[1]:.9g}; runtime={runtime:.2f}s')
    with (Path(__file__).with_name('surface_beta.csv')).open('w') as f:
        w=csv.writer(f); w.writerow(['beta','time','boundary_mean','boundary_energy_beta_scaled','runtime_seconds']); w.writerows(rows)
    wed_rows=[]
    wed=Path(__file__).with_name('wed_allencahn.py')
    for beta in (.5,1.,2.,5.):
        label=str(beta).replace('.','p'); folder=OUT/f'beta_{label}'; inp=folder/'case.in'; archive=folder/'wed.npz'
        cmd=['python3',str(wed),str(inp),'--eps','0.05','--Smax','4','--ds','0.1','--output',str(archive),'--maxiter','100']
        start=time.perf_counter(); proc=subprocess.run(cmd,cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); rt=time.perf_counter()-start
        (folder/'wed.log').write_text(proc.stdout)
        if proc.returncode: raise RuntimeError(f'WED beta={beta} failed; inspect {folder}/wed.log')
        z=np.load(archive)
        for t,u in zip(z['time'],z['solution']):
            mean,_=surface_stats(u,nodes,edges,beta); wed_rows.append((beta,float(t),mean,int(z['iterations']),rt,float(z['gradient_inf'])))
        print(f'WED beta={beta:g}: iterations={int(z["iterations"])}, runtime={rt:.2f}s, grad={float(z["gradient_inf"]):.3e}')
    with (Path(__file__).with_name('wed_surface_beta.csv')).open('w') as f:
        w=csv.writer(f); w.writerow(['beta','physical_time','boundary_mean','newton_iterations','runtime_seconds','gradient_inf']); w.writerows(wed_rows)
if __name__=='__main__': main()
