#!/usr/bin/env python3
"""Natural Allen-Cahn/dynamic-boundary Q1 mesh refinement, compared to fine OOFEM."""
from pathlib import Path
import math, subprocess, xml.etree.ElementTree as ET, csv
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=Path('/tmp/stvowu_literature_spatial'); OUT.mkdir(exist_ok=True)
OOFEM=ROOT/'build/debug_magneto/oofem'; rows=[]; solutions={}
for n in (10,20,40):
 d=OUT/f'{n}x{n}'; d.mkdir(exist_ok=True)
 nodes=(n+1)**2; quads=n*n; edges=4*n; elems=quads+edges; lines=[]
 lines += [f'case_{n}.out','Natural Allen-Cahn relaxation for dynamic-boundary FEM refinement',
           'nltransienttransportproblem nsteps 40 deltat 0.0025 alpha 1.0 rtol 1.e-10 nsmax 40 nmodules 1',
           'vtkxml tstep_all domain_all primvars 1 6','domain heattransfer','OutputManager tstep_all dofman_all element_all',
           f'ndofman {nodes} nelem {elems} ncrosssect 1 nmat 1 nbc 0 nic {n+1} nltf 1 nset {n+2}']
 for j in range(n+1):
  for i in range(n+1):
   k=j*(n+1)+i+1; lines.append(f'node {k} coords 3 {i/n:.10f} {j/n:.10f} 0.0')
 eid=1
 for j in range(n):
  for i in range(n):
   a=j*(n+1)+i+1; b=a+1; dd=a+n+1; c=dd+1
   lines.append(f'quad1ht {eid} nodes 4 {a} {b} {c} {dd} nip 9'); eid+=1
 for i in range(n): lines.append(f'line1ht {eid} nodes 2 {i+1} {i+2} nip 3'); eid+=1
 for j in range(n): lines.append(f'line1ht {eid} nodes 2 {j*(n+1)+n+1} {(j+1)*(n+1)+n+1} nip 3'); eid+=1
 for i in range(n,0,-1): lines.append(f'line1ht {eid} nodes 2 {n*(n+1)+i+1} {n*(n+1)+i} nip 3'); eid+=1
 for j in range(n,0,-1): lines.append(f'line1ht {eid} nodes 2 {j*(n+1)+1} {(j-1)*(n+1)+1} nip 3'); eid+=1
 lines += ['SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1','allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0']
 for i in range(n+1):
  val=.1+.05*math.cos(2*math.pi*i/n)
  lines.append(f'InitialCondition {i+1} dofs 1 10 Conditions 1 u {val:.15g} set {i+2}')
 lines.append('ConstantFunction 1 f(t) 1.0')
 lines.append(f'Set 1 elementranges {{(1 {elems})}}')
 for i in range(n+1):
  ids=[j*(n+1)+i+1 for j in range(n+1)]
  lines.append(f'Set {i+2} nodes {n+1} '+' '.join(map(str,ids)))
 inp=d/'case.in'; inp.write_text('\n'.join(lines)+'\n')
 p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=d,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (d/'run.log').write_text(p.stdout)
 if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout: raise SystemExit(f'OOFEM failed for {n}x{n}; see {d}/run.log')
 f=d/f'case_{n}.out.m0.40.vtu'; root=ET.parse(f).getroot(); piece=root.find('.//Piece')
 u=np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))
 solutions[n]=u
 u0=np.array([.1+.05*math.cos(2*math.pi*i/n) for j in range(n+1) for i in range(n+1)])
 rows.append({'n':n,'h':1/n,'nodes':nodes,'elements':elems,'T':.1,'dt':.0025,'initial_rms':float(np.sqrt(np.mean(u0*u0)))})
ref=solutions[40]
for row,n in zip(rows,(10,20,40)):
 stride=40//n; coarse=solutions[n]; sample=np.array([ref[j*41*stride+i*stride] for j in range(n+1) for i in range(n+1)])
 delta=coarse-sample; row['nodal_rms_vs_40']=float(np.sqrt(np.mean(delta**2))); row['max_abs_vs_40']=float(np.max(np.abs(delta)))
with (Path(__file__).with_name('spatial_convergence.csv')).open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
