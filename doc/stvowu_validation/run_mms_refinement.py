#!/usr/bin/env python3
"""Run the coupled bulk/boundary MMS on nested Q1 meshes and tabulate nodal errors."""
from pathlib import Path
import math, subprocess, csv, time, xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OOFEM=ROOT/'build/debug_magneto/oofem'
WORK=Path('/tmp/stvowu_mesh_refinement'); WORK.mkdir(exist_ok=True)
T=0.02; DT=0.0005; Lx=2.0; Ly=1.0

def exact(x,y,t=0):
    c=math.cos(math.pi*x/Lx)*math.cos(math.pi*y/Ly)+0.25*math.sin(math.pi*x/Lx)*math.sin(math.pi*y/Ly)
    return 0.2+0.1*math.exp(-t)*c

def build(nx,ny):
    nodes=[(j*(nx+1)+i+1, Lx*i/nx, Ly*j/ny) for j in range(ny+1) for i in range(nx+1)]
    elems=[]
    for j in range(ny):
        for i in range(nx):
            n=j*(nx+1)+i+1
            elems.append(('quad1ht',[n,n+1,n+nx+2,n+nx+1]))
    # Boundary line elements use the same scalar field and material.
    elems += [('line1ht',[i+1,i+2]) for i in range(nx)]
    elems += [('line1ht',[j*(nx+1)+nx+1,(j+1)*(nx+1)+nx+1]) for j in range(ny)]
    elems += [('line1ht',[(ny)*(nx+1)+i+1,(ny)*(nx+1)+i]) for i in range(nx,0,-1)]
    elems += [('line1ht',[j*(nx+1)+1,(j-1)*(nx+1)+1]) for j in range(ny,0,-1)]
    nnode=len(nodes); nelem=len(elems); nsteps=round(T/DT)
    out=f'mms_q1_{nx}x{ny}.out'
    lines=[out,'Stefanelli-Voso-Wu bulk plus dynamic-boundary MMS',
        f'nltransienttransportproblem nsteps {nsteps} deltat {DT} alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1',
        'vtkxml tstep_all domain_all primvars 1 6','domain heattransfer','OutputManager tstep_all dofman_all element_all',
        f'ndofman {nnode} nelem {nelem} ncrosssect 1 nmat 1 nbc 0 nic {nnode} nltf 1 nset {nnode+1}']
    lines += [f'node {k} coords 3 {x:.15g} {y:.15g} 0.0' for k,x,y in nodes]
    # Integrate Q1 Allen-Cahn cubic reactions and quartic potential with a
    # 3x3 rule; the element default 2x2 rule is under-integrated for them.
    lines += [f'{typ} {k} nodes {len(conn)} '+' '.join(map(str,conn)) + (' nip 9' if typ == 'quad1ht' else '') for k,(typ,conn) in enumerate(elems,1)]
    lines += ['SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
        'allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 manufactured 1 mmsbase 0.2 mmsamplitude 0.1 mmsrate 1.0 mmsbeta 0.25 mmslx 2.0 mmsly 1.0']
    for k,x,y in nodes:
        lines.append(f'InitialCondition {k} dofs 1 10 Conditions 1 u {exact(x,y):.16g} set {k+1}')
    lines += ['ConstantFunction 1 f(t) 1.0',f'Set 1 elementranges {{(1 {nelem})}}']
    lines += [f'Set {k+1} nodes 1 {k}' for k,_,_ in nodes]
    inp=WORK/f'mms_{nx}x{ny}.in'; inp.write_text('\n'.join(lines)+'\n')
    start=time.perf_counter()
    p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=WORK,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    runtime=time.perf_counter()-start
    (WORK/f'mms_{nx}x{ny}.log').write_text(p.stdout)
    if p.returncode or 'Error' in p.stdout or 'Warning' in p.stdout:
        raise RuntimeError(f'{inp.name} failed; see log')
    f=WORK/f'{out}.m0.{nsteps}.vtu'; root=ET.parse(f).getroot(); piece=root.find('.//Piece')
    xyz=np.array(list(map(float,piece.find('./Points/DataArray').text.split()))).reshape(-1,3)
    vals=np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))
    err=np.array([v-exact(x,y,T) for (x,y,_),v in zip(xyz,vals)])
    return nx,ny,len(vals),math.sqrt(float(err@err)/len(err)),float(np.max(np.abs(err))),runtime

rows=[build(nx,ny) for nx,ny in [(8,4),(16,8),(32,16)]]
orders=[]
for a,b in zip(rows,rows[1:]): orders.append(math.log(a[3]/b[3],2))
out=ROOT/'doc/stvowu_validation/spatial_convergence.csv'
with out.open('w') as f:
    w=csv.writer(f); w.writerow(['nx','ny','nodes','nodal_rms','max_abs','runtime_seconds','observed_order_from_previous']);
    for i,row in enumerate(rows): w.writerow([*row, '' if i==0 else orders[i-1]])
for i,row in enumerate(rows): print(row,'p=',None if i==0 else orders[i-1])
