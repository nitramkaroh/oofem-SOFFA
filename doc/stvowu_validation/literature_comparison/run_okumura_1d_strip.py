#!/usr/bin/env python3
"""Reproduce Okumura dissertation Ch. 4, Computation Example 2, in OOFEM.

The 1D interval is extruded to a unit-depth 2D strip. Only the two x-end edges
are assigned dynamic-boundary elements; top and bottom remain natural Neumann
edges. A y-independent Q1 field then reduces to the 1D bulk/surface problem.
Separate bulk and end-edge materials set the published boundary capacity 10.
"""
from pathlib import Path
import math, subprocess

ROOT=Path(__file__).resolve().parents[3]
OUT=Path('/tmp/stvowu_okumura_ch4_case2')
OUT.mkdir(exist_ok=True)
NX=100; NY=1; DT=0.001; T=0.7; NSTEPS=round(T/DT)
NN=(NX+1)*(NY+1); NE=NX*NY+2
lines=[
    'okumura_ch4_case2.out','Okumura thesis Chapter 4, Computation Example 2',
    f'nltransienttransportproblem nsteps {NSTEPS} deltat {DT} alpha 1.0 rtol 1.e-10 nsmax 30 nmodules 1',
    'vtkxml tstep_all domain_all primvars 1 6','domain heattransfer',
    'OutputManager tstep_all dofman_all element_all',
    f'ndofman {NN} nelem {NE} ncrosssect 2 nmat 2 nbc 0 nic {NX+1} nltf 1 nset {NX+3}'
]
for j in range(NY+1):
    for i in range(NX+1):
        node=j*(NX+1)+i+1
        lines.append(f'node {node} coords 3 {i/NX:.12g} {j/NY:.12g} 0.0')
eid=1
for j in range(NY):
    for i in range(NX):
        a=j*(NX+1)+i+1; b=a+1; d=a+NX+1; c=d+1
        lines.append(f'quad1ht {eid} nodes 4 {a} {b} {c} {d} nip 9'); eid+=1
# Dynamic 0D endpoints of the 1D interval, represented as unit-height edges.
lines.append(f'line1ht {eid} nodes 2 1 {NX+2} nip 3'); eid+=1
lines.append(f'line1ht {eid} nodes 2 {NX+1} {NN} nip 3')
lines += [
    'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
    'SimpleTransportCS 2 thickness 1.0 area 1.0 mat 2 set 2',
    'allencahn 1 d 1.0 k 1.0 c 1.0 mobility 100.0',
    'allencahn 2 d 1.0 k 1.0 c 10.0 mobility 100.0 surfacebeta 1.0'
]
for i in range(NX+1):
    value=math.exp(-500.0*(i/NX-.5)**2)
    n0=i+1; n1=NX+1+i+1
    lines.append(f'InitialCondition {i+1} dofs 1 10 Conditions 1 u {value:.16g} set {i+3}')
lines += [
    'ConstantFunction 1 f(t) 1.0',
    f'Set 1 elementranges {{(1 {NX})}}',
    f'Set 2 elements 2 {NX+1} {NX+2}'
]
for i in range(NX+1):
    n0=i+1; n1=NX+1+i+1
    lines.append(f'Set {i+3} nodes 2 {n0} {n1}')
inp=OUT/'case2.in';inp.write_text('\n'.join(lines)+'\n')
exe=ROOT/'build/debug_magneto/oofem'
p=subprocess.run([str(exe),'-f',str(inp)],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(OUT/'run.log').write_text(p.stdout)
print(p.stdout[-4000:])
if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout:
    raise SystemExit(f'OOFEM failed; inspect {OUT}/run.log')
print(f'Completed Okumura Ch.4 example 2: Nx={NX}, dt={DT}, T={T}; outputs {OUT}')
