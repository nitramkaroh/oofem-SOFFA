#!/usr/bin/env python3
"""Qualitative reproduction of Xiao & Chen (2025), Case 4, unpenalized branch.

Their paper does not state epsilon for this case. We therefore label epsilon=0.1
as an explicit exploratory assumption and do not treat this run as a parameter-
matched benchmark.
"""
from pathlib import Path
import math, subprocess, argparse

ROOT = Path(__file__).resolve().parents[3]
parser=argparse.ArgumentParser()
parser.add_argument('--n',type=int,default=20)
parser.add_argument('--dt',type=float,default=0.0001)
parser.add_argument('--eps',type=float,default=0.1)
parser.add_argument('--end-time',type=float,default=0.1)
args=parser.parse_args()
N=args.n  # exploratory mesh; paper uses h=1/300
EPS=args.eps  # undocumented in the paper; qualitative assumption only
DT=args.dt
T=args.end_time
NSTEPS = round(T / DT)
OUT = Path(f'/tmp/stvowu_xiao_chen_case4_n{N}_dt{DT:g}_eps{EPS:g}')
OUT.mkdir(exist_ok=True)
nodes = (N + 1) ** 2
elems = N * N + 4 * N
lines = [
    'xiao_chen_case4.out', 'Xiao-Chen Case 4 unpenalized, exploratory epsilon=0.1',
    f'nltransienttransportproblem nsteps {NSTEPS} deltat {DT} alpha 1.0 rtol 1.e-9 nsmax 30 nmodules 1',
    'vtkxml tstep_all domain_all primvars 1 6', 'domain heattransfer',
    'OutputManager tstep_all dofman_all element_all',
    f'ndofman {nodes} nelem {elems} ncrosssect 1 nmat 1 nbc 0 nic {nodes} nltf 1 nset {nodes+1}'
]
for j in range(N+1):
    for i in range(N+1):
        k=j*(N+1)+i+1
        lines.append(f'node {k} coords 3 {i/N:.12g} {j/N:.12g} 0')
eid=1
for j in range(N):
    for i in range(N):
        a=j*(N+1)+i+1; b=a+1; d=a+N+1; c=d+1
        lines.append(f'quad1ht {eid} nodes 4 {a} {b} {c} {d} nip 9'); eid+=1
for i in range(N): lines.append(f'line1ht {eid} nodes 2 {i+1} {i+2} nip 3'); eid+=1
for j in range(N): lines.append(f'line1ht {eid} nodes 2 {j*(N+1)+N+1} {(j+1)*(N+1)+N+1} nip 3'); eid+=1
for i in range(N,0,-1): lines.append(f'line1ht {eid} nodes 2 {N*(N+1)+i+1} {N*(N+1)+i} nip 3'); eid+=1
for j in range(N,0,-1): lines.append(f'line1ht {eid} nodes 2 {j*(N+1)+1} {(j-1)*(N+1)+1} nip 3'); eid+=1
lines += [
    'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
    f'allencahn 1 d 1.0 k 0.01 c 1.0 mobility {0.01/EPS**2:.12g} surfacebeta 1.0'
]
# The published figure depicts +1 in the square and -1 outside; square touches y=0.
for i in range(N+1):
    x=i/N
    for j in range(N+1):
        y=j/N
        value=1.0 if (0.25 <= x <= 0.75 and 0.0 <= y <= 0.5) else -1.0
        node=j*(N+1)+i+1
        lines.append(f'InitialCondition {node} dofs 1 10 Conditions 1 u {value:g} set {node+1}')
lines.append('ConstantFunction 1 f(t) 1.0')
lines.append(f'Set 1 elementranges {{(1 {elems})}}')
for node in range(1,nodes+1): lines.append(f'Set {node+1} nodes 1 {node}')
inp=OUT/'case4.in'; inp.write_text('\n'.join(lines)+'\n')
exe=ROOT/'build/debug_magneto/oofem'
run=subprocess.run([str(exe),'-f',str(inp)],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(OUT/'run.log').write_text(run.stdout)
print(run.stdout[-5000:])
if run.returncode or '0 error(s) and 0 warning(s)' not in run.stdout:
    raise SystemExit(f'OOFEM failed; inspect {OUT}/run.log')
print(f'Completed Case 4 at T={T}; EPS={EPS} is assumed, not reported by paper. Outputs: {OUT}')
