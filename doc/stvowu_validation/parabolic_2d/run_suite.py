#!/usr/bin/env python3
"""Reproducible mesh/time and boundary orientation runs for parabolic 2D AC."""
import subprocess, sys
from pathlib import Path
from run_case import HERE, make_input, OOFEM

ROOT_RUNS = HERE / 'runs'
CASES = [
    # Spatial sequence with identical physical parameters and time step.
    (20,.0005,.1,200,.1,'space_n20_dt0p0005',False,False),
    (30,.0005,.1,200,.1,'space_n30_dt0p0005',False,False),
    (40,.0005,.1,200,.1,'space_n40_dt0p0005',False,False),
    # A finer space sequence at a shared final time and shared time step.
    (20,.000625,.0625,100,.1,'space_T0p0625_n20_dt0p000625',False,False),
    (30,.000625,.0625,100,.1,'space_T0p0625_n30_dt0p000625',False,False),
    (40,.000625,.0625,100,.1,'space_T0p0625_n40_dt0p000625',False,False),
    (60,.000625,.0625,100,.1,'space_T0p0625_n60_dt0p000625',False,False),
    # Time sequence at fixed n=20; dt=1e-4 is also run with output every step.
    (20,.0004,.1,250,.1,'time_n20_dt0p0004',False,False),
    (20,.0002,.1,500,.1,'time_n20_dt0p0002',False,False),
    (20,.0001,.1,1,.1,'energy_n20_dt0p0001_every_step',False,False),
    # Small natural case isolates orientation/order invariance without a long run.
    (8,.001,.01,10,.1,'orientation_base_n8_dt0p001',False,False),
    (8,.001,.01,10,.1,'orientation_reverse_n8_dt0p001',True,True),
]
for n,dt,end_time,out_every,eps,name,reverse,permute in CASES:
    d=(ROOT_RUNS/name).resolve()
    steps=round(end_time/dt)
    make_input(n,dt,end_time,d,out_every,eps,reverse,permute)
    # Parameterize the output cadence to keep only snapshots for the spatial
    # trajectory and final state for the temporal/orientation controls.
    if name.startswith('time_') or name.startswith('orientation_'):
        lines=(d/'case.in').read_text().replace(f'vtkxml tstep_step {out_every}',f'vtkxml tstep_step {steps}')
        lines=lines.replace(f'OutputManager tstep_step {out_every}',f'OutputManager tstep_step {steps}')
        (d/'case.in').write_text(lines)
    final=d/f'case.out.m0.{steps}.vtu'
    oldlog=d/'run.log'
    if final.exists() and oldlog.exists() and 'ANALYSIS FINISHED' in oldlog.read_text():
        print(f'SKIP completed {name}',flush=True)
        continue
    print(f'RUN {name}: n={n}, dt={dt:g}, steps={steps}',flush=True)
    result=subprocess.run([str(OOFEM),'-f',str(d/'case.in')],cwd=d,text=True,
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (d/'run.log').write_text(result.stdout)
    print(result.stdout[-1000:],flush=True)
    if result.returncode or '0 error(s) and 0 warning(s)' not in result.stdout:
        raise SystemExit(f'OOFEM failed for {name}; inspect {d}/run.log')
print('ALL RUNS COMPLETED')
