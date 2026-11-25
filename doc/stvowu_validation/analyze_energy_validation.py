#!/usr/bin/env python3
"""Recompute energy histories from high-precision OOFEM text outputs.

No OOFEM solve is launched. The full precision dof values in the text output
are used so the reconstructed Newmark velocity is not contaminated by VTU's
shorter ASCII formatting.
"""
import csv, importlib.util, re
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('energy_runner',HERE/'run_energy_validation.py')
runner=importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)

CASES=[
    ('parabolic',16,8,.002,'energy_validation_runs_highprecision'),
    ('parabolic',16,8,.001,'energy_validation_runs_highprecision'),
    ('parabolic',16,8,.0005,'energy_validation_runs_highprecision'),
    ('parabolic',8,4,.0005,'energy_validation_runs_highprecision'),
    ('parabolic',32,16,.0005,'energy_validation_runs'),
    ('hyperbolic',16,8,.002,'energy_validation_runs_highprecision'),
    ('hyperbolic',16,8,.001,'energy_validation_runs_highprecision'),
    ('hyperbolic',16,8,.0005,'energy_validation_runs_highprecision'),
    ('hyperbolic',8,4,.0005,'energy_validation_runs_highprecision'),
    ('hyperbolic',32,16,.0005,'energy_validation_runs_highprecision'),
]

def evaluate(kind,nx,ny,dt,folder):
    tag=f'{kind}_n{nx}x{ny}_dt{dt:g}'.replace('.','p')
    run=HERE/folder/tag
    txt=run/f'{tag}.out'
    if not txt.exists(): raise FileNotFoundError(txt)
    xy,quads,edges,boundary=runner.mesh(nx,ny)
    nsteps=round(runner.T/dt)
    hist=runner.read_text_history(txt,nsteps,len(xy),dt)
    U=np.vstack([runner.fields(xy)]+[item[1] for item in hist])
    V=np.zeros_like(U)
    if kind=='hyperbolic':
        for j in range(nsteps): V[j+1]=2*(U[j+1]-U[j])/dt-V[j]
    else:
        V[1:]=(U[1:]-U[:-1])/dt
    Mb,Ms=runner.matrices(xy,quads,edges); Db=Ds=0.; E0=None; rows=[]
    for j,u in enumerate(U):
        bg,bp,kin,sg,sp=runner.energy_parts(u,V[j],xy,quads,edges)
        if kind=='parabolic': kin=0.
        Et=bg+bp+kin+sg+sp
        if E0 is None: E0=Et
        if j:
            if kind=='parabolic':
                du=U[j]-U[j-1]
                db=runner.C*(du@Mb@du)/dt
                ds=runner.CS*(du@Ms@du)/dt
            else:
                v0,v1=V[j-1],V[j]
                db=runner.RHO*runner.DELTA*dt/3*(v0@Mb@v0+v0@Mb@v1+v1@Mb@v1)
                ds=runner.CS*dt/3*(v0@Ms@v0+v0@Ms@v1+v1@Ms@v1)
            Db+=db; Ds+=ds
        residual=Et-E0+Db+Ds
        rows.append({'kind':kind,'mesh':f'{nx}x{ny}','dt':dt,'time':j*dt,
            'bulk_gradient':bg,'bulk_potential':bp,'bulk_kinetic':kin,
            'surface_gradient':sg,'surface_potential':sp,'energy_total':Et,
            'dissipation_bulk_cumulative':Db,'dissipation_surface_cumulative':Ds,
            'balance_residual':residual,'normalized_balance_residual':residual/E0})
    return rows

allrows=[]; summary=[]
for case in CASES:
    rows=evaluate(*case); allrows.extend(rows); end=rows[-1]
    summary.append({'kind':case[0],'mesh':f'{case[1]}x{case[2]}','dt':case[3],'T':runner.T,
       'E0':rows[0]['energy_total'],'ET':end['energy_total'],
       'bulk_dissipation':end['dissipation_bulk_cumulative'],
       'surface_dissipation':end['dissipation_surface_cumulative'],
       'balance_residual':end['balance_residual'],
       'relative_balance_residual':end['normalized_balance_residual'],
       'max_abs_residual':max(abs(r['balance_residual']) for r in rows),
       'positive_energy_steps':sum(rows[i]['energy_total']>rows[i-1]['energy_total']+1e-12 for i in range(1,len(rows)))})

for name,rows in [('energy_validation_hp_time_series.csv',allrows),('energy_validation_hp_summary.csv',summary)]:
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
for r in summary: print(r)
