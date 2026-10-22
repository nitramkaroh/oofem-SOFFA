#!/usr/bin/env python3
"""Compare WED horizon truncations for epsilon=.05 on a common physical window."""
from pathlib import Path
import csv, re, subprocess
import numpy as np
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]; OOFEM=ROOT/'build/debug_magneto/oofem'
WED=Path(__file__).with_name('wed_allencahn.py'); BASE=Path(__file__).with_name('allencahn_rectangle_dynamic.in')
OUT=Path('/tmp/stvowu_wed_horizon'); EPS=.05; DS=.1; TCOMPARE=.1

def main():
    OUT.mkdir(exist_ok=True); lines=BASE.read_text().splitlines(); lines[0]='horizon_reference.out'
    lines[2]='nltransienttransportproblem nsteps 200 deltat 0.001 alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1'
    lines=[re.sub(r'manufactured 1','manufactured 0',s) for s in lines]
    inp=OUT/'direct.in'; inp.write_text('\n'.join(lines)+'\n')
    p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/'direct.log').write_text(p.stdout)
    if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout: raise RuntimeError('direct reference failed')
    paths={}
    for smax in (6.,12.,18.):
        archive=OUT/f'wed_S{smax:g}.npz'; cmd=['python3',str(WED),str(inp),'--eps',str(EPS),'--Smax',str(smax),'--ds',str(DS),'--output',str(archive),'--maxiter','100']
        p=subprocess.run(cmd,cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); (OUT/f'wed_S{smax:g}.log').write_text(p.stdout)
        if p.returncode: raise RuntimeError(f'WED Smax={smax:g} failed')
        z=np.load(archive); times=np.arange(0,TCOMPARE+1.e-12,.001)
        U=np.column_stack([np.interp(times,z['time'],z['solution'][:,j]) for j in range(z['solution'].shape[1])]); paths[smax]=(times,U,z)
    ref=paths[18.][1]; rows=[]
    for smax in (6.,12.,18.):
        times,U,z=paths[smax]; diff=U-ref
        rows.append((smax,EPS,DS,EPS*DS,TCOMPARE,float(np.sqrt(np.trapz(np.mean(diff*diff,axis=1),times))),float(np.max(np.abs(diff))),int(z['iterations']),float(z['gradient_inf'])))
    with (Path(__file__).with_name('wed_horizon_sensitivity.csv')).open('w') as f:
        w=csv.writer(f); w.writerow(['Smax','epsilon','ds','physical_dt','comparison_window','space_time_RMS_vs_S18','max_nodal_difference_vs_S18','Newton_iterations','gradient_inf']); w.writerows(rows)
    for r in rows: print(f'Smax={r[0]:g}: RMS={r[5]:.6e}, max={r[6]:.6e}, Newton={r[7]}')
if __name__=='__main__': main()
