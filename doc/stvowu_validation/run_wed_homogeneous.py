#!/usr/bin/env python3
"""Reproduce WED epsilon-limit check on the spatially homogeneous invariant solution."""
from pathlib import Path
import argparse, math, subprocess, csv, time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OOFEM=ROOT/'build/debug_magneto/oofem'
WED=Path(__file__).with_name('wed_allencahn.py')
INPUT=Path(__file__).with_name('wed_homogeneous.in')

def analytic(t,u0=.3):
    y0=u0*u0
    return math.sqrt(1.0/(1.0+(1.0/y0-1.0)*math.exp(-2*t)))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output-dir',default='/tmp/stvowu_wed_homogeneous'); a=p.parse_args()
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    direct=out/'reference.in'; lines=INPUT.read_text().splitlines(); lines[0]='homogeneous_direct.out'; direct.write_text('\n'.join(lines)+'\n')
    run=subprocess.run([str(OOFEM),'-f',str(direct)],cwd=out,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'oofem.log').write_text(run.stdout)
    if run.returncode or '0 error(s) and 0 warning(s)' not in run.stdout: raise RuntimeError('OOFEM reference transient failed; see oofem.log')
    rows=[]
    integrated=[]
    import xml.etree.ElementTree as ET
    ref_times=np.arange(0,1001)*0.001
    ref_means=[0.3]
    for nstep in range(1,201):
        root=ET.parse(out/f'homogeneous_direct.out.m0.{nstep}.vtu').getroot()
        ref_means.append(float(np.mean(np.fromstring(root.find('.//PointData/DataArray').text,sep=' '))))
    ref_times=ref_times[:len(ref_means)]
    for eps in (.1,.05,.025):
        label=f'{eps:g}'.replace('.','p')
        archive=out/f'wed_eps_{label}.npz'
        cmd=['python3',str(WED),str(direct),'--eps',str(eps),'--Smax','12','--ds','0.1','--output',str(archive),'--maxiter','100']
        start=time.perf_counter(); proc=subprocess.run(cmd,cwd=out,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); runtime=time.perf_counter()-start
        (out/f'wed_eps_{label}.log').write_text(proc.stdout)
        if proc.returncode: raise RuntimeError(f'WED failed for epsilon={eps}; see its log')
        z=np.load(archive); tvals=z['time']; solution=z['solution']
        mask=tvals <= .2 + 1.e-12; tgrid=tvals[mask]; exact_grid=np.array([analytic(float(t)) for t in tgrid]); wed_grid=solution[mask].mean(axis=1); ref_grid=np.interp(tgrid,ref_times,ref_means)
        err_ref=wed_grid-ref_grid; err_exact=wed_grid-exact_grid
        ebulk=math.sqrt(2.0*float(np.trapz(err_ref**2,tgrid))); egamma=math.sqrt(6.0*float(np.trapz(err_ref**2,tgrid)))
        ebulk_exact=math.sqrt(2.0*float(np.trapz(err_exact**2,tgrid))); egamma_exact=math.sqrt(6.0*float(np.trapz(err_exact**2,tgrid)))
        integrated.append((eps,ebulk,egamma,ebulk_exact,egamma_exact,int(z['iterations']),runtime,float(z['gradient_inf'])))
        for t in (.05,.1,.2):
            j=round(t/(tvals[1]-tvals[0])); value=float(solution[j].mean()); exact=analytic(t)
            nstep=round(t/.001); root=__import__('xml.etree.ElementTree',fromlist=['ElementTree']).parse(out/f'homogeneous_direct.out.m0.{nstep}.vtu').getroot()
            direct_values=np.array(list(map(float,root.find('.//PointData/DataArray').text.split())))
            rows.append((eps,t,value,exact,abs(value-exact),float(direct_values.mean()),float(z['gradient_inf']),float(z['gradient_check_relative']),float(z['hessian_check_relative']),bool(z['success'])))
    with (out/'wed_integrated.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['epsilon','bulk_L2_vs_OOFEM','boundary_L2_vs_OOFEM','bulk_L2_vs_analytic','boundary_L2_vs_analytic','newton_iterations','runtime_seconds','gradient_inf']); w.writerows(integrated)
    with (out/'wed_homogeneous.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['epsilon','time','wed_value','analytic_value','absolute_error','oofem_mean','newton_gradient_inf','gradient_check_relative','hessian_check_relative','optimizer_success']); w.writerows(rows)
    print('Integrated space-time errors on t=[0,0.2], |Omega|=2, |Gamma|=6')
    for row in integrated: print(f'epsilon={row[0]:g}, Ebulk={row[1]:.7g}, Egamma={row[2]:.7g}, Newton={row[5]}, seconds={row[6]:.3f}')
    print('epsilon,time,WED error vs analytic, OOFEM error vs analytic')
    for row in rows:
        eps,t,val,exact,err,direct_value,*_=row
        print(f'{eps:g},{t:g},{err:.9g},{abs(direct_value-exact):.9g}')

if __name__=='__main__': main()
