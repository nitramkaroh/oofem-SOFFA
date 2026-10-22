#!/usr/bin/env python3
"""Refine WED s-step for epsilon=0.0125 and compare on t in [0,0.1]."""
from pathlib import Path
import csv, subprocess, time, xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent
WED=HERE/'wed_allencahn.py'
OUT=Path('/tmp/stvowu_wed_robustness')
INPUT=OUT/'direct.in'
OLD=OUT/'epsilon_0.0125/path.npz'
def vtu(path):
    root=ET.parse(path).getroot()
    return np.fromstring(root.find('.//Piece/PointData/DataArray').text,sep=' ')
def main():
    import sys
    sys.path.insert(0,str(HERE))
    from run_energy_balance import geometry,assembled
    folder=OUT/'epsilon_0.0125_ds_0.2';folder.mkdir(exist_ok=True)
    output=folder/'path.npz';log=folder/'wed.log'
    cmd=['python3',str(WED),str(INPUT),'--eps','0.0125','--Smax','12','--ds','0.2','--tol','2e-5','--maxiter','100','--output',str(output)]
    start=time.perf_counter();p=subprocess.run(cmd,cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);runtime=time.perf_counter()-start;log.write_text(p.stdout)
    if p.returncode:raise RuntimeError(f'WED failed (exit {p.returncode}): {p.stdout}')
    z=np.load(output);old=np.load(OLD)
    nodes,elems,u0=geometry(INPUT);Mb,Mg=assembled(nodes,elems)
    t=np.arange(101)*.001
    ref=np.asarray([u0]+[vtu(OUT/f'direct_reference.out.m0.{i}.vtu') for i in range(1,101)])
    def interpolate(path,times):
        return np.column_stack([np.interp(t,path['time'],path['solution'][:,j]) for j in range(path['solution'].shape[1])])
    u=interpolate(z,t);uold=interpolate(old,t)
    def errs(a):
        d=a-ref
        eb=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',d,Mb,d),t)))
        eg=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',d,Mg,d),t)))
        return eb,eg
    eb,eg=errs(u);oldb,oldg=errs(uold)
    du=u-uold
    db=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',du,Mb,du),t)))
    dg=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',du,Mg,du),t)))
    rows=[('ds=0.4',.4,.005,oldb,oldg,0.,0.,int(old['iterations']),float(old['gradient_inf'])),
          ('ds=0.2',.2,.0025,eb,eg,db,dg,int(z['iterations']),float(z['gradient_inf']))]
    with (HERE/'wed_small_epsilon_time_check.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['run','ds','physical_dt','bulk_L2_space_time_error_vs_OOFEM','boundary_L2_space_time_error_vs_OOFEM','bulk_L2_space_time_difference_vs_ds0p4','boundary_L2_space_time_difference_vs_ds0p4','newton_iterations','gradient_inf']);w.writerows(rows)
    print(p.stdout,end='')
    print(f'error OOFEM: bulk {oldb:.9g} -> {eb:.9g}; boundary {oldg:.9g} -> {eg:.9g}')
    print(f'difference from ds=.4 path: bulk {db:.9g}; boundary {dg:.9g}; runtime={runtime:.3f}s')
if __name__=='__main__':main()
