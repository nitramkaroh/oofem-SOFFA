#!/usr/bin/env python3
"""Compare WED paths with direct OOFEM for a nonuniform 2D initial field."""
from pathlib import Path
import csv,re,subprocess,time,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[2]; OOFEM=ROOT/'build/debug_magneto/oofem'
WED=Path(__file__).with_name('wed_allencahn.py'); BASE=Path(__file__).with_name('allencahn_rectangle_dynamic.in')
OUT=Path('/tmp/stvowu_wed_nonuniform'); T=.1; DT=.001

def readvtu(f):
 root=ET.parse(f).getroot(); p=root.find('.//Piece')
 return np.array(list(map(float,p.find('./PointData/DataArray').text.split())))

def main():
 OUT.mkdir(exist_ok=True); lines=BASE.read_text().splitlines(); lines[0]='direct_nonuniform.out'
 lines[2]=f'nltransienttransportproblem nsteps {round(.2/DT)} deltat {DT} alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1'
 lines=[re.sub(r'manufactured 1','manufactured 0',s) for s in lines]
 direct=OUT/'direct.in'; direct.write_text('\n'.join(lines)+'\n')
 start=time.perf_counter(); p=subprocess.run([str(OOFEM),'-f',str(direct)],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); direct_runtime=time.perf_counter()-start
 (OUT/'direct.log').write_text(p.stdout)
 if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout: raise RuntimeError('direct OOFEM reference failed')
 import sys
 sys.path.insert(0,str(Path(__file__).parent)); from run_energy_balance import geometry,assembled
 nodes,elems,u0=geometry(direct); Mb,Mg=assembled(nodes,elems)
 ref=[u0]
 for i in range(1,round(T/DT)+1): ref.append(readvtu(OUT/f'direct_nonuniform.out.m0.{i}.vtu'))
 ref=np.asarray(ref); ref_t=np.arange(ref.shape[0])*DT; target=ref_t
 rows=[]
 for eps in (.1,.05,.025,.0125):
  label=str(eps).replace('.','p'); archive=OUT/f'wed_eps_{label}.npz'
  cmd=['python3',str(WED),str(direct),'--eps',str(eps),'--Smax','12','--ds','0.1','--output',str(archive),'--maxiter','100']
  start=time.perf_counter(); p=subprocess.run(cmd,cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); runtime=time.perf_counter()-start
  (OUT/f'wed_{label}.log').write_text(p.stdout)
  if p.returncode: raise RuntimeError(f'WED epsilon={eps} failed')
  z=np.load(archive); wt=z['time']; wu=z['solution']; interp=np.column_stack([np.interp(target,wt,wu[:,j]) for j in range(wu.shape[1])]); diff=interp-ref
  qb=np.einsum('ti,ij,tj->t',diff,Mb,diff); qg=np.einsum('ti,ij,tj->t',diff,Mg,diff)
  eb=float(np.sqrt(np.trapz(qb,target))); eg=float(np.sqrt(np.trapz(qg,target)))
  rows.append((eps,eb,eg,int(z['iterations']),runtime,float(z['gradient_inf']),float(z['gradient_check_relative']),float(z['hessian_check_relative']),direct_runtime))
  print(f'epsilon={eps:g}: Ebulk={eb:.8g}, Egamma={eg:.8g}, Newton={int(z["iterations"])}, runtime={runtime:.3f}s, grad={float(z["gradient_inf"]):.2e}')
 time_rows=[]
 for ds in (.05,.025):
  archive=OUT/f'wed_eps_0p05_ds_{str(ds).replace(".","p")}.npz'
  cmd=['python3',str(WED),str(direct),'--eps','0.05','--Smax','12','--ds',str(ds),'--output',str(archive),'--maxiter','100']
  start=time.perf_counter(); p=subprocess.run(cmd,cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT); runtime=time.perf_counter()-start
  (OUT/f'wed_ds_{str(ds).replace(".","p")}.log').write_text(p.stdout)
  if p.returncode: raise RuntimeError(f'WED temporal refinement ds={ds} failed')
  z=np.load(archive); wt=z['time']; wu=z['solution']; interp=np.column_stack([np.interp(target,wt,wu[:,j]) for j in range(wu.shape[1])]); diff=interp-ref
  qb=np.einsum('ti,ij,tj->t',diff,Mb,diff); qg=np.einsum('ti,ij,tj->t',diff,Mg,diff)
  eb=float(np.sqrt(np.trapz(qb,target))); eg=float(np.sqrt(np.trapz(qg,target)))
  time_rows.append((0.05,ds,0.05*ds,eb,eg,int(z['iterations']),runtime,float(z['gradient_inf'])))
  print(f'epsilon=0.05 ds={ds:g}: Ebulk={eb:.8g}, Egamma={eg:.8g}, Newton={int(z["iterations"])}, runtime={runtime:.3f}s')
 with (Path(__file__).parent/'wed_nonuniform.csv').open('w') as f:
  w=csv.writer(f); w.writerow(['epsilon','bulk_L2_space_time','boundary_L2_space_time','newton_iterations','runtime_seconds','gradient_inf','gradient_check_relative','hessian_check_relative','oofem_reference_runtime_seconds']); w.writerows(rows)
 with (Path(__file__).parent/'wed_nonuniform_time_refinement.csv').open('w') as f:
  w=csv.writer(f); w.writerow(['epsilon','ds','physical_dt','bulk_L2_space_time','boundary_L2_space_time','newton_iterations','runtime_seconds','gradient_inf']); w.writerows(time_rows)
if __name__=='__main__': main()
