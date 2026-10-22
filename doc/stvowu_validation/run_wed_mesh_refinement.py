#!/usr/bin/env python3
"""Compare WED meshes at fixed epsilon against a fine WED trajectory."""
from pathlib import Path
import csv,re,subprocess,time,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[2]; OOFEM=ROOT/'build/debug_magneto/oofem'
WED=Path(__file__).with_name('wed_allencahn.py'); MESHES=Path('/tmp/stvowu_mesh_refinement'); OUT=Path('/tmp/stvowu_wed_mesh'); DT=.001; T=.1

def readvtu(f):
 r=ET.parse(f).getroot(); p=r.find('.//Piece'); return np.array(list(map(float,p.find('./PointData/DataArray').text.split())))

def main():
 import sys
 sys.path.insert(0,str(Path(__file__).parent)); from run_energy_balance import geometry,assembled
 OUT.mkdir(exist_ok=True); cases={}
 for tag in ('8x4','16x8','32x16'):
  source=MESHES/f'mms_{tag}.in'; lines=source.read_text().splitlines(); nx,ny=map(int,tag.split('x'))
  lines[0]=f'direct_{tag}.out'; lines[2]=f'nltransienttransportproblem nsteps {round(T/DT)} deltat {DT} alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1'
  lines=[re.sub(r'manufactured 1','manufactured 0',s) for s in lines]
  folder=OUT/tag;folder.mkdir(exist_ok=True);inp=folder/'direct.in';inp.write_text('\n'.join(lines)+'\n');cases[tag]=(folder,inp,nx,ny)
  if tag=='32x16':
   p=subprocess.run([str(OOFEM),'-f',str(inp)],cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(folder/'direct.log').write_text(p.stdout)
   if p.returncode or '0 error(s) and 0 warning(s)' not in p.stdout:raise RuntimeError('fine direct OOFEM reference failed')
 fine_folder,fine_inp,_,_=cases['32x16'];fine_nodes,_,fine_u0=geometry(fine_inp);fine_ref=[fine_u0]
 for i in range(1,round(T/DT)+1):fine_ref.append(readvtu(fine_folder/f'direct_32x16.out.m0.{i}.vtu'))
 fine_ref=np.asarray(fine_ref);times=np.arange(len(fine_ref))*DT;fine_ids={tuple(np.round(x,12)):j for j,x in enumerate(fine_nodes.values())}
 rows=[]; wed_paths={}
 for tag in ('8x4','16x8','32x16'):
  folder,inp,nx,ny=cases[tag]; nodes,elems,_=geometry(inp);Mb,Mg=assembled(nodes,elems)
  zpath=folder/'wed.npz';cmd=['python3',str(WED),str(inp),'--eps','0.05','--Smax','12','--ds','0.1','--output',str(zpath),'--maxiter','100']
  start=time.perf_counter();p=subprocess.run(cmd,cwd=folder,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);runtime=time.perf_counter()-start;(folder/'wed.log').write_text(p.stdout)
  if p.returncode:raise RuntimeError(f'WED mesh {tag} failed')
  z=np.load(zpath);path=np.column_stack([np.interp(times,z['time'],z['solution'][:,j]) for j in range(z['solution'].shape[1])]);
  wed_paths[tag]=(path,int(z['iterations']),runtime,float(z['gradient_inf']))
 fine_wed=wed_paths['32x16'][0]
 for tag in ('8x4','16x8','32x16'):
  folder,inp,nx,ny=cases[tag]; nodes,elems,_=geometry(inp);Mb,Mg=assembled(nodes,elems);wed,nit,runtime,gn=wed_paths[tag]
  # A mesh study must keep epsilon fixed: compare WED with WED. The separate
  # OOFEM comparison contains the finite-epsilon bias and is reported apart.
  ix=np.array([fine_ids[tuple(np.round(x,12))] for x in nodes.values()]);dw=wed-fine_wed[:,ix]
  ew_b=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',dw,Mb,dw),times)));ew_g=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',dw,Mg,dw),times)))
  do=wed-fine_ref[:,ix]
  eo_b=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',do,Mb,do),times)));eo_g=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',do,Mg,do),times)))
  h=2.0/nx; rows.append((tag,h,len(nodes),ew_b,ew_g,eo_b,eo_g,nit,runtime,gn,'',''))
  print(f'h={h:g} ({tag}, {len(nodes)} nodes): WED/WED Ebulk={ew_b:.8g}, Egamma={ew_g:.8g}; WED/OOFEM Ebulk={eo_b:.8g}, Egamma={eo_g:.8g}; Newton={nit}, runtime={runtime:.2f}s')
 pbulk=np.log(rows[0][3]/rows[1][3])/np.log(2.0); pgamma=np.log(rows[0][4]/rows[1][4])/np.log(2.0)
 rows[0]=(*rows[0][:10],pbulk,pgamma); rows[1]=(*rows[1][:10],pbulk,pgamma)
 with (Path(__file__).parent/'wed_mesh_refinement.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['mesh','h','nodes','bulk_L2_space_time_vs_fine_WED_same_epsilon','boundary_L2_space_time_vs_fine_WED_same_epsilon','bulk_L2_space_time_vs_fine_OOFEM_includes_epsilon_bias','boundary_L2_space_time_vs_fine_OOFEM_includes_epsilon_bias','newton_iterations','runtime_seconds','gradient_inf','observed_order_bulk','observed_order_boundary']);w.writerows(rows)
if __name__=='__main__':main()
