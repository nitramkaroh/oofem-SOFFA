#!/usr/bin/env python3
"""Targeted WED time-step and mesh stability checks for the validation report."""
from pathlib import Path
import csv, re, subprocess, time, xml.etree.ElementTree as ET
import os
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
WED=HERE/'wed_allencahn.py'
INPUT=HERE/'allencahn_rectangle_dynamic.in'
MESH_ROOT=Path('/tmp/stvowu_mesh_refinement')
OUT=Path('/tmp/stvowu_wed_robustness')
S=12.0
OOFEM=ROOT/'build/debug_magneto/oofem'

def run(input_file, eps, ds, tol, output, log):
    cmd=['python3',str(WED),str(input_file),'--eps',str(eps),'--Smax',str(S),
         '--ds',str(ds),'--tol',str(tol),'--maxiter','100','--output',str(output)]
    start=time.perf_counter()
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    proc=subprocess.run(cmd,cwd=output.parent,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    runtime=time.perf_counter()-start
    log.write_text(proc.stdout)
    if proc.returncode:
        raise RuntimeError(f'WED failed (exit {proc.returncode}) for {input_file}, epsilon={eps}, ds={ds}: {proc.stdout}')
    z=np.load(output)
    return z, runtime

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    # Use the same non-forced OOFEM trajectory as the epsilon study as the
    # direct reference, with a fine physical time step.
    direct_lines=INPUT.read_text().splitlines()
    direct_lines[0]='direct_reference.out'
    direct_lines[2]='nltransienttransportproblem nsteps 100 deltat 0.001 alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1'
    direct_lines=[re.sub(r'manufactured 1','manufactured 0',line) for line in direct_lines]
    direct=OUT/'direct.in';direct.write_text('\n'.join(direct_lines)+'\n')
    proc=subprocess.run([str(OOFEM),'-f',str(direct)],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/'direct.log').write_text(proc.stdout)
    if proc.returncode or '0 error(s) and 0 warning(s)' not in proc.stdout: raise RuntimeError('direct OOFEM reference failed')
    def readvtu(path):
        root=ET.parse(path).getroot();piece=root.find('.//Piece')
        return np.fromstring(piece.find('./PointData/DataArray').text,sep=' ')
    import sys
    sys.path.insert(0,str(HERE))
    from run_energy_balance import geometry, assembled
    nodes,elems,u0=geometry(direct);Mb,Mg=assembled(nodes,elems)
    target=np.arange(101)*.001
    ref=np.asarray([u0]+[readvtu(OUT/f'direct_reference.out.m0.{i}.vtu') for i in range(1,101)])
    # Keep the physical step fixed at 0.005 while varying epsilon. Smax=12
    # leaves all four trajectories longer than the common comparison window.
    eps_rows=[]
    for eps in (.1,.05,.025,.0125):
        ds=.005/eps
        folder=OUT/f'epsilon_{eps:g}';folder.mkdir(exist_ok=True)
        z,runtime=run(INPUT,eps,ds,2.e-5,folder/'path.npz',folder/'wed.log')
        interp=np.column_stack([np.interp(target,z['time'],z['solution'][:,j]) for j in range(z['solution'].shape[1])])
        diff=interp-ref
        eb=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',diff,Mb,diff),target)))
        eg=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',diff,Mg,diff),target)))
        eps_rows.append((eps,ds,eps*ds,float(z['time'][-1]),len(z['time'])-1,int(z['iterations']),runtime,float(z['gradient_inf']),float(z['objective']),eb,eg))
        print(f'epsilon={eps:g}, ds={ds:g}, dt={eps*ds:g}, Ebulk={eb:.4e}, Egamma={eg:.4e}, Newton={int(z["iterations"])}, grad={float(z["gradient_inf"]):.3e}')
    with (HERE/'wed_epsilon_fixed_dt.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['epsilon','ds','physical_dt','physical_horizon','time_intervals','newton_iterations','runtime_seconds','gradient_inf','objective','bulk_L2_space_time_vs_OOFEM','boundary_L2_space_time_vs_OOFEM']);w.writerows(eps_rows)

    # Recompute all three mesh levels with half ds and a ten-times tighter
    # scaled-gradient tolerance; compare coarse WED paths to the fine WED path.
    mesh_paths={}; mesh_meta={}
    for tag in ('8x4','16x8','32x16'):
        inp=MESH_ROOT/f'mms_{tag}.in';folder=OUT/'meshes'/tag;folder.mkdir(parents=True,exist_ok=True)
        z,runtime=run(inp,.05,.05,2.e-6,folder/'path.npz',folder/'wed.log')
        mesh_paths[tag]=z['solution'];mesh_meta[tag]=(runtime,int(z['iterations']),float(z['gradient_inf']),z['time'])
        print(f'mesh={tag}, dt=0.0025, tol=2e-6, Newton={int(z["iterations"])}, runtime={runtime:.2f}s, grad={float(z["gradient_inf"]):.3e}')
    # The mesh trajectories use common physical times and nested rectangular
    # nodes. Compare on the same physical window [0,0.1] as the baseline mesh
    # study; the fine-grid trajectory is a numerical reference, not exact.
    fine_times=mesh_meta['32x16'][3]
    tf=fine_times[fine_times<=.1+1.e-12]
    fine=np.column_stack([np.interp(tf,fine_times,mesh_paths['32x16'][:,j]) for j in range(mesh_paths['32x16'].shape[1])])
    fine_inp=MESH_ROOT/'mms_32x16.in'
    fine_nodes,_,_=geometry(fine_inp)
    fine_ids={tuple(np.round(x,12)):j for j,x in enumerate(fine_nodes.values())}
    rows=[]
    for tag in ('8x4','16x8','32x16'):
        inp=MESH_ROOT/f'mms_{tag}.in';nodes,elems,_=geometry(inp);Mb,Mg=assembled(nodes,elems)
        ids=np.array([fine_ids[tuple(np.round(x,12))] for x in nodes.values()])
        local=np.column_stack([np.interp(tf,mesh_meta[tag][3],mesh_paths[tag][:,j]) for j in range(mesh_paths[tag].shape[1])])
        diff=local-fine[:,ids]
        eb=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',diff,Mb,diff),tf)))
        eg=float(np.sqrt(np.trapz(np.einsum('ti,ij,tj->t',diff,Mg,diff),tf)))
        h=2/int(tag.split('x')[0]);runtime,nit,gn,_=mesh_meta[tag]
        rows.append((tag,h,len(nodes),eb,eg,runtime,nit,gn))
    pb=np.log(rows[0][3]/rows[1][3])/np.log(2);pg=np.log(rows[0][4]/rows[1][4])/np.log(2)
    rows[0]=(*rows[0],pb,pg);rows[1]=(*rows[1],pb,pg);rows[2]=(*rows[2],'','')
    with (HERE/'wed_mesh_robustness.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['mesh','h','nodes','bulk_L2_space_time_vs_fine_WED','boundary_L2_space_time_vs_fine_WED','runtime_seconds','newton_iterations','gradient_inf','observed_order_bulk','observed_order_boundary']);w.writerows(rows)

if __name__=='__main__': main()
