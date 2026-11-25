#!/usr/bin/env python3
"""Reproduce the Nizovtseva et al. hyperbolic Allen-Cahn traveling wave in OOFEM.

All outputs are placed in new ``nizovtseva_bulk_v2`` folders; the earlier runs
are left untouched. Requires numpy and the repository's debug_magneto/oofem.
"""
from __future__ import annotations
import argparse, csv, math, subprocess, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OOFEM = ROOT / 'build/debug_magneto/oofem'
TAU = 0.2
T = 0.2
D = MOBILITY = DAMPING = 1.0
RHO = TAU * MOBILITY
BETA, GAMMA = 0.25, 0.5
LAMBDA = math.sqrt(2.0 + 9.0 * TAU) / 2.0
SPEED = 1.5 / LAMBDA


def exact(x, t):
    z = np.clip(LAMBDA * np.asarray(x) - 1.5 * MOBILITY * t, -700., 700.)
    return 1. / (1. + np.exp(z))


def exact_velocity(x, t):
    u = exact(x, t)
    return 1.5 * MOBILITY * u * (1. - u)


def make_input(nx, dt, xmin, xmax, folder, velocity_scale=1.0):
    ny = 1
    nn, ne = (nx + 1) * 2, nx
    nsteps = round(T / dt)
    if abs(nsteps * dt - T) > 1.e-12:
        raise ValueError('dt must divide final time exactly')
    node = lambda i, j: j * (nx + 1) + i + 1
    lines = [
        'case.out', 'Nizovtseva et al. (2016), Eq. (29), bulk traveling wave',
        f'nltransienttransportproblem nsteps {nsteps} deltat {dt:.16g} alpha 1.0 rtol 1.e-10 nsmax 100 nmodules 1 wave 1 newmarkbeta {BETA} newmarkgamma {GAMMA} initialvelocity 0.0 wavedamping {DAMPING/RHO:.16g}',
        'vtkxml tstep_all domain_all primvars 1 6',
        'domain heattransfer', 'OutputManager tstep_all dofman_all element_all',
        f'ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}',
    ]
    values_u, values_v = [], []
    for j in range(2):
        for i in range(nx + 1):
            x = xmin + (xmax - xmin) * i / nx
            lines.append(f'node {node(i,j)} coords 3 {x:.17g} {float(j):.1f} 0.0')
            values_u.append(float(exact(x, 0.)))
            values_v.append(float(velocity_scale * exact_velocity(x, 0.)))
    for i in range(nx):
        lines.append(f'quad1ht {i+1} nodes 4 {node(i,0)} {node(i+1,0)} {node(i+1,1)} {node(i,1)} nip 9')
    lines += [
        'SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1',
        f'allencahn 1 d 1.0 k {D} c 1.0 mobility {MOBILITY} hyperbolic 1 rho {RHO} surfacecapacity 1.0',
    ]
    for n, u in enumerate(values_u, 1):
        lines.append(f'InitialCondition {n} dofs 1 10 Conditions 1 u {u:.17g} set {n+1}')
    lines.append('ConstantFunction 1 f(t) 1.0')
    lines.append(f'Set 1 elementranges {{(1 {ne})}}')
    lines += [f'Set {n+1} nodes 1 {n}' for n in range(1, nn+1)]
    lines[2] += f' initialvelocityvalues {nn} ' + ' '.join(f'{v:.17g}' for v in values_v)
    return '\n'.join(lines) + '\n', np.asarray(values_u), np.asarray(values_v)


def read_vtu(path):
    root = ET.parse(path).getroot()
    points = np.fromstring(root.find('.//Points/DataArray').text, sep=' ').reshape(-1, 3)
    pd = root.find('.//PointData')
    arr = next(a for a in pd.findall('DataArray') if a.attrib.get('Name') == 'Temperature')
    u = np.fromstring(arr.text or '', sep=' ')
    cells = root.find('.//Cells')
    arrays = {a.attrib.get('Name'): np.fromstring(a.text or '', sep=' ', dtype=int) for a in cells.findall('DataArray')}
    conn = arrays['connectivity'].reshape(-1, 4)
    return points, u, conn


def load_history(folder):
    pvd = ET.parse(folder / 'case.out.m0.pvd').getroot()
    frames=[]
    for ds in pvd.findall('.//DataSet'):
        pts,u,conn=read_vtu(folder/ds.attrib['file'])
        frames.append((float(ds.attrib['timestep']),pts,u,conn))
    return frames


def metrics(x, u, xmin, xmax, tag, nx, dt):
    # Two identical y-layers are intentionally collapsed: this is the unique
    # uniform x-grid of the exact y-independent strip solution.
    xu=np.unique(x)
    un=np.array([np.mean(u[np.isclose(x,xx)]) for xx in xu])
    ue=exact(xu,T)
    err=un-ue
    inner=np.abs(xu) <= 6.0 + 1.e-12
    return {'tag':tag,'nx':nx,'h':(xmax-xmin)/nx,'dt':dt,'T':T,
            'rms_full':float(np.sqrt(np.mean(err**2))),
            'max_full':float(np.max(np.abs(err))),
            'rms_inner_absx_le6':float(np.sqrt(np.mean(err[inner]**2))),
            'max_inner_absx_le6':float(np.max(np.abs(err[inner]))),
            'u_min':float(np.min(un)),'u_max':float(np.max(un)),
            'x0_exact':SPEED*T,'x0_numerical':zero_crossing(xu,un),
            'speed_numerical':zero_crossing(xu,un)/T}


def zero_crossing(x,u):
    # exact interface defined by u=1/2, linearly interpolated between nodes
    for i in range(len(x)-1):
        if (u[i]-.5)*(u[i+1]-.5) <= 0 and u[i] != u[i+1]:
            return float(x[i]+(.5-u[i])*(x[i+1]-x[i])/(u[i+1]-u[i]))
    return float('nan')


def write_csv(path, rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)


def run(nx,dt,xmin,xmax,tag,velocity_scale=1.0):
    folder=HERE/f'nizovtseva_bulk_v2_{tag}'
    folder.mkdir(exist_ok=True)
    inp,u0,v0=make_input(nx,dt,xmin,xmax,folder,velocity_scale)
    (folder/'case.in').write_text(inp)
    # Check the nodal velocity data against d/dt of Eq. (29), independent of OOFEM.
    x0=np.tile(np.linspace(xmin,xmax,nx+1),2)
    vcheck=velocity_scale*exact_velocity(x0,0.)
    v_init_error=float(np.max(np.abs(v0-vcheck)))
    p=subprocess.run([str(OOFEM),'-f',str(folder/'case.in')],cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (folder/'run.log').write_text(p.stdout)
    if p.returncode or 'ANALYSIS FINISHED' not in p.stdout:
        raise RuntimeError(f'OOFEM failed: {folder}/run.log')
    hist=load_history(folder)
    time,pts,u,conn=hist[-1]
    if abs(time-T)>1.e-10: raise RuntimeError(f'final output time {time} != {T}')
    row=metrics(pts[:,0],u,xmin,xmax,tag,nx,dt)
    row['initial_velocity_max_mismatch']=v_init_error
    row['initial_velocity_max_abs']=float(np.max(np.abs(v0)))
    print(row,flush=True)
    return row,hist


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--suite',choices=['all','space','time','domain','energy','initial'],default='all')
    args=ap.parse_args()
    if not OOFEM.exists(): raise FileNotFoundError(OOFEM)
    spatial=[]; temporal=[]; domains=[]
    if args.suite in ('all','space'):
        for n in (80,160,320): spatial.append(run(n,.001,-12.,12.,f'space_n{n}')[0])
        spatial.sort(key=lambda r:r['nx'])
        for a,b in zip(spatial,spatial[1:]):
            b['order_full']=math.log(a['rms_full']/b['rms_full'],2.)
            b['order_inner']=math.log(a['rms_inner_absx_le6']/b['rms_inner_absx_le6'],2.)
        write_csv(HERE/'nizovtseva_bulk_v2_spatial.csv',spatial)
    if args.suite in ('all','time'):
        # Larger dt makes the O(dt^2) Newmark error resolvable above the h=.075
        # spatial error; compare exact errors and adjacent-dt differences.
        for dt in (.04,.02,.01,.005): temporal.append(run(320,dt,-12.,12.,f'time_n320_dt{dt:g}')[0])
        temporal.sort(key=lambda r:r['dt'],reverse=True)
        for a,b in zip(temporal,temporal[1:]):
            a['rms_to_next_finer']=float('nan')
            a['inner_rms_to_next_finer']=float('nan')
        # Fields read back from each final output for pairwise temporal differences.
        fields=[]
        for r in temporal:
            f=HERE/f"nizovtseva_bulk_v2_{r['tag']}"
            fr=load_history(f)[-1]
            x=fr[1][:,0]; u=fr[2]; xu=np.unique(x)
            fields.append((xu,np.array([np.mean(u[np.isclose(x,z)]) for z in xu])))
        for i in range(len(temporal)-1):
            diff=fields[i][1]-fields[i+1][1]
            inner=np.abs(fields[i][0])<=6.
            temporal[i]['rms_to_next_finer']=float(np.sqrt(np.mean(diff**2)))
            temporal[i]['inner_rms_to_next_finer']=float(np.sqrt(np.mean(diff[inner]**2)))
        if all(temporal[i]['rms_to_next_finer']>0 for i in range(2)):
            temporal[0]['order_pairwise']=math.log(temporal[0]['rms_to_next_finer']/temporal[1]['rms_to_next_finer'],2.)
        write_csv(HERE/'nizovtseva_bulk_v2_temporal.csv',temporal)
    if args.suite in ('all','domain'):
        # Same h=.1 and dt=.001, compare truncated [-8,8] against [-12,12].
        # nx 160/240 respectively; report whole domain and common interior.
        domains.append(run(160,.001,-8.,8.,'domain_m8_p8')[0])
        domains.append(run(240,.001,-12.,12.,'domain_m12_p12')[0])
        write_csv(HERE/'nizovtseva_bulk_v2_domain.csv',domains)
    if args.suite in ('all','initial'):
        row,_=run(160,.001,-12.,12.,'initial_velocity_plus1pct',1.01)
        write_csv(HERE/'nizovtseva_bulk_v2_initial_sensitivity.csv',[row])
    if args.suite in ('all','energy'):
        # Energy history from an independent postprocessor over each saved frame.
        energy_history(HERE/'nizovtseva_bulk_v2_space_n160')


def energy_history(folder):
    frames=load_history(folder)
    # The OOFEM PVD starts at t=dt; t=0 is reconstructed from the exact
    # initial condition and the exact nodal velocity written in case.in.
    pts0,_,conn0=frames[0][1],frames[0][2],frames[0][3]
    frames.insert(0,(0.0,pts0,exact(pts0[:,0],0.0),conn0))
    times=np.array([f[0] for f in frames]); U=np.array([f[2] for f in frames])
    pts,_,conn=frames[0][1],frames[0][2],frames[0][3]
    x=pts[:,0]; y=pts[:,1]
    # Uniform Q1 rectangles; 3x3 Gauss integrates the quartic Q1 potential exactly.
    gp=[-math.sqrt(3/5),0.,math.sqrt(3/5)]; gw=[5/9,8/9,5/9]; records=[]
    vel=np.gradient(U,times,axis=0,edge_order=2)
    vel[0]=exact_velocity(pts[:,0],0.0)
    for it,t in enumerate(times):
        ek=ep=kin=dis=0.
        for ids in conn:
            xy=pts[ids,:2]; u=U[it,ids]; v=vel[it,ids]
            for ixi,xi in enumerate(gp):
                for jeta,eta in enumerate(gp):
                    N=.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)])
                    dN=.25*np.array([[-(1-eta),-(1-xi)],[(1-eta),-(1+xi)],[(1+eta),(1+xi)],[-(1+eta),(1-xi)]])
                    J=xy.T@dN; det=np.linalg.det(J)
                    gradN=dN@np.linalg.inv(J)
                    uq=N@u; vq=N@v; g=u@gradN
                    w=det*gw[ixi]*gw[jeta]
                    kin += .5*RHO*vq*vq*w
                    ek += .5*D*(g@g)*w
                    ep += MOBILITY*.25*(uq*uq-1.)**2*w
                    dis += DAMPING*vq*vq*w
        records.append({'time':t,'kinetic':kin,'gradient':ek,'potential':ep,'energy':kin+ek+ep,'dissipation_rate':dis})
    integral=np.trapz([r['dissipation_rate'] for r in records],times)
    e0,ef=records[0]['energy'],records[-1]['energy']
    cumulative=np.zeros(len(times))
    cumulative[1:]=np.cumsum(.5*np.diff(times)*(np.asarray([r['dissipation_rate'] for r in records[:-1]])+np.asarray([r['dissipation_rate'] for r in records[1:]])))
    for i,r in enumerate(records): r['energy_balance_residual']=r['energy']-e0+cumulative[i]
    write_csv(HERE/'nizovtseva_bulk_v2_energy.csv',records)
    print({'energy_initial':e0,'energy_final':ef,'dissipated_integral':integral,'balance_residual_final':ef-e0+integral,'relative_balance_residual':(ef-e0+integral)/max(abs(e0),1.e-30)})

if __name__=='__main__': main()
