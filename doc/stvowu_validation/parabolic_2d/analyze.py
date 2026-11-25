#!/usr/bin/env python3
"""Postprocess structured Q1 outputs: energy, traces, fields, and interface metrics."""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def load_vtu(path: Path, n: int) -> np.ndarray:
    piece = ET.parse(path).getroot().find('.//Piece')
    xyz = np.fromstring(piece.find('./Points/DataArray').text, sep=' ').reshape(-1, 3)
    arr = piece.find('./PointData/DataArray')
    if arr is None: raise ValueError(f'No point field in {path}')
    vals = np.fromstring(arr.text, sep=' ')
    u = np.empty((n+1,n+1))
    for (x,y,_), v in zip(xyz, vals): u[int(round(y*n)),int(round(x*n))]=v
    return u

def q1_at(u, r, s):
    N=np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])/4
    dr=np.array([-(1-s),(1-s),(1+s),-(1+s)])/4
    ds=np.array([-(1-r),-(1+r),(1+r),(1-r)])/4
    return N@u, dr@u, ds@u

def energy(u, eps):
    # Three-point Gauss exactly integrates quartic F(u_h) for bilinear Q1.
    z,w=np.polynomial.legendre.leggauss(3); eb=0.
    n=u.shape[0]-1; h=1/n
    for j in range(n):
        for i in range(n):
            v=np.array([u[j,i],u[j,i+1],u[j+1,i+1],u[j+1,i]])
            for r,wr in zip(z,w):
                for s,ws in zip(z,w):
                    val,dr,ds=q1_at(v,r,s)
                    gu2=4*(dr*dr+ds*ds)/(h*h)
                    pot=(val*val-1)**2/(4*eps*eps)
                    eb+=(.5*gu2+pot)*h*h*wr*ws/4
    eg=0.
    # Boundary integral: G=F, tangential gradient from linear Q1 trace.
    z,w=np.polynomial.legendre.leggauss(3)
    sides=[u[0,:],u[-1,:],u[:,0],u[:,-1]]
    for side in sides:
        for i in range(n):
            a,b=side[i],side[i+1]
            for r,wr in zip(z,w):
                N0,N1=(1-r)/2,(1+r)/2
                val=N0*a+N1*b; grad=(b-a)/h
                eg+=(.5*grad*grad+(val*val-1)**2/(4*eps*eps))*h*wr/2
    return eb,eg

def mass_square(v, boundary=False):
    """Integrate v_h^2 with exact-enough 3x3/3-point Q1 Gauss rules."""
    n=v.shape[0]-1; h=1/n; z,w=np.polynomial.legendre.leggauss(3); total=0.
    if boundary:
        for side in (v[0,:],v[-1,:],v[:,0],v[:,-1]):
            for i in range(n):
                a,b=side[i],side[i+1]
                for r,wr in zip(z,w): total+=((1-r)*a/2+(1+r)*b/2)**2*h*wr/2
    else:
        for j in range(n):
            for i in range(n):
                vals=np.array([v[j,i],v[j,i+1],v[j+1,i+1],v[j+1,i]])
                for r,wr in zip(z,w):
                    for s,ws in zip(z,w):
                        val=q1_at(vals,r,s)[0]
                        total+=val*val*h*h*wr*ws/4
    return total

def crossings(values, coords, level=0.):
    out=[]
    for i in range(len(values)-1):
        a,b=values[i]-level,values[i+1]-level
        if a*b<0:
            out.append(coords[i]+(-a)/(b-a)*(coords[i+1]-coords[i]))
    return out

def interface(u):
    n=u.shape[0]-1; y=np.linspace(0,1,n+1); x=np.linspace(0,1,n+1)
    mid=n//2
    # choose the upper crossing of u=0 on x=1/2
    top=crossings(u[:,mid],y)
    bottom=crossings(u[0,:],x)
    return (max(top) if top else float('nan'),
            min(bottom) if bottom else float('nan'),
            max(bottom) if bottom else float('nan'))

def interp_grid(u, xx, yy):
    n=u.shape[0]-1; x=np.asarray(xx); y=np.asarray(yy)
    ix=np.minimum((x*n).astype(int),n-1); iy=np.minimum((y*n).astype(int),n-1)
    r=2*(x*n-ix)-1; s=2*(y*n-iy)-1
    out=np.empty(np.broadcast(x,y).shape)
    for k in np.ndindex(out.shape):
        v=np.array([u[iy[k],ix[k]],u[iy[k],ix[k]+1],u[iy[k]+1,ix[k]+1],u[iy[k]+1,ix[k]]])
        out[k]=q1_at(v,r[k],s[k])[0]
    return out

def rms(a): return float(np.sqrt(np.mean(np.asarray(a)**2)))

def vtu_path(folder,n,dt,step): return folder/f'case.out.m0.{step}.vtu'

def analyze_case(folder,n,dt,eps,steps):
    files=sorted(folder.glob('case.out.m0.*.vtu'),key=lambda p:int(p.stem.split('.')[-1]))
    rows=[]
    for f in files:
        step=int(f.stem.split('.')[-1]); u=load_vtu(f,n); eb,es=energy(u,eps)
        top,left,right=interface(u)
        rows.append({'step':step,'time':step*dt,'E_bulk':eb,'E_surface':es,'E_total':eb+es,
                     'u_min':float(u.min()),'u_max':float(u.max()),'top_y_xmid_u0':top,
                     'bottom_left_u0':left,'bottom_right_u0':right})
    if n==20 and folder.name in ('space_n20_dt0p0001','energy_n20_dt0p0001_every_step'):
        init=np.full((n+1,n+1),-1.)
        init[:n//2+1,n//4:3*n//4+1]=1.
        eb,es=energy(init,eps); rows.insert(0,{'step':0,'time':0.,'E_bulk':eb,'E_surface':es,'E_total':eb+es,
            'u_min':-1.,'u_max':1.,'top_y_xmid_u0':.5,'bottom_left_u0':.25,'bottom_right_u0':.75})
        # Secant-rate dissipation diagnostic over each sampled interval.
        prev=init; prevE=eb+es; budget=[]
        for row in rows[1:]:
            cur=load_vtu(vtu_path(folder,n,dt,int(row['step'])),n)
            dtime=row['time']-(budget[-1]['time'] if budget else 0.)
            change=cur-prev
            db=mass_square(change)/(0.01*dtime)
            ds=mass_square(change,True)/(0.01*dtime)
            drop=prevE-row['E_total']
            budget.append({'time':row['time'],'delta_E':-drop,'bulk_dissipation_secant':db,
                           'surface_dissipation_secant':ds,'dissipation_secant_sum':db+ds,
                           'balance_residual_secant':-drop+db+ds})
            prev=cur;prevE=row['E_total']
        with (HERE/'energy_balance_secant.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=budget[0].keys());w.writeheader();w.writerows(budget)
    with (folder/'energy_interface.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    final=load_vtu(vtu_path(folder,n,dt,steps),n)
    return final,rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=HERE/'runs');p.add_argument('--eps',type=float,default=.1);a=p.parse_args()
    cases={}
    erun=a.root/'energy_n20_dt0p0001_every_step'
    ecsv=erun/'energy_interface.csv'
    if not ecsv.exists() or sum(1 for _ in ecsv.open()) != 1002:
        analyze_case(erun,20,.0001,a.eps,1000)
    # Spatial sequence: all cases end at T=.1 and share dt=5e-4.
    for n in (20,30,40):
        d=a.root/f'space_n{n}_dt0p0005'
        cases[('space',n,.0005)]=analyze_case(d,n,.0005,a.eps,200)
    refrow=cases[('space',40,.0005)][1][-1]
    m=4*40+1; coords=np.linspace(0,1,m); X,Y=np.meshgrid(coords,coords)
    refv=interp_grid(cases[('space',40,.0005)][0],X,Y)
    space=[]
    for n in (20,30,40):
        u=cases[('space',n,.0005)][0];v=interp_grid(u,X,Y);d=v-refv
        trace=np.concatenate([d[0,:],d[-1,:],d[:,0],d[:,-1]])
        bulk_l2=float(np.sqrt(np.trapz(np.trapz(d*d,coords,axis=1),coords)))
        trace_l2=float(np.sqrt(sum(np.trapz(z*z,coords) for z in (d[0,:],d[-1,:],d[:,0],d[:,-1]))))
        erow=cases[('space',n,.0005)][1][-1]
        space.append({'n':n,'h':1/n,'dt':.0005,'bulk_L2_vs_n40':bulk_l2,'bulk_RMS_grid_vs_n40':rms(d),'bulk_max_vs_n40':float(np.max(np.abs(d))),
            'trace_L2_vs_n40':trace_l2,'trace_RMS_nodes_vs_n40':rms(trace),'trace_max_vs_n40':float(np.max(np.abs(trace))),
            'E_bulk':erow['E_bulk'],'E_surface':erow['E_surface'],'E_total':erow['E_total'],
            'E_bulk_rel_diff_pct_vs_n40':100*(erow['E_bulk']-refrow['E_bulk'])/refrow['E_bulk'],
            'E_surface_rel_diff_pct_vs_n40':100*(erow['E_surface']-refrow['E_surface'])/refrow['E_surface'],
            'E_total_rel_diff_pct_vs_n40':100*(erow['E_total']-refrow['E_total'])/refrow['E_total'],
            'u_min':erow['u_min'],'u_max':erow['u_max'],
            'top_y_xmid_u0':erow['top_y_xmid_u0'],'bottom_left_u0':erow['bottom_left_u0'],
            'bottom_right_u0':erow['bottom_right_u0']})
    with (HERE/'spatial_refinement.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=space[0].keys());w.writeheader();w.writerows(space)
    # A finer independent spatial sequence at T=.0625 permits n=60 within the
    # available runtime while keeping h and dt identical across meshes.
    short={}
    for n in (20,30,40,60):
        d=a.root/f'space_T0p0625_n{n}_dt0p000625'
        short[n]=analyze_case(d,n,.000625,a.eps,100)
    shortref=short[60][1][-1]
    ms=4*60+1;cs=np.linspace(0,1,ms);XS,YS=np.meshgrid(cs,cs);shortrefv=interp_grid(short[60][0],XS,YS)
    shortrows=[]
    for n in (20,30,40,60):
        us=short[n][0]; vv=interp_grid(us,XS,YS); dd=vv-shortrefv
        tr=(dd[0,:],dd[-1,:],dd[:,0],dd[:,-1]); sr=short[n][1][-1]
        shortrows.append({'n':n,'h':1/n,'dt':.000625,'T':.0625,
            'bulk_L2_vs_n60':float(np.sqrt(np.trapz(np.trapz(dd*dd,cs,axis=1),cs))),
            'bulk_max_vs_n60':float(np.max(np.abs(dd))),
            'trace_L2_vs_n60':float(np.sqrt(sum(np.trapz(z*z,cs) for z in tr))),
            'trace_max_vs_n60':float(max(np.max(np.abs(z)) for z in tr)),
            'E_bulk':sr['E_bulk'],'E_surface':sr['E_surface'],'E_total':sr['E_total'],
            'E_total_rel_diff_pct_vs_n60':100*(sr['E_total']-shortref['E_total'])/shortref['E_total'],
            'u_min':sr['u_min'],'u_max':sr['u_max'],
            'top_y_xmid_u0':sr['top_y_xmid_u0'],'bottom_left_u0':sr['bottom_left_u0'],
            'bottom_right_u0':sr['bottom_right_u0']})
    with (HERE/'spatial_refinement_T0p0625.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=shortrows[0].keys());w.writeheader();w.writerows(shortrows)
    # Temporal sequence at fixed n=20, with the existing dt=1e-4 run as reference.
    temporal=[]; runs=[]
    for dt,folder,step in ((.0004,a.root/'time_n20_dt0p0004',250),(.0002,a.root/'time_n20_dt0p0002',500),(.0001,a.root/'energy_n20_dt0p0001_every_step',1000)):
        u=load_vtu(vtu_path(folder,20,dt,step),20);runs.append((dt,u,energy(u,a.eps),interface(u)))
    refu=runs[-1][1]
    for dt,u,(eb,es),ints in runs:
        diff=u-refu;trace=np.concatenate([diff[0,:],diff[-1,:],diff[:,0],diff[:,-1]])
        trace_l2=math.sqrt(mass_square(diff,True))
        temporal.append({'dt':dt,'bulk_nodal_RMS_vs_dt1e-4':rms(diff),'bulk_nodal_max_vs_dt1e-4':float(np.max(np.abs(diff))),
            'bulk_L2_vs_dt1e-4':math.sqrt(mass_square(diff)),'trace_L2_vs_dt1e-4':trace_l2,
            'trace_nodal_RMS_vs_dt1e-4':rms(trace),'trace_nodal_max_vs_dt1e-4':float(np.max(np.abs(trace))),
            'E_bulk':eb,'E_surface':es,'E_total':eb+es,'top_y_xmid_u0':ints[0],
            'bottom_left_u0':ints[1],'bottom_right_u0':ints[2]})
    adjacent=[]; adjacent_trace=[]
    for (_,u0,_,_),(_,u1,_,_) in zip(runs[:-1],runs[1:]):
        adjacent.append(rms(u0-u1))
        dd=u0-u1; adjacent_trace.append(rms(np.concatenate([dd[0,:],dd[-1,:],dd[:,0],dd[:,-1]])))
    ratio=adjacent[0]/adjacent[1] if adjacent[1] else float('nan')
    order=math.log(ratio,2) if ratio>0 else float('nan')
    for i,row in enumerate(temporal):
        row['adjacent_field_RMS_to_next_finer']=adjacent[i] if i<len(adjacent) else float('nan')
        row['adjacent_trace_RMS_to_next_finer']=adjacent_trace[i] if i<len(adjacent_trace) else float('nan')
        row['observed_order_from_adjacent_differences']=order
    uc=load_vtu(vtu_path(a.root/'space_n20_dt0p0005',20,.0005,200),20)
    ebc,esc=energy(uc,a.eps); ic=interface(uc); dc=uc-refu
    tc=np.concatenate([dc[0,:],dc[-1,:],dc[:,0],dc[:,-1]])
    tcl2=math.sqrt(mass_square(dc,True))
    temporal.append({'dt':.0005,'bulk_nodal_RMS_vs_dt1e-4':rms(dc),'bulk_nodal_max_vs_dt1e-4':float(np.max(np.abs(dc))),
        'bulk_L2_vs_dt1e-4':math.sqrt(mass_square(dc)),'trace_L2_vs_dt1e-4':tcl2,
        'trace_nodal_RMS_vs_dt1e-4':rms(tc),'trace_nodal_max_vs_dt1e-4':float(np.max(np.abs(tc))),
        'E_bulk':ebc,'E_surface':esc,'E_total':ebc+esc,'top_y_xmid_u0':ic[0],
        'bottom_left_u0':ic[1],'bottom_right_u0':ic[2],
        'adjacent_field_RMS_to_next_finer':float('nan'),'adjacent_trace_RMS_to_next_finer':float('nan'),
        'observed_order_from_adjacent_differences':order})
    with (HERE/'temporal_refinement.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=temporal[0].keys());w.writeheader();w.writerows(temporal)
    # Orientation/permutation invariance on a small square mesh.
    base=load_vtu(vtu_path(a.root/'orientation_base_n8_dt0p001',8,.001,10),8)
    other=load_vtu(vtu_path(a.root/'orientation_reverse_n8_dt0p001',8,.001,10),8)
    (HERE/'geometry_sensitivity.csv').write_text('case,n,dt,bulk_nodal_RMS_vs_standard,bulk_nodal_max_vs_standard\n'
        f'reversed_and_permuted_boundary,8,0.001,{rms(other-base):.16g},{np.max(np.abs(other-base)):.16g}\n')
    # Figures.
    fig,axs=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
    for ax,n in zip(axs,(20,30,40)):
        u=cases[('space',n,.0005)][0];xy=np.linspace(0,1,n+1)
        im=ax.contourf(xy,xy,u,levels=np.linspace(-1.05,1.05,22),cmap='viridis',vmin=-1,vmax=1,extend='both')
        ax.contour(xy,xy,u,levels=[0],colors='k',linewidths=1)
        ax.set_title(f'Q1 {n}×{n}, T=0.1');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_aspect('equal')
    fig.colorbar(im,ax=axs,label='u');fig.savefig(HERE/'spatial_fields_T0p1.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(9,8),layout='constrained')
    for ax,n in zip(axs.flat,(20,30,40,60)):
        u=short[n][0];xy=np.linspace(0,1,n+1)
        im=ax.contourf(xy,xy,u,levels=np.linspace(-1.05,1.05,22),cmap='viridis',vmin=-1,vmax=1,extend='both')
        ax.contour(xy,xy,u,levels=[0],colors='white',linewidths=1)
        ax.set_title(f'Q1 {n}×{n}, T=0.0625');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_aspect('equal')
    fig.colorbar(im,ax=axs,label='u');fig.savefig(HERE/'spatial_fields_T0p0625.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    ax.loglog([r['h'] for r in shortrows[:-1]],[r['bulk_L2_vs_n60'] for r in shortrows[:-1]],'o-',label='bulk L2')
    ax.loglog([r['h'] for r in shortrows[:-1]],[r['trace_L2_vs_n60'] for r in shortrows[:-1]],'s-',label='trace L2')
    ax.invert_xaxis();ax.set(xlabel='h',ylabel='difference to n=60',title='Spatial refinement, T=0.0625');ax.grid(True,which='both');ax.legend()
    fig.savefig(HERE/'spatial_convergence_T0p0625.png',dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    axs[0].loglog([r['h'] for r in space[:-1]],[r['bulk_L2_vs_n40'] for r in space[:-1]],'o-',label='bulk L2')
    axs[0].loglog([r['h'] for r in space[:-1]],[r['trace_L2_vs_n40'] for r in space[:-1]],'s-',label='trace L2')
    axs[0].invert_xaxis();axs[0].set_xlabel('h');axs[0].set_ylabel('difference to n=40');axs[0].grid(True,which='both');axs[0].legend()
    for n in (20,30,40):
        r=cases[('space',n,.0005)][1][-1];axs[1].scatter([n],[r['E_bulk']],label=f'n={n}')
    axs[1].set(xlabel='N elements per side',ylabel='bulk energy, T=0.1',title='Spatial energy sensitivity');axs[1].grid(True);axs[1].legend()
    fig.savefig(HERE/'refinement_metrics.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.2,4),layout='constrained')
    energy_rows=list(csv.DictReader((HERE/'runs/energy_n20_dt0p0001_every_step/energy_interface.csv').open()))
    times=[float(r['time']) for r in energy_rows]
    for key,label in (('E_bulk','bulk'),('E_surface','surface'),('E_total','total')):
        ax.plot(times,[float(r[key]) for r in energy_rows],label=label)
    ax.set(xlabel='t',ylabel='energy',title='Energy components, n=20, Δt=10⁻⁴; ε=0.1 assumed');ax.grid(True);ax.legend()
    fig.savefig(HERE/'energy_time.png',dpi=180);plt.close(fig)
    budget=list(csv.DictReader((HERE/'energy_balance_secant.csv').open()))
    bt=np.array([float(r['time']) for r in budget]);dstep=1.e-4
    drop_rate=np.array([-float(r['delta_E'])/dstep for r in budget])
    db=np.array([float(r['bulk_dissipation_secant']) for r in budget])
    ds=np.array([float(r['surface_dissipation_secant']) for r in budget])
    resid=np.array([float(r['balance_residual_secant']) for r in budget])
    relative=np.abs(resid)/np.maximum(np.abs(np.array([float(r['delta_E']) for r in budget])),1.e-30)
    fig,axs=plt.subplots(2,1,figsize=(7,6),layout='constrained',sharex=True)
    axs[0].plot(bt,drop_rate,label='−ΔE/Δt');axs[0].plot(bt,db,label='bulk secant dissipation');axs[0].plot(bt,ds,label='surface secant dissipation')
    axs[0].plot(bt,db+ds,'--',label='sum of secant dissipation');axs[0].set_ylabel('rate');axs[0].grid(True);axs[0].legend(fontsize=8,ncol=2)
    axs[1].semilogy(bt,np.maximum(relative,1.e-12));axs[1].set(xlabel='t',ylabel='|balance residual| / |ΔE|',title='Secant energy-balance diagnostic');axs[1].grid(True,which='both')
    fig.savefig(HERE/'energy_balance_secant.png',dpi=180);plt.close(fig)
    print('Spatial:');[print(r) for r in space]
    print('Spatial T=.0625:');[print(r) for r in shortrows]
    print('Temporal:');[print(r) for r in temporal]
    print('Geometry:',(HERE/'geometry_sensitivity.csv').read_text().strip())

if __name__=='__main__': main()
