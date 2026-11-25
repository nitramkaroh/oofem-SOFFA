#!/usr/bin/env python3
"""Postprocess the three OOFEM temporal-refinement runs for literature comparison."""
from pathlib import Path
import csv, importlib.util
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'doc/stvowu_validation/run_energy_balance.py'
TEMP=Path('/tmp/stvowu_energy_balance')
spec=importlib.util.spec_from_file_location('energy_balance',SOURCE)
eb=importlib.util.module_from_spec(spec); spec.loader.exec_module(eb)
rows=[]; trajectories={}; endpoints={}
for dt,label in [(0.005,'dt_0p005'),(0.0025,'dt_0p0025'),(0.00125,'dt_0p00125')]:
    folder=TEMP/label; inp=folder/'case.in'
    nodes,elems,u0=eb.geometry(inp); mb,mg=eb.assembled(nodes,elems); mass=mb+mg
    energy0=sum(eb.energy(u0,nodes,elems))
    files=sorted(folder.glob('energy_balance_*.out.m0.*.vtu'),key=lambda f:int(f.stem.split('.')[-1]))
    times=[0.0]; es=[energy0]; ebvals=[eb.energy(u0,nodes,elems)[0]]; egvals=[eb.energy(u0,nodes,elems)[1]]
    for i,f in enumerate(files,1):
        u=eb.load_values(f); a,b=eb.energy(u,nodes,elems)
        times.append(i*dt); es.append(a+b); ebvals.append(a); egvals.append(b)
    trajectories[dt]=(np.array(times),np.array(es),np.array(ebvals),np.array(egvals)); endpoints[dt]=eb.load_values(files[-1])
    mass0=float(np.ones(len(u0))@mass@u0); massT=float(np.ones(len(u0))@mass@endpoints[dt])
    rows.append({'dt':dt,'steps':len(files),'E0':energy0,'ET':es[-1],'energy_drop':energy0-es[-1],
                 'max_energy_increase':float(np.max(np.diff(es))), 'mass_initial':mass0,'mass_final':massT})
# Consistent weighted norms for common spatial mesh and fixed final time.
for row,dt in zip(rows,(.005,.0025,.00125)):
    d=endpoints[dt]-endpoints[.00125]
    n=len(d); nodes,elems,_=eb.geometry(TEMP/'dt_0p00125'/'case.in'); mb,mg=eb.assembled(nodes,elems); M=mb+mg
    row['final_error_vs_dt_0p00125_L2']=float(np.sqrt(max(0,d@M@d)))
# Empirical time order from the ratio of successive endpoint differences.
d1=endpoints[.005]-endpoints[.0025]; d2=endpoints[.0025]-endpoints[.00125]
nodes,elems,_=eb.geometry(TEMP/'dt_0p00125'/'case.in'); mb,mg=eb.assembled(nodes,elems); M=mb+mg
ratio=float(np.sqrt(d1@M@d1)/np.sqrt(d2@M@d2)); observed=float(np.log2(ratio))
with (HERE/'metrics.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with (HERE/'energy_trajectory.csv').open('w',newline='') as f:
    w=csv.writer(f); w.writerow(['time','E_dt_0.005','E_dt_0.0025','E_dt_0.00125'])
    # Shared common times only, interpolated for visualization.
    ts=np.linspace(0,.1,81)
    vals=[np.interp(ts,trajectories[dt][0],trajectories[dt][1]) for dt in (.005,.0025,.00125)]
    for i,t in enumerate(ts): w.writerow([t,*[v[i] for v in vals]])
plt.figure(figsize=(7,4.5))
for dt in (.005,.0025,.00125):
    t,e,_,_=trajectories[dt]; plt.plot(t,e,label=fr'$\Delta t={dt:g}$')
plt.xlabel('t'); plt.ylabel('total Allen–Cahn energy'); plt.grid(True,alpha=.25); plt.legend(); plt.tight_layout()
plt.savefig(HERE/'energy_decay.png',dpi=180); plt.close()
(HERE/'observed_order.txt').write_text(f'endpoint difference ratio = {ratio:.8g}\nobserved temporal order = {observed:.8g}\n')
print(f'endpoint difference ratio {ratio:.6f}; observed temporal order {observed:.4f}')
for row in rows: print(row)
