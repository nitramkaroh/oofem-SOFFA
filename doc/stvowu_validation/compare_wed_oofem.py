#!/usr/bin/env python3
"""Compare WED minimizer paths with the standard OOFEM transient solution."""
from pathlib import Path
import xml.etree.ElementTree as ET, numpy as np, csv, math
p=Path('/tmp/stvowu_validation'); work=Path(__file__).parent

def read_vtu(step):
    root=ET.parse(p/f'allen_cahn_direct.out.m0.{step}.vtu').getroot(); piece=root.find('.//Piece')
    xyz=np.array(list(map(float,piece.find('./Points/DataArray').text.split()))).reshape(-1,3)
    u=np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))
    return xyz,u
rows=[]
for eps,filename in [(0.1,'wed_newton_0.1.npz'),(.05,'wed_newton_0.05.npz'),(.025,'wed_newton_0.025.npz')]:
    z=np.load(p/filename); times=z['time']; U=z['solution']; final=float(times[-1])
    for t in (.05,.1,.2):
        if t>final+1.e-12: continue
        j=int(round(t/float(times[1])))
        _,u_ref=read_vtu(round(t/.001))
        diff=U[j]-u_ref
        rows.append((eps,t,math.sqrt(float(diff@diff)/len(diff)),float(np.max(np.abs(diff))),float(z['gradient_inf']),str(z['success'])))
with (work/'wed_vs_oofem.csv').open('w') as f:
    w=csv.writer(f); w.writerow(['epsilon','physical_time','nodal_rms_difference','max_abs_difference','optimizer_gradient_inf','optimizer_success']); w.writerows(rows)
for r in rows: print(r)
