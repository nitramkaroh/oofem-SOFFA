#!/usr/bin/env python3
"""Spatial check of the Allen-Cahn bulk against its stationary tanh front."""
from __future__ import annotations

import csv
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OOFEM=ROOT/"build/debug_magneto/oofem"
T,DT=0.2,0.001


def make_input(n:int)->str:
    ny=1; nn=(n+1)*(ny+1); ne=n*ny
    nid=lambda i,j:j*(n+1)+i+1
    lines=["case.out","Stationary hyperbolic Allen-Cahn tanh front",
           f"nltransienttransportproblem nsteps {round(T/DT)} deltat {DT} alpha 1.0 rtol 1.e-8 nsmax 100 nmodules 1 wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0",
           "vtkxml tstep_all domain_all primvars 1 6","domain heattransfer",
           "OutputManager tstep_all dofman_all element_all",
           f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}"]
    for j in range(2):
        for i in range(n+1):
            x=-8+16*i/n
            lines.append(f"node {nid(i,j)} coords 3 {x:.16g} {float(j):.1f} 0.0")
    e=1
    for i in range(n):
        lines.append(f"quad1ht {e} nodes 4 {nid(i,0)} {nid(i+1,0)} {nid(i+1,1)} {nid(i,1)} nip 9")
        e+=1
    lines += ["SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
              "allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 hyperbolic 1 rho 1.0 surfacecapacity 1.0"]
    for j in range(2):
        for i in range(n+1):
            node=nid(i,j); x=-8+16*i/n; u=math.tanh(x/math.sqrt(2))
            lines.append(f"InitialCondition {node} dofs 1 10 Conditions 1 u {u:.17g} set {node+1}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {ne})}}")
    lines += [f"Set {node+1} nodes 1 {node}" for node in range(1,nn+1)]
    return "\n".join(lines)+"\n"


def read_field(folder:Path)->np.ndarray:
    pvd=ET.parse(folder/"case.out.m0.pvd").getroot()
    vtu=folder/pvd.findall(".//DataSet")[-1].attrib["file"]
    tree=ET.parse(vtu)
    arr=next(a for a in tree.findall(".//PointData/DataArray") if a.attrib.get("Name")=="Temperature")
    return np.fromstring(arr.text or "",sep=" ")


def main()->None:
    rows=[]
    for n in (40,80,160):
        folder=HERE/f"kink_n{n}"
        folder.mkdir(parents=True,exist_ok=True)
        (folder/"case.in").write_text(make_input(n))
        run=subprocess.run([str(OOFEM),"-f",str(folder/"case.in")],cwd=folder,
                           text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (folder/"run.log").write_text(run.stdout)
        if run.returncode or "ANALYSIS FINISHED" not in run.stdout:
            raise RuntimeError(f"OOFEM failed for N={n}; see {folder/'run.log'}")
        u=read_field(folder)
        x=np.array([-8+16*i/n for j in range(2) for i in range(n+1)])
        ref=np.tanh(x/math.sqrt(2))
        err=u-ref
        rows.append({"nx":n,"h":16/n,"T":T,"rms_error":float(np.sqrt(np.mean(err**2))),
                     "max_error":float(np.max(np.abs(err))),"min_u":float(u.min()),"max_u":float(u.max())})
    for a,b in zip(rows,rows[1:]):
        a["observed_order"]=math.log(a["rms_error"]/b["rms_error"],2)
    rows[-1]["observed_order"]=""
    with (HERE/"stationary_kink_space_refinement.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["nx","h","T","rms_error","max_error","min_u","max_u","observed_order"])
        w.writeheader();w.writerows(rows)
    for row in rows: print(row)


if __name__=="__main__": main()
