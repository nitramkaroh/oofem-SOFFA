#!/usr/bin/env python3
"""Check the Hochbruck--Leibold exact profile on a polygonal disk in OOFEM.

Forcing terms are rebuilt for OOFEM's bulk-inertial/first-order-surface
operator. This checks the exact profile and Newmark coupling, not the paper's
second-order acoustic boundary equation or its IMEX discretization.
"""
from __future__ import annotations

import csv
import math
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OOFEM = ROOT / "build/debug_magneto/oofem"
GMSH = "/usr/bin/gmsh"
T = 0.2


def gmsh_mesh(folder: Path, h: float) -> tuple[list[tuple[float,float]], list[tuple[int,int,int]], list[tuple[int,int]]]:
    geo = folder / "disk.geo"
    msh = folder / "disk.msh"
    geo.write_text(f"""SetFactory(\"OpenCASCADE\");
Disk(1) = {{0, 0, 0, 1, 1}};
Mesh.CharacteristicLengthMin = {h};
Mesh.CharacteristicLengthMax = {h};
Mesh.MshFileVersion = 2.2;
""")
    proc = subprocess.run([GMSH, "-2", str(geo), "-format", "msh2", "-o", str(msh)],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode or not msh.exists():
        raise RuntimeError(f"gmsh failed for h={h}: {proc.stdout}")
    lines = msh.read_text().splitlines()
    i = lines.index("$Nodes")
    count = int(lines[i+1]); raw_nodes = lines[i+2:i+2+count]
    tags = {}
    coords = []
    for row in raw_nodes:
        vals = row.split(); tags[int(vals[0])] = len(coords)+1
        coords.append((float(vals[1]), float(vals[2])))
    j = lines.index("$Elements")
    count = int(lines[j+1]); raw_elements = lines[j+2:j+2+count]
    tris, edges = [], []
    for row in raw_elements:
        vals = list(map(int, row.split())); etype, ntags = vals[1], vals[2]
        conn = vals[3+ntags:]
        if etype == 2:
            tris.append(tuple(tags[v] for v in conn))
        elif etype == 1:
            edges.append(tuple(tags[v] for v in conn))
    if not tris or not edges:
        raise RuntimeError("gmsh mesh did not contain first-order triangles and boundary lines")
    return coords, tris, edges


def make_input(folder: Path, h: float, dt: float) -> tuple[str, list[tuple[float,float]]]:
    coords, tris, edges = gmsh_mesh(folder, h)
    nn, ne = len(coords), len(tris)+len(edges)
    steps = round(T/dt)
    velocities = [2*math.pi*x*y for x,y in coords]
    lines = [
        "case.out", "Hochbruck-Leibold profile on a gmsh disk; OOFEM boundary operator",
        f"nltransienttransportproblem nsteps {steps} deltat {dt:.16g} alpha 1.0 rtol 1.e-9 nsmax 60 nmodules 1 wave 1 newmarkbeta 0.25 newmarkgamma 0.5 initialvelocity 0.0 wavedamping 1.0",
        "vtkxml tstep_all domain_all primvars 1 6", "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc 0 nic {nn} nltf 1 nset {nn+1}",
    ]
    for n,(x,y) in enumerate(coords,1):
        lines.append(f"node {n} coords 3 {x:.16g} {y:.16g} 0.0")
    eid=1
    for tri in tris:
        lines.append(f"tr1ht {eid} nodes 3 {tri[0]} {tri[1]} {tri[2]}"); eid+=1
    for edge in edges:
        lines.append(f"line1ht {eid} nodes 2 {edge[0]} {edge[1]} nip 3"); eid+=1
    lines += [
        "SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
        "allencahn 1 d 1.0 k 1.0 c 1.0 mobility 1.0 hyperbolic 1 rho 1.0 surfacecapacity 1.0 manufactured 1 mmsbase 0.0 mmsamplitude 1.0 mmsrate 6.283185307179586 mmsbeta 0.0 mmslx 1.0 mmsly 1.0 mmstemporal 1 mmsspatial 3",
    ]
    for n in range(1,nn+1):
        lines.append(f"InitialCondition {n} dofs 1 10 Conditions 1 u 0.0 set {n+1}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {ne})}}")
    lines += [f"Set {n+1} nodes 1 {n}" for n in range(1,nn+1)]
    lines[2] += f" initialvelocityvalues {nn} " + " ".join(f"{v:.17g}" for v in velocities)
    return "\n".join(lines)+"\n", coords


def final_field(folder: Path) -> tuple[np.ndarray,np.ndarray]:
    pvd=folder/"case.out.m0.pvd"; root=ET.parse(pvd).getroot()
    vtu=folder/root.findall(".//DataSet")[-1].attrib["file"]
    tree=ET.parse(vtu).getroot()
    coords=np.fromstring(tree.find(".//Points/DataArray").text,sep=" ").reshape(-1,3)
    arr=next(a for a in tree.findall(".//PointData/DataArray") if a.attrib.get("Name")=="Temperature")
    return coords,np.fromstring(arr.text or "",sep=" ")


def run(h: float, dt: float, tag: str) -> tuple[np.ndarray,np.ndarray,dict]:
    folder=HERE/f"hl_disk_{tag}";folder.mkdir(parents=True,exist_ok=True)
    inp, _=make_input(folder,h,dt);(folder/"case.in").write_text(inp)
    res=subprocess.run([str(OOFEM),"-f",str(folder/"case.in")],cwd=folder,
                       stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (folder/"run.log").write_text(res.stdout)
    if res.returncode or "ANALYSIS FINISHED" not in res.stdout:
        raise RuntimeError(f"OOFEM failed for {tag}; see {folder/'run.log'}")
    coords,u=final_field(folder);ref=math.sin(2*math.pi*T)*coords[:,0]*coords[:,1]
    row={"tag":tag,"h_target":h,"dt":dt,"T":T,
         "nodes":len(u),"rms_error_vs_exact_profile":float(np.sqrt(np.mean((u-ref)**2))),
         "max_error_vs_exact_profile":float(np.max(np.abs(u-ref)))}
    print(row,flush=True);return coords,u,row


def main() -> None:
    rows=[];temporal=[]
    for dt in (0.01,0.005,0.0025):
        _,u,r=run(0.16,dt,f"time_h016_dt{dt:g}");rows.append(r);temporal.append((dt,u))
    for (d0,u0),(d1,u1) in zip(temporal,temporal[1:]):
        next(r for r in rows if r["dt"]==d0)["rms_to_next_finer"]=float(np.sqrt(np.mean((u0-u1)**2)))
    rows[0]["observed_temporal_order"]=math.log(rows[0]["rms_to_next_finer"]/rows[1]["rms_to_next_finer"],2)
    for h in (0.3,0.16,0.09):
        _,_,r=run(h,0.0025,f"space_h{h:g}_dt0p0025");rows.append(r)
    space=[r for r in rows if r["tag"].startswith("space_")];space.sort(key=lambda r:r["h_target"],reverse=True)
    for a,b in zip(space,space[1:]):b["observed_spatial_order"]=math.log(a["rms_error_vs_exact_profile"]/b["rms_error_vs_exact_profile"],2)
    cols=sorted({k for r in rows for k in r})
    with (HERE/"hochbruck_leibold_disk_metrics.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)


if __name__=="__main__":main()
