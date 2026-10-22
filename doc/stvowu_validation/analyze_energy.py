from pathlib import Path
import xml.etree.ElementTree as ET, math, csv
import numpy as np
p=Path('/tmp/stvowu_validation')
def energy(file):
 root=ET.parse(file).getroot(); piece=root.find('.//Piece')
 xyz=np.array(list(map(float,piece.find('./Points/DataArray').text.split()))).reshape(-1,3)[:,:2]
 vals=np.array(list(map(float,piece.find('./PointData/DataArray').text.split())))
 cells=piece.find('./Cells'); conn=np.array(list(map(int,cells.find("./DataArray[@Name='connectivity']").text.split()))); offsets=np.array(list(map(int,cells.find("./DataArray[@Name='offsets']").text.split()))); types=list(map(int,cells.find("./DataArray[@Name='types']").text.split()))
 F=lambda u: .25*(u*u-1)**2
 total=0.
 a=1/math.sqrt(3)
 for i,typ in enumerate(types):
  ids=conn[offsets[i-1] if i else 0:offsets[i]]; X=xyz[ids]; U=vals[ids]
  if typ==9:
   for r in [-a,a]:
    for s in [-a,a]:
     N=np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])*.25
     dN=np.array([[-(1-s),-(1-r)],[(1-s),-(1+r)],[(1+s),(1+r)],[-(1+s),(1-r)]])*.25
     J=X.T@dN; det=np.linalg.det(J); grad=np.linalg.solve(J.T,dN.T).T.T@U
     u=N@U; total+=(.5*grad@grad+F(u))*det
  elif typ==3:
   length=np.linalg.norm(X[1]-X[0]); tangent=(X[1]-X[0])/length
   for r in [-a,a]:
    N=np.array([(1-r)/2,(1+r)/2]); dUdr=np.array([-1.,1.])@U; grad=dUdr/length
    u=N@U; total+=(.5*grad*grad+F(u))*length/2
  else: raise ValueError(typ)
 return total
rows=[]
for step in range(1,41):
 f=p/f'energy_decay.out.m0.{step}.vtu'; e=energy(f); rows.append((step,step*.0125,e))
print('first last minmax monotonically decreasing',rows[0],rows[-1],min(b[2]-a[2] for a,b in zip(rows,rows[1:])),sum(b[2]>a[2]+1e-11 for a,b in zip(rows,rows[1:])))
with open('doc/stvowu_validation/energy.csv','w') as f:
 w=csv.writer(f); w.writerow(['step','time','discrete_total_energy']); w.writerows(rows)
