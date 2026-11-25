#!/usr/bin/env python3
"""Digitize Okumura thesis Fig. 4.4 and compare it to the OOFEM energy curve.

The figure is printed on dissertation page 71 (PDF page 80). Pixel-to-data
values are approximate and are retained as a separate CSV for transparency.
"""
from pathlib import Path
import argparse, subprocess, tempfile, csv
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

BASE=Path(__file__).parent
parser=argparse.ArgumentParser()
parser.add_argument('--pdf',default='/home/marco/Documents/Work/Riccardo/refs/31954_Dissertation.pdf')
args=parser.parse_args()
with tempfile.TemporaryDirectory(prefix='okumura_fig44_') as tmp:
    prefix=str(Path(tmp)/'page')
    subprocess.run(['pdftoppm','-f','80','-l','80','-png','-r','300',args.pdf,prefix],check=True,stdout=subprocess.DEVNULL)
    page=Image.open(prefix+'-080.png').convert('L')
    # Crop containing only Fig. 4.4; coordinates are for A4 at 300 dpi.
    graph=np.asarray(page.crop((1350,2380,2130,2835)))
left,right,top,bottom=38,739,20,431
x_tick_px=[138,238,338,438,538,638]
samples=[]; previous_y=88
for x in range(42,735):
    if any(abs(x-z)<=2 for z in x_tick_px): continue
    ys=np.flatnonzero(graph[top+8:bottom,x]<130)+top+8
    ys=ys[ys<bottom-1]
    # Ignore the legend in the upper right; the energy curve is below it.
    if x>550: ys=ys[ys>75]
    if not len(ys): continue
    y=int(ys[np.argmin(np.abs(ys-previous_y))]);previous_y=y
    t=(x-left)/(right-left)*.7
    energy=(bottom-y)/(bottom-top)*90
    samples.append((t,energy))
with (BASE/'okumura_ch4_example2_paper_fig44_digitized.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['time','paper_energy_digitized_from_Fig4_4']);w.writerows(samples)

paper_t=np.array([r[0] for r in samples]);paper_e=np.array([r[1] for r in samples])
oof=np.genfromtxt(BASE/'okumura_ch4_example2_energy.csv',delimiter=',',names=True)
mask=(oof['time']>=.01)&(oof['time']<=.4)
ref=np.interp(oof['time'][mask],paper_t,paper_e)
err=oof['OOFEM_Q1_total_energy'][mask]-ref
summary={'samples':int(mask.sum()),'comparison_interval':'0.01 <= t <= 0.4',
         'paper_initial_energy_digitized':float(paper_e[0]),
         'oofem_energy_at_0.01':float(oof['OOFEM_Q1_total_energy'][np.argmin(abs(oof['time']-.01))]),
         'oofem_energy_at_0.1':float(oof['OOFEM_Q1_total_energy'][np.argmin(abs(oof['time']-.1))]),
         'oofem_energy_at_0.2':float(oof['OOFEM_Q1_total_energy'][np.argmin(abs(oof['time']-.2))]),
         'oofem_energy_at_0.3':float(oof['OOFEM_Q1_total_energy'][np.argmin(abs(oof['time']-.3))]),
         'oofem_energy_at_0.4':float(oof['OOFEM_Q1_total_energy'][np.argmin(abs(oof['time']-.4))]),
         'digitized_curve_RMSE':float(np.sqrt(np.mean(err**2))),
         'digitized_curve_MAE':float(np.mean(np.abs(err))),
         'digitized_curve_max_abs_error':float(np.max(np.abs(err)))}
with (BASE/'okumura_ch4_example2_literature_comparison.txt').open('w') as f:
    for k,v in summary.items():f.write(f'{k}: {v}\n')
fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
ax.plot(oof['time'],oof['OOFEM_Q1_total_energy'],label='OOFEM Q1 FEM',lw=2)
ax.plot(paper_t,paper_e,label='Okumura Fig. 4.4, digitized',lw=1.5,alpha=.8)
ax.set(xlim=(0,.7),ylim=(0,90),xlabel='t',ylabel='Total energy',title='Okumura dissertation, Ch. 4 Example 2')
ax.grid(alpha=.25);ax.legend();fig.savefig(BASE/'okumura_ch4_example2_energy_overlay.png',dpi=180)
print(summary)
