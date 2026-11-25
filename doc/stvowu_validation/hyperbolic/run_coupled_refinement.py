#!/usr/bin/env python3
"""Temporal refinement of the existing bulk-inertial/dynamic-boundary case."""
from __future__ import annotations

import csv
import math
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "doc/stvowu_validation/allencahn_wave_dynamic.in"
OOFEM = ROOT / "build/debug_magneto/oofem"
T = 0.02
DTS = (0.002, 0.001, 0.0005)


def read_final(path: Path) -> np.ndarray:
    pvd = path.parent / "case.out.m0.pvd"
    root = ET.parse(pvd).getroot()
    fname = root.findall(".//DataSet")[-1].attrib["file"]
    tree = ET.parse(pvd.parent / fname)
    arr = next(a for a in tree.findall(".//PointData/DataArray")
               if a.attrib.get("Name") == "Temperature")
    return np.fromstring(arr.text or "", sep=" ")


def main() -> None:
    source = SOURCE.read_text()
    fields = []
    for dt in DTS:
        nsteps = round(T / dt)
        folder = HERE / f"coupled_dt{dt:g}"
        folder.mkdir(parents=True, exist_ok=True)
        text = source.replace("allencahn_wave_dynamic_result.out", "case.out", 1)
        text = re.sub(r"nsteps\s+\d+\s+deltat\s+\S+", f"nsteps {nsteps} deltat {dt:.16g}", text, count=1)
        (folder / "case.in").write_text(text)
        run = subprocess.run([str(OOFEM), "-f", str(folder / "case.in")], cwd=folder,
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (folder / "run.log").write_text(run.stdout)
        if run.returncode or "ANALYSIS FINISHED" not in run.stdout:
            raise RuntimeError(f"OOFEM failed for dt={dt}; see {folder/'run.log'}")
        fields.append((dt, read_final(folder / "case.out")))

    rows=[]
    for dt, u in fields:
        diff=u-fields[-1][1]
        rows.append({"dt":dt,"T":T,"rms_vs_finest":float(np.sqrt(np.mean(diff**2))),
                     "max_vs_finest":float(np.max(np.abs(diff))),"pairwise_rms_to_next":"",
                     "observed_order":""})
    for i in range(len(rows)-1):
        pairwise=float(np.sqrt(np.mean((fields[i][1]-fields[i+1][1])**2)))
        rows[i]["pairwise_rms_to_next"] = pairwise
    if len(rows) >= 3:
        a,b=rows[0]["pairwise_rms_to_next"],rows[1]["pairwise_rms_to_next"]
        rows[0]["observed_order"] = math.log(a/b,2) if a and b else ""
    with (HERE/"coupled_time_refinement.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    for row in rows: print(row)


if __name__ == "__main__":
    main()
