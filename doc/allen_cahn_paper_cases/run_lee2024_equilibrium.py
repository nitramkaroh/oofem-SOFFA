#!/usr/bin/env python3
"""Generate and run the 1D equilibrium test from Lee et al. (2024), Fig. 7."""
from pathlib import Path
import math
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "doc/allen_cahn_paper_cases/lee2024_equilibrium"
OOFEM = ROOT / "build/opt_magneto/oofem"
NX = 200
H = 2.0 / NX
DT = 0.1 * H * H
NSTEPS = 2000
EPS = 4.0 * H / (2.0 * math.sqrt(2.0) * math.atanh(0.9))


def generate(out=OUT, dt=DT, nsteps=NSTEPS, stride=200):
    out.mkdir(parents=True, exist_ok=True)
    title = out.name
    lines = [
        f"{title}.out",
        "Lee et al. (2024), Fig. 7 equilibrium profile, OOFEM P1 line elements",
        f"nltransienttransportproblem nsteps {nsteps} deltat {dt:.16g} alpha 1.0 rtol 1.e-11 nsmax 40 nmodules 1",
        f"vtkxml tstep_step {stride} domain_all primvars 1 6",
        "domain heattransfer",
        "OutputManager tstep_all dofman_all element_all",
        f"ndofman {NX+1} nelem {NX} ncrosssect 1 nmat 1 nbc 0 nic {NX+1} nltf 1 nset {NX+2}",
    ]
    values = []
    for i in range(NX+1):
        x = -1.0 + i*H
        val = 0.5*math.tanh(x/(math.sqrt(2.0)*EPS))
        values.append(val)
        lines.append(f"node {i+1} coords 3 {x:.16g} 0 0")
    for i in range(NX):
        lines.append(f"line1ht {i+1} nodes 2 {i+1} {i+2}")
    lines += [
        "SimpleTransportCS 1 area 1.0 mat 1 set 1",
        f"allencahn 1 d 1.0 k 1.0 c 1.0 mobility {1.0/(EPS*EPS):.16g}",
    ]
    for i, val in enumerate(values, 1):
        lines.append(f"InitialCondition {i} dofs 1 10 Conditions 1 u {val:.16g} set {i+1}")
    lines.append("ConstantFunction 1 f(t) 1.0")
    lines.append(f"Set 1 elementranges {{(1 {NX})}}")
    for i in range(1, NX+2):
        lines.append(f"Set {i+1} nodes 1 {i}")
    lines[5] = f"OutputManager tstep_step {stride} dofman_all element_all"
    (out / "case.in").write_text("\n".join(lines)+"\n")


def run(out=OUT):
    with (out / "run.log").open("w") as log:
        subprocess.run([str(OOFEM), "-f", "case.in"], cwd=out,
                       stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "t5":
        out = ROOT / "doc/allen_cahn_paper_cases/lee2024_equilibrium_t5"
        generate(out, DT, 500000, 50000)
        run(out)
        print(f"finished Lee et al. equilibrium: epsilon={EPS:.8g}, dt={DT:.8g}, T=5")
    else:
        generate()
        run()
        print(f"finished Lee et al. equilibrium: epsilon={EPS:.8g}, dt={DT:.8g}, T={DT*NSTEPS:g}")
