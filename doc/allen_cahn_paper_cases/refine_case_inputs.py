#!/usr/bin/env python3
"""Generate uniformly refined P1 versions of the six paper-comparison cases.

The ellipse initial data use the signed Euclidean distance to the ellipse.
Outputs go to a fresh run root so the original runs remain available.
"""
from pathlib import Path
import argparse
import math

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "reproduction_inputs"


def signed_ellipse_distance(x, y, a=0.6, b=0.2):
    # Closest point on the ellipse, found by minimizing over its parameter.
    # A dense bracket followed by golden-section refinement is robust at axes.
    def f(t):
        return (a * math.cos(t) - x) ** 2 + (b * math.sin(t) - y) ** 2
    m = 256
    k = min(range(m), key=lambda i: f(2 * math.pi * i / m))
    lo, hi = 2 * math.pi * (k - 1) / m, 2 * math.pi * (k + 1) / m
    gr = (math.sqrt(5) - 1) / 2
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    for _ in range(45):
        if f(c) < f(d):
            hi, d = d, c
            c = hi - gr * (hi - lo)
        else:
            lo, c = c, d
            d = lo + gr * (hi - lo)
    dist = math.sqrt(f((lo + hi) / 2))
    inside = (x / a) ** 2 + (y / b) ** 2 < 1
    return -dist if inside else dist


def make_case(template: Path, dest: Path, n: int, dt_override=None, nsteps_override=None,
              bdf_order=1, full_newton=False):
    lines = template.read_text().splitlines()
    name = template.stem
    is_wang = name.startswith("wang")
    eps_digits = name.split("eps0p")[1]
    eps = int(eps_digits) / (10 if len(eps_digits) == 1 else 100)
    dtline = lines[2]
    nsteps = int(dtline.split("nsteps ")[1].split()[0]) if nsteps_override is None else nsteps_override
    dt = float(dtline.split("deltat ")[1].split()[0]) if dt_override is None else dt_override
    oldout = lines[0]
    prefix = name.replace("fig1_", "fig1_").replace("fig3_", "fig3_")
    outname = f"{prefix}_refined_n{n}.out"
    domain_min, domain_max = ((0.0, 1.0) if is_wang else (-1.0, 1.0))
    h = (domain_max - domain_min) / n
    nn, ne = (n + 1) ** 2, 2 * n * n
    nset = nn + (1 if is_wang else 2)
    text = [outname,
            f"{name} uniformly refined mesh n={n}, h={h:.8g}, eps={eps:g}",
            f"nltransienttransportproblem nsteps {nsteps} deltat {dt:.12g} alpha 1.0"
            + (f" bdforder {bdf_order}" if bdf_order > 1 else "")
            + (" manrmsteps 1" if bdf_order > 1 or full_newton else "")
            + " rtol 1.e-10 nsmax 40 nmodules 1",
            (f"vtkxml tstep_step {max(1, nsteps // 100)} domain_all primvars 1 6"
             if bdf_order > 1 else "vtkxml tstep_all domain_all primvars 1 6"),
            "domain heattransfer",
            f"OutputManager tstep_step {nsteps} dofman_all element_all",
            f"ndofman {nn} nelem {ne} ncrosssect 1 nmat 1 nbc {0 if is_wang else 1} nic {nn} nltf 1 nset {nset}"]
    def node_id(i, j):
        return j * (n + 1) + i + 1
    coords = {}
    for j in range(n + 1):
        y = domain_min + h * j
        for i in range(n + 1):
            x = domain_min + h * i
            idx = node_id(i, j)
            coords[idx] = (x, y)
            text.append(f"node {idx} coords 3 {x:.16g} {y:.16g} 0")
    e = 0
    for j in range(n):
        for i in range(n):
            bl, br = node_id(i, j), node_id(i + 1, j)
            tl, tr = node_id(i, j + 1), node_id(i + 1, j + 1)
            e += 1; text.append(f"tr1ht {e} nodes 3 {bl} {br} {tr}")
            e += 1; text.append(f"tr1ht {e} nodes 3 {bl} {tr} {tl}")
    text += ["SimpleTransportCS 1 thickness 1.0 area 1.0 mat 1 set 1",
             f"allencahn 1 d 1.0 k 1.0 c 1.0 mobility {1/eps**2:.16g}"]
    if not is_wang:
        boundary = [node_id(i, j) for j in range(n + 1) for i in range(n + 1)
                    if i in (0, n) or j in (0, n)]
        text.append(f"BoundaryCondition 1 loadTimeFunction 1 dofs 1 10 values 1 1.0 set {nn + 2}")
    for idx, (x, y) in coords.items():
        if is_wang:
            u = math.cos(2 * math.pi * x) * math.cos(2 * math.pi * y)
        else:
            lam = signed_ellipse_distance(x, y)
            u = math.tanh(lam / (math.sqrt(2) * eps))
        text.append(f"InitialCondition {idx} dofs 1 10 Conditions 1 u {u:.16g} set {idx + 1}")
    text.append("ConstantFunction 1 f(t) 1.0")
    text.append(f"Set 1 elementranges {{(1 {ne})}}")
    for idx in range(1, nn + 1):
        text.append(f"Set {idx + 1} nodes 1 {idx}")
    if not is_wang:
        text.append(f"Set {nn + 2} nodes {len(boundary)} " + " ".join(map(str, boundary)))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(text) + "\n")
    return outname, h, nsteps, dt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=40, help="number of grid divisions per side")
    parser.add_argument("--output", type=Path, default=ROOT / "refined_inputs_n40")
    parser.add_argument("--akrivis-paper", action="store_true",
                        help="use Fig. 3 h=epsilon/2 and tau=epsilon^2 with BDF5 (BDF startup, not Radau IIA)")
    parser.add_argument("--wang-paper", action="store_true",
                        help="use Fig. 1 mesh h=1/20 and delta-t=0.01 (bulk Allen-Cahn BE, not the paper's coupled scheme)")
    args = parser.parse_args()
    for template in sorted(INPUTS.glob("*.in")):
        group = "wang" if template.name.startswith("wang") else "akrivis"
        if args.akrivis_paper:
            if group != "akrivis":
                continue
            eps_digits = template.stem.split("eps0p")[1]
            eps = int(eps_digits) / (10 if len(eps_digits) == 1 else 100)
            n = round(4.0 / eps)  # domain [-1,1], h=epsilon/2
            dt = eps * eps
            nsteps = math.ceil(0.1 / dt)
            result = make_case(template, args.output / group / template.stem / "case.in",
                               n, dt, nsteps, bdf_order=5)
        elif args.wang_paper:
            if group != "wang":
                continue
            result = make_case(template, args.output / group / template.stem / "case.in",
                               20, 0.01, 100, full_newton=True)
        else:
            result = make_case(template, args.output / group / template.stem / "case.in", args.n)
        print(group, template.stem, "h=", result[1], "steps=", result[2], "dt=", result[3])


if __name__ == "__main__":
    main()
