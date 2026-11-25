#!/usr/bin/env python3
"""Post-process Allen--Cahn energy from OOFEM ASCII VTU/PVD output.

Initial scope: one AllenCahn material/cross-section, affine or bilinear
quad1ht cells (VTK_QUAD) and dynamic line1ht cells (VTK_LINE), ASCII VTU.
No solver run is started. If nodal Newmark velocity is absent, kinetic and
total hyperbolic energy are left blank rather than estimated from snapshots.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def data_array(parent: ET.Element, name: str | None = None) -> tuple[list[float], int]:
    for arr in parent.findall("DataArray"):
        if name is None or arr.attrib.get("Name") == name:
            if arr.attrib.get("format", "ascii").lower() != "ascii":
                raise ValueError("Only ASCII VTU arrays are supported in this first version")
            ncomp = int(arr.attrib.get("NumberOfComponents", "1"))
            return [float(v) for v in (arr.text or "").split()], ncomp
    raise ValueError(f"VTU DataArray {name!r} is missing")


def read_vtu(path: Path, field: str, velocity_field: str) -> dict:
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError(f"No VTU Piece found in {path}")
    xyz, ncoord = data_array(piece.find("./Points"))
    if ncoord != 3 or len(xyz) % 3:
        raise ValueError(f"Unexpected point coordinate array in {path}")
    points = [xyz[i:i + 3] for i in range(0, len(xyz), 3)]
    pdata = piece.find("./PointData")
    if pdata is None:
        raise ValueError(f"No PointData in {path}")
    vals, nc = data_array(pdata, field)
    if nc != 1 or len(vals) != len(points):
        raise ValueError(f"{field} must be a scalar at each node in {path}")
    velocity = None
    try:
        velocity, nv = data_array(pdata, velocity_field)
        if nv != 1 or len(velocity) != len(points):
            raise ValueError(f"{velocity_field} must be a scalar at each node")
    except ValueError as exc:
        if f"DataArray {velocity_field!r} is missing" not in str(exc):
            raise

    cells = piece.find("./Cells")
    conn, _ = data_array(cells, "connectivity")
    offsets, _ = data_array(cells, "offsets")
    types, _ = data_array(cells, "types")
    conn = [int(x) for x in conn]
    offsets = [int(x) for x in offsets]
    types = [int(x) for x in types]
    bulk, boundary = [], []
    start = 0
    for end, typ in zip(offsets, types):
        nodes = conn[start:end]
        start = end
        if typ == 9 and len(nodes) == 4:  # VTK_QUAD / quad1ht
            bulk.append(nodes)
        elif typ == 3 and len(nodes) == 2:  # VTK_LINE / line1ht
            boundary.append(nodes)
        else:
            raise ValueError(f"Unsupported VTK cell type {typ} in {path}; expected quad and line")
    if not bulk:
        raise ValueError(f"No quad1ht (VTK_QUAD) bulk cells in {path}")
    return {"points": points, "u": vals, "v": velocity, "bulk": bulk, "boundary": boundary}


def read_frames(input_path: Path, field: str, velocity_field: str) -> list[tuple[float, dict]]:
    first_line = next((line.strip() for line in input_path.read_text(errors="replace").splitlines() if line.strip()), "")
    if not first_line:
        raise ValueError(f"Cannot read OOFEM output base name from {input_path}")
    result_base = Path(first_line)
    if not result_base.is_absolute():
        result_base = input_path.parent / result_base
    pvd = Path(str(result_base) + ".m0.pvd")
    if not pvd.exists():
        # Single-frame fallback: input.out.m0.1.vtu is also common.
        candidates = sorted(result_base.parent.glob(result_base.name + ".m0.*.vtu"))
        if len(candidates) != 1:
            raise FileNotFoundError(f"Expected PVD {pvd} or exactly one matching VTU")
        match = re.search(r"\.m0\.([0-9.eE+-]+)\.vtu$", candidates[0].name)
        t = float(match.group(1)) if match else 0.0
        return [(t, read_vtu(candidates[0], field, velocity_field))]
    root = ET.parse(pvd).getroot()
    datasets = root.findall(".//DataSet")
    if not datasets:
        raise ValueError(f"No DataSet entries in {pvd}")
    frames = []
    for ds in datasets:
        vtu = (pvd.parent / ds.attrib["file"]).resolve()
        frames.append((float(ds.attrib["timestep"]), read_vtu(vtu, field, velocity_field)))
    return frames


def parse_record(line: str) -> dict[str, float]:
    # Parse scalar key/value pairs, ignoring the leading OOFEM record name/id.
    tokens = line.split()
    result = {}
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    for key, value in re.findall(rf"([A-Za-z_]\w*)\s+({number})(?=\s|$)", " ".join(tokens[2:])):
        result[key.lower()] = float(value)
    return result


def read_parameters(input_path: Path) -> dict:
    lines = input_path.read_text(errors="replace").splitlines()
    ac = [line for line in lines if re.match(r"^\s*allencahn\s+\d+\s", line, re.I)]
    cs = [line for line in lines if re.match(r"^\s*SimpleTransportCS\s+\d+\s", line, re.I)]
    problem = next((line for line in lines if re.match(r"^\s*nltransienttransportproblem\s", line, re.I)), "")
    if len(ac) != 1 or len(cs) != 1:
        raise ValueError("First version requires exactly one allencahn material and one SimpleTransportCS")
    mat = parse_record(ac[0])
    section = parse_record(cs[0])
    cfg = parse_record(problem)
    # Defaults mirror AllenCahnMaterial and NLTransientTransportProblem.
    material_hyperbolic = int(mat.get("hyperbolic", 0)) == 1
    wave_mode = int(cfg.get("wave", 0)) == 1
    if material_hyperbolic != wave_mode:
        raise ValueError("Inconsistent hyperbolic setup: material hyperbolic and problem wave flags must both be 1 or both be 0")
    if "thickness" not in section:
        raise ValueError("SimpleTransportCS must explicitly provide thickness for bulk energy integration")
    return {
        "k": mat.get("k", 1.0), "mobility": mat.get("mobility", 1.0),
        "surfacebeta": mat.get("surfacebeta", 1.0), "c": mat.get("c", 1.0),
        "rho": mat.get("rho", 1.0), "surfacecapacity": mat.get("surfacecapacity", 1.0),
        "hyperbolic": material_hyperbolic and wave_mode,
        "manufactured": int(mat.get("manufactured", 0)) == 1,
        "thickness": section.get("thickness", 1.0), "area": section.get("area", 0.0),
        "wave_damping": cfg.get("wavedamping", 0.0), "input_line": ac[0],
        "section_line": cs[0],
    }


def shape_quad(xi: float, eta: float) -> tuple[list[float], list[list[float]]]:
    n = [0.25 * (1-xi) * (1-eta), 0.25 * (1+xi) * (1-eta),
         0.25 * (1+xi) * (1+eta), 0.25 * (1-xi) * (1+eta)]
    dn = [[-0.25*(1-eta), -0.25*(1-xi)], [0.25*(1-eta), -0.25*(1+xi)],
          [0.25*(1+eta), 0.25*(1+xi)], [-0.25*(1+eta), 0.25*(1-xi)]]
    return n, dn


def integrate_energy(frame: dict, p: dict) -> dict[str, float | None]:
    # Three-point Gauss-Legendre is exact for the quartic potential on affine
    # Q1 quads and linear lines; for distorted quads it is a robust starting rule.
    qx = [-math.sqrt(3/5), 0.0, math.sqrt(3/5)]
    qw = [5/9, 8/9, 5/9]
    out = {"bulk_gradient": 0.0, "bulk_potential": 0.0,
           "surface_gradient": 0.0, "surface_potential": 0.0, "kinetic": None}
    pts, u = frame["points"], frame["u"]
    for ids in frame["bulk"]:
        xy = [[pts[i][0], pts[i][1]] for i in ids]
        ue = [u[i] for i in ids]
        for i, xi in enumerate(qx):
            for j, eta in enumerate(qx):
                n, dn = shape_quad(xi, eta)
                # A = dx/d(xi,eta), stored physical rows by natural columns.
                a00 = sum(xy[a][0] * dn[a][0] for a in range(4))
                a01 = sum(xy[a][0] * dn[a][1] for a in range(4))
                a10 = sum(xy[a][1] * dn[a][0] for a in range(4))
                a11 = sum(xy[a][1] * dn[a][1] for a in range(4))
                det = a00*a11 - a01*a10
                if abs(det) < 1e-15:
                    raise ValueError("Degenerate bulk element")
                inv = [[a11/det, -a01/det], [-a10/det, a00/det]]
                grad = [0.0, 0.0]
                uq = sum(n[a]*ue[a] for a in range(4))
                for a in range(4):
                    # grad_x = grad_nat * inverse(dx/dnat)
                    gx = dn[a][0]*inv[0][0] + dn[a][1]*inv[1][0]
                    gy = dn[a][0]*inv[0][1] + dn[a][1]*inv[1][1]
                    grad[0] += gx*ue[a]
                    grad[1] += gy*ue[a]
                w = abs(det) * qw[i]*qw[j] * p["thickness"]
                out["bulk_gradient"] += 0.5*p["k"]*(grad[0]**2+grad[1]**2)*w
                out["bulk_potential"] += 0.25*p["mobility"]*(uq*uq-1.0)**2*w

    for ids in frame["boundary"]:
        x0, x1 = pts[ids[0]], pts[ids[1]]
        dx, dy = x1[0]-x0[0], x1[1]-x0[1]
        jac = 0.5*math.hypot(dx, dy)
        if jac <= 0:
            raise ValueError("Degenerate surface line")
        ue = [u[ids[0]], u[ids[1]]]
        for xi, wi in zip(qx, qw):
            n = [0.5*(1-xi), 0.5*(1+xi)]
            uq = n[0]*ue[0] + n[1]*ue[1]
            grad_s = (ue[1]-ue[0])/(2*jac)
            w = jac*wi*p["area"]
            out["surface_gradient"] += 0.5*p["k"]*grad_s*grad_s*w
            out["surface_potential"] += 0.25*p["surfacebeta"]*p["mobility"]*(uq*uq-1.0)**2*w

    if p["hyperbolic"] and frame["v"] is not None:
        kinetic = 0.0
        v = frame["v"]
        for ids in frame["bulk"]:
            xy = [[pts[i][0], pts[i][1]] for i in ids]
            ve = [v[i] for i in ids]
            for i, xi in enumerate(qx):
                for j, eta in enumerate(qx):
                    n, dn = shape_quad(xi, eta)
                    a00 = sum(xy[a][0]*dn[a][0] for a in range(4))
                    a01 = sum(xy[a][0]*dn[a][1] for a in range(4))
                    a10 = sum(xy[a][1]*dn[a][0] for a in range(4))
                    a11 = sum(xy[a][1]*dn[a][1] for a in range(4))
                    det = a00*a11-a01*a10
                    vq = sum(n[a]*ve[a] for a in range(4))
                    kinetic += 0.5*p["rho"]*vq*vq*abs(det)*qw[i]*qw[j]*p["thickness"]
        out["kinetic"] = kinetic
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path, help="OOFEM .in file corresponding to exported VTU/PVD")
    ap.add_argument("--output", type=Path, help="CSV destination (must not already exist)")
    ap.add_argument("--field", default="Temperature", help="scalar primary PointData array")
    ap.add_argument("--velocity-field", default="Velocity", help="scalar nodal velocity PointData array")
    args = ap.parse_args()
    inp = args.input.resolve()
    dest = args.output.resolve() if args.output else inp.with_name(inp.stem + "_energy_postprocessed.csv")
    if dest.exists():
        ap.error(f"Refusing to overwrite {dest}; choose a new --output path")
    p = read_parameters(inp)
    if p["manufactured"]:
        print("WARNING: manufactured source is enabled; the unforced energy law in energy_theory.md does not apply.", file=sys.stderr)
    frames = read_frames(inp, args.field, args.velocity_field)
    first = frames[0][1]
    rows = []
    for time, frame in frames:
        if frame["points"] != first["points"] or frame["bulk"] != first["bulk"] or frame["boundary"] != first["boundary"]:
            raise ValueError("Mesh/connectivity changed between output frames")
        e = integrate_energy(frame, p)
        pot = e["bulk_gradient"] + e["bulk_potential"] + e["surface_gradient"] + e["surface_potential"]
        kin = e["kinetic"] if p["hyperbolic"] else 0.0
        total = pot + kin if kin is not None else None
        status = ("not_applicable" if not p["hyperbolic"] else
                  ("available" if kin is not None else "missing_velocity"))
        rows.append({"time": time, **e, "potential_energy": pot, "total_energy": total,
                     "kinetic_status": status})
    dest.parent.mkdir(parents=True, exist_ok=True)
    fields = ["time", "bulk_gradient", "bulk_potential", "surface_gradient", "surface_potential",
              "kinetic", "potential_energy", "total_energy", "kinetic_status"]
    with dest.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} energy frames to {dest}")
    print(f"Coefficients: k={p['k']:g}, mobility={p['mobility']:g}, surfacebeta={p['surfacebeta']:g}, "
          f"rho={p['rho']:g}, thickness={p['thickness']:g}, line area={p['area']:g}")
    if p["hyperbolic"] and first["v"] is None:
        print("Kinetic and total energies are blank: no scalar nodal Velocity array was exported.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
