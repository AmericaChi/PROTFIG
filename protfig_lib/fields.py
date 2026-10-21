"""
fields.py — 3D fields for coloring surfaces without external programs (no APBS):

  coulombic(sel)  ChimeraX-style Coulomb potential: phi = sum 332 q / (eps r),
                  eps = 4r  ->  phi = 83 q / r^2   [kcal/(mol·e)]
  lipophilic(sel) molecular lipophilicity potential (MLP, Fauchère-type):
                  sum f_i exp(-r), f_i per atom type

They are computed on a grid, written as OpenDX, and PyMOL colors the surface with a ramp
(surface_color = ramp). colorbar() draws the color bar separately with matplotlib.
"""
import os

import numpy as np
from pymol import cmd

from .core import ONE_LETTER

ONE_LETTER_SET = set(ONE_LETTER)

# simple formal charges per (resn, atom); neutral HIS (like ChimeraX without protonation)
CHARGES = {
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("LYS", "NZ"): 1.0,
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5,
}
ION_CHARGES = {"ZN": 2, "MG": 2, "CA": 2, "MN": 2, "FE": 2, "FE2": 2, "NI": 2, "CU": 2, "CO": 2, "CD": 2,
               "NA": 1, "K": 1, "CL": -1, "BR": -1, "IOD": -1}


def charges(sel):
    """(xyz, q) of the charged atoms in sel (simple charges + termini + ions + formal_charge)."""
    m = cmd.get_model(sel)
    xyz = np.array([a.coord for a in m.atom]) if m.atom else np.zeros((0, 3))
    q = np.zeros(len(m.atom))
    first, last = {}, {}  # per protein chain: resv of the first/last residue
    for a in m.atom:
        if a.resn in ONE_LETTER_SET:
            v = int(a.resi_number)
            first[a.chain] = min(first.get(a.chain, v), v)
            last[a.chain] = max(last.get(a.chain, v), v)
    has_oxt = {(a.chain, int(a.resi_number)) for a in m.atom if a.name == "OXT"}
    for i, a in enumerate(m.atom):
        if (a.resn, a.name) in CHARGES:
            q[i] = CHARGES[(a.resn, a.name)]
        elif a.resn in ION_CHARGES:
            q[i] = ION_CHARGES[a.resn]
        elif a.formal_charge:
            q[i] = a.formal_charge
        if a.resn in ONE_LETTER_SET:
            v = int(a.resi_number)
            if v == first[a.chain] and a.name == "N":
                q[i] += 1.0
            if v == last[a.chain]:
                if (a.chain, v) in has_oxt and a.name in ("O", "OXT"):
                    q[i] += -0.5
                elif (a.chain, v) not in has_oxt and a.name == "O":
                    q[i] += -1.0
    keep = q != 0
    return xyz[keep], q[keep]


def _grid(xyz, pad, spacing, max_points=4.0e6):
    lo = xyz.min(0) - pad
    hi = xyz.max(0) + pad
    while np.prod(np.ceil((hi - lo) / spacing) + 1) > max_points:
        spacing *= 1.15
    n = (np.ceil((hi - lo) / spacing) + 1).astype(int)
    return lo, spacing, n


def _splat(lo, h, n, xyz, w, kernel, cutoff):
    """Sum w_i * kernel(d) over grid points within cutoff of each atom."""
    field = np.zeros(int(np.prod(n)))
    k = int(np.ceil(cutoff / h))
    off = np.stack(np.meshgrid(*[np.arange(-k, k + 1)] * 3, indexing="ij"), -1).reshape(-1, 3)
    for start in range(0, len(xyz), 256):
        X = xyz[start:start + 256]
        W = w[start:start + 256]
        base = np.round((X - lo) / h).astype(int)  # (m,3)
        idx = base[:, None, :] + off[None, :, :]  # (m,K,3)
        pts = lo + idx * h
        d = np.linalg.norm(pts - X[:, None, :], axis=2)
        ok = (d <= cutoff) & np.all((idx >= 0) & (idx < n), axis=2)
        flat = np.ravel_multi_index(tuple(idx[ok].T), tuple(n))
        vals = (W[:, None] * kernel(d))[ok]
        field += np.bincount(flat, weights=vals, minlength=field.size)
    return field.reshape(n)


def write_dx(path, field, lo, h):
    nx, ny, nz = field.shape
    with open(path, "w") as fh:
        fh.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        fh.write(f"origin {lo[0]:.4f} {lo[1]:.4f} {lo[2]:.4f}\n")
        fh.write(f"delta {h:.4f} 0 0\ndelta 0 {h:.4f} 0\ndelta 0 0 {h:.4f}\n")
        fh.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        fh.write(f"object 3 class array type double rank 0 items {field.size} data follows\n")
        flat = field.ravel()  # x slow, z fast (OpenDX order)
        n3 = len(flat) // 3 * 3
        np.savetxt(fh, flat[:n3].reshape(-1, 3), fmt="%.4e")
        if n3 < len(flat):
            fh.write(" ".join(f"{v:.4e}" for v in flat[n3:]) + "\n")
        fh.write('attribute "dep" string "positions"\n')
        fh.write('object "regular positions regular connections" class field\n')
        fh.write('component "positions" value 1\ncomponent "connections" value 2\ncomponent "data" value 3\n')


def coulombic(sel, region_sel, path, spacing=1.0, cutoff=15.0):
    """Coulomb potential (eps = 4r) on a grid around region_sel; writes a .dx."""
    xyz, q = charges(sel)
    lo, h, n = _grid(cmd.get_coords(region_sel), 8.0, spacing)
    field = _splat(lo, h, n, xyz, q, lambda d: 83.0 / np.maximum(d, 2.0) ** 2, cutoff)
    write_dx(path, field, lo, h)
    return field


def lipophilic_weights(sel):
    """f_i per atom (simplified Fauchère/Broto-type scale): apolar C +, polar N/O -."""
    groups = [
        (f"({sel}) and elem C and not (neighbor elem N+O)", 0.45),
        (f"({sel}) and elem C and (neighbor elem N+O)", -0.10),
        (f"({sel}) and elem S", 0.55),
        (f"({sel}) and elem Cl+Br+I+F", 0.60),
        (f"({sel}) and elem N", -0.75),
        (f"({sel}) and elem O", -0.65),
        (f"({sel}) and ((resn LYS and name NZ) or (resn ARG and name NH1+NH2+NE) or (resn ASP and name OD1+OD2) "
         f"or (resn GLU and name OE1+OE2))", -0.40),  # extra for charged groups
    ]
    xyz, w = [], []
    for s, f in groups:
        c = cmd.get_coords(s)
        if c is not None:
            xyz.append(c)
            w.append(np.full(len(c), f))
    return np.vstack(xyz), np.concatenate(w)


def surface_probe_points(surf_obj, tmp_obj, probe=1.4, max_pts=60000):
    """Points where PyMOL evaluates the map to color the surface (surface_ramp_above_mode=1):
    surface vertices shifted probe Å along their normal. Obtained by exporting only surf_obj
    to .obj (the only way to read a surface's vertices through the API)."""
    on = cmd.get_names("objects", enabled_only=1)
    for o in on:
        if o != surf_obj:
            cmd.disable(o)
    cmd.save(tmp_obj)
    for o in on:
        cmd.enable(o)
    v, n = [], []
    with open(tmp_obj) as fh:
        for line in fh:
            if line.startswith("v "):
                v.append(line.split()[1:4])
            elif line.startswith("vn "):
                n.append(line.split()[1:4])
    os.remove(tmp_obj)
    v, n = np.array(v, float), np.array(n, float)
    step = max(1, len(v) // max_pts)
    return v[::step] + probe * n[::step]


def sample(field, lo, h, pts):
    idx = np.round((pts - lo) / h).astype(int)
    ok = np.all((idx >= 0) & (idx < field.shape), axis=1)
    return field[tuple(idx[ok].T)]


def lipophilic(sel, region_sel, path, spacing=1.0, cutoff=5.0):
    xyz, w = lipophilic_weights(sel)
    lo, h, n = _grid(cmd.get_coords(region_sel), 8.0, spacing)
    field = _splat(lo, h, n, xyz, w, lambda d: np.exp(-d), cutoff)
    write_dx(path, field, lo, h)
    return field, lo, h


def auto_levels(field, lo, h, surf_obj, tmp_obj):
    """(p10, median, p90) of the field over the surface: the ramp is centered on the median (otherwise
    a mostly polar protein comes out a single color) and uses the range that is actually visible."""
    vals = sample(field, lo, h, surface_probe_points(surf_obj, tmp_obj))
    if not len(vals):
        return (-0.5, 0.0, 0.5)
    return tuple(float(np.percentile(vals, q)) for q in (10, 50, 90))


def surface_ramp(surf_obj, map_path, map_name, ramp_name, levels, colors):
    """Color the surface of surf_obj by the map value evaluated just outside it."""
    cmd.load(map_path, map_name)
    os.remove(map_path)  # the map is kept inside the .pse
    # ramp_new with a list of RGB colors segfaults in PyMOL 3.1: use named colors instead
    names = []
    for i, c in enumerate(colors):
        names.append(f"{ramp_name}_c{i}")
        cmd.set_color(names[-1], [float(x) for x in c])
    cmd.ramp_new(ramp_name, map_name, [float(v) for v in levels], names)
    cmd.set("surface_ramp_above_mode", 1)
    cmd.set("surface_color", ramp_name, surf_obj)
    cmd.disable(ramp_name)
    cmd.disable(map_name)


def colorbar(path, levels, colors, label, ticks=None, ticklabels=None):
    """Horizontal color bar (png + svg) for assembling the figure; same colors as the ramp."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    lo, hi = levels[0], levels[-1]
    pos = [(v - lo) / (hi - lo) for v in levels]
    cm = LinearSegmentedColormap.from_list("mf", list(zip(pos, colors)))
    fig, ax = plt.subplots(figsize=(3.2, 0.62), dpi=300)
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.55, top=0.92)
    ax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=cm, extent=[lo, hi, 0, 1])
    ax.set_yticks([])
    ax.set_xticks(ticks if ticks is not None else [lo, hi])
    ax.set_xticklabels(ticklabels or [f"{t:.2g}" for t in (ticks if ticks is not None else [lo, hi])])
    ax.tick_params(labelsize=8, length=2, pad=1)
    ax.set_xlabel(label, fontsize=8, labelpad=1)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
    for ext in ("png", "svg"):
        fig.savefig(f"{path}.{ext}", transparent=True)
    plt.close(fig)
