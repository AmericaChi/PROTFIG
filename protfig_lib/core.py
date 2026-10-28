"""
core.py — shared PROTFIG pieces: loading (file, PDB ID or user .pse), cleanup,
binding-partner detection (ligand / peptide / protein), per-chain colors, automatic views
and labels.
"""
import math
import os
import re
import sys

import numpy as np
from pymol import cmd

class ProtFigError(Exception):
    """User error (missing file, empty selection, style without a partner...): shown
    as 'ERROR: message', without a traceback."""


ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "SEP": "S", "TPO": "T", "PTR": "Y",
}
# crystallization / buffer additives and simulation solvent/ions (GROMACS, CHARMM, AMBER):
# removed unless they are the chosen partner
JUNK_HET = ("DMS GOL EDO PEG PG4 PGE 1PE P6G SO4 PO4 ACT ACY FMT TRS EPE MES BME MPD IPA EOH NO3 "
            "SCN NH4 IMD CIT TLA MLI PGO BU3 DIO MRD SIN PEU 12P 15P P33 2PE 7PE XPE IOD BR CL NA K "
            "SOL TIP3 TIP4 TP3 T3P SPC WAT HOH DOD SOD CLA POT NA+ CL- K+").split()
METALS = "ZN MG CA MN FE FE2 NI CU CO CD".split()
# lipids (CHARMM-GUI, MemProtMD, atomistic Martini...). Only the membrane preset keeps them.
LIPIDS = ("DPP DPPC POPC POPE POPG POPS POPA DOPC DOPE DOPG DOPS DMPC DMPE DLPC DLPE DSPC PLPC "
          "SDPC SOPC CHL1 CHOL CLR CHS PSM SSM PIP2 PI2 SAPI POPI CDL CDL2 LMPG DPC LDA OLC").split()
# chains with up to this many residues are treated as peptides (e.g. GLP-1 analogs, ~30-40 aa)
PEPTIDE_MAX = 50
# glycans: never picked as the automatic "ligand"
GLYCANS = "NAG NDG BMA MAN FUC GAL GLC BGC SIA NGA FUL XYS A2G".split()

# (disabled) object holding the cleaned structure; saved in the .pse so you can re-run
# with --session <prefix>.pse after reorienting or editing selections by hand in PyMOL
BASE = "protfig_input"
# selections the user can define in their .pse that PROTFIG recognizes
USER_SELECTIONS = ("partner", "ligand", "site", "highlight", "key", "mutation", "focus")

BASE_SETTINGS = {
    "bg_rgb": "white",
    "ambient": 0.5,
    "spec_reflect": 0.0,
    "specular": 0.0,
    "direct": 0.45,
    "reflect": 0.45,
    "light_count": 2,
    "ray_shadows": 0,
    "depth_cue": 1,
    "fog_start": 0.55,
    "transparency_mode": 2,
    "ray_trace_mode": 0,
    "ray_trace_gain": 0.1,
    "ray_trace_color": "black",
    "orthoscopic": 0,
    "field_of_view": 20,
    "surface_quality": 1,
    "antialias": 2,
    "cartoon_smooth_loops": 1,
    "cartoon_flat_sheets": 1,
    "cartoon_fancy_helices": 0,
    "cartoon_side_chain_helper": 1,
    "cartoon_loop_radius": 0.25,
    "cartoon_transparency": 0.0,
    "stick_radius": 0.22,
    "sphere_scale": 1.0,
    "valence": 0,
    "label_font_id": 5,
    "label_color": "black",
    "float_labels": 0,
    "label_position": (0.0, 0.0, 3.0),  # in front of the atom: not hidden by its own stick
    "dash_color": "black",
    "dash_gap": 0.35,
    "dash_length": 0.35,
    "dash_radius": 0.07,
    "dash_round_ends": 1,
    "mesh_width": 1.0,
    "mesh_radius": 0.016,
    "ray_interior_color": "grey85",  # faces cut by the clipping plane: light gray, not black
}

# per-chain palettes (in --chains order)
PALETTES = {
    "vivid": [(0.18, 0.82, 0.18), (0.62, 0.05, 0.62), (0.07, 0.12, 0.85), (0.95, 0.55, 0.10),
              (0.10, 0.70, 0.80), (0.85, 0.20, 0.20), (0.95, 0.85, 0.15), (0.50, 0.50, 0.50)],
    # soft ChimeraX-style tones (sage green, teal, pink, sky blue, sand)
    "chimerax": [(0.40, 0.60, 0.40), (0.20, 0.50, 0.53), (0.78, 0.58, 0.58), (0.65, 0.85, 0.90),
                 (0.88, 0.82, 0.58), (0.80, 0.60, 0.40), (0.60, 0.60, 0.80), (0.70, 0.70, 0.70)],
    # muted (lavender, dusty rose, sage green, orange)
    "muted": [(0.53, 0.53, 0.83), (0.76, 0.45, 0.45), (0.55, 0.76, 0.55), (0.93, 0.74, 0.45),
              (0.55, 0.75, 0.85), (0.80, 0.65, 0.80), (0.75, 0.75, 0.55), (0.70, 0.70, 0.70)],
    "pastel": [(0.70, 0.88, 0.75), (0.95, 0.85, 0.70), (0.85, 0.80, 0.92), (0.75, 0.82, 0.92),
               (0.95, 0.80, 0.82), (0.90, 0.92, 0.75), (0.80, 0.90, 0.90), (0.85, 0.85, 0.85)],
    "illustrative": [(1.00, 0.94, 0.62), (0.70, 0.84, 0.60), (0.95, 0.75, 0.85), (0.70, 0.85, 0.95),
                     (0.95, 0.80, 0.60), (0.85, 0.85, 0.85)],
    # alternating light/dark grays (oligomers)
    "grays": [(0.86, 0.86, 0.86), (0.42, 0.42, 0.42)],
    # membrane protein: blue, then tones that contrast with the gray bilayer
    "membrane": [(0.45, 0.58, 0.88), (0.55, 0.78, 0.55), (0.85, 0.55, 0.75), (0.9, 0.75, 0.35)],
}


# ---------------------------------------------------------------- colors
def rgb_of(c):
    """'red', '0.2 0.3 0.4', '#33aa55' or tuple -> (r, g, b) in 0..1."""
    if isinstance(c, (tuple, list)):
        return tuple(float(x) for x in c)
    c = c.strip()
    if c.startswith("#") and len(c) == 7:
        return tuple(int(c[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
    parts = re.split(r"[ ,]+", c)
    if len(parts) == 3:
        return tuple(float(x) for x in parts)
    return tuple(cmd.get_color_tuple(cmd.get_color_index(c)))


def mix(c, other, t):
    """blend c with other: t=0 -> c, t=1 -> other."""
    c, o = rgb_of(c), rgb_of(other)
    return tuple((1 - t) * x + t * y for x, y in zip(c, o))


def define(name, c):
    cmd.set_color(name, list(rgb_of(c)))
    return name


def color_by_element(sel, carbon):
    """carbons = carbon, the rest by element (N blue, O red, S yellow, ...)."""
    cmd.util.cbag(sel)
    cmd.color(carbon, f"({sel}) and elem C")


# ---------------------------------------------------------------- loading
def _escape_resi(resi):
    return resi.replace("-", "\\-")


def residue_sel(chain, resi):
    return f"chain {chain} and resi {_escape_resi(str(resi))}"


def residues(sel):
    """Sorted [(chain, resi, resn)] with a CA (or any atom for non-protein residues)."""
    out = {}
    cmd.iterate(sel, "out[(chain, resi, resv)] = resn", space={"out": out})
    return [(c, r, n) for (c, r, v), n in sorted(out.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][1]))]


def load_input(a):
    """Load a.input (file / PDB ID / .pse) and optionally a.session (.pse).

    Returns (session_view, user_sels) and leaves the cleaned structure in the BASE object.
    user_sels: {name: selection on BASE} for the USER_SELECTIONS present in the .pse.
    """
    session = a.session or (a.input if a.input.lower().endswith(".pse") else None)
    structure = None if (a.input.lower().endswith(".pse")) else a.input
    view, user = None, {}

    if session:
        cmd.load(session)
        view = cmd.get_view()
        names = cmd.get_names("objects")
        if structure:
            _load_structure(structure, "_mf_new", os.path.dirname(os.path.abspath(a.prefix)))
            src = "_mf_new"
        elif BASE in names:
            src = BASE
        else:
            mols = [n for n in cmd.get_names("objects", enabled_only=1) if cmd.get_type(n) == "object:molecule"]
            if not mols:
                raise ProtFigError(f"{session} has no enabled molecular objects")
            src = " or ".join(mols)
        cmd.create("_mf_tmp", f"({src})")
        sels = cmd.get_names("selections")
        # chain order from a previous PROTFIG run (empty selection protfig_chains_B_C_A)
        tag = [n for n in sels if n.startswith("protfig_chains_")]
        chain_tag = tag[0][len("protfig_chains_"):].replace("_", ",") if tag else None
        for s in USER_SELECTIONS:
            if s in sels and cmd.count_atoms(s):
                # 'in' compares identifiers (chain/resi/name...), so it works even if the
                # selection was made on another object in the .pse (e.g. LIGAND or SITE)
                cmd.select(f"_mf_user_{s}", f"_mf_tmp in ({s})")
        for n in cmd.get_names("objects"):
            if n != "_mf_tmp":
                cmd.delete(n)
        for n in cmd.get_names("selections"):
            if not n.startswith("_mf_user_"):
                cmd.delete(n)
        cmd.set_name("_mf_tmp", BASE)
        for s in USER_SELECTIONS:
            if f"_mf_user_{s}" in cmd.get_names("selections"):
                user[s] = f"_mf_user_{s}"
        if chain_tag:
            user["_chains"] = chain_tag
        # settings from the user's .pse are not inherited: the preset defines the style
        view_keep = view
        _reset_settings_keep_objects()
        view = view_keep
    else:
        _load_structure(structure, BASE, os.path.dirname(os.path.abspath(a.prefix)))
    return view, user


def _reset_settings_keep_objects():
    # cmd.reinitialize("settings") resets all global settings to their defaults without deleting objects
    cmd.reinitialize("settings")


def _load_structure(src, name, fetch_dir="."):
    if os.path.exists(src):
        cmd.load(src, name)
    elif re.fullmatch(r"[0-9][A-Za-z0-9]{3}", src):
        cmd.set("fetch_path", os.path.abspath(fetch_dir))  # the downloaded .cif goes next to the outputs
        cmd.fetch(src, name, type="cif")
    else:
        raise ProtFigError(f"{src} does not exist (and is not a PDB ID)")


def clean(a, user, keep_lipids=False):
    """Remove H, altlocs, water (unless --waters), additives and MD solvent/ions (unless
    --keep_het), lipids (unless keep_lipids) and unrequested chains. Assign chain 'A' to
    chainless polymers (MD outputs)."""
    if not getattr(a, "keep_het", False):
        cmd.remove(f"{BASE} and resn {'+'.join(r for r in JUNK_HET if r not in ('HOH', 'WAT', 'DOD'))}"
                   + _keep_partner(a, user))
    cmd.remove(f"{BASE} and hydro")
    cmd.remove(f"{BASE} and not alt ''+A")
    cmd.alter(BASE, "alt=''")
    if cmd.count_atoms(f'{BASE} and polymer and chain ""'):
        cmd.alter(f'{BASE} and chain ""', "chain = segi.strip()[:1] or 'A'")
    if a.chains:
        chains = [c.strip() for c in a.chains.split(",") if c.strip()]
        present = cmd.get_chains(BASE)
        missing = [c for c in chains if c not in present]
        if missing:
            raise ProtFigError(f"--chains: chain {', '.join(missing)} is not in the structure. " + chain_help(missing))
        cmd.remove(f"{BASE} and not chain {'+'.join(chains)} and not resn {'+'.join(LIPIDS)}")
    if not a.waters:
        cmd.remove(f"{BASE} and solvent")
    if not keep_lipids:
        cmd.remove(f"{BASE} and resn {'+'.join(LIPIDS)}" + _keep_partner(a, user))
    if cmd.count_atoms(f"{BASE} and polymer.protein and ss H+S") == 0:
        cmd.dss(BASE)
    if cmd.count_atoms(f"{BASE} and polymer") == 0:
        raise ProtFigError("no protein/nucleic acid left in the structure after cleanup (is --chains correct?)")


def chain_help(missing):
    """Lists the chains present; for mmCIF files, also says which author chain a missing ID is the
    label (PDB website) name of, e.g. 7NRE: label chain B = author chain C."""
    msg = f"Available chains: {', '.join(cmd.get_chains(BASE))}."
    for c in missing:
        auth = cmd.get_chains(f"{BASE} and polymer and segi {c}")
        if auth and c not in auth:
            msg += (f"\n  Label chain {c} (the name some PDB pages show) is author chain {auth[0]} in this"
                    f" file; PROTFIG uses author chains, so use {auth[0]}.")
    return msg


def _keep_partner(a, user):
    keep = "".join(f" and not ({user[k]})" for k in ("partner", "ligand") if k in user)
    if a.partner:
        keep += f" and not ({BASE} and ({a.partner}))"
    return keep


def chain_order(a, user=None):
    present = cmd.get_chains(f"{BASE} and polymer")
    given = a.chains or (user or {}).get("_chains")
    if given:
        order = [c.strip() for c in given.split(",") if c.strip() in present]
        return order + [c for c in present if c not in order]
    return present


def partner_kind(sel):
    """'small' (small molecule), 'peptide' (<= PEPTIDE_MAX residues) or 'protein'."""
    n = cmd.count_atoms(f"({sel}) and polymer and name CA+P")
    if n == 0:
        return "small"
    return "peptide" if n <= PEPTIDE_MAX else "protein"


def find_partner(a, user, mode="required"):
    """Binding partner on BASE -> ('ligand' selection, kind) or (None, None).

    Priority: --partner > 'partner'/'ligand' selection in the .pse > automatic. Automatic:
      1) the largest organic residue (excluding additives, lipids and glycans) -> 'small'
      2) if mode == 'required': a short chain (<= PEPTIDE_MAX residues) -> 'peptide';
         otherwise, if there are >= 2 chains, the smallest chain -> 'protein'
    mode 'optional' only looks for small molecules; 'none' does no automatic search.
    """
    if a.partner:
        sel = f"{BASE} and ({a.partner})"
    elif "partner" in user:
        sel = user["partner"]
    elif "ligand" in user:
        sel = user["ligand"]
    elif mode == "none":
        return None, None
    else:
        sel = None
        counts = {}
        extra = f" and chain {a.partner_chain}" if a.partner_chain else ""
        skip = "+".join(JUNK_HET + LIPIDS + GLYCANS)
        cmd.iterate(f"{BASE} and organic and not resn {skip}{extra}",
                    "counts[(chain, resi, resn)] = counts.get((chain, resi, resn), 0) + 1",
                    space={"counts": counts})
        if counts:
            (c, r, n), _ = max(counts.items(), key=lambda kv: (kv[1], -ord(kv[0][0][:1] or "A")))
            sel = f"{BASE} and {residue_sel(c, r)} and resn {n}"
        elif mode == "required":
            sizes = []
            for c in cmd.get_chains(f"{BASE} and polymer"):
                n = cmd.count_atoms(f"{BASE} and chain {c} and polymer and name CA+P")
                if n:
                    sizes.append((n, c))
            if len(sizes) >= 2 or (sizes and min(sizes)[0] <= PEPTIDE_MAX):
                sel = f"{BASE} and chain {min(sizes)[1]} and polymer"
        if sel is None:
            return None, None
    if cmd.count_atoms(sel) == 0:
        what = a.partner or "'partner' selection in the .pse"
        named = re.findall(r"\bchain\s+([\w+]+)", a.partner or "")
        missing = [c for tok in named for c in tok.split("+") if not cmd.count_atoms(f"{BASE} and chain {c}")]
        hint = chain_help(missing) if missing else "Is the residue name or chain correct?"
        raise ProtFigError(f'--partner "{what}" selects no atoms. {hint}')
    cmd.select("ligand", sel)
    kind = partner_kind("ligand")
    if cmd.count_atoms(f"{BASE} and polymer and not ligand") == 0:
        raise ProtFigError("the partner covers the whole structure; specify another one with --partner")
    # interface core: the whole ligand/peptide, or the partner protein residues within
    # --site_cutoff of the receptor (same criterion as the receptor-side site)
    if kind == "protein":
        cmd.select("lig_core", f"byres (ligand and (all within {a.site_cutoff} of ({BASE} and polymer and not ligand)))")
        if cmd.count_atoms("lig_core") == 0:
            raise ProtFigError("the protein partner does not touch the rest of the structure (no interface)")
    else:
        cmd.select("lig_core", "ligand")
    return "ligand", kind


# ---------------------------------------------------------------- views
def _rot(z, up):
    z = z / np.linalg.norm(z)
    up = up - np.dot(up, z) * z
    if np.linalg.norm(up) < 1e-6:
        up = np.cross(z, [1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.cross(z, [0.0, 1.0, 0.0])
    up = up / np.linalg.norm(up)
    x = np.cross(up, z)
    x /= np.linalg.norm(x)
    return np.vstack([x, np.cross(z, x), z])


def _pca(xyz):
    c = xyz.mean(axis=0)
    w, v = np.linalg.eigh(np.cov((xyz - c).T))
    return c, v[:, ::-1]  # columns: long, middle, short axis


def adjust(R, a):
    """--flip / --roll / --turn 'x:30,y:-20' applied in camera axes before framing."""
    turns = []
    if a.flip:
        turns.append(("z", 180.0))
    if a.roll:
        turns.append(("z", a.roll))
    for tok in (a.turn or "").split(","):
        if ":" in tok:
            ax, ang = tok.split(":")
            turns.append((ax.strip().lower(), float(ang)))
    for ax, ang in turns:
        t = math.radians(ang)
        c, s = math.cos(t), math.sin(t)
        M = {"x": [[1, 0, 0], [0, c, -s], [0, s, c]],
             "y": [[c, 0, s], [0, 1, 0], [-s, 0, c]],
             "z": [[c, -s, 0], [s, c, 0], [0, 0, 1]]}[ax]
        R = np.array(M) @ R
    return R


def frame(R, center, frame_xyz, aspect, margin=1.06, front_z=None, fit_box=False):
    """PyMOL view with rotation R (rows = camera axes), centered on center, framing frame_xyz.

    front_z: if given, front clipping plane at that camera z coordinate (relative to center);
    anything closer to the camera is clipped (to look inside a pocket).
    fit_box: center on the on-screen bounding box instead (long flexible tails no longer push
    the molecule off-center).
    """
    cam = (frame_xyz - center) @ R.T
    if fit_box:
        mid = (cam.min(0) + cam.max(0)) / 2
        mid[2] = 0.0
        center = center + mid @ R
        cam = cam - mid
    half = math.radians(cmd.get_setting_float("field_of_view") / 2.0)
    zmax = np.max(cam[:, 2]) if front_z is None else min(front_z, np.max(cam[:, 2]))
    need = max(np.max(np.abs(cam[:, 1])), np.max(np.abs(cam[:, 0])) / aspect)
    dist = margin * need / math.tan(half) + zmax
    front = max(1.0, dist - (front_z if front_z is not None else np.max(cam[:, 2]) + 5))
    back = dist - np.min(cam[:, 2]) + 5
    return tuple(R.T.flatten().tolist()) + (0.0, 0.0, -float(dist)) + tuple(center.tolist()) + (front, back, -20.0)


def overview_R(sel, aspect, toward=None, up_sel=None):
    """Short axis toward the camera (most spread-out view); long axis along the image's long side.

    toward: selection that must face the camera (e.g. the ligand).
    up_sel: selection that must end up on top (or on the left if the image is landscape).
    """
    xyz = cmd.get_coords(sel)
    c, V = _pca(xyz)
    z = V[:, 2]
    tw = cmd.get_coords(toward) if toward else None
    if tw is not None and np.dot(tw.mean(0) - c, z) < 0:
        z = -z
    up = V[:, 0] if aspect < 1.15 else V[:, 1]
    ups = cmd.get_coords(up_sel) if up_sel else None
    if ups is not None:
        ref = ups.mean(0) - c
        if aspect >= 1.15:  # landscape image: the long axis is x; up_sel on the left
            if np.dot(ref, np.cross(up, z)) > 0:
                up = -up
        elif np.dot(ref, up) < 0:
            up = -up
    return c, _rot(z, up)


def site_R(lig, env_sel, aspect, up_sel=None, side=False):
    """Camera looking at the ligand from the pocket mouth: z = ligand centroid - environment centroid.

    side=True (protein-protein interfaces): side view of the interface, with the interface normal
    (receptor -> partner) vertical and the camera on the widest side of the interface.
    """
    L = cmd.get_coords(lig)
    lc = L.mean(axis=0)
    env = cmd.get_coords(env_sel)
    if side:
        n = lc - env.mean(axis=0)
        n /= np.linalg.norm(n)
        both = np.vstack([L, env])
        _, V = _pca(both)
        x = V[:, 0] - np.dot(V[:, 0], n) * n
        z = np.cross(x / np.linalg.norm(x), n)
        if aspect < 1.0:
            return (lc + env.mean(axis=0)) / 2, _rot(z, x)
        return (lc + env.mean(axis=0)) / 2, _rot(z, n)
    z = lc - env.mean(axis=0)
    if np.linalg.norm(z) < 1.0:
        z = _pca(env)[1][:, 2]
    z /= np.linalg.norm(z)
    # ligand long axis along the image's long side; the camera is placed perpendicular to that
    # axis (otherwise the ligand is foreshortened), keeping the pocket direction as far as possible
    lax = _pca(L)[1][:, 0] if len(L) > 2 else np.array([1.0, 0.0, 0.0])
    zp = z - np.dot(z, lax) * lax
    if np.linalg.norm(zp) > 0.2:
        z = zp / np.linalg.norm(zp)
    lax = lax - np.dot(lax, z) * z
    if aspect >= 1.0:
        up = np.cross(z, lax)
    else:
        up = lax
    ups = cmd.get_coords(up_sel) if up_sel else None
    if ups is not None and np.dot(ups.mean(0) - lc, up) < 0:
        up = -up
    return lc, _rot(z, up)


def facing_R(rec_sel, pocket_sel, aspect, up_sel=None):
    """Pocket facing the camera (camera along receptor center -> pocket centroid), receptor long
    axis vertical and up_sel (default: first 20 residues) on top."""
    xyz = cmd.get_coords(rec_sel)
    c = xyz.mean(axis=0)
    z = cmd.get_coords(pocket_sel).mean(axis=0) - c
    z /= np.linalg.norm(z)
    _, V = _pca(xyz)
    if abs(np.dot(z, V[:, 0])) > 0.7:
        # pocket at one end of the long axis (e.g. a peptide entering a GPCR from the top):
        # side view, with the long axis vertical and the pocket at top/bottom, not end-on
        long_ = V[:, 0] * np.sign(np.dot(z, V[:, 0]))
        side = z - np.dot(z, long_) * long_
        z = side / np.linalg.norm(side) if np.linalg.norm(side) > 0.15 else V[:, 1]
        return c, _rot(z, long_)
    up = next((V[:, i] for i in range(3) if np.linalg.norm(V[:, i] - np.dot(V[:, i], z) * z) > 0.3), V[:, 1])
    if not up_sel:
        # first 20 residues (explicit list: with negative numbering "resi -8-11" is ambiguous)
        first = residues(f"({rec_sel}) and name CA")[:20]
        up_sel = " or ".join(f"(({rec_sel}) and {residue_sel(c, r)})" for c, r, _ in first) or None
    top = cmd.get_coords(up_sel) if up_sel else None
    if top is not None and np.dot(top.mean(0) - c, up) < 0:
        up = -up
    return c, _rot(z, up)


def membrane_R(protein_sel, lipid_sel, up_sel=None):
    """Side view of a bilayer: membrane normal (lowest-variance axis of the phosphorus atoms)
    vertical, camera looking at the protein from its widest side; up_sel on top."""
    heads = f"({lipid_sel}) and elem P"
    L = cmd.get_coords(heads) if cmd.count_atoms(heads) else cmd.get_coords(lipid_sel)
    _, V = _pca(L)
    normal = V[:, 2]
    P = cmd.get_coords(protein_sel)
    c = P.mean(axis=0)
    inplane = P - c - np.outer((P - c) @ normal, normal)
    _, W = _pca(inplane + c)
    x = W[:, 0] - np.dot(W[:, 0], normal) * normal
    z = np.cross(x / np.linalg.norm(x), normal)
    ups = cmd.get_coords(up_sel) if up_sel else None
    if ups is not None and np.dot(ups.mean(0) - c, normal) < 0:
        normal = -normal
        z = -z
    return c, _rot(z, normal)


def vis_height(view):
    """visible height (Å) at the plane of the origin."""
    return 2.0 * abs(view[11]) * math.tan(math.radians(cmd.get_setting_float("field_of_view") / 2.0))


# ---------------------------------------------------------------- labels
def _declutter(P, texts, R, size, avoid=None, iters=300):
    """Separate labels that overlap in the image plane (boxes ~ text width) and move them off
    the atoms in avoid (e.g. the ligand), pushing them radially outward."""
    S = (P @ R.T)[:, :2].copy()
    w = np.array([0.62 * size * len(t) + 0.3 * size for t in texts])
    h = 1.25 * size
    A = (avoid @ R.T)[:, :2] if avoid is not None and len(avoid) else None
    ac = A.mean(0) if A is not None else None
    for _ in range(iters):
        moved = False
        if A is not None:
            for i in range(len(S)):
                d = np.abs(A - S[i])
                if np.any((d[:, 0] < w[i] / 2 + 0.9) & (d[:, 1] < h / 2 + 0.9)):
                    v = S[i] - ac
                    S[i] += 0.4 * v / (np.linalg.norm(v) or 1.0)
                    moved = True
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                d = S[j] - S[i]
                ox = (w[i] + w[j]) / 2 - abs(d[0])
                oy = h - abs(d[1])
                if ox > 0 and oy > 0:
                    moved = True
                    if oy < ox:
                        s = (oy / 2 + 0.05) * (1 if d[1] >= 0 else -1)
                        S[i, 1] -= s
                        S[j, 1] += s
                    else:
                        s = (ox / 2 + 0.05) * (1 if d[0] >= 0 else -1)
                        S[i, 0] -= s
                        S[j, 0] += s
        if not moved:
            break
    cam = P @ R.T
    cam[:, :2] = S
    return cam @ R


# width/height of the final image (set by protfig.py); the headless viewport doesn't reflect it
ASPECT = 1.0
# user label options (set by protfig.py): size scale, font, single color
# (None = the preset's) and residue label text format
LABEL_OPTS = {"scale": 1.0, "font": None, "color": None, "style": "one", "dark_bg": False}
FONTS = {"sans": 5, "sans-italic": 6, "sans-bold": 7, "sans-bold-italic": 8, "serif": 9, "serif-bold": 10,
         "mono": 11, "mono-bold": 13}
LABEL_FORMATS = {"one": "{one}{resi}", "three": "{resn}{resi}", "chain": "{chain}:{one}{resi}",
                 "chain-three": "{chain}:{resn}{resi}"}

# negative label_size = height in Å, but the visible capital-letter height is ~1/3 of that
LABEL_SCALE = 2.6
# SVG font-size per Å of label_size (calibrated against the PyMOL render)
SVG_EM = 0.5


def add_labels(obj, items, R, size, font=7, outline=None, view=None, avoid=None, aspect=None):
    """items: [(text, pos_xyz, color)]. Labels as pseudoatoms (one per item) in object obj.

    size: approximate text height in Å (at the site plane). With view, labels are moved
    forward to just behind the front clipping plane (never hidden by sticks or surfaces),
    correcting x,y for perspective so they project onto the same spot.
    """
    if not items:
        return
    size *= LABEL_OPTS["scale"]
    font = FONTS.get(LABEL_OPTS["font"], font) if LABEL_OPTS["font"] else font
    if LABEL_OPTS["color"]:
        items = [(t, p, rgb_of(LABEL_OPTS["color"])) for t, p, _ in items]
    elif LABEL_OPTS.get("dark_bg"):
        # dark background: the preset's dark label colors are lightened so they stay readable
        items = [(t, p, mix(c, "white", 0.75) if sum(rgb_of(c)) / 3 < 0.45 else c) for t, p, c in items]
    P = np.array([it[1] for it in items], dtype=float)
    if view is not None and avoid is not None:
        avoid = avoid - np.array(view[12:15])  # same relative coords as P below
    origin0 = np.array(view[12:15]) if view is not None else np.zeros(3)
    P = _declutter(P - origin0, [it[0] for it in items], R, size, avoid=avoid) + origin0
    if view is not None:
        dist, origin = -view[11], np.array(view[12:15])
        cam = (P - origin) @ R.T
        # label_position pushes the text 3 Å further toward the camera: leave a margin to the clipping plane
        zl = dist - view[15] - 4.2
        cam[:, :2] *= ((dist - zl) / (dist - cam[:, 2]))[:, None]
        cam[:, 2] = zl
        # inside the frame: the text (centered on the point) must not spill over the edges
        half_h = (dist - zl) * math.tan(math.radians(cmd.get_setting_float("field_of_view") / 2.0))
        half_w = half_h * (aspect or ASPECT)
        w = np.array([0.62 * size * len(it[0]) + 0.3 * size for it in items]) / 2
        cam[:, 0] = np.clip(cam[:, 0], -half_w + w + 0.3 * size, half_w - w - 0.3 * size)
        cam[:, 1] = np.clip(cam[:, 1], -half_h + size, half_h - size)
        P = cam @ R + origin
    for i, ((text, _, col), pos) in enumerate(zip(items, P)):
        name = f"L{i}"
        cmd.pseudoatom(obj, name=name, pos=[float(x) for x in pos], label=text)
        cname = define(f"mf_lab_{obj}_{i}", col)
        cmd.set("label_color", cname, f"{obj} and name {name}")
    cmd.hide("everything", obj)
    cmd.show("labels", obj)
    cmd.set("label_size", -size * LABEL_SCALE, obj)
    cmd.set("label_font_id", font, obj)
    if outline:
        cmd.set("label_outline_color", outline, obj)


# PyMOL fonts (label_font_id) -> family, weight and style in the SVG
SVG_FONTS = {5: ("sans", "normal", "normal"), 6: ("sans", "normal", "italic"), 7: ("sans", "bold", "normal"),
             8: ("sans", "bold", "italic"), 9: ("serif", "normal", "normal"), 10: ("serif", "bold", "normal"),
             11: ("mono", "normal", "normal"), 13: ("mono", "bold", "normal")}
SVG_FAMILIES = {"sans": "Arial, Helvetica, sans-serif", "serif": "'Times New Roman', Times, serif",
                "mono": "'Courier New', Courier, monospace"}


def labels_svg(path, objs, width, height, dpi, halo_rgb=None, halo_px=0):
    """Labels of objs as editable text in an SVG with the same physical size as the PNG
    (width x height px at dpi), to overlay on <name>_nolabels.png in Illustrator/Inkscape.
    Projects each pseudoatom with the current view (perspective, vertical field_of_view). halo_rgb:
    background-colored outline halo_px wide (the same halo as the PNG)."""
    from xml.sax.saxutils import escape
    view = cmd.get_view()
    R = np.array(view[:9]).reshape(3, 3)
    origin, shift = np.array(view[12:15]), np.array(view[9:12])
    tan_h = math.tan(math.radians(cmd.get_setting_float("field_of_view") / 2.0))
    aspect = width / height
    out = []
    for obj in objs:
        atoms = []
        cmd.iterate_state(1, f"{obj} and not label ''", "atoms.append((label, x, y, z, s.label_color))",
                          space={"atoms": atoms})
        size = abs(cmd.get_setting_float("label_size", obj))
        font = SVG_FONTS.get(cmd.get_setting_int("label_font_id", obj), SVG_FONTS[5])
        for text, x, y, z, col in atoms:
            cam = (np.array([x, y, z]) - origin) @ R + shift
            depth = -cam[2]
            px = (cam[0] / (depth * tan_h * aspect) + 1) / 2 * width
            py = (1 - cam[1] / (depth * tan_h)) / 2 * height
            fs = size / (2 * depth * tan_h) * height * SVG_EM
            rgb = cmd.get_color_tuple(col) if isinstance(col, int) and col >= 0 else \
                cmd.get_color_tuple(cmd.get_color_index(cmd.get_setting_text("label_color", obj)) or 0)
            out.append((text, px, py, fs, rgb, font))
    pt = 72.0 / dpi  # px -> pt: same physical size as the PNG at dpi
    hexc = lambda c: "#%02x%02x%02x" % tuple(int(round(255 * v)) for v in c)
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width * pt:.2f}pt" height="{height * pt:.2f}pt" '
             f'viewBox="0 0 {width} {height}">']
    for text, px, py, fs, rgb, (fam, weight, style) in out:
        halo = (f' stroke="{hexc(halo_rgb)}" stroke-width="{2 * halo_px}" stroke-linejoin="round" '
                f'paint-order="stroke"') if halo_rgb is not None else ""
        lines.append(f'<text x="{px:.1f}" y="{py + 0.36 * fs:.1f}" font-family="{SVG_FAMILIES[fam]}" '
                     f'font-size="{fs:.1f}" font-weight="{weight}" font-style="{style}" text-anchor="middle" '
                     f'fill="{hexc(rgb)}"{halo}>{escape(text)}</text>')
    lines.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return out


def residue_label_items(res_list, obj, center, color_of, offset=2.2, fmt=None):
    """Label position: side-chain tip, pushed away from center."""
    items = []
    for chain, resi, resn in res_list:
        rs = f"{obj} and {residue_sel(chain, resi)}"
        xyz = cmd.get_coords(f"{rs} and not name N+C+O")
        if xyz is None:
            xyz = cmd.get_coords(rs)
        if xyz is None:
            continue
        ca = cmd.get_coords(f"{rs} and name CA")
        base = ca[0] if ca is not None else xyz.mean(0)
        tip = xyz[np.argmax(np.linalg.norm(xyz - base, axis=1))]
        d = tip - center
        d = d / (np.linalg.norm(d) or 1.0)
        text = (fmt or LABEL_FORMATS.get(LABEL_OPTS["style"], "{one}{resi}")).format(one=ONE_LETTER.get(resn, resn), resn=resn.capitalize(), resi=resi, chain=chain)
        items.append((text, tip + offset * d, color_of(chain)))
    return items


# ---------------------------------------------------------------- contacts
def polar_contacts(name, s1, s2, cutoff=3.6, mode="geometric"):
    """Partner-receptor H-bonds: N/O donor/acceptor pairs within cutoff. mode 'geometric' uses
    PyMOL's distance + angle criterion (find_pairs mode=1, works without H); 'distance' keeps
    every donor-acceptor pair within cutoff (catches e.g. Lys NZ or Ser OG, whose estimated H
    position can fail the angle test). Each pair is drawn as a dash in the object name, so the
    .png and the returned list always match.

    Returns [(at1, at2, dist)] with at = (chain, resn, resi, name).
    """
    cmd.set("h_bond_cutoff_center", cutoff)
    cmd.set("h_bond_cutoff_edge", cutoff)
    out = []
    if mode == "distance":
        pairs = cmd.find_pairs(f"({s1}) and donor", f"({s2}) and acceptor", mode=0, cutoff=cutoff)
        pairs += [p for p in cmd.find_pairs(f"({s1}) and acceptor", f"({s2}) and donor", mode=0, cutoff=cutoff)
                  if p not in pairs]
    else:
        pairs = cmd.find_pairs(f"({s1}) and (donor or acceptor)", f"({s2}) and (donor or acceptor)",
                               mode=1, cutoff=cutoff)
    for (m1, i1), (m2, i2) in pairs:
        a1, a2 = f"{m1} and index {i1}", f"{m2} and index {i2}"
        at = []
        for s in (a1, a2):
            cmd.iterate(s, "at.append((chain, resn, resi, name))", space={"at": at})
        out.append((at[0], at[1], cmd.get_distance(a1, a2)))
        cmd.distance(name, a1, a2)
    if out:
        cmd.hide("labels", name)
    return sorted(out, key=lambda t: (t[0][2], t[2]))


def write_contacts(path, header, pairs):
    with open(path, "w") as fh:
        fh.write(header)
        fh.write(f"# n = {len(pairs)}\n# ligand               receptor               dist(A)\n")
        for (c1, n1, r1, a1), (c2, n2, r2, a2), d in pairs:
            fh.write(f"{c1}:{n1}{r1:<6s} {a1:<5s}   {c2}:{n2}{r2:<6s} {a2:<5s}   {d:.2f}\n")
