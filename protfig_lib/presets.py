"""
presets.py — figure styles. A preset = a function that builds the scene and returns the view.

Each preset is registered with @preset(...), which declares:
  partner  'required' | 'optional' | 'none'  whether it needs a binding partner (ligand,
           peptide or protein) and whether PROTFIG should detect it automatically
  kinds    partner types it makes sense for: 'small', 'peptide', 'protein'
  size     default image size (px at 300 dpi)

All colors and details go through ctx.opt(name, default): if the user gave the option
(CLI or .toml) theirs is used, otherwise the preset's. Options are documented in protfig.py.
"""
import os
import sys

import numpy as np
from pymol import cmd

from . import fields
from .core import (BASE, ONE_LETTER, WATERS, readable, _load_structure, GLYCANS, LIPIDS, METALS, ProtFigError, PALETTES, add_labels, adjust, color_by_element, define, facing_R, frame,
                   front_R, membrane_R, mix, overview_R, polar_contacts, unoccluded_R, residue_label_items, residue_sel, residues, rgb_of, site_R,
                   vis_height, write_contacts)

HB_YELLOW = (1.0, 0.88, 0.0)  # H-bonds: yellow, readable on white and on most surfaces
HB_GOLD = (0.914, 0.635, 0.231)  # H-bonds in site (lavender palette)
# close-ups: no smoothed loops or flattened sheets, so the cartoon passes through the CA atoms the
# side-chain sticks start from (otherwise the sticks float next to the ribbon)
CLOSEUP_CARTOON = {"cartoon_smooth_loops": 0, "cartoon_flat_sheets": 0}
LONG_PEPTIDE = 15  # site: peptides longer than this show only their H-bonded residues in sticks, like proteins
PRESETS = {}
# atoms drawn as sticks for a residue: PyMOL's 'sidechain' skips non-standard residues (MSE, SEP,
# staples...), so --nonstandard_sticks on switches to "everything but the backbone N, C, O"
SIDECHAIN = "(sidechain or name CA)"


def set_sidechain_mode(nonstandard):
    global SIDECHAIN
    SIDECHAIN = "(not name N+C+O+OXT)" if nonstandard else "(sidechain or name CA)"


def preset(partner, kinds, size, glycans=False):
    """glycans=True: the style accepts --glycans on (cartoon and surface styles, not close-ups)."""
    def deco(fn):
        PRESETS[fn.__name__.rstrip("_")] = {"fn": fn, "partner": partner, "kinds": kinds, "size": size,
                                            "glycans": glycans, "doc": " ".join((fn.__doc__ or "").split())}
        return fn
    return deco


class Ctx:
    def __init__(self, a, lig, kind, chains, user, session_view):
        self.a = a
        self.lig = lig  # 'ligand' or None
        self.kind = kind  # 'small' | 'peptide' | 'protein' | None
        self.chains = chains
        self.user = user
        self.session_view = session_view
        self.aspect = a.width / a.height
        self.colors = {}
        self.receptor = f"{BASE} and polymer" + (" and not ligand" if lig else "")

    def opt(self, name, default):
        v = getattr(self.a, name, None)
        return default if v is None else v

    def rgb(self, name, default):
        return rgb_of(self.opt(name, default))

    def labels_on(self, default=True):
        return {"on": True, "off": False}.get(self.a.labels, default)

    # -- per-chain colors
    def chain_colors(self, palette):
        pal = PALETTES[self.a.palette or palette]
        over, pending = {}, []
        for tok in (self.a.colors or "").replace(",", ";").split(";"):
            if ":" in tok:  # "A,B:palecyan" gives A and B the same color
                k, v = tok.split(":", 1)
                for c in pending + [k.strip()]:
                    over[c] = rgb_of(v)
                pending = []
            elif tok.strip():
                pending.append(tok.strip())
        for i, c in enumerate(self.chains):
            self.colors[c] = over.get(c, pal[i % len(pal)])
            define(f"mf_chain_{c}", self.colors[c])
            define(f"mf_light_{c}", mix(self.colors[c], "white", 0.45))
            define(f"mf_dark_{c}", mix(self.colors[c], "black", 0.25))
        return self.colors

    def partner_rgb(self, default_small):
        """Partner color: --partner_color > its chain color (peptide/protein) > default."""
        if self.a.partner_color:
            return rgb_of(self.a.partner_color)
        if self.kind != "small" and self.colors:
            ch = cmd.get_chains("ligand")
            if ch and ch[0] in self.colors:
                return self.colors[ch[0]]
        return rgb_of(default_small)

    def chain_names(self):
        out = {}
        for tok in (self.a.chain_names or "").replace(";", ",").split(","):
            if ":" in tok:
                k, v = tok.split(":", 1)
                out[k.strip()] = v.strip()
        return out

    # -- selections with priority CLI > .pse > auto
    def sel(self, key, default=None):
        val = getattr(self.a, key, None)
        if val:
            return f"{BASE} and ({val})"
        if key in self.user:
            return self.user[key]
        return default

    def site(self, cutoff):
        s = self.sel("site")
        if s:
            return f"byres ({s}) and polymer"
        if not self.lig:
            return None
        return f"byres ({self.receptor} within {cutoff} of lig_core)"

    # -- views
    def view(self, kind, center_sel=None, frame_sel=None, margin=1.06, clip=None, fit_box=False, unocclude=None,
             front=None, upright=False):
        """unocclude: (group selection, environment selection) -> turn the automatic view so that the
        group is not hidden behind the environment (covalent preset). front: (partner, receptor) ->
        partner in front of the receptor, running across the image (outline --style comic)."""
        """kind: 'overview' | 'site'. If --view or a .pse view is given, that one is used."""
        a = self.a
        if a.view:
            vals = [float(t) for t in open(a.view).read().replace(",", " ").replace("(", " ").replace(")", " ").split()]
            return tuple(vals[:18])
        if self.session_view and not a.auto_view:
            return self.session_view
        if kind == "site" and self.lig:
            env = f"{self.receptor} within 12 of lig_core"
            c, R = site_R("lig_core", env, self.aspect, up_sel=self.sel("focus"), side=self.kind == "protein")
            if unocclude:
                R = unoccluded_R(R, cmd.get_coords(unocclude[0]) - c, cmd.get_coords(unocclude[1]) - c,
                                 prefer=unocclude[2] if len(unocclude) > 2 else "face")
            R = adjust(R, a)
            front = None
            if clip is not None:
                front = float(np.max((cmd.get_coords("lig_core") - c) @ R[2])) + clip
            return frame(R, c, cmd.get_coords(frame_sel or "lig_core"), self.aspect, margin=margin, front_z=front,
                         fit_box=fit_box)
        sel = center_sel or f"{BASE} and polymer"
        up = self.sel("focus") or (f"{BASE} and chain {self.chains[0]}" if len(self.chains) > 1 else None)
        if front:
            c, R = front_R(*front)
        else:
            # upright: long axis vertical even in a landscape image (e.g. a beta barrel standing up)
            c, R = overview_R(sel, 1.0 if upright else self.aspect, toward="ligand" if self.lig else None, up_sel=up)
        if unocclude:
            R = unoccluded_R(R, cmd.get_coords(unocclude[0]) - c, cmd.get_coords(unocclude[1]) - c)
        R = adjust(R, a)
        return frame(R, c, cmd.get_coords(frame_sel or sel), self.aspect, margin=margin, fit_box=fit_box)


def R_of(view):
    return np.array(view[:9]).reshape(3, 3).T


def label_size(view, frac):
    return frac * vis_height(view)


def finish_lighting(**kw):
    for k, v in kw.items():
        cmd.set(k, v)


# ======================================================================= shared building blocks
def show_partner(ctx, carbon, style="sticks", radius=0.25, cartoon=None):
    """PARTNER_MOL object: small molecule in sticks (or ball-and-stick); peptide/protein in cartoon
    with interface side chains in sticks. carbon = carbon color; cartoon = cartoon color
    (default: lightened carbon)."""
    if not ctx.lig:
        return None
    cmd.create("PARTNER_MOL", "ligand")
    cmd.hide("everything", "PARTNER_MOL")
    # the side-chain helper hides atoms named like nucleic-acid backbone (P, O1P...): phosphotyrosine
    # and other phosphorylated residues would lose their phosphate oxygens
    cmd.set("cartoon_side_chain_helper", 0, "PARTNER_MOL")
    c_name = define("mf_partner_c", carbon)
    if ctx.kind == "small":
        cmd.show("sticks", "PARTNER_MOL")
        if style == "ballstick":
            cmd.show("spheres", "PARTNER_MOL")
            cmd.set("stick_radius", radius * 0.45, "PARTNER_MOL")
            cmd.set("sphere_scale", radius * 0.62, "PARTNER_MOL")
        else:
            cmd.set("stick_radius", radius, "PARTNER_MOL")
        color_by_element("PARTNER_MOL", c_name)
    else:
        cmd.show("cartoon", "PARTNER_MOL")
        if ctx.opt("partner_sticks", "on") == "on":  # --partner_sticks off: cartoon only
            cmd.show("sticks", f"byres (PARTNER_MOL in lig_core) and {SIDECHAIN}")
        cmd.set("stick_radius", radius, "PARTNER_MOL")
        color_by_element("PARTNER_MOL", c_name)
        cmd.set("cartoon_color", define("mf_partner_cartoon", cartoon or mix(carbon, "white", 0.3)), "PARTNER_MOL")
    return "PARTNER_MOL"


def partner_core_obj():
    return "(PARTNER_MOL in lig_core)"


def show_metals(near="ligand", obj="METAL_IONS", scale=0.35):
    """Metal ions as gray spheres; near=None: all of them, scale: sphere size (larger for whole-molecule views)."""
    s = f"{BASE} and resn {'+'.join(METALS)}" + (f" and (all within 6 of {near})" if near else "")
    if cmd.count_atoms(s):
        cmd.create(obj, s)
        cmd.show("spheres", obj)
        cmd.set("sphere_scale", scale, obj)
        cmd.color("grey50", obj)


def hbonds(ctx, partner_sel=None, dash_radius=0.08, name="HBONDS", default_color=HB_YELLOW):
    """Partner <-> receptor H-bonds (yellow by default; --hbond_color / --hbonds off)."""
    if not ctx.lig or ctx.a.hbonds == "off":
        return []
    pairs = polar_contacts(name, "PARTNER_MOL", partner_sel or ctx.receptor, cutoff=ctx.a.hbond_cutoff,
                           mode=ctx.a.hbond_mode)
    if pairs:
        cmd.set("dash_color", define("mf_hbond", ctx.rgb("hbond_color", default_color)), name)
        cmd.set("dash_radius", dash_radius, name)
        write_contacts(f"{ctx.a.prefix}_hbonds.txt",
                       f"# partner-receptor H-bonds (N/O, <= {ctx.a.hbond_cutoff} A, {ctx.a.hbond_mode} criterion)\n",
                       pairs)
    return pairs


def distance_labels(ctx, pairs, view, size, obj="DIST_LABELS", add=True):
    """'3.0 Å' labels at the midpoint of each H-bond (add=False: return the items)."""
    items = []
    R = R_of(view)
    col = ctx.rgb("label_color", (0, 0, 0)) if ctx.a.label_color else (0, 0, 0)
    for (c1, n1, r1, a1), (c2, n2, r2, a2), d in pairs:
        p1 = cmd.get_coords(f"PARTNER_MOL and chain {c1} and resi \\{r1} and name {a1}")
        p2 = cmd.get_coords(f"{BASE} and chain {c2} and resi \\{r2} and name {a2}")
        if p1 is None or p2 is None:
            continue
        items.append((f"{d:.1f} Å", (p1[0] + p2[0]) / 2 + R[0] * size * 1.6, col))
    if not add:
        return items
    add_labels(obj, items, R, size, font=7, view=view)


def write_site(ctx, site_obj, cutoff):
    res = residues(f"({site_obj}) and name CA")
    with open(f"{ctx.a.prefix}_site.txt", "w") as fh:
        fh.write(f"# site = receptor residues within {cutoff} A of the partner (heavy atoms); n = {len(res)}\n")
        fh.write(", ".join(f"{c}:{n}{r}" for c, r, n in res) + "\n")
        if ctx.kind in ("peptide", "protein"):
            pres = residues(f"{partner_core_obj()} and name CA")
            fh.write(f"\n# partner interface (residues within {cutoff} A of the receptor); n = {len(pres)}\n")
            fh.write(", ".join(f"{c}:{n}{r}" for c, r, n in pres) + "\n")
    return res


def site_labels(ctx, view, site_obj, color_of, frac=0.03, font=7, partner_color=None, extra=(), pairs=(),
                only_hbonded=False, add=True):
    """Labels for the site residues and, if the partner is a peptide/protein, for its side too.

    Ligand: all site residues are labeled. Peptide/protein (large interfaces): only residues
    forming H-bonds, on both sides (if there are none, the site residues). --label_sel picks
    the receptor residues by hand. only_hbonded: H-bonded residues only, for a ligand too.
    add=False: return (items, avoid) instead of drawing them.
    """
    if not ctx.labels_on(True):
        return ([], None) if not add else None
    hb_rec = {(c, r) for _, (c, _, r, _), _ in pairs}
    hb_par = {(c, r) for (c, _, r, _), _, _ in pairs}
    res = residues(f"({site_obj}) and name CA")
    if ctx.a.label_sel:
        res = residues(f"({site_obj}) and name CA and ({ctx.a.label_sel})")
    elif (only_hbonded or ctx.kind in ("peptide", "protein")) and hb_rec:
        res = [x for x in res if (x[0], x[1]) in hb_rec]
    res += [x for x in extra if (x[0], x[1]) not in {(y[0], y[1]) for y in res}]
    cen = cmd.get_coords("lig_core").mean(0)
    items = residue_label_items(res, site_obj, cen, color_of)
    avoid = cmd.get_coords("PARTNER_MOL")
    if ctx.kind in ("peptide", "protein"):
        pres = residues(f"{partner_core_obj()} and name CA")
        if hb_par:
            pres = [x for x in pres if (x[0], x[1]) in hb_par]
        rec_c = cmd.get_coords(site_obj).mean(0)
        items += residue_label_items(pres, "PARTNER_MOL", rec_c, lambda c: partner_color or (0.1, 0.4, 0.5))
        avoid = None
    if not add:
        return items, avoid
    add_labels("LABELS", items, R_of(view), label_size(view, frac), font=font, view=view, avoid=avoid)


def site_sticks(ctx, cutoff, obj="SITE_STICKS", radius=0.22):
    cmd.create(obj, ctx.site(cutoff))
    cmd.hide("everything", obj)
    cmd.show("sticks", f"{obj} and {SIDECHAIN}")
    cmd.show("sticks", f"{obj} and resn GLY+PRO and name CA+N+C")
    cmd.set("stick_radius", radius, obj)
    return obj


def load_ref(a, name):
    """--ref: a file, or a PDB code (downloaded next to the outputs, like the main input)."""
    _load_structure(a.ref, name, os.path.dirname(os.path.abspath(a.prefix)))


def show_mutation(ctx):
    """--mutation (any preset): the mutated residue as sticks inside translucent spheres. With
    --mutation_original on and --ref (the wild type, superposed around the mutation), the original
    residue too, as translucent gray sticks."""
    a = ctx.a
    mut = ctx.sel("mutation")
    if not mut or "MUT_STICKS" in cmd.get_names("objects"):
        return
    cmd.create("MUT_STICKS", f"byres ({mut}) and polymer")
    if not cmd.count_atoms("MUT_STICKS"):
        print(f"NOTE: --mutation {a.mutation} selects no residue")
        return
    cmd.hide("everything", "MUT_STICKS")
    cmd.show("sticks", f"MUT_STICKS and not name N+C+O+OXT")
    cmd.show("spheres", "MUT_STICKS and not name N+C+O+OXT+CA")
    color_by_element("MUT_STICKS", define("mf_mut_c", (0.45, 0.72, 0.35)))
    cmd.set("sphere_color", define("mf_mut_sph", ctx.rgb("mutation_color", (0.3, 0.6, 0.95))), "MUT_STICKS")
    cmd.set("sphere_transparency", 0.45, "MUT_STICKS")
    cmd.set("stick_radius", 0.2, "MUT_STICKS")
    if "SITE_STICKS" in cmd.get_names("objects"):
        cmd.remove("SITE_STICKS and (SITE_STICKS in MUT_STICKS)")
    if ctx.opt("mutation_original", "off") != "on":
        return
    if not a.ref:
        print("NOTE: --mutation_original needs --ref (the wild-type structure)")
        return
    load_ref(a, "REF_MUT")
    cmd.remove("REF_MUT and (hydro or solvent)")
    cmd.remove("REF_MUT and not alt ''+A")
    cmd.super("REF_MUT and polymer and name CA", f"{BASE} and polymer and name CA and (all within 15 of MUT_STICKS)",
              cycles=5)
    orig = " or ".join(f"(REF_MUT and {residue_sel(c, r)} and polymer)" for c, r, _ in residues("MUT_STICKS and name CA"))
    if cmd.count_atoms(orig):
        cmd.create("MUT_ORIGINAL", orig)
        cmd.hide("everything", "MUT_ORIGINAL")
        cmd.show("sticks", "MUT_ORIGINAL and not name N+C+O+OXT")
        color_by_element("MUT_ORIGINAL", define("mf_mut_orig", (0.35, 0.35, 0.35)))
        cmd.set("stick_radius", 0.17, "MUT_ORIGINAL")
        cmd.set("sphere_transparency", 0.62, "MUT_STICKS")  # see the original residue through the spheres
    cmd.delete("REF_MUT")


MUT_FONT = 10  # serif bold: in sans fonts the I of isoleucine (M49I) reads as a 1


def mutation_labels(ctx, view):
    """Label of each --mutation residue (wild type, number, mutant: e.g. M49I when the original is
    known with --mutation_original, else I49), unless the preset already labeled it."""
    if getattr(ctx, "mutation_labeled", False):
        return
    items = mutation_label_items(ctx)
    if items:
        add_labels("MUT_LABELS", items, R_of(view), label_size(view, 0.03), font=MUT_FONT, view=view)


def mutation_label_items(ctx):
    """[(text, position, color)] for the --mutation residues: at the side-chain tip, in a dark tone of
    the sphere color."""
    if "MUT_STICKS" not in cmd.get_names("objects") or not ctx.labels_on(True):
        return []
    orig = {}
    if "MUT_ORIGINAL" in cmd.get_names("objects"):
        cmd.iterate("MUT_ORIGINAL and name CA", "orig[(chain, resi)] = resn", space={"orig": orig})
    col = mix(ctx.rgb("mutation_color", (0.3, 0.6, 0.95)), "black", 0.45)
    center = cmd.get_coords(f"{BASE} and ligand") if ctx.lig else None
    center = center.mean(0) if center is not None else cmd.get_coords("MUT_STICKS").mean(0)
    res = residues("MUT_STICKS and name CA")
    items = []
    for c, r, n in res:
        fmt = None
        if (c, r) in orig:
            fmt = f"{ONE_LETTER.get(orig[(c, r)], orig[(c, r)])}{{resi}}{{one}}"
        items += residue_label_items([(c, r, n)], "MUT_STICKS", center, lambda ch: col, fmt=fmt)
    return items


def show_ref_ligand(ctx):
    """--ref: the ligand of a second structure (same residue name, superposed on the site) as a
    translucent gray ghost, to compare binding modes (e.g. wild type vs mutant)."""
    load_ref(ctx.a, "REF_FULL")
    cmd.remove("REF_FULL and hydro")
    cmd.remove("REF_FULL and not alt ''+A")
    # superpose on the CAs around the site (not the whole protein)
    cmd.super("REF_FULL and polymer and name CA", f"{ctx.receptor} and name CA and (all within 15 of ligand)", cycles=5)
    lig_resn = set()
    cmd.iterate("ligand", "lig_resn.add(resn)", space={"lig_resn": lig_resn})
    ref_lig = f"REF_FULL and resn {'+'.join(lig_resn)} and (all within 4 of ligand)"
    if cmd.count_atoms(ref_lig):
        cmd.create("REF_PARTNER", f"byres ({ref_lig})")
        cmd.hide("everything", "REF_PARTNER")
        cmd.show("sticks", "REF_PARTNER")
        cmd.set("stick_radius", 0.2, "REF_PARTNER")
        cmd.color("grey50", "REF_PARTNER")
        cmd.set("stick_transparency", 0.55, "REF_PARTNER")
    cmd.delete("REF_FULL")


def trim_copies(ctx):
    """Close-ups without --chains: crystal copies away from the partner are left out, and a
    peptide/protein partner takes the second color of the palette (contrast with the receptor,
    whatever its chain ID)."""
    a = ctx.a
    if a.chains:
        return
    touching = cmd.get_chains(f"({ctx.receptor}) within {a.site_cutoff} of ligand")
    if touching:
        cmd.remove(f"({ctx.receptor}) and not chain {'+'.join(touching)}")
    present = cmd.get_chains(f"{BASE} and polymer")
    par = [c for c in ctx.chains if c in cmd.get_chains("ligand and polymer")]
    rec = [c for c in ctx.chains if c in present and c not in par]
    ctx.chains = rec[:1] + par + rec[1:] if par else rec


def refine_sticks(ctx, pairs, carbon_of, radius=0.22):
    """Sticks of a close-up after the H-bonds are known (site, all styles): H-bonded residues
    always shown; modified residues whole; large interfaces (protein, peptide > LONG_PEPTIDE) only the
    H-bonded residues; partner backbone hidden; backbone H-bonds show their C=O / N-CA piece.
    carbon_of(chain) -> carbon color name of the receptor sticks."""
    for _, (c2, _, r2, _), _ in pairs:  # H-bonded residues are shown even if outside the cutoff
        cmd.create("SITE_STICKS", f"SITE_STICKS or ({BASE} and chain {c2} and resi \\{r2})")
    # (re)style SITE_STICKS: re-creating it above adds residues without their sticks or colors
    cmd.hide("everything", "SITE_STICKS")
    cmd.show("sticks", f"SITE_STICKS and ({SIDECHAIN} or (resn GLY+PRO and name CA+N+C))")
    cmd.set("stick_radius", radius, "SITE_STICKS")
    for c in cmd.get_chains("SITE_STICKS"):
        color_by_element(f"SITE_STICKS and chain {c}", carbon_of(c))
    # modified residues (pTyr, pSer, MSE...): PyMOL's "sidechain" leaves out their phosphate etc.
    std = "ALA+ARG+ASN+ASP+CYS+GLN+GLU+GLY+HIS+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL"
    for obj in ("SITE_STICKS", "PARTNER_MOL"):
        if obj in cmd.get_names("objects"):
            cmd.show("sticks", f"{obj} and polymer and not resn {std} and not name N+C+O+OXT"
                     + (" and (PARTNER_MOL in lig_core)" if obj == "PARTNER_MOL" else ""))
    n_par = cmd.count_atoms("ligand and polymer and name CA")
    if ctx.kind == "protein" or (ctx.kind == "peptide" and n_par > LONG_PEPTIDE):
        # large interfaces (proteins, long peptides): only the H-bonded residues in sticks (the labeled ones)
        hb_rec = " or ".join(f"(chain {c} and resi \\{r})" for _, (c, _, r, _), _ in pairs) or "none"
        hb_par = " or ".join(f"(chain {c} and resi \\{r})" for (c, _, r, _), _, _ in pairs) or "none"
        cmd.hide("sticks", f"SITE_STICKS and not byres ({hb_rec})")
        cmd.hide("sticks", f"PARTNER_MOL and not byres ({hb_par})")
    if "PARTNER_MOL" in cmd.get_names("objects") and ctx.kind != "small":
        cmd.hide("sticks", "PARTNER_MOL and name N+C+O+OXT and not (resn GLY+PRO and name N+C)")
    # H-bonds to backbone atoms: show that piece of backbone (C=O or N-CA), or the dash ends in the air
    for obj, (c, _, r, at) in [(o, x) for p_, r_, _ in pairs for o, x in (("PARTNER_MOL", p_), ("SITE_STICKS", r_))]:
        piece = {"O": "C+O+CA", "OXT": "C+OXT+CA", "N": "N+CA"}.get(at)
        if piece:
            cmd.show("sticks", f"{obj} and chain {c} and resi \\{r} and name {piece}")


def palette_pair(ctx, receptor_default, partner_default):
    """(receptor, partner) colors of the surface presets: their own defaults, or with --palette the
    light tones of its first and second colors (lavender, nordic, sage, candy, vivid...)."""
    if not ctx.a.palette:
        return rgb_of(receptor_default), rgb_of(partner_default)
    pal = PALETTES[ctx.a.palette]
    return mix(pal[0], "white", 0.45), mix(pal[1 % len(pal)], "white", 0.45)


def need_partner(ctx, name):
    if not ctx.lig:
        raise ProtFigError(f"preset '{name}' needs a binding partner: no ligand or second chain was found. "
                 f"Specify it with --partner (e.g. --partner \"resn ATP\" or --partner \"chain B\").")


def surface_behind(ctx, view, sel, ahead, name="SURFACE", radius=22):
    """Object with the atoms of sel lying behind the plane z = partner center + ahead."""
    R, origin = R_of(view), np.array(view[12:15])
    zc = float(((cmd.get_coords("lig_core") - origin) @ R[2]).mean())
    cmd.create("_mf_depth", f"({sel}) and (all within {radius} of lig_core)")
    cmd.alter_state(1, "_mf_depth", "b = (x - ox) * r0 + (y - oy) * r1 + (z - oz) * r2",
                    space={"ox": origin[0], "oy": origin[1], "oz": origin[2],
                           "r0": R[2][0], "r1": R[2][1], "r2": R[2][2]})
    cmd.create(name, f"_mf_depth and b < {zc + ahead}")
    cmd.delete("_mf_depth")


# ======================================================================= presets
def show_glycans(ctx):
    """--glycans on: glycans (NAG, MAN, FUC...) as thick sticks, C in --glycan_color. Leaves the
    partner alone (if the glycan is the partner, the style already draws it)."""
    sel = f"{BASE} and resn {'+'.join(GLYCANS)}" + (" and not ligand" if ctx.lig else "")
    cmd.create("GLYCANS", sel)
    if cmd.count_atoms("GLYCANS") == 0:
        print("WARNING: --glycans on, but the structure has no glycans (NAG, MAN, FUC...)")
        return
    cmd.hide("everything", "GLYCANS")
    cmd.show("sticks", "GLYCANS")
    cmd.set("stick_radius", 0.7, "GLYCANS")  # visible at whole-molecule scale
    color_by_element("GLYCANS", define("mf_glycan", ctx.rgb("glycan_color", (0.25, 0.55, 0.85))))


def show_het(ctx, scale=0.35):
    """--keep_het: the remaining non-polymer residues (additives, buffers, ions, unassigned
    ligands) as gray sticks, single atoms as spheres. Skips the partner, waters, glycans and
    metals already drawn by --metals / site."""
    sel = (f"{BASE} and not polymer and not solvent and not resn {'+'.join(GLYCANS)}"
           + (" and not ligand" if ctx.lig else "")
           + (f" and not resn {'+'.join(METALS)}" if "METAL_IONS" in cmd.get_names("objects") else ""))
    cmd.create("HET_RES", sel)
    if cmd.count_atoms("HET_RES") == 0:
        print("WARNING: --keep_het, but the structure has no additives or ions left")
        return
    cmd.hide("everything", "HET_RES")
    cmd.show("sticks", "HET_RES")
    cmd.set("stick_radius", 0.2 if scale < 0.5 else 0.7, "HET_RES")
    color_by_element("HET_RES", define("mf_het", ctx.rgb("het_color", (0.6, 0.6, 0.6))))
    lone = cmd.get_model("HET_RES")
    by_res = {}
    for at in lone.atom:
        by_res.setdefault((at.chain, at.resi, at.resn), []).append(at)
    singles = [f"(chain {c} and resi \\{r} and resn {n})" for (c, r, n), ats in by_res.items() if len(ats) == 1]
    if singles:
        cmd.show("spheres", "HET_RES and (" + " or ".join(singles) + ")")
        cmd.set("sphere_scale", scale, "HET_RES")


def _overview_plain(ctx):
    """overview, default style: cartoon per chain + translucent surface of the same color, ligand in sticks."""
    a = ctx.a
    ctx.chain_colors("vivid")
    cmd.create("CARTOON", f"{BASE} and polymer")
    cmd.create("SURFACE", f"{BASE} and polymer")
    cmd.hide("everything", "CARTOON or SURFACE")
    cmd.show("cartoon", "CARTOON")
    cmd.show("surface", "SURFACE")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.72), "SURFACE")
    cmd.set("cartoon_loop_radius", 0.3)
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"CARTOON and chain {c}")
        cmd.color(f"mf_light_{c}", f"SURFACE and chain {c}")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.15, 0.75, 0.85)), radius=0.45)
        hbonds(ctx, dash_radius=0.1)
    return ctx.view("overview", margin=1.0)


OVERVIEW_STYLES = ("outline", "comic", "illustrative", "pastel", "assembly")


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1500, 2000), glycans=True)
def overview(ctx):
    """Whole molecule. --style: (default) cartoon per chain + translucent surface; outline (outlined
    cartoon, --highlight region, key residues), comic (white cartoon, green partner in front),
    illustrative (flat shading, helices as cylinders; --repr surface), pastel (N->C gradient),
    assembly (surface, one shade per chain, side view). --highlight / --key and --chain_names work
    in every style."""
    a = ctx.a
    style_fn = {"outline": _overview_outline, "comic": _overview_outline, "illustrative": _overview_illustrative,
                "pastel": _overview_pastel, "assembly": _overview_assembly}.get(a.style, _overview_plain)
    view = style_fn(ctx)
    if a.style not in ("outline", "comic"):
        overview_highlight_key(ctx)
    chain_name_labels(ctx, view)
    return view


def overview_highlight_key(ctx):
    """--highlight (a region, in --highlight_color) and --key (residues as sticks) for the overview
    styles that do not draw them themselves."""
    a = ctx.a
    objs = [o for o in ("CARTOON", "SURFACE") if o in cmd.get_names("objects")]
    hl = ctx.sel("highlight") if (a.highlight or "").strip().lower() != "interface" else ctx.site(a.site_cutoff)
    if hl and objs:
        for o in objs:
            cmd.color(define("mf_hl", ctx.rgb("highlight_color", (0.98, 0.52, 0.1))), f"{o} and ({o} in ({hl}))")
    key = ctx.sel("key")
    if key:
        cmd.create("KEY", f"byres ({key}) and polymer")
        cmd.hide("everything", "KEY")
        cmd.show("sticks", f"KEY and {SIDECHAIN}")
        cmd.set("stick_radius", 0.25, "KEY")
        color_by_element("KEY", define("mf_key", ctx.rgb("key_color", (1.0, 1.0, 0.1))))


def chain_name_labels(ctx, view):
    """--chain_names "A:Kinase,B:Cyclin" (and --partner_name for a small ligand): each name beside its
    chain, on the outer side (left for chains left of the center, right for the others), level with
    the chain's center."""
    a = ctx.a
    names = ctx.chain_names()
    if not (names or a.partner_name) or not ctx.labels_on(True):
        return
    R = R_of(view)
    size = label_size(view, 0.05)
    poly = cmd.get_coords(f"{BASE} and polymer")
    zf = (poly @ R.T)[:, 2].max() - 4  # +3 from label_position: inside the clip
    mid_x = np.median((poly @ R.T)[:, 0])
    items = []
    for c in ctx.chains:
        if c not in names or not cmd.count_atoms(f"{BASE} and polymer and chain {c}"):
            continue
        cam = cmd.get_coords(f"{BASE} and polymer and chain {c}") @ R.T
        cy = np.median(cam[:, 1])
        band = cam[np.abs(cam[:, 1] - cy) < 6]
        band = band if len(band) else cam
        half = 0.33 * size * len(names[c]) + 2.0
        x = band[:, 0].min() - half if np.median(cam[:, 0]) <= mid_x else band[:, 0].max() + half
        items.append((names[c], np.array([x, cy, zf]) @ R, readable(mix(ctx.colors.get(c, (0.2, 0.2, 0.2)), "black", 0.1))))
    if ctx.kind == "small" and a.partner_name and "PARTNER_MOL" in cmd.get_names("objects"):
        L = cmd.get_coords("PARTNER_MOL") @ R.T
        p = np.array([L[:, 0].min() - 0.33 * size * len(a.partner_name) - 3, L[:, 1].max() + 4, zf])
        items.append((a.partner_name, p @ R, readable(ctx.partner_rgb((0.1, 0.7, 0.8)))))
    add_labels("CHAIN_LABELS", items, R, size, font=5, view=view)


@preset(partner="required", kinds=("small", "peptide", "protein"), size=(2000, 1700))
def site(ctx):
    """Binding-site close-up: light cartoon per chain, site residues in sticks colored by
    chain, partner in sticks (ligand in yellow with a translucent surface) or in cartoon + interface
    sticks (peptide/protein), H-bonds and labels on both sides. --highlight interface (peptide/protein):
    the cartoon of the interface residues of both sides in a deeper tone. --style faded: the protein as a
    ghost-gray cartoon, only the sticks in the palette colors. --distances on: H-bond lengths in Å.
    --ref (small-molecule partner): the ligand of a 2nd structure as a translucent gray ghost.
    A covalently bound partner (retinal-Lys, drug-Cys...) is drawn bonded to its residue, labeled.
    --style outline: illustration look, opaque pastel cartoon and sticks (light palette colors) with a
    black outline."""
    a = ctx.a
    need_partner(ctx, "site")
    trim_copies(ctx)
    ctx.chain_colors("lavender")
    finish_lighting(**CLOSEUP_CARTOON)
    cut = a.site_cutoff
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.35), "CARTOON")
    site_sticks(ctx, cut)
    faded = a.style == "faded"
    outlined = a.style == "outline"
    for c in ctx.chains:
        cmd.color(f"mf_light_{c}", f"CARTOON and chain {c}")
        color_by_element(f"SITE_STICKS and chain {c}", f"mf_chain_{c}")
    if faded:  # protein as a ghost; the backbone under the site less faded, so the sticks sit on it
        cmd.color(define("mf_ghost", ctx.rgb("receptor_color", "grey70")), "CARTOON")
        cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.72), "CARTOON")
        cmd.set("cartoon_transparency", 0.35, "CARTOON and byres (CARTOON in SITE_STICKS)")
    if outlined:  # opaque pastel cartoon, sticks in the same light tone (illustration look)
        cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.0), "CARTOON")
    pc = ctx.partner_rgb((0.91, 0.76, 0.11))  # small ligand: yellow (peptide/protein: its chain color)
    show_partner(ctx, pc, radius=0.24)
    if ctx.kind != "small":
        cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.0 if outlined else 0.35), "PARTNER_MOL")
        if outlined:
            light_pc = mix(pc, "white", 0.45)
            cmd.set("cartoon_color", define("mf_partner_light", light_pc), "PARTNER_MOL")
            color_by_element("PARTNER_MOL", "mf_partner_light")
    if a.ref and ctx.kind == "small":
        show_ref_ligand(ctx)
    if ctx.opt("metals", "on") == "on":
        show_metals()
    if ctx.kind == "small" and ctx.opt("partner_surface", "on") == "on":
        cmd.create("PARTNER_SURF", "ligand")
        cmd.hide("everything", "PARTNER_SURF")
        cmd.show("surface", "PARTNER_SURF")
        cmd.color(define("mf_ligsurf", ctx.rgb("partner_surface_color", (0.97, 0.90, 0.60))), "PARTNER_SURF")
        cmd.set("transparency", 0.55, "PARTNER_SURF")
    # --highlight interface (peptide/protein partners): the cartoon of the interface residues of both
    # sides in a tone between the light cartoon and the sticks
    if (a.highlight or "").strip().lower() == "interface":
        if ctx.kind == "small":
            print("NOTE: --highlight interface in site applies to peptide/protein partners only")
        else:
            rec_hl = f"CARTOON and byres (CARTOON in ({ctx.site(cut)}))"
            for c in ctx.chains:
                cmd.color(define(f"mf_hl_{c}", mix(ctx.colors[c], "white", 0.28)), f"({rec_hl}) and chain {c}")
            cmd.set("cartoon_transparency", 0.15, rec_hl)
            par_hl = f"PARTNER_MOL and ({partner_core_obj()})"
            cmd.set("cartoon_color", define("mf_par_hl", mix(pc, "white", 0.15)), par_hl)
            cmd.set("cartoon_transparency", 0.15, par_hl)
    pairs = hbonds(ctx, dash_radius=0.08, default_color=HB_GOLD)
    refine_sticks(ctx, pairs, (lambda c: f"mf_light_{c}") if outlined else (lambda c: f"mf_chain_{c}"),
                  radius=0.25 if outlined else 0.22)
    if outlined:
        finish_lighting(ray_trace_mode=1, ray_trace_gain=0.15, ambient=0.5, direct=0.45, specular=0,
                        spec_reflect=0, cartoon_loop_radius=0.3, cartoon_fancy_helices=0, ray_shadows=0, fog_start=0.6)
    # covalent partner (retinal-Lys, drug-Cys, PLP; GFP chromophore built into the chain): the bond
    # drawn between the partner and its residue, which is labeled; chain neighbors of a built-in group hidden
    link_res = []
    if ctx.kind == "small":
        links, side, link_res, _ = partner_links(ctx)
        if side:
            atoms = " or ".join(f"index {x[4]}" for p, r, _, _ in side for x in (p, r))
            cmd.create("COV_BOND", f"{BASE} and ({atoms})")
            cmd.hide("everything", "COV_BOND")
            cmd.show("sticks", "COV_BOND")
            cmd.set("stick_radius", 0.22, "COV_BOND")
            color_by_element("COV_BOND", define("mf_cov_c", pc))
        elif links:
            nb = " or ".join(f"(chain {r[0]} and resi \\{r[1]})" for _, r, _, _ in links)
            cmd.remove(f"SITE_STICKS and byres ({nb})")  # built-in group: its chain neighbors are not "site"
    write_site(ctx, "SITE_STICKS", cut)
    show_mutation(ctx)  # before framing, so the mutated residue is in view
    frame_sel = f"SITE_STICKS or {partner_core_obj()}" + (" or MUT_STICKS" if "MUT_STICKS" in cmd.get_names("objects") else "")
    view = ctx.view("site", frame_sel=frame_sel, margin=1.08, clip=a.clip, fit_box=True,
                    unocclude=("lig_core", "CARTOON and name N+CA+C+O", "turn"))
    # residue labels, H-bond distances and the mutation, separated together (no overlaps)
    items, avoid = site_labels(ctx, view, "SITE_STICKS", lambda c: readable(mix(ctx.colors.get(c, (0, 0, 0)), "black", 0.15)),
                               partner_color=readable(mix(pc, "black", 0.25)), pairs=pairs, add=False,
                               only_hbonded=faded,  # faded: the H-bonds (and their lengths) are the point
                               extra=link_res)
    size = label_size(view, 0.03)
    if pairs and ctx.labels_on(True) and ctx.opt("distances", "off") == "on":
        items += distance_labels(ctx, pairs, view, size, add=False)
    if items:
        add_labels("LABELS", items, R_of(view), size, font=7, view=view, avoid=avoid)
    mut_items = mutation_label_items(ctx)
    if mut_items:
        add_labels("MUT_LABELS", mut_items, R_of(view), size, font=MUT_FONT, view=view)
        ctx.mutation_labeled = True
    return view


def covalent_links(partner_sel, receptor_sel, cutoff):
    """Partner-receptor atom pairs within cutoff Å -> [(partner_atom, receptor_atom, dist, bonded)],
    each atom as (chain, resi, resn, name, index). Bonds already in the file are within any sensible
    cutoff, so the distance search finds them too; bonded says whether the file had the bond."""
    def atom(obj, idx):
        out = []
        cmd.iterate(f"{obj} and index {idx}", "out.append((chain, resi, resn, name, index))", space={"out": out})
        return out[0]
    links = []
    for (o1, i1), (o2, i2) in cmd.find_pairs(partner_sel, receptor_sel, cutoff=cutoff, mode=0):
        p, r = atom(o1, i1), atom(o2, i2)
        d = cmd.get_distance(f"{o1} and index {i1}", f"{o2} and index {i2}")
        bonded = cmd.count_atoms(f"({o1} and index {i1}) and bound_to ({o2} and index {i2})") > 0
        links.append((p, r, d, bonded))
    return sorted(links, key=lambda x: x[2])


def partner_links(ctx, required=False):
    """Covalent links between the partner and the receptor (bond in the file, or atoms within
    --covalent_cutoff Å, default 2.0: the bond is added). Writes <name>_covalent.txt. Returns
    (links, side, rec_res, par_res): side = links to receptor side chains (retinal-Lys, drug-Cys);
    links only to backbone N/C mean a group built into the chain (GFP chromophore)."""
    a = ctx.a
    cut = ctx.opt("covalent_cutoff", 2.0)
    rec = f"({ctx.receptor}) and not resn {'+'.join(METALS)} and not hydro"
    links = covalent_links(f"{BASE} and ligand and not hydro", rec, cut)
    if not links:
        if required:
            raise ProtFigError(f"covalent: no partner atom within {cut} Å of the receptor; "
                               "raise --covalent_cutoff or check --partner")
        return [], [], [], []
    for p, r, d, bonded in links:
        if not bonded:
            cmd.bond(f"{BASE} and index {p[4]}", f"{BASE} and index {r[4]}")
        print(f"covalent link: {p[0]}:{p[2]}{p[1]} {p[3]} - {r[0]}:{r[2]}{r[1]} {r[3]}  {d:.2f} Å"
              + ("" if bonded else "  (bond not in the file: added)"))
    with open(f"{a.prefix}_covalent.txt", "w") as fh:
        fh.write("partner_atom\treceptor_atom\tdistance_A\tbond_in_file\n")
        for p, r, d, bonded in links:
            fh.write(f"{p[0]}:{p[2]}{p[1]}:{p[3]}\t{r[0]}:{r[2]}{r[1]}:{r[3]}\t{d:.2f}\t{'yes' if bonded else 'no'}\n")
    side = [x for x in links if x[1][3] not in ("N", "C")]
    rec_res = sorted({(r[0], r[1], r[2]) for _, r, _, _ in side})
    par_res = sorted({(p[0], p[1], p[2]) for p, _, _, _ in links})
    return links, side, rec_res, par_res


@preset(partner="required", kinds=("small", "peptide"), size=(1800, 1800))
def covalent(ctx):
    """Whole protein with a group covalently bound to it or built into it (GFP chromophore,
    retinal-Lys in opsins, covalent drugs, PLP...): cartoon, the group (and the residue it is bonded to)
    in sticks, long axis vertical, view turned so the group is not hidden. The close-up of the link is
    site (it detects covalent links). Bonds missing from the file are added for partner-receptor atoms
    within --covalent_cutoff (default 2.0 Å)."""
    a = ctx.a
    need_partner(ctx, "covalent")
    ctx.chain_colors("lavender")
    links, side, rec_res, par_res = partner_links(ctx, required=True)
    rec_sel = " or ".join(f"({BASE} and {residue_sel(c, i)})" for c, i, _ in rec_res) or "none"
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.set("cartoon_loop_radius", 0.3)
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"CARTOON and chain {c}")
    if a.cartoon_transparency is not None:
        cmd.set("cartoon_transparency", a.cartoon_transparency, "CARTOON")
    # partner + bonded receptor residue(s) in one object: PyMOL only draws bonds within an object
    pc = ctx.partner_rgb((0.91, 0.76, 0.11))  # yellow, as a ligand in site
    cmd.create("LINKED", f"({BASE} and ligand and not hydro) or ({rec_sel})")
    cmd.hide("everything", "LINKED")
    cmd.show("sticks", "LINKED")
    cmd.set("stick_radius", 0.45, "LINKED")
    color_by_element("LINKED", define("mf_partner_c", pc))
    for c, i, _ in rec_res:
        color_by_element(f"LINKED and {residue_sel(c, i)}", f"mf_chain_{c}")
    if ctx.opt("metals", "off") == "on":
        show_metals()
    finish_lighting(ambient=0.5, direct=0.55, specular=0.2)
    view = ctx.view("overview", margin=1.05, fit_box=True, unocclude=("LINKED", "CARTOON and name N+CA+C+O"),
                    upright=True)
    # labels (off by default): the group and its bonded residue, at their far ends
    if ctx.labels_on(False):
        center = cmd.get_coords("LINKED").mean(0)
        items = []
        for res_list, color_of in ((par_res, lambda c: mix(pc, "black", 0.3)),
                                   (rec_res, lambda c: mix(ctx.colors.get(c, (0, 0, 0)), "black", 0.2))):
            for c, i, n in res_list:
                xyz = cmd.get_coords(f"LINKED and {residue_sel(c, i)}")
                tip = xyz[np.argmax(np.linalg.norm(xyz - center, axis=1))]
                d = (tip - center) / (np.linalg.norm(tip - center) or 1.0)
                text = residue_label_items([(c, i, n)], "LINKED", center, color_of)[0][0]
                items.append((text, tip + 6.0 * d, readable(color_of(c))))
        add_labels("LABELS", items, R_of(view), label_size(view, 0.03), font=7, view=view)
    return view


@preset(partner="required", kinds=("small", "peptide"), size=(1500, 1700))
def pocket(ctx):
    """Surface pocket (per chain) in a "cutaway" view: only the surface behind the partner, which
    lies in its groove; site sticks. For ligands and peptides (as cartoon). --palette lavender etc. to
    change the default vivid colors."""
    a = ctx.a
    need_partner(ctx, "pocket")
    ctx.chain_colors("vivid")
    if ctx.kind == "peptide" and a.partner_sticks is None:
        a.partner_sticks = "off"  # a peptide as cartoon only: the groove is the point (--partner_sticks on: side chains)
    show_partner(ctx, ctx.partner_rgb((0.05, 0.62, 0.72)), radius=0.3)
    view = ctx.view("site", frame_sel="PARTNER_MOL", margin=1.7 if ctx.kind == "small" else 1.25)
    surface_behind(ctx, view, ctx.receptor, a.clip if a.clip is not None else 1.5)
    cmd.hide("everything", "SURFACE")
    cmd.show("surface", "SURFACE")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.1), "SURFACE")
    cmd.create("SITE_STICKS", f"byres (({ctx.site(a.site_cutoff)}) and ({BASE} in SURFACE))")
    cmd.hide("everything", "SITE_STICKS")
    cmd.show("sticks", f"SITE_STICKS and {SIDECHAIN}")
    cmd.set("stick_radius", 0.28, "SITE_STICKS")
    for c in ctx.chains:
        cmd.color(f"mf_light_{c}", f"SURFACE and chain {c}")
        cmd.color(f"mf_chain_{c}", f"SITE_STICKS and chain {c}")
    write_site(ctx, "SITE_STICKS", a.site_cutoff)
    finish_lighting(ray_shadows=0, ambient=0.5, direct=0.5)
    return view


def _overview_outline(ctx):
    """Cartoon with black outline: base color (--receptor_color) + highlighted region (--highlight,
    --highlight_color), partner in cyan sticks (peptide/protein in cartoon), key residues
    (those forming H-bonds, plus --key) in yellow sticks. --highlight interface shades the
    interface on both sides (--partner_highlight_color) and shows key residues of both sides.
    --style comic: white cartoon, peptide/protein partner as a green cartoon, no sticks or H-bonds.
    A peptide/protein partner is turned to face the camera."""
    a = ctx.a
    comic = a.style == "comic"
    if comic:
        a.partner_sticks = a.partner_sticks or "off"
        a.hbonds = a.hbonds or "off"
    want_iface = (a.highlight or "").strip().lower() == "interface"
    iface = want_iface and ctx.kind in ("peptide", "protein")  # both sides; a small ligand: receptor side only
    if comic and ctx.lig and not a.chains:  # crystal copies away from the partner are left out
        touching = cmd.get_chains(f"({ctx.receptor}) within {a.site_cutoff} of ligand")
        if touching:
            cmd.remove(f"({ctx.receptor}) and not chain {'+'.join(touching)}")
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.color(define("mf_base", ctx.rgb("receptor_color", (1, 1, 1) if comic else (0.25, 0.42, 0.75))), "CARTOON")
    hl_c = ctx.rgb("highlight_color", (0.98, 0.52, 0.1))
    hl = ctx.site(a.site_cutoff) if want_iface else ctx.sel("highlight")
    if hl:
        cmd.color(define("mf_hl", hl_c), f"CARTOON in ({hl})")
    finish_lighting(cartoon_fancy_helices=1, cartoon_highlight_color=-1, cartoon_loop_radius=0.3,
                    cartoon_oval_length=1.35)
    par_c = rgb_of(a.partner_color) if comic and a.partner_color else \
        rgb_of("#5CCB86") if comic else ctx.partner_rgb((0.35, 0.85, 0.9))
    if ctx.lig:
        if iface and a.partner_sticks is None:
            a.partner_sticks = "off"  # only the key residues in sticks
        show_partner(ctx, par_c, radius=0.22, cartoon=par_c if iface or comic else None)
        if ctx.kind == "peptide" and not iface and a.partner_sticks != "off":
            cmd.hide("cartoon", "PARTNER_MOL")
            cmd.show("sticks", "PARTNER_MOL")
    par_hl_c = ctx.rgb("partner_highlight_color", mix(par_c, "black", 0.35))
    if iface:
        cmd.color(define("mf_par_hl", par_hl_c), f"PARTNER_MOL and ({partner_core_obj()})")
        cmd.set("cartoon_color", "mf_par_hl", f"PARTNER_MOL and ({partner_core_obj()})")
    pairs = hbonds(ctx, "CARTOON", dash_radius=0.09)
    hb_res = {(c, r) for _, (c, _, r, _), _ in pairs}
    if iface:
        hb_res |= {(c, r) for (c, _, r, _), _, _ in pairs}
    key = " or ".join([f"({ctx.sel('key')})"] * bool(ctx.sel("key"))
                      + [f"({BASE} and chain {c} and resi \\{r})" for c, r in sorted(hb_res)])
    if key:
        cmd.create("KEY", f"byres ({key}) and polymer")
        cmd.hide("everything", "KEY")
        cmd.show("sticks", f"KEY and {SIDECHAIN}")
        cmd.set("stick_radius", 0.22, "KEY")
        if iface and ctx.lig and not a.key_color:  # each side in a darker shade of its interface color
            color_by_element("KEY and not (KEY in PARTNER_MOL)", define("mf_key", mix(hl_c, "black", 0.2)))
            color_by_element("KEY and (KEY in PARTNER_MOL)", define("mf_key_par", mix(par_hl_c, "black", 0.2)))
        else:
            color_by_element("KEY", define("mf_key", ctx.rgb("key_color", (1.0, 1.0, 0.1))))
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.05, ambient=0.45, direct=0.55, specular=0.15,
                    spec_reflect=0.15, fog_start=0.4, fog=0.6)
    pp = ctx.kind in ("peptide", "protein")
    if comic and pp:  # partner in front, across the image, tightly framed
        view = ctx.view("overview", margin=1.02, fit_box=True, front=("PARTNER_MOL", "CARTOON"))
    else:
        view = ctx.view("overview", margin=1.04,
                        unocclude=("PARTNER_MOL", "CARTOON and name N+CA+C+O") if pp else None)
    if key and ctx.labels_on(False) and ctx.lig:
        rec = residues("KEY and not (KEY in PARTNER_MOL) and name CA")
        par = residues("KEY and (KEY in PARTNER_MOL) and name CA")
        items = residue_label_items(rec, "KEY", cmd.get_coords("KEY").mean(0) if not par else
                                    cmd.get_coords("KEY and (KEY in PARTNER_MOL)").mean(0), lambda c: mix(hl_c, "black", 0.45))
        if par:
            items += residue_label_items(par, "KEY", cmd.get_coords("KEY and not (KEY in PARTNER_MOL)").mean(0) if rec else
                                         cmd.get_coords("KEY").mean(0), lambda c: mix(par_hl_c, "black", 0.45))
        add_labels("LABELS", items, R_of(view), label_size(view, 0.03), font=7, view=view)
    return view


def _overview_illustrative(ctx):
    """Illustration style (flat shading + outline): helices as pastel cylinders per chain, or a
    flat surface with --repr surface. Works equally for a protein alone, protein-ligand or
    protein-protein complexes (each chain its own color)."""
    a = ctx.a
    ctx.chain_colors("illustrative")
    obj = "SURFACE" if a.repr == "surface" else "CARTOON"
    cmd.create(obj, ctx.receptor if ctx.kind == "small" else f"{BASE} and polymer")
    cmd.hide("everything", obj)
    if a.repr == "surface":
        cmd.show("surface", obj)
    else:
        cmd.show("cartoon", obj)
        finish_lighting(cartoon_cylindrical_helices=1, cartoon_loop_radius=0.35, cartoon_flat_sheets=1)
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"{obj} and chain {c}")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.92, 0.82, 0.9)), radius=0.22)
    finish_lighting(ray_trace_mode=3, ray_trace_gain=0.0, ambient=0.8, direct=0.3, reflect=0.2, specular=0,
                    depth_cue=0, ray_shadows=0)
    return ctx.view("overview", margin=1.05)


def _overview_pastel(ctx):
    """Pastel cartoon with an N->C gradient within each chain and soft lighting.
    --waters shows crystallographic waters as spheres."""
    a = ctx.a
    pal = PALETTES[a.palette or "pastel"]
    cmd.create("CARTOON", ctx.receptor if ctx.kind == "small" else f"{BASE} and polymer")
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    finish_lighting(cartoon_loop_radius=0.35, cartoon_oval_length=1.3, cartoon_oval_width=0.35, cartoon_rect_width=0.5)
    k = 0
    for c in ctx.chains:
        res = residues(f"CARTOON and chain {c} and name CA")
        stops = [pal[(k + j) % len(pal)] for j in range(4)]
        k += 4
        nb = 24
        for i, (_, resi, _) in enumerate(res):
            t = i / max(1, len(res) - 1) * (len(stops) - 1)
            j = min(int(t), len(stops) - 2)
            b = int(i / max(1, len(res)) * nb)
            name = f"mf_p_{c}_{b}"
            if i == 0 or int((i - 1) / max(1, len(res)) * nb) != b:
                define(name, mix(stops[j], stops[j + 1], t - j))
            cmd.color(name, f"CARTOON and chain {c} and resi \\{resi}")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.6, 0.6, 0.6)), radius=0.25)
    if a.waters:
        cmd.create("WATERS", f"{BASE} and solvent")
        cmd.show("spheres", "WATERS")
        cmd.set("sphere_scale", 0.4, "WATERS")
        cmd.color(define("mf_water", (0.95, 0.4, 0.4)), "WATERS")
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.0, ray_trace_color=define("mf_outline", (0.35, 0.35, 0.35)),
                    ambient=0.4, direct=0.6, reflect=0.3, specular=0.35, spec_reflect=0.35, shininess=25,
                    light_count=3, ray_shadows=0, depth_cue=0)
    return ctx.view("overview", margin=1.05)


def _field_surface(ctx, kind):
    a = ctx.a
    region = a.region or ("site" if ctx.lig else "all")
    if region == "site":
        need_partner(ctx, "electrostatic/lipophilic --region site")
    surf_sel = ctx.receptor + (" and (all within 20 of lig_core)" if region == "site" else "")
    cmd.create("SURFACE", ctx.receptor)
    cmd.hide("everything", "SURFACE")
    cmd.show("surface", "SURFACE")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.22 if region == "site" else 0.0), "SURFACE")
    path = f"{a.prefix}_{kind}.dx"
    if kind == "esp":
        fields.coulombic(f"{BASE} and not solvent", surf_sel, path)
        lv = (-a.range, 0, a.range) if a.range else (-10, 0, 10)
        cols, lab = [(0.85, 0.1, 0.1), (1, 1, 1), (0.1, 0.1, 0.85)], "kcal per (mol·e)"
    else:
        fld, lo, h = fields.lipophilic(ctx.receptor, surf_sel, path)
        lv = (-a.range, 0, a.range) if a.range else fields.auto_levels(fld, lo, h, "SURFACE", f"{a.prefix}_tmp.obj")
        cols, lab = [(0.1, 0.55, 0.55), (1, 1, 1), (0.7, 0.55, 0.1)], "Molecular lipophilicity potential"
    fields.surface_ramp("SURFACE", path, f"MAP_{kind.upper()}", f"RAMP_{kind.upper()}", lv, cols)
    fields.colorbar(f"{a.prefix}_colorbar", list(lv), cols, lab)
    if region == "site":
        cmd.create("CARTOON", f"{ctx.receptor} and (all within 20 of lig_core)")
        cmd.hide("everything", "CARTOON")
        cmd.show("cartoon", "CARTOON")
        cmd.set("cartoon_loop_radius", 0.2)
        cmd.set("cartoon_transparency", 0.3, "CARTOON")
        cmd.color("grey60", "CARTOON")
        site_sticks(ctx, a.site_cutoff, radius=0.16)
        cmd.hide("sticks", "SITE_STICKS and not sidechain")
        color_by_element("SITE_STICKS", "grey40")
    if ctx.lig:
        show_partner(ctx, ctx.partner_rgb((0.2, 0.55, 0.55) if kind == "lipo" else (0.3, 0.3, 0.3)), radius=0.25)
    finish_lighting(ambient=0.55, direct=0.4, specular=0.0, ray_trace_mode=1, ray_trace_gain=0.0,
                    ray_trace_color="grey20", ray_shadows=0)
    if region == "site":
        view = ctx.view("site", frame_sel=f"SITE_STICKS or {partner_core_obj()}", margin=1.05, clip=a.clip)
        write_site(ctx, "SITE_STICKS", a.site_cutoff)
        site_labels(ctx, view, "SITE_STICKS", lambda c: (0.05, 0.05, 0.05))
        return view
    return ctx.view("overview", margin=1.05)


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1800, 1800), glycans=True)
def electrostatic(ctx):
    """Surface colored by Coulomb potential (ε = 4r), red -10 / white / blue +10
    kcal/(mol·e) (--range). Around the partner (--region site) or the whole molecule (--region all)."""
    return _field_surface(ctx, "esp")


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1800, 1800), glycans=True)
def lipophilic(ctx):
    """Surface colored by lipophilicity (MLP): hydrophilic teal / white / hydrophobic gold.
    Range set automatically from the visible surface (--range to fix it)."""
    return _field_surface(ctx, "lipo")


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1600, 1600), glycans=True)
def property_(ctx):
    """Cartoon colored by a per-residue property: B-factor (default) or a --values file
    'chain resi value' (e.g. ConSurf conservation), teal-white-maroon palette + color bar."""
    a = ctx.a
    cmd.create("CARTOON", ctx.receptor if ctx.kind == "small" else f"{BASE} and polymer")
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    if a.values:
        cmd.alter("CARTOON", "b = -999")
        for line in open(a.values):
            t = line.split()
            if len(t) >= 3 and not line.startswith("#"):
                cmd.alter(f"CARTOON and chain {t[0]} and resi \\{t[1]}", f"b = {float(t[2])}")
        known = "CARTOON and not b = -999"
        prop_label = a.prop_label or "Value"
    else:
        bs = {}  # per-residue B-factor = mean over its atoms
        cmd.iterate("CARTOON", "bs.setdefault((chain, resi), []).append(b)", space={"bs": bs})
        for (c, r), v in bs.items():
            cmd.alter(f"CARTOON and chain {c} and resi \\{r}", f"b = {sum(v) / len(v)}")
        known = "CARTOON"
        prop_label = a.prop_label or "B-factor (Å²)"
    vals = []
    cmd.iterate(f"({known}) and name CA", "vals.append(b)", space={"vals": vals})
    lo, hi = a.prop_range if a.prop_range else (float(np.percentile(vals, 5)), float(np.percentile(vals, 95)))
    cols = [(0.05, 0.42, 0.45), (0.55, 0.78, 0.8), (1, 1, 1), (0.9, 0.62, 0.75), (0.5, 0.1, 0.3)]
    n = 32
    for i in range(n):
        x = i / (n - 1) * (len(cols) - 1)
        j = min(int(x), len(cols) - 2)
        define(f"mf_prop_{i}", mix(cols[j], cols[j + 1], x - j))
    step = (hi - lo) / n
    for i in range(n):
        v0 = -1e9 if i == 0 else lo + i * step
        v1 = 1e9 if i == n - 1 else lo + (i + 1) * step
        cmd.color(f"mf_prop_{i}", f"({known}) and b > {v0} and b < {v1} or (({known}) and b = {v0})")
    if a.values:
        cmd.color("grey80", "CARTOON and b = -999")
    key = ctx.sel("key")
    if key:
        cmd.show("sticks", f"byres (CARTOON in ({key})) and {SIDECHAIN}")
        cmd.util.cnc(f"byres (CARTOON in ({key}))")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.1, 0.8, 0.85)), radius=0.22)
    finish_lighting(cartoon_loop_radius=0.3, ray_trace_mode=1, ray_trace_gain=0.0, ray_trace_color="grey25", ambient=0.55)
    fields.colorbar(f"{a.prefix}_colorbar", list(np.linspace(lo, hi, len(cols))), cols, prop_label, ticks=[lo, hi],
                    ticklabels=a.prop_ticks.split(",") if a.prop_ticks else [f"{lo:.0f}", f"{hi:.0f}"])
    view = ctx.view("overview", margin=1.04, fit_box=True)
    if key and ctx.labels_on(True):
        lab = residues(f"byres (CARTOON in ({key})) and name CA")
        add_labels("LABELS", residue_label_items(lab, "CARTOON", cmd.get_coords("CARTOON").mean(0), lambda c: (0, 0, 0)),
                   R_of(view), label_size(view, 0.03), view=view)
    return view


# AlphaFold confidence bands (pLDDT), with the AlphaFold Database colors
PLDDT_BANDS = [(90, "Very high", "#0053D6"), (70, "Confident", "#65CBF3"), (50, "Low", "#FFDB13"),
               (0, "Very low", "#FF7D45")]


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1600, 1600), glycans=True)
def plddt(ctx):
    """AlphaFold model colored by confidence (pLDDT, read from the B-factor column) with the
    AlphaFold Database colors: >90 dark blue, 70-90 light blue, 50-70 yellow, <50 orange.
    Discrete color bar; --key residues in sticks."""
    a = ctx.a
    cmd.create("CARTOON", ctx.receptor if ctx.kind == "small" else f"{BASE} and polymer")
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    vals = []
    cmd.iterate("CARTOON and name CA", "vals.append(b)", space={"vals": vals})
    if vals and (max(vals) > 100 or min(vals) < 0 or max(vals) <= 1):
        print(f"[PROTFIG] warning: B-factors range {min(vals):.1f}-{max(vals):.1f}; this does not look like "
              "pLDDT (0-100). Is this an AlphaFold model?", file=sys.stderr)
    bs = {}  # per-residue pLDDT = mean over its atoms (AlphaFold 3 writes per-atom values)
    cmd.iterate("CARTOON", "bs.setdefault((chain, resi), []).append(b)", space={"bs": bs})
    for (c, r), v in bs.items():
        cmd.alter(f"CARTOON and chain {c} and resi \\{r}", f"b = {sum(v) / len(v)}")
    for i, (cut, _, hexc) in enumerate(PLDDT_BANDS):
        upper = PLDDT_BANDS[i - 1][0] if i else None
        sel = f"CARTOON and (b > {cut} or b = {cut})" + (f" and b < {upper}" if upper is not None else "")
        cmd.color(define(f"mf_plddt_{i}", rgb_of(hexc)), sel)
    key = ctx.sel("key")
    if key:
        cmd.show("sticks", f"byres (CARTOON in ({key})) and {SIDECHAIN}")
        cmd.util.cnc(f"byres (CARTOON in ({key}))")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.6, 0.6, 0.6)), radius=0.22)
    finish_lighting(cartoon_loop_radius=0.3, ray_trace_mode=1, ray_trace_gain=0.0, ray_trace_color="grey25", ambient=0.55)
    cols = [rgb_of(h) for _, _, h in reversed(PLDDT_BANDS)]
    levels, stops = [], []
    for i, col in enumerate(cols):  # equal-width blocks; repeated stops -> sharp steps between bands
        levels += [i, i + 1]
        stops += [col, col]
    fields.colorbar(f"{a.prefix}_colorbar", levels, stops, a.prop_label or "pLDDT", ticks=list(range(5)),
                    ticklabels=["0", "50", "70", "90", "100"])
    view = ctx.view("overview", margin=1.05, fit_box=True)
    if key and ctx.labels_on(True):
        lab = residues(f"byres (CARTOON in ({key})) and name CA")
        add_labels("LABELS", residue_label_items(lab, "CARTOON", cmd.get_coords("CARTOON").mean(0), lambda c: (0, 0, 0)),
                   R_of(view), label_size(view, 0.03), view=view)
    return view


@preset(partner="required", kinds=("small", "peptide"), size=(1800, 1800))
def density(ctx):
    """Site with an electron density mesh (--map 2Fo-Fc .ccp4, e.g. from PDBe; --map_level σ)
    over the partner and the site residues; gold H-bonds; --labels on: H-bonded residues labeled."""
    a = ctx.a
    finish_lighting(**CLOSEUP_CARTOON)
    need_partner(ctx, "density")
    if not a.map:
        raise ProtFigError("preset 'density' needs --map (2Fo-Fc map, e.g. https://www.ebi.ac.uk/pdbe/entry-files/download/<id>.ccp4)")
    trim_copies(ctx)
    ctx.chain_colors("muted")
    cut = a.site_cutoff
    cmd.create("CARTOON", f"{ctx.receptor} and (all within 20 of lig_core)")
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.6), "CARTOON")
    site_sticks(ctx, cut, radius=0.25)
    for c in ctx.chains:
        cmd.color(f"mf_light_{c}", f"CARTOON and chain {c}")
        color_by_element(f"SITE_STICKS and chain {c}", f"mf_chain_{c}")
    show_partner(ctx, ctx.partner_rgb((0.95, 0.75, 0.4)), radius=0.25)
    cmd.load(a.map, "MAP")
    cmd.map_double("MAP")  # PDBe maps come on a coarse grid (~1 Å): smooths the mesh
    # the region is taken from the original object: copies (cmd.create) lack the crystallographic
    # symmetry PyMOL uses to expand a unit-cell map out to the model
    cmd.isomesh("MESH", "MAP", a.map_level, f"ligand or ({BASE} in (SITE_STICKS and (sidechain or name CA)))", carve=1.8)
    cmd.show("mesh", "MESH")
    cmd.refresh()  # in headless mode the mesh is computed lazily
    cmd.color(define("mf_mesh", ctx.rgb("mesh_color", (0.62, 0.62, 0.62))), "MESH")
    cmd.disable("MAP")
    pairs = hbonds(ctx, dash_radius=0.07, default_color=HB_GOLD)
    refine_sticks(ctx, pairs, lambda c: f"mf_chain_{c}", radius=0.25)
    if pairs:
        cmd.set("dash_gap", 0.3, "HBONDS")
        cmd.set("dash_length", 0.15, "HBONDS")
    write_site(ctx, "SITE_STICKS", cut)
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.0, ambient=0.55, specular=0.1, fog_start=0.35)
    view = ctx.view("site", frame_sel=f"SITE_STICKS or {partner_core_obj()}", margin=1.0, clip=a.clip, fit_box=True,
                    unocclude=("lig_core", "CARTOON and name N+CA+C+O", "turn"))
    if ctx.labels_on(False):  # the mesh is busy: labels only on request
        site_labels(ctx, view, "SITE_STICKS", lambda c: readable(mix(ctx.colors.get(c, (0, 0, 0)), "black", 0.15)),
                    pairs=pairs, only_hbonded=True)
    return view


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1800, 1800), glycans=True)
def superpose(ctx):
    """Superposition: the input structure (reference) colored by chain + --ref (model)
    in light gray; RMSD in _rmsd.txt."""
    a = ctx.a
    if not a.ref:
        raise ProtFigError("preset 'superpose' needs --ref (structure to superpose)")
    ctx.chain_colors("muted")
    cmd.create("REFERENCE", f"{BASE} and polymer")
    cmd.hide("everything", "REFERENCE")
    cmd.show("cartoon", "REFERENCE")
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"REFERENCE and chain {c}")
    load_ref(a, "REF_MODEL")
    cmd.remove("REF_MODEL and (hydro or solvent)")
    cmd.remove("REF_MODEL and not alt ''+A")
    if a.ref_chains:  # only these chains of --ref are aligned and drawn
        keep = [c.strip() for c in a.ref_chains.split(",") if c.strip()]
        missing = [c for c in keep if c not in cmd.get_chains("REF_MODEL")]
        if missing:
            raise ProtFigError(f"--ref_chains: chain {', '.join(missing)} is not in {a.ref}. "
                               f"Available chains: {', '.join(cmd.get_chains('REF_MODEL'))}.")
        cmd.remove(f"REF_MODEL and not chain {'+'.join(keep)}")
    rms = cmd.super("REF_MODEL and polymer and name CA", "REFERENCE and name CA")
    cmd.hide("everything", "REF_MODEL")
    cmd.show("cartoon", "REF_MODEL and polymer")
    cmd.color(define("mf_model", ctx.rgb("ref_color", (0.88, 0.88, 0.88))), "REF_MODEL")
    # thinner, translucent model cartoon: where both overlap the reference shows without mottling
    for k, v in (("cartoon_oval_length", 1.0), ("cartoon_oval_width", 0.2), ("cartoon_rect_length", 1.2),
                 ("cartoon_rect_width", 0.35), ("cartoon_loop_radius", 0.18), ("cartoon_transparency", 0.3)):
        cmd.set(k, v, "REF_MODEL")
    if ctx.kind == "small":
        show_partner(ctx, ctx.partner_rgb((0.95, 0.75, 0.4)), radius=0.3)
    with open(f"{a.prefix}_rmsd.txt", "w") as fh:
        fh.write(f"# super(REF_MODEL={a.ref} CA, REFERENCE={a.input} CA)\nRMSD = {rms[0]:.3f} A over {rms[1]} CA\n")
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.0, ray_trace_color="grey30", cartoon_fancy_helices=1,
                    ambient=0.5, specular=0.2)
    return ctx.view("overview", margin=1.04, fit_box=True)


# ----------------------------------------------------------------------- new presets
@preset(partner="required", kinds=("small", "peptide", "protein"), size=(1600, 1350), glycans=True)
def footprint(ctx):
    """Binding footprint on the receptor surface: receptor as a gray surface, pocket (residues
    within --pocket_cutoff Å of the partner, default 4.5) with surface colored by element (pink C;
    --palette: its first color)
    + outlined sticks, partner hidden (--show_partner on to show it). Camera facing the pocket,
    whole receptor in view, N-terminus up (--focus)."""
    a = ctx.a
    need_partner(ctx, "footprint")
    cut = ctx.opt("pocket_cutoff", 4.5)
    pocket_sel = f"byres ({ctx.receptor} within {cut} of lig_core)"
    pocket_rgb, partner_rgb = palette_pair(ctx, "lightpink", "palecyan")
    c_pocket = define("mf_pocket_c", ctx.rgb("site_color", pocket_rgb))
    cmd.create("RECEPTOR_REST", f"{ctx.receptor} and not ({pocket_sel})")
    cmd.create("POCKET", pocket_sel)
    cmd.create("POCKET_STICKS", pocket_sel)
    cmd.hide("everything", "RECEPTOR_REST or POCKET or POCKET_STICKS")
    cmd.show("surface", "RECEPTOR_REST")
    cmd.color(define("mf_rest", ctx.rgb("receptor_color", (0.7, 0.7, 0.7))), "RECEPTOR_REST")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.25), "RECEPTOR_REST")
    cmd.show("surface", "POCKET")
    color_by_element("POCKET", c_pocket)
    cmd.set("transparency", 0.45, "POCKET")
    cmd.show("sticks", "POCKET_STICKS")
    cmd.set("stick_radius", 0.3, "POCKET_STICKS")
    color_by_element("POCKET_STICKS", c_pocket)
    show_partner(ctx, ctx.partner_rgb(partner_rgb), radius=0.3)
    if ctx.opt("show_partner", "off") != "on":
        cmd.disable("PARTNER_MOL")
    res = residues("POCKET and name CA")
    with open(f"{a.prefix}_site.txt", "w") as fh:
        fh.write(f"# pocket = byres(receptor within {cut} A of the partner), heavy atoms; n = {len(res)}\n")
        fh.write(", ".join(f"{c}:{n}{r}" for c, r, n in res) + "\n")
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.3, ambient=0.5, direct=0.45, specular=0, spec_reflect=0,
                    ray_shadows=1, fog_start=0.45)
    if a.view or (ctx.session_view and not a.auto_view):
        return ctx.view("overview")
    c, R = facing_R(ctx.receptor, "POCKET", ctx.aspect, up_sel=ctx.sel("focus"))
    R = adjust(R, a)
    fr = ctx.receptor + (" or ligand" if ctx.opt("show_partner", "off") == "on" else "")
    return frame(R, c, cmd.get_coords(fr), ctx.aspect, margin=1.04, fit_box=True)


@preset(partner="required", kinds=("peptide", "protein"), size=(1800, 1200), glycans=True)
def profile(ctx):
    """Side view of the whole complex: receptor as a solid light-pink surface, peptide/protein partner
    as a pale-cyan cartoon on top, black outline. The camera looks along the
    binding interface, with the partner above the receptor and the long axis of the complex
    horizontal. --palette (light tones of its first two colors) or --receptor_color / --partner_color to
    recolor; --partner_sticks on adds side chains."""
    a = ctx.a
    need_partner(ctx, "profile")
    if a.partner_sticks is None:
        a.partner_sticks = "off"
    cmd.create("SURFACE", ctx.receptor)
    cmd.hide("everything", "SURFACE")
    cmd.show("surface", "SURFACE")
    rec_rgb, par_rgb = palette_pair(ctx, "lightpink", "palecyan")
    cmd.color(define("mf_rec", ctx.rgb("receptor_color", rec_rgb)), "SURFACE")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.0), "SURFACE")
    par_c = rgb_of(a.partner_color) if a.partner_color else par_rgb
    show_partner(ctx, par_c, radius=0.3, cartoon=par_c)
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.1, ambient=0.5, direct=0.5, specular=0, spec_reflect=0,
                    cartoon_transparency=0, cartoon_loop_radius=0.3, cartoon_fancy_helices=0, ray_shadows=0)
    if a.view or (ctx.session_view and not a.auto_view):
        return ctx.view("overview")
    # partner above the receptor (up = receptor centroid -> partner), partner long axis horizontal,
    # camera raised 30 degrees so a partner lying in a groove is not hidden by its walls
    L, P = cmd.get_coords("lig_core"), cmd.get_coords(ctx.receptor)
    up = L.mean(0) - P.mean(0)
    up /= np.linalg.norm(up)
    both = np.vstack([L, P])
    ax = np.linalg.svd(L - L.mean(0))[2][0] if len(L) > 2 else np.linalg.svd(both - both.mean(0))[2][0]
    x = ax - np.dot(ax, up) * up
    x /= np.linalg.norm(x)
    t = np.radians(30.0)
    z = np.cos(t) * np.cross(x, up) + np.sin(t) * up
    y = up - np.dot(up, z) * z
    y /= np.linalg.norm(y)
    R = adjust(np.vstack([np.cross(y, z), y, z]), a)
    xyz = np.vstack([cmd.get_coords(ctx.receptor), cmd.get_coords("PARTNER_MOL")])
    cam = xyz @ R.T
    center = ((cam.max(axis=0) + cam.min(axis=0)) / 2) @ R
    return frame(R, center, xyz, ctx.aspect, margin=1.05)


def _overview_assembly(ctx):
    """Oligomeric assembly as a surface, one color per chain (default: alternating light/dark
    grays) with black outline. --palette / --colors to recolor; --chains to select chains."""
    ctx.chain_colors("grays")
    cmd.create("SURFACE", f"{BASE} and polymer")
    cmd.hide("everything", "SURFACE")
    cmd.show("surface", "SURFACE")
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"SURFACE and chain {c}")
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.1, ambient=0.35, direct=0.65, reflect=0.3, specular=0.15,
                    spec_reflect=0.1, ray_shadows=0, depth_cue=0, surface_quality=0)
    a = ctx.a
    if a.view or (ctx.session_view and not a.auto_view):
        return ctx.view("overview")
    # side view: the symmetry axis (normal to the plane of the chain centers; with < 3 chains, the
    # short axis) vertical, the assembly's widest extent horizontal
    xyz = cmd.get_coords("SURFACE")
    c = xyz.mean(0)
    cen = np.array([cmd.get_coords(f"SURFACE and chain {ch}").mean(0) for ch in cmd.get_chains("SURFACE")])
    V = np.linalg.svd(xyz - c)[2]
    up = np.linalg.svd(cen - cen.mean(0))[2][2] if len(cen) >= 3 else V[2]
    x = V[0] - np.dot(V[0], up) * up
    if np.linalg.norm(x) < 0.2:
        x = V[1] - np.dot(V[1], up) * up
    x /= np.linalg.norm(x)
    R = adjust(np.vstack([x, up, np.cross(x, up)]), a)
    return frame(R, c, xyz, ctx.aspect, margin=1.03, fit_box=True)


@preset(partner="none", kinds=("peptide", "protein"), size=(1800, 1000))
def membrane(ctx):
    """Membrane protein in its bilayer (MD system, CHARMM-GUI or MemProtMD): protein as cartoon
    colored by chain (--repr spheres / surface), bilayer as spheres (C gray, O red, P orange, N blue) trimmed to
    --membrane_radius Å around the protein with the front lipids removed, partner (--partner,
    e.g. a nanobody) in cartoon. Side view with the membrane normal vertical. --waters: the MD water
    as a translucent pale-blue surface on both faces."""
    a = ctx.a
    lip = f"{BASE} and resn {'+'.join(LIPIDS)}" + (f" or ({BASE} and ({a.lipids}))" if a.lipids else "")
    if cmd.count_atoms(lip) == 0:
        raise ProtFigError("preset 'membrane' needs lipids in the structure (MD system: CHARMM-GUI, "
                 "MemProtMD...). If your lipids have a different name, use --lipids \"resn XXX\".")
    ctx.chain_colors("membrane")
    cmd.create("PROTEIN", ctx.receptor)
    cmd.hide("everything", "PROTEIN")
    rep = a.repr or "cartoon"  # --repr spheres / surface
    cmd.show(rep, "PROTEIN")
    if rep == "cartoon":
        finish_lighting(cartoon_loop_radius=0.35, cartoon_oval_length=1.3)
    for c in ctx.chains:
        cmd.color(f"mf_chain_{c}", f"PROTEIN and chain {c}")
    if ctx.lig:
        orange = rgb_of(a.partner_color) if a.partner_color else (0.95, 0.6, 0.15)
        show_partner(ctx, orange, radius=0.25, cartoon=orange)
        cmd.hide("sticks", "PARTNER_MOL")
    c, R = membrane_R("PROTEIN", lip, up_sel="ligand" if ctx.lig else ctx.sel("focus"))
    R = adjust(R, a)
    # lipids: within the lateral radius and behind the protein's mid-plane
    rad = ctx.opt("membrane_radius", 40.0)
    cmd.create("_mf_lip", f"byres (({lip}) and (all within {rad} of PROTEIN))")
    cmd.alter_state(1, "_mf_lip", "b = (x - ox) * r0 + (y - oy) * r1 + (z - oz) * r2",
                    space={"ox": c[0], "oy": c[1], "oz": c[2], "r0": R[2][0], "r1": R[2][1], "r2": R[2][2]})
    cmd.create("LIPIDS", "byres (_mf_lip and b < 2)")
    cmd.delete("_mf_lip")
    cmd.hide("everything", "LIPIDS")
    cmd.show("spheres", "LIPIDS")
    cmd.color(define("mf_lip_c", ctx.rgb("lipid_color", (0.62, 0.62, 0.62))), "LIPIDS and elem C")
    cmd.color(define("mf_lip_o", (0.9, 0.1, 0.1)), "LIPIDS and elem O")
    cmd.color(define("mf_lip_p", (1.0, 0.6, 0.0)), "LIPIDS and elem P")
    cmd.color(define("mf_lip_n", (0.2, 0.2, 0.9)), "LIPIDS and elem N")
    finish_lighting(ray_trace_mode=0, ambient=0.45, direct=0.55, specular=0.25, spec_reflect=0.2, shininess=20,
                    ray_shadows=0, depth_cue=0, sphere_quality=2)
    if a.waters:  # MD water: a thin layer on both faces, as a translucent pale-blue surface, behind the protein
        water = f"{BASE} and resn {'+'.join(WATERS)}"
        layer = str(ctx.opt("water_layer", "10")).strip().lower()  # Å from the bilayer, or "all" (the whole box)
        near = "" if layer == "all" else f" and ((all within {float(layer)} of LIPIDS) or (all within {float(layer)} of PROTEIN))"
        cmd.create("_mf_wat", f"byres (({water}){near})")
        if cmd.count_atoms("_mf_wat"):
            cmd.alter_state(1, "_mf_wat", "b = (x - ox) * r0 + (y - oy) * r1 + (z - oz) * r2",
                            space={"ox": c[0], "oy": c[1], "oz": c[2], "r0": R[2][0], "r1": R[2][1], "r2": R[2][2]})
            cmd.create("WATER", "byres (_mf_wat and b < 2)")
            cmd.hide("everything", "WATER")
            cmd.show("surface", "WATER")
            cmd.color(define("mf_water", ctx.rgb("water_color", (0.62, 0.78, 1.0))), "WATER")
            cmd.set("transparency", 0.6, "WATER")
        cmd.delete("_mf_wat")
    fr = "PROTEIN or LIPIDS" + (" or PARTNER_MOL" if ctx.lig else "") + (" or WATER" if "WATER" in cmd.get_names("objects") else "")
    return frame(R, c, cmd.get_coords(fr), ctx.aspect, margin=1.02, fit_box=True)
