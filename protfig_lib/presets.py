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
import sys

import numpy as np
from pymol import cmd

from . import fields
from .core import (BASE, GLYCANS, LIPIDS, METALS, ProtFigError, PALETTES, add_labels, adjust, color_by_element, define, facing_R, frame,
                   membrane_R, mix, overview_R, polar_contacts, residue_label_items, residue_sel, residues, rgb_of, site_R,
                   vis_height, write_contacts)

HB_YELLOW = (1.0, 0.88, 0.0)  # H-bonds: yellow, readable on white and on most surfaces
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
    def view(self, kind, center_sel=None, frame_sel=None, margin=1.06, clip=None, fit_box=False):
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
            R = adjust(R, a)
            front = None
            if clip is not None:
                front = float(np.max((cmd.get_coords("lig_core") - c) @ R[2])) + clip
            return frame(R, c, cmd.get_coords(frame_sel or "lig_core"), self.aspect, margin=margin, front_z=front)
        sel = center_sel or f"{BASE} and polymer"
        up = self.sel("focus") or (f"{BASE} and chain {self.chains[0]}" if len(self.chains) > 1 else None)
        c, R = overview_R(sel, self.aspect, toward="ligand" if self.lig else None, up_sel=up)
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


def hbonds(ctx, partner_sel=None, dash_radius=0.08, name="HBONDS"):
    """Partner <-> receptor H-bonds (yellow by default; --hbond_color / --hbonds off)."""
    if not ctx.lig or ctx.a.hbonds == "off":
        return []
    pairs = polar_contacts(name, "PARTNER_MOL", partner_sel or ctx.receptor, cutoff=ctx.a.hbond_cutoff,
                           mode=ctx.a.hbond_mode)
    if pairs:
        cmd.set("dash_color", define("mf_hbond", ctx.rgb("hbond_color", HB_YELLOW)), name)
        cmd.set("dash_radius", dash_radius, name)
        write_contacts(f"{ctx.a.prefix}_hbonds.txt",
                       f"# partner-receptor H-bonds (N/O, <= {ctx.a.hbond_cutoff} A, {ctx.a.hbond_mode} criterion)\n",
                       pairs)
    return pairs


def distance_labels(ctx, pairs, view, size, obj="DIST_LABELS"):
    """'3.0 Å' labels at the midpoint of each H-bond."""
    items = []
    R = R_of(view)
    col = ctx.rgb("label_color", (0, 0, 0)) if ctx.a.label_color else (0, 0, 0)
    for (c1, n1, r1, a1), (c2, n2, r2, a2), d in pairs:
        p1 = cmd.get_coords(f"PARTNER_MOL and chain {c1} and resi \\{r1} and name {a1}")
        p2 = cmd.get_coords(f"{BASE} and chain {c2} and resi \\{r2} and name {a2}")
        if p1 is None or p2 is None:
            continue
        items.append((f"{d:.1f} Å", (p1[0] + p2[0]) / 2 + R[0] * size * 1.6, col))
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


def site_labels(ctx, view, site_obj, color_of, frac=0.03, font=7, partner_color=None, extra=(), pairs=()):
    """Labels for the site residues and, if the partner is a peptide/protein, for its side too.

    Ligand: all site residues are labeled. Peptide/protein (large interfaces): only residues
    forming H-bonds, on both sides (if there are none, the site residues). --label_sel picks
    the receptor residues by hand.
    """
    if not ctx.labels_on(True):
        return
    hb_rec = {(c, r) for _, (c, _, r, _), _ in pairs}
    hb_par = {(c, r) for (c, _, r, _), _, _ in pairs}
    res = residues(f"({site_obj}) and name CA")
    if ctx.a.label_sel:
        res = residues(f"({site_obj}) and name CA and ({ctx.a.label_sel})")
    elif ctx.kind in ("peptide", "protein") and hb_rec:
        res = [x for x in res if (x[0], x[1]) in hb_rec]
    res += list(extra)
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
    add_labels("LABELS", items, R_of(view), label_size(view, frac), font=font, view=view, avoid=avoid)


def site_sticks(ctx, cutoff, obj="SITE_STICKS", radius=0.22):
    cmd.create(obj, ctx.site(cutoff))
    cmd.hide("everything", obj)
    cmd.show("sticks", f"{obj} and {SIDECHAIN}")
    cmd.show("sticks", f"{obj} and resn GLY+PRO and name CA+N+C")
    cmd.set("stick_radius", radius, obj)
    return obj


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


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1500, 2000), glycans=True)
def overview(ctx):
    """Whole complex: cartoon per chain + translucent surface of the same color, ligand in
    sticks and chain names (--chain_names). Works for any complex."""
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
    view = ctx.view("overview", margin=1.0)
    names = ctx.chain_names()
    if (names or a.partner_name) and ctx.labels_on(True):
        # each chain's name to the left of its on-screen edge, level with its center
        R = R_of(view)
        size = label_size(view, 0.05)
        zf = (cmd.get_coords("CARTOON") @ R.T)[:, 2].max() - 4  # +3 from label_position: inside the clip
        items = []
        for c in ctx.chains:
            if c not in names:
                continue
            cam = cmd.get_coords(f"CARTOON and chain {c}") @ R.T
            cy = np.median(cam[:, 1])
            band = cam[np.abs(cam[:, 1] - cy) < 6]
            x = (band if len(band) else cam)[:, 0].min() - 0.33 * size * len(names[c]) - 2.0
            items.append((names[c], np.array([x, cy, zf]) @ R, mix(ctx.colors[c], "black", 0.1)))
        if ctx.kind == "small" and a.partner_name:
            L = cmd.get_coords("PARTNER_MOL") @ R.T
            p = np.array([L[:, 0].min() - 0.33 * size * len(a.partner_name) - 3, L[:, 1].max() + 4, zf])
            items.append((a.partner_name, p @ R, ctx.partner_rgb((0.1, 0.7, 0.8))))
        add_labels("CHAIN_LABELS", items, R, size, font=5, view=view)
    return view


@preset(partner="required", kinds=("small", "peptide", "protein"), size=(2000, 1700))
def site(ctx):
    """Binding-site close-up: light cartoon per chain, site residues in sticks colored by
    chain, partner in sticks (ligand with a translucent blue surface) or in cartoon + interface
    sticks (peptide/protein), H-bonds and labels on both sides."""
    a = ctx.a
    need_partner(ctx, "site")
    ctx.chain_colors("vivid")
    cut = a.site_cutoff
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.35), "CARTOON")
    site_sticks(ctx, cut)
    for c in ctx.chains:
        cmd.color(f"mf_light_{c}", f"CARTOON and chain {c}")
        color_by_element(f"SITE_STICKS and chain {c}", f"mf_chain_{c}")
    pc = ctx.partner_rgb((0.0, 0.8, 0.82))
    show_partner(ctx, pc, radius=0.24)
    if ctx.kind != "small":
        cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.35), "PARTNER_MOL")
    if ctx.opt("metals", "on") == "on":
        show_metals()
    if ctx.kind == "small" and ctx.opt("partner_surface", "on") == "on":
        cmd.create("PARTNER_SURF", "ligand")
        cmd.hide("everything", "PARTNER_SURF")
        cmd.show("surface", "PARTNER_SURF")
        cmd.color(define("mf_ligsurf", ctx.rgb("partner_surface_color", (0.25, 0.35, 0.95))), "PARTNER_SURF")
        cmd.set("transparency", 0.55, "PARTNER_SURF")
    pairs = hbonds(ctx, dash_radius=0.08)
    for _, (c2, _, r2, _), _ in pairs:  # H-bonded residues are shown even if outside the cutoff
        cmd.create("SITE_STICKS", f"SITE_STICKS or ({BASE} and chain {c2} and resi \\{r2})")
    write_site(ctx, "SITE_STICKS", cut)
    frame_sel = f"SITE_STICKS or {partner_core_obj()}"
    view = ctx.view("site", frame_sel=frame_sel, margin=1.12, clip=a.clip)
    site_labels(ctx, view, "SITE_STICKS", lambda c: mix(ctx.colors.get(c, (0, 0, 0)), "black", 0.15),
                partner_color=mix(pc, "black", 0.25), pairs=pairs)
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


@preset(partner="required", kinds=("small", "peptide"), size=(2000, 1700))
def covalent(ctx):
    """Covalently bound ligand or peptide (retinal-Lys in opsins, peptide acyl-enzymes...): the
    partner and the receptor residue it is bonded to in sticks, in different colors, with the bond
    between them; the rest in cartoon (--interface_sticks on: interface residues in sticks too);
    labels on the two linked residues. The bond is taken from the file or, if missing, added for
    partner-receptor atoms within --covalent_cutoff (default 2.0 Å)."""
    a = ctx.a
    need_partner(ctx, "covalent")
    ctx.chain_colors("vivid")
    cut = ctx.opt("covalent_cutoff", 2.0)
    rec = f"({ctx.receptor}) and not resn {'+'.join(METALS)} and not hydro"
    links = covalent_links(f"{BASE} and ligand and not hydro", rec, cut)
    if not links:
        raise ProtFigError(f"covalent: no partner atom within {cut} Å of the receptor; "
                           "raise --covalent_cutoff or check --partner")
    for p, r, d, bonded in links:
        if not bonded:
            cmd.bond(f"{BASE} and index {p[4]}", f"{BASE} and index {r[4]}")
        print(f"covalent link: {p[0]}:{p[2]}{p[1]} {p[3]} - {r[0]}:{r[2]}{r[1]} {r[3]}  {d:.2f} Å"
              + ("" if bonded else "  (bond not in the file: added)"))
    rec_res = sorted({(r[0], r[1], r[2]) for _, r, _, _ in links})
    par_res = sorted({(p[0], p[1], p[2]) for p, _, _, _ in links})
    rec_sel = " or ".join(f"({BASE} and {residue_sel(c, i)})" for c, i, _ in rec_res)

    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.5), "CARTOON")  # the link is often inside a helix
    for c in ctx.chains:
        cmd.color(f"mf_light_{c}", f"CARTOON and chain {c}")
    frame_sel = "LINKED"
    if ctx.opt("interface_sticks", "off") == "on":
        site_sticks(ctx, a.site_cutoff)
        cmd.remove(f"SITE_STICKS and ({rec_sel})")
        for c in ctx.chains:
            color_by_element(f"SITE_STICKS and chain {c}", define(f"mf_site_{c}", mix(ctx.colors[c], "white", 0.55)))
        cmd.set("stick_radius", 0.16, "SITE_STICKS")  # thinner and lighter than the linked residues
        write_site(ctx, "SITE_STICKS", a.site_cutoff)
        frame_sel = "LINKED or SITE_STICKS"

    # partner + bonded receptor residue(s) in one object: PyMOL only draws bonds within an object
    pc = ctx.partner_rgb((0.0, 0.8, 0.82))
    cmd.create("LINKED", f"({BASE} and ligand and not hydro) or ({rec_sel})")
    cmd.hide("everything", "LINKED")
    cmd.show("sticks", "LINKED")
    cmd.set("stick_radius", 0.27, "LINKED")
    color_by_element("LINKED", define("mf_partner_c", pc))
    for c, i, _ in rec_res:
        color_by_element(f"LINKED and {residue_sel(c, i)}", f"mf_chain_{c}")
    if ctx.opt("metals", "off") == "on":
        show_metals()
    with open(f"{a.prefix}_covalent.txt", "w") as fh:
        fh.write("partner_atom\treceptor_atom\tdistance_A\tbond_in_file\n")
        for p, r, d, bonded in links:
            fh.write(f"{p[0]}:{p[2]}{p[1]}:{p[3]}\t{r[0]}:{r[2]}{r[1]}:{r[3]}\t{d:.2f}\t{'yes' if bonded else 'no'}\n")

    # labels at the end of each linked residue farthest from the bond, so the two sides do not
    # overlap; their anchors are framed too (temporary pseudoatoms) so no label falls off the image
    items = []
    if ctx.labels_on(True):
        center = cmd.get_coords(" or ".join(f"({BASE} and index {x[4]})" for p, r, _, _ in links for x in (p, r))).mean(0)
        for res_list, color_of in ((par_res, lambda c: mix(pc, "black", 0.3)),
                                   (rec_res, lambda c: mix(ctx.colors.get(c, (0, 0, 0)), "black", 0.2))):
            for c, i, n in res_list:
                xyz = cmd.get_coords(f"LINKED and {residue_sel(c, i)}")
                tip = xyz[np.argmax(np.linalg.norm(xyz - center, axis=1))]
                d = (tip - center) / (np.linalg.norm(tip - center) or 1.0)
                text = residue_label_items([(c, i, n)], "LINKED", center, color_of)[0][0]
                items.append((text, tip + 3.5 * d, color_of(c)))
        for k, (_, pos, _) in enumerate(items):
            cmd.pseudoatom("_mf_label_pts", pos=[float(v) for v in pos], name=f"L{k}")
        frame_sel += " or _mf_label_pts"
    # default clip: helices or strands in front of the link would hide it
    view = ctx.view("site", frame_sel=frame_sel, margin=1.3, clip=a.clip if a.clip is not None else 6.0)
    cmd.delete("_mf_label_pts")
    if items:
        add_labels("LABELS", items, R_of(view), label_size(view, 0.035), font=7, view=view)
    return view


@preset(partner="required", kinds=("small", "peptide"), size=(1500, 1700))
def pocket(ctx):
    """Surface pocket (per chain) in a "cutaway" view: only the surface behind the partner, which
    lies in its groove; site sticks. For ligands and peptides."""
    a = ctx.a
    need_partner(ctx, "pocket")
    ctx.chain_colors("vivid")
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


@preset(partner="required", kinds=("small", "peptide", "protein"), size=(2000, 2000))
def site_faded(ctx):
    """Site on a ghost-gray cartoon: residues in yellow sticks, ligand in ball-and-stick
    (peptide/protein in cartoon + sticks), H-bonds with distances in Å. Optional --ref
    (ligand from a 2nd structure in translucent gray) and --mutation (translucent spheres)."""
    a = ctx.a
    need_partner(ctx, "site_faded")
    cut = a.site_cutoff
    cmd.create("CARTOON", f"{ctx.receptor} and (all within 22 of lig_core)")
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.color(define("mf_ghost", ctx.rgb("receptor_color", "grey70")), "CARTOON")
    cmd.set("cartoon_transparency", ctx.opt("cartoon_transparency", 0.72), "CARTOON")
    site_sticks(ctx, cut, radius=0.17)
    color_by_element("SITE_STICKS", define("mf_site_c", ctx.rgb("site_color", (1.0, 0.95, 0.45))))
    cmd.color(define("mf_N_soft", (0.55, 0.55, 0.95)), "SITE_STICKS and elem N")
    cmd.color(define("mf_O_soft", (0.95, 0.45, 0.45)), "SITE_STICKS and elem O")
    show_partner(ctx, ctx.partner_rgb((0.2, 0.9, 0.95)), style="ballstick", radius=0.38)
    if ctx.kind == "small":
        cmd.color(define("mf_N_lig", (0.05, 0.05, 0.9)), "PARTNER_MOL and elem N")
        cmd.color(define("mf_O_lig", (0.95, 0.4, 0.4)), "PARTNER_MOL and elem O")
        cmd.color(define("mf_F_lig", (0.75, 0.95, 0.95)), "PARTNER_MOL and elem F")
        cmd.color(define("mf_Cl_lig", (0.1, 0.75, 0.1)), "PARTNER_MOL and elem Cl")
    else:
        cmd.set("cartoon_transparency", 0.4, "PARTNER_MOL")
    mut = ctx.sel("mutation")
    if mut:
        cmd.create("MUT_STICKS", f"byres ({mut})")
        cmd.hide("everything", "MUT_STICKS")
        cmd.show("sticks", f"MUT_STICKS and {SIDECHAIN}")
        cmd.show("spheres", "MUT_STICKS and sidechain")
        color_by_element("MUT_STICKS", define("mf_mut_c", (0.45, 0.72, 0.35)))
        cmd.set("sphere_color", define("mf_mut_sph", ctx.rgb("mutation_color", (0.3, 0.6, 0.95))), "MUT_STICKS")
        cmd.set("sphere_transparency", 0.45, "MUT_STICKS")
        cmd.set("stick_radius", 0.2, "MUT_STICKS")
        cmd.remove("SITE_STICKS and (SITE_STICKS in MUT_STICKS)")
    if a.ref and ctx.kind == "small":
        cmd.load(a.ref, "REF_FULL")
        cmd.remove("REF_FULL and hydro")
        cmd.remove("REF_FULL and not alt ''+A")
        # superpose on the CAs around the site (not the whole protein)
        cmd.super("REF_FULL and polymer and name CA",
                  f"{ctx.receptor} and name CA and (all within 15 of ligand)", cycles=5)
        lig_resn = set()
        cmd.iterate("PARTNER_MOL", "lig_resn.add(resn)", space={"lig_resn": lig_resn})
        ref_lig = f"REF_FULL and resn {'+'.join(lig_resn)} and (all within 4 of PARTNER_MOL)"
        if cmd.count_atoms(ref_lig):
            cmd.create("REF_PARTNER", f"byres ({ref_lig})")
            cmd.hide("everything", "REF_PARTNER")
            cmd.show("sticks", "REF_PARTNER")
            cmd.show("spheres", "REF_PARTNER")
            cmd.set("sphere_scale", 0.42, "REF_PARTNER")
            cmd.set("stick_radius", 0.24, "REF_PARTNER")
            cmd.color("grey50", "REF_PARTNER")
            cmd.set("sphere_transparency", 0.6, "REF_PARTNER")
            cmd.set("stick_transparency", 0.6, "REF_PARTNER")
        cmd.delete("REF_FULL")
    pairs = hbonds(ctx, dash_radius=0.11)
    if pairs:
        cmd.set("dash_gap", 0.25, "HBONDS")
        cmd.set("dash_length", 0.3, "HBONDS")
    write_site(ctx, "SITE_STICKS", cut)
    finish_lighting(ambient=0.62, direct=0.35, spec_reflect=0.25, specular=0.4, shininess=30, ray_shadows=0)
    frame_sel = f"SITE_STICKS or {partner_core_obj()}" + (" or MUT_STICKS" if mut else "")
    view = ctx.view("site", frame_sel=frame_sel, margin=0.98, clip=a.clip)
    gold = (0.78, 0.6, 0.1)
    site_labels(ctx, view, "SITE_STICKS", lambda c: gold, frac=0.04, font=7, partner_color=(0.1, 0.5, 0.6),
                extra=residues("MUT_STICKS and name CA") if mut else (), pairs=pairs)
    if pairs and ctx.labels_on(True) and ctx.opt("distances", "on") == "on":
        distance_labels(ctx, pairs, view, label_size(view, 0.032))
    return view


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(2400, 1500), glycans=True)
def outline(ctx):
    """Cartoon with black outline: base color (--receptor_color) + highlighted region (--highlight,
    --highlight_color), partner in cyan sticks (peptide/protein in cartoon), key residues
    (those forming H-bonds, plus --key) in yellow sticks. --highlight interface shades the
    interface on both sides (--partner_highlight_color) and shows key residues of both sides."""
    a = ctx.a
    want_iface = (a.highlight or "").strip().lower() == "interface"
    iface = want_iface and ctx.kind in ("peptide", "protein")  # both sides; a small ligand: receptor side only
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.color(define("mf_base", ctx.rgb("receptor_color", (0.25, 0.42, 0.75))), "CARTOON")
    hl_c = ctx.rgb("highlight_color", (0.98, 0.52, 0.1))
    hl = ctx.site(a.site_cutoff) if want_iface else ctx.sel("highlight")
    if hl:
        cmd.color(define("mf_hl", hl_c), f"CARTOON in ({hl})")
    finish_lighting(cartoon_fancy_helices=1, cartoon_highlight_color=-1, cartoon_loop_radius=0.3,
                    cartoon_oval_length=1.35)
    par_c = ctx.partner_rgb((0.35, 0.85, 0.9))
    if ctx.lig:
        if iface and a.partner_sticks is None:
            a.partner_sticks = "off"  # only the key residues in sticks
        show_partner(ctx, par_c, radius=0.22, cartoon=par_c if iface else None)
        if ctx.kind == "peptide" and not iface:
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
    view = ctx.view("overview", margin=1.04)
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


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1600, 1600), glycans=True)
def illustrative(ctx):
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


@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(2000, 1500), glycans=True)
def pastel(ctx):
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
    view = ctx.view("overview", margin=1.05)
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
    over the partner and the site residues."""
    a = ctx.a
    need_partner(ctx, "density")
    if not a.map:
        raise ProtFigError("preset 'density' needs --map (2Fo-Fc map, e.g. https://www.ebi.ac.uk/pdbe/entry-files/download/<id>.ccp4)")
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
    pairs = hbonds(ctx, dash_radius=0.07)
    if pairs:
        cmd.set("dash_gap", 0.3, "HBONDS")
        cmd.set("dash_length", 0.15, "HBONDS")
    write_site(ctx, "SITE_STICKS", cut)
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.0, ambient=0.55, specular=0.1, fog_start=0.35)
    return ctx.view("site", frame_sel=f"SITE_STICKS or {partner_core_obj()}", margin=1.0, clip=a.clip)


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
    cmd.load(a.ref, "REF_MODEL")
    cmd.remove("REF_MODEL and (hydro or solvent)")
    cmd.remove("REF_MODEL and not alt ''+A")
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
    return ctx.view("overview", margin=1.05)


# ----------------------------------------------------------------------- new presets
@preset(partner="required", kinds=("small", "peptide", "protein"), size=(1600, 1350), glycans=True)
def binding_surface(ctx):
    """Binding footprint on the receptor surface: receptor as a gray surface, pocket (residues
    within --pocket_cutoff Å of the partner, default 8) with surface colored by element (pink C)
    + outlined sticks, partner hidden (--show_partner on to show it). Camera facing the pocket,
    whole receptor in view, N-terminus up (--focus)."""
    a = ctx.a
    need_partner(ctx, "binding_surface")
    cut = ctx.opt("pocket_cutoff", 8.0)
    pocket_sel = f"byres ({ctx.receptor} within {cut} of lig_core)"
    c_pocket = define("mf_pocket_c", ctx.rgb("site_color", "lightpink"))
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
    show_partner(ctx, ctx.partner_rgb("palecyan"), radius=0.3)
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
    return frame(R, c, cmd.get_coords(fr), ctx.aspect, margin=1.08)


@preset(partner="required", kinds=("small", "peptide", "protein"), size=(1280, 1072))
def interface_closeup(ctx):
    """Outlined interface close-up: receptor in pink cartoon + interface sticks, partner in
    cyan cartoon + sticks, H-bonds in yellow. No labels by default (--labels on to add them)."""
    a = ctx.a
    need_partner(ctx, "interface_closeup")
    rec_c = ctx.rgb("receptor_color", "lightpink")
    par_c = rgb_of(a.partner_color) if a.partner_color else rgb_of("palecyan")
    cmd.create("CARTOON", ctx.receptor)
    cmd.hide("everything", "CARTOON")
    cmd.show("cartoon", "CARTOON")
    cmd.color(define("mf_rec", rec_c), "CARTOON")
    site_sticks(ctx, a.site_cutoff, radius=0.25)
    color_by_element("SITE_STICKS", "mf_rec")
    show_partner(ctx, par_c, radius=0.25, cartoon=par_c)
    pairs = hbonds(ctx, dash_radius=0.07)
    if pairs:
        cmd.set("dash_gap", 0.35, "HBONDS")
        cmd.set("dash_length", 0.45, "HBONDS")
    write_site(ctx, "SITE_STICKS", a.site_cutoff)
    finish_lighting(ray_trace_mode=1, ray_trace_gain=0.15, ambient=0.5, direct=0.45, specular=0, spec_reflect=0,
                    cartoon_transparency=0, cartoon_loop_radius=0.3, cartoon_fancy_helices=0, ray_shadows=0,
                    fog_start=0.6)
    view = ctx.view("site", frame_sel=f"SITE_STICKS or {partner_core_obj()}", margin=0.8, clip=a.clip)
    if ctx.labels_on(False):
        site_labels(ctx, view, "SITE_STICKS", lambda c: mix(rec_c, "black", 0.55),
                    partner_color=mix(par_c, "black", 0.55), pairs=pairs)
    return view


@preset(partner="required", kinds=("peptide", "protein"), size=(1800, 1200), glycans=True)
def profile(ctx):
    """Side view of the whole complex: receptor as a solid light-pink surface, peptide/protein partner
    as a pale-cyan cartoon on top, black outline. The camera looks along the
    binding interface, with the partner above the receptor and the long axis of the complex
    horizontal. --receptor_color / --partner_color to recolor; --partner_sticks on adds side chains."""
    a = ctx.a
    need_partner(ctx, "profile")
    if a.partner_sticks is None:
        a.partner_sticks = "off"
    cmd.create("SURFACE", ctx.receptor)
    cmd.hide("everything", "SURFACE")
    cmd.show("surface", "SURFACE")
    cmd.color(define("mf_rec", ctx.rgb("receptor_color", "lightpink")), "SURFACE")
    cmd.set("transparency", ctx.opt("surface_transparency", 0.0), "SURFACE")
    par_c = rgb_of(a.partner_color) if a.partner_color else rgb_of("palecyan")
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


@preset(partner="none", kinds=(), size=(1280, 960), glycans=True)
def oligomer(ctx):
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
    return ctx.view("overview", margin=1.03)


@preset(partner="none", kinds=("peptide", "protein"), size=(1800, 1000))
def membrane(ctx):
    """Membrane protein in its bilayer (MD system, CHARMM-GUI or MemProtMD): protein as spheres
    colored by chain, bilayer as spheres (C gray, O red, P orange, N blue) trimmed to
    --membrane_radius Å around the protein with the front lipids removed, partner (--partner,
    e.g. a nanobody) in cartoon. Side view with the membrane normal vertical."""
    a = ctx.a
    lip = f"{BASE} and resn {'+'.join(LIPIDS)}" + (f" or ({BASE} and ({a.lipids}))" if a.lipids else "")
    if cmd.count_atoms(lip) == 0:
        raise ProtFigError("preset 'membrane' needs lipids in the structure (MD system: CHARMM-GUI, "
                 "MemProtMD...). If your lipids have a different name, use --lipids \"resn XXX\".")
    ctx.chain_colors("membrane")
    cmd.create("PROTEIN", ctx.receptor)
    cmd.hide("everything", "PROTEIN")
    cmd.show("spheres", "PROTEIN")
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
    fr = "PROTEIN or LIPIDS" + (" or PARTNER_MOL" if ctx.lig else "")
    return frame(R, c, cmd.get_coords(fr), ctx.aspect, margin=1.0)
