"""
PROTFIG — publication-quality molecular figures from style presets (PyMOL, headless).

USAGE
  protfig <structure> --preset <style> [options]           (or: pymol -cq protfig.py -- ...)
  protfig --list                                           available presets
  protfig --config figures.toml                            one or more figures from a file

  <structure>: .pdb / .cif file, a PDB code (downloaded) or a .pse session.

EXAMPLES
  protfig 8G66.cif --preset site --chains B,C,A
  protfig 1BRS --preset interface_closeup --chains A,D --labels on
  protfig my_protein.pdb --preset illustrative --palette pastel --bg transparent
  protfig figures/site/8G66/8G66.pse --preset site_faded       (reuse view and selections)

OUTPUTS  (one folder per figure, grouped by style)
  figures/<preset>/<name>/<name>.png            final figure (300 dpi)
                         /<name>.pse            PyMOL session (open, adjust and reuse)
                         /<name>_run.txt        command, date, partner, chains and files
                         /<name>_view.txt       camera view (reusable with --view)
                         /<name>_*              extras: site residues, H-bonds, color bar,
                                                editable label layer...

ADJUSTING BY HAND IN PYMOL
  1) run automatically; 2) open <name>.pse in PyMOL, rotate/zoom and, optionally, define
  selections with these names: partner, site, highlight, key, mutation, focus; save;
  3) run again with that .pse as the structure (with the same or another style).
"""
import argparse
import datetime
import glob
import os
import shlex
import sys
import textwrap

# PyMOL defines __script__ when running `pymol -cq protfig.py`; plain python defines __file__
SCRIPT = os.path.abspath(globals().get("__script__") or globals().get("__file__") or sys.argv[0])
HERE = os.path.dirname(SCRIPT)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from pymol import cmd  # noqa: E402

from protfig_lib import core  # noqa: E402
from protfig_lib.core import ProtFigError  # noqa: E402
from protfig_lib.presets import PRESETS, Ctx, show_glycans  # noqa: E402

KIND_NAMES = {"small": "ligand", "peptide": "peptide", "protein": "protein"}
ONOFF = ("on", "off")


def build_parser():
    p = argparse.ArgumentParser(
        prog="protfig", usage="protfig <structure> --preset <style> [options]   |   protfig --list   |   "
                             "protfig --config figures.toml",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)
    g = p.add_argument_group("basic")
    g.add_argument("input", nargs="?", help="structure: .pdb, .cif, PDB code or .pse")
    g.add_argument("prefix", nargs="?", help=argparse.SUPPRESS)  # compatibility: explicit output prefix
    g.add_argument("-p", "--preset", default="overview", choices=sorted(PRESETS), metavar="STYLE",
                   help="figure style (see --list). Default: overview")
    g.add_argument("--out", default="figures", help="root output folder (default: figures)")
    g.add_argument("--name", help="figure name (default: input file name)")
    g.add_argument("--config", help=".toml file with one or more figures (see examples/figures.toml)")
    g.add_argument("--jobs", type=int, default=1, help="with --config: figures in parallel (default 1)")
    g.add_argument("--draft", action="store_true", help="quick draft: small image, no antialiasing")
    g.add_argument("--list", action="store_true", help="list the presets and which systems they work with")
    g.add_argument("-h", "--help", action="help", help="show this help")

    g = p.add_argument_group("what to show (selections in PyMOL syntax)")
    g.add_argument("--chains", help="chains to include and color order, e.g. B,C,A (default: all)")
    g.add_argument("--partner", "--ligand", dest="partner",
                   help='binding partner: ligand, peptide or protein, e.g. "resn ATP" or "chain B". '
                        "Default: automatic (largest organic molecule, otherwise the smallest chain)")
    g.add_argument("--partner_chain", "--ligand_chain", dest="partner_chain",
                   help="look for the automatic ligand only in this chain")
    g.add_argument("--site", help="site residues (default: those within --site_cutoff of the partner)")
    g.add_argument("--site_cutoff", type=float, default=4.0,
                   help="site/interface residues: within N Å of the partner (default 4)")
    g.add_argument("--pocket_cutoff", type=float,
                   help="binding_surface: pocket residues within N Å of the partner (default 8)")
    g.add_argument("--highlight", help="region to highlight (outline)")
    g.add_argument("--key", help="key residues as sticks (outline, property)")
    g.add_argument("--mutation", help="mutated residue, as translucent spheres (site_faded)")
    g.add_argument("--focus", help="selection to place at the top in the automatic view")
    g.add_argument("--show_partner", choices=ONOFF, help="binding_surface: show the partner (default off)")
    g.add_argument("--partner_sticks", choices=ONOFF,
                   help="peptide/protein partner: side chains as sticks (default on); off = cartoon only")
    g.add_argument("--ref", help="second structure (superpose; site_faded: its ligand as a ghost)")
    g.add_argument("--map", help="density map .ccp4/.map/.mrc (density)")
    g.add_argument("--map_level", type=float, default=1.0, help="mesh contour level in σ (default 1.0)")
    g.add_argument("--values", help="property: text file of 'chain resi value' lines (e.g. conservation)")
    g.add_argument("--lipids", help='membrane: lipid selection if not recognized, e.g. "resn POPX"')
    g.add_argument("--membrane_radius", type=float, help="membrane: bilayer radius around the protein in Å (default 40)")
    g.add_argument("--region", choices=("site", "all"), help="electrostatic/lipophilic: site or whole molecule")
    g.add_argument("--repr", choices=("cartoon", "surface"), default="cartoon", help="illustrative: cartoon or surface")
    g.add_argument("--glycans", choices=ONOFF, help="glycans (NAG, MAN...) as sticks (default off); cartoon and surface styles")
    g.add_argument("--waters", action="store_true", help="keep crystallographic waters (pastel shows them)")

    g = p.add_argument_group("colors (PyMOL name, '#rrggbb' or 'r g b' from 0 to 1)")
    g.add_argument("--palette", choices=sorted(core.PALETTES), metavar="PALETTE",
                   help="per-chain palette: " + ", ".join(sorted(core.PALETTES)))
    g.add_argument("--colors", help='colors for specific chains: "A:salmon;B:#88cc44;C:0.2 0.5 0.8"')
    g.add_argument("--partner_color", help="partner color (carbons / cartoon)")
    g.add_argument("--receptor_color", "--base_color", dest="receptor_color",
                   help="receptor color (outline, site_faded, binding_surface, interface_closeup)")
    g.add_argument("--site_color", help="site residue color (site_faded, binding_surface)")
    g.add_argument("--highlight_color", help="--highlight color (default orange)")
    g.add_argument("--key_color", help="--key color (default yellow)")
    g.add_argument("--hbond_color", help="H-bond color (default yellow)")
    g.add_argument("--mutation_color", help="color of the --mutation spheres")
    g.add_argument("--ref_color", help="color of the --ref structure in superpose")
    g.add_argument("--mesh_color", help="density mesh color")
    g.add_argument("--glycan_color", help="glycan carbon color (default blue)")
    g.add_argument("--lipid_color", help="lipid carbon color")
    g.add_argument("--partner_surface", choices=ONOFF, help="site: blue ligand surface (default on)")
    g.add_argument("--partner_surface_color", help="site: ligand surface color")
    g.add_argument("--bg", default="white", help="background: white (default), black, transparent or a color")
    g.add_argument("--transparent", dest="bg", action="store_const", const="transparent", help=argparse.SUPPRESS)
    g.add_argument("--outline", choices=ONOFF, help="black outline around shapes (default: depends on the style)")
    g.add_argument("--outline_color", help="outline color")
    g.add_argument("--surface_transparency", type=float, help="surface transparency (0-1)")
    g.add_argument("--cartoon_transparency", type=float, help="cartoon transparency (0-1)")

    g = p.add_argument_group("labels")
    g.add_argument("--labels", choices=ONOFF, help="residue and name labels (default: depends on the style)")
    g.add_argument("--no_labels", dest="labels", action="store_const", const="off", help=argparse.SUPPRESS)
    g.add_argument("--label_size", type=float, default=1.0, help="relative label size (default 1.0)")
    g.add_argument("--font", choices=sorted(core.FONTS), metavar="FONT",
                   help="font: " + ", ".join(sorted(core.FONTS)) + " (default: depends on the style)")
    g.add_argument("--label_color", help="a single color for all labels (default: the color of their chain)")
    g.add_argument("--label_style", choices=sorted(core.LABEL_FORMATS), default="one",
                   help="format: one=H353, three=HIS353, chain=A:H353, chain-three=A:HIS353")
    g.add_argument("--label_sel", help="label only these site residues")
    g.add_argument("--label_halo", choices=ONOFF, default="on", help="background-colored halo around the text")
    g.add_argument("--chain_names", help='chain names (overview): "B:CRBN,C:CK1α,A:DDB1"')
    g.add_argument("--partner_name", "--ligand_name", dest="partner_name", help="ligand name (overview)")
    g.add_argument("--distances", choices=ONOFF, help="site_faded: H-bond distances (default on)")

    g = p.add_argument_group("H-bonds")
    g.add_argument("--hbonds", choices=ONOFF, default="on", help="draw H-bonds (default on)")
    g.add_argument("--no_hbonds", dest="hbonds", action="store_const", const="off", help=argparse.SUPPRESS)
    g.add_argument("--hbond_cutoff", type=float, default=3.5, help="maximum N/O distance in Å (default 3.5)")
    g.add_argument("--hbond_mode", choices=("geometric", "distance"), default="geometric",
                   help="geometric = distance + angle (default); distance = every donor-acceptor pair within "
                        "--hbond_cutoff (catches Lys/Ser/Thr/Tyr H-bonds that the angle test can miss)")

    g = p.add_argument_group("color scales (electrostatic, lipophilic, property)")
    g.add_argument("--range", type=float, help="symmetric scale range (electrostatic: default ±10)")
    g.add_argument("--prop_label", help="property: color bar title")
    g.add_argument("--prop_range", type=float, nargs=2, metavar=("MIN", "MAX"), help="property: scale range")
    g.add_argument("--prop_ticks", help='property: text at the ends of the color bar, e.g. "Least,Most"')

    g = p.add_argument_group("camera")
    g.add_argument("--view", help="file with the view (18 numbers from get_view, e.g. <name>_view.txt)")
    g.add_argument("--session", help=".pse whose view and selections are applied to <structure>")
    g.add_argument("--auto_view", action="store_true", help="ignore the .pse view and compute a new one")
    g.add_argument("--roll", type=float, default=0.0, help="rotate the image in-plane (degrees)")
    g.add_argument("--turn", help='extra rotations: "y:30,x:-15" (y = horizontal, x = vertical)')
    g.add_argument("--flip", action="store_true", help="rotate 180° in-plane")
    g.add_argument("--clip", type=float, help="clip anything more than N Å in front of the partner")

    g = p.add_argument_group("image")
    g.add_argument("--width", type=int, help="width in px (default: depends on the style)")
    g.add_argument("--height", type=int, help="height in px (default: depends on the style)")
    g.add_argument("--dpi", type=int, default=300)
    g.add_argument("--format", default="png",
                   help="figure format(s), comma-separated: " + ", ".join(IMAGE_FORMATS) +
                        " (default png). jpg and pdf lack transparency: flattened onto white")
    g.add_argument("--no_png", action="store_true", help="only the .pse (fast, for orienting by hand)")
    g.add_argument("--set", action="append", default=[], metavar="SETTING=VALUE",
                   help="any PyMOL setting, applied last (repeatable)")
    return p


def parse(argv):
    p = build_parser()
    a = p.parse_args(argv)
    if a.list or a.config:
        return a
    if not a.input:
        p.error("missing input structure (or use --config / --list)")
    a.formats = [f.strip().lower().lstrip(".") for f in a.format.split(",") if f.strip()]
    a.formats = ["jpg" if f == "jpeg" else "tiff" if f == "tif" else f for f in a.formats]
    if not a.formats or any(f not in IMAGE_FORMATS for f in a.formats):
        p.error(f"--format {a.format}: valid formats are {', '.join(IMAGE_FORMATS)}")
    meta = PRESETS[a.preset]
    w, h = meta["size"]
    a.width, a.height = a.width or w, a.height or h
    if a.draft:
        a.width, a.height = int(a.width * 0.4), int(a.height * 0.4)
    if not a.prefix:
        name = a.name or os.path.splitext(os.path.basename(a.input))[0]
        a.prefix = os.path.join(a.out, a.preset, name, name)
    return a


def list_presets():
    print("PROTFIG presets (--preset). 'Works with' = type of binding partner accepted;"
          " all also work with the protein alone, except those marked with *.\n")
    for name, m in sorted(PRESETS.items()):
        kinds = ", ".join(KIND_NAMES[k] for k in m["kinds"]) or "protein(s) without partner"
        star = " *" if m["partner"] == "required" else ""
        print(f"{name}{star}\n    works with: {kinds}")
        print(textwrap.fill(m["doc"], width=96, initial_indent="    ", subsequent_indent="    ") + "\n")
    print("* requires a binding partner (automatic or --partner).")


IMAGE_FORMATS = ("png", "tiff", "jpg", "pdf")
LABEL_OBJECTS = ("LABELS", "CHAIN_LABELS", "DIST_LABELS")


def render(a, bg_rgb):
    """Final PNG. If there are labels, the scene and the labels are rendered separately and the
    labels are composited on top (when ray tracing, a label behind a translucent cartoon or a
    surface gets hidden). Also writes <name>_nolabels.png and <name>_labels.svg (labels as
    editable text, same size) for editing in Illustrator/Inkscape."""
    labels = [o for o in LABEL_OBJECTS if o in cmd.get_names("objects", enabled_only=1)]
    if not labels:
        cmd.png(f"{a.prefix}.png", width=a.width, height=a.height, dpi=a.dpi, ray=1)
        return
    from PIL import Image, ImageFilter

    enabled = cmd.get_names("objects", enabled_only=1)
    for o in labels:
        cmd.disable(o)
    cmd.png(f"{a.prefix}_nolabels.png", width=a.width, height=a.height, dpi=a.dpi, ray=1)
    for o in enabled:
        cmd.disable(o)
    for o in labels:
        cmd.enable(o)
    opaque = cmd.get_setting_int("ray_opaque_background")
    cmd.set("ray_opaque_background", 0)
    cmd.png(f"{a.prefix}_labels.png", width=a.width, height=a.height, dpi=a.dpi, ray=1)
    cmd.set("ray_opaque_background", opaque)
    for o in enabled:
        cmd.enable(o)
    base = Image.open(f"{a.prefix}_nolabels.png").convert("RGBA")
    top = Image.open(f"{a.prefix}_labels.png").convert("RGBA")
    if a.label_halo == "on":
        # background-colored halo around the text (readable over sticks and surfaces)
        r = max(3, int(round(a.width / 600)) * 2 + 1)
        halo_a = top.getchannel("A").filter(ImageFilter.MaxFilter(r)).point(lambda v: int(v * 0.85))
        halo = Image.new("RGBA", top.size, tuple(int(255 * x) for x in bg_rgb) + (0,))
        halo.putalpha(halo_a)
        base = Image.alpha_composite(base, halo)
    Image.alpha_composite(base, top).save(f"{a.prefix}.png", dpi=(a.dpi, a.dpi))
    os.remove(f"{a.prefix}_labels.png")
    halo_on = a.label_halo == "on"
    core.labels_svg(f"{a.prefix}_labels.svg", labels, a.width, a.height, a.dpi,
                    halo_rgb=bg_rgb if halo_on else None, halo_px=(r - 1) // 2 if halo_on else 0)


def convert_formats(a):
    """Copies <name>.png (and _nolabels.png) to the --format formats, at the same dpi; if png
    is not in the list, the .png files are deleted. TIFF uses LZW compression (lossless)."""
    fmts = a.formats
    from PIL import Image

    for png in (f"{a.prefix}.png", f"{a.prefix}_nolabels.png"):
        if not os.path.exists(png):
            continue
        im = Image.open(png)
        im.load()
        stem = png[:-4]
        for f in fmts:
            if f == "png":
                continue
            out = im
            if f in ("jpg", "pdf") and im.mode in ("RGBA", "LA", "P"):
                out = Image.new("RGB", im.size, "white")
                out.paste(im.convert("RGBA"), mask=im.convert("RGBA").getchannel("A"))
            if f == "tiff":
                out.save(f"{stem}.tif", dpi=(a.dpi, a.dpi), compression="tiff_lzw")
            elif f == "jpg":
                out.convert("RGB").save(f"{stem}.jpg", dpi=(a.dpi, a.dpi), quality=95, subsampling=0)
            elif f == "pdf":
                out.convert("RGB").save(f"{stem}.pdf", resolution=a.dpi)
        if "png" not in fmts:
            os.remove(png)


def write_run(a, argv, chains, partner_txt):
    """<name>_run.txt: where the figure came from and how to regenerate it."""
    files = sorted(os.path.basename(f) for f in glob.glob(f"{glob.escape(a.prefix)}*"))
    with open(f"{a.prefix}_run.txt", "w") as fh:
        fh.write(f"# PROTFIG {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
        fh.write(f"folder:    {os.getcwd()}\n")
        fh.write("command:   pymol -cq " + shlex.quote(SCRIPT) + " -- " + " ".join(shlex.quote(x) for x in argv) + "\n")
        fh.write(f"input:     {a.input}" + (f"  (session: {a.session})" if a.session else "") + "\n")
        fh.write(f"style:     {a.preset}\nchains:    {','.join(chains)}\npartner:   {partner_txt}\n")
        fh.write(f"image:     {a.width}x{a.height} px, {a.dpi} dpi" + (" (draft)" if a.draft else "") + "\n")
        fh.write("files:\n" + "".join(f"  {f}\n" for f in files + [os.path.basename(a.prefix) + "_run.txt"]))


def apply_background(a):
    """Background, plus a matching default label color (light text on a dark background)."""
    if a.bg == "transparent":
        cmd.set("bg_rgb", "white")
        cmd.set("ray_opaque_background", 0)
        return (1.0, 1.0, 1.0)
    rgb = core.rgb_of(a.bg)
    cmd.set("bg_rgb", list(rgb))
    cmd.set("ray_opaque_background", 1)
    if sum(rgb) / 3 < 0.4:
        core.LABEL_OPTS["dark_bg"] = True
    return rgb


def run_one(argv):
    a = parse(argv)
    if a.list:
        list_presets()
        return
    if a.config:
        run_config(a.config, a.jobs, extra=strip_opts(argv, ("--config", "--jobs")))
        return
    meta = PRESETS[a.preset]
    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)

    cmd.reinitialize()
    session_view, user = core.load_input(a)
    for k, v in core.BASE_SETTINGS.items():
        cmd.set(k, v)
    core.LABEL_OPTS.update(scale=a.label_size, font=a.font, color=a.label_color, style=a.label_style, dark_bg=False)
    bg_rgb = apply_background(a)
    core.clean(a, user, keep_lipids=a.preset == "membrane")
    chains = core.chain_order(a, user)
    lig, kind = core.find_partner(a, user, meta["partner"])
    if kind and meta["kinds"] and kind not in meta["kinds"]:
        print(f"WARNING: style '{a.preset}' is not designed for a {KIND_NAMES[kind]} partner; "
              f"works with: {', '.join(KIND_NAMES[k] for k in meta['kinds'])}.")

    core.ASPECT = a.width / a.height
    ctx = Ctx(a, lig, kind, chains, user, session_view)
    view = meta["fn"](ctx)
    if a.glycans == "on":
        if meta["glycans"]:
            show_glycans(ctx)
        else:
            print(f"WARNING: --glycans does not apply to style '{a.preset}' (only cartoon and surface styles: "
                  f"{', '.join(n for n, m in sorted(PRESETS.items()) if m['glycans'])}).")

    if a.outline == "on":
        if cmd.get_setting_int("ray_trace_mode") not in (1, 3):
            cmd.set("ray_trace_mode", 1)
    elif a.outline == "off":
        cmd.set("ray_trace_mode", 0)
    if a.outline_color:
        cmd.set("ray_trace_color", core.define("mf_outline_user", a.outline_color))
    if a.draft:
        cmd.set("antialias", 0)
    for s in a.set:
        k, v = s.split("=", 1)
        cmd.set(k.strip(), v.strip())
    cmd.set_view(view)
    with open(f"{a.prefix}_view.txt", "w") as fh:
        fh.write(", ".join(f"{v:.6f}" for v in view) + "\n")

    # editable selections in the .pse (reused when running again from the .pse)
    if lig:
        cmd.select("partner", f"{core.BASE} and ligand")
    if "SITE_STICKS" in cmd.get_names("objects"):
        cmd.select("site", f"{core.BASE} in SITE_STICKS")
    cmd.select("protfig_chains_" + "_".join(chains), "none")  # color order, for reusing the .pse
    for s in ("ligand", "lig_core"):
        cmd.delete(s) if s in cmd.get_names("selections") else None
    cmd.disable(core.BASE)
    cmd.deselect()
    cmd.save(f"{a.prefix}.pse")
    if not a.no_png:
        render(a, bg_rgb)
        convert_formats(a)
    partner_txt = "-"
    if lig:
        r = core.residues(f"{core.BASE} and partner")
        partner_txt = (f"{KIND_NAMES[kind]} {r[0][2]} {r[0][0]}:{r[0][1]}" if len(r) == 1
                       else f"{KIND_NAMES[kind]}, {len(r)} residues (chain {r[0][0]})")
    write_run(a, argv, chains, partner_txt)
    ext = {"tiff": "tif"}.get(a.formats[0], a.formats[0])
    print(f"[PROTFIG:{a.preset}] {a.prefix}.{ext}  chains {','.join(chains)}; partner: {partner_txt}")


def config_argvs(path):
    """Converts a .toml into a list of argv. Keys = option names (without --).

        [defaults]            # optional: applies to all figures
        out = "figures"
        [[figure]]
        input = "8G66.cif"
        preset = "site"
        chains = "B,C,A"
        labels = "on"
    A file without [[figure]] describes a single figure with top-level keys.
    """
    import tomllib

    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    defaults = data.get("defaults", {})
    figures = data.get("figure") or [{k: v for k, v in data.items() if k != "defaults"}]
    out = []
    for fig in figures:
        opts = {**defaults, **fig}
        if "input" not in opts:
            raise ProtFigError(f"in {path}: a figure has no 'input'")
        argv = [str(opts.pop("input"))]
        for k, v in opts.items():
            flag = "--" + k
            if v is True:
                argv.append(flag)
            elif v is False or v is None:
                continue
            elif isinstance(v, list):
                if k == "set":
                    for item in v:
                        argv += [flag, str(item)]
                else:
                    argv += [flag] + [str(x) for x in v]
            else:
                argv += [flag, str(v)]
        out.append(argv)
    return out


def strip_opts(argv, names):
    """argv without the options in names (or their values): the rest applies to every figure."""
    out, skip = [], False
    for x in argv:
        if skip:
            skip = False
        elif x in names:
            skip = True
        elif not any(x.startswith(n + "=") for n in names):
            out.append(x)
    return out


def run_config(path, jobs=1, extra=()):
    """Runs every figure in a .toml (paths relative to the .toml). One error does not stop the rest.

    jobs > 1: figures in parallel, each in its own PyMOL process with its log in the figure folder.
    Figures that use a .pse as input run last (they may depend on the others).
    extra: command-line options added to every figure (e.g. --draft).
    """
    here = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(path)))
    argvs = [argv + list(extra) for argv in config_argvs(os.path.basename(path))]
    failed = []
    try:
        if jobs <= 1:
            for argv in argvs:
                try:
                    run_one(argv)
                except (ProtFigError, SystemExit) as e:
                    failed.append(argv)
                    print(f"ERROR in {' '.join(argv[:3])}: {e}")
        else:
            first = [v for v in argvs if not v[0].lower().endswith(".pse")]
            later = [v for v in argvs if v[0].lower().endswith(".pse")]
            for wave in (first, later):
                failed += _run_parallel(wave, jobs)
    finally:
        os.chdir(here)
    print(f"{len(argvs) - len(failed)}/{len(argvs)} figures done" + (f"; {len(failed)} failed" if failed else ""))


def _run_parallel(argvs, jobs):
    import subprocess
    from concurrent.futures import ThreadPoolExecutor

    pymol = os.environ.get("PYMOL") or sys.executable.replace("python", "pymol")
    if not os.path.exists(pymol):
        import shutil
        pymol = shutil.which("pymol") or "pymol"

    def one(argv):
        prefix = parse(argv).prefix
        os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
        with open(prefix + ".log", "w") as log:
            r = subprocess.run([pymol, "-cq", SCRIPT, "--"] + argv, stdout=log, stderr=subprocess.STDOUT)
        lines = open(prefix + ".log").read().splitlines()
        ok = r.returncode == 0 and any(l.startswith("[PROTFIG:") for l in lines)
        msg = next((l for l in lines if l.startswith("[PROTFIG:")), None) if ok else \
            f"ERROR in {prefix}.png (see {prefix}.log): " + next((l for l in reversed(lines) if l.strip()), "")
        print(msg, flush=True)
        return None if ok else argv

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return [v for v in ex.map(one, argvs) if v]


def main(argv):
    """Entry point: user errors as 'ERROR: ...' (no traceback) and exit code 1."""
    try:
        run_one(argv)
    except ProtFigError as e:
        print(f"ERROR: {e}", file=sys.stderr, flush=True)
        sys.stdout.flush()
        os._exit(1)  # inside PyMOL, sys.exit prints a traceback and exits with code 0
    except SystemExit as e:  # argparse (--help, invalid option): same problem
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(e.code if isinstance(e.code, int) else 1)


if __name__ in ("__main__", "pymol"):
    main(sys.argv[1:])
