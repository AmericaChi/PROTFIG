"""
PROTFIG — publication-quality molecular figures from style presets (PyMOL, headless).

USAGE
  protfig <structure> --preset <style> [options]           (or: pymol -cq protfig.py -- ...)
  protfig --list                                           available presets
  protfig --config figures.toml                            one or more figures from a file
  protfig <structure> --preset site [options] --save_preset mySite    save options as a preset
  protfig <structure> --preset mySite                      reuse them (options given here win)

  <structure>: .pdb / .cif file, a PDB code (downloaded) or a .pse session.

EXAMPLES
  protfig 8G66.cif --preset site --chains B,C,A
  protfig 1BRS --preset site --style outline --chains A,D
  protfig my_protein.pdb --preset overview --style illustrative --palette pastel --bg transparent
  protfig figures/site/8G66/8G66.pse --preset site --style faded   (reuse view and selections)

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
import json
import os
import re
import shlex
import sys
import textwrap
import time

# PyMOL defines __script__ when running `pymol -cq protfig.py`; plain python defines __file__
SCRIPT = os.path.abspath(globals().get("__script__") or globals().get("__file__") or sys.argv[0])
HERE = os.path.dirname(SCRIPT)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from pymol import cmd  # noqa: E402

from protfig_lib import core  # noqa: E402
from protfig_lib.core import ProtFigError  # noqa: E402
from protfig_lib.presets import (PRESETS, Ctx, mutation_labels, set_sidechain_mode, show_glycans, show_het,
                                 show_metals, show_mutation)  # noqa: E402

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
    g.add_argument("-p", "--preset", default="overview", choices=sorted(list(PRESETS) + list(ALIASES)), metavar="STYLE",
                   help="figure style (see --list), or a preset saved with --save_preset. Default: overview")
    g.add_argument("--save_preset", metavar="NAME",
                   help="save the options of this command (all but the structure and output names) as a "
                        "preset NAME, used as --preset NAME (folder: $PROTFIG_PRESETS or ~/.protfig/presets)")
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
                   help="footprint: receptor residues within N Å of the partner (default 4.5)")
    g.add_argument("--highlight", help="region to highlight (overview, site); \"interface\" shades both sides of a protein-protein/peptide interface")
    g.add_argument("--key", help="key residues as sticks (overview, property)")
    g.add_argument("--mutation", help="mutated residue(s), as sticks in translucent spheres with a label (any style)")
    g.add_argument("--mutation_original", choices=ONOFF,
                   help="with --mutation and --ref (the wild type): also the original residue, in translucent gray")
    g.add_argument("--covalent_cutoff", type=float,
                   help="covalent: max partner-receptor distance in Å for a bond missing from the file (default 2.0)")
    g.add_argument("--style", choices=("faded", "outline", "comic", "illustrative", "pastel", "assembly"),
                   help="site: faded (gray protein, colored sticks), outline (pastel with black outline); "
                        "overview: outline, comic, illustrative, pastel, assembly")
    g.add_argument("--scope", choices=("whole", "site"),
                   help="covalent: whole protein (default); site = the close-up, now --preset site")
    g.add_argument("--interface_sticks", choices=ONOFF,
                   help=argparse.SUPPRESS)  # former covalent --scope site option, accepted for old .toml files
    g.add_argument("--focus", help="selection to place at the top in the automatic view")
    g.add_argument("--show_partner", choices=ONOFF, help="footprint: show the partner (default off)")
    g.add_argument("--partner_sticks", choices=ONOFF,
                   help="peptide/protein partner: side chains as sticks (default on); off = cartoon only")
    g.add_argument("--ref", help="second structure: file or PDB code (superpose; site: its ligand as a ghost)")
    g.add_argument("--ref_chains", help="superpose: chains of --ref to align and draw (default: all)")
    g.add_argument("--map", help="density map .ccp4/.map/.mrc (density)")
    g.add_argument("--map_level", type=float, default=1.0, help="mesh contour level in σ (default 1.0)")
    g.add_argument("--values", help="property: text file of 'chain resi value' lines (e.g. conservation)")
    g.add_argument("--lipids", help='membrane: lipid selection if not recognized, e.g. "resn POPX"')
    g.add_argument("--membrane_radius", type=float, help="membrane: bilayer radius around the protein in Å (default 40)")
    g.add_argument("--region", choices=("site", "all"), help="electrostatic/lipophilic: site or whole molecule")
    g.add_argument("--repr", choices=("cartoon", "surface", "spheres"),
                   help="overview --style illustrative: cartoon (default) or surface; membrane: cartoon (default), spheres or surface")
    g.add_argument("--glycans", choices=ONOFF, help="glycans (NAG, MAN...) as sticks (default off); cartoon and surface styles")
    g.add_argument("--waters", action="store_true",
                   help="keep waters (overview --style pastel: crystallographic waters as spheres; membrane: MD water as a translucent surface)")
    g.add_argument("--water_color", help="membrane --waters: water surface color (default pale blue)")
    g.add_argument("--water_layer", help="membrane --waters: water within N Å of the bilayer (default 10), or 'all' for the whole box")
    g.add_argument("--keep_het", action="store_true",
                   help="keep and show crystallization additives and ions (GOL, SO4, NA, CL...), removed by default")
    g.add_argument("--nonstandard_sticks", choices=ONOFF,
                   help="sticks also show side chains of non-standard residues (MSE, SEP, staples...); default off")
    g.add_argument("--metals", choices=ONOFF,
                   help="metal ions (ZN, MG, CA...) as spheres (default: on in site, off elsewhere)")

    g = p.add_argument_group("colors (PyMOL name, '#rrggbb' or 'r g b' from 0 to 1)")
    g.add_argument("--palette", choices=sorted(core.PALETTES), metavar="PALETTE",
                   help="per-chain palette: " + ", ".join(sorted(core.PALETTES)))
    g.add_argument("--colors", help='colors for specific chains: "A:salmon;B:#88cc44;C:0.2 0.5 0.8"')
    g.add_argument("--partner_color", help="partner color (carbons / cartoon)")
    g.add_argument("--receptor_color", "--base_color", dest="receptor_color",
                   help="receptor color (overview --style outline, site --style faded, footprint)")
    g.add_argument("--site_color", help="footprint color (footprint)")
    g.add_argument("--highlight_color", help="--highlight color (default orange)")
    g.add_argument("--partner_highlight_color", help="partner-side interface color with --highlight interface (default: darker partner color)")
    g.add_argument("--key_color", help="--key color (default yellow)")
    g.add_argument("--hbond_color", help="H-bond color (default yellow)")
    g.add_argument("--mutation_color", help="color of the --mutation spheres")
    g.add_argument("--ref_color", help="color of the --ref structure in superpose")
    g.add_argument("--mesh_color", help="density mesh color")
    g.add_argument("--glycan_color", help="glycan carbon color (default blue)")
    g.add_argument("--het_color", help="--keep_het carbon color (default gray)")
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
    g.add_argument("--distances", choices=ONOFF, help="site: H-bond lengths in Å (default off)")

    g = p.add_argument_group("H-bonds")
    g.add_argument("--hbonds", choices=ONOFF, help="draw H-bonds (default on; off in overview --style comic)")
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
    g.add_argument("--width", type=int, help="width in px (default: 180 mm at --dpi, i.e. 2126 px at 300 dpi)")
    g.add_argument("--height", type=int, help="height in px (default: 3/4 of the width)")
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
    argv, user_preset = apply_user_preset(p, argv)
    a = p.parse_args(argv)
    a.user_preset, a.saved_only = user_preset, False
    if a.save_preset:
        save_user_preset(a)
        a.saved_only = not (a.input or a.list or a.config)  # only saving: no figure
    if a.list or a.config or a.saved_only:
        return a
    if not a.input:
        p.error("missing input structure (or use --config / --list)")
    if a.preset == "covalent" and a.scope == "site":  # the close-up of a covalent link is now site
        print("NOTE: 'covalent --scope site' is now '--preset site' (it detects covalent links)")
        a.preset = "site"
    if a.preset in ALIASES:  # former style -> its replacement, with the options that reproduce it
        new, opts = ALIASES[a.preset]
        print(f"NOTE: '{a.preset}' is now '" + " ".join([f"--preset {new}"] + [f"--{k} {v}" for k, v in opts.items()]) + "'")
        a.preset = new
        for k, v in opts.items():
            if getattr(a, k, None) is None:
                setattr(a, k, v)
    styles = {"site": ("faded", "outline"), "overview": ("outline", "comic", "illustrative", "pastel", "assembly")}
    if a.style and a.style not in styles.get(a.preset, ()):
        ok = [p for p, st in styles.items() if a.style in st]
        print(f"WARNING: --style {a.style} applies to {', '.join(ok)} only")
    a.formats = [f.strip().lower().lstrip(".") for f in a.format.split(",") if f.strip()]
    a.formats = ["jpg" if f == "jpeg" else "tiff" if f == "tif" else f for f in a.formats]
    if not a.formats or any(f not in IMAGE_FORMATS for f in a.formats):
        p.error(f"--format {a.format}: valid formats are {', '.join(IMAGE_FORMATS)}")
    meta = PRESETS[a.preset]
    # standard image: 180 mm wide (double column in Nature, Science, PNAS), 4:3, at --dpi
    if not a.width and not a.height:
        a.width = int(round(STANDARD_WIDTH_MM / 25.4 * a.dpi))
    a.width = a.width or int(round(a.height * 4 / 3))
    a.height = a.height or int(round(a.width * 3 / 4))
    if a.draft:
        a.width, a.height = int(a.width * 0.4), int(a.height * 0.4)
    if not a.prefix:
        name = a.name or os.path.splitext(os.path.basename(a.input))[0]
        a.prefix = os.path.join(a.out, a.user_preset or a.preset, name, name)
    return a


# --- saved presets (--save_preset NAME, then --preset NAME): a preset plus fixed options, in a .toml
USER_PRESETS = os.path.expanduser(os.environ.get("PROTFIG_PRESETS") or "~/.protfig/presets")
NOT_SAVED = ("input", "prefix", "out", "name", "config", "jobs", "list", "save_preset", "draft", "no_png")


def user_preset_path(name):
    return os.path.join(USER_PRESETS, name + ".toml")


def is_user_preset(name):
    return name not in PRESETS and name not in ALIASES and os.path.isfile(user_preset_path(name))


def read_user_preset(name, seen=()):
    """Options of a saved preset {option: value}, checked as if given on the command line.
    Its base preset may be another saved preset (options here win)."""
    import tomllib

    path = user_preset_path(name)
    if name in seen:
        raise ProtFigError(f"saved preset '{name}' uses itself as base ({path})")
    try:
        with open(path, "rb") as fh:
            opts = tomllib.load(fh)
    except tomllib.TOMLDecodeError as e:
        raise ProtFigError(f"saved preset '{name}': {path} is not valid TOML ({e})")
    base = str(opts.get("preset", "overview"))
    if is_user_preset(base):
        opts = {**read_user_preset(base, seen + (name,)), **{k: v for k, v in opts.items() if k != "preset"}}
    bad = [k for k in opts if k in NOT_SAVED]
    if bad:
        raise ProtFigError(f"saved preset '{name}' ({path}): '{bad[0]}' cannot be part of a preset")
    try:  # same checks as on the command line (option names, choices, numbers)
        ns = build_parser().parse_args(["_"] + opts_to_argv(opts))
    except SystemExit:
        raise ProtFigError(f"saved preset '{name}': invalid option in {path} (see the message above); "
                           "edit or delete the file")
    return {k: getattr(ns, k) for k in opts}


def apply_user_preset(p, argv):
    """--preset NAME of a saved preset: its options become the defaults (so the command line wins)
    and its base preset the --preset. Returns (argv without --preset NAME, NAME or None)."""
    for i, x in enumerate(argv):
        if x in ("-p", "--preset") and i + 1 < len(argv):
            name, span = argv[i + 1], 2
        elif x.startswith("--preset="):
            name, span = x.split("=", 1)[1], 1
        else:
            continue
        if not is_user_preset(name):
            return argv, None  # a built-in preset, or an unknown name (argparse reports it)
        p.set_defaults(**read_user_preset(name))
        return argv[:i] + argv[i + span:], name
    return argv, None


def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    return json.dumps(str(v), ensure_ascii=False)  # a JSON string is a valid TOML basic string


def save_user_preset(a):
    """Writes the options that differ from the defaults (plus --preset) as <NAME>.toml."""
    name = a.save_preset
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        raise ProtFigError(f"--save_preset {name}: use only letters, digits, '_' and '-'")
    if name in PRESETS or name in ALIASES:
        raise ProtFigError(f"--save_preset {name}: '{name}' is a built-in preset; choose another name")
    p = build_parser()
    defaults = p.parse_args(["_"])
    flags = {o[2:] for act in p._actions for o in act.option_strings if o.startswith("--")}
    opts = {"preset": a.preset}
    for k, v in vars(a).items():
        if k in NOT_SAVED or k == "preset" or k not in flags or not hasattr(defaults, k):
            continue
        if v != getattr(defaults, k) and v is not None and v is not False and v != []:
            opts[k] = list(v) if isinstance(v, tuple) else v
    path = user_preset_path(name)
    existed = os.path.exists(path)
    os.makedirs(USER_PRESETS, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"# PROTFIG saved preset '{name}' ({datetime.datetime.now():%Y-%m-%d %H:%M}).\n"
                 f"# Use: protfig <structure> --preset {name} [options]   (options given there win)\n"
                 "# Keys are the command-line options without '--'; edit or delete freely.\n")
        fh.write("".join(f"{k} = {toml_value(v)}\n" for k, v in opts.items()))
    print(f"{'Replaced' if existed else 'Saved'} preset '{name}': {path}\n    {cli_text(opts)}")


def cli_text(opts):
    return " ".join(shlex.quote(x) for x in opts_to_argv(opts))


def user_presets():
    return sorted(os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(USER_PRESETS, "*.toml")))


def list_presets():
    print("PROTFIG presets (--preset). 'Works with' = type of binding partner accepted;"
          " all also work with the protein alone, except those marked with *.\n")
    for name, m in sorted(PRESETS.items()):
        kinds = ", ".join(KIND_NAMES[k] for k in m["kinds"]) or "protein(s) without partner"
        star = " *" if m["partner"] == "required" else ""
        print(f"{name}{star}\n    works with: {kinds}")
        print(textwrap.fill(m["doc"], width=96, initial_indent="    ", subsequent_indent="    ") + "\n")
    print("* requires a binding partner (automatic or --partner).")
    names = user_presets()
    if names:
        import tomllib
        print(f"\nSaved presets ({USER_PRESETS}):")
        for n in names:
            try:
                with open(user_preset_path(n), "rb") as fh:
                    opts = tomllib.load(fh)
                print(f"{n}\n    {cli_text(opts)}")
            except (OSError, tomllib.TOMLDecodeError) as e:
                print(f"{n}\n    (cannot read: {e})")


# former styles, now options of another one: name -> (style, option defaults)
ALIASES = {"site_faded": ("site", {"style": "faded"}),
           "interface_closeup": ("site", {"style": "outline", "palette": "candy", "labels": "off"}),
           "binding_surface": ("footprint", {}),
           "outline": ("overview", {"style": "outline"}),
           "illustrative": ("overview", {"style": "illustrative"}),
           "pastel": ("overview", {"style": "pastel"}),
           "oligomer": ("overview", {"style": "assembly"})}
STANDARD_WIDTH_MM = 180  # default image width; height 3/4 of it
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


def write_run(a, argv, chains, partner_txt, times=None):
    """<name>_run.txt: where the figure came from and how to regenerate it."""
    files = sorted(os.path.basename(f) for f in glob.glob(f"{glob.escape(a.prefix)}*"))
    with open(f"{a.prefix}_run.txt", "w") as fh:
        fh.write(f"# PROTFIG {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
        fh.write(f"folder:    {os.getcwd()}\n")
        fh.write("command:   pymol -cq " + shlex.quote(SCRIPT) + " -- " + " ".join(shlex.quote(x) for x in argv) + "\n")
        fh.write(f"input:     {a.input}" + (f"  (session: {a.session})" if a.session else "") + "\n")
        if a.user_preset:
            fh.write(f"saved preset: {a.user_preset} ({user_preset_path(a.user_preset)})\n")
        fh.write(f"style:     {a.preset}\nchains:    {','.join(chains)}\npartner:   {partner_txt}\n")
        fh.write(f"image:     {a.width}x{a.height} px, {a.dpi} dpi" + (" (draft)" if a.draft else "") + "\n")
        if times:  # wall clock inside PROTFIG; PyMOL start-up is not included
            fh.write(f"time:      setup {times[0]:.1f} s (load, analysis, view), render {times[1]:.1f} s (ray tracing, files)\n")
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


def is_closeup(preset):
    return preset in ("site", "pocket", "density")


def run_one(argv):
    a = parse(argv)
    if a.saved_only:
        return
    if a.list:
        list_presets()
        return
    if a.config:
        run_config(a.config, a.jobs, extra=strip_opts(argv, ("--config", "--jobs")))
        return
    meta = PRESETS[a.preset]
    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)
    t_start = time.perf_counter()

    cmd.reinitialize()
    session_view, user = core.load_input(a)
    for k, v in core.BASE_SETTINGS.items():
        cmd.set(k, v)
    core.LABEL_OPTS.update(scale=a.label_size, font=a.font, color=a.label_color, style=a.label_style, dark_bg=False)
    bg_rgb = apply_background(a)
    core.clean(a, user, keep_lipids=a.preset == "membrane")
    chains = core.chain_order(a, user)
    mode = meta["partner"]
    if ((a.highlight or "").strip().lower() == "interface" or a.style == "comic") and mode == "optional":
        mode = "required"  # --highlight interface, --style comic: also look for a partner chain
    lig, kind = core.find_partner(a, user, mode)
    if kind and meta["kinds"] and kind not in meta["kinds"]:
        print(f"WARNING: style '{a.preset}' is not designed for a {KIND_NAMES[kind]} partner; "
              f"works with: {', '.join(KIND_NAMES[k] for k in meta['kinds'])}.")

    core.ASPECT = a.width / a.height
    core.IMAGE.update(height=a.height, dpi=a.dpi, scale=0.4 if a.draft else 1.0)  # draft: same look, smaller
    set_sidechain_mode(a.nonstandard_sticks == "on")
    ctx = Ctx(a, lig, kind, chains, user, session_view)
    view = meta["fn"](ctx)
    if a.mutation:  # any style (site draws it itself, before framing)
        show_mutation(ctx)
        mutation_labels(ctx, view)
    if a.glycans == "on":
        if meta["glycans"]:
            show_glycans(ctx)
        else:
            print(f"WARNING: --glycans does not apply to style '{a.preset}' (only cartoon and surface styles: "
                  f"{', '.join(n for n, m in sorted(PRESETS.items()) if m['glycans'])}).")
    # close-ups (site...) keep small spheres; whole-molecule views need them larger to be seen
    scale = 0.35 if is_closeup(a.preset) else 1.2
    if a.metals == "on" and "METAL_IONS" not in cmd.get_names("objects"):
        show_metals(near=None, scale=scale)
    if a.keep_het:
        show_het(ctx, scale=scale)

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
    t_setup = time.perf_counter()
    if not a.no_png:
        render(a, bg_rgb)
        convert_formats(a)
    times = (t_setup - t_start, time.perf_counter() - t_setup)
    partner_txt = "-"
    if lig:
        r = core.residues(f"{core.BASE} and partner")
        partner_txt = (f"{KIND_NAMES[kind]} {r[0][2]} {r[0][0]}:{r[0][1]}" if len(r) == 1
                       else f"{KIND_NAMES[kind]}, {len(r)} residues (chain {r[0][0]})")
    chains = ctx.chains  # the preset may leave out chains (e.g. crystal copies away from the partner)
    write_run(a, argv, chains, partner_txt, times)
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
        out.append([str(opts.pop("input"))] + opts_to_argv(opts))
    return out


def opts_to_argv(opts):
    """{option: value} -> command-line arguments (true -> flag; lists -> values, or repeated --set)."""
    argv = []
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
    return argv


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
