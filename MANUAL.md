# PROTFIG manual

1. [Basic use](#1-basic-use)
2. [The binding partner and the site](#2-the-binding-partner-and-the-site)
3. [Options](#3-options)
4. [Adjusting a figure in PyMOL](#4-adjusting-a-figure-in-pymol)
5. [Many figures from one file](#5-many-figures-from-one-file)
6. [Outputs](#6-outputs)
7. [Membrane proteins](#7-membrane-proteins)
8. [Notes and limitations](#8-notes-and-limitations)
9. [Adding a preset](#9-adding-a-preset)

## 1. Basic use

```sh
protfig <structure> --preset <preset> [options]
```

- `<structure>` is a `.pdb` or `.cif` file, a PDB code (it is downloaded) or a PyMOL `.pse` session.
- `protfig --list` describes every preset; `protfig -h` lists every option.
- `--draft` renders a small preview in seconds.
- `protfig` is a shortcut for `pymol -cq protfig.py -- <options>`. If PyMOL is not on your
  `PATH`, run `PYMOL=/path/to/pymol ./protfig ...`.

PROTFIG automatically removes hydrogens, waters, crystallization additives (glycerol, PEG,
sulfate...) and simulation solvent and ions. Keep crystallographic waters with `--waters`.

`--chains B,C,A` keeps only those chains and sets the color order: the first chain gets the first
color of the palette. Chain IDs are the author chain IDs (`auth` on the PDB website), which can
differ from the label IDs some pages show (in 7NRE, label chain B is author chain C). If a chain is
not found, PROTFIG lists the chains present and suggests the right ID.

## 2. The binding partner and the site

Site and interface presets are built around a **binding partner**: a ligand, a peptide or a
protein. If you don't specify one, PROTFIG chooses it this way:

1. The largest organic non-polymer residue, counted in heavy atoms. Crystallization additives,
   lipids and glycans are skipped.
2. If there is none, and the preset needs a partner, the smallest chain is used, provided the
   structure has at least two chains.

A polymer partner with ≤ 50 residues is treated as a **peptide**; a larger one is treated as a
**protein**. PROTFIG prints its choice (for example `partner: ligand YOT B:502`) and writes it to
`_run.txt`. To choose the partner yourself:

```sh
protfig 1ATP --preset site --partner "resn ATP"
protfig complex.pdb --preset site --partner "chain C"
protfig complex.pdb --preset site --partner_chain B      # automatic search, only in chain B
```

**Site residues** are the receptor residues with any heavy atom within `--site_cutoff` (default
4 Å) of the partner. Residues hydrogen-bonded to the partner are always included. For protein
partners, the same cutoff selects the partner side of the interface; peptides are shown in full.
`binding_surface` uses a wider cutoff, `--pocket_cutoff` (default 8 Å), so that the footprint is a
continuous patch. You can also give the site by hand with `--site "resi 41+145+166"`.

## 3. Options

Colors accept PyMOL names (`salmon`), `#rrggbb` or `"r g b"` (0-1). Defaults depend on the preset.

**What to show**

| Option | Effect |
|---|---|
| `--chains B,C,A` | chains to keep and their color order |
| `--partner`, `--partner_chain` | choose the binding partner |
| `--site`, `--site_cutoff`, `--pocket_cutoff` | site residues |
| `--highlight "resi 1-140"`, `--key "resi 41+145"` | highlighted region and key residues (`outline`); key residues are added to those forming H-bonds |
| `--highlight interface` | shade the interface on both sides of a protein-protein or peptide complex, with H-bond and `--key` residues of both sides in sticks (`outline`) |
| `--mutation "resi 49"`, `--ref other.cif` | mutated residue; second structure (`site_faded`, `superpose`) |
| `--map map.ccp4`, `--map_level 1.0` | electron density (`density`) |
| `--values file.txt` | per-residue values, lines `chain resi value` (`property`) |
| `--prop_label`, `--prop_range`, `--prop_ticks` | color bar title, scale range and end labels (`property`; `--prop_label` also in `plddt`) |
| `--region site\|all` | site or whole molecule (`electrostatic`, `lipophilic`) |
| `--glycans on` | glycans as sticks (cartoon and surface presets) |
| `--show_partner on` | show the partner in `binding_surface` |
| `--partner_sticks off` | peptide or protein partner as cartoon only, without side-chain sticks |
| `--waters` | keep crystallographic waters |

**Colors**

| Option | Effect |
|---|---|
| `--palette` | `vivid`, `chimerax`, `muted`, `pastel`, `illustrative`, `grays`, `membrane` |
| `--colors "A:salmon;B:#88cc44"` | color of specific chains |
| `--partner_color`, `--receptor_color`, `--site_color`, `--hbond_color` | main elements |
| `--highlight_color`, `--partner_highlight_color`, `--key_color`, `--mutation_color`, `--glycan_color`, `--lipid_color`, `--mesh_color`, `--ref_color` | preset-specific elements |
| `--bg white\|black\|transparent\|<color>` | background (default white) |
| `--outline on\|off`, `--outline_color` | black outline of the shapes |
| `--surface_transparency`, `--cartoon_transparency` | 0 to 1 |

**Labels and hydrogen bonds**

| Option | Effect |
|---|---|
| `--labels on\|off` | show or hide labels |
| `--label_size 1.3` | relative size (1.0 = default) |
| `--font sans-bold` | `sans`, `sans-bold`, `serif`, `serif-bold`, `mono`, ... |
| `--label_color black` | one color for all labels (default: the chain color) |
| `--label_style one\|three\|chain\|chain-three` | `H353`, `HIS353`, `A:H353`, `A:HIS353` |
| `--label_sel "resi 41+145"` | label only these residues |
| `--chain_names "A:Kinase,B:Inhibitor"`, `--partner_name` | names drawn in `overview` |
| `--hbonds on\|off`, `--hbond_cutoff 3.5`, `--hbond_mode geometric\|distance` | hydrogen bonds |
| `--distances on\|off` | distances in Å (`site_faded`) |

**Camera and image**

| Option | Effect |
|---|---|
| `--roll 30`, `--turn "y:30,x:-15"`, `--flip` | rotate the automatic view |
| `--clip 5` | remove what lies more than 5 Å in front of the partner |
| `--view file.txt`, `--session view.pse`, `--auto_view` | reuse a saved view, or ignore the one in a `.pse` |
| `--width`, `--height`, `--dpi` | image size (default: depends on the preset, 300 dpi) |
| `--format png,tiff,jpg,pdf` | one or several image formats (default `png`) |
| `--out folder`, `--name name` | output location (default `figures/<preset>/<input name>/`) |
| `--set ray_trace_mode=1` | any PyMOL setting (repeatable) |

## 4. Adjusting a figure in PyMOL

Every figure saves its PyMOL session, `<name>.pse`. It contains the cleaned structure and, when
they apply, the selections `partner` and `site`.

1. Open the `.pse` in PyMOL. Rotate or zoom, and edit `partner` or `site` if you wish. You can also
   create selections named `highlight`, `key`, `mutation` or `focus` (this part goes on top in
   automatic views).
2. Save the session.
3. Run PROTFIG with the session as input, with the same or another preset:

```sh
protfig figures/site/8G66/8G66.pse --preset site          # same figure, your view
protfig figures/site/8G66/8G66.pse --preset site_faded    # another preset, same view
```

The saved view, selections and chain colors are kept. Your own `.pse` made from scratch works too.

## 5. Many figures from one file

Each `[[figure]]` block is one figure. Keys are the command-line options without `--`, and paths
are relative to the file.

```toml
[defaults]            # optional, applies to all figures
out = "figures"

[[figure]]
input = "inputs/8G66.cif"
preset = "site"
chains = "B,C,A"

[[figure]]
input = "inputs/1BRS.cif"
preset = "interface_closeup"
chains = "A,D"
```

```sh
protfig --config figures.toml --jobs 4      # 4 figures in parallel
```

[`examples/figures.toml`](examples/figures.toml) has one example of every preset; run
`examples/run_examples.sh` to download the inputs and render them all.

## 6. Outputs

```
figures/<preset>/<name>/
├── <name>.png            final figure
├── <name>.pse            PyMOL session
├── <name>_nolabels.png   the figure without labels   } only when the
├── <name>_labels.svg     the labels as editable text } figure has labels
├── <name>_site.txt       site residues
├── <name>_hbonds.txt     hydrogen bonds (the same ones that are drawn)
├── <name>_colorbar.svg   color bar (electrostatic, lipophilic, property)
├── <name>_view.txt       camera view, reusable with --view
└── <name>_run.txt        command, date, partner, chains and file list
```

`_labels.svg` has the same physical size as the image. Place it over `_nolabels.png` in
Illustrator or Inkscape to move or edit labels as text. `--format tiff,pdf` also saves those
formats (TIFF is lossless). JPG and PDF have no transparency, so with `--bg transparent` they come
out on white.

## 7. Membrane proteins

The `membrane` preset needs a system with explicit lipids:

- **A simulation frame** (GROMACS, CHARMM-GUI, MemProtMD...) can be used directly:
  `protfig frame.pdb --preset membrane`. Water and ions are removed automatically. Export the
  frame with whole molecules (for example `gmx trjconv -pbc mol -center`), and name unusual lipids
  with `--lipids "resn XXXX"`.
- **An experimental structure** without a membrane can be placed into a simulated bilayer with
  `tools/embed_in_membrane.py`, which superposes it onto the simulated protein:

```sh
curl -o 6xox_dppc.pdb https://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/6xox_default_dppc/files/structures/at.pdb
pymol -cq tools/embed_in_membrane.py -- 6xox_dppc.pdb 6XOX.cif 6XOX_in_membrane.pdb --chains R,A,B,G,N
protfig 6XOX_in_membrane.pdb --preset membrane --partner "chain A+B+G+N"
```

The membrane protein is drawn as spheres; a partner given with `--partner` is drawn as cartoon.

## 8. Notes and limitations

- **Hydrogen bonds** are N/O donor-acceptor pairs within 3.5 Å that pass PyMOL's distance and
  angle test, estimated without hydrogens. The test can miss bonds from groups whose hydrogens
  rotate (Lys, Ser, Thr, Tyr). `--hbond_mode distance` keeps every donor-acceptor pair within the
  cutoff instead. Check `_hbonds.txt` before reporting them.
- **Electrostatics**: a Coulomb potential with a distance-dependent dielectric (ε = 4r, r in Å) and
  formal charges (Asp, Glu, Lys, Arg and termini at neutral pH, His neutral, metal ions, and charges
  annotated in the file). It is a fast approximation; use Poisson-Boltzmann (e.g. APBS) for
  quantitative work.
- **pLDDT**: the `plddt` preset reads AlphaFold confidence from the B-factor column (AlphaFold 2
  and 3 models, `.pdb` or `.cif`). On an experimental structure it would color by B-factor, so use
  it only for predicted models.
- **Lipophilicity**: molecular lipophilicity potential with simplified atomic contributions.
- **Density**: PDBe 2Fo-Fc maps can be downloaded from
  `https://www.ebi.ac.uk/pdbe/entry-files/download/<pdb id>.ccp4`.
- PROTFIG makes **individual panels**. Assemble multi-panel figures in a graphics editor.

## 9. Adding a preset

Presets live in [`protfig_lib/presets.py`](protfig_lib/presets.py). Each one is a function that
builds PyMOL objects and returns the camera view:

```python
@preset(partner="optional", kinds=("small", "peptide", "protein"), size=(1600, 1600))
def my_style(ctx):
    """Description shown by protfig --list."""
    ctx.chain_colors("pastel")                          # honors --palette and --colors
    cmd.create("CARTOON", ctx.receptor)
    cmd.show("cartoon", "CARTOON")
    color = ctx.rgb("receptor_color", (0.3, 0.5, 0.8))  # user option, or this default
    ...
    return ctx.view("overview")                         # or ctx.view("site", ...)
```

Read every visual choice through `ctx.opt` or `ctx.rgb`, so that users can change it from the
command line.
