"""
embed_in_membrane.py — places an experimental complex (PDB/cryo-EM) inside a lipid bilayer
taken from a simulation (MemProtMD, CHARMM-GUI...), so it can be used with `protfig --preset membrane`.

The simulated protein is used only as a guide: the experimental complex (with all its chains
and ligands) is superposed onto it and replaces the simulated protein. Water and ions are discarded.

  pymol -cq embed_in_membrane.py -- <membrane.pdb> <complex.cif> <output.pdb> [--chains R,A,B,G,N]

Example (GLP-1R–Gs with orforglipron in DPPC):
  curl -o 6xox_dppc.pdb https://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/6xox_default_dppc/files/structures/at.pdb
  pymol -cq embed_in_membrane.py -- 6xox_dppc.pdb 6XOX.cif 6XOX_in_membrane.pdb --chains R,A,B,G,N
"""
import argparse
import sys

from pymol import cmd

MD_SOLVENT = "SOL+TIP3+TIP4+TP3+T3P+SPC+WAT+HOH+NA+CL+K+SOD+CLA+POT"


def main(argv):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("membrane", help="simulated system with lipids (PDB)")
    p.add_argument("complex", help="experimental structure to insert (.pdb/.cif)")
    p.add_argument("output", help="output PDB: complex + lipids")
    p.add_argument("--chains", help="chains of the complex to keep (default: all)")
    a = p.parse_args(argv)

    cmd.load(a.membrane, "md")
    cmd.remove(f"md and (resn {MD_SOLVENT} or hydro)")
    cmd.load(a.complex, "xtal")
    cmd.remove("xtal and (hydro or solvent)")
    keep = f"xtal and chain {'+'.join(a.chains.split(','))}" if a.chains else "xtal"
    # structural superposition (independent of residue numbering) of the whole complex onto the
    # simulated protein; super ignores the parts that have no equivalent
    rms = cmd.super(f"({keep}) and polymer and name CA", "md and polymer and name CA")
    print(f"[embed] superposition RMSD: {rms[0]:.2f} Å over {rms[1]} CA")
    if rms[1] < 50:
        print("[embed] WARNING: few atoms superposed; does the simulation correspond to this protein?")
    # MemProtMD lipids are spread over chains B..Z, which clash with those of the complex:
    # move them all to a free chain
    used = set(cmd.get_chains(keep))
    free = next(c for c in "LMWXYZ9876" if c not in used)
    cmd.alter("md and not polymer", f"chain = '{free}'; segi = ''")
    cmd.create("system", f"({keep}) or (md and not polymer)")
    cmd.save(a.output, "system")
    print(f"[embed] {a.output}: complex ({keep}) + lipids from {a.membrane}")


if __name__ in ("__main__", "pymol"):
    main(sys.argv[1:])
