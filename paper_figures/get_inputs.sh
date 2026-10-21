#!/bin/sh
# Downloads the manuscript inputs (PDB + MemProtMD simulation) and builds the membrane system.
cd "$(dirname "$0")"
mkdir -p inputs
for id in 6XOX 7KI0 8YW3 9SV3 9DTY 9ZTJ 9JN0; do
  [ -s inputs/$id.cif ] || curl -fsSL -o inputs/$id.cif https://files.rcsb.org/download/$id.cif
done
[ -s inputs/6XOX_membrane_dppc.pdb ] || curl -fsSL -o inputs/6XOX_membrane_dppc.pdb \
  https://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/6xox_default_dppc/files/structures/at.pdb
[ -s inputs/6XOX_in_membrane.pdb ] || ${PYMOL:-pymol} -cq ../tools/embed_in_membrane.py -- \
  inputs/6XOX_membrane_dppc.pdb inputs/6XOX.cif inputs/6XOX_in_membrane.pdb --chains R,A,B,G,N
