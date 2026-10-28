#!/bin/sh
# Downloads the public structures used by the examples (only missing ones) into inputs/.
cd "$(dirname "$0")"
mkdir -p inputs
get() {  # get <file> <url>
  [ -s "inputs/$1" ] && return
  echo "downloading $1"
  curl -fsSL -o "inputs/$1" "$2" || { echo "ERROR downloading $2" >&2; rm -f "inputs/$1"; }
}
for id in 8G66 8DZ0 8DZ1 1BRS 1YCR 3P0G 1U19; do get "$id.cif" "https://files.rcsb.org/download/$id.cif"; done
get 9P3Y.cif "https://files.rcsb.org/download/9P3Y.cif"
get 8dz0_2fofc.ccp4 "https://www.ebi.ac.uk/pdbe/entry-files/download/8dz0.ccp4"
# β2-adrenergic receptor in a DPPC bilayer (MemProtMD simulation) + nanobody from 3P0G
get 3P0G_membrane_dppc.pdb \
  "https://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/3p0g_default_dppc/files/structures/at.pdb"
[ -s inputs/3P0G_membrane_nanobody.pdb ] || ${PYMOL:-pymol} -cq ../tools/embed_in_membrane.py -- \
  inputs/3P0G_membrane_dppc.pdb inputs/3P0G.cif inputs/3P0G_membrane_nanobody.pdb --chains A,B
