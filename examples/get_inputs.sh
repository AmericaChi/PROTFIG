#!/bin/sh
# Downloads the public structures used by the examples (only missing ones) into inputs/.
cd "$(dirname "$0")"
mkdir -p inputs
get() {  # get <file> <url>
  [ -s "inputs/$1" ] && return
  echo "downloading $1"
  curl -fsSL -o "inputs/$1" "$2" || { echo "ERROR downloading $2" >&2; rm -f "inputs/$1"; }
}
for id in 8G66 8DZ0 8DZ1 1BRS 1YCR 1U19 1EMA 4HT6 6M0J; do get "$id.cif" "https://files.rcsb.org/download/$id.cif"; done
get 9P3Y.cif "https://files.rcsb.org/download/9P3Y.cif"
get 8dz0_2fofc.ccp4 "https://www.ebi.ac.uk/pdbe/entry-files/download/8dz0.ccp4"
