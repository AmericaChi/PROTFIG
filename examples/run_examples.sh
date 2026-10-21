#!/bin/sh
# Downloads the example structures and renders every figure in figures.toml + the gallery.
cd "$(dirname "$0")"
./get_inputs.sh
../protfig --config figures.toml --jobs "${JOBS:-4}"
${PYMOL:-pymol} -cq make_gallery.py -- figures ../docs
