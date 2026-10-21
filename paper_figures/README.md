# Paper figures

The panels of the PROTFIG article, and the file that generates them.

- `figures.toml`: one `[[figure]]` per panel, with the exact options used.
- `get_inputs.sh`: downloads the public inputs (PDB entries 6XOX, 7KI0, 8YW3, 9SV3, 9DTY, 9ZTJ
  and 9JN0, plus a MemProtMD simulation of 6XOX) and builds the membrane system.
- `panels/<preset>/<name>/<name>.png`: the rendered panels.

To render them again (outputs overwrite `panels/`):

```sh
./get_inputs.sh
../protfig --config figures.toml --jobs 8
```
