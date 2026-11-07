# PROTFIG benchmark

Every preset is run on PDB entries that are not used anywhere else in PROTFIG (examples, paper
figures, development tests), and failures and times are reported.

```sh
cd benchmark
python3 run_benchmark.py                    # all presets (PYMOL=/path/to/pymol if pymol is not on PATH)
python3 run_benchmark.py --preset density   # one preset
```

## Cases

[`cases.csv`](cases.csv): one row per figure, 3 per preset (for `site`, one per partner type:
protein, peptide, ligand, the hardest of each).

| Column | Meaning |
|---|---|
| `preset` | preset to run |
| `id` | PDB entry (downloaded from the RCSB PDB) or a model in [`data/`](data) |
| `type` | partner type: `protein`, `peptide`, `ligand`, or `-` for presets without a partner |
| `expected_partner` | the partner PROTFIG should pick: `chain X` or `resn XXX` (`-`: not checked) |
| `options` | extra options; `{in}` is the inputs folder (maps and second structures are downloaded too) |

Cases include harder situations on purpose: cofactors next to the ligand (4COX), lipids or
cholesterol in the crystal (2RH1, 3EML, 1C3W, 5C1M), glycans (1A4G, 1FSS, 1IGT), a peptide with an
acetyl cap (1EVH), several copies in the crystal (1HHK, 2HYY), an enzyme split in several chains
(1CHO), a peptidomimetic covalent inhibitor (6LU7).

Inputs are downloaded into `inputs/` (not timed, not in git): structures from the RCSB PDB, 2Fo-Fc
maps from PDBe. The AlphaFold models for `plddt` are in [`data/`](data) (AlphaFold Protein Structure
Database, CC BY 4.0; Jumper et al. 2021, Varadi et al. 2024), because the database does not serve
scripted downloads. Figures go to `figures/` (not in git).

`membrane` is not benchmarked yet: it needs simulated bilayers (MemProtMD), whose server was not
reachable when the benchmark was run.

## Errors and improvements

- **Automatic** (`results/<preset>_results.csv`): the run failed; the detected partner is not the
  expected one (type and chain / residue name); blank image. The number of H-bonds is reported
  (`n_hbonds`) but zero is not an error: some ligands bind without direct H-bonds (by aromatic
  stacking or through waters).
- **Visual** (`results/<preset>_review.csv`): checked by a person looking at each figure. Two
  columns, a number each (`0` if none; blank = not reviewed yet; a description typed here counts
  as 1 and is moved to `visual_notes`):
  - `visual_errors`: the figure is wrong or misleading (e.g. a crystal copy drawn as if it were part
    of the complex).
  - `visual_improvements`: the figure is right and readable but could look better (framing,
    orientation, a buried ligand hard to see).

  Then update the summary without rendering again: `python3 run_benchmark.py --summary_only`. A
  full run keeps the reviews.

[`results/summary.csv`](results/summary.csv) gives, per preset and in total, the figures with errors
(automatic or visual; a figure counts once) and the error rate, the figures without errors that could
be improved, the figures with no remark (`figures_ok`), and the time per figure. Until the review is
done these columns read `pending`. [`benchmark_errors.txt`](benchmark_errors.txt) describes every
error and improvement: what is seen, the cause, a workaround and the fix planned.

## Times

Wall clock of the whole `protfig <file> --preset <preset>` command (PyMOL start-up, analysis, ray
tracing), default options (2126 × 1594 px, 300 dpi), one figure at a time. For each figure the
first run is discarded (warm-up) and the median of the next 3 is reported, with the range.
`setup_s` (loading, analysis, view) and `render_s` (ray tracing, writing files) are PROTFIG's own
split of the last run. [`results/environment.txt`](results/environment.txt) records the machine and
the software versions.

## Columns of `<preset>_results.csv`

| Column | Meaning |
|---|---|
| `status` | `ok` or `failed` |
| `partner_expected`, `partner_detected`, `partner_ok` | partner check |
| `chains_shown` | chains drawn (crystal copies away from the partner are left out) |
| `n_hbonds` | partner–receptor H-bonds drawn (same as `<name>_hbonds.txt`) |
| `auto_errors`, `auto_error_notes` | number and description of automatic errors |
| `time_median_s`, `time_min_s`, `time_max_s` | wall-clock time per figure |
| `setup_s`, `render_s` | PROTFIG's split of the last run |
