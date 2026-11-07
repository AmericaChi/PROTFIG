"""
run_benchmark.py — reproducible test of the PROTFIG presets on PDB entries not used elsewhere.

  python3 run_benchmark.py [cases.csv] [--preset site] [--repeats 3] [--warmup 1] [--summary_only]

Each row of cases.csv is one figure: preset, PDB id (or a model in data/), partner type, expected
partner (`chain X`, `resn XXX`, or `-` for presets without a partner) and extra options ({in} = the
inputs folder). Structures and maps are downloaded first (not timed); then each figure is made with
`protfig <file> --preset <preset> [options]`, one at a time, and the script writes:

  results/<preset>_results.csv   one row per figure: automatic checks and times
  results/<preset>_review.csv    visual review, filled in by hand (kept between runs)
  results/summary.csv            errors, possible improvements and times per preset and in total
  results/environment.txt        machine, OS, PyMOL, Python and PROTFIG version

Times: wall clock of the whole command (PyMOL start-up, analysis, ray tracing), one process at a
time; the first --warmup runs of each figure are discarded and the median of --repeats runs is
reported. setup/render are PROTFIG's own split of the last run (from <name>_run.txt).
Uses only the Python standard library (Pillow, if present, adds the blank-image check).
"""
import argparse
import csv
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PROTFIG = os.path.join(HERE, "..", "protfig")
INPUTS, DATA, FIGURES, RESULTS = (os.path.join(HERE, d) for d in ("inputs", "data", "figures", "results"))


def sh(command):
    try:
        return subprocess.run(command, shell=True, capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception:
        return ""


def environment():
    pymol = os.environ.get("PYMOL", "pymol")
    mem = sh("sysctl -n hw.memsize")
    lines = [
        f"date:      {time.strftime('%Y-%m-%d %H:%M')}",
        f"os:        {platform.platform()}",
        f"cpu:       {sh('sysctl -n machdep.cpu.brand_string') or sh('lscpu | grep -m1 \"Model name\" | cut -d: -f2') or platform.processor()}",
        f"cores:     {os.cpu_count()}",
        f"memory:    {str(round(int(mem) / 2**30)) + ' GB' if mem else sh('grep MemTotal /proc/meminfo')}",
        f"python:    {platform.python_version()}",
        f"pymol:     {(sh(pymol + ' -cq -d \"print(cmd.get_version()[0])\"').splitlines() or [''])[-1]}",
        f"protfig:   {sh('git -C ' + HERE + ' rev-parse --short HEAD') or 'unknown'}"
        + (" (with local changes)" if sh('git -C ' + HERE + ' status --porcelain -- ..') else ""),
    ]
    return "\n".join(lines) + "\n"


def fetch(name):
    """inputs/<name>: XXXX.cif from the RCSB PDB, xxxx.ccp4 (2Fo-Fc map) from PDBe."""
    path = os.path.join(INPUTS, name)
    if os.path.exists(path) and os.path.getsize(path):
        return path
    stem, ext = os.path.splitext(name)
    url = {".cif": f"https://files.rcsb.org/download/{stem.upper()}.cif",
           ".ccp4": f"https://www.ebi.ac.uk/pdbe/entry-files/download/{stem.lower()}.ccp4"}[ext]
    urllib.request.urlretrieve(url, path)
    return path


def case_input(c):
    """The structure of a case: a model shipped in data/ (e.g. AlphaFold DB), else <id>.cif."""
    local = os.path.join(DATA, f"{c['id']}.pdb")
    return local if os.path.exists(local) else fetch(f"{c['id']}.cif")


def case_options(c):
    opts = (c.get("options") or "").replace("{in}", INPUTS).split()
    for o in opts:  # files named in the options (maps, second structures) are downloaded too
        if o.startswith(INPUTS):
            fetch(os.path.basename(o))
    return opts


def read_run(prefix):
    """partner, chains and PROTFIG's setup/render times from <name>_run.txt."""
    out = {"partner": "", "chains": "", "setup": "", "render": ""}
    try:
        txt = open(prefix + "_run.txt").read()
    except OSError:
        return out
    for key in ("partner", "chains"):
        m = re.search(rf"^{key}:\s+(.*)$", txt, re.M)
        out[key] = m.group(1).strip() if m else ""
    m = re.search(r"setup ([\d.]+) s.*render ([\d.]+) s", txt)
    if m:
        out["setup"], out["render"] = m.group(1), m.group(2)
    return out


def n_hbonds(prefix):
    try:
        m = re.search(r"# n = (\d+)", open(prefix + "_hbonds.txt").read())
        return int(m.group(1)) if m else 0
    except OSError:
        return 0


def blank_image(png):
    """True if (almost) nothing was drawn; None if Pillow is not available."""
    try:
        from PIL import Image
    except ImportError:
        return None
    px = list(Image.open(png).convert("RGBA").resize((200, 150)).getdata())
    return sum(1 for r, g, b, a in px if a > 0 and (r, g, b) != (255, 255, 255)) / len(px) < 0.01


def partner_ok(case_type, expected, detected):
    """The detected partner (from _run.txt) has the expected type and chain / residue name."""
    if not detected.startswith(case_type):
        return False
    what, value = expected.split()
    return f"(chain {value})" in detected if what == "chain" else f" {value} " in f" {detected} "


def run_case(c, repeats, warmup):
    pdb, preset = c["id"], c["preset"]
    command = [PROTFIG, case_input(c), "--preset", preset, "--name", pdb, "--out", FIGURES] + case_options(c)
    times, log, code = [], "", 0
    for k in range(warmup + repeats):
        t = time.perf_counter()
        r = subprocess.run(command, capture_output=True, text=True)
        dt = time.perf_counter() - t
        log, code = r.stdout + r.stderr, r.returncode
        if k >= warmup:
            times.append(dt)
    prefix = os.path.join(FIGURES, preset, pdb, pdb)
    run = read_run(prefix)
    failed = code != 0 or "ERROR" in log or "Traceback" in log or not os.path.exists(prefix + ".png")
    errors = []
    if failed:
        last = [ln for ln in log.splitlines() if "ERROR" in ln or "Error" in ln]
        errors.append("run failed: " + (last[-1].strip() if last else f"exit code {code}"))
    else:
        if c["expected_partner"] != "-" and not partner_ok(c["type"], c["expected_partner"], run["partner"]):
            errors.append(f"partner: expected {c['type']} {c['expected_partner']}, got {run['partner']}")
        if blank_image(prefix + ".png"):
            errors.append("blank image")
    return {
        "preset": preset, "id": pdb, "type": c["type"], "description": c["description"],
        "options": c.get("options", ""), "status": "failed" if failed else "ok",
        "partner_expected": c["expected_partner"], "partner_detected": run["partner"],
        "partner_ok": "" if failed or c["expected_partner"] == "-" else ("no" if any(e.startswith("partner") for e in errors) else "yes"),
        "chains_shown": run["chains"], "n_hbonds": n_hbonds(prefix),
        "auto_errors": len(errors), "auto_error_notes": "; ".join(errors),
        "time_median_s": f"{statistics.median(times):.1f}", "time_min_s": f"{min(times):.1f}",
        "time_max_s": f"{max(times):.1f}", "setup_s": run["setup"], "render_s": run["render"], "runs": repeats,
    }


def visual_count(entry, column="visual_errors"):
    """visual_errors / visual_improvements: a number (0 = reviewed, none); free text = 1 (the text
    is a note); blank = not reviewed yet (None)."""
    v = (entry.get(column) or "").strip()
    return None if not v else (int(v) if v.isdigit() else 1)


def write_review(rows, preset):
    """<preset>_review.csv, one row per figure, kept between runs. A description typed in
    visual_errors or visual_improvements is moved to visual_notes and counted as 1."""
    path = os.path.join(RESULTS, f"{preset}_review.csv")
    old = {r["id"]: r for r in csv.DictReader(open(path, newline=""))} if os.path.exists(path) else {}
    review = {}
    for r in rows:
        e = old.get(r["id"], {})
        notes = (e.get("visual_notes") or "").strip()
        entry = {"id": r["id"], "type": r["type"]}
        for col in ("visual_errors", "visual_improvements"):
            n = (e.get(col) or "").strip()
            if n and not n.isdigit():
                n, notes = "1", "; ".join(x for x in (n, notes) if x)
            entry[col] = n
        entry["visual_notes"] = notes
        review[r["id"]] = entry
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "type", "visual_errors", "visual_improvements", "visual_notes"])
        w.writeheader()
        w.writerows(review.values())
    return review


def summary_row(group, rows, review):
    """Errors: the figure is wrong (run failed, wrong partner, misleading figure). Improvements: the
    figure is right and readable but could look better. A figure counts once in each column."""
    visual = [visual_count(review[r["id"]]) for r in rows]
    better = [visual_count(review[r["id"]], "visual_improvements") for r in rows]
    reviewed = sum(v is not None for v in visual)
    errors = [int(r["auto_errors"]) > 0 or (v or 0) > 0 for r, v in zip(rows, visual)]
    improve = [not e and (b or 0) > 0 for e, b in zip(errors, better)]
    t = [float(r["time_median_s"]) for r in rows if r["status"] == "ok"]
    pct = lambda k: f"{100 * k / len(rows):.0f}"
    return {
        "group": group, "figures": len(rows),
        "failed_runs": sum(r["status"] == "failed" for r in rows),
        "figures_auto_errors": sum(int(r["auto_errors"]) > 0 for r in rows),
        "figures_reviewed": reviewed,
        "figures_visual_errors": sum((v or 0) > 0 for v in visual) if reviewed else "pending",
        "figures_with_errors": sum(errors),
        "error_rate_pct": pct(sum(errors)),
        "figures_to_improve": sum(improve) if reviewed else "pending",
        "improve_rate_pct": pct(sum(improve)) if reviewed else "pending",
        "figures_ok": len(rows) - sum(errors) - sum(improve) if reviewed else "pending",
        "time_median_s": f"{statistics.median(t):.1f}" if t else "",
        "time_range_s": f"{min(t):.1f}-{max(t):.1f}" if t else "",
    }


def write_summary():
    """results/summary.csv from every <preset>_results.csv and its review."""
    rows_all, review_all, out = [], {}, []
    for f in sorted(os.listdir(RESULTS)):
        if not f.endswith("_results.csv"):
            continue
        preset = f[: -len("_results.csv")]
        rows = list(csv.DictReader(open(os.path.join(RESULTS, f))))
        review = write_review(rows, preset)
        out.append(summary_row(preset, rows, review))
        rows_all += rows
        review_all.update({(preset, k): v for k, v in review.items()})
    if not rows_all:
        return
    keyed = [dict(r, id=(r["preset"], r["id"])) for r in rows_all]
    out.append(summary_row("all", keyed, review_all))
    with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cases", nargs="?", default=os.path.join(HERE, "cases.csv"))
    p.add_argument("--preset", help="only the cases of this preset")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--summary_only", action="store_true",
                   help="only update the summary from the existing results and the visual reviews (no rendering)")
    a = p.parse_args()
    for d in (INPUTS, FIGURES, RESULTS):
        os.makedirs(d, exist_ok=True)
    if a.summary_only:
        write_summary()
        print(os.path.join(RESULTS, "summary.csv"))
        return
    with open(os.path.join(RESULTS, "environment.txt"), "w") as fh:
        fh.write(environment())
    cases = [c for c in csv.DictReader(open(a.cases)) if not a.preset or c["preset"] == a.preset]
    for c in cases:  # downloads are not timed
        case_input(c)
        case_options(c)
    by_preset = {}
    for c in cases:
        r = run_case(c, a.repeats, a.warmup)
        by_preset.setdefault(c["preset"], []).append(r)
        print(f"{r['preset']:13s} {r['id']:10s} {r['status']:6s} {r['time_median_s']:>6s} s  {r['auto_error_notes']}", flush=True)
    for preset, rows in by_preset.items():
        with open(os.path.join(RESULTS, f"{preset}_results.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    write_summary()
    print(f"\n{RESULTS}/summary.csv  (visual review: fill in <preset>_review.csv, then --summary_only)")


if __name__ == "__main__":
    sys.exit(main())
