"""
protfig_web — web interface for PROTFIG (works locally and on HPC nodes through an SSH tunnel).

    ./protfig-web                         # http://127.0.0.1:8050/?token=...
    ./protfig-web --port 9000 --workdir ~/protfig_runs
    ./protfig-web --host 0.0.0.0          # compute node reached through a tunnel from the login node

Every figure runs `pymol -cq protfig.py -- <options>` in a subprocess, exactly as on the command line,
with its own folder <workdir>/jobs/<id>/ (log, settings, outputs). The form is built from
protfig.build_parser(), so new command-line options appear in the interface automatically.
"""
import argparse
import datetime
import glob
import io
import json
import os
import secrets
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import uuid
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # repository root (protfig.py, protfig_lib/, docs/)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from flask import Flask, abort, jsonify, redirect, request, send_file, send_from_directory  # noqa: E402
from werkzeug.utils import secure_filename  # noqa: E402

import protfig  # noqa: E402  (imports pymol.cmd, but nothing is run in this process)

SCRIPT = os.path.join(ROOT, "protfig.py")
THUMBS = os.path.join(ROOT, "docs", "presets")
# options managed by the interface itself, never shown as form fields
HIDDEN = {"help", "input", "prefix", "out", "config", "jobs", "list", "draft", "preset", "style"}
# options whose value is a file: the form offers Upload and Browse for them
FILE_OPTS = {"input", "ref", "map", "values", "view", "session"}
BROWSE_EXT = (".pdb", ".ent", ".cif", ".mmcif", ".pse", ".ccp4", ".map", ".mrc", ".txt", ".gz")
MARKER = "[PROTFIG:"

app = Flask(__name__, static_folder=os.path.join(HERE, "static"), static_url_path="/static")
CFG = {"token": None, "workdir": None, "launch_dir": os.getcwd(), "pymol": None}
JOBS = {}
JOBS_LOCK = threading.Lock()
SLOTS = threading.Semaphore(2)


# --- access token (like Jupyter): other users on a shared HPC node can reach localhost ports

@app.before_request
def check_token():
    tok = CFG["token"]
    if not tok or request.path.startswith("/static/"):
        return None
    given = request.args.get("token")
    if given and secrets.compare_digest(given, tok):
        resp = redirect(request.path)
        resp.set_cookie("protfig_token", tok, httponly=True, samesite="Strict")
        return resp
    sent = request.cookies.get("protfig_token") or request.headers.get("X-Protfig-Token") or ""
    if secrets.compare_digest(sent, tok):
        return None
    return ("PROTFIG web: open the link with ?token=... printed when the server started.", 403)


# --- option schema, from the command-line parser

def option_schema():
    """Groups of form fields built from protfig.build_parser()."""
    p = protfig.build_parser()
    groups = []
    for g in p._action_groups:
        fields = []
        for act in g._group_actions:
            if act.help == argparse.SUPPRESS or act.dest in HIDDEN or not act.option_strings:
                continue
            flag = next((o for o in act.option_strings if o.startswith("--")), act.option_strings[0])
            if flag[2:] != act.dest:
                continue  # alternative spellings that store into another option (e.g. --transparent)
            f = {"dest": act.dest, "flag": flag, "help": act.help or "", "default": act.default,
                 "file": act.dest in FILE_OPTS}
            if isinstance(act, argparse._StoreTrueAction):
                f["kind"], f["default"] = "flag", False
            elif isinstance(act, argparse._AppendAction):
                f["kind"], f["default"] = "list", ""
            elif act.choices:
                f["kind"], f["choices"] = "choice", list(act.choices)
            elif act.nargs == 2:
                f["kind"], f["names"] = "pair", list(act.metavar or ("min", "max"))
            elif act.type in (int, float):
                f["kind"], f["numtype"] = "number", act.type.__name__
            else:
                f["kind"] = "text"
            fields.append(f)
        if fields:
            groups.append({"title": g.title, "fields": fields})
    return groups


def schema_index():
    return {f["dest"]: f for g in option_schema() for f in g["fields"]}


def thumbnail(preset, style=None):
    if not os.path.isdir(THUMBS):
        return None
    files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(THUMBS, f"{preset}__*.jpg")))
    all_styles = {s for st in protfig.STYLES.values() for s in st}
    if style:
        pick = [f for f in files if f[:-4].endswith("_" + style)]
    else:
        pick = [f for f in files if not any(f[:-4].endswith("_" + s) for s in all_styles)]
    pick = pick or files
    return f"/thumbs/{pick[0]}" if pick else None


def presets_info():
    import tomllib

    built_in = []
    for name, m in sorted(protfig.PRESETS.items()):
        styles = [{"name": "", "thumb": thumbnail(name)}]
        styles += [{"name": s, "thumb": thumbnail(name, s)} for s in protfig.STYLES.get(name, ())]
        built_in.append({"name": name, "doc": m["doc"], "partner": m["partner"],
                         "kinds": [protfig.KIND_NAMES[k] for k in m["kinds"]], "styles": styles})
    saved = []
    for n in protfig.user_presets():
        try:
            with open(protfig.user_preset_path(n), "rb") as fh:
                opts = tomllib.load(fh)
            saved.append({"name": n, "base": str(opts.get("preset", "overview")), "command": protfig.cli_text(opts)})
        except (OSError, tomllib.TOMLDecodeError) as e:
            saved.append({"name": n, "base": "?", "command": f"(cannot read: {e})"})
    return {"built_in": built_in, "saved": saved, "folder": protfig.USER_PRESETS}


# --- form -> argv

class FormError(ValueError):
    pass


def form_to_argv(form):
    """{'input':..., 'preset':..., 'style':..., 'draft':bool, 'options':{dest: value}} -> argv."""
    inp = str(form.get("input") or "").strip()
    if not inp:
        raise FormError("choose a structure: upload a file, browse the server, or type a PDB code")
    argv = [inp, "--preset", str(form.get("preset") or "overview")]
    if form.get("style"):
        argv += ["--style", str(form["style"])]
    index = schema_index()
    opts = {}
    for dest, v in (form.get("options") or {}).items():
        f = index.get(dest)
        if f is None:
            raise FormError(f"unknown option '{dest}'")
        if f["kind"] == "flag":
            if v is True:
                opts[dest] = True
        elif f["kind"] == "list":
            items = [x.strip() for x in str(v or "").splitlines() if x.strip()]
            if items:
                opts[dest] = items
        elif f["kind"] == "pair":
            vals = [str(x).strip() for x in (v or []) if str(x).strip()]
            if len(vals) == 2:
                opts[dest] = vals
            elif vals:
                raise FormError(f"{f['flag']} needs two values")
        else:
            s = str(v if v is not None else "").strip()
            if s and s != str(f["default"] if f["default"] is not None else ""):
                opts[dest] = s
    argv += protfig.opts_to_argv(opts)
    if form.get("draft"):
        argv.append("--draft")
    return argv


def command_text(argv):
    return "protfig " + " ".join(shlex.quote(x) for x in argv)


# --- jobs

class Job:
    def __init__(self, jid, argv, form, folder, preview=False):
        self.id, self.argv, self.form, self.folder = jid, argv, form, folder
        self.preview = preview  # draft preview: temporary, only the latest one is kept
        self.discard = False  # delete its folder once it stops
        self.status, self.returncode, self.result = "queued", None, None
        self.created = time.time()
        self.started = self.finished = None
        self.log = []
        self.proc = None
        self.cancelled = False

    def meta(self):
        return {"id": self.id, "preview": self.preview, "argv": self.argv, "form": self.form,
                "command": command_text(self.argv),
                "status": self.status, "returncode": self.returncode, "result": self.result,
                "created": self.created, "started": self.started, "finished": self.finished}

    def save(self):
        try:
            with open(os.path.join(self.folder, "job.json"), "w") as fh:
                json.dump(self.meta(), fh, indent=1)
        except OSError:
            pass  # a discarded preview whose folder is already gone


def jobs_dir():
    return os.path.join(CFG["workdir"], "jobs")


def previews_dir():
    return os.path.join(CFG["workdir"], "previews")


def stop(job):
    job.cancelled = True
    if job.proc and job.proc.poll() is None:
        try:
            os.killpg(job.proc.pid, signal.SIGTERM)
        except OSError:
            job.proc.terminate()


def remove(job):
    shutil.rmtree(job.folder, ignore_errors=True)
    with JOBS_LOCK:
        JOBS.pop(job.id, None)


def discard_previews():
    """A new preview replaces the previous ones: stop them and delete their files."""
    with JOBS_LOCK:
        old = [j for j in JOBS.values() if j.preview]
    for job in old:
        job.discard = True
        if job.status in ("queued", "running"):
            stop(job)  # run_job removes it when the process ends
        else:
            remove(job)


def run_job(job):
    try:
        _run_job(job)
    finally:
        if job.discard:
            remove(job)


def _run_job(job):
    with SLOTS:
        if job.cancelled:
            job.status = "cancelled"
            return
        job.status, job.started = "running", time.time()
        job.save()
        cmdline = [CFG["pymol"], "-cq", SCRIPT, "--"] + job.argv + ["--out", os.path.join(job.folder, "figures")]
        logf = open(os.path.join(job.folder, "log.txt"), "w")
        try:
            job.proc = subprocess.Popen(cmdline, cwd=CFG["launch_dir"], stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, text=True, bufsize=1,
                                        start_new_session=True)  # own process group: cancel stops all of it
            for line in job.proc.stdout:
                line = line.rstrip("\n")
                job.log.append(line)
                logf.write(line + "\n")
                logf.flush()
                if line.startswith(MARKER):
                    path = line.split("] ", 1)[1].split("  ", 1)[0].strip()
                    job.result = os.path.relpath(path, job.folder)
            job.returncode = job.proc.wait()
        except OSError as e:
            job.log.append(f"ERROR: cannot start PyMOL ({CFG['pymol']}): {e}")
            job.returncode = -1
        finally:
            logf.close()
        ok = job.returncode == 0 and (job.result or "--save_preset" in job.argv or "--no_png" in job.argv)
        job.status = "cancelled" if job.cancelled else ("done" if ok else "failed")
        job.finished = time.time()
        job.save()


def load_jobs():
    for path in glob.glob(os.path.join(jobs_dir(), "*", "job.json")):
        try:
            with open(path) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            continue
        job = Job(m["id"], m["argv"], m.get("form", {}), os.path.dirname(path))
        for k in ("status", "returncode", "result", "created", "started", "finished"):
            setattr(job, k, m.get(k))
        if job.status in ("queued", "running"):
            job.status = "failed"  # the server stopped while it was running
            job.log = ["(interrupted: the server stopped during this job)"]
        try:
            with open(os.path.join(job.folder, "log.txt")) as fh:
                job.log = fh.read().splitlines() or job.log
        except OSError:
            pass
        JOBS[job.id] = job


def get_job(jid):
    job = JOBS.get(jid)
    if job is None:
        abort(404)
    return job


def output_files(job):
    """Files of the figure folder (where the result is), or every file of the job's outputs."""
    figs = os.path.join(job.folder, "figures")
    base = os.path.join(job.folder, os.path.dirname(job.result)) if job.result else figs
    out = []
    for dirpath, _, names in os.walk(base):
        for n in sorted(names):
            full = os.path.join(dirpath, n)
            rel = os.path.relpath(full, job.folder)
            out.append({"name": os.path.relpath(full, base), "path": rel, "size": os.path.getsize(full),
                        "url": f"/api/jobs/{job.id}/files/{rel}"})
    return out


# --- routes

@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.get("/thumbs/<path:name>")
def thumbs(name):
    return send_from_directory(THUMBS, name)


@app.get("/api/schema")
def api_schema():
    return jsonify({"groups": option_schema(), "presets": presets_info(),
                    "launch_dir": CFG["launch_dir"], "host": socket.gethostname()})


@app.get("/api/presets")
def api_presets():
    return jsonify(presets_info())


@app.post("/api/preview")
def api_preview():
    try:
        return jsonify({"command": command_text(form_to_argv(request.get_json(force=True)))})
    except FormError as e:
        return jsonify({"command": None, "error": str(e)})


@app.post("/api/upload")
def api_upload():
    f = request.files.get("file")
    if not f or not f.filename:
        return jsonify({"error": "no file"}), 400
    name = secure_filename(f.filename) or "upload"
    folder = os.path.join(CFG["workdir"], "uploads", uuid.uuid4().hex[:8])
    os.makedirs(folder)
    path = os.path.join(folder, name)  # own folder: the figure is named after the original file
    f.save(path)
    return jsonify({"path": path})


@app.get("/api/browse")
def api_browse():
    path = os.path.abspath(os.path.expanduser(request.args.get("path") or CFG["launch_dir"]))
    if os.path.isfile(path):
        path = os.path.dirname(path)
    try:
        entries = sorted(os.scandir(path), key=lambda e: e.name.lower())
    except OSError as e:
        return jsonify({"error": f"{path}: {e.strerror}"}), 400
    dirs, files = [], []
    for e in entries:
        if e.name.startswith("."):
            continue
        try:
            if e.is_dir():
                dirs.append(e.name)
            elif e.name.lower().endswith(BROWSE_EXT):
                files.append({"name": e.name, "size": e.stat().st_size})
        except OSError:
            continue
    return jsonify({"path": path, "parent": os.path.dirname(path), "dirs": dirs, "files": files})


@app.post("/api/jobs")
def api_submit():
    form = request.get_json(force=True)
    try:
        argv = form_to_argv(form)
    except FormError as e:
        return jsonify({"error": str(e)}), 400
    preview = bool(form.get("draft"))
    if preview:
        discard_previews()
    jid = datetime.datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:4]
    folder = os.path.join(previews_dir() if preview else jobs_dir(), jid)
    os.makedirs(folder)
    job = Job(jid, argv, form, folder, preview)
    job.save()
    with JOBS_LOCK:
        JOBS[jid] = job
    threading.Thread(target=run_job, args=(job,), daemon=True).start()
    return jsonify(job.meta())


@app.get("/api/jobs")
def api_jobs():
    jobs = sorted((j for j in JOBS.values() if not j.preview), key=lambda j: j.created, reverse=True)
    return jsonify([j.meta() for j in jobs])


@app.get("/api/jobs/<jid>")
def api_job(jid):
    job = get_job(jid)
    since = request.args.get("since", 0, type=int)
    m = job.meta()
    m["log"], m["log_len"] = job.log[since:], len(job.log)
    if job.status in ("done", "failed", "cancelled"):
        m["files"] = output_files(job)
        if job.result and job.result.lower().endswith((".png", ".jpg")):
            m["image"] = f"/api/jobs/{job.id}/files/{job.result}"
        elif job.result and job.result.lower().endswith((".tif", ".tiff")):
            m["image"] = f"/api/jobs/{job.id}/preview"  # browsers cannot show TIFF
    return jsonify(m)


@app.get("/api/jobs/<jid>/preview")
def api_preview_image(jid):
    """PNG copy of a TIFF-only result, for display in the browser."""
    from PIL import Image

    job = get_job(jid)
    if not job.result:
        abort(404)
    im = Image.open(os.path.join(job.folder, job.result))
    im.thumbnail((2400, 2400))
    buf = io.BytesIO()
    im.save(buf, "PNG")
    buf.seek(0)
    return send_file(buf, mimetype="image/png")


@app.post("/api/jobs/<jid>/cancel")
def api_cancel(jid):
    job = get_job(jid)
    stop(job)
    if job.status == "queued":
        job.status = "cancelled"
        job.save()
    return jsonify(job.meta())


@app.delete("/api/jobs/<jid>")
def api_delete(jid):
    job = get_job(jid)
    if job.status in ("queued", "running"):
        return jsonify({"error": "cancel the job first"}), 400
    remove(job)
    return jsonify({"ok": True})


@app.get("/api/jobs/<jid>/files/<path:rel>")
def api_file(jid, rel):
    job = get_job(jid)
    return send_from_directory(job.folder, rel, as_attachment=request.args.get("download") == "1")


@app.get("/api/jobs/<jid>/zip")
def api_zip(jid):
    job = get_job(jid)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in output_files(job):
            z.write(os.path.join(job.folder, f["path"]), f["name"])
        z.writestr("command.txt", command_text(job.argv) + "\n")
    buf.seek(0)
    name = os.path.splitext(os.path.basename(job.result))[0] if job.result else job.id
    return send_file(buf, mimetype="application/zip", as_attachment=True, download_name=f"{name}.zip")


# --- start-up

def find_pymol():
    env = os.environ.get("PYMOL")
    if env:
        return shutil.which(env) or env
    sibling = os.path.join(os.path.dirname(sys.executable), "pymol")
    return sibling if os.path.exists(sibling) else shutil.which("pymol")


def main(argv=None):
    p = argparse.ArgumentParser(prog="protfig-web", description="PROTFIG web interface.")
    p.add_argument("--host", default="127.0.0.1",
                   help="address to listen on (default 127.0.0.1; 0.0.0.0 on a compute node reached from the login node)")
    p.add_argument("--port", type=int, default=8050)
    p.add_argument("--workdir", default="protfig_web_runs",
                   help="folder for uploads and job outputs (default ./protfig_web_runs)")
    p.add_argument("--max_jobs", type=int, default=2, help="figures rendered at the same time (default 2)")
    p.add_argument("--max_upload_mb", type=int, default=2048, help="largest upload in MB (default 2048)")
    p.add_argument("--no_token", action="store_true", help="no access token (only on a single-user machine)")
    a = p.parse_args(argv)

    CFG["pymol"] = find_pymol()
    if not CFG["pymol"]:
        sys.exit("protfig-web: PyMOL not found. Activate the protfig environment or set PYMOL=/path/to/pymol")
    CFG["workdir"] = os.path.abspath(os.path.expanduser(a.workdir))
    os.makedirs(jobs_dir(), exist_ok=True)
    shutil.rmtree(previews_dir(), ignore_errors=True)  # previews are temporary
    CFG["token"] = None if a.no_token else secrets.token_urlsafe(24)
    app.config["MAX_CONTENT_LENGTH"] = a.max_upload_mb * 1024 * 1024
    global SLOTS
    SLOTS = threading.Semaphore(max(1, a.max_jobs))
    load_jobs()

    host = socket.gethostname()
    url = f"http://localhost:{a.port}/" + (f"?token={CFG['token']}" if CFG["token"] else "")
    print(f"PROTFIG web on {host}:{a.port}   (PyMOL: {CFG['pymol']}; outputs: {CFG['workdir']})")
    print(f"\n  Open:  {url}\n")
    print("  From your computer, if this runs on a remote machine (HPC), first open a tunnel:")
    if a.host in ("127.0.0.1", "localhost"):
        print(f"    ssh -N -L {a.port}:localhost:{a.port} <user>@{host}")
    else:
        print(f"    ssh -N -L {a.port}:{host}:{a.port} <user>@<login-node>")
    print("  Stop with Ctrl+C.\n", flush=True)
    app.run(host=a.host, port=a.port, threaded=True, debug=False)


if __name__ == "__main__":
    main()
