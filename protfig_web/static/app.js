// PROTFIG web client: builds the form from /api/schema, submits jobs and polls their progress.
"use strict";

const $ = (s) => document.querySelector(s);
const el = (tag, attrs = {}, ...kids) => {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined && v !== false) e.setAttribute(k, v === true ? "" : v);
  }
  for (const k of kids.flat()) if (k !== null && k !== undefined) e.append(k);
  return e;
};
const api = async (url, opts = {}) => {
  const r = await fetch(url, opts);
  const data = await r.json().catch(() => ({ error: `HTTP ${r.status}` }));
  if (!r.ok && !data.error) data.error = `HTTP ${r.status}`;
  return data;
};
const postJSON = (url, body) =>
  api(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
const fmtSize = (n) => (n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(0)} KB` : `${(n / 1048576).toFixed(1)} MB`);
const fmtTime = (t) => (t ? new Date(t * 1000).toLocaleString() : "");

const state = { schema: null, fields: {}, form: { input: "", preset: "overview", style: "", options: {} },
                current: null, logLen: 0, poll: null };
const STORE_KEY = "protfig-web-form";

// ---------- form state ----------

function saveForm() {
  try { localStorage.setItem(STORE_KEY, JSON.stringify(state.form)); } catch (e) { /* storage unavailable */ }
}
function loadSavedForm() {
  try { return JSON.parse(localStorage.getItem(STORE_KEY) || "null"); } catch (e) { return null; }
}

function setForm(form) {
  state.form = { input: form.input || "", preset: form.preset || "overview", style: form.style || "",
                 options: { ...(form.options || {}) } };
  for (const k of Object.keys(state.form.options)) if (!state.fields[k]) delete state.form.options[k];
  $("#input").value = state.form.input;
  renderPresets();
  renderOptions();
  changed();
}

let previewTimer = null;
function changed() {
  saveForm();
  clearTimeout(previewTimer);
  previewTimer = setTimeout(async () => {
    const r = await postJSON("/api/preview", state.form);
    $("#command").textContent = r.command || `protfig …   (${r.error})`;
  }, 200);
}

// ---------- presets and styles ----------

function allPresets() { return state.schema.presets.built_in; }
function currentBuiltIn() {
  const p = state.form.preset;
  const saved = state.schema.presets.saved.find((s) => s.name === p);
  return allPresets().find((b) => b.name === (saved ? saved.base : p));
}

function card(label, thumb, selected, onclick, title) {
  return el("button", { class: "card" + (selected ? " sel" : ""), type: "button", onclick, title },
    thumb ? el("img", { src: thumb, alt: "", loading: "lazy" }) : el("div", { class: "noimg" }, "no preview"),
    el("span", {}, label));
}

function renderPresets() {
  const grid = $("#presets");
  grid.replaceChildren(...allPresets().map((p) =>
    card(p.name, p.styles[0].thumb, p.name === state.form.preset, () => {
      state.form.preset = p.name; state.form.style = ""; renderPresets(); changed();
    }, p.doc)));

  const saved = state.schema.presets.saved;
  $("#saved-presets").replaceChildren(...(saved.length ? [el("span", { class: "small muted" }, "Saved:")] : []),
    ...saved.map((s) => el("button", { class: "chip" + (s.name === state.form.preset ? " sel" : ""), type: "button",
      title: s.command, onclick: () => { state.form.preset = s.name; state.form.style = ""; renderPresets(); changed(); } },
      s.name)));

  const p = currentBuiltIn();
  const isSaved = p && p.name !== state.form.preset;
  const styles = $("#styles");
  styles.replaceChildren();
  if (p && p.styles.length > 1 && !isSaved) {
    styles.append(el("div", { class: "label" }, `Style of ${p.name}:`), ...p.styles.map((s) => card(s.name || "default", s.thumb, s.name === state.form.style, () => {
      state.form.style = s.name; renderPresets(); changed();
    })));
  }
  let doc = p ? p.doc : "";
  if (p && p.partner === "required") doc += " Needs a binding partner (automatic or --partner).";
  if (p && p.kinds.length) doc += ` Works with: ${p.kinds.join(", ")}.`;
  if (isSaved) {
    const s = saved.find((x) => x.name === state.form.preset);
    doc = `Saved preset based on '${s.base}': ${s.command}. Options set below override it.`;
  }
  $("#preset-doc").textContent = doc;
}

// ---------- options ----------

function fieldInput(f) {
  const opts = state.form.options;
  const set = (v) => {
    const empty = v === "" || v === false || v === null || (Array.isArray(v) && v.every((x) => x === ""));
    if (empty || (f.default !== null && String(v) === String(f.default))) delete opts[f.dest];
    else opts[f.dest] = v;
    markChanged(f.dest);
    changed();
  };
  const val = opts[f.dest];
  const ph = f.default !== null && f.default !== undefined && f.default !== false && f.default !== "" ? `default: ${f.default}` : "";
  if (f.kind === "flag") {
    return el("input", { type: "checkbox", id: "o-" + f.dest, checked: !!val, onchange: (e) => set(e.target.checked) });
  }
  if (f.kind === "choice") {
    const s = el("select", { id: "o-" + f.dest, onchange: (e) => set(e.target.value) },
      el("option", { value: "" }, f.default !== null && f.default !== undefined ? `(default: ${f.default})` : "(default)"),
      ...f.choices.map((c) => el("option", { value: c, selected: val === c }, c)));
    return s;
  }
  if (f.kind === "pair") {
    const v = Array.isArray(val) ? val : ["", ""];
    const a = el("input", { type: "number", step: "any", placeholder: f.names[0], value: v[0] });
    const b = el("input", { type: "number", step: "any", placeholder: f.names[1], value: v[1] });
    const upd = () => set([a.value, b.value]);
    a.oninput = upd; b.oninput = upd;
    return el("div", { class: "pair" }, a, b);
  }
  if (f.kind === "list") {
    return el("textarea", { id: "o-" + f.dest, placeholder: "one per line, e.g. ray_shadows=0",
      oninput: (e) => set(e.target.value) }, val || "");
  }
  const input = el("input", { type: f.kind === "number" ? "number" : "text", step: f.numtype === "int" ? "1" : "any",
    id: "o-" + f.dest, placeholder: ph, value: val ?? "", spellcheck: "false", oninput: (e) => set(e.target.value) });
  if (f.file) {
    return el("div", { class: "field-row" }, input,
      el("button", { class: "btn small", type: "button", "data-upload": "o-" + f.dest }, "Upload"),
      el("button", { class: "btn small", type: "button", "data-browse": "o-" + f.dest }, "Browse"));
  }
  return input;
}

function renderOptions() {
  const box = $("#options");
  const open = new Set([...box.querySelectorAll("details[open]")].map((d) => d.dataset.title));
  box.replaceChildren(...state.schema.groups.map((g) => {
    const d = el("details", { "data-title": g.title, open: open.has(g.title) || g.title === "basic" },
      el("summary", {}, g.title, el("span", { class: "badge", hidden: true })),
      ...g.fields.map((f) => el("div", { class: "opt", "data-dest": f.dest,
        "data-search": `${f.flag} ${f.help} ${g.title}`.toLowerCase() },
        el("label", { for: "o-" + f.dest }, f.flag), fieldInput(f), el("div", { class: "help" }, f.help))));
    return d;
  }));
  for (const k of Object.keys(state.fields)) markChanged(k);
  filterOptions();
}

function markChanged(dest) {
  const row = document.querySelector(`.opt[data-dest="${dest}"]`);
  if (row) row.classList.toggle("changed", dest in state.form.options);
  const n = Object.keys(state.form.options).length;
  $("#changed-count").textContent = n ? `${n} changed` : "";
  for (const d of document.querySelectorAll("#options details")) {
    const c = d.querySelectorAll(".opt.changed").length;
    const b = d.querySelector(".badge");
    b.hidden = !c; b.textContent = c;
  }
}

function filterOptions() {
  const q = $("#opt-search").value.trim().toLowerCase();
  for (const d of document.querySelectorAll("#options details")) {
    let any = false;
    for (const row of d.querySelectorAll(".opt")) {
      const show = !q || q.split(/\s+/).every((w) => row.dataset.search.includes(w));
      row.hidden = !show; any ||= show;
    }
    d.hidden = !any;
    if (q && any) d.open = true;
  }
}

// ---------- files: upload and server browser ----------

let uploadTarget = null;
function setTarget(id, value) {
  const input = document.getElementById(id);
  input.value = value;
  input.dispatchEvent(new Event("input"));
}

async function uploadFile(file) {
  const fd = new FormData();
  fd.append("file", file);
  showError(`Uploading ${file.name}…`);
  const r = await api("/api/upload", { method: "POST", body: fd });
  if (r.error) return showError(`Upload failed: ${r.error}`);
  showError(null);
  setTarget(uploadTarget, r.path);
}

let browseTarget = null;
async function browse(path) {
  const r = await api("/api/browse?path=" + encodeURIComponent(path || ""));
  const list = $("#br-list");
  if (r.error) { list.replaceChildren(el("li", {}, r.error)); return; }
  $("#br-path").value = r.path;
  list.replaceChildren(
    el("li", { class: "dir", onclick: () => browse(r.parent) }, ".."),
    ...r.dirs.map((d) => el("li", { class: "dir", onclick: () => browse(r.path + "/" + d) }, d)),
    ...r.files.map((f) => el("li", { onclick: () => { setTarget(browseTarget, r.path + "/" + f.name); $("#browser").close(); } },
      el("span", {}, f.name), el("span", { class: "muted small" }, fmtSize(f.size)))));
}

// ---------- jobs ----------

function showError(msg) {
  const e = $("#form-error");
  e.hidden = !msg; e.textContent = msg || "";
}

async function submit(draft, form = state.form) {
  showError(null);
  const r = await postJSON("/api/jobs", { ...form, draft });
  if (r.error) return showError(r.error);
  watch(r.id);
  refreshHistory();
}

function watch(id) {
  clearTimeout(state.poll);
  state.current = id; state.logLen = 0;
  $("#current").hidden = false;
  $("#cur-log").textContent = "";
  $("#cur-image").replaceChildren(el("div", { class: "spinner" }, "Rendering…"));
  $("#cur-files").replaceChildren();
  $("#cur-zip").hidden = true;
  $("#cur-full").hidden = true;
  $("#cur-note").textContent = "";
  $("#cur-log-box").open = false;
  poll();
}

async function poll() {
  const id = state.current;
  const j = await api(`/api/jobs/${id}?since=${state.logLen}`);
  if (id !== state.current || j.error) return;
  if (j.log.length) {
    $("#cur-log").textContent += j.log.join("\n") + "\n";
    $("#cur-log").scrollTop = 1e9;
  }
  state.logLen = j.log_len;
  const pill = $("#cur-status");
  pill.className = "pill " + j.status; pill.textContent = j.status;
  const t = j.finished && j.started ? `${(j.finished - j.started).toFixed(1)} s` : fmtTime(j.created);
  $("#cur-time").textContent = t;
  const active = j.status === "queued" || j.status === "running";
  $("#cur-cancel").hidden = !active;
  // saved figures get a link (#job=...); previews are temporary
  history.replaceState(null, "", j.preview ? location.pathname : "#job=" + id);
  $("#cur-note").textContent = j.preview ? "Draft preview, not saved: the next preview replaces it." : "";
  if (active) { state.poll = setTimeout(poll, 700); return; }

  $("#cur-image").replaceChildren(j.image ? el("a", { href: j.image, target: "_blank" },
    el("img", { src: j.image + "?t=" + Date.now(), alt: "figure" })) : el("div", { class: "spinner" },
    j.status === "done" ? "No image (see the files below)."
      : ($("#cur-log").textContent.split("\n").reverse().find((l) => l.startsWith("ERROR")) || "No figure: see the log.")));
  const files = j.preview ? [] : (j.files || []);
  $("#cur-files").replaceChildren(...files.map((f) => el("li", {},
    el("a", { href: f.url + "?download=1" }, f.name), " ", el("span", { class: "muted" }, fmtSize(f.size)))));
  const zip = $("#cur-zip");
  zip.hidden = !files.length; zip.href = `/api/jobs/${id}/zip`;
  const full = $("#cur-full");
  full.hidden = !(j.preview && j.status === "done");
  full.onclick = () => submit(false, j.form);  // exactly the settings of this preview
  if (j.status !== "done") $("#cur-log-box").open = true;
  if (j.argv.includes("--save_preset")) reloadPresets();
  refreshHistory();
}

async function refreshHistory() {
  const jobs = await api("/api/jobs");
  if (!Array.isArray(jobs)) return;
  $("#history").replaceChildren(...jobs.map((j) => el("li", {},
    el("span", { class: "pill " + j.status }, j.status),
    el("span", { class: "cmd", title: j.command + "\n" + fmtTime(j.created), onclick: () => watch(j.id) },
      j.command.replace(/^protfig /, "")),
    el("button", { class: "btn small", type: "button", title: "Put these settings back in the form",
      onclick: () => { if (j.form && j.form.preset) setForm(j.form); } }, "Reuse"),
    el("button", { class: "btn small", type: "button", title: "Delete this job and its files",
      disabled: j.status === "queued" || j.status === "running",
      onclick: async () => {
        if (!confirm("Delete this job and its files?")) return;
        await api(`/api/jobs/${j.id}`, { method: "DELETE" });
        if (state.current === j.id) $("#current").hidden = true;
        refreshHistory();
      } }, "Delete"))));
  if (!jobs.length) $("#history").append(el("li", { class: "muted" }, "No figures yet."));
}

async function reloadPresets() {
  state.schema.presets = await api("/api/presets");
  renderPresets();
}

// ---------- start ----------

async function init() {
  state.schema = await api("/api/schema");
  if (state.schema.error) { document.body.textContent = state.schema.error; return; }
  for (const g of state.schema.groups) for (const f of g.fields) state.fields[f.dest] = f;
  $("#hostinfo").textContent = `${state.schema.host} · relative paths from ${state.schema.launch_dir}`;
  setForm(loadSavedForm() || state.form);

  $("#input").addEventListener("input", (e) => { state.form.input = e.target.value; changed(); });
  $("#opt-search").addEventListener("input", filterOptions);
  $("#reset-opts").addEventListener("click", () => { state.form.options = {}; renderOptions(); changed(); });
  $("#run-draft").addEventListener("click", () => submit(true));
  $("#run-full").addEventListener("click", () => submit(false));
  $("#cur-cancel").addEventListener("click", () => api(`/api/jobs/${state.current}/cancel`, { method: "POST" }));
  $("#copy-cmd").addEventListener("click", () => navigator.clipboard?.writeText($("#command").textContent));
  $("#file-picker").addEventListener("change", (e) => { if (e.target.files[0]) uploadFile(e.target.files[0]); e.target.value = ""; });
  $("#br-go").addEventListener("click", () => browse($("#br-path").value));
  $("#br-path").addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); browse(e.target.value); } });
  document.addEventListener("click", (e) => {
    const up = e.target.closest("[data-upload]");
    if (up) { uploadTarget = up.dataset.upload; $("#file-picker").click(); }
    const br = e.target.closest("[data-browse]");
    if (br) {
      browseTarget = br.dataset.browse;
      const cur = document.getElementById(browseTarget).value;
      browse(cur.includes("/") ? cur : "");
      $("#browser").showModal();
    }
  });
  refreshHistory();
  const m = location.hash.match(/job=([\w-]+)/);
  if (m) watch(m[1]);
}
init();
