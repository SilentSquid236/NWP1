"""
The per-run HTML viewer: one self-contained index.html beside the PNGs.

No external libraries and no web server. The product list is embedded as
JSON, and the per-hour data are JavaScript files (maps/data/hHHH.js,
sHHH.js) loaded with <script> tags, which browsers allow from file://.

  * left: products grouped Surface / Upper air / Analysis / Verification
  * top: forecast hours, play, and the arrow keys (left/right = hour,
    up/down = product)
  * hover: the values under the cursor for the product shown
  * click: the model sounding at that grid point, laid out like SHARPpy
    (Blumberg et al. 2017): skew-T
    log-p, wind barbs, hodograph and a few indices. The model is dry, so
    forecast soundings have no dewpoint and no CAPE; hour 0 (the analysis)
    shows a dewpoint from the analysed humidity.
"""

import json

PAGE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>NWP1 __RUN__</title>
<style>
 body{margin:0;font:14px/1.35 system-ui,Segoe UI,Helvetica,Arial,sans-serif;background:#eef0f3;color:#1d2330}
 header{background:#1d2a44;color:#fff;padding:8px 14px;display:flex;gap:18px;align-items:baseline}
 header b{font-size:16px} header span{opacity:.85;font-size:13px}
 #wrap{display:flex;min-height:calc(100vh - 40px)}
 nav{width:250px;flex:none;background:#fff;border-right:1px solid #d5d9e0;padding:8px 0;overflow-y:auto}
 nav h3{margin:12px 14px 4px;font-size:12px;text-transform:uppercase;letter-spacing:.06em;color:#6a7285}
 nav a{display:block;padding:5px 14px;color:#1d2330;text-decoration:none;cursor:pointer}
 nav a:hover{background:#eef3fb} nav a.on{background:#1d62d6;color:#fff}
 main{flex:1;padding:10px 14px;min-width:0}
 #hours{display:flex;flex-wrap:wrap;gap:3px;margin-bottom:8px;align-items:center}
 #hours button{border:1px solid #c6ccd6;background:#fff;border-radius:3px;padding:3px 6px;font-size:12px;cursor:pointer;min-width:34px}
 #hours button.on{background:#1d62d6;color:#fff;border-color:#1d62d6}
 #hours button:disabled{color:#c2c6ce;background:#f4f5f7;cursor:default}
 #play{font-weight:bold;margin-right:8px}
 #mapbox{position:relative;display:inline-block;max-width:100%}
 #map{max-width:100%;height:auto;border:1px solid #c6ccd6;background:#fff;display:block;cursor:crosshair}
 #tip{position:absolute;pointer-events:none;background:rgba(20,26,40,.92);color:#fff;font:12px/1.35 Consolas,monospace;
      padding:6px 8px;border-radius:4px;white-space:pre;display:none;z-index:5}
 #note{color:#6a7285;font-size:12px;margin-top:6px}
 #snd{position:fixed;inset:0;background:rgba(0,0,0,.55);display:none;z-index:10;align-items:center;justify-content:center}
 #sndbox{background:#000;border:1px solid #555;padding:8px;position:relative}
 #sndclose{position:absolute;right:10px;top:6px;color:#ddd;cursor:pointer;font:18px sans-serif}
</style></head><body>
<header><b>NWP1</b><span>__RUN__</span><span id="status"></span></header>
<div id="wrap"><nav id="menu"></nav>
<main><div id="hours"></div>
<div id="mapbox"><img id="map" alt="map"><div id="tip"></div></div>
<div id="note"></div></main></div>
<div id="snd"><div id="sndbox"><span id="sndclose">&#x2715;</span><canvas id="skew" width="1000" height="640"></canvas></div></div>
<script>
const M = __MANIFEST__;
const groups = ["Surface","Upper air","Analysis","Verification"];
let prod = M.default_product, hour = M.products[prod].hours[0], timer = null;
window.NWPH = window.NWPH || {}; window.NWPS = window.NWPS || {};
const cache = {h:{}, s:{}};
function fmt(h){return "F"+String(h).padStart(3,"0");}
function pad3(h){return String(h).padStart(3,"0");}

/* ---------- menu, hours ---------- */
function menu(){
  const nav = document.getElementById("menu"); nav.innerHTML = "";
  for (const g of groups){
    const keys = Object.keys(M.products).filter(k => M.products[k].group === g);
    if (!keys.length) continue;
    const h = document.createElement("h3"); h.textContent = g; nav.appendChild(h);
    for (const k of keys){
      const a = document.createElement("a"); a.textContent = M.products[k].label;
      a.className = (k === prod ? "on" : "");
      a.onclick = () => { prod = k; const hs = M.products[k].hours;
        if (!hs.includes(hour)) hour = hs.reduce((b,x)=>Math.abs(x-hour)<Math.abs(b-hour)?x:b, hs[0]);
        draw(); };
      nav.appendChild(a);
    }
  }
}
function strip(){
  const div = document.getElementById("hours"); div.innerHTML = "";
  const b = document.createElement("button"); b.id = "play"; b.textContent = timer ? "\u275a\u275a" : "\u25b6";
  b.onclick = toggle; div.appendChild(b);
  const hs = M.products[prod].hours;
  for (const h of M.all_hours){
    const x = document.createElement("button"); x.textContent = fmt(h);
    x.disabled = !hs.includes(h); x.className = (h === hour ? "on" : "");
    x.onclick = () => { hour = h; draw(); }; div.appendChild(x);
  }
}
function draw(){
  document.getElementById("map").src = prod + "_" + fmt(hour).toLowerCase() + ".png";
  document.getElementById("note").textContent = M.note;
  document.getElementById("status").textContent = M.status;
  menu(); strip(); loadData("h", hour, ()=>{});
}
function step(d){ const hs = M.products[prod].hours; let i = hs.indexOf(hour) + d;
  if (i < 0) i = hs.length - 1; if (i >= hs.length) i = 0; hour = hs[i]; draw(); }
function toggle(){ if (timer){ clearInterval(timer); timer = null; } else { timer = setInterval(()=>step(1), 700); } strip(); }
document.addEventListener("keydown", e => {
  if (e.key === "Escape") { document.getElementById("snd").style.display = "none"; return; }
  if (e.key === "ArrowRight") step(1); else if (e.key === "ArrowLeft") step(-1);
  else if (e.key === "ArrowUp" || e.key === "ArrowDown"){
    const keys = Object.keys(M.products); let i = keys.indexOf(prod) + (e.key === "ArrowDown" ? 1 : -1);
    i = (i + keys.length) % keys.length; prod = keys[i];
    const hs = M.products[prod].hours; if (!hs.includes(hour)) hour = hs[0]; draw(); e.preventDefault(); }
});

/* ---------- data files ---------- */
function unpack(o){
  const bin = atob(o.b64), buf = new ArrayBuffer(bin.length), u8 = new Uint8Array(buf);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  const q = new Int16Array(buf), n = q.length / o.names.length, out = {};
  o.names.forEach((name, c) => { out[name] = {q: q.subarray(c*n, (c+1)*n), s: o.scale[c], o: o.offset[c]}; });
  return out;
}
function val(ch, idx){ const v = ch.q[idx]; return v === -32768 ? NaN : ch.o + v * ch.s; }
function loadData(kind, h, done){
  const store = kind === "h" ? window.NWPH : window.NWPS;
  if (cache[kind][h]) return done(cache[kind][h]);
  if (!M.data_hours.includes(h)) return done(null);
  if (store[h]) { cache[kind][h] = {meta: store[h], ch: unpack(store[h])}; return done(cache[kind][h]); }
  const sc = document.createElement("script"); sc.src = "data/" + kind + pad3(h) + ".js";
  sc.onload = () => { cache[kind][h] = {meta: store[h], ch: unpack(store[h])}; done(cache[kind][h]); };
  sc.onerror = () => done(null);
  document.body.appendChild(sc);
}

/* ---------- pixel -> lat/lon -> grid point ---------- */
function toGrid(ev){
  const img = document.getElementById("map"), r = img.getBoundingClientRect(), F = M.frame;
  const sx = img.naturalWidth / r.width, sy = img.naturalHeight / r.height;
  const px = (ev.clientX - r.left) * sx, py = (ev.clientY - r.top) * sy;
  const fx = (px - F.ax_left) / F.ax_w, fy = (py - F.ax_top) / F.ax_h;
  if (fx < 0 || fx > 1 || fy < 0 || fy > 1) return null;
  const x = F.extent[0] + fx * (F.extent[1] - F.extent[0]);
  const y = F.extent[3] - fy * (F.extent[3] - F.extent[2]);
  const P = F.proj, dy = P.rho0 - y, rho = Math.sign(P.n) * Math.hypot(x, dy);
  const th = Math.atan2(x, dy);
  const lon = P.lon0 + th / P.n * 180 / Math.PI;
  const lat = (2 * Math.atan(Math.pow(P.R * P.F / rho, 1 / P.n)) - Math.PI / 2) * 180 / Math.PI;
  const G = M.grid, j = Math.round((lat - G.lat0) / G.dlat), i = Math.round((lon - G.lon0) / G.dlon);
  if (j < 0 || j >= G.ny || i < 0 || i >= G.nx) return null;
  return {lat, lon, j, i, idx: j * G.nx + i, px: ev.clientX - r.left, py: ev.clientY - r.top};
}
function wind(ch, name, idx){
  const u = val(ch[name + "_u"], idx), v = val(ch[name + "_v"], idx);
  if (!isFinite(u)) return "below ground";
  const spd = Math.hypot(u, v); let dir = (Math.atan2(-u, -v) * 180 / Math.PI + 360) % 360;
  return String(Math.round(dir / 5) * 5).padStart(3, "0") + "\u00b0 " + Math.round(spd) + " kt";
}
const img = document.getElementById("map"), tip = document.getElementById("tip");
img.addEventListener("mousemove", ev => {
  const g = toGrid(ev); if (!g) { tip.style.display = "none"; return; }
  const d = cache.h[hour];
  let lines = [g.lat.toFixed(2) + "\u00b0N " + (-g.lon).toFixed(2) + "\u00b0W"];
  if (d){
    const ch = d.ch, want = (M.products[prod].hover || []).concat(["terrain"]);
    for (const name of want){
      const meta = M.hover_meta[name]; if (!meta) continue;
      if (ch[name + "_u"]) lines.push(meta[0] + ": " + wind(ch, name, g.idx));
      else if (ch[name]) { const v = val(ch[name], g.idx);
        lines.push(meta[0] + ": " + (isFinite(v) ? v.toFixed(meta[2]) + " " + meta[1] : "below ground")); }
    }
  }
  lines.push("click: sounding");
  tip.textContent = lines.join("\n"); tip.style.display = "block";
  tip.style.left = (g.px + 14) + "px"; tip.style.top = (g.py + 10) + "px";
});
img.addEventListener("mouseleave", () => { tip.style.display = "none"; });
img.addEventListener("click", ev => {
  const g = toGrid(ev); if (!g) return;
  loadData("h", hour, hd => loadData("s", hour, sd => { if (sd) sounding(g, sd, hd); }));
});
document.getElementById("sndclose").onclick = () => { document.getElementById("snd").style.display = "none"; };
document.getElementById("snd").onclick = e => { if (e.target.id === "snd") e.currentTarget.style.display = "none"; };

/* ---------- the sounding ---------- */
const RD = 287.05, G0 = 9.80665, KAP = 0.2857;
function column(g, sd, hd){
  const ch = sd.ch, m = sd.meta, st = m.stride || 1, nxs = m.nx || M.grid.nx;
  const nys = Math.ceil(M.grid.ny / st);
  const js = Math.min(Math.round(g.j / st), nys - 1), is = Math.min(Math.round(g.i / st), nxs - 1);
  const idx = js * nxs + is, full = (js * st) * M.grid.nx + is * st;
  const terr = hd ? val(hd.ch.terrain, full) : 0;
  g.sj = js * st; g.si = is * st;
  let L = [];
  if (m.kind === "sigma"){
    const ps = val(ch.ps, idx);
    for (let k = m.sigma.length - 1; k >= 0; k--){
      L.push({p: m.p_top + m.sigma[k] * (ps - m.p_top), T: val(ch["T" + k], idx),
              u: val(ch["u" + k], idx), v: val(ch["v" + k], idx), rh: NaN});
    }
    let z = terr + RD * L[0].T / G0 * Math.log(ps / L[0].p);
    L[0].z = z;
    for (let k = 1; k < L.length; k++){ z += RD * 0.5 * (L[k].T + L[k-1].T) / G0 * Math.log(L[k-1].p / L[k].p); L[k].z = z; }
    L.unshift({p: ps, T: L[0].T + 0.0065 * (L[0].z - terr), u: L[0].u, v: L[0].v, rh: NaN, z: terr});
  } else {
    m.levels_hPa.forEach((ph, k) => {
      const z = val(ch["HGT" + k], idx);
      if (z >= terr) L.push({p: ph * 100, T: val(ch["TMP" + k], idx), u: val(ch["UGRD" + k], idx),
                             v: val(ch["VGRD" + k], idx), rh: val(ch["RH" + k], idx), z: z});
    });
    L.sort((a, b) => b.p - a.p);
  }
  return {L, terr};
}
function dewpoint(Tk, rh){ const T = Tk - 273.15, a = 17.625, b = 243.04;
  const gam = Math.log(Math.max(rh, 1) / 100) + a * T / (b + T); return b * gam / (a - gam); }
function interpZ(L, key, zq){
  for (let k = 1; k < L.length; k++) if (L[k].z >= zq){ const a = L[k-1], b = L[k], w = (zq - a.z) / (b.z - a.z); return a[key] + w * (b[key] - a[key]); }
  return NaN;
}
function interpP(L, key, pq){
  for (let k = 1; k < L.length; k++) if (L[k].p <= pq){ const a = L[k-1], b = L[k], w = Math.log(a.p / pq) / Math.log(a.p / b.p); return a[key] + w * (b[key] - a[key]); }
  return NaN;
}
function barb(ctx, x, y, u, v, color){
  const spd = Math.hypot(u, v); if (!isFinite(spd)) return;
  const ang = Math.atan2(-u, -v);           /* direction the wind comes from */
  ctx.save(); ctx.translate(x, y); ctx.rotate(ang); ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 1.3;
  const len = 30; ctx.beginPath(); ctx.moveTo(0, 0); ctx.lineTo(0, -len); ctx.stroke();
  if (spd < 2.5) { ctx.beginPath(); ctx.arc(0, 0, 3, 0, 2 * Math.PI); ctx.stroke(); ctx.restore(); return; }
  let r = Math.round(spd / 5) * 5, pos = -len;
  while (r >= 50){ ctx.beginPath(); ctx.moveTo(0, pos); ctx.lineTo(10, pos + 3); ctx.lineTo(0, pos + 7); ctx.fill(); pos += 9; r -= 50; }
  while (r >= 10){ ctx.beginPath(); ctx.moveTo(0, pos); ctx.lineTo(11, pos - 5); ctx.stroke(); pos += 5; r -= 10; }
  if (r >= 5){ if (pos === -len) pos += 4; ctx.beginPath(); ctx.moveTo(0, pos); ctx.lineTo(6, pos - 3); ctx.stroke(); }
  ctx.restore();
}
function sounding(g, sd, hd){
  const {L, terr} = column(g, sd, hd), c = document.getElementById("skew"), ctx = c.getContext("2d");
  document.getElementById("snd").style.display = "flex";
  ctx.fillStyle = "#000"; ctx.fillRect(0, 0, c.width, c.height);
  const X0 = 55, X1 = 600, Y0 = 20, Y1 = 600, PB = 105000, PT = 20000;
  const yOf = p => Y1 - (Math.log(PB / p) / Math.log(PB / PT)) * (Y1 - Y0);
  const xOf = (Tc, p) => X0 + (Tc + 40) / 90 * (X1 - X0) + (Y1 - yOf(p)) * 0.95;
  ctx.save(); ctx.beginPath(); ctx.rect(X0, Y0, X1 - X0, Y1 - Y0); ctx.clip();
  ctx.lineWidth = 1; ctx.font = "11px sans-serif";
  for (let p = 100000; p >= PT; p -= 10000){ ctx.strokeStyle = "#444"; ctx.beginPath(); ctx.moveTo(X0, yOf(p)); ctx.lineTo(X1, yOf(p)); ctx.stroke(); }
  for (let T = -120; T <= 50; T += 10){ ctx.strokeStyle = T === 0 ? "#3a6fd8" : "#333"; ctx.setLineDash(T === 0 ? [] : [4, 4]);
    ctx.beginPath(); ctx.moveTo(xOf(T, PB), yOf(PB)); ctx.lineTo(xOf(T, PT), yOf(PT)); ctx.stroke(); }
  ctx.setLineDash([]); ctx.strokeStyle = "#5a3a1a";
  for (let th = 250; th <= 450; th += 10){ ctx.beginPath();
    for (let p = PB; p >= PT; p -= 2500){ const T = th * Math.pow(p / 100000, KAP) - 273.15;
      if (p === PB) ctx.moveTo(xOf(T, p), yOf(p)); else ctx.lineTo(xOf(T, p), yOf(p)); } ctx.stroke(); }
  const line = (key, color, w) => { ctx.strokeStyle = color; ctx.lineWidth = w; ctx.beginPath(); let started = false;
    for (const l of L){ const v = key(l); if (!isFinite(v)) continue;
      if (!started){ ctx.moveTo(xOf(v, l.p), yOf(l.p)); started = true; } else ctx.lineTo(xOf(v, l.p), yOf(l.p)); } ctx.stroke(); };
  line(l => l.T - 273.15, "#ff3030", 2.5);
  const hasTd = L.some(l => isFinite(l.rh));
  if (hasTd) line(l => dewpoint(l.T, l.rh), "#30d030", 2.5);
  ctx.restore();
  ctx.fillStyle = "#bbb"; ctx.font = "11px sans-serif";
  for (let p = 100000; p >= PT; p -= 10000) ctx.fillText(String(p / 100), 12, yOf(p) + 4);
  for (let T = -40; T <= 40; T += 10) ctx.fillText(String(T), xOf(T, PB) - 8, Y1 + 14);
  ctx.strokeStyle = "#777"; ctx.strokeRect(X0, Y0, X1 - X0, Y1 - Y0);
  for (const l of L){ if (l.p < PT || !isFinite(l.u)) continue; barb(ctx, 640, yOf(l.p), l.u / 0.514444, l.v / 0.514444, "#ddd"); }
  /* hodograph */
  const HX = 820, HY = 175, HR = 145, KMAX = 80;
  ctx.strokeStyle = "#444"; ctx.lineWidth = 1;
  for (let k = 20; k <= KMAX; k += 20){ ctx.beginPath(); ctx.arc(HX, HY, HR * k / KMAX, 0, 2 * Math.PI); ctx.stroke(); ctx.fillStyle = "#777"; ctx.fillText(k + "", HX + HR * k / KMAX + 2, HY - 2); }
  ctx.beginPath(); ctx.moveTo(HX - HR, HY); ctx.lineTo(HX + HR, HY); ctx.moveTo(HX, HY - HR); ctx.lineTo(HX, HY + HR); ctx.stroke();
  const segs = [[0, 3000, "#ff4040"], [3000, 6000, "#40d040"], [6000, 9000, "#e0e040"], [9000, 12000, "#40d0e0"]];
  for (const [a, b, col] of segs){ ctx.strokeStyle = col; ctx.lineWidth = 2; ctx.beginPath(); let st = false;
    for (let z = a; z <= b; z += 250){ const u = interpZ(L, "u", terr + z) / 0.514444, v = interpZ(L, "v", terr + z) / 0.514444;
      if (!isFinite(u)) continue; const x = HX + u / KMAX * HR, y = HY - v / KMAX * HR;
      if (!st){ ctx.moveTo(x, y); st = true; } else ctx.lineTo(x, y); } ctx.stroke(); }
  /* indices */
  const T3 = interpZ(L, "T", terr + 3000), T0 = L[0].T;
  const lr03 = (T0 - T3) / 3, t700 = interpP(L, "T", 70000), t500 = interpP(L, "T", 50000);
  const z700 = interpP(L, "z", 70000), z500 = interpP(L, "z", 50000), lr75 = (t700 - t500) / ((z500 - z700) / 1000);
  const shear = zt => { const u = interpZ(L, "u", terr + zt) - L[0].u, v = interpZ(L, "v", terr + zt) - L[0].v; return Math.hypot(u, v) / 0.514444; };
  let fzl = NaN; for (let k = 1; k < L.length; k++) if (L[k-1].T >= 273.15 && L[k].T < 273.15){ fzl = L[k-1].z + (L[k-1].T - 273.15) / (L[k-1].T - L[k].T) * (L[k].z - L[k-1].z) - terr; break; }
  const valid = new Date(new Date(M.init + "Z").getTime() + hour * 3600e3);
  const txt = [
    "NWP1 model sounding   " + fmt(hour) + "   valid " + valid.toISOString().slice(0, 13).replace("T", " ") + "Z",
    g.lat.toFixed(2) + "\u00b0N " + (-g.lon).toFixed(2) + "\u00b0W   column (" + g.sj + ", " + g.si + ")   terrain " + Math.round(terr) + " m",
    "Surface pressure " + (L[0].p / 100).toFixed(1) + " mb   T " + (T0 - 273.15).toFixed(1) + " \u00b0C (" + ((T0 - 273.15) * 1.8 + 32).toFixed(0) + " \u00b0F)",
    "", "Lapse rate sfc-3 km   " + lr03.toFixed(1) + " \u00b0C/km",
    "Lapse rate 700-500 mb " + (isFinite(lr75) ? lr75.toFixed(1) + " \u00b0C/km" : "n/a"),
    "Freezing level        " + (isFinite(fzl) ? Math.round(fzl) + " m AGL" : (T0 < 273.15 ? "at the surface" : "above the top")),
    "Bulk shear 0-1 km     " + Math.round(shear(1000)) + " kt",
    "Bulk shear 0-6 km     " + (isFinite(shear(6000)) ? Math.round(shear(6000)) + " kt" : "n/a"),
    "", hasTd ? "Dewpoint (green) from the analysed humidity." : "The model is dry: no dewpoint, no CAPE.",
    "Hodograph: 0-3 km red, 3-6 green, 6-9 yellow, 9-12 cyan (kt).", "Esc or click outside to close."];
  ctx.fillStyle = "#eee"; ctx.font = "13px Consolas, monospace";
  txt.forEach((t, k) => ctx.fillText(t, 675, 350 + k * 19));
}
draw();
</script></body></html>
"""


def write_viewer(path, run_label, products, status, note, extra=None):
    """products: {key: {"group", "label", "hours": [int, ...], "hover": [...]}} in menu order."""
    all_hours = sorted({h for p in products.values() for h in p["hours"]})
    default = next((k for k in ("mslp", "an_sfc") if k in products), next(iter(products)))
    manifest = {"products": products, "all_hours": all_hours, "default_product": default,
                "status": status, "note": note, "data_hours": [], "frame": None,
                "grid": None, "hover_meta": {}, "init": ""}
    manifest.update(extra or {})
    html = (PAGE.replace("__RUN__", run_label)
                .replace("__MANIFEST__", json.dumps(manifest)))
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
