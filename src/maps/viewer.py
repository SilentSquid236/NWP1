"""
The per-run HTML viewer: one self-contained index.html beside the PNGs.

No external libraries, no web server: the product list and hours are
embedded as JSON, so the page works from a copied folder (file://) or through
an SSH tunnel. Left: products grouped as Surface / Upper air / Analysis /
Verification. Top: the forecast-hour strip, play button, and the arrow keys
(left/right = hour, up/down = product).
"""

import json

PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>NWP1 __RUN__</title>
<style>
 body{margin:0;font:14px/1.35 system-ui,Segoe UI,Helvetica,Arial,sans-serif;background:#eef0f3;color:#1d2330}
 header{background:#1d2a44;color:#fff;padding:8px 14px;display:flex;gap:18px;align-items:baseline}
 header b{font-size:16px} header span{opacity:.8;font-size:13px}
 #wrap{display:flex;min-height:calc(100vh - 40px)}
 nav{width:250px;background:#fff;border-right:1px solid #d5d9e0;padding:8px 0;overflow-y:auto}
 nav h3{margin:12px 14px 4px;font-size:12px;text-transform:uppercase;letter-spacing:.06em;color:#6a7285}
 nav a{display:block;padding:5px 14px;color:#1d2330;text-decoration:none;cursor:pointer}
 nav a:hover{background:#eef3fb} nav a.on{background:#1d62d6;color:#fff}
 nav a.off{color:#b0b5bf;cursor:default;background:none}
 main{flex:1;padding:10px 14px}
 #hours{display:flex;flex-wrap:wrap;gap:3px;margin-bottom:8px;align-items:center}
 #hours button{border:1px solid #c6ccd6;background:#fff;border-radius:3px;padding:3px 6px;font-size:12px;cursor:pointer;min-width:34px}
 #hours button.on{background:#1d62d6;color:#fff;border-color:#1d62d6}
 #hours button:disabled{color:#c2c6ce;background:#f4f5f7;cursor:default}
 #play{font-weight:bold;margin-right:8px}
 img{max-width:100%;height:auto;border:1px solid #c6ccd6;background:#fff;display:block}
 #note{color:#6a7285;font-size:12px;margin-top:6px}
</style></head><body>
<header><b>NWP1</b><span>__RUN__</span><span id="status"></span></header>
<div id="wrap"><nav id="menu"></nav>
<main><div id="hours"></div><img id="map" alt="map"><div id="note"></div></main></div>
<script>
const M = __MANIFEST__;
const groups = ["Surface","Upper air","Analysis","Verification"];
let prod = M.default_product, hour = M.products[prod].hours[0], timer = null;
function fmt(h){return "F"+String(h).padStart(3,"0");}
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
  menu(); strip();
}
function step(d){ const hs = M.products[prod].hours; let i = hs.indexOf(hour) + d;
  if (i < 0) i = hs.length - 1; if (i >= hs.length) i = 0; hour = hs[i]; draw(); }
function toggle(){ if (timer){ clearInterval(timer); timer = null; } else { timer = setInterval(()=>step(1), 700); } strip(); }
document.addEventListener("keydown", e => {
  if (e.key === "ArrowRight") step(1); else if (e.key === "ArrowLeft") step(-1);
  else if (e.key === "ArrowUp" || e.key === "ArrowDown"){
    const keys = Object.keys(M.products); let i = keys.indexOf(prod) + (e.key === "ArrowDown" ? 1 : -1);
    i = (i + keys.length) % keys.length; prod = keys[i];
    const hs = M.products[prod].hours; if (!hs.includes(hour)) hour = hs[0]; draw(); e.preventDefault(); }
});
draw();
</script></body></html>
"""


def write_viewer(path, run_label, products, status, note):
    """products: {key: {"group", "label", "hours": [int, ...]}} in menu order."""
    all_hours = sorted({h for p in products.values() for h in p["hours"]})
    default = next((k for k in ("mslp", "an_sfc") if k in products), next(iter(products)))
    manifest = {"products": products, "all_hours": all_hours,
                "default_product": default, "status": status, "note": note}
    html = (PAGE.replace("__RUN__", run_label)
                .replace("__MANIFEST__", json.dumps(manifest)))
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
