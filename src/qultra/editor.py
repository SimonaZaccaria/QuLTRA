"""
Graphical circuit editor for QuLTRA.

Draw a circuit in your browser, set every parameter of every component,
and get the netlist back in Python as real qultra objects.

    import qultra as qu

    design = qu.open_editor()          # opens the editor, waits for "Send to Python"
    circuit = design.circuit()         # same as qu.QCircuit(design.net, f_min, f_max)

Save and reopen a drawing:

    design = qu.open_editor("my_circuit.json")   # loads the file if it exists, saves it on send
    design = qu.load_design("my_circuit.json")   # rebuilds the netlist without opening the editor

The editor runs on a small local web server (127.0.0.1 only) that stops
as soon as the circuit has been sent. No internet connection is needed.
"""

from __future__ import annotations

import json
import os
import threading
import webbrowser
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

__all__ = ["open_editor", "load_design", "CircuitDesign"]

# Classes the editor can produce and the arguments that are not parameters.
_CLASSES = ("C", "L", "J", "R", "CPW", "CPW_coupler", "new_CPW", "new_CPW_coupler")
# Fallback location of these classes, used only if qultra/__init__.py does not export them.
_MODULES = {
    "new_CPW": "qultra.components.new_cpw",
    "new_CPW_coupler": "qultra.components.new_cpw_coupler",
}
_MULTI = ("CPW_coupler", "new_CPW_coupler")  # components with a list of 4 nodes
_NODE_KEYS = ("name", "class", "node_minus", "node_plus", "nodes")
_FORMAT_VERSION = 1


@dataclass
class CircuitDesign:
    """A circuit drawn in the editor.

    Attributes
    ----------
    net : list
        The qultra components (``qu.C``, ``qu.J``, ``qu.CPW_coupler``, ...),
        ready to be passed to ``qu.QCircuit``.
    f_min, f_max : float or None
        ``f_starting_point`` and ``f_end_point`` set in the editor [GHz].
    names : list of str
        The name of each component in ``net``, in the same order.
    components : list of dict
        The same netlist as plain data (nodes and parameter values).
    script : str
        An equivalent standalone Python script.
    design : dict
        The drawing itself (positions and wires), used to reopen it.
    """

    net: list
    f_min: float | None
    f_max: float | None
    names: list = field(default_factory=list)
    components: list = field(default_factory=list)
    script: str = ""
    design: dict = field(default_factory=dict, repr=False)

    def circuit(self, f_min: float | None = None, f_max: float | None = None):
        """Build ``qu.QCircuit(net, f_min, f_max)``. Arguments override the editor values."""
        from .core import QCircuit

        f1 = self.f_min if f_min is None else f_min
        f2 = self.f_max if f_max is None else f_max
        if f1 is None or f2 is None:
            raise ValueError("Set f_min and f_max, in the editor or as arguments of circuit().")
        return QCircuit(self.net, f1, f2)

    def to_python(self) -> str:
        """Return a standalone Python script that builds this circuit."""
        return self.script

    def save(self, path: str | os.PathLike) -> None:
        """Save the drawing and the netlist to a JSON file."""
        data = {
            "qultra_editor": _FORMAT_VERSION,
            "f_starting_point": self.f_min,
            "f_end_point": self.f_max,
            "components": self.components,
            "design": self.design,
        }
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, ensure_ascii=False)

    def __iter__(self):
        # Allows: net, f_min, f_max = qu.open_editor()
        return iter((self.net, self.f_min, self.f_max))


def _get_class(cls_name: str):
    """Find a qultra component class, also when it is not exported by qultra/__init__.py."""
    import importlib
    import qultra

    cls = getattr(qultra, cls_name, None)
    if cls is None and cls_name in _MODULES:
        cls = getattr(importlib.import_module(_MODULES[cls_name]), cls_name)
    if cls is None:
        raise ValueError(f"qultra has no component {cls_name!r}")
    return cls


def _build_net(components: list) -> tuple[list, list]:
    """Turn the editor's component list into qultra objects."""

    net, names = [], []
    for c in components:
        name = c.get("name", "?")
        cls_name = c.get("class")
        if cls_name not in _CLASSES:
            raise ValueError(f"{name}: unknown component class {cls_name!r}")
        try:
            cls = _get_class(cls_name)
        except ImportError as exc:
            raise ValueError(f"{name}: cannot import {cls_name} ({exc})") from exc
        kwargs = {k: v for k, v in c.items() if k not in _NODE_KEYS}
        missing = [k for k, v in kwargs.items() if v is None]
        if missing:
            raise ValueError(f"{name}: invalid value for {', '.join(missing)}")
        # phi_ext only exists in recent qultra versions: pass it only when it is used.
        if "phi_ext" in kwargs and not isinstance(kwargs["phi_ext"], list) and kwargs["phi_ext"] == 0:
            kwargs.pop("phi_ext")
        if "N" in kwargs:
            kwargs["N"] = int(kwargs["N"])
        try:
            if cls_name in _MULTI:
                obj = cls(list(c["nodes"]), **kwargs)
            else:
                obj = cls(c["node_minus"], c["node_plus"], **kwargs)
        except TypeError as exc:
            if "phi_ext" in kwargs and "phi_ext" in str(exc):
                raise ValueError(
                    f"{name}: this qultra version does not support phi_ext; set it to 0 or update qultra"
                ) from exc
            raise ValueError(f"{name}: {exc}") from exc
        net.append(obj)
        names.append(name)
    return net, names


def _design_from_payload(payload: dict) -> CircuitDesign:
    result = payload.get("result") or {}
    components = result.get("components") or []
    if not components:
        raise ValueError("the circuit has no components")
    net, names = _build_net(components)
    return CircuitDesign(
        net=net,
        f_min=result.get("f_starting_point"),
        f_max=result.get("f_end_point"),
        names=names,
        components=components,
        script=payload.get("script", ""),
        design=payload.get("design") or {},
    )


def load_design(path: str | os.PathLike) -> CircuitDesign:
    """Rebuild a circuit saved by the editor, without opening it."""
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    components = data.get("components") or []
    net, names = _build_net(components)
    return CircuitDesign(
        net=net,
        f_min=data.get("f_starting_point"),
        f_max=data.get("f_end_point"),
        names=names,
        components=components,
        design=data.get("design") or {},
        script=_script_from_components(components, data.get("f_starting_point"), data.get("f_end_point")),
    )


def _script_from_components(components, f_min, f_max) -> str:
    lines, imports = [], []
    for c in components:
        kw = {k: v for k, v in c.items() if k not in _NODE_KEYS and v is not None}
        if "phi_ext" in kw and not isinstance(kw["phi_ext"], list) and kw["phi_ext"] == 0:
            kw.pop("phi_ext")
        args = ", ".join(f"{k}={v!r}" for k, v in kw.items())
        cls = c["class"]
        head = repr(list(c["nodes"])) if cls in _MULTI else f"{c['node_minus']}, {c['node_plus']}"
        call = f"qu.{cls}"
        lines.append(f"    {call}({head}, {args}),  # {c.get('name', '')}")
    body = "\n".join(lines)
    return (
        "import qultra as qu\n" + "".join(imports) + "\nnet = [\n" + body + "\n]\n\n"
        f"f_min = {f_min!r}  # GHz\nf_max = {f_max!r}  # GHz\n\n"
        "circuit = qu.QCircuit(net, f_min, f_max)\n"
    )


def open_editor(
    path: str | os.PathLike | None = None,
    f_min: float | None = None,
    f_max: float | None = None,
    *,
    open_browser: bool = True,
    port: int = 0,
    timeout: float | None = None,
) -> CircuitDesign:
    """Open the graphical editor and wait until the circuit is sent back.

    Parameters
    ----------
    path : str, optional
        JSON file for the drawing. If it exists it is loaded in the editor;
        when you press "Send to Python" the drawing is saved there.
    f_min, f_max : float, optional
        Initial ``f_starting_point`` and ``f_end_point`` [GHz].
    open_browser : bool
        Open the editor in the default browser (otherwise just print the address).
    port : int
        Local port to use; 0 picks a free one.
    timeout : float, optional
        Seconds to wait before giving up. ``None`` waits until the circuit is sent
        (stop with Ctrl+C or the notebook's interrupt button).

    Returns
    -------
    CircuitDesign
        ``design.net`` is the list of qultra components, ``design.circuit()``
        builds the ``QCircuit``. It can also be unpacked:
        ``net, f_min, f_max = qu.open_editor()``.
    """
    initial = None
    if path is not None and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            saved = json.load(fh)
        initial = saved.get("design")
        if f_min is None:
            f_min = saved.get("f_starting_point")
        if f_max is None:
            f_max = saved.get("f_end_point")

    done = threading.Event()
    box: dict = {}
    page = _HTML.encode("utf-8")
    init = json.dumps({
        "qultra_editor": _FORMAT_VERSION,
        "design": initial,
        "f_min": f_min,
        "f_max": f_max,
        "path": os.fspath(path) if path is not None else None,
    }).encode("utf-8")

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):  # keep the console quiet
            pass

        def _send(self, code, body, ctype):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path in ("/", "/index.html"):
                self._send(200, page, "text/html; charset=utf-8")
            elif self.path == "/api/init":
                self._send(200, init, "application/json")
            else:
                self._send(404, b"not found", "text/plain")

        def do_POST(self):
            if self.path != "/api/submit":
                self._send(404, b"not found", "text/plain")
                return
            try:
                length = int(self.headers.get("Content-Length", 0))
                payload = json.loads(self.rfile.read(length).decode("utf-8"))
                design = _design_from_payload(payload)
                if path is not None:
                    design.save(path)
            except Exception as exc:  # report the problem in the browser, keep waiting
                self._send(200, json.dumps({"ok": False, "error": str(exc)}).encode(), "application/json")
                return
            box["design"] = design
            self._send(200, b'{"ok": true}', "application/json")
            done.set()

    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    print(f"QuLTRA editor running at {url}  (waiting for 'Send to Python')")
    if open_browser:
        webbrowser.open(url)
    try:
        if not done.wait(timeout):
            raise TimeoutError("No circuit was sent from the editor.")
    finally:
        server.shutdown()
        server.server_close()

    design = box["design"]
    print(f"Received {len(design.net)} components from the editor.")
    return design


_HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>QuLTRA circuit editor</title>
<style>
:root{
  --paper:#eaf1f5; --grid:#d2dee6; --grid2:#bccdd8; --ink:#1c3550; --ink-soft:#5a7085;
  --panel:#f7fafc; --line:#c9d6df; --accent:#d9480f; --pin:#2b7a78; --code:#0f2236; --code-bg:#e3ebf1; --bad:#c92a2a;
  color-scheme:light;
  box-sizing:border-box;
  padding-top:env(safe-area-inset-top,0px);
  padding-bottom:env(safe-area-inset-bottom,0px);
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --paper:#111c24; --grid:#1a2833; --grid2:#243644; --ink:#cfe1ee; --ink-soft:#8ba3b5;
    --panel:#0d161d; --line:#24394a; --accent:#ff7a3d; --pin:#4fc2bd; --code:#d6e6f2; --code-bg:#15232e; --bad:#ff6b6b;
    color-scheme:dark;
  }
}
:root[data-theme="dark"]{
  --paper:#111c24; --grid:#1a2833; --grid2:#243644; --ink:#cfe1ee; --ink-soft:#8ba3b5;
  --panel:#0d161d; --line:#24394a; --accent:#ff7a3d; --pin:#4fc2bd; --code:#d6e6f2; --code-bg:#15232e; --bad:#ff6b6b;
  color-scheme:dark;
}
*,*::before,*::after{box-sizing:inherit}
html{height:100%;scroll-padding-top:env(safe-area-inset-top,0px)}
body{height:100%;margin:0;background:var(--panel);color:var(--ink);font:15px/1.45 "Barlow",system-ui,-apple-system,"Segoe UI",sans-serif}
.app{height:100%;display:grid;grid-template-columns:176px 1fr 380px;min-height:0}
aside,section{min-height:0}
.palette{border-right:1px solid var(--line);padding:16px 12px;overflow:auto}
.palette h1{font-size:19px;font-weight:600;margin:0 0 4px}
.palette p{margin:0 0 14px;color:var(--ink-soft);font-size:13px}
.part{display:flex;align-items:center;gap:8px;width:100%;padding:6px 8px;margin-bottom:6px;border:1px solid var(--line);border-radius:6px;background:transparent;color:var(--ink);font:inherit;font-size:14px;cursor:grab;text-align:left}
.part:hover{border-color:var(--ink-soft)}
.part code{font:11px "JetBrains Mono",ui-monospace,monospace;color:var(--ink-soft);display:block}
.part:focus-visible,button:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.part svg{width:44px;height:30px;flex:none;stroke:var(--ink);fill:none;stroke-width:2}
.tools{margin-top:14px;display:flex;flex-direction:column;gap:6px}
.btn{font:inherit;font-size:14px;padding:6px 10px;border:1px solid var(--line);border-radius:6px;background:transparent;color:var(--ink);cursor:pointer}
.btn:hover{border-color:var(--ink-soft)}
.btn.warn{color:var(--accent)}
.btn.primary{background:var(--ink);color:var(--panel);border-color:var(--ink);font-weight:600}
.btn.primary:hover{opacity:.9}
.sendmsg{font-size:13px;margin:8px 0 0;color:var(--ink-soft)}
.sendmsg.bad{color:var(--bad)}
.stage{position:relative;background:var(--paper);overflow:hidden}
#cv{width:100%;height:100%;display:block;touch-action:none;user-select:none;-webkit-user-select:none}
.hint{position:absolute;left:12px;bottom:10px;font-size:12.5px;color:var(--ink-soft);pointer-events:none;max-width:calc(100% - 24px)}
.side{border-left:1px solid var(--line);display:flex;flex-direction:column;overflow:auto}
.props{padding:16px;border-bottom:1px solid var(--line)}
.props h2,.out h2{font-size:15px;font-weight:600;margin:0 0 10px}
.props h2 code{font:13px "JetBrains Mono",ui-monospace,monospace;color:var(--ink-soft);font-weight:400;margin-left:6px}
.props .empty{color:var(--ink-soft);font-size:14px;margin:0}
.field{display:flex;flex-direction:column;gap:3px;margin-bottom:10px}
.field label{font-size:13px;color:var(--ink-soft)}
.field label code{font:12.5px "JetBrains Mono",ui-monospace,monospace;color:var(--ink)}
.field input{font:inherit;padding:6px 8px;border:1px solid var(--line);border-radius:6px;background:var(--paper);color:var(--ink);width:100%}
.field input.bad{border-color:var(--bad)}
.field small{font-size:12.5px;color:var(--ink-soft)}
.field small.bad{color:var(--bad)}
.conn{font:12.5px/1.6 "JetBrains Mono",ui-monospace,monospace;color:var(--ink-soft);margin:0 0 10px;padding:8px 10px;background:var(--code-bg);border-radius:6px}
.row{display:flex;gap:8px;flex-wrap:wrap}
.out{padding:16px;flex:1;display:flex;flex-direction:column;min-height:320px}
.freq{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-bottom:6px}
.tabs{display:flex;gap:4px;margin:4px 0 8px}
.tab{font:inherit;font-size:13.5px;padding:4px 10px;border:1px solid transparent;border-radius:6px;background:transparent;color:var(--ink-soft);cursor:pointer}
.tab[aria-selected="true"]{border-color:var(--line);color:var(--ink);background:var(--code-bg)}
pre{flex:1;margin:0;padding:12px;background:var(--code-bg);color:var(--code);border-radius:6px;font:12.5px/1.5 "JetBrains Mono",ui-monospace,Menlo,Consolas,monospace;overflow:auto;white-space:pre;min-height:180px}
.warns{margin:8px 0 0;padding:0;list-style:none;font-size:13px;color:var(--bad)}
.warns li{margin-top:2px}
.outbar{display:flex;justify-content:space-between;align-items:center;margin-top:8px;gap:8px}
.status{font-size:13px;color:var(--ink-soft)}
/* svg */
.wire{stroke:var(--ink);stroke-width:2;fill:none}
.wire.sel{stroke:var(--accent);stroke-width:3}
.hit{stroke:transparent;stroke-width:14;fill:none;cursor:pointer}
.body{stroke:var(--ink);stroke-width:2;fill:none}
.comp.sel .body{stroke:var(--accent)}
.comp.err .body{stroke:var(--bad)}
.box{fill:transparent;cursor:move}
.pin{fill:var(--paper);stroke:var(--pin);stroke-width:2;cursor:crosshair}
.pin:hover,.pin.on{fill:var(--pin)}
.lbl{fill:var(--ink);font:600 12px "Barlow",system-ui,sans-serif;pointer-events:none}
.val{fill:var(--ink-soft);font:11.5px "Barlow",system-ui,sans-serif;pointer-events:none}
.netl{fill:var(--pin);font:10.5px "JetBrains Mono",ui-monospace,monospace;pointer-events:none}
.dot{fill:var(--ink)}
#rubber{stroke:var(--pin);stroke-width:2;stroke-dasharray:5 4}
@media (max-width:900px){
  .app{grid-template-columns:1fr;grid-template-rows:auto 62vh auto;height:auto}
  .palette{border-right:0;border-bottom:1px solid var(--line)}
  .parts{display:flex;flex-wrap:wrap;gap:6px}
  .part{width:auto;margin:0}
  .tools{flex-direction:row;flex-wrap:wrap}
  .side{border-left:0;border-top:1px solid var(--line)}
}
</style>
</head>
<body>
<div class="app">
  <aside class="palette">
    <h1>QuLTRA circuit</h1>
    <p>Drag or tap a component to add it.</p>
    <div class="parts" id="parts"></div>
    <div class="tools">
      <button class="btn" id="btnFit">Fit to view</button>
      <button class="btn" id="btnExample">Load example</button>
      <button class="btn warn" id="btnClear">Clear schematic</button>
    </div>
  </aside>
  <section class="stage">
    <svg id="cv" xmlns="http://www.w3.org/2000/svg">
      <defs>
        <pattern id="g1" width="20" height="20" patternUnits="userSpaceOnUse"><path d="M20 0H0V20" fill="none" stroke="var(--grid)" stroke-width="1"/></pattern>
        <pattern id="g5" width="100" height="100" patternUnits="userSpaceOnUse"><rect width="100" height="100" fill="url(#g1)"/><path d="M100 0H0V100" fill="none" stroke="var(--grid2)" stroke-width="1"/></pattern>
      </defs>
      <rect x="-20000" y="-20000" width="40000" height="40000" fill="url(#g5)"/>
      <g id="wires"></g>
      <g id="comps"></g>
      <line id="rubber" visibility="hidden"/>
    </svg>
    <div class="hint">Tap a pin, then another pin, to connect them. The numbers next to the pins are the QuLTRA nodes. Drag the background to pan, scroll to zoom. R rotates, Delete removes.</div>
  </section>
  <aside class="side">
    <div class="props" id="props"></div>
    <div class="out">
      <h2>Analysis</h2>
      <div class="freq">
        <div class="field"><label for="fMin"><code>f_starting_point</code> [GHz]</label><input id="fMin" autocomplete="off"></div>
        <div class="field"><label for="fMax"><code>f_end_point</code> [GHz]</label><input id="fMax" autocomplete="off"></div>
      </div>
      <div class="tabs" role="tablist" id="tabs">
        <button class="tab" role="tab" data-fmt="script" aria-selected="true">Script</button>
        <button class="tab" role="tab" data-fmt="list" aria-selected="false">Netlist only</button>
        <button class="tab" role="tab" data-fmt="json" aria-selected="false">JSON</button>
      </div>
      <pre id="code"></pre>
      <ul class="warns" id="warns"></ul>
      <p class="sendmsg" id="sendMsg" hidden></p>
      <div class="outbar">
        <span class="status" id="status"></span>
        <span class="row"><button class="btn" id="btnCopy">Copy</button><button class="btn primary" id="btnSend" hidden>Send to Python</button></span>
      </div>
    </div>
  </aside>
</div>

<script>
(() => {
const GRID = 20;
const TWO = [[-40,0],[40,0]];
/* Parameters taken from the qultra constructors (src/qultra/components) */
const TYPES = {
  C: {label:"Capacitor", cls:"C", prefix:"C", pins:TWO, hh:22,
      d:"M-40 0H-5M-5 -14V14M5 -14V14M5 0H40",
      params:[{k:"C_value", u:"F", d:"90e-15", t:"pos", h:"Capacitance"}]},
  L: {label:"Inductor", cls:"L", prefix:"L", pins:TWO, hh:22,
      d:"M-40 0H-24a6 6 0 0 1 12 0a6 6 0 0 1 12 0a6 6 0 0 1 12 0a6 6 0 0 1 12 0H40",
      params:[{k:"L_value", u:"H", d:"10e-9", t:"pos", h:"Inductance"},
              {k:"phi_ext", u:"rad", d:"0", t:"num", opt:true, h:"External flux as reduced phase (needs QuLTRA from GitHub; omitted when 0)"}]},
  J: {label:"Josephson junction", cls:"J", prefix:"J", pins:TWO, hh:22,
      d:"M-40 0H-10M10 0H40M-10 -10L10 10M-10 10L10 -10",
      params:[{k:"J_value", u:"H", d:"8e-9", t:"pos", h:"Total linear inductance (of the whole array if N > 1)"},
              {k:"N", u:"", d:"1", t:"int", h:"Number of junctions in series"},
              {k:"phi_ext", u:"rad", d:"0", t:"num", opt:true, h:"External flux as reduced phase (needs QuLTRA from GitHub; omitted when 0)"}]},
  R: {label:"Resistor", cls:"R", prefix:"R", pins:TWO, hh:22,
      d:"M-40 0H-24L-20 -8L-12 8L-4 -8L4 8L12 -8L20 8L24 0H40",
      params:[{k:"R_value", u:"Ω", d:"50", t:"pos", h:"Resistance"}]},
  CPW: {label:"CPW line", cls:"CPW", prefix:"T", pins:TWO, hh:22,
      d:"M-40 0H-28M28 0H40M-28 -8H28V8H-28ZM-28 0H28",
      params:[{k:"l", u:"m", d:"4.5e-3", t:"pos", h:"Line length"},
              {k:"Z0", u:"Ω", d:"50", t:"pos", h:"Characteristic impedance"},
              {k:"phi_ext", u:"rad", d:"", t:"pair", empty:true, h:"External phase at the two ends, 2 values (empty = None)"}]},
  CPW_coupler: {label:"CPW coupler", cls:"CPW_coupler", multi:true, prefix:"K", pins:[[-40,-20],[40,-20],[-40,20],[40,20]], hh:34,
      d:"M-40 -20H40M-40 20H40M-28 -12H28V12H-28Z",
      params:[{k:"gap", u:"µm", d:"10, 10, 10, 10", t:"list", h:"Gaps between strips, ground to ground"},
              {k:"cpw", u:"µm", d:"15, 10, 15", t:"list", h:"Strip widths (3 values = central ground plane)"},
              {k:"l", u:"m", d:"0.7e-3", t:"pos", h:"Coupler length"}]},
  new_CPW: {label:"FlipChip CPW", cls:"new_CPW", prefix:"W", pins:TWO, hh:22,
      d:"M-40 0H-28M28 0H40M-28 -8H28V8H-28ZM-28 -3H28M-28 3H28",
      params:[{k:"l", u:"m", d:"4.5e-3", t:"pos", h:"Line length"},
              {k:"w", u:"µm", d:"10", t:"pos", h:"Width of the center conductor"},
              {k:"g", u:"µm", d:"6", t:"pos", h:"Gap to the ground planes"},
              {k:"hs", u:"µm", d:"", t:"pos", empty:true, h:"Chip-to-chip distance in flip-chip (empty = single chip)"},
              {k:"phi_ext", u:"rad", d:"", t:"pair", empty:true, h:"External phase at the two ends, 2 values (empty = None)"}]},
  new_CPW_coupler: {label:"FlipChip Coupler", cls:"new_CPW_coupler", multi:true, prefix:"X", pins:[[-40,-20],[40,-20],[-40,20],[40,20]], hh:34,
      d:"M-40 -20H40M-40 20H40M-28 -12H28V12H-28ZM-28 0H28",
      params:[{k:"gap", u:"µm", d:"10, 10, 10, 10", t:"list", h:"Gaps between strips, ground to ground"},
              {k:"cpw", u:"µm", d:"15, 10, 15", t:"list", h:"Strip widths (3 values = central ground plane)"},
              {k:"l", u:"m", d:"0.7e-3", t:"pos", h:"Coupler length"},
              {k:"hs", u:"µm", d:"", t:"pos", empty:true, h:"Chip-to-chip distance in flip-chip (empty = single chip)"}]},
  GND: {label:"Ground (node 0)", cls:null, prefix:"GND", pins:[[0,-20]], hh:14,
      d:"M0 -20V0M-14 0H14M-9 5H9M-4 10H4", params:[]},
};

const $ = id => document.getElementById(id);
const svg = $("cv"), gW = $("wires"), gC = $("comps"), rubber = $("rubber");
const KEY = "editor-qultra";
let state, sel = null, pending = null, drag = null, pan = null, fmt = "script";
let view = {x:-300, y:-200, w:600, h:400};

/* ---------- numbers ---------- */
const SI = {f:1e-15, p:1e-12, n:1e-9, u:1e-6, "µ":1e-6, m:1e-3, k:1e3, M:1e6, G:1e9};
function parseNum(s){
  const m = String(s).trim().replace(",", ".").match(/^([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([fpnuµmkMG])?\s*[A-Za-zΩ]*$/);
  if(!m) return NaN;
  return parseFloat(m[1]) * (m[2] ? SI[m[2]] : 1);
}
const pyNum = v => String(Number(v.toPrecision(12)));
function parseList(s){
  const parts = String(s).split(/[\s,;]+/).filter(Boolean);
  if(!parts.length) return null;
  const v = parts.map(Number);
  return v.every(x => isFinite(x) && x > 0) ? v : null;
}
/* returns {ok, value, msg} for one parameter */
function checkParam(p, raw){
  if(p.empty && String(raw ?? "").trim() === "") return {ok:true, value:null};
  if(p.t === "pair"){
    const parts = String(raw).split(/[\s,;]+/).filter(Boolean).map(parseNum);
    return parts.length === 2 && parts.every(isFinite) ? {ok:true, value:parts} : {ok:false, msg:"Enter exactly 2 numbers separated by a comma, or leave empty"};
  }
  if(p.t === "list"){
    const v = parseList(raw);
    return v ? {ok:true, value:v} : {ok:false, msg:"Enter positive numbers separated by commas"};
  }
  if(p.t === "int"){
    const n = Number(String(raw).trim());
    return Number.isInteger(n) && n >= 1 ? {ok:true, value:n} : {ok:false, msg:"Must be an integer of at least 1"};
  }
  const v = parseNum(raw);
  if(!isFinite(v)) return {ok:false, msg:"Invalid value. Examples: 90e-15, 90f, 4.5m"};
  if(p.t === "pos" && v <= 0) return {ok:false, msg:"Must be greater than zero"};
  return {ok:true, value:v};
}
function checkComp(c){
  const T = TYPES[c.type], res = {}, errs = [];
  T.params.forEach(p => { const r = checkParam(p, c.params[p.k]); res[p.k] = r; if(!r.ok) errs.push(`${p.k}: ${r.msg}`); });
  if(T.multi && res.gap.ok && res.cpw.ok){
    if(res.gap.value.length !== res.cpw.value.length + 1){
      const m = `gap needs ${res.cpw.value.length + 1} values (one more than cpw), it has ${res.gap.value.length}`;
      res.gap = {ok:false, msg:m}; errs.push(m);
    }
    if(res.cpw.value.length < 2){ res.cpw = {ok:false, msg:"At least 2 widths are needed"}; errs.push("cpw: at least 2 widths are needed"); }
  }
  return {res, errs};
}
function short(c){
  const p = TYPES[c.params ? c.type : c.type].params[0]; if(!p) return "";
  const raw = c.params[p.k];
  return p.t === "list" ? `l = ${c.params.l} m` : `${raw}${p.u && !/[A-Za-zΩ]$/.test(String(raw).trim()) ? " " + p.u : ""}`;
}

/* ---------- storage ---------- */
function save(){ try{ localStorage.setItem(KEY, JSON.stringify(state)); }catch(e){} }
function load(){ try{ const s = JSON.parse(localStorage.getItem(KEY)); if(s && Array.isArray(s.comps)) return s; }catch(e){} return null; }
function example(){
  /* The lossy circuit from the QuLTRA Basics tutorial */
  const P = (o) => o;
  return {seq:20, fmin:"3", fmax:"12", comps:[
    {id:"c1", type:"C", x:-320, y:40, rot:270, name:"C1", params:P({C_value:"90e-15"})},
    {id:"c2", type:"J", x:-240, y:40, rot:270, name:"J1", params:P({J_value:"8e-9", N:"1", phi_ext:"0"})},
    {id:"c3", type:"C", x:-140, y:0, rot:0, name:"C2", params:P({C_value:"5e-15"})},
    {id:"c4", type:"CPW", x:0, y:0, rot:180, name:"T1", params:P({l:"3.8e-3", Z0:"50"})},
    {id:"c5", type:"CPW_coupler", x:140, y:20, rot:0, name:"K1", params:P({gap:"10, 10, 10, 10", cpw:"15, 10, 15", l:"0.7e-3"})},
    {id:"c6", type:"R", x:100, y:120, rot:270, name:"R1", params:P({R_value:"50"})},
    {id:"c7", type:"R", x:180, y:120, rot:270, name:"R2", params:P({R_value:"50"})},
    {id:"g1", type:"GND", x:-320, y:100, rot:0, name:"GND", params:{}},
    {id:"g2", type:"GND", x:-240, y:100, rot:0, name:"GND", params:{}},
    {id:"g3", type:"GND", x:100, y:180, rot:0, name:"GND", params:{}},
    {id:"g4", type:"GND", x:180, y:180, rot:0, name:"GND", params:{}},
    {id:"g5", type:"GND", x:240, y:20, rot:0, name:"GND", params:{}},
  ], wires:[
    {id:"w1", a:{c:"c1",p:0}, b:{c:"g1",p:0}},
    {id:"w2", a:{c:"c2",p:0}, b:{c:"g2",p:0}},
    {id:"w3", a:{c:"c1",p:1}, b:{c:"c2",p:1}},
    {id:"w4", a:{c:"c2",p:1}, b:{c:"c3",p:0}},
    {id:"w5", a:{c:"c3",p:1}, b:{c:"c4",p:1}},
    {id:"w6", a:{c:"c4",p:0}, b:{c:"c5",p:0}},
    {id:"w7", a:{c:"c5",p:1}, b:{c:"g5",p:0}},
    {id:"w8", a:{c:"c5",p:2}, b:{c:"c6",p:1}},
    {id:"w9", a:{c:"c5",p:3}, b:{c:"c7",p:1}},
    {id:"w10", a:{c:"c6",p:0}, b:{c:"g3",p:0}},
    {id:"w11", a:{c:"c7",p:0}, b:{c:"g4",p:0}},
  ]};
}
state = load() || example();
if(!state.seq) state.seq = 1;
if(state.fmin === undefined) state.fmin = "3";
if(state.fmax === undefined) state.fmax = "12";

/* ---------- geometry ---------- */
const snap = v => Math.round(v / GRID) * GRID;
const comp = id => state.comps.find(c => c.id === id);
function rotPt(c, px, py){
  const r = c.rot * Math.PI / 180, cs = Math.round(Math.cos(r)), sn = Math.round(Math.sin(r));
  return {x: c.x + px*cs - py*sn, y: c.y + px*sn + py*cs};
}
const pinPos = (c, i) => { const [px,py] = TYPES[c.type].pins[i]; return rotPt(c, px, py); };
function toSvg(e){
  const p = svg.createSVGPoint(); p.x = e.clientX; p.y = e.clientY;
  return p.matrixTransform(svg.getScreenCTM().inverse());
}
function setView(){ svg.setAttribute("viewBox", `${view.x} ${view.y} ${view.w} ${view.h}`); }
function fitAspect(){
  const w = svg.clientWidth || 1, h = svg.clientHeight || 1;
  const cx = view.x + view.w/2, cy = view.y + view.h/2;
  view.w = view.h * (w/h); view.x = cx - view.w/2; view.y = cy - view.h/2; setView();
}
function fitAll(){
  if(!state.comps.length){ view = {x:-300,y:-200,w:600,h:400}; fitAspect(); return; }
  let x0=Infinity,y0=Infinity,x1=-Infinity,y1=-Infinity;
  state.comps.forEach(c => { x0=Math.min(x0,c.x-70); y0=Math.min(y0,c.y-70); x1=Math.max(x1,c.x+70); y1=Math.max(y1,c.y+70); });
  const w = svg.clientWidth || 1, h = svg.clientHeight || 1, ar = w/h;
  let vw = x1-x0, vh = y1-y0;
  if(vw/vh > ar) vh = vw/ar; else vw = vh*ar;
  vw = Math.max(vw, 300); vh = vw/ar;
  view = {x:(x0+x1)/2 - vw/2, y:(y0+y1)/2 - vh/2, w:vw, h:vh}; setView();
}

/* ---------- nodes ---------- */
function computeNets(){
  const parent = {}, key = (c,p) => c + ":" + p;
  const find = k => { while(parent[k] !== k){ parent[k] = parent[parent[k]]; k = parent[k]; } return k; };
  const union = (a,b) => { a = find(a); b = find(b); if(a !== b) parent[b] = a; };
  state.comps.forEach(c => TYPES[c.type].pins.forEach((_,i) => parent[key(c.id,i)] = key(c.id,i)));
  state.wires.forEach(w => { const a = key(w.a.c,w.a.p), b = key(w.b.c,w.b.p); if(parent[a] && parent[b]) union(a, b); });
  const grounds = state.comps.filter(c => c.type === "GND").map(c => key(c.id,0));
  for(let i = 1; i < grounds.length; i++) union(grounds[0], grounds[i]);
  const gRoot = grounds.length ? find(grounds[0]) : null;
  const names = {}, pinNet = {}, count = {}; let n = 1;
  state.comps.forEach(c => { if(c.type === "GND") return;
    TYPES[c.type].pins.forEach((_,i) => {
      const r = find(key(c.id,i));
      if(!(r in names)) names[r] = r === gRoot ? 0 : n++;
      pinNet[key(c.id,i)] = names[r];
      count[names[r]] = (count[names[r]] || 0) + 1;
    });
  });
  // ground pins
  grounds.forEach(g => pinNet[g] = 0);
  // wires per pin (to find unconnected pins)
  const linked = {};
  state.wires.forEach(w => { linked[key(w.a.c,w.a.p)] = 1; linked[key(w.b.c,w.b.p)] = 1; });
  return {pinNet, count, linked, nNodes:n};
}

/* ---------- output ---------- */
function buildOutput(){
  const {pinNet, count, linked, nNodes} = computeNets();
  const parts = state.comps.filter(c => c.type !== "GND");
  const warns = [], bad = new Set();
  const items = [];
  parts.forEach(c => {
    const T = TYPES[c.type], {res, errs} = checkComp(c);
    const nodes = T.pins.map((_,i) => pinNet[c.id + ":" + i]);
    errs.forEach(e => warns.push(`${c.name}: ${e}`));
    if(errs.length) bad.add(c.id);
    T.pins.forEach((_,i) => { if(!linked[c.id + ":" + i]) warns.push(`${c.name}: pin ${i+1} is not connected`); });
    if(nodes.length === 2 && nodes[0] === nodes[1]){ warns.push(`${c.name}: both pins are on the same node (QuLTRA rejects this as a short)`); bad.add(c.id); }
    if(nodes.length === 4 && nodes.every(n => n === nodes[0])){ warns.push(`${c.name}: all pins are on the same node`); bad.add(c.id); }
    items.push({c, T, res, nodes});
  });
  if(parts.length && !state.comps.some(c => c.type === "GND")) warns.push("No ground: QuLTRA needs node 0.");
  const fmin = parseNum(state.fmin), fmax = parseNum(state.fmax);
  if(!isFinite(fmin) || !isFinite(fmax) || fmin >= fmax) warns.push("Invalid frequency range: f_starting_point must be lower than f_end_point.");

  const val = (p, r, raw) => !r.ok ? "None" : (p.t === "list" || p.t === "pair") ? `[${r.value.map(pyNum).join(", ")}]` : p.t === "int" ? String(r.value) : pyNum(r.value);
  const pyLine = ({c, T, res, nodes}) => {
    const args = T.params.filter(p => !(res[p.k].ok && (res[p.k].value === null || (p.opt && res[p.k].value === 0)))).map(p => `${p.k}=${val(p, res[p.k])}`);
    const head = T.multi ? `[${nodes.join(", ")}]` : nodes.join(", ");
    return `${T.mod ? "" : "qu."}${T.cls}(${head}, ${args.join(", ")})`;
  };
  const width = Math.max(0, ...items.map(it => pyLine(it).length + 1));
  let list = "net = [\n";
  items.forEach(it => { const l = pyLine(it) + ","; list += `    ${l.padEnd(width + 2)}# ${it.c.name}\n`; });
  list += "]\n";

  const f1 = isFinite(fmin) ? pyNum(fmin) : "None", f2 = isFinite(fmax) ? pyNum(fmax) : "None";
  const mods = [...new Set(items.filter(it => it.T.mod).map(it => `from ${it.T.mod} import ${it.T.cls}\n`))].join("");
  const script = `import qultra as qu\n${mods}\n${list}\nf_min = ${f1}  # GHz\nf_max = ${f2}  # GHz\n\ncircuit = qu.QCircuit(net, f_min, f_max)\n\nprint("Mode frequencies (GHz):", circuit.mode_frequencies())\ncircuit.show_all()\n`;

  const json = JSON.stringify({
    f_starting_point: isFinite(fmin) ? fmin : null, f_end_point: isFinite(fmax) ? fmax : null,
    components: items.map(({c, T, res, nodes}) => {
      const o = {name:c.name, class:T.cls};
      if(T.multi) o.nodes = nodes; else { o.node_minus = nodes[0]; o.node_plus = nodes[1]; }
      T.params.forEach(p => { const r = res[p.k]; if(r.ok && r.value === null) return; o[p.k] = r.ok ? r.value : null; });
      return o;
    })
  }, null, 2);

  return {script, list, json, warns, bad, pinNet, nParts:parts.length, nNodes:parts.length ? nNodes : 0};
}

/* ---------- drawing ---------- */
const esc = s => String(s).replace(/[&<>"]/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[ch]));
let lastOut = null;

function renderCanvas(){
  const out = lastOut = buildOutput();
  const counts = {};
  state.wires.forEach(w => [w.a,w.b].forEach(e => { const k = e.c+":"+e.p; counts[k] = (counts[k]||0)+1; }));
  let wh = "";
  state.wires.forEach(w => {
    const A = comp(w.a.c), B = comp(w.b.c); if(!A || !B) return;
    const a = pinPos(A, w.a.p), b = pinPos(B, w.b.p);
    const d = `M${a.x} ${a.y}L${b.x} ${a.y}L${b.x} ${b.y}`;
    const isSel = sel && sel.kind === "wire" && sel.id === w.id;
    wh += `<g data-wire="${w.id}"><path class="hit" d="${d}"/><path class="wire${isSel?" sel":""}" d="${d}"/></g>`;
  });
  Object.keys(counts).forEach(k => { if(counts[k] > 1){ const [c,p] = k.split(":"); const C = comp(c); if(C){ const q = pinPos(C,+p); wh += `<circle class="dot" cx="${q.x}" cy="${q.y}" r="3.5"/>`; } } });
  gW.innerHTML = wh;

  let ch = "";
  state.comps.forEach(c => {
    const T = TYPES[c.type], isSel = sel && sel.kind === "comp" && sel.id === c.id;
    const bh = T.multi ? 28 : 18;
    ch += `<g class="comp${isSel?" sel":""}${out.bad.has(c.id)?" err":""}" data-comp="${c.id}">`;
    ch += `<g transform="translate(${c.x} ${c.y}) rotate(${c.rot})"><rect class="box" x="-32" y="-${bh}" width="64" height="${bh*2}"/><path class="body" d="${T.d}"/></g>`;
    if(c.type !== "GND"){
      const vert = c.rot % 180 !== 0;
      const lx = vert ? c.x + T.hh + 4 : c.x, ly = vert ? c.y - 2 : c.y - T.hh - 4, anchor = vert ? "start" : "middle";
      const cpl = T.multi && !vert;
      const ny = cpl ? c.y - 50 : ly, vy = cpl ? c.y - 36 : vert ? ly + 14 : c.y + T.hh + 12;
      ch += `<text class="lbl" x="${lx}" y="${ny}" text-anchor="${anchor}">${esc(c.name)}</text>`;
      ch += `<text class="val" x="${lx}" y="${vy}" text-anchor="${anchor}">${esc(short(c))}</text>`;
    }
    T.pins.forEach((pp,i) => {
      const q = pinPos(c,i), on = pending && pending.c === c.id && pending.p === i;
      ch += `<circle class="pin${on?" on":""}" data-pin="${i}" data-c="${c.id}" cx="${q.x}" cy="${q.y}" r="5"/>`;
      if(c.type !== "GND"){
        const o = rotPt({x:0,y:0,rot:c.rot}, pp[0] === 0 ? 8 : Math.sign(pp[0])*8, -9);
        ch += `<text class="netl" x="${q.x + o.x}" y="${q.y + o.y + 3}" text-anchor="middle">${out.pinNet[c.id+":"+i]}</text>`;
      }
    });
    ch += `</g>`;
  });
  gC.innerHTML = ch;
  renderOutput();
}

function renderOutput(){
  const out = lastOut || buildOutput();
  $("code").textContent = out[fmt];
  $("warns").innerHTML = out.warns.map(w => `<li>${esc(w)}</li>`).join("");
  $("status").textContent = `${out.nParts} components, ${out.nNodes} nodes (including 0)`;
}

function connText(c){
  const T = TYPES[c.type], out = lastOut || buildOutput();
  const n = T.pins.map((_,i) => out.pinNet[c.id + ":" + i]);
  if(T.multi) return `nodes = [${n.join(", ")}]<br>1──2 on top, 3──4 at the bottom`;
  return `node_minus = ${n[0]}<br>node_plus  = ${n[1]}`;
}
function renderProps(){
  const p = $("props");
  if(sel && sel.kind === "comp"){
    const c = comp(sel.id); if(!c){ sel = null; return renderProps(); }
    const T = TYPES[c.type];
    if(c.type === "GND"){
      p.innerHTML = `<h2>Ground</h2><p class="empty">Every pin connected here goes to node 0. You can use several ground symbols.</p><div class="row" style="margin-top:10px"><button class="btn warn" id="pDel">Delete</button></div>`;
      $("pDel").addEventListener("click", deleteSel);
      return;
    }
    const {res} = checkComp(c);
    let h = `<h2>${T.label}<code>${T.mod ? "" : "qu."}${T.cls}</code></h2>`;
    h += `<div class="field"><label for="pName">Name (only used to identify it in the script)</label><input id="pName" value="${esc(c.name)}" autocomplete="off"></div>`;
    T.params.forEach(pr => {
      const r = res[pr.k];
      h += `<div class="field"><label for="pp_${pr.k}"><code>${pr.k}</code>${pr.u ? ` [${pr.u}]` : ""}</label>
        <input id="pp_${pr.k}" data-k="${pr.k}" value="${esc(c.params[pr.k])}" autocomplete="off" class="${r.ok ? "" : "bad"}">
        <small id="ph_${pr.k}" class="${r.ok ? "" : "bad"}">${esc(r.ok ? hintFor(pr, r) : r.msg)}</small></div>`;
    });
    h += `<p class="conn" id="pConn">${connText(c)}</p>`;
    h += `<div class="row"><button class="btn" id="pRot">Rotate 90°</button><button class="btn warn" id="pDel">Delete</button></div>`;
    p.innerHTML = h;
    $("pName").addEventListener("input", e => { c.name = e.target.value.trim() || T.prefix; renderCanvas(); save(); });
    p.querySelectorAll("input[data-k]").forEach(inp => inp.addEventListener("input", e => {
      c.params[e.target.dataset.k] = e.target.value;
      renderCanvas(); save();
      const {res} = checkComp(c);
      T.params.forEach(pr => {
        const r = res[pr.k], inEl = $("pp_" + pr.k), hEl = $("ph_" + pr.k);
        inEl.classList.toggle("bad", !r.ok); hEl.classList.toggle("bad", !r.ok);
        hEl.textContent = r.ok ? hintFor(pr, r) : r.msg;
      });
    }));
    $("pRot").addEventListener("click", rotateSel);
    $("pDel").addEventListener("click", deleteSel);
  } else if(sel && sel.kind === "wire"){
    p.innerHTML = `<h2>Wire</h2><div class="row"><button class="btn warn" id="pDel">Delete wire</button></div>`;
    $("pDel").addEventListener("click", deleteSel);
  } else {
    p.innerHTML = `<h2>Parameters</h2><p class="empty">Select a component to set all of its parameters. Enter values in scientific notation (90e-15) or with SI prefixes (90f, 8n, 4.5m).</p>`;
  }
}
function hintFor(pr, r){
  if(r.value === null) return `${pr.h}: not set`;
  if(pr.t === "list") return `${pr.h}: ${r.value.length} values`;
  if(pr.t === "pair") return `${pr.h} = [${r.value.map(pyNum).join(", ")}]`;
  if(pr.t === "int") return pr.h;
  return `${pr.h} = ${pyNum(r.value)}${pr.u ? " " + pr.u : ""}`;
}
function refreshConn(){ const el = $("pConn"); if(el && sel && sel.kind === "comp"){ const c = comp(sel.id); if(c) el.innerHTML = connText(c); } }
function render(){ renderCanvas(); renderProps(); }

/* ---------- actions ---------- */
function nextName(type){
  const pre = TYPES[type].prefix; if(type === "GND") return "GND";
  let n = 1; const used = new Set(state.comps.map(c => c.name));
  while(used.has(pre + n)) n++;
  return pre + n;
}
function addComp(type, x, y){
  const params = {}; TYPES[type].params.forEach(p => params[p.k] = p.d);
  const c = {id:"c" + (state.seq++), type, x:snap(x), y:snap(y), rot:0, name:nextName(type), params};
  while(state.comps.some(o => Math.abs(o.x - c.x) < 100 && Math.abs(o.y - c.y) < 80)){ c.y += GRID*4; }
  state.comps.push(c); sel = {kind:"comp", id:c.id}; render(); save();
}
function addWire(a, b){
  if(a.c === b.c && a.p === b.p) return;
  const dup = state.wires.some(w => (w.a.c===a.c&&w.a.p===a.p&&w.b.c===b.c&&w.b.p===b.p) || (w.a.c===b.c&&w.a.p===b.p&&w.b.c===a.c&&w.b.p===a.p));
  if(!dup) state.wires.push({id:"w" + (state.seq++), a:{...a}, b:{...b}});
  save();
}
function deleteSel(){
  if(!sel) return;
  if(sel.kind === "comp"){ state.comps = state.comps.filter(c => c.id !== sel.id); state.wires = state.wires.filter(w => w.a.c !== sel.id && w.b.c !== sel.id); }
  else state.wires = state.wires.filter(w => w.id !== sel.id);
  sel = null; pending = null; rubber.setAttribute("visibility","hidden"); render(); save();
}
function rotateSel(){
  if(!sel || sel.kind !== "comp") return;
  const c = comp(sel.id); c.rot = (c.rot + 90) % 360; renderCanvas(); save();
}

/* ---------- palette ---------- */
const partsEl = $("parts");
Object.entries(TYPES).forEach(([k,T]) => {
  const b = document.createElement("button");
  b.className = "part"; b.draggable = true;
  b.innerHTML = `<svg viewBox="-46 -26 92 52"><path d="${T.d}"/></svg><span>${T.label}${T.cls ? `<code>${T.mod ? "" : "qu."}${T.cls}</code>` : ""}</span>`;
  b.addEventListener("dragstart", e => { e.dataTransfer.setData("text/plain", k); e.dataTransfer.effectAllowed = "copy"; });
  b.addEventListener("click", () => addComp(k, view.x + view.w/2, view.y + view.h/2));
  partsEl.appendChild(b);
});
svg.addEventListener("dragover", e => { e.preventDefault(); e.dataTransfer.dropEffect = "copy"; });
svg.addEventListener("drop", e => {
  e.preventDefault(); const k = e.dataTransfer.getData("text/plain");
  if(TYPES[k]){ const p = toSvg(e); addComp(k, p.x, p.y); }
});

/* ---------- interaction ---------- */
svg.addEventListener("pointerdown", e => {
  const t = e.target, p = toSvg(e);
  const pinEl = t.closest("[data-pin]");
  if(pinEl){
    const here = {c:pinEl.dataset.c, p:+pinEl.dataset.pin};
    if(!pending) pending = here; else { addWire(pending, here); pending = null; }
    rubber.setAttribute("visibility", pending ? "visible" : "hidden");
    if(pending){ const q = pinPos(comp(here.c), here.p); ["x1","x2"].forEach(a => rubber.setAttribute(a, q.x)); ["y1","y2"].forEach(a => rubber.setAttribute(a, q.y)); }
    renderCanvas(); refreshConn(); return;
  }
  const cEl = t.closest("[data-comp]");
  if(cEl){
    const c = comp(cEl.dataset.comp);
    const same = sel && sel.kind === "comp" && sel.id === c.id;
    sel = {kind:"comp", id:c.id}; drag = {id:c.id, dx:p.x - c.x, dy:p.y - c.y};
    svg.setPointerCapture(e.pointerId);
    if(same) renderCanvas(); else render();
    return;
  }
  const wEl = t.closest("[data-wire]");
  if(wEl){ sel = {kind:"wire", id:wEl.dataset.wire}; render(); return; }
  const had = sel !== null;
  sel = null; pending = null; rubber.setAttribute("visibility","hidden");
  pan = {sx:e.clientX, sy:e.clientY, vx:view.x, vy:view.y};
  svg.setPointerCapture(e.pointerId);
  if(had) render(); else renderCanvas();
});
svg.addEventListener("pointermove", e => {
  if(drag){
    const c = comp(drag.id), p = toSvg(e);
    const nx = snap(p.x - drag.dx), ny = snap(p.y - drag.dy);
    if(nx !== c.x || ny !== c.y){ c.x = nx; c.y = ny; renderCanvas(); }
  } else if(pan){
    const k = view.w / svg.clientWidth;
    view.x = pan.vx - (e.clientX - pan.sx)*k; view.y = pan.vy - (e.clientY - pan.sy)*k; setView();
  } else if(pending){
    const p = toSvg(e); rubber.setAttribute("x2", p.x); rubber.setAttribute("y2", p.y);
  }
});
const endPointer = () => { if(drag) save(); drag = null; pan = null; };
svg.addEventListener("pointerup", endPointer);
svg.addEventListener("pointercancel", endPointer);
svg.addEventListener("wheel", e => {
  e.preventDefault();
  const p = toSvg(e), f = Math.exp(e.deltaY * 0.0015);
  const nw = Math.min(Math.max(view.w * f, 150), 6000), r = nw / view.w;
  view.x = p.x - (p.x - view.x)*r; view.y = p.y - (p.y - view.y)*r; view.w = nw; view.h *= r; setView();
}, {passive:false});
document.addEventListener("keydown", e => {
  if(e.target.tagName === "INPUT") return;
  if(e.key === "Delete" || e.key === "Backspace"){ e.preventDefault(); deleteSel(); }
  else if(e.key === "r" || e.key === "R") rotateSel();
  else if(e.key === "Escape"){ pending = null; sel = null; rubber.setAttribute("visibility","hidden"); render(); }
});

/* ---------- analysis panel ---------- */
$("fMin").value = state.fmin; $("fMax").value = state.fmax;
$("fMin").addEventListener("input", e => { state.fmin = e.target.value; lastOut = buildOutput(); renderOutput(); save(); });
$("fMax").addEventListener("input", e => { state.fmax = e.target.value; lastOut = buildOutput(); renderOutput(); save(); });
$("tabs").addEventListener("click", e => {
  const b = e.target.closest(".tab"); if(!b) return;
  fmt = b.dataset.fmt;
  document.querySelectorAll(".tab").forEach(t => t.setAttribute("aria-selected", t === b ? "true" : "false"));
  renderOutput();
});
$("btnCopy").addEventListener("click", async () => {
  const txt = $("code").textContent;
  try{ await navigator.clipboard.writeText(txt); }
  catch(_){ const r = document.createRange(); r.selectNodeContents($("code")); const s = getSelection(); s.removeAllRanges(); s.addRange(r); try{ document.execCommand("copy"); }catch(__){} s.removeAllRanges(); }
  $("btnCopy").textContent = "Copied";
  setTimeout(() => $("btnCopy").textContent = "Copy", 1400);
});
$("btnClear").addEventListener("click", () => { if(confirm("Clear the schematic?")){ state = {seq:1, fmin:state.fmin, fmax:state.fmax, comps:[], wires:[]}; sel = pending = null; render(); save(); } });
$("btnExample").addEventListener("click", () => { state = example(); $("fMin").value = state.fmin; $("fMax").value = state.fmax; sel = pending = null; render(); fitAll(); save(); });
$("btnFit").addEventListener("click", fitAll);

/* ---------- link with Python (qultra.open_editor) ---------- */
function sendMsg(text, bad){ const m = $("sendMsg"); m.hidden = false; m.textContent = text; m.classList.toggle("bad", !!bad); }
fetch("/api/init").then(r => r.ok ? r.json() : null).then(d => {
  if(!d || !d.qultra_editor) return;
  $("btnSend").hidden = false;
  if(d.design && Array.isArray(d.design.comps)) state = d.design;
  else state = {seq:1, fmin:"3", fmax:"12", comps:[], wires:[]};
  if(!state.seq) state.seq = 1;
  if(d.f_min !== null && d.f_min !== undefined) state.fmin = String(d.f_min);
  if(d.f_max !== null && d.f_max !== undefined) state.fmax = String(d.f_max);
  if(state.fmin === undefined) state.fmin = "3";
  if(state.fmax === undefined) state.fmax = "12";
  $("fMin").value = state.fmin; $("fMax").value = state.fmax;
  sel = pending = null; render(); fitAll();
  sendMsg(d.path ? `Editing ${d.path}. Press “Send to Python” when you are done.` : "Press “Send to Python” when you are done.");
}).catch(() => {});
$("btnSend").addEventListener("click", async () => {
  const out = buildOutput();
  if(!state.comps.some(c => c.type !== "GND")){ sendMsg("Add at least one component first.", true); return; }
  if(out.warns.length && !confirm("The circuit has warnings:\n\n" + out.warns.join("\n") + "\n\nSend it to Python anyway?")) return;
  $("btnSend").disabled = true;
  try{
    const r = await fetch("/api/submit", {method:"POST", headers:{"Content-Type":"application/json"},
      body: JSON.stringify({design: state, result: JSON.parse(out.json), script: out.script})});
    const d = await r.json();
    if(d.ok){ sendMsg("Sent to Python. You can close this tab.", false); }
    else { sendMsg("Python could not build the circuit: " + d.error, true); $("btnSend").disabled = false; }
  }catch(e){
    sendMsg("Python is no longer listening. Run open_editor() again.", true); $("btnSend").disabled = false;
  }
});

new ResizeObserver(fitAspect).observe(svg);
render(); fitAll();
})();
</script>
</body>
</html>
'''