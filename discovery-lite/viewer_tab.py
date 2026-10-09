"""
3D structure viewer tab, built on 3Dmol.js inside a QWebEngineView.

This replaces the old standalone main.py window with a tab that lives inside
redock_gui.py and can be driven two ways:

  1. Automatically: RedockValidationTab emits `runCompleted(work_dir, target,
     result)` when a run finishes, and MainWindow wires that straight into
     ViewerTab.load_run_result(). The receptor, native (crystal) ligand, and
     the docked mode-1 pose all get loaded and overlaid, colored differently,
     so a PASS/FAIL is something you can *see*, not just read off an RMSD
     number.
  2. Manually: separate "Load Receptor / Load Ligand / Load Pose" buttons for
     browsing to any .pdb/.pdbqt file, e.g. from the Pocket Detection or
     Docking Only tabs, or files that never went through this app at all.

Unlike main.py's original viewer (which rebuilt the whole HTML page - and
re-ran the 3Dmol.js download - on every single load), the 3Dmol page here is
built once. Everything after that - loading a new structure, restyling,
toggling spin/surface, highlighting contacts - goes through
QWebEnginePage.runJavaScript() calls to a small set of JS functions defined
on the page. Structure text is passed through json.dumps() rather than
hand-rolled string escaping, which is what main.py did (replace("\\", "\\\\")
etc.) and is easy to get subtly wrong on real PDB files.
"""

import base64
import json
import os

from PyQt6.QtCore import Qt
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QPushButton, QLabel,
    QTextEdit, QComboBox, QCheckBox, QFileDialog, QMessageBox
)


PAGE_HTML = """
<html>
<head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  html, body { margin: 0; padding: 0; height: 100%; background: #ffffff; }
  #viewer { width: 100%; height: 100%; position: relative; }
</style>
</head>
<body>
<div id="viewer"></div>
<script>
  let viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
  let models = {};    // id -> GLModel
  let surfaces = {};  // id -> [surface handles]

  function clearModel(id) {
    if (surfaces[id]) {
      surfaces[id].forEach(function(s) { viewer.removeSurface(s); });
      delete surfaces[id];
    }
    if (models[id]) {
      viewer.removeModel(models[id]);
      delete models[id];
    }
  }

  function clearAll() {
    Object.keys(models).forEach(clearModel);
    viewer.render();
  }

  function loadStructure(id, data, format, styleJson, doZoom) {
    clearModel(id);
    let m = viewer.addModel(data, format);
    m.setStyle({}, styleJson);
    models[id] = m;
    if (doZoom) { viewer.zoomTo(); }
    viewer.render();
  }

  function restyle(id, styleJson) {
    if (models[id]) {
      models[id].setStyle({}, styleJson);
      viewer.render();
    }
  }

  function setSurface(id, show, opacity, color) {
    if (surfaces[id]) {
      surfaces[id].forEach(function(s) { viewer.removeSurface(s); });
      delete surfaces[id];
    }
    if (show && models[id]) {
      let handle = viewer.addSurface($3Dmol.SurfaceType.VDW,
        { opacity: opacity, color: color }, { model: models[id] });
      surfaces[id] = [handle];
    }
    viewer.render();
  }

  function highlightContacts(receptorId, aroundId, on, distance, baseStyleJson) {
    if (!models[receptorId]) { return; }
    if (on && models[aroundId]) {
      models[receptorId].setStyle({}, baseStyleJson);
      models[receptorId].setStyle(
        { within: { distance: distance, sel: { model: models[aroundId] } } },
        { stick: { colorscheme: "yellowCarbon", radius: 0.18 } }
      );
    } else {
      models[receptorId].setStyle({}, baseStyleJson);
    }
    viewer.render();
  }

  function setSpin(on) { viewer.spin(on ? "y" : false); }

  function setBackground(color) { viewer.setBackgroundColor(color); viewer.render(); }

  function zoomToAll() { viewer.zoomTo(); viewer.render(); }

  function captureImage() { return viewer.pngURI(); }
</script>
</body>
</html>
"""

# Style presets, keyed by what the receptor style dropdown shows.
RECEPTOR_STYLES = {
    "Cartoon": {"cartoon": {"colorscheme": "spectrum"}},
    "Cartoon (by chain)": {"cartoon": {"colorscheme": "chain"}},
    "Cartoon (white)": {"cartoon": {"color": "white"}},
    "Sticks": {"stick": {"colorscheme": "grayCarbon", "radius": 0.15}},
    "Lines": {"line": {}},
    "Surface only": {"cartoon": {"hidden": True}},  # surface handled separately
}

LIGAND_STYLES = {
    "Sticks": {"stick": {"radius": 0.22}},
    "Spheres": {"sphere": {"scale": 0.35}},
    "Sticks + Spheres": {"stick": {"radius": 0.18}, "sphere": {"scale": 0.18}},
}

NATIVE_LIGAND_COLOR = "limeCarbon"
POSE_LIGAND_COLOR = "magentaCarbon"

# The Viewer tab loads raw file content straight into 3Dmol.js - there's no
# obabel conversion step here like there is in the docking pipeline scripts,
# so only formats 3Dmol.js can parse natively belong in these dialogs/maps.
# (.smi/.cml/.mrv/.mol are valid *docking* inputs elsewhere in this app, but
# aren't directly renderable here without converting them first.)
_VIEWER_FORMAT_BY_EXT = {
    "pdb": "pdb", "ent": "pdb", "pdbqt": "pdbqt", "sdf": "sdf", "mol2": "mol2", "ml2": "mol2",
}
RECEPTOR_FILE_FILTER = "Receptor Files (*.pdb *.pdbqt *.ent);;All Files (*)"
LIGAND_FILE_FILTER = "Ligand Files (*.pdb *.pdbqt *.sdf *.mol2 *.ml2 *.ent);;All Files (*)"
POSE_FILE_FILTER = "Docked Pose Files (*.pdbqt *.pdb *.sdf *.mol2);;All Files (*)"


def _guess_3dmol_format(path):
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    return _VIEWER_FORMAT_BY_EXT.get(ext, "pdb")


def _js_str(text):
    """Safely embed arbitrary text (PDB data, paths) as a JS string literal."""
    return json.dumps(text)


def _read_text(path):
    with open(path, "r", errors="replace") as f:
        return f.read()


def extract_first_pose_as_pdb(pdbqt_path):
    """Pull just the first MODEL...ENDMDL block out of a multi-pose Vina
    out.pdbqt (mode 1 = best pose) and return it as PDB-compatible text.
    Falls back to reading the whole file if there's no MODEL wrapper at all
    (a single-pose pdbqt). AutoDock's extra trailing columns (partial charge,
    atom type) are trimmed off since they sit past where a normal PDB parser
    stops reading - harmless to drop them, we only need coordinates here."""
    lines = []
    in_model = False
    seen_model = False
    with open(pdbqt_path, "r", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                if seen_model:
                    break
                in_model = True
                seen_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                if seen_model and not in_model:
                    continue
                lines.append(line[:66].rstrip("\n") + "\n")
    if not lines:
        with open(pdbqt_path, "r", errors="replace") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    lines.append(line[:66].rstrip("\n") + "\n")
    lines.append("END\n")
    return "".join(lines)


def summarize_structure(path, label):
    """Same chain/ATOM/HETATM tally main.py showed, generalized to any file."""
    chains = set()
    atom_count = 0
    hetatm_count = 0
    with open(path, "r", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM"):
                if len(line) > 21:
                    chains.add(line[21])
                atom_count += 1
            elif line.startswith("HETATM"):
                hetatm_count += 1
    lines = [f"{label}: {os.path.basename(path)}"]
    if chains:
        lines.append(f"Chains: {', '.join(sorted(chains))}")
    lines.append(f"ATOM lines: {atom_count}")
    lines.append(f"HETATM lines: {hetatm_count}")
    return "\n".join(lines)


class ViewerTab(QWidget):
    """Tab: 3Dmol structure viewer with receptor/ligand/pose overlay,
    style controls, and an auto-load hook fed by RedockValidationTab."""

    def __init__(self):
        super().__init__()
        self._last_run = None  # (work_dir, target, result) from the last auto-load
        self._loaded = {"receptor": None, "ligand": None, "pose": None}  # id -> source path/desc

        outer = QHBoxLayout()
        self.setLayout(outer)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter)

        # ---- left panel: load buttons, style controls, info box ----
        left = QWidget()
        left.setMaximumWidth(320)
        left_layout = QVBoxLayout()
        left.setLayout(left_layout)

        left_layout.addWidget(QLabel("<b>Load</b>"))
        self.reload_latest_btn = QPushButton("Reload Latest Redocking Run")
        self.reload_latest_btn.clicked.connect(self.reload_latest_run)
        self.reload_latest_btn.setEnabled(False)
        left_layout.addWidget(self.reload_latest_btn)

        load_receptor_btn = QPushButton("Load Receptor...")
        load_receptor_btn.clicked.connect(self.browse_receptor)
        left_layout.addWidget(load_receptor_btn)

        load_ligand_btn = QPushButton("Load Ligand (native/reference)...")
        load_ligand_btn.clicked.connect(self.browse_ligand)
        left_layout.addWidget(load_ligand_btn)

        load_pose_btn = QPushButton("Load Docked Pose...")
        load_pose_btn.clicked.connect(self.browse_pose)
        left_layout.addWidget(load_pose_btn)

        clear_btn = QPushButton("Clear Viewer")
        clear_btn.clicked.connect(self.clear_viewer)
        left_layout.addWidget(clear_btn)

        left_layout.addWidget(QLabel("<b>Style</b>"))

        style_row1 = QHBoxLayout()
        style_row1.addWidget(QLabel("Receptor:"))
        self.receptor_style = QComboBox()
        self.receptor_style.addItems(list(RECEPTOR_STYLES.keys()))
        self.receptor_style.currentTextChanged.connect(self.apply_receptor_style)
        style_row1.addWidget(self.receptor_style)
        left_layout.addLayout(style_row1)

        style_row2 = QHBoxLayout()
        style_row2.addWidget(QLabel("Ligands:"))
        self.ligand_style = QComboBox()
        self.ligand_style.addItems(list(LIGAND_STYLES.keys()))
        self.ligand_style.currentTextChanged.connect(self.apply_ligand_style)
        style_row2.addWidget(self.ligand_style)
        left_layout.addLayout(style_row2)

        self.surface_check = QCheckBox("Receptor surface (transparent)")
        self.surface_check.stateChanged.connect(self.apply_surface)
        left_layout.addWidget(self.surface_check)

        self.contacts_check = QCheckBox("Highlight residues near native ligand (4.5 A)")
        self.contacts_check.stateChanged.connect(self.apply_contacts)
        left_layout.addWidget(self.contacts_check)

        self.spin_check = QCheckBox("Spin")
        self.spin_check.stateChanged.connect(self.apply_spin)
        left_layout.addWidget(self.spin_check)

        bg_row = QHBoxLayout()
        bg_row.addWidget(QLabel("Background:"))
        self.bg_combo = QComboBox()
        self.bg_combo.addItems(["White", "Black"])
        self.bg_combo.currentTextChanged.connect(self.apply_background)
        bg_row.addWidget(self.bg_combo)
        left_layout.addLayout(bg_row)

        view_row = QHBoxLayout()
        reset_btn = QPushButton("Reset View")
        reset_btn.clicked.connect(self.reset_view)
        save_img_btn = QPushButton("Save Image...")
        save_img_btn.clicked.connect(self.save_image)
        view_row.addWidget(reset_btn)
        view_row.addWidget(save_img_btn)
        left_layout.addLayout(view_row)

        self.status_label = QLabel("Nothing loaded yet.")
        self.status_label.setWordWrap(True)
        left_layout.addWidget(self.status_label)

        left_layout.addWidget(QLabel("<b>Structure Info</b>"))
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        left_layout.addWidget(self.info_box)

        left_layout.addStretch()
        splitter.addWidget(left)

        # ---- right panel: the 3Dmol view itself ----
        self.view = QWebEngineView()
        self.view.setHtml(PAGE_HTML)
        splitter.addWidget(self.view)
        splitter.setSizes([320, 700])

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------

    def _run_js(self, script, callback=None):
        if callback is not None:
            self.view.page().runJavaScript(script, callback)
        else:
            self.view.page().runJavaScript(script)

    def _current_receptor_style(self):
        return RECEPTOR_STYLES[self.receptor_style.currentText()]

    def _current_ligand_style(self, color):
        style = dict(LIGAND_STYLES[self.ligand_style.currentText()])
        # stamp the color onto whichever rep(s) this preset uses
        for rep in style:
            style[rep] = dict(style[rep])
            style[rep]["colorscheme"] = color
        return style

    def _load_model(self, slot_id, path, fmt, style_json, do_zoom=True):
        data = _read_text(path)
        script = f"loadStructure({_js_str(slot_id)}, {_js_str(data)}, {_js_str(fmt)}, {json.dumps(style_json)}, {json.dumps(do_zoom)});"
        self._run_js(script)
        self._loaded[slot_id] = path

    def _refresh_info_box(self):
        blocks = []
        if self._loaded.get("receptor"):
            blocks.append(summarize_structure(self._loaded["receptor"], "Receptor"))
        if self._loaded.get("ligand"):
            blocks.append(summarize_structure(self._loaded["ligand"], "Ligand"))
        self.info_box.setPlainText("\n\n".join(blocks) if blocks else "")

    # ------------------------------------------------------------------
    # Auto-load hook (connected to RedockValidationTab.runCompleted)
    # ------------------------------------------------------------------

    def load_run_result(self, work_dir, target, result):
        """Load receptor.pdb + native_ligand.pdb + mode-1 of out.pdbqt from a
        just-finished Redocking Validation run, overlaid so PASS/FAIL is
        something you can see, not just a number."""
        self._last_run = (work_dir, target, result)
        self.reload_latest_btn.setEnabled(True)

        receptor_path = os.path.join(work_dir, "receptor.pdb")
        native_path = os.path.join(work_dir, "native_ligand.pdb")
        pose_path = os.path.join(work_dir, "out.pdbqt")

        if not os.path.isfile(receptor_path):
            self.status_label.setText(f"Run for {target} finished, but receptor.pdb wasn't found in {work_dir}.")
            return

        self._load_model("receptor", receptor_path, "pdb", self._current_receptor_style(), do_zoom=False)

        if os.path.isfile(native_path):
            self._load_model("ligand", native_path, "pdb", self._current_ligand_style(NATIVE_LIGAND_COLOR), do_zoom=False)
        else:
            self._loaded["ligand"] = None

        if os.path.isfile(pose_path):
            pose_text = extract_first_pose_as_pdb(pose_path)
            script = (
                f"loadStructure('pose', {_js_str(pose_text)}, 'pdb', "
                f"{json.dumps(self._current_ligand_style(POSE_LIGAND_COLOR))}, true);"
            )
            self._run_js(script)
            self._loaded["pose"] = pose_path
        else:
            self._loaded["pose"] = None
            self._run_js("zoomToAll();")

        self._refresh_info_box()
        legend = "green = native crystal pose, magenta = docked mode 1" if self._loaded["pose"] and self._loaded["ligand"] else ""
        self.status_label.setText(f"Loaded {target} - {result}. {legend}".strip())

        # re-apply any active surface/contacts toggles to the freshly loaded models
        self.apply_surface()
        self.apply_contacts()

    def reload_latest_run(self):
        if self._last_run:
            self.load_run_result(*self._last_run)

    # ------------------------------------------------------------------
    # Manual load (any file, any tab)
    # ------------------------------------------------------------------

    def load_receptor_from_path(self, path):
        if not path or not os.path.isfile(path):
            return
        fmt = _guess_3dmol_format(path)
        self._load_model("receptor", path, fmt, self._current_receptor_style(), do_zoom=True)
        self._refresh_info_box()
        self.status_label.setText(f"Receptor loaded: {os.path.basename(path)}")
        self.apply_surface()
        self.apply_contacts()

    def load_ligand_from_path(self, path):
        if not path or not os.path.isfile(path):
            return
        fmt = _guess_3dmol_format(path)
        self._load_model("ligand", path, fmt, self._current_ligand_style(NATIVE_LIGAND_COLOR), do_zoom=False)
        self._run_js("zoomToAll();")
        self._refresh_info_box()
        self.status_label.setText(f"Ligand loaded: {os.path.basename(path)}")
        self.apply_contacts()

    def load_pose_from_path(self, path):
        if not path or not os.path.isfile(path):
            return
        if path.lower().endswith(".pdbqt"):
            text = extract_first_pose_as_pdb(path)
            fmt = "pdb"
        else:
            text = _read_text(path)
            fmt = _guess_3dmol_format(path)
        script = (
            f"loadStructure('pose', {_js_str(text)}, {_js_str(fmt)}, "
            f"{json.dumps(self._current_ligand_style(POSE_LIGAND_COLOR))}, true);"
        )
        self._run_js(script)
        self._loaded["pose"] = path
        self.status_label.setText(f"Docked pose loaded: {os.path.basename(path)}")

    def browse_receptor(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Receptor File", "", RECEPTOR_FILE_FILTER)
        if path:
            self.load_receptor_from_path(path)

    def browse_ligand(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Ligand File", "", LIGAND_FILE_FILTER)
        if path:
            self.load_ligand_from_path(path)

    def browse_pose(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Docked Pose File", "", POSE_FILE_FILTER)
        if path:
            self.load_pose_from_path(path)

    def clear_viewer(self):
        self._run_js("clearAll();")
        self._loaded = {"receptor": None, "ligand": None, "pose": None}
        self.info_box.clear()
        self.status_label.setText("Nothing loaded yet.")

    # ------------------------------------------------------------------
    # Style controls
    # ------------------------------------------------------------------

    def apply_receptor_style(self, *_args):
        if self._loaded.get("receptor"):
            self._run_js(f"restyle('receptor', {json.dumps(self._current_receptor_style())});")
            self.apply_surface()
            self.apply_contacts()

    def apply_ligand_style(self, *_args):
        if self._loaded.get("ligand"):
            self._run_js(f"restyle('ligand', {json.dumps(self._current_ligand_style(NATIVE_LIGAND_COLOR))});")
        if self._loaded.get("pose"):
            self._run_js(f"restyle('pose', {json.dumps(self._current_ligand_style(POSE_LIGAND_COLOR))});")

    def apply_surface(self, *_args):
        show = self.surface_check.isChecked() and bool(self._loaded.get("receptor"))
        self._run_js(f"setSurface('receptor', {json.dumps(show)}, 0.5, 'white');")

    def apply_contacts(self, *_args):
        show = self.contacts_check.isChecked() and bool(self._loaded.get("receptor")) and bool(self._loaded.get("ligand"))
        base_style = json.dumps(self._current_receptor_style())
        self._run_js(f"highlightContacts('receptor', 'ligand', {json.dumps(show)}, 4.5, {base_style});")

    def apply_spin(self, *_args):
        self._run_js(f"setSpin({json.dumps(self.spin_check.isChecked())});")

    def apply_background(self, text):
        color = "black" if text == "Black" else "white"
        self._run_js(f"setBackground({_js_str(color)});")

    def reset_view(self):
        self._run_js("zoomToAll();")

    def save_image(self):
        def handle_result(data_uri):
            if not data_uri or "," not in data_uri:
                QMessageBox.warning(self, "Save failed", "Nothing is loaded in the viewer yet.")
                return
            _header, b64data = data_uri.split(",", 1)
            try:
                raw = base64.b64decode(b64data)
            except Exception:
                QMessageBox.warning(self, "Save failed", "Could not decode the image data.")
                return
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Viewer Image", "structure.png", "PNG Files (*.png)")
            if not file_path:
                return
            try:
                with open(file_path, "wb") as f:
                    f.write(raw)
            except OSError as exc:
                QMessageBox.warning(self, "Save failed", f"Could not write file:\n{exc}")
                return
            QMessageBox.information(self, "Saved", f"Saved image to:\n{file_path}")

        self._run_js("captureImage();", handle_result)
