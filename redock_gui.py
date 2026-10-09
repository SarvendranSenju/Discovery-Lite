import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QLineEdit, QComboBox, QPushButton, QSpinBox, QDoubleSpinBox, QTextEdit, QFileDialog, QLabel,
    QTabWidget, QCheckBox
)
from PyQt6.QtCore import QProcess, pyqtSignal
from PyQt6.QtGui import QFont, QTextCursor

from session_memory import SessionMemory
from results_log import ResultsLogStore, ResultsLogTab
from viewer_tab import ViewerTab


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Shared Qt file-dialog filters, kept in sync with what the backend scripts
# actually accept (see _OBABEL_FORMAT_BY_EXT in redock_and_validate.py /
# prepare_structures.py for ligands; receptor/complex parsing reads raw
# PDB-style ATOM/HETATM records, so .pdbqt and the legacy .ent PDB extension
# work there too, but not formats like mmCIF that use a different layout).
RECEPTOR_FILE_FILTER = "Receptor Files (*.pdb *.pdbqt *.ent);;All Files (*)"
COMPLEX_FILE_FILTER = "PDB Complex Files (*.pdb *.pdbqt *.ent);;All Files (*)"
LIGAND_FILE_FILTER = (
    "Ligand Files (*.sdf *.mol *.mol2 *.ml2 *.pdb *.pdbqt *.smi *.cml *.mrv);;All Files (*)"
)

# Shared across every tab: recent receptor/ligand/PDB paths (session_memory.json)
# and the aggregated affinity+RMSD log (results_log.json). Both are plain JSON
# files next to the scripts here - delete either to reset it.
memory = SessionMemory(os.path.join(SCRIPT_DIR, "session_memory.json"))
results_store = ResultsLogStore(os.path.join(SCRIPT_DIR, "results_log.json"))


def make_recent_combo(key, placeholder):
    """An editable combo box pre-filled with the last-used value for `key`
    and a dropdown of recent values, backed by session_memory.json. Behaves
    like a QLineEdit (free typing works) with a shortcut for reused paths -
    call .currentText() / .setCurrentText() on it exactly as you would
    .text() / .setText() on a QLineEdit."""
    combo = QComboBox()
    combo.setEditable(True)
    combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
    combo.lineEdit().setPlaceholderText(placeholder)
    combo.addItems(memory.get_recent(key))
    combo.setCurrentText(memory.get_last(key, ""))
    return combo


class RedockValidationTab(QWidget):
    """Tab 1: full pipeline with auto-detected native ligand + RMSD validation."""

    # (work_dir, target_label, PASS/FAIL) - emitted once a run finishes with a
    # parseable RESULT_ROW line, so the Viewer tab can auto-load it.
    runCompleted = pyqtSignal(str, str, str)

    def __init__(self):
        super().__init__()
        self.script_path = os.path.join(SCRIPT_DIR, "redock_and_validate.py")
        self.process = None

        layout = QVBoxLayout()
        self.setLayout(layout)

        form = QFormLayout()

        pdb_row = QHBoxLayout()
        self.pdb_input = make_recent_combo("redock_pdb", "e.g. 8FK4, or path to a .pdb file")
        browse_pdb_btn = QPushButton("Browse")
        browse_pdb_btn.clicked.connect(self.browse_pdb)
        pdb_row.addWidget(self.pdb_input)
        pdb_row.addWidget(browse_pdb_btn)
        form.addRow("PDB ID / File:", pdb_row)

        ligand_row = QHBoxLayout()
        self.ligand_input = make_recent_combo("redock_ligand", "Optional: SMILES, file path, or URL")
        browse_ligand_btn = QPushButton("Browse")
        browse_ligand_btn.clicked.connect(self.browse_ligand)
        ligand_row.addWidget(self.ligand_input)
        ligand_row.addWidget(browse_ligand_btn)
        form.addRow("Ligand Override:", ligand_row)

        self.exhaustiveness = QSpinBox()
        self.exhaustiveness.setRange(1, 64)
        self.exhaustiveness.setValue(8)
        form.addRow("Exhaustiveness:", self.exhaustiveness)

        self.num_modes = QSpinBox()
        self.num_modes.setRange(1, 20)
        self.num_modes.setValue(9)
        form.addRow("Num Modes:", self.num_modes)

        self.rmsd_threshold = QDoubleSpinBox()
        self.rmsd_threshold.setRange(0.1, 10.0)
        self.rmsd_threshold.setSingleStep(0.1)
        self.rmsd_threshold.setValue(2.0)
        form.addRow("RMSD Pass Threshold (A):", self.rmsd_threshold)

        layout.addLayout(form)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run Redocking Validation")
        self.run_btn.clicked.connect(self.run_pipeline)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_pipeline)
        self.stop_btn.setEnabled(False)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.stop_btn)
        layout.addLayout(btn_row)

        self.result_label = QLabel("")
        self.result_label.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        layout.addWidget(self.result_label)

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Consolas", 10))
        layout.addWidget(self.console)

    def browse_pdb(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PDB File", "", COMPLEX_FILE_FILTER)
        if file_path:
            self.pdb_input.setCurrentText(file_path)

    def browse_ligand(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Ligand File", "", LIGAND_FILE_FILTER)
        if file_path:
            self.ligand_input.setCurrentText(file_path)

    def run_pipeline(self):
        pdb_val = self.pdb_input.currentText().strip()
        if not pdb_val:
            self.console.append("ERROR: Enter a PDB ID or select a PDB file first.\n")
            return

        args = ["-u", self.script_path, pdb_val]

        ligand_val = self.ligand_input.currentText().strip()
        if ligand_val:
            args += ["--ligand-input", ligand_val]

        args += ["--exhaustiveness", str(self.exhaustiveness.value())]
        args += ["--num-modes", str(self.num_modes.value())]
        args += ["--rmsd-threshold", str(self.rmsd_threshold.value())]

        memory.remember("redock_pdb", pdb_val)
        if ligand_val:
            memory.remember("redock_ligand", ligand_val)

        self.console.clear()
        self.result_label.setText("")
        self.result_label.setStyleSheet("")
        self.console.append(f"Running: python3 {' '.join(args)}\n")

        self.process = QProcess(self)
        self.process.setWorkingDirectory(SCRIPT_DIR)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.finished.connect(self.handle_finished)
        self.process.start("python3", args)

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

    def handle_output(self):
        data = self.process.readAllStandardOutput().data().decode(errors="replace")
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(data)
        self.console.moveCursor(QTextCursor.MoveOperation.End)

    def handle_finished(self):
        full_text = self.console.toPlainText()
        if "PASS: Mode 1 RMSD" in full_text:
            self.result_label.setText("PASS - redocking validated")
            self.result_label.setStyleSheet("color: green;")
        elif "FAIL" in full_text:
            self.result_label.setText("FAIL - check grid box / results")
            self.result_label.setStyleSheet("color: red;")
        else:
            self.result_label.setText("Run finished - check log")
            self.result_label.setStyleSheet("color: orange;")

        self._log_result(full_text)

        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _log_result(self, full_text):
        """Parse the RESULT_ROW,... line redock_and_validate.py prints at the
        end of a run and feed it to the shared Results Log tab."""
        for line in full_text.splitlines():
            if not line.startswith("RESULT_ROW,"):
                continue
            parts = line.split(",")
            if len(parts) < 9:
                continue
            try:
                affinity = float(parts[3]) if parts[3] else None
                rmsd = float(parts[4]) if parts[4] else None
                work_dir = ",".join(parts[8:])
                results_store.add_result({
                    "source": "Redocking Validation",
                    "target": parts[1],
                    "ligand": parts[2],
                    "affinity": affinity,
                    "rmsd": rmsd,
                    "result": parts[5],
                    "exhaustiveness": parts[6],
                    "num_modes": parts[7],
                    "work_dir": work_dir,
                })
                self.runCompleted.emit(work_dir, parts[1], parts[5])
            except ValueError:
                pass
            return  # only one RESULT_ROW per run

    def stop_pipeline(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.process.kill()
            self.console.append("\n--- Stopped by user ---\n")
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)


class DockOnlyTab(QWidget):
    """Tab 2: straight docking run, no native ligand / no RMSD needed.
    Box can be auto-computed from a small reference ligand file, typed in
    manually, or filled in automatically from the Pocket Detection tab."""

    def __init__(self, viewer_tab=None):
        super().__init__()
        self.script_path = os.path.join(SCRIPT_DIR, "dock_only.py")
        self.process = None
        self.viewer_tab = viewer_tab

        layout = QVBoxLayout()
        self.setLayout(layout)

        form = QFormLayout()

        receptor_row = QHBoxLayout()
        self.receptor_input = make_recent_combo("dock_receptor", "Path to receptor .pdb (or already-prepared .pdbqt)")
        browse_receptor_btn = QPushButton("Browse")
        browse_receptor_btn.clicked.connect(self.browse_receptor)
        receptor_row.addWidget(self.receptor_input)
        receptor_row.addWidget(browse_receptor_btn)
        form.addRow("Receptor:", receptor_row)

        ligand_row = QHBoxLayout()
        self.ligand_input = make_recent_combo("dock_ligand", "SMILES string, or path to .sdf/.mol/.mol2/.pdb")
        browse_ligand_btn = QPushButton("Browse")
        browse_ligand_btn.clicked.connect(self.browse_ligand)
        ligand_row.addWidget(self.ligand_input)
        ligand_row.addWidget(browse_ligand_btn)
        form.addRow("Ligand to Dock:", ligand_row)

        ref_row = QHBoxLayout()
        self.reference_input = make_recent_combo(
            "dock_reference", "Optional: small file to auto-compute the box from (e.g. native_ligand.pdb)"
        )
        browse_ref_btn = QPushButton("Browse")
        browse_ref_btn.clicked.connect(self.browse_reference)
        ref_row.addWidget(self.reference_input)
        ref_row.addWidget(browse_ref_btn)
        form.addRow("Reference Ligand (box):", ref_row)

        box_note = QLabel(
            "Box source: use a Reference Ligand file above, OR send a predicted pocket from the "
            "Pocket Detection tab, OR fill in the coordinates manually below."
        )
        box_note.setWordWrap(True)
        layout.addLayout(form)
        layout.addWidget(box_note)

        if self.viewer_tab is not None:
            view_inputs_btn = QPushButton("Preview Receptor + Ligand in 3D Viewer")
            view_inputs_btn.clicked.connect(self.view_inputs_in_3d)
            layout.addWidget(view_inputs_btn)

        box_form = QFormLayout()
        center_row = QHBoxLayout()
        self.center_x = QDoubleSpinBox(); self.center_x.setRange(-999, 999); self.center_x.setDecimals(3)
        self.center_y = QDoubleSpinBox(); self.center_y.setRange(-999, 999); self.center_y.setDecimals(3)
        self.center_z = QDoubleSpinBox(); self.center_z.setRange(-999, 999); self.center_z.setDecimals(3)
        center_row.addWidget(QLabel("X:")); center_row.addWidget(self.center_x)
        center_row.addWidget(QLabel("Y:")); center_row.addWidget(self.center_y)
        center_row.addWidget(QLabel("Z:")); center_row.addWidget(self.center_z)
        box_form.addRow("Box Center:", center_row)

        size_row = QHBoxLayout()
        self.size_x = QDoubleSpinBox(); self.size_x.setRange(1, 200); self.size_x.setValue(20)
        self.size_y = QDoubleSpinBox(); self.size_y.setRange(1, 200); self.size_y.setValue(20)
        self.size_z = QDoubleSpinBox(); self.size_z.setRange(1, 200); self.size_z.setValue(20)
        size_row.addWidget(QLabel("X:")); size_row.addWidget(self.size_x)
        size_row.addWidget(QLabel("Y:")); size_row.addWidget(self.size_y)
        size_row.addWidget(QLabel("Z:")); size_row.addWidget(self.size_z)
        box_form.addRow("Box Size:", size_row)

        self.exhaustiveness = QSpinBox()
        self.exhaustiveness.setRange(1, 64)
        self.exhaustiveness.setValue(8)
        box_form.addRow("Exhaustiveness:", self.exhaustiveness)

        self.num_modes = QSpinBox()
        self.num_modes.setRange(1, 20)
        self.num_modes.setValue(9)
        box_form.addRow("Num Modes:", self.num_modes)

        layout.addLayout(box_form)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run Docking")
        self.run_btn.clicked.connect(self.run_pipeline)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_pipeline)
        self.stop_btn.setEnabled(False)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.stop_btn)
        layout.addLayout(btn_row)

        self.result_label = QLabel("")
        self.result_label.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        layout.addWidget(self.result_label)

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Consolas", 10))
        layout.addWidget(self.console)

    def set_box_from_pocket(self, center_x, center_y, center_z, size=22.0):
        """Called by the Pocket Detection tab when the user sends a chosen pocket here."""
        self.reference_input.clearEditText()
        self.center_x.setValue(center_x)
        self.center_y.setValue(center_y)
        self.center_z.setValue(center_z)
        self.size_x.setValue(size)
        self.size_y.setValue(size)
        self.size_z.setValue(size)

    def view_inputs_in_3d(self):
        """Sends the currently-entered receptor/ligand *inputs* (not a docked
        result - dock_only.py isn't part of this bundle so its output layout
        isn't known here) to the Viewer tab, for a quick sanity-check look
        before running."""
        receptor_val = self.receptor_input.currentText().strip()
        ligand_val = self.ligand_input.currentText().strip()
        if not receptor_val and not ligand_val:
            self.console.append("ERROR: Enter a receptor and/or ligand first.\n")
            return
        if receptor_val and os.path.isfile(receptor_val):
            self.viewer_tab.load_receptor_from_path(receptor_val)
        if ligand_val and os.path.isfile(ligand_val):
            self.viewer_tab.load_ligand_from_path(ligand_val)
        if ligand_val and not os.path.isfile(ligand_val):
            self.console.append(
                f"Note: '{ligand_val}' isn't a local file (looks like a SMILES/URL), so it can't be "
                "previewed directly - the Viewer tab only opens structure files.\n"
            )

    def browse_receptor(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Receptor File", "", RECEPTOR_FILE_FILTER)
        if file_path:
            self.receptor_input.setCurrentText(file_path)

    def browse_ligand(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Ligand File", "", LIGAND_FILE_FILTER)
        if file_path:
            self.ligand_input.setCurrentText(file_path)

    def browse_reference(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Reference Ligand File", "", LIGAND_FILE_FILTER)
        if file_path:
            self.reference_input.setCurrentText(file_path)

    def run_pipeline(self):
        receptor_val = self.receptor_input.currentText().strip()
        ligand_val = self.ligand_input.currentText().strip()
        reference_val = self.reference_input.currentText().strip()

        if not receptor_val:
            self.console.append("ERROR: Select a receptor file first.\n")
            return
        if not ligand_val:
            self.console.append("ERROR: Provide a ligand (SMILES or file) first.\n")
            return
        if not reference_val and (self.center_x.value() == 0 and self.center_y.value() == 0 and self.center_z.value() == 0):
            self.console.append(
                "ERROR: Either provide a Reference Ligand file, send a pocket from the Pocket Detection "
                "tab, or set Box Center X/Y/Z manually (all three are currently 0).\n"
            )
            return

        args = ["-u", self.script_path, receptor_val, "--ligand-input", ligand_val]

        if reference_val:
            args += ["--reference-ligand", reference_val]
        else:
            args += [
                "--center-x", str(self.center_x.value()),
                "--center-y", str(self.center_y.value()),
                "--center-z", str(self.center_z.value()),
                "--size-x", str(self.size_x.value()),
                "--size-y", str(self.size_y.value()),
                "--size-z", str(self.size_z.value()),
            ]

        args += ["--exhaustiveness", str(self.exhaustiveness.value())]
        args += ["--num-modes", str(self.num_modes.value())]

        memory.remember("dock_receptor", receptor_val)
        memory.remember("dock_ligand", ligand_val)
        if reference_val:
            memory.remember("dock_reference", reference_val)

        self.console.clear()
        self.result_label.setText("")
        self.result_label.setStyleSheet("")
        self.console.append(f"Running: python3 {' '.join(args)}\n")

        self.process = QProcess(self)
        self.process.setWorkingDirectory(SCRIPT_DIR)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.finished.connect(self.handle_finished)
        self.process.start("python3", args)

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

    def handle_output(self):
        data = self.process.readAllStandardOutput().data().decode(errors="replace")
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(data)
        self.console.moveCursor(QTextCursor.MoveOperation.End)

    def handle_finished(self):
        full_text = self.console.toPlainText()
        if "DOCKING_COMPLETE: best affinity" in full_text:
            best = full_text.split("DOCKING_COMPLETE: best affinity =")[-1].split("kcal/mol")[0].strip()
            self.result_label.setText(f"Done - best affinity: {best} kcal/mol")
            self.result_label.setStyleSheet("color: green;")
            self._log_result(best)
        else:
            self.result_label.setText("Run finished - check log for issues")
            self.result_label.setStyleSheet("color: orange;")

        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _log_result(self, best_affinity_str):
        """dock_only.py isn't part of this bundle, so there's no RESULT_ROW
        marker to parse here (unlike the Redocking Validation tab). This logs
        what the GUI already has: the best affinity plus the inputs used.
        RMSD is left blank - there's no crystal pose to validate an arbitrary
        docked ligand against. If dock_only.py gets a RESULT_ROW line added
        later (matching redock_and_validate.py's format), swap this out for
        the same parsing _log_result uses on the Redocking Validation tab."""
        try:
            affinity = float(best_affinity_str)
        except ValueError:
            affinity = None
        results_store.add_result({
            "source": "Docking Only",
            "target": self.receptor_input.currentText().strip(),
            "ligand": self.ligand_input.currentText().strip(),
            "affinity": affinity,
            "rmsd": None,
            "result": "N/A",
            "exhaustiveness": self.exhaustiveness.value(),
            "num_modes": self.num_modes.value(),
            "work_dir": "",
        })

    def stop_pipeline(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.process.kill()
            self.console.append("\n--- Stopped by user ---\n")
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)


class PocketDetectionTab(QWidget):
    """Tab 3: runs P2Rank on a receptor to auto-detect candidate binding pockets,
    then lets the user send a chosen pocket's coordinates to the Docking Only tab."""

    def __init__(self, dock_only_tab: DockOnlyTab):
        super().__init__()
        self.dock_only_tab = dock_only_tab
        self.script_path = os.path.join(SCRIPT_DIR, "find_pockets.py")
        self.process = None
        self.pockets = []  # list of dicts parsed from POCKET_ROW lines

        layout = QVBoxLayout()
        self.setLayout(layout)

        form = QFormLayout()
        receptor_row = QHBoxLayout()
        self.receptor_input = make_recent_combo(
            "pocket_receptor", "Path to receptor .pdb (protein only, e.g. a cleaned structure)"
        )
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_receptor)
        receptor_row.addWidget(self.receptor_input)
        receptor_row.addWidget(browse_btn)
        form.addRow("Receptor:", receptor_row)
        layout.addLayout(form)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Find Pockets")
        self.run_btn.clicked.connect(self.run_pipeline)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_pipeline)
        self.stop_btn.setEnabled(False)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.stop_btn)
        layout.addLayout(btn_row)

        send_row = QHBoxLayout()
        send_row.addWidget(QLabel("Pocket rank to send:"))
        self.pocket_rank = QSpinBox()
        self.pocket_rank.setRange(1, 50)
        self.pocket_rank.setValue(1)
        send_row.addWidget(self.pocket_rank)
        self.send_btn = QPushButton("Send to Docking Only tab")
        self.send_btn.clicked.connect(self.send_pocket)
        send_row.addWidget(self.send_btn)
        send_row.addStretch()
        layout.addLayout(send_row)

        self.status_label = QLabel("")
        self.status_label.setFont(QFont("Arial", 11, QFont.Weight.Bold))
        layout.addWidget(self.status_label)

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Consolas", 10))
        layout.addWidget(self.console)

    def browse_receptor(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Receptor File", "", RECEPTOR_FILE_FILTER)
        if file_path:
            self.receptor_input.setCurrentText(file_path)

    def run_pipeline(self):
        receptor_val = self.receptor_input.currentText().strip()
        if not receptor_val:
            self.console.append("ERROR: Select a receptor PDB file first.\n")
            return

        args = ["-u", self.script_path, receptor_val]

        memory.remember("pocket_receptor", receptor_val)

        self.console.clear()
        self.status_label.setText("")
        self.status_label.setStyleSheet("")
        self.pockets = []
        self.console.append(f"Running: python3 {' '.join(args)}\n")

        self.process = QProcess(self)
        self.process.setWorkingDirectory(SCRIPT_DIR)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.finished.connect(self.handle_finished)
        self.process.start("python3", args)

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

    def handle_output(self):
        data = self.process.readAllStandardOutput().data().decode(errors="replace")
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(data)
        self.console.moveCursor(QTextCursor.MoveOperation.End)

    def handle_finished(self):
        full_text = self.console.toPlainText()
        self.pockets = []
        for line in full_text.splitlines():
            if line.startswith("POCKET_ROW,"):
                parts = line.split(",")
                # POCKET_ROW,rank,score,probability,center_x,center_y,center_z,n_residues
                try:
                    self.pockets.append({
                        "rank": int(parts[1]),
                        "score": float(parts[2]),
                        "probability": float(parts[3]) if parts[3] else None,
                        "center_x": float(parts[4]),
                        "center_y": float(parts[5]),
                        "center_z": float(parts[6]),
                        "n_residues": int(parts[7]) if len(parts) > 7 and parts[7] else None,
                    })
                except (ValueError, IndexError):
                    continue

        if "POCKETS_COMPLETE" in full_text and self.pockets:
            self.status_label.setText(f"Found {len(self.pockets)} candidate pocket(s) - pick a rank above and send it.")
            self.status_label.setStyleSheet("color: green;")
            self.pocket_rank.setMaximum(max(p["rank"] for p in self.pockets))
        elif "POCKETS_COMPLETE" in full_text:
            self.status_label.setText("Finished, but no pockets were predicted for this structure.")
            self.status_label.setStyleSheet("color: orange;")
        else:
            self.status_label.setText("Run finished - check log for issues")
            self.status_label.setStyleSheet("color: red;")

        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def send_pocket(self):
        target_rank = self.pocket_rank.value()
        match = next((p for p in self.pockets if p["rank"] == target_rank), None)
        if not match:
            self.console.append(f"\nERROR: no pocket with rank {target_rank} found in the last run.\n")
            return
        self.dock_only_tab.set_box_from_pocket(match["center_x"], match["center_y"], match["center_z"])
        self.console.append(
            f"\nSent pocket rank {target_rank} (center {match['center_x']:.2f}, "
            f"{match['center_y']:.2f}, {match['center_z']:.2f}) to the Docking Only tab.\n"
        )

    def stop_pipeline(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.process.kill()
            self.console.append("\n--- Stopped by user ---\n")
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)


class PrepareStructuresTab(QWidget):
    """Tab 4: prepares a receptor and/or ligand into PDBQT form and saves the
    result(s) to a folder of the user's choosing - no docking run involved."""

    def __init__(self):
        super().__init__()
        self.script_path = os.path.join(SCRIPT_DIR, "prepare_structures.py")
        self.process = None

        layout = QVBoxLayout()
        self.setLayout(layout)

        form = QFormLayout()

        receptor_row = QHBoxLayout()
        self.receptor_input = make_recent_combo("prepare_receptor", "Optional: path to receptor .pdb")
        browse_receptor_btn = QPushButton("Browse")
        browse_receptor_btn.clicked.connect(self.browse_receptor)
        receptor_row.addWidget(self.receptor_input)
        receptor_row.addWidget(browse_receptor_btn)
        form.addRow("Receptor:", receptor_row)

        ligand_row = QHBoxLayout()
        self.ligand_input = make_recent_combo("prepare_ligand", "Optional: SMILES string, or path to .sdf/.mol/.mol2/.pdb")
        browse_ligand_btn = QPushButton("Browse")
        browse_ligand_btn.clicked.connect(self.browse_ligand)
        ligand_row.addWidget(self.ligand_input)
        ligand_row.addWidget(browse_ligand_btn)
        form.addRow("Ligand:", ligand_row)

        ligand_name_row = QHBoxLayout()
        self.ligand_name = QLineEdit()
        self.ligand_name.setPlaceholderText("Optional: output filename for the ligand (e.g. syringic_acid)")
        ligand_name_row.addWidget(self.ligand_name)
        form.addRow("Ligand Output Name:", ligand_name_row)

        output_row = QHBoxLayout()
        self.output_dir = QLineEdit()
        self.output_dir.setPlaceholderText("Folder to save prepared .pdbqt file(s) into")
        self.output_dir.setText(memory.get_last("prepare_output_dir", ""))
        browse_output_btn = QPushButton("Choose Folder")
        browse_output_btn.clicked.connect(self.browse_output_dir)
        output_row.addWidget(self.output_dir)
        output_row.addWidget(browse_output_btn)
        form.addRow("Save To:", output_row)

        self.ai_advise_checkbox = QCheckBox("Get AI structure advisory (report-only)")
        self.ai_advise_checkbox.toggled.connect(self.toggle_ai_advisory)
        form.addRow("Optional Add-on:", self.ai_advise_checkbox)

        self.ai_backend_combo = QComboBox()
        self.ai_backend_combo.addItems(["anthropic", "gemini"])
        self.ai_backend_combo.setEnabled(False)
        form.addRow("AI Backend:", self.ai_backend_combo)

        self.ai_pdb_id_input = make_recent_combo(
            "prepare_ai_pdb_id", "Optional but recommended: PDB code (e.g. 8FK4) for RCSB/PubMed context"
        )
        self.ai_pdb_id_input.setEnabled(False)
        form.addRow("PDB ID (for AI context):", self.ai_pdb_id_input)

        self.ai_model_input = QLineEdit()
        self.ai_model_input.setPlaceholderText("Optional: override model (leave blank for the backend's default)")
        self.ai_model_input.setEnabled(False)
        form.addRow("AI Model Override:", self.ai_model_input)

        self.ai_api_key_input = QLineEdit()
        self.ai_api_key_input.setPlaceholderText("Optional: override API key (leave blank to use env var)")
        self.ai_api_key_input.setEchoMode(QLineEdit.EchoMode.Password)
        self.ai_api_key_input.setEnabled(False)
        form.addRow("AI API Key Override:", self.ai_api_key_input)

        ai_note = QLabel(
            "AI Structure Advisory is report-only - it never edits the receptor, it just writes a "
            "*_ai_advisory.md/.json alongside the prepared .pdbqt with protein identification and "
            "prep guidance grounded in RCSB/PubMed data. Requires a Receptor above, internet access, "
            "and an ANTHROPIC_API_KEY or GEMINI_API_KEY environment variable (or the override field above)."
        )
        ai_note.setWordWrap(True)

        note = QLabel(
            "Provide a receptor, a ligand, or both. Each prepared file is saved directly into "
            "the chosen folder as <name>.pdbqt - handy for building a reusable library of "
            "prepared receptors/ligands ahead of a docking run."
        )
        note.setWordWrap(True)

        layout.addLayout(form)
        layout.addWidget(note)
        layout.addWidget(ai_note)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Prepare && Save")
        self.run_btn.clicked.connect(self.run_pipeline)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_pipeline)
        self.stop_btn.setEnabled(False)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.stop_btn)
        layout.addLayout(btn_row)

        self.result_label = QLabel("")
        self.result_label.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        layout.addWidget(self.result_label)

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setFont(QFont("Consolas", 10))
        layout.addWidget(self.console)

    def browse_receptor(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Receptor File", "", RECEPTOR_FILE_FILTER)
        if file_path:
            self.receptor_input.setCurrentText(file_path)

    def browse_ligand(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Ligand File", "", LIGAND_FILE_FILTER)
        if file_path:
            self.ligand_input.setCurrentText(file_path)

    def browse_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Folder")
        if dir_path:
            self.output_dir.setText(dir_path)

    def toggle_ai_advisory(self, checked):
        self.ai_backend_combo.setEnabled(checked)
        self.ai_pdb_id_input.setEnabled(checked)
        self.ai_model_input.setEnabled(checked)
        self.ai_api_key_input.setEnabled(checked)

    def run_pipeline(self):
        receptor_val = self.receptor_input.currentText().strip()
        ligand_val = self.ligand_input.currentText().strip()
        ligand_name_val = self.ligand_name.text().strip()
        output_dir_val = self.output_dir.text().strip()

        if not receptor_val and not ligand_val:
            self.console.append("ERROR: Provide a receptor, a ligand, or both.\n")
            return
        if not output_dir_val:
            self.console.append("ERROR: Choose a folder to save the prepared file(s) into.\n")
            return

        args = ["-u", self.script_path, "--output-dir", output_dir_val]
        if receptor_val:
            args += ["--receptor", receptor_val]
        if ligand_val:
            args += ["--ligand", ligand_val]
        if ligand_name_val:
            args += ["--ligand-name", ligand_name_val]

        if self.ai_advise_checkbox.isChecked():
            if not receptor_val:
                self.console.append(
                    "NOTE: AI Structure Advisory is checked but no Receptor was given - it only runs "
                    "alongside a receptor, so this run will skip the advisory step.\n"
                )
            args += ["--ai-advise", "--backend", self.ai_backend_combo.currentText()]
            pdb_id_val = self.ai_pdb_id_input.currentText().strip()
            if pdb_id_val:
                args += ["--pdb-id", pdb_id_val]
                memory.remember("prepare_ai_pdb_id", pdb_id_val)
            model_val = self.ai_model_input.text().strip()
            if model_val:
                args += ["--model", model_val]
            api_key_val = self.ai_api_key_input.text().strip()
            if api_key_val:
                args += ["--api-key", api_key_val]

        if receptor_val:
            memory.remember("prepare_receptor", receptor_val)
        if ligand_val:
            memory.remember("prepare_ligand", ligand_val)
        memory.remember("prepare_output_dir", output_dir_val)

        self.console.clear()
        self.result_label.setText("")
        self.result_label.setStyleSheet("")
        display_args = list(args)
        if "--api-key" in display_args:
            idx = display_args.index("--api-key")
            if idx + 1 < len(display_args):
                display_args[idx + 1] = "***"
        self.console.append(f"Running: python3 {' '.join(display_args)}\n")

        self.process = QProcess(self)
        self.process.setWorkingDirectory(SCRIPT_DIR)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.finished.connect(self.handle_finished)
        self.process.start("python3", args)

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

    def handle_output(self):
        data = self.process.readAllStandardOutput().data().decode(errors="replace")
        self.console.moveCursor(QTextCursor.MoveOperation.End)
        self.console.insertPlainText(data)
        self.console.moveCursor(QTextCursor.MoveOperation.End)

    def handle_finished(self):
        full_text = self.console.toPlainText()
        saved = []
        for line in full_text.splitlines():
            if line.startswith("RECEPTOR_SAVED,"):
                saved.append(("Receptor", line.split(",", 1)[1]))
            elif line.startswith("LIGAND_SAVED,"):
                saved.append(("Ligand", line.split(",", 1)[1]))
            elif line.startswith("AI_ADVISORY_SAVED,"):
                saved.append(("AI Advisory", line.split(",", 1)[1]))

        if "PREP_COMPLETE" in full_text and saved:
            summary = "; ".join(f"{kind} -> {path}" for kind, path in saved)
            if self.ai_advise_checkbox.isChecked() and "AI advisory FAILED" in full_text:
                summary += " (AI advisory failed - see log; prep itself still succeeded)"
                self.result_label.setStyleSheet("color: orange;")
            else:
                self.result_label.setStyleSheet("color: green;")
            self.result_label.setText(f"Saved: {summary}")
        else:
            self.result_label.setText("Run finished - check log for issues")
            self.result_label.setStyleSheet("color: orange;")

        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def stop_pipeline(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.process.kill()
            self.console.append("\n--- Stopped by user ---\n")
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Discovery Lite - Docking Tools")
        self.setGeometry(150, 150, 1100, 750)

        tabs = QTabWidget()
        self.setCentralWidget(tabs)

        viewer_tab = ViewerTab()
        redock_tab = RedockValidationTab()
        dock_only_tab = DockOnlyTab(viewer_tab)
        pocket_tab = PocketDetectionTab(dock_only_tab)
        prepare_tab = PrepareStructuresTab()
        results_tab = ResultsLogTab(results_store)

        # Auto-load: as soon as a Redocking Validation run finishes, the
        # receptor + native ligand + docked mode-1 pose show up overlaid in
        # the Viewer tab, whether or not it's the one currently open.
        redock_tab.runCompleted.connect(viewer_tab.load_run_result)

        tabs.addTab(redock_tab, "Redocking Validation")
        tabs.addTab(dock_only_tab, "Docking Only")
        tabs.addTab(pocket_tab, "Pocket Detection")
        tabs.addTab(prepare_tab, "Prepare Structures")
        tabs.addTab(results_tab, "Results Log")
        tabs.addTab(viewer_tab, "Viewer")


app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec())
