"""
Aggregates affinity + RMSD across docking runs into one persisted, sortable
log (results_log.json, next to the GUI scripts) so a whole set of runs -
e.g. every PDB ID in a validation sweep - lands in one table you can sort,
review, and export straight to CSV for a manuscript table.

ResultsLogStore is the on-disk store (list of row dicts, one per completed
run) plus a Qt signal so any open ResultsLogTab updates live as new runs
finish. ResultsLogTab is the QWidget that displays it.
"""

import csv
import json
import os
import threading
from datetime import datetime

from PyQt6.QtCore import QObject, pyqtSignal, Qt, QUrl
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QPushButton, QLabel, QFileDialog, QMessageBox, QAbstractItemView, QHeaderView
)

# (row dict key, column header). Order here = column order in the table and CSV export.
COLUMNS = [
    ("timestamp", "Timestamp"),
    ("source", "Source"),
    ("target", "PDB / Receptor"),
    ("ligand", "Ligand"),
    ("affinity", "Affinity (kcal/mol)"),
    ("rmsd", "RMSD (A)"),
    ("result", "Result"),
    ("exhaustiveness", "Exh."),
    ("num_modes", "Modes"),
    ("work_dir", "Work Dir"),
]

_NUMERIC_KEYS = {"affinity", "rmsd", "exhaustiveness", "num_modes"}


class ResultsLogStore(QObject):
    """Owns results_log.json. Shared by every tab that produces results and
    by the Results Log tab that displays them."""

    resultAdded = pyqtSignal(dict)
    logCleared = pyqtSignal()

    def __init__(self, path):
        super().__init__()
        self.path = path
        self._lock = threading.Lock()

    def load_all(self):
        if not os.path.isfile(self.path):
            return []
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                data = json.load(f)
            return data if isinstance(data, list) else []
        except (json.JSONDecodeError, OSError):
            return []

    def add_result(self, row):
        """row keys: source, target, ligand, affinity, rmsd, result,
        exhaustiveness, num_modes, work_dir. affinity/rmsd may be None."""
        row = dict(row)
        row.setdefault("timestamp", datetime.now().isoformat(timespec="seconds"))
        with self._lock:
            rows = self.load_all()
            rows.append(row)
            try:
                tmp_path = self.path + ".tmp"
                with open(tmp_path, "w", encoding="utf-8") as f:
                    json.dump(rows, f, indent=2)
                os.replace(tmp_path, self.path)
            except OSError:
                return
        self.resultAdded.emit(row)

    def clear(self):
        with self._lock:
            try:
                if os.path.isfile(self.path):
                    os.remove(self.path)
            except OSError:
                return
        self.logCleared.emit()


class ResultsLogTab(QWidget):
    """Tab: sortable log of affinity + RMSD across every completed run."""

    def __init__(self, store: ResultsLogStore):
        super().__init__()
        self.store = store
        self.store.resultAdded.connect(lambda _row: self.refresh())
        self.store.logCleared.connect(self.refresh)

        layout = QVBoxLayout()
        self.setLayout(layout)

        btn_row = QHBoxLayout()
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.refresh)
        self.export_btn = QPushButton("Export CSV...")
        self.export_btn.clicked.connect(self.export_csv)
        self.open_folder_btn = QPushButton("Open Run Folder")
        self.open_folder_btn.clicked.connect(self.open_selected_folder)
        self.clear_btn = QPushButton("Clear Log")
        self.clear_btn.clicked.connect(self.clear_log)
        btn_row.addWidget(self.refresh_btn)
        btn_row.addWidget(self.export_btn)
        btn_row.addWidget(self.open_folder_btn)
        btn_row.addWidget(self.clear_btn)
        btn_row.addStretch()
        self.count_label = QLabel("")
        btn_row.addWidget(self.count_label)
        layout.addLayout(btn_row)

        self.table = QTableWidget()
        self.table.setColumnCount(len(COLUMNS))
        self.table.setHorizontalHeaderLabels([label for _key, label in COLUMNS])
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table)

        note = QLabel(
            "Every completed Redocking Validation run is logged here automatically, with "
            "affinity + RMSD pulled straight from that run's output. Docking Only runs are "
            "logged too, with RMSD left blank (there's no crystal pose to validate an "
            "arbitrary docked ligand against). Click a column header to sort; "
            "'Export CSV' gives a table you can drop straight into a manuscript."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        self.refresh()

    def refresh(self):
        rows = self.store.load_all()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            work_dir = row.get("work_dir", "")
            for c, (key, _label) in enumerate(COLUMNS):
                value = row.get(key, "")
                item = QTableWidgetItem()
                if key in _NUMERIC_KEYS and value not in (None, ""):
                    try:
                        item.setData(Qt.ItemDataRole.DisplayRole, round(float(value), 3))
                    except (TypeError, ValueError):
                        item.setText(str(value))
                else:
                    item.setText("" if value is None else str(value))
                item.setData(Qt.ItemDataRole.UserRole, work_dir)
                self.table.setItem(r, c, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self.count_label.setText(f"{len(rows)} run(s) logged")

    def export_csv(self):
        rows = self.store.load_all()
        if not rows:
            QMessageBox.information(self, "Nothing to export", "The results log is empty.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results Log", "results_log.csv", "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow([label for _key, label in COLUMNS])
                for row in rows:
                    writer.writerow([row.get(key, "") for key, _label in COLUMNS])
        except OSError as exc:
            QMessageBox.warning(self, "Export failed", f"Could not write CSV:\n{exc}")
            return
        QMessageBox.information(self, "Exported", f"Saved {len(rows)} row(s) to:\n{file_path}")

    def open_selected_folder(self):
        selected = self.table.selectedItems()
        if not selected:
            QMessageBox.information(self, "No row selected", "Select a row first.")
            return
        work_dir = selected[0].data(Qt.ItemDataRole.UserRole)
        if not work_dir or not os.path.isdir(work_dir):
            QMessageBox.warning(self, "Folder not found", "This run's work folder isn't available on disk.")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(work_dir))

    def clear_log(self):
        if QMessageBox.question(
            self, "Clear results log",
            "This permanently deletes results_log.json. Continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        ) == QMessageBox.StandardButton.Yes:
            self.store.clear()
