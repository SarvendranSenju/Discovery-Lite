import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QFileDialog, QTextEdit, QSplitter
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtCore import Qt


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Discovery Lite")
        self.setGeometry(100, 100, 1000, 700)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout()
        central_widget.setLayout(layout)

        # Button to open PDB file
        self.open_button = QPushButton("Open PDB File")
        self.open_button.clicked.connect(self.open_pdb_file)
        layout.addWidget(self.open_button)

        # Splitter: info box on left, 3D viewer on right
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)

        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setMaximumWidth(300)
        splitter.addWidget(self.info_box)

        self.viewer = QWebEngineView()
        splitter.addWidget(self.viewer)

        splitter.setSizes([300, 700])

    def open_pdb_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB Files (*.pdb)")
        if file_path:
            self.display_pdb_info(file_path)
            self.render_structure(file_path)

    def display_pdb_info(self, file_path):
        chains = set()
        residue_count = 0
        hetatm_count = 0

        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    chains.add(line[21])
                    residue_count += 1
                elif line.startswith("HETATM"):
                    hetatm_count += 1

        info_text = f"File: {file_path}\n"
        info_text += f"Chains found: {', '.join(sorted(chains))}\n"
        info_text += f"ATOM lines: {residue_count}\n"
        info_text += f"HETATM lines: {hetatm_count}\n"

        self.info_box.setPlainText(info_text)

    def render_structure(self, file_path):
        with open(file_path, "r") as f:
            pdb_data = f.read()

        # Escape for embedding inside JS template literal
        pdb_data_escaped = pdb_data.replace("\\", "\\\\").replace("`", "\\`")

        html = f"""
        <html>
        <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style> body {{ margin: 0; }} #viewer {{ width: 100vw; height: 100vh; }} </style>
        </head>
        <body>
        <div id="viewer"></div>
        <script>
            let viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
            let pdbData = `{pdb_data_escaped}`;
            viewer.addModel(pdbData, "pdb");
            viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
            viewer.zoomTo();
            viewer.render();
        </script>
        </body>
        </html>
        """

        self.viewer.setHtml(html)


app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec())