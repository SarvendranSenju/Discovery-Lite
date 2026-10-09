# Discovery Lite

A desktop GUI that wraps a molecular-docking workflow (structure preparation, pocket detection, docking and redocking validation) around existing command-line tools, so a full run can be done without hand-editing config files.

> **AI-assistance disclosure.** The code in this repository was developed with AI assistance (Claude Code, Anthropic). The author specified the workflow and requirements, ran the tools on real structures, tested the application and reported bugs. The author is not the primary writer of the code. Docking results from this tool should be checked before use in any publication.

## What it does

The GUI (`redock_gui.py`) has six tabs:

| Tab | Purpose |
|---|---|
| Redocking Validation | Splits a complex PDB into receptor and native ligand, prepares both, docks with AutoDock Vina and reports redocking RMSD (PASS/FAIL against a threshold) |
| Docking Only | Docks a chosen ligand into a prepared receptor with a user-defined grid box |
| Pocket Detection | Runs P2Rank to predict binding pockets |
| Prepare Structures | Cleans a receptor and converts it to PDBQT; optional AI advisory report (report-only) |
| Results Log | Records each run (affinity, RMSD, parameters) |
| Viewer | 3D view of structures and poses |

Each tab calls a command-line worker script (`redock_and_validate.py`, `dock_only.py`, `find_pockets.py`, `prepare_structures.py`), which can also be run on its own.

## Requirements

Python 3.10+ and:

```
pip install -r requirements.txt
```

The app shells out to these tools, which are **not** included and must be installed and on your `PATH`:

- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
- [Open Babel](https://openbabel.org) (`obabel`)
- [ADFRsuite](https://ccsb.scripps.edu/adfr/) (`prepare_receptor`)
- [P2Rank](https://github.com/rdk/p2rank) (`prank`; needs a Java runtime)

Developed on Windows with WSL2.

## Run

```
python redock_gui.py
```

Or run a worker directly, e.g. redocking of PDB entry 8FK4:

```
python redock_and_validate.py 8FK4
```

## Optional AI structure advisor

`prepare_structures.py --ai-advise` can produce a report-only advisory for a receptor (heteroatom handling, protonation notes), using the Anthropic or Gemini API. It never edits your structure files. It needs `ANTHROPIC_API_KEY` or `GEMINI_API_KEY` set as an **environment variable**. Never commit an API key. See [docs/AI_ADVISOR.md](docs/AI_ADVISOR.md).

## Building a Windows executable

`build_windows.bat` freezes the GUI and workers with PyInstaller. It bundles only the Python/Qt side, not the external tools above.

## Status and limitations

- No benchmark of redocking accuracy across many complexes is included yet.
- Results depend on correct receptor and ligand preparation; always inspect poses.
- `legacy/main_prototype.py` is an early viewer-only prototype kept for history.

## Licence

See `LICENSE`. This project uses PyQt6, which is GPL-licensed, so a GPL-compatible licence is required.
