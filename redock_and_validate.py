"""
All-in-one redocking validation pipeline - FULLY AUTOMATIC VERSION
Complex PDB (receptor + native ligand together)  ->  auto-split  ->  auto-prepared
receptor.pdbqt + ligand.pdbqt  ->  auto-centered conf.txt  ->  run AutoDock Vina
->  RMSD validation

No PyMOL step needed. Just point PDB_INPUT at:
  - a 4-character PDB ID, e.g. "8FK4"   -> auto-downloaded from rcsb.org
  - a local .pdb file path              -> used as-is
  - a URL to a .pdb file                -> auto-downloaded

The script will:
  1. Fetch the structure (download if it's an ID/URL)
  2. Auto-detect the native ligand (the largest HETATM group that isn't water,
     a crystallization additive, or a lone ion) - used to define the grid box
     and the RMSD reference pose
  3. Split it into receptor.pdb + native_ligand.pdb
  4. Auto-center a docking grid box around the native ligand
  5. Build the ligand actually used for docking. By default this is just the
     extracted native_ligand.pdb, upgraded automatically with RCSB's "ideal"
     SDF for that ligand code when available (correct bond orders/protonation -
     PDB coordinate extraction alone can't know these). You can also override
     the docking ligand entirely with LIGAND_INPUT below, supplying:
       - a local .sdf / .mol / .mol2 / .pdb file
       - a SMILES string, e.g. "CC(=O)Oc1ccccc1C(=O)O"
       - a URL to a downloadable .sdf/.mol/.mol2/.pdb file
     (3D coordinates are generated automatically for SMILES/2D inputs)
  6. Auto-prepare receptor.pdbqt (uses ADFRsuite's `prepare_receptor` if
     available, otherwise falls back to obabel) and ligand_for_docking.pdbqt
  7. Write conf.txt
  8. Run AutoDock Vina
  9. Compute redocking RMSD per pose and PASS/FAIL against a threshold
     (only meaningful when the docked ligand is chemically the native one)

Requires: obabel and vina available on PATH (or set full paths below).
Optional: ADFRsuite's `prepare_receptor` on PATH for higher-quality receptor prep.

Run with:
    python redock_and_validate.py 8FK4
    python redock_and_validate.py C:\\path\\to\\complex.pdb
    python redock_and_validate.py                      (uses PDB_INPUT below)

GUI note: this script also accepts command-line flags so redock_gui.py can
override the constants below without editing this file:
    --ligand-input, --exhaustiveness, --num-modes, --rmsd-threshold,
    --force-ligand-resname
"""

import subprocess
import argparse
import shutil
import sys
import os
import re
import urllib.request
from math import sqrt
from collections import Counter, OrderedDict

# ============================================================
# EDIT THIS SECTION (or just pass a PDB ID / path as a command-line arg)
# ============================================================
PDB_INPUT = "8FK4"           # PDB ID (auto-downloaded), local .pdb path, or URL to a .pdb
                              # that contains BOTH the receptor and its native ligand.
                              # Overridable from the command line: python redock_and_validate.py 8FK4

WORK_DIR = None              # None = auto: ./docking_runs/<PDB_ID>/   (or set an explicit path)

VINA_EXE             = "vina"     # change to full path if not on PATH
OBABEL_EXE           = "obabel"   # change to full path if not on PATH
PREPARE_RECEPTOR_EXE = "prepare_receptor"   # optional (ADFRsuite); auto-falls back to obabel if missing

RECEPTOR_PDBQT_OVERRIDE = None    # set a path here if you already have a prepared receptor.pdbqt
                                   # and want to skip auto receptor-prep entirely

# --- what to actually dock ---
LIGAND_INPUT = None
    # None (default) = dock the native ligand auto-extracted from the complex above.
    # Or set this to any ONE of:
    #   - a local file path:  r"C:\path\to\ligand.sdf"   (.sdf/.mol/.mol2/.pdb/.pdbqt)
    #   - a SMILES string:    "CC(=O)Oc1ccccc1C(=O)O"
    #   - a URL:              "https://.../compound.sdf"
    # The grid box and the RMSD reference pose ALWAYS come from the real crystallographic
    # native-ligand coordinates regardless of this setting - only the molecule that gets
    # docked changes. If you dock something other than the native ligand, treat the RMSD
    # column as not meaningful (there's no crystal pose for a molecule that wasn't solved).
LIGAND_INPUT_FORMAT = None
    # Force the obabel input format code (e.g. "sdf", "mol2", "smi") if auto-detection from
    # the file extension / URL ever guesses wrong. Usually leave this as None.
AUTO_FETCH_IDEAL_LIGAND_SDF = True
    # Only applies when LIGAND_INPUT is None. Bond orders/protonation can't be recovered from
    # PDB coordinates alone, so before falling back to the raw extracted PDB, the script tries
    # to download RCSB's "ideal" SDF for the auto-detected ligand code (proper chemistry) and
    # docks that instead. Silently falls back to the extracted PDB if the download fails.

# --- ligand auto-detection ---
# The script picks the largest "drug-like" HETATM group as the native ligand automatically.
# If it ever picks the wrong group (rare - can happen with unusual cofactors), force it here:
FORCE_LIGAND_RESNAME = None    # e.g. "ACR"
FORCE_LIGAND_CHAIN   = None    # e.g. "A"
FORCE_LIGAND_RESSEQ  = None    # e.g. "501"
MIN_LIGAND_ATOMS     = 5       # heavy-atom floor for something to even be considered a ligand candidate

BOX_PADDING = 5.0    # Angstroms added around the ligand's bounding box on each side
MIN_BOX_SIZE = 20.0  # Angstroms - floor so the box isn't too tight even for small ligands
EXHAUSTIVENESS = 8
NUM_MODES = 9
SEED = 42            # fixed seed = reproducible runs

RMSD_PASS_THRESHOLD = 2.0
# ============================================================

RCSB_URL_TEMPLATE = "https://files.rcsb.org/download/{id}.pdb"

WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}

ION_RESNAMES = {
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO",
    "CD", "HG", "LI", "CS", "BA", "SR", "AG", "AU", "PT", "AL", "GA", "IN", "PB", "SN", "TL",
}

CRYO_BUFFER_RESNAMES = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "DMS", "ACT", "FMT",
    "TRS", "BME", "MPD", "IPA", "MOH", "EOH", "CO3", "NO3", "UNK", "UNX", "UNL", "IOD",
    "BR", "AZI", "ACY", "CIT", "TAR", "MRD", "BOG", "LDA", "P6G", "P4G", "OGA", "BU3",
    "SIN", "MES", "IMD", "EPE", "BEZ", "PLM", "OLA", "MYR", "PLC", "GSH", "NH4", "OXL",
}

AD_TYPE_TO_ELEMENT = {
    "C": "C", "A": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
    "S": "S", "SA": "S", "P": "P", "F": "F", "CL": "CL", "BR": "BR",
    "I": "I", "MG": "MG", "CA": "CA", "MN": "MN", "FE": "FE", "ZN": "ZN",
    "H": None, "HD": None,
}


def safe_local_copy(path, work_dir):
    """Some older command-line tools (e.g. ADFRsuite's prepare_receptor) mishandle
    spaces in file paths even when properly quoted by the caller. If the given path
    is a local file with a space in it, copy it to work_dir under a space-free name
    and return that new path; otherwise return the path unchanged."""
    if os.path.isfile(path) and " " in path:
        safe_name = os.path.basename(path).replace(" ", "_")
        safe_path = os.path.join(work_dir, safe_name)
        shutil.copyfile(path, safe_path)
        print(f"NOTE: input path contains spaces, which breaks some older CLI tools. "
              f"Copied to a space-free path: {safe_path}")
        return safe_path
    return path


def check_tools():
    missing = []
    for exe, name in [(VINA_EXE, "vina"), (OBABEL_EXE, "obabel")]:
        if shutil.which(exe) is None and not os.path.isfile(exe):
            missing.append((exe, name))
    if missing:
        for exe, name in missing:
            print(f"ERROR: could not find '{exe}' on PATH or as a file.")
            print(f"  -> Run 'where {name}' (Windows) or 'which {name}' (Mac/Linux) to locate it, "
                  f"then set {name.upper()}_EXE to the full path at the top of this script.")
        sys.exit(1)
    if shutil.which(PREPARE_RECEPTOR_EXE) is None and not os.path.isfile(PREPARE_RECEPTOR_EXE):
        print(f"NOTE: optional tool '{PREPARE_RECEPTOR_EXE}' (ADFRsuite) not found - receptor prep will "
              f"fall back to obabel, which is a reasonable but slightly less rigorous option.")


# ---------------------------------------------------------------------------
# STEP 1: resolve input (download if it's a PDB ID or URL)
# ---------------------------------------------------------------------------

def resolve_input_pdb(pdb_input, work_dir):
    if os.path.isfile(pdb_input):
        return safe_local_copy(pdb_input, work_dir)

    if pdb_input.lower().startswith(("http://", "https://")):
        url = pdb_input
        out_path = os.path.join(work_dir, os.path.basename(url.split("?")[0]) or "downloaded.pdb")
    elif re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_input):
        pdb_id = pdb_input.upper()
        url = RCSB_URL_TEMPLATE.format(id=pdb_id)
        out_path = os.path.join(work_dir, f"{pdb_id}.pdb")
    else:
        print(f"ERROR: '{pdb_input}' is not an existing file, a URL, or a 4-character PDB ID.")
        sys.exit(1)

    print(f"Downloading structure from {url} ...")
    try:
        return download_file(url, work_dir, out_path=out_path)
    except Exception as e:
        print(f"ERROR: failed to download '{url}': {e}")
        print("  -> Check your internet connection, or download the PDB manually and point PDB_INPUT "
              "at the local file instead.")
        sys.exit(1)


def download_file(url, work_dir, out_path=None):
    """Generic file download helper, reused for the complex PDB, the optional ideal-SDF
    fetch, and any URL passed via LIGAND_INPUT. Raises on failure - callers decide how to
    handle that (hard exit vs. silent fallback)."""
    if out_path is None:
        out_path = os.path.join(work_dir, os.path.basename(url.split("?")[0]) or "downloaded_file")
    req = urllib.request.Request(url, headers={"User-Agent": "redock-pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp, open(out_path, "wb") as f:
        f.write(resp.read())
    print(f"Saved to {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# STEP 2: parse the complex, auto-detect the native ligand, split into
#          receptor.pdb + native_ligand.pdb
# ---------------------------------------------------------------------------

def parse_pdb_records(path):
    """Parse ATOM/HETATM lines from the first MODEL only (handles NMR ensembles),
    skipping alternate conformers other than the primary one."""
    records = []
    model_count = 0
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                model_count += 1
                if model_count >= 1:
                    break
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54:
                continue
            alt_loc = line[16].strip()
            if alt_loc not in ("", "A"):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            elem = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not elem:
                elem = "".join(c for c in atom_name if c.isalpha())[:1].upper()
            records.append({
                "record_type": line[:6].strip(),
                "atom_name": atom_name,
                "res_name": line[17:20].strip(),
                "chain_id": line[21].strip() or " ",
                "res_seq": line[22:26].strip(),
                "icode": line[26].strip(),
                "element": elem,
                "x": x, "y": y, "z": z,
                "raw_line": line if line.endswith("\n") else line + "\n",
            })
    return records


def classify_and_select_ligand(records, force_resname=None, force_chain=None,
                                force_resseq=None, min_atoms=5):
    groups = OrderedDict()
    for r in records:
        if r["record_type"] != "HETATM" or r["element"] == "H":
            continue
        if r["res_name"] in WATER_RESNAMES:
            continue
        key = (r["chain_id"], r["res_seq"], r["icode"], r["res_name"])
        groups.setdefault(key, []).append(r)

    print("HETATM groups found (waters excluded):")
    for key, atoms in groups.items():
        chain, resseq, icode, resname = key
        tag = ""
        if resname in ION_RESNAMES:
            tag = "  [ion - kept in receptor, not a ligand candidate]"
        elif resname in CRYO_BUFFER_RESNAMES:
            tag = "  [buffer/cryoprotectant - excluded]"
        elif len(atoms) < min_atoms:
            tag = f"  [only {len(atoms)} heavy atoms - excluded]"
        print(f"  {resname:<5} chain {chain} resSeq {resseq}{icode:<2} - {len(atoms):>3} heavy atoms{tag}")

    if force_resname:
        for key, atoms in groups.items():
            chain, resseq, icode, resname = key
            if resname != force_resname:
                continue
            if force_chain and chain != force_chain:
                continue
            if force_resseq and str(resseq) != str(force_resseq):
                continue
            return key, atoms
        print(f"ERROR: forced ligand resname={force_resname} chain={force_chain} "
              f"resSeq={force_resseq} not found among the HETATM groups above.")
        sys.exit(1)

    candidates = {
        key: atoms for key, atoms in groups.items()
        if key[3] not in ION_RESNAMES and key[3] not in CRYO_BUFFER_RESNAMES and len(atoms) >= min_atoms
    }
    if not candidates:
        print("ERROR: no ligand candidate found automatically.")
        print("  -> The ligand may be unusually small, or excluded by the buffer/ion list above.")
        print("  -> Set FORCE_LIGAND_RESNAME (and optionally FORCE_LIGAND_CHAIN/FORCE_LIGAND_RESSEQ) "
              "at the top of this script and re-run.")
        sys.exit(1)

    best_key = max(candidates, key=lambda k: len(candidates[k]))
    if len(candidates) > 1:
        print(f"\nMultiple ligand candidates found - auto-selected the largest: "
              f"{best_key[3]} chain {best_key[0]} resSeq {best_key[1]}{best_key[2]} "
              f"({len(candidates[best_key])} atoms).")
        print("  If that's wrong, set FORCE_LIGAND_RESNAME / FORCE_LIGAND_CHAIN / FORCE_LIGAND_RESSEQ "
              "at the top of this script.")
    return best_key, candidates[best_key]


def write_native_ligand_pdb(ligand_atoms, out_path):
    with open(out_path, "w") as f:
        for a in ligand_atoms:
            f.write(a["raw_line"])
        f.write("END\n")


def write_receptor_pdb(records, ligand_key, out_path):
    """Protein + any cofactors/ions, minus water and minus the chosen native ligand."""
    with open(out_path, "w") as f:
        for r in records:
            if r["record_type"] == "HETATM":
                if r["res_name"] in WATER_RESNAMES:
                    continue
                key = (r["chain_id"], r["res_seq"], r["icode"], r["res_name"])
                if key == ligand_key:
                    continue
            f.write(r["raw_line"])
        f.write("END\n")


# ---------------------------------------------------------------------------
# Ligand heavy-atom parsing (used for the box + RMSD, reads back native_ligand.pdb)
# ---------------------------------------------------------------------------

def parse_ligand_heavy_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")):
                name = line[12:16].strip()
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError:
                    continue
                elem = line[76:78].strip().upper() if len(line) >= 78 else ""
                if not elem:
                    elem = "".join(c for c in name if c.isalpha())[:1].upper()
                if elem == "H":
                    continue
                atoms.append({"element": elem, "x": x, "y": y, "z": z})
    return atoms


def compute_box(atoms, padding, min_size):
    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]
    center = ((max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, (max(zs) + min(zs)) / 2)
    size = (
        max(max(xs) - min(xs) + 2 * padding, min_size),
        max(max(ys) - min(ys) + 2 * padding, min_size),
        max(max(zs) - min(zs) + 2 * padding, min_size),
    )
    return center, size


# ---------------------------------------------------------------------------
# STEP 4/5: resolve + convert the docking ligand (native / sdf / smiles / url),
# then auto receptor prep
# ---------------------------------------------------------------------------

IDEAL_SDF_URL_TEMPLATE = "https://files.rcsb.org/ligands/download/{code}_ideal.sdf"

_OBABEL_FORMAT_BY_EXT = {
    "sdf": "sdf", "mol": "mol", "mol2": "mol2", "ml2": "mol2",
    "pdb": "pdb", "pdbqt": "pdbqt", "smi": "smi", "cml": "cml", "mrv": "mrv",
}
_FORMATS_NEEDING_GEN3D = {"sdf", "mol", "mol2", "smi", "cml", "mrv"}


def _looks_like_sdf(path):
    try:
        with open(path, "r", errors="ignore") as f:
            head = f.read(2000)
        return ("$$$$" in head) or ("V2000" in head) or ("V3000" in head)
    except OSError:
        return False


def guess_obabel_informat(path, forced_format=None):
    if forced_format:
        return forced_format
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    return _OBABEL_FORMAT_BY_EXT.get(ext, ext or "sdf")


def resolve_ligand_source(work_dir, ligand_key, native_ligand_pdb):
    """Decide what to actually dock. Returns (source, kind, is_custom) where
    kind is 'file' or 'smiles', and is_custom flags whether the user explicitly
    overrode the native ligand (used later to caveat the RMSD column)."""
    if LIGAND_INPUT:
        if os.path.isfile(LIGAND_INPUT):
            safe_path = safe_local_copy(LIGAND_INPUT, work_dir)
            print(f"Using local ligand file: {safe_path}")
            return safe_path, "file", True
        if LIGAND_INPUT.lower().startswith(("http://", "https://")):
            print(f"Downloading ligand from {LIGAND_INPUT} ...")
            try:
                downloaded = download_file(LIGAND_INPUT, work_dir)
            except Exception as e:
                print(f"ERROR: failed to download LIGAND_INPUT '{LIGAND_INPUT}': {e}")
                sys.exit(1)
            return downloaded, "file", True
        # Not a file, not a URL -> treat as a raw SMILES string
        print(f"Using SMILES ligand: {LIGAND_INPUT}")
        return LIGAND_INPUT, "smiles", True

    # Default: dock the native ligand extracted from the complex, optionally upgraded
    # with RCSB's ideal SDF (proper bond orders/protonation) for the same ligand code.
    if AUTO_FETCH_IDEAL_LIGAND_SDF and ligand_key:
        resname = ligand_key[3]
        url = IDEAL_SDF_URL_TEMPLATE.format(code=resname)
        ideal_path = os.path.join(work_dir, f"{resname}_ideal.sdf")
        print(f"Attempting to fetch RCSB ideal SDF for ligand code '{resname}' "
              f"(gives correct bond orders/protonation vs. PDB-coordinate guessing)...")
        try:
            download_file(url, work_dir, out_path=ideal_path)
            if os.path.getsize(ideal_path) > 0 and _looks_like_sdf(ideal_path):
                print(f"Using RCSB ideal SDF for docking: {ideal_path}")
                return ideal_path, "file", False
            print("Downloaded file doesn't look like a valid SDF - falling back to the extracted PDB.")
        except Exception as e:
            print(f"Ideal SDF not available for '{resname}' ({e}) - falling back to the extracted PDB.")

    return native_ligand_pdb, "file", False


def convert_ligand_to_pdbqt(source, kind, out_path, forced_format=None):
    if kind == "smiles":
        cmd = [OBABEL_EXE, f"-:{source}", "-O", out_path, "--gen3d", "-h", "--partialcharge", "gasteiger"]
    else:
        informat = guess_obabel_informat(source, forced_format)
        cmd = [OBABEL_EXE, f"-i{informat}", source, "-O", out_path, "-h", "--partialcharge", "gasteiger"]
        if informat in _FORMATS_NEEDING_GEN3D:
            cmd.insert(-3, "--gen3d")  # before -h/--partialcharge/gasteiger, ensures real 3D coords
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0 or not os.path.isfile(out_path):
        print("obabel STDERR:", result.stderr)
        if kind == "smiles":
            print("  -> Check that the SMILES string is valid.")
        sys.exit(1)


def prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt):
    prep_exe = PREPARE_RECEPTOR_EXE
    if shutil.which(prep_exe) or os.path.isfile(prep_exe):
        cmd = [prep_exe, "-r", receptor_pdb, "-o", receptor_pdbqt]
        print("Running:", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode == 0 and os.path.isfile(receptor_pdbqt):
            print(f"Receptor prepared with {prep_exe} -> {receptor_pdbqt}")
            return
        print(f"{prep_exe} failed or produced no output, falling back to obabel.")
        print("STDERR:", result.stderr)

    cmd = [OBABEL_EXE, receptor_pdb, "-O", receptor_pdbqt, "-xr", "--partialcharge", "gasteiger"]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0 or not os.path.isfile(receptor_pdbqt):
        print("obabel STDERR:", result.stderr)
        print("ERROR: could not prepare receptor.pdbqt automatically.")
        print("  -> Install ADFRsuite and make sure 'prepare_receptor' is on PATH for more reliable "
              "results, or prepare the receptor yourself and set RECEPTOR_PDBQT_OVERRIDE at the top "
              "of this script.")
        sys.exit(1)
    print(f"Receptor prepared with obabel (fallback) -> {receptor_pdbqt}")
    print("NOTE: obabel receptor prep usually works fine, but ADFRsuite's prepare_receptor gives more "
          "reliable AutoDock atom typing - worth installing if redocking results look off.")


def write_conf(conf_path, receptor, ligand, center, size, out_pdbqt):
    with open(conf_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"ligand = {ligand}\n\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n\n")
        f.write(f"exhaustiveness = {EXHAUSTIVENESS}\n")
        f.write(f"num_modes = {NUM_MODES}\n")
        f.write(f"seed = {SEED}\n\n")
        f.write(f"out = {out_pdbqt}\n")


def run_vina(conf_path, log_path):
    cmd = [VINA_EXE, "--config", conf_path]
    print("Running:", " ".join(cmd))
    with open(log_path, "w") as logf:
        result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        print(f"Vina exited with an error - check {log_path}")
        sys.exit(1)
    print(f"Vina finished. Log written to {log_path}")


def parse_vina_affinities(pdbqt_path):
    affinities = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                affinities.append(float(parts[3]))
    return affinities


def parse_docked_models(path):
    models, current = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current = []
            elif line.startswith("ENDMDL"):
                models.append(current)
            elif line.startswith(("ATOM", "HETATM")):
                tokens = line.split()
                if len(tokens) < 9:
                    continue
                adtype = tokens[-1].upper()
                elem = AD_TYPE_TO_ELEMENT.get(adtype)
                if elem is None:
                    continue
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError:
                    continue
                current.append({"element": elem, "x": x, "y": y, "z": z})
    if current and current not in models:
        models.append(current)
    return models


def nearest_neighbor_rmsd(native, docked, max_dist=10.0):
    used = [False] * len(docked)
    sq_dists, unmatched = [], 0
    for na in native:
        best_j, best_d2 = None, None
        for j, da in enumerate(docked):
            if used[j] or da["element"] != na["element"]:
                continue
            d2 = (na["x"] - da["x"]) ** 2 + (na["y"] - da["y"]) ** 2 + (na["z"] - da["z"]) ** 2
            if best_d2 is None or d2 < best_d2:
                best_d2, best_j = d2, j
        if best_j is not None and sqrt(best_d2) <= max_dist:
            used[best_j] = True
            sq_dists.append(best_d2)
        else:
            unmatched += 1
    if not sq_dists:
        return None, 0, unmatched
    return sqrt(sum(sq_dists) / len(sq_dists)), len(sq_dists), unmatched


def parse_args():
    parser = argparse.ArgumentParser(description="Automated redocking validation pipeline")
    parser.add_argument("pdb_input", nargs="?", default=None,
                         help="PDB ID, local .pdb path, or URL")
    parser.add_argument("--ligand-input", default=None,
                         help="Override docking ligand: file path, SMILES, or URL")
    parser.add_argument("--exhaustiveness", type=int, default=None)
    parser.add_argument("--num-modes", type=int, default=None)
    parser.add_argument("--rmsd-threshold", type=float, default=None)
    parser.add_argument("--force-ligand-resname", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    pdb_input = args.pdb_input or PDB_INPUT

    global LIGAND_INPUT, EXHAUSTIVENESS, NUM_MODES, RMSD_PASS_THRESHOLD, FORCE_LIGAND_RESNAME
    if args.ligand_input is not None:
        LIGAND_INPUT = args.ligand_input
    if args.exhaustiveness is not None:
        EXHAUSTIVENESS = args.exhaustiveness
    if args.num_modes is not None:
        NUM_MODES = args.num_modes
    if args.rmsd_threshold is not None:
        RMSD_PASS_THRESHOLD = args.rmsd_threshold
    if args.force_ligand_resname is not None:
        FORCE_LIGAND_RESNAME = args.force_ligand_resname

    work_dir = WORK_DIR
    if not work_dir:
        stem = (os.path.splitext(os.path.basename(pdb_input))[0]
                if os.path.isfile(pdb_input) else re.sub(r"\W+", "_", pdb_input))
        work_dir = os.path.join(os.getcwd(), "docking_runs", stem)
    os.makedirs(work_dir, exist_ok=True)

    check_tools()

    print("=" * 60)
    print("STEP 1: Resolving input structure")
    print("=" * 60)
    pdb_path = resolve_input_pdb(pdb_input, work_dir)

    print("\n" + "=" * 60)
    print("STEP 2: Parsing structure & auto-detecting the native ligand")
    print("=" * 60)
    records = parse_pdb_records(pdb_path)
    if not records:
        print(f"ERROR: no ATOM/HETATM records parsed from {pdb_path}")
        sys.exit(1)
    ligand_key, ligand_group_atoms = classify_and_select_ligand(
        records, FORCE_LIGAND_RESNAME, FORCE_LIGAND_CHAIN, FORCE_LIGAND_RESSEQ, MIN_LIGAND_ATOMS
    )
    chain, resseq, icode, resname = ligand_key
    print(f"\nSelected native ligand: {resname} chain {chain} resSeq {resseq}{icode} "
          f"({len(ligand_group_atoms)} heavy atoms)")

    native_ligand_pdb = os.path.join(work_dir, "native_ligand.pdb")
    receptor_pdb = os.path.join(work_dir, "receptor.pdb")
    write_native_ligand_pdb(ligand_group_atoms, native_ligand_pdb)
    write_receptor_pdb(records, ligand_key, receptor_pdb)
    print(f"Wrote {native_ligand_pdb}")
    print(f"Wrote {receptor_pdb}")

    print("\n" + "=" * 60)
    print("STEP 3: Parsing native ligand heavy atoms & computing grid box")
    print("=" * 60)
    native_atoms = parse_ligand_heavy_atoms(native_ligand_pdb)
    if not native_atoms:
        print(f"ERROR: no heavy atoms parsed from {native_ligand_pdb}")
        sys.exit(1)
    print(f"Native heavy atoms: {len(native_atoms)}")
    print("Element counts:", dict(Counter(a["element"] for a in native_atoms)))

    center, size = compute_box(native_atoms, BOX_PADDING, MIN_BOX_SIZE)
    print(f"Box center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"Box size:   ({size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f})")

    print("\n" + "=" * 60)
    print("STEP 4: Resolving & converting the docking ligand to PDBQT")
    print("=" * 60)
    ligand_source, ligand_kind, is_custom_ligand = resolve_ligand_source(
        work_dir, ligand_key, native_ligand_pdb
    )
    if is_custom_ligand:
        print("NOTE: docking a user-supplied ligand, not (necessarily) the native crystallographic one.")
        print("      The RMSD column in STEP 8 only means something if this is chemically the same")
        print("      molecule as the co-crystallized ligand - ignore it if you're docking something else.")
    ligand_pdbqt = os.path.join(work_dir, "ligand_for_docking.pdbqt")
    convert_ligand_to_pdbqt(ligand_source, ligand_kind, ligand_pdbqt, LIGAND_INPUT_FORMAT)

    print("\n" + "=" * 60)
    print("STEP 5: Preparing receptor PDBQT")
    print("=" * 60)
    if RECEPTOR_PDBQT_OVERRIDE:
        receptor_pdbqt = safe_local_copy(RECEPTOR_PDBQT_OVERRIDE, work_dir)
        print(f"Using pre-prepared receptor override: {receptor_pdbqt}")
        if not os.path.isfile(receptor_pdbqt):
            print(f"ERROR: RECEPTOR_PDBQT_OVERRIDE is set but the file was not found: {receptor_pdbqt}")
            sys.exit(1)
    else:
        receptor_pdbqt = os.path.join(work_dir, "receptor.pdbqt")
        prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

    print("\n" + "=" * 60)
    print("STEP 6: Writing conf.txt")
    print("=" * 60)
    conf_path = os.path.join(work_dir, "conf.txt")
    out_pdbqt = os.path.join(work_dir, "out.pdbqt")
    write_conf(conf_path, receptor_pdbqt, ligand_pdbqt, center, size, out_pdbqt)
    print(f"conf.txt written to {conf_path}")

    print("\n" + "=" * 60)
    print("STEP 7: Running Vina")
    print("=" * 60)
    log_path = os.path.join(work_dir, "vina_log.txt")
    run_vina(conf_path, log_path)

    print("\n" + "=" * 60)
    print("STEP 8: Validating redocking RMSD")
    print("=" * 60)
    affinities = parse_vina_affinities(out_pdbqt)
    models = parse_docked_models(out_pdbqt)

    print(f"{'Mode':<6}{'Affinity':<12}{'HeavyAtoms':<12}{'Matched':<10}{'Unmatched':<12}{'RMSD (A)':<10}")
    mode1_rmsd = None
    mode1_affinity = affinities[0] if affinities else None
    for i, model in enumerate(models, start=1):
        rmsd, n_matched, n_unmatched = nearest_neighbor_rmsd(native_atoms, model)
        rmsd_str = f"{rmsd:.3f}" if rmsd is not None else "N/A"
        aff = affinities[i - 1] if i - 1 < len(affinities) else float("nan")
        print(f"{i:<6}{aff:<12.2f}{len(model):<12}{n_matched:<10}{n_unmatched:<12}{rmsd_str:<10}")
        if i == 1:
            mode1_rmsd = rmsd

    print("\n" + "=" * 60)
    if mode1_rmsd is not None and mode1_rmsd <= RMSD_PASS_THRESHOLD:
        pass_fail = "PASS"
        print(f"PASS: Mode 1 RMSD = {mode1_rmsd:.3f} A (<= {RMSD_PASS_THRESHOLD} A threshold)")
        print("Grid box is validated - safe to proceed with your ligand panel using this conf.txt.")
    elif mode1_rmsd is not None:
        pass_fail = "FAIL"
        print(f"FAIL: Mode 1 RMSD = {mode1_rmsd:.3f} A (> {RMSD_PASS_THRESHOLD} A threshold)")
        print("Grid box likely needs re-centering - check the center/size values above against")
        print("the true binding site before trusting any comparative docking results.")
    else:
        pass_fail = "FAIL"
        print("FAIL: could not compute RMSD for mode 1 (no atoms matched).")
    print("=" * 60)

    # Machine-readable summary line for redock_gui.py's Results Log tab (same
    # convention as POCKET_ROW / RECEPTOR_SAVED elsewhere in this project).
    # Free-text fields are comma-sanitized since this line is comma-delimited.
    ligand_label = (LIGAND_INPUT if is_custom_ligand else f"native ({resname})")
    result_fields = [
        "RESULT_ROW",
        str(pdb_input).replace(",", ";"),
        str(ligand_label).replace(",", ";"),
        f"{mode1_affinity:.3f}" if mode1_affinity is not None else "",
        f"{mode1_rmsd:.3f}" if mode1_rmsd is not None else "",
        pass_fail,
        str(EXHAUSTIVENESS),
        str(NUM_MODES),
        str(work_dir).replace(",", ";"),
    ]
    print(",".join(result_fields))


if __name__ == "__main__":
    main()
