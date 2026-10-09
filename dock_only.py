"""
Docking-only pipeline (no RMSD validation).

Use this when you already know your binding site (e.g. from a validated
redocking run) and just want to dock a ligand against a receptor.

Two ways to define the search box:
  1. --reference-ligand   Give a small .pdb/.sdf/.mol2 file containing just the
                          ligand/pocket-defining molecule (e.g. the native_ligand.pdb
                          produced by redock_and_validate.py). The box is auto-computed
                          around it, same logic as the validation script.
  2. --center-x/y/z + --size-x/y/z    Type the box coordinates directly (e.g. from a
                          conf.txt you already trust, or from literature).

Run with:
    python dock_only.py receptor.pdb --ligand-input ligand.sdf --reference-ligand native_ligand.pdb
    python dock_only.py receptor.pdb --ligand-input "CC(=O)Oc1ccccc1C(=O)O" \
        --center-x 12.3 --center-y 4.5 --size-x 20 --size-y 20 --size-z 20 --center-z -8.1
"""

import subprocess
import argparse
import shutil
import sys
import os


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

VINA_EXE = "vina"
OBABEL_EXE = "obabel"
PREPARE_RECEPTOR_EXE = "prepare_receptor"

BOX_PADDING = 5.0
MIN_BOX_SIZE = 20.0

_OBABEL_FORMAT_BY_EXT = {
    "sdf": "sdf", "mol": "mol", "mol2": "mol2", "ml2": "mol2",
    "pdb": "pdb", "pdbqt": "pdbqt", "smi": "smi", "cml": "cml", "mrv": "mrv",
}
_FORMATS_NEEDING_GEN3D = {"sdf", "mol", "mol2", "smi", "cml", "mrv"}


def check_tools():
    missing = []
    for exe, name in [(VINA_EXE, "vina"), (OBABEL_EXE, "obabel")]:
        if shutil.which(exe) is None and not os.path.isfile(exe):
            missing.append(name)
    if missing:
        for name in missing:
            print(f"ERROR: could not find '{name}' on PATH.")
        sys.exit(1)


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


def compute_box_from_reference(ref_path, padding, min_size):
    atoms = parse_ligand_heavy_atoms(ref_path)
    if not atoms:
        print(f"ERROR: no heavy atoms parsed from reference ligand file {ref_path}")
        sys.exit(1)
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


def guess_obabel_informat(path):
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    return _OBABEL_FORMAT_BY_EXT.get(ext, ext or "sdf")


def convert_ligand_to_pdbqt(source, out_path):
    is_smiles = not os.path.isfile(source)
    if is_smiles:
        cmd = [OBABEL_EXE, f"-:{source}", "-O", out_path, "--gen3d", "-h", "--partialcharge", "gasteiger"]
    else:
        informat = guess_obabel_informat(source)
        cmd = [OBABEL_EXE, f"-i{informat}", source, "-O", out_path, "-h", "--partialcharge", "gasteiger"]
        if informat in _FORMATS_NEEDING_GEN3D:
            cmd.insert(-3, "--gen3d")
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0 or not os.path.isfile(out_path):
        print("obabel STDERR:", result.stderr)
        sys.exit(1)


def _file_has_atoms(path):
    """A file existing isn't enough - prepare_receptor/obabel can exit 0 while
    writing an empty or atom-less file. Check it actually has ATOM/HETATM lines."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                return True
    return False


def prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt):
    if not _file_has_atoms(receptor_pdb):
        print(f"ERROR: input receptor file has no ATOM/HETATM records or doesn't exist: {receptor_pdb}")
        print("  -> Check the file path and that the file itself isn't empty/corrupted.")
        sys.exit(1)

    if shutil.which(PREPARE_RECEPTOR_EXE) or os.path.isfile(PREPARE_RECEPTOR_EXE):
        cmd = [PREPARE_RECEPTOR_EXE, "-r", receptor_pdb, "-o", receptor_pdbqt]
        print("Running:", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode == 0 and _file_has_atoms(receptor_pdbqt):
            print(f"Receptor prepared with {PREPARE_RECEPTOR_EXE} -> {receptor_pdbqt}")
            return
        if result.returncode == 0:
            print(f"WARNING: {PREPARE_RECEPTOR_EXE} exited with code 0 but wrote an empty/atom-less "
                  f"file - treating this as a failure and falling back to obabel.")
        else:
            print(f"{PREPARE_RECEPTOR_EXE} failed, falling back to obabel.")
        print("STDERR:", result.stderr)

    cmd = [OBABEL_EXE, receptor_pdb, "-O", receptor_pdbqt, "-xr", "--partialcharge", "gasteiger"]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0 or not _file_has_atoms(receptor_pdbqt):
        print("obabel STDERR:", result.stderr)
        print("ERROR: could not prepare a valid (non-empty) receptor.pdbqt automatically.")
        print("  -> Both prepare_receptor and obabel failed to produce a usable receptor file.")
        print("  -> Inspect the input PDB manually - it may have a formatting issue neither tool handled.")
        sys.exit(1)
    print(f"Receptor prepared with obabel (fallback) -> {receptor_pdbqt}")


def write_conf(conf_path, receptor, ligand, center, size, exhaustiveness, num_modes, out_pdbqt):
    with open(conf_path, "w") as f:
        f.write(f"receptor = {receptor}\n")
        f.write(f"ligand = {ligand}\n\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n\n")
        f.write(f"exhaustiveness = {exhaustiveness}\n")
        f.write(f"num_modes = {num_modes}\n\n")
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


def parse_args():
    parser = argparse.ArgumentParser(description="Docking-only pipeline (no RMSD validation)")
    parser.add_argument("receptor_input", help="Local receptor .pdb (or already-prepared .pdbqt) file")
    parser.add_argument("--ligand-input", required=True, help="Ligand file path, or a SMILES string")
    parser.add_argument("--reference-ligand", default=None,
                         help="Optional: small file (native_ligand.pdb etc.) used only to auto-compute the box")
    parser.add_argument("--center-x", type=float, default=None)
    parser.add_argument("--center-y", type=float, default=None)
    parser.add_argument("--center-z", type=float, default=None)
    parser.add_argument("--size-x", type=float, default=20.0)
    parser.add_argument("--size-y", type=float, default=20.0)
    parser.add_argument("--size-z", type=float, default=20.0)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    return parser.parse_args()


def main():
    args = parse_args()
    check_tools()

    stem = os.path.splitext(os.path.basename(args.receptor_input))[0]
    work_dir = os.path.join(os.getcwd(), "docking_runs", f"{stem}_dockonly")
    os.makedirs(work_dir, exist_ok=True)

    receptor_input = safe_local_copy(args.receptor_input, work_dir)
    ligand_input = safe_local_copy(args.ligand_input, work_dir)
    reference_ligand = safe_local_copy(args.reference_ligand, work_dir) if args.reference_ligand else None

    print("=" * 60)
    print("STEP 1: Determining search box")
    print("=" * 60)
    if reference_ligand:
        center, size = compute_box_from_reference(reference_ligand, BOX_PADDING, MIN_BOX_SIZE)
        print(f"Box auto-computed from reference ligand: {reference_ligand}")
    elif args.center_x is not None and args.center_y is not None and args.center_z is not None:
        center = (args.center_x, args.center_y, args.center_z)
        size = (args.size_x, args.size_y, args.size_z)
        print("Using manually supplied box coordinates.")
    else:
        print("ERROR: provide either --reference-ligand, or all three of --center-x/--center-y/--center-z.")
        sys.exit(1)
    print(f"Box center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"Box size:   ({size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f})")

    print("\n" + "=" * 60)
    print("STEP 2: Preparing receptor PDBQT")
    print("=" * 60)
    if receptor_input.lower().endswith(".pdbqt"):
        receptor_pdbqt = receptor_input
        print(f"Receptor already a .pdbqt, using as-is: {receptor_pdbqt}")
    else:
        receptor_pdbqt = os.path.join(work_dir, "receptor.pdbqt")
        prepare_receptor_pdbqt(receptor_input, receptor_pdbqt)

    print("\n" + "=" * 60)
    print("STEP 3: Converting ligand to PDBQT")
    print("=" * 60)
    ligand_pdbqt = os.path.join(work_dir, "ligand_for_docking.pdbqt")
    convert_ligand_to_pdbqt(ligand_input, ligand_pdbqt)

    print("\n" + "=" * 60)
    print("STEP 4: Writing conf.txt")
    print("=" * 60)
    conf_path = os.path.join(work_dir, "conf.txt")
    out_pdbqt = os.path.join(work_dir, "out.pdbqt")
    write_conf(conf_path, receptor_pdbqt, ligand_pdbqt, center, size,
               args.exhaustiveness, args.num_modes, out_pdbqt)
    print(f"conf.txt written to {conf_path}")

    print("\n" + "=" * 60)
    print("STEP 5: Running Vina")
    print("=" * 60)
    log_path = os.path.join(work_dir, "vina_log.txt")
    run_vina(conf_path, log_path)

    print("\n" + "=" * 60)
    print("STEP 6: Results")
    print("=" * 60)
    affinities = parse_vina_affinities(out_pdbqt)
    print(f"{'Mode':<6}{'Affinity (kcal/mol)':<22}")
    for i, aff in enumerate(affinities, start=1):
        print(f"{i:<6}{aff:<22.2f}")

    if affinities:
        print(f"\nDOCKING_COMPLETE: best affinity = {affinities[0]:.2f} kcal/mol")
    else:
        print("\nDOCKING_COMPLETE: no poses found - check vina_log.txt")
    print("=" * 60)


if __name__ == "__main__":
    main()