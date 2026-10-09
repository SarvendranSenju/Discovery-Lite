"""
Prepare a receptor and/or ligand as PDBQT files, independent of running any
docking. Useful for building up a library of pre-prepared structures (e.g. one
receptor.pdbqt reused across many ligands), or for prepping a ligand panel
ahead of time.

At least one of --receptor / --ligand must be given. Each is optional so this
script can be used for receptor-only or ligand-only prep.

Run with:
    python prepare_structures.py --receptor receptor.pdb --output-dir ./prepared
    python prepare_structures.py --ligand "CC(=O)Oc1ccccc1C(=O)O" --ligand-name aspirin --output-dir ./prepared
    python prepare_structures.py --receptor receptor.pdb --ligand ligand.sdf --output-dir ./prepared

Optional AI structure advisory (report-only, needs ANTHROPIC_API_KEY + internet):
    python prepare_structures.py --receptor receptor.pdb --output-dir ./prepared \\
        --ai-advise --pdb-id 8FK4
See structure_advisor.py for what this does and does not do (it never
auto-edits the receptor - it writes a report alongside the prepared files).
"""

import argparse
import os
import shutil
import subprocess
import sys

try:
    from structure_advisor import run_advisor, format_report
    _ADVISOR_AVAILABLE = True
except ImportError:
    _ADVISOR_AVAILABLE = False

VINA_EXE = "vina"
OBABEL_EXE = "obabel"
PREPARE_RECEPTOR_EXE = "prepare_receptor"

_OBABEL_FORMAT_BY_EXT = {
    "sdf": "sdf", "mol": "mol", "mol2": "mol2", "ml2": "mol2",
    "pdb": "pdb", "pdbqt": "pdbqt", "smi": "smi", "cml": "cml", "mrv": "mrv",
}
_FORMATS_NEEDING_GEN3D = {"sdf", "mol", "mol2", "smi", "cml", "mrv"}


def safe_local_copy(path, work_dir):
    """Some older command-line tools (e.g. ADFRsuite's prepare_receptor) mishandle
    spaces in file paths even when properly quoted by the caller. If the given path
    is a local file with a space in it, copy it to work_dir under a space-free name
    and return that new path; otherwise return the path unchanged."""
    if path and os.path.isfile(path) and " " in path:
        safe_name = os.path.basename(path).replace(" ", "_")
        safe_path = os.path.join(work_dir, safe_name)
        shutil.copyfile(path, safe_path)
        print(f"NOTE: input path contains spaces, which breaks some older CLI tools. "
              f"Copied to a space-free path: {safe_path}")
        return safe_path
    return path


def check_tools(need_receptor, need_ligand):
    if (need_receptor or need_ligand) and shutil.which(OBABEL_EXE) is None and not os.path.isfile(OBABEL_EXE):
        print("ERROR: could not find 'obabel' on PATH.")
        sys.exit(1)


def _file_has_atoms(path):
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
        sys.exit(1)
    print(f"Receptor prepared with obabel (fallback) -> {receptor_pdbqt}")


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
        if is_smiles:
            print("  -> Check that the SMILES string is valid.")
        sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare a receptor and/or ligand as PDBQT files")
    parser.add_argument("--receptor", default=None, help="Local receptor .pdb file to prepare")
    parser.add_argument("--receptor-name", default=None,
                         help="Base filename to use for the prepared receptor (default: derived from --receptor)")
    parser.add_argument("--ligand", default=None, help="Ligand file path, or a SMILES string")
    parser.add_argument("--ligand-name", default=None,
                         help="Base filename to use for the prepared ligand (default: derived from --ligand, "
                              "or 'ligand' for a raw SMILES string)")
    parser.add_argument("--output-dir", required=True, help="Folder to save the prepared file(s) into")
    parser.add_argument("--ai-advise", action="store_true",
                         help="Fetch RCSB/PubMed context and ask an LLM for protein-aware prep guidance "
                              "(report-only - written alongside the prepared receptor, never auto-applied). "
                              "Needs an API key (see --backend) and internet access.")
    parser.add_argument("--backend", choices=["anthropic", "gemini"], default="anthropic",
                         help="Which API to use for --ai-advise. 'anthropic' (default) needs "
                              "ANTHROPIC_API_KEY. 'gemini' needs GEMINI_API_KEY and has an ongoing "
                              "free tier - get a key free at https://aistudio.google.com/apikey")
    parser.add_argument("--pdb-id", default=None,
                         help="4-character PDB ID for --ai-advise to pull RCSB/PubMed metadata for. "
                              "If omitted, --ai-advise still runs but only on your local receptor file's "
                              "own HETATM content, which is a much weaker basis for identification.")
    parser.add_argument("--model", default=None,
                         help="Model string for --ai-advise (default set per-backend in structure_advisor.py)")
    parser.add_argument("--api-key", default=None,
                         help="API key for --ai-advise (default: ANTHROPIC_API_KEY or GEMINI_API_KEY "
                              "env var, per --backend)")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.receptor and not args.ligand:
        print("ERROR: provide at least one of --receptor or --ligand.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    check_tools(need_receptor=bool(args.receptor), need_ligand=bool(args.ligand))

    # Scratch space for space-free copies of inputs; kept separate from the
    # user's chosen output folder so we never leave temp copies mixed in there.
    scratch_dir = os.path.join(args.output_dir, ".prep_scratch")
    os.makedirs(scratch_dir, exist_ok=True)

    if args.receptor:
        print("=" * 60)
        print("Preparing receptor")
        print("=" * 60)
        receptor_input = safe_local_copy(args.receptor, scratch_dir)
        if receptor_input.lower().endswith(".pdbqt"):
            print("Input is already a .pdbqt - copying as-is instead of re-preparing.")
            stem = args.receptor_name or os.path.splitext(os.path.basename(args.receptor))[0]
            out_path = os.path.join(args.output_dir, f"{stem}.pdbqt")
            shutil.copyfile(receptor_input, out_path)
        else:
            stem = args.receptor_name or os.path.splitext(os.path.basename(args.receptor))[0]
            out_path = os.path.join(args.output_dir, f"{stem}.pdbqt")
            prepare_receptor_pdbqt(receptor_input, out_path)
        print(f"RECEPTOR_SAVED,{out_path}")

        if args.ai_advise:
            print("\n" + "=" * 60)
            print("AI structure advisory (report-only)")
            print("=" * 60)
            if not _ADVISOR_AVAILABLE:
                print("NOTE: structure_advisor.py not found next to this script - skipping.")
            elif not args.pdb_id:
                print("NOTE: no --pdb-id given, so this can only reason from the local file's own "
                      "HETATM content, not RCSB/PubMed context. Pass --pdb-id for a much stronger "
                      "identification if this receptor came from a known PDB entry.")
            advisor_kwargs = {}
            if args.model:
                advisor_kwargs["model"] = args.model
            if _ADVISOR_AVAILABLE:
                result = run_advisor(
                    pdb_id=args.pdb_id, pdb_path=receptor_input, api_key=args.api_key,
                    backend=args.backend, **advisor_kwargs
                )
                report_text = format_report(result)
                print(report_text)
                if result["status"] == "ok":
                    report_path = os.path.join(args.output_dir, f"{stem}_ai_advisory.md")
                    json_path = os.path.join(args.output_dir, f"{stem}_ai_advisory.json")
                    with open(report_path, "w") as f:
                        f.write(report_text)
                    with open(json_path, "w") as f:
                        import json as _json
                        _json.dump(result, f, indent=2)
                    print(f"\nAI_ADVISORY_SAVED,{report_path}")

    if args.ligand:
        print("\n" + "=" * 60)
        print("Preparing ligand")
        print("=" * 60)
        ligand_input = safe_local_copy(args.ligand, scratch_dir)
        is_smiles = not os.path.isfile(ligand_input)
        if args.ligand_name:
            stem = args.ligand_name
        elif is_smiles:
            stem = "ligand"
        else:
            stem = os.path.splitext(os.path.basename(args.ligand))[0]
        out_path = os.path.join(args.output_dir, f"{stem}.pdbqt")
        convert_ligand_to_pdbqt(ligand_input, out_path)
        print(f"LIGAND_SAVED,{out_path}")

    shutil.rmtree(scratch_dir, ignore_errors=True)
    print("\nPREP_COMPLETE")


if __name__ == "__main__":
    main()
