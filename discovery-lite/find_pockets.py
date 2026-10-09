"""
Runs P2Rank on a receptor PDB file and prints predicted binding pockets,
ranked by score, with their center coordinates - ready to feed into
dock_only.py's --center-x/--center-y/--center-z box definition.

Requires: 'prank' (P2Rank) on PATH.

Run with:
    python find_pockets.py receptor.pdb
"""

import argparse
import csv
import os
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Predict ligand binding pockets with P2Rank")
    parser.add_argument("pdb_path", help="Path to a receptor .pdb file")
    parser.add_argument("--output-dir", default=None,
                         help="Where P2Rank writes its output (default: docking_runs/<name>_pockets)")
    return parser.parse_args()


def run_p2rank(pdb_path, output_dir):
    cmd = ["prank", "predict", "-f", pdb_path, "-o", output_dir, "-visualizations", "0"]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("P2Rank STDERR:", result.stderr)
        print("ERROR: P2Rank failed to run. Check that 'prank' is on PATH and the input PDB is valid.")
        sys.exit(1)


def find_predictions_csv(output_dir, pdb_path):
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    expected = os.path.join(output_dir, f"{stem}.pdb_predictions.csv")
    if os.path.isfile(expected):
        return expected
    # Fallback: search the output dir for any *_predictions.csv, in case the
    # naming convention differs slightly between P2Rank versions.
    for fname in os.listdir(output_dir):
        if fname.endswith("_predictions.csv"):
            return os.path.join(output_dir, fname)
    print(f"ERROR: could not find a *_predictions.csv file in {output_dir}")
    sys.exit(1)


def print_pockets(csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        # P2Rank's CSV headers sometimes have leading spaces - normalize them.
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        rows = list(reader)

    if not rows:
        print("No pockets predicted for this structure.")
        return

    print(f"\n{'Rank':<6}{'Score':<10}{'Probability':<14}{'Center (x, y, z)':<32}{'Residues':<10}")
    for row in rows:
        rank = row.get("rank", "").strip()
        score = row.get("score", "").strip()
        prob = row.get("probability", "").strip()
        cx = row.get("center_x", "").strip()
        cy = row.get("center_y", "").strip()
        cz = row.get("center_z", "").strip()
        residue_ids = row.get("residue_ids", "").strip()
        n_residues = len(residue_ids.split()) if residue_ids else 0

        center_str = f"({float(cx):.2f}, {float(cy):.2f}, {float(cz):.2f})"
        print(f"{rank:<6}{score:<10}{prob:<14}{center_str:<32}{n_residues:<10}")

        # Machine-parseable line for the GUI to pick up reliably
        print(f"POCKET_ROW,{rank},{score},{prob},{cx},{cy},{cz},{n_residues}")


def main():
    args = parse_args()
    if not os.path.isfile(args.pdb_path):
        print(f"ERROR: file not found: {args.pdb_path}")
        sys.exit(1)

    output_dir = args.output_dir
    if not output_dir:
        stem = os.path.splitext(os.path.basename(args.pdb_path))[0]
        output_dir = os.path.join(os.getcwd(), "docking_runs", f"{stem}_pockets")
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Running P2Rank pocket prediction")
    print("=" * 60)
    run_p2rank(args.pdb_path, output_dir)

    csv_path = find_predictions_csv(output_dir, args.pdb_path)
    print(f"\nReading predictions from {csv_path}")
    print_pockets(csv_path)
    print("\nPOCKETS_COMPLETE")


if __name__ == "__main__":
    main()
