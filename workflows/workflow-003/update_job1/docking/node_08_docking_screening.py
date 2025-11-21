#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 8: In-silico screening using SMINA
Input: 
  - docking/output/{pdb_id}_A_NAD_fixed_with_NAD.pdbqt
  - preparation/output/config.txt
  - ligand/output/ligand_library/ (SDF files)
Output: docking/output/ (docking result files)
"""

import os
import subprocess

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PREPARATION_DIR = os.path.join(SCRIPT_DIR, "..", "preparation", "output")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "..", "ligand", "output")
LIGAND_LIBRARY_DIR = os.path.join(LIGAND_DIR, "ligand_library")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
CONFIG_FILE = os.path.join(PREPARATION_DIR, "config.txt")

def check_smina():
    """Ensure smina command is available."""
    print("\n--- Checking smina installation ---")
    try:
        result = subprocess.run(
            ["smina", "--version"], capture_output=True, text=True, timeout=10
        )
        print("✔ smina found:")
        print(result.stdout.strip() or result.stderr.strip())
    except FileNotFoundError:
        print("✗ Error: smina not found in PATH. Please install Smina and try again.")
        exit(1)
    except Exception as e:
        print(f"✗ Unexpected error checking smina: {e}")
        exit(1)

def extract_affinity(log_file):
    """Extract top binding energy from a Smina log file."""
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1])
    except Exception as e:
        print(f"Could not read {log_file}: {e}")
    return None

def run_docking(ligand_name, receptor_pdbqt, ligand_path, output_dir):
    """Dock a single ligand using Smina."""
    print(f"\n--- Docking {ligand_name} ---")

    out_sdf = os.path.join(output_dir, f"{ligand_name}_docked.sdf")
    out_log = os.path.join(output_dir, f"{ligand_name}_docking.log")

    smina_cmd = [
        "smina",
        "-r",
        receptor_pdbqt,
        "-l",
        ligand_path,
        "--config",
        CONFIG_FILE,
        "-o",
        out_sdf,
        "--log",
        out_log,
        "--scoring",
        "vina",
        "--num_modes",
        "1",
    ]

    print("Command:")
    print(" ".join(smina_cmd))

    try:
        subprocess.run(
            smina_cmd, check=True, capture_output=True, text=True, timeout=600
        )
        print(f"✔ Docking complete: {ligand_name}")
    except subprocess.TimeoutExpired:
        print(f"⚠ Docking timed out for {ligand_name}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Docking failed for {ligand_name}: {e}")
    except Exception as e:
        print(f"✗ Unexpected error during docking: {e}")

    affinity = extract_affinity(out_log)
    if affinity is not None:
        print(f"Predicted binding energy: {affinity:.2f} kcal/mol")
    else:
        print("⚠ Could not extract binding energy.")

def main():
    """Main execution function"""
    print("=== Node 8: In-silico screening using SMINA ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Check smina
    check_smina()
    
    # Get PDB ID from environment variable
    pdb_id = os.getenv("PARAM_PDB_ID")
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    
    # Prepare receptor
    receptor_pdbqt = os.path.join(OUTPUT_DIR, f"{pdb_id}_A_NAD_fixed_with_NAD.pdbqt")
    if not os.path.exists(receptor_pdbqt):
        print(f"❌ Error: {receptor_pdbqt} not found.")
        print("   Please run Node 7 (node_07_prepare_receptor.py) first.")
        exit(1)
    
    # Check config file
    if not os.path.exists(CONFIG_FILE):
        print(f"❌ Error: {CONFIG_FILE} not found.")
        print("   Please run Node 3 (node_03_ligand_center_config.py) first.")
        exit(1)
    
    # Check ligand library
    if not os.path.isdir(LIGAND_LIBRARY_DIR):
        print(f"✗ Error: ligand folder '{LIGAND_LIBRARY_DIR}' not found.")
        print("   Please run Node 6 (node_06_generate_variants.py) first.")
        exit(1)

    # Dock all ligands in ligand_library/
    ligands = [f for f in os.listdir(LIGAND_LIBRARY_DIR) if f.endswith(".sdf")]
    if not ligands:
        print(f"✗ No .sdf files found in {LIGAND_LIBRARY_DIR}/")
        exit(1)

    print(f"\nFound {len(ligands)} ligand(s) to dock:\n" + "\n".join(ligands))

    for ligand_file in ligands:
        ligand_path = os.path.join(LIGAND_LIBRARY_DIR, ligand_file)
        ligand_name = os.path.splitext(ligand_file)[0]
        run_docking(ligand_name, receptor_pdbqt, ligand_path, OUTPUT_DIR)

    print("\n=== Docking complete ===")
    print(f"Results saved in: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()

