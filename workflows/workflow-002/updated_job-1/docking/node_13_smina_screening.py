#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test
Node 11: SMINA - In-silico Screening - In-silico screening
Input: config.txt, cleaned PDB file, selected_compounds/ (individual SDF files)
Output: Docking score (docking result files)
"""

import glob
import os
import subprocess

DEFAULT_PDB_ID = "5Y7J"
# Default chain ID (should match Node 9 protein extraction)
DEFAULT_CHAIN_ID = "A"
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PREPARATION_DIR = os.path.join(SCRIPT_DIR, "..", "preparation", "output")
LIGAND_DIR = os.path.join(SCRIPT_DIR, "..", "ligand", "output")
SELECTED_COMPOUNDS_DIR = os.path.join(LIGAND_DIR, "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
CONFIG_FILE = os.path.join(PREPARATION_DIR, "config.txt")

def main():
    """Main execution function"""
    print("=== Node 11: SMINA - In-silico screening ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID and chain ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Find receptor file (cleaned PDB from Node 9)
    # Priority: chain-specific cleaned PDB, then general cleaned PDB
    receptor_file = None
    if chain_id:
        # Try chain-specific file first (e.g., 5Y7J_chain_A_clean.pdb)
        receptor_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_chain_{chain_id}_clean.pdb")
        if not os.path.exists(receptor_file):
            # Try multiple chains format (e.g., 5Y7J_chain_AB_clean.pdb)
            receptor_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_chain_{chain_id}_clean.pdb")
    
    # Fallback to general cleaned PDB
    if not receptor_file or not os.path.exists(receptor_file):
        receptor_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_clean.pdb")
    
    if not os.path.exists(receptor_file):
        print(f"❌ Error: Receptor file not found.")
        print(f"   Expected: {receptor_file}")
        print("   Please run Node 9 (node_09_protein_extraction.py) first.")
        exit(1)
    
    if not os.path.exists(CONFIG_FILE):
        print(f"❌ Error: {CONFIG_FILE} not found.")
        print("   Please run Node 9 (node_09_ligand_center_identification.py) first.")
        exit(1)
    
    if not os.path.exists(SELECTED_COMPOUNDS_DIR):
        print(f"❌ Error: {SELECTED_COMPOUNDS_DIR} directory not found.")
        print("   Please run Node 6 (node_06_real_ligand_addition.py) first.")
        exit(1)
    
    # Find all SDF files in selected_compounds directory
    sdf_pattern = os.path.join(SELECTED_COMPOUNDS_DIR, "*.sdf")
    ligand_files = sorted(glob.glob(sdf_pattern))
    
    if not ligand_files:
        print(f"❌ Error: No SDF files found in {SELECTED_COMPOUNDS_DIR}.")
        print("   Please run Node 6 (node_06_real_ligand_addition.py) first.")
        exit(1)
    
    print(f"Receptor file: {receptor_file}")
    print(f"Config file: {CONFIG_FILE}")
    print(f"Ligand directory: {SELECTED_COMPOUNDS_DIR}")
    print(f"Number of ligands to dock: {len(ligand_files)}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Check smina command
    check_smina()
    smina_path = get_smina_path()
    
    # Process each ligand file
    successful_count = 0
    failed_count = 0
    
    for i, ligand_file in enumerate(ligand_files, 1):
        # Get base filename without extension for output naming
        base_name = os.path.splitext(os.path.basename(ligand_file))[0]
        out_sdf = os.path.join(OUTPUT_DIR, f"{base_name}_docked.sdf")
        out_log = os.path.join(OUTPUT_DIR, f"{base_name}_log.txt")
        
        print(f"\n[{i}/{len(ligand_files)}] Docking {base_name}...")
        
        try:
            result = subprocess.run(
                [
                    smina_path,
                    "-r",
                    receptor_file,
                    "-l",
                    ligand_file,
                    "--config",
                    CONFIG_FILE,
                    "-o",
                    out_sdf,
                    "--log",
                    out_log,
                    "--scoring",
                    "vina",
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            print(f"✓ Docking complete for {base_name}.")
            
            affinity = extract_affinity_from_log(out_log)
            if affinity is not None:
                print(f"  Binding energy: {affinity:.2f} kcal/mol")
            
            successful_count += 1
        
        except subprocess.TimeoutExpired:
            print(f"✗ Docking timeout for {base_name}.")
            failed_count += 1
        except subprocess.CalledProcessError as e:
            print(f"✗ Error occurred during docking for {base_name}: {e}")
            if e.stderr:
                print(f"  Error details: {e.stderr[:200]}...")
            failed_count += 1
        except Exception as e:
            print(f"✗ Unexpected error occurred during docking for {base_name}: {e}")
            failed_count += 1
    
    print(f"\n✅ Docking screening complete.")
    print(f"   Successful: {successful_count}/{len(ligand_files)}")
    if failed_count > 0:
        print(f"   Failed: {failed_count}/{len(ligand_files)}")
    print(f"Results saved in {OUTPUT_DIR}/ directory.")


def get_smina_path():
    """Get smina executable path (same directory as this script)"""
    smina_path = os.path.join(SCRIPT_DIR, "smina.osx.12")
    
    # Check if smina.osx.12 exists in the same directory
    if os.path.exists(smina_path):
        return smina_path
    
    # Fallback to system smina command
    return "smina"


def check_smina():
    """Check smina command"""
    smina_path = get_smina_path()
    
    try:
        result = subprocess.run(
            [smina_path, "--help"], capture_output=True, text=True, timeout=10
        )
        print("smina command is available.")
    except subprocess.TimeoutExpired:
        print("smina command check timed out.")
    except FileNotFoundError:
        print(f"❌ Error: {smina_path} not found.")
        print("Please verify that smina is correctly installed.")
        exit(1)
    except Exception as e:
        print(f"Error occurred while checking smina command: {e}")


def extract_affinity_from_log(log_file):
    """Extract binding energy from log file"""
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):  # Get mode1 line
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            affinity = float(parts[1])  # Affinity value (kcal/mol)
                            return affinity
                        except ValueError:
                            pass
    except Exception as e:
        print(f"  Log file read error: {e}")
    return None


if __name__ == "__main__":
    main()

