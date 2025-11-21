#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 7: Prepare receptor for docking (convert to PDBQT)
Input: preparation/output/{pdb_id}_A_NAD_fixed_with_NAD.pdb
Output: docking/output/{pdb_id}_A_NAD_fixed_with_NAD.pdbqt
"""

import os
import subprocess

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PREPARATION_DIR = os.path.join(SCRIPT_DIR, "..", "preparation", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def prepare_receptor(pdb_id):
    """Prepare receptor PDBQT (preserve NAD)."""
    print("\n--- Preparing receptor ---")

    receptor_pdb = os.path.join(PREPARATION_DIR, f"{pdb_id}_A_NAD_fixed_with_NAD.pdb")
    receptor_pdbqt = os.path.join(OUTPUT_DIR, f"{pdb_id}_A_NAD_fixed_with_NAD.pdbqt")

    if os.path.exists(receptor_pdbqt):
        print(f"✔ Using existing receptor file: {receptor_pdbqt}")
        return receptor_pdbqt

    if not os.path.exists(receptor_pdb):
        print(f"❌ Error: {receptor_pdb} not found.")
        print("   Please run Node 4 (node_04_fix_structure.py) first.")
        exit(1)

    cmd = [
        "prepare_receptor4.py",
        "-r",
        receptor_pdb,
        "-o",
        receptor_pdbqt,
        "-A",
        "hydrogens",
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"Receptor prepared successfully: {receptor_pdbqt}")
    except FileNotFoundError:
        print("✗ prepare_receptor4.py not found. Install MGLTools/AutoDockTools.")
        exit(1)
    except subprocess.CalledProcessError as e:
        print(f"✗ Error preparing receptor:\n{e}")
        exit(1)

    return receptor_pdbqt

def main():
    """Main execution function"""
    print("=== Node 7: Prepare receptor for docking ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID from environment variable
    pdb_id = os.getenv("PARAM_PDB_ID")
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    
    prepare_receptor(pdb_id)

if __name__ == "__main__":
    main()

