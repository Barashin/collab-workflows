#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 6: Real ligand addition - Combine prepared ligands with real ligand from PDB
Input: ligands_3D.sdf (from Node 4), real_ligand.sdf (from Node 8)
Output: real_ligands_3D_sdf (combined ligands)
"""

import os
import glob
from rdkit import Chem

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Silva mounts inputs from depends_on nodes to the current directory
# inputs = ["*.pdb", "real_ligand.sdf"] means:
# - *.pdb files from 1-protein_preparation are mounted at ./outputs/*.pdb
# - real_ligand.sdf from 3-docking_preparation is mounted at ./outputs/real_ligand.sdf
INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
# real_ligand.sdf is mounted from 3-docking_preparation at ./outputs/real_ligand.sdf
REAL_LIGAND_SDF = os.path.join(SCRIPT_DIR, "outputs", "real_ligand.sdf")

def main():
    """Main execution function"""
    print("=== Node 6: Real ligand addition ===")
    
    # Path to real ligand from Node 8 (mounted by Silva from 3-docking_preparation)
    real_ligand_sdf = REAL_LIGAND_SDF
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 4 (node_04_prepare_ligands.py) first.")
        exit(1)
    
    if not os.path.exists(real_ligand_sdf):
        print(f"⚠ Warning: {real_ligand_sdf} not found.")
        print("   Skipping real ligand addition. Using only prepared ligands.")
        return
    
    # Find all prepared ligand SDF files
    sdf_pattern = os.path.join(INPUT_DIR, "*.sdf")
    prepared_ligand_files = sorted(glob.glob(sdf_pattern))
    
    # Exclude real_ligand if it's already in the directory
    prepared_ligand_files = [f for f in prepared_ligand_files 
                            if not os.path.basename(f).startswith("real_ligand")]
    
    if not prepared_ligand_files:
        print(f"⚠ Warning: No prepared ligand files found in {INPUT_DIR}")
        return
    
    print(f"Found {len(prepared_ligand_files)} prepared ligand(s)")
    print(f"Adding real ligand from: {os.path.basename(real_ligand_sdf)}")
    
    # Copy real ligand to selected_compounds directory
    real_ligand_dest = os.path.join(OUTPUT_DIR, "real_ligand.sdf")
    try:
        import shutil
        shutil.copy2(real_ligand_sdf, real_ligand_dest)
        print(f"✅ Real ligand added to: {os.path.basename(real_ligand_dest)}")
    except Exception as e:
        print(f"⚠ Warning: Error copying real ligand: {e}")
        return
    
    print(f"✅ Real ligand addition complete. Total ligands: {len(prepared_ligand_files) + 1}")


if __name__ == "__main__":
    main()

