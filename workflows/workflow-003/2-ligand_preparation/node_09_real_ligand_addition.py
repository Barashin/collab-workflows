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
# - *.pdb files from 1-protein_preparation may be mounted at root or ./outputs/*.pdb
# - real_ligand.sdf from 3-docking_preparation may be mounted at root or ./outputs/real_ligand.sdf
INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
# real_ligand.sdf may be mounted at multiple locations
REAL_LIGAND_ROOT = os.path.join(SCRIPT_DIR, "real_ligand.sdf")
REAL_LIGAND_OUTPUTS = os.path.join(SCRIPT_DIR, "outputs", "real_ligand.sdf")

def main():
    """Main execution function"""
    print("=== Node 6: Real ligand addition ===")
    
    # Ensure selected_compounds directory exists
    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find real_ligand.sdf from 3-docking_preparation (mounted by Silva)
    # Check both root and outputs directories
    real_ligand_source = None
    if os.path.exists(REAL_LIGAND_ROOT):
        real_ligand_source = REAL_LIGAND_ROOT
        print(f"✓ Found real_ligand.sdf at root: {REAL_LIGAND_ROOT}")
    elif os.path.exists(REAL_LIGAND_OUTPUTS):
        real_ligand_source = REAL_LIGAND_OUTPUTS
        print(f"✓ Found real_ligand.sdf at outputs: {REAL_LIGAND_OUTPUTS}")
    
    if not real_ligand_source:
        print(f"❌ Error: real_ligand.sdf not found.")
        print(f"   Searched in: {REAL_LIGAND_ROOT}")
        print(f"   Searched in: {REAL_LIGAND_OUTPUTS}")
        print("   Please ensure 3-docking_preparation has completed and real_ligand.sdf is available.")
        exit(1)
    
    # Find all prepared ligand SDF files
    sdf_pattern = os.path.join(INPUT_DIR, "*.sdf")
    prepared_ligand_files = sorted(glob.glob(sdf_pattern))
    
    # Exclude real_ligand if it's already in the directory
    prepared_ligand_files = [f for f in prepared_ligand_files 
                            if not os.path.basename(f).startswith("real_ligand")]
    
    print(f"Found {len(prepared_ligand_files)} prepared ligand(s) in {INPUT_DIR}")
    print(f"Adding real ligand from: {os.path.basename(real_ligand_source)}")
    
    # Copy real ligand to selected_compounds directory
    real_ligand_dest = os.path.join(OUTPUT_DIR, "real_ligand.sdf")
    
    # Check if real_ligand.sdf already exists in destination
    if os.path.exists(real_ligand_dest):
        print(f"⚠ real_ligand.sdf already exists in {OUTPUT_DIR}, overwriting...")
    
    try:
        import shutil
        shutil.copy2(real_ligand_source, real_ligand_dest)
        print(f"✅ Real ligand copied to: {real_ligand_dest}")
        
        # Verify the file was copied successfully
        if os.path.exists(real_ligand_dest):
            file_size = os.path.getsize(real_ligand_dest)
            print(f"   File size: {file_size} bytes")
            
            # Count total ligands after addition
            all_ligand_files = sorted(glob.glob(sdf_pattern))
            total_ligands = len(all_ligand_files)
            print(f"✅ Real ligand addition complete. Total ligands in selected_compounds: {total_ligands}")
            
            # Verify real_ligand.sdf is in the list
            ligand_basenames = [os.path.basename(f) for f in all_ligand_files]
            if "real_ligand.sdf" in ligand_basenames:
                print(f"✓ Verified: real_ligand.sdf is in selected_compounds directory")
            else:
                print(f"⚠ Warning: real_ligand.sdf was copied but not found in the file list")
        else:
            print(f"❌ Error: real_ligand.sdf was not copied successfully")
            exit(1)
    except Exception as e:
        print(f"❌ Error copying real ligand: {e}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()

