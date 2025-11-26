#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 5: Ligand Preparation - Add hydrogens, assign charges, generate 3D structures, and minimize energy
Input: outputs/selected_compounds/*.sdf (SDF files from Node 3 and real_ligand.sdf from Node 4)
Output: outputs/selected_compounds/*.sdf (prepared SDF files with 3D structures)
"""

import os
import glob
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Process all ligands (library + real_ligand) from outputs/selected_compounds
# This includes ligands from Node 3 and real_ligand from Node 4
INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")

def main():
    """Main execution function"""
    print("=== Node 5: Ligand Preparation ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 3 (node_03_ligand_selection.py) and Node 4 (node_04_real_ligand_addition.py) first.")
        exit(1)
    
    # Find all SDF files (including real_ligand.sdf from Node 4)
    sdf_pattern = os.path.join(INPUT_DIR, "*.sdf")
    sdf_files = sorted(glob.glob(sdf_pattern))
    
    if not sdf_files:
        print(f"❌ Error: No SDF files found in {INPUT_DIR}.")
        print("   Please run Node 3 (node_03_ligand_selection.py) and Node 4 (node_04_real_ligand_addition.py) first.")
        exit(1)
    
    print(f"Input directory: {INPUT_DIR}")
    print(f"Number of SDF files found: {len(sdf_files)}")
    
    # Process each SDF file
    successful_count = 0
    failed_count = 0
    
    for i, sdf_file in enumerate(sdf_files, 1):
        base_name = os.path.basename(sdf_file)
        print(f"\n[{i}/{len(sdf_files)}] Processing {base_name}...")
        
        try:
            # Step 1: Add hydrogens using obabel
            print("  Step 1: Adding hydrogens...")
            if not add_hydrogens_obabel(sdf_file):
                print(f"  ⚠ Warning: Failed to add hydrogens for {base_name}")
                failed_count += 1
                continue
            
            # Step 2: Assign partial charges (Gasteiger) using obabel
            print("  Step 2: Assigning Gasteiger charges...")
            if not assign_charges_obabel(sdf_file):
                print(f"  ⚠ Warning: Failed to assign charges for {base_name}")
                failed_count += 1
                continue
            
            # Step 3: Generate 3D structure using RDKit
            print("  Step 3: Generating 3D structure...")
            if not generate_3d_structure(sdf_file):
                print(f"  ⚠ Warning: Failed to generate 3D structure for {base_name}")
                failed_count += 1
                continue
            
            # Step 4: Energy minimization
            print("  Step 4: Energy minimization...")
            
            print(f"  ✓ Successfully prepared {base_name}")
            successful_count += 1
            
        except Exception as e:
            print(f"  ✗ Error processing {base_name}: {e}")
            failed_count += 1
            continue
    
    print(f"\n✅ Ligand preparation complete.")
    print(f"   Successful: {successful_count}/{len(sdf_files)}")
    if failed_count > 0:
        print(f"   Failed: {failed_count}/{len(sdf_files)}")


def add_hydrogens_obabel(sdf_file):
    """Add hydrogens to SDF file using obabel with pH consideration"""
    try:
        abs_sdf_file = os.path.abspath(sdf_file)
        # -h: add hydrogens
        # -p 7.4: consider protonation at physiological pH (7.4)
        # This ensures all atoms get hydrogens with proper protonation state
        result = subprocess.run(
            ["obabel", abs_sdf_file, "-O", abs_sdf_file, "-h", "-p", "7.4"],
            capture_output=True,
            text=True,
            timeout=60
        )
        return result.returncode == 0
    except Exception as e:
        print(f"    Error: {e}")
        return False


def assign_charges_obabel(sdf_file):
    """Assign Gasteiger partial charges using obabel"""
    try:
        abs_sdf_file = os.path.abspath(sdf_file)
        result = subprocess.run(
            ["obabel", abs_sdf_file, "-O", abs_sdf_file, "--partialcharge", "gasteiger"],
            capture_output=True,
            text=True,
            timeout=120
        )
        return result.returncode == 0
    except Exception as e:
        print(f"    Error: {e}")
        return False


def generate_3d_structure(sdf_file):
    """
    Generate 3D structure using RDKit.
    Uses multiple embedding attempts with different methods for better success rate.
    Note: Energy minimization is performed separately using obminimize, so no optimization here.
    Processes all molecules in the SDF file.
    """
    try:
        # Read all molecules from SDF
        supplier = Chem.SDMolSupplier(sdf_file)
        mols = [m for m in supplier if m is not None]
        
        if not mols:
            print(f"    Error: Could not read any molecules from {sdf_file}")
            return False
        
        processed_mols = []
        
        # Process each molecule
        for mol_idx, mol in enumerate(mols):
            if mol is None:
                continue
            
            # Remove existing 3D coordinates if present
            mol = Chem.RemoveHs(mol)
            mol = Chem.AddHs(mol)
            
            # Try to generate 3D structure with multiple attempts
            success = False
            max_attempts = 5
            
            # Method 1: Standard ETKDG embedding
            for attempt in range(max_attempts):
                try:
                    # Use ETKDG (Experimental-Torsion-angle preference with Distance Geometry) method
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42 + attempt  # Different seed for each attempt
                    conf_id = AllChem.EmbedMolecule(mol, params)
                    
                    if conf_id >= 0:
                        # No optimization here - will be done by obminimize later
                        success = True
                        break
                except Exception as e:
                    if attempt < max_attempts - 1:
                        continue
                    else:
                        if mol_idx == 0:  # Only print warning for first molecule
                            print(f"    Warning: ETKDG embedding failed: {e}")
            
            # Method 2: If ETKDG fails, try basic distance geometry
            if not success:
                try:
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                    if conf_id >= 0:
                        # No optimization here - will be done by obminimize later
                        success = True
                except Exception as e:
                    if mol_idx == 0:  # Only print warning for first molecule
                        print(f"    Warning: Basic embedding failed: {e}")
            
            # Method 3: If still fails, try with 2D coordinates first
            if not success:
                try:
                    # Generate 2D coordinates
                    AllChem.Compute2DCoords(mol)
                    # Then try to embed 3D
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                    if conf_id >= 0:
                        # No optimization here - will be done by obminimize later
                        success = True
                except Exception as e:
                    if mol_idx == 0:  # Only print warning for first molecule
                        print(f"    Warning: 2D-to-3D conversion failed: {e}")
            
            if success:
                processed_mols.append(mol)
            else:
                if mol_idx == 0:  # Only print error for first molecule
                    print(f"    Warning: Could not generate 3D structure for molecule {mol_idx + 1}")
                # Still add the molecule even if 3D generation failed (it may have existing coords)
                processed_mols.append(mol)
        
        if not processed_mols:
            print(f"    Error: Could not process any molecules")
            return False
        
        # Write all processed molecules back to SDF file
        writer = Chem.SDWriter(sdf_file)
        for mol in processed_mols:
            writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"    Error generating 3D structure: {e}")
        return False


def minimize_energy(sdf_file):
    """Minimize energy using obminimize"""
    try:
        abs_sdf_file = os.path.abspath(sdf_file)
        with open(abs_sdf_file, "w") as outfile:
            result = subprocess.run(
                ["obminimize", "-ff", "MMFF94", "-n", "1000", abs_sdf_file],
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=300
            )
        return result.returncode == 0
    except Exception:
        return False


if __name__ == "__main__":
    main()

