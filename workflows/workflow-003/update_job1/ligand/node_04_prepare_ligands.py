#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 4: Ligand Preparation - Add hydrogens, assign charges, generate 3D structures, and minimize energy
Input: selected_compounds/*.sdf (SDF files from Node 5)
Output: selected_compounds/*.sdf (prepared SDF files with 3D structures)
"""

import os
import glob
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "output", "selected_compounds")

def main():
    """Main execution function"""
    print("=== Node 4: Ligand Preparation ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 3 (node_03_ligand_selection.py) first.")
        exit(1)
    
    # Find all SDF files
    sdf_pattern = os.path.join(INPUT_DIR, "*.sdf")
    sdf_files = sorted(glob.glob(sdf_pattern))
    
    if not sdf_files:
        print(f"❌ Error: No SDF files found in {INPUT_DIR}.")
        print("   Please run Node 3 (node_03_ligand_selection.py) first.")
        exit(1)
    
    print(f"Input directory: {INPUT_DIR}")
    print(f"Number of SDF files found: {len(sdf_files)}")
    
    # Check if obabel is available
    if not check_obabel():
        print("❌ Error: Open Babel (obabel) is not available.")
        print("   Please install Open Babel: https://openbabel.org/wiki/Category:Installation")
        exit(1)
    
    # Check if obminimize is available
    if not check_obminimize():
        print("⚠ Warning: obminimize is not available. Energy minimization will be skipped.")
        print("   Please install Open Babel tools for energy minimization.")
        use_obminimize = False
    else:
        use_obminimize = True
    
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
            if use_obminimize:
                print("  Step 4: Energy minimization...")
                if not minimize_energy(sdf_file):
                    print(f"  ⚠ Warning: Failed to minimize energy for {base_name}")
                    # Continue even if minimization fails
            
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


def check_obabel():
    """Check if obabel command is available"""
    try:
        result = subprocess.run(
            ["obabel", "-V"], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_obminimize():
    """Check if obminimize command is available"""
    try:
        result = subprocess.run(
            ["obminimize", "-h"], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def add_hydrogens_obabel(sdf_file):
    """Add hydrogens to SDF file using obabel"""
    try:
        # Create temporary file
        temp_file = sdf_file + ".tmp"
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        abs_temp_file = os.path.abspath(temp_file)
        
        # obabel requires explicit output format specification
        cmd = ["obabel", abs_sdf_file, "-h", "-osdf", "-O", abs_temp_file]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # obabel outputs "N molecules converted" to stderr even on success
        if result.returncode == 0:
            if os.path.exists(abs_temp_file) and os.path.getsize(abs_temp_file) > 0:
                # Replace original file
                os.replace(abs_temp_file, abs_sdf_file)
                return True
            else:
                print(f"    Error: Output file was not created or is empty")
                return False
        else:
            # Print error details for debugging
            if result.stderr:
                print(f"    obabel stderr: {result.stderr[:200]}")
            if result.stdout:
                print(f"    obabel stdout: {result.stdout[:200]}")
            if os.path.exists(abs_temp_file):
                os.remove(abs_temp_file)
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout while adding hydrogens")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def assign_charges_obabel(sdf_file):
    """Assign Gasteiger partial charges using obabel"""
    try:
        # Create temporary file
        temp_file = sdf_file + ".tmp"
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        abs_temp_file = os.path.abspath(temp_file)
        
        result = subprocess.run(
            ["obabel", abs_sdf_file, "--partialcharge", "gasteiger", "-osdf", "-O", abs_temp_file],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        # obabel outputs "N molecules converted" to stderr even on success
        if result.returncode == 0:
            if os.path.exists(abs_temp_file) and os.path.getsize(abs_temp_file) > 0:
                # Replace original file
                os.replace(abs_temp_file, abs_sdf_file)
                return True
            else:
                print(f"    Error: Output file was not created or is empty")
                return False
        else:
            if os.path.exists(abs_temp_file):
                os.remove(abs_temp_file)
            if result.stderr:
                print(f"    obabel stderr: {result.stderr[:200]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout while assigning charges")
        return False
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
        # Create temporary file
        temp_file = sdf_file + ".tmp"
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        abs_temp_file = os.path.abspath(temp_file)
        
        with open(abs_temp_file, "w") as outfile:
            result = subprocess.run(
                ["obminimize", "-ff", "MMFF94", "-n", "1000", abs_sdf_file],
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=300
            )
        
        if result.returncode == 0 and os.path.exists(abs_temp_file) and os.path.getsize(abs_temp_file) > 0:
            # Replace original file
            os.replace(abs_temp_file, abs_sdf_file)
            return True
        else:
            if os.path.exists(abs_temp_file):
                os.remove(abs_temp_file)
            if result.stderr:
                print(f"    obminimize stderr: {result.stderr[:200]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout during energy minimization (skipping)")
        return False
    except Exception as e:
        print(f"    Error during energy minimization: {e}")
        return False


if __name__ == "__main__":
    main()

