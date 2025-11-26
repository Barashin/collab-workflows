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
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input", "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
# Fallback to outputs for backward compatibility
if not os.path.exists(INPUT_DIR):
    INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")

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
    obabel_cmd = find_obabel()
    if not obabel_cmd:
        print("❌ Error: Open Babel (obabel) is not available.")
        print("   Please install Open Babel: https://openbabel.org/wiki/Category:Installation")
        exit(1)
    else:
        print(f"✓ Found obabel at: {obabel_cmd}")
    
    # Check if obminimize is available
    # Try to find obminimize in the same directory as obabel first
    obminimize_cmd = None
    if obabel_cmd:
        obabel_dir = os.path.dirname(obabel_cmd)
        obminimize_in_same_dir = os.path.join(obabel_dir, "obminimize")
        if os.path.exists(obminimize_in_same_dir):
            try:
                result = subprocess.run(
                    [obminimize_in_same_dir, "-h"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    obminimize_cmd = obminimize_in_same_dir
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass
    
    # If not found in same directory, try general search
    if not obminimize_cmd:
        obminimize_cmd = find_obminimize()
    
    if not obminimize_cmd:
        print("⚠ Warning: obminimize is not available. Energy minimization will be skipped.")
        print(f"   Searched in: {obabel_dir if obabel_cmd else 'N/A'}")
        print("   Please install Open Babel tools for energy minimization.")
        use_obminimize = False
    else:
        print(f"✓ Found obminimize at: {obminimize_cmd}")
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


def find_obabel():
    """Find obabel command in common locations"""
    # First, try to find in PATH using shutil.which
    obabel_path = shutil.which("obabel")
    if obabel_path:
        try:
            result = subprocess.run(
                [obabel_path, "-V"], 
                capture_output=True, 
                text=True, 
                timeout=5
            )
            if result.returncode == 0:
                return obabel_path
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    
    # Try to find conda environment paths dynamically
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, "bin", "obabel")
        if os.path.exists(conda_bin):
            try:
                result = subprocess.run(
                    [conda_bin, "-V"], 
                    capture_output=True, 
                    text=True, 
                    timeout=5
                )
                if result.returncode == 0:
                    return conda_bin
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass
    
    # Check common conda installation paths
    possible_paths = [
        "/opt/conda/envs/in_silico/bin/obabel",
        "/opt/conda/bin/obabel",
        "/usr/local/conda/bin/obabel",
        "/usr/local/bin/obabel",
        "/usr/bin/obabel",
    ]
    for path in possible_paths:
        if os.path.exists(path):
            try:
                result = subprocess.run(
                    [path, "-V"], 
                    capture_output=True, 
                    text=True, 
                    timeout=5
                )
                if result.returncode == 0:
                    return path
            except (FileNotFoundError, subprocess.TimeoutExpired):
                continue
    return None


def find_obminimize():
    """Find obminimize command in common locations"""
    # First, try to find in PATH using shutil.which
    obminimize_path = shutil.which("obminimize")
    if obminimize_path:
        try:
            result = subprocess.run(
                [obminimize_path, "-h"], 
                capture_output=True, 
                text=True, 
                timeout=5
            )
            if result.returncode == 0:
                return obminimize_path
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    
    # Try to find conda environment paths dynamically
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, "bin", "obminimize")
        if os.path.exists(conda_bin):
            try:
                result = subprocess.run(
                    [conda_bin, "-h"], 
                    capture_output=True, 
                    text=True, 
                    timeout=5
                )
                if result.returncode == 0:
                    return conda_bin
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass
    
    # Check common conda installation paths
    possible_paths = [
        "/opt/conda/envs/in_silico/bin/obminimize",
        "/opt/conda/bin/obminimize",
        "/usr/local/conda/bin/obminimize",
        "/usr/local/bin/obminimize",
        "/usr/bin/obminimize",
    ]
    for path in possible_paths:
        if os.path.exists(path):
            try:
                result = subprocess.run(
                    [path, "-h"], 
                    capture_output=True, 
                    text=True, 
                    timeout=5
                )
                if result.returncode == 0:
                    return path
            except (FileNotFoundError, subprocess.TimeoutExpired):
                continue
    return None


def check_obabel():
    """Check if obabel command is available"""
    return find_obabel() is not None


def check_obminimize():
    """Check if obminimize command is available"""
    return find_obminimize() is not None


def add_hydrogens_obabel(sdf_file):
    """Add hydrogens to SDF file using obabel"""
    try:
        # Find obabel command
        obabel_cmd = find_obabel()
        if not obabel_cmd:
            print(f"    Error: obabel command not found")
            return False
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        
        # Check if input file exists and is readable
        if not os.path.exists(abs_sdf_file):
            print(f"    Error: Input file does not exist: {abs_sdf_file}")
            return False
        
        if os.path.getsize(abs_sdf_file) == 0:
            print(f"    Error: Input file is empty: {abs_sdf_file}")
            return False
        
        # Create temporary file for output (obabel cannot directly overwrite input file)
        temp_file = abs_sdf_file + ".tmp"
        
        # Try different command formats - some versions of obabel have different syntax
        # Format 1: obabel input.sdf -h -osdf -O output.sdf
        cmd = [obabel_cmd, abs_sdf_file, "-h", "-osdf", "-O", temp_file]
        
        # Debug: print command being executed
        print(f"    Running: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # If format 1 fails with "cannot write output format", try format 2
        if result.returncode != 0 and "cannot write output format" in result.stderr.lower():
            print(f"    Trying alternative command format...")
            # Format 2: obabel input.sdf -O output.sdf -h (let obabel detect format from extension)
            cmd2 = [obabel_cmd, abs_sdf_file, "-O", temp_file, "-h"]
            print(f"    Running: {' '.join(cmd2)}")
            result = subprocess.run(
                cmd2,
                capture_output=True,
                text=True,
                timeout=60
            )
            if result.returncode == 0:
                result = result  # Use the successful result
            else:
                # Format 3: obabel -isdf input.sdf -osdf output.sdf -h
                print(f"    Trying format 3...")
                cmd3 = [obabel_cmd, "-isdf", abs_sdf_file, "-osdf", temp_file, "-h"]
                print(f"    Running: {' '.join(cmd3)}")
                result = subprocess.run(
                    cmd3,
                    capture_output=True,
                    text=True,
                    timeout=60
                )
        
        # Check if output file was created successfully
        if result.returncode == 0:
            if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
                # Replace original file with temporary file
                os.replace(temp_file, abs_sdf_file)
                if result.stderr:
                    # Print success message from stderr (e.g., "1 molecule converted")
                    print(f"    {result.stderr.strip()}")
                return True
            else:
                print(f"    Error: Output file was not created or is empty")
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                if result.stderr:
                    print(f"    obabel stderr: {result.stderr[:300]}")
                if result.stdout:
                    print(f"    obabel stdout: {result.stdout[:300]}")
                return False
        else:
            # Clean up temp file if it exists
            if os.path.exists(temp_file):
                os.remove(temp_file)
            # Print error details for debugging
            print(f"    Return code: {result.returncode}")
            if result.stderr:
                print(f"    Full stderr: {result.stderr}")
                stderr_lower = result.stderr.lower()
                if "error" in stderr_lower or "cannot" in stderr_lower or "failed" in stderr_lower:
                    print(f"    obabel error: {result.stderr[:500]}")
                else:
                    print(f"    obabel stderr: {result.stderr[:500]}")
            if result.stdout:
                print(f"    obabel stdout: {result.stdout[:500]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout while adding hydrogens")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def assign_charges_obabel(sdf_file):
    """Assign Gasteiger partial charges using obabel"""
    try:
        # Find obabel command
        obabel_cmd = find_obabel()
        if not obabel_cmd:
            print(f"    Error: obabel command not found")
            return False
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        
        # Check if input file exists and is readable
        if not os.path.exists(abs_sdf_file):
            print(f"    Error: Input file does not exist: {abs_sdf_file}")
            return False
        
        if os.path.getsize(abs_sdf_file) == 0:
            print(f"    Error: Input file is empty: {abs_sdf_file}")
            return False
        
        # Create temporary file for output (obabel cannot directly overwrite input file)
        temp_file = abs_sdf_file + ".tmp"
        
        # Try different command formats
        # Format 1: obabel input.sdf --partialcharge gasteiger -osdf -O output.sdf
        cmd = [obabel_cmd, abs_sdf_file, "--partialcharge", "gasteiger", "-osdf", "-O", temp_file]
        
        # Debug: print command being executed
        print(f"    Running: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        # If format 1 fails with "cannot write output format", try format 2
        if result.returncode != 0 and "cannot write output format" in result.stderr.lower():
            print(f"    Trying alternative command format...")
            # Format 2: obabel input.sdf -O output.sdf --partialcharge gasteiger
            cmd2 = [obabel_cmd, abs_sdf_file, "-O", temp_file, "--partialcharge", "gasteiger"]
            print(f"    Running: {' '.join(cmd2)}")
            result = subprocess.run(
                cmd2,
                capture_output=True,
                text=True,
                timeout=120
            )
            if result.returncode != 0:
                # Format 3: obabel -isdf input.sdf -osdf output.sdf --partialcharge gasteiger
                print(f"    Trying format 3...")
                cmd3 = [obabel_cmd, "-isdf", abs_sdf_file, "-osdf", temp_file, "--partialcharge", "gasteiger"]
                print(f"    Running: {' '.join(cmd3)}")
                result = subprocess.run(
                    cmd3,
                    capture_output=True,
                    text=True,
                    timeout=120
                )
        
        # Check if output file was created successfully
        if result.returncode == 0:
            if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
                # Replace original file with temporary file
                os.replace(temp_file, abs_sdf_file)
                if result.stderr:
                    # Print success message from stderr (e.g., "1 molecule converted")
                    print(f"    {result.stderr.strip()}")
                return True
            else:
                print(f"    Error: Output file was not created or is empty")
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                if result.stderr:
                    print(f"    obabel stderr: {result.stderr[:300]}")
                if result.stdout:
                    print(f"    obabel stdout: {result.stdout[:300]}")
                return False
        else:
            # Clean up temp file if it exists
            if os.path.exists(temp_file):
                os.remove(temp_file)
            # Print error details for debugging
            print(f"    Return code: {result.returncode}")
            if result.stderr:
                print(f"    Full stderr: {result.stderr}")
                stderr_lower = result.stderr.lower()
                if "error" in stderr_lower or "cannot" in stderr_lower or "failed" in stderr_lower:
                    print(f"    obabel error: {result.stderr[:500]}")
                else:
                    print(f"    obabel stderr: {result.stderr[:500]}")
            if result.stdout:
                print(f"    obabel stdout: {result.stdout[:500]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout while assigning charges")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        import traceback
        traceback.print_exc()
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
        # Find obminimize command
        obminimize_cmd = find_obminimize()
        if not obminimize_cmd:
            print(f"    Error: obminimize command not found")
            return False
        
        # Use absolute path to avoid path issues
        abs_sdf_file = os.path.abspath(sdf_file)
        
        # Check if input file exists and is readable
        if not os.path.exists(abs_sdf_file):
            print(f"    Error: Input file does not exist: {abs_sdf_file}")
            return False
        
        if os.path.getsize(abs_sdf_file) == 0:
            print(f"    Error: Input file is empty: {abs_sdf_file}")
            return False
        
        # Create temporary file for output
        # obminimize outputs to stdout, so we redirect to temp file
        temp_file = sdf_file + ".tmp"
        abs_temp_file = os.path.abspath(temp_file)
        
        # obminimize command: obminimize -ff MMFF94 -n 1000 $i
        # Output goes to stdout, so we redirect to temp file
        with open(abs_temp_file, "w") as outfile:
            result = subprocess.run(
                [obminimize_cmd, "-ff", "MMFF94", "-n", "1000", abs_sdf_file],
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=300
            )
        
        # Check if output file was created and is not empty
        if result.returncode == 0:
            if os.path.exists(abs_temp_file) and os.path.getsize(abs_temp_file) > 0:
                # Replace original file
                os.replace(abs_temp_file, abs_sdf_file)
                return True
            else:
                print(f"    Error: Output file was not created or is empty")
                if os.path.exists(abs_temp_file):
                    os.remove(abs_temp_file)
                if result.stderr:
                    print(f"    obminimize stderr: {result.stderr[:300]}")
                return False
        else:
            if os.path.exists(abs_temp_file):
                os.remove(abs_temp_file)
            if result.stderr:
                print(f"    obminimize stderr: {result.stderr[:300]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"    Timeout during energy minimization (skipping)")
        return False
    except Exception as e:
        print(f"    Error during energy minimization: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    main()

