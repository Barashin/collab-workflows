#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 4: Real ligand addition - Download real ligand from PDB and add to selected compounds
Input: global_params.json (pdb_id, ligand_name), outputs/selected_compounds/ (from Node 3)
Output: outputs/selected_compounds/real_ligand.sdf (real ligand downloaded from PDB and added to library)
Note: This runs before Node 5 (prepare_ligands) so that real_ligand is also prepared
"""

import os
import json
import glob
import urllib.request
from rdkit import Chem

# Default values
DEFAULT_PDB_ID = "5Y7J"
DEFAULT_LIGAND_NAME = "8OL"

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
SELECTED_COMPOUNDS_DIR = os.path.join(OUTPUT_DIR, "selected_compounds")
# real_ligand.sdf will be downloaded from PDB
REAL_LIGAND_OUTPUT = os.path.join(SELECTED_COMPOUNDS_DIR, "real_ligand.sdf")

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", DEFAULT_LIGAND_NAME)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, DEFAULT_LIGAND_NAME

def download_ligand_from_pdb(ligand_name, output_file):
    """
    Download ligand SDF file from RCSB PDB
    
    Args:
        ligand_name: Three-letter ligand name (e.g., "RMZ")
        output_file: Path to save the downloaded SDF file
    
    Returns:
        bool: True if download successful, False otherwise
    """
    # RCSB PDB ligand download URL
    # Try ideal SDF first (most common format)
    url = f"https://files.rcsb.org/ligands/download/{ligand_name.upper()}_ideal.sdf"
    
    print(f"Downloading ligand '{ligand_name}' from RCSB PDB...")
    print(f"URL: {url}")
    
    try:
        # Download the file
        urllib.request.urlretrieve(url, output_file)
        
        # Verify the file was downloaded and is not empty
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            # Try to read with RDKit to verify it's a valid SDF
            try:
                suppl = Chem.SDMolSupplier(output_file)
                mol = next(suppl, None)
                if mol is not None:
                    print(f"✅ Successfully downloaded ligand '{ligand_name}' from PDB")
                    print(f"   File saved to: {output_file}")
                    print(f"   File size: {os.path.getsize(output_file)} bytes")
                    return True
                else:
                    print(f"⚠ Warning: Downloaded file exists but could not be parsed as SDF")
                    # Try alternative URL format
                    return download_ligand_alternative(ligand_name, output_file)
            except Exception as e:
                print(f"⚠ Warning: Error reading SDF file: {e}")
                # Try alternative URL format
                return download_ligand_alternative(ligand_name, output_file)
        else:
            print(f"⚠ Warning: Downloaded file is empty, trying alternative URL...")
            return download_ligand_alternative(ligand_name, output_file)
            
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"⚠ Ligand '{ligand_name}' not found at ideal SDF URL, trying alternative format...")
            return download_ligand_alternative(ligand_name, output_file)
        else:
            print(f"❌ HTTP Error {e.code}: {e.reason}")
            return False
    except Exception as e:
        print(f"❌ Error downloading ligand: {e}")
        return False

def download_ligand_alternative(ligand_name, output_file):
    """
    Try alternative URL format for downloading ligand from PDB
    
    Args:
        ligand_name: Three-letter ligand name
        output_file: Path to save the downloaded SDF file
    
    Returns:
        bool: True if download successful, False otherwise
    """
    # Try alternative URL without "_ideal" suffix
    alt_url = f"https://files.rcsb.org/ligands/download/{ligand_name.upper()}.sdf"
    
    print(f"Trying alternative URL: {alt_url}")
    
    try:
        urllib.request.urlretrieve(alt_url, output_file)
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            try:
                suppl = Chem.SDMolSupplier(output_file)
                mol = next(suppl, None)
                if mol is not None:
                    print(f"✅ Successfully downloaded ligand '{ligand_name}' from PDB (alternative URL)")
                    print(f"   File saved to: {output_file}")
                    print(f"   File size: {os.path.getsize(output_file)} bytes")
                    return True
            except Exception:
                pass
        
        print(f"❌ Error: Could not download ligand '{ligand_name}' from PDB")
        print(f"   Tried URLs:")
        print(f"   - https://files.rcsb.org/ligands/download/{ligand_name.upper()}_ideal.sdf")
        print(f"   - https://files.rcsb.org/ligands/download/{ligand_name.upper()}.sdf")
        return False
        
    except Exception as e:
        print(f"❌ Error downloading ligand from alternative URL: {e}")
        return False

def main():
    """Main execution function"""
    print("=== Node 4: Real ligand addition ===")
    
    # Ensure selected_compounds directory exists
    os.makedirs(SELECTED_COMPOUNDS_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or global_ligand_name
    
    print(f"PDB ID: {pdb_id}")
    print(f"Ligand name: {ligand_name}")
    
    # Find all prepared ligand SDF files in selected_compounds
    sdf_pattern = os.path.join(SELECTED_COMPOUNDS_DIR, "*.sdf")
    prepared_ligand_files = sorted(glob.glob(sdf_pattern))
    
    # Exclude real_ligand if it's already in the directory
    prepared_ligand_files = [f for f in prepared_ligand_files 
                            if not os.path.basename(f).startswith("real_ligand")]
    
    print(f"Found {len(prepared_ligand_files)} prepared ligand(s) in {SELECTED_COMPOUNDS_DIR}")
    
    # Check if real_ligand.sdf already exists in destination
    if os.path.exists(REAL_LIGAND_OUTPUT):
        print(f"⚠ real_ligand.sdf already exists in {SELECTED_COMPOUNDS_DIR}, overwriting...")
    
    # Download real ligand from PDB
    print(f"\nDownloading real ligand '{ligand_name}' from RCSB PDB...")
    success = download_ligand_from_pdb(ligand_name, REAL_LIGAND_OUTPUT)
    
    if not success:
        print(f"❌ Error: Failed to download ligand '{ligand_name}' from PDB")
        print(f"   Please check that the ligand name '{ligand_name}' is correct in global_params.json")
        exit(1)
    
    # Verify the file was downloaded successfully
    if os.path.exists(REAL_LIGAND_OUTPUT):
        # Count total ligands after addition
        all_ligand_files = sorted(glob.glob(sdf_pattern))
        total_ligands = len(all_ligand_files)
        print(f"\n✅ Real ligand addition complete. Total ligands in selected_compounds: {total_ligands}")
        
        # Verify real_ligand.sdf is in the list
        ligand_basenames = [os.path.basename(f) for f in all_ligand_files]
        if "real_ligand.sdf" in ligand_basenames:
            print(f"✓ Verified: real_ligand.sdf is in selected_compounds directory")
        else:
            print(f"⚠ Warning: real_ligand.sdf was downloaded but not found in the file list")
    else:
        print(f"❌ Error: real_ligand.sdf was not downloaded successfully")
        exit(1)


if __name__ == "__main__":
    main()

