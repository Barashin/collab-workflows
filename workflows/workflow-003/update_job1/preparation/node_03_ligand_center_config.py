#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 3: Calculate ligand center and generate docking config
Input: protein/output/{pdb_id}.pdb, PARAM_LIGAND_NAME
Output: preparation/output/config.txt
"""

import os

import numpy as np
from Bio.PDB import PDBParser

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROTEIN_DIR = os.path.join(SCRIPT_DIR, "..", "protein", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def calculate_ligand_center(pdb_file, ligand_name):
    """Compute ligand center coordinates"""
    print(f"\n=== Calculating ligand ({ligand_name}) center coordinates ===")

    extracted_ligand_coords = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target_protein", pdb_file)

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_name:
                    coords = [atom.get_coord() for atom in residue]
                    if coords:
                        extracted_ligand_coords.append(coords)

    if not extracted_ligand_coords:
        print(f"❌ Ligand '{ligand_name}' not found in {pdb_file}")
        exit(1)

    coords_array = np.array(extracted_ligand_coords[0])
    center = np.mean(coords_array, axis=0)

    print(
        f"Center of {ligand_name}: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})"
    )
    return center, coords_array

def generate_docking_config(center_data):
    """Generate docking config file"""
    print("\n=== Generating docking configuration file ===")

    center, coords_array = center_data
    lig_min = coords_array.min(axis=0)
    lig_max = coords_array.max(axis=0)
    extent = lig_max - lig_min
    padding = 8.0
    min_size = 20.0
    size_vec = np.maximum(extent + padding, min_size)

    exhaustiveness = os.getenv("PARAM_EXHAUSTIVENESS", "8")
    num_modes = os.getenv("PARAM_NUM_MODES", "1")
    energy_range = os.getenv("PARAM_ENERGY_RANGE", "3")

    config_path = os.path.join(OUTPUT_DIR, "config.txt")
    config_lines = [
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {size_vec[0]:.1f}",
        f"size_y   = {size_vec[1]:.1f}",
        f"size_z   = {size_vec[2]:.1f}",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {num_modes}",
        f"energy_range = {energy_range}",
    ]

    with open(config_path, "w") as f:
        f.write("\n".join(config_lines) + "\n")

    print(f"✅ Wrote docking config: {config_path}")

def main():
    """Main execution function"""
    print("=== Node 3: Calculate ligand center and generate docking config ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get parameters from environment variables
    pdb_id = os.getenv("PARAM_PDB_ID")
    ligand_name = os.getenv("PARAM_LIGAND_NAME")
    
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    if not ligand_name:
        print("❌ Error: PARAM_LIGAND_NAME environment variable is not set.")
        exit(1)
    
    pdb_file = os.path.join(PROTEIN_DIR, f"{pdb_id}.pdb")
    
    if not os.path.exists(pdb_file):
        print(f"❌ Error: {pdb_file} not found.")
        print("   Please run Node 1 (node_01_download_pdb.py) first.")
        exit(1)
    
    # Calculate ligand center
    center_data = calculate_ligand_center(pdb_file, ligand_name)
    
    # Generate docking config
    generate_docking_config(center_data)

if __name__ == "__main__":
    main()

