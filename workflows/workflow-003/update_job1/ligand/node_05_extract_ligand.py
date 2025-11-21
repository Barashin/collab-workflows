#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 5: Extract ligand from PDB file
Input: preparation/output/{pdb_id}_A_NAD.pdb, PARAM_LIGAND_NAME
Output: ligand/output/{ligand_name}.sdf
"""

import os

from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PREPARATION_DIR = os.path.join(SCRIPT_DIR, "..", "preparation", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

class LigandSelect(Select):
    """Select only residues matching ligand name"""
    
    def __init__(self, ligand_name):
        self.ligand_name = ligand_name
    
    def accept_residue(self, residue):
        return residue.get_resname() == self.ligand_name

def extract_ligand_from_pdb(pdb_path, ligand_name, output_sdf):
    """Extract the ligand from a PDB file and save as SDF (keeping original conformation)"""
    if not os.path.exists(pdb_path):
        print(f"File not found: {pdb_path}")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Save only the ligand residue to a temporary PDB
    temp_pdb = f"{ligand_name}_temp.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb, LigandSelect(ligand_name))

    # Read molecule from the PDB without generating new coordinates
    mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False)
    if mol is None:
        print(f"RDKit could not parse the extracted {ligand_name} PDB.")
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)
        return

    # Save to SDF (preserves PDB coordinates)
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()

    os.remove(temp_pdb)
    print(f"Extracted {ligand_name} ligand and saved as: {output_sdf}")

def main():
    """Main execution function"""
    print("=== Node 5: Extract ligand from PDB file ===")
    
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
    
    pdb_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_A_NAD.pdb")
    output_sdf = os.path.join(OUTPUT_DIR, f"{ligand_name}.sdf")
    
    if not os.path.exists(pdb_file):
        print(f"❌ Error: {pdb_file} not found.")
        print("   Please run Node 2 (node_02_extract_chain.py) first.")
        exit(1)
    
    extract_ligand_from_pdb(pdb_file, ligand_name, output_sdf)

if __name__ == "__main__":
    main()

