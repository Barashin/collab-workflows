#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 2: Extract Chain A with NAD cofactor
Input: protein/output/{pdb_id}.pdb
Output: preparation/output/{pdb_id}_A_NAD.pdb
"""

import os

from Bio.PDB import PDBIO, PDBParser, Select

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROTEIN_DIR = os.path.join(SCRIPT_DIR, "..", "protein", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def extract_chain_a_with_nad(pdb_file, output_pdb_file):
    """Extract Chain A with NAD cofactor (exclude ligand)"""
    print("\n=== Extracting Chain A with NAD cofactor ===")

    cofactor_name = "NAD"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(set(chains))
    print(f"Chains present: {', '.join(unique_chains)}")

    if "A" not in unique_chains:
        print("⚠ Chain A not found. Using original file.")
        import shutil
        shutil.copy(pdb_file, output_pdb_file)
        return output_pdb_file

    class AChainAndNADSelect(Select):
        def accept_chain(self, chain):
            return chain.id == "A"

        def accept_residue(self, residue):
            if residue.get_parent().id == "A":
                return True
            return residue.get_resname() == cofactor_name

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, AChainAndNADSelect())
    print(f"✅ Extracted: {output_pdb_file}")

    return output_pdb_file

def main():
    """Main execution function"""
    print("=== Node 2: Extract Chain A with NAD cofactor ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID from environment variable
    pdb_id = os.getenv("PARAM_PDB_ID")
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    
    pdb_file = os.path.join(PROTEIN_DIR, f"{pdb_id}.pdb")
    output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_A_NAD.pdb")
    
    if not os.path.exists(pdb_file):
        print(f"❌ Error: {pdb_file} not found.")
        print("   Please run Node 1 (node_01_download_pdb.py) first.")
        exit(1)
    
    extract_chain_a_with_nad(pdb_file, output_pdb_file)

if __name__ == "__main__":
    main()

