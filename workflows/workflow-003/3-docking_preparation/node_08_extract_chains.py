#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 8: Chain extraction - Extract chains and real ligand
Input: 5Y7J.pdb
Output: 5Y7J_chain.pdb (PDB file with extracted chains), real_ligand.sdf (extracted ligand)

Environment variables:
    PDB_ID: PDB ID (default: "5Y7J")
    LIGAND_NAME: Ligand name to extract (default: "8OL")
    CHAIN_ID: Chain ID(s) to extract, comma-separated for multiple chains (e.g., "A" or "A,B")
              Default is "A". The ligand will be extracted only from the selected chains.
"""

import os
import json
import tempfile
from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem

DEFAULT_PDB_ID = "5Y7J"
LIGAND_NAME = "8OL"
# Default chain ID (can be overridden by environment variable CHAIN_ID)
# Can be a single chain (e.g., "A") or multiple chains separated by commas (e.g., "A,B")
DEFAULT_CHAIN_ID = "A"  # Default is chain A

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json first, then workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if not os.path.exists(global_params_file):
        global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", LIGAND_NAME)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, LIGAND_NAME
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Silva mounts inputs from depends_on nodes to the current directory
# inputs = ["*.pdb"] means PDB files are mounted at the working directory root
# Try both root directory and outputs directory
INPUT_DIR_ROOT = SCRIPT_DIR  # Root of working directory
INPUT_DIR_OUTPUTS = os.path.join(SCRIPT_DIR, "outputs")  # outputs subdirectory
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

def main():
    """Main execution function"""
    print("=== Node 8: Chain extraction ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or global_ligand_name
    
    # Try to find PDB file in root directory first, then outputs directory
    input_pdb_file = os.path.join(INPUT_DIR_ROOT, f"{pdb_id}.pdb")
    if not os.path.exists(input_pdb_file):
        input_pdb_file = os.path.join(INPUT_DIR_OUTPUTS, f"{pdb_id}.pdb")
    
    output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_chain.pdb")
    output_ligand_sdf = os.path.join(OUTPUT_DIR, "real_ligand.sdf")
    
    print(f"Input file: {input_pdb_file}")
    print(f"Output file: {output_pdb_file}")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 5 (node_05_download_pdb.py) first.")
        exit(1)
    
    # Read PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb_file)
    
    # Check chains
    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(list(set(chains)))
    print(f"\nChains present in PDB file: {', '.join(unique_chains)}")
    
    # Get chain IDs to extract from environment variable
    chain_id_env = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Parse chain IDs (can be comma-separated like "A,B" or single like "A")
    if chain_id_env:
        selected_chains = [c.strip().upper() for c in chain_id_env.split(",")]
        selected_chains = [c for c in selected_chains if c]  # Remove empty strings
    else:
        # If not specified, extract all chains
        selected_chains = unique_chains
        print(f"No CHAIN_ID specified. Will extract all chains: {', '.join(selected_chains)}")
    
    # Validate that all selected chains exist in the PDB file
    missing_chains = [c for c in selected_chains if c not in unique_chains]
    if missing_chains:
        print(f"❌ Error: The following chains are not present in the PDB file: {', '.join(missing_chains)}")
        print(f"   Available chains: {', '.join(unique_chains)}")
        print(f"   Please set CHAIN_ID environment variable to one or more of the available chains.")
        print(f"   Example: export CHAIN_ID=A or export CHAIN_ID=A,B")
        exit(1)
    
    print(f"Selected chains to extract: {', '.join(selected_chains)}")
    
    # Class to select specified chains and ligand
    class ChainAndLigandSelect(Select):
        def __init__(self, chain_ids, ligand_name):
            self.chain_ids = set(chain_ids)
            self.ligand_name = ligand_name
        
        def accept_chain(self, chain):
            return chain.id in self.chain_ids
        
        def accept_residue(self, residue):
            return (
                residue.get_parent().id in self.chain_ids
                or residue.get_resname().strip().upper() == self.ligand_name.upper()
            )
    
    # Extract selected chains
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, ChainAndLigandSelect(selected_chains, ligand_name))
    print(f"✅ Chain extraction complete: {output_pdb_file}")
    
    # Extract real ligand from selected chains and save as SDF
    extract_ligand_to_sdf(input_pdb_file, ligand_name, output_ligand_sdf, selected_chains)


def extract_ligand_to_sdf(pdb_file, ligand_name, output_sdf_file, selected_chains):
    """
    Extract a specific ligand from selected chains in PDB file and save as SDF.
    
    Args:
        pdb_file: Path to input PDB file
        ligand_name: Name of the ligand residue
        output_sdf_file: Path to output SDF file
        selected_chains: List of chain IDs to extract ligand from
    """
    try:
        # Parse PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        
        # Create a Select class to extract only the specified ligand from selected chains
        class LigandSelect(Select):
            def __init__(self, chain_ids, ligand_name):
                self.chain_ids = set(chain_ids)
                self.ligand_name = ligand_name.upper()
            
            def accept_residue(self, residue):
                residue_name = residue.get_resname().strip().upper()
                chain_id = residue.get_parent().id
                # Only extract ligand if it's in one of the selected chains
                return (residue_name == self.ligand_name and 
                        chain_id in self.chain_ids)
        
        # Save ligand to temporary PDB file
        temp_fd, temp_pdb_file = tempfile.mkstemp(suffix=".pdb", prefix="ligand_")
        os.close(temp_fd)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(temp_pdb_file, LigandSelect(selected_chains, ligand_name))
        
        # Check if temporary file has any content (ligand was found)
        if os.path.getsize(temp_pdb_file) == 0:
            print(f"⚠ Warning: Ligand '{ligand_name}' not found in selected chains: {', '.join(selected_chains)}")
            os.remove(temp_pdb_file)
            return
        
        # Read molecule from temporary PDB using RDKit
        mol = Chem.MolFromPDBFile(temp_pdb_file, removeHs=False)
        
        if mol is None:
            print(f"⚠ Warning: Could not parse ligand '{ligand_name}' using RDKit.")
            os.remove(temp_pdb_file)
            return
        
        # Save to SDF
        writer = Chem.SDWriter(output_sdf_file)
        writer.write(mol)
        writer.close()
        
        # Clean up temporary file
        os.remove(temp_pdb_file)
        
        print(f"✅ Extracted ligand '{ligand_name}' from chain(s) {', '.join(selected_chains)} and saved to: {os.path.basename(output_sdf_file)}")
        
    except Exception as e:
        print(f"⚠ Warning: Error extracting ligand '{ligand_name}': {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

