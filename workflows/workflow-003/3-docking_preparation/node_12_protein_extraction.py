#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 9: OpenMM PDBFixer - Clean protein - Clean up protein structure
Input: 5Y7J.pdb (full protein structure)
Output: 5Y7J_clean.pdb (cleaned up PDB file with selected chains)
"""

import os
import json
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from Bio.PDB import PDBIO, PDBParser, Select

DEFAULT_PDB_ID = "5Y7J"
# Default chain ID (can be overridden by environment variable, None means all chains)
DEFAULT_CHAIN_ID =  "A"

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

def main():
    """Main execution function"""
    print("=== Node 9: OpenMM PDBFixer - Clean up protein ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    selected_chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Try to use chain-extracted PDB from node_08 first (in outputs/), fallback to original PDB (in input/)
    # node_08 outputs to outputs/ in the same node, so check outputs/ first
    chain_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_chain.pdb")
    original_pdb_file = os.path.join(INPUT_DIR, f"{pdb_id}.pdb")
    
    input_pdb_file = None
    if os.path.exists(chain_pdb_file):
        input_pdb_file = chain_pdb_file
        print(f"Using chain-extracted PDB from node_08: {input_pdb_file}")
    elif os.path.exists(original_pdb_file):
        input_pdb_file = original_pdb_file
        print(f"Using original PDB file from input: {input_pdb_file}")
    
    if not input_pdb_file:
        print(f"❌ Error: PDB file not found.")
        print(f"   Searched for chain-extracted: {chain_pdb_file}")
        print(f"   Searched for original: {original_pdb_file}")
        if os.path.exists(OUTPUT_DIR):
            output_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.pdb')]
            if output_files:
                print(f"   Available PDB files in outputs: {output_files}")
        if os.path.exists(INPUT_DIR):
            input_files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
            if input_files:
                print(f"   Available PDB files in input: {input_files}")
        print("   Please ensure 1-protein_preparation has completed and Node 8 (node_08_extract_chains.py) has run.")
        exit(1)
    
    print(f"Input file: {input_pdb_file}")
    
    try:
        # First, check available chains and select if needed
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_pdb_file)
        
        # Get all chains in the structure
        all_chains = []
        for model in structure:
            for chain in model:
                if chain.id not in all_chains:
                    all_chains.append(chain.id)
        
        all_chains = sorted(all_chains)
        print(f"\nAvailable chains in PDB file: {', '.join(all_chains)}")
        
        # Determine which chains to keep
        if selected_chain_id is None:
            chains_to_keep = all_chains
            print(f"Using all chains: {', '.join(chains_to_keep)}")
        else:
            if selected_chain_id not in all_chains:
                print(f"❌ Error: Chain '{selected_chain_id}' not found in PDB file.")
                print(f"   Available chains: {', '.join(all_chains)}")
                exit(1)
            chains_to_keep = [selected_chain_id]
            print(f"Selected chain: {selected_chain_id}")
        
        # Generate output file name based on selected chains
        if len(chains_to_keep) == len(all_chains):
            # All chains selected
            chain_suffix = "clean"
        elif len(chains_to_keep) == 1:
            # Single chain selected
            chain_suffix = f"chain_{chains_to_keep[0]}_clean"
        else:
            # Multiple chains selected
            chain_suffix = f"chain_{''.join(chains_to_keep)}_clean"
        
        output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_{chain_suffix}.pdb")
        
        print(f"Output file: {output_pdb_file}")
        
        # Extract selected chains if needed
        # If using chain-extracted PDB from node_08, skip chain extraction
        temp_pdb_file = input_pdb_file
        is_chain_extracted = input_pdb_file.endswith("_chain.pdb")
        
        if selected_chain_id is not None and not is_chain_extracted:
            # Create a temporary file with only selected chains
            import tempfile
            temp_fd, temp_pdb_file = tempfile.mkstemp(suffix=".pdb", prefix="selected_chains_")
            os.close(temp_fd)
            
            class ChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in chains_to_keep
            
            io = PDBIO()
            io.set_structure(structure)
            io.save(temp_pdb_file, ChainSelect())
            print(f"Extracted chains {', '.join(chains_to_keep)} to temporary file")
        elif is_chain_extracted:
            print(f"Using already chain-extracted PDB from node_08, skipping chain extraction")
        
        # Process with PDBFixer
        print(f"\nProcessing {temp_pdb_file} with PDBFixer...")
        fixer = PDBFixer(filename=temp_pdb_file)
        
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)  # OK if receptor only
        fixer.addMissingHydrogens()
        
        # Write fixed structure
        with open(output_pdb_file, "w") as fout:
            PDBFile.writeFile(fixer.topology, fixer.positions, fout)
        
        # Clean up temporary file if created
        if temp_pdb_file != input_pdb_file and os.path.exists(temp_pdb_file):
            os.remove(temp_pdb_file)
        
        print(f"\n✅ Cleanup complete: {output_pdb_file}")
        if selected_chain_id:
            print(f"   Contains chain(s): {', '.join(chains_to_keep)}")
    except Exception as e:
        print(f"❌ Error occurred during PDBFixer processing: {e}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()

