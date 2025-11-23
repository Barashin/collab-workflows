#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 6: Protein input - Verify protein input file
Input: pdb.id (environment variable, global_params.json, or default value)
Output: {pdb_id}.pdb (verify already downloaded file)
"""

import json
import os

DEFAULT_PDB_ID = "5Y7J"
# Get script directory and set input directory relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", None)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, None

def main():
    """Main execution function"""
    print("=== Node 6: Verify protein input file ===")
    
    # Get PDB ID (priority: environment variable > global_params.json > default)
    pdb_id = os.environ.get("PARAM_PDB_ID") or os.environ.get("PDB_ID")
    if not pdb_id:
        pdb_id, _ = load_global_params()
        if pdb_id:
            print(f"📋 Loaded PDB ID from global_params.json: {pdb_id}")
    
    if not pdb_id:
        pdb_id = DEFAULT_PDB_ID
        print(f"⚠ Using default PDB ID: {pdb_id}")
    
    pdb_file = os.path.join(INPUT_DIR, f"{pdb_id}.pdb")
    
    print(f"PDB ID: {pdb_id}")
    print(f"File to verify: {pdb_file}")
    
    if os.path.exists(pdb_file):
        file_size = os.path.getsize(pdb_file)
        print(f"✅ File exists: {pdb_file} ({file_size / 1024:.2f} KB)")
    else:
        print(f"❌ Error: {pdb_file} not found.")
        print("   Please run Node 5 (node_05_download_pdb.py) first.")
        exit(1)


if __name__ == "__main__":
    main()

