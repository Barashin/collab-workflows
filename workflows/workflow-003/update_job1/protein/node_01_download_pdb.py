#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 1: Download PDB file
Input: PARAM_PDB_ID (environment variable)
Output: protein/output/{pdb_id}.pdb
"""

import os
import urllib.request

# Get script directory and set output directory relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def main():
    """Main execution function"""
    print("=== Node 1: Download PDB file ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID from environment variable
    pdb_id = os.getenv("PARAM_PDB_ID")
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    
    pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    print(f"PDB ID: {pdb_id}")
    print(f"Download URL: {url}")
    print(f"Output file: {pdb_file}")
    
    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"✅ Download complete: {pdb_file}")
    except Exception as e:
        print(f"❌ Error occurred during download: {e}")
        exit(1)


if __name__ == "__main__":
    main()

