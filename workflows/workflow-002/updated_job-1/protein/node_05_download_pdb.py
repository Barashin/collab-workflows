#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 5: Download File - Download PDB file
Input: pdb.id (environment variable or default value)
Output: 5Y7J.pdb (based on PDB ID)
"""

import os
import urllib.request

# Default PDB ID
DEFAULT_PDB_ID = "5Y7J"
# Get script directory and set output directory relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def main():
    """Main execution function"""
    print("=== Node 5: Download PDB file ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID (from environment variable or default value)
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
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

