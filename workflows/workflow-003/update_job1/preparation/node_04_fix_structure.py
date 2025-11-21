#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 4: Fix structure with PDBFixer and reattach NAD
Input: preparation/output/{pdb_id}_A_NAD.pdb
Output: preparation/output/{pdb_id}_A_NAD_fixed_with_NAD.pdb
"""

import os

from openmm.app import PDBFile
from pdbfixer import PDBFixer

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def fix_pdb_structure(input_pdb_file, output_pdb_file):
    """Fix structure using PDBFixer, keeping NAD"""
    print(f"\n=== Fixing {input_pdb_file} with PDBFixer ===")

    fixed_pdb_file = output_pdb_file.replace("_with_NAD.pdb", "_fixed.pdb")
    fixer = PDBFixer(filename=input_pdb_file)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Keep cofactors (remove only water)
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens()

    with open(fixed_pdb_file, "w") as fout:
        PDBFile.writeFile(fixer.topology, fixer.positions, fout)

    print(f"✅ Fixed structure saved: {fixed_pdb_file}")
    return fixed_pdb_file

def reattach_nad(original_pdb, fixed_pdb, output_pdb):
    """Reattach NAD residues from original PDB if missing after fixing"""
    with (
        open(original_pdb) as orig,
        open(fixed_pdb) as fixed,
        open(output_pdb, "w") as out,
    ):
        fixed_lines = fixed.readlines()
        nad_lines = [
            line
            for line in orig
            if "NAD" in line and line.startswith(("HETATM", "ATOM"))
        ]

        # Only append NAD lines if not already present
        if not any("NAD" in line for line in fixed_lines):
            fixed_lines += nad_lines
            print("🔁 NAD was missing — reattached to fixed structure.")
        else:
            print("✅ NAD already present, no reattachment needed.")

        out.writelines(fixed_lines)

    print(f"✅ Final receptor with NAD: {output_pdb}")
    return output_pdb

def main():
    """Main execution function"""
    print("=== Node 4: Fix structure with PDBFixer and reattach NAD ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get PDB ID from environment variable
    pdb_id = os.getenv("PARAM_PDB_ID")
    if not pdb_id:
        print("❌ Error: PARAM_PDB_ID environment variable is not set.")
        exit(1)
    
    input_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_A_NAD.pdb")
    output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_A_NAD_fixed_with_NAD.pdb")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 2 (node_02_extract_chain.py) first.")
        exit(1)
    
    # Fix structure
    fixed_pdb_file = fix_pdb_structure(input_pdb_file, output_pdb_file)
    
    # Reattach NAD if needed
    reattach_nad(input_pdb_file, fixed_pdb_file, output_pdb_file)

if __name__ == "__main__":
    main()

