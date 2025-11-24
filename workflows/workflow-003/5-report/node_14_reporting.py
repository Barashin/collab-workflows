#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 12: Reporting - Generate report
Input: Docking score (docking result files from Node 11)
Output: Score report (score report)
"""

import glob
import json
import os
import re
import shutil
from pathlib import Path

DEFAULT_PDB_ID = "5Y7J"
DEFAULT_CHAIN_ID = "A"

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
                return params.get("pdb_id", DEFAULT_PDB_ID)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Silva mounts inputs from depends_on nodes to the current directory
# inputs = ["*.pdb", "docking_results"] means:
# - *.pdb files from 3-docking_preparation may be at root or ./outputs/*.pdb
# - docking_results directory from 4-docking may be at root or ./outputs/docking_results
# Try both root directory and outputs directory
INPUT_DIR_ROOT = os.path.join(SCRIPT_DIR, "docking_results")
INPUT_DIR_OUTPUTS = os.path.join(SCRIPT_DIR, "outputs", "docking_results")
PREPARATION_DIR_ROOT = SCRIPT_DIR  # Root of working directory
PREPARATION_DIR_OUTPUTS = os.path.join(SCRIPT_DIR, "outputs")  # outputs subdirectory
# Silva expects outputs in the directory specified in outputs = ["results"]
# Output should be in "results" directory (not "outputs/results")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "results")
RANKING_FILE = os.path.join(OUTPUT_DIR, "docking_ranking.txt")

def main():
    """Main execution function"""
    print("=== Node 12: Generate report ===")
    
    # Find docking_results directory in both root and outputs directories
    INPUT_DIR = None
    if os.path.exists(INPUT_DIR_ROOT):
        INPUT_DIR = INPUT_DIR_ROOT
        print(f"✓ Found docking_results in root directory")
    elif os.path.exists(INPUT_DIR_OUTPUTS):
        INPUT_DIR = INPUT_DIR_OUTPUTS
        print(f"✓ Found docking_results in outputs directory")
    
    if not INPUT_DIR or not os.path.exists(INPUT_DIR):
        print(f"❌ Error: docking_results directory not found.")
        print(f"   Searched in: {INPUT_DIR_ROOT}")
        print(f"   Searched in: {INPUT_DIR_OUTPUTS}")
        print("   Please run Node 11 (node_11_smina_screening.py) first.")
        exit(1)
    
    # Find preparation directory (for PDB files)
    PREPARATION_DIR = None
    if os.path.exists(PREPARATION_DIR_ROOT):
        # Check if there are PDB files in root
        root_pdb_files = glob.glob(os.path.join(PREPARATION_DIR_ROOT, "*.pdb"))
        if root_pdb_files:
            PREPARATION_DIR = PREPARATION_DIR_ROOT
            print(f"✓ Found PDB files in root directory")
    if not PREPARATION_DIR and os.path.exists(PREPARATION_DIR_OUTPUTS):
        outputs_pdb_files = glob.glob(os.path.join(PREPARATION_DIR_OUTPUTS, "*.pdb"))
        if outputs_pdb_files:
            PREPARATION_DIR = PREPARATION_DIR_OUTPUTS
            print(f"✓ Found PDB files in outputs directory")
    
    if not PREPARATION_DIR:
        PREPARATION_DIR = PREPARATION_DIR_OUTPUTS  # Default to outputs
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Input directory (docking_results): {INPUT_DIR}")
    print(f"Preparation directory (PDB files): {PREPARATION_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Generate docking result ranking
    results = generate_docking_ranking(INPUT_DIR)
    
    # Copy top 3 compound files and generate complexes
    copy_top_compounds(results, INPUT_DIR, PREPARATION_DIR)
    
    print("✅ Report generation complete.")


def parse_smina_log(log_file):
    """
    Extract mode1 affinity value from Smina log file.
    
    Args:
        log_file (str): Path to log file
    
    Returns:
        float or None: Binding energy value (kcal/mol) or None
    """
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):  # Get mode1 line
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            affinity = float(parts[1])  # Affinity value (kcal/mol)
                            return affinity
                        except ValueError:
                            pass
    except Exception as e:
        print(f"Warning: Error occurred while reading {log_file}: {e}")
    return None  # mode1 data not found


def generate_docking_ranking(input_dir):
    """Generate docking result ranking"""
    print("\n--- Generate docking result ranking ---")
    
    # Get docking log files
    log_files = glob.glob(os.path.join(input_dir, "*_log.txt"))
    
    if not log_files:
        print(f"❌ Error: No log files found in {input_dir}.")
        print(f"   Searched pattern: {os.path.join(input_dir, '*_log.txt')}")
        if os.path.exists(input_dir):
            available_files = os.listdir(input_dir)
            print(f"   Available files: {available_files[:20]}...")  # Show first 20 files
        exit(1)
    
    # Extract affinity from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))
    
    if not results:
        print("❌ Error: No valid docking results found.")
        exit(1)
    
    # Sort by affinity (binding energy) in ascending order (lower = stronger binding)
    results.sort(key=lambda x: x[1])
    
    # Write results to text file
    with open(RANKING_FILE, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (sorted by binding strength)\n")
        f.write("=" * 50 + "\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            f.write(
                f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol\n"
            )
    
    # Display results
    print(f"Ranked docking results for {len(results)} compounds")
    print(f"Results saved to {RANKING_FILE}")
    
    # Display top 10 results
    print("\n--- Top 10 compounds ---")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol")
    
    return results


def copy_top_compounds(results, input_dir, preparation_dir):
    """Copy top 3 compound files and receptor PDB file separately"""
    print("\n--- Copy top 3 compound files and receptor ---")
    
    if not results:
        print("Warning: No results available for copying.")
        return
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Find receptor file with multiple fallback options
    receptor_file = None
    
    # Priority 1: Chain-specific cleaned PDB (e.g., 4OHU_chain_A_clean.pdb)
    if chain_id:
        receptor_file = os.path.join(preparation_dir, f"{pdb_id}_chain_{chain_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 2: General cleaned PDB (e.g., 4OHU_clean.pdb)
    if not receptor_file or not os.path.exists(receptor_file):
        receptor_file = os.path.join(preparation_dir, f"{pdb_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 3: Search for any PDB file with _clean in the name
    if not receptor_file or not os.path.exists(receptor_file):
        clean_files = glob.glob(os.path.join(preparation_dir, f"{pdb_id}*_clean.pdb"))
        if clean_files:
            receptor_file = clean_files[0]
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 4: Search for any PDB file starting with pdb_id
    if not receptor_file or not os.path.exists(receptor_file):
        pdb_files = glob.glob(os.path.join(preparation_dir, f"{pdb_id}*.pdb"))
        if pdb_files:
            # Prefer files with "clean" in the name
            clean_files = [f for f in pdb_files if "clean" in os.path.basename(f).lower()]
            if clean_files:
                receptor_file = clean_files[0]
            else:
                receptor_file = pdb_files[0]
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Get top 3 compounds
    top_n = min(3, len(results))
    top_compounds = results[:top_n]
    
    # Destination directory
    dst_dir = Path(OUTPUT_DIR)
    dst_dir.mkdir(exist_ok=True)
    
    # Copy receptor PDB file (used for docking)
    if receptor_file and os.path.exists(receptor_file):
        receptor_dst = dst_dir / f"receptor_{os.path.basename(receptor_file)}"
        shutil.copy(receptor_file, receptor_dst)
        print(f"✓ Copied receptor PDB: {receptor_dst.name}")
    else:
        print(f"⚠ Warning: Receptor file not found in {preparation_dir}")
        print(f"   Searched for: {pdb_id}_chain_{chain_id}_clean.pdb, {pdb_id}_clean.pdb, {pdb_id}*_clean.pdb")
        if os.path.exists(preparation_dir):
            available_files = [f for f in os.listdir(preparation_dir) if f.endswith('.pdb')]
            print(f"   Available PDB files: {available_files[:10]}...")
    
    copied_count = 0
    
    for rank, (ligand_name, affinity) in enumerate(top_compounds, 1):
        print(f"\nRank {rank}: {ligand_name} (Binding energy: {affinity:.2f} kcal/mol)")
        
        # Copy docked ligand SDF file
        src_sdf = Path(input_dir) / f"{ligand_name}_docked.sdf"
        if src_sdf.exists():
            dst_sdf = dst_dir / f"top{rank}_{ligand_name}_docked.sdf"
            shutil.copy(src_sdf, dst_sdf)
            print(f"  ✓ Copied docked ligand: {dst_sdf.name}")
            copied_count += 1
        else:
            print(f"  ⚠ Warning: Docked SDF file not found: {src_sdf.name}")
    
    print(f"\n✅ Copied receptor PDB and {copied_count} top {top_n} ligand file(s) to {OUTPUT_DIR}/")




if __name__ == "__main__":
    main()

