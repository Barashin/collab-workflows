#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 12: Reporting - Generate report
Input: Docking score (docking result files from Node 11)
Output: Score report (score report)
"""

import glob
import os
import re
import shutil
from pathlib import Path

DEFAULT_PDB_ID = "5Y7J"
DEFAULT_CHAIN_ID = "A"

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "docking", "output")
PREPARATION_DIR = os.path.join(SCRIPT_DIR, "..", "preparation", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
RANKING_FILE = os.path.join(OUTPUT_DIR, "docking_ranking.txt")

def main():
    """Main execution function"""
    print("=== Node 12: Generate report ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 11 (node_11_smina_screening.py) first.")
        exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Input directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Generate docking result ranking
    results = generate_docking_ranking()
    
    # Copy top 3 compound files and generate complexes
    copy_top_compounds(results)
    
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


def generate_docking_ranking():
    """Generate docking result ranking"""
    print("\n--- Generate docking result ranking ---")
    
    # Get docking log files
    log_files = glob.glob(os.path.join(INPUT_DIR, "*_log.txt"))
    
    if not log_files:
        print(f"❌ Error: No log files found in {INPUT_DIR}.")
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


def copy_top_compounds(results):
    """Copy top 3 compound files and receptor PDB file separately"""
    print("\n--- Copy top 3 compound files and receptor ---")
    
    if not results:
        print("Warning: No results available for copying.")
        return
    
    # Get receptor file path
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    receptor_file = None
    if chain_id:
        receptor_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_chain_{chain_id}_clean.pdb")
    
    if not receptor_file or not os.path.exists(receptor_file):
        receptor_file = os.path.join(PREPARATION_DIR, f"{pdb_id}_clean.pdb")
    
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
        print(f"⚠ Warning: Receptor file not found: {receptor_file}")
    
    copied_count = 0
    
    for rank, (ligand_name, affinity) in enumerate(top_compounds, 1):
        print(f"\nRank {rank}: {ligand_name} (Binding energy: {affinity:.2f} kcal/mol)")
        
        # Copy docked ligand SDF file
        src_sdf = Path(INPUT_DIR) / f"{ligand_name}_docked.sdf"
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

