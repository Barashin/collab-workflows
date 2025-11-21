#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 9: Generate docking results report
Input: docking/output/ (docking result files)
Output: report/output/docking_ranking.txt, report/output/*_docked.sdf
"""

import glob
import os
import re
import shutil
from pathlib import Path

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCKING_DIR = os.path.join(SCRIPT_DIR, "..", "docking", "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

def parse_smina_log(log_file):
    """
    Extract mode 1 affinity value from a Smina log file.

    Args:
        log_file (str): Path to the log file

    Returns:
        float or None: Binding energy value (kcal/mol), or None if not found
    """
    with open(log_file, "r") as f:
        for line in f:
            if line.strip().startswith("1 "):  # Get the line for mode 1
                parts = line.split()
                if len(parts) > 1:
                    try:
                        affinity = float(parts[1])  # Affinity value (kcal/mol)
                        return affinity
                    except ValueError:
                        pass
    return None  # No mode 1 data found

def generate_docking_ranking():
    """Generate docking result ranking"""
    print("\n=== Generating docking result ranking ===")

    # Get docking log files
    log_files = glob.glob(os.path.join(DOCKING_DIR, "*_docking.log"))

    # Extract affinity values from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_docking.log", "")
            results.append((compound_name, affinity))

    # Sort by affinity (lower = stronger binding)
    results.sort(key=lambda x: x[1])

    # Output results to a text file
    output_file = os.path.join(OUTPUT_DIR, "docking_ranking.txt")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (Strongest Binding First)\n")
        f.write("----------------------------------------\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            f.write(
                f"Rank {rank}: Compound {compound}, Binding Energy: {affinity:.2f} kcal/mol\n"
            )

    # Display summary
    print(f"Ranked docking results for {len(results)} compounds.")
    print(f"Results saved to {output_file}")

    # Display top 10 results
    print("\n=== Top 10 Compounds ===")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        print(f"Rank {rank}: Compound {compound}, Binding Energy: {affinity:.2f} kcal/mol")

def copy_top_compound():
    """Copy the top compound's file"""
    print("\n=== Copying Top Compound File ===")

    ranking_file = os.path.join(OUTPUT_DIR, "docking_ranking.txt")
    if not os.path.exists(ranking_file):
        print(f"❌ Error: {ranking_file} not found.")
        return

    # Read ranking file
    with open(ranking_file, encoding="utf-8") as f:
        text = f.read()

    # Extract top compound name
    m = re.search(r"Rank\s*1:.*Compound\s+([^\s,，]+)", text)
    if not m:
        print("⚠ Top compound name not found in ranking file.")
        return
    top_ligand = m.group(1)
    print(f"Top ligand: {top_ligand}")

    # Source SDF file path
    src = os.path.join(DOCKING_DIR, f"{top_ligand}_docked.sdf")

    if not os.path.exists(src):
        print(f"⚠ Top compound file not found: {src}")
        return

    # Destination directory (output)
    dst_dir = OUTPUT_DIR
    os.makedirs(dst_dir, exist_ok=True)

    # Destination path
    dst = os.path.join(dst_dir, os.path.basename(src))

    # Copy file
    shutil.copy(src, dst)

    print(f"Copied {src} → {dst}")

def main():
    """Main execution function"""
    print("=== Node 9: Generate docking results report ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Check docking results directory
    if not os.path.exists(DOCKING_DIR):
        print(f"❌ Error: {DOCKING_DIR} not found.")
        print("   Please run Node 8 (node_08_docking_screening.py) first.")
        exit(1)
    
    # Analyze docking results and generate ranking
    generate_docking_ranking()

    # Copy top-ranked compound file
    copy_top_compound()

    print("\n✅ Report generation complete.")

if __name__ == "__main__":
    main()

