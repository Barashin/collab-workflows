#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 9: Ligand Center Identification - Generate docking configuration file
Input: 5Y7J_ligand.txt (from Node 7), 5Y7J_chain.pdb (from Node 8)
Output: config.txt (docking configuration file)
"""

import os

DEFAULT_PDB_ID = "5Y7J"
# Default ligand name (can be overridden by environment variable)
DEFAULT_LIGAND_NAME = "8OL"
# Default chain ID (can be overridden by environment variable, should match Node 9)
DEFAULT_CHAIN_ID = "A"
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "output")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
OUTPUT_CONFIG = os.path.join(OUTPUT_DIR, "config.txt")

def main():
    """Main execution function"""
    print("=== Node 9: Ligand Center Identification ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get chain selection (should match Node 9)
    selected_chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    input_ligand_list_file = os.path.join(INPUT_DIR, "ligand_list.txt")
    
    print(f"\nInput ligand list file: {input_ligand_list_file}")
    print(f"Output config file: {OUTPUT_CONFIG}")
    
    if not os.path.exists(input_ligand_list_file):
        print(f"\n❌ Error: {input_ligand_list_file} not found.")
        print("   Please run Node 7 (node_07_ligand_loading.py) first.")
        exit(1)
    
    try:
        # Generate config.txt from ligand_list.txt
        generate_config_file(input_ligand_list_file, OUTPUT_CONFIG, selected_chain_id)
        print(f"\n✅ Successfully generated docking configuration file: {OUTPUT_CONFIG}")
        
    except Exception as e:
        print(f"\n❌ Error occurred during config file generation: {e}")
        exit(1)


def generate_config_file(ligand_list_file, output_config_file, selected_chain_id):
    """
    Generate docking configuration file from ligand list file.
    Reads the ligand center coordinates from ligand_list.txt that binds to the selected chain.
    
    Args:
        ligand_list_file: Path to ligand_list.txt file
        output_config_file: Path to output config.txt file
        selected_chain_id: Chain ID selected in Node 9 protein extraction (should match the chain in the input PDB file)
    """
    # Get ligand selection from environment variables or use defaults
    selected_ligand_name = os.environ.get("LIGAND_NAME", DEFAULT_LIGAND_NAME)
    
    # Read ligand center coordinates from ligand_list.txt
    center_x = None
    center_y = None
    center_z = None
    current_chain = None
    current_ligand = None
    
    print(f"\nSearching for ligand in chain '{selected_chain_id}' (matching Node 9 protein extraction selection)...")
    
    # Priority 1: Find specified ligand in the selected chain
    if selected_chain_id:
        with open(ligand_list_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Check for chain header
                if line.startswith("Chain "):
                    current_chain = line.split()[1].rstrip(":")
                    continue
                
                # Check for ligand name
                if line.startswith("Ligand: "):
                    current_ligand = line.split(":", 1)[1].strip()
                    continue
                
                # Check for center coordinates line
                if line.startswith("Center coordinates"):
                    # Parse: "Center coordinates (x, y, z): 42.286, 13.914, 31.628"
                    if ":" in line and current_chain == selected_chain_id:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            ligand_match = (current_ligand and 
                                          current_ligand.upper() == selected_ligand_name.upper())
                            if ligand_match:
                                center_x, center_y, center_z = coords
                                print(f"✓ Found ligand '{current_ligand}' in chain '{current_chain}'")
                                break
    
    # Priority 2: If no matching ligand found, use any ligand in the selected chain
    if center_x is None and selected_chain_id:
        print(f"⚠ Ligand '{selected_ligand_name}' not found in chain '{selected_chain_id}'")
        print(f"   Searching for any ligand in chain '{selected_chain_id}'...")
        
        with open(ligand_list_file, "r", encoding="utf-8") as f:
            current_chain = None
            for line in f:
                line = line.strip()
                
                if line.startswith("Chain "):
                    current_chain = line.split()[1].rstrip(":")
                    continue
                
                if line.startswith("Ligand: "):
                    current_ligand = line.split(":", 1)[1].strip()
                    continue
                
                if line.startswith("Center coordinates") and current_chain == selected_chain_id:
                    if ":" in line:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            center_x, center_y, center_z = coords
                            print(f"✓ Using ligand '{current_ligand}' from chain '{current_chain}'")
                            break
    
    # Priority 3: If still no ligand found, use first ligand in the file (fallback)
    if center_x is None:
        print(f"⚠ No ligand found in chain '{selected_chain_id}'. Using first ligand from file as fallback...")
        with open(ligand_list_file, "r", encoding="utf-8") as f:
            current_chain = None
            for line in f:
                line = line.strip()
                
                if line.startswith("Chain "):
                    current_chain = line.split()[1].rstrip(":")
                    continue
                
                if line.startswith("Ligand: "):
                    current_ligand = line.split(":", 1)[1].strip()
                    continue
                
                if line.startswith("Center coordinates"):
                    if ":" in line:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            center_x, center_y, center_z = coords
                            print(f"✓ Using first ligand '{current_ligand}' from chain '{current_chain}' (fallback)")
                            break
    
    if center_x is None or center_y is None or center_z is None:
        raise ValueError("Could not parse center coordinates from ligand_list.txt")
    
    print(f"   Center coordinates (x, y, z): {center_x:.3f}, {center_y:.3f}, {center_z:.3f}")
    
    # Generate config file
    config_lines = [
        f"center_x = {center_x:.3f}",
        f"center_y = {center_y:.3f}",
        f"center_z = {center_z:.3f}",
        "size_x   = 15",  # Use fixed value
        "size_y   = 15",
        "size_z   = 15",
        "exhaustiveness = 8",  # Search intensity
        "num_modes = 5",  # Number of output poses
        "energy_range = 4",  # Maximum energy difference between output poses
    ]
    
    with open(output_config_file, "w", encoding="utf-8") as f:
        f.write("\n".join(config_lines) + "\n")


if __name__ == "__main__":
    main()

