#!/bin/bash
set -e

# ============================================================================
# Pre-run stage: Data preparation
# ============================================================================

# ----------------------------------------------------------------------------
# Protein processing
# ----------------------------------------------------------------------------

# Node 1: Download PDB file
echo "=== Node 1: Downloading PDB file ==="
python3 protein/node_01_download_pdb.py

# ----------------------------------------------------------------------------
# Preparation processing
# ----------------------------------------------------------------------------

# Node 2: Extract Chain A with NAD cofactor
echo "=== Node 2: Extracting Chain A with NAD cofactor ==="
python3 preparation/node_02_extract_chain.py

# Node 3: Calculate ligand center and generate docking config
echo "=== Node 3: Calculating ligand center and generating docking config ==="
python3 preparation/node_03_ligand_center_config.py

# Node 4: Fix structure with PDBFixer and reattach NAD
echo "=== Node 4: Fixing structure with PDBFixer and reattaching NAD ==="
python3 preparation/node_04_fix_structure.py

# ----------------------------------------------------------------------------
# Ligand processing
# ----------------------------------------------------------------------------

# Node 5: Extract ligand from PDB file
echo "=== Node 5: Extracting ligand from PDB file ==="
python3 ligand/node_05_extract_ligand.py

# Node 6: Generate ligand variants
echo "=== Node 6: Generating ligand variants ==="
python3 ligand/node_06_generate_variants.py

echo "✅ Pre-run stage completed successfully!"

