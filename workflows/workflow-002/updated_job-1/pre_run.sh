#!/bin/bash
set -e

# ============================================================================
# Pre-run stage: Data preparation
# ============================================================================

# ----------------------------------------------------------------------------
# Ligand processing
# ----------------------------------------------------------------------------

# Node 1: Download ligand files
echo "=== Node 1: Downloading ligand files ==="
python ligand/node_01_download_ligands.py

# Node 3: Unpack ligand ZIP file
echo "=== Node 3: Unpacking ligand ZIP file ==="
python ligand/node_02_unpack_ligands.py

# Node 5: Select ligands
echo "=== Node 5: Selecting ligands ==="
python ligand/node_03_ligand_selection.py

# Node 4: Prepare ligands (add hydrogens, assign charges, generate 3D structures)
echo "=== Node 4: Preparing ligands ==="
python ligand/node_04_prepare_ligands.py

# Node 6: Ligand view
echo "=== Node 6: Generating ligand view ==="
python ligand/node_10_ligand_view.py

# ----------------------------------------------------------------------------
# Protein processing
# ----------------------------------------------------------------------------

# Node 2: Download PDB file
echo "=== Node 2: Downloading PDB file ==="
python protein/node_05_download_pdb.py

# Node 4: Verify protein input file
echo "=== Node 4: Verifying protein input file ==="
python protein/node_06_protein_input.py

# ----------------------------------------------------------------------------
# Preparation processing
# ----------------------------------------------------------------------------

# Node 7: Ligand loading (extract ligand information from PDB)
echo "=== Node 7: Loading ligand information ==="
python preparation/node_07_ligand_loading.py

# Node 7: Extract chains
echo "=== Node 7: Extracting chains ==="
python preparation/node_08_extract_chains.py

# Node 6: Real ligand addition (combine prepared ligands with real ligand from PDB)
echo "=== Node 6: Adding real ligand ==="
python ligand/node_09_real_ligand_addition.py

# Node 8: Ligand center identification
echo "=== Node 8: Identifying ligand center ==="
python preparation/node_11_ligand_center_identification.py

# Node 9: Clean up protein structure
echo "=== Node 9: Cleaning up protein structure ==="
python preparation/node_12_protein_extraction.py

echo "✅ Pre-run stage completed successfully!"
