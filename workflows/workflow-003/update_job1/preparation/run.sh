#!/bin/bash
set -e

python3 node_07_ligand_loading.py
python3 node_08_extract_chains.py
python3 node_11_ligand_center_identification.py
python3 node_12_protein_extraction.py

