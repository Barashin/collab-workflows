#!/bin/bash
set -e

python3 node_02_extract_chain.py
python3 node_03_ligand_center_config.py
python3 node_04_fix_structure.py

