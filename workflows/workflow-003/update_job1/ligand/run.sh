#!/bin/bash
set -e

python3 node_01_download_ligands.py
python3 node_02_unpack_ligands.py
python3 node_03_ligand_selection.py
python3 node_04_prepare_ligands.py
python3 node_09_real_ligand_addition.py
python3 node_10_ligand_view.py

