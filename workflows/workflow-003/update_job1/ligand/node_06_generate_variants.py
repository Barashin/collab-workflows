#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 6: Generate ligand variants
Input: ligand/output/{ligand_name}.sdf
Output: ligand/output/ligand_library/ (multiple SDF files)
"""

import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")
LIGAND_LIBRARY_DIR = os.path.join(OUTPUT_DIR, "ligand_library")

def generate_variants(smiles):
    """Generate functional group variants of a molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    substitutions = {
        "hydroxyl_to_amine": ("[OX2H]", "N"),
        "hydroxyl_to_thiol": ("[OX2H]", "S"),
        "hydroxyl_to_halogen": ("[OX2H]", "Cl"),
        "amine_to_hydroxyl": ("[NX3;H2]", "O"),
        "carboxyl_to_amide": ("C(=O)[OH]", "C(=O)N"),
        "carboxyl_to_ester": ("C(=O)[OH]", "C(=O)OC"),
        "halogen_to_hydroxyl": ("[F,Cl,Br,I]", "O"),
        "halogen_to_amine": ("[F,Cl,Br,I]", "N"),
        "nitro_to_amine": ("[NX3](=O)=O", "N"),
    }

    results = {}
    for name, (pattern, replacement) in substitutions.items():
        patt = Chem.MolFromSmarts(pattern)
        repl = Chem.MolFromSmiles(replacement)
        if patt and repl:
            replaced = Chem.ReplaceSubstructs(mol, patt, repl, replaceAll=True)
            variants = [Chem.MolToSmiles(r, canonical=True) for r in replaced]
            results[name] = list(set(variants))
    return results

def save_variant_figure(base_smiles, variants, filename):
    """Save all variant molecules as an SVG figure."""
    mols = [Chem.MolFromSmiles(base_smiles)]
    labels = ["original"]
    for name, smiles_list in variants.items():
        for s in smiles_list:
            mols.append(Chem.MolFromSmiles(s))
            labels.append(name)

    n_cols = 3
    n_rows = (len(mols) + n_cols - 1) // n_cols

    drawer = rdMolDraw2D.MolDraw2DSVG(300 * n_cols, 300 * n_rows, 300, 300)
    drawer.DrawMolecules(mols, legends=labels)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")

    with open(filename, "w") as f:
        f.write(svg)
    print(f"Saved SVG figure as '{filename}'")

def generate_3d_structures(base_smiles, variants, out_dir):
    """Generate 3D coordinates for all unique variants and the original molecule."""
    os.makedirs(out_dir, exist_ok=True)
    base_mol = Chem.MolFromSmiles(base_smiles)
    base_smiles_canonical = Chem.MolToSmiles(base_mol, canonical=True)

    # --- Generate 3D for original molecule ---
    base_3d = Chem.AddHs(base_mol)
    print("\nGenerating 3D structure for original molecule...")
    if AllChem.EmbedMolecule(base_3d, AllChem.ETKDGv3()) == 0:
        AllChem.MMFFOptimizeMolecule(base_3d)
        base_path = os.path.join(out_dir, "original.sdf")
        writer = Chem.SDWriter(base_path)
        writer.write(base_3d)
        writer.close()
        print(f"Original molecule saved: {base_path}")
    else:
        print("Failed to embed original molecule.")

    # --- Generate 3D for variants ---
    generated = []
    for name, smiles_list in variants.items():
        for s in smiles_list:
            if s == base_smiles_canonical:
                continue  # skip identical
            mol = Chem.MolFromSmiles(s)
            mol = Chem.AddHs(mol)
            success = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if success == 0:
                AllChem.MMFFOptimizeMolecule(mol)
                out_path = os.path.join(out_dir, f"{name}.sdf")
                writer = Chem.SDWriter(out_path)
                writer.write(mol)
                writer.close()
                generated.append(out_path)
                print(f"3D variant saved: {out_path}")
            else:
                print(f"Failed to embed: {name} ({s})")

    if not generated:
        print("\nℹNo new variants required 3D generation.")
    else:
        print(f"\nGenerated {len(generated)} 3D structure(s) in '{out_dir}'.")

def sdf_to_smiles(sdf_filepath):
    """Extract SMILES from SDF file"""
    smiles_list = []
    with Chem.SDMolSupplier(sdf_filepath) as supplier:
        for mol in supplier:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol, canonical=True)
                smiles_list.append(smiles)
    return smiles_list

def main():
    """Main execution function"""
    print("=== Node 6: Generate ligand variants ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get parameters from environment variables
    ligand_name = os.getenv("PARAM_LIGAND_NAME")
    
    if not ligand_name:
        print("❌ Error: PARAM_LIGAND_NAME environment variable is not set.")
        exit(1)
    
    input_sdf = os.path.join(OUTPUT_DIR, f"{ligand_name}.sdf")
    
    if not os.path.exists(input_sdf):
        print(f"❌ Error: {input_sdf} not found.")
        print("   Please run Node 5 (node_05_extract_ligand.py) first.")
        exit(1)
    
    # Extract SMILES from SDF
    smiles_output = sdf_to_smiles(input_sdf)
    if len(smiles_output) == 0:
        print(f"❌ Error: Could not extract SMILES from {input_sdf}")
        exit(1)
    
    base_smiles = smiles_output[0]
    print(f"Base SMILES: {base_smiles}")
    
    # Generate variants
    variants = generate_variants(base_smiles)

    print("Generated variants:")
    for name, smiles_list in variants.items():
        print(f"\n{name}:")
        for s in smiles_list:
            print("  ", s)

    # Save 2D figure and 3D models
    variant_figure_path = os.path.join(OUTPUT_DIR, "variants.svg")
    save_variant_figure(base_smiles, variants, variant_figure_path)
    generate_3d_structures(base_smiles, variants, LIGAND_LIBRARY_DIR)

if __name__ == "__main__":
    main()

