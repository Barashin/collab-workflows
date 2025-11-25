# Workflow-003: In-Silico Virtual Screening

## Overview

This workflow automates **virtual screening** of ligand compounds against a protein target using **SMINA**, an enhanced fork of AutoDock Vina.  
It integrates molecular structure extraction, ligand preparation, receptor preparation, docking, and result ranking into a reproducible computational pipeline managed by **Silva**.

The target system modeled here is based on **PDB ID: 4OHU** (default: 5Y7J), with the ligand of interest **2TK** (default: 8OL).

---

## Pipeline Structure

The workflow is organized into **5 main nodes**, each containing multiple Python scripts executed sequentially:

### Node 1: Protein Preparation (`1-protein_preparation`)

Prepares the target protein structure for docking.

**Scripts:**
- **`node_05_download_pdb.py`**: Downloads the PDB file from RCSB PDB database
- **`node_06_protein_input.py`**: Verifies the downloaded protein input file

**Process:**
1. Downloads the PDB structure file based on `pdb_id` parameter
2. Verifies the downloaded file exists and is valid

**Input:**
- `pdb_id` (from `global_params.json` or environment variable `PARAM_PDB_ID`)

**Output:**
```
outputs/
└── {pdb_id}.pdb
```

**Dependencies:** None (entry point)

---

### Node 2: Ligand Preparation (`2-ligand_preparation`)

Downloads, unpacks, selects, and prepares ligand compounds for docking.

**Scripts:**
- **`node_01_download_ligands.py`**: Downloads ligand library ZIP file
- **`node_02_unpack_ligands.py`**: Unpacks the ligand library
- **`node_03_ligand_selection.py`**: Selects specific ligands from the library
- **`node_04_prepare_ligands.py`**: Prepares ligands (adds hydrogens, assigns charges, minimizes energy)
- **`node_09_real_ligand_addition.py`**: Adds the real ligand from PDB to the selected compounds
- **`node_10_ligand_view.py`**: Generates visualization of ligands

**Process:**
1. Downloads a ligand library ZIP file from a specified URL
2. Unpacks the library to extract SDF files
3. Selects specific compounds based on criteria
4. Prepares each ligand using Open Babel (adds hydrogens, assigns charges, minimizes energy)
5. Adds the real ligand extracted from the PDB structure
6. Generates visualization files

**Input:**
- `*.pdb` from `1-protein_preparation`
- `real_ligand.sdf` from `3-docking_preparation`

**Output:**
```
outputs/
├── ligands.zip
├── constructed_library/
│   └── *.sdf
├── selected_compounds/
│   ├── *.sdf
│   └── real_ligand.sdf
├── ligand.csv
└── variants.svg
```

**Dependencies:** `1-protein_preparation`, `3-docking_preparation`

---

### Node 3: Docking Preparation (`3-docking_preparation`)

Prepares the docking configuration and cleans the protein structure.

**Scripts:**
- **`node_07_ligand_loading.py`**: Extracts ligand information from PDB file and generates `ligand_list.txt`
- **`node_08_extract_chains.py`**: Extracts specific protein chains and the real ligand from PDB
- **`node_11_ligand_center_identification.py`**: Generates `config.txt` for SMINA docking based on ligand center coordinates
- **`node_12_protein_extraction.py`**: Cleans up the protein structure using PDBFixer

**Process:**
1. Loads ligand information from the PDB file
2. Extracts specified protein chains (default: chain A) and the real ligand
3. Calculates the binding site center using the real ligand's coordinates
4. Generates docking configuration file (`config.txt`) with center coordinates and grid size
5. Cleans the protein structure using PDBFixer (adds missing residues, adds hydrogens, fixes structural issues)

**Input:**
- `*.pdb` from `1-protein_preparation`

**Output:**
```
outputs/
├── {pdb_id}_chain.pdb
├── {pdb_id}_chain_{chain_id}_clean.pdb
├── real_ligand.sdf
├── ligand_list.txt
└── config.txt
```

**Dependencies:** `1-protein_preparation`

**Configuration File (`config.txt`):**
```
center_x = <x_coordinate>
center_y = <y_coordinate>
center_z = <z_coordinate>
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 8
num_modes = 1
energy_range = 4.0
```

---

### Node 4: Docking (`4-docking`)

Performs molecular docking of each ligand against the prepared receptor using **SMINA**.

**Scripts:**
- **`node_13_smina_screening.py`**: Executes SMINA docking for each ligand in `selected_compounds/`

**Process:**
1. Verifies SMINA installation
2. Locates the receptor PDB file (cleaned protein structure)
3. Reads docking configuration from `config.txt`
4. Iterates through each ligand SDF file in `selected_compounds/`
5. Runs SMINA docking with parameters from `config.txt`
6. Extracts predicted binding affinities from docking logs

**Configuration:**
- Scoring: Vina
- Number of modes: 1 (configurable via `config.txt`)
- Exhaustiveness: 8 (configurable via `config.txt`)
- Energy range: 4.0 kcal/mol (configurable via `config.txt`)
- Grid center and size: From `config.txt`

**Input:**
- `*.pdb` from `3-docking_preparation` (cleaned protein structure)
- `selected_compounds/` from `2-ligand_preparation` (ligand SDF files)
- `config.txt` from `3-docking_preparation` (docking configuration)

**Output:**
```
outputs/
└── docking_results/
    ├── {ligand_name}_docked.sdf
    └── {ligand_name}_log.txt
```

**Dependencies:** `2-ligand_preparation`, `3-docking_preparation`

---

### Node 5: Report Generation (`5-report`)

Aggregates docking results and generates a ranked summary of ligand performance.

**Scripts:**
- **`node_14_reporting.py`**: Generates docking result ranking and copies top compounds

**Process:**
1. Parses each `*_log.txt` file to extract the top (mode 1) binding affinity value
2. Sorts ligands by binding energy (ascending order = stronger binding)
3. Generates a human-readable report in `results/docking_ranking.txt`
4. Copies the top 3 ranked docked ligand structures and receptor PDB to `results/`

**Input:**
- `*.pdb` from `3-docking_preparation` (receptor structure)
- `docking_results/` from `4-docking` (docking result files)

**Output:**
```
results/
├── docking_ranking.txt
├── receptor_{pdb_id}_chain_{chain_id}_clean.pdb
├── top1_{ligand_name}_docked.sdf
├── top2_{ligand_name}_docked.sdf
└── top3_{ligand_name}_docked.sdf
```

**Dependencies:** `3-docking_preparation`, `4-docking`

---

## Execution Flow

The workflow is executed using **Silva**, which automatically manages dependencies and execution order:

```
1-protein_preparation (Node 5-6)
    ↓
3-docking_preparation (Node 7-8, 11-12) ← 1-protein_preparation
    ↓                                    ↓
2-ligand_preparation (Node 1-4, 9-10) ← 1-protein_preparation, 3-docking_preparation
    ↓
4-docking (Node 13) ← 2-ligand_preparation, 3-docking_preparation
    ↓
5-report (Node 14) ← 3-docking_preparation, 4-docking
```

**Execution Order:**
1. **1-protein_preparation** runs first (no dependencies)
2. **3-docking_preparation** runs after 1-protein_preparation completes
3. **2-ligand_preparation** runs after both 1-protein_preparation and 3-docking_preparation complete
4. **4-docking** runs after both 2-ligand_preparation and 3-docking_preparation complete
5. **5-report** runs last, after 3-docking_preparation and 4-docking complete

---

## Configuration

### Global Parameters (`global_params.json`)

Located at `workflows/workflow-003/global_params.json`:

```json
{
  "pdb_id": "4OHU",
  "ligand_name": "2TK"
}
```

**Parameters:**
- `pdb_id`: The PDB ID of the target protein structure (default: "5Y7J")
- `ligand_name`: The three-letter residue name of the ligand to extract (default: "8OL")

### Environment Variables

Parameters can also be set via environment variables (takes precedence over `global_params.json`):
- `PARAM_PDB_ID` or `PDB_ID`: PDB ID
- `PARAM_LIGAND_NAME` or `LIGAND_NAME`: Ligand name
- `CHAIN_ID`: Chain ID(s) to extract (default: "A")

---

## Expected Runtime

| Stage | Duration (Approx.) |
|-------|--------------------|
| Protein Preparation | ~1-2 min |
| Ligand Preparation | ~5-10 min (depends on number of ligands) |
| Docking Preparation | ~2-3 min |
| Docking (per ligand) | ~2-5 min |
| Full Library (20-30 ligands) | ~10 min |
| Report Generation | <1 min |

**Total Runtime:** ~20 min for a typical screening of 20-30 ligands

---

## Output Structure

```
workflow-003/
├── global_params.json
├── 1-protein_preparation/
│   └── outputs/
│       └── {pdb_id}.pdb
├── 2-ligand_preparation/
│   └── outputs/
│       ├── ligands.zip
│       ├── constructed_library/
│       ├── selected_compounds/
│       │   ├── *.sdf
│       │   └── real_ligand.sdf
│       ├── ligand.csv
│       └── variants.svg
├── 3-docking_preparation/
│   └── outputs/
│       ├── {pdb_id}_chain.pdb
│       ├── {pdb_id}_chain_{chain_id}_clean.pdb
│       ├── real_ligand.sdf
│       ├── ligand_list.txt
│       └── config.txt
├── 4-docking/
│   └── outputs/
│       └── docking_results/
│           ├── {ligand_name}_docked.sdf
│           └── {ligand_name}_log.txt
└── 5-report/
    └── outputs/
        ├── docking_ranking.txt
        └── results/
            ├── receptor_{pdb_id}_chain_{chain_id}_clean.pdb
            ├── top1_{ligand_name}_docked.sdf
            ├── top2_{ligand_name}_docked.sdf
            └── top3_{ligand_name}_docked.sdf
```

---

## Docking Result Ranking

The final report (`outputs/results/docking_ranking.txt`) contains a ranked list of ligands sorted by binding energy (stronger binding = more negative energy):

```
Docking Result Ranking (sorted by binding strength)
==================================================
Rank 1: Compound {ligand_name}, Binding energy: -10.10 kcal/mol
Rank 2: Compound {ligand_name}, Binding energy: -9.90 kcal/mol
Rank 3: Compound {ligand_name}, Binding energy: -9.80 kcal/mol
...
```

**Interpretation:**
- Ligands with more negative binding energies exhibit stronger predicted affinity
- The top-ranked compounds are the most promising drug candidates
- Typical binding energies range from -5 to -12 kcal/mol

---

## Dependencies

All dependencies are installed inside the Docker container (`chiral.sakuracr.jp/smina:2025_11_06`):

- **Python 3.12**
- **RDKit** - Cheminformatics library for ligand manipulation
- **SMINA** - Enhanced AutoDock Vina for molecular docking
- **BioPython** - Biological computation library for PDB parsing
- **OpenMM / PDBFixer** - Protein structure cleaning and preparation
- **Open Babel** - Ligand preparation (hydrogen addition, charge assignment, energy minimization)
- **NumPy** - Numerical computing

---

## Running the Workflow

### Using Silva (Recommended)

The workflow is designed to run in **Silva**, which automatically:
- Manages Docker containers
- Resolves node dependencies
- Mounts input/output files between nodes
- Executes nodes in the correct order

1. Configure `global_params.json` with your target PDB ID and ligand name
2. Run the workflow in Silva TUI
3. Monitor execution progress
4. Review results in `5-report/outputs/results/`

### Manual Execution (Testing)

For testing individual nodes, you can run them manually:

```bash
# Activate conda environment
conda activate in_silico

# Run individual nodes
cd 1-protein_preparation
bash run.sh

cd ../3-docking_preparation
bash run.sh

cd ../2-ligand_preparation
bash run.sh

cd ../4-docking
bash run.sh

cd ../5-report
bash run.sh
```

---

## Notes

1. **File Mounting**: Silva mounts outputs from upstream nodes to the current node's working directory. Scripts check both the root directory and `outputs/` subdirectory to handle different mounting patterns.

2. **Ligand Preparation**: The `node_04_prepare_ligands.py` script uses Open Babel to prepare ligands. Failed preparations are logged but don't stop the workflow.

3. **Protein Cleaning**: PDBFixer may remove cofactors (like NAD) during cleaning. The workflow handles this by preserving important cofactors when possible.

4. **Docking Configuration**: The `config.txt` file is generated automatically based on the real ligand's coordinates. You can manually edit it to adjust the docking search space.

5. **SMINA Path**: In Docker containers, SMINA is located at `/usr/local/bin/smina`. For local testing on macOS, a local `smina.osx.12` binary is used if available.

---

## References

- **SMINA:** https://sourceforge.net/projects/smina/  
- **AutoDock Vina:** Trott & Olson, *J. Comput. Chem.* 31, 455–461 (2010).  
- **PDBFixer:** https://github.com/openmm/pdbfixer  
- **Open Babel:** http://openbabel.org/  
- **RDKit:** https://www.rdkit.org/  
- **BioPython:** https://biopython.org/
