#!/bin/bash
set -e

# ============================================================================
# Load global parameters from global_params.json
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLOBAL_PARAMS_FILE="${SCRIPT_DIR}/../global_params.json"

if [ -f "$GLOBAL_PARAMS_FILE" ]; then
    # Use Python to parse JSON and set environment variables
    export PARAM_PDB_ID=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['pdb_id'])")
    export PARAM_LIGAND_NAME=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['ligand_name'])")
    # Also set for compatibility
    export PDB_ID="${PARAM_PDB_ID}"
    export LIGAND_NAME="${PARAM_LIGAND_NAME}"
    echo "✅ Loaded parameters from global_params.json:"
    echo "   PARAM_PDB_ID=${PARAM_PDB_ID}"
    echo "   PARAM_LIGAND_NAME=${PARAM_LIGAND_NAME}"
else
    echo "⚠️  Warning: global_params.json not found at $GLOBAL_PARAMS_FILE"
    echo "   Please set PARAM_PDB_ID and PARAM_LIGAND_NAME environment variables manually."
fi

# ============================================================================
# Create input and outputs directories
# ============================================================================
INPUT_DIR="${SCRIPT_DIR}/input"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Clean previous outputs to avoid nested directories when Silva copies results
# This is critical to prevent outputs/outputs creation on re-runs
if [ -d "${OUTPUT_DIR}" ]; then
    echo "Cleaning previous outputs in ${OUTPUT_DIR}"
    rm -rf "${OUTPUT_DIR:?}/"*
fi

# ============================================================================
# Move files from Silva mounts to input/ directory
# ============================================================================
# Move PDB files from root directory (Silva mounts) to input/
if [ -d "${SCRIPT_DIR}" ]; then
    for pdb_file in "${SCRIPT_DIR}"/*.pdb; do
        if [ -f "${pdb_file}" ]; then
            mv "${pdb_file}" "${INPUT_DIR}/"
            echo "✓ Moved $(basename "${pdb_file}") to input/"
        fi
    done
fi

# Move docking_results directory from root directory (Silva mounts) to input/
if [ -d "${SCRIPT_DIR}/docking_results" ]; then
    mv "${SCRIPT_DIR}/docking_results" "${INPUT_DIR}/"
    echo "✓ Moved docking_results directory to input/"
fi

# ============================================================================
# Run nodes
# ============================================================================
python3 node_14_reporting.py

# ============================================================================
# Post-run cleanup: Fix nested outputs if present (outputs/outputs -> outputs)
# ============================================================================
if [ -d "${OUTPUT_DIR}/outputs" ]; then
    echo "⚠️ Detected nested outputs directory. Flattening structure..."
    # Move contents if not empty
    if [ "$(ls -A "${OUTPUT_DIR}/outputs")" ]; then
        mv "${OUTPUT_DIR}/outputs"/* "${OUTPUT_DIR}/"
    fi
    # Remove empty nested directory
    rm -rf "${OUTPUT_DIR}/outputs"
    echo "✓ Flattened outputs directory"
fi

