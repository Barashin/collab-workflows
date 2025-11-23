#!/bin/bash
set -e

# ============================================================================
# Load global parameters from global_params.json
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Try update_job1/global_params.json first, then workflow-003/global_params.json
GLOBAL_PARAMS_FILE="${SCRIPT_DIR}/global_params.json"
if [ ! -f "$GLOBAL_PARAMS_FILE" ]; then
    GLOBAL_PARAMS_FILE="${SCRIPT_DIR}/../global_params.json"
fi

if [ -f "$GLOBAL_PARAMS_FILE" ]; then
    # Use Python to parse JSON and set environment variables
    export PARAM_PDB_ID=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['pdb_id'])")
    export PARAM_LIGAND_NAME=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['ligand_name'])")
    # Also set LIGAND_NAME for compatibility with some scripts
    export LIGAND_NAME="${PARAM_LIGAND_NAME}"
    export PDB_ID="${PARAM_PDB_ID}"
    echo "✅ Loaded parameters from global_params.json:"
    echo "   PARAM_PDB_ID=${PARAM_PDB_ID}"
    echo "   PARAM_LIGAND_NAME=${PARAM_LIGAND_NAME}"
    echo "   LIGAND_NAME=${LIGAND_NAME}"
    echo "   PDB_ID=${PDB_ID}"
else
    echo "⚠️  Warning: global_params.json not found at $GLOBAL_PARAMS_FILE"
    echo "   Please set PARAM_PDB_ID and PARAM_LIGAND_NAME environment variables manually."
fi

# ============================================================================
# Pre-run stage: Data preparation
# ============================================================================

# ----------------------------------------------------------------------------
# Protein processing
# ----------------------------------------------------------------------------
echo "=== Running protein processing ==="
cd "${SCRIPT_DIR}/protein"
bash run.sh
cd "${SCRIPT_DIR}"

# ----------------------------------------------------------------------------
# Preparation processing
# ----------------------------------------------------------------------------
echo "=== Running preparation processing ==="
cd "${SCRIPT_DIR}/preparation"
bash run.sh
cd "${SCRIPT_DIR}"

# ----------------------------------------------------------------------------
# Ligand processing
# ----------------------------------------------------------------------------
echo "=== Running ligand processing ==="
cd "${SCRIPT_DIR}/ligand"
bash run.sh
cd "${SCRIPT_DIR}"

echo "✅ Pre-run stage completed successfully!"

