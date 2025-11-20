#!/bin/bash
set -e

# ============================================================================
# Run stage: In-silico screening
# ============================================================================

# ----------------------------------------------------------------------------
# Docking processing
# ----------------------------------------------------------------------------

# Node 13: SMINA - In-silico screening
echo "=== Node 13: Running SMINA in-silico screening ==="
python docking/node_13_smina_screening.py

echo "✅ Run stage completed successfully!"
