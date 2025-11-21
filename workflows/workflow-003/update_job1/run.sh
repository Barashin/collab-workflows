#!/bin/bash
set -e

# ============================================================================
# Run stage: In-silico screening
# ============================================================================

# ----------------------------------------------------------------------------
# Docking processing
# ----------------------------------------------------------------------------

# Node 7: Prepare receptor for docking
echo "=== Node 7: Preparing receptor for docking ==="
python3 docking/node_07_prepare_receptor.py

# Node 8: In-silico screening using SMINA
echo "=== Node 8: Running in-silico screening using SMINA ==="
python3 docking/node_08_docking_screening.py

echo "✅ Run stage completed successfully!"

