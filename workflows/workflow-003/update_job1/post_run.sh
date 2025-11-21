#!/bin/bash
set -e

# ============================================================================
# Post-run stage: Report generation
# ============================================================================

# ----------------------------------------------------------------------------
# Report processing
# ----------------------------------------------------------------------------

# Node 9: Generate docking results report
echo "=== Node 9: Generating docking results report ==="
python3 report/node_09_reporting.py

echo "✅ Post-run stage completed successfully!"

