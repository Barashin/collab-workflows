#!/bin/bash
set -e

# ============================================================================
# Post-run stage: Report generation
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------------------------------------------------------
# Report processing
# ----------------------------------------------------------------------------

# Node 14: Generate docking results report
echo "=== Node 14: Generating docking results report ==="
python3 "${SCRIPT_DIR}/report/node_14_reporting.py"

echo "✅ Post-run stage completed successfully!"

