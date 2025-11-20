#!/bin/bash
set -e

# ============================================================================
# Post-run stage: Report generation
# ============================================================================

# ----------------------------------------------------------------------------
# Report processing
# ----------------------------------------------------------------------------

# Node 14: Generate report
echo "=== Node 14: Generating report ==="
python report/node_14_reporting.py

echo "✅ Post-run stage completed successfully!"
