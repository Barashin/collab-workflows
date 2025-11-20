#!/bin/bash
set -e

# ============================================================================
# Post-run stage: Report generation
# ============================================================================

# ----------------------------------------------------------------------------
# Report processing
# ----------------------------------------------------------------------------

# Node 12: Generate report
echo "=== Node 12: Generating report ==="
python report/node_14_reporting.py

echo "✅ Post-run stage completed successfully!"
