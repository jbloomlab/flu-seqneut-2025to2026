#!/bin/bash
# Run all titer boxplot notebooks and save figures to results/
# Usage: bash run_notebooks.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NOTEBOOKS_DIR="$SCRIPT_DIR/notebooks"

echo "Running titer boxplot analyses..."

python "$NOTEBOOKS_DIR/subclade_titers.py"
echo "  subclade_titers.py done"

python "$NOTEBOOKS_DIR/mutation_titers.py"
echo "  mutation_titers.py done"

echo "Figures saved to $SCRIPT_DIR/results/"
