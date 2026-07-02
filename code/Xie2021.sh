#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
echo "This dataset-specific entry script has been merged into ${SCRIPT_DIR}/main.sh."
echo "Please edit PROJECT_ROOT and DATASETS in main.sh, then run: bash main.sh"
