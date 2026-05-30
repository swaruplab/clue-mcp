#!/bin/bash
#SBATCH --job-name=cmap_duckdb
#SBATCH --account=vswarup_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=scripts/logs/step3_%j.out
#SBATCH --error=scripts/logs/step3_%j.err

# Step 3: Build the DuckDB metadata index.
# Submit from the repo root:  sbatch scripts/slurm_step3.sh

set -euo pipefail

echo "Job started: $(date)"
echo "Node: $(hostname)"

set +u
source ~/.bashrc
conda activate "${CONDA_ENV:-model-ad_env}"
set -u

export CMAP_PROJECT_ROOT="${CMAP_PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
mkdir -p "$CMAP_PROJECT_ROOT/scripts/logs"

cd "$CMAP_PROJECT_ROOT"
python scripts/03_build_duckdb.py

echo "Job finished: $(date)"
