#!/bin/bash
#SBATCH --job-name=cmap_agg
#SBATCH --account=vswarup_lab
#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=scripts/logs/step2_%j.out
#SBATCH --error=scripts/logs/step2_%j.err

# Step 2: Aggregate replicate signatures and pre-rank gene lists for one class.
# Submit from the repo root, selecting the class via CMAP_PERT_CLASS:
#   sbatch --export=ALL,CMAP_PERT_CLASS=knockdown scripts/slurm_step2.sh

set -euo pipefail

export CMAP_PERT_CLASS="${CMAP_PERT_CLASS:-drug}"

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Perturbation class: $CMAP_PERT_CLASS"

set +u
source ~/.bashrc
conda activate "${CONDA_ENV:-model-ad_env}"
set -u

export CMAP_PROJECT_ROOT="${CMAP_PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
mkdir -p "$CMAP_PROJECT_ROOT/scripts/logs"

cd "$CMAP_PROJECT_ROOT"
python scripts/02_aggregate_signatures.py

echo "Job finished: $(date)"
