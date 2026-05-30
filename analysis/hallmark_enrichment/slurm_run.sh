#!/bin/bash
#SBATCH --job-name=hallmark_enrichment
#SBATCH --account=vswarup_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=analysis/hallmark_enrichment/slurm_%j.out
#SBATCH --error=analysis/hallmark_enrichment/slurm_%j.err

# Run the cross-Hallmark CMap enrichment showcase (~12 min on 8 cores).
# Submit from the repo root:  sbatch analysis/hallmark_enrichment/slurm_run.sh

set -euo pipefail

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Memory requested: 64G"
echo "CPUs: 8"

source ~/.bashrc
conda activate "${CONDA_ENV:-model-ad_env}"

export CMAP_PROJECT_ROOT="${CMAP_PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
cd "$CMAP_PROJECT_ROOT"

python analysis/hallmark_enrichment/run_hallmark_enrichment.py

echo "Job finished: $(date)"
