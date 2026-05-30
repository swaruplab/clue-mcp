#!/bin/bash
#SBATCH --job-name=cmap_enrich
#SBATCH --account=vswarup_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=analysis/slurm_%j.out
#SBATCH --error=analysis/slurm_%j.err

# Run the per-gene-list CMap enrichment analysis on SLURM.
#
# Usage:
#   GENE_FILE=/path/to/genes.csv sbatch analysis/slurm_run.sh
#   (omit GENE_FILE to use the default test_genes.txt at the repo root)

set -euo pipefail

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Memory requested: 64G"

source ~/.bashrc
conda activate "${CONDA_ENV:-model-ad_env}"

export CMAP_PROJECT_ROOT="${CMAP_PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
cd "$CMAP_PROJECT_ROOT"

GENE_FILE="${GENE_FILE:-}"
if [ -n "$GENE_FILE" ]; then
    echo "Gene file: $GENE_FILE"
    python analysis/run_enrichment.py "$GENE_FILE"
else
    echo "No GENE_FILE specified; using default test_genes.txt"
    python analysis/run_enrichment.py
fi

echo "Job finished: $(date)"
