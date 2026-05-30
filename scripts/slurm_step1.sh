#!/bin/bash
#SBATCH --job-name=cmap_parse
#SBATCH --account=vswarup_lab
#SBATCH --partition=hugemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=scripts/logs/step1_%j.out
#SBATCH --error=scripts/logs/step1_%j.err

# Step 1: Parse a LINCS 2020 Level-5 GCTX into per-signature parquet files for
# one perturbation class (drug | knockdown | overexpression).
#
# Submit from the repo root, selecting the class via CMAP_PERT_CLASS:
#   sbatch --export=ALL,CMAP_PERT_CLASS=drug          scripts/slurm_step1.sh
#   sbatch --export=ALL,CMAP_PERT_CLASS=knockdown     scripts/slurm_step1.sh
#   sbatch --export=ALL,CMAP_PERT_CLASS=overexpression scripts/slurm_step1.sh
#
# knockdown parses trt_sh + trt_xpr together and may need --partition=maxmem.

set -euo pipefail

export CMAP_PERT_CLASS="${CMAP_PERT_CLASS:-drug}"

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Perturbation class: $CMAP_PERT_CLASS"

# Activate conda env (override CONDA_ENV to use a different one).
# Disable nounset while sourcing bashrc: some site /etc/bashrc files reference
# unset vars (e.g. BASHRCSOURCED) and would abort the job under `set -u`.
set +u
source ~/.bashrc
conda activate "${CONDA_ENV:-model-ad_env}"
set -u

# Project root: dir where sbatch was invoked from (the repo root).
export CMAP_PROJECT_ROOT="${CMAP_PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
mkdir -p "$CMAP_PROJECT_ROOT/scripts/logs"

cd "$CMAP_PROJECT_ROOT"
python scripts/01_parse_gctx.py

echo "Job finished: $(date)"
