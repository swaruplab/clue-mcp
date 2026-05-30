#!/bin/bash
# Submit the full build pipeline for all three perturbation classes, with SLURM
# dependencies so each class runs step1 -> step2, and a single step3 (DuckDB)
# runs once every class has finished.
#
# Run from the repo root (the directory holding scripts/, data/, the *.gctx):
#   bash scripts/build_all.sh
#
# Override the classes to build:  CLASSES="drug overexpression" bash scripts/build_all.sh
# knockdown is submitted to maxmem (trt_sh + trt_xpr together are large).

set -euo pipefail

CLASSES="${CLASSES:-drug knockdown overexpression}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p scripts/logs

step2_ids=()
for cls in $CLASSES; do
    part="hugemem"
    [ "$cls" = "knockdown" ] && part="maxmem"

    j1=$(sbatch --parsable --partition="$part" \
            --export=ALL,CMAP_PERT_CLASS="$cls" scripts/slurm_step1.sh)
    echo "[$cls] step1 -> job $j1 (partition $part)"

    j2=$(sbatch --parsable --partition="$part" --dependency=afterok:"$j1" \
            --export=ALL,CMAP_PERT_CLASS="$cls" scripts/slurm_step2.sh)
    echo "[$cls] step2 -> job $j2 (after $j1)"
    step2_ids+=("$j2")
done

dep=$(IFS=: ; echo "${step2_ids[*]}")
j3=$(sbatch --parsable --dependency=afterok:"$dep" scripts/slurm_step3.sh)
echo "step3 (DuckDB, all classes) -> job $j3 (after ${dep})"
echo "Done submitting. Track with: squeue -u \"$USER\""
