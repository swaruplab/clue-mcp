#!/usr/bin/env python3
"""
Step 2: Aggregate replicate signatures and pre-compute ranked gene lists for one
perturbation class (CMAP_PERT_CLASS; default 'drug').

Consensus key = (cmap_name, cell_iname, pert_type). Keeping pert_type in the key
keeps shRNA (trt_sh) and CRISPR (trt_xpr) signatures of the same gene distinct,
so the knockdown tool can still filter by method. Singletons (n_reps=1) are kept
as-is (they are already exemplar signatures); groups with >=2 reps take the
per-gene median z-score.

Outputs to data/processed/<class>/:
    zscore_matrix.npy    (n_sigs, n_genes) float32
    rank_matrix.npy      (n_sigs, n_genes) int16
    metadata.parquet     n_sigs rows (generic perturbation schema)
and the shared data/processed/gene_names.txt.
"""

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_config import CLASS_HAS_MOA, selected_class

DATA_DIR = os.environ.get(
    "CMAP_PROJECT_ROOT",
    str(Path(__file__).resolve().parent.parent),
)
PERT_CLASS = selected_class()
IN_DIR = f"{DATA_DIR}/data/intermediate/{PERT_CLASS}"
OUT_DIR = f"{DATA_DIR}/data/processed/{PERT_CLASS}"
SHARED_DIR = f"{DATA_DIR}/data/processed"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(f"Step 2: Aggregate & Pre-Rank  (class = {PERT_CLASS})")
print("=" * 60)

t0 = time.time()
print("\nLoading intermediate data...")
zscores = pd.read_parquet(f"{IN_DIR}/zscore_matrix.parquet")
print(f"  zscore_matrix: {zscores.shape}  ({time.time()-t0:.1f}s)")
sig_meta = pd.read_parquet(f"{IN_DIR}/sig_metadata.parquet")
print(f"  sig_metadata: {len(sig_meta):,} rows")

print("\nAggregating by (cmap_name, cell_iname, pert_type)...")
t1 = time.time()
records = []
groups = sig_meta.groupby(["cmap_name", "cell_iname", "pert_type"], dropna=False)
n_groups = len(groups)
n_singleton = n_multi = n_missing = 0

for i, ((pert_name, cell, pert_type), group) in enumerate(groups):
    if (i + 1) % 20000 == 0:
        print(f"  {i+1:,}/{n_groups:,} groups ({time.time()-t1:.1f}s)")

    available = [s for s in group["sig_id"].tolist() if s in zscores.columns]
    if not available:
        n_missing += 1
        continue
    if len(available) == 1:
        consensus = zscores[available[0]].values
        n_singleton += 1
    else:
        consensus = zscores[available].median(axis=1).values
        n_multi += 1

    method = group["method"].mode()
    records.append({
        "pert_name": pert_name,
        "cell_iname": cell,
        "pert_type": pert_type,
        "method": method.iloc[0] if len(method) else None,
        "n_reps": len(available),
        "dose_mode": group["pert_idose"].mode().iloc[0] if "pert_idose" in group and not group["pert_idose"].mode().empty else "NA",
        "time_mode": group["pert_itime"].mode().iloc[0] if "pert_itime" in group and not group["pert_itime"].mode().empty else "NA",
        "zscores": consensus,
    })

consensus_df = pd.DataFrame(records)
gene_names = zscores.index.tolist()
print(f"\nAggregation done in {time.time()-t1:.1f}s")
print(f"  consensus signatures: {len(consensus_df):,} "
      f"(singletons {n_singleton:,}, aggregated {n_multi:,}, missing {n_missing:,})")
print(f"  unique perturbagens: {consensus_df['pert_name'].nunique():,}")
print(f"  unique cell lines:   {consensus_df['cell_iname'].nunique():,}")

print("\nPre-computing rank matrix...")
t2 = time.time()
all_zscores = np.vstack(consensus_df["zscores"].values)
all_ranks = np.apply_along_axis(rankdata, 1, all_zscores).astype(np.int16)
print(f"  {all_ranks.shape} in {time.time()-t2:.1f}s")

# -- generic metadata table --
sig_table = consensus_df[
    ["pert_name", "cell_iname", "pert_type", "method", "n_reps", "dose_mode", "time_mode"]
].copy()
sig_table.insert(0, "pert_class", PERT_CLASS)
sig_table["sig_idx"] = range(len(sig_table))
# cmap_name kept as alias for backward compatibility.
sig_table["cmap_name"] = sig_table["pert_name"]

if CLASS_HAS_MOA.get(PERT_CLASS):
    compoundinfo = pd.read_csv(f"{DATA_DIR}/compoundinfo_beta.txt", sep="\t")
    cdedup = compoundinfo.drop_duplicates(subset=["cmap_name"])[["cmap_name", "moa", "target"]]
    sig_table = sig_table.merge(cdedup, left_on="pert_name", right_on="cmap_name", how="left",
                                suffixes=("", "_drop"))
    sig_table = sig_table.drop(columns=[c for c in sig_table.columns if c.endswith("_drop")])
else:
    sig_table["moa"] = pd.NA
    sig_table["target"] = pd.NA

print(f"\nSaving to {OUT_DIR}/ ...")
sig_table.to_parquet(f"{OUT_DIR}/metadata.parquet", index=False)
np.save(f"{OUT_DIR}/zscore_matrix.npy", all_zscores.astype(np.float32))
np.save(f"{OUT_DIR}/rank_matrix.npy", all_ranks)
with open(f"{SHARED_DIR}/gene_names.txt", "w") as f:
    f.write("\n".join(gene_names))
print(f"  metadata.parquet: {len(sig_table):,} rows")
print(f"  zscore_matrix.npy: {all_zscores.shape} (~{all_zscores.astype(np.float32).nbytes/1e9:.1f} GB)")
print(f"  rank_matrix.npy:   {all_ranks.shape} (~{all_ranks.nbytes/1e9:.1f} GB)")
print(f"  gene_names.txt:    {len(gene_names):,} genes (shared)")
print(f"\nStep 2 complete in {(time.time()-t0)/60:.1f} min")
