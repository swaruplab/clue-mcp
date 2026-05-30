#!/usr/bin/env python3
"""
Step 1: Parse LINCS 2020 Level-5 GCTX into annotated DataFrames for one
perturbation class.

Select the class with CMAP_PERT_CLASS (drug | knockdown | overexpression);
default 'drug'. For multi-source classes (knockdown = trt_sh + trt_xpr) each
source GCTX is parsed and the exemplar signatures are concatenated, tagged with
pert_type and a sub-method label (shRNA / CRISPR).

Filters to exemplar signatures only (is_exemplar_sig == 1). For genetic
perturbations the perturbed entity is a gene; its symbol lives in cmap_name,
exactly like a compound name, so downstream aggregation is identical.

Outputs (per class) to data/intermediate/<class>/:
    zscore_matrix.parquet   genes x signatures
    sig_metadata.parquet    one row per signature
"""

import os
import time
from pathlib import Path

import pandas as pd
from cmapPy.pandasGEXpress import parse_gctx

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_config import (
    CLASS_PERT_TYPES,
    PERT_TYPE_SOURCES,
    selected_class,
)

DATA_DIR = os.environ.get(
    "CMAP_PROJECT_ROOT",
    str(Path(__file__).resolve().parent.parent),
)

PERT_CLASS = selected_class()
OUT_DIR = f"{DATA_DIR}/data/intermediate/{PERT_CLASS}"
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(f"Step 1: Parse GCTX  (class = {PERT_CLASS})")
print("=" * 60)

t0 = time.time()
print("\nLoading metadata...")
siginfo = pd.read_csv(f"{DATA_DIR}/siginfo_beta.txt", sep="\t", low_memory=False)
print(f"  siginfo_beta.txt: {len(siginfo):,} rows ({time.time()-t0:.1f}s)")

geneinfo = pd.read_csv(f"{DATA_DIR}/geneinfo_beta.txt", sep="\t")
gene_map = geneinfo.set_index("gene_id")["gene_symbol"].to_dict()

# Columns we carry forward (only those present in siginfo).
KEEP = [c for c in
        ["sig_id", "cmap_name", "cell_iname", "pert_type",
         "pert_idose", "pert_itime", "is_exemplar_sig"]
        if c in siginfo.columns]

zscore_frames = []
meta_frames = []

for pert_type in CLASS_PERT_TYPES[PERT_CLASS]:
    gctx_name, method = PERT_TYPE_SOURCES[pert_type]
    gctx_path = f"{DATA_DIR}/{gctx_name}"

    sub = siginfo[
        (siginfo["pert_type"] == pert_type) &
        (siginfo["is_exemplar_sig"] == 1)
    ].copy()
    sub["method"] = method
    print(f"\n[{pert_type}] exemplar signatures: {len(sub):,}  "
          f"(perturbagens: {sub['cmap_name'].nunique():,}, "
          f"cell lines: {sub['cell_iname'].nunique():,})")

    print(f"  Parsing {gctx_path} ...")
    t1 = time.time()
    gctoo = parse_gctx.parse(gctx_path, cid=sub["sig_id"].tolist())
    print(f"  parsed {gctoo.data_df.shape} in {time.time()-t1:.1f}s")

    # Keep the unique integer gene_id index here; map to symbols after concat.
    # Multiple gene_ids can share a symbol, so a symbol index is non-unique and
    # breaks the axis=1 alignment when >1 source is concatenated (sh + xpr).
    zscore_frames.append(gctoo.data_df)
    meta_frames.append(sub[KEEP + ["method"]])

# Concatenate sources (knockdown = sh + xpr). Columns = signatures.
zscores = pd.concat(zscore_frames, axis=1)
zscores.index = zscores.index.map(lambda x: gene_map.get(int(x), str(x)))
sig_meta = pd.concat(meta_frames, axis=0, ignore_index=True)
print(f"\nCombined: {zscores.shape[1]:,} signatures x {zscores.shape[0]:,} genes")

print(f"\nSaving to {OUT_DIR}/ ...")
t2 = time.time()
zscores.to_parquet(f"{OUT_DIR}/zscore_matrix.parquet")
sig_meta.to_parquet(f"{OUT_DIR}/sig_metadata.parquet")
print(f"  saved ({time.time()-t2:.1f}s)")

print(f"\nStep 1 complete in {(time.time()-t0)/60:.1f} min")
