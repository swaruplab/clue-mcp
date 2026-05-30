#!/usr/bin/env python3
"""
Step 3: Build the DuckDB metadata index from every per-class metadata.parquet
that exists under data/processed/<class>/.

Produces one indexed `signatures` table spanning all built classes (with a
pert_class column) plus a `genes` lookup table. Safe to re-run after building
each class; it picks up whichever classes are present.
"""

import os
import time
from pathlib import Path

import duckdb
import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_config import CLASS_PERT_TYPES

DATA_DIR = os.environ.get(
    "CMAP_PROJECT_ROOT",
    str(Path(__file__).resolve().parent.parent),
)
PROC_DIR = f"{DATA_DIR}/data/processed"

print("=" * 60)
print("Step 3: Build DuckDB (all available classes)")
print("=" * 60)
t0 = time.time()

frames = []
for cls in CLASS_PERT_TYPES:
    meta_path = f"{PROC_DIR}/{cls}/metadata.parquet"
    if os.path.exists(meta_path):
        df = pd.read_parquet(meta_path)
        df["pert_class"] = cls
        frames.append(df)
        print(f"  + {cls}: {len(df):,} signatures")
    else:
        print(f"  - {cls}: not built (skipped)")

if not frames:
    raise SystemExit("No per-class metadata.parquet found; run steps 1-2 first.")

signatures = pd.concat(frames, axis=0, ignore_index=True)

db_path = f"{PROC_DIR}/cmap.duckdb"
if os.path.exists(db_path):
    os.remove(db_path)
con = duckdb.connect(db_path)

con.execute("CREATE TABLE signatures AS SELECT * FROM signatures")
con.execute("CREATE INDEX idx_class ON signatures(pert_class)")
con.execute("CREATE INDEX idx_pert ON signatures(pert_name)")
con.execute("CREATE INDEX idx_cell ON signatures(cell_iname)")
if "moa" in signatures.columns:
    con.execute("CREATE INDEX idx_moa ON signatures(moa)")
print(f"\n  Table 'signatures': {len(signatures):,} rows across {len(frames)} class(es)")

with open(f"{PROC_DIR}/gene_names.txt") as f:
    genes = [g.strip() for g in f if g.strip()]
gene_df = pd.DataFrame({"gene_idx": range(len(genes)), "gene_symbol": genes})
con.execute("CREATE TABLE genes AS SELECT * FROM gene_df")
con.execute("CREATE INDEX idx_gene ON genes(gene_symbol)")
print(f"  Table 'genes': {len(genes):,} genes")

summary = con.execute("""
    SELECT pert_class,
           COUNT(*) AS n_sigs,
           COUNT(DISTINCT pert_name) AS n_perturbagens,
           COUNT(DISTINCT cell_iname) AS n_cells
    FROM signatures GROUP BY pert_class ORDER BY pert_class
""").fetchdf()
print("\nPer-class summary:")
print(summary.to_string(index=False))

con.close()
print(f"\nStep 3 complete in {time.time()-t0:.1f}s -> {db_path} "
      f"({os.path.getsize(db_path)/1e6:.1f} MB)")
