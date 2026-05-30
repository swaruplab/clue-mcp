# CMap Drug Enrichment Pipeline: GCTX → API → MCP

## End-to-End Guide — LINCS 2020 Build

---

## Data Inventory (What We Have)

All files are located in `/dfs7/swaruplab/shared_lab/Clue_database/`.

### Level 5 GCTX Files (Replicate-Collapsed Signatures)

Downloaded from `s3://macchiato.clue.io/builds/LINCS2020/level5/`.

| File | Pert Type | Signatures x Genes | Size |
|------|-----------|---------------------|------|
| `level5_beta_trt_cp_n720216x12328.gctx` | Compound treatments | 720,216 x 12,328 | 34 GB |
| `level5_beta_trt_sh_n238351x12328.gctx` | shRNA knockdowns | 238,351 x 12,328 | 11 GB |
| `level5_beta_trt_xpr_n142901x12328.gctx` | CRISPR perturbations | 142,901 x 12,328 | 6.1 GB |
| `level5_beta_ctl_n58022x12328.gctx` | Vehicle controls | 58,022 x 12,328 | 2.7 GB |
| `level5_beta_trt_oe_n34171x12328.gctx` | Overexpression | 34,171 x 12,328 | 1.6 GB |
| `level5_beta_trt_misc_n8283x12328.gctx` | Miscellaneous | 8,283 x 12,328 | 390 MB |

**Total: ~1,201,944 signatures across 12,328 genes (~54 GB)**

### Metadata Files

Downloaded from `s3://macchiato.clue.io/builds/LINCS2020/`.

| File | Rows | Description |
|------|------|-------------|
| `siginfo_beta.txt` (444 MB) | 1,201,944 | Signature metadata -- maps `sig_id` to compound, dose, cell line, time, QC metrics |
| `geneinfo_beta.txt` (1.1 MB) | 12,328 | Gene info -- maps `gene_id` to `gene_symbol`, `ensembl_id`, `feature_space` (landmark vs inferred) |
| `compoundinfo_beta.txt` (4.5 MB) | 39,321 | Compound metadata -- `pert_id`, `cmap_name`, `moa`, `target`, `canonical_smiles`, `inchi_key` |
| `cellinfo_beta.txt` (38 KB) | 240 | Cell line metadata -- lineage, disease, donor info |
| `README.txt` (2.2 KB) | -- | Official dataset documentation |

### Helper Scripts

| File | Description |
|------|-------------|
| `cmap_downloads` | Bash script used to download the Level 5 GCTX files |
| `links2020_level5_urls.txt` | S3 URLs for all six Level 5 GCTX files |

---

## Metadata Column Schemas

### siginfo_beta.txt (Key Columns)

| Column | Description |
|--------|-------------|
| `sig_id` | Unique signature identifier |
| `pert_id` | CMap perturbagen BRD ID |
| `pert_type` | Perturbation type: `trt_cp`, `trt_sh`, `trt_oe`, `trt_xpr`, `trt_misc`, `ctl` |
| `cmap_name` | Perturbagen name (drug name for `trt_cp`) |
| `cell_iname` | Cell line name |
| `pert_dose` / `pert_idose` | Dose (numeric / with unit) |
| `pert_time` / `pert_itime` | Treatment duration (numeric / with unit) |
| `nsample` | Number of replicates collapsed into this signature |
| `tas` | Transcriptional Activity Score |
| `cc_q75` | 75th percentile replicate correlation |
| `is_hiq` | High-quality signature flag |
| `qc_pass` | QC pass flag |
| `is_exemplar_sig` | 1 = best replicate per condition (use these for analysis) |

### geneinfo_beta.txt

| Column | Description |
|--------|-------------|
| `gene_id` | Entrez gene ID |
| `gene_symbol` | HGNC gene symbol |
| `ensembl_id` | Ensembl gene ID |
| `gene_title` | Full gene name |
| `gene_type` | e.g. protein-coding, ncRNA |
| `feature_space` | `lm` (978 landmark genes) or `inferred` (~11,350 computationally inferred) |

### compoundinfo_beta.txt

| Column | Description |
|--------|-------------|
| `pert_id` | CMap BRD perturbagen ID |
| `cmap_name` | Compound name |
| `target` | Protein target(s) |
| `moa` | Mechanism of action |
| `canonical_smiles` | SMILES structure |
| `inchi_key` | InChI key |
| `compound_aliases` | Alternative names |

**Note:** One compound may have multiple rows (one per MOA/target/structure combination).

### cellinfo_beta.txt (Key Columns)

| Column | Description |
|--------|-------------|
| `cell_iname` | Cell line name (primary key) |
| `cell_type` | Cell type |
| `cell_lineage` | Tissue lineage |
| `primary_disease` | Associated disease |
| `donor_sex` | Donor sex |
| `growth_medium` | Culture medium |

---

## Architecture Overview

```
LINCS 2020 Level 5 GCTX files (already downloaded)
    |
[Step 1] Parse trt_cp GCTX with cmapPy -> pandas DataFrames
    |
[Step 2] Filter exemplars, aggregate, pre-rank -> Parquet + numpy
    |
[Step 3] Index into DuckDB
    |
[Step 4] Build enrichment engine (connectivity scoring)
    |
[Step 5] Wrap in FastAPI
    |
[Step 6] Wrap in MCP server
    |
[Step 7] Deploy: GitHub + Zenodo + HuggingFace Space
```

---

## Step 0: Data Already Downloaded

All raw data is in `/dfs7/swaruplab/shared_lab/Clue_database/`.

The Level 5 GCTX files and metadata were downloaded from the LINCS 2020 S3 bucket. **No GEO download needed** -- this is the newer, larger dataset (Dec 2020 release, ~1.2M signatures vs the older GSE92742's ~473K).

If you need to re-download metadata or additional levels:

```bash
cd /dfs7/swaruplab/shared_lab/Clue_database

# Metadata (already present)
wget -c https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt
wget -c https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt
wget -c https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt
wget -c https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt

# Optional: field definitions reference
wget -c "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/LINCS2020%20Release%20Metadata%20Field%20Definitions.xlsx"

# Optional: instance-level metadata (Level 3/4 profiles, large file)
wget -c https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/instinfo_beta.txt
```

---

## Step 1: Parse GCTX Files

```bash
pip install cmapPy pandas pyarrow numpy scipy
```

```python
# scripts/01_parse_gctx.py
"""
Parse LINCS 2020 Level 5 GCTX files into annotated pandas DataFrames.
GCTX = HDF5 format. Rows = 12,328 genes, Cols = signatures.
Level 5 = moderated z-scores (replicate-collapsed) -- what we want for enrichment.

Key difference from older GSE92742 data:
  - LINCS 2020 files are already split by pert_type (no filtering needed)
  - Column names changed: gene_id (not pr_gene_id), gene_symbol (not pr_gene_symbol)
  - Use siginfo_beta.txt (not sig_info.txt)
  - Filter on is_exemplar_sig == 1 to get best replicate per condition
"""

import pandas as pd
import numpy as np
from cmapPy.pandasGEOplus import parse_gctx

DATA_DIR = "/dfs7/swaruplab/shared_lab/Clue_database"

# -- Load metadata --
siginfo = pd.read_csv(f"{DATA_DIR}/siginfo_beta.txt", sep="\t", low_memory=False)
geneinfo = pd.read_csv(f"{DATA_DIR}/geneinfo_beta.txt", sep="\t")
compoundinfo = pd.read_csv(f"{DATA_DIR}/compoundinfo_beta.txt", sep="\t")

print(f"Total signatures in metadata: {len(siginfo)}")
print(f"Perturbation types: {siginfo['pert_type'].value_counts().to_dict()}")

# -- Filter for drug (trt_cp) exemplar signatures --
# is_exemplar_sig == 1 gives the single best replicate per treatment condition
drug_sigs = siginfo[
    (siginfo["pert_type"] == "trt_cp") &
    (siginfo["is_exemplar_sig"] == 1)
].copy()
print(f"Drug exemplar signatures: {len(drug_sigs)}")

# -- Parse the trt_cp GCTX (this is the big one -- 34 GB, may need ~16+ GB RAM) --
# The file already contains ONLY compound treatment signatures
gctx_path = f"{DATA_DIR}/level5_beta_trt_cp_n720216x12328.gctx"

# Option A: Load only exemplar signatures (recommended -- much smaller)
sig_ids = drug_sigs["sig_id"].tolist()
gctoo = parse_gctx.parse(gctx_path, cid=sig_ids)

# Option B: If memory allows, load all trt_cp signatures
# gctoo = parse_gctx.parse(gctx_path)

# gctoo.data_df: rows = genes (gene_id), cols = signatures (sig_id)
# Values = moderated z-scores
print(f"Matrix shape: {gctoo.data_df.shape}")

# -- Map gene IDs to symbols --
# LINCS 2020 uses gene_id / gene_symbol (not pr_gene_id / pr_gene_symbol)
gene_map = geneinfo.set_index("gene_id")["gene_symbol"].to_dict()
gctoo.data_df.index = gctoo.data_df.index.map(lambda x: gene_map.get(int(x), str(x)))

# Save intermediate
gctoo.data_df.to_parquet("data/intermediate/zscore_matrix.parquet")
drug_sigs.to_parquet("data/intermediate/drug_sig_metadata.parquet")

print("Step 1 complete: z-score matrix saved")
```

---

## Step 2: Aggregate and Pre-Rank Signatures

```python
# scripts/02_aggregate_signatures.py
"""
Aggregate replicate signatures per compound and pre-compute ranked gene lists.
This is the key preprocessing step that makes queries fast.

LINCS 2020 note: if using exemplar signatures (is_exemplar_sig == 1),
these are already the best single replicate per condition. You may still
want to aggregate across cell lines or doses.
"""

import pandas as pd
import numpy as np

# -- Load --
zscores = pd.read_parquet("data/intermediate/zscore_matrix.parquet")
sig_meta = pd.read_parquet("data/intermediate/drug_sig_metadata.parquet")

# LINCS 2020 compound metadata
compoundinfo = pd.read_csv(
    "/dfs7/swaruplab/shared_lab/Clue_database/compoundinfo_beta.txt", sep="\t"
)

# -- Aggregate: median z-score per compound x cell line --
# Group by (cmap_name, cell_iname) to get one consensus signature per drug-cell combo
# Note: LINCS 2020 uses cmap_name (not pert_iname) and cell_iname (not cell_id)

sig_meta_indexed = sig_meta.set_index("sig_id")

records = []
groups = sig_meta.groupby(["cmap_name", "cell_iname"])

for (drug, cell), group in groups:
    sig_ids = group["sig_id"].tolist()
    available = [s for s in sig_ids if s in zscores.columns]
    if len(available) < 2:
        continue  # skip singletons -- unreliable

    consensus = zscores[available].median(axis=1)

    records.append({
        "cmap_name": drug,
        "cell_iname": cell,
        "n_reps": len(available),
        "dose_mode": group["pert_idose"].mode().iloc[0] if "pert_idose" in group else "NA",
        "time_mode": group["pert_itime"].mode().iloc[0] if "pert_itime" in group else "NA",
        "zscores": consensus.values,
    })

consensus_df = pd.DataFrame(records)
gene_names = zscores.index.tolist()

print(f"Consensus signatures: {len(consensus_df)}")
print(f"Unique compounds: {consensus_df['cmap_name'].nunique()}")
print(f"Unique cell lines: {consensus_df['cell_iname'].nunique()}")

# -- Pre-compute ranks (ascending: most downregulated = rank 1) --
all_zscores = np.vstack(consensus_df["zscores"].values)  # (n_sigs, 12328)
from scipy.stats import rankdata
all_ranks = np.apply_along_axis(rankdata, 1, all_zscores).astype(np.int16)

# -- Build final tables --

# Signature metadata table
sig_table = consensus_df[["cmap_name", "cell_iname", "n_reps", "dose_mode", "time_mode"]].copy()
sig_table["sig_idx"] = range(len(sig_table))

# Merge in MOA and target info from compoundinfo_beta
# Note: compoundinfo may have multiple rows per compound (one per MOA/target combo)
compoundinfo_dedup = compoundinfo.drop_duplicates(subset=["cmap_name"])
sig_table = sig_table.merge(
    compoundinfo_dedup[["cmap_name", "moa", "target"]],
    on="cmap_name",
    how="left"
)

# Save everything
sig_table.to_parquet("data/processed/signature_metadata.parquet", index=False)
np.save("data/processed/zscore_matrix.npy", all_zscores.astype(np.float32))
np.save("data/processed/rank_matrix.npy", all_ranks)

with open("data/processed/gene_names.txt", "w") as f:
    f.write("\n".join(gene_names))

print("Step 2 complete: aggregated signatures and ranks saved")
print(f"  zscore_matrix.npy: {all_zscores.shape}")
print(f"  rank_matrix.npy: {all_ranks.shape}")
print(f"  signature_metadata.parquet: {len(sig_table)} rows")
```

---

## Step 3: Index into DuckDB

```python
# scripts/03_build_duckdb.py
"""
Build a DuckDB database for fast metadata queries.
The actual enrichment runs against numpy arrays (faster),
but DuckDB handles all the filtering and lookup queries.
"""

import duckdb
import pandas as pd

con = duckdb.connect("data/processed/cmap.duckdb")

# -- Signature metadata --
sig_meta = pd.read_parquet("data/processed/signature_metadata.parquet")
con.execute("CREATE TABLE signatures AS SELECT * FROM sig_meta")

# -- Create indexes for common queries --
con.execute("CREATE INDEX idx_pert ON signatures(cmap_name)")
con.execute("CREATE INDEX idx_cell ON signatures(cell_iname)")
con.execute("CREATE INDEX idx_moa ON signatures(moa)")

# -- Gene lookup table --
with open("data/processed/gene_names.txt") as f:
    genes = [g.strip() for g in f.readlines()]

gene_df = pd.DataFrame({"gene_idx": range(len(genes)), "gene_symbol": genes})
con.execute("CREATE TABLE genes AS SELECT * FROM gene_df")
con.execute("CREATE INDEX idx_gene ON genes(gene_symbol)")

# -- Summary stats --
result = con.execute("""
    SELECT 
        COUNT(*) as n_signatures,
        COUNT(DISTINCT cmap_name) as n_compounds,
        COUNT(DISTINCT cell_iname) as n_cell_lines,
        COUNT(DISTINCT moa) as n_moas
    FROM signatures
""").fetchone()

print(f"Database built: {result[0]} signatures, {result[1]} compounds, "
      f"{result[2]} cell lines, {result[3]} MOAs")

con.close()
```

---

## Step 4: Build the Enrichment Engine

```python
# cmap_enrichment/engine.py
"""
Core enrichment engine: given up/down gene sets, find drugs that 
reverse the signature (negative connectivity = therapeutic candidates).

Implements the Weighted Connectivity Score (WTCS) from Subramanian et al. 2017.
"""

import numpy as np
import pandas as pd
import duckdb
from pathlib import Path
from typing import Optional

DATA_DIR = Path(__file__).parent.parent / "data" / "processed"


class CMapEngine:
    def __init__(self, data_dir: str = None):
        self.data_dir = Path(data_dir) if data_dir else DATA_DIR
        
        # Load pre-computed arrays into memory
        # For LINCS 2020 trt_cp exemplars, expect rank_matrix to be larger than the old dataset
        self.rank_matrix = np.load(self.data_dir / "rank_matrix.npy")    # (n_sigs, n_genes)
        self.zscore_matrix = np.load(self.data_dir / "zscore_matrix.npy")
        
        with open(self.data_dir / "gene_names.txt") as f:
            self.gene_names = [g.strip() for g in f.readlines()]
        self.gene_to_idx = {g: i for i, g in enumerate(self.gene_names)}
        
        self.n_sigs, self.n_genes = self.rank_matrix.shape
        self.con = duckdb.connect(str(self.data_dir / "cmap.duckdb"), read_only=True)
    
    def _ks_enrichment(self, ranks_col: np.ndarray, gene_indices: np.ndarray) -> float:
        """
        KS-like enrichment score for a single signature against a gene set.
        Adapted from the CMap WTCS method.
        
        ranks_col: rank vector for one signature (length = n_genes)
        gene_indices: indices of genes in the query set
        """
        n = len(ranks_col)
        n_set = len(gene_indices)
        if n_set == 0:
            return 0.0
        
        # Sort gene set by their rank in this signature
        set_ranks = np.sort(ranks_col[gene_indices])
        
        # Vectorized ES computation
        positions = np.arange(1, n_set + 1)
        hit_scores = set_ranks / n  # normalized rank positions
        
        # ES = max deviation from expected uniform distribution
        expected = positions / n_set
        deviations_up = expected - hit_scores
        
        es = deviations_up[np.argmax(np.abs(deviations_up))]
        return float(es)
    
    def _wtcs(self, ranks_row: np.ndarray, up_idx: np.ndarray, down_idx: np.ndarray) -> float:
        """
        Weighted Connectivity Score.
        Positive WTCS = signature mimics the query.
        Negative WTCS = signature REVERSES the query (what we want for drug discovery).
        """
        es_up = self._ks_enrichment(ranks_row, up_idx) if len(up_idx) > 0 else 0.0
        es_down = self._ks_enrichment(ranks_row, down_idx) if len(down_idx) > 0 else 0.0
        
        if es_up == 0 and es_down == 0:
            return 0.0
        if np.sign(es_up) != np.sign(es_down):
            return 0.0  # discordant -- not interpretable
        
        return (es_up - es_down) / 2
    
    def query_enrichment(
        self,
        genes_up: list[str],
        genes_down: list[str] = None,
        cell_line: str = None,
        top_n: int = 50,
    ) -> pd.DataFrame:
        """
        Main query: find drugs that reverse your gene signature.
        
        Parameters:
            genes_up: upregulated genes in your disease/condition
            genes_down: downregulated genes (optional but improves specificity)
            cell_line: filter to specific cell line (e.g., "A549", "MCF7")
            top_n: number of top results to return
        
        Returns:
            DataFrame with columns: cmap_name, cell_iname, wtcs, moa, target, n_reps
            Sorted by WTCS ascending (most negative = strongest reversal)
        """
        # Map gene symbols to indices
        up_idx = np.array([self.gene_to_idx[g] for g in genes_up if g in self.gene_to_idx])
        down_idx = np.array([self.gene_to_idx[g] for g in (genes_down or []) if g in self.gene_to_idx])
        
        # Report unmapped genes
        unmapped_up = [g for g in genes_up if g not in self.gene_to_idx]
        unmapped_down = [g for g in (genes_down or []) if g not in self.gene_to_idx]
        
        if len(up_idx) == 0 and len(down_idx) == 0:
            raise ValueError(f"No query genes found in L1000 landmark genes. "
                           f"Unmapped: {unmapped_up + unmapped_down}")
        
        # Optional cell line filter
        if cell_line:
            mask = self.con.execute(
                "SELECT sig_idx FROM signatures WHERE cell_iname = ?", [cell_line]
            ).fetchnumpy()["sig_idx"]
        else:
            mask = np.arange(self.n_sigs)
        
        # Compute WTCS for all (filtered) signatures
        scores = np.zeros(len(mask))
        for i, sig_idx in enumerate(mask):
            if len(down_idx) > 0:
                scores[i] = self._wtcs(self.rank_matrix[sig_idx], up_idx, down_idx)
            else:
                scores[i] = self._ks_enrichment(self.rank_matrix[sig_idx], up_idx)
        
        # Get top reversals (most negative scores)
        top_indices = np.argsort(scores)[:top_n]
        
        # Build results
        result_indices = mask[top_indices]
        results = self.con.execute(
            f"SELECT * FROM signatures WHERE sig_idx IN ({','.join(map(str, result_indices))})"
        ).fetchdf()
        
        # Add scores
        score_map = dict(zip(result_indices, scores[top_indices]))
        results["wtcs"] = results["sig_idx"].map(score_map)
        results = results.sort_values("wtcs", ascending=True)
        
        # Add query stats
        results.attrs["n_up_mapped"] = len(up_idx)
        results.attrs["n_down_mapped"] = len(down_idx)
        results.attrs["n_up_unmapped"] = len(unmapped_up)
        results.attrs["n_down_unmapped"] = len(unmapped_down)
        results.attrs["unmapped_genes"] = unmapped_up + unmapped_down
        
        return results
    
    def get_compound_info(self, compound_name: str) -> pd.DataFrame:
        """Get all signatures and metadata for a compound."""
        return self.con.execute(
            "SELECT * FROM signatures WHERE cmap_name = ?", [compound_name]
        ).fetchdf()
    
    def get_signature(self, sig_idx: int) -> pd.Series:
        """Get z-score vector for a specific signature."""
        return pd.Series(
            self.zscore_matrix[sig_idx],
            index=self.gene_names,
            name=f"sig_{sig_idx}"
        )
    
    def list_cell_lines(self) -> list[str]:
        """List available cell lines."""
        return self.con.execute(
            "SELECT DISTINCT cell_iname FROM signatures ORDER BY cell_iname"
        ).fetchdf()["cell_iname"].tolist()
    
    def list_moas(self) -> list[str]:
        """List available mechanisms of action."""
        return self.con.execute(
            "SELECT DISTINCT moa FROM signatures WHERE moa IS NOT NULL ORDER BY moa"
        ).fetchdf()["moa"].tolist()
    
    def search_compounds(self, query: str) -> pd.DataFrame:
        """Fuzzy search for compound names."""
        return self.con.execute(
            "SELECT DISTINCT cmap_name, moa, target FROM signatures "
            "WHERE cmap_name ILIKE ? ORDER BY cmap_name",
            [f"%{query}%"]
        ).fetchdf()
```

---

## Step 5: FastAPI Server

```python
# api/main.py
"""
REST API for CMap drug enrichment queries.
Deploy on HuggingFace Spaces or any container host.
"""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from typing import Optional
import os

from cmap_enrichment.engine import CMapEngine

app = FastAPI(
    title="CMap Drug Enrichment API",
    description=(
        "Query the Connectivity Map (CMap/L1000) LINCS 2020 dataset for drugs that "
        "reverse a gene expression signature. Community-hosted preservation of CLUE.io data."
    ),
    version="1.0.0",
)

DATA_DIR = os.environ.get("CMAP_DATA_DIR", "data/processed")
engine = CMapEngine(data_dir=DATA_DIR)


# -- Request/Response Models --

class EnrichmentRequest(BaseModel):
    genes_up: list[str] = Field(
        ...,
        description="Upregulated genes (HGNC symbols) in your condition",
        examples=[["APOE", "CLU", "TREM2", "C1QB", "CD68"]],
    )
    genes_down: list[str] = Field(
        default=[],
        description="Downregulated genes (optional, improves specificity)",
        examples=[["SYN1", "SNAP25", "SLC17A7"]],
    )
    cell_line: Optional[str] = Field(
        default=None,
        description="Filter to specific cell line (e.g., 'A549', 'MCF7', 'HT29')",
    )
    top_n: int = Field(
        default=50,
        ge=1,
        le=500,
        description="Number of top results to return",
    )

class EnrichmentResult(BaseModel):
    cmap_name: str
    cell_iname: str
    wtcs: float
    moa: Optional[str]
    target: Optional[str]
    n_reps: int

class EnrichmentResponse(BaseModel):
    results: list[EnrichmentResult]
    query_stats: dict


# -- Endpoints --

@app.post("/enrichment", response_model=EnrichmentResponse)
async def query_enrichment(request: EnrichmentRequest):
    """
    Given up/down-regulated genes, find drugs that REVERSE the signature.
    
    Most negative WTCS = strongest reversal = top therapeutic candidates.
    """
    try:
        df = engine.query_enrichment(
            genes_up=request.genes_up,
            genes_down=request.genes_down,
            cell_line=request.cell_line,
            top_n=request.top_n,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
    return EnrichmentResponse(
        results=df.to_dict(orient="records"),
        query_stats={
            "n_up_mapped": df.attrs.get("n_up_mapped", 0),
            "n_down_mapped": df.attrs.get("n_down_mapped", 0),
            "unmapped_genes": df.attrs.get("unmapped_genes", []),
        },
    )

@app.get("/compound/{compound_name}")
async def get_compound(compound_name: str):
    """Get all signatures and metadata for a compound."""
    df = engine.get_compound_info(compound_name)
    if df.empty:
        raise HTTPException(status_code=404, detail=f"Compound '{compound_name}' not found")
    return df.to_dict(orient="records")

@app.get("/signature/{sig_idx}")
async def get_signature(sig_idx: int):
    """Get the z-score vector for a specific signature."""
    try:
        series = engine.get_signature(sig_idx)
    except IndexError:
        raise HTTPException(status_code=404, detail=f"Signature index {sig_idx} out of range")
    return {"sig_idx": sig_idx, "zscores": series.to_dict()}

@app.get("/search")
async def search_compounds(q: str):
    """Fuzzy search for compounds by name."""
    return engine.search_compounds(q).to_dict(orient="records")

@app.get("/cell_lines")
async def list_cell_lines():
    """List all available cell lines."""
    return engine.list_cell_lines()

@app.get("/moas")
async def list_moas():
    """List all mechanisms of action."""
    return engine.list_moas()

@app.get("/health")
async def health():
    return {
        "status": "ok",
        "n_signatures": engine.n_sigs,
        "n_genes": engine.n_genes,
    }
```

---

## Step 6: MCP Server

```python
# mcp_server/server.py
"""
MCP server for CMap drug enrichment.
Runs locally -- reads DuckDB + numpy from local data directory.

Install: pip install clue-mcp
Run:     python -m cmap_enrichment.mcp_server
"""

import json
import sys
from typing import Any
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

from cmap_enrichment.engine import CMapEngine

server = Server("clue-mcp")
engine = None


def get_engine():
    global engine
    if engine is None:
        engine = CMapEngine()
    return engine


@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="cmap_drug_enrichment",
            description=(
                "Find drugs that REVERSE a gene expression signature using the "
                "Connectivity Map (CMap/L1000) LINCS 2020 database (~720K compound "
                "signatures, ~39K compounds). Provide upregulated and optionally "
                "downregulated genes from your condition. Returns compounds ranked by "
                "connectivity score (most negative = strongest reversal = top "
                "therapeutic candidate)."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "genes_up": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Upregulated gene symbols (HGNC)",
                    },
                    "genes_down": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Downregulated gene symbols (optional)",
                        "default": [],
                    },
                    "cell_line": {
                        "type": "string",
                        "description": "Filter by cell line (e.g. A549, MCF7)",
                    },
                    "top_n": {
                        "type": "integer",
                        "description": "Number of top results",
                        "default": 25,
                    },
                },
                "required": ["genes_up"],
            },
        ),
        Tool(
            name="cmap_compound_info",
            description="Get CMap signature metadata for a compound (MOA, targets, cell lines tested).",
            inputSchema={
                "type": "object",
                "properties": {
                    "compound_name": {
                        "type": "string",
                        "description": "Compound name (e.g. 'sirolimus', 'vorinostat')",
                    },
                },
                "required": ["compound_name"],
            },
        ),
        Tool(
            name="cmap_search_compounds",
            description="Search CMap for compounds by partial name match.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Partial compound name to search",
                    },
                },
                "required": ["query"],
            },
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[TextContent]:
    eng = get_engine()
    
    if name == "cmap_drug_enrichment":
        try:
            df = eng.query_enrichment(
                genes_up=arguments["genes_up"],
                genes_down=arguments.get("genes_down", []),
                cell_line=arguments.get("cell_line"),
                top_n=arguments.get("top_n", 25),
            )
            
            # Format results for LLM consumption
            lines = [
                f"CMap Drug Enrichment Results (LINCS 2020)",
                f"Query: {len(arguments['genes_up'])} up genes, "
                f"{len(arguments.get('genes_down', []))} down genes",
                f"Mapped: {df.attrs.get('n_up_mapped', '?')} up, "
                f"{df.attrs.get('n_down_mapped', '?')} down",
                "",
            ]
            
            unmapped = df.attrs.get("unmapped_genes", [])
            if unmapped:
                lines.append(f"Unmapped (not in L1000 landmarks): {', '.join(unmapped)}")
                lines.append("")
            
            lines.append(f"Top {len(df)} reversals (negative WTCS = reversal):")
            lines.append("-" * 70)
            
            for _, row in df.iterrows():
                lines.append(
                    f"  {row['cmap_name']:<25} | WTCS: {row['wtcs']:+.4f} | "
                    f"Cell: {row['cell_iname']:<6} | MOA: {row.get('moa', 'N/A')}"
                )
            
            return [TextContent(type="text", text="\n".join(lines))]
        
        except ValueError as e:
            return [TextContent(type="text", text=f"Error: {str(e)}")]
    
    elif name == "cmap_compound_info":
        df = eng.get_compound_info(arguments["compound_name"])
        if df.empty:
            return [TextContent(type="text", text=f"Compound '{arguments['compound_name']}' not found")]
        return [TextContent(type="text", text=df.to_string())]
    
    elif name == "cmap_search_compounds":
        df = eng.search_compounds(arguments["query"])
        if df.empty:
            return [TextContent(type="text", text=f"No compounds matching '{arguments['query']}'")]
        return [TextContent(type="text", text=df.to_string())]
    
    return [TextContent(type="text", text=f"Unknown tool: {name}")]


async def main():
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())
```

**Claude Code config** (`~/.claude/claude_desktop_config.json` or CLAUDE.md):

```json
{
  "mcpServers": {
    "cmap": {
      "command": "python",
      "args": ["-m", "cmap_enrichment.mcp_server"],
      "env": {
        "CMAP_DATA_DIR": "/path/to/data/processed"
      }
    }
  }
}
```

---

## Step 7: Deploy

### 7a. GitHub Repository Structure

```
clue-mcp/
├── README.md
├── LICENSE                    # MIT
├── pyproject.toml
├── Dockerfile
├── scripts/
│   ├── 01_parse_gctx.py
│   ├── 02_aggregate_signatures.py
│   └── 03_build_duckdb.py
├── cmap_enrichment/
│   ├── __init__.py
│   ├── engine.py
│   └── mcp_server.py
├── api/
│   ├── main.py
│   └── requirements.txt
└── data/
    ├── raw/                   # Symlink or mount to /dfs7/swaruplab/shared_lab/Clue_database/
    └── processed/             # .gitignore'd; downloaded from Zenodo
        ├── cmap.duckdb
        ├── rank_matrix.npy
        ├── zscore_matrix.npy
        ├── signature_metadata.parquet
        └── gene_names.txt
```

### 7b. pyproject.toml

```toml
[project]
name = "clue-mcp"
version = "1.0.0"
description = "CMap/L1000 LINCS 2020 drug enrichment engine -- community preservation of CLUE.io"
readme = "README.md"
license = {text = "MIT"}
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.24",
    "pandas>=2.0",
    "scipy>=1.10",
    "duckdb>=0.9",
    "mcp>=1.0",
]

[project.optional-dependencies]
api = ["fastapi>=0.100", "uvicorn>=0.20"]
build = ["cmapPy>=4.0", "pyarrow>=12.0"]

[project.scripts]
cmap-mcp = "cmap_enrichment.mcp_server:main"
```

### 7c. Dockerfile

```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY pyproject.toml .
COPY cmap_enrichment/ cmap_enrichment/
COPY api/ api/

RUN pip install ".[api]"

# Data mounted or downloaded at runtime
ENV CMAP_DATA_DIR=/app/data/processed

EXPOSE 8000
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### 7d. HuggingFace Space (app.py)

```python
# For HuggingFace Spaces deployment
# Space type: Docker or Gradio

import gradio as gr
from cmap_enrichment.engine import CMapEngine

engine = CMapEngine(data_dir="data/processed")

def run_enrichment(genes_up_text, genes_down_text, cell_line, top_n):
    genes_up = [g.strip() for g in genes_up_text.split(",") if g.strip()]
    genes_down = [g.strip() for g in genes_down_text.split(",") if g.strip()]
    
    if not genes_up:
        return "Please enter at least one upregulated gene."
    
    df = engine.query_enrichment(
        genes_up=genes_up,
        genes_down=genes_down,
        cell_line=cell_line if cell_line else None,
        top_n=int(top_n),
    )
    
    # Format for display
    display_cols = ["cmap_name", "wtcs", "cell_iname", "moa", "target", "n_reps"]
    return df[display_cols].to_markdown(index=False)

demo = gr.Interface(
    fn=run_enrichment,
    inputs=[
        gr.Textbox(
            label="Upregulated Genes (comma-separated)",
            placeholder="APOE, CLU, TREM2, C1QB, CD68, TYROBP",
        ),
        gr.Textbox(
            label="Downregulated Genes (optional, comma-separated)",
            placeholder="SYN1, SNAP25, SLC17A7, CAMK2A",
        ),
        gr.Dropdown(
            choices=[""] + engine.list_cell_lines(),
            label="Cell Line Filter (optional)",
        ),
        gr.Slider(minimum=10, maximum=200, value=50, step=10, label="Top N Results"),
    ],
    outputs=gr.Markdown(label="Drug Enrichment Results"),
    title="CMap Drug Enrichment (LINCS 2020)",
    description=(
        "Community-hosted CMap/L1000 LINCS 2020 drug enrichment engine. "
        "Enter upregulated (and optionally downregulated) genes from your condition. "
        "Returns drugs ranked by connectivity score -- most negative = strongest reversal."
    ),
)

demo.launch()
```

### 7e. Zenodo Upload Checklist

> **See the dedicated, up-to-date guide: [Releasing the database](data-release.md).** It
> pins the exact archive name (`clue_mcp_processed.tar.gz`) and layout that
> `scripts/download_data.sh` expects, plus the post-publish steps to wire the DOI back
> into the repo. The sketch below is kept for historical context.

```bash
# Files to deposit on Zenodo:
data/processed/
├── cmap.duckdb                      # ~100-200 MB
├── rank_matrix.npy                  # ~2-4 GB (larger with LINCS 2020)
├── zscore_matrix.npy                # ~4-8 GB
├── signature_metadata.parquet       # ~20 MB
└── gene_names.txt                   # <1 MB

# Zenodo metadata:
# Title: "CMap/L1000 LINCS 2020 Level 5 Drug Signatures -- Processed for Enrichment Analysis"
# Description: "Pre-processed Connectivity Map drug perturbation signatures
#   from the LINCS 2020 data release (Dec 2020). ~720K compound treatment
#   signatures across 12,328 genes and ~39K compounds. Aggregated consensus
#   signatures with pre-computed gene rankings for fast enrichment queries.
#   Community preservation of CLUE.io functionality."
# License: CC-BY-4.0
# Related identifiers:
#   - LINCS 2020 S3: s3://macchiato.clue.io/builds/LINCS2020/
#   - doi:10.1016/j.cell.2017.10.049 (Subramanian et al. 2017)
#   - GitHub repo URL
```

---

## Quick Start for Other Labs

```bash
# Option 1: Use the hosted API (no install)
curl -X POST https://huggingface.co/spaces/YOUR_SPACE/enrichment \
  -H "Content-Type: application/json" \
  -d '{"genes_up": ["APOE","CLU","TREM2"], "genes_down": ["SYN1","SNAP25"]}'

# Option 2: Run locally with Docker
git clone https://github.com/swaruplab/clue-mcp
cd clue-mcp
# Download data from Zenodo (DOI link in README)
wget -P data/processed/ https://zenodo.org/records/XXXXX/files/cmap_processed.tar.gz
tar xzf data/processed/cmap_processed.tar.gz -C data/processed/
docker build -t cmap . && docker run -p 8000:8000 -v ./data:/app/data cmap

# Option 3: Python package + MCP server
pip install clue-mcp
# Download data...
cmap-mcp  # starts MCP server for Claude Code / Operon
```

---

## Memory / Performance Estimates

| Component | Disk | RAM at Runtime |
|-----------|------|----------------|
| rank_matrix.npy (~720K sigs x 12328 genes, int16) | ~17 GB | ~17 GB |
| zscore_matrix.npy (float32) | ~34 GB | ~34 GB (lazy-loadable) |
| cmap.duckdb | ~200 MB | ~100 MB |
| **API server (ranks only, exemplars)** | **~4-8 GB** | **~8-16 GB** |

**Optimization**: With LINCS 2020's larger dataset, memory management is critical:

1. **Use exemplar signatures only** (`is_exemplar_sig == 1`) to reduce from 720K to a more manageable subset
2. **Memory-mapped numpy** for constrained environments:
   ```python
   self.rank_matrix = np.load(path / "rank_matrix.npy", mmap_mode="r")
   ```
3. **For HuggingFace free tier** (16 GB RAM): load only ranks, use mmap, and consider further filtering by QC metrics (`is_hiq == 1`, `qc_pass == 1`)

---

## Key Differences from GSE92742 (Old Guide)

| Aspect | Old (GSE92742) | Current (LINCS 2020) |
|--------|---------------|----------------------|
| Source | GEO (NCBI FTP) | S3 `macchiato.clue.io` |
| Total signatures | ~473K + 118K | ~1.2M |
| Drug (trt_cp) signatures | ~200K | **720,216** |
| Files | Single monolithic GCTX | Split by pert_type |
| Drug name column | `pert_iname` | `cmap_name` |
| Cell line column | `cell_id` | `cell_iname` |
| Gene ID column | `pr_gene_id` | `gene_id` |
| Gene symbol column | `pr_gene_symbol` | `gene_symbol` |
| Landmark flag | `pr_is_lm` | `feature_space == "lm"` |
| Metadata files | `sig_info.txt`, `pert_info.txt`, `gene_info.txt` | `siginfo_beta.txt`, `compoundinfo_beta.txt`, `geneinfo_beta.txt`, `cellinfo_beta.txt` |
| Exemplar flag | N/A | `is_exemplar_sig` (filter to 1 for best replicate) |
| QC metrics | Limited | `tas`, `cc_q75`, `is_hiq`, `qc_pass` |
| Compounds | ~20K | **~39K** |

---

## Additional Data Available (Not Yet Downloaded)

The LINCS 2020 S3 bucket also provides:

```
# Level 4 (replicate-normalized profiles -- individual replicates before collapsing)
level4/level4_beta_trt_cp_n1805898x12328.gctx   # 1.8M individual profiles

# Level 3 (landmark-only, replicate profiles)
level3/level3_beta_trt_cp_n1805898x12328.gctx

# Instance-level metadata (for Level 3/4)
instinfo_beta.txt

# Google BigQuery access for arbitrary subsets:
# https://console.cloud.google.com/bigquery?p=cmap-big-table&d=cmap_lincs_public_views&page=dataset
```
