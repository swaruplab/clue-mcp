# Python Library Tutorial

A walkthrough of using the `cmap_enrichment` library from Python — querying the LINCS 2020 database for drugs that reverse a gene-expression signature.

> **Prerequisite:** the processed data must exist under `data/processed/` — shared `gene_names.txt` + `cmap.duckdb`, plus a per-class subfolder (`drug/`, `knockdown/`, `overexpression/`) each holding `rank_matrix.npy`, `zscore_matrix.npy`, and `metadata.parquet`. Pull it from Zenodo with `scripts/download_data.sh`, or build it with the [pipeline](pipeline.md). Only `drug/` is needed to follow this page.

---

## 1. Install

```bash
git clone https://github.com/swaruplab/clue-mcp.git
cd clue-mcp
pip install -e ".[viz]"      # core + matplotlib + networkx
```

## 2. Load the engine

```python
from cmap_enrichment import CMapEngine

# Reads data/processed/ from CWD by default — pass data_dir to override.
# pert_class defaults to "drug"; pass "knockdown" or "overexpression" for the others.
engine = CMapEngine(data_dir="data/processed", pert_class="drug")

print(f"Loaded {engine.n_sigs:,} {engine.pert_class} signatures over {engine.n_genes:,} genes")
# Loaded 3,493 drug signatures over 12,328 genes
```

Loading takes ~5–10 seconds and memory-maps the rank + z-score matrices (they are
never read fully into RAM), so the resident footprint stays small.

### Querying across all three perturbation classes

`EngineRegistry` lazily manages one engine per class and reports which datasets are
installed — the same component the MCP server and REST API use:

```python
from cmap_enrichment import EngineRegistry

reg = EngineRegistry(data_dir="data/processed")
print(reg.available())   # {'drug': True, 'knockdown': False, 'overexpression': False}

drug_eng = reg.get("drug")                 # or reg.get("cmap_drug_enrichment")
res = drug_eng.rank_perturbations(genes_up=["APOE", "CLU"], top_n=10)
res["reversing"]   # candidate therapeutics (most-negative WTCS)
res["mimicking"]   # signature-inducers   (most-positive WTCS)
res["stats"]       # mapping diagnostics
```

## 3. Run an enrichment query

```python
# Disease signature: microglial activation in Alzheimer's
genes_up   = ["APOE", "CLU", "TREM2", "CD68", "C1QB", "C1QA", "CTSB"]
genes_down = ["SYN1", "SNAP25", "SLC17A7", "GAP43"]

hits = engine.query_enrichment(
    genes_up=genes_up,
    genes_down=genes_down,
    top_n=20,
)

# Sorted ascending by WTCS — most negative = strongest reversal
print(hits[["cmap_name", "cell_iname", "wtcs", "moa", "target"]].head(10))
```

### Filtering to a specific cell line

```python
hits_a549 = engine.query_enrichment(
    genes_up=genes_up,
    cell_line="A549",      # restricts the search space
    top_n=10,
)
```

### Single-set queries

Pass only `genes_up` (no `genes_down`) and the engine drops back to a one-sided KS enrichment.

## 4. Inspect mapping diagnostics

The result DataFrame carries metadata about your query in `.attrs`:

```python
print("Up genes mapped to L1000:",   hits.attrs["n_up_mapped"])
print("Down genes mapped to L1000:", hits.attrs["n_down_mapped"])
print("Unmapped:",                   hits.attrs["unmapped_genes"])
```

L1000 only profiles 12,328 genes (978 landmark + 11,350 inferred), so some query genes will always be missing.

## 5. Explore the database

```python
# What cell lines are available?
engine.list_cell_lines()
# ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2', ...]  (41 total)

# What mechanisms of action have signatures?
engine.list_moas()[:5]
# ['Acetylcholine receptor antagonist', 'Adrenergic receptor agonist', ...]

# Fuzzy compound search
engine.search_compounds("statin")
#   cmap_name      moa                            target
# 0 atorvastatin   HMG-CoA reductase inhibitor    HMGCR
# 1 simvastatin    HMG-CoA reductase inhibitor    HMGCR
# ...

# Get every signature for a single compound
engine.get_compound_info("sirolimus")

# Get the full z-score vector for a single signature
zscores = engine.get_signature(42)   # pandas.Series indexed by gene symbol
```

## 6. Permutation-based significance

The library returns raw WTCS scores. For FDR-corrected drug rankings, use the higher-level analysis driver:

```bash
python analysis/run_enrichment.py path/to/your_genes.csv
```

It runs 10,000 permutations of random gene sets to compute empirical p-values and Benjamini-Hochberg FDR, plus generates 11 publication-quality figures. See [Running an analysis](usage.md) for full output details.

---

## Next steps

- [REST API tutorial](api-tutorial.md) — same engine, exposed over HTTP
- [MCP server setup](mcp-setup.md) — make the engine callable from Claude
- [Hallmark cross-pathway showcase](showcase.md) — the engine applied to all 50 MSigDB Hallmark gene sets
