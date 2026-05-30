# Examples

Minimal end-to-end examples for each of the three interfaces.

| File | What it shows |
|------|---------------|
| [`quickstart_library.py`](quickstart_library.py) | Querying all three perturbation classes (drug / knockdown / overexpression) from Python via `EngineRegistry` |
| [`quickstart_api.sh`](quickstart_api.sh) | `curl` against the FastAPI server |
| [`sample_genes_ad_microglia.csv`](sample_genes_ad_microglia.csv) | Two-column (gene, log2FC) signature — drop into `analysis/run_enrichment.py` |
| [`sample_genes_mtor_pathway.txt`](sample_genes_mtor_pathway.txt) | Single-column gene list — same idea, no fold-change values |
| [`hallmark_signatures_walkthrough.md`](hallmark_signatures_walkthrough.md) | Step-by-step recap of the bundled cross-MSigDB-Hallmark showcase |

## Prerequisites

The processed CMap database must exist at `../data/processed/`. Pull it from Zenodo with `scripts/download_data.sh`, or rebuild it with the [3-step pipeline](../README.md#3-run-the-processing-pipeline).

```bash
# Verify — shared files plus one subfolder per installed perturbation class
ls ../data/processed/
# Expect: gene_names.txt  cmap.duckdb  drug/  [knockdown/]  [overexpression/]
ls ../data/processed/drug/
# Expect: rank_matrix.npy  zscore_matrix.npy  metadata.parquet
```

Only the `drug/` class is required for the examples to run; `knockdown/` and
`overexpression/` are optional and the scripts skip any class that isn't installed.

## Run them

```bash
# 1. Python library
python examples/quickstart_library.py

# 2. REST API (in two terminals)
#    Terminal A:
uvicorn api.main:app --port 8000
#    Terminal B:
bash examples/quickstart_api.sh

# 3. Run the bundled analysis driver on a sample signature
python analysis/run_enrichment.py examples/sample_genes_ad_microglia.csv
```
