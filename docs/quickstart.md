# Quickstart

This page takes you from a fresh clone to your first drug-enrichment result. Budget about 10 minutes (plus the one-time database build, which is separate).

---

## 1. Install

```bash
git clone https://github.com/swaruplab/clue-mcp.git
cd clue-mcp

# Core engine only
pip install -e .

# …or pick the extras you need
pip install -e ".[viz]"     # + matplotlib, networkx  (analysis figures)
pip install -e ".[api]"     # + fastapi, uvicorn       (REST API)
pip install -e ".[mcp]"     # + mcp SDK                (Claude integration)
pip install -e ".[all]"     # everything
```

`clue-mcp` requires Python ≥ 3.10. The importable package is `cmap_enrichment`.

---

## 2. Get the database

The engine queries a processed copy of the **LINCS 2020** Connectivity Map, shipped on Zenodo as **per-class bundles** (drug ≈ 9 GB, plus optional knockdown and overexpression). It's *not* bundled with the repo (too large for git, and it's a derived product of the Broad's data — see [Method & data](method.md)). The repo ships only **code** plus the small **precomputed [showcase](showcase.md)** outputs.

Pick whichever path fits you:

=== "Download the prebuilt DB (recommended)"

    Download the already-processed database from Zenodo — no raw files, no HPC
    node, no rebuild. The **drug** bundle alone is enough to get started:

    ```bash
    bash scripts/download_data.sh            # base + every available class
    bash scripts/download_data.sh drug       # …or just one class
    # or, if the record id isn't wired into your checkout yet:
    ZENODO_RECORD=<record-id> bash scripts/download_data.sh
    ```

    It downloads into `data/processed/` and verifies each archive's checksum.
    The Zenodo DOI is listed on the [repo README](https://github.com/swaruplab/clue-mcp#getting-the-data).

=== "Already have it"

    Anyone in a lab that already built the database (or mounted a shared copy)
    just points the engine at it:

    ```bash
    export CMAP_DATA_DIR=/path/to/data/processed
    ```

    The directory must contain `gene_names.txt`, `cmap.duckdb`, and at least the
    `drug/` subfolder (`rank_matrix.npy`, `zscore_matrix.npy`, `metadata.parquet`).

=== "Build it from scratch"

    Download the raw LINCS 2020 GCTX files (~34 GB) and run the 3-step pipeline
    (Steps 1–2 need ~256 GB RAM, so run them on an HPC node). Full instructions:

    [Building the database →](pipeline.md)

!!! tip "No data yet? You can still explore — 0 GB."
    The bundled **[50-Hallmark showcase](showcase.md)** ships pre-computed results
    (tables + figures), so you can see exactly what the engine produces before
    downloading or building anything.

---

## 3. Your first query (Python)

```python
from cmap_enrichment import CMapEngine

engine = CMapEngine(data_dir="data/processed")
print(f"{engine.n_sigs:,} signatures × {engine.n_genes:,} genes")

hits = engine.query_enrichment(
    genes_up=["APOE", "CLU", "TREM2", "CD68", "C1QB"],   # disease-up genes
    genes_down=["SYN1", "SNAP25", "SLC17A7"],            # disease-down genes
    top_n=10,
)

# WTCS sorted ascending — most negative = strongest reversal
print(hits[["cmap_name", "cell_iname", "wtcs", "moa", "target"]])
```

How to read it:

| WTCS | Meaning |
|------|---------|
| **Negative** | Drug **reverses** your signature → therapeutic candidate |
| **Positive** | Drug **mimics** your signature |
| **~ 0** | No meaningful connectivity |

Full walkthrough, including cell-line filtering and database exploration: **[Python tutorial](python-tutorial.md)**.

---

## 4. The other three interfaces

=== "MCP (Claude)"

    Register the server once, then just *ask*:

    ```bash
    pip install -e ".[mcp]"
    CMAP_DATA_DIR=data/processed python mcp_server/server.py   # smoke test (Ctrl-C to quit)
    ```

    > **You:** I have an AD microglia signature — up: APOE, CLU, TREM2; down: SYN1, SNAP25. What drugs might reverse it?
    >
    > **Claude:** *(calls `cmap_drug_enrichment`)* The top reversal candidates are…

    Full setup for Claude Desktop / Claude Code / Cursor: **[MCP server setup](mcp-setup.md)**.

=== "REST API"

    ```bash
    pip install -e ".[api]"
    uvicorn api.main:app --port 8000      # interactive docs at /docs
    ```

    ```bash
    curl -X POST http://localhost:8000/enrich/drug \
      -H "Content-Type: application/json" \
      -d '{"genes_up": ["APOE","CLU","TREM2","CD68","C1QB"],
           "genes_down": ["SYN1","SNAP25","SLC17A7"], "top_n": 10}'
    ```

    There's also `/enrich/knockdown` and `/enrich/overexpression` with the same body.

    Full reference: **[REST API tutorial](api-tutorial.md)**.

=== "Analysis driver"

    For permutation FDR, MOA enrichment, and 11 figure panels from a gene file:

    ```bash
    python analysis/run_enrichment.py examples/sample_genes_ad_microglia.csv
    ```

    What every output table and figure means: **[Running an analysis](usage.md)**.

---

## Ready-made example inputs

The [`examples/`](https://github.com/swaruplab/clue-mcp/tree/main/examples) folder has drop-in gene lists:

| File | Format |
|------|--------|
| `sample_genes_ad_microglia.csv` | 2-column (gene, log2FC) — computes a full up/down WTCS |
| `sample_genes_mtor_pathway.txt` | 1-column gene list — single-set KS enrichment |

---

## Where to next

- **[Showcase — 50 Hallmarks](showcase.md)** — a complete, pre-computed worked example.
- **[Method & data](method.md)** — the science behind the score.
- **[Building the database](pipeline.md)** — generate the engine from raw GCTX files.
