# clue-mcp

> Turn a disease gene-expression signature into ranked lists of candidate drugs, gene-knockdown targets, and overexpression hits — from Python, a REST API, or directly inside Claude via the **Model Context Protocol (MCP)**.

**clue-mcp** queries the **Connectivity Map (CMap) LINCS 2020** database from the Broad Institute. Give it a list of up- and/or down-regulated genes; it scores perturbation signatures across **12,328 genes** with the Weighted Connectivity Score (WTCS) and returns the perturbations most likely to *reverse* (or *mimic*) your signature.

The catalog is split into **three perturbation classes**, each exposed as its own framed tool/endpoint so an LLM picks the right scientific interpretation:

| Class | Perturbation type(s) | Reverses your signature ⇒ |
|-------|----------------------|---------------------------|
| **Drug** | small molecules (`trt_cp`) | candidate **therapeutic** to give |
| **Knockdown** | shRNA + CRISPR (`trt_sh`, `trt_xpr`) | candidate **driver gene** to silence |
| **Overexpression** | ORF (`trt_oe`) | gene whose **forced expression** opposes the state |

The drug class also runs permutation-based FDR via the analysis driver.

📖 **Documentation & tutorials:** **https://swaruplab.github.io/clue-mcp/**

Built on the [LINCS 2020 Level 5](https://clue.io/data/CMap2020#LINCS2020) dataset.

**Four ways to use it:**

| Interface | For | Entry point |
|-----------|-----|-------------|
| 🤖 **MCP server** | Asking Claude/LLMs "what drugs reverse this signature, and what genes should I knock down?" in plain English | [`mcp_server/`](mcp_server/) · [setup](docs/mcp-setup.md) |
| 🐍 **Python library** | Notebooks, batch jobs, custom pipelines | [`cmap_enrichment/`](cmap_enrichment/) · [tutorial](docs/python-tutorial.md) |
| 🌐 **REST API** | Web apps, other languages, hosted service | [`api/`](api/) · [tutorial](docs/api-tutorial.md) |
| 📊 **Analysis driver** | Publication figures + permutation FDR | [`analysis/run_enrichment.py`](analysis/run_enrichment.py) |

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Setup](#setup)
  - [Install](#1-install)
  - [Download LINCS 2020 Data](#2-download-lincs-2020-data)
  - [Run the Processing Pipeline](#3-run-the-processing-pipeline)
- [Quickstart](#quickstart)
- [Running an Enrichment Analysis](#running-an-enrichment-analysis)
  - [Input Formats](#input-gene-file-formats)
  - [Run](#run-the-analysis)
  - [Output](#output)
  - [Parameters](#configurable-parameters)
- [Python Library](#using-the-enrichment-engine-as-a-python-library)
- [REST API](#rest-api)
- [MCP Server](#mcp-server-claude-integration)
- [Worked Example: 50 MSigDB Hallmarks](#worked-example-50-msigdb-hallmarks)
- [Documentation](#documentation)
- [Method](#method)
- [Data Source](#data-source)

---

## Overview

```
                        ┌─────────────────────────┐
                        │   Your Gene Signature    │
                        │  (up/down-regulated genes)│
                        └────────────┬────────────┘
                                     │
                                     v
┌────────────────────────────────────────────────────────────────────┐
│                    CMap Enrichment Engine                          │
│                                                                    │
│   rank_matrix.npy (3.1 GB)  ─┐                                    │
│   zscore_matrix.npy (6.1 GB) ├──> WTCS scoring ──> FDR correction │
│   cmap.duckdb (9 MB)        ─┘    (KS-like)       (permutation)  │
│                                                                    │
│   3,493 consensus signatures | 511 compounds | 41 cell lines      │
└────────────────────────────────────┬───────────────────────────────┘
                                     │
                                     v
                        ┌─────────────────────────┐
                        │    Ranked Drug List      │
                        │  + FDR, MOA, targets,    │
                        │  heatmaps, networks,     │
                        │  volcano plots           │
                        └─────────────────────────┘
```

---

## Repository Structure

```
.
├── scripts/                        # Data processing pipeline (per perturbation class)
│   ├── build_config.py             #   pert_type -> GCTX file / class / method map
│   ├── 01_parse_gctx.py            #   Parse GCTX -> Parquet (exemplar signatures)
│   ├── 02_aggregate_signatures.py  #   Aggregate per (perturbagen x cell line x method), pre-rank
│   ├── 03_build_duckdb.py          #   Index all classes' metadata into one DuckDB
│   ├── slurm_step1.sh              #   SLURM Step 1 (forwards CMAP_PERT_CLASS)
│   ├── slurm_step2.sh              #   SLURM Step 2 (forwards CMAP_PERT_CLASS)
│   ├── slurm_step3.sh              #   SLURM Step 3 (32 GB RAM)
│   ├── build_all.sh                #   Submit step1->step2 per class + one step3 (deps)
│   ├── download_data.sh            #   Download + verify per-class bundles from Zenodo
│   └── make_zenodo_archives.sh     #   Build the Zenodo tarballs from data/processed
│
├── cmap_enrichment/                # Core enrichment engine (importable library)
│   ├── __init__.py
│   ├── perturbation_classes.py     #   Single source of truth: the 3 classes + framing
│   ├── engine.py                   #   CMapEngine class (WTCS algorithm, per class)
│   ├── registry.py                 #   EngineRegistry: lazy per-class engines
│   └── zenodo_manifest.json        #   DOI + per-class archive checksums
│
├── analysis/                       # Enrichment analyses on user gene lists
│   ├── run_enrichment.py           #   Drug analysis driver (permutation FDR + figures)
│   ├── slurm_run.sh                #   SLURM job wrapper (64 GB RAM)
│   └── hallmark_enrichment/        #   Worked example: 50 MSigDB Hallmarks
│
├── api/                            # REST API server (FastAPI) — 3 /enrich/* endpoints
│   └── main.py
│
├── mcp_server/                     # MCP server for Claude/LLM integration (3 tools)
│   ├── __init__.py
│   └── server.py
│
├── tests/                          # pytest suite (synthetic bundle, no big download)
│   ├── conftest.py
│   ├── test_perturbation_classes.py
│   ├── test_engine.py
│   └── test_mcp_tools.py
│
├── examples/                       # Quickstart scripts + sample gene lists
│   ├── quickstart_library.py
│   ├── quickstart_api.sh
│   └── sample_genes_*.{csv,txt}
│
├── docs/                           # Long-form tutorials + website source (MkDocs)
│   ├── index.md  quickstart.md  usage.md
│   ├── python-tutorial.md  api-tutorial.md  mcp-setup.md
│   ├── showcase.md  pipeline.md  method.md
│   └── claude_desktop_config.example.json
│
├── data/                           # Generated/downloaded (not tracked in git)
│   ├── intermediate/               #   Step 1 outputs (per class)
│   └── processed/                  #   gene_names.txt, cmap.duckdb, + drug/ knockdown/ overexpression/
│
├── mkdocs.yml                      # Website configuration
├── pyproject.toml                  # pip-installable package config (+ clue-mcp-server)
├── requirements.txt                # legacy pinned dependencies
├── LICENSE                         # MIT
└── README.md
```

---

## Setup

### 1. Install

Clone and install as an editable package — this exposes the `cmap_enrichment` library on your `PYTHONPATH` and lets you opt into the API/MCP/viz extras.

```bash
git clone https://github.com/swaruplab/clue-mcp.git
cd clue-mcp

# Core engine only
pip install -e .

# Or pick the extras you need:
pip install -e ".[viz]"          # + matplotlib, networkx (for analysis/ scripts)
pip install -e ".[api]"          # + fastapi, uvicorn
pip install -e ".[mcp]"          # + mcp SDK
pip install -e ".[pipeline]"     # + cmapPy (only needed to (re)build the database)
pip install -e ".[all]"          # everything
```

The legacy `requirements.txt` is still provided for `pip install -r` installs.

**Key packages:** `numpy`, `pandas`, `scipy`, `duckdb`, `cmapPy`, `matplotlib`, `networkx`, `fastapi`, `mcp`

### 2. Getting the data

The engine queries a **processed database** under `data/processed/`. **This repository does not ship that database** — it's too large for git and is a derived product of the Broad's LINCS data (see [data policy](https://clue.io/connectopedia/data_use_policy)). The repo contains only **code** plus the small **precomputed [showcase](analysis/hallmark_enrichment/)** outputs.

The data ships on Zenodo as **per-class bundles** plus a small shared base (`gene_names.txt`, `cmap.duckdb`). You can download only the classes you need — the **drug** bundle is enough to get started, and the knockdown/overexpression tools simply report "not available" until their bundles are present.

Choose the path that fits you:

| Tier | You get | Cost | How |
|------|---------|------|-----|
| **Look** | Browse precomputed results; nothing to install | **0 GB** | Open [`analysis/hallmark_enrichment/`](analysis/hallmark_enrichment/) or the [website](https://swaruplab.github.io/clue-mcp/showcase/) |
| **Run** *(recommended)* | The queryable database | **~9 GB (drug) +** | `bash scripts/download_data.sh` (prebuilt bundles from Zenodo) |
| **Rebuild** | Regenerate the DB from raw files | 34 GB + HPC | Section 3 below |

```bash
# Download base + all available per-class bundles (no raw files, no HPC needed)
bash scripts/download_data.sh
#   …or a specific class:   bash scripts/download_data.sh drug
#   …or pin a record:       ZENODO_RECORD=<record-id> bash scripts/download_data.sh
```

> **Prebuilt database DOI:** [`10.5281/zenodo.20465969`](https://doi.org/10.5281/zenodo.20465969). `scripts/download_data.sh` resolves this automatically.

> **Already have a copy?** Point the engine at it and skip the rest of setup:
> ```bash
> export CMAP_DATA_DIR=/path/to/data/processed
> ```

---

### 3. (Advanced) Build the database from raw data

*Only needed if you're regenerating the database instead of downloading it.*

Download the GCTX file(s) for the class(es) you want, plus the shared metadata, from [LINCS 2020](https://clue.io/data/CMap2020#LINCS2020) and place them in the project root:

| File | Size | Class |
|------|------|-------|
| `level5_beta_trt_cp_n720216x12328.gctx` | 34 GB | drug |
| `level5_beta_trt_sh_n238351x12328.gctx` | 11 GB | knockdown (shRNA) |
| `level5_beta_trt_xpr_n142901x12328.gctx` | 6 GB | knockdown (CRISPR) |
| `level5_beta_trt_oe_n34171x12328.gctx` | 1.6 GB | overexpression |
| `siginfo_beta.txt` | 444 MB | shared — signature metadata |
| `geneinfo_beta.txt` | 1.1 MB | shared — gene metadata |
| `compoundinfo_beta.txt` | 4.5 MB | shared — MoA/target (drug only) |
| `cellinfo_beta.txt` | 38 KB | shared — cell line metadata |

The download URLs are listed in `links2020_level5_urls.txt`.

Then run the pipeline. Steps 1–2 run **per class** (selected via the `CMAP_PERT_CLASS` env var, default `drug`); Step 3 indexes whichever classes are present into one DuckDB.

> **Note:** Steps 1 and 2 require **~256 GB RAM** and should be run on an HPC cluster (the knockdown class is the heaviest — use the `maxmem` partition).

```bash
# Build one class at a time …
CMAP_PERT_CLASS=drug sbatch scripts/slurm_step1.sh
CMAP_PERT_CLASS=drug sbatch scripts/slurm_step2.sh
sbatch scripts/slurm_step3.sh

# … or submit every class with the right SLURM dependencies in one shot:
bash scripts/build_all.sh
```

**Pipeline output** (saved to `data/processed/`):

| Path | Description |
|------|-------------|
| `gene_names.txt` | 12,328 gene symbols (shared row index) |
| `cmap.duckdb` | Indexed DuckDB over all built classes |
| `<class>/rank_matrix.npy` | Pre-ranked gene lists per signature (int16) |
| `<class>/zscore_matrix.npy` | Z-score matrix, n_sigs × 12,328 (float32) |
| `<class>/metadata.parquet` | Perturbagen, cell line, pert_type, method, MoA/target |

where `<class>` ∈ {`drug`, `knockdown`, `overexpression`}.

**Drug database:** 3,493 consensus signatures from 511 compounds across 41 cell lines.

---

## Quickstart

Once the pipeline has produced `data/processed/`, the fastest way to see something working is:

```bash
# 1. Python library: one query against an AD-microglia signature
python examples/quickstart_library.py

# 2. REST API (two terminals)
uvicorn api.main:app --port 8000          # terminal A
bash examples/quickstart_api.sh           # terminal B

# 3. Full analysis with permutation FDR + 11 figure panels
python analysis/run_enrichment.py examples/sample_genes_ad_microglia.csv
```

More sample inputs and a step-by-step walkthrough are in [`examples/`](examples/).

---

## Running an Enrichment Analysis

### Input Gene File Formats

The analysis script (`analysis/run_enrichment.py`) auto-detects the input format.

**Format A: Gene symbols only** (1-column, treated as upregulated)

```
# One gene per line (.txt, .csv, .tsv)
APOE
CLU
TREM2
CD68
C1QB
```

**Format B: Gene symbols with values** (2-column, log2FC / kME / etc.)

```
gene_symbol,log2FC
APOE,2.31
CLU,1.87
SYN1,-1.54
SNAP25,-2.01
```

- Positive values = upregulated, negative = downregulated
- With 2-column input, the script computes **WTCS** (up vs down gene sets)
- With 1-column input, it uses a single-set **KS enrichment** score
- Supports `.txt`, `.csv`, `.tsv` with auto-detected delimiters
- Optional header row is auto-detected

### Run the Analysis

```bash
# Direct execution (requires ~10 GB RAM to load matrices)
python analysis/run_enrichment.py /path/to/your_genes.csv

# Via SLURM on HPC (recommended)
GENE_FILE=/path/to/your_genes.csv sbatch analysis/slurm_run.sh
```

Output files are written to the same directory as the script.

### Output

The script produces publication-quality **CSV tables** and **figures**.

#### CSV Tables

| File | Description |
|------|-------------|
| `drug_level_summary.csv` | All compounds ranked by median enrichment score, with permutation p-values and FDR |
| `top50_reversing_signatures.csv` | Top 50 individual signatures with strongest reversal |
| `top50_mimicking_signatures.csv` | Top 50 individual signatures with strongest mimicry |
| `all_signature_scores.csv` | Complete signature-level enrichment scores |
| `moa_enrichment.csv` | Mechanism of action enrichment (hypergeometric test + FDR) |
| `cell_line_stats.csv` | Per-cell-line score statistics |
| `summary.txt` | Human-readable results summary |

#### Figures (11 panels, all with dual views: all drugs vs named drugs)

| Figure | Description |
|--------|-------------|
| `plot1_score_distribution.png` | Enrichment score distribution with 5th/95th percentile cutoffs |
| `plot2_top25_reversing_signature_level.png` | Top 25 reversing signatures (bar chart) |
| `plot3_top25_mimicking_signature_level.png` | Top 25 mimicking signatures (bar chart) |
| `plot4_top25_reversing_drug_level.png` | Top 25 reversing drugs, median across cell lines, with FDR annotations |
| `plot5_top25_mimicking_drug_level.png` | Top 25 mimicking drugs, median across cell lines, with FDR annotations |
| `plot6_moa_enrichment.png` | MOA fold-enrichment among top 5% reversing drugs |
| `plot7_cellline_comparison.png` | Median enrichment score by cell line |
| `plot8a/b_heatmap_*.png` | Drug x cell line enrichment score heatmaps |
| `plot9a/b_zscore_heatmap_*.png` | Z-score heatmaps of query genes in top reversing signatures |
| `plot10_volcano.png` | Volcano plot (median score vs -log10 FDR) |
| `plot11a/b_network_*.png` | Drug-gene bipartite network (edges where \|z-score\| > threshold) |

> Each multi-panel figure shows **(A) all drugs** and **(B) named drugs only** (excluding BRD-coded compound identifiers).

### Configurable Parameters

Set at the top of `analysis/run_enrichment.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_PERM` | 10,000 | Number of permutations for FDR calculation |
| `TOP_N` | 50 | Number of top results to save per direction |
| `ZSCORE_THRESHOLD` | 1.5 | \|z-score\| cutoff for drug-gene network edges |
| `MIN_REPS` | 2 | Minimum replicates per consensus signature |

---

## Using the Enrichment Engine as a Python Library

```python
from cmap_enrichment import EngineRegistry

# Lazily manages one engine per class; reports which datasets are installed.
reg = EngineRegistry(data_dir="data/processed")
print(reg.available())          # {'drug': True, 'knockdown': False, 'overexpression': False}

# Find DRUGS that reverse / mimic your disease signature
res = reg.get("drug").rank_perturbations(
    genes_up=["APOE", "CLU", "TREM2", "CD68"],
    genes_down=["SYN1", "SNAP25", "SLC17A7"],
    cell_line="A549",   # optional: filter to a specific cell line
    top_n=25,
)
res["reversing"]   # candidate therapeutics (most-negative WTCS)
res["mimicking"]   # signature-inducers   (most-positive WTCS)

# Find gene KNOCKDOWNS that reverse it — optionally restrict to one technology
reg.get("knockdown").rank_perturbations(
    genes_up=["APOE", "CLU"], method="CRISPR", top_n=25,
)
```

Or instantiate a single class directly:

```python
from cmap_enrichment import CMapEngine

engine = CMapEngine(data_dir="data/processed", pert_class="drug")
engine.search_perturbagens("statin")     # Fuzzy perturbagen search
engine.get_perturbagen_info("sirolimus")  # All signatures for a perturbagen
engine.get_signature(42)                  # Z-score vector for signature index 42
engine.list_cell_lines()                  # Available cell lines
engine.list_moas()                        # MoAs (drug class only)
```

---

## REST API

Start the FastAPI server:

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8000

# With custom data directory:
CMAP_DATA_DIR=/path/to/data/processed uvicorn api.main:app
```

### Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| `POST` | `/enrich/drug` | Drug (`trt_cp`) enrichment from up/down gene lists |
| `POST` | `/enrich/knockdown` | Knockdown (`trt_sh` + `trt_xpr`) enrichment; accepts `method` |
| `POST` | `/enrich/overexpression` | Overexpression (`trt_oe`) enrichment |
| `GET` | `/classes` | The three classes, their framing, and which are installed |
| `GET` | `/{pert_class}/cell_lines` | Cell lines available for a class |
| `GET` | `/{pert_class}/search?q=...` | Fuzzy perturbagen search within a class |
| `GET` | `/health` | Health check + per-class availability |

Each `/enrich/*` call returns both `reversing` and `mimicking` lists plus `query_stats`. A class whose data isn't installed returns **503**.

### Example

```bash
curl -X POST http://localhost:8000/enrich/drug \
  -H "Content-Type: application/json" \
  -d '{
    "genes_up": ["APOE", "CLU", "TREM2", "CD68", "C1QB"],
    "genes_down": ["SYN1", "SNAP25", "SLC17A7"],
    "top_n": 10
  }'
```

Interactive API docs available at `http://localhost:8000/docs`.

---

## MCP Server (Claude Integration)

The MCP server wraps the enrichment engine for use with Claude and other LLM agents via the [Model Context Protocol](https://modelcontextprotocol.io/).

```bash
python mcp_server/server.py
```

**Available tools** — three framed enrichment tools plus helpers:

| Tool | Description |
|------|-------------|
| `cmap_drug_enrichment` | Small-molecule **drugs** that reverse/mimic the signature (with MoA/target) |
| `cmap_target_knockdown` | Gene **loss-of-function** (shRNA + CRISPR) that reverse/mimic; optional `method` |
| `cmap_target_overexpression` | Gene **gain-of-function** (ORF) that reverse/mimic |
| `cmap_list_perturbation_classes` | Describe the three classes + which datasets are installed |
| `cmap_search_perturbagens` | Search perturbagens (drugs/genes) by partial name |
| `cmap_list_cell_lines` | List available cell lines for a class |

Each tool takes `genes_up` (required), optional `genes_down`, `cell_line`, `direction` (`reversing`/`mimicking`/`both`), and `top_n`; only `cmap_target_knockdown` adds `method` (`shRNA`/`CRISPR`). The result text is narrated with class-specific wording so the LLM frames drugs, knockdowns, and overexpression correctly.

For Claude Desktop / Claude Code / Cursor registration, copy [`docs/claude_desktop_config.example.json`](docs/claude_desktop_config.example.json) into your client's MCP config and replace the absolute paths. Full walkthrough in [`docs/mcp-setup.md`](docs/mcp-setup.md).

---

## Worked Example: 50 MSigDB Hallmarks

The repo ships a full cross-pathway showcase under [`analysis/hallmark_enrichment/`](analysis/hallmark_enrichment/) that runs the enrichment engine against all 50 [MSigDB Hallmark](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H) gene sets, aggregates the results, and produces a master heatmap and bipartite drug ↔ hallmark network.

```bash
sbatch analysis/hallmark_enrichment/slurm_run.sh    # ~12 min on 8 cores
```

| Output | Description |
|--------|-------------|
| `aggregated/cross_hallmark_drug_matrix.csv` | Drugs × 50 Hallmarks score matrix |
| `aggregated/top5_drugs_summary.csv`         | Best broad-spectrum reversal candidates |
| `plots/master_heatmap.png`                  | Top 5 drugs × 50 Hallmarks, grouped by biology |
| `plots/drug_hallmark_network.png`           | Bipartite drug ↔ Hallmark network |

Full walkthrough: [`examples/hallmark_signatures_walkthrough.md`](examples/hallmark_signatures_walkthrough.md).

---

## Documentation

| Doc | What it covers |
|-----|----------------|
| [`docs/python-tutorial.md`](docs/python-tutorial.md)   | Library usage with worked queries |
| [`docs/api-tutorial.md`](docs/api-tutorial.md)         | REST API endpoints, curl + Python client examples |
| [`docs/mcp-setup.md`](docs/mcp-setup.md)               | Wiring the MCP server into Claude Desktop / Claude Code / Cursor |
| [`docs/pipeline.md`](docs/pipeline.md)                 | End-to-end data engineering: GCTX → DuckDB |
| [`examples/`](examples/)                                | Runnable quickstarts and sample gene lists |

---

## Method

This pipeline implements the **Weighted Connectivity Score (WTCS)** from [Subramanian et al., *Cell* (2017)](https://doi.org/10.1016/j.cell.2017.10.049).

### Scoring

1. For each consensus drug signature, all 12,328 genes are ranked by z-score (rank 1 = most downregulated).
2. A **KS-like enrichment score** is computed separately for the up-regulated gene set (ES_up) and the down-regulated gene set (ES_down).
3. The **Weighted Connectivity Score** combines both:

```
WTCS = (ES_up - ES_down) / 2
```

| WTCS | Interpretation |
|------|----------------|
| **Negative** | Drug **reverses** the query signature (therapeutic candidate) |
| **Positive** | Drug **mimics** the query signature |
| **Zero** | No significant connectivity |

### Statistical Significance

- **Drug-level p-values:** 10,000 permutation tests (random gene sets of matching size), then Benjamini-Hochberg FDR correction.
- **MOA enrichment:** Hypergeometric test on the top 5% of reversing drugs, with BH-FDR correction.

### Processing Pipeline

```
Raw GCTX (720K sigs x 12K genes, 34 GB)
    │
    ├── Step 1: Filter to exemplar signatures (is_exemplar_sig == 1)
    │           Map gene IDs to HGNC symbols
    │           Output: Parquet files
    │
    ├── Step 2: Aggregate by (compound, cell line) using median z-scores
    │           Pre-compute rank matrix (scipy.stats.rankdata)
    │           Merge MOA/target annotations from compoundinfo
    │           Output: numpy arrays + metadata
    │
    └── Step 3: Build indexed DuckDB database
                Output: cmap.duckdb
```

---

## Data Source

**[LINCS 2020 Level 5](https://clue.io/data/CMap2020#LINCS2020)** — Replicate-collapsed z-score signatures from the Broad Institute Connectivity Map.

| Statistic | Value |
|-----------|-------|
| Total signatures | 1,201,944 |
| Genes per signature | 12,328 (978 landmark + 11,350 inferred) |
| Compound treatment signatures | 720,216 |
| Consensus signatures (after aggregation) | 3,493 |
| Unique compounds | 511 |
| Cell lines | 41 |

The underlying data are also queryable via **Google BigQuery** (extract arbitrary subsets without downloading the static files): [`cmap-big-table` public views](https://console.cloud.google.com/bigquery?p=cmap-big-table&d=cmap_lincs_public_views&page=dataset).

**Upstream reference:** [GCTX format](https://clue.io/connectopedia/gctx_format) · [Data levels](https://clue.io/connectopedia/data_levels)

---

## Citation

If you use this pipeline, please cite the underlying CMap methodology:

> Subramanian, A., et al. "A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles." *Cell* 171.6 (2017): 1437-1452. [doi:10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)

---

## License

This pipeline is open source. The underlying LINCS 2020 data is provided by the Broad Institute under their [data use policy](https://clue.io/connectopedia/data_use_policy).
