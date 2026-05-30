---
hide:
  - navigation
---

<div class="clue-hero" markdown>

# clue-mcp

<p class="tagline">
Turn a disease gene-expression signature into a ranked list of candidate drugs,
gene-knockdown targets, and overexpression hits —
from Python, a REST API, or directly inside Claude via the Model Context Protocol.
</p>

<p class="clue-badges">
<a href="https://github.com/swaruplab/clue-mcp"><img alt="GitHub" src="https://img.shields.io/badge/GitHub-clue--mcp-181717?logo=github"></a>
<img alt="License" src="https://img.shields.io/badge/license-MIT-blue">
<img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white">
<img alt="Status" src="https://img.shields.io/badge/status-beta-orange">
</p>

[Quickstart :material-rocket-launch:](quickstart.md){ .md-button .md-button--primary }
[See the showcase :material-chart-box:](showcase.md){ .md-button }
[Use it in Claude :material-robot:](mcp-setup.md){ .md-button }

</div>

---

## What it does

You have a gene signature — a set of genes **up-** and **down-regulated** in your disease, cell state, or perturbation. **clue-mcp** asks which perturbations push expression in the **opposite** direction (reverse it) or the **same** direction (mimic it), across three biologically distinct classes:

> - **Drugs** — which small molecules might *reverse* this state? → drug-repurposing leads.
> - **Knockdowns** — which gene *loss-of-function* perturbations reverse it? → candidate driver targets to silence (shRNA + CRISPR).
> - **Overexpression** — which gene *gain-of-function* perturbations reverse it? → genes whose forced expression opposes the state.

It answers by scoring your signature against the **Connectivity Map (CMap) LINCS 2020** perturbation database from the Broad Institute, using the **Weighted Connectivity Score (WTCS)** from [Subramanian *et al.*, *Cell* 2017](https://doi.org/10.1016/j.cell.2017.10.049). Each class is exposed as its own tool/endpoint with tailored scientific framing, so the most-negative-scoring perturbations are your top reversal candidates *of that kind*.

<div class="clue-stats" markdown>
<span class="stat"><span class="num">3,493</span><span class="label">Consensus signatures</span></span>
<span class="stat"><span class="num">511</span><span class="label">Compounds</span></span>
<span class="stat"><span class="num">41</span><span class="label">Cell lines</span></span>
<span class="stat"><span class="num">12,328</span><span class="label">Genes scored</span></span>
</div>

The consensus signatures are aggregated from **720,216** LINCS 2020 Level 5 compound treatment profiles. See [Method & data](method.md) for the full provenance.

---

## Four ways to use it

<div class="clue-grid" markdown>

<div class="clue-card" markdown>
### :material-robot: MCP server
Ask Claude, in plain English, *"what drugs reverse this signature, and what genes should I knock down?"* and let it call the engine for you. Three framed enrichment tools (drug / knockdown / overexpression) plus helpers, over the [Model Context Protocol](https://modelcontextprotocol.io/).

[Set up the MCP server →](mcp-setup.md)
</div>

<div class="clue-card" markdown>
### :material-language-python: Python library
`from cmap_enrichment import CMapEngine`. Drop it into notebooks, batch jobs, and custom pipelines.

[Python tutorial →](python-tutorial.md)
</div>

<div class="clue-card" markdown>
### :material-api: REST API
A FastAPI server exposing the same engine over HTTP/JSON, with interactive Swagger docs.

[REST API tutorial →](api-tutorial.md)
</div>

<div class="clue-card" markdown>
### :material-chart-bell-curve: Analysis driver
One command → permutation FDR, MOA enrichment, and 11 publication-quality figure panels.

[Running an analysis →](usage.md)
</div>

</div>

---

## A 30-second taste

```python
from cmap_enrichment import CMapEngine

engine = CMapEngine(data_dir="data/processed")

hits = engine.query_enrichment(
    genes_up=["APOE", "CLU", "TREM2", "CD68", "C1QB"],   # e.g. AD microglia
    genes_down=["SYN1", "SNAP25", "SLC17A7"],
    top_n=10,
)
print(hits[["cmap_name", "wtcs", "moa", "target"]])
# Most-negative WTCS = strongest reversal = top therapeutic candidate
```

Inside Claude (with the [MCP server](mcp-setup.md) registered) the same query is just a sentence:

> *I have an Alzheimer's microglia signature: up = APOE, CLU, TREM2, CD68, C1QB; down = SYN1, SNAP25, SLC17A7. What drugs might reverse it?*

---

## See it on real biology

The repo ships a fully worked, pre-computed **[showcase](showcase.md)**: clue-mcp run against **all 50 MSigDB Hallmark gene sets**, with a cross-pathway drug heatmap and a drug ↔ hallmark network.

[![Master heatmap of top drugs across 50 Hallmarks](assets/master_heatmap.png){ width="640" }](showcase.md)

---

## Get started

1. **[Quickstart](quickstart.md)** — install and run your first query.
2. **[Building the database](pipeline.md)** — how the LINCS 2020 GCTX files become the queryable engine (only needed once, on an HPC node).
3. **[Method & data](method.md)** — the WTCS scoring, statistics, and data provenance.

!!! note "Do I need the 9 GB database?"
    The repo ships **code + the small precomputed [showcase](showcase.md)** — not the
    ~9 GB processed database (too big for git, and a derived product of the Broad's data).
    Three tiers:

    - **Look (0 GB):** browse the [showcase](showcase.md) results — no download at all.
    - **Run (download ~9 GB):** `bash scripts/download_data.sh` pulls the prebuilt database from Zenodo — no HPC needed. *(Recommended.)*
    - **Rebuild (34 GB + HPC):** generate it yourself via [Building the database](pipeline.md).

    Already have a copy? Just `export CMAP_DATA_DIR=/path/to/data/processed`.

---

<small>
clue-mcp is developed by the [Swarup Lab](https://swaruplab.bio.uci.edu/) at UC Irvine and released under the MIT License.
The underlying LINCS 2020 / Connectivity Map data is © the Broad Institute, used under its
[data use policy](https://clue.io/connectopedia/data_use_policy) and **not** redistributed here.
</small>
