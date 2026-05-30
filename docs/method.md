# Method & data

clue-mcp implements the **Weighted Connectivity Score (WTCS)** from
[Subramanian *et al.*, *Cell* 2017](https://doi.org/10.1016/j.cell.2017.10.049),
the method behind the Broad Institute's Connectivity Map.

---

## Scoring

1. For each consensus drug signature, all **12,328 genes** are ranked by z-score (rank 1 = most downregulated by the drug).
2. A **KS-like enrichment score** is computed separately for your up-regulated gene set (`ES_up`) and your down-regulated gene set (`ES_down`) — each measures how strongly that set clusters toward one end of the drug's ranking.
3. The two are combined into the **Weighted Connectivity Score**:

    ```
    WTCS = (ES_up − ES_down) / 2
    ```

| WTCS | Interpretation |
|------|----------------|
| **Negative** | Drug **reverses** the query signature → therapeutic candidate |
| **Positive** | Drug **mimics** the query signature |
| **≈ 0** | No significant connectivity (or discordant up/down enrichment) |

If you supply only `genes_up`, the engine falls back to a single-set, one-sided KS enrichment.

---

## Statistical significance

The [analysis driver](usage.md) adds two layers of significance testing on top of the raw scores:

- **Drug-level p-values** — 10,000 permutation tests using random gene sets of matching size, followed by Benjamini–Hochberg FDR correction.
- **MOA enrichment** — a hypergeometric test for mechanisms of action over-represented among the top 5% of reversing drugs, again with BH-FDR correction.

---

## Processing pipeline

Raw GCTX → queryable engine, in three steps (detailed in [Building the database](pipeline.md)):

```text
Raw GCTX (720K signatures × 12K genes, 34 GB)
    │
    ├── Step 1  Filter to exemplar signatures (is_exemplar_sig == 1);
    │           map gene IDs → HGNC symbols.                      → Parquet
    │
    ├── Step 2  Aggregate by (compound, cell line) via median z-scores;
    │           pre-compute the rank matrix; merge MOA/target.    → numpy + metadata
    │
    └── Step 3  Index metadata into DuckDB.                       → cmap.duckdb
```

The enrichment itself runs against the pre-ranked numpy matrices (fast); DuckDB handles
all metadata filtering and lookups.

---

## Data source

**[LINCS 2020 Level 5](https://clue.io/data/CMap2020#LINCS2020)** — replicate-collapsed
z-score signatures from the Broad Institute Connectivity Map.

| Statistic | Value |
|-----------|-------|
| Total LINCS 2020 signatures | 1,201,944 |
| Genes per signature | 12,328 (978 landmark + 11,350 inferred) |
| Compound-treatment signatures | 720,216 |
| Consensus signatures (after aggregation) | 3,493 |
| Unique compounds | 511 |
| Cell lines | 41 |

!!! warning "The CMap data is not redistributed here"
    clue-mcp ships only code and small, derived example outputs. The underlying
    LINCS 2020 data is © the Broad Institute and must be downloaded directly from
    [clue.io](https://clue.io/data/CMap2020#LINCS2020) under its
    [data use policy](https://clue.io/connectopedia/data_use_policy).

---

## Citation

If you use clue-mcp, please cite the underlying CMap methodology:

> Subramanian, A., *et al.* "A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles." *Cell* 171.6 (2017): 1437–1452. [doi:10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)
