# Hallmark Signatures Walkthrough

This walks through the bundled showcase that runs the CMap enrichment engine against **all 50 MSigDB Hallmark gene sets**, aggregates the results, and produces a cross-pathway heatmap and network plot.

It is both a useful biological screen *and* a demonstration of how to scale the engine to many parallel queries.

---

## What it does

For each of the 50 Hallmark gene sets:

1. Looks up the gene list from MSigDB (cached in `data/msigdb/h.all.v2026.1.Hs.json`).
2. Maps gene symbols onto the L1000 12,328-gene panel.
3. Calls `CMapEngine.query_enrichment(genes_up=...)` for every drug signature.
4. Runs 10,000 permutations of size-matched random gene sets to compute empirical p-values.
5. Applies Benjamini–Hochberg FDR correction across all drug-hallmark pairs.

After all 50 hallmarks complete:

- Aggregates into a 511-drug × 50-hallmark matrix of mean enrichment scores.
- Ranks drugs by a composite of effect size and breadth (how many hallmarks they significantly reverse).
- Renders a master heatmap of the top 5 drugs across all 50 hallmarks, color-coded by biological category (immune, proliferation, signaling, metabolic, development, stress/DNA damage, other).
- Renders a bipartite drug ↔ hallmark network, with edge weight = significance.

## Run it

```bash
# On SLURM (recommended — ~12 minutes on 8 cores, 64 GB RAM)
sbatch analysis/hallmark_enrichment/slurm_run.sh

# Or locally
python analysis/hallmark_enrichment/run_hallmark_enrichment.py
```

## Outputs

```
analysis/hallmark_enrichment/
├── per_hallmark/                            # one folder per Hallmark
│   └── HALLMARK_E2F_TARGETS/
│       ├── all_signature_scores.csv
│       ├── drug_level_summary.csv
│       ├── top50_reversing_signatures.csv
│       └── top50_mimicking_signatures.csv
├── aggregated/
│   ├── cross_hallmark_drug_matrix.csv       # drugs x 50 hallmarks score matrix
│   ├── top5_drugs_summary.csv
│   └── hallmark_statistics.csv
├── plots/
│   ├── master_heatmap.png                   # top 5 drugs x 50 hallmarks
│   ├── expanded_heatmap_top20.png
│   └── drug_hallmark_network.png
└── summary.txt
```

## Headline result (from the bundled run)

The 5 most broadly drug-reversing compounds across all 50 Hallmarks were:

| Rank | Compound | MOA | # significant hallmarks |
|------|----------|-----|--------------------------|
| 1 | staurosporine  | PKC / CDK inhibitor      | 18 |
| 2 | dactinomycin   | RNA polymerase inhibitor | 14 |
| 3 | A-443654       | AKT inhibitor            | 12 |
| 4 | emetine        | Protein synthesis inhibitor | 11 |
| 5 | PIK-75         | PI3K inhibitor           | 11 |

Proliferation-associated Hallmarks (E2F Targets, MYC Targets, G2/M Checkpoint) were the most consistently reversed; immune Hallmarks were predominantly *mimicked* (positive WTCS).

See the full write-up in [the bundled report](../report_2026-04-05_2110.md) (if present in your checkout).

## Adapt it

To run the same workflow on a different gene-set collection (KEGG, GO BP, Reactome, your own), edit `run_hallmark_enrichment.py` and:

1. Replace `MSIGDB_URL` / `MSIGDB_JSON` with your gene-set source.
2. Update the `HALLMARK_CATEGORIES` dict if you want biological column groupings.
3. Adjust `TOP_DRUGS` if you want a wider master heatmap.
