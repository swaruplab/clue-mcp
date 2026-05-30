# Running an analysis

The library returns raw WTCS scores. For a **full, publication-ready analysis** — permutation FDR, MOA enrichment, and 11 figure panels — use the analysis driver:

```bash
python analysis/run_enrichment.py /path/to/your_genes.csv
```

This page documents its inputs, outputs, and parameters.

---

## Input gene-file formats

`analysis/run_enrichment.py` auto-detects the format of `.txt`, `.csv`, and `.tsv` files (delimiter and optional header are detected automatically).

=== "Format A — symbols only"

    One gene per line. Treated as a single **upregulated** set (one-sided KS enrichment).

    ```text
    APOE
    CLU
    TREM2
    CD68
    C1QB
    ```

=== "Format B — symbols + values"

    Two columns: gene symbol and a numeric value (log2FC, kME, …).
    **Positive = upregulated, negative = downregulated.** The driver computes the
    two-sided **WTCS** (up vs. down gene sets).

    ```csv
    gene_symbol,log2FC
    APOE,2.31
    CLU,1.87
    SYN1,-1.54
    SNAP25,-2.01
    ```

---

## Run it

```bash
# Direct (loads the matrices into ~10 GB RAM)
python analysis/run_enrichment.py /path/to/your_genes.csv

# On an HPC scheduler (recommended)
GENE_FILE=/path/to/your_genes.csv sbatch analysis/slurm_run.sh
```

Outputs are written next to the script.

---

## Output tables

| File | Description |
|------|-------------|
| `drug_level_summary.csv` | Every compound ranked by median enrichment score, with permutation p-values and FDR |
| `top50_reversing_signatures.csv` | 50 individual signatures with the strongest reversal |
| `top50_mimicking_signatures.csv` | 50 individual signatures with the strongest mimicry |
| `all_signature_scores.csv` | Complete signature-level enrichment scores |
| `moa_enrichment.csv` | Mechanism-of-action enrichment (hypergeometric test + FDR) |
| `cell_line_stats.csv` | Per-cell-line score statistics |
| `summary.txt` | Human-readable results summary |

## Output figures

Each multi-panel figure renders **(A) all drugs** and **(B) named drugs only** (excluding BRD-coded identifiers).

| Figure | Description |
|--------|-------------|
| `plot1_score_distribution.png` | Score distribution with 5th/95th-percentile cutoffs |
| `plot2/3_top25_*_signature_level.png` | Top 25 reversing / mimicking signatures |
| `plot4/5_top25_*_drug_level.png` | Top 25 reversing / mimicking drugs (median across cell lines, FDR-annotated) |
| `plot6_moa_enrichment.png` | MOA fold-enrichment among the top 5% reversing drugs |
| `plot7_cellline_comparison.png` | Median enrichment score by cell line |
| `plot8a/b_heatmap_*.png` | Drug × cell-line enrichment heatmaps |
| `plot9a/b_zscore_heatmap_*.png` | Z-score heatmaps of query genes in top reversing signatures |
| `plot10_volcano.png` | Volcano plot (median score vs. −log₁₀ FDR) |
| `plot11a/b_network_*.png` | Drug ↔ gene bipartite network (edges where \|z\| > threshold) |

---

## Configurable parameters

Set at the top of `analysis/run_enrichment.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_PERM` | 10,000 | Permutations for the FDR null distribution |
| `TOP_N` | 50 | Top results saved per direction |
| `ZSCORE_THRESHOLD` | 1.5 | \|z-score\| cutoff for drug–gene network edges |
| `MIN_REPS` | 2 | Minimum replicates per consensus signature (quality filter) |

---

## Scaling to many gene sets

To run the driver across a whole gene-set collection at once (and aggregate into a
cross-pathway heatmap + network), see the **[50-Hallmark showcase](showcase.md)**, which
does exactly this for all 50 MSigDB Hallmark sets.

## Next steps

- [Method & data](method.md) — how the WTCS and the statistics are computed.
- [Python library](python-tutorial.md) — call the engine directly without the driver.
