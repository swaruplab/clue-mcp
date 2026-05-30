#!/usr/bin/env python3
"""
CMap Enrichment Analysis — v3

Generic input:
  - Accepts any .txt, .csv, .tsv gene list
  - 1-column: gene symbols only (treated as upregulated set)
  - 2-column: gene_symbol + numeric value (log2FC, kME, etc.)
    → positive values = upregulated, negative = downregulated
    → magnitude used for z-score weighted ranking

Each figure panel has two subplots:
  (A) All drugs   (B) Named drugs only (excluding BRD- coded IDs)

Usage:
  python run_enrichment.py                      # uses default test_genes.txt
  python run_enrichment.py /path/to/genes.csv   # custom gene file
"""

import sys, os, re
import numpy as np
import duckdb
import csv
from pathlib import Path
from collections import Counter, defaultdict
from scipy.stats import rankdata, hypergeom, zscore as scipy_zscore
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import networkx as nx

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════
BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed"
OUT_DIR = Path(__file__).resolve().parent

# Accept gene file as CLI argument, fallback to test_genes.txt
GENE_FILE = Path(sys.argv[1]) if len(sys.argv) > 1 else BASE_DIR / "test_genes.txt"

N_PERM = 10000
TOP_N = 50
ZSCORE_THRESHOLD = 1.5   # for network edges
MIN_REPS = 2             # minimum replicates per signature (quality filter)

# Color palette
C_REV = '#D7263D'     # reversing / negative — vivid red
C_MIM = '#1B998B'     # mimicking / positive — teal
C_BLUE = '#3A86FF'    # neutral / distribution
C_PURPLE = '#7B2D8E'  # MOA
C_GOLD = '#F4A261'    # cell lines
C_GREY = '#BFBFBF'    # non-significant / coded IDs
C_BG = '#FAFAFA'      # background

# ═══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def safe_str(val):
    """Convert NaN/None to empty string."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ''
    return str(val)


def is_named_drug(name):
    """True if NOT a BRD- coded identifier."""
    return not name.startswith('BRD-')


def parse_gene_file(filepath):
    """
    Parse a gene file in any common format.

    Supported formats:
      - .txt, .csv, .tsv (auto-detected delimiter)
      - 1-column: gene symbols (one per line, or comma-separated on one line)
      - 2-column: gene_symbol, numeric_value (log2FC, kME, etc.)
      - Optional header row (auto-detected)

    Returns:
      genes: list of gene symbol strings
      values: list of float values (or None if 1-column)
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    with open(filepath) as f:
        raw = f.read().strip()

    if not raw:
        raise ValueError(f"Gene file is empty: {filepath}")

    lines = raw.split('\n')

    # Detect delimiter
    if ext == '.csv':
        delim = ','
    elif ext == '.tsv':
        delim = '\t'
    else:
        # Auto-detect: check first data line for tabs, then commas
        first_line = lines[0]
        if '\t' in first_line:
            delim = '\t'
        elif ',' in first_line:
            delim = ','
        else:
            delim = None  # single column, whitespace or one-per-line

    # Special case: single line with many comma-separated genes
    # (like the original test_genes.txt format: HEADER\tGENE1,GENE2,...)
    if len(lines) <= 2:
        parts = lines[0].split('\t')
        if len(parts) == 2 and ',' in parts[1]:
            # Old format: HEADER\tGENE1,GENE2,...
            gene_list = [g.strip() for g in parts[1].split(',') if g.strip()]
            return gene_list, None
        elif len(parts) == 1 and ',' in parts[0]:
            gene_list = [g.strip() for g in parts[0].split(',') if g.strip()]
            return gene_list, None

    # Parse as rows
    rows = []
    for line in lines:
        if delim:
            fields = [f.strip() for f in line.split(delim)]
        else:
            fields = line.strip().split()
        if fields and fields[0]:
            rows.append(fields)

    if not rows:
        raise ValueError(f"No data found in {filepath}")

    # Detect header: if first row's second field (if exists) is not numeric
    has_header = False
    n_cols = max(len(r) for r in rows)
    if n_cols >= 2:
        try:
            float(rows[0][1])
        except (ValueError, IndexError):
            has_header = True
    elif n_cols == 1:
        # Check if first entry looks like a header keyword
        if rows[0][0].upper() in ('GENE', 'GENE_SYMBOL', 'SYMBOL', 'GENES',
                                   'GENE_SYMBOLS', 'GENE_NAME', 'NAME'):
            has_header = True

    data_rows = rows[1:] if has_header else rows

    genes = []
    values = []
    has_values = False

    for row in data_rows:
        gene = row[0].strip()
        if not gene:
            continue
        genes.append(gene)
        if len(row) >= 2:
            try:
                val = float(row[1])
                values.append(val)
                has_values = True
            except ValueError:
                values.append(None)
        else:
            values.append(None)

    if has_values and any(v is not None for v in values):
        # Fill None values with 0
        values = [v if v is not None else 0.0 for v in values]
        return genes, values
    else:
        return genes, None


def make_dual_barh(data_all, data_named, val_key, label_fn, color_all, color_named,
                   title_all, title_named, xlabel, filename, n=25, fdr_key=None):
    """Create side-by-side horizontal bar charts (all vs named drugs)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9), sharey=False)

    for ax, data, color, title in [
        (ax1, data_all[:n], color_all, title_all),
        (ax2, data_named[:n], color_named, title_named),
    ]:
        if not data:
            ax.set_visible(False)
            continue
        labels = [label_fn(d) for d in data]
        vals = [d[val_key] for d in data]
        y_pos = np.arange(len(labels))

        if fdr_key:
            bar_colors = [color if d.get(fdr_key, 1) < 0.05 else C_GREY for d in data]
        else:
            bar_colors = color

        bars = ax.barh(y_pos, vals, color=bar_colors, edgecolor='white',
                       height=0.72, zorder=3)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.axvline(0, color='black', linewidth=0.5, zorder=2)
        ax.grid(axis='x', alpha=0.2, zorder=1)
        ax.set_facecolor('#FCFCFC')

        # FDR annotations
        if fdr_key:
            for i, d in enumerate(data):
                fdr_val = d.get(fdr_key, 1)
                fdr_txt = f"FDR={fdr_val:.3f}" if fdr_val >= 0.001 else "FDR<0.001"
                if vals[i] < 0:
                    ax.text(vals[i] + 0.001, i, fdr_txt, va='center', ha='left',
                            fontsize=6.5, color='#333333')
                else:
                    ax.text(vals[i] - 0.001, i, fdr_txt, va='center', ha='right',
                            fontsize=6.5, color='#333333')

    fig.tight_layout(w_pad=3)
    fig.savefig(OUT_DIR / filename, dpi=200, bbox_inches='tight')
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("CMap Enrichment Analysis v3")
print("=" * 65)

print(f"\nGene file: {GENE_FILE}")
print("Loading CMap database...")
t0 = time.time()

rank_matrix = np.load(DATA_DIR / "rank_matrix.npy", mmap_mode='r')
zscore_matrix = np.load(DATA_DIR / "zscore_matrix.npy", mmap_mode='r')
with open(DATA_DIR / "gene_names.txt") as f:
    gene_names = [g.strip() for g in f.readlines()]
gene_to_idx = {g: i for i, g in enumerate(gene_names)}
n_sigs, n_genes = rank_matrix.shape

con = duckdb.connect(str(DATA_DIR / "cmap.duckdb"), read_only=True)
print(f"  {n_sigs:,} signatures x {n_genes:,} genes ({time.time()-t0:.1f}s)")

# ═══════════════════════════════════════════════════════════════════════════
# PARSE INPUT GENES
# ═══════════════════════════════════════════════════════════════════════════
print("\nParsing gene file...")
raw_genes, gene_values = parse_gene_file(GENE_FILE)

mapped_genes = [g for g in raw_genes if g in gene_to_idx]
unmapped_genes = [g for g in raw_genes if g not in gene_to_idx]
gene_idx_list = np.array([gene_to_idx[g] for g in mapped_genes])
n_set = len(gene_idx_list)

has_values = gene_values is not None
if has_values:
    # Map values to mapped genes only
    gene_val_map = {g: v for g, v in zip(raw_genes, gene_values) if g in gene_to_idx}
    mapped_values = np.array([gene_val_map[g] for g in mapped_genes], dtype=np.float32)
    # Convert to z-scores for weighting
    if np.std(mapped_values) > 0:
        gene_weights = scipy_zscore(mapped_values)
    else:
        gene_weights = mapped_values
    # Split into up/down based on sign of original values
    up_mask = mapped_values >= 0
    down_mask = mapped_values < 0
    up_idx = gene_idx_list[up_mask]
    down_idx = gene_idx_list[down_mask]
    print(f"  Input: {len(raw_genes)} genes with values")
    print(f"  Mapped: {n_set} (up: {up_mask.sum()}, down: {down_mask.sum()})")
    print(f"  Value range: [{mapped_values.min():.3f}, {mapped_values.max():.3f}]")
else:
    gene_weights = None
    up_idx = gene_idx_list
    down_idx = np.array([], dtype=int)
    print(f"  Input: {len(raw_genes)} genes (no values — treating all as upregulated)")
    print(f"  Mapped: {n_set}")

if unmapped_genes:
    print(f"  Unmapped ({len(unmapped_genes)}): {', '.join(unmapped_genes[:20])}"
          + (f"... +{len(unmapped_genes)-20} more" if len(unmapped_genes) > 20 else ""))

if n_set == 0:
    print("ERROR: No query genes mapped to CMap genes. Exiting.")
    sys.exit(1)

# ═══════════════════════════════════════════════════════════════════════════
# QUALITY FILTER — restrict to signatures with sufficient replicates
# ═══════════════════════════════════════════════════════════════════════════
if MIN_REPS > 1:
    quality_idx = con.execute(
        f"SELECT sig_idx FROM signatures WHERE n_reps >= {MIN_REPS} ORDER BY sig_idx"
    ).fetchnumpy()["sig_idx"]
    print(f"\nQuality filter: n_reps >= {MIN_REPS}")
    print(f"  {len(quality_idx):,} / {n_sigs:,} signatures pass ({100*len(quality_idx)/n_sigs:.1f}%)")
else:
    quality_idx = np.arange(n_sigs)
    print(f"\nNo quality filter applied ({n_sigs:,} signatures)")

# ═══════════════════════════════════════════════════════════════════════════
# COMPUTE ENRICHMENT SCORES
# ═══════════════════════════════════════════════════════════════════════════
print(f"\nComputing enrichment scores across {len(quality_idx):,} signatures...")
t1 = time.time()

if has_values and len(down_idx) > 0:
    # WTCS: separate up and down gene sets
    def vectorized_ks(rank_mat, gene_indices):
        """Vectorized KS enrichment for a gene set across all signatures."""
        n_g = len(gene_indices)
        if n_g == 0:
            return np.zeros(rank_mat.shape[0], dtype=np.float32)
        qr = np.array(rank_mat[:, gene_indices], dtype=np.float32)
        qr_sorted = np.sort(qr, axis=1)
        expected = (np.arange(1, n_g + 1) / n_g).astype(np.float32)
        hit = qr_sorted / n_genes
        dev = hit - expected[None, :]
        max_idx = np.argmax(np.abs(dev), axis=1)
        return dev[np.arange(dev.shape[0]), max_idx]

    es_up = vectorized_ks(rank_matrix, up_idx)
    es_down = vectorized_ks(rank_matrix, down_idx)

    # WTCS: (es_up - es_down) / 2 when concordant, else 0
    concordant = np.sign(es_up) != np.sign(es_down)
    scores = np.where(concordant, (es_up - es_down) / 2, 0.0).astype(np.float32)
    # Where one is zero, fall back to the non-zero one
    both_zero = (es_up == 0) & (es_down == 0)
    scores[both_zero] = 0.0
    scoring_method = "WTCS (up + down gene sets)"

else:
    # Single gene set — standard KS
    query_ranks = np.array(rank_matrix[:, gene_idx_list], dtype=np.float32)
    query_ranks_sorted = np.sort(query_ranks, axis=1)
    expected = (np.arange(1, n_set + 1) / n_set).astype(np.float32)
    hit_scores = query_ranks_sorted / n_genes
    deviations = hit_scores - expected[None, :]
    max_dev_idx = np.argmax(np.abs(deviations), axis=1)
    scores = deviations[np.arange(n_sigs), max_dev_idx]
    scoring_method = "KS (single gene set)"

print(f"  Method: {scoring_method}")
print(f"  Done in {time.time()-t1:.1f}s")
print(f"  Score range: [{scores.min():.4f}, {scores.max():.4f}]")
print(f"  Mean: {scores.mean():.4f}, Std: {scores.std():.4f}")

# ═══════════════════════════════════════════════════════════════════════════
# METADATA & DRUG-LEVEL AGGREGATION
# ═══════════════════════════════════════════════════════════════════════════
print("\nLoading metadata & aggregating...")
t2 = time.time()

sig_meta_df = con.execute("SELECT * FROM signatures ORDER BY sig_idx").fetchdf()
meta_dict = {}
for _, row in sig_meta_df.iterrows():
    meta_dict[row['sig_idx']] = row.to_dict()

quality_set = set(quality_idx)

drug_scores = defaultdict(list)
drug_cells = defaultdict(set)
drug_moa = {}
drug_target = {}

for idx in range(n_sigs):
    if idx not in quality_set:
        continue
    m = meta_dict[idx]
    drug = m['cmap_name']
    drug_scores[drug].append(float(scores[idx]))
    drug_cells[drug].add(m['cell_iname'])
    moa_val = m.get('moa')
    if drug not in drug_moa and safe_str(moa_val):
        drug_moa[drug] = safe_str(moa_val)
        drug_target[drug] = safe_str(m.get('target', ''))

drug_summary = []
for drug in drug_scores:
    s = np.array(drug_scores[drug])
    drug_summary.append({
        'drug': drug,
        'median_score': float(np.median(s)),
        'mean_score': float(np.mean(s)),
        'min_score': float(np.min(s)),
        'max_score': float(np.max(s)),
        'std_score': float(np.std(s)) if len(s) > 1 else 0.0,
        'n_signatures': len(s),
        'n_cell_lines': len(drug_cells[drug]),
        'frac_negative': float(np.mean(s < 0)),
        'moa': drug_moa.get(drug, ''),
        'target': drug_target.get(drug, ''),
    })

drug_summary.sort(key=lambda x: x['median_score'])
print(f"  {len(drug_summary):,} compounds ({time.time()-t2:.1f}s)")

# ── Permutation p-values ─────────────────────────────────────────────────
print("Computing permutation p-values...")
t3 = time.time()
rng = np.random.default_rng(42)
all_scores_arr = scores[quality_idx].copy()

drug_pvals = {}
n_sig_groups = defaultdict(list)
for d in drug_summary:
    n_sig_groups[d['n_signatures']].append(d)

for n_s, grp in n_sig_groups.items():
    null_medians = np.array([
        np.median(rng.choice(all_scores_arr, size=n_s, replace=False))
        for _ in range(N_PERM)
    ], dtype=np.float32)

    for d in grp:
        p_rev = float(max(np.mean(null_medians <= d['median_score']), 1/N_PERM))
        p_mim = float(max(np.mean(null_medians >= d['median_score']), 1/N_PERM))
        drug_pvals[d['drug']] = {'pval_reversal': p_rev, 'pval_mimicking': p_mim}

# BH FDR correction
def bh_fdr(items, pval_key, fdr_key):
    ranked = sorted(items, key=lambda x: x[1][pval_key])
    n = len(ranked)
    for i, (drug, pv) in enumerate(ranked, 1):
        pv[fdr_key] = min(pv[pval_key] * n / i, 1.0)
    prev = 1.0
    for drug, pv in reversed(ranked):
        pv[fdr_key] = min(pv[fdr_key], prev)
        prev = pv[fdr_key]

bh_fdr(list(drug_pvals.items()), 'pval_reversal', 'fdr_reversal')
bh_fdr(list(drug_pvals.items()), 'pval_mimicking', 'fdr_mimicking')

for d in drug_summary:
    d.update(drug_pvals[d['drug']])

n_sig_rev = sum(1 for d in drug_summary if d['fdr_reversal'] < 0.05)
n_sig_mim = sum(1 for d in drug_summary if d['fdr_mimicking'] < 0.05)
print(f"  Done ({time.time()-t3:.1f}s). Sig reversals: {n_sig_rev}, Sig mimickers: {n_sig_mim}")

# ── MOA enrichment (top 5%) ─────────────────────────────────────────────
print("Computing MOA enrichment...")
top_pct = 0.05
n_top = max(1, int(len(drug_summary) * top_pct))
top_rev_set = {d['drug'] for d in drug_summary[:n_top]}

all_moa_drugs = {d['drug'] for d in drug_summary if safe_str(d['moa'])}
sig_rev_moa_drugs = {d['drug'] for d in drug_summary if d['drug'] in top_rev_set and safe_str(d['moa'])}

moa_to_all = defaultdict(set)
moa_to_sig = defaultdict(set)
for d in drug_summary:
    ms = safe_str(d['moa'])
    if not ms:
        continue
    for moa in ms.split('|'):
        moa = moa.strip()
        if not moa:
            continue
        moa_to_all[moa].add(d['drug'])
        if d['drug'] in top_rev_set:
            moa_to_sig[moa].add(d['drug'])

N_total = len(all_moa_drugs)
K_sig = len(sig_rev_moa_drugs)

moa_enrichment = []
for moa in moa_to_all:
    n_moa = len(moa_to_all[moa])
    k_hit = len(moa_to_sig.get(moa, set()))
    if k_hit < 2:
        continue
    pval = hypergeom.sf(k_hit - 1, N_total, n_moa, K_sig)
    moa_enrichment.append({
        'moa': moa, 'n_drugs_total': n_moa, 'n_drugs_significant': k_hit,
        'fold_enrichment': (k_hit / K_sig) / (n_moa / N_total) if K_sig > 0 and N_total > 0 else 0,
        'pval': pval,
    })
moa_enrichment.sort(key=lambda x: x['pval'])
n_moa_tests = len(moa_enrichment)
for i, m in enumerate(moa_enrichment, 1):
    m['fdr'] = min(m['pval'] * n_moa_tests / i, 1.0)
prev = 1.0
for m in reversed(moa_enrichment):
    m['fdr'] = min(m['fdr'], prev)
    prev = m['fdr']

n_sig_moa = sum(1 for m in moa_enrichment if m['fdr'] < 0.05)
print(f"  Enriched MOAs (FDR<0.05): {n_sig_moa}")

# ── Per-cell-line stats ──────────────────────────────────────────────────
cell_scores = defaultdict(list)
for idx in quality_idx:
    cell_scores[meta_dict[idx]['cell_iname']].append(float(scores[idx]))

cell_stats = []
for cell, s in cell_scores.items():
    s_arr = np.array(s)
    cell_stats.append({
        'cell_line': cell, 'n_signatures': len(s),
        'mean_score': float(np.mean(s_arr)),
        'median_score': float(np.median(s_arr)),
        'frac_negative': float(np.mean(s_arr < 0)),
    })
cell_stats.sort(key=lambda x: x['median_score'])

# ── Top signature results (quality-filtered) ────────────────────────────
quality_scores = scores[quality_idx]
quality_order = np.argsort(quality_scores)
sorted_idx = quality_idx[quality_order]
reverse_idx = sorted_idx[:TOP_N]
mimic_idx = sorted_idx[-TOP_N:][::-1]

def build_results(indices):
    return [{
        'sig_idx': int(idx),
        'cmap_name': meta_dict[idx]['cmap_name'],
        'cell_iname': meta_dict[idx]['cell_iname'],
        'moa': safe_str(meta_dict[idx].get('moa')),
        'target': safe_str(meta_dict[idx].get('target')),
        'dose_mode': meta_dict[idx].get('dose_mode') or '',
        'time_mode': meta_dict[idx].get('time_mode') or '',
        'n_reps': meta_dict[idx]['n_reps'],
        'enrichment_score': float(scores[idx]),
    } for idx in indices]

results_reverse = build_results(reverse_idx)
results_mimic = build_results(mimic_idx)

# ── Named-drug subsets ───────────────────────────────────────────────────
named_drug_summary = [d for d in drug_summary if is_named_drug(d['drug'])]
named_drug_summary_rev = [d for d in named_drug_summary if d['n_cell_lines'] >= 2]
named_drug_summary_mim = sorted(
    [d for d in named_drug_summary if d['n_cell_lines'] >= 2],
    key=lambda x: -x['median_score']
)
all_drug_summary_rev = [d for d in drug_summary if d['n_cell_lines'] >= 2]
all_drug_summary_mim = sorted(
    [d for d in drug_summary if d['n_cell_lines'] >= 2],
    key=lambda x: -x['median_score']
)

named_set = {d['drug'] for d in named_drug_summary}
named_reverse_idx = [i for i in sorted_idx if meta_dict[i]['cmap_name'] in named_set][:TOP_N]
named_mimic_idx = [i for i in sorted_idx[::-1] if meta_dict[i]['cmap_name'] in named_set][:TOP_N]

n_named = len(named_drug_summary)
n_coded = len(drug_summary) - n_named
print(f"\nNamed drugs: {n_named:,}  |  BRD- coded: {n_coded:,}")

# ═══════════════════════════════════════════════════════════════════════════
# SAVE CSVs
# ═══════════════════════════════════════════════════════════════════════════
print("\nSaving CSVs...")

sig_fields = ['sig_idx', 'cmap_name', 'cell_iname', 'moa', 'target',
              'dose_mode', 'time_mode', 'n_reps', 'enrichment_score']
drug_fields = ['drug', 'median_score', 'mean_score', 'min_score', 'max_score', 'std_score',
               'n_signatures', 'n_cell_lines', 'frac_negative', 'moa', 'target',
               'pval_reversal', 'fdr_reversal', 'pval_mimicking', 'fdr_mimicking']

def write_csv(fn, data, fields):
    with open(OUT_DIR / fn, 'w', newline='') as f:
        csv.DictWriter(f, fieldnames=fields).writeheader()
        csv.DictWriter(f, fieldnames=fields).writerows(data)

# Fix: write header and rows with same writer
def write_csv(fn, data, fields):
    with open(OUT_DIR / fn, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(data)

write_csv('top50_reversing_signatures.csv', results_reverse, sig_fields)
write_csv('top50_mimicking_signatures.csv', results_mimic, sig_fields)
write_csv('drug_level_summary.csv', drug_summary, drug_fields)
write_csv('moa_enrichment.csv', moa_enrichment,
          ['moa', 'n_drugs_total', 'n_drugs_significant', 'fold_enrichment', 'pval', 'fdr'])
write_csv('cell_line_stats.csv', cell_stats,
          ['cell_line', 'n_signatures', 'mean_score', 'median_score', 'frac_negative'])

# All scores (large)
with open(OUT_DIR / 'all_signature_scores.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sig_idx', 'cmap_name', 'cell_iname', 'moa', 'enrichment_score'])
    for idx in sorted_idx:
        m = meta_dict[idx]
        w.writerow([idx, m['cmap_name'], m['cell_iname'], safe_str(m.get('moa')), f"{scores[idx]:.6f}"])

print(f"  CSVs saved to {OUT_DIR}")

# ═══════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════
plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'figure.dpi': 200, 'savefig.bbox': 'tight', 'axes.spines.top': False,
    'axes.spines.right': False,
})

q_scores = scores[quality_idx]
pct5 = np.percentile(q_scores, 5)
pct95 = np.percentile(q_scores, 95)

# ── PLOT 1: Score distribution ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 5))
ax.hist(q_scores, bins=120, color=C_BLUE, edgecolor='white', linewidth=0.2, alpha=0.85)
ax.axvline(0, color='black', linestyle='--', linewidth=0.8)
ax.axvline(pct5, color=C_REV, linestyle='--', linewidth=1.2,
           label=f'5th pctl ({pct5:.3f})')
ax.axvline(pct95, color=C_MIM, linestyle='--', linewidth=1.2,
           label=f'95th pctl ({pct95:.3f})')
ax.set_xlabel('Enrichment Score')
ax.set_ylabel('Number of Signatures')
ax.set_title(f'CMap Enrichment Score Distribution\n'
             f'{n_set} query genes | {len(quality_idx):,} signatures (n_reps>={MIN_REPS}) | {len(drug_summary):,} compounds | Method: {scoring_method}')
ax.legend(fontsize=9, framealpha=0.9)
ax.grid(axis='y', alpha=0.15)
fig.savefig(OUT_DIR / 'plot1_score_distribution.png', dpi=200)
plt.close()

# ── PLOT 2: Top 25 reversing SIGNATURES (dual) ──────────────────────────
def sig_label(idx):
    return f"{meta_dict[idx]['cmap_name']} ({meta_dict[idx]['cell_iname']})"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
for ax, idxs, title_suffix in [
    (ax1, reverse_idx[:25], "All Drugs"),
    (ax2, named_reverse_idx[:25], "Named Drugs (excl. BRD-)"),
]:
    labels = [sig_label(i) for i in idxs]
    vals = [float(scores[i]) for i in idxs]
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, vals, color=C_REV, edgecolor='white', height=0.72, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Enrichment Score')
    ax.set_title(f'Top 25 Reversing — Signature-Level (Single Experiments) — {title_suffix}', fontweight='bold')
    ax.axvline(0, color='black', linewidth=0.5, zorder=2)
    ax.grid(axis='x', alpha=0.15, zorder=1)
    ax.set_facecolor('#FCFCFC')
fig.tight_layout(w_pad=3)
fig.savefig(OUT_DIR / 'plot2_top25_reversing_signature_level.png', dpi=200)
plt.close()

# ── PLOT 3: Top 25 mimicking SIGNATURES (dual) ──────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
for ax, idxs, title_suffix in [
    (ax1, mimic_idx[:25], "All Drugs"),
    (ax2, named_mimic_idx[:25], "Named Drugs (excl. BRD-)"),
]:
    labels = [sig_label(i) for i in idxs]
    vals = [float(scores[i]) for i in idxs]
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, vals, color=C_MIM, edgecolor='white', height=0.72, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Enrichment Score')
    ax.set_title(f'Top 25 Mimicking — Signature-Level (Single Experiments) — {title_suffix}', fontweight='bold')
    ax.axvline(0, color='black', linewidth=0.5, zorder=2)
    ax.grid(axis='x', alpha=0.15, zorder=1)
    ax.set_facecolor('#FCFCFC')
fig.tight_layout(w_pad=3)
fig.savefig(OUT_DIR / 'plot3_top25_mimicking_signature_level.png', dpi=200)
plt.close()

# ── PLOT 4: Top 25 reversing DRUGS (dual) ────────────────────────────────
make_dual_barh(
    all_drug_summary_rev, named_drug_summary_rev, 'median_score',
    lambda d: f"{d['drug']} (n={d['n_cell_lines']})", C_REV, C_REV,
    'Top 25 Reversing — Drug-Level (Median Across Cell Lines) — All',
    'Top 25 Reversing — Drug-Level (Median Across Cell Lines) — Named (excl. BRD-)',
    'Median Enrichment Score', 'plot4_top25_reversing_drug_level.png',
    fdr_key='fdr_reversal'
)

# ── PLOT 5: Top 25 mimicking DRUGS (dual) ────────────────────────────────
make_dual_barh(
    all_drug_summary_mim, named_drug_summary_mim, 'median_score',
    lambda d: f"{d['drug']} (n={d['n_cell_lines']})", C_MIM, C_MIM,
    'Top 25 Mimicking — Drug-Level (Median Across Cell Lines) — All',
    'Top 25 Mimicking — Drug-Level (Median Across Cell Lines) — Named (excl. BRD-)',
    'Median Enrichment Score', 'plot5_top25_mimicking_drug_level.png',
    fdr_key='fdr_mimicking'
)

# ── PLOT 6: MOA enrichment (dual: top 5% threshold, top 10% threshold) ──
sig_moas_5 = [m for m in moa_enrichment if m['fdr'] < 0.25][:20]  # relaxed for display
if sig_moas_5:
    fig, ax = plt.subplots(figsize=(11, max(4, len(sig_moas_5) * 0.38 + 1)))
    labels = [f"{m['moa'][:45]} ({m['n_drugs_significant']}/{m['n_drugs_total']})" for m in sig_moas_5]
    vals = [m['fold_enrichment'] for m in sig_moas_5]
    bar_colors = [C_PURPLE if m['fdr'] < 0.05 else C_GREY for m in sig_moas_5]
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, vals, color=bar_colors, edgecolor='white', height=0.72, zorder=3)
    for i, m in enumerate(sig_moas_5):
        ax.text(vals[i] + 0.05, i, f"FDR={m['fdr']:.3f}", va='center', fontsize=7, color='#333')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Fold Enrichment')
    ax.set_title(f'MOA Enrichment Among Top {top_pct*100:.0f}% Reversing Drugs\n'
                 f'(Purple = FDR<0.05, Grey = FDR 0.05-0.25)', fontweight='bold')
    ax.axvline(1, color='black', linestyle='--', linewidth=0.5)
    ax.grid(axis='x', alpha=0.15, zorder=1)
    ax.set_facecolor('#FCFCFC')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'plot6_moa_enrichment.png', dpi=200)
    plt.close()
else:
    print("  No MOAs enriched at FDR<0.25; skipping plot6")

# ── PLOT 7: Cell line comparison ─────────────────────────────────────────
cell_stats_top = sorted(cell_stats, key=lambda x: -x['n_signatures'])[:20]
cell_stats_top.sort(key=lambda x: x['median_score'])

fig, ax = plt.subplots(figsize=(11, 6))
labels = [f"{c['cell_line']} (n={c['n_signatures']:,})" for c in cell_stats_top]
vals = [c['median_score'] for c in cell_stats_top]
colors = [C_REV if v < 0 else C_MIM for v in vals]
y_pos = np.arange(len(labels))
ax.barh(y_pos, vals, color=colors, edgecolor='white', height=0.72, zorder=3)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('Median Enrichment Score')
ax.set_title('Median Enrichment Score by Cell Line (top 20 by coverage)', fontweight='bold')
ax.axvline(0, color='black', linewidth=0.5, zorder=2)
ax.grid(axis='x', alpha=0.15, zorder=1)
ax.set_facecolor('#FCFCFC')
fig.tight_layout()
fig.savefig(OUT_DIR / 'plot7_cellline_comparison.png', dpi=200)
plt.close()

# ── PLOT 8: Drug x Cell Line heatmap (dual) ─────────────────────────────
top_cells_list = [c['cell_line'] for c in sorted(cell_stats, key=lambda x: -x['n_signatures'])[:10]]

def make_heatmap(drug_list, cells, title, filename):
    if not drug_list or not cells:
        return
    drug_cell_map = defaultdict(dict)
    drug_set = set(drug_list)
    cell_set = set(cells)
    for idx in range(n_sigs):
        m = meta_dict[idx]
        if m['cmap_name'] in drug_set and m['cell_iname'] in cell_set:
            drug_cell_map[m['cmap_name']][m['cell_iname']] = float(scores[idx])

    mat = np.full((len(drug_list), len(cells)), np.nan)
    for i, drug in enumerate(drug_list):
        for j, cell in enumerate(cells):
            if cell in drug_cell_map.get(drug, {}):
                mat[i, j] = drug_cell_map[drug][cell]

    valid = mat[~np.isnan(mat)]
    if len(valid) == 0:
        return
    vmax = np.percentile(np.abs(valid), 95)

    fig, ax = plt.subplots(figsize=(12, max(5, len(drug_list) * 0.4 + 1)))
    masked = np.ma.masked_invalid(mat)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad('#F0F0F0')
    im = ax.imshow(masked, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_xticks(range(len(cells)))
    ax.set_xticklabels(cells, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(drug_list)))
    ax.set_yticklabels(drug_list, fontsize=9)
    ax.set_title(title, fontweight='bold', fontsize=12)
    plt.colorbar(im, ax=ax, label='Enrichment Score', shrink=0.8, pad=0.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, dpi=200)
    plt.close()

all_rev_drugs_heatmap = [d['drug'] for d in drug_summary if d['n_cell_lines'] >= 3][:20]
named_rev_drugs_heatmap = [d['drug'] for d in named_drug_summary if d['n_cell_lines'] >= 3][:20]

make_heatmap(all_rev_drugs_heatmap, top_cells_list,
             'Top 20 Reversing Drugs x Top 10 Cell Lines — All',
             'plot8a_heatmap_all.png')
make_heatmap(named_rev_drugs_heatmap, top_cells_list,
             'Top 20 Reversing Drugs x Top 10 Cell Lines — Named (excl. BRD-)',
             'plot8b_heatmap_named.png')

# ── PLOT 9: Z-score heatmap (dual) ──────────────────────────────────────
def make_zscore_heatmap(sig_indices, title, filename, n=15):
    idxs = sig_indices[:n]
    if not len(idxs):
        return
    zmat = np.zeros((len(idxs), len(gene_idx_list)))
    for i, si in enumerate(idxs):
        zmat[i] = zscore_matrix[si, gene_idx_list]

    fig, ax = plt.subplots(figsize=(18, max(4, len(idxs) * 0.45 + 1)))
    vmax = np.percentile(np.abs(zmat), 95)
    im = ax.imshow(zmat, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_yticks(range(len(idxs)))
    ax.set_yticklabels([sig_label(s) for s in idxs], fontsize=8)
    tick_interval = max(1, len(mapped_genes) // 25)
    xtick_pos = list(range(0, len(mapped_genes), tick_interval))
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels([mapped_genes[i] for i in xtick_pos], fontsize=6, rotation=90)
    ax.set_title(title, fontweight='bold', fontsize=12)
    plt.colorbar(im, ax=ax, label='Z-score', shrink=0.7, pad=0.01)
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, dpi=200)
    plt.close()

make_zscore_heatmap(list(reverse_idx), 'Query Gene Z-scores in Top 15 Reversing Signatures — All',
                    'plot9a_zscore_heatmap_all.png')
make_zscore_heatmap(named_reverse_idx, 'Query Gene Z-scores in Top 15 Reversing Signatures — Named (excl. BRD-)',
                    'plot9b_zscore_heatmap_named.png')

# ── PLOT 10: Volcano (dual) ─────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

multi_drugs = [d for d in drug_summary if d['n_cell_lines'] >= 2]

for ax, title_suffix, highlight_only_named in [
    (ax1, "All Drugs", False),
    (ax2, "Named Drugs (excl. BRD-)", True),
]:
    x_vals, y_vals, c_vals, sizes = [], [], [], []
    for d in multi_drugs:
        x = d['median_score']
        fdr = d['fdr_reversal'] if x < 0 else d['fdr_mimicking']
        y = -np.log10(max(fdr, 1e-10))
        x_vals.append(x)
        y_vals.append(y)

        if highlight_only_named and not is_named_drug(d['drug']):
            c_vals.append('#EEEEEE')
            sizes.append(6)
        elif x < -0.10:
            c_vals.append(C_REV)
            sizes.append(14)
        elif x > 0.10:
            c_vals.append(C_MIM)
            sizes.append(14)
        else:
            c_vals.append(C_GREY)
            sizes.append(8)

    ax.scatter(x_vals, y_vals, s=sizes, c=c_vals, alpha=0.6, edgecolors='none', zorder=3)
    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.8, label='FDR=0.05')
    ax.axvline(0, color='black', linestyle='--', linewidth=0.5)

    # Label top named drugs
    label_src = named_drug_summary_rev[:10] if highlight_only_named else all_drug_summary_rev[:10]
    for d in label_src:
        fdr = max(d['fdr_reversal'], 1e-10)
        ax.annotate(d['drug'], (d['median_score'], -np.log10(fdr)),
                    fontsize=6, color=C_REV, fontweight='bold',
                    xytext=(5, 3), textcoords='offset points')
    label_src_m = named_drug_summary_mim[:10] if highlight_only_named else all_drug_summary_mim[:10]
    for d in label_src_m:
        fdr = max(d['fdr_mimicking'], 1e-10)
        ax.annotate(d['drug'], (d['median_score'], -np.log10(fdr)),
                    fontsize=6, color=C_MIM, fontweight='bold',
                    xytext=(5, 3), textcoords='offset points')

    ax.set_xlabel('Median Enrichment Score')
    ax.set_ylabel('-log10(FDR)')
    ax.set_title(f'Volcano — {title_suffix}', fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.1, zorder=1)
    ax.set_facecolor('#FCFCFC')

fig.tight_layout(w_pad=3)
fig.savefig(OUT_DIR / 'plot10_volcano.png', dpi=200)
plt.close()

# ── PLOT 11: Drug–Gene Network (dual) ───────────────────────────────────
def make_network(sig_indices, title, filename, n_drugs=25):
    """Build and draw a drug-gene bipartite network."""
    G = nx.Graph()
    drug_nodes = set()
    gene_nodes = set()
    edges = []

    for idx in sig_indices[:n_drugs]:
        m = meta_dict[idx]
        drug_label = f"{m['cmap_name']}\n({m['cell_iname']})"
        drug_nodes.add(drug_label)
        zscores_vec = zscore_matrix[idx, gene_idx_list]

        for j, gene in enumerate(mapped_genes):
            z = float(zscores_vec[j])
            if abs(z) >= ZSCORE_THRESHOLD:
                gene_nodes.add(gene)
                edges.append((drug_label, gene, z))

    for d in drug_nodes:
        G.add_node(d, node_type='drug')
    for g in gene_nodes:
        G.add_node(g, node_type='gene')
    for d, g, z in edges:
        if G.has_edge(d, g):
            if abs(z) > abs(G[d][g]['weight']):
                G[d][g]['weight'] = z
        else:
            G.add_edge(d, g, weight=z)

    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    if G.number_of_nodes() < 3:
        print(f"  Network too sparse for {filename}; skipping")
        return

    pos = nx.spring_layout(G, k=2.5, iterations=150, seed=42)

    fig, ax = plt.subplots(figsize=(20, 16))
    ax.set_facecolor(C_BG)
    fig.patch.set_facecolor(C_BG)

    drugs_g = [n for n in G.nodes if G.nodes[n]['node_type'] == 'drug']
    genes_g = [n for n in G.nodes if G.nodes[n]['node_type'] == 'gene']

    drug_sizes = [400 + G.degree(n) * 60 for n in drugs_g]
    gene_deg = [G.degree(n) for n in genes_g]
    max_gd = max(gene_deg) if gene_deg else 1
    gene_sizes = [120 + (d / max_gd) * 400 for d in gene_deg]

    edge_list = list(G.edges(data=True))
    edge_w = [abs(e[2]['weight']) for e in edge_list]
    edge_c = ['#2166AC' if e[2]['weight'] < 0 else '#B2182B' for e in edge_list]
    max_ew = max(edge_w) if edge_w else 1
    edge_widths = [0.5 + (w / max_ew) * 2.5 for w in edge_w]

    nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1]) for e in edge_list],
                           edge_color=edge_c, width=edge_widths, alpha=0.35, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=drugs_g, node_color=C_REV,
                           node_size=drug_sizes, edgecolors='#8B0000',
                           linewidths=1.0, alpha=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=genes_g, node_color=C_BLUE,
                           node_size=gene_sizes, edgecolors='#1A3A5C',
                           linewidths=0.7, alpha=0.85, ax=ax)

    drug_lpos = {n: (pos[n][0], pos[n][1] + 0.04) for n in drugs_g}
    gene_lpos = {n: (pos[n][0], pos[n][1] - 0.03) for n in genes_g}

    nx.draw_networkx_labels(G, drug_lpos,
                            labels={n: n.replace('\n', ' ') for n in drugs_g},
                            font_size=7, font_weight='bold', font_color='#4A0000', ax=ax)
    nx.draw_networkx_labels(G, gene_lpos,
                            labels={n: n for n in genes_g},
                            font_size=7, font_weight='bold', font_color='#0A2A4A', ax=ax)

    legend_el = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_REV,
               markersize=12, markeredgecolor='#8B0000', label='Drug'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_BLUE,
               markersize=10, markeredgecolor='#1A3A5C', label='Input Gene'),
        Line2D([0], [0], color='#2166AC', linewidth=2, alpha=0.6, label='Downreg (z<0)'),
        Line2D([0], [0], color='#B2182B', linewidth=2, alpha=0.6, label='Upreg (z>0)'),
    ]
    ax.legend(handles=legend_el, loc='upper left', fontsize=9, framealpha=0.9, edgecolor='gray')
    ax.set_title(f'{title}\n(edges: |z-score| > {ZSCORE_THRESHOLD})',
                 fontsize=14, fontweight='bold')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(OUT_DIR / filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Network saved: {filename} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")

print("\nGenerating network plots...")
make_network(list(reverse_idx), 'Drug-Gene Network: Top 25 Reversing — All',
             'plot11a_network_all.png')
make_network(named_reverse_idx, 'Drug-Gene Network: Top 25 Reversing — Named (excl. BRD-)',
             'plot11b_network_named.png')

# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
with open(OUT_DIR / 'summary.txt', 'w') as f:
    f.write("CMap Enrichment Analysis v3\n")
    f.write("=" * 65 + "\n\n")
    f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Gene file: {GENE_FILE.name}\n")
    f.write(f"Scoring method: {scoring_method}\n")
    f.write(f"Database: {n_sigs:,} total signatures, {len(quality_idx):,} pass quality filter "
            f"(n_reps>={MIN_REPS})\n")
    f.write(f"  {len(drug_summary):,} compounds, {len(cell_stats)} cell lines\n\n")
    f.write(f"Input genes: {len(raw_genes)}\n")
    f.write(f"Mapped to L1000: {n_set}\n")
    if has_values:
        f.write(f"  Upregulated: {int(sum(mapped_values >= 0))}\n")
        f.write(f"  Downregulated: {int(sum(mapped_values < 0))}\n")
    f.write(f"Unmapped: {len(unmapped_genes)}")
    if unmapped_genes:
        f.write(f" ({', '.join(unmapped_genes)})")
    f.write("\n\n")
    f.write(f"Score statistics (quality-filtered):\n")
    f.write(f"  Range: [{q_scores.min():.4f}, {q_scores.max():.4f}]\n")
    f.write(f"  Mean: {q_scores.mean():.4f}, Std: {q_scores.std():.4f}\n")
    f.write(f"  5th pctl: {pct5:.4f}, 95th: {pct95:.4f}\n\n")
    f.write(f"Drug-level results ({N_PERM:,} permutations):\n")
    f.write(f"  Significant reversals (FDR<0.05): {n_sig_rev:,}\n")
    f.write(f"  Significant mimickers (FDR<0.05): {n_sig_mim:,}\n")
    f.write(f"  Named drugs: {n_named:,} / BRD- coded: {n_coded:,}\n\n")

    f.write("Top 15 Reversing Named Drugs (median score, n_cell_lines >= 2):\n")
    f.write("-" * 65 + "\n")
    for d in named_drug_summary_rev[:15]:
        fdr_s = f"{d['fdr_reversal']:.4f}" if d['fdr_reversal'] >= 0.0001 else "<0.0001"
        f.write(f"  {d['drug']:25s} median={d['median_score']:.4f}  "
                f"n_cells={d['n_cell_lines']:2d}  FDR={fdr_s}")
        if safe_str(d['moa']):
            f.write(f"  MOA: {d['moa'][:40]}")
        f.write("\n")

    f.write("\nTop 15 Mimicking Named Drugs (median score, n_cell_lines >= 2):\n")
    f.write("-" * 65 + "\n")
    for d in named_drug_summary_mim[:15]:
        fdr_s = f"{d['fdr_mimicking']:.4f}" if d['fdr_mimicking'] >= 0.0001 else "<0.0001"
        f.write(f"  {d['drug']:25s} median={d['median_score']:+.4f}  "
                f"n_cells={d['n_cell_lines']:2d}  FDR={fdr_s}")
        if safe_str(d['moa']):
            f.write(f"  MOA: {d['moa'][:40]}")
        f.write("\n")

    sig_moas_txt = [m for m in moa_enrichment if m['fdr'] < 0.25][:15]
    if sig_moas_txt:
        f.write(f"\nMOA Enrichment (top {top_pct*100:.0f}% reversing drugs, FDR<0.25):\n")
        f.write("-" * 65 + "\n")
        for m in sig_moas_txt:
            f.write(f"  {m['moa'][:45]:45s} fold={m['fold_enrichment']:.2f}  "
                    f"({m['n_drugs_significant']}/{m['n_drugs_total']})  FDR={m['fdr']:.4f}\n")

total_time = time.time() - t0
print(f"\nAnalysis complete in {total_time/60:.1f} minutes")
print(f"Outputs: {OUT_DIR}")
print("\nFiles:")
for p in sorted(OUT_DIR.glob('*')):
    if p.suffix in ('.py', '.sh'):
        continue
    if p.name.startswith('slurm_'):
        continue
    sz = p.stat().st_size
    if sz > 1e6:
        print(f"  {p.name} ({sz/1e6:.1f} MB)")
    elif sz > 1e3:
        print(f"  {p.name} ({sz/1e3:.1f} KB)")
    else:
        print(f"  {p.name} ({sz} B)")
