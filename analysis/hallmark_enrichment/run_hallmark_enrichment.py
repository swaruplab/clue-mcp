#!/usr/bin/env python3
"""
MSigDB Hallmark Enrichment Analysis

Downloads all 50 Hallmark gene sets from MSigDB, runs CMap enrichment
for each hallmark separately, aggregates results, and produces:
  - Per-hallmark enrichment CSVs
  - Cross-hallmark drug score matrix
  - Master heatmap (top 5 drugs x 50 hallmarks)
  - Biology-focused network plot
  - Summary report
"""

import sys, os, json, csv, time, urllib.request, ssl
import numpy as np
import duckdb
from pathlib import Path
from collections import defaultdict
from scipy.stats import rankdata

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import networkx as nx

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════
BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed"
MSIGDB_DIR = BASE_DIR / "data" / "msigdb"
OUT_DIR = Path(__file__).resolve().parent
PER_HALLMARK_DIR = OUT_DIR / "per_hallmark"
AGG_DIR = OUT_DIR / "aggregated"
PLOT_DIR = OUT_DIR / "plots"

MSIGDB_URL = "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.json"
MSIGDB_JSON = MSIGDB_DIR / "h.all.v2026.1.Hs.json"

N_PERM = 10000
MIN_REPS = 2
TOP_N = 50
TOP_DRUGS = 5  # for master heatmap / network

# Color palette
C_REV = '#D7263D'
C_MIM = '#1B998B'
C_BLUE = '#3A86FF'
C_PURPLE = '#7B2D8E'
C_GOLD = '#F4A261'
C_GREY = '#BFBFBF'
C_BG = '#FAFAFA'

# Biological grouping of hallmark gene sets
HALLMARK_CATEGORIES = {
    'Immune': [
        'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB',
        'HALLMARK_IL_6_JAK_STAT3_SIGNALING', 'HALLMARK_IL_2_STAT5_SIGNALING',
        'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
        'HALLMARK_COMPLEMENT', 'HALLMARK_ALLOGRAFT_REJECTION',
        'HALLMARK_COAGULATION',
    ],
    'Proliferation': [
        'HALLMARK_E2F_TARGETS', 'HALLMARK_G2_M_CHECKPOINT', 'HALLMARK_MITOTIC_SPINDLE',
        'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2',
    ],
    'Signaling': [
        'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_MTORC1_SIGNALING',
        'HALLMARK_HEDGEHOG_SIGNALING', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
        'HALLMARK_NOTCH_SIGNALING', 'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_KRAS_SIGNALING_UP', 'HALLMARK_KRAS_SIGNALING_DN',
    ],
    'Metabolic': [
        'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_FATTY_ACID_METABOLISM',
        'HALLMARK_GLYCOLYSIS', 'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
        'HALLMARK_BILE_ACID_METABOLISM', 'HALLMARK_ADIPOGENESIS',
        'HALLMARK_XENOBIOTIC_METABOLISM', 'HALLMARK_PEROXISOME',
        'HALLMARK_HEME_METABOLISM', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
    ],
    'Stress & DNA Damage': [
        'HALLMARK_P53_PATHWAY', 'HALLMARK_APOPTOSIS', 'HALLMARK_HYPOXIA',
        'HALLMARK_UV_RESPONSE_UP', 'HALLMARK_UV_RESPONSE_DN',
        'HALLMARK_DNA_REPAIR', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
    ],
    'Development': [
        'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_ANGIOGENESIS',
        'HALLMARK_MYOGENESIS', 'HALLMARK_SPERMATOGENESIS', 'HALLMARK_PANCREAS_BETA_CELLS',
        'HALLMARK_APICAL_JUNCTION', 'HALLMARK_APICAL_SURFACE',
    ],
    'Other': [
        'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_ESTROGEN_RESPONSE_LATE',
        'HALLMARK_ANDROGEN_RESPONSE', 'HALLMARK_PROTEIN_SECRETION',
    ],
}

# Reverse lookup: hallmark -> category
HALLMARK_TO_CATEGORY = {}
for cat, members in HALLMARK_CATEGORIES.items():
    for h in members:
        HALLMARK_TO_CATEGORY[h] = cat

CATEGORY_COLORS = {
    'Immune': '#E63946',
    'Proliferation': '#457B9D',
    'Signaling': '#2A9D8F',
    'Metabolic': '#E9C46A',
    'Stress & DNA Damage': '#F4A261',
    'Development': '#264653',
    'Other': '#BFBFBF',
}


def safe_str(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ''
    return str(val)


def is_named_drug(name):
    return not name.startswith('BRD-')


def short_hallmark_name(name):
    """HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB -> TNFa Signaling via NFkB"""
    s = name.replace('HALLMARK_', '')
    s = s.replace('_', ' ').title()
    # Fix common abbreviations after title-casing
    replacements = [
        ('Tnf Alpha', 'TNFa'), ('Nf Kb', 'NFkB'),
        ('Il 6', 'IL-6'), ('Il 2', 'IL-2'),
        ('Jak', 'JAK'), ('Stat3', 'STAT3'), ('Stat5', 'STAT5'),
        ('Mtor', 'mTOR'), ('Mtorc1', 'mTORC1'),
        ('Pi3K', 'PI3K'), ('Akt', 'AKT'), ('Kras', 'KRAS'),
        ('E2F', 'E2F'), ('G2 M', 'G2/M'), ('Myc', 'MYC'),
        ('P53', 'p53'), ('Dna', 'DNA'), ('Uv', 'UV'),
        ('Tgf', 'TGF'), ('Wnt', 'WNT'),
    ]
    for old, new in replacements:
        s = s.replace(old, new)
    return s


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction on a list of p-values."""
    n = len(pvals)
    if n == 0:
        return []
    ranked = np.argsort(pvals)
    fdr = np.ones(n)
    for i, idx in enumerate(ranked):
        fdr[idx] = pvals[idx] * n / (i + 1)
    # Enforce monotonicity
    fdr = np.minimum(fdr, 1.0)
    prev = 1.0
    for idx in reversed(ranked):
        fdr[idx] = min(fdr[idx], prev)
        prev = fdr[idx]
    return fdr


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: DOWNLOAD MSigDB HALLMARK GENE SETS
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("MSigDB Hallmark Enrichment Analysis")
print("=" * 65)
t0 = time.time()

if MSIGDB_JSON.exists():
    print(f"\nMSigDB JSON already exists: {MSIGDB_JSON}")
else:
    print(f"\nDownloading MSigDB Hallmark gene sets...")
    MSIGDB_DIR.mkdir(parents=True, exist_ok=True)
    # Try direct download; some servers need SSL context relaxation
    try:
        ctx = ssl.create_default_context()
        req = urllib.request.Request(MSIGDB_URL, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, context=ctx, timeout=120) as resp:
            data = resp.read()
        with open(MSIGDB_JSON, 'wb') as f:
            f.write(data)
        print(f"  Downloaded {len(data)} bytes -> {MSIGDB_JSON}")
    except Exception as e:
        print(f"  Direct download failed: {e}")
        print("  Trying with relaxed SSL...")
        try:
            ctx = ssl._create_unverified_context()
            req = urllib.request.Request(MSIGDB_URL, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, context=ctx, timeout=120) as resp:
                data = resp.read()
            with open(MSIGDB_JSON, 'wb') as f:
                f.write(data)
            print(f"  Downloaded {len(data)} bytes -> {MSIGDB_JSON}")
        except Exception as e2:
            print(f"  ERROR: Could not download MSigDB JSON: {e2}")
            print("  Please download manually from:")
            print(f"  {MSIGDB_URL}")
            print(f"  Save to: {MSIGDB_JSON}")
            sys.exit(1)

# Parse JSON
with open(MSIGDB_JSON) as f:
    msigdb_data = json.load(f)

# MSigDB JSON format: list of objects or dict of gene set name -> info
# Handle both formats
if isinstance(msigdb_data, dict):
    hallmark_sets = {}
    for name, info in msigdb_data.items():
        if isinstance(info, dict):
            genes = info.get('geneSymbols', info.get('genes', []))
        elif isinstance(info, list):
            genes = info
        else:
            continue
        hallmark_sets[name] = genes
elif isinstance(msigdb_data, list):
    hallmark_sets = {}
    for entry in msigdb_data:
        if isinstance(entry, dict):
            name = entry.get('systematicName', entry.get('name', entry.get('geneset_name', '')))
            if not name:
                # Try other keys
                for k in entry:
                    if 'HALLMARK' in str(entry[k]):
                        name = entry[k]
                        break
            genes = entry.get('geneSymbols', entry.get('genes', []))
            if name and genes:
                hallmark_sets[name] = genes
else:
    print(f"ERROR: Unexpected JSON format (type: {type(msigdb_data)})")
    sys.exit(1)

print(f"\nLoaded {len(hallmark_sets)} hallmark gene sets:")
for name in sorted(hallmark_sets.keys()):
    print(f"  {name}: {len(hallmark_sets[name])} genes")

if len(hallmark_sets) == 0:
    print("ERROR: No gene sets found in MSigDB JSON. Check file format.")
    print(f"  JSON top-level type: {type(msigdb_data)}")
    if isinstance(msigdb_data, dict):
        print(f"  Keys (first 5): {list(msigdb_data.keys())[:5]}")
    elif isinstance(msigdb_data, list) and len(msigdb_data) > 0:
        print(f"  First entry keys: {list(msigdb_data[0].keys()) if isinstance(msigdb_data[0], dict) else type(msigdb_data[0])}")
    sys.exit(1)

# ═══════════════════════════════════════════════════════════════════════════
# LOAD CMAP DATA
# ═══════════════════════════════════════════════════════════════════════════
print("\nLoading CMap database...")
t1 = time.time()

rank_matrix = np.load(DATA_DIR / "rank_matrix.npy", mmap_mode='r')
zscore_matrix = np.load(DATA_DIR / "zscore_matrix.npy", mmap_mode='r')
with open(DATA_DIR / "gene_names.txt") as f:
    gene_names = [g.strip() for g in f.readlines()]
gene_to_idx = {g: i for i, g in enumerate(gene_names)}
n_sigs, n_genes = rank_matrix.shape

con = duckdb.connect(str(DATA_DIR / "cmap.duckdb"), read_only=True)
print(f"  {n_sigs:,} signatures x {n_genes:,} genes ({time.time()-t1:.1f}s)")

# Quality filter
quality_idx = con.execute(
    f"SELECT sig_idx FROM signatures WHERE n_reps >= {MIN_REPS} ORDER BY sig_idx"
).fetchnumpy()["sig_idx"]
print(f"  Quality filter (n_reps >= {MIN_REPS}): {len(quality_idx):,} / {n_sigs:,} signatures")

# Load all metadata
sig_meta_df = con.execute("SELECT * FROM signatures ORDER BY sig_idx").fetchdf()
meta_dict = {}
for _, row in sig_meta_df.iterrows():
    meta_dict[row['sig_idx']] = row.to_dict()
quality_set = set(quality_idx)


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


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: RUN ENRICHMENT FOR EACH HALLMARK
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("Running enrichment for each hallmark gene set...")
print("=" * 65)

hallmark_results = {}  # hallmark_name -> dict with scores, drug_summary, etc.

for h_idx, (hallmark_name, hallmark_genes) in enumerate(sorted(hallmark_sets.items())):
    t_h = time.time()
    print(f"\n[{h_idx+1}/{len(hallmark_sets)}] {hallmark_name}")

    # Map genes
    mapped = [g for g in hallmark_genes if g in gene_to_idx]
    unmapped = [g for g in hallmark_genes if g not in gene_to_idx]
    gene_idx_list = np.array([gene_to_idx[g] for g in mapped])
    n_set = len(gene_idx_list)

    print(f"  Genes: {len(hallmark_genes)} total, {n_set} mapped, {len(unmapped)} unmapped")

    if n_set < 5:
        print(f"  SKIPPING: too few mapped genes ({n_set})")
        continue

    # Compute enrichment scores (single gene set = KS)
    scores = vectorized_ks(rank_matrix, gene_idx_list)

    # Drug-level aggregation
    drug_scores_map = defaultdict(list)
    drug_cells = defaultdict(set)
    drug_moa = {}
    drug_target = {}

    for idx in quality_idx:
        m = meta_dict[idx]
        drug = m['cmap_name']
        drug_scores_map[drug].append(float(scores[idx]))
        drug_cells[drug].add(m['cell_iname'])
        moa_val = m.get('moa')
        if drug not in drug_moa and safe_str(moa_val):
            drug_moa[drug] = safe_str(moa_val)
            drug_target[drug] = safe_str(m.get('target', ''))

    drug_summary = []
    for drug in drug_scores_map:
        s = np.array(drug_scores_map[drug])
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

    # Permutation p-values for drug-level scores
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

    # BH FDR
    drugs_list = [d['drug'] for d in drug_summary]
    rev_pvals = np.array([drug_pvals[d]['pval_reversal'] for d in drugs_list])
    mim_pvals = np.array([drug_pvals[d]['pval_mimicking'] for d in drugs_list])
    rev_fdr = bh_fdr(rev_pvals)
    mim_fdr = bh_fdr(mim_pvals)

    for i, d in enumerate(drug_summary):
        d['pval_reversal'] = drug_pvals[d['drug']]['pval_reversal']
        d['fdr_reversal'] = float(rev_fdr[i])
        d['pval_mimicking'] = drug_pvals[d['drug']]['pval_mimicking']
        d['fdr_mimicking'] = float(mim_fdr[i])

    # Save per-hallmark results
    h_dir = PER_HALLMARK_DIR / hallmark_name
    h_dir.mkdir(parents=True, exist_ok=True)

    drug_fields = ['drug', 'median_score', 'mean_score', 'min_score', 'max_score', 'std_score',
                   'n_signatures', 'n_cell_lines', 'frac_negative', 'moa', 'target',
                   'pval_reversal', 'fdr_reversal', 'pval_mimicking', 'fdr_mimicking']

    with open(h_dir / 'drug_level_summary.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=drug_fields)
        w.writeheader()
        w.writerows(drug_summary)

    # Top signatures
    quality_scores = scores[quality_idx]
    quality_order = np.argsort(quality_scores)
    sorted_q_idx = quality_idx[quality_order]

    sig_fields = ['sig_idx', 'cmap_name', 'cell_iname', 'moa', 'target', 'enrichment_score']
    for fname, idxs in [('top50_reversing.csv', sorted_q_idx[:TOP_N]),
                         ('top50_mimicking.csv', sorted_q_idx[-TOP_N:][::-1])]:
        rows = []
        for idx in idxs:
            m = meta_dict[idx]
            rows.append({
                'sig_idx': int(idx), 'cmap_name': m['cmap_name'],
                'cell_iname': m['cell_iname'], 'moa': safe_str(m.get('moa')),
                'target': safe_str(m.get('target')), 'enrichment_score': f"{scores[idx]:.6f}"
            })
        with open(h_dir / fname, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=sig_fields)
            w.writeheader()
            w.writerows(rows)

    # Store results for aggregation
    hallmark_results[hallmark_name] = {
        'n_genes_input': len(hallmark_genes),
        'n_genes_mapped': n_set,
        'n_unmapped': len(unmapped),
        'drug_summary': drug_summary,
        'score_mean': float(np.mean(quality_scores)),
        'score_std': float(np.std(quality_scores)),
    }

    n_sig_rev = sum(1 for d in drug_summary if d['fdr_reversal'] < 0.05)
    n_sig_mim = sum(1 for d in drug_summary if d['fdr_mimicking'] < 0.05)
    elapsed = time.time() - t_h
    print(f"  Scores: mean={np.mean(quality_scores):.4f}, std={np.std(quality_scores):.4f}")
    print(f"  Drugs: {len(drug_summary):,} total, {n_sig_rev} sig reversals, {n_sig_mim} sig mimickers")
    print(f"  Done in {elapsed:.1f}s")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: AGGREGATE ACROSS ALL HALLMARKS
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("Aggregating results across all hallmarks...")
print("=" * 65)

hallmark_names_sorted = sorted(hallmark_results.keys())
n_hallmarks = len(hallmark_names_sorted)

# Build drug -> hallmark score matrix
# For each hallmark, get the drug-level median score (named drugs, n_cell_lines >= 2)
all_drugs = set()
for h_name in hallmark_names_sorted:
    for d in hallmark_results[h_name]['drug_summary']:
        if is_named_drug(d['drug']) and d['n_cell_lines'] >= 2:
            all_drugs.add(d['drug'])

all_drugs_sorted = sorted(all_drugs)
drug_to_row = {d: i for i, d in enumerate(all_drugs_sorted)}

# Build the matrix
score_matrix = np.full((len(all_drugs_sorted), n_hallmarks), np.nan, dtype=np.float32)
fdr_matrix = np.full((len(all_drugs_sorted), n_hallmarks), 1.0, dtype=np.float32)

for j, h_name in enumerate(hallmark_names_sorted):
    for d in hallmark_results[h_name]['drug_summary']:
        if d['drug'] in drug_to_row:
            score_matrix[drug_to_row[d['drug']], j] = d['median_score']
            fdr_matrix[drug_to_row[d['drug']], j] = d['fdr_reversal']

# Composite score for ranking: combine effect size (mean median across hallmarks)
# with breadth (number of hallmarks with FDR < 0.05 reversal)
drug_composite = []
for i, drug in enumerate(all_drugs_sorted):
    row = score_matrix[i]
    valid = row[~np.isnan(row)]
    if len(valid) == 0:
        continue
    mean_score = float(np.nanmean(row))
    n_sig_hallmarks = int(np.sum(fdr_matrix[i] < 0.05))
    # Composite: weighted combination of mean reversal score and breadth
    # More negative mean = better reverser; more hallmarks significant = broader effect
    composite = mean_score - 0.01 * n_sig_hallmarks  # lower is better
    drug_composite.append({
        'drug': drug,
        'mean_median_score': mean_score,
        'n_sig_hallmarks': n_sig_hallmarks,
        'composite_score': composite,
        'moa': '',
        'target': '',
    })

# Get MOA/target info from any hallmark's drug_summary
for dc in drug_composite:
    for h_name in hallmark_names_sorted:
        for d in hallmark_results[h_name]['drug_summary']:
            if d['drug'] == dc['drug'] and safe_str(d['moa']):
                dc['moa'] = d['moa']
                dc['target'] = d['target']
                break
        if dc['moa']:
            break

drug_composite.sort(key=lambda x: x['composite_score'])

# Top 5 drugs
top5 = drug_composite[:TOP_DRUGS]
print(f"\nTop {TOP_DRUGS} reversing drugs (composite score):")
for i, d in enumerate(top5):
    print(f"  {i+1}. {d['drug']} (mean_median={d['mean_median_score']:.4f}, "
          f"n_sig_hallmarks={d['n_sig_hallmarks']}, MOA: {d['moa'][:50]})")

# Save aggregated CSVs
AGG_DIR.mkdir(parents=True, exist_ok=True)

# Cross-hallmark matrix (full)
with open(AGG_DIR / 'cross_hallmark_drug_matrix.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['drug'] + [short_hallmark_name(h) for h in hallmark_names_sorted])
    for i, drug in enumerate(all_drugs_sorted):
        row = [drug] + [f"{score_matrix[i, j]:.6f}" if not np.isnan(score_matrix[i, j]) else ''
               for j in range(n_hallmarks)]
        w.writerow(row)

# Top 5 drugs summary
with open(AGG_DIR / 'top5_drugs_summary.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['rank', 'drug', 'mean_median_score', 'n_sig_hallmarks',
                                       'composite_score', 'moa', 'target'])
    w.writeheader()
    for i, d in enumerate(top5):
        w.writerow({'rank': i+1, **d})

# Hallmark statistics
with open(AGG_DIR / 'hallmark_statistics.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['hallmark', 'short_name', 'category',
                                       'n_genes_input', 'n_genes_mapped', 'n_unmapped',
                                       'score_mean', 'score_std'])
    w.writeheader()
    for h_name in hallmark_names_sorted:
        hr = hallmark_results[h_name]
        w.writerow({
            'hallmark': h_name,
            'short_name': short_hallmark_name(h_name),
            'category': HALLMARK_TO_CATEGORY.get(h_name, 'Other'),
            'n_genes_input': hr['n_genes_input'],
            'n_genes_mapped': hr['n_genes_mapped'],
            'n_unmapped': hr['n_unmapped'],
            'score_mean': f"{hr['score_mean']:.6f}",
            'score_std': f"{hr['score_std']:.6f}",
        })

print(f"\nAggregated data saved to {AGG_DIR}")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: MASTER HEATMAP
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("Generating master heatmap...")
print("=" * 65)

PLOT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 13, 'axes.labelsize': 11,
    'figure.dpi': 200, 'savefig.bbox': 'tight', 'axes.spines.top': False,
    'axes.spines.right': False,
})

# Order hallmarks by biological category
ordered_hallmarks = []
category_boundaries = []
category_labels = []
for cat in ['Immune', 'Proliferation', 'Signaling', 'Metabolic', 'Stress & DNA Damage', 'Development', 'Other']:
    members = [h for h in HALLMARK_CATEGORIES.get(cat, []) if h in hallmark_results]
    # Add any hallmarks not in our predefined categories
    if cat == 'Other':
        categorized = set()
        for c, ms in HALLMARK_CATEGORIES.items():
            categorized.update(ms)
        uncategorized = [h for h in hallmark_names_sorted if h not in categorized and h in hallmark_results]
        members = members + uncategorized
    if members:
        start = len(ordered_hallmarks)
        ordered_hallmarks.extend(members)
        category_boundaries.append((start, start + len(members)))
        category_labels.append(cat)

hallmark_col_idx = {h: i for i, h in enumerate(ordered_hallmarks)}

# Extract top5 drug rows in the biologically-ordered hallmark columns
top5_names = [d['drug'] for d in top5]
top5_moas = [d['moa'][:40] if d['moa'] else 'Unknown' for d in top5]

heatmap_data = np.full((len(top5_names), len(ordered_hallmarks)), np.nan)
heatmap_fdr = np.full((len(top5_names), len(ordered_hallmarks)), 1.0)

for i, drug in enumerate(top5_names):
    drug_row = drug_to_row[drug]
    for j_orig, h_name in enumerate(hallmark_names_sorted):
        if h_name in hallmark_col_idx:
            j_new = hallmark_col_idx[h_name]
            heatmap_data[i, j_new] = score_matrix[drug_row, j_orig]
            heatmap_fdr[i, j_new] = fdr_matrix[drug_row, j_orig]

# Create heatmap
fig = plt.figure(figsize=(24, 8))
gs = gridspec.GridSpec(2, 2, height_ratios=[0.3, 10], width_ratios=[10, 0.5],
                       hspace=0.02, wspace=0.03)

# Category color bar (top)
ax_cat = fig.add_subplot(gs[0, 0])
for start, end in category_boundaries:
    cat = category_labels[category_boundaries.index((start, end))]
    color = CATEGORY_COLORS.get(cat, '#BFBFBF')
    ax_cat.barh(0, end - start, left=start, color=color, edgecolor='white', linewidth=0.5, height=0.8)
    mid = (start + end) / 2
    if (end - start) >= 3:
        ax_cat.text(mid, 0, cat, ha='center', va='center', fontsize=7, fontweight='bold')
ax_cat.set_xlim(-0.5, len(ordered_hallmarks) - 0.5)
ax_cat.set_ylim(-0.5, 0.5)
ax_cat.axis('off')

# Main heatmap
ax_heat = fig.add_subplot(gs[1, 0])
valid_vals = heatmap_data[~np.isnan(heatmap_data)]
if len(valid_vals) > 0:
    vmax = max(abs(np.nanmin(heatmap_data)), abs(np.nanmax(heatmap_data)))
    vmax = max(vmax, 0.01)  # Prevent zero vmax
else:
    vmax = 0.1

cmap = plt.cm.RdBu_r.copy()
cmap.set_bad('#F0F0F0')
masked = np.ma.masked_invalid(heatmap_data)
im = ax_heat.imshow(masked, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax,
                     interpolation='nearest')

# Significance annotations
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        fdr = heatmap_fdr[i, j]
        if fdr < 0.001:
            ax_heat.text(j, i, '***', ha='center', va='center', fontsize=6, fontweight='bold')
        elif fdr < 0.01:
            ax_heat.text(j, i, '**', ha='center', va='center', fontsize=6, fontweight='bold')
        elif fdr < 0.05:
            ax_heat.text(j, i, '*', ha='center', va='center', fontsize=7)

# Labels
ax_heat.set_xticks(range(len(ordered_hallmarks)))
ax_heat.set_xticklabels([short_hallmark_name(h) for h in ordered_hallmarks],
                         rotation=65, ha='right', fontsize=7)
ax_heat.set_yticks(range(len(top5_names)))
ax_heat.set_yticklabels([f"{name}  [{moa}]" for name, moa in zip(top5_names, top5_moas)],
                         fontsize=9)

# Category boundary lines
for start, end in category_boundaries:
    if start > 0:
        ax_heat.axvline(start - 0.5, color='black', linewidth=1.0, linestyle='-')

ax_heat.set_title(f'Top {TOP_DRUGS} Reversing Drugs Across {n_hallmarks} MSigDB Hallmark Gene Sets\n'
                   f'(* FDR<0.05, ** FDR<0.01, *** FDR<0.001 | Blue=Reversing, Red=Mimicking)',
                   fontsize=13, fontweight='bold', pad=10)

# Colorbar
ax_cb = fig.add_subplot(gs[1, 1])
plt.colorbar(im, cax=ax_cb, label='Median Enrichment Score')

fig.savefig(PLOT_DIR / 'master_heatmap.png', dpi=250, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"  Heatmap saved: {PLOT_DIR / 'master_heatmap.png'}")

# ── Also make an expanded heatmap with top 20 drugs ────────────────────
top20 = drug_composite[:20]
top20_names = [d['drug'] for d in top20]
top20_moas = [d['moa'][:35] if d['moa'] else '?' for d in top20]

heatmap20 = np.full((len(top20_names), len(ordered_hallmarks)), np.nan)
for i, drug in enumerate(top20_names):
    drug_row = drug_to_row[drug]
    for j_orig, h_name in enumerate(hallmark_names_sorted):
        if h_name in hallmark_col_idx:
            j_new = hallmark_col_idx[h_name]
            heatmap20[i, j_new] = score_matrix[drug_row, j_orig]

fig, ax = plt.subplots(figsize=(24, 14))
masked20 = np.ma.masked_invalid(heatmap20)
v20 = max(abs(np.nanmin(heatmap20)), abs(np.nanmax(heatmap20)), 0.01)
im20 = ax.imshow(masked20, aspect='auto', cmap=cmap, vmin=-v20, vmax=v20, interpolation='nearest')
ax.set_xticks(range(len(ordered_hallmarks)))
ax.set_xticklabels([short_hallmark_name(h) for h in ordered_hallmarks], rotation=65, ha='right', fontsize=7)
ax.set_yticks(range(len(top20_names)))
ax.set_yticklabels([f"{n} [{m}]" for n, m in zip(top20_names, top20_moas)], fontsize=8)
for start, end in category_boundaries:
    if start > 0:
        ax.axvline(start - 0.5, color='black', linewidth=1.0)
ax.set_title(f'Top 20 Reversing Drugs Across {n_hallmarks} Hallmarks\n(Blue=Reversing, Red=Mimicking)',
             fontsize=13, fontweight='bold')
plt.colorbar(im20, ax=ax, label='Median Enrichment Score', shrink=0.6, pad=0.02)
fig.tight_layout()
fig.savefig(PLOT_DIR / 'expanded_heatmap_top20.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Expanded heatmap saved: {PLOT_DIR / 'expanded_heatmap_top20.png'}")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: NETWORK PLOT (Biology-Focused)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("Generating network plot...")
print("=" * 65)

G = nx.Graph()

# Add drug nodes (top 5)
for d in top5:
    G.add_node(d['drug'], node_type='drug',
               moa=d['moa'][:40] if d['moa'] else 'Unknown',
               composite=d['composite_score'])

# Add hallmark nodes and edges
for j, h_name in enumerate(ordered_hallmarks):
    cat = HALLMARK_TO_CATEGORY.get(h_name, 'Other')
    short_name = short_hallmark_name(h_name)

    # Check if any top5 drug has a significant connection to this hallmark
    has_connection = False
    for d in top5:
        drug_row = drug_to_row[d['drug']]
        j_orig = hallmark_names_sorted.index(h_name)
        score = score_matrix[drug_row, j_orig]
        fdr = fdr_matrix[drug_row, j_orig]
        if not np.isnan(score) and fdr < 0.1:  # relaxed threshold for network
            has_connection = True
            break

    if not has_connection:
        continue

    G.add_node(short_name, node_type='hallmark', category=cat, full_name=h_name)

    for d in top5:
        drug_row = drug_to_row[d['drug']]
        j_orig = hallmark_names_sorted.index(h_name)
        score = score_matrix[drug_row, j_orig]
        fdr = fdr_matrix[drug_row, j_orig]
        if not np.isnan(score) and fdr < 0.1:
            G.add_edge(d['drug'], short_name, weight=float(score), fdr=float(fdr))

if G.number_of_nodes() < 3:
    print("  Network too sparse; skipping network plot")
else:
    # Remove isolated nodes
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    # Layout: use spring layout with category-based initial positions
    pos = nx.spring_layout(G, k=3.0, iterations=200, seed=42)

    fig, ax = plt.subplots(figsize=(22, 18))
    ax.set_facecolor(C_BG)
    fig.patch.set_facecolor(C_BG)

    drugs_g = [n for n in G.nodes if G.nodes[n]['node_type'] == 'drug']
    hallmarks_g = [n for n in G.nodes if G.nodes[n]['node_type'] == 'hallmark']

    # Drug node sizes based on degree
    drug_sizes = [800 + G.degree(n) * 80 for n in drugs_g]

    # Hallmark node sizes based on degree
    hall_deg = [G.degree(n) for n in hallmarks_g]
    max_hd = max(hall_deg) if hall_deg else 1
    hall_sizes = [200 + (d / max_hd) * 500 for d in hall_deg]

    # Hallmark node colors by category
    hall_colors = [CATEGORY_COLORS.get(G.nodes[n].get('category', 'Other'), '#BFBFBF')
                   for n in hallmarks_g]

    # Edge styling
    edge_list = list(G.edges(data=True))
    edge_colors = ['#2166AC' if e[2]['weight'] < 0 else '#B2182B' for e in edge_list]
    edge_widths = [0.5 + abs(e[2]['weight']) * 20 for e in edge_list]  # scale by score
    edge_alphas = [0.8 if e[2]['fdr'] < 0.05 else 0.3 for e in edge_list]

    # Draw edges
    for e, c, w, a in zip(edge_list, edge_colors, edge_widths, edge_alphas):
        nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1])],
                               edge_color=[c], width=min(w, 5.0), alpha=a, ax=ax)

    # Draw hallmark nodes
    nx.draw_networkx_nodes(G, pos, nodelist=hallmarks_g, node_color=hall_colors,
                           node_size=hall_sizes, edgecolors='#333333',
                           linewidths=0.8, alpha=0.85, ax=ax)

    # Draw drug nodes (stars or circles)
    nx.draw_networkx_nodes(G, pos, nodelist=drugs_g, node_color=C_REV,
                           node_size=drug_sizes, edgecolors='#8B0000',
                           linewidths=1.5, alpha=0.95, ax=ax, node_shape='s')

    # Labels
    drug_labels = {}
    for n in drugs_g:
        moa = G.nodes[n].get('moa', '')
        drug_labels[n] = f"{n}\n({moa})" if moa and moa != 'Unknown' else n

    hall_labels = {n: n for n in hallmarks_g}

    # Offset labels slightly
    drug_lpos = {n: (pos[n][0], pos[n][1] + 0.06) for n in drugs_g}
    hall_lpos = {n: (pos[n][0], pos[n][1] - 0.04) for n in hallmarks_g}

    nx.draw_networkx_labels(G, drug_lpos, labels=drug_labels,
                            font_size=8, font_weight='bold', font_color='#4A0000', ax=ax)
    nx.draw_networkx_labels(G, hall_lpos, labels=hall_labels,
                            font_size=7, font_color='#1A1A1A', ax=ax)

    # Legend
    legend_elements = [
        Patch(facecolor=C_REV, edgecolor='#8B0000', label='Drug (top 5 reversing)'),
    ]
    for cat in ['Immune', 'Proliferation', 'Signaling', 'Metabolic', 'Stress & DNA Damage', 'Development', 'Other']:
        if any(G.nodes[n].get('category') == cat for n in hallmarks_g):
            legend_elements.append(
                Patch(facecolor=CATEGORY_COLORS[cat], edgecolor='#333', label=f'Hallmark: {cat}'))
    legend_elements.extend([
        Line2D([0], [0], color='#2166AC', linewidth=3, alpha=0.7, label='Reversing (neg score)'),
        Line2D([0], [0], color='#B2182B', linewidth=3, alpha=0.7, label='Mimicking (pos score)'),
        Line2D([0], [0], color='black', linewidth=2, alpha=0.8, label='FDR < 0.05 (solid)'),
        Line2D([0], [0], color='black', linewidth=1, alpha=0.3, label='FDR 0.05-0.10 (faint)'),
    ])

    ax.legend(handles=legend_elements, loc='upper left', fontsize=9, framealpha=0.92,
              edgecolor='gray', ncol=2)
    ax.set_title(f'Drug-Hallmark Pathway Network\n'
                 f'Top {TOP_DRUGS} Reversing Drugs x Significantly Affected Hallmark Pathways (FDR<0.10)\n'
                 f'Edge thickness = enrichment score magnitude, Color = direction',
                 fontsize=14, fontweight='bold')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(PLOT_DIR / 'drug_hallmark_network.png', dpi=200, bbox_inches='tight',
                facecolor=C_BG)
    plt.close()
    print(f"  Network saved: {PLOT_DIR / 'drug_hallmark_network.png'}")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6: SUMMARY REPORT
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("Writing summary report...")
print("=" * 65)

total_time = time.time() - t0

with open(OUT_DIR / 'summary.txt', 'w') as f:
    f.write("MSigDB Hallmark Enrichment Analysis\n")
    f.write("=" * 65 + "\n\n")
    f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Total runtime: {total_time/60:.1f} minutes\n\n")

    f.write(f"Database: {n_sigs:,} total signatures, {len(quality_idx):,} pass quality filter "
            f"(n_reps>={MIN_REPS})\n")
    f.write(f"Permutations: {N_PERM:,}\n\n")

    f.write(f"Hallmark Gene Sets Analyzed: {n_hallmarks}\n")
    f.write("-" * 65 + "\n")
    for h_name in hallmark_names_sorted:
        hr = hallmark_results[h_name]
        cat = HALLMARK_TO_CATEGORY.get(h_name, 'Other')
        f.write(f"  {short_hallmark_name(h_name):45s}  "
                f"genes={hr['n_genes_mapped']:3d}/{hr['n_genes_input']:3d}  "
                f"cat={cat}\n")

    f.write(f"\n\nTop {TOP_DRUGS} Reversing Drugs (Composite Score)\n")
    f.write("=" * 65 + "\n")
    for i, d in enumerate(top5):
        f.write(f"\n  {i+1}. {d['drug']}\n")
        f.write(f"     Mean median score: {d['mean_median_score']:.4f}\n")
        f.write(f"     Significant hallmarks (FDR<0.05): {d['n_sig_hallmarks']}\n")
        f.write(f"     MOA: {d['moa'] or 'Unknown'}\n")
        f.write(f"     Target: {d['target'] or 'Unknown'}\n")

        # List which hallmarks this drug most strongly reverses
        drug_row = drug_to_row[d['drug']]
        h_scores = []
        for j, h_name in enumerate(hallmark_names_sorted):
            s = score_matrix[drug_row, j]
            fdr = fdr_matrix[drug_row, j]
            if not np.isnan(s):
                h_scores.append((h_name, s, fdr))
        h_scores.sort(key=lambda x: x[1])
        f.write(f"     Top 5 reversed hallmarks:\n")
        for h_name, s, fdr in h_scores[:5]:
            fdr_s = f"{fdr:.4f}" if fdr >= 0.0001 else "<0.0001"
            f.write(f"       {short_hallmark_name(h_name):40s}  score={s:.4f}  FDR={fdr_s}\n")

    f.write(f"\n\nBiological Interpretation\n")
    f.write("=" * 65 + "\n")
    f.write("The top reversing drugs were identified by a composite score that\n")
    f.write("combines the mean enrichment score across all hallmark pathways\n")
    f.write("with the breadth of significant reversals (FDR<0.05). Drugs that\n")
    f.write("strongly and broadly reverse hallmark pathway signatures are\n")
    f.write("prioritized as potential therapeutic candidates.\n\n")

    # Pathway category summary
    f.write("Pathway Category Summary:\n")
    f.write("-" * 65 + "\n")
    for cat in ['Immune', 'Proliferation', 'Signaling', 'Metabolic', 'Stress & DNA Damage', 'Development']:
        members = [h for h in HALLMARK_CATEGORIES.get(cat, []) if h in hallmark_results]
        if not members:
            continue
        # Average of top 5 drugs' scores in this category
        cat_scores = []
        for d in top5:
            drug_row = drug_to_row[d['drug']]
            for h_name in members:
                j = hallmark_names_sorted.index(h_name)
                s = score_matrix[drug_row, j]
                if not np.isnan(s):
                    cat_scores.append(s)
        avg = np.mean(cat_scores) if cat_scores else 0
        f.write(f"  {cat:25s}: {len(members)} hallmarks, avg score for top5 drugs = {avg:.4f}\n")

    f.write(f"\n\nOutput Files\n")
    f.write("=" * 65 + "\n")
    f.write(f"  Per-hallmark results: {PER_HALLMARK_DIR}/\n")
    f.write(f"  Aggregated data:      {AGG_DIR}/\n")
    f.write(f"  Plots:                {PLOT_DIR}/\n")
    f.write(f"  This summary:         {OUT_DIR / 'summary.txt'}\n")

print(f"\nSummary saved: {OUT_DIR / 'summary.txt'}")
print(f"\nAnalysis complete in {total_time/60:.1f} minutes")
print(f"Outputs: {OUT_DIR}")

# List output files
print("\nKey output files:")
for p in sorted(PLOT_DIR.glob('*.png')):
    sz = p.stat().st_size
    print(f"  {p.relative_to(OUT_DIR)} ({sz/1e3:.1f} KB)")
for p in sorted(AGG_DIR.glob('*.csv')):
    sz = p.stat().st_size
    print(f"  {p.relative_to(OUT_DIR)} ({sz/1e3:.1f} KB)")
