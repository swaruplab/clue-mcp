# Task: upload the clue-mcp database to Zenodo

**For the Claude Code agent running on the server that holds the processed data.**

Goal: deposit the database as **per-class archives** so anyone can fetch exactly
the classes they need with `scripts/download_data.sh`. There are up to **four**
tarballs — a small shared `base` plus one archive per perturbation class
(`drug`, `knockdown`, `overexpression`). The `base` archive is always required;
the class archives are optional and uploaded only for classes that have been built.

The archive names and contents below are the source of truth consumed by
`cmap_enrichment/zenodo_manifest.json` and `scripts/download_data.sh` — keep them
exactly as written.

---

## 1. Archives to build (and their contents)

Find the finished `data/processed/` directory (on the Swarup cluster this is
typically `/dfs7/swaruplab/shared_lab/Clue_database/data/processed/`). It holds the
shared files plus one subfolder per built class:

```
data/processed/
├── gene_names.txt                    # shared
├── cmap.duckdb                       # shared (index over all built classes)
├── drug/{zscore_matrix.npy, rank_matrix.npy, metadata.parquet}
├── knockdown/{zscore_matrix.npy, rank_matrix.npy, metadata.parquet}
└── overexpression/{zscore_matrix.npy, rank_matrix.npy, metadata.parquet}
```

| Archive | Contents | Approx size | Required |
|---------|----------|-------------|:--------:|
| `clue_mcp_base.tar.gz` | `gene_names.txt`, `cmap.duckdb` | 9 MB | ✅ always |
| `clue_mcp_drug.tar.gz` | `drug/zscore_matrix.npy`, `drug/rank_matrix.npy`, `drug/metadata.parquet` | ~9 GB | if built |
| `clue_mcp_knockdown.tar.gz` | `knockdown/…` (3 files) | larger (sh+xpr) | if built |
| `clue_mcp_overexpression.tar.gz` | `overexpression/…` (3 files) | ~2 GB | if built |

Each class archive keeps its **leading class folder** (`drug/…`) so it extracts
straight into `data/processed/`.

### Do NOT upload
- Raw GCTX files (`*.gctx`) — large and © Broad; users get those from clue.io.
- Raw metadata txt (`siginfo_beta.txt`, `compoundinfo_beta.txt`, `geneinfo_beta.txt`, `cellinfo_beta.txt`).
- Anything in `data/intermediate/`.
- Analysis outputs (`analysis/…`) — those already live in the GitHub repo.

---

## 2. Build the archives

A helper builds every archive that has inputs present and writes a `.sha256`
next to each:

```bash
bash scripts/make_zenodo_archives.sh            # reads ./data/processed, writes to $HOME
#   …or:  bash scripts/make_zenodo_archives.sh /path/to/data/processed /path/to/out
```

It only builds archives whose inputs exist, so running it after only the drug
class is built produces `clue_mcp_base.tar.gz` + `clue_mcp_drug.tar.gz`.

Sanity-check a class archive keeps its leading folder:

```bash
tar tzf ~/clue_mcp_drug.tar.gz | head   # expect: drug/zscore_matrix.npy, …
```

---

## 3. Deposit on Zenodo

Upload all of the produced `.tar.gz` files (and their `.sha256` files) to a
**single Zenodo record**. Use these deposit settings so the repo's docs stay accurate:

| Field | Value |
|-------|-------|
| **Files** | every `clue_mcp_*.tar.gz` + matching `.sha256` |
| **Upload type** | Dataset |
| **Title** | CMap / L1000 LINCS 2020 — processed perturbation signatures for clue-mcp |
| **Authors** | Swarup Lab, UC Irvine |
| **License** | CC-BY-4.0 |
| **Version** | `2020.1` |

**Description:**

> Pre-processed Connectivity Map perturbation signatures from the LINCS 2020
> release (Dec 2020), split into three classes — small-molecule **drugs**
> (`trt_cp`), gene **knockdown** (`trt_sh` shRNA + `trt_xpr` CRISPR), and gene
> **overexpression** (`trt_oe`) — over 12,328 genes, each with a pre-computed
> gene-rank matrix for fast WTCS enrichment. Built for the clue-mcp engine
> (github.com/swaruplab/clue-mcp). Derived from data © the Broad Institute; see
> https://clue.io/connectopedia/data_use_policy.

**Related identifiers:** *is derived from* `s3://macchiato.clue.io/builds/LINCS2020/`;
*cites* `doi:10.1016/j.cell.2017.10.049`; *is supplement to* the clue-mcp GitHub repo.

### Optional: upload via the Zenodo API instead of the web UI

```bash
# Needs a Zenodo personal access token with deposit:write scope
ZENODO_TOKEN=...                          # set this
BASE=https://zenodo.org/api               # use https://sandbox.zenodo.org/api to test first

# Create a new empty deposition; capture its id and bucket URL
resp=$(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" \
  -H "Content-Type: application/json" -X POST "$BASE/deposit/depositions" -d '{}')
echo "$resp" | python3 -c 'import sys,json;d=json.load(sys.stdin);print("id=",d["id"]);print("bucket=",d["links"]["bucket"])'

# Upload each archive to the bucket URL printed above
for f in ~/clue_mcp_*.tar.gz ~/clue_mcp_*.tar.gz.sha256; do
  curl -s --upload-file "$f" \
    -H "Authorization: Bearer $ZENODO_TOKEN" \
    "<bucket-url>/$(basename "$f")"
done
```

Then add the metadata (title/description/license/etc.) and **publish** — either in
the web UI or via `PUT $BASE/deposit/depositions/<id>` followed by
`POST $BASE/deposit/depositions/<id>/actions/publish`. Publishing mints the DOI.

---

## 4. Wire the deposit back into the repo

After publishing, edit `cmap_enrichment/zenodo_manifest.json`:

1. Set `record_id` (the number in `https://zenodo.org/records/<record-id>`) and `doi`.
2. Paste each archive's **sha256** (from the `.sha256` files) into its `sha256` field,
   replacing `FILL_AFTER_DEPOSIT`. The downloader verifies these.
3. Replace the `10.5281/zenodo.XXXXXXX` placeholder in `README.md`.

The full post-publish checklist is in [`docs/data-release.md`](docs/data-release.md).

---

## 5. Verify it worked

```bash
# Downloads base + every available class, verifying each archive's sha256
ZENODO_RECORD=<record-id> bash scripts/download_data.sh
#   …or one class:  ZENODO_RECORD=<record-id> bash scripts/download_data.sh drug
```
