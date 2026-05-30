# Releasing the database (maintainers)

This page is for **maintainers** depositing the processed database so that
users can fetch it with `scripts/download_data.sh` instead of rebuilding it.

The database ships as **per-class archives** on one Zenodo record: a small shared
`clue_mcp_base.tar.gz` (always required) plus one archive per perturbation class
(`clue_mcp_drug`, `clue_mcp_knockdown`, `clue_mcp_overexpression`) — upload only the
classes you have built. The archive names and contents are defined in
[`cmap_enrichment/zenodo_manifest.json`](https://github.com/swaruplab/clue-mcp/blob/main/cmap_enrichment/zenodo_manifest.json),
which the downloader reads.

---

## 1. Assemble the archives

From a machine that has a finished `data/processed/` (with `gene_names.txt`,
`cmap.duckdb`, and one `<class>/` subfolder per built class):

```bash
bash scripts/make_zenodo_archives.sh        # reads ./data/processed, writes to $HOME
```

It builds every archive whose inputs are present and writes a `.sha256` beside each:

| Archive | Contents |
|---------|----------|
| `clue_mcp_base.tar.gz` | `gene_names.txt`, `cmap.duckdb` |
| `clue_mcp_drug.tar.gz` | `drug/{zscore_matrix.npy, rank_matrix.npy, metadata.parquet}` |
| `clue_mcp_knockdown.tar.gz` | `knockdown/…` |
| `clue_mcp_overexpression.tar.gz` | `overexpression/…` |

!!! note "Class archives keep their leading folder"
    Each class archive contains `drug/…`, `knockdown/…`, etc. so it extracts
    straight into `data/processed/`. Only the `base` archive is flat.

---

## 2. Create the Zenodo deposit

At [zenodo.org/deposit](https://zenodo.org/deposit), **New upload**, attach every
`clue_mcp_*.tar.gz` (+ its `.sha256`) to a single record, and fill in:

| Field | Value |
|-------|-------|
| **Files** | every `clue_mcp_*.tar.gz` (+ the `.sha256` files) |
| **Upload type** | Dataset |
| **Title** | *CMap / L1000 LINCS 2020 — processed perturbation signatures for clue-mcp* |
| **Authors** | Swarup Lab, UC Irvine |
| **License** | CC-BY-4.0 |
| **Version** | Match the build (e.g. `2020.1`) |

**Description** (suggested):

> Pre-processed Connectivity Map perturbation signatures from the LINCS 2020
> release (Dec 2020), split into drug (`trt_cp`), knockdown (`trt_sh`+`trt_xpr`),
> and overexpression (`trt_oe`) classes over 12,328 genes, each with a pre-computed
> gene-rank matrix for fast WTCS enrichment. Built for the clue-mcp engine
> (github.com/swaruplab/clue-mcp). Derived from data © the Broad Institute; see
> https://clue.io/connectopedia/data_use_policy.

**Related/alternate identifiers:**

- *is derived from* — `s3://macchiato.clue.io/builds/LINCS2020/`
- *cites* — `doi:10.1016/j.cell.2017.10.049` (Subramanian et al. 2017)
- *is supplement to* — the clue-mcp GitHub repository URL

Publish to mint the DOI.

---

## 3. Wire the deposit back into the repo

After publishing, note the numeric **record id** (the number in
`https://zenodo.org/records/<record-id>`) and the **DOI**, then edit
[`cmap_enrichment/zenodo_manifest.json`](https://github.com/swaruplab/clue-mcp/blob/main/cmap_enrichment/zenodo_manifest.json):

- [ ] Set `record_id` and `doi`.
- [ ] Paste each archive's **sha256** (from the `.sha256` files) into its `sha256` field, replacing `FILL_AFTER_DEPOSIT` — the downloader verifies these.
- [ ] Replace the `10.5281/zenodo.XXXXXXX` placeholder in `README.md` ("Getting the data") with the real DOI.
- [ ] Replace the DOI mention on the [Quickstart](quickstart.md) page if you pinned a specific record.
- [ ] (Optional) add a Zenodo DOI badge to the README.
- [ ] Commit, and verify end-to-end in a clean clone:

```bash
git clone https://github.com/swaruplab/clue-mcp.git && cd clue-mcp
pip install -e .
bash scripts/download_data.sh          # downloads + verifies base + each class
python examples/quickstart_library.py  # confirms the engine loads and queries
```

---

## Re-releasing a new build

Zenodo versions deposits: use **New version** on the existing record so the
**concept DOI** keeps resolving to the latest. Bump the `Version` field and update
`ZENODO_RECORD` to the new record id.
