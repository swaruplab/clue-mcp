#!/usr/bin/env bash
# Build the Zenodo upload archives from a finished data/processed/ tree.
#
# Produces (in $OUT_DIR, default ~):
#   clue_mcp_base.tar.gz            gene_names.txt + cmap.duckdb
#   clue_mcp_drug.tar.gz           drug/{zscore,rank}_matrix.npy + drug/metadata.parquet
#   clue_mcp_knockdown.tar.gz      knockdown/...
#   clue_mcp_overexpression.tar.gz overexpression/...
# plus a .sha256 next to each. Only builds archives whose inputs exist.
#
# Usage:  bash scripts/make_zenodo_archives.sh [data/processed dir] [out dir]
# After publishing, paste the record id + sha256 values into
# cmap_enrichment/zenodo_manifest.json (see ZENODO_UPLOAD.md).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROC="${1:-$ROOT/data/processed}"
OUT="${2:-$HOME}"
cd "$PROC"

sha_of() {
    if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi
}

# Use pigz (parallel gzip) when present, else plain gzip. -h follows symlinks so
# classes whose matrices are symlinked (e.g. drug -> shared store) archive their
# real contents, not a dangling link.
if command -v pigz >/dev/null; then COMPRESS=(--use-compress-program=pigz); else COMPRESS=(-z); fi

make_archive() {
    local name="$1"; shift
    local tarball="$OUT/$name"
    echo ">> building $tarball"
    tar -c "${COMPRESS[@]}" -h -f "$tarball" "$@"
    sha_of "$tarball" | tee "$tarball.sha256"
}

# Shared base
if [[ -f gene_names.txt && -f cmap.duckdb ]]; then
    make_archive clue_mcp_base.tar.gz gene_names.txt cmap.duckdb
else
    echo "-- skip base (need gene_names.txt + cmap.duckdb)"
fi

# Per-class
for cls in drug knockdown overexpression; do
    if [[ -f "$cls/rank_matrix.npy" ]]; then
        make_archive "clue_mcp_${cls}.tar.gz" \
            "$cls/zscore_matrix.npy" "$cls/rank_matrix.npy" "$cls/metadata.parquet"
    else
        echo "-- skip $cls (not built)"
    fi
done

echo "Archives written to $OUT"
