#!/usr/bin/env bash
# Download and verify the processed clue-mcp data bundles from Zenodo.
#
# Usage:
#   bash scripts/download_data.sh                 # all classes + shared base
#   bash scripts/download_data.sh drug            # just the drug class (+ base)
#   bash scripts/download_data.sh knockdown overexpression
#
# Config:
#   CMAP_DATA_DIR   where to extract (default: <repo>/data/processed)
#   ZENODO_RECORD   Zenodo record id (overrides record_id in zenodo_manifest.json)
#
# The 'base' archive (gene_names.txt + cmap.duckdb) is always fetched.
# Each class archive extracts into <CMAP_DATA_DIR>/<class>/.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$ROOT/cmap_enrichment/zenodo_manifest.json"
DATA_DIR="${CMAP_DATA_DIR:-$ROOT/data/processed}"

command -v python3 >/dev/null || { echo "python3 required" >&2; exit 1; }

# Read a value out of the manifest with python (no jq dependency).
mq() { python3 -c "import json,sys; d=json.load(open('$MANIFEST')); print(eval(sys.argv[1]))" "$1"; }

RECORD="${ZENODO_RECORD:-$(mq 'd["record_id"]')}"
BASE_URL="$(mq 'd["base_url"]')"
if [[ "$RECORD" == *X* || -z "$RECORD" ]]; then
    echo "ERROR: Zenodo record id is not set." >&2
    echo "  Set ZENODO_RECORD=<id>, or fill record_id in cmap_enrichment/zenodo_manifest.json" >&2
    echo "  (the deposit + DOI workflow is described in ZENODO_UPLOAD.md)." >&2
    exit 1
fi

CLASSES=("$@")
if [ "${#CLASSES[@]}" -eq 0 ]; then
    CLASSES=(drug knockdown overexpression)
fi
# Always include the shared base archive, first.
ARCHIVES=(base "${CLASSES[@]}")

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

fetch_one() {
    local key="$1"
    local file sha url
    file="$(mq "d['archives']['$key']['file']")"
    sha="$(mq "d['archives']['$key']['sha256']")"
    url="$BASE_URL/$RECORD/files/$file"

    echo ">> $key : $file"
    if command -v wget >/dev/null; then
        wget -c -O "$file" "$url"
    else
        curl -fL -C - -o "$file" "$url"
    fi

    if [[ "$sha" != FILL_* && -n "$sha" ]]; then
        echo "   verifying sha256..."
        if command -v sha256sum >/dev/null; then
            echo "$sha  $file" | sha256sum -c -
        else
            test "$(shasum -a 256 "$file" | awk '{print $1}')" = "$sha" \
                || { echo "   checksum MISMATCH for $file" >&2; exit 1; }
        fi
    else
        echo "   (no checksum in manifest yet — skipping verification)"
    fi

    echo "   extracting..."
    tar xzf "$file"
    rm -f "$file"
}

for key in "${ARCHIVES[@]}"; do
    fetch_one "$key"
done

echo "Done. Data in: $DATA_DIR"
ls -1 "$DATA_DIR"
