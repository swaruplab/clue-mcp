#!/usr/bin/env bash
# Run a handful of example requests against a locally-running CMap API.
#
# Start the server first (in another terminal):
#   uvicorn api.main:app --port 8000

set -euo pipefail

HOST="${CMAP_API:-http://localhost:8000}"

echo "==> Health check (per-class availability)"
curl -s "$HOST/health" | python -m json.tool

echo
echo "==> Describe the three perturbation classes"
curl -s "$HOST/classes" | python -m json.tool

echo
echo "==> Drug cell lines (first 5)"
curl -s "$HOST/drug/cell_lines" | python -c 'import sys, json; print(json.load(sys.stdin)[:5])'

echo
echo "==> Search drugs matching 'statin'"
curl -s "$HOST/drug/search?q=statin" | python -m json.tool

echo
echo "==> Drug enrichment (AD-like microglial signature)"
curl -s -X POST "$HOST/enrich/drug" \
    -H "Content-Type: application/json" \
    -d '{
      "genes_up":   ["APOE", "CLU", "TREM2", "CD68", "C1QB"],
      "genes_down": ["SYN1", "SNAP25", "SLC17A7"],
      "top_n": 5
    }' | python -m json.tool

echo
echo "==> Target knockdown enrichment, CRISPR only (skips with 503 if not installed)"
curl -s -X POST "$HOST/enrich/knockdown" \
    -H "Content-Type: application/json" \
    -d '{
      "genes_up":   ["APOE", "CLU", "TREM2", "CD68", "C1QB"],
      "genes_down": ["SYN1", "SNAP25", "SLC17A7"],
      "method": "CRISPR",
      "top_n": 5
    }' | python -m json.tool
