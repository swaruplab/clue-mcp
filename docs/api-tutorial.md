# REST API Tutorial

The repo ships with a [FastAPI](https://fastapi.tiangolo.com/) wrapper around the enrichment engine. It exposes the same operations as the Python library over HTTP/JSON.

---

## 1. Install + start

```bash
pip install -e ".[api]"

# Default: reads data/processed/ from CWD
uvicorn api.main:app --host 0.0.0.0 --port 8000

# Or point at a non-default data directory
CMAP_DATA_DIR=/path/to/data/processed uvicorn api.main:app --port 8000
```

The first request to each class pays a ~10s cold-start (memory-mapping that class's matrices). Subsequent requests are millisecond-scale. Classes are loaded lazily — a server that only ever serves `/enrich/drug` never pays for the knockdown matrices.

Interactive Swagger docs: **http://localhost:8000/docs**

## 2. Endpoints

One enrichment endpoint per perturbation class — same request body, class-specific framing.

| Method | Path | Purpose |
|--------|------|---------|
| `POST` | `/enrich/drug`               | Small-molecule **drugs** (`trt_cp`) that reverse / mimic the signature |
| `POST` | `/enrich/knockdown`          | Gene **loss-of-function** (`trt_sh` + `trt_xpr`) that reverse / mimic — accepts `method` |
| `POST` | `/enrich/overexpression`     | Gene **gain-of-function** (`trt_oe`) that reverse / mimic |
| `GET`  | `/classes`                   | The three classes, their `pert_types`, framing, and which are installed |
| `GET`  | `/{pert_class}/cell_lines`   | Cell lines available for a class (`drug` / `knockdown` / `overexpression`) |
| `GET`  | `/{pert_class}/search?q=...` | Fuzzy perturbagen-name search within a class |
| `GET`  | `/health`                    | Readiness + per-class availability map |

A class whose dataset isn't installed returns **503** from its endpoints; `GET /classes`
and `GET /health` always work and report which classes are available.

## 3. Enrichment request

The request body is identical across the three `/enrich/*` endpoints (only `method` is
honored by `/enrich/knockdown`):

```bash
curl -X POST http://localhost:8000/enrich/drug \
  -H "Content-Type: application/json" \
  -d '{
    "genes_up":   ["APOE", "CLU", "TREM2", "CD68", "C1QB"],
    "genes_down": ["SYN1", "SNAP25", "SLC17A7"],
    "cell_line":  "A549",
    "top_n":      10
  }'
```

Response — both directions are returned in one call:

```json
{
  "pert_class": "drug",
  "reversing": [
    {
      "pert_name": "sirolimus",
      "cell_iname": "A549",
      "wtcs": -0.4823,
      "moa": "mTOR inhibitor",
      "target": "MTOR|FKBP1A",
      "n_reps": 6,
      "method": null
    }
  ],
  "mimicking": [ ... ],
  "query_stats": {
    "n_up_mapped": 5,
    "n_down_mapped": 3,
    "unmapped_genes": []
  }
}
```

`wtcs` is the Weighted Connectivity Score from Subramanian *et al.* 2017. The `reversing`
list is sorted most-negative-first (strongest reversal); `mimicking` is most-positive-first.

To restrict a knockdown query to one technology:

```bash
curl -X POST http://localhost:8000/enrich/knockdown \
  -H "Content-Type: application/json" \
  -d '{"genes_up": ["APOE", "CLU"], "method": "CRISPR", "top_n": 10}'
```

## 4. Python client

```python
import requests

r = requests.post(
    "http://localhost:8000/enrich/drug",
    json={
        "genes_up":   ["APOE", "CLU", "TREM2"],
        "genes_down": ["SYN1", "SNAP25"],
        "top_n":      25,
    },
    timeout=120,
)
r.raise_for_status()
body = r.json()
reversing = body["reversing"]   # candidate therapeutics
mimicking = body["mimicking"]   # signature-inducers
```

## 5. Production deployment

For multi-worker / public-facing deployment:

```bash
gunicorn api.main:app \
  -k uvicorn.workers.UvicornWorker \
  -w 1 \
  --timeout 120 \
  --bind 0.0.0.0:8000
```

> **Prefer few workers.** The matrices are memory-mapped, so the OS page cache is shared across workers on one host — but each class still backs a large file. Start with `-w 1` and scale horizontally across machines rather than piling workers onto one node.

## 6. Behind a reverse proxy

If you're putting the API behind nginx, allow at least 60 s for the first request (cold start):

```nginx
location /cmap/ {
    proxy_pass http://127.0.0.1:8000/;
    proxy_read_timeout 120s;
    proxy_connect_timeout 30s;
}
```

---

## Next steps

- [MCP server](mcp-setup.md) — same engine, but Claude/LLMs can call it directly
- [Python library](python-tutorial.md) — skip the HTTP layer for batch jobs
