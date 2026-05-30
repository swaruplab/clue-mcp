"""
REST API for CMap/L1000 connectivity enrichment.

Three perturbation classes, each with its own endpoint and scientific framing:

    POST /enrich/drug            small-molecule drugs            (trt_cp)
    POST /enrich/knockdown       gene loss-of-function (sh+xpr)  (trt_sh/trt_xpr)
    POST /enrich/overexpression  gene gain-of-function (oe)      (trt_oe)

All three share one scoring engine (cmap_enrichment). Community-hosted
preservation of CLUE.io / LINCS 2020 data.
"""

from __future__ import annotations

import os
import sys
from typing import Optional

from fastapi import FastAPI, HTTPException, Path as PathParam
from pydantic import BaseModel, Field

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cmap_enrichment.perturbation_classes import PERTURBATION_CLASSES, get_class
from cmap_enrichment.registry import EngineRegistry

app = FastAPI(
    title="CMap Connectivity Enrichment API",
    description=(
        "Query the Connectivity Map (CMap/L1000) LINCS 2020 dataset for perturbations "
        "that reverse or mimic a gene expression signature — across drugs (trt_cp), "
        "gene knockdowns (trt_sh/trt_xpr), and gene overexpression (trt_oe)."
    ),
    version="2.0.0",
)

registry = EngineRegistry()


# -- models ----------------------------------------------------------------

class EnrichmentRequest(BaseModel):
    genes_up: list[str] = Field(
        ...,
        description="Upregulated genes (HGNC symbols) in your condition",
        examples=[["APOE", "CLU", "TREM2", "C1QB", "CD68"]],
    )
    genes_down: list[str] = Field(
        default=[],
        description="Downregulated genes (optional, improves specificity)",
        examples=[["SYN1", "SNAP25", "SLC17A7"]],
    )
    cell_line: Optional[str] = Field(default=None, description="Filter to one cell line")
    method: Optional[str] = Field(
        default=None,
        description="Knockdown only: restrict to 'shRNA' or 'CRISPR' (default = pool both)",
    )
    top_n: int = Field(default=50, ge=1, le=500, description="Hits per direction")


class Hit(BaseModel):
    pert_name: str
    cell_iname: Optional[str] = None
    wtcs: float
    n_reps: Optional[int] = None
    moa: Optional[str] = None
    target: Optional[str] = None
    method: Optional[str] = None


class EnrichmentResponse(BaseModel):
    pert_class: str
    reversing: list[Hit]
    mimicking: list[Hit]
    query_stats: dict


def _hits(df) -> list[Hit]:
    import pandas as pd

    out = []
    for _, r in df.iterrows():
        def g(col):
            v = r.get(col) if col in df.columns else None
            return None if v is None or (isinstance(v, float) and pd.isna(v)) or str(v) in ("nan", "<NA>") else v

        out.append(Hit(
            pert_name=str(r["pert_name"]),
            cell_iname=g("cell_iname"),
            wtcs=float(r["wtcs"]),
            n_reps=int(r["n_reps"]) if "n_reps" in df.columns and pd.notna(r.get("n_reps")) else None,
            moa=g("moa"),
            target=g("target"),
            method=g("method"),
        ))
    return out


def _enrich(pert_class: str, req: EnrichmentRequest) -> EnrichmentResponse:
    key = get_class(pert_class).key
    if not registry.available().get(key):
        raise HTTPException(
            status_code=503,
            detail=f"Dataset '{key}' is not available on this server. "
                   f"Build or download it first (scripts/download_data.sh {key}).",
        )
    try:
        res = registry.get(key).rank_perturbations(
            genes_up=req.genes_up,
            genes_down=req.genes_down,
            cell_line=req.cell_line,
            method=req.method,
            top_n=req.top_n,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    return EnrichmentResponse(
        pert_class=key,
        reversing=_hits(res["reversing"]),
        mimicking=_hits(res["mimicking"]),
        query_stats=res["stats"],
    )


# -- endpoints -------------------------------------------------------------

@app.post("/enrich/drug", response_model=EnrichmentResponse)
async def enrich_drug(request: EnrichmentRequest):
    """Find small-molecule DRUGS (trt_cp) that reverse or mimic the signature."""
    return _enrich("drug", request)


@app.post("/enrich/knockdown", response_model=EnrichmentResponse)
async def enrich_knockdown(request: EnrichmentRequest):
    """Find gene LOSS-of-function (shRNA + CRISPR) perturbations that reverse or mimic."""
    return _enrich("knockdown", request)


@app.post("/enrich/overexpression", response_model=EnrichmentResponse)
async def enrich_overexpression(request: EnrichmentRequest):
    """Find gene GAIN-of-function (ORF over-expression) perturbations that reverse or mimic."""
    return _enrich("overexpression", request)


@app.get("/classes")
async def list_classes():
    """Describe the three perturbation classes and their availability."""
    avail = registry.available()
    return {
        key: {
            "tool_name": spec.tool_name,
            "pert_types": list(spec.pert_types),
            "perturbagen": spec.perturbagen_plural,
            "reversing_meaning": spec.reversing_meaning,
            "mimicking_meaning": spec.mimicking_meaning,
            "available": avail.get(key, False),
        }
        for key, spec in PERTURBATION_CLASSES.items()
    }


@app.get("/{pert_class}/cell_lines")
async def list_cell_lines(pert_class: str = PathParam(...)):
    key = get_class(pert_class).key
    if not registry.available().get(key):
        raise HTTPException(status_code=503, detail=f"Dataset '{key}' not available.")
    return registry.get(key).list_cell_lines()


@app.get("/{pert_class}/search")
async def search(pert_class: str, q: str):
    key = get_class(pert_class).key
    if not registry.available().get(key):
        raise HTTPException(status_code=503, detail=f"Dataset '{key}' not available.")
    return registry.get(key).search_perturbagens(q).to_dict(orient="records")


@app.get("/health")
async def health():
    return {"status": "ok", "classes": registry.available()}
