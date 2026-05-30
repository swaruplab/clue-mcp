"""
MCP server for CMap/L1000 connectivity enrichment.

Exposes three scientifically-framed enrichment tools so an LLM picks the right
interpretation for the user's question:

    cmap_drug_enrichment          small-molecule DRUGS              (trt_cp)
    cmap_target_knockdown         gene LOSS-of-function (sh + xpr)  (trt_sh/trt_xpr)
    cmap_target_overexpression    gene GAIN-of-function (oe)        (trt_oe)

Runs locally; reads the per-class numpy/parquet bundles from CMAP_DATA_DIR.

Run:  python mcp_server/server.py
"""

from __future__ import annotations

import os
import sys
from typing import Any

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cmap_enrichment.perturbation_classes import (
    PERTURBATION_CLASSES,
    TOOL_TO_CLASS,
    get_class,
)
from cmap_enrichment.registry import EngineRegistry

server = Server("clue-mcp")
registry = EngineRegistry()


# -- tool catalog ----------------------------------------------------------

def _enrichment_input_schema(spec) -> dict:
    props: dict[str, Any] = {
        "genes_up": {
            "type": "array",
            "items": {"type": "string"},
            "description": "Upregulated gene symbols (HGNC) in your condition",
        },
        "genes_down": {
            "type": "array",
            "items": {"type": "string"},
            "description": "Downregulated gene symbols (optional, improves specificity)",
            "default": [],
        },
        "cell_line": {
            "type": "string",
            "description": "Filter to a single cell line (e.g. A549, MCF7, PC3)",
        },
        "direction": {
            "type": "string",
            "enum": ["reversing", "mimicking", "both"],
            "description": (
                "Which end to report: 'reversing' (opposes the signature), "
                "'mimicking' (reproduces it), or 'both'"
            ),
            "default": "both",
        },
        "top_n": {
            "type": "integer",
            "description": "Number of hits per direction",
            "default": 20,
        },
    }
    if spec.has_method:
        props["method"] = {
            "type": "string",
            "enum": list(spec.methods),
            "description": "Restrict to one knockdown method (shRNA or CRISPR); default = pool both",
        }
    return {"type": "object", "properties": props, "required": ["genes_up"]}


@server.list_tools()
async def list_tools() -> list[Tool]:
    tools: list[Tool] = []
    for spec in PERTURBATION_CLASSES.values():
        tools.append(
            Tool(
                name=spec.tool_name,
                description=spec.description,
                inputSchema=_enrichment_input_schema(spec),
            )
        )

    tools.append(
        Tool(
            name="cmap_list_perturbation_classes",
            description=(
                "List the three CMap perturbation classes (drug / knockdown / "
                "overexpression), their scientific framing, and whether each dataset "
                "is currently available. Call this first if unsure which enrichment "
                "tool fits the question."
            ),
            inputSchema={"type": "object", "properties": {}},
        )
    )
    tools.append(
        Tool(
            name="cmap_search_perturbagens",
            description=(
                "Search perturbagens (drug names or perturbed gene symbols) by partial "
                "match within one perturbation class."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {"type": "string", "description": "Partial name to search"},
                    "pert_class": {
                        "type": "string",
                        "enum": list(PERTURBATION_CLASSES),
                        "description": "Which class to search",
                        "default": "drug",
                    },
                },
                "required": ["query"],
            },
        )
    )
    tools.append(
        Tool(
            name="cmap_list_cell_lines",
            description="List available cell lines for a perturbation class.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pert_class": {
                        "type": "string",
                        "enum": list(PERTURBATION_CLASSES),
                        "default": "drug",
                    },
                },
            },
        )
    )
    return tools


# -- result narration ------------------------------------------------------

def _fmt(val: Any, default: str = "N/A") -> str:
    if val is None:
        return default
    s = str(val)
    if s == "nan" or s == "<NA>" or s == "":
        return default
    return s


def _narrate_table(df, spec, heading: str, meaning: str) -> list[str]:
    lines = [heading, f"  ({meaning})", "-" * 78]
    for _, row in df.iterrows():
        parts = [f"  {_fmt(row['pert_name']):<24} | WTCS: {row['wtcs']:+.4f}",
                 f"Cell: {_fmt(row.get('cell_iname')):<7}"]
        if spec.has_moa:
            parts.append(f"MOA: {_fmt(row.get('moa'))}")
            parts.append(f"Target: {_fmt(row.get('target'))}")
        if spec.has_method:
            parts.append(f"Method: {_fmt(row.get('method'))}")
        lines.append(" | ".join(parts))
    return lines


def _run_enrichment(tool_name: str, args: dict[str, Any]) -> str:
    spec = get_class(tool_name)
    if not registry.available().get(spec.key, False):
        return (
            f"The '{spec.key}' dataset ({'/'.join(spec.pert_types)}) is not available "
            f"in CMAP_DATA_DIR yet. Build it with scripts/01-03 for these pert_types, "
            f"or download the bundle (scripts/download_data.sh {spec.key})."
        )

    engine = registry.get(spec.key)
    direction = args.get("direction", "both")
    top_n = args.get("top_n", 20)

    res = engine.rank_perturbations(
        genes_up=args["genes_up"],
        genes_down=args.get("genes_down", []),
        cell_line=args.get("cell_line"),
        method=args.get("method"),
        top_n=top_n,
    )
    stats = res["stats"]

    lines = [
        f"CMap {spec.perturbagen_singular.title()} Enrichment (LINCS 2020, {'/'.join(spec.pert_types)})",
        f"Query: {len(args['genes_up'])} up, {len(args.get('genes_down', []))} down  |  "
        f"mapped {stats['n_up_mapped']} up / {stats['n_down_mapped']} down  |  "
        f"{stats['n_candidates']:,} candidate {spec.perturbagen_plural}",
    ]
    if stats["unmapped_genes"]:
        lines.append(f"Unmapped (not in L1000): {', '.join(stats['unmapped_genes'])}")
    lines.append("")

    if direction in ("reversing", "both"):
        lines += _narrate_table(
            res["reversing"], spec,
            f"TOP REVERSING {spec.perturbagen_plural.upper()} (negative WTCS):",
            spec.reversing_meaning,
        )
        lines.append("")
    if direction in ("mimicking", "both"):
        lines += _narrate_table(
            res["mimicking"], spec,
            f"TOP MIMICKING {spec.perturbagen_plural.upper()} (positive WTCS):",
            spec.mimicking_meaning,
        )
    return "\n".join(lines)


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[TextContent]:
    try:
        if name in TOOL_TO_CLASS:
            return [TextContent(type="text", text=_run_enrichment(name, arguments))]

        if name == "cmap_list_perturbation_classes":
            avail = registry.available()
            lines = ["CMap perturbation classes:", ""]
            for spec in PERTURBATION_CLASSES.values():
                status = "available" if avail.get(spec.key) else "NOT built/downloaded"
                lines += [
                    f"## {spec.tool_name}  [{status}]",
                    f"  pert_types : {', '.join(spec.pert_types)}",
                    f"  perturbagen: {spec.perturbagen_plural} ({spec.unit})",
                    f"  reversing  : {spec.reversing_meaning}",
                    f"  mimicking  : {spec.mimicking_meaning}",
                    "",
                ]
            return [TextContent(type="text", text="\n".join(lines))]

        if name == "cmap_search_perturbagens":
            key = get_class(arguments.get("pert_class", "drug")).key
            if not registry.available().get(key):
                return [TextContent(type="text", text=f"Dataset '{key}' not available.")]
            df = registry.get(key).search_perturbagens(arguments["query"])
            if df.empty:
                return [TextContent(type="text", text=f"No {key} perturbagens matching '{arguments['query']}'.")]
            return [TextContent(type="text", text=df.to_string(index=False))]

        if name == "cmap_list_cell_lines":
            key = get_class(arguments.get("pert_class", "drug")).key
            if not registry.available().get(key):
                return [TextContent(type="text", text=f"Dataset '{key}' not available.")]
            cells = registry.get(key).list_cell_lines()
            return [TextContent(type="text", text=f"Cell lines for {key} ({len(cells)}):\n" + ", ".join(cells))]

        return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except (ValueError, FileNotFoundError) as e:
        return [TextContent(type="text", text=f"Error: {e}")]


async def main():
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


def run() -> None:
    """Console-script entry point (clue-mcp-server)."""
    import asyncio

    asyncio.run(main())


if __name__ == "__main__":
    run()
