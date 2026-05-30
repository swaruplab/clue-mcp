"""MCP tool catalog + result-narration framing (no live MCP transport needed)."""

import pytest

pytest.importorskip("mcp", reason="mcp SDK not installed; skipping MCP server tests")

import mcp_server.server as srv  # noqa: E402
from cmap_enrichment.registry import EngineRegistry  # noqa: E402


def test_three_enrichment_tools_plus_helpers(data_dir, monkeypatch):
    monkeypatch.setattr(srv, "registry", EngineRegistry(data_dir))
    import asyncio

    tools = asyncio.run(srv.list_tools())
    names = {t.name for t in tools}
    assert {"cmap_drug_enrichment", "cmap_target_knockdown",
            "cmap_target_overexpression"} <= names
    assert "cmap_list_perturbation_classes" in names


def test_knockdown_tool_exposes_method_enum(data_dir, monkeypatch):
    monkeypatch.setattr(srv, "registry", EngineRegistry(data_dir))
    import asyncio

    tools = {t.name: t for t in asyncio.run(srv.list_tools())}
    kd = tools["cmap_target_knockdown"].inputSchema["properties"]
    assert kd["method"]["enum"] == ["shRNA", "CRISPR"]
    assert "method" not in tools["cmap_drug_enrichment"].inputSchema["properties"]


def test_result_narration_is_class_specific(data_dir, monkeypatch):
    monkeypatch.setattr(srv, "registry", EngineRegistry(data_dir))

    drug = srv._run_enrichment("cmap_drug_enrichment",
                               {"genes_up": ["AAA", "BBB"], "top_n": 2})
    kd = srv._run_enrichment("cmap_target_knockdown",
                             {"genes_up": ["AAA", "BBB"], "top_n": 2})
    oe = srv._run_enrichment("cmap_target_overexpression",
                             {"genes_up": ["AAA", "BBB"], "top_n": 2})

    assert "DRUG" in drug.upper() and "MOA" in drug.upper()
    assert "KNOCKDOWN" in kd.upper()
    assert "OVEREXPRES" in oe.upper() or "OVEREXPRESSED" in oe.upper()


def test_missing_dataset_reports_cleanly(tmp_path, monkeypatch):
    # Empty data dir -> no classes available.
    monkeypatch.setattr(srv, "registry", EngineRegistry(tmp_path))
    out = srv._run_enrichment("cmap_target_knockdown", {"genes_up": ["AAA"]})
    assert "not available" in out.lower()
