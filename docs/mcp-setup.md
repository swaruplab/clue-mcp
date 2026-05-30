# MCP Server Setup

The repo ships an [MCP](https://modelcontextprotocol.io/) server that exposes the CMap enrichment engine as a set of tools any MCP-compatible LLM client (Claude Desktop, Claude Code, Cursor, etc.) can call directly.

---

## 1. Install

```bash
pip install -e ".[mcp]"
```

This pulls in the `mcp` Python SDK on top of the core engine dependencies.

## 2. Verify the server runs

```bash
CMAP_DATA_DIR=/abs/path/to/data/processed python mcp_server/server.py
```

The server speaks MCP over stdio, so it'll appear to hang — that's expected. Press `Ctrl-C` to quit. If you see a Python traceback, fix that before continuing.

## 3. Register with Claude Desktop

Edit your Claude Desktop config (`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS, `%APPDATA%\Claude\claude_desktop_config.json` on Windows) and add the server:

```json
{
  "mcpServers": {
    "clue-mcp": {
      "command": "/abs/path/to/conda/envs/model-ad_env/bin/python",
      "args": ["/abs/path/to/clue-mcp/mcp_server/server.py"],
      "env": {
        "CMAP_DATA_DIR": "/abs/path/to/clue-mcp/data/processed"
      }
    }
  }
}
```

A ready-to-edit template is in [`claude_desktop_config.example.json`](claude_desktop_config.example.json).

Restart Claude Desktop. The CMap tools will appear in the "tools" menu of any conversation.

## 4. Register with Claude Code

Add to `.claude/settings.json` (project) or `~/.claude/settings.json` (user):

```json
{
  "mcpServers": {
    "clue-mcp": {
      "command": "python",
      "args": ["mcp_server/server.py"],
      "env": {
        "CMAP_DATA_DIR": "data/processed"
      }
    }
  }
}
```

## 5. Available tools

The engine splits the LINCS 2020 catalog into three **perturbation classes**, each
exposed as its own enrichment tool with scientific framing tailored to that class.
This lets the LLM pick the right interpretation — a drug that *reverses* a disease
signature is a therapeutic lead, whereas a gene whose *knockdown* reverses it is a
candidate driver target.

| Tool | Perturbation class | Arguments | Returns |
|------|--------------------|-----------|---------|
| `cmap_drug_enrichment` | Small molecules (`trt_cp`) | `genes_up: list[str]`, `genes_down?`, `cell_line?`, `direction?`, `top_n?` | Ranked drugs; reversing = candidate therapeutics, mimicking = signature-inducers. Includes MoA / target. |
| `cmap_target_knockdown` | Loss-of-function (`trt_sh` shRNA + `trt_xpr` CRISPR) | `genes_up`, `genes_down?`, `cell_line?`, `method?` (`shRNA`/`CRISPR`), `direction?`, `top_n?` | Ranked knocked-down genes; reversing = candidate driver targets to suppress. |
| `cmap_target_overexpression` | Gain-of-function (`trt_oe`) | `genes_up`, `genes_down?`, `cell_line?`, `direction?`, `top_n?` | Ranked overexpressed genes; reversing = genes whose forced expression opposes the signature. |
| `cmap_list_perturbation_classes` | — | — | The three classes, their pert_types, and which datasets are installed. |
| `cmap_search_perturbagens` | any | `query: str`, `pert_class?` | Perturbagens (drugs/genes) whose name matches the partial query. |
| `cmap_list_cell_lines` | any | `pert_class?` | Cell lines available for a given class. |

`direction` accepts `reversing` (most-negative WTCS), `mimicking` (most-positive WTCS),
or `both` (default). Only `cmap_target_knockdown` accepts `method`, because shRNA vs.
CRISPR is meaningful only for loss-of-function.

## 6. Example conversation

> **You:** I have a gene signature from an Alzheimer's microglia dataset: up = APOE, CLU, TREM2, CD68, C1QB; down = SYN1, SNAP25, SLC17A7. What drugs might reverse it, and are there any candidate driver genes I should consider silencing?
>
> **Claude:** *(calls `cmap_drug_enrichment` with those gene lists)*
> The top drug reversal candidates in LINCS 2020 are: ...
> *(then calls `cmap_target_knockdown` with the same lists)*
> And the knockdowns that most strongly reverse the signature — i.e. candidate driver genes — are: ...

## 7. Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ModuleNotFoundError: cmap_enrichment` | `pip install -e .` not run, or wrong Python | Use an absolute path to the env's Python in the config |
| Server starts but tools never appear in Claude | JSON syntax error in the config | Validate with `python -m json.tool < claude_desktop_config.json` |
| `FileNotFoundError: rank_matrix.npy` | `CMAP_DATA_DIR` not set or wrong | Point it at the `data/processed` directory containing `gene_names.txt` and the per-class `drug/`, `knockdown/`, `overexpression/` subfolders |
| Tool replies "dataset not available" | That class's tarball isn't downloaded | Run `scripts/download_data.sh`; the drug class works first, knockdown/overexpression need their archives |
| Cold-start takes 15+ s | Matrices loading from disk on first call | Expected — the server stays warm after the first tool call |

---

## Next steps

- [Python library](python-tutorial.md) — same operations, direct API
- [REST API](api-tutorial.md) — same operations, HTTP
