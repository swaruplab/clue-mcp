"""
Lazy, per-class engine registry shared by the MCP server and REST API.

Each perturbation class (drug / knockdown / overexpression) is backed by its own
data bundle and its own CMapEngine. Engines are instantiated on first use so a
process that only ever queries drugs never pays to memory-map the knockdown
matrices. Missing bundles are reported clearly rather than crashing the server.
"""

from __future__ import annotations

import os
from pathlib import Path

from .engine import CMapEngine, DATA_DIR
from .perturbation_classes import PERTURBATION_CLASSES, get_class


def default_data_dir() -> Path:
    return Path(os.environ.get("CMAP_DATA_DIR", str(DATA_DIR)))


class EngineRegistry:
    def __init__(self, data_dir: str | Path | None = None):
        self.data_dir = Path(data_dir) if data_dir else default_data_dir()
        self._engines: dict[str, CMapEngine] = {}

    def available(self) -> dict[str, bool]:
        """Map each class key -> whether its processed data is present on disk."""
        out: dict[str, bool] = {}
        for key in PERTURBATION_CLASSES:
            class_dir = self.data_dir / key
            flat = self.data_dir  # legacy drug-only flat layout
            out[key] = (class_dir / "rank_matrix.npy").exists() or (
                key == "drug" and (flat / "rank_matrix.npy").exists()
            )
        return out

    def get(self, key_or_tool: str) -> CMapEngine:
        key = get_class(key_or_tool).key
        if key not in self._engines:
            self._engines[key] = CMapEngine(data_dir=self.data_dir, pert_class=key)
        return self._engines[key]
