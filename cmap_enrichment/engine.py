"""
Core enrichment engine.

Given an up/down gene signature, score every CMap/L1000 perturbation by how
strongly it REVERSES (negative connectivity) or MIMICS (positive connectivity)
that signature, using the Weighted Connectivity Score (WTCS) from Subramanian
et al. 2017.

A single CMapEngine instance is bound to ONE perturbation class
(drug | knockdown | overexpression). The same scoring math is shared across all
three classes; only the candidate set (filtered by pert_type via the per-class
data bundle) and the result framing differ. See ``perturbation_classes.py``.

Data layout (per-class, preferred)::

    <data_dir>/
        gene_names.txt                      # shared L1000 gene universe
        cmap.duckdb                         # optional, all classes
        drug/           {zscore_matrix.npy, rank_matrix.npy, metadata.parquet}
        knockdown/      {zscore_matrix.npy, rank_matrix.npy, metadata.parquet}
        overexpression/ {zscore_matrix.npy, rank_matrix.npy, metadata.parquet}

Legacy layout (flat, drug-only) is still supported: if
``<data_dir>/<class>/`` does not exist, the engine falls back to reading
``zscore_matrix.npy`` / ``rank_matrix.npy`` / ``signature_metadata.parquet``
directly from ``<data_dir>`` and treats them as the drug class.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .perturbation_classes import PERTURBATION_CLASSES, get_class

# Default: <repo_root>/data/processed
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"


class CMapEngine:
    """WTCS enrichment engine for one perturbation class."""

    def __init__(self, data_dir: str | Path | None = None, pert_class: str = "drug"):
        self.data_dir = Path(data_dir) if data_dir else DATA_DIR
        self.pert_class = get_class(pert_class).key
        self._spec = PERTURBATION_CLASSES[self.pert_class]

        class_dir = self.data_dir / self.pert_class
        legacy = not (class_dir / "rank_matrix.npy").exists()
        mat_dir = self.data_dir if legacy else class_dir

        if not (mat_dir / "rank_matrix.npy").exists():
            raise FileNotFoundError(
                f"No processed data for class '{self.pert_class}' under {self.data_dir}. "
                f"Looked in {class_dir}/ and {self.data_dir}/. "
                f"Run the build pipeline (scripts/01-03) for this perturbation type, "
                f"or download the bundle (scripts/download_data.sh)."
            )

        # Memory-map the big matrices -- never pull 3-6 GB into RAM.
        self.rank_matrix = np.load(mat_dir / "rank_matrix.npy", mmap_mode="r")
        self.zscore_matrix = np.load(mat_dir / "zscore_matrix.npy", mmap_mode="r")

        # Gene universe is shared; tolerate it living alongside the matrices.
        gene_file = self.data_dir / "gene_names.txt"
        if not gene_file.exists():
            gene_file = mat_dir / "gene_names.txt"
        with open(gene_file) as f:
            self.gene_names = [g.strip() for g in f if g.strip()]
        self.gene_to_idx = {g: i for i, g in enumerate(self.gene_names)}

        # Metadata: new per-class metadata.parquet or legacy signature_metadata.parquet
        meta_path = mat_dir / "metadata.parquet"
        if not meta_path.exists():
            meta_path = mat_dir / "signature_metadata.parquet"
        self.meta = self._normalize_meta(pd.read_parquet(meta_path))

        self.n_sigs, self.n_genes = self.rank_matrix.shape
        if len(self.meta) != self.n_sigs:
            raise ValueError(
                f"Metadata rows ({len(self.meta)}) != matrix rows ({self.n_sigs}) "
                f"for class '{self.pert_class}'."
            )

    # -- metadata normalization ------------------------------------------------

    def _normalize_meta(self, meta: pd.DataFrame) -> pd.DataFrame:
        """
        Guarantee the generic columns the engine relies on, deriving them from a
        legacy drug-only table when necessary so existing bundles keep working.
        """
        meta = meta.copy()

        # Canonical perturbagen label (compound name OR gene symbol).
        if "pert_name" not in meta.columns:
            if "cmap_name" in meta.columns:
                meta["pert_name"] = meta["cmap_name"]
            else:
                raise ValueError("Metadata lacks both 'pert_name' and 'cmap_name'.")
        # Keep cmap_name as an alias for backward compatibility.
        if "cmap_name" not in meta.columns:
            meta["cmap_name"] = meta["pert_name"]

        if "pert_class" not in meta.columns:
            meta["pert_class"] = self.pert_class
        if "pert_type" not in meta.columns:
            meta["pert_type"] = self._spec.pert_types[0]
        if "method" not in meta.columns:
            meta["method"] = pd.NA
        for col in ("moa", "target", "n_reps", "cell_iname"):
            if col not in meta.columns:
                meta[col] = pd.NA
        if "sig_idx" not in meta.columns:
            meta["sig_idx"] = np.arange(len(meta))

        return meta.reset_index(drop=True)

    # -- WTCS scoring (unchanged math) ----------------------------------------

    def _ks_enrichment(self, ranks_row: np.ndarray, gene_indices: np.ndarray) -> float:
        """KS-like enrichment score for one signature against one gene set."""
        n = len(ranks_row)
        n_set = len(gene_indices)
        if n_set == 0:
            return 0.0

        set_ranks = np.sort(ranks_row[gene_indices])
        positions = np.arange(1, n_set + 1)
        hit_scores = set_ranks / n
        expected = positions / n_set
        deviations_up = hit_scores - expected
        es = deviations_up[np.argmax(np.abs(deviations_up))]
        return float(es)

    def _wtcs(self, ranks_row: np.ndarray, up_idx: np.ndarray, down_idx: np.ndarray) -> float:
        """
        Weighted Connectivity Score.
        Negative = signature REVERSES the query; positive = signature MIMICS it.
        """
        es_up = self._ks_enrichment(ranks_row, up_idx) if len(up_idx) > 0 else 0.0
        es_down = self._ks_enrichment(ranks_row, down_idx) if len(down_idx) > 0 else 0.0

        if es_up == 0 and es_down == 0:
            return 0.0
        if np.sign(es_up) == np.sign(es_down):
            return 0.0
        return (es_up - es_down) / 2

    def _map_genes(self, genes: list[str]) -> tuple[np.ndarray, list[str]]:
        idx = np.array([self.gene_to_idx[g] for g in genes if g in self.gene_to_idx], dtype=int)
        unmapped = [g for g in genes if g not in self.gene_to_idx]
        return idx, unmapped

    def _candidate_mask(self, cell_line: Optional[str], method: Optional[str]) -> np.ndarray:
        """Row indices (into this class's matrices) that pass cell/method filters."""
        sel = pd.Series(True, index=self.meta.index)
        if cell_line:
            sel &= self.meta["cell_iname"] == cell_line
        if method:
            if not self._spec.has_method:
                raise ValueError(
                    f"Class '{self.pert_class}' has no sub-method to filter on."
                )
            sel &= self.meta["method"] == method
        mask = self.meta.index[sel].to_numpy()
        if len(mask) == 0:
            raise ValueError(
                f"No '{self.pert_class}' signatures match "
                f"cell_line={cell_line!r}, method={method!r}."
            )
        return mask

    def _score(self, up_idx: np.ndarray, down_idx: np.ndarray, mask: np.ndarray) -> np.ndarray:
        scores = np.zeros(len(mask), dtype=np.float64)
        use_down = len(down_idx) > 0
        for i, sig_idx in enumerate(mask):
            row = self.rank_matrix[sig_idx]
            if use_down:
                scores[i] = self._wtcs(row, up_idx, down_idx)
            else:
                scores[i] = self._ks_enrichment(row, up_idx)
        return scores

    # -- public query API ------------------------------------------------------

    def rank_perturbations(
        self,
        genes_up: list[str],
        genes_down: list[str] | None = None,
        cell_line: str | None = None,
        method: str | None = None,
        top_n: int = 50,
    ) -> dict:
        """
        Score every candidate perturbation in this class against the query and
        return the strongest reversers and mimickers in a single scoring pass.

        Returns a dict with:
            reversing : DataFrame, top_n most-negative WTCS (strongest reversal)
            mimicking : DataFrame, top_n most-positive WTCS (strongest mimic)
            stats     : dict of mapping diagnostics
        """
        up_idx, unmapped_up = self._map_genes(genes_up)
        down_idx, unmapped_down = self._map_genes(genes_down or [])

        if len(up_idx) == 0 and len(down_idx) == 0:
            raise ValueError(
                "No query genes found in the L1000 gene universe. "
                f"Unmapped: {unmapped_up + unmapped_down}"
            )

        mask = self._candidate_mask(cell_line, method)
        scores = self._score(up_idx, down_idx, mask)

        order = np.argsort(scores)  # ascending: most-negative (reversing) first
        rev_sel = mask[order[:top_n]]
        rev_scores = scores[order[:top_n]]
        mim_sel = mask[order[::-1][:top_n]]
        mim_scores = scores[order[::-1][:top_n]]

        stats = {
            "pert_class": self.pert_class,
            "n_up_mapped": len(up_idx),
            "n_down_mapped": len(down_idx),
            "n_up_unmapped": len(unmapped_up),
            "n_down_unmapped": len(unmapped_down),
            "unmapped_genes": unmapped_up + unmapped_down,
            "n_candidates": len(mask),
        }

        return {
            "reversing": self._assemble(rev_sel, rev_scores, stats, ascending=True),
            "mimicking": self._assemble(mim_sel, mim_scores, stats, ascending=False),
            "stats": stats,
        }

    def _assemble(self, sel: np.ndarray, scores: np.ndarray, stats: dict,
                  ascending: bool = True) -> pd.DataFrame:
        cols = ["pert_name", "cmap_name", "cell_iname", "n_reps", "pert_class", "pert_type"]
        if self._spec.has_moa:
            cols += ["moa", "target"]
        if self._spec.has_method:
            cols += ["method"]
        cols = [c for c in cols if c in self.meta.columns]

        out = self.meta.loc[sel, cols].copy()
        out["wtcs"] = scores
        out = out.sort_values("wtcs", ascending=ascending).reset_index(drop=True)
        for k, v in stats.items():
            out.attrs[k] = v
        return out

    def query_enrichment(
        self,
        genes_up: list[str],
        genes_down: list[str] | None = None,
        cell_line: str | None = None,
        method: str | None = None,
        top_n: int = 50,
    ) -> pd.DataFrame:
        """
        Backward-compatible helper: return only the top_n REVERSING perturbations
        (most-negative WTCS), with mapping diagnostics on ``df.attrs``.
        """
        res = self.rank_perturbations(
            genes_up, genes_down, cell_line=cell_line, method=method, top_n=top_n
        )
        return res["reversing"]

    # -- metadata helpers (parquet-driven; DuckDB not required) ----------------

    def get_perturbagen_info(self, name: str) -> pd.DataFrame:
        return self.meta[self.meta["pert_name"].str.casefold() == name.casefold()].copy()

    # Legacy alias.
    def get_compound_info(self, compound_name: str) -> pd.DataFrame:
        return self.get_perturbagen_info(compound_name)

    def get_signature(self, sig_idx: int) -> pd.Series:
        if sig_idx < 0 or sig_idx >= self.n_sigs:
            raise IndexError(sig_idx)
        return pd.Series(
            np.asarray(self.zscore_matrix[sig_idx]),
            index=self.gene_names,
            name=f"{self.pert_class}_sig_{sig_idx}",
        )

    def list_cell_lines(self) -> list[str]:
        return sorted(self.meta["cell_iname"].dropna().unique().tolist())

    def list_moas(self) -> list[str]:
        if "moa" not in self.meta.columns:
            return []
        return sorted(self.meta["moa"].dropna().unique().tolist())

    def search_perturbagens(self, query: str) -> pd.DataFrame:
        hits = self.meta[self.meta["pert_name"].str.contains(query, case=False, na=False)]
        cols = [c for c in ("pert_name", "moa", "target", "method") if c in hits.columns]
        return hits[cols].drop_duplicates().sort_values("pert_name").reset_index(drop=True)

    # Legacy alias.
    def search_compounds(self, query: str) -> pd.DataFrame:
        return self.search_perturbagens(query)
