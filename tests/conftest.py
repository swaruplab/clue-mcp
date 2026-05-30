"""Shared fixtures: build a tiny synthetic per-class data bundle on disk so the
engine can be exercised without the multi-GB Zenodo download."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

# Make the package importable when running pytest from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

GENES = ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF", "GGG", "HHH"]


def _write_class(root: Path, cls: str, rows: list[dict], zmat: np.ndarray) -> None:
    cdir = root / cls
    cdir.mkdir(parents=True, exist_ok=True)
    np.save(cdir / "zscore_matrix.npy", zmat.astype(np.float32))
    ranks = np.apply_along_axis(rankdata, 1, zmat).astype(np.int16)
    np.save(cdir / "rank_matrix.npy", ranks)
    meta = pd.DataFrame(rows)
    meta["sig_idx"] = range(len(meta))
    meta.to_parquet(cdir / "metadata.parquet", index=False)


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    """A synthetic data/processed tree with all three classes."""
    (tmp_path / "gene_names.txt").write_text("\n".join(GENES))

    n = len(GENES)
    # Signature A strongly up in AAA/BBB (so it mimics an AAA/BBB-up query);
    # signature B is its mirror (reverses that query).
    up_profile = np.array([3, 3, 0, 0, -1, -1, 0, 0], float)
    rng = np.random.default_rng(0)

    def mat(profiles):
        return np.vstack(profiles)

    _write_class(
        tmp_path, "drug",
        [
            {"pert_name": "mimic_drug", "cmap_name": "mimic_drug", "cell_iname": "A549",
             "pert_type": "trt_cp", "method": None, "n_reps": 2, "moa": "HDAC inhibitor", "target": "HDAC1"},
            {"pert_name": "reverse_drug", "cmap_name": "reverse_drug", "cell_iname": "MCF7",
             "pert_type": "trt_cp", "method": None, "n_reps": 1, "moa": "mTOR inhibitor", "target": "MTOR"},
            {"pert_name": "noise_drug", "cmap_name": "noise_drug", "cell_iname": "A549",
             "pert_type": "trt_cp", "method": None, "n_reps": 3, "moa": None, "target": None},
        ],
        mat([up_profile, -up_profile, rng.normal(size=n)]),
    )

    _write_class(
        tmp_path, "knockdown",
        [
            {"pert_name": "TP53", "cmap_name": "TP53", "cell_iname": "A549",
             "pert_type": "trt_sh", "method": "shRNA", "n_reps": 1, "moa": None, "target": None},
            {"pert_name": "TP53", "cmap_name": "TP53", "cell_iname": "A549",
             "pert_type": "trt_xpr", "method": "CRISPR", "n_reps": 1, "moa": None, "target": None},
            {"pert_name": "MYC", "cmap_name": "MYC", "cell_iname": "PC3",
             "pert_type": "trt_sh", "method": "shRNA", "n_reps": 2, "moa": None, "target": None},
        ],
        mat([-up_profile, -up_profile * 0.8, up_profile]),
    )

    _write_class(
        tmp_path, "overexpression",
        [
            {"pert_name": "EGFR", "cmap_name": "EGFR", "cell_iname": "A549",
             "pert_type": "trt_oe", "method": None, "n_reps": 1, "moa": None, "target": None},
            {"pert_name": "PTEN", "cmap_name": "PTEN", "cell_iname": "MCF7",
             "pert_type": "trt_oe", "method": None, "n_reps": 1, "moa": None, "target": None},
        ],
        mat([up_profile, -up_profile]),
    )

    return tmp_path


@pytest.fixture
def legacy_data_dir(tmp_path: Path) -> Path:
    """A flat, drug-only bundle using the OLD column/file names, to prove the
    engine still loads pre-three-class data."""
    (tmp_path / "gene_names.txt").write_text("\n".join(GENES))
    n = len(GENES)
    z = np.vstack([np.array([3, 3, 0, 0, -1, -1, 0, 0], float),
                   np.array([-3, -3, 0, 0, 1, 1, 0, 0], float)])
    np.save(tmp_path / "zscore_matrix.npy", z.astype(np.float32))
    np.save(tmp_path / "rank_matrix.npy",
            np.apply_along_axis(rankdata, 1, z).astype(np.int16))
    meta = pd.DataFrame([
        {"cmap_name": "old_drug_a", "cell_iname": "A549", "n_reps": 1,
         "moa": "x", "target": "y", "sig_idx": 0},
        {"cmap_name": "old_drug_b", "cell_iname": "MCF7", "n_reps": 1,
         "moa": None, "target": None, "sig_idx": 1},
    ])
    meta.to_parquet(tmp_path / "signature_metadata.parquet", index=False)
    return tmp_path
