"""Engine behavior on the synthetic per-class bundle."""

import numpy as np
import pytest

from cmap_enrichment.engine import CMapEngine


def test_loads_each_class(data_dir):
    for cls, n in [("drug", 3), ("knockdown", 3), ("overexpression", 2)]:
        eng = CMapEngine(data_dir=data_dir, pert_class=cls)
        assert eng.n_sigs == n
        assert eng.pert_class == cls
        # Matrices are memory-mapped, never loaded fully into RAM.
        assert isinstance(eng.rank_matrix, np.memmap)


def test_class_only_sees_its_own_signatures(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="knockdown")
    assert set(eng.meta["pert_type"]) <= {"trt_sh", "trt_xpr"}
    assert "trt_cp" not in set(eng.meta["pert_type"])


def test_reversing_vs_mimicking_direction(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    res = eng.rank_perturbations(genes_up=["AAA", "BBB"], genes_down=["EEE", "FFF"], top_n=3)
    # reverse_drug has the mirror profile -> strongest reverser (most negative).
    assert res["reversing"].iloc[0]["pert_name"] == "reverse_drug"
    assert res["reversing"].iloc[0]["wtcs"] < 0
    # mimic_drug shares the query profile -> strongest mimic (most positive).
    assert res["mimicking"].iloc[0]["pert_name"] == "mimic_drug"
    assert res["mimicking"].iloc[0]["wtcs"] > 0


def test_method_filter_on_knockdown(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="knockdown")
    res = eng.rank_perturbations(genes_up=["AAA", "BBB"], method="CRISPR", top_n=5)
    assert set(res["reversing"]["method"]) == {"CRISPR"}
    assert res["stats"]["n_candidates"] == 1


def test_method_filter_rejected_for_drug(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    with pytest.raises(ValueError):
        eng.rank_perturbations(genes_up=["AAA"], method="shRNA")


def test_cell_line_filter(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    res = eng.rank_perturbations(genes_up=["AAA", "BBB"], cell_line="MCF7", top_n=5)
    assert set(res["reversing"]["cell_iname"]) == {"MCF7"}


def test_unmapped_genes_reported(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    res = eng.rank_perturbations(genes_up=["AAA", "NOT_A_GENE"], top_n=3)
    assert res["stats"]["n_up_mapped"] == 1
    assert "NOT_A_GENE" in res["stats"]["unmapped_genes"]


def test_all_unmapped_raises(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    with pytest.raises(ValueError):
        eng.rank_perturbations(genes_up=["ZZZ"], genes_down=["QQQ"])


def test_drug_result_has_moa_knockdown_does_not(data_dir):
    drug = CMapEngine(data_dir=data_dir, pert_class="drug")
    kd = CMapEngine(data_dir=data_dir, pert_class="knockdown")
    drug_cols = drug.rank_perturbations(genes_up=["AAA"], top_n=2)["reversing"].columns
    kd_cols = kd.rank_perturbations(genes_up=["AAA"], top_n=2)["reversing"].columns
    assert "moa" in drug_cols and "method" not in drug_cols
    assert "method" in kd_cols and "moa" not in kd_cols


def test_legacy_flat_drug_bundle_still_loads(legacy_data_dir):
    eng = CMapEngine(data_dir=legacy_data_dir, pert_class="drug")
    assert eng.n_sigs == 2
    # Normalized columns are synthesized from the legacy schema.
    assert "pert_name" in eng.meta.columns
    assert "pert_class" in eng.meta.columns
    res = eng.rank_perturbations(genes_up=["AAA", "BBB"], top_n=2)
    assert res["reversing"].iloc[0]["pert_name"] == "old_drug_b"


def test_backward_compatible_query_enrichment(data_dir):
    eng = CMapEngine(data_dir=data_dir, pert_class="drug")
    df = eng.query_enrichment(genes_up=["AAA", "BBB"], top_n=3)
    assert df.iloc[0]["wtcs"] <= df.iloc[-1]["wtcs"]  # ascending (reversing first)
