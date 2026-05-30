"""Integrity of the perturbation-class config — the single source of truth."""

from cmap_enrichment.perturbation_classes import (
    PERT_TYPE_TO_CLASS,
    PERTURBATION_CLASSES,
    TOOL_TO_CLASS,
    get_class,
)


def test_three_classes_present():
    assert set(PERTURBATION_CLASSES) == {"drug", "knockdown", "overexpression"}


def test_pert_types_are_disjoint_and_complete():
    seen = []
    for spec in PERTURBATION_CLASSES.values():
        seen += list(spec.pert_types)
    assert sorted(seen) == ["trt_cp", "trt_oe", "trt_sh", "trt_xpr"]
    assert len(seen) == len(set(seen)), "pert_types must map to exactly one class"


def test_class_mapping():
    assert PERT_TYPE_TO_CLASS["trt_cp"] == "drug"
    assert PERT_TYPE_TO_CLASS["trt_sh"] == "knockdown"
    assert PERT_TYPE_TO_CLASS["trt_xpr"] == "knockdown"
    assert PERT_TYPE_TO_CLASS["trt_oe"] == "overexpression"


def test_tool_names_unique_and_resolvable():
    assert set(TOOL_TO_CLASS) == {
        "cmap_drug_enrichment",
        "cmap_target_knockdown",
        "cmap_target_overexpression",
    }
    for tool, key in TOOL_TO_CLASS.items():
        assert get_class(tool).key == key


def test_only_drug_has_moa_only_knockdown_has_method():
    assert PERTURBATION_CLASSES["drug"].has_moa is True
    assert PERTURBATION_CLASSES["knockdown"].has_moa is False
    assert PERTURBATION_CLASSES["overexpression"].has_moa is False
    assert PERTURBATION_CLASSES["knockdown"].has_method is True
    assert PERTURBATION_CLASSES["knockdown"].methods == ("shRNA", "CRISPR")


def test_descriptions_carry_class_specific_framing():
    d = PERTURBATION_CLASSES["drug"].description.lower()
    k = PERTURBATION_CLASSES["knockdown"].description.lower()
    o = PERTURBATION_CLASSES["overexpression"].description.lower()
    assert "drug" in d and "loss-of-function" not in d
    assert "loss-of-function" in k and ("shrna" in k or "crispr" in k)
    assert "gain-of-function" in o


def test_get_class_rejects_unknown():
    import pytest

    with pytest.raises(KeyError):
        get_class("nonsense")
