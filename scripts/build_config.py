"""
Shared configuration for the build pipeline (scripts 01-03).

Maps each perturbation class to the LINCS 2020 Level-5 GCTX file(s) it is built
from and the metadata needed to tag its signatures. Select the class to build
with the CMAP_PERT_CLASS env var (drug | knockdown | overexpression); default
is 'drug'.
"""

from __future__ import annotations

import os

# pert_type -> (gctx filename, sub-method label or None)
PERT_TYPE_SOURCES: dict[str, tuple[str, str | None]] = {
    "trt_cp": ("level5_beta_trt_cp_n720216x12328.gctx", None),
    "trt_sh": ("level5_beta_trt_sh_n238351x12328.gctx", "shRNA"),
    "trt_xpr": ("level5_beta_trt_xpr_n142901x12328.gctx", "CRISPR"),
    "trt_oe": ("level5_beta_trt_oe_n34171x12328.gctx", None),
}

# class -> ordered list of pert_types it includes
CLASS_PERT_TYPES: dict[str, list[str]] = {
    "drug": ["trt_cp"],
    "knockdown": ["trt_sh", "trt_xpr"],
    "overexpression": ["trt_oe"],
}

# Which classes carry compound MoA/target annotation (drugs only).
CLASS_HAS_MOA = {"drug": True, "knockdown": False, "overexpression": False}


def selected_class() -> str:
    cls = os.environ.get("CMAP_PERT_CLASS", "drug").strip()
    if cls not in CLASS_PERT_TYPES:
        raise SystemExit(
            f"CMAP_PERT_CLASS={cls!r} invalid. Choose one of {list(CLASS_PERT_TYPES)}."
        )
    return cls
