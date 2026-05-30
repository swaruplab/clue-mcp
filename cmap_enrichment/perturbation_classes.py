"""
Single source of truth for the three CMap/L1000 perturbation classes exposed by
clue-mcp. Each class wraps the *same* connectivity-scoring engine but carries its
own scientific framing so that an LLM (via the MCP tools) and a human reader pick
the correct biological interpretation:

    drug            trt_cp          small-molecule compounds      (pharmacology)
    knockdown       trt_sh, trt_xpr gene loss-of-function         (shRNA + CRISPR KO)
    overexpression  trt_oe          gene gain-of-function         (ORF over-expression)

The engine, MCP server, REST API, build scripts, and docs all read these
definitions so wording and pert_type filtering stay consistent in one place.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class PerturbationClass:
    # Identity
    key: str                      # internal key: "drug" | "knockdown" | "overexpression"
    tool_name: str                # MCP tool name / API path segment
    pert_types: tuple[str, ...]   # LINCS pert_type values that belong to this class

    # Vocabulary used to narrate results
    unit: str                     # what a single perturbagen is ("compound", "gene")
    perturbagen_singular: str     # "drug", "gene knockdown", "gene overexpression"
    perturbagen_plural: str       # "drugs", "gene knockdowns", "overexpressed genes"

    # Connectivity-score interpretation (negative WTCS = reversal)
    reversing_meaning: str        # what a top reversing hit implies biologically
    mimicking_meaning: str        # what a top mimicking hit implies biologically

    # LLM-facing tool description (when to pick this tool)
    description: str

    # Schema flags
    has_moa: bool = False         # MoA/target columns only meaningful for drugs
    has_method: bool = False      # sub-method (shRNA vs CRISPR) only for knockdown
    methods: tuple[str, ...] = field(default_factory=tuple)


PERTURBATION_CLASSES: dict[str, PerturbationClass] = {
    "drug": PerturbationClass(
        key="drug",
        tool_name="cmap_drug_enrichment",
        pert_types=("trt_cp",),
        unit="compound",
        perturbagen_singular="drug",
        perturbagen_plural="drugs",
        reversing_meaning=(
            "a small molecule whose transcriptional effect OPPOSES the query "
            "signature — a candidate to therapeutically reverse the condition "
            "(drug-repurposing lead)"
        ),
        mimicking_meaning=(
            "a small molecule that REPRODUCES the query signature — useful as a "
            "positive control or to flag a drug that may worsen/induce the state"
        ),
        description=(
            "Find SMALL-MOLECULE DRUGS whose transcriptional signature reverses (or "
            "mimics) a query gene signature, using the Connectivity Map (CMap/L1000) "
            "LINCS 2020 compound perturbations (trt_cp). Use this when the question is "
            "pharmacological: 'what drug could reverse this disease state?', 'drug "
            "repurposing candidates', 'compounds that phenocopy this signature'. "
            "Provide upregulated and (optionally) downregulated gene symbols. Returns "
            "compounds ranked by weighted connectivity score (most negative = strongest "
            "reversal = top therapeutic candidate), with mechanism-of-action and target "
            "annotation. NOT for genetic perturbations — use cmap_target_knockdown or "
            "cmap_target_overexpression for gene loss-/gain-of-function."
        ),
        has_moa=True,
    ),
    "knockdown": PerturbationClass(
        key="knockdown",
        tool_name="cmap_target_knockdown",
        pert_types=("trt_sh", "trt_xpr"),
        unit="gene",
        perturbagen_singular="gene knockdown",
        perturbagen_plural="gene knockdowns",
        reversing_meaning=(
            "knocking down this gene OPPOSES the query signature — the gene's activity "
            "may DRIVE/SUSTAIN the condition, making it a candidate therapeutic target "
            "(inhibition would be beneficial)"
        ),
        mimicking_meaning=(
            "knocking down this gene REPRODUCES the query signature — loss of this "
            "gene's function may CAUSE the state; the gene is likely protective/required"
        ),
        description=(
            "Find GENE LOSS-OF-FUNCTION perturbations whose signature reverses (or "
            "mimics) a query gene signature, using CMap/L1000 LINCS 2020 genetic "
            "knockdowns: shRNA (trt_sh) and CRISPR knockout (trt_xpr), pooled as "
            "loss-of-function. Use this when the question is about CAUSAL TARGETS or "
            "DEPENDENCIES: 'which gene, if inhibited/knocked out, would reverse this "
            "signature?', 'what target drives this state?'. Provide up/down gene "
            "symbols. Returns knocked-down GENES ranked by connectivity score (most "
            "negative = knockdown most strongly reverses the query = candidate target "
            "to inhibit). Optionally restrict to a single method (shRNA or CRISPR). "
            "This is GENE perturbation, not drugs — use cmap_drug_enrichment for "
            "compounds, cmap_target_overexpression for gain-of-function."
        ),
        has_method=True,
        methods=("shRNA", "CRISPR"),
    ),
    "overexpression": PerturbationClass(
        key="overexpression",
        tool_name="cmap_target_overexpression",
        pert_types=("trt_oe",),
        unit="gene",
        perturbagen_singular="gene overexpression",
        perturbagen_plural="overexpressed genes",
        reversing_meaning=(
            "over-expressing this gene OPPOSES the query signature — restoring/boosting "
            "this gene's activity may reverse the condition (the gene is a candidate to "
            "ACTIVATE/agonize, or is lost in the disease state)"
        ),
        mimicking_meaning=(
            "over-expressing this gene REPRODUCES the query signature — excess of this "
            "gene's activity may be SUFFICIENT to induce the state (candidate oncogene/"
            "driver to inhibit)"
        ),
        description=(
            "Find GENE GAIN-OF-FUNCTION perturbations whose signature reverses (or "
            "mimics) a query gene signature, using CMap/L1000 LINCS 2020 ORF "
            "over-expression (trt_oe). Use this when the question is about SUFFICIENCY "
            "or ACTIVATION: 'which gene, if over-expressed/activated, would reverse this "
            "signature?', 'what gain-of-function phenocopies this state?'. Provide "
            "up/down gene symbols. Returns over-expressed GENES ranked by connectivity "
            "score (most negative = over-expression most strongly reverses the query = "
            "candidate to activate; most positive = over-expression mimics/induces the "
            "state, candidate driver). This is GENE gain-of-function, not drugs and not "
            "knockdown — use cmap_drug_enrichment or cmap_target_knockdown otherwise."
        ),
    ),
}

# Convenience lookups -------------------------------------------------------

# tool_name -> class key
TOOL_TO_CLASS: dict[str, str] = {
    cls.tool_name: key for key, cls in PERTURBATION_CLASSES.items()
}

# pert_type -> class key  (e.g. "trt_sh" -> "knockdown")
PERT_TYPE_TO_CLASS: dict[str, str] = {
    pt: key for key, cls in PERTURBATION_CLASSES.items() for pt in cls.pert_types
}

# pert_type -> sub-method label (only knockdown distinguishes methods)
PERT_TYPE_TO_METHOD: dict[str, str] = {
    "trt_sh": "shRNA",
    "trt_xpr": "CRISPR",
}


def get_class(key_or_tool: str) -> PerturbationClass:
    """Resolve a PerturbationClass from either its internal key or its tool name."""
    if key_or_tool in PERTURBATION_CLASSES:
        return PERTURBATION_CLASSES[key_or_tool]
    if key_or_tool in TOOL_TO_CLASS:
        return PERTURBATION_CLASSES[TOOL_TO_CLASS[key_or_tool]]
    raise KeyError(
        f"Unknown perturbation class or tool: {key_or_tool!r}. "
        f"Valid keys: {list(PERTURBATION_CLASSES)}; "
        f"valid tools: {list(TOOL_TO_CLASS)}"
    )
