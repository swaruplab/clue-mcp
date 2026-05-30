"""
Minimal example: query the CMap engine across all three perturbation classes
to interpret a disease-like gene signature.

  - drugs           -> candidate therapeutics that reverse the signature
  - knockdowns      -> candidate driver genes to silence
  - overexpression  -> genes whose forced expression opposes the signature

Run from the repo root:
    python examples/quickstart_library.py
"""

from pathlib import Path

from cmap_enrichment import EngineRegistry

# Point at data/processed/ relative to the repo root.
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "processed"

# Toy AD microglia signature:
#   up   = neuroinflammation / lipid handling markers
#   down = synaptic / glutamatergic markers
GENES_UP = ["APOE", "CLU", "TREM2", "CD68", "C1QB", "C1QA", "CTSB"]
GENES_DOWN = ["SYN1", "SNAP25", "SLC17A7", "GAP43"]


def main() -> None:
    registry = EngineRegistry(data_dir=DATA_DIR)
    available = registry.available()
    print(f"Data dir: {DATA_DIR}")
    print(f"Installed classes: {[k for k, ok in available.items() if ok] or 'none'}\n")

    print(f"Query: {len(GENES_UP)} up genes, {len(GENES_DOWN)} down genes")

    for pert_class in ("drug", "knockdown", "overexpression"):
        print("\n" + "=" * 90)
        print(f"CLASS: {pert_class}")
        if not available.get(pert_class):
            print("  dataset not installed — run scripts/download_data.sh")
            continue

        engine = registry.get(pert_class)
        res = engine.rank_perturbations(
            genes_up=GENES_UP,
            genes_down=GENES_DOWN,
            top_n=10,
        )
        stats = res["stats"]
        print(
            f"  {engine.n_sigs:,} signatures; mapped {stats['n_up_mapped']} up / "
            f"{stats['n_down_mapped']} down to L1000"
        )
        if stats["unmapped_genes"]:
            print(f"  Unmapped: {stats['unmapped_genes']}")

        reversing = res["reversing"]
        cols = [c for c in ("pert_name", "cell_iname", "method", "wtcs", "moa", "target")
                if c in reversing.columns]
        print("\n  Top reversers (most-negative WTCS first):")
        print(reversing[cols].head(10).to_string(index=False))

    print("\nNext steps:")
    print("  - For permutation FDR + plots:  python analysis/run_enrichment.py <gene_file>")
    print("  - For REST API:                 uvicorn api.main:app --port 8000")
    print("  - For Claude/LLM access:        see docs/mcp-setup.md")


if __name__ == "__main__":
    main()
