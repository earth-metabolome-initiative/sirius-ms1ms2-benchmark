"""
Summarize Sirius annotation results across positive and negative ionization
directories, and write a Markdown report with rich per-run and aggregate
statistics.

Expected directory layout (per ionization mode):
  <ionization>/
    <sample>/
      sirius_normal/sirius/annotation_results.csv.gz
      sirius_ms2_only/sirius/annotation_results.csv.gz

Usage:
  python summarize_results.py
  python summarize_results.py --positive-dir positive --negative-dir negative \
                              --output summary.md
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
import uuid

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

IONIZATION_DIRS: dict[str, str] = {
    "positive": "positive",
    "negative": "negative",
}
SIRIUS_MODES: list[str] = ["sirius_normal", "sirius_ms2_only"]
ANNOTATION_FILE: str = "sirius/annotation_results.csv.gz"

# Columns we always expect
REQUIRED_COLUMNS: list[str] = [
    "inchiKey",
    "smiles",
    "xlogP",
    "rank",
    "csiScore",
    "tanimotoSimilarity",
    "mcesDistToTopHit",
    "molecularFormula",
    "adduct",
    "formulaId",
    "feature_id",
    "npc_pathway",
    "npc_superclass",
    "npc_class",
]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def load_annotation_file(path: Path) -> pd.DataFrame | None:
    """Read a gzipped CSV; return None if missing or empty."""
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            df = pd.read_csv(fh)
        if df.empty:
            return None
        # Keep only expected columns that are actually present
        present = [c for c in REQUIRED_COLUMNS if c in df.columns]
        return df[present]
    except Exception as exc:  # noqa: BLE001
        print(f"  WARNING: could not read {path}: {exc}", file=sys.stderr)
        return None


def collect_mode_data(
    ionization_dir: Path,
    mode: str,
) -> pd.DataFrame:
    """
    Walk *ionization_dir* and concatenate every annotation CSV for *mode*.

    Adds columns ``sample`` and ``ionization`` for traceability.
    """
    frames: list[pd.DataFrame] = []
    for sample_dir in sorted(ionization_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        annotation_path = sample_dir / mode / ANNOTATION_FILE
        df = load_annotation_file(annotation_path)
        if df is None:
            continue
        df = df.copy()
        df["sample"] = sample_dir.name
        df["ionization"] = ionization_dir.name

        # we generate a uuid per feature_id
        feature_to_uuid = {fid: str(uuid.uuid4()) for fid in df["feature_id"].unique()}
        df["uuid"] = df["feature_id"].map(feature_to_uuid)
        frames.append(df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def _safe_top1(df: pd.DataFrame) -> pd.DataFrame:
    """Return only rank-1 hits (best candidate per spectrum)."""
    if "rank" not in df.columns or df.empty:
        return df
    return df[df["rank"] == 1]


def compute_stats(df: pd.DataFrame, label: str) -> dict[str, object]:
    """
    Derive a rich set of statistics from an annotation DataFrame.

    Parameters
    ----------
    df:
        Full (all-ranks) annotation table for one mode + ionization.
    label:
        Human-readable label used only for the ``label`` key in the result.
    """
    stats: dict[str, object] = {"label": label}

    if df.empty:
        stats["n_samples"] = 0
        return stats

    top1 = _safe_top1(df)

    # --- Coverage ---
    stats["n_samples"] = int(df["sample"].nunique()) if "sample" in df.columns else "—"

    # Spectra that received at least one annotation
    annotated_spectra = (
        int(df[df["inchiKey"].notna()]["feature_id"].nunique())
        if "inchiKey" in df.columns
        else "—"
    )
    stats["n_annotated_spectra"] = annotated_spectra


    # --- Candidates per spectrum ---
    if "feature_id" in df.columns:
        cands_per_spectrum = df.groupby("uuid").size()
        stats["avg_candidates_per_spectrum"] = round(cands_per_spectrum.mean(), 2)
        stats["median_candidates_per_spectrum"] = round(cands_per_spectrum.median(), 2)
        stats["max_candidates_per_spectrum"] = int(cands_per_spectrum.max())
        stats["min_candidates_per_spectrum"] = int(cands_per_spectrum.min())

    # --- csiScore (all ranks) ---
    if "csiScore" in df.columns:
        scores = df["csiScore"].dropna()
        stats["avg_csiScore_all"] = round(scores.mean(), 4) if not scores.empty else "—"
        stats["median_csiScore_all"] = round(scores.median(), 4) if not scores.empty else "—"
        stats["std_csiScore_all"] = round(scores.std(), 4) if not scores.empty else "—"

    # --- csiScore (top-1 only) ---
    if "csiScore" in top1.columns:
        scores_top1 = top1["csiScore"].dropna()
        stats["avg_csiScore_top1"] = round(scores_top1.mean(), 4) if not scores_top1.empty else "—"
        stats["median_csiScore_top1"] = round(scores_top1.median(), 4) if not scores_top1.empty else "—"
        stats["std_csiScore_top1"] = round(scores_top1.std(), 4) if not scores_top1.empty else "—"

    # --- Tanimoto similarity (top-1) ---
    if "tanimotoSimilarity" in top1.columns:
        tan = top1["tanimotoSimilarity"].dropna()
        stats["avg_tanimoto_top1"] = round(tan.mean(), 4) if not tan.empty else "—"
        stats["median_tanimoto_top1"] = round(tan.median(), 4) if not tan.empty else "—"
        # Fraction with tanimoto ≥ 0.5 (good structural similarity)
        if not tan.empty:
            stats["pct_tanimoto_ge_0.5"] = round(100 * (tan >= 0.5).mean(), 1)
        else:
            stats["pct_tanimoto_ge_0.5"] = "—"

    # --- MCES distance to top hit (top-1) ---
    if "mcesDistToTopHit" in top1.columns:
        mces = top1["mcesDistToTopHit"].dropna()
        stats["avg_mcesDistToTopHit_top1"] = round(mces.mean(), 4) if not mces.empty else "—"
        stats["median_mcesDistToTopHit_top1"] = round(mces.median(), 4) if not mces.empty else "—"

    # --- xlogP (top-1) ---
    if "xlogP" in top1.columns:
        logp = top1["xlogP"].dropna()
        stats["avg_xlogP_top1"] = round(logp.mean(), 2) if not logp.empty else "—"
        stats["median_xlogP_top1"] = round(logp.median(), 2) if not logp.empty else "—"
        stats["std_xlogP_top1"] = round(logp.std(), 2) if not logp.empty else "—"

    # --- Unique structures ---
    if "inchiKey" in top1.columns:
        stats["n_unique_inchikeys_top1"] = int(top1["inchiKey"].dropna().nunique())
    if "molecularFormula" in top1.columns:
        stats["n_unique_formulas_top1"] = int(top1["molecularFormula"].dropna().nunique())
    if "smiles" in top1.columns:
        stats["n_unique_smiles_top1"] = int(top1["smiles"].dropna().nunique())

    # --- Adduct diversity ---
    if "adduct" in top1.columns:
        stats["n_unique_adducts"] = int(top1["adduct"].dropna().nunique())
        most_common_adduct = top1["adduct"].dropna().mode()
        stats["most_common_adduct"] = most_common_adduct.iloc[0] if not most_common_adduct.empty else "—"

    # --- NPC classifications (top-1, unique spectra) ---
    for npc_col in ("npc_pathway", "npc_superclass", "npc_class"):
        if npc_col in top1.columns:
            n_classified = int(top1[npc_col].dropna().shape[0])
            stats[f"n_classified_{npc_col}"] = n_classified
            n_unique = int(top1[npc_col].dropna().nunique())
            stats[f"n_unique_{npc_col}"] = n_unique
            top_val = top1[npc_col].dropna().mode()
            stats[f"top_{npc_col}"] = top_val.iloc[0] if not top_val.empty else "—"

    # --- Per-sample uniformity: std of annotation counts across samples ---
    if "sample" in df.columns and "feature_id" in df.columns:
        per_sample_counts = df.groupby("sample")["feature_id"].nunique()
        stats["avg_spectra_per_sample"] = round(per_sample_counts.mean(), 1)
        stats["std_spectra_per_sample"] = round(per_sample_counts.std(), 1)

    return stats


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

# Human-readable labels for every stat key
STAT_LABELS: dict[str, str] = {
    "n_samples": "Number of samples",
    "n_annotated_spectra": "Annotated spectra",
    "annotation_rate_pct": "Annotation rate (%)",
    "avg_candidates_per_spectrum": "Avg candidates / spectrum",
    "median_candidates_per_spectrum": "Median candidates / spectrum",
    "max_candidates_per_spectrum": "Max candidates / spectrum",
    "min_candidates_per_spectrum": "Min candidates / spectrum",
    "avg_csiScore_all": "Avg CSI:FingerID score (all ranks)",
    "median_csiScore_all": "Median CSI:FingerID score (all ranks)",
    "std_csiScore_all": "Std CSI:FingerID score (all ranks)",
    "avg_csiScore_top1": "Avg CSI:FingerID score (rank 1)",
    "median_csiScore_top1": "Median CSI:FingerID score (rank 1)",
    "std_csiScore_top1": "Std CSI:FingerID score (rank 1)",
    "avg_tanimoto_top1": "Avg Tanimoto similarity (rank 1)",
    "median_tanimoto_top1": "Median Tanimoto similarity (rank 1)",
    "pct_tanimoto_ge_0.5": "% hits with Tanimoto ≥ 0.5",
    "avg_mcesDistToTopHit_top1": "Avg MCES dist to top hit (rank 1)",
    "median_mcesDistToTopHit_top1": "Median MCES dist to top hit (rank 1)",
    "avg_xlogP_top1": "Avg xlogP (rank 1)",
    "median_xlogP_top1": "Median xlogP (rank 1)",
    "std_xlogP_top1": "Std xlogP (rank 1)",
    "n_unique_inchikeys_top1": "Unique InChIKeys (rank 1)",
    "n_unique_formulas_top1": "Unique molecular formulas (rank 1)",
    "n_unique_smiles_top1": "Unique SMILES (rank 1)",
    "n_unique_adducts": "Unique adducts",
    "most_common_adduct": "Most common adduct",
    "n_classified_npc_pathway": "Spectra with NPC pathway",
    "n_unique_npc_pathway": "Unique NPC pathways",
    "top_npc_pathway": "Most frequent NPC pathway",
    "n_classified_npc_superclass": "Spectra with NPC superclass",
    "n_unique_npc_superclass": "Unique NPC superclasses",
    "top_npc_superclass": "Most frequent NPC superclass",
    "n_classified_npc_class": "Spectra with NPC class",
    "n_unique_npc_class": "Unique NPC classes",
    "top_npc_class": "Most frequent NPC class",
    "avg_spectra_per_sample": "Avg spectra per sample",
    "std_spectra_per_sample": "Std spectra per sample",
}

# Logical groupings for the table sections
STAT_GROUPS: list[tuple[str, list[str]]] = [
    (
        "Coverage",
        [
            "n_samples",
            "n_annotated_spectra",
            "annotation_rate_pct",
            "avg_spectra_per_sample",
            "std_spectra_per_sample",
        ],
    ),
    (
        "Candidates per spectrum",
        [
            "avg_candidates_per_spectrum",
            "median_candidates_per_spectrum",
            "max_candidates_per_spectrum",
            "min_candidates_per_spectrum",
        ],
    ),
    (
        "CSI:FingerID score",
        [
            "avg_csiScore_top1",
            "median_csiScore_top1",
            "std_csiScore_top1",
            "avg_csiScore_all",
            "median_csiScore_all",
            "std_csiScore_all",
        ],
    ),
    (
        "Structural similarity",
        [
            "avg_tanimoto_top1",
            "median_tanimoto_top1",
            "pct_tanimoto_ge_0.5",
            "avg_mcesDistToTopHit_top1",
            "median_mcesDistToTopHit_top1",
        ],
    ),
    (
        "Physicochemical properties",
        [
            "avg_xlogP_top1",
            "median_xlogP_top1",
            "std_xlogP_top1",
        ],
    ),
    (
        "Chemical diversity",
        [
            "n_unique_inchikeys_top1",
            "n_unique_formulas_top1",
            "n_unique_smiles_top1",
            "n_unique_adducts",
            "most_common_adduct",
        ],
    ),
    (
        "NPC chemical classification (rank 1)",
        [
            "n_classified_npc_pathway",
            "n_unique_npc_pathway",
            "top_npc_pathway",
            "n_classified_npc_superclass",
            "n_unique_npc_superclass",
            "top_npc_superclass",
            "n_classified_npc_class",
            "n_unique_npc_class",
            "top_npc_class",
        ],
    ),
]


def build_markdown(all_stats: list[dict[str, object]]) -> str:
    """
    Render *all_stats* (one dict per combination of ionization × mode) as a
    grouped Markdown document.
    """
    lines: list[str] = []

    lines.append("# Sirius Annotation Summary\n")
    lines.append(
        "Statistics are computed across **all samples** in each ionization directory.  \n"
        "**Rank 1** = best candidate per spectrum; **all ranks** = every candidate returned.\n"
    )

    # Column headers derived from the labels of each stats dict
    col_labels: list[str] = [s["label"] for s in all_stats]  # type: ignore[index]
    header = "| Metric | " + " | ".join(col_labels) + " |"
    separator = "| :--- | " + " | ".join([":---:"] * len(all_stats)) + " |"

    for group_name, stat_keys in STAT_GROUPS:
        # Check at least one stat key is present in any column
        available = [k for k in stat_keys if any(k in s for s in all_stats)]
        if not available:
            continue

        lines.append(f"\n## {group_name}\n")
        lines.append(header)
        lines.append(separator)

        for key in available:
            row_label = STAT_LABELS.get(key, key)
            values: list[str] = []
            for stats in all_stats:
                val = stats.get(key, "—")
                values.append("—" if val is None else str(val))
            lines.append(f"| {row_label} | " + " | ".join(values) + " |")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize Sirius annotation results and write a Markdown report."
    )
    parser.add_argument(
        "--positive-dir",
        type=Path,
        default=Path("positive"),
        help="Root directory containing positive-mode sample folders (default: positive/).",
    )
    parser.add_argument(
        "--negative-dir",
        type=Path,
        default=Path("negative"),
        help="Root directory containing negative-mode sample folders (default: negative/).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("summary.md"),
        help="Output Markdown file path (default: summary.md).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    ionization_dirs: dict[str, Path] = {
        "positive": args.positive_dir,
        "negative": args.negative_dir,
    }

    all_stats: list[dict[str, object]] = []

    for ion_name, ion_dir in ionization_dirs.items():
        if not ion_dir.exists():
            print(f"WARNING: {ion_dir} does not exist — skipping.", file=sys.stderr)
            continue

        for mode in SIRIUS_MODES:
            mode_label = mode.replace("sirius_", "").replace("_", " ")
            label = f"{ion_name.capitalize()} — {mode_label}"
            print(f"Loading data for: {label} …")

            df = collect_mode_data(ion_dir, mode)
            stats = compute_stats(df, label=label)
            all_stats.append(stats)
            print(
                f"  → {stats.get('n_samples', 0)} sample(s), "
                f"{stats.get('n_annotated_spectra', '—')} annotated."
            )

    if not all_stats:
        print("ERROR: no data loaded — nothing to summarize.", file=sys.stderr)
        sys.exit(1)

    md = build_markdown(all_stats)
    args.output.write_text(md, encoding="utf-8")
    print(f"\nSummary written to {args.output}")


if __name__ == "__main__":
    main()