"""Deeper RBP-centered APA and skipped-exon mechanism summaries.

This module builds on outputs from ``rbp_rna_processing_links.py``.  It avoids
recomputing CLIP overlaps and instead adds event-level annotations, phase-module
summaries, positional SE binding summaries, and candidate prioritization.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log10
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from scripts import rbp_rna_processing_links as rbp
except ModuleNotFoundError:  # pragma: no cover - supports direct script execution.
    import rbp_rna_processing_links as rbp


ANALYSIS = "rbp_rna_processing_links"
FIGURE_FORMATS = ("pdf", "png")
PALETTE = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "gold": "#E69F00",
    "sky": "#56B4E9",
    "gray": "#7A7A7A",
    "light_gray": "#E6E6E6",
    "dark": "#222222",
}


@dataclass(frozen=True)
class DeepDiveOutputs:
    annotated_apa: Path
    apa_utr_bias: Path
    phase_modules: Path
    se_region_distribution: Path
    se_region_bias: Path
    candidate_priority: Path
    manifest: Path


def resolve_project_root(project_root: str | os.PathLike[str] | None = None) -> Path:
    if project_root is not None:
        return Path(project_root).expanduser().resolve()
    workspace_root = Path(os.environ.get("WORKSPACE_ROOT", "/home/wdeng3/workspace/Codespace"))
    candidates = [
        workspace_root / "PDAC_Rhythmicity",
        Path.cwd() / "PDAC_Rhythmicity",
        Path.cwd(),
        Path(__file__).resolve().parents[1],
    ]
    for candidate in candidates:
        if candidate.name == "PDAC_Rhythmicity" and (candidate / "data").exists():
            return candidate.resolve()
    return Path(__file__).resolve().parents[1]


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt", ".bed"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _write_tsv(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def _relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _minus_log10(value: object, cap: float = 20.0) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric) or numeric <= 0:
        return 0.0
    return float(min(-log10(float(numeric)), cap))


def _safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _event_key(frame: pd.DataFrame) -> pd.Series:
    if "apa_event_id" in frame.columns and frame["apa_event_id"].notna().any():
        apa = frame["apa_event_id"].fillna("").astype(str)
    else:
        apa = pd.Series("", index=frame.index)
    if "event_id" in frame.columns:
        se = frame["event_id"].fillna("").astype(str)
    else:
        se = pd.Series("", index=frame.index)
    out = apa.where(apa.ne(""), se)
    return out.map(rbp.standardize_event_id)


def _bh_fdr(p_values: Iterable[object]) -> np.ndarray:
    values = pd.to_numeric(pd.Series(list(p_values)), errors="coerce").to_numpy(dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    mask = np.isfinite(values)
    if not mask.any():
        return out
    finite = values[mask]
    order = np.argsort(finite)
    ranked = finite[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q = np.empty_like(finite)
    q[order] = np.clip(adjusted, 0.0, 1.0)
    out[mask] = q
    return out


def _circular_mean_lag(lags: Iterable[object]) -> dict[str, float]:
    values = pd.to_numeric(pd.Series(list(lags)), errors="coerce").dropna()
    values = values[values.between(-12.0, 12.0)]
    if values.empty:
        return {"n": 0, "mean_lag_hours": np.nan, "median_lag_hours": np.nan, "resultant_length": np.nan}
    theta = ((values.to_numpy(dtype=float) % 24.0) / 24.0) * 2.0 * np.pi
    z = np.mean(np.exp(1j * theta))
    mean = (np.angle(z) / (2.0 * np.pi) * 24.0) % 24.0
    if mean > 12.0:
        mean -= 24.0
    return {
        "n": int(len(values)),
        "mean_lag_hours": float(mean),
        "median_lag_hours": float(values.median()),
        "resultant_length": float(abs(z)),
    }


def _lag_timing_class(mean_lag: float) -> str:
    if pd.isna(mean_lag):
        return "insufficient"
    abs_lag = abs(float(mean_lag))
    if abs_lag <= 2.0:
        return "near_synchronous"
    if abs_lag <= 6.0:
        return "short_lag"
    if abs_lag <= 10.0:
        return "long_lag"
    return "anti_phase"


def _lag_order(mean_lag: float) -> str:
    if pd.isna(mean_lag):
        return "insufficient"
    if abs(float(mean_lag)) <= 2.0:
        return "synchronous"
    return "event_after_rbp" if float(mean_lag) > 0 else "event_before_rbp"


def load_utr_annotation(project_root: Path) -> pd.DataFrame:
    path = project_root / "data" / "Ref" / "hg38_extracted_3UTR.bed"
    columns = [
        "transcript_id",
        "transcript_core",
        "target_gene_symbol",
        "chrom",
        "strand",
        "utr_length_bp",
        "utr_span_bp",
        "n_utr_intervals",
        "n_gene_utr_isoforms",
        "utr_length_rank",
        "utr_length_quantile",
        "utr_length_label",
    ]
    if not path.exists():
        return pd.DataFrame(columns=columns)
    bed = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
    parts = bed["name"].astype(str).str.split("|", expand=True)
    bed["transcript_id"] = parts[0].fillna("")
    bed["transcript_core"] = bed["transcript_id"].str.replace(r"\.\d+$", "", regex=True)
    bed["target_gene_symbol"] = parts[1].fillna("").map(rbp.standardize_symbol).replace("NA", "")
    bed["length"] = _safe_numeric(bed["end"]) - _safe_numeric(bed["start"])
    bed = bed[bed["transcript_id"].ne("") & bed["target_gene_symbol"].ne("") & bed["length"].gt(0)].copy()
    if bed.empty:
        return pd.DataFrame(columns=columns)
    agg = (
        bed.groupby(["transcript_id", "transcript_core", "target_gene_symbol", "chrom", "strand"], dropna=False)
        .agg(
            utr_length_bp=("length", "sum"),
            utr_start=("start", "min"),
            utr_end=("end", "max"),
            n_utr_intervals=("length", "size"),
        )
        .reset_index()
    )
    agg["utr_span_bp"] = _safe_numeric(agg["utr_end"]) - _safe_numeric(agg["utr_start"])
    agg = agg.drop(columns=["utr_start", "utr_end"])
    agg["n_gene_utr_isoforms"] = agg.groupby("target_gene_symbol")["transcript_id"].transform("nunique")
    agg["utr_length_rank"] = agg.groupby("target_gene_symbol")["utr_length_bp"].rank(method="average", pct=False)
    denom = (agg["n_gene_utr_isoforms"] - 1).replace(0, np.nan)
    agg["utr_length_quantile"] = (agg["utr_length_rank"] - 1) / denom
    agg.loc[agg["n_gene_utr_isoforms"].eq(1), "utr_length_quantile"] = 0.5
    agg["utr_length_label"] = "middle_utr"
    agg.loc[agg["n_gene_utr_isoforms"].eq(1), "utr_length_label"] = "single_annotated_utr"
    agg.loc[agg["n_gene_utr_isoforms"].gt(1) & agg["utr_length_quantile"].le(0.25), "utr_length_label"] = "short_utr_isoform"
    agg.loc[agg["n_gene_utr_isoforms"].gt(1) & agg["utr_length_quantile"].ge(0.75), "utr_length_label"] = "long_utr_isoform"
    return agg[columns].copy()


def annotate_apa_utr_rank(apa_phase: pd.DataFrame, utr: pd.DataFrame) -> pd.DataFrame:
    out = apa_phase.copy()
    if out.empty:
        return out
    out["transcript_id_clean"] = out.get("transcript_id", out.get("apa_event_id", "")).fillna("").astype(str).str.split("|").str[0]
    out["transcript_core"] = out["transcript_id_clean"].str.replace(r"\.\d+$", "", regex=True)
    if utr.empty:
        for col in [
            "utr_length_bp",
            "utr_span_bp",
            "n_utr_intervals",
            "n_gene_utr_isoforms",
            "utr_length_rank",
            "utr_length_quantile",
            "utr_length_label",
        ]:
            out[col] = np.nan if col != "utr_length_label" else "missing_utr_annotation"
        return out

    full_map = utr.drop_duplicates("transcript_id").set_index("transcript_id")
    core_map = utr.sort_values("utr_length_bp").drop_duplicates("transcript_core").set_index("transcript_core")
    merge_cols = [
        "target_gene_symbol",
        "chrom",
        "strand",
        "utr_length_bp",
        "utr_span_bp",
        "n_utr_intervals",
        "n_gene_utr_isoforms",
        "utr_length_rank",
        "utr_length_quantile",
        "utr_length_label",
    ]
    full = out["transcript_id_clean"].map(full_map[merge_cols].to_dict("index"))
    core = out["transcript_core"].map(core_map[merge_cols].to_dict("index"))
    rows = []
    for full_item, core_item in zip(full, core):
        rows.append(full_item if isinstance(full_item, dict) else (core_item if isinstance(core_item, dict) else {}))
    annot = pd.DataFrame(rows, index=out.index)
    for col in merge_cols:
        out[col] = annot[col] if col in annot else np.nan
    out["utr_annotation_available"] = out["utr_length_bp"].notna()
    out["utr_length_label"] = out["utr_length_label"].fillna("missing_utr_annotation")
    return out


def compute_apa_utr_bias(
    annotated_apa: pd.DataFrame,
    apa_lags: pd.DataFrame,
    n_permutations: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    if annotated_apa.empty or apa_lags.empty:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    apa = annotated_apa[annotated_apa["utr_annotation_available"].astype(bool)].copy()
    apa["event_key"] = apa["apa_event_id"].map(rbp.standardize_event_id)
    apa = apa.drop_duplicates(["comparison", "event_key"])
    lag = apa_lags.copy()
    lag["event_key"] = lag["apa_event_id"].map(rbp.standardize_event_id)
    lag = lag.drop_duplicates(["comparison", "rbp_symbol", "event_key"])

    rows = []
    for (comparison, rbp_symbol), grp in lag.groupby(["comparison", "rbp_symbol"], dropna=False):
        background = apa[apa["comparison"].eq(comparison)].copy()
        if background.empty:
            continue
        bound_keys = set(grp["event_key"])
        bound = background[background["event_key"].isin(bound_keys)].copy()
        if len(bound) < 8 or len(background) < 20:
            continue
        bound_q = _safe_numeric(bound["utr_length_quantile"]).dropna().to_numpy(dtype=float)
        bg_q = _safe_numeric(background["utr_length_quantile"]).dropna().to_numpy(dtype=float)
        bound_len = _safe_numeric(bound["utr_length_bp"]).dropna().to_numpy(dtype=float)
        bg_len = _safe_numeric(background["utr_length_bp"]).dropna().to_numpy(dtype=float)
        if len(bound_q) < 8 or len(bg_q) < 20:
            continue
        observed_delta = float(np.nanmedian(bound_q) - np.nanmedian(bg_q))
        if len(bg_q) >= len(bound_q):
            perm_idx = rng.integers(0, len(bg_q), size=(int(n_permutations), len(bound_q)))
            perm_delta = np.nanmedian(bg_q[perm_idx], axis=1) - np.nanmedian(bg_q)
            empirical_p = float((np.sum(np.abs(perm_delta) >= abs(observed_delta)) + 1) / (len(perm_delta) + 1))
        else:
            empirical_p = np.nan
        long_bound = float(np.mean(bound_q >= 0.75))
        long_bg = float(np.mean(bg_q >= 0.75))
        short_bound = float(np.mean(bound_q <= 0.25))
        short_bg = float(np.mean(bg_q <= 0.25))
        examples = (
            bound.sort_values(["utr_length_quantile", "fdr"], ascending=[False, True])
            .get("host_gene_symbol", pd.Series(dtype=str))
            .dropna()
            .astype(str)
            .head(8)
            .tolist()
        )
        rows.append(
            {
                "comparison": comparison,
                "rbp_symbol": rbp_symbol,
                "n_bound_apa_events": int(len(bound)),
                "n_background_apa_events": int(len(background)),
                "median_bound_utr_length_bp": float(np.nanmedian(bound_len)) if len(bound_len) else np.nan,
                "median_background_utr_length_bp": float(np.nanmedian(bg_len)) if len(bg_len) else np.nan,
                "median_bound_utr_quantile": float(np.nanmedian(bound_q)),
                "median_background_utr_quantile": float(np.nanmedian(bg_q)),
                "delta_median_utr_quantile": observed_delta,
                "long_utr_fraction_bound": long_bound,
                "long_utr_fraction_background": long_bg,
                "delta_long_utr_fraction": long_bound - long_bg,
                "short_utr_fraction_bound": short_bound,
                "short_utr_fraction_background": short_bg,
                "delta_short_utr_fraction": short_bound - short_bg,
                "empirical_p_value": empirical_p,
                "top_long_utr_examples": ",".join(dict.fromkeys(examples)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["empirical_fdr"] = _bh_fdr(out["empirical_p_value"])
    out["utr_bias_direction"] = "mixed_utr_rank"
    out.loc[out["delta_median_utr_quantile"].ge(0.10), "utr_bias_direction"] = "longer_utr_rank"
    out.loc[out["delta_median_utr_quantile"].le(-0.10), "utr_bias_direction"] = "shorter_utr_rank"
    return out.sort_values(
        ["empirical_fdr", "n_bound_apa_events", "delta_median_utr_quantile"],
        ascending=[True, False, False],
    )


def summarize_phase_modules(triplets: pd.DataFrame) -> pd.DataFrame:
    if triplets.empty:
        return pd.DataFrame()
    df = triplets.copy()
    df["event_key"] = _event_key(df)
    dedup_cols = ["comparison", "event_kind", "rbp_symbol", "event_key", "target_gene_id"]
    df = df.sort_values(["max_score", "n_accessions"], ascending=[False, False]).drop_duplicates(dedup_cols)
    rows = []
    for keys, grp in df.groupby(["comparison", "event_kind", "rbp_symbol"], dropna=False):
        comparison, event_kind, rbp_symbol = keys
        rbp_event = _circular_mean_lag(grp["rbp_to_event_lag_hours"])
        event_mrna = _circular_mean_lag(grp["event_to_mrna_lag_hours"])
        rbp_event_lag = _safe_numeric(grp["rbp_to_event_lag_hours"])
        event_mrna_lag = _safe_numeric(grp["event_to_mrna_lag_hours"])
        rows.append(
            {
                "comparison": comparison,
                "event_kind": event_kind,
                "rbp_symbol": rbp_symbol,
                "n_triplet_rows_deduplicated": int(len(grp)),
                "n_events": int(grp["event_key"].nunique()),
                "n_targets": int(grp["target_gene_id"].fillna("").replace("", np.nan).nunique()),
                "n_validation_replicated_events": int(
                    grp.loc[grp.get("validation_replicated", False).astype(bool), "event_key"].nunique()
                ),
                "rbp_to_event_mean_lag_hours": rbp_event["mean_lag_hours"],
                "rbp_to_event_median_lag_hours": rbp_event["median_lag_hours"],
                "rbp_to_event_resultant_length": rbp_event["resultant_length"],
                "event_to_mrna_mean_lag_hours": event_mrna["mean_lag_hours"],
                "event_to_mrna_median_lag_hours": event_mrna["median_lag_hours"],
                "event_to_mrna_resultant_length": event_mrna["resultant_length"],
                "rbp_event_same_phase_fraction": float(rbp_event_lag.abs().le(2.0).mean()),
                "rbp_event_antiphase_fraction": float((rbp_event_lag.abs() - 12.0).abs().le(2.0).mean()),
                "event_mrna_same_phase_fraction": float(event_mrna_lag.abs().le(2.0).mean()),
                "event_mrna_antiphase_fraction": float((event_mrna_lag.abs() - 12.0).abs().le(2.0).mean()),
                "median_clip_max_score": float(_safe_numeric(grp["max_score"]).median()) if "max_score" in grp else np.nan,
                "median_clip_n_accessions": float(_safe_numeric(grp["n_accessions"]).median()) if "n_accessions" in grp else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    out["rbp_event_timing_class"] = out["rbp_to_event_mean_lag_hours"].map(_lag_timing_class)
    out["rbp_event_order"] = out["rbp_to_event_mean_lag_hours"].map(_lag_order)
    out["event_mrna_timing_class"] = out["event_to_mrna_mean_lag_hours"].map(_lag_timing_class)
    out["event_mrna_order"] = out["event_to_mrna_mean_lag_hours"].map(
        lambda x: "mrna_after_event" if pd.notna(x) and x > 2 else ("mrna_before_event" if pd.notna(x) and x < -2 else "synchronous")
    )
    return out.sort_values(["event_kind", "n_events", "rbp_to_event_resultant_length"], ascending=[True, False, False])


def summarize_se_region_bias(se_lags: pd.DataFrame, clip_links: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if se_lags.empty:
        return pd.DataFrame(), pd.DataFrame()
    se = se_lags.copy()
    se["event_key"] = se["event_id"].astype(str).map(rbp.standardize_event_id)
    se["region_group"] = se["region"].astype(str).str.replace("SE_", "", regex=False)
    se["region_group"] = se["region_group"].replace(
        {
            "upstream_flank": "upstream_flank",
            "downstream_flank": "downstream_flank",
            "exon_body": "exon_body",
        }
    )
    se = se.drop_duplicates(["comparison", "rbp_symbol", "event_key", "region_group"])

    rows = []
    for keys, grp in se.groupby(["comparison", "rbp_symbol"], dropna=False):
        comparison, rbp_symbol = keys
        total = len(grp)
        for region in ["upstream_flank", "exon_body", "downstream_flank"]:
            count = int(grp["region_group"].eq(region).sum())
            rows.append(
                {
                    "comparison": comparison,
                    "rbp_symbol": rbp_symbol,
                    "region_group": region,
                    "n_region_links": count,
                    "n_total_region_links": total,
                    "region_fraction": float(count / total) if total else np.nan,
                }
            )
    dist = pd.DataFrame(rows)
    if dist.empty:
        return dist, pd.DataFrame()
    wide = dist.pivot_table(
        index=["comparison", "rbp_symbol"],
        columns="region_group",
        values="region_fraction",
        fill_value=0,
    ).reset_index()
    for col in ["upstream_flank", "exon_body", "downstream_flank"]:
        if col not in wide.columns:
            wide[col] = 0.0
    count_wide = dist.pivot_table(
        index=["comparison", "rbp_symbol"],
        columns="region_group",
        values="n_region_links",
        fill_value=0,
        aggfunc="sum",
    ).reset_index()
    count_wide["n_total_region_links"] = count_wide[[c for c in count_wide.columns if c not in {"comparison", "rbp_symbol"}]].sum(axis=1)
    summary = wide.merge(count_wide[["comparison", "rbp_symbol", "n_total_region_links"]], on=["comparison", "rbp_symbol"], how="left")
    summary["flank_fraction"] = summary["upstream_flank"] + summary["downstream_flank"]
    summary["dominant_region"] = summary[["upstream_flank", "exon_body", "downstream_flank"]].idxmax(axis=1)
    summary["dominant_region_fraction"] = summary[["upstream_flank", "exon_body", "downstream_flank"]].max(axis=1)
    summary = summary.sort_values(["n_total_region_links", "dominant_region_fraction"], ascending=[False, False])
    return dist, summary


def build_candidate_priority(
    candidate_summary: pd.DataFrame,
    phase_modules: pd.DataFrame,
    apa_utr_bias: pd.DataFrame,
    se_region_bias: pd.DataFrame,
    target_contrast: pd.DataFrame,
) -> pd.DataFrame:
    if phase_modules.empty:
        return pd.DataFrame()
    out = phase_modules.copy()
    if not candidate_summary.empty:
        out = out.merge(candidate_summary, on=["event_kind", "rbp_symbol"], how="left", suffixes=("", "_global"))
    if not target_contrast.empty:
        contrast = (
            target_contrast.groupby(["comparison", "event_kind", "rbp_symbol"], dropna=False)
            .agg(
                best_contrast_fdr=("fdr", "min"),
                max_rate_difference=("rate_difference", "max"),
                max_odds_ratio=("odds_ratio", "max"),
                max_phase_resultant_delta=("phase_resultant_delta", "max"),
            )
            .reset_index()
        )
        out = out.merge(contrast, on=["comparison", "event_kind", "rbp_symbol"], how="left")
    if not apa_utr_bias.empty:
        apa_bias = apa_utr_bias.copy()
        apa_bias["event_kind"] = "APA"
        out = out.merge(
            apa_bias[
                [
                    "comparison",
                    "event_kind",
                    "rbp_symbol",
                    "n_bound_apa_events",
                    "median_bound_utr_quantile",
                    "delta_median_utr_quantile",
                    "delta_long_utr_fraction",
                    "empirical_p_value",
                    "empirical_fdr",
                    "utr_bias_direction",
                    "top_long_utr_examples",
                ]
            ],
            on=["comparison", "event_kind", "rbp_symbol"],
            how="left",
        )
    if not se_region_bias.empty:
        se_bias = se_region_bias.copy()
        se_bias["event_kind"] = "SE"
        out = out.merge(
            se_bias[
                [
                    "comparison",
                    "event_kind",
                    "rbp_symbol",
                    "flank_fraction",
                    "dominant_region",
                    "dominant_region_fraction",
                    "n_total_region_links",
                ]
            ],
            on=["comparison", "event_kind", "rbp_symbol"],
            how="left",
        )

    out["enrichment_score"] = out.get("best_enrichment_fdr", np.nan).map(_minus_log10)
    out["contrast_score"] = out.get("best_contrast_fdr", np.nan).map(_minus_log10) if "best_contrast_fdr" in out else 0.0
    out["lag_coherence_score"] = out.get("best_lag_coherence_fdr", np.nan).map(_minus_log10)
    out["phase_module_score"] = _safe_numeric(out["rbp_to_event_resultant_length"]).fillna(0) + _safe_numeric(
        out["event_to_mrna_resultant_length"]
    ).fillna(0)
    out["replication_score"] = np.log10(_safe_numeric(out.get("n_validation_replicated_events", 0)).fillna(0) + 1.0)
    out["triplet_scale_score"] = np.log10(_safe_numeric(out["n_events"]).fillna(0) + 1.0)
    out["utr_or_region_score"] = 0.0
    if "delta_median_utr_quantile" in out:
        out["utr_or_region_score"] = out["utr_or_region_score"] + _safe_numeric(out["delta_median_utr_quantile"]).abs().fillna(0)
    if "dominant_region_fraction" in out:
        out["utr_or_region_score"] = out["utr_or_region_score"] + _safe_numeric(out["dominant_region_fraction"]).fillna(0)
    out["deep_priority_score"] = (
        out["enrichment_score"].clip(0, 10)
        + out["contrast_score"].clip(0, 10)
        + out["lag_coherence_score"].clip(0, 8)
        + 2.0 * out["phase_module_score"].clip(0, 2)
        + out["replication_score"].clip(0, 3)
        + out["triplet_scale_score"].clip(0, 4)
        + 2.0 * out["utr_or_region_score"].clip(0, 2)
    )
    out["mechanistic_axis"] = ""
    if "utr_bias_direction" in out:
        out.loc[out["event_kind"].eq("APA"), "mechanistic_axis"] = out.loc[out["event_kind"].eq("APA"), "utr_bias_direction"].fillna(
            "no_utr_rank_bias"
        )
    if "dominant_region" in out:
        out.loc[out["event_kind"].eq("SE"), "mechanistic_axis"] = (
            "SE_" + out.loc[out["event_kind"].eq("SE"), "dominant_region"].fillna("no_region_bias")
        )
    out["suggested_panel_use"] = "backup"
    out.loc[
        out.get("evidence_tier", "").astype(str).eq("main_candidate")
        | (out["deep_priority_score"].ge(out["deep_priority_score"].quantile(0.85))),
        "suggested_panel_use",
    ] = "main_or_supplement"
    out.loc[
        out["event_kind"].eq("SE") & out["deep_priority_score"].ge(out["deep_priority_score"].quantile(0.75)),
        "suggested_panel_use",
    ] = "splicing_supplement"
    return out.sort_values(["deep_priority_score", "n_events"], ascending=[False, False])


def _load_font(size: int):
    from PIL import ImageFont

    for path in [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial.ttf",
    ]:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def _save_png_pdf(image, figures_dir: Path, name: str) -> list[Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    png_path = figures_dir / f"{name}.png"
    pdf_path = figures_dir / f"{name}.pdf"
    image.save(png_path)
    try:
        from reportlab.lib.utils import ImageReader
        from reportlab.pdfgen import canvas

        width_px, height_px = image.size
        width_pt = width_px * 0.72
        height_pt = height_px * 0.72
        c = canvas.Canvas(str(pdf_path), pagesize=(width_pt, height_pt))
        c.drawImage(ImageReader(str(png_path)), 0, 0, width=width_pt, height=height_pt)
        c.showPage()
        c.save()
    except ModuleNotFoundError:
        image.save(pdf_path, "PDF", resolution=300.0)
    return [pdf_path, png_path]


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    color = color.lstrip("#")
    return tuple(int(color[i : i + 2], 16) for i in (0, 2, 4))


def _blend(c1: str, c2: str, t: float) -> tuple[int, int, int]:
    a = np.array(_hex_to_rgb(c1), dtype=float)
    b = np.array(_hex_to_rgb(c2), dtype=float)
    return tuple(np.clip(a * (1 - t) + b * t, 0, 255).astype(int))


def _draw_text(draw, xy, text: str, font, fill=None, anchor=None):
    draw.text(xy, str(text), font=font, fill=fill or _hex_to_rgb(PALETTE["dark"]), anchor=anchor)


def plot_apa_utr_bias(apa_utr_bias: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if apa_utr_bias.empty:
        return []
    from PIL import Image, ImageDraw

    plot = apa_utr_bias.copy()
    plot["abs_delta"] = _safe_numeric(plot["delta_median_utr_quantile"]).abs()
    plot = plot.sort_values(["empirical_fdr", "abs_delta", "n_bound_apa_events"], ascending=[True, False, False]).head(18)
    row_h = 34
    width = 1200
    height = 150 + row_h * len(plot)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = _load_font(18)
    small = _load_font(14)
    bold = _load_font(22)
    left = 310
    right = width - 85
    top = 95
    axis_y = top + row_h * len(plot) + 15
    _draw_text(draw, (30, 28), "APA events bound by rhythmic RBPs are stratified by annotated 3'UTR isoform rank", bold)
    _draw_text(draw, (30, 58), "Points show median within-gene 3'UTR length quantile: gray=all rhythmic APA, color=RBP-bound APA", small, PALETTE["gray"])
    draw.line((left, axis_y, right, axis_y), fill=_hex_to_rgb(PALETTE["dark"]), width=2)
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = int(left + tick * (right - left))
        draw.line((x, axis_y - 6, x, axis_y + 6), fill=_hex_to_rgb(PALETTE["dark"]), width=1)
        _draw_text(draw, (x, axis_y + 12), f"{tick:.2g}", small, anchor="ma")
    _draw_text(draw, ((left + right) // 2, axis_y + 44), "Within-gene 3'UTR length quantile", font, anchor="ma")
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        label = f"{row.comparison.replace('_vs_', ' vs ')}  {row.rbp_symbol}  n={row.n_bound_apa_events}"
        _draw_text(draw, (30, y - 7), label, small)
        bg = float(row.median_background_utr_quantile)
        bd = float(row.median_bound_utr_quantile)
        x_bg = int(left + np.clip(bg, 0, 1) * (right - left))
        x_bd = int(left + np.clip(bd, 0, 1) * (right - left))
        color = PALETTE["orange"] if bd >= bg else PALETTE["blue"]
        draw.line((x_bg, y, x_bd, y), fill=_hex_to_rgb(PALETTE["light_gray"]), width=5)
        draw.ellipse((x_bg - 6, y - 6, x_bg + 6, y + 6), outline=_hex_to_rgb(PALETTE["gray"]), width=2)
        draw.ellipse((x_bd - 8, y - 8, x_bd + 8, y + 8), fill=_hex_to_rgb(color))
        fdr = row.empirical_fdr
        fdr_label = f"FDR={fdr:.3g}" if pd.notna(fdr) else "FDR=NA"
        _draw_text(draw, (right + 15, y - 7), fdr_label, small, PALETTE["gray"])
    return _save_png_pdf(img, figures_dir, "rbp_apa_utr_rank_bias")


def plot_phase_modules(phase_modules: pd.DataFrame, candidate_priority: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if phase_modules.empty:
        return []
    from PIL import Image, ImageDraw

    if not candidate_priority.empty:
        keys = candidate_priority.head(24)[["comparison", "event_kind", "rbp_symbol"]]
        plot = phase_modules.merge(keys, on=["comparison", "event_kind", "rbp_symbol"], how="inner")
    else:
        plot = phase_modules.sort_values(["n_events", "rbp_to_event_resultant_length"], ascending=[False, False]).head(24)
    plot = plot.drop_duplicates(["comparison", "event_kind", "rbp_symbol"]).head(24)
    row_h = 34
    width = 1250
    height = 155 + row_h * len(plot)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = _load_font(18)
    small = _load_font(14)
    bold = _load_font(22)
    left = 370
    right = width - 125
    top = 98
    axis_y = top + row_h * len(plot) + 16
    _draw_text(draw, (30, 28), "RBP-event-mRNA phase modules", bold)
    _draw_text(draw, (30, 58), "Positive lag means the downstream layer peaks after the upstream layer", small, PALETTE["gray"])
    draw.line((left, axis_y, right, axis_y), fill=_hex_to_rgb(PALETTE["dark"]), width=2)
    for tick in [-12, -6, 0, 6, 12]:
        x = int(left + (tick + 12) / 24 * (right - left))
        draw.line((x, top - 18, x, axis_y + 6), fill=_hex_to_rgb(PALETTE["light_gray"] if tick else PALETTE["gray"]), width=1)
        _draw_text(draw, (x, axis_y + 12), str(tick), small, anchor="ma")
    _draw_text(draw, ((left + right) // 2, axis_y + 44), "Mean circular lag (hours)", font, anchor="ma")
    legend_x = right - 265
    draw.ellipse((legend_x, 28, legend_x + 14, 42), fill=_hex_to_rgb(PALETTE["orange"]))
    _draw_text(draw, (legend_x + 22, 25), "RBP -> event", small)
    draw.ellipse((legend_x + 135, 28, legend_x + 149, 42), fill=_hex_to_rgb(PALETTE["green"]))
    _draw_text(draw, (legend_x + 157, 25), "event -> mRNA", small)
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        label = f"{row.event_kind}  {row.rbp_symbol}  {row.comparison.replace('_vs_', ' vs ')}  events={row.n_events}"
        _draw_text(draw, (30, y - 7), label, small)
        lag1 = row.rbp_to_event_mean_lag_hours
        lag2 = row.event_to_mrna_mean_lag_hours
        x1 = int(left + (np.clip(float(lag1), -12, 12) + 12) / 24 * (right - left)) if pd.notna(lag1) else None
        x2 = int(left + (np.clip(float(lag2), -12, 12) + 12) / 24 * (right - left)) if pd.notna(lag2) else None
        if x1 is not None and x2 is not None:
            draw.line((x1, y, x2, y), fill=_hex_to_rgb(PALETTE["light_gray"]), width=3)
        if x1 is not None:
            draw.ellipse((x1 - 8, y - 8, x1 + 8, y + 8), fill=_hex_to_rgb(PALETTE["orange"]))
        if x2 is not None:
            draw.rectangle((x2 - 7, y - 7, x2 + 7, y + 7), fill=_hex_to_rgb(PALETTE["green"]))
    return _save_png_pdf(img, figures_dir, "rbp_phase_module_cascade")


def plot_se_region_distribution(se_region_distribution: pd.DataFrame, se_region_bias: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if se_region_distribution.empty or se_region_bias.empty:
        return []
    from PIL import Image, ImageDraw

    top_keys = se_region_bias.head(16)[["comparison", "rbp_symbol"]]
    plot = se_region_distribution.merge(top_keys, on=["comparison", "rbp_symbol"], how="inner")
    row_h = 34
    width = 1100
    height = 140 + row_h * len(top_keys)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    small = _load_font(14)
    font = _load_font(18)
    bold = _load_font(22)
    left = 310
    right = width - 80
    top = 88
    colors = {
        "upstream_flank": PALETTE["blue"],
        "exon_body": PALETTE["gold"],
        "downstream_flank": PALETTE["purple"],
    }
    _draw_text(draw, (30, 28), "CLIP positional support around rhythmic skipped exons", bold)
    legend_x = left
    for region in ["upstream_flank", "exon_body", "downstream_flank"]:
        draw.rectangle((legend_x, 56, legend_x + 18, 70), fill=_hex_to_rgb(colors[region]))
        _draw_text(draw, (legend_x + 25, 53), region.replace("_", " "), small)
        legend_x += 185
    for i, key in enumerate(top_keys.itertuples(index=False)):
        y = top + i * row_h
        sub = plot[plot["comparison"].eq(key.comparison) & plot["rbp_symbol"].eq(key.rbp_symbol)].set_index("region_group")
        label = f"{key.comparison.replace('_vs_', ' vs ')}  {key.rbp_symbol}"
        _draw_text(draw, (30, y - 7), label, small)
        x0 = left
        for region in ["upstream_flank", "exon_body", "downstream_flank"]:
            frac = float(sub.loc[region, "region_fraction"]) if region in sub.index else 0.0
            x1 = x0 + int(frac * (right - left))
            draw.rectangle((x0, y - 10, x1, y + 10), fill=_hex_to_rgb(colors[region]))
            x0 = x1
        draw.rectangle((left, y - 10, right, y + 10), outline=_hex_to_rgb(PALETTE["dark"]), width=1)
    _draw_text(draw, ((left + right) // 2, height - 35), "Fraction of rhythmic SE CLIP links", font, anchor="ma")
    return _save_png_pdf(img, figures_dir, "rbp_se_region_distribution")


def plot_candidate_matrix(candidate_priority: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if candidate_priority.empty:
        return []
    from PIL import Image, ImageDraw

    features = [
        "enrichment_score",
        "contrast_score",
        "lag_coherence_score",
        "phase_module_score",
        "replication_score",
        "triplet_scale_score",
        "utr_or_region_score",
    ]
    plot = candidate_priority.head(22).copy()
    width = 1180
    cell_w = 105
    row_h = 31
    left = 355
    top = 105
    height = top + row_h * len(plot) + 95
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    small = _load_font(13)
    font = _load_font(16)
    bold = _load_font(22)
    _draw_text(draw, (30, 28), "Deep RBP-processing candidate evidence matrix", bold)
    _draw_text(draw, (30, 58), "Rows are ranked by integrated enrichment, phase coherence, replication, and event-annotation evidence", small, PALETTE["gray"])
    labels = {
        "enrichment_score": "enrich",
        "contrast_score": "target bg",
        "lag_coherence_score": "lag coh",
        "phase_module_score": "phase R",
        "replication_score": "replicate",
        "triplet_scale_score": "events",
        "utr_or_region_score": "UTR/SE pos",
    }
    maxima = {}
    for feat in features:
        values = _safe_numeric(plot[feat]).fillna(0)
        maxima[feat] = max(float(values.quantile(0.95)), 1e-9)
        x = left + features.index(feat) * cell_w + cell_w // 2
        _draw_text(draw, (x, top - 28), labels[feat], small, anchor="ma")
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        label = f"{row.event_kind} {row.rbp_symbol} {row.comparison.replace('_vs_', ' vs ')}"
        _draw_text(draw, (30, y + 3), label, small)
        for j, feat in enumerate(features):
            value = getattr(row, feat)
            value = 0.0 if pd.isna(value) else float(value)
            t = float(np.clip(value / maxima[feat], 0, 1))
            color = _blend("#F3F3F3", PALETTE["blue"], t)
            x0 = left + j * cell_w
            draw.rectangle((x0, y - 10, x0 + cell_w - 5, y + 15), fill=color, outline=(255, 255, 255))
        score_x = left + len(features) * cell_w + 12
        _draw_text(draw, (score_x, y + 3), f"{row.deep_priority_score:.1f}", small, PALETTE["dark"])
    _draw_text(draw, (left + len(features) * cell_w + 12, top - 28), "score", small)
    return _save_png_pdf(img, figures_dir, "rbp_deep_candidate_evidence_matrix")


def make_figures(
    apa_utr_bias: pd.DataFrame,
    phase_modules: pd.DataFrame,
    se_region_distribution: pd.DataFrame,
    se_region_bias: pd.DataFrame,
    candidate_priority: pd.DataFrame,
    figures_dir: Path,
) -> list[Path]:
    saved: list[Path] = []
    try:
        saved.extend(plot_apa_utr_bias(apa_utr_bias, figures_dir))
        saved.extend(plot_phase_modules(phase_modules, candidate_priority, figures_dir))
        saved.extend(plot_se_region_distribution(se_region_distribution, se_region_bias, figures_dir))
        saved.extend(plot_candidate_matrix(candidate_priority, figures_dir))
    except Exception as exc:  # Keep table generation usable in lean notebook runtimes.
        print(f"Figure generation skipped: {type(exc).__name__}: {exc}")
    return saved


def run_deep_dive(
    project_root: str | os.PathLike[str] | None = None,
    analysis: str = ANALYSIS,
    n_permutations: int = 2000,
    make_plots: bool = True,
) -> dict[str, pd.DataFrame | list[Path] | DeepDiveOutputs]:
    project_root = resolve_project_root(project_root)
    output_dir = project_root / "data" / analysis
    figures_dir = project_root / "data" / "figures" / analysis
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    apa_phase = _read_table(output_dir / "apa_expression_phase_lag_inputs.tsv")
    apa_lags = _read_table(output_dir / "rbp_apa_phase_lags.tsv")
    se_lags = _read_table(output_dir / "rbp_se_phase_lags.tsv")
    triplets = _read_table(output_dir / "rbp_event_gene_triplets.tsv")
    clip_links = _read_table(output_dir / "clip_supported_rbp_target_links_consensus.tsv")
    candidate_summary = _read_table(output_dir / "rbp_main_figure_candidate_summary.tsv")
    target_contrast = _read_table(output_dir / "rbp_target_background_contrast.tsv")

    utr = load_utr_annotation(project_root)
    annotated_apa = annotate_apa_utr_rank(apa_phase, utr)
    apa_utr_bias = compute_apa_utr_bias(annotated_apa, apa_lags, n_permutations=n_permutations)
    phase_modules = summarize_phase_modules(triplets)
    se_region_distribution, se_region_bias = summarize_se_region_bias(se_lags, clip_links)
    candidate_priority = build_candidate_priority(
        candidate_summary=candidate_summary,
        phase_modules=phase_modules,
        apa_utr_bias=apa_utr_bias,
        se_region_bias=se_region_bias,
        target_contrast=target_contrast,
    )

    paths = DeepDiveOutputs(
        annotated_apa=_write_tsv(annotated_apa, output_dir / "rbp_apa_utr_rank_annotation.tsv"),
        apa_utr_bias=_write_tsv(apa_utr_bias, output_dir / "rbp_apa_utr_rank_bias.tsv"),
        phase_modules=_write_tsv(phase_modules, output_dir / "rbp_event_phase_module_summary.tsv"),
        se_region_distribution=_write_tsv(se_region_distribution, output_dir / "rbp_se_region_distribution.tsv"),
        se_region_bias=_write_tsv(se_region_bias, output_dir / "rbp_se_region_bias.tsv"),
        candidate_priority=_write_tsv(candidate_priority, output_dir / "rbp_deep_candidate_priority.tsv"),
        manifest=output_dir / "rbp_deep_dive_manifest.tsv",
    )

    saved_figures = []
    if make_plots:
        saved_figures = make_figures(
            apa_utr_bias=apa_utr_bias,
            phase_modules=phase_modules,
            se_region_distribution=se_region_distribution,
            se_region_bias=se_region_bias,
            candidate_priority=candidate_priority,
            figures_dir=figures_dir,
        )

    manifest = pd.DataFrame(
        [
            {"metric": "annotated_apa_rows", "value": len(annotated_apa), "detail": _relative_path(paths.annotated_apa, project_root)},
            {
                "metric": "annotated_apa_utr_available_fraction",
                "value": float(annotated_apa.get("utr_annotation_available", pd.Series(dtype=bool)).mean())
                if len(annotated_apa)
                else np.nan,
                "detail": "Fraction of APA phase rows with matching transcript-level 3'UTR annotation.",
            },
            {"metric": "apa_utr_bias_rows", "value": len(apa_utr_bias), "detail": _relative_path(paths.apa_utr_bias, project_root)},
            {"metric": "phase_module_rows", "value": len(phase_modules), "detail": _relative_path(paths.phase_modules, project_root)},
            {"metric": "se_region_bias_rows", "value": len(se_region_bias), "detail": _relative_path(paths.se_region_bias, project_root)},
            {"metric": "candidate_priority_rows", "value": len(candidate_priority), "detail": _relative_path(paths.candidate_priority, project_root)},
            {
                "metric": "saved_deep_dive_figures",
                "value": len(saved_figures),
                "detail": ";".join(_relative_path(Path(path), project_root) for path in saved_figures),
            },
            {
                "metric": "apa_utr_rank_caveat",
                "value": "annotation_proxy",
                "detail": "3'UTR directionality is inferred from annotated transcript-level 3'UTR length rank, not direct PAS coordinates.",
            },
        ]
    )
    _write_tsv(manifest, paths.manifest)

    return {
        "paths": paths,
        "annotated_apa": annotated_apa,
        "apa_utr_bias": apa_utr_bias,
        "phase_modules": phase_modules,
        "se_region_distribution": se_region_distribution,
        "se_region_bias": se_region_bias,
        "candidate_priority": candidate_priority,
        "manifest": manifest,
        "saved_figures": saved_figures,
    }


if __name__ == "__main__":
    results = run_deep_dive()
    manifest = results["manifest"]
    print(manifest.to_string(index=False))
