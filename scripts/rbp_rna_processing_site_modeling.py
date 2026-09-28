"""Regulatory-site availability and sample-level usage models for RBP modules.

This module extends the RBP-centered APA/SE analysis with two mechanistic
checks:

1. APA regulatory-site availability proxies from annotated 3'UTR length rank,
   TargetScan miRNA site burden, and broad 3'UTR CLIP burden.
2. Exploratory sample-level PDUI/PSI models asking whether RBP expression adds
   explanatory power beyond host-gene expression in a within-event framework.

The sample-level models are prioritization screens. Observations share samples
and events, so p-values should not be read as fully independent tests.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

try:
    from scripts import rbp_rna_processing_links as rbp
    from scripts import rbp_rna_processing_deep_dive as rbp_deep
except ModuleNotFoundError:  # pragma: no cover - supports direct script execution.
    import rbp_rna_processing_links as rbp
    import rbp_rna_processing_deep_dive as rbp_deep


ANALYSIS = "rbp_rna_processing_links"
FOCUS_RBPS = {"DDX3X", "MBNL2", "RBFOX2", "PRPF8", "NUDT21", "CPSF6", "FIP1L1"}

APA_USAGE_FILES = {
    "TCGA_Tumor": "data/polyAPA/TCGA/filtered_apa/TCGA_variable_apa_events.csv",
    "CPTAC_Tumor": "data/polyAPA/CPTAC/filtered_apa/CPTAC_variable_apa_events.csv",
    "GTEx_Normal": "data/polyAPA/GTEx/filtered_apa/GTEx_variable_apa_events.csv",
}
SE_USAGE_FILES = {
    "TCGA_Tumor": "data/AS/TCGA/filtered_events/TCGA_SE_variable_AS.csv",
    "CPTAC_Tumor": "data/AS/CPTAC/filtered_events/CPTAC_SE_variable_AS.csv",
    "GTEx_Normal": "data/AS/GTEx/post/SE.MATS.JC.txt",
}
SE_SAMPLE_FILES = {
    "TCGA_Tumor": "data/AS/TCGA/all_samples.txt",
    "CPTAC_Tumor": "data/AS/CPTAC/all_samples.txt",
    "GTEx_Normal": "data/AS/GTEx/all_samples.txt",
}
EXPRESSION_FILES = {
    "TCGA_Tumor": "data/Merged/Group_TPM_TCGA-PAAD_Tumor.csv",
    "CPTAC_Tumor": "data/Merged/Group_TPM_CPTAC-3_Tumor.csv",
    "GTEx_Normal": "data/Merged/Group_TPM_GTEx_Normal.csv",
}

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
    "red": "#B2182B",
    "navy": "#2166AC",
}


@dataclass(frozen=True)
class SiteAvailabilityOutputs:
    site_availability: Path
    site_bias: Path
    site_manifest: Path


@dataclass(frozen=True)
class SampleModelOutputs:
    model_summary: Path
    event_effects: Path
    model_manifest: Path


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


def _read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt", ".bed"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False, **kwargs)


def _write_tsv(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def _relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _safe_numeric(series: pd.Series | Iterable[object]) -> pd.Series:
    return pd.to_numeric(pd.Series(series), errors="coerce")


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


def _zscore(values: pd.Series | np.ndarray) -> pd.Series:
    series = _safe_numeric(values)
    mean = series.mean(skipna=True)
    sd = series.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - mean) / sd


def _event_key(value: object) -> str:
    return rbp.standardize_event_id(value)


def _transcript_core(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    transcript = text.split("|")[0] if text else ""
    return re.sub(r"\.\d+$", "", transcript)


def _transcript_id(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return text.split("|")[0] if text else ""


def _standardize_sample_id(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().strip('"')
    if not text:
        return ""
    text = Path(text).name
    for suffix in [
        ".Aligned.sortedByCoord.out.patched.md.bam",
        ".Aligned.sortedByCoord.out.bam",
        ".patched.md.bam",
        ".bam",
    ]:
        if text.endswith(suffix):
            text = text[: -len(suffix)]
    text = text.replace(", ", "_").replace(",", "_").replace(" ", "")
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


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


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    color = color.lstrip("#")
    return tuple(int(color[i : i + 2], 16) for i in (0, 2, 4))


def _blend(c1: str, c2: str, t: float) -> tuple[int, int, int]:
    a = np.array(_hex_to_rgb(c1), dtype=float)
    b = np.array(_hex_to_rgb(c2), dtype=float)
    return tuple(np.clip(a * (1 - t) + b * t, 0, 255).astype(int))


def _draw_text(draw, xy, text: str, font, fill=None, anchor=None):
    draw.text(xy, str(text), font=font, fill=fill or _hex_to_rgb(PALETTE["dark"]), anchor=anchor)


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


def load_mirna_site_burden(project_root: Path) -> pd.DataFrame:
    path = project_root / "data" / "apa_mirna_rhythmic_connection" / "apa_mirna_gene_triads.tsv"
    columns = [
        "apa_event_key",
        "transcript_id",
        "transcript_core",
        "host_gene_symbol",
        "n_mirna_ids",
        "n_seed_families",
        "n_targetscan_transcripts",
        "total_conserved_sites",
        "total_nonconserved_sites",
        "total_mirna_sites",
    ]
    if not path.exists():
        return pd.DataFrame(columns=columns)
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    wanted = [
        "mirna_id_target_db",
        "mirna_id",
        "seed_family",
        "n_targetscan_transcripts",
        "targetscan_total_conserved_sites",
        "targetscan_total_nonconserved_sites",
        "transcript_id",
        "apa_event_id",
        "host_gene_symbol",
        "target_gene_symbol",
        "apa_gene_symbol",
    ]
    usecols = [col for col in wanted if col in header]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    if df.empty or "apa_event_id" not in df.columns:
        return pd.DataFrame(columns=columns)
    df["apa_event_key"] = df["apa_event_id"].map(_event_key)
    df["transcript_id"] = df.get("transcript_id", df["apa_event_id"]).fillna(df["apa_event_id"]).map(_transcript_id)
    df["transcript_core"] = df["transcript_id"].map(_transcript_core)
    if "host_gene_symbol" not in df.columns:
        df["host_gene_symbol"] = ""
    fallback_symbol = df.get("apa_gene_symbol", df.get("target_gene_symbol", pd.Series("", index=df.index)))
    df["host_gene_symbol"] = df["host_gene_symbol"].fillna(fallback_symbol).map(rbp.standardize_symbol)
    df["mirna_key"] = df.get("mirna_id_target_db", df.get("mirna_id", pd.Series("", index=df.index))).fillna("").astype(str)
    df["seed_family"] = df.get("seed_family", pd.Series("", index=df.index)).fillna("").astype(str)
    for col in [
        "n_targetscan_transcripts",
        "targetscan_total_conserved_sites",
        "targetscan_total_nonconserved_sites",
    ]:
        if col not in df.columns:
            df[col] = 0
        df[col] = _safe_numeric(df[col]).fillna(0)
    df = df[df["apa_event_key"].ne("")].copy()
    if df.empty:
        return pd.DataFrame(columns=columns)

    dedup = df.drop_duplicates(["apa_event_key", "mirna_key", "seed_family"]).copy()
    out = (
        dedup.groupby(["apa_event_key", "transcript_id", "transcript_core", "host_gene_symbol"], dropna=False)
        .agg(
            n_mirna_ids=("mirna_key", lambda x: int(pd.Series(x).replace("", np.nan).nunique())),
            n_seed_families=("seed_family", lambda x: int(pd.Series(x).replace("", np.nan).nunique())),
            n_targetscan_transcripts=("n_targetscan_transcripts", "max"),
            total_conserved_sites=("targetscan_total_conserved_sites", "sum"),
            total_nonconserved_sites=("targetscan_total_nonconserved_sites", "sum"),
        )
        .reset_index()
    )
    out["total_mirna_sites"] = out["total_conserved_sites"] + out["total_nonconserved_sites"]
    return out[columns].copy()


def load_clip_utr_site_burden(project_root: Path, analysis: str = ANALYSIS) -> pd.DataFrame:
    path = project_root / "data" / analysis / "clip_supported_rbp_target_links_consensus.tsv"
    columns = [
        "host_gene_symbol",
        "n_3utr_clip_rbps",
        "n_3utr_clip_links",
        "n_3utr_clip_accessions",
        "n_3utr_clip_methods",
        "max_3utr_clip_score",
        "median_3utr_clip_score",
    ]
    if not path.exists():
        return pd.DataFrame(columns=columns)
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    wanted = [
        "target_gene_symbol",
        "region",
        "rbp_symbol",
        "n_accessions",
        "n_methods",
        "max_score",
        "median_score",
    ]
    usecols = [col for col in wanted if col in header]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    if df.empty:
        return pd.DataFrame(columns=columns)
    df["region"] = df.get("region", pd.Series("", index=df.index)).astype(str)
    df = df[df["region"].str.contains("3UTR", case=False, na=False)].copy()
    if df.empty:
        return pd.DataFrame(columns=columns)
    df["host_gene_symbol"] = df.get("target_gene_symbol", pd.Series("", index=df.index)).map(rbp.standardize_symbol)
    df = df[df["host_gene_symbol"].ne("")]
    for col in ["n_accessions", "n_methods", "max_score", "median_score"]:
        if col not in df.columns:
            df[col] = 0
        df[col] = _safe_numeric(df[col]).fillna(0)
    out = (
        df.groupby("host_gene_symbol", dropna=False)
        .agg(
            n_3utr_clip_rbps=("rbp_symbol", lambda x: int(pd.Series(x).map(rbp.standardize_symbol).replace("", np.nan).nunique())),
            n_3utr_clip_links=("host_gene_symbol", "size"),
            n_3utr_clip_accessions=("n_accessions", "sum"),
            n_3utr_clip_methods=("n_methods", "sum"),
            max_3utr_clip_score=("max_score", "max"),
            median_3utr_clip_score=("median_score", "median"),
        )
        .reset_index()
    )
    return out[columns].copy()


def build_site_availability(project_root: Path, analysis: str = ANALYSIS) -> pd.DataFrame:
    output_dir = project_root / "data" / analysis
    apa_phase = _read_table(output_dir / "apa_expression_phase_lag_inputs.tsv")
    if apa_phase.empty:
        apa_phase = rbp.load_apa_phase_tables(project_root)
    if apa_phase.empty:
        return pd.DataFrame()

    utr = rbp_deep.load_utr_annotation(project_root)
    site = rbp_deep.annotate_apa_utr_rank(apa_phase, utr)
    site["apa_event_key"] = site["apa_event_id"].map(_event_key)
    site["transcript_id"] = site.get("transcript_id", site["apa_event_id"]).fillna(site["apa_event_id"]).map(_transcript_id)
    site["transcript_core"] = site["transcript_id"].map(_transcript_core)
    site["host_gene_symbol"] = site.get("host_gene_symbol", pd.Series("", index=site.index)).map(rbp.standardize_symbol)

    mirna = load_mirna_site_burden(project_root)
    if not mirna.empty:
        merge_cols = [
            "apa_event_key",
            "n_mirna_ids",
            "n_seed_families",
            "n_targetscan_transcripts",
            "total_conserved_sites",
            "total_nonconserved_sites",
            "total_mirna_sites",
        ]
        site = site.merge(
            mirna[merge_cols].drop_duplicates("apa_event_key"),
            on="apa_event_key",
            how="left",
        )
    for col in [
        "n_mirna_ids",
        "n_seed_families",
        "n_targetscan_transcripts",
        "total_conserved_sites",
        "total_nonconserved_sites",
        "total_mirna_sites",
    ]:
        if col not in site.columns:
            site[col] = 0
        site[col] = _safe_numeric(site[col]).fillna(0)

    clip_burden = load_clip_utr_site_burden(project_root, analysis)
    if not clip_burden.empty:
        site = site.merge(clip_burden, on="host_gene_symbol", how="left")
    for col in [
        "n_3utr_clip_rbps",
        "n_3utr_clip_links",
        "n_3utr_clip_accessions",
        "n_3utr_clip_methods",
        "max_3utr_clip_score",
        "median_3utr_clip_score",
    ]:
        if col not in site.columns:
            site[col] = 0
        site[col] = _safe_numeric(site[col]).fillna(0)

    components = pd.DataFrame(
        {
            "mirna_sites": np.log1p(site["total_mirna_sites"]),
            "seed_families": np.log1p(site["n_seed_families"]),
            "clip_rbps": np.log1p(site["n_3utr_clip_rbps"]),
            "clip_links": np.log1p(site["n_3utr_clip_links"]),
            "utr_length": np.log1p(_safe_numeric(site.get("utr_length_bp", 0)).fillna(0)),
        },
        index=site.index,
    )
    z_components = components.apply(_zscore)
    site["regulatory_site_availability_score"] = z_components.mean(axis=1)
    site["has_mirna_site_proxy"] = site["total_mirna_sites"].gt(0)
    site["has_3utr_clip_site_proxy"] = site["n_3utr_clip_rbps"].gt(0)
    site["site_availability_proxy_note"] = (
        "TargetScan transcript/gene miRNA burden plus broad 3'UTR CLIP burden; "
        "not direct PAS-gain/loss interval counting."
    )
    return site


def compute_site_availability_bias(site: pd.DataFrame, apa_lags: pd.DataFrame) -> pd.DataFrame:
    if site.empty or apa_lags.empty:
        return pd.DataFrame()
    bg = site.copy()
    bg["apa_event_key"] = bg["apa_event_key"].map(_event_key)
    bg = bg.drop_duplicates(["comparison", "apa_event_key"])
    lag = apa_lags.copy()
    lag["apa_event_key"] = lag["apa_event_id"].map(_event_key)
    lag = lag.drop_duplicates(["comparison", "rbp_symbol", "apa_event_key"])

    metrics = [
        "n_seed_families",
        "total_mirna_sites",
        "n_3utr_clip_rbps",
        "n_3utr_clip_links",
        "regulatory_site_availability_score",
        "utr_length_quantile",
    ]
    rows = []
    for (comparison, rbp_symbol), grp in lag.groupby(["comparison", "rbp_symbol"], dropna=False):
        background = bg[bg["comparison"].eq(comparison)].copy()
        if len(background) < 20:
            continue
        bound_keys = set(grp["apa_event_key"])
        bound = background[background["apa_event_key"].isin(bound_keys)].copy()
        if len(bound) < 8:
            continue
        row = {
            "comparison": comparison,
            "rbp_symbol": rbp_symbol,
            "n_bound_apa_events": int(len(bound)),
            "n_background_apa_events": int(len(background)),
        }
        for metric in metrics:
            b = _safe_numeric(bound.get(metric, pd.Series(dtype=float))).dropna()
            a = _safe_numeric(background.get(metric, pd.Series(dtype=float))).dropna()
            if len(b) and len(a):
                row[f"median_bound_{metric}"] = float(b.median())
                row[f"median_background_{metric}"] = float(a.median())
                row[f"delta_{metric}"] = float(b.median() - a.median())
                try:
                    row[f"p_value_{metric}"] = float(stats.mannwhitneyu(b, a, alternative="two-sided").pvalue)
                except ValueError:
                    row[f"p_value_{metric}"] = np.nan
            else:
                row[f"median_bound_{metric}"] = np.nan
                row[f"median_background_{metric}"] = np.nan
                row[f"delta_{metric}"] = np.nan
                row[f"p_value_{metric}"] = np.nan
        examples = (
            bound.assign(
                example_score=lambda x: _safe_numeric(x["total_mirna_sites"]).fillna(0)
                + _safe_numeric(x["n_3utr_clip_links"]).fillna(0)
                + _safe_numeric(x["utr_length_quantile"]).fillna(0)
            )
            .sort_values(["example_score", "fdr"], ascending=[False, True])
            .get("host_gene_symbol", pd.Series(dtype=str))
            .dropna()
            .astype(str)
            .head(8)
            .tolist()
        )
        row["top_site_available_examples"] = ",".join(dict.fromkeys(examples))
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["availability_bias_p_value"] = out["p_value_regulatory_site_availability_score"]
    out["availability_bias_fdr"] = _bh_fdr(out["availability_bias_p_value"])
    out["availability_bias_direction"] = "mixed_or_flat"
    out.loc[out["delta_regulatory_site_availability_score"].ge(0.15), "availability_bias_direction"] = "higher_site_availability"
    out.loc[out["delta_regulatory_site_availability_score"].le(-0.15), "availability_bias_direction"] = "lower_site_availability"
    return out.sort_values(
        ["availability_bias_fdr", "n_bound_apa_events", "delta_regulatory_site_availability_score"],
        ascending=[True, False, False],
    )


def _delta_color(value: float, scale: float) -> tuple[int, int, int]:
    if pd.isna(value):
        return _hex_to_rgb("#F3F3F3")
    t = float(np.clip(abs(value) / max(scale, 1e-9), 0, 1))
    if value >= 0:
        return _blend("#F7F7F7", PALETTE["red"], t)
    return _blend("#F7F7F7", PALETTE["navy"], t)


def plot_site_availability_bias(site_bias: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if site_bias.empty:
        return []
    from PIL import Image, ImageDraw

    metrics = [
        ("delta_n_seed_families", "seed fam"),
        ("delta_total_mirna_sites", "miRNA sites"),
        ("delta_n_3utr_clip_rbps", "3'UTR RBPs"),
        ("delta_n_3utr_clip_links", "3'UTR links"),
        ("delta_regulatory_site_availability_score", "site score"),
    ]
    plot = site_bias.copy()
    plot["rank_score"] = _safe_numeric(plot["delta_regulatory_site_availability_score"]).abs().fillna(0)
    plot = plot.sort_values(["availability_bias_fdr", "rank_score", "n_bound_apa_events"], ascending=[True, False, False]).head(22)
    row_h = 31
    cell_w = 108
    left = 335
    top = 108
    width = 1120
    height = top + row_h * len(plot) + 105
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    small = _load_font(13)
    font = _load_font(16)
    bold = _load_font(22)

    _draw_text(draw, (30, 28), "APA regulatory-site availability around RBP-bound isoforms", bold)
    _draw_text(
        draw,
        (30, 58),
        "Cells show bound minus background medians; red=higher availability in bound APA, blue=lower",
        small,
        PALETTE["gray"],
    )
    for j, (_, label) in enumerate(metrics):
        x = left + j * cell_w + cell_w // 2
        _draw_text(draw, (x, top - 28), label, small, anchor="ma")

    scales = {}
    for metric, _ in metrics:
        values = _safe_numeric(plot[metric]).abs().dropna()
        scales[metric] = max(float(values.quantile(0.9)) if len(values) else 1.0, 1e-6)

    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        label = f"{row.comparison.replace('_vs_', ' vs ')}  {row.rbp_symbol}  n={row.n_bound_apa_events}"
        _draw_text(draw, (30, y + 2), label, small)
        for j, (metric, _) in enumerate(metrics):
            value = getattr(row, metric)
            x0 = left + j * cell_w
            draw.rectangle(
                (x0, y - 11, x0 + cell_w - 5, y + 15),
                fill=_delta_color(value, scales[metric]),
                outline=(255, 255, 255),
            )
            if pd.notna(value):
                _draw_text(draw, (x0 + cell_w // 2, y + 2), f"{float(value):.2g}", small, anchor="ma")
        fdr = row.availability_bias_fdr
        fdr_label = f"FDR={fdr:.2g}" if pd.notna(fdr) else "FDR=NA"
        _draw_text(draw, (left + len(metrics) * cell_w + 20, y + 2), fdr_label, small, PALETTE["gray"])
    _draw_text(draw, (left + len(metrics) * cell_w + 20, top - 28), "site score", small)
    _draw_text(
        draw,
        (left + 230, height - 35),
        "Availability proxy: TargetScan miRNA burden + broad 3'UTR CLIP burden + annotated 3'UTR rank",
        font,
        PALETTE["gray"],
        anchor="ma",
    )
    return _save_png_pdf(img, figures_dir, "rbp_apa_regulatory_site_availability")


def run_site_availability(
    project_root: str | os.PathLike[str] | None = None,
    analysis: str = ANALYSIS,
    make_plots: bool = True,
) -> dict[str, pd.DataFrame | list[Path]]:
    project_root = resolve_project_root(project_root)
    output_dir = project_root / "data" / analysis
    figures_dir = project_root / "data" / "figures" / analysis
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    site = build_site_availability(project_root, analysis)
    apa_lags = _read_table(output_dir / "rbp_apa_phase_lags.tsv")
    site_bias = compute_site_availability_bias(site, apa_lags)

    site_path = _write_tsv(site, output_dir / "rbp_apa_regulatory_site_availability.tsv")
    site_bias_path = _write_tsv(site_bias, output_dir / "rbp_apa_regulatory_site_bias.tsv")
    figures = plot_site_availability_bias(site_bias, figures_dir) if make_plots else []

    manifest = pd.DataFrame(
        [
            {"key": "site_availability_rows", "value": len(site), "detail": _relative_path(site_path, project_root)},
            {"key": "site_bias_rows", "value": len(site_bias), "detail": _relative_path(site_bias_path, project_root)},
            {
                "key": "availability_proxy_caveat",
                "value": "proxy",
                "detail": "Uses transcript-level TargetScan/CLIP burden and UTR rank; direct PAS-specific gained/lost motif intervals remain future work.",
            },
            {"key": "figures", "value": len(figures), "detail": ";".join(_relative_path(p, project_root) for p in figures)},
        ]
    )
    manifest_path = _write_tsv(manifest, output_dir / "rbp_apa_regulatory_site_manifest.tsv")
    return {
        "site_availability": site,
        "site_bias": site_bias,
        "site_manifest": manifest,
        "figures": figures,
        "paths": SiteAvailabilityOutputs(site_path, site_bias_path, manifest_path),
    }


def load_expression_matrix(project_root: Path, cohort_key: str) -> pd.DataFrame:
    rel = EXPRESSION_FILES.get(cohort_key)
    if rel is None:
        return pd.DataFrame()
    path = project_root / rel
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.empty or "Gene" not in df.columns:
        return pd.DataFrame()
    df["gene_symbol"] = df["Gene"].map(rbp.standardize_symbol)
    df = df[df["gene_symbol"].ne("")].drop(columns=["Gene"])
    sample_cols = [col for col in df.columns if col != "gene_symbol"]
    rename = {col: _standardize_sample_id(col) for col in sample_cols}
    df = df.rename(columns=rename)
    value_cols = [col for col in df.columns if col != "gene_symbol" and col]
    values = df[value_cols].apply(pd.to_numeric, errors="coerce")
    values.insert(0, "gene_symbol", df["gene_symbol"].values)
    values = values.groupby("gene_symbol", dropna=False).mean()
    return np.log2(values + 1.0)


def load_apa_usage_matrix(project_root: Path, cohort_key: str, event_ids: set[str] | None = None) -> pd.DataFrame:
    rel = APA_USAGE_FILES.get(cohort_key)
    if rel is None:
        return pd.DataFrame()
    path = project_root / rel
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, low_memory=False)
    if df.empty:
        return pd.DataFrame()
    first = df.columns[0]
    df = df.rename(columns={first: "event_id"})
    df["event_id"] = df["event_id"].map(_event_key)
    df = df[df["event_id"].ne("")]
    if event_ids:
        df = df[df["event_id"].isin(event_ids)]
    if df.empty:
        return pd.DataFrame()
    sample_cols = [col for col in df.columns if col != "event_id"]
    rename = {col: _standardize_sample_id(col) for col in sample_cols}
    df = df.rename(columns=rename)
    value_cols = [col for col in df.columns if col != "event_id" and col]
    values = df[value_cols].apply(pd.to_numeric, errors="coerce")
    values.insert(0, "event_id", df["event_id"].values)
    return values.drop_duplicates("event_id").set_index("event_id")


def _read_rmats_samples(project_root: Path, cohort_key: str) -> list[str]:
    rel = SE_SAMPLE_FILES.get(cohort_key)
    if rel is None:
        return []
    path = project_root / rel
    if not path.exists():
        return []
    samples = []
    for item in path.read_text().strip().split(","):
        sample = _standardize_sample_id(item)
        if sample:
            samples.append(sample)
    return samples


def load_se_usage_matrix(
    project_root: Path,
    cohort_key: str,
    event_ids: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rel = SE_USAGE_FILES.get(cohort_key)
    if rel is None:
        return pd.DataFrame(), pd.DataFrame()
    path = project_root / rel
    if not path.exists():
        return pd.DataFrame(), pd.DataFrame()
    header = pd.read_csv(path, nrows=0).columns
    id_col = "ID" if "ID" in header else header[0]
    gene_col = "geneSymbol" if "geneSymbol" in header else ("GeneID" if "GeneID" in header else id_col)
    inc_col = "IncLevel1" if "IncLevel1" in header else ""
    if not inc_col:
        return pd.DataFrame(), pd.DataFrame()
    df = pd.read_csv(path, usecols=[id_col, gene_col, inc_col], low_memory=False)
    df["event_id"] = df[id_col].map(_event_key)
    df["gene_symbol"] = df[gene_col].map(rbp.standardize_symbol)
    if event_ids:
        df = df[df["event_id"].isin(event_ids)]
    df = df[df["event_id"].ne("")]
    samples = _read_rmats_samples(project_root, cohort_key)
    if df.empty or not samples:
        return pd.DataFrame(), pd.DataFrame()

    rows = []
    index = []
    genes = []
    for row in df.itertuples(index=False):
        event_id = getattr(row, "event_id")
        values = str(getattr(row, inc_col)).split(",")
        if len(values) != len(samples):
            continue
        numeric = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
        rows.append(numeric)
        index.append(event_id)
        genes.append(getattr(row, "gene_symbol"))
    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    usage = pd.DataFrame(rows, index=index, columns=samples)
    usage = usage[~usage.index.duplicated(keep="first")]
    meta = pd.DataFrame({"event_id": index, "gene_symbol": genes}).drop_duplicates("event_id").set_index("event_id")
    return usage, meta


def select_model_candidates(project_root: Path, analysis: str = ANALYSIS, max_models: int = 48) -> pd.DataFrame:
    path = project_root / "data" / analysis / "rbp_deep_candidate_priority.tsv"
    candidates = _read_table(path)
    if candidates.empty:
        rows = []
        for rbp_symbol in sorted(FOCUS_RBPS):
            for comparison, cfg in rbp.COMPARISONS.items():
                rows.append(
                    {
                        "comparison": comparison,
                        "cohort_key": cfg["tumor_key"],
                        "event_kind": "APA",
                        "rbp_symbol": rbp_symbol,
                        "deep_priority_score": 0,
                        "candidate_reason": "focus_rbp_fallback",
                    }
                )
        return pd.DataFrame(rows)
    candidates["rbp_symbol"] = candidates["rbp_symbol"].map(rbp.standardize_symbol)
    candidates["deep_priority_score"] = _safe_numeric(candidates.get("deep_priority_score", 0)).fillna(0)
    candidates["candidate_reason"] = "deep_priority"
    keep = candidates[
        candidates["suggested_panel_use"].astype(str).ne("backup")
        | candidates["rbp_symbol"].isin(FOCUS_RBPS)
        | candidates["deep_priority_score"].ge(candidates["deep_priority_score"].quantile(0.70))
    ].copy()
    keep = keep.sort_values(["deep_priority_score", "n_events"], ascending=[False, False])
    if len(keep) < max_models:
        keep = pd.concat([keep, candidates.sort_values("deep_priority_score", ascending=False).head(max_models)], ignore_index=True)
    keep = keep.drop_duplicates(["comparison", "event_kind", "rbp_symbol"]).head(max_models)
    cohort_map = {comparison: cfg["tumor_key"] for comparison, cfg in rbp.COMPARISONS.items()}
    keep["cohort_key"] = keep["comparison"].map(cohort_map)
    return keep[["comparison", "cohort_key", "event_kind", "rbp_symbol", "deep_priority_score", "candidate_reason"]].copy()


def _aggregate_lag_events(
    lag_group: pd.DataFrame,
    event_kind: str,
    site_availability: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if lag_group.empty:
        return pd.DataFrame()
    lag = lag_group.copy()
    event_col = "apa_event_id" if event_kind == "APA" else "event_id"
    gene_col = "host_gene_symbol" if event_kind == "APA" else "gene_symbol"
    lag["event_id"] = lag[event_col].map(_event_key)
    lag["gene_symbol_model"] = lag[gene_col].map(rbp.standardize_symbol)
    for col in ["n_merged_peaks", "n_accessions", "n_methods", "n_samples", "max_score", "median_score"]:
        if col not in lag.columns:
            lag[col] = 0
        lag[col] = _safe_numeric(lag[col]).fillna(0)
    agg = (
        lag.groupby("event_id", dropna=False)
        .agg(
            gene_symbol=("gene_symbol_model", "first"),
            n_clip_rows=("event_id", "size"),
            n_merged_peaks=("n_merged_peaks", "max"),
            n_accessions=("n_accessions", "max"),
            n_methods=("n_methods", "max"),
            n_samples=("n_samples", "max"),
            max_score=("max_score", "max"),
            median_score=("median_score", "median"),
        )
        .reset_index()
    )
    agg = agg[agg["event_id"].ne("") & agg["gene_symbol"].ne("")]
    if event_kind == "APA" and site_availability is not None and not site_availability.empty:
        site_cols = [
            "comparison",
            "apa_event_key",
            "n_seed_families",
            "total_mirna_sites",
            "n_3utr_clip_rbps",
            "n_3utr_clip_links",
            "utr_length_quantile",
            "regulatory_site_availability_score",
        ]
        available_cols = [col for col in site_cols if col in site_availability.columns]
        site = site_availability[available_cols].drop_duplicates(["comparison", "apa_event_key"])
        comparison = lag_group["comparison"].iloc[0] if "comparison" in lag_group else ""
        site = site[site["comparison"].eq(comparison)].rename(columns={"apa_event_key": "event_id"})
        agg = agg.merge(site.drop(columns=["comparison"]), on="event_id", how="left")
    for col in ["n_seed_families", "total_mirna_sites", "n_3utr_clip_rbps", "n_3utr_clip_links", "utr_length_quantile"]:
        if col not in agg.columns:
            agg[col] = 0
        agg[col] = _safe_numeric(agg[col]).fillna(0)
    return agg.sort_values(["event_id", "gene_symbol"], ascending=[True, True])


def _fit_ols(y: pd.Series, x: pd.DataFrame) -> dict[str, float]:
    data = pd.concat([y.rename("y"), x], axis=1).dropna()
    if len(data) <= x.shape[1] + 2:
        return {"n": len(data), "r2": np.nan, "adj_r2": np.nan, "coef_last": np.nan, "p_last": np.nan}
    yv = data["y"].to_numpy(dtype=float)
    xv = data.drop(columns=["y"]).to_numpy(dtype=float)
    xv = np.column_stack([np.ones(len(xv)), xv])
    coef = np.linalg.pinv(xv) @ yv
    fitted = xv @ coef
    resid = yv - fitted
    sse = float(np.sum(resid**2))
    sst = float(np.sum((yv - np.mean(yv)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    p = xv.shape[1]
    dof = len(yv) - p
    adj_r2 = 1.0 - (1.0 - r2) * (len(yv) - 1) / dof if pd.notna(r2) and dof > 0 else np.nan
    p_last = np.nan
    if dof > 0:
        sigma2 = sse / dof
        cov = sigma2 * np.linalg.pinv(xv.T @ xv)
        se = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
        if se[-1] > 0:
            t_stat = coef[-1] / se[-1]
            p_last = float(2 * stats.t.sf(abs(t_stat), dof))
    return {
        "n": int(len(data)),
        "r2": float(r2),
        "adj_r2": float(adj_r2) if pd.notna(adj_r2) else np.nan,
        "coef_last": float(coef[-1]),
        "p_last": p_last,
    }


def _fit_model_table(model_df: pd.DataFrame) -> dict[str, float]:
    df = model_df.copy()
    required = ["event_id", "sample_id", "usage", "host_expr", "rbp_expr"]
    df = df.dropna(subset=required)
    if len(df) < 100 or df["event_id"].nunique() < 3 or df["sample_id"].nunique() < 10:
        return {
            "n_observations": int(len(df)),
            "n_events_modeled": int(df["event_id"].nunique()),
            "n_samples_overlap": int(df["sample_id"].nunique()),
            "model_type": "within_event_centered_ols",
            "base_r2": np.nan,
            "full_r2": np.nan,
            "delta_r2": np.nan,
            "beta_rbp_expr_z": np.nan,
            "p_value": np.nan,
        }

    event_groups = df.groupby("event_id", dropna=False)
    df["usage_within_event"] = df["usage"] - event_groups["usage"].transform("mean")
    df["host_expr_within_event_z"] = _zscore(df["host_expr"] - event_groups["host_expr"].transform("mean"))
    df["rbp_expr_within_event_z"] = _zscore(df["rbp_expr"] - event_groups["rbp_expr"].transform("mean"))
    df = df.dropna(subset=["usage_within_event", "host_expr_within_event_z", "rbp_expr_within_event_z"])
    if len(df) < 100 or df["event_id"].nunique() < 3 or df["sample_id"].nunique() < 10:
        return {
            "n_observations": int(len(df)),
            "n_events_modeled": int(df["event_id"].nunique()),
            "n_samples_overlap": int(df["sample_id"].nunique()),
            "model_type": "within_event_centered_ols",
            "base_r2": np.nan,
            "full_r2": np.nan,
            "delta_r2": np.nan,
            "beta_rbp_expr_z": np.nan,
            "p_value": np.nan,
        }

    base = _fit_ols(df["usage_within_event"], df[["host_expr_within_event_z"]])
    full = _fit_ols(df["usage_within_event"], df[["host_expr_within_event_z", "rbp_expr_within_event_z"]])
    return {
        "n_observations": int(len(df)),
        "n_events_modeled": int(df["event_id"].nunique()),
        "n_samples_overlap": int(df["sample_id"].nunique()),
        "model_type": "within_event_centered_ols",
        "base_r2": base["r2"],
        "full_r2": full["r2"],
        "delta_r2": float(full["r2"] - base["r2"]) if pd.notna(base["r2"]) and pd.notna(full["r2"]) else np.nan,
        "beta_rbp_expr_z": full["coef_last"],
        "p_value": full["p_last"],
    }


def _event_level_effects(model_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for event_id, grp in model_df.groupby("event_id", dropna=False):
        df = grp.copy()
        if len(df) < 20:
            continue
        df["host_expr_z"] = _zscore(df["host_expr"])
        df["rbp_expr_z"] = _zscore(df["rbp_expr"])
        df = df.dropna(subset=["usage", "host_expr_z", "rbp_expr_z"])
        if len(df) < 20:
            continue
        base = _fit_ols(df["usage"], df[["host_expr_z"]])
        full = _fit_ols(df["usage"], df[["host_expr_z", "rbp_expr_z"]])
        rows.append(
            {
                "event_id": event_id,
                "host_gene_symbol": df["host_gene_symbol"].iloc[0],
                "n_samples": int(len(df)),
                "event_mean_usage": float(df["usage"].mean()),
                "event_base_r2": base["r2"],
                "event_full_r2": full["r2"],
                "event_delta_r2": float(full["r2"] - base["r2"]) if pd.notna(base["r2"]) and pd.notna(full["r2"]) else np.nan,
                "event_beta_rbp_expr_z": full["coef_last"],
                "event_p_value": full["p_last"],
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["event_fdr"] = _bh_fdr(out["event_p_value"])
    return out


def _build_model_dataframe(
    usage: pd.DataFrame,
    expression: pd.DataFrame,
    event_info: pd.DataFrame,
    rbp_symbol: str,
    event_kind: str,
    max_events: int | None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    if usage.empty or expression.empty or event_info.empty or rbp_symbol not in expression.index:
        return pd.DataFrame(), {
            "candidate_events": int(len(event_info)),
            "events_with_usage": 0,
            "events_with_host_expression": 0,
            "overlap_samples": 0,
        }
    samples = [sample for sample in usage.columns if sample in expression.columns]
    if not samples:
        return pd.DataFrame(), {
            "candidate_events": int(len(event_info)),
            "events_with_usage": 0,
            "events_with_host_expression": 0,
            "overlap_samples": 0,
        }
    event_info = event_info[event_info["event_id"].isin(usage.index)].copy()
    if max_events is not None and max_events > 0 and len(event_info) > max_events:
        # Deterministic fallback for constrained reruns. Main analyses should use
        # all available events so event inclusion is not driven by a composite score.
        event_info = event_info.sort_values(["event_id", "gene_symbol"]).head(max_events).copy()
    event_info = event_info[event_info["gene_symbol"].isin(expression.index)].copy()
    if event_info.empty:
        return pd.DataFrame(), {
            "candidate_events": int(len(event_info)),
            "events_with_usage": 0,
            "events_with_host_expression": 0,
            "overlap_samples": int(len(samples)),
        }
    rbp_expr = expression.loc[rbp_symbol, samples].astype(float)
    rows = []
    for event in event_info.itertuples(index=False):
        event_id = event.event_id
        host = event.gene_symbol
        usage_values = usage.loc[event_id, samples].astype(float)
        host_expr = expression.loc[host, samples].astype(float)
        tmp = pd.DataFrame(
            {
                "event_id": event_id,
                "event_kind": event_kind,
                "host_gene_symbol": host,
                "sample_id": samples,
                "usage": usage_values.to_numpy(dtype=float),
                "host_expr": host_expr.to_numpy(dtype=float),
                "rbp_expr": rbp_expr.to_numpy(dtype=float),
            }
        )
        tmp["event_mean_usage"] = float(tmp["usage"].mean(skipna=True))
        rows.append(tmp)
    if not rows:
        return pd.DataFrame(), {
            "candidate_events": int(len(event_info)),
            "events_with_usage": 0,
            "events_with_host_expression": 0,
            "overlap_samples": int(len(samples)),
        }
    model_df = pd.concat(rows, ignore_index=True, sort=False)
    model_df = model_df.dropna(subset=["usage", "host_expr", "rbp_expr"])
    return model_df, {
        "candidate_events": int(len(event_info)),
        "events_with_usage": int(event_info["event_id"].nunique()),
        "events_with_host_expression": int(event_info["gene_symbol"].nunique()),
        "overlap_samples": int(len(samples)),
    }


def plot_sample_model_delta_r2(model_summary: pd.DataFrame, figures_dir: Path) -> list[Path]:
    if model_summary.empty:
        return []
    from PIL import Image, ImageDraw

    plot = model_summary.copy()
    plot["delta_r2_rank"] = _safe_numeric(plot["delta_r2"]).fillna(-np.inf)
    plot = plot.sort_values(["delta_r2_rank", "n_events_modeled"], ascending=[False, False]).head(24)
    if plot.empty:
        return []
    row_h = 32
    left = 405
    right = 1040
    top = 100
    width = 1220
    height = top + row_h * len(plot) + 105
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    small = _load_font(13)
    font = _load_font(16)
    bold = _load_font(22)
    _draw_text(draw, (30, 28), "Within-event PDUI/PSI models: added explanatory power of RBP expression", bold)
    _draw_text(
        draw,
        (30, 58),
        "Delta R2 compares full model to within-event host-expression baseline; p-values are exploratory",
        small,
        PALETTE["gray"],
    )
    max_x = max(float(_safe_numeric(plot["delta_r2"]).max()), 0.01)
    max_x = min(max_x * 1.15, max(max_x, 0.05))
    draw.line((left, top - 10, left, top + row_h * len(plot) + 12), fill=_hex_to_rgb(PALETTE["gray"]), width=1)
    for tick in np.linspace(0, max_x, 5):
        x = int(left + tick / max_x * (right - left))
        draw.line((x, top + row_h * len(plot) + 12, x, top + row_h * len(plot) + 18), fill=_hex_to_rgb(PALETTE["dark"]), width=1)
        _draw_text(draw, (x, top + row_h * len(plot) + 25), f"{tick:.3f}", small, anchor="ma")
    for i, row in enumerate(plot.itertuples(index=False)):
        y = top + i * row_h
        label = f"{row.event_kind}  {row.rbp_symbol}  {row.comparison.replace('_vs_', ' vs ')}  events={row.n_events_modeled}"
        _draw_text(draw, (30, y - 4), label, small)
        delta = 0.0 if pd.isna(row.delta_r2) else max(float(row.delta_r2), 0.0)
        x = int(left + min(delta, max_x) / max_x * (right - left))
        color = PALETTE["blue"] if row.event_kind == "APA" else PALETTE["orange"]
        draw.rectangle((left, y - 9, x, y + 11), fill=_hex_to_rgb(color))
        direction = "+" if pd.notna(row.beta_rbp_expr_z) and row.beta_rbp_expr_z >= 0 else "-"
        fdr = row.model_fdr
        fdr_label = f"FDR={fdr:.2g}" if pd.notna(fdr) else "FDR=NA"
        _draw_text(draw, (right + 18, y - 4), f"{direction} {fdr_label}", small, PALETTE["gray"])
    _draw_text(draw, ((left + right) // 2, height - 35), "Incremental R2 from sample-level RBP expression", font, anchor="ma")
    return _save_png_pdf(img, figures_dir, "rbp_sample_level_usage_model_delta_r2")


def run_sample_level_models(
    project_root: str | os.PathLike[str] | None = None,
    analysis: str = ANALYSIS,
    make_plots: bool = True,
    max_models: int = 48,
    max_events_per_model: int | None = None,
) -> dict[str, pd.DataFrame | list[Path]]:
    project_root = resolve_project_root(project_root)
    output_dir = project_root / "data" / analysis
    figures_dir = project_root / "data" / "figures" / analysis
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    site_path = output_dir / "rbp_apa_regulatory_site_availability.tsv"
    site_availability = _read_table(site_path) if site_path.exists() else build_site_availability(project_root, analysis)
    apa_lags = _read_table(output_dir / "rbp_apa_phase_lags.tsv")
    se_lags = _read_table(output_dir / "rbp_se_phase_lags.tsv")
    candidates = select_model_candidates(project_root, analysis, max_models=max_models)

    expression_cache: dict[str, pd.DataFrame] = {}
    usage_cache: dict[tuple[str, str], tuple[pd.DataFrame, pd.DataFrame]] = {}
    summary_rows = []
    event_effect_frames = []
    manifest_rows = []

    for candidate in candidates.itertuples(index=False):
        comparison = candidate.comparison
        cohort_key = candidate.cohort_key
        event_kind = candidate.event_kind
        rbp_symbol = rbp.standardize_symbol(candidate.rbp_symbol)
        lag_source = apa_lags if event_kind == "APA" else se_lags
        if lag_source.empty:
            continue
        lag_group = lag_source[
            lag_source["comparison"].eq(comparison)
            & lag_source["cohort_key"].eq(cohort_key)
            & lag_source["rbp_symbol"].map(rbp.standardize_symbol).eq(rbp_symbol)
        ].copy()
        if lag_group.empty:
            continue
        event_info = _aggregate_lag_events(lag_group, event_kind, site_availability)
        if event_info.empty:
            continue
        event_ids = set(event_info["event_id"].astype(str))

        if cohort_key not in expression_cache:
            expression_cache[cohort_key] = load_expression_matrix(project_root, cohort_key)
        expression = expression_cache[cohort_key]
        usage_key = (cohort_key, event_kind)
        if usage_key not in usage_cache:
            if event_kind == "APA":
                usage_cache[usage_key] = (load_apa_usage_matrix(project_root, cohort_key, event_ids=None), pd.DataFrame())
            else:
                usage_cache[usage_key] = load_se_usage_matrix(project_root, cohort_key, event_ids=None)
        usage, se_meta = usage_cache[usage_key]
        if event_kind == "SE" and not se_meta.empty:
            event_info = event_info.merge(
                se_meta[["gene_symbol"]].rename(columns={"gene_symbol": "se_usage_gene_symbol"}),
                left_on="event_id",
                right_index=True,
                how="left",
            )
            event_info["gene_symbol"] = event_info["gene_symbol"].where(
                event_info["gene_symbol"].ne(""),
                event_info["se_usage_gene_symbol"].fillna(""),
            )
            event_info = event_info.drop(columns=["se_usage_gene_symbol"])

        model_df, counts = _build_model_dataframe(
            usage=usage,
            expression=expression,
            event_info=event_info,
            rbp_symbol=rbp_symbol,
            event_kind=event_kind,
            max_events=max_events_per_model,
        )
        fit = _fit_model_table(model_df)
        summary_row = {
            "comparison": comparison,
            "cohort_key": cohort_key,
            "event_kind": event_kind,
            "usage_metric": "PDUI" if event_kind == "APA" else "PSI",
            "rbp_symbol": rbp_symbol,
            "candidate_reason": candidate.candidate_reason,
            "deep_priority_score": float(candidate.deep_priority_score),
            "n_candidate_events": counts["candidate_events"],
            "n_events_with_usage": counts["events_with_usage"],
            "n_events_with_host_expression": counts["events_with_host_expression"],
            "n_overlap_samples_available": counts["overlap_samples"],
            **fit,
        }
        summary_rows.append(summary_row)

        event_effects = _event_level_effects(model_df)
        if not event_effects.empty:
            event_effects.insert(0, "rbp_symbol", rbp_symbol)
            event_effects.insert(0, "usage_metric", "PDUI" if event_kind == "APA" else "PSI")
            event_effects.insert(0, "event_kind", event_kind)
            event_effects.insert(0, "cohort_key", cohort_key)
            event_effects.insert(0, "comparison", comparison)
            event_effect_frames.append(event_effects)

        manifest_rows.append(
            {
                "comparison": comparison,
                "cohort_key": cohort_key,
                "event_kind": event_kind,
                "rbp_symbol": rbp_symbol,
                "candidate_events": counts["candidate_events"],
                "events_modeled": fit["n_events_modeled"],
                "samples_modeled": fit["n_samples_overlap"],
                "observations_modeled": fit["n_observations"],
            }
        )

    model_summary = pd.DataFrame(summary_rows)
    if not model_summary.empty:
        model_summary["model_fdr"] = _bh_fdr(model_summary["p_value"])
        model_summary["rbp_effect_direction"] = "flat_or_unstable"
        model_summary.loc[model_summary["beta_rbp_expr_z"].ge(0), "rbp_effect_direction"] = "higher_rbp_higher_usage"
        model_summary.loc[model_summary["beta_rbp_expr_z"].lt(0), "rbp_effect_direction"] = "higher_rbp_lower_usage"
        model_summary = model_summary.sort_values(["delta_r2", "n_events_modeled"], ascending=[False, False])
    event_effects = pd.concat(event_effect_frames, ignore_index=True, sort=False) if event_effect_frames else pd.DataFrame()
    if not event_effects.empty:
        event_effects["global_event_fdr"] = _bh_fdr(event_effects["event_p_value"])
        event_effects = event_effects.sort_values(["event_delta_r2", "n_samples"], ascending=[False, False])

    summary_path = _write_tsv(model_summary, output_dir / "rbp_sample_level_usage_model_summary.tsv")
    event_path = _write_tsv(event_effects, output_dir / "rbp_sample_level_usage_model_event_effects.tsv")
    figures = plot_sample_model_delta_r2(model_summary, figures_dir) if make_plots else []

    manifest = pd.DataFrame(manifest_rows)
    caveat = pd.DataFrame(
        [
            {
                "comparison": "all",
                "cohort_key": "tumor_cohorts_only",
                "event_kind": "APA_SE",
                "rbp_symbol": "model_caveat",
                "candidate_events": len(candidates),
                "events_modeled": int(model_summary["n_events_modeled"].sum()) if not model_summary.empty else 0,
                "samples_modeled": int(model_summary["n_samples_overlap"].max()) if not model_summary.empty else 0,
                "observations_modeled": int(model_summary["n_observations"].sum()) if not model_summary.empty else 0,
                "detail": "Exploratory within-event centered OLS; event-sample observations are not independent. Event-level targetability/mean-usage effects are absorbed by within-event centering. Normal-cohort sample-level modeling is not run because current CLIP-supported phase-lag candidates are tumor-cohort keyed.",
            }
        ]
    )
    manifest = pd.concat([manifest, caveat], ignore_index=True, sort=False)
    manifest_path = _write_tsv(manifest, output_dir / "rbp_sample_level_usage_model_manifest.tsv")

    return {
        "model_summary": model_summary,
        "event_effects": event_effects,
        "model_manifest": manifest,
        "figures": figures,
        "paths": SampleModelOutputs(summary_path, event_path, manifest_path),
    }
