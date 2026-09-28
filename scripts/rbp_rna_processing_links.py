"""RBP-centered links among APA, skipped exons, and rhythmic expression.

This module supports the PDAC_Rhythmicity RBP extension notebook. It keeps the
binding-evidence logic explicit: primary RBP-target claims should use
CLIP-supported gene links or peak overlaps, while motif-only evidence remains
outside this core analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

try:
    from scipy.stats import fisher_exact
except Exception:  # pragma: no cover - scipy is expected in analysis envs.
    fisher_exact = None


ANALYSIS_SLUG = "rbp_rna_processing_links"

EXPRESSION_COHORTS = {
    "TCGA_Tumor": "TCGA-PAAD_Tumor",
    "CPTAC_Tumor": "CPTAC-3_Tumor",
    "GTEx_Normal": "GTEx_Normal",
}

COMPARISONS = {
    "TCGA_vs_GTEx": {
        "tumor_key": "TCGA_Tumor",
        "normal_key": "GTEx_Normal",
        "tumor_label": "TCGA-PAAD Tumor",
        "normal_label": "GTEx Normal",
    },
    "CPTAC_vs_GTEx": {
        "tumor_key": "CPTAC_Tumor",
        "normal_key": "GTEx_Normal",
        "tumor_label": "CPTAC-3 Tumor",
        "normal_label": "GTEx Normal",
    },
}

RBP_COLUMNS = [
    "gene_symbol",
    "ensembl_id",
    "rbp_family",
    "functional_class",
    "clip_available",
]

CLIP_LINK_COLUMNS = [
    "rbp_symbol",
    "target_gene_symbol",
    "target_ensg",
    "event_id",
    "region",
    "source",
    "evidence_type",
    "genome_build",
    "chrom",
    "peak_start",
    "peak_end",
    "method",
    "sample",
    "accession",
    "score",
    "clip_file",
]

CLIP_SUPPORT_COLUMNS = [
    "n_merged_peaks",
    "n_accessions",
    "n_methods",
    "n_samples",
    "max_score",
    "median_score",
    "peak_width_bp",
    "consensus_filter",
]

CLIP_ALL_LINK_COLUMNS = CLIP_LINK_COLUMNS + CLIP_SUPPORT_COLUMNS

CONSENSUS_SCORE_THRESHOLD = 20.0
CONSENSUS_MIN_ACCESSIONS = 2
CONSENSUS_MIN_PEAKS = 2
CONSENSUS_MERGE_DISTANCE_BP = 20


@dataclass(frozen=True)
class RbpRecord:
    gene_symbol: str
    rbp_family: str
    functional_class: str


def strip_version(value: object) -> str:
    """Return an Ensembl-like identifier without version suffix."""
    if pd.isna(value):
        return ""
    text = str(value).strip().strip('"')
    if not text:
        return ""
    return text.split(".")[0]


def standardize_symbol(value: object) -> str:
    """Normalize a gene or RBP symbol for joins."""
    if pd.isna(value):
        return ""
    return str(value).strip().strip('"').upper()


def standardize_event_id(value: object) -> str:
    """Normalize event IDs that may be read as ints, floats, or strings."""
    if pd.isna(value):
        return ""
    text = str(value).strip().strip('"')
    if text.lower() in {"", "nan", "none", "<na>"}:
        return ""
    if text.endswith(".0"):
        prefix = text[:-2]
        if prefix.isdigit():
            return prefix
    return text


def ensure_output_dirs(project_root: Path, analysis_slug: str = ANALYSIS_SLUG) -> dict[str, Path]:
    """Create and return analysis output, figure, and log directories."""
    project_root = Path(project_root)
    dirs = {
        "output": project_root / "data" / analysis_slug,
        "figures": project_root / "data" / "figures" / analysis_slug,
        "logs": project_root / "data" / "logs" / analysis_slug,
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    """Write a TSV with parent directory creation and return its path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def _records_to_catalog(records: Sequence[RbpRecord]) -> pd.DataFrame:
    rows = [
        {
            "gene_symbol": standardize_symbol(r.gene_symbol),
            "ensembl_id": "",
            "rbp_family": r.rbp_family,
            "functional_class": r.functional_class,
            "clip_available": False,
        }
        for r in records
    ]
    catalog = pd.DataFrame(rows).drop_duplicates("gene_symbol")
    return catalog.sort_values(["functional_class", "gene_symbol"]).reset_index(drop=True)


def default_rbp_catalog() -> pd.DataFrame:
    """Return a curated RBP catalog covering splicing, APA, 3'UTR, and miRNA layers."""
    records = [
        # Alternative splicing and spliceosome-associated RBPs.
        RbpRecord("SRSF1", "SR protein", "splicing"),
        RbpRecord("SRSF2", "SR protein", "splicing"),
        RbpRecord("SRSF3", "SR protein", "splicing"),
        RbpRecord("SRSF4", "SR protein", "splicing"),
        RbpRecord("SRSF5", "SR protein", "splicing"),
        RbpRecord("SRSF6", "SR protein", "splicing"),
        RbpRecord("SRSF7", "SR protein", "splicing"),
        RbpRecord("SRSF8", "SR protein", "splicing"),
        RbpRecord("SRSF9", "SR protein", "splicing"),
        RbpRecord("SRSF10", "SR protein", "splicing"),
        RbpRecord("SRSF11", "SR protein", "splicing"),
        RbpRecord("SRSF12", "SR protein", "splicing"),
        RbpRecord("HNRNPA0", "hnRNP", "splicing"),
        RbpRecord("HNRNPA1", "hnRNP", "splicing"),
        RbpRecord("HNRNPA2B1", "hnRNP", "splicing"),
        RbpRecord("HNRNPA3", "hnRNP", "splicing"),
        RbpRecord("HNRNPAB", "hnRNP", "splicing"),
        RbpRecord("HNRNPC", "hnRNP", "splicing"),
        RbpRecord("HNRNPD", "hnRNP", "splicing"),
        RbpRecord("HNRNPDL", "hnRNP", "splicing"),
        RbpRecord("HNRNPF", "hnRNP", "splicing"),
        RbpRecord("HNRNPH1", "hnRNP", "splicing"),
        RbpRecord("HNRNPH2", "hnRNP", "splicing"),
        RbpRecord("HNRNPH3", "hnRNP", "splicing"),
        RbpRecord("HNRNPK", "hnRNP", "splicing"),
        RbpRecord("HNRNPL", "hnRNP", "splicing"),
        RbpRecord("HNRNPLL", "hnRNP", "splicing"),
        RbpRecord("HNRNPM", "hnRNP", "splicing"),
        RbpRecord("HNRNPR", "hnRNP", "splicing"),
        RbpRecord("HNRNPU", "hnRNP", "splicing"),
        RbpRecord("RBFOX1", "RRM", "splicing"),
        RbpRecord("RBFOX2", "RRM", "splicing"),
        RbpRecord("RBFOX3", "RRM", "splicing"),
        RbpRecord("PTBP1", "RRM", "splicing"),
        RbpRecord("PTBP2", "RRM", "splicing"),
        RbpRecord("PTBP3", "RRM", "splicing"),
        RbpRecord("QKI", "STAR/KH", "splicing"),
        RbpRecord("NOVA1", "KH", "splicing"),
        RbpRecord("NOVA2", "KH", "splicing"),
        RbpRecord("MBNL1", "zinc finger", "splicing"),
        RbpRecord("MBNL2", "zinc finger", "splicing"),
        RbpRecord("MBNL3", "zinc finger", "splicing"),
        RbpRecord("CELF1", "CELF", "splicing"),
        RbpRecord("CELF2", "CELF", "splicing"),
        RbpRecord("CELF3", "CELF", "splicing"),
        RbpRecord("CELF4", "CELF", "splicing"),
        RbpRecord("CELF5", "CELF", "splicing"),
        RbpRecord("CELF6", "CELF", "splicing"),
        RbpRecord("ESRP1", "RRM", "splicing"),
        RbpRecord("ESRP2", "RRM", "splicing"),
        RbpRecord("TRA2A", "SR-related", "splicing"),
        RbpRecord("TRA2B", "SR-related", "splicing"),
        RbpRecord("U2AF1", "spliceosome", "splicing"),
        RbpRecord("U2AF2", "spliceosome", "splicing"),
        RbpRecord("SF3B1", "spliceosome", "splicing"),
        RbpRecord("SF3B2", "spliceosome", "splicing"),
        RbpRecord("SF3B3", "spliceosome", "splicing"),
        RbpRecord("SF3B4", "spliceosome", "splicing"),
        RbpRecord("PRPF3", "spliceosome", "splicing"),
        RbpRecord("PRPF4", "spliceosome", "splicing"),
        RbpRecord("PRPF6", "spliceosome", "splicing"),
        RbpRecord("PRPF8", "spliceosome", "splicing"),
        RbpRecord("PRPF18", "spliceosome", "splicing"),
        RbpRecord("PRPF19", "spliceosome", "splicing"),
        RbpRecord("PRPF31", "spliceosome", "splicing"),
        RbpRecord("PRPF38A", "spliceosome", "splicing"),
        RbpRecord("PRPF38B", "spliceosome", "splicing"),
        RbpRecord("SNRPA", "snRNP", "splicing"),
        RbpRecord("SNRPB", "snRNP", "splicing"),
        RbpRecord("SNRPC", "snRNP", "splicing"),
        RbpRecord("SNRPE", "snRNP", "splicing"),
        RbpRecord("SNRPF", "snRNP", "splicing"),
        RbpRecord("SNRPG", "snRNP", "splicing"),
        # 3'UTR stability, localization, translation, and miRNA-related RBPs.
        RbpRecord("ELAVL1", "ELAV/Hu", "3utr_stability_translation"),
        RbpRecord("ELAVL2", "ELAV/Hu", "3utr_stability_translation"),
        RbpRecord("ELAVL3", "ELAV/Hu", "3utr_stability_translation"),
        RbpRecord("ELAVL4", "ELAV/Hu", "3utr_stability_translation"),
        RbpRecord("IGF2BP1", "IGF2BP", "3utr_stability_translation"),
        RbpRecord("IGF2BP2", "IGF2BP", "3utr_stability_translation"),
        RbpRecord("IGF2BP3", "IGF2BP", "3utr_stability_translation"),
        RbpRecord("PUM1", "Pumilio", "3utr_stability_translation"),
        RbpRecord("PUM2", "Pumilio", "3utr_stability_translation"),
        RbpRecord("ZFP36", "TTP", "3utr_stability_translation"),
        RbpRecord("ZFP36L1", "TTP", "3utr_stability_translation"),
        RbpRecord("ZFP36L2", "TTP", "3utr_stability_translation"),
        RbpRecord("TIA1", "TIA/TIAR", "3utr_stability_translation"),
        RbpRecord("TIAL1", "TIA/TIAR", "3utr_stability_translation"),
        RbpRecord("CPEB1", "CPEB", "3utr_stability_translation"),
        RbpRecord("CPEB2", "CPEB", "3utr_stability_translation"),
        RbpRecord("CPEB3", "CPEB", "3utr_stability_translation"),
        RbpRecord("CPEB4", "CPEB", "3utr_stability_translation"),
        RbpRecord("MSI1", "Musashi", "3utr_stability_translation"),
        RbpRecord("MSI2", "Musashi", "3utr_stability_translation"),
        RbpRecord("LIN28A", "LIN28", "3utr_stability_translation"),
        RbpRecord("LIN28B", "LIN28", "3utr_stability_translation"),
        RbpRecord("SYNCRIP", "hnRNP", "3utr_stability_translation"),
        RbpRecord("PCBP1", "PCBP", "3utr_stability_translation"),
        RbpRecord("PCBP2", "PCBP", "3utr_stability_translation"),
        RbpRecord("PCBP3", "PCBP", "3utr_stability_translation"),
        RbpRecord("PCBP4", "PCBP", "3utr_stability_translation"),
        RbpRecord("AGO1", "Argonaute", "mirna_effector"),
        RbpRecord("AGO2", "Argonaute", "mirna_effector"),
        RbpRecord("AGO3", "Argonaute", "mirna_effector"),
        RbpRecord("AGO4", "Argonaute", "mirna_effector"),
        RbpRecord("TNRC6A", "GW182", "mirna_effector"),
        RbpRecord("TNRC6B", "GW182", "mirna_effector"),
        RbpRecord("TNRC6C", "GW182", "mirna_effector"),
        RbpRecord("DICER1", "miRNA processing", "mirna_effector"),
        RbpRecord("DROSHA", "miRNA processing", "mirna_effector"),
        RbpRecord("DGCR8", "miRNA processing", "mirna_effector"),
        RbpRecord("TARBP1", "dsRBD", "mirna_effector"),
        RbpRecord("TARBP2", "dsRBD", "mirna_effector"),
        # APA and cleavage/polyadenylation machinery.
        RbpRecord("CPSF1", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CPSF2", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CPSF3", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CPSF4", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CPSF6", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CPSF7", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CSTF1", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CSTF2", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CSTF3", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("NUDT21", "CFIm", "apa_machinery"),
        RbpRecord("FIP1L1", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("PAPOLA", "poly(A) polymerase", "apa_machinery"),
        RbpRecord("PAPOLG", "poly(A) polymerase", "apa_machinery"),
        RbpRecord("PCF11", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("CLP1", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("SYMPK", "cleavage/polyadenylation", "apa_machinery"),
        RbpRecord("PABPN1", "poly(A)-binding", "apa_machinery"),
        RbpRecord("PABPC1", "poly(A)-binding", "3utr_stability_translation"),
        RbpRecord("PABPC4", "poly(A)-binding", "3utr_stability_translation"),
        # Broad RNA processing and stress granule RBPs often available in CLIP resources.
        RbpRecord("DDX1", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX3X", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX5", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX6", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX17", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX21", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX39A", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX39B", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX42", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX46", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX47", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX50", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX54", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DDX56", "DEAD-box helicase", "rna_processing"),
        RbpRecord("DHX9", "DEAH-box helicase", "rna_processing"),
        RbpRecord("DHX15", "DEAH-box helicase", "rna_processing"),
        RbpRecord("DHX30", "DEAH-box helicase", "rna_processing"),
        RbpRecord("DHX36", "DEAH-box helicase", "rna_processing"),
        RbpRecord("LARP1", "La-related", "3utr_stability_translation"),
        RbpRecord("LARP4", "La-related", "3utr_stability_translation"),
        RbpRecord("LARP6", "La-related", "3utr_stability_translation"),
        RbpRecord("LARP7", "La-related", "rna_processing"),
        RbpRecord("YBX1", "cold-shock", "rna_processing"),
        RbpRecord("RBM3", "RBM", "rna_processing"),
        RbpRecord("RBM4", "RBM", "splicing"),
        RbpRecord("RBM5", "RBM", "splicing"),
        RbpRecord("RBM6", "RBM", "splicing"),
        RbpRecord("RBM7", "RBM", "rna_processing"),
        RbpRecord("RBM8A", "RBM", "rna_processing"),
        RbpRecord("RBM10", "RBM", "splicing"),
        RbpRecord("RBM14", "RBM", "rna_processing"),
        RbpRecord("RBM15", "RBM", "rna_processing"),
        RbpRecord("RBM17", "RBM", "splicing"),
        RbpRecord("RBM20", "RBM", "splicing"),
        RbpRecord("RBM22", "RBM", "splicing"),
        RbpRecord("RBM24", "RBM", "splicing"),
        RbpRecord("RBM25", "RBM", "splicing"),
        RbpRecord("RBM26", "RBM", "splicing"),
        RbpRecord("RBM27", "RBM", "splicing"),
        RbpRecord("RBM28", "RBM", "rna_processing"),
        RbpRecord("RBM38", "RBM", "3utr_stability_translation"),
        RbpRecord("RBM39", "RBM", "splicing"),
        RbpRecord("RBM41", "RBM", "rna_processing"),
        RbpRecord("RBM42", "RBM", "splicing"),
        RbpRecord("RBM45", "RBM", "rna_processing"),
        RbpRecord("RBM47", "RBM", "splicing"),
        RbpRecord("FUS", "FET", "rna_processing"),
        RbpRecord("TAF15", "FET", "rna_processing"),
        RbpRecord("EWSR1", "FET", "rna_processing"),
        RbpRecord("SFPQ", "DBHS", "rna_processing"),
        RbpRecord("NONO", "DBHS", "rna_processing"),
        RbpRecord("PSPC1", "DBHS", "rna_processing"),
        RbpRecord("MATR3", "zinc finger", "rna_processing"),
        RbpRecord("KHDRBS1", "STAR/KH", "splicing"),
        RbpRecord("KHDRBS2", "STAR/KH", "splicing"),
        RbpRecord("KHDRBS3", "STAR/KH", "splicing"),
        RbpRecord("FMR1", "FMRP", "3utr_stability_translation"),
        RbpRecord("FXR1", "FMRP", "3utr_stability_translation"),
        RbpRecord("FXR2", "FMRP", "3utr_stability_translation"),
        RbpRecord("G3BP1", "stress granule", "3utr_stability_translation"),
        RbpRecord("G3BP2", "stress granule", "3utr_stability_translation"),
        RbpRecord("NCL", "nucleolin", "rna_processing"),
        RbpRecord("NPM1", "nucleophosmin", "rna_processing"),
        RbpRecord("CIRBP", "cold-shock", "3utr_stability_translation"),
    ]
    return _records_to_catalog(records)


def load_expression_tables(
    project_root: Path,
    cohorts: Mapping[str, str] = EXPRESSION_COHORTS,
) -> pd.DataFrame:
    """Load harmonic regression outputs for the requested cohorts."""
    rows = []
    base = Path(project_root) / "data" / "Merged" / "regression_windowed"
    for cohort_key, cohort_token in cohorts.items():
        path = base / f"all_genes_regression_periodic_{cohort_token}.txt"
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep=r"\s+", engine="python")
        frame.columns = [str(c).strip().strip('"') for c in frame.columns]
        rename = {
            "Gene": "gene_symbol",
            "ENSG": "ensg",
            "Amplitude": "amplitude",
            "p": "p_value",
            "q": "q_value",
            "acrophase": "acrophase_hours",
            "Group": "group",
        }
        frame = frame.rename(columns=rename)
        for col in ["gene_symbol", "ensg", "amplitude", "p_value", "q_value", "acrophase_hours", "group"]:
            if col not in frame.columns:
                frame[col] = np.nan
        frame = frame[["group", "gene_symbol", "ensg", "amplitude", "p_value", "q_value", "acrophase_hours"]].copy()
        frame["gene_symbol"] = frame["gene_symbol"].map(standardize_symbol)
        frame["ensg"] = frame["ensg"].map(strip_version)
        for col in ["amplitude", "p_value", "q_value", "acrophase_hours"]:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
        frame["cohort_key"] = cohort_key
        frame["cohort_label"] = cohort_token.replace("_", " ")
        frame["source_path"] = str(path)
        rows.append(frame)
    if not rows:
        return pd.DataFrame(
            columns=[
                "group",
                "gene_symbol",
                "ensg",
                "amplitude",
                "p_value",
                "q_value",
                "acrophase_hours",
                "cohort_key",
                "cohort_label",
                "source_path",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def map_catalog_to_expression(catalog: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    """Fill catalog Ensembl IDs from expression tables and keep CLIP availability boolean."""
    catalog = catalog.copy()
    catalog["gene_symbol"] = catalog["gene_symbol"].map(standardize_symbol)
    expr_map = (
        expression.dropna(subset=["gene_symbol", "ensg"])
        .query("gene_symbol != '' and ensg != ''")
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")["ensg"]
        .to_dict()
    )
    existing_ensg = catalog["ensembl_id"] if "ensembl_id" in catalog.columns else pd.Series("", index=catalog.index)
    catalog["ensembl_id"] = existing_ensg.map(strip_version)
    catalog["ensembl_id"] = catalog.apply(
        lambda r: strip_version(r["ensembl_id"]) or expr_map.get(r["gene_symbol"], ""),
        axis=1,
    )
    if "clip_available" not in catalog.columns:
        catalog["clip_available"] = False
    catalog["clip_available"] = catalog["clip_available"].fillna(False).astype(bool)
    return catalog[RBP_COLUMNS].sort_values("gene_symbol").reset_index(drop=True)


def annotate_rbp_rhythmicity(
    catalog: pd.DataFrame,
    expression: pd.DataFrame,
    p_threshold: float = 0.05,
    q_threshold: float | None = 0.05,
) -> pd.DataFrame:
    """Join RBP catalog to expression rhythmicity results by gene symbol."""
    catalog = map_catalog_to_expression(catalog, expression)
    rbp_symbols = set(catalog["gene_symbol"])
    frame = expression[expression["gene_symbol"].isin(rbp_symbols)].copy()
    frame = frame.merge(catalog, on="gene_symbol", how="left", suffixes=("", "_catalog"))
    if "ensembl_id" in frame.columns:
        frame["ensembl_id"] = frame["ensembl_id"].where(frame["ensembl_id"].ne(""), frame["ensg"])
    else:
        frame["ensembl_id"] = frame["ensg"]
    frame["p_rhythmic"] = frame["p_value"] < p_threshold
    if q_threshold is None:
        frame["q_rhythmic"] = pd.NA
        frame["rhythmic"] = frame["p_rhythmic"]
        threshold_label = f"p<{p_threshold:g}"
    else:
        frame["q_rhythmic"] = frame["q_value"] < q_threshold
        frame["rhythmic"] = frame["p_rhythmic"] & frame["q_rhythmic"]
        threshold_label = f"p<{p_threshold:g};q<{q_threshold:g}"
    frame["rhythmicity_threshold"] = threshold_label
    keep = [
        "cohort_key",
        "cohort_label",
        "gene_symbol",
        "ensembl_id",
        "rbp_family",
        "functional_class",
        "clip_available",
        "amplitude",
        "p_value",
        "q_value",
        "acrophase_hours",
        "p_rhythmic",
        "q_rhythmic",
        "rhythmic",
        "rhythmicity_threshold",
    ]
    return frame[keep].sort_values(["cohort_key", "gene_symbol"]).reset_index(drop=True)


def summarize_rbp_rhythmicity(rbp_expression: pd.DataFrame) -> pd.DataFrame:
    """Summarize RBP rhythmicity by cohort and by requested tumor-normal comparisons."""
    rows = []
    for cohort_key, grp in rbp_expression.groupby("cohort_key", dropna=False):
        rows.append(
            {
                "summary_type": "cohort",
                "comparison": "",
                "cohort_key": cohort_key,
                "n_rbps_tested": int(grp["gene_symbol"].nunique()),
                "n_rhythmic": int(grp.loc[grp["rhythmic"], "gene_symbol"].nunique()),
                "rhythmic_fraction": float(grp.loc[grp["rhythmic"], "gene_symbol"].nunique() / max(grp["gene_symbol"].nunique(), 1)),
                "rbp_set": ",".join(sorted(grp.loc[grp["rhythmic"], "gene_symbol"].unique())),
            }
        )
    for comparison, cfg in COMPARISONS.items():
        tumor = rbp_expression[rbp_expression["cohort_key"].eq(cfg["tumor_key"])]
        normal = rbp_expression[rbp_expression["cohort_key"].eq(cfg["normal_key"])]
        tumor_set = set(tumor.loc[tumor["rhythmic"], "gene_symbol"])
        normal_set = set(normal.loc[normal["rhythmic"], "gene_symbol"])
        rows.extend(
            [
                {
                    "summary_type": "comparison_set",
                    "comparison": comparison,
                    "cohort_key": "tumor_only",
                    "n_rbps_tested": int(len(set(tumor["gene_symbol"]) | set(normal["gene_symbol"]))),
                    "n_rhythmic": int(len(tumor_set - normal_set)),
                    "rhythmic_fraction": np.nan,
                    "rbp_set": ",".join(sorted(tumor_set - normal_set)),
                },
                {
                    "summary_type": "comparison_set",
                    "comparison": comparison,
                    "cohort_key": "normal_only",
                    "n_rbps_tested": int(len(set(tumor["gene_symbol"]) | set(normal["gene_symbol"]))),
                    "n_rhythmic": int(len(normal_set - tumor_set)),
                    "rhythmic_fraction": np.nan,
                    "rbp_set": ",".join(sorted(normal_set - tumor_set)),
                },
                {
                    "summary_type": "comparison_set",
                    "comparison": comparison,
                    "cohort_key": "shared",
                    "n_rbps_tested": int(len(set(tumor["gene_symbol"]) | set(normal["gene_symbol"]))),
                    "n_rhythmic": int(len(tumor_set & normal_set)),
                    "rhythmic_fraction": np.nan,
                    "rbp_set": ",".join(sorted(tumor_set & normal_set)),
                },
            ]
        )
    primary = set(
        rbp_expression.loc[
            rbp_expression["cohort_key"].eq("TCGA_Tumor") & rbp_expression["rhythmic"],
            "gene_symbol",
        ]
    )
    validation = set(
        rbp_expression.loc[
            rbp_expression["cohort_key"].eq("CPTAC_Tumor") & rbp_expression["rhythmic"],
            "gene_symbol",
        ]
    )
    rows.append(
        {
            "summary_type": "validation",
            "comparison": "TCGA_Tumor_and_CPTAC_Tumor",
            "cohort_key": "validation_replicated",
            "n_rbps_tested": int(
                rbp_expression.loc[rbp_expression["cohort_key"].isin(["TCGA_Tumor", "CPTAC_Tumor"]), "gene_symbol"].nunique()
            ),
            "n_rhythmic": int(len(primary & validation)),
            "rhythmic_fraction": np.nan,
            "rbp_set": ",".join(sorted(primary & validation)),
        }
    )
    return pd.DataFrame(rows)


def phase_to_radians(hours: object) -> np.ndarray:
    """Convert circadian phase hours to radians on [0, 2*pi)."""
    return (pd.to_numeric(hours, errors="coerce").astype(float).to_numpy() % 24.0) / 24.0 * 2.0 * np.pi


def circular_lag_hours(target_phase: object, source_phase: object) -> np.ndarray:
    """Return target-source circular lag in hours, wrapped to [-12, 12]."""
    target = pd.to_numeric(target_phase, errors="coerce").astype(float)
    source = pd.to_numeric(source_phase, errors="coerce").astype(float)
    lag = ((target - source + 12.0) % 24.0) - 12.0
    return lag.to_numpy()


def phase_summary(phases: Iterable[object]) -> dict[str, float]:
    """Summarize circular concentration of 0-24h phases."""
    values = pd.to_numeric(pd.Series(list(phases)), errors="coerce").dropna()
    values = values[(values >= 0) & (values <= 24)]
    n = int(values.size)
    if n == 0:
        return {
            "n": 0,
            "mean_phase_hours": np.nan,
            "resultant_length": np.nan,
            "axial_r2": np.nan,
            "axis_hours": np.nan,
            "rayleigh_p_approx": np.nan,
        }
    theta = values.to_numpy() / 24.0 * 2.0 * np.pi
    c = np.cos(theta).mean()
    s = np.sin(theta).mean()
    mean_angle = np.arctan2(s, c) % (2.0 * np.pi)
    resultant = float(np.hypot(c, s))
    axial_theta = 2.0 * theta
    ac = np.cos(axial_theta).mean()
    ass = np.sin(axial_theta).mean()
    axial_r2 = float(ac * ac + ass * ass)
    axis_angle = (np.arctan2(ass, ac) / 2.0) % np.pi
    # Simple large-sample Rayleigh approximation, sufficient for QC summaries.
    rayleigh_p = float(np.exp(-n * resultant * resultant))
    return {
        "n": n,
        "mean_phase_hours": float(mean_angle / (2.0 * np.pi) * 24.0),
        "resultant_length": resultant,
        "axial_r2": axial_r2,
        "axis_hours": float(axis_angle / (2.0 * np.pi) * 24.0),
        "rayleigh_p_approx": rayleigh_p,
    }


def summarize_phase_by_group(df: pd.DataFrame, group_cols: Sequence[str], phase_col: str = "acrophase_hours") -> pd.DataFrame:
    """Return circular summaries for phase distributions by group."""
    rows = []
    if df.empty:
        return pd.DataFrame(columns=list(group_cols) + list(phase_summary([]).keys()))
    for keys, grp in df.groupby(list(group_cols), dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row.update(phase_summary(grp[phase_col]))
        rows.append(row)
    return pd.DataFrame(rows)


def _read_tsv_if_exists(path: Path) -> pd.DataFrame:
    if not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def load_apa_phase_tables(project_root: Path) -> pd.DataFrame:
    """Load APA-mRNA phase-lag tables for TCGA-vs-GTEx and CPTAC-vs-GTEx."""
    base = Path(project_root) / "data" / "polyAPA" / "comparison"
    paths = {
        "TCGA_vs_GTEx": base / "apa_expression_phase_lag_tcga_vs_gtex.tsv",
        "CPTAC_vs_GTEx": base / "apa_expression_phase_lag_cptac_vs_gtex.tsv",
    }
    frames = []
    for comparison, path in paths.items():
        frame = _read_tsv_if_exists(path)
        if frame.empty:
            continue
        frame["comparison"] = comparison
        frame["cohort_key"] = COMPARISONS[comparison]["tumor_key"]
        frame["normal_key"] = "GTEx_Normal"
        frame["host_ensg"] = frame.get("host_ensg", "").map(strip_version)
        for col in [
            "host_gene_symbol",
            "expr_gene_symbol",
            "apa_gene_symbol",
        ]:
            if col in frame.columns:
                frame[col] = frame[col].map(standardize_symbol)
        for col in [
            "apa_acrophase_hours",
            "expr_acrophase_hours",
            "phase_lag_hours",
            "amplitude",
            "relative_amplitude",
            "p_value",
            "fdr",
            "expr_amplitude",
            "expr_p",
            "expr_q",
        ]:
            if col in frame.columns:
                frame[col] = pd.to_numeric(frame[col], errors="coerce")
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_apa_membership(project_root: Path) -> pd.DataFrame:
    """Load APA host-gene set membership tables for matched gene backgrounds."""
    base = Path(project_root) / "data" / "polyAPA" / "comparison"
    paths = {
        "TCGA_vs_GTEx": base / "apa_host_gene_set_membership_tcga_vs_gtex.tsv",
        "CPTAC_vs_GTEx": base / "apa_host_gene_set_membership_cptac_vs_gtex.tsv",
    }
    frames = []
    for comparison, path in paths.items():
        frame = _read_tsv_if_exists(path)
        if frame.empty:
            continue
        frame["comparison"] = comparison
        frame["cohort_key"] = COMPARISONS[comparison]["tumor_key"]
        frame["normal_key"] = "GTEx_Normal"
        frame["host_ensg"] = frame.get("host_ensg", "").map(strip_version)
        for col in frame.columns:
            if col.endswith("host_gene_symbol") or col in {"host_gene_symbol"}:
                frame[col] = frame[col].map(standardize_symbol)
        if "host_gene_symbol" in frame.columns:
            frame["host_gene_symbol"] = frame["host_gene_symbol"].map(standardize_symbol)
        for col in frame.columns:
            if col.endswith("_apa_rhythmic"):
                frame[col] = frame[col].map(lambda x: str(x).strip().lower() in {"true", "1", "yes"})
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _load_se_universe_for_cohort(project_root: Path, cohort_dir: str) -> pd.DataFrame:
    path = Path(project_root) / "data" / "AS" / cohort_dir / "filtered_events" / f"{cohort_dir}_SE_variable_AS.csv"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    rename = {
        "ID": "event_id",
        "GeneID": "gene_id",
        "geneSymbol": "gene_symbol",
        "exonStart_0base": "exon_start",
        "exonEnd": "exon_end",
        "upstreamES": "upstream_start",
        "upstreamEE": "upstream_end",
        "downstreamES": "downstream_start",
        "downstreamEE": "downstream_end",
    }
    frame = frame.rename(columns=rename)
    for col in ["event_id", "gene_id", "gene_symbol", "chr", "strand"]:
        if col not in frame.columns:
            frame[col] = ""
    for col in ["exon_start", "exon_end", "upstream_start", "upstream_end", "downstream_start", "downstream_end"]:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame["event_id"] = frame["event_id"].astype(str)
    frame["gene_id"] = frame["gene_id"].map(strip_version)
    frame["gene_symbol"] = frame["gene_symbol"].map(standardize_symbol)
    frame["event_type"] = "SE"
    frame["cohort_dir"] = cohort_dir
    return frame


SE_CANONICAL_COHORTS = [
    ("TCGA", "TCGA_Tumor", "TCGA_vs_GTEx"),
    ("CPTAC", "CPTAC_Tumor", "CPTAC_vs_GTEx"),
]


def _read_se_canonical_rhythmicity(
    project_root: Path, cohort_dir: str, cohort_key: str
) -> pd.DataFrame:
    """Read the canonical (notebook 3.1.2) per-cohort SE harmonic-regression results.

    Source: data/AS/<cohort_dir>/rhythmic_results/<cohort_key>/SE/all_results.csv, the
    same SE rhythmicity profile used elsewhere in the manuscript. Holds every tested SE
    event with its rhythmicity p-value, BH-FDR, and amplitude. This replaces the older
    genome-wide AS table, whose pooled-FDR correction collapsed the CPTAC SE set to ~1.
    """
    path = (
        Path(project_root)
        / "data"
        / "AS"
        / cohort_dir
        / "rhythmic_results"
        / cohort_key
        / "SE"
        / "all_results.csv"
    )
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    if "event_id" not in frame.columns:
        return pd.DataFrame()
    frame["event_id"] = frame["event_id"].astype(str).map(standardize_event_id)
    for col in ["acrophase_hours", "p_value", "fdr", "qvalue", "amplitude", "relative_amplitude"]:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    if "fdr" not in frame.columns and "qvalue" in frame.columns:
        frame["fdr"] = frame["qvalue"]
    return frame


def load_se_universe(project_root: Path) -> pd.DataFrame:
    """Load the SE event background restricted to the canonical 3.1.2 tested set.

    Coordinates come from the filtered variable-SE table; the universe is limited to the
    events notebook 3.1.2 actually tested for rhythmicity, so the enrichment background
    and the rhythmic positive set share one SE event profile per cohort.
    """
    frames = []
    for cohort_dir, cohort_key, comparison in SE_CANONICAL_COHORTS:
        frame = _load_se_universe_for_cohort(project_root, cohort_dir)
        if frame.empty:
            continue
        frame["event_id"] = frame["event_id"].astype(str).map(standardize_event_id)
        canonical = _read_se_canonical_rhythmicity(project_root, cohort_dir, cohort_key)
        if not canonical.empty:
            tested = set(canonical["event_id"])
            frame = frame[frame["event_id"].isin(tested)].copy()
        frame["cohort_key"] = cohort_key
        frame["comparison"] = comparison
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_se_events(
    project_root: Path,
    p_threshold: float = 0.05,
    fdr_threshold: float | None = None,
) -> pd.DataFrame:
    """Load rhythmic SE events from the canonical 3.1.2 per-cohort SE results.

    The rhythmic positive set matches notebook 3.1.2 (nominal p < ``p_threshold``) unless
    ``fdr_threshold`` is supplied, in which case BH-FDR within SE is used instead.
    Coordinates are merged from the filtered SE universe by event_id.
    """
    universe = load_se_universe(project_root)
    if universe.empty:
        return pd.DataFrame()
    coord_cols = [
        "event_id",
        "gene_id",
        "gene_symbol",
        "chr",
        "strand",
        "exon_start",
        "exon_end",
        "upstream_start",
        "upstream_end",
        "downstream_start",
        "downstream_end",
    ]
    frames = []
    for cohort_dir, cohort_key, comparison in SE_CANONICAL_COHORTS:
        canonical = _read_se_canonical_rhythmicity(project_root, cohort_dir, cohort_key)
        if canonical.empty:
            continue
        if fdr_threshold is not None and "fdr" in canonical.columns:
            rhythmic = canonical[canonical["fdr"].le(fdr_threshold)].copy()
        else:
            rhythmic = canonical[canonical["p_value"].lt(p_threshold)].copy()
        coords = universe.loc[
            universe["cohort_dir"].eq(cohort_dir),
            [c for c in coord_cols if c in universe.columns],
        ].drop_duplicates("event_id")
        rhythmic = rhythmic.drop(
            columns=[c for c in ["gene_id", "gene_symbol", "chr", "strand"] if c in rhythmic.columns]
        )
        rhythmic = rhythmic.merge(coords, on="event_id", how="left")
        rhythmic["event_type"] = "SE"
        rhythmic["cohort_key"] = cohort_key
        rhythmic["comparison"] = comparison
        frames.append(rhythmic)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


CPTAC_SE_THRESHOLD_SOURCES = {
    "genome_wide": {
        "path": ("data", "AS", "CPTAC", "rhythmic_results", "genome_wide_rhythmic_AS_all.csv"),
        "label": "Genome-wide AS FDR",
        "detail": "Current main SE loader source; FDR adjusted across genome-wide AS calls.",
    },
    "se_specific": {
        "path": ("data", "AS", "CPTAC", "rhythmic_results", "SE", "CPTAC_SE_rhythmic_AS_all.csv"),
        "label": "SE-specific FDR",
        "detail": "Event-type-specific SE rhythmicity table; FDR adjusted within SE calls.",
    },
    "cptac_tumor_se": {
        "path": ("data", "AS", "CPTAC", "rhythmic_results", "CPTAC_Tumor", "SE", "all_results.csv"),
        "label": "CPTAC tumor SE",
        "detail": "CPTAC_Tumor SE rhythmicity table; CPTAC normal is not read.",
    },
}


CPTAC_SE_THRESHOLD_SPECS = [
    {"threshold_label": "FDR<=0.10", "threshold_family": "fdr", "fdr_le": 0.10},
    {"threshold_label": "FDR<=0.15", "threshold_family": "fdr", "fdr_le": 0.15},
    {"threshold_label": "FDR<=0.20", "threshold_family": "fdr", "fdr_le": 0.20},
    {"threshold_label": "FDR<=0.25", "threshold_family": "fdr", "fdr_le": 0.25},
    {"threshold_label": "p<0.05", "threshold_family": "p_value_only", "p_lt": 0.05},
    {"threshold_label": "p<0.10", "threshold_family": "p_value_only", "p_lt": 0.10},
    {"threshold_label": "p<0.20", "threshold_family": "p_value_only", "p_lt": 0.20},
    {"threshold_label": "p<0.10+rAMP>=0.15", "threshold_family": "p_value_amplitude", "p_lt": 0.10, "relative_amplitude_ge": 0.15},
    {"threshold_label": "p<0.10+amp>=0.05", "threshold_family": "p_value_amplitude", "p_lt": 0.10, "amplitude_ge": 0.05},
    {"threshold_label": "p<0.10+PSIrange>=0.20", "threshold_family": "p_value_psi", "p_lt": 0.10, "psi_range_ge": 0.20},
    {
        "threshold_label": "p<0.10+amp>=0.05+PSI",
        "threshold_family": "p_value_amplitude_psi",
        "p_lt": 0.10,
        "amplitude_ge": 0.05,
        "psi_range_ge": 0.20,
        "mean_psi_between": (0.05, 0.95),
    },
    {
        "threshold_label": "p<0.20+rAMP>=0.15+PSI",
        "threshold_family": "p_value_amplitude_psi",
        "p_lt": 0.20,
        "relative_amplitude_ge": 0.15,
        "psi_range_ge": 0.20,
        "mean_psi_between": (0.05, 0.95),
    },
]


def _read_cptac_se_threshold_source(project_root: Path, source_key: str, universe: pd.DataFrame) -> pd.DataFrame:
    """Read one CPTAC tumor SE rhythmicity source and align it to the SE universe."""
    source = CPTAC_SE_THRESHOLD_SOURCES[source_key]
    path = Path(project_root).joinpath(*source["path"])
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    event_col = "event_id" if "event_id" in frame.columns else ("ID" if "ID" in frame.columns else "")
    if not event_col:
        return pd.DataFrame()
    frame["event_id"] = frame[event_col].astype(str).map(standardize_event_id)
    for col in ["gene_id", "gene_symbol"]:
        if col not in frame.columns:
            frame[col] = ""
    frame["gene_id"] = frame["gene_id"].map(strip_version)
    frame["gene_symbol"] = frame["gene_symbol"].map(standardize_symbol)
    numeric_cols = [
        "acrophase_hours",
        "p_value",
        "fdr",
        "qvalue",
        "amplitude",
        "relative_amplitude",
        "mean_psi",
        "sd_psi",
        "min_psi",
        "max_psi",
        "psi_range",
        "r_squared",
    ]
    for col in numeric_cols:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    if "fdr" not in frame.columns and "qvalue" in frame.columns:
        frame["fdr"] = frame["qvalue"]
    if "qvalue" not in frame.columns and "fdr" in frame.columns:
        frame["qvalue"] = frame["fdr"]

    coords = universe[universe["cohort_dir"].eq("CPTAC")][
        [
            "event_id",
            "gene_id",
            "gene_symbol",
            "chr",
            "strand",
            "exon_start",
            "exon_end",
            "upstream_start",
            "upstream_end",
            "downstream_start",
            "downstream_end",
        ]
    ].copy()
    coords["event_id"] = coords["event_id"].astype(str).map(standardize_event_id)
    coords = coords.drop_duplicates("event_id")
    se_ids = set(coords["event_id"])
    frame = frame[frame["event_id"].isin(se_ids)].copy()
    drop_cols = [c for c in ["gene_id", "gene_symbol", "chr", "strand", "exonStart", "exonEnd"] if c in frame.columns]
    frame = frame.drop(columns=drop_cols).merge(coords, on="event_id", how="left")
    frame["event_type"] = "SE"
    frame["cohort_dir"] = "CPTAC"
    frame["cohort_key"] = "CPTAC_Tumor"
    frame["comparison"] = "CPTAC_vs_GTEx"
    frame["threshold_source"] = source_key
    frame["threshold_source_label"] = source["label"]
    frame["threshold_source_detail"] = source["detail"]
    frame["threshold_source_path"] = str(path)
    return frame


def _threshold_mask(frame: pd.DataFrame, spec: Mapping[str, object]) -> pd.Series:
    """Return a boolean mask for one rhythmicity threshold specification."""
    mask = pd.Series(True, index=frame.index)
    if "fdr_le" in spec:
        mask &= pd.to_numeric(frame.get("fdr", pd.Series(np.nan, index=frame.index)), errors="coerce").le(float(spec["fdr_le"]))
    if "p_lt" in spec:
        mask &= pd.to_numeric(frame.get("p_value", pd.Series(np.nan, index=frame.index)), errors="coerce").lt(float(spec["p_lt"]))
    if "relative_amplitude_ge" in spec:
        mask &= pd.to_numeric(frame.get("relative_amplitude", pd.Series(np.nan, index=frame.index)), errors="coerce").ge(float(spec["relative_amplitude_ge"]))
    if "amplitude_ge" in spec:
        mask &= pd.to_numeric(frame.get("amplitude", pd.Series(np.nan, index=frame.index)), errors="coerce").ge(float(spec["amplitude_ge"]))
    if "psi_range_ge" in spec:
        mask &= pd.to_numeric(frame.get("psi_range", pd.Series(np.nan, index=frame.index)), errors="coerce").ge(float(spec["psi_range_ge"]))
    if "mean_psi_between" in spec:
        lower, upper = spec["mean_psi_between"]
        values = pd.to_numeric(frame.get("mean_psi", pd.Series(np.nan, index=frame.index)), errors="coerce")
        mask &= values.between(float(lower), float(upper))
    return mask.fillna(False)


def load_cptac_se_threshold_sensitivity(
    project_root: Path,
    sources: Sequence[str] = ("genome_wide", "se_specific", "cptac_tumor_se"),
    threshold_specs: Sequence[Mapping[str, object]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load CPTAC tumor SE event sets across relaxed rhythmicity thresholds.

    The returned event table contains one row per selected event and threshold.
    The summary table records source-specific event counts and filter metadata.
    """
    universe = load_se_universe(project_root)
    universe = universe[universe["cohort_dir"].eq("CPTAC")].copy()
    if universe.empty:
        return pd.DataFrame(), pd.DataFrame()
    threshold_specs = list(threshold_specs or CPTAC_SE_THRESHOLD_SPECS)
    event_frames = []
    summary_rows = []
    for source_key in sources:
        if source_key not in CPTAC_SE_THRESHOLD_SOURCES:
            continue
        frame = _read_cptac_se_threshold_source(project_root, source_key, universe)
        if frame.empty:
            continue
        for rank, spec in enumerate(threshold_specs, start=1):
            label = str(spec["threshold_label"])
            mask = _threshold_mask(frame, spec)
            selected = frame[mask].copy()
            selected["threshold_rank"] = rank
            selected["threshold_label"] = label
            selected["threshold_family"] = str(spec.get("threshold_family", "custom"))
            selected["threshold_p_lt"] = spec.get("p_lt", np.nan)
            selected["threshold_fdr_le"] = spec.get("fdr_le", np.nan)
            selected["threshold_relative_amplitude_ge"] = spec.get("relative_amplitude_ge", np.nan)
            selected["threshold_amplitude_ge"] = spec.get("amplitude_ge", np.nan)
            selected["threshold_psi_range_ge"] = spec.get("psi_range_ge", np.nan)
            selected["threshold_mean_psi_between"] = (
                f"{spec['mean_psi_between'][0]}-{spec['mean_psi_between'][1]}"
                if "mean_psi_between" in spec
                else ""
            )
            if not selected.empty:
                event_frames.append(selected)
            summary_rows.append(
                {
                    "threshold_source": source_key,
                    "threshold_source_label": CPTAC_SE_THRESHOLD_SOURCES[source_key]["label"],
                    "threshold_source_detail": CPTAC_SE_THRESHOLD_SOURCES[source_key]["detail"],
                    "threshold_rank": rank,
                    "threshold_label": label,
                    "threshold_family": str(spec.get("threshold_family", "custom")),
                    "n_universe_events": int(frame["event_id"].nunique()),
                    "n_events": int(selected["event_id"].nunique()),
                    "n_genes": int(selected["gene_symbol"].nunique()) if "gene_symbol" in selected else 0,
                    "median_p_value": float(selected["p_value"].median()) if "p_value" in selected and selected["p_value"].notna().any() else np.nan,
                    "median_fdr": float(selected["fdr"].median()) if "fdr" in selected and selected["fdr"].notna().any() else np.nan,
                    "median_amplitude": float(selected["amplitude"].median()) if "amplitude" in selected and selected["amplitude"].notna().any() else np.nan,
                    "median_relative_amplitude": float(selected["relative_amplitude"].median()) if "relative_amplitude" in selected and selected["relative_amplitude"].notna().any() else np.nan,
                    "median_psi_range": float(selected["psi_range"].median()) if "psi_range" in selected and selected["psi_range"].notna().any() else np.nan,
                    "source_path": str(Path(project_root).joinpath(*CPTAC_SE_THRESHOLD_SOURCES[source_key]["path"])),
                }
            )
    events = pd.concat(event_frames, ignore_index=True, sort=False) if event_frames else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = summary.sort_values(["threshold_source", "threshold_rank"]).reset_index(drop=True)
    return events, summary


def discover_clip_files(project_root: Path, clip_dir: Path | None = None) -> pd.DataFrame:
    """Return a manifest for CLIP/eCLIP/POSTAR files under data/Ref/clip_rbp."""
    clip_dir = Path(clip_dir) if clip_dir is not None else Path(project_root) / "data" / "Ref" / "clip_rbp"
    rows = []
    if not clip_dir.exists():
        return pd.DataFrame(
            [
                {
                    "clip_dir": str(clip_dir),
                    "clip_file": "",
                    "file_type": "missing_directory",
                    "source_guess": "",
                    "genome_build_guess": "",
                    "n_rows": 0,
                    "status": "CLIP directory not found; CLIP-supported analyses are gated off.",
                }
            ]
        )
    suffixes = {".bed", ".gz", ".tsv", ".csv", ".txt", ".narrowpeak", ".broadpeak"}
    for path in sorted(clip_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in suffixes and not any(str(path).lower().endswith(x) for x in [".bed.gz", ".tsv.gz", ".txt.gz"]):
            continue
        lower = path.name.lower()
        source_guess = "ENCODE_eCLIP" if "encode" in lower or "eclip" in lower else ("POSTAR3" if "postar" in lower else "CLIP")
        build_guess = "hg38" if "hg38" in lower or "grch38" in lower else ("hg19" if "hg19" in lower or "grch37" in lower else "")
        rows.append(
            {
                "clip_dir": str(clip_dir),
                "clip_file": str(path),
                "file_type": _guess_clip_file_type(path),
                "source_guess": source_guess,
                "genome_build_guess": build_guess,
                "n_rows": np.nan,
                "status": "discovered",
            }
        )
    if not rows:
        rows.append(
            {
                "clip_dir": str(clip_dir),
                "clip_file": "",
                "file_type": "no_supported_files",
                "source_guess": "",
                "genome_build_guess": "",
                "n_rows": 0,
                "status": "No supported CLIP files found; add POSTAR3/ENCODE files here.",
            }
        )
    return pd.DataFrame(rows)


def _guess_clip_file_type(path: Path) -> str:
    lower = str(path).lower()
    if "postar3" in lower and lower.endswith((".txt", ".txt.gz")):
        return "peak_bed"
    if lower.endswith((".bed", ".bed.gz", ".narrowpeak", ".broadpeak")):
        return "peak_bed"
    return "interaction_table"


def _read_any_table(path: Path, max_rows: int | None = None) -> pd.DataFrame:
    lower = str(path).lower()
    compression = "gzip" if lower.endswith(".gz") else None
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return pd.read_csv(path, nrows=max_rows, compression=compression)
    sep = "\t" if lower.endswith((".tsv", ".tsv.gz", ".bed", ".bed.gz", ".txt", ".txt.gz", ".narrowpeak", ".broadpeak")) else None
    return pd.read_csv(path, sep=sep, nrows=max_rows, compression=compression, comment="#", engine="python")


def _infer_rbp_from_filename(path: Path, catalog_symbols: set[str]) -> str:
    stem = Path(str(path).replace(".gz", "")).stem.upper()
    parts = stem.replace("-", "_").replace(".", "_").split("_")
    for part in parts:
        if part in catalog_symbols:
            return part
    for symbol in sorted(catalog_symbols, key=len, reverse=True):
        if symbol in stem:
            return symbol
    return ""


def load_clip_interactions(
    project_root: Path,
    catalog: pd.DataFrame,
    expression: pd.DataFrame,
    clip_dir: Path | None = None,
) -> pd.DataFrame:
    """Load CLIP-derived RBP-target interaction tables if present."""
    clip_dir = Path(clip_dir) if clip_dir is not None else Path(project_root) / "data" / "Ref" / "clip_rbp"
    if not clip_dir.exists():
        return pd.DataFrame(columns=CLIP_LINK_COLUMNS)
    catalog_symbols = set(catalog["gene_symbol"].map(standardize_symbol))
    symbol_to_ensg = (
        expression.dropna(subset=["gene_symbol", "ensg"])
        .assign(gene_symbol=lambda x: x["gene_symbol"].map(standardize_symbol), ensg=lambda x: x["ensg"].map(strip_version))
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")["ensg"]
        .to_dict()
    )
    rows = []
    for path in sorted(clip_dir.rglob("*")):
        if not path.is_file() or _guess_clip_file_type(path) == "peak_bed":
            continue
        try:
            table = _read_any_table(path)
        except Exception:
            continue
        table.columns = [str(c).strip() for c in table.columns]
        lower_cols = {c.lower(): c for c in table.columns}
        rbp_col = _pick_column(lower_cols, ["rbp", "rbp_symbol", "protein", "clip_rbp", "rbp_gene"])
        target_col = _pick_column(lower_cols, ["target", "target_gene", "target_symbol", "gene", "gene_symbol", "gene_name", "genename", "host_gene_symbol"])
        target_ensg_col = _pick_column(lower_cols, ["target_ensg", "target_ensembl", "ensembl", "gene_id", "ensg"])
        region_col = _pick_column(lower_cols, ["region", "annotation", "feature", "binding_region"])
        build_col = _pick_column(lower_cols, ["genome_build", "genome", "assembly"])
        if target_col is None and target_ensg_col is None:
            continue
        rbp_from_file = _infer_rbp_from_filename(path, catalog_symbols)
        source_guess = "POSTAR3" if "postar" in path.name.lower() else ("ENCODE_eCLIP" if "eclip" in path.name.lower() else "CLIP")
        for _, row in table.iterrows():
            rbp = standardize_symbol(row[rbp_col]) if rbp_col else rbp_from_file
            if rbp not in catalog_symbols:
                continue
            target_symbol = standardize_symbol(row[target_col]) if target_col else ""
            target_ensg = strip_version(row[target_ensg_col]) if target_ensg_col else symbol_to_ensg.get(target_symbol, "")
            if not target_symbol and target_ensg:
                target_symbol = _lookup_symbol_by_ensg(expression, target_ensg)
            if not target_symbol and not target_ensg:
                continue
            region = str(row[region_col]).strip() if region_col else "gene_level_clip"
            genome_build = str(row[build_col]).strip() if build_col else _guess_genome_build_from_name(path.name)
            rows.append(
                {
                    "rbp_symbol": rbp,
                    "target_gene_symbol": target_symbol,
                    "target_ensg": target_ensg,
                    "event_id": "",
                    "region": region,
                    "source": source_guess,
                    "evidence_type": "clip_interaction_table",
                    "genome_build": genome_build,
                    "clip_file": str(path),
                }
            )
    if not rows:
        return pd.DataFrame(columns=CLIP_LINK_COLUMNS)
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def _pick_column(lower_cols: Mapping[str, str], aliases: Sequence[str]) -> str | None:
    for alias in aliases:
        if alias.lower() in lower_cols:
            return lower_cols[alias.lower()]
    return None


def _lookup_symbol_by_ensg(expression: pd.DataFrame, ensg: str) -> str:
    hit = expression[expression["ensg"].map(strip_version).eq(strip_version(ensg))]
    if hit.empty:
        return ""
    return standardize_symbol(hit.iloc[0]["gene_symbol"])


def _guess_genome_build_from_name(name: str) -> str:
    lower = name.lower()
    if "hg38" in lower or "grch38" in lower:
        return "hg38"
    if "hg19" in lower or "grch37" in lower:
        return "hg19"
    return ""


def _iter_bed_like_clip_peak_chunks(
    path: Path,
    catalog_symbols: set[str],
    chunksize: int = 1_000_000,
):
    """Yield normalized BED-like CLIP peak chunks."""
    lower = str(path).lower()
    compression = "gzip" if lower.endswith(".gz") else None
    is_postar3 = "postar3" in lower
    source = "POSTAR3" if is_postar3 else ("ENCODE_eCLIP" if "eclip" in lower or "encode" in lower else "CLIP_peak")
    genome_build = _guess_genome_build_from_name(path.name)
    if is_postar3 and not genome_build:
        genome_build = "hg38"
    rbp_from_file = _infer_rbp_from_filename(path, catalog_symbols)
    reader = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        compression=compression,
        chunksize=chunksize,
        low_memory=False,
    )
    for chunk in reader:
        if chunk.shape[1] < 3:
            continue
        peak = pd.DataFrame(
            {
                "chrom": chunk.iloc[:, 0].astype(str),
                "start": pd.to_numeric(chunk.iloc[:, 1], errors="coerce"),
                "end": pd.to_numeric(chunk.iloc[:, 2], errors="coerce"),
            }
        )
        peak["strand"] = "."
        peak["rbp_symbol"] = rbp_from_file
        peak["method"] = ""
        peak["sample"] = ""
        peak["accession"] = ""
        peak["score"] = np.nan
        if is_postar3 and chunk.shape[1] >= 10:
            peak["strand"] = chunk.iloc[:, 4].astype(str).where(chunk.iloc[:, 4].astype(str).isin(["+", "-"]), ".")
            peak["rbp_symbol"] = chunk.iloc[:, 5].map(standardize_symbol)
            peak["method"] = chunk.iloc[:, 6].astype(str)
            peak["sample"] = chunk.iloc[:, 7].astype(str)
            peak["accession"] = chunk.iloc[:, 8].astype(str)
            peak["score"] = pd.to_numeric(chunk.iloc[:, 9], errors="coerce")
        else:
            if chunk.shape[1] >= 4:
                maybe_rbp = chunk.iloc[:, 3].map(standardize_symbol)
                peak["rbp_symbol"] = maybe_rbp.where(maybe_rbp.isin(catalog_symbols), peak["rbp_symbol"])
            if chunk.shape[1] >= 6:
                peak["strand"] = chunk.iloc[:, 5].astype(str).where(chunk.iloc[:, 5].astype(str).isin(["+", "-"]), ".")
        peak["source"] = source
        peak["genome_build"] = genome_build
        peak["clip_file"] = str(path)
        peak = peak.dropna(subset=["start", "end"])
        peak = peak[peak["rbp_symbol"].isin(catalog_symbols)]
        if not peak.empty:
            yield peak


def _read_bed_like_clip_peaks(
    path: Path,
    catalog_symbols: set[str],
    chunksize: int = 1_000_000,
) -> pd.DataFrame:
    """Read BED-like CLIP peaks, including POSTAR3's 10-column peak table."""
    frames = list(_iter_bed_like_clip_peak_chunks(path, catalog_symbols, chunksize=chunksize))
    if not frames:
        return pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "rbp_symbol",
                "strand",
                "method",
                "sample",
                "accession",
                "score",
                "source",
                "genome_build",
                "clip_file",
            ]
        )
    return pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)


def load_clip_peaks(project_root: Path, catalog: pd.DataFrame, clip_dir: Path | None = None) -> pd.DataFrame:
    """Load BED-like CLIP peaks with normalized RBP symbols."""
    clip_dir = Path(clip_dir) if clip_dir is not None else Path(project_root) / "data" / "Ref" / "clip_rbp"
    columns = [
        "chrom",
        "start",
        "end",
        "rbp_symbol",
        "strand",
        "method",
        "sample",
        "accession",
        "score",
        "source",
        "genome_build",
        "clip_file",
    ]
    if not clip_dir.exists():
        return pd.DataFrame(columns=columns)
    catalog_symbols = set(catalog["gene_symbol"].map(standardize_symbol))
    frames = []
    for path in sorted(clip_dir.rglob("*")):
        if not path.is_file() or _guess_clip_file_type(path) != "peak_bed":
            continue
        try:
            peak = _read_bed_like_clip_peaks(path, catalog_symbols)
        except Exception:
            continue
        if not peak.empty:
            frames.append(peak)
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)


def load_utr_bed(project_root: Path) -> pd.DataFrame:
    """Load extracted hg38 3'UTR intervals used for CLIP peak-to-gene links."""
    path = Path(project_root) / "data" / "Ref" / "hg38_extracted_3UTR.bed"
    if not path.exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "target_gene_symbol", "strand", "region", "target_ensg"])
    bed = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
    parts = bed["name"].astype(str).str.split("|", expand=True)
    bed["transcript_id"] = parts[0].fillna("")
    bed["target_gene_symbol"] = parts[1].fillna("").map(standardize_symbol)
    bed["target_gene_symbol"] = bed["target_gene_symbol"].replace("NA", "")
    bed["target_ensg"] = ""
    bed["region"] = "3UTR"
    return bed[["chrom", "start", "end", "target_gene_symbol", "target_ensg", "strand", "region", "transcript_id"]]


def make_se_overlap_intervals(se_events: pd.DataFrame, intron_window_bp: int = 250) -> pd.DataFrame:
    """Create skipped-exon body and flanking intronic windows for CLIP overlap."""
    if se_events.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "target_gene_symbol", "target_ensg", "event_id", "strand", "region"])
    rows = []
    for _, row in se_events.iterrows():
        chrom = row.get("chr", "")
        start = pd.to_numeric(row.get("exon_start"), errors="coerce")
        end = pd.to_numeric(row.get("exon_end"), errors="coerce")
        if pd.isna(start) or pd.isna(end) or end <= start:
            continue
        event_id = str(row.get("event_id", ""))
        target_symbol = standardize_symbol(row.get("gene_symbol", ""))
        target_ensg = strip_version(row.get("gene_id", ""))
        strand = str(row.get("strand", "."))
        starts = [
            max(0, int(start) - intron_window_bp),
            int(start),
            int(end),
        ]
        ends = [
            int(start),
            int(end),
            int(end) + intron_window_bp,
        ]
        labels = ["SE_upstream_flank", "SE_exon_body", "SE_downstream_flank"]
        for interval_start, interval_end, label in zip(starts, ends, labels):
            if interval_end <= interval_start:
                continue
            rows.append(
                {
                    "chrom": chrom,
                    "start": interval_start,
                    "end": interval_end,
                    "target_gene_symbol": target_symbol,
                    "target_ensg": target_ensg,
                    "event_id": event_id,
                    "strand": strand,
                    "region": label,
                    "comparison": row.get("comparison", ""),
                    "cohort_key": row.get("cohort_key", ""),
                }
            )
    return pd.DataFrame(rows)


def overlap_clip_peaks_to_intervals(peaks: pd.DataFrame, intervals: pd.DataFrame) -> pd.DataFrame:
    """Overlap BED-like CLIP peaks with genomic intervals without requiring pybedtools."""
    if peaks.empty or intervals.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    rows = []
    peak = peaks.copy()
    ints = intervals.copy()
    for frame in [peak, ints]:
        frame["chrom"] = frame["chrom"].astype(str)
        frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
        frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    peak = peak.dropna(subset=["start", "end"])
    ints = ints.dropna(subset=["start", "end"])
    for chrom, chrom_peaks in peak.groupby("chrom"):
        chrom_intervals = ints[ints["chrom"].eq(chrom)]
        if chrom_intervals.empty:
            continue
        chrom_intervals = chrom_intervals.sort_values("start").reset_index(drop=True)
        interval_starts = chrom_intervals["start"].to_numpy()
        for _, p in chrom_peaks.iterrows():
            candidates = chrom_intervals.iloc[: np.searchsorted(interval_starts, p["end"], side="left")]
            candidates = candidates[candidates["end"].gt(p["start"])]
            if candidates.empty:
                continue
            for _, interval in candidates.iterrows():
                rows.append(
                    {
                        "rbp_symbol": p["rbp_symbol"],
                        "target_gene_symbol": interval.get("target_gene_symbol", ""),
                        "target_ensg": strip_version(interval.get("target_ensg", "")),
                        "event_id": str(interval.get("event_id", "")),
                        "region": interval.get("region", ""),
                        "source": p.get("source", "CLIP_peak"),
                        "evidence_type": "clip_peak_overlap",
                        "genome_build": p.get("genome_build", ""),
                        "chrom": p.get("chrom", ""),
                        "peak_start": p.get("start", np.nan),
                        "peak_end": p.get("end", np.nan),
                        "method": p.get("method", ""),
                        "sample": p.get("sample", ""),
                        "accession": p.get("accession", ""),
                        "score": p.get("score", np.nan),
                        "clip_file": p.get("clip_file", ""),
                        "cohort_key": interval.get("cohort_key", ""),
                        "comparison": interval.get("comparison", ""),
                    }
                )
    if not rows:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def _expand_intervals_to_bins(intervals: pd.DataFrame, bin_size: int = 50_000) -> pd.DataFrame:
    """Expand target intervals into genomic bins for fast chunked peak overlap."""
    if intervals.empty:
        return pd.DataFrame()
    frame = intervals.copy()
    frame["chrom"] = frame["chrom"].astype(str)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["start", "end"])
    frame = frame[frame["end"].gt(frame["start"])].copy()
    rows = []
    for idx, row in frame.reset_index(drop=True).iterrows():
        start_bin = int(row["start"]) // bin_size
        end_bin = max(start_bin, (int(row["end"]) - 1) // bin_size)
        for bin_id in range(start_bin, end_bin + 1):
            item = row.to_dict()
            item["interval_id"] = idx
            item["bin"] = bin_id
            rows.append(item)
    return pd.DataFrame(rows)


def _expand_peaks_to_bins(peaks: pd.DataFrame, bin_size: int = 50_000) -> pd.DataFrame:
    """Expand peak chunks into bins; most CLIP peaks stay in one bin."""
    if peaks.empty:
        return pd.DataFrame()
    frame = peaks.copy()
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["start", "end"])
    frame = frame[frame["end"].gt(frame["start"])].copy()
    start_bins = (frame["start"].astype(int) // bin_size).to_numpy()
    end_bins = ((frame["end"].astype(int) - 1) // bin_size).to_numpy()
    if np.array_equal(start_bins, end_bins):
        frame["bin"] = start_bins
        return frame
    frame["start_bin"] = start_bins
    frame["end_bin"] = end_bins
    single = frame[frame["start_bin"].eq(frame["end_bin"])].copy()
    single["bin"] = single["start_bin"]
    multi = frame[frame["start_bin"].ne(frame["end_bin"])].copy()
    if multi.empty:
        return single.drop(columns=["start_bin", "end_bin"])
    lengths = (multi["end_bin"] - multi["start_bin"] + 1).astype(int).to_numpy()
    repeated = multi.iloc[np.repeat(np.arange(len(multi)), lengths)].copy()
    repeated["bin"] = np.concatenate(
        [np.arange(start, end + 1, dtype=int) for start, end in zip(multi["start_bin"].astype(int), multi["end_bin"].astype(int))]
    )
    return pd.concat([single, repeated], ignore_index=True, sort=False).drop(columns=["start_bin", "end_bin"])


def _deduplicate_clip_links(links: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate peak-derived links while retaining the maximum available score."""
    if links.empty:
        return links
    for col in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]:
        if col not in links.columns:
            links[col] = "" if col != "score" else np.nan
    numeric_cols = [
        "score",
        "peak_start",
        "peak_end",
        "n_merged_peaks",
        "n_accessions",
        "n_methods",
        "n_samples",
        "max_score",
        "median_score",
        "peak_width_bp",
    ]
    for col in numeric_cols:
        if col in links.columns:
            links[col] = pd.to_numeric(links[col], errors="coerce")
    group_cols = [c for c in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"] if c != "score"]
    out = (
        links.groupby(group_cols, dropna=False, as_index=False)["score"]
        .max()
        .sort_values(["rbp_symbol", "region", "target_gene_symbol", "event_id"])
        .reset_index(drop=True)
    )
    return out


def _join_unique(values: Iterable[object], max_items: int = 12) -> str:
    """Join unique, non-missing values for compact evidence summaries."""
    cleaned = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if not text or text in {"nan", "None", "<NA>"}:
            continue
        cleaned.append(text)
    unique = sorted(set(cleaned))
    if len(unique) > max_items:
        return ",".join(unique[:max_items]) + f",...(+{len(unique) - max_items})"
    return ",".join(unique)


def _finalize_consensus_peak(
    base: Mapping[str, object],
    segment: pd.DataFrame,
    score_threshold: float,
    min_accessions: int,
    min_peak_support: int,
) -> dict | None:
    """Return one merged peak row if it passes the consensus evidence rule."""
    scores = pd.to_numeric(segment["score"], errors="coerce")
    accessions = [
        str(x).strip()
        for x in segment.get("accession", pd.Series(dtype=str))
        if not pd.isna(x) and str(x).strip() and str(x).strip() not in {"nan", "None", "<NA>"}
    ]
    methods = [
        str(x).strip()
        for x in segment.get("method", pd.Series(dtype=str))
        if not pd.isna(x) and str(x).strip() and str(x).strip() not in {"nan", "None", "<NA>"}
    ]
    samples = [
        str(x).strip()
        for x in segment.get("sample", pd.Series(dtype=str))
        if not pd.isna(x) and str(x).strip() and str(x).strip() not in {"nan", "None", "<NA>"}
    ]
    n_accessions = len(set(accessions))
    n_methods = len(set(methods))
    n_samples = len(set(samples))
    max_score = float(scores.max()) if scores.notna().any() else np.nan
    median_score = float(scores.median()) if scores.notna().any() else np.nan
    n_peaks = int(len(segment))

    filters = []
    if n_accessions >= min_accessions:
        filters.append(f"accessions_ge_{min_accessions}")
    if np.isfinite(max_score) and max_score >= score_threshold:
        filters.append(f"score_ge_{score_threshold:g}")
    if n_accessions == 0 and n_peaks >= min_peak_support:
        filters.append(f"peaks_ge_{min_peak_support}")
    if not filters:
        return None

    start = int(pd.to_numeric(segment["peak_start"], errors="coerce").min())
    end = int(pd.to_numeric(segment["peak_end"], errors="coerce").max())
    row = dict(base)
    row.update(
        {
            "evidence_type": "clip_peak_consensus_overlap",
            "peak_start": start,
            "peak_end": end,
            "method": _join_unique(methods),
            "sample": _join_unique(samples),
            "accession": _join_unique(accessions),
            "score": max_score,
            "clip_file": _join_unique(segment.get("clip_file", pd.Series(dtype=str)), max_items=4),
            "n_merged_peaks": n_peaks,
            "n_accessions": n_accessions,
            "n_methods": n_methods,
            "n_samples": n_samples,
            "max_score": max_score,
            "median_score": median_score,
            "peak_width_bp": end - start,
            "consensus_filter": ";".join(filters),
        }
    )
    return row


def build_consensus_clip_peak_links(
    raw_peak_links: pd.DataFrame,
    merge_distance_bp: int = CONSENSUS_MERGE_DISTANCE_BP,
    score_threshold: float = CONSENSUS_SCORE_THRESHOLD,
    min_accessions: int = CONSENSUS_MIN_ACCESSIONS,
    min_peak_support: int = CONSENSUS_MIN_PEAKS,
) -> pd.DataFrame:
    """Merge raw CLIP peak overlaps into high-confidence target-level binding segments."""
    if raw_peak_links.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    links = raw_peak_links.copy()
    for col in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]:
        if col not in links.columns:
            links[col] = "" if col not in {"score", "peak_start", "peak_end"} else np.nan
    links["peak_start"] = pd.to_numeric(links["peak_start"], errors="coerce")
    links["peak_end"] = pd.to_numeric(links["peak_end"], errors="coerce")
    links["score"] = pd.to_numeric(links["score"], errors="coerce")
    links = links.dropna(subset=["peak_start", "peak_end"])
    links = links[links["peak_end"].gt(links["peak_start"])].copy()
    if links.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])

    text_cols = [
        "rbp_symbol",
        "target_gene_symbol",
        "target_ensg",
        "event_id",
        "region",
        "source",
        "genome_build",
        "chrom",
        "cohort_key",
        "comparison",
    ]
    for col in text_cols:
        links[col] = links[col].fillna("").astype(str)
    links["rbp_symbol"] = links["rbp_symbol"].map(standardize_symbol)
    links["target_gene_symbol"] = links["target_gene_symbol"].map(standardize_symbol)
    links["target_ensg"] = links["target_ensg"].map(strip_version)
    links["event_id"] = links["event_id"].map(standardize_event_id)

    group_cols = [
        "rbp_symbol",
        "target_gene_symbol",
        "target_ensg",
        "event_id",
        "region",
        "source",
        "genome_build",
        "chrom",
        "cohort_key",
        "comparison",
    ]
    links = links.sort_values(group_cols + ["peak_start", "peak_end"], kind="mergesort").reset_index(drop=True)
    group_change = links[group_cols].ne(links[group_cols].shift()).any(axis=1)
    group_id = group_change.cumsum()
    previous_group_max_end = links.groupby(group_id, sort=False)["peak_end"].cummax().shift()
    previous_group_max_end[group_change] = np.nan
    segment_change = group_change | links["peak_start"].gt(previous_group_max_end + merge_distance_bp)
    links["_segment_id"] = segment_change.cumsum()

    grouped = links.groupby(group_cols + ["_segment_id"], sort=False, dropna=False)
    out = grouped.agg(
        peak_start=("peak_start", "min"),
        peak_end=("peak_end", "max"),
        n_merged_peaks=("peak_start", "size"),
        max_score=("score", "max"),
        median_score=("score", "median"),
        n_accessions=("accession", lambda x: len(set(v for v in x.astype(str) if v and v not in {"nan", "None", "<NA>"}))),
        n_methods=("method", lambda x: len(set(v for v in x.astype(str) if v and v not in {"nan", "None", "<NA>"}))),
        n_samples=("sample", lambda x: len(set(v for v in x.astype(str) if v and v not in {"nan", "None", "<NA>"}))),
        method=("method", _join_unique),
        sample=("sample", _join_unique),
        accession=("accession", _join_unique),
        clip_file=("clip_file", lambda x: _join_unique(x, max_items=4)),
    ).reset_index()
    if out.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])

    out["max_score"] = pd.to_numeric(out["max_score"], errors="coerce")
    out["median_score"] = pd.to_numeric(out["median_score"], errors="coerce")
    out["n_accessions"] = pd.to_numeric(out["n_accessions"], errors="coerce").fillna(0).astype(int)
    out["n_merged_peaks"] = pd.to_numeric(out["n_merged_peaks"], errors="coerce").fillna(0).astype(int)
    score_pass = out["max_score"].ge(score_threshold).fillna(False)
    accession_pass = out["n_accessions"].ge(min_accessions)
    peak_support_pass = out["n_accessions"].eq(0) & out["n_merged_peaks"].ge(min_peak_support)
    out = out[score_pass | accession_pass | peak_support_pass].copy()
    if out.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])

    filters = []
    for _, row in out.iterrows():
        row_filters = []
        if row["n_accessions"] >= min_accessions:
            row_filters.append(f"accessions_ge_{min_accessions}")
        if pd.notna(row["max_score"]) and row["max_score"] >= score_threshold:
            row_filters.append(f"score_ge_{score_threshold:g}")
        if row["n_accessions"] == 0 and row["n_merged_peaks"] >= min_peak_support:
            row_filters.append(f"peaks_ge_{min_peak_support}")
        filters.append(";".join(row_filters))
    out["consensus_filter"] = filters
    out["evidence_type"] = "clip_peak_consensus_overlap"
    out["score"] = out["max_score"]
    out["peak_width_bp"] = pd.to_numeric(out["peak_end"], errors="coerce") - pd.to_numeric(out["peak_start"], errors="coerce")
    out = out.drop(columns=["_segment_id"], errors="ignore")
    for col in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]:
        if col not in out.columns:
            out[col] = "" if col not in {"score", "peak_start", "peak_end"} else np.nan
    return _deduplicate_clip_links(out[CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]])


def summarize_clip_filtering_qc(
    links: pd.DataFrame,
    raw_peak_links: pd.DataFrame | None = None,
    score_threshold: float = CONSENSUS_SCORE_THRESHOLD,
    min_accessions: int = CONSENSUS_MIN_ACCESSIONS,
    min_peak_support: int = CONSENSUS_MIN_PEAKS,
    merge_distance_bp: int = CONSENSUS_MERGE_DISTANCE_BP,
    cache_status: str = "",
) -> pd.DataFrame:
    """Build a compact manifest for the CLIP consensus filtering strategy."""
    raw_peak_links = raw_peak_links if raw_peak_links is not None else pd.DataFrame()
    rows = [
        {"metric": "binding_map_mode", "value": "target_interval_consensus_peak_overlap", "detail": "Raw CLIP peaks are overlapped to 3'UTR/SE intervals, merged, and filtered before target claims."},
        {"metric": "merge_distance_bp", "value": merge_distance_bp, "detail": "Nearby peaks within this distance are merged per RBP-target/event-region."},
        {"metric": "score_threshold", "value": score_threshold, "detail": "Consensus peaks pass when max POSTAR3/confidence score is at least this value."},
        {"metric": "min_accessions", "value": min_accessions, "detail": "Consensus peaks pass when this many distinct accessions support the merged segment."},
        {"metric": "min_peak_support_no_accession", "value": min_peak_support, "detail": "Fallback support rule for peak files without accession IDs."},
        {"metric": "cache_status", "value": cache_status, "detail": "Whether consensus links were loaded from cache or rebuilt."},
        {"metric": "raw_peak_overlap_links", "value": int(len(raw_peak_links)), "detail": "Raw peak-overlap rows before consensus filtering."},
        {"metric": "final_clip_links", "value": int(len(links)), "detail": "CLIP links used by downstream APA/SE/mRNA analyses."},
    ]
    if not links.empty:
        evidence_counts = links["evidence_type"].fillna("").astype(str).value_counts()
        for evidence, count in evidence_counts.items():
            rows.append({"metric": f"evidence_type:{evidence}", "value": int(count), "detail": "Final link count by evidence type."})
        if "consensus_filter" in links.columns:
            exploded = (
                links["consensus_filter"]
                .fillna("")
                .astype(str)
                .str.split(";")
                .explode()
                .replace("", np.nan)
                .dropna()
                .value_counts()
            )
            for filt, count in exploded.items():
                rows.append({"metric": f"consensus_filter:{filt}", "value": int(count), "detail": "Final consensus link count by passing rule."})
        for source, grp in links.groupby("source"):
            rows.append({"metric": f"source:{source}:links", "value": int(len(grp)), "detail": "Final link count by CLIP source."})
            rows.append({"metric": f"source:{source}:rbps", "value": int(grp["rbp_symbol"].nunique()), "detail": "RBP count by CLIP source."})
    return pd.DataFrame(rows)


def filter_clip_links_by_support(clip_links: pd.DataFrame, support_mode: str = "consensus_all") -> pd.DataFrame:
    """Filter consensus CLIP links for support-threshold sensitivity analyses."""
    if clip_links.empty:
        return clip_links.copy()
    support_mode = support_mode.lower()
    links = clip_links.copy()
    for col in ["n_accessions", "n_merged_peaks", "max_score"]:
        if col not in links.columns:
            links[col] = np.nan
        links[col] = pd.to_numeric(links[col], errors="coerce")
    if support_mode in {"consensus_all", "all"}:
        mask = pd.Series(True, index=links.index)
    elif support_mode in {"multi_accession", "accessions_ge_2"}:
        mask = links["n_accessions"].ge(CONSENSUS_MIN_ACCESSIONS)
    elif support_mode in {"high_score", "score_ge_20"}:
        mask = links["max_score"].ge(CONSENSUS_SCORE_THRESHOLD)
    elif support_mode in {"strict", "strict_multi_accession_and_score"}:
        mask = links["n_accessions"].ge(CONSENSUS_MIN_ACCESSIONS) & links["max_score"].ge(CONSENSUS_SCORE_THRESHOLD)
    elif support_mode in {"multi_peak", "peaks_ge_2"}:
        mask = links["n_merged_peaks"].ge(CONSENSUS_MIN_PEAKS)
    else:
        raise ValueError(
            "support_mode must be one of consensus_all, multi_accession, high_score, strict_multi_accession_and_score, or multi_peak."
        )
    return links[mask.fillna(False)].copy().reset_index(drop=True)


def overlap_clip_peak_chunks_to_intervals(
    project_root: Path,
    catalog: pd.DataFrame,
    intervals: pd.DataFrame,
    clip_dir: Path | None = None,
    bin_size: int = 50_000,
    chunksize: int = 250_000,
) -> pd.DataFrame:
    """Stream CLIP peaks and return deduplicated overlaps to target intervals."""
    clip_dir = Path(clip_dir) if clip_dir is not None else Path(project_root) / "data" / "Ref" / "clip_rbp"
    if not clip_dir.exists() or intervals.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    interval_bins = _expand_intervals_to_bins(intervals, bin_size=bin_size)
    if interval_bins.empty:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    interval_bounds = (
        intervals.assign(
            chrom=lambda x: x["chrom"].astype(str),
            start=lambda x: pd.to_numeric(x["start"], errors="coerce"),
            end=lambda x: pd.to_numeric(x["end"], errors="coerce"),
        )
        .dropna(subset=["start", "end"])
        .groupby("chrom")
        .agg(min_start=("start", "min"), max_end=("end", "max"))
    )
    interval_chroms = set(interval_bounds.index)
    catalog_symbols = set(catalog["gene_symbol"].map(standardize_symbol))
    link_frames = []
    for path in sorted(clip_dir.rglob("*")):
        if not path.is_file() or _guess_clip_file_type(path) != "peak_bed":
            continue
        try:
            peak_iter = _iter_bed_like_clip_peak_chunks(path, catalog_symbols, chunksize=chunksize)
            for peak_chunk in peak_iter:
                peak_chunk = peak_chunk[peak_chunk["chrom"].astype(str).isin(interval_chroms)].copy()
                if peak_chunk.empty:
                    continue
                min_start = peak_chunk["chrom"].astype(str).map(interval_bounds["min_start"])
                max_end = peak_chunk["chrom"].astype(str).map(interval_bounds["max_end"])
                peak_chunk = peak_chunk[peak_chunk["end"].gt(min_start) & peak_chunk["start"].lt(max_end)].copy()
                if peak_chunk.empty:
                    continue
                peak_bins = _expand_peaks_to_bins(peak_chunk, bin_size=bin_size)
                if peak_bins.empty:
                    continue
                merged = peak_bins.merge(
                    interval_bins,
                    on=["chrom", "bin"],
                    how="inner",
                    suffixes=("_peak", "_interval"),
                )
                if merged.empty:
                    continue
                overlap = merged[
                    merged["start_peak"].lt(merged["end_interval"])
                    & merged["end_peak"].gt(merged["start_interval"])
                ].copy()
                if overlap.empty:
                    continue
                if "strand_interval" in overlap.columns and "strand_peak" in overlap.columns:
                    interval_strand = overlap["strand_interval"].astype(str)
                    peak_strand = overlap["strand_peak"].astype(str)
                    overlap = overlap[
                        interval_strand.isin(["+", "-"]).eq(False)
                        | peak_strand.isin(["+", "-"]).eq(False)
                        | interval_strand.eq(peak_strand)
                    ].copy()
                if overlap.empty:
                    continue
                links = pd.DataFrame(
                    {
                        "rbp_symbol": overlap["rbp_symbol"],
                        "target_gene_symbol": overlap.get("target_gene_symbol", ""),
                        "target_ensg": overlap.get("target_ensg", "").map(strip_version)
                        if "target_ensg" in overlap.columns
                        else "",
                        "event_id": overlap.get("event_id", "").astype(str)
                        if "event_id" in overlap.columns
                        else "",
                        "region": overlap.get("region", ""),
                        "source": overlap.get("source", ""),
                        "evidence_type": "clip_peak_overlap",
                        "genome_build": overlap.get("genome_build", ""),
                        "chrom": overlap.get("chrom", ""),
                        "peak_start": pd.to_numeric(overlap.get("start_peak", np.nan), errors="coerce"),
                        "peak_end": pd.to_numeric(overlap.get("end_peak", np.nan), errors="coerce"),
                        "method": overlap.get("method", ""),
                        "sample": overlap.get("sample", ""),
                        "accession": overlap.get("accession", ""),
                        "score": pd.to_numeric(overlap.get("score", np.nan), errors="coerce"),
                        "clip_file": overlap.get("clip_file", ""),
                        "cohort_key": overlap.get("cohort_key", ""),
                        "comparison": overlap.get("comparison", ""),
                    }
                )
                link_frames.append(_deduplicate_clip_links(links))
        except Exception:
            continue
    if not link_frames:
        return pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    return _deduplicate_clip_links(pd.concat(link_frames, ignore_index=True, sort=False))


def build_clip_links(
    project_root: Path,
    catalog: pd.DataFrame,
    expression: pd.DataFrame,
    se_events: pd.DataFrame | None = None,
    apa_membership: pd.DataFrame | None = None,
    clip_dir: Path | None = None,
    use_cache: bool = True,
    peak_evidence_mode: str = "consensus",
    consensus_score_threshold: float = CONSENSUS_SCORE_THRESHOLD,
    consensus_min_accessions: int = CONSENSUS_MIN_ACCESSIONS,
    consensus_min_peak_support: int = CONSENSUS_MIN_PEAKS,
    consensus_merge_distance_bp: int = CONSENSUS_MERGE_DISTANCE_BP,
    return_qc: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame] | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build all currently available CLIP-supported RBP-target links."""
    project_root = Path(project_root)
    peak_evidence_mode = peak_evidence_mode.lower()
    if peak_evidence_mode not in {"consensus", "raw"}:
        raise ValueError("peak_evidence_mode must be 'consensus' or 'raw'.")
    cache_suffix = "consensus" if peak_evidence_mode == "consensus" else "raw"
    cache_path = project_root / "data" / ANALYSIS_SLUG / f"clip_supported_rbp_target_links_{cache_suffix}.tsv"
    if use_cache and cache_path.exists() and cache_path.stat().st_size > 0:
        try:
            cached = pd.read_csv(cache_path, sep="\t", low_memory=False)
            if not cached.empty:
                for col in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]:
                    if col not in cached.columns:
                        cached[col] = "" if col != "score" else np.nan
                numeric_cols = {
                    "score",
                    "peak_start",
                    "peak_end",
                    "n_merged_peaks",
                    "n_accessions",
                    "n_methods",
                    "n_samples",
                    "max_score",
                    "median_score",
                    "peak_width_bp",
                }
                for col in [c for c in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"] if c not in numeric_cols]:
                    cached[col] = cached[col].fillna("").astype(str).replace({"nan": "", "None": "", "<NA>": ""})
                cached["event_id"] = cached["event_id"].map(standardize_event_id)
                for col in numeric_cols:
                    if col in cached.columns:
                        cached[col] = pd.to_numeric(cached[col], errors="coerce")
                catalog_for_coverage = map_catalog_to_expression(catalog, expression)
                coverage = summarize_clip_coverage(cached, catalog_for_coverage)
                qc = summarize_clip_filtering_qc(
                    cached,
                    raw_peak_links=pd.DataFrame(),
                    score_threshold=consensus_score_threshold,
                    min_accessions=consensus_min_accessions,
                    min_peak_support=consensus_min_peak_support,
                    merge_distance_bp=consensus_merge_distance_bp,
                    cache_status=f"loaded:{cache_path.name}",
                )
                if return_qc:
                    return cached, coverage, qc
                return cached, coverage
        except Exception:
            pass
    catalog = map_catalog_to_expression(catalog, expression)
    interactions = load_clip_interactions(project_root, catalog, expression, clip_dir=clip_dir)
    utr_intervals = load_utr_bed(project_root)
    symbol_to_ensg = (
        expression.dropna(subset=["gene_symbol", "ensg"])
        .assign(gene_symbol=lambda x: x["gene_symbol"].map(standardize_symbol), ensg=lambda x: x["ensg"].map(strip_version))
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")["ensg"]
        .to_dict()
    )
    if not utr_intervals.empty:
        utr_intervals["target_ensg"] = utr_intervals["target_gene_symbol"].map(symbol_to_ensg).fillna("")
        if apa_membership is not None and not apa_membership.empty:
            apa_symbols = set()
            for col in apa_membership.columns:
                if col.endswith("host_gene_symbol") or col == "host_gene_symbol":
                    apa_symbols.update(apa_membership[col].dropna().map(standardize_symbol))
            apa_symbols.discard("")
            if apa_symbols:
                utr_intervals = utr_intervals[utr_intervals["target_gene_symbol"].isin(apa_symbols)].copy()
    peak_intervals = []
    if not utr_intervals.empty:
        peak_intervals.append(utr_intervals)
    se_intervals = make_se_overlap_intervals(se_events if se_events is not None else pd.DataFrame())
    if not se_intervals.empty:
        peak_intervals.append(se_intervals)
    if peak_intervals:
        raw_cache_path = project_root / "data" / ANALYSIS_SLUG / "clip_supported_rbp_target_links_raw_peak_overlap.tsv"
        if use_cache and raw_cache_path.exists() and raw_cache_path.stat().st_size > 0:
            try:
                raw_peak_links = pd.read_csv(raw_cache_path, sep="\t", low_memory=False)
            except Exception:
                raw_peak_links = pd.DataFrame()
        else:
            raw_peak_links = pd.DataFrame()
        if raw_peak_links.empty:
            raw_peak_links = overlap_clip_peak_chunks_to_intervals(
                project_root,
                catalog,
                pd.concat(peak_intervals, ignore_index=True, sort=False),
                clip_dir=clip_dir,
            )
            try:
                write_tsv(raw_peak_links, raw_cache_path)
            except Exception:
                pass
        if peak_evidence_mode == "consensus":
            peak_links = build_consensus_clip_peak_links(
                raw_peak_links,
                merge_distance_bp=consensus_merge_distance_bp,
                score_threshold=consensus_score_threshold,
                min_accessions=consensus_min_accessions,
                min_peak_support=consensus_min_peak_support,
            )
        else:
            peak_links = raw_peak_links
    else:
        raw_peak_links = pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
        peak_links = pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    links = pd.concat([interactions, peak_links], ignore_index=True, sort=False)
    if links.empty:
        links = pd.DataFrame(columns=CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"])
    for col in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"]:
        if col not in links.columns:
            links[col] = "" if col != "score" else np.nan
    numeric_cols = {
        "score",
        "peak_start",
        "peak_end",
        "n_merged_peaks",
        "n_accessions",
        "n_methods",
        "n_samples",
        "max_score",
        "median_score",
        "peak_width_bp",
    }
    for col in [c for c in CLIP_ALL_LINK_COLUMNS + ["cohort_key", "comparison"] if c not in numeric_cols]:
        links[col] = links[col].fillna("").astype(str).replace({"nan": "", "None": "", "<NA>": ""})
    for col in numeric_cols:
        if col in links.columns:
            links[col] = pd.to_numeric(links[col], errors="coerce")
    links["event_id"] = links["event_id"].map(standardize_event_id)
    links["rbp_symbol"] = links["rbp_symbol"].map(standardize_symbol)
    links["target_gene_symbol"] = links["target_gene_symbol"].map(standardize_symbol)
    links["target_ensg"] = links["target_ensg"].map(strip_version)
    links = links.drop_duplicates().reset_index(drop=True)
    coverage = summarize_clip_coverage(links, catalog)
    qc = summarize_clip_filtering_qc(
        links,
        raw_peak_links=raw_peak_links,
        score_threshold=consensus_score_threshold,
        min_accessions=consensus_min_accessions,
        min_peak_support=consensus_min_peak_support,
        merge_distance_bp=consensus_merge_distance_bp,
        cache_status=f"rebuilt:{cache_path.name}",
    )
    try:
        write_tsv(links, cache_path)
    except Exception:
        pass
    if return_qc:
        return links, coverage, qc
    return links, coverage


def summarize_clip_coverage(links: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    """Summarize CLIP target counts by RBP and evidence region."""
    if links.empty:
        return pd.DataFrame(
            columns=[
                "rbp_symbol",
                "functional_class",
                "rbp_family",
                "n_clip_links",
                "n_target_genes",
                "n_target_events",
                "regions",
                "sources",
            ]
        )
    catalog_small = catalog[["gene_symbol", "functional_class", "rbp_family"]].rename(columns={"gene_symbol": "rbp_symbol"})
    rows = []
    for rbp, grp in links.groupby("rbp_symbol"):
        rows.append(
            {
                "rbp_symbol": rbp,
                "n_clip_links": int(len(grp)),
                "n_target_genes": int(grp.loc[grp["target_gene_symbol"].ne("") | grp["target_ensg"].ne(""), ["target_gene_symbol", "target_ensg"]].drop_duplicates().shape[0]),
                "n_target_events": int(grp.loc[grp.get("event_id", "").astype(str).ne(""), "event_id"].nunique()) if "event_id" in grp else 0,
                "regions": ",".join(sorted(grp["region"].dropna().astype(str).unique())),
                "sources": ",".join(sorted(grp["source"].dropna().astype(str).unique())),
            }
        )
    out = pd.DataFrame(rows)
    return out.merge(catalog_small, on="rbp_symbol", how="left").sort_values(["n_clip_links", "rbp_symbol"], ascending=[False, True])


def bh_fdr(p_values: Sequence[object]) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    p = pd.to_numeric(pd.Series(p_values), errors="coerce").to_numpy(dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[valid])]
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q[order] = np.clip(adjusted, 0, 1)
    return q


def _gene_key_series(df: pd.DataFrame, ensg_col: str, symbol_col: str) -> pd.Series:
    ensg = df.get(ensg_col, pd.Series("", index=df.index)).map(strip_version)
    symbol = df.get(symbol_col, pd.Series("", index=df.index)).map(standardize_symbol)
    return ensg.where(ensg.ne(""), symbol)


def compute_gene_target_enrichment(
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    universe_genes: Iterable[str],
    positive_genes: Iterable[str],
    comparison: str,
    cohort_key: str,
    region_keywords: Sequence[str] = ("3UTR", "gene_level_clip"),
    min_clip_targets: int = 3,
) -> pd.DataFrame:
    """Test whether positive genes are enriched among CLIP targets of rhythmic RBPs."""
    universe = {strip_version(g) if str(g).startswith("ENS") else standardize_symbol(g) for g in universe_genes if str(g)}
    positive = {strip_version(g) if str(g).startswith("ENS") else standardize_symbol(g) for g in positive_genes if str(g)}
    positive &= universe
    if not universe or clip_links.empty:
        return pd.DataFrame()
    links = clip_links.copy()
    links["target_key"] = _gene_key_series(links, "target_ensg", "target_gene_symbol")
    links = links[links["target_key"].isin(universe)]
    if region_keywords:
        pattern = "|".join([str(x) for x in region_keywords])
        links = links[links["region"].astype(str).str.contains(pattern, case=False, na=False)]
    rhythmic_rbps = set(
        rbp_expression.loc[
            rbp_expression["cohort_key"].eq(cohort_key) & rbp_expression["rhythmic"],
            "gene_symbol",
        ].map(standardize_symbol)
    )
    rows = []
    for rbp, grp in links.groupby("rbp_symbol"):
        if rbp not in rhythmic_rbps:
            continue
        targets = set(grp["target_key"]) & universe
        if len(targets) < min_clip_targets:
            continue
        a = len(targets & positive)
        b = len(targets - positive)
        c = len(positive - targets)
        d = len(universe - positive - targets)
        if fisher_exact is None:
            odds_ratio, p_value = np.nan, np.nan
        else:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "comparison": comparison,
                "cohort_key": cohort_key,
                "rbp_symbol": rbp,
                "n_universe_genes": len(universe),
                "n_positive_genes": len(positive),
                "n_clip_targets_in_universe": len(targets),
                "n_positive_clip_targets": a,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "regions_used": ",".join(region_keywords),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = bh_fdr(out["p_value"])
        out = out.sort_values(["fdr", "p_value", "n_positive_clip_targets"], ascending=[True, True, False])
    return out


def compute_event_target_enrichment(
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    universe_events: Iterable[str],
    positive_events: Iterable[str],
    comparison: str,
    cohort_key: str,
    region_keywords: Sequence[str] = ("SE_",),
    min_clip_events: int = 3,
) -> pd.DataFrame:
    """Test whether rhythmic SE events are enriched among CLIP-overlapped events."""
    universe = {str(x) for x in universe_events if str(x)}
    positive = {str(x) for x in positive_events if str(x)} & universe
    if not universe or clip_links.empty:
        return pd.DataFrame()
    links = clip_links.copy()
    links["event_id"] = links.get("event_id", "").astype(str)
    links = links[links["event_id"].isin(universe)]
    if region_keywords:
        pattern = "|".join([str(x) for x in region_keywords])
        links = links[links["region"].astype(str).str.contains(pattern, case=False, na=False)]
    rhythmic_rbps = set(
        rbp_expression.loc[
            rbp_expression["cohort_key"].eq(cohort_key) & rbp_expression["rhythmic"],
            "gene_symbol",
        ].map(standardize_symbol)
    )
    rows = []
    for rbp, grp in links.groupby("rbp_symbol"):
        if rbp not in rhythmic_rbps:
            continue
        targets = set(grp["event_id"]) & universe
        if len(targets) < min_clip_events:
            continue
        a = len(targets & positive)
        b = len(targets - positive)
        c = len(positive - targets)
        d = len(universe - positive - targets)
        if fisher_exact is None:
            odds_ratio, p_value = np.nan, np.nan
        else:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "comparison": comparison,
                "cohort_key": cohort_key,
                "rbp_symbol": rbp,
                "n_universe_events": len(universe),
                "n_positive_events": len(positive),
                "n_clip_events_in_universe": len(targets),
                "n_positive_clip_events": a,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "regions_used": ",".join(region_keywords),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = bh_fdr(out["p_value"])
        out = out.sort_values(["fdr", "p_value", "n_positive_clip_events"], ascending=[True, True, False])
    return out


def _normalize_target_ids(values: Iterable[object], id_type: str) -> pd.Series:
    """Normalize gene or event identifiers for set operations."""
    series = pd.Series(list(values), dtype="object")
    if id_type == "gene":
        return series.map(lambda x: strip_version(x) if str(x).startswith("ENS") else standardize_symbol(x))
    if id_type == "event":
        return series.map(standardize_event_id)
    raise ValueError("id_type must be 'gene' or 'event'.")


def _clip_target_key(links: pd.DataFrame, id_type: str) -> pd.Series:
    """Return the CLIP target key used for gene-level APA or event-level SE tests."""
    if id_type == "gene":
        return _gene_key_series(links, "target_ensg", "target_gene_symbol")
    if id_type == "event":
        return links.get("event_id", pd.Series("", index=links.index)).map(standardize_event_id)
    raise ValueError("id_type must be 'gene' or 'event'.")


def _phase_resultant_length(hours: Iterable[object]) -> float:
    """Return circular mean resultant length for phase-like hours."""
    return float(phase_summary(hours).get("resultant_length", np.nan))


def compute_target_background_contrast(
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    universe_ids: Iterable[object],
    positive_ids: Iterable[object],
    phase_table: pd.DataFrame,
    phase_id_col: str,
    phase_col: str,
    comparison: str,
    cohort_key: str,
    event_kind: str,
    id_type: str,
    region_keywords: Sequence[str],
    min_clip_targets: int = 3,
) -> pd.DataFrame:
    """Compare rhythmic event rates and phase concentration in CLIP targets vs non-target background."""
    universe = set(_normalize_target_ids(universe_ids, id_type))
    universe.discard("")
    positive = set(_normalize_target_ids(positive_ids, id_type)) & universe
    if not universe or clip_links.empty:
        return pd.DataFrame()

    links = clip_links.copy()
    links["target_key"] = _clip_target_key(links, id_type)
    links = links[links["target_key"].isin(universe)].copy()
    if region_keywords:
        pattern = "|".join([str(x) for x in region_keywords])
        links = links[links["region"].astype(str).str.contains(pattern, case=False, na=False)].copy()
    if links.empty:
        return pd.DataFrame()

    rhythmic_rbps = set(
        rbp_expression.loc[
            rbp_expression["cohort_key"].eq(cohort_key) & rbp_expression["rhythmic"],
            "gene_symbol",
        ].map(standardize_symbol)
    )
    phase = phase_table.copy()
    if phase.empty or phase_id_col not in phase.columns or phase_col not in phase.columns:
        phase = pd.DataFrame(columns=["target_key", "phase_hours"])
    else:
        phase["target_key"] = _normalize_target_ids(phase[phase_id_col], id_type)
        phase["phase_hours"] = pd.to_numeric(phase[phase_col], errors="coerce")
        phase = phase[phase["target_key"].isin(positive)].dropna(subset=["phase_hours"]).copy()

    rows = []
    for rbp, grp in links.groupby("rbp_symbol"):
        rbp = standardize_symbol(rbp)
        if rbp not in rhythmic_rbps:
            continue
        targets = set(grp["target_key"]) & universe
        if len(targets) < min_clip_targets:
            continue
        background = universe - targets
        target_positive = targets & positive
        background_positive = background & positive
        a = len(target_positive)
        b = len(targets - positive)
        c = len(background_positive)
        d = len(background - positive)
        if fisher_exact is None:
            odds_ratio, p_value = np.nan, np.nan
        else:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        target_phases = phase.loc[phase["target_key"].isin(target_positive), "phase_hours"]
        background_phases = phase.loc[phase["target_key"].isin(background_positive), "phase_hours"]
        target_phase_summary = phase_summary(target_phases)
        background_phase_summary = phase_summary(background_phases)
        rows.append(
            {
                "comparison": comparison,
                "cohort_key": cohort_key,
                "event_kind": event_kind,
                "rbp_symbol": rbp,
                "n_universe": len(universe),
                "n_positive": len(positive),
                "n_clip_targets": len(targets),
                "n_background": len(background),
                "n_positive_clip_targets": a,
                "n_positive_background": c,
                "target_rhythmic_rate": a / len(targets) if targets else np.nan,
                "background_rhythmic_rate": c / len(background) if background else np.nan,
                "rate_difference": (a / len(targets) if targets else np.nan) - (c / len(background) if background else np.nan),
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "target_phase_n": target_phase_summary["n"],
                "target_phase_resultant_length": target_phase_summary["resultant_length"],
                "target_phase_mean_hours": target_phase_summary["mean_phase_hours"],
                "background_phase_n": background_phase_summary["n"],
                "background_phase_resultant_length": background_phase_summary["resultant_length"],
                "phase_resultant_delta": target_phase_summary["resultant_length"] - background_phase_summary["resultant_length"]
                if pd.notna(target_phase_summary["resultant_length"]) and pd.notna(background_phase_summary["resultant_length"])
                else np.nan,
                "regions_used": ",".join(region_keywords),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = bh_fdr(out["p_value"])
        out = out.sort_values(["fdr", "p_value", "rate_difference"], ascending=[True, True, False])
    return out


def summarize_phase_lag_coherence(
    lag_df: pd.DataFrame,
    lag_col: str = "rbp_to_event_lag_hours",
    min_lags: int = 5,
    n_permutations: int = 500,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Summarize RBP-event lag concentration and compare it with a pooled phase-lag null."""
    if lag_df.empty or lag_col not in lag_df.columns:
        return pd.DataFrame()
    rng = np.random.default_rng(random_seed)
    frame = lag_df.copy()
    frame[lag_col] = pd.to_numeric(frame[lag_col], errors="coerce")
    frame = frame.dropna(subset=[lag_col])
    if frame.empty:
        return pd.DataFrame()

    rows = []
    for (comparison, event_kind), base in frame.groupby(["comparison", "event_kind"], dropna=False):
        pool = base[lag_col].dropna().to_numpy(dtype=float)
        if pool.size < min_lags:
            continue
        for rbp, grp in base.groupby("rbp_symbol", dropna=False):
            values = grp[lag_col].dropna().to_numpy(dtype=float)
            if values.size < min_lags:
                continue
            phase_values = values % 24.0
            observed = phase_summary(phase_values)
            null_r = []
            for _ in range(n_permutations):
                sampled = rng.choice(pool, size=values.size, replace=True) % 24.0
                null_r.append(_phase_resultant_length(sampled))
            null_r = np.asarray(null_r, dtype=float)
            p_value = (float((null_r >= observed["resultant_length"]).sum()) + 1.0) / (len(null_r) + 1.0)
            mean_lag = ((observed["mean_phase_hours"] + 12.0) % 24.0) - 12.0 if pd.notna(observed["mean_phase_hours"]) else np.nan
            rows.append(
                {
                    "comparison": comparison,
                    "event_kind": event_kind,
                    "rbp_symbol": standardize_symbol(rbp),
                    "lag_col": lag_col,
                    "n_lags": int(values.size),
                    "mean_lag_hours": mean_lag,
                    "median_lag_hours": float(np.median(values)),
                    "lag_resultant_length": observed["resultant_length"],
                    "null_mean_resultant_length": float(np.nanmean(null_r)) if null_r.size else np.nan,
                    "empirical_p_value": p_value,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["empirical_fdr"] = bh_fdr(out["empirical_p_value"])
        out = out.sort_values(["empirical_fdr", "empirical_p_value", "lag_resultant_length"], ascending=[True, True, False])
    return out


def compute_apa_phase_lags(
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    apa_phase: pd.DataFrame,
    comparison: str,
    cohort_key: str,
) -> pd.DataFrame:
    """Build rhythmic RBP -> CLIP-supported APA target -> rhythmic mRNA lag rows."""
    if clip_links.empty or apa_phase.empty:
        return pd.DataFrame()
    rbp_phase = rbp_expression[
        rbp_expression["cohort_key"].eq(cohort_key) & rbp_expression["rhythmic"]
    ][["gene_symbol", "acrophase_hours", "amplitude", "p_value", "q_value"]].copy()
    rbp_phase = rbp_phase.rename(
        columns={
            "gene_symbol": "rbp_symbol",
            "acrophase_hours": "rbp_acrophase_hours",
            "amplitude": "rbp_amplitude",
            "p_value": "rbp_p_value",
            "q_value": "rbp_q_value",
        }
    )
    events = apa_phase[apa_phase["comparison"].eq(comparison)].copy()
    events["target_key"] = _gene_key_series(events, "host_ensg", "host_gene_symbol")
    links = clip_links.copy()
    links["target_key"] = _gene_key_series(links, "target_ensg", "target_gene_symbol")
    links = links[
        links["region"].astype(str).str.contains("3UTR|gene_level_clip", case=False, na=False)
    ].copy()
    for col in ["n_accessions", "n_merged_peaks", "max_score", "median_score"]:
        if col not in links.columns:
            links[col] = np.nan
        links[col] = pd.to_numeric(links[col], errors="coerce")
    links = links.sort_values(["n_accessions", "n_merged_peaks", "max_score"], ascending=[False, False, False])
    links = links.drop_duplicates(["rbp_symbol", "target_key", "region", "source", "evidence_type", "genome_build"])
    joined = links.merge(events, on="target_key", how="inner", suffixes=("", "_apa"))
    joined = joined.merge(rbp_phase, on="rbp_symbol", how="inner")
    if joined.empty:
        return pd.DataFrame()
    joined["event_kind"] = "APA"
    joined["comparison"] = comparison
    joined["cohort_key"] = cohort_key
    joined["event_phase_hours"] = joined["apa_acrophase_hours"]
    joined["target_mrna_phase_hours"] = joined["expr_acrophase_hours"]
    joined["rbp_to_event_lag_hours"] = circular_lag_hours(joined["event_phase_hours"], joined["rbp_acrophase_hours"])
    joined["event_to_mrna_lag_hours"] = circular_lag_hours(joined["target_mrna_phase_hours"], joined["event_phase_hours"])
    keep = [
        "comparison",
        "cohort_key",
        "event_kind",
        "rbp_symbol",
        "target_key",
        "host_ensg",
        "host_gene_symbol",
        "apa_event_id",
        "region",
        "source",
        "evidence_type",
        "genome_build",
        "rbp_acrophase_hours",
        "event_phase_hours",
        "target_mrna_phase_hours",
        "rbp_to_event_lag_hours",
        "event_to_mrna_lag_hours",
        "rbp_amplitude",
        "amplitude",
        "expr_amplitude",
        "fdr",
        "expr_q",
        "clip_file",
        "n_merged_peaks",
        "n_accessions",
        "n_methods",
        "n_samples",
        "max_score",
        "median_score",
        "peak_width_bp",
        "consensus_filter",
    ]
    for col in keep:
        if col not in joined.columns:
            joined[col] = ""
    return joined[keep].drop_duplicates().reset_index(drop=True)


def compute_se_phase_lags(
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    se_events: pd.DataFrame,
    expression: pd.DataFrame,
    comparison: str,
    cohort_key: str,
) -> pd.DataFrame:
    """Build rhythmic RBP -> CLIP-supported SE target -> rhythmic mRNA lag rows."""
    if clip_links.empty or se_events.empty:
        return pd.DataFrame()
    rbp_phase = rbp_expression[
        rbp_expression["cohort_key"].eq(cohort_key) & rbp_expression["rhythmic"]
    ][["gene_symbol", "acrophase_hours", "amplitude", "p_value", "q_value"]].copy()
    rbp_phase = rbp_phase.rename(
        columns={
            "gene_symbol": "rbp_symbol",
            "acrophase_hours": "rbp_acrophase_hours",
            "amplitude": "rbp_amplitude",
            "p_value": "rbp_p_value",
            "q_value": "rbp_q_value",
        }
    )
    expr_phase = expression[
        expression["cohort_key"].eq(cohort_key) & expression["p_value"].lt(0.05)
    ][["gene_symbol", "ensg", "acrophase_hours", "amplitude", "q_value"]].copy()
    expr_phase = expr_phase.rename(
        columns={
            "gene_symbol": "expr_gene_symbol",
            "ensg": "gene_id",
            "acrophase_hours": "target_mrna_phase_hours",
            "amplitude": "expr_amplitude",
            "q_value": "expr_q",
        }
    )
    expr_phase["gene_id"] = expr_phase["gene_id"].map(strip_version)
    events = se_events[se_events["comparison"].eq(comparison)].copy()
    events["event_id"] = events["event_id"].astype(str)
    links = clip_links.copy()
    links["event_id"] = links.get("event_id", "").astype(str)
    links = links[links["region"].astype(str).str.contains("SE_", case=False, na=False)].copy()
    for col in ["n_accessions", "n_merged_peaks", "max_score", "median_score"]:
        if col not in links.columns:
            links[col] = np.nan
        links[col] = pd.to_numeric(links[col], errors="coerce")
    links = links.sort_values(["n_accessions", "n_merged_peaks", "max_score"], ascending=[False, False, False])
    links = links.drop_duplicates(["rbp_symbol", "event_id", "region", "source", "evidence_type", "genome_build"])
    joined = links.merge(events, on="event_id", how="inner", suffixes=("", "_se"))
    joined = joined.merge(rbp_phase, on="rbp_symbol", how="inner")
    joined = joined.merge(expr_phase, on="gene_id", how="left")
    if joined.empty:
        return pd.DataFrame()
    joined["event_kind"] = "SE"
    joined["comparison"] = comparison
    joined["cohort_key"] = cohort_key
    joined["event_phase_hours"] = joined["acrophase_hours"]
    joined["rbp_to_event_lag_hours"] = circular_lag_hours(joined["event_phase_hours"], joined["rbp_acrophase_hours"])
    joined["event_to_mrna_lag_hours"] = circular_lag_hours(joined["target_mrna_phase_hours"], joined["event_phase_hours"])
    keep = [
        "comparison",
        "cohort_key",
        "event_kind",
        "rbp_symbol",
        "event_id",
        "gene_id",
        "gene_symbol",
        "region",
        "source",
        "evidence_type",
        "genome_build",
        "rbp_acrophase_hours",
        "event_phase_hours",
        "target_mrna_phase_hours",
        "rbp_to_event_lag_hours",
        "event_to_mrna_lag_hours",
        "rbp_amplitude",
        "amplitude",
        "expr_amplitude",
        "fdr",
        "expr_q",
        "clip_file",
        "n_merged_peaks",
        "n_accessions",
        "n_methods",
        "n_samples",
        "max_score",
        "median_score",
        "peak_width_bp",
        "consensus_filter",
    ]
    for col in keep:
        if col not in joined.columns:
            joined[col] = ""
    return joined[keep].drop_duplicates().reset_index(drop=True)


def build_triplets(apa_lags: pd.DataFrame, se_lags: pd.DataFrame) -> pd.DataFrame:
    """Combine APA and SE lag rows and annotate validation-replicated triplets."""
    frames = []
    if apa_lags is not None and not apa_lags.empty:
        apa = apa_lags.copy()
        apa["target_gene_id"] = apa["host_ensg"].map(strip_version)
        apa["target_gene_symbol"] = apa["host_gene_symbol"].map(standardize_symbol)
        apa["event_id"] = apa["apa_event_id"].astype(str)
        frames.append(apa)
    if se_lags is not None and not se_lags.empty:
        se = se_lags.copy()
        se["target_gene_id"] = se["gene_id"].map(strip_version)
        se["target_gene_symbol"] = se["gene_symbol"].map(standardize_symbol)
        frames.append(se)
    if not frames:
        return pd.DataFrame()
    triplets = pd.concat(frames, ignore_index=True, sort=False)
    triplets["triplet_key"] = (
        triplets["rbp_symbol"].map(standardize_symbol)
        + "|"
        + triplets["event_kind"].astype(str)
        + "|"
        + triplets["target_gene_id"].fillna("").astype(str)
    )
    replicated = (
        triplets.groupby("triplet_key")["comparison"]
        .nunique()
        .rename("n_comparisons")
        .reset_index()
    )
    triplets = triplets.merge(replicated, on="triplet_key", how="left")
    triplets["validation_replicated"] = triplets["n_comparisons"].ge(2)
    lag_abs = triplets["rbp_to_event_lag_hours"].abs().fillna(12) + triplets["event_to_mrna_lag_hours"].abs().fillna(12)
    triplets["phase_lag_priority_score"] = triplets["validation_replicated"].astype(int) * 100 - lag_abs
    return triplets.sort_values(
        ["validation_replicated", "phase_lag_priority_score"],
        ascending=[False, False],
    ).reset_index(drop=True)


def empty_table(columns: Sequence[str]) -> pd.DataFrame:
    """Return an empty table with a stable schema."""
    return pd.DataFrame(columns=list(columns))
