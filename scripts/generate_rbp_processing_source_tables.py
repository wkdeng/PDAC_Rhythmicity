#!/usr/bin/env python3
"""Generate source tables for the current 3.4.3 RBP-processing figure section."""

from __future__ import annotations

from pathlib import Path
import argparse

import numpy as np
import pandas as pd

try:
    from scripts import rbp_rna_processing_links as rbp
except ModuleNotFoundError:  # pragma: no cover - supports direct script execution.
    import rbp_rna_processing_links as rbp


ANALYSIS = "rbp_rna_processing_links"

CURRENT_PANEL_SOURCE_TABLES = [
    "rbp_rhythmicity_summary.tsv",
    "rbp_expression_rhythmicity.tsv",
    "rbp_phase_circular_density.tsv",
    "rbp_phase_polarization_summary.tsv",
    "rbp_tumor_normal_rhythmic_overlap.tsv",
    "rbp_overlap_reactome_enrichment_top.tsv",
    "rbp_apa_clip_enrichment.tsv",
    "rbp_se_clip_enrichment.tsv",
    "rbp_target_background_contrast.tsv",
    "rbp_phase_lag_coherence.tsv",
    "rbp_event_gene_triplets.tsv",
    "clip_supported_rbp_target_links_consensus.tsv",
    "rbp_main_figure_candidate_summary.tsv",
    "rbp_apa_phase_lags.tsv",
    "rbp_se_phase_lags.tsv",
    "apa_expression_phase_lag_inputs.tsv",
    "apa_host_gene_membership_inputs.tsv",
    "se_variable_event_universe.tsv",
    "se_rhythmic_events_p0.05.tsv",
]


def _cohort_prefix(comparison: str) -> str:
    return "TCGA" if comparison == "TCGA_vs_GTEx" else "CPTAC"


def _rhythmic_mask(frame: pd.DataFrame, rhythmic_col: str) -> pd.Series:
    if rhythmic_col not in frame.columns:
        return pd.Series(False, index=frame.index)
    series = frame[rhythmic_col]
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce").fillna(0).ne(0)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def _save_table(
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    output_dir: Path,
    stage: str,
    filename: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    df = pd.DataFrame(df)
    path = rbp.write_tsv(df, output_dir / filename)
    tables[filename] = df
    manifest_rows.append(
        {
            "stage": stage,
            "table": filename,
            "path": str(path),
            "rows": int(df.shape[0]),
            "columns": int(df.shape[1]),
            "bytes": int(path.stat().st_size) if path.exists() else 0,
        }
    )
    return df


def circular_density(phases_hours: pd.Series, n_grid: int = 360, bw: float = 0.3) -> pd.DataFrame:
    """Von-Mises-like circular kernel density on a 0-24h phase scale."""
    values = pd.to_numeric(pd.Series(phases_hours), errors="coerce").dropna()
    values = values[(values >= 0) & (values <= 24)]
    if len(values) < 3:
        return pd.DataFrame(columns=["theta", "phase_hours", "density"])
    theta = (values.to_numpy(dtype=float) % 24.0) / 24.0 * 2.0 * np.pi
    grid = np.linspace(0, 2.0 * np.pi, n_grid + 1)[:-1]
    kappa = 1.0 / (float(bw) * float(bw))
    density = np.array([np.exp(kappa * np.cos(theta - g)).mean() for g in grid])
    density = density / (density.sum() * (2.0 * np.pi / n_grid))
    out = pd.DataFrame(
        {
            "theta": grid,
            "phase_hours": grid / (2.0 * np.pi) * 24.0,
            "density": density,
        }
    )
    return pd.concat(
        [
            out,
            pd.DataFrame(
                {
                    "theta": [2.0 * np.pi],
                    "phase_hours": [24.0],
                    "density": [out["density"].iloc[0]],
                }
            ),
        ],
        ignore_index=True,
    )


def axial_phase_summary(phases_hours: pd.Series) -> dict[str, float]:
    """Summarize antipodal phase orientation with the second trigonometric moment."""
    values = pd.to_numeric(pd.Series(phases_hours), errors="coerce").dropna()
    values = values[(values >= 0) & (values <= 24)]
    n = int(len(values))
    if n == 0:
        return {"n": 0, "axial_r2": np.nan, "axis_hours": np.nan, "opposite_axis_hours": np.nan}
    theta = (values.to_numpy(dtype=float) % 24.0) / 24.0 * 2.0 * np.pi
    z2 = np.mean(np.exp(1j * 2.0 * theta))
    axis_hours = ((np.angle(z2) / 2.0) % np.pi) * 24.0 / (2.0 * np.pi)
    return {
        "n": n,
        "axial_r2": float(np.abs(z2)),
        "axis_hours": float(axis_hours),
        "opposite_axis_hours": float((axis_hours + 12.0) % 24.0),
    }


def build_catalog_expression_tables(
    project_root: Path,
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    catalog = rbp.default_rbp_catalog()
    expression = rbp.load_expression_tables(project_root)
    catalog = rbp.map_catalog_to_expression(catalog, expression)
    rbp_expression = rbp.annotate_rbp_rhythmicity(catalog, expression, p_threshold=0.05, q_threshold=0.05)
    rbp_summary = rbp.summarize_rbp_rhythmicity(rbp_expression)
    rbp_phase_summary = rbp.summarize_phase_by_group(
        rbp_expression[rbp_expression["rhythmic"]],
        ["cohort_key", "functional_class"],
        phase_col="acrophase_hours",
    )

    forbidden = sorted([x for x in expression["cohort_key"].dropna().unique() if "CPTAC" in x and "Normal" in x])
    if forbidden:
        raise ValueError(f"CPTAC normal should not be present in this analysis: {forbidden}")

    _save_table(tables, manifest_rows, output_dir, "catalog_expression", "rbp_catalog.tsv", catalog)
    _save_table(tables, manifest_rows, output_dir, "catalog_expression", "expression_rhythmicity_inputs.tsv", expression)
    _save_table(tables, manifest_rows, output_dir, "catalog_expression", "rbp_expression_rhythmicity.tsv", rbp_expression)
    _save_table(tables, manifest_rows, output_dir, "catalog_expression", "rbp_rhythmicity_summary.tsv", rbp_summary)
    _save_table(tables, manifest_rows, output_dir, "catalog_expression", "rbp_phase_summary.tsv", rbp_phase_summary)
    return {
        "catalog": catalog,
        "expression": expression,
        "rbp_expression": rbp_expression,
        "rbp_summary": rbp_summary,
        "rbp_phase_summary": rbp_phase_summary,
    }


def build_event_input_tables(
    project_root: Path,
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
) -> dict[str, pd.DataFrame]:
    apa_phase = rbp.load_apa_phase_tables(project_root)
    apa_membership = rbp.load_apa_membership(project_root)
    se_universe = rbp.load_se_universe(project_root)
    # Rhythmic SE positive set now matches notebook 3.1.2 (canonical per-cohort SE,
    # nominal p < 0.05) rather than the older genome-wide-FDR table.
    se_events = rbp.load_se_events(project_root)

    input_summary = pd.DataFrame(
        [
            {
                "input": "APA phase lags",
                "n_rows": len(apa_phase),
                "n_genes_or_events": apa_phase.get("host_ensg", pd.Series(dtype=str)).nunique(),
            },
            {
                "input": "APA membership",
                "n_rows": len(apa_membership),
                "n_genes_or_events": apa_membership.get("host_ensg", pd.Series(dtype=str)).nunique(),
            },
            {
                "input": "SE universe",
                "n_rows": len(se_universe),
                "n_genes_or_events": se_universe.get("event_id", pd.Series(dtype=str)).nunique(),
            },
            {
                "input": "Rhythmic SE p0.05 (3.1.2 canonical)",
                "n_rows": len(se_events),
                "n_genes_or_events": se_events.get("event_id", pd.Series(dtype=str)).nunique(),
            },
        ]
    )

    _save_table(tables, manifest_rows, output_dir, "event_inputs", "apa_expression_phase_lag_inputs.tsv", apa_phase)
    _save_table(tables, manifest_rows, output_dir, "event_inputs", "apa_host_gene_membership_inputs.tsv", apa_membership)
    _save_table(tables, manifest_rows, output_dir, "event_inputs", "se_variable_event_universe.tsv", se_universe)
    _save_table(tables, manifest_rows, output_dir, "event_inputs", "se_rhythmic_events_p0.05.tsv", se_events)
    _save_table(tables, manifest_rows, output_dir, "event_inputs", "analysis_input_summary.tsv", input_summary)
    return {
        "apa_phase": apa_phase,
        "apa_membership": apa_membership,
        "se_universe": se_universe,
        "se_events": se_events,
        "input_summary": input_summary,
    }


def build_clip_source_tables(
    project_root: Path,
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    catalog: pd.DataFrame,
    expression: pd.DataFrame,
    se_universe: pd.DataFrame,
    apa_membership: pd.DataFrame,
    use_clip_cache: bool = True,
) -> dict[str, pd.DataFrame]:
    clip_manifest = rbp.discover_clip_files(project_root)
    clip_links, clip_coverage, clip_filtering_qc = rbp.build_clip_links(
        project_root,
        catalog=catalog,
        expression=expression,
        se_events=se_universe,
        apa_membership=apa_membership,
        use_cache=use_clip_cache,
        peak_evidence_mode="consensus",
        consensus_score_threshold=rbp.CONSENSUS_SCORE_THRESHOLD,
        consensus_min_accessions=rbp.CONSENSUS_MIN_ACCESSIONS,
        consensus_min_peak_support=rbp.CONSENSUS_MIN_PEAKS,
        consensus_merge_distance_bp=rbp.CONSENSUS_MERGE_DISTANCE_BP,
        return_qc=True,
    )
    catalog = catalog.copy()
    catalog["clip_available"] = catalog["gene_symbol"].isin(set(clip_coverage.get("rbp_symbol", [])))

    _save_table(tables, manifest_rows, output_dir, "clip_links", "rbp_catalog.tsv", catalog)
    _save_table(tables, manifest_rows, output_dir, "clip_links", "clip_source_manifest.tsv", clip_manifest)
    _save_table(tables, manifest_rows, output_dir, "clip_links", "clip_supported_rbp_target_links.tsv", clip_links)
    _save_table(tables, manifest_rows, output_dir, "clip_links", "clip_supported_rbp_target_links_consensus.tsv", clip_links)
    _save_table(tables, manifest_rows, output_dir, "clip_links", "clip_coverage_summary.tsv", clip_coverage)
    _save_table(tables, manifest_rows, output_dir, "clip_links", "clip_consensus_filtering_manifest.tsv", clip_filtering_qc)
    return {
        "catalog": catalog,
        "clip_manifest": clip_manifest,
        "clip_links": clip_links,
        "clip_coverage": clip_coverage,
        "clip_filtering_qc": clip_filtering_qc,
    }


def build_rbp_phase_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    rbp_expression: pd.DataFrame,
    bandwidth: float = 0.3,
) -> dict[str, pd.DataFrame]:
    cohort_order = ["TCGA_Tumor", "GTEx_Normal", "CPTAC_Tumor"]
    phase_df = rbp_expression[rbp_expression["rhythmic"]].dropna(subset=["acrophase_hours"]).copy()
    phase_density_frames = []
    phase_polarization_rows = []
    for cohort in cohort_order:
        values = phase_df.loc[phase_df["cohort_key"].eq(cohort), "acrophase_hours"]
        density = circular_density(values, bw=bandwidth)
        if not density.empty:
            density["cohort_key"] = cohort
            phase_density_frames.append(density)
        phase_polarization_rows.append({"cohort_key": cohort, **axial_phase_summary(values)})

    phase_density = (
        pd.concat(phase_density_frames, ignore_index=True, sort=False)
        if phase_density_frames
        else rbp.empty_table(["theta", "phase_hours", "density", "cohort_key"])
    )
    phase_polarization = pd.DataFrame(phase_polarization_rows)
    _save_table(tables, manifest_rows, output_dir, "rbp_phase", "rbp_phase_circular_density.tsv", phase_density)
    _save_table(tables, manifest_rows, output_dir, "rbp_phase", "rbp_phase_polarization_summary.tsv", phase_polarization)
    return {"rbp_phase_density": phase_density, "rbp_phase_polarization": phase_polarization}


def build_overlap_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    rbp_expression: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    overlap_rows = []
    overlap_gene_rows = []
    for comparison, cfg in rbp.COMPARISONS.items():
        tumor_key = cfg["tumor_key"]
        normal_key = cfg["normal_key"]
        tested = set(
            rbp_expression.loc[
                rbp_expression["cohort_key"].isin([tumor_key, normal_key]),
                "gene_symbol",
            ].map(rbp.standardize_symbol)
        )
        tumor_set = set(
            rbp_expression.loc[
                rbp_expression["cohort_key"].eq(tumor_key) & rbp_expression["rhythmic"],
                "gene_symbol",
            ].map(rbp.standardize_symbol)
        )
        normal_set = set(
            rbp_expression.loc[
                rbp_expression["cohort_key"].eq(normal_key) & rbp_expression["rhythmic"],
                "gene_symbol",
            ].map(rbp.standardize_symbol)
        )
        categories = {
            "shared": tumor_set & normal_set,
            "tumor_only": tumor_set - normal_set,
            "normal_only": normal_set - tumor_set,
            "not_rhythmic": tested - tumor_set - normal_set,
        }
        for category, genes in categories.items():
            overlap_rows.append(
                {
                    "comparison": comparison,
                    "tumor_key": tumor_key,
                    "normal_key": normal_key,
                    "category": category,
                    "n_rbps": len(genes),
                    "rbp_symbols": ",".join(sorted(genes)),
                }
            )
            for gene in sorted(genes):
                overlap_gene_rows.append(
                    {
                        "comparison": comparison,
                        "tumor_key": tumor_key,
                        "normal_key": normal_key,
                        "category": category,
                        "gene_symbol": gene,
                    }
                )

    overlap_summary = pd.DataFrame(overlap_rows)
    overlap_genes = pd.DataFrame(overlap_gene_rows)
    _save_table(tables, manifest_rows, output_dir, "overlap", "rbp_tumor_normal_rhythmic_overlap.tsv", overlap_summary)
    _save_table(tables, manifest_rows, output_dir, "overlap", "rbp_tumor_normal_rhythmic_overlap_genes.tsv", overlap_genes)
    return {"rbp_overlap_summary": overlap_summary, "rbp_overlap_genes": overlap_genes}


def build_reactome_overlap_tables(
    workspace_root: Path,
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    catalog: pd.DataFrame,
    rbp_overlap_genes: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    reactome_path = workspace_root / "SingleCellRhythmicity" / "data" / "Ref" / "Reactome_anno_HS.txt"
    reactome_cols = ["ensembl_id", "pathway_id", "url", "pathway_name", "evidence", "species"]
    category_focus = ["tumor_only", "shared", "normal_only"]
    category_label_map = {"tumor_only": "Tumor only", "shared": "Shared", "normal_only": "GTEx only"}

    if reactome_path.exists():
        reactome = pd.read_csv(reactome_path, sep="\t", names=reactome_cols, dtype=str)
        reactome["ensembl_id"] = reactome["ensembl_id"].map(rbp.strip_version)
        reactome = reactome[
            reactome["species"].eq("Homo sapiens")
            & reactome["ensembl_id"].ne("")
            & reactome["pathway_id"].notna()
            & reactome["pathway_name"].notna()
        ].drop_duplicates(["ensembl_id", "pathway_id"])
    else:
        reactome = rbp.empty_table(reactome_cols)

    rbp_ensembl_map = catalog[["gene_symbol", "ensembl_id"]].copy()
    rbp_ensembl_map["gene_symbol"] = rbp_ensembl_map["gene_symbol"].map(rbp.standardize_symbol)
    rbp_ensembl_map["ensembl_id"] = rbp_ensembl_map["ensembl_id"].map(rbp.strip_version)
    rbp_ensembl_map = rbp_ensembl_map[
        rbp_ensembl_map["gene_symbol"].ne("") & rbp_ensembl_map["ensembl_id"].ne("")
    ].drop_duplicates()
    ensembl_to_symbols = rbp_ensembl_map.groupby("ensembl_id")["gene_symbol"].apply(lambda x: sorted(set(x))).to_dict()
    reactome_background = set(rbp_ensembl_map["ensembl_id"]) & set(reactome.get("ensembl_id", pd.Series(dtype=str)))

    pathway_records = []
    if len(reactome) and reactome_background:
        for (pathway_id, pathway_name, url), grp in reactome.groupby(["pathway_id", "pathway_name", "url"], dropna=False):
            pathway_genes = set(grp["ensembl_id"]) & reactome_background
            if len(pathway_genes) < 2:
                continue
            pathway_records.append(
                {
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                    "url": url,
                    "pathway_background_genes": pathway_genes,
                    "pathway_background_size": len(pathway_genes),
                }
            )

    def symbols_for_ensembl(ids: set[str]) -> list[str]:
        symbols = []
        for ensg in sorted(set(ids)):
            symbols.extend(ensembl_to_symbols.get(ensg, []))
        return sorted(set(symbols))

    query_rows = []
    enrichment_frames = []
    for comparison in ["TCGA_vs_GTEx", "CPTAC_vs_GTEx"]:
        for category in category_focus:
            query_symbols = sorted(
                set(
                    rbp_overlap_genes.loc[
                        rbp_overlap_genes["comparison"].eq(comparison)
                        & rbp_overlap_genes["category"].eq(category),
                        "gene_symbol",
                    ].map(rbp.standardize_symbol)
                )
            )
            query_map = pd.DataFrame({"gene_symbol": query_symbols}).merge(rbp_ensembl_map, on="gene_symbol", how="left")
            query_ensembl = set(query_map["ensembl_id"].dropna().map(rbp.strip_version))
            query_ensembl = {x for x in query_ensembl if x} & reactome_background
            query_rows.append(
                {
                    "comparison": comparison,
                    "category": category,
                    "category_label": category_label_map[category],
                    "n_query_symbols": len(query_symbols),
                    "n_query_symbols_with_ensembl": int(
                        query_map["ensembl_id"].dropna().map(rbp.strip_version).ne("").sum()
                    ),
                    "n_query_rbps_in_reactome_background": len(query_ensembl),
                    "query_symbols": ",".join(query_symbols),
                    "query_symbols_in_reactome": ",".join(symbols_for_ensembl(query_ensembl)),
                }
            )

            rows = []
            universe = set(reactome_background)
            for pathway in pathway_records:
                pathway_genes = set(pathway["pathway_background_genes"])
                overlap_genes = query_ensembl & pathway_genes
                a = len(overlap_genes)
                b = len(query_ensembl - pathway_genes)
                c = len(pathway_genes - query_ensembl)
                d = len(universe - (query_ensembl | pathway_genes))
                if rbp.fisher_exact is None or len(query_ensembl) == 0:
                    odds_ratio = np.nan
                    p_value = np.nan
                else:
                    odds_ratio, p_value = rbp.fisher_exact([[a, b], [c, d]], alternative="greater")
                rows.append(
                    {
                        "comparison": comparison,
                        "category": category,
                        "category_label": category_label_map[category],
                        "pathway_id": pathway["pathway_id"],
                        "pathway_name": pathway["pathway_name"],
                        "url": pathway["url"],
                        "n_background_rbps": len(universe),
                        "n_query_rbps": len(query_ensembl),
                        "pathway_background_size": pathway["pathway_background_size"],
                        "overlap_count": a,
                        "odds_ratio": odds_ratio,
                        "p_value": p_value,
                        "overlap_ensembl": ",".join(sorted(overlap_genes)),
                        "overlap_symbols": ",".join(symbols_for_ensembl(overlap_genes)),
                    }
                )
            group_df = pd.DataFrame(rows)
            if len(group_df):
                group_df["fdr"] = rbp.bh_fdr(group_df["p_value"])
                enrichment_frames.append(group_df)

    query_sets = pd.DataFrame(query_rows)
    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True, sort=False)
        if enrichment_frames
        else rbp.empty_table(
            [
                "comparison",
                "category",
                "category_label",
                "pathway_id",
                "pathway_name",
                "url",
                "n_background_rbps",
                "n_query_rbps",
                "pathway_background_size",
                "overlap_count",
                "odds_ratio",
                "p_value",
                "fdr",
                "overlap_ensembl",
                "overlap_symbols",
            ]
        )
    )
    if len(enrichment):
        enrichment = enrichment.sort_values(
            ["comparison", "category", "p_value", "fdr", "overlap_count"],
            ascending=[True, True, True, True, False],
        )
        top = (
            enrichment[enrichment["overlap_count"].gt(0)]
            .sort_values(["comparison", "category", "p_value", "fdr", "overlap_count"], ascending=[True, True, True, True, False])
            .groupby(["comparison", "category"], dropna=False)
            .head(8)
            .reset_index(drop=True)
        )
    else:
        top = rbp.empty_table(enrichment.columns)

    manifest = pd.DataFrame(
        [
            {
                "source": "Reactome_anno_HS",
                "source_path": str(reactome_path),
                "source_exists": bool(reactome_path.exists()),
                "background": "RBP catalog genes mapped to Ensembl IDs and present in Reactome",
                "n_reactome_rows": int(len(reactome)),
                "n_reactome_pathways_tested": int(len(pathway_records)),
                "n_background_rbps": int(len(reactome_background)),
                "min_pathway_background_rbps": 2,
                "test": "one-sided Fisher exact test with BH FDR per comparison/category query set",
            }
        ]
    )

    _save_table(tables, manifest_rows, output_dir, "reactome_overlap", "rbp_overlap_reactome_query_sets.tsv", query_sets)
    _save_table(tables, manifest_rows, output_dir, "reactome_overlap", "rbp_overlap_reactome_enrichment.tsv", enrichment)
    _save_table(tables, manifest_rows, output_dir, "reactome_overlap", "rbp_overlap_reactome_enrichment_top.tsv", top)
    _save_table(tables, manifest_rows, output_dir, "reactome_overlap", "rbp_overlap_reactome_manifest.tsv", manifest)
    return {
        "rbp_overlap_pathway_query_sets": query_sets,
        "rbp_overlap_pathway_enrichment": enrichment,
        "rbp_overlap_pathway_top": top,
        "rbp_overlap_pathway_manifest": manifest,
    }


def build_apa_link_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    apa_membership: pd.DataFrame,
    apa_phase: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    enrichment_frames = []
    lag_frames = []
    for comparison, cfg in rbp.COMPARISONS.items():
        cohort_key = cfg["tumor_key"]
        membership = apa_membership[apa_membership["comparison"].eq(comparison)].copy()
        if membership.empty:
            continue
        rhythmic_col = f"{_cohort_prefix(comparison)}_apa_rhythmic"
        rhythmic = _rhythmic_mask(membership, rhythmic_col)
        universe = membership["host_ensg"].map(rbp.strip_version)
        positive = membership.loc[rhythmic, "host_ensg"].map(rbp.strip_version)
        enr = rbp.compute_gene_target_enrichment(
            clip_links=clip_links,
            rbp_expression=rbp_expression,
            universe_genes=universe,
            positive_genes=positive,
            comparison=comparison,
            cohort_key=cohort_key,
            region_keywords=("3UTR", "gene_level_clip"),
            min_clip_targets=3,
        )
        if not enr.empty:
            enrichment_frames.append(enr)
        lags = rbp.compute_apa_phase_lags(
            clip_links=clip_links,
            rbp_expression=rbp_expression,
            apa_phase=apa_phase,
            comparison=comparison,
            cohort_key=cohort_key,
        )
        if not lags.empty:
            lag_frames.append(lags)

    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True)
        if enrichment_frames
        else rbp.empty_table(
            [
                "comparison",
                "cohort_key",
                "rbp_symbol",
                "n_universe_genes",
                "n_positive_genes",
                "n_clip_targets_in_universe",
                "n_positive_clip_targets",
                "odds_ratio",
                "p_value",
                "fdr",
                "regions_used",
            ]
        )
    )
    lags = (
        pd.concat(lag_frames, ignore_index=True)
        if lag_frames
        else rbp.empty_table(
            [
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
            ]
        )
    )

    _save_table(tables, manifest_rows, output_dir, "apa_links", "rbp_apa_clip_enrichment.tsv", enrichment)
    _save_table(tables, manifest_rows, output_dir, "apa_links", "rbp_apa_phase_lags.tsv", lags)
    return {"apa_enrichment": enrichment, "apa_lags": lags}


def build_se_link_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    expression: pd.DataFrame,
    se_universe: pd.DataFrame,
    se_events: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    enrichment_frames = []
    lag_frames = []
    for comparison, cfg in rbp.COMPARISONS.items():
        cohort_key = cfg["tumor_key"]
        universe = se_universe.loc[se_universe["comparison"].eq(comparison), "event_id"].astype(str)
        positive = se_events.loc[se_events["comparison"].eq(comparison), "event_id"].astype(str)
        enr = rbp.compute_event_target_enrichment(
            clip_links=clip_links,
            rbp_expression=rbp_expression,
            universe_events=universe,
            positive_events=positive,
            comparison=comparison,
            cohort_key=cohort_key,
            region_keywords=("SE_",),
            min_clip_events=3,
        )
        if not enr.empty:
            enrichment_frames.append(enr)
        lags = rbp.compute_se_phase_lags(
            clip_links=clip_links,
            rbp_expression=rbp_expression,
            se_events=se_events,
            expression=expression,
            comparison=comparison,
            cohort_key=cohort_key,
        )
        if not lags.empty:
            lag_frames.append(lags)

    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True)
        if enrichment_frames
        else rbp.empty_table(
            [
                "comparison",
                "cohort_key",
                "rbp_symbol",
                "n_universe_events",
                "n_positive_events",
                "n_clip_events_in_universe",
                "n_positive_clip_events",
                "odds_ratio",
                "p_value",
                "fdr",
                "regions_used",
            ]
        )
    )
    lags = (
        pd.concat(lag_frames, ignore_index=True)
        if lag_frames
        else rbp.empty_table(
            [
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
            ]
        )
    )

    _save_table(tables, manifest_rows, output_dir, "se_links", "rbp_se_clip_enrichment.tsv", enrichment)
    _save_table(tables, manifest_rows, output_dir, "se_links", "rbp_se_phase_lags.tsv", lags)
    return {"se_enrichment": enrichment, "se_lags": lags}


def build_target_contrast_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    clip_links: pd.DataFrame,
    rbp_expression: pd.DataFrame,
    apa_membership: pd.DataFrame,
    apa_phase: pd.DataFrame,
    se_universe: pd.DataFrame,
    se_events: pd.DataFrame,
    apa_lags: pd.DataFrame,
    se_lags: pd.DataFrame,
    n_lag_permutations: int = 500,
) -> dict[str, pd.DataFrame]:
    apa_contrast_frames = []
    se_contrast_frames = []
    for comparison, cfg in rbp.COMPARISONS.items():
        cohort_key = cfg["tumor_key"]
        membership = apa_membership[apa_membership["comparison"].eq(comparison)].copy()
        if not membership.empty:
            rhythmic_col = f"{_cohort_prefix(comparison)}_apa_rhythmic"
            rhythmic = _rhythmic_mask(membership, rhythmic_col)
            apa_contrast = rbp.compute_target_background_contrast(
                clip_links=clip_links,
                rbp_expression=rbp_expression,
                universe_ids=membership["host_ensg"].map(rbp.strip_version),
                positive_ids=membership.loc[rhythmic, "host_ensg"].map(rbp.strip_version),
                phase_table=apa_phase[apa_phase["comparison"].eq(comparison)],
                phase_id_col="host_ensg",
                phase_col="apa_acrophase_hours",
                comparison=comparison,
                cohort_key=cohort_key,
                event_kind="APA",
                id_type="gene",
                region_keywords=("3UTR", "gene_level_clip"),
                min_clip_targets=3,
            )
            if not apa_contrast.empty:
                apa_contrast_frames.append(apa_contrast)

        se_contrast = rbp.compute_target_background_contrast(
            clip_links=clip_links,
            rbp_expression=rbp_expression,
            universe_ids=se_universe.loc[se_universe["comparison"].eq(comparison), "event_id"].astype(str),
            positive_ids=se_events.loc[se_events["comparison"].eq(comparison), "event_id"].astype(str),
            phase_table=se_events[se_events["comparison"].eq(comparison)],
            phase_id_col="event_id",
            phase_col="acrophase_hours",
            comparison=comparison,
            cohort_key=cohort_key,
            event_kind="SE",
            id_type="event",
            region_keywords=("SE_",),
            min_clip_targets=3,
        )
        if not se_contrast.empty:
            se_contrast_frames.append(se_contrast)

    apa_target_contrast = pd.concat(apa_contrast_frames, ignore_index=True) if apa_contrast_frames else rbp.empty_table([])
    se_target_contrast = pd.concat(se_contrast_frames, ignore_index=True) if se_contrast_frames else rbp.empty_table([])
    target_contrast = (
        pd.concat([apa_target_contrast, se_target_contrast], ignore_index=True, sort=False)
        if len(apa_target_contrast) or len(se_target_contrast)
        else rbp.empty_table([])
    )

    lag_input = pd.concat([apa_lags, se_lags], ignore_index=True, sort=False)
    lag_coherence = rbp.summarize_phase_lag_coherence(
        lag_input,
        lag_col="rbp_to_event_lag_hours",
        min_lags=5,
        n_permutations=n_lag_permutations,
        random_seed=42,
    )

    _save_table(tables, manifest_rows, output_dir, "target_contrast", "rbp_apa_target_background_contrast.tsv", apa_target_contrast)
    _save_table(tables, manifest_rows, output_dir, "target_contrast", "rbp_se_target_background_contrast.tsv", se_target_contrast)
    _save_table(tables, manifest_rows, output_dir, "target_contrast", "rbp_target_background_contrast.tsv", target_contrast)
    _save_table(tables, manifest_rows, output_dir, "target_contrast", "rbp_phase_lag_coherence.tsv", lag_coherence)
    return {
        "apa_target_contrast": apa_target_contrast,
        "se_target_contrast": se_target_contrast,
        "target_contrast": target_contrast,
        "lag_coherence": lag_coherence,
    }


def build_triplet_candidate_tables(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    manifest_rows: list[dict[str, object]],
    apa_enrichment: pd.DataFrame,
    se_enrichment: pd.DataFrame,
    apa_lags: pd.DataFrame,
    se_lags: pd.DataFrame,
    target_contrast: pd.DataFrame,
    lag_coherence: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    triplets = rbp.build_triplets(apa_lags, se_lags)
    enrichment_for_candidates = pd.concat(
        [apa_enrichment.assign(event_kind="APA"), se_enrichment.assign(event_kind="SE")],
        ignore_index=True,
        sort=False,
    )
    if enrichment_for_candidates.empty:
        candidates = rbp.empty_table([])
    else:
        enrichment_for_candidates["sig_fdr05"] = pd.to_numeric(
            enrichment_for_candidates["fdr"], errors="coerce"
        ).lt(0.05)
        candidate_rows = []
        for (event_kind, rbp_symbol), grp in enrichment_for_candidates.groupby(["event_kind", "rbp_symbol"], dropna=False):
            sig = grp[grp["sig_fdr05"]].copy()
            finite_or = pd.to_numeric(grp["odds_ratio"], errors="coerce").replace(np.inf, np.nan)
            candidate_rows.append(
                {
                    "event_kind": event_kind,
                    "rbp_symbol": rbp_symbol,
                    "n_sig_enrichment_comparisons": int(sig["comparison"].nunique()),
                    "sig_enrichment_comparisons": ",".join(sorted(sig["comparison"].dropna().astype(str).unique())),
                    "best_enrichment_fdr": float(pd.to_numeric(grp["fdr"], errors="coerce").min())
                    if pd.to_numeric(grp["fdr"], errors="coerce").notna().any()
                    else np.nan,
                    "max_enrichment_odds_ratio": float(finite_or.max()) if finite_or.notna().any() else np.nan,
                }
            )
        candidates = pd.DataFrame(candidate_rows)

        if not target_contrast.empty:
            contrast_summary = (
                target_contrast.groupby(["event_kind", "rbp_symbol"], dropna=False)
                .agg(
                    max_rate_difference=("rate_difference", "max"),
                    best_target_background_fdr=("fdr", "min"),
                    max_phase_resultant_delta=("phase_resultant_delta", "max"),
                )
                .reset_index()
            )
            candidates = candidates.merge(contrast_summary, on=["event_kind", "rbp_symbol"], how="left")

        if not lag_coherence.empty:
            coherence_summary = (
                lag_coherence.groupby(["event_kind", "rbp_symbol"], dropna=False)
                .agg(
                    best_lag_coherence_fdr=("empirical_fdr", "min"),
                    max_lag_resultant_length=("lag_resultant_length", "max"),
                    max_lag_count=("n_lags", "max"),
                )
                .reset_index()
            )
            candidates = candidates.merge(coherence_summary, on=["event_kind", "rbp_symbol"], how="left")

        if not triplets.empty:
            triplet_summary = (
                triplets.groupby(["event_kind", "rbp_symbol"], dropna=False)
                .agg(
                    n_triplet_rows=("triplet_key", "size"),
                    n_unique_triplets=("triplet_key", "nunique"),
                    n_replicated_rows=("validation_replicated", "sum"),
                    n_unique_targets=("target_gene_symbol", "nunique"),
                )
                .reset_index()
            )
            candidates = candidates.merge(triplet_summary, on=["event_kind", "rbp_symbol"], how="left")

        for col in ["n_triplet_rows", "n_unique_triplets", "n_replicated_rows", "n_unique_targets"]:
            if col not in candidates.columns:
                candidates[col] = 0
            candidates[col] = candidates[col].fillna(0).astype(int)
        candidates["evidence_tier"] = "exploratory"
        candidates.loc[
            candidates["event_kind"].eq("APA")
            & candidates["n_sig_enrichment_comparisons"].ge(2)
            & candidates["n_replicated_rows"].gt(0),
            "evidence_tier",
        ] = "main_candidate"
        candidates.loc[
            candidates["event_kind"].eq("SE") & candidates["n_sig_enrichment_comparisons"].ge(1),
            "evidence_tier",
        ] = "exploratory_se"
        tier_order = {"main_candidate": 0, "exploratory_se": 1, "exploratory": 2}
        candidates["evidence_tier_order"] = candidates["evidence_tier"].map(tier_order).fillna(9).astype(int)
        candidates = candidates.sort_values(
            ["evidence_tier_order", "n_sig_enrichment_comparisons", "n_replicated_rows", "best_enrichment_fdr"],
            ascending=[True, False, False, True],
        ).drop(columns=["evidence_tier_order"])

    _save_table(tables, manifest_rows, output_dir, "triplets_candidates", "rbp_event_gene_triplets.tsv", triplets)
    _save_table(tables, manifest_rows, output_dir, "triplets_candidates", "rbp_main_figure_candidate_summary.tsv", candidates)
    return {"triplets": triplets, "main_figure_candidates": candidates}


def generate_current_source_tables(
    project_root: str | Path,
    analysis: str = ANALYSIS,
    n_lag_permutations: int = 500,
    use_clip_cache: bool = True,
) -> dict[str, object]:
    """Generate current upstream source tables used by the RBP-processing panels."""
    project_root = Path(project_root)
    workspace_root = project_root.parent
    dirs = rbp.ensure_output_dirs(project_root, analysis)
    output_dir = dirs["output"]

    tables: dict[str, pd.DataFrame] = {}
    manifest_rows: list[dict[str, object]] = []

    core = build_catalog_expression_tables(project_root, output_dir, tables, manifest_rows)
    events = build_event_input_tables(project_root, output_dir, tables, manifest_rows)
    clip = build_clip_source_tables(
        project_root,
        output_dir,
        tables,
        manifest_rows,
        catalog=core["catalog"],
        expression=core["expression"],
        se_universe=events["se_universe"],
        apa_membership=events["apa_membership"],
        use_clip_cache=use_clip_cache,
    )
    build_rbp_phase_tables(output_dir, tables, manifest_rows, core["rbp_expression"])
    overlap = build_overlap_tables(output_dir, tables, manifest_rows, core["rbp_expression"])
    reactome = build_reactome_overlap_tables(
        workspace_root,
        output_dir,
        tables,
        manifest_rows,
        clip["catalog"],
        overlap["rbp_overlap_genes"],
    )
    apa = build_apa_link_tables(
        output_dir,
        tables,
        manifest_rows,
        clip["clip_links"],
        core["rbp_expression"],
        events["apa_membership"],
        events["apa_phase"],
    )
    se = build_se_link_tables(
        output_dir,
        tables,
        manifest_rows,
        clip["clip_links"],
        core["rbp_expression"],
        core["expression"],
        events["se_universe"],
        events["se_events"],
    )
    contrast = build_target_contrast_tables(
        output_dir,
        tables,
        manifest_rows,
        clip["clip_links"],
        core["rbp_expression"],
        events["apa_membership"],
        events["apa_phase"],
        events["se_universe"],
        events["se_events"],
        apa["apa_lags"],
        se["se_lags"],
        n_lag_permutations=n_lag_permutations,
    )
    triplets = build_triplet_candidate_tables(
        output_dir,
        tables,
        manifest_rows,
        apa["apa_enrichment"],
        se["se_enrichment"],
        apa["apa_lags"],
        se["se_lags"],
        contrast["target_contrast"],
        contrast["lag_coherence"],
    )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = rbp.write_tsv(manifest, output_dir / "rbp_processing_source_table_generation_manifest.tsv")

    missing_current = []
    for filename in CURRENT_PANEL_SOURCE_TABLES:
        path = output_dir / filename
        if not path.exists():
            missing_current.append(filename)
    if missing_current:
        raise FileNotFoundError(f"Missing current panel source tables after generation: {missing_current}")

    return {
        "manifest": manifest,
        "manifest_path": manifest_path,
        "tables": tables,
        "core": core,
        "events": events,
        "clip": clip,
        "overlap": overlap,
        "reactome": reactome,
        "apa": apa,
        "se": se,
        "contrast": contrast,
        "triplets": triplets,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--analysis", default=ANALYSIS)
    parser.add_argument("--n-lag-permutations", type=int, default=500)
    parser.add_argument("--no-clip-cache", action="store_true")
    args = parser.parse_args()
    results = generate_current_source_tables(
        project_root=args.project_root,
        analysis=args.analysis,
        n_lag_permutations=args.n_lag_permutations,
        use_clip_cache=not args.no_clip_cache,
    )
    manifest = results["manifest"]
    print(f"Wrote {len(manifest)} source-table records to {results['manifest_path']}")
    print(manifest[["stage", "table", "rows", "columns"]].to_string(index=False))


if __name__ == "__main__":
    main()
