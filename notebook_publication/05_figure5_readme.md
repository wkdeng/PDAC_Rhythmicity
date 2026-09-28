# Figure 5 and Supplementary Figures S14–S15

Run `05_figure5_rbp_processing.ipynb` from top to bottom. It embeds the RBP analysis and plotting implementation and does not import project scripts or execute the legacy notebook. The original source is `3.4.3_apa_se_rbp_rhythmic_connection.ipynb`; the detailed panel mapping is in `figure5_s14_s15_provenance.md`. The reference screenshots' S12/S13 are numbered S14/S15 here.

## Environment and inputs

Use Python 3.10 or later with NumPy, pandas, SciPy, matplotlib, statsmodels, IPython, Jupyter, nbformat, and nbclient. Pillow and ReportLab are needed only if invoking the included legacy image-export helpers. The notebook was run in the existing `jptsvr` environment. Exact versions are exported to `data/publication/fig05_rbp_processing/session_info.json`.

Keep the repository-relative input layout or set `CHRONO_PROJECT_ROOT` to its root. The notebook finds the root from the working directory or its parent; it contains no macOS-specific paths. `FIG05_REACTOME_FILE` can override the project-local Reactome annotation path.

Required inputs include:

- Expression rhythmicity and phases: `data/Merged/regression_windowed/all_genes_regression_periodic_{TCGA-PAAD_Tumor,CPTAC-3_Tumor,GTEx_Normal}.txt`.
- APA phase and membership results: `data/polyAPA/comparison/apa_expression_phase_lag_{tcga,cptac}_vs_gtex.tsv` and corresponding `apa_host_gene_set_membership_*` tables.
- Filtered SE events and existing SE rhythmicity results under `data/AS/{TCGA,CPTAC}/`.
- Sample-level expression, APA usage, and SE usage: `data/Merged/Group_TPM_*.csv`, `data/polyAPA/{TCGA,CPTAC}/filtered_apa/*_variable_apa_events.csv`, and SE filtered-event tables with `all_samples.txt`.
- External annotations: `data/Ref/hg38_extracted_3UTR.bed` and `data/Ref/Reactome_anno_HS.txt`. Raw rebuilding additionally needs `data/Ref/clip_rbp/POSTAR3_human_hg38.txt.gz`.
- Existing TargetScan-derived miRNA burden: `data/apa_mirna_rhythmic_connection/apa_mirna_gene_triads.tsv`.
- Default processed CLIP snapshot: `data/rbp_rna_processing_links/clip_supported_rbp_target_links_consensus.tsv`.

The starting expression/APA/SE inputs are processed analysis products downstream of RNA, DaPars2, and rMATS processing. `process.smk` does not directly produce every required file. Upstream CHIRAL inference and general APA/SE preprocessing are outside this notebook's boundary. A public repository must supply or document access to these inputs and external annotations; uploading the notebook alone does not supply the data.

`USE_REVIEWED_CLIP_SNAPSHOT=True` preserves the historical consensus-link input and records its fingerprint. This does not prove that the historical snapshot matches the present raw POSTAR3 reference. Setting it to `False` uses the included raw-overlap implementation; that substantially longer route was not used for the delivered run.

## Outputs and layout

- Notebook: `notebook_publication/05_figure5_rbp_processing.ipynb`.
- Regenerated analysis tables, panel source tables, input manifest, parameters, and software versions: `data/publication/fig05_rbp_processing/`.
- Figures: `data/figures/publication/fig05_rbp_processing/`.
- Execution timing and logs: `data/logs/publication/fig05_rbp_processing/`.

There are 24 individual panels and three assemblies, each exported as PNG and PDF and displayed inline. PNGs omit panel letters; PDFs include them. Filenames follow `Figure5A_rhythmic_rbp_landscape.png/pdf` and analogous content-based stems. One module is 52.5 mm; font roles are configured globally with a 7-pt minimum. Dense panels use wider modules to retain readable labels.

Figure 5 uses a four-module-wide grid: A–D; E/F/G; H/I; J/K; L/K; M/N. K spans two rows. S14 contains two full-width Reactome panels, each 210 × 52.5 mm (four modules), assembled as a 210 × 105 mm figure. S15 retains its eight original panel identities, including the repeated CPTAC phase distribution from Figure 5D.

## Differences from historical screenshots

The historical renderer preferred legacy aliases over newer canonical analysis tables. The old TCGA evidence matrix can be reconstructed from stale candidate, site-availability, and sample-model aliases. Its sample-model alias contains 23 models, versus 18 in the current canonical TCGA table; its site-availability alias has 10 significant RBPs, versus 21 in the canonical table. The old positional CLIP panel contains 10 rows, whereas the current canonical TCGA region-bias table contains eight.

This notebook uses newly regenerated canonical tables with the original analysis and six-criterion selection logic. Consequently, evidence-selected rows in Figure 5K/L and S15F/G, Figure 5M's positional support, and site-availability displays can differ from the screenshots. These are refreshed analysis inputs, not changes made to force a visual match. The current upstream APA phase input also differs in 39 of 9,236 rows from the older stored snapshot, propagating to phase-module summaries while preserving the top candidate ordering. Original files are preserved.

## Interpretation

CLIP supplies external binding evidence, not cohort-specific binding measurements. Phase links are associations. Sample-level OLS models preserve the original exploratory procedure, including shared samples and events; their p-values should not be interpreted as independent-replicate confirmation. Regulatory-site availability is an annotation-derived miRNA/CLIP burden and UTR-rank proxy, not measured accessibility. Reactome plots show nominal enrichment context; no displayed term passes FDR < 0.05. CPTAC and TCGA comparisons share GTEx and annotation inputs.

## Validation of the delivered run

The regenerated rhythmicity counts and CLIP enrichment tables match the previous canonical results. All 48 sample-level model identities (38 APA, 10 SE) and the top-ten order match; maximum absolute change in incremental R² is 1.04 × 10⁻¹⁶. The run evaluated about 5.13 million event–sample observations and 32,218 event-level effects. Analysis-stage timing is recorded separately from export timing.

`data/logs/publication/fig05_rbp_processing/final_validation.json` records the final notebook and export checks: successful code-cell execution, 27 inline PNGs, 27 PNG/PDF pairs, physical dimensions, PDF font sizes, panel letters, and page clipping.

## Annotation refinements (2026-09-23)

Phase panels show `n` and subscript `R₂` beneath the title, without axis-hour text. Enrichment threshold subtitles use regular weight. Figure 5K and S15F retain the criterion names but omit the explanatory axis sentence. These presentation edits use the existing canonical publication tables and preserve the original analysis input hashes.

## Statistical notation (2026-09-25)

Figure 5 and S14–S15 use lowercase italic *p* for nominal p-values and lowercase italic *q* for FDR-adjusted values, including transformed scales and thresholds. The displayed values and underlying calculations are unchanged.
