# Figure 4 and Supplementary Figures S10–S13

Run `04_figure4_apa_mirna_coordination.ipynb` with the Python 3 kernel. All statistical functions, reference parsing, panel drawing, layout, and export code are embedded. No legacy notebook or project script is imported at runtime.

## Thresholds and model

- APA events and miRNAs: **nominal P < 0.05** after the documented assay-specific QC.
- mRNA: **BH q < 0.05 and stored peak-to-trough `Amplitude > 0.1`** on **TMM-normalized log2 CPM**, using the exact versionless ENSG memberships from Figures 2/3. The notebook checks equality with each saved significant-gene table. This additional amplitude filter is part of the prior-figure contract.
- APA/miRNA use intercept + cosine + sine regression with a 24 h period, at least ten finite observations, and nonconstant rows. Raw P and BH q are both exported; selection uses P.
- APA uses untransformed PDUI; TCGA/CPTAC miRNA uses log2(CPM+1), and GTEx miRNA uses log2(TPM+1). Gene expression uses the supplied TMM logCPM matrix without another log transformation.
- Correlation, other category comparisons, and phase-lag polarization retain BH q < 0.05 within their stated test families. Panels D and I use raw two-sided Fisher P values. Panel H displays the lowest-P GO terms with raw P < 0.01, capped at 20 per cohort; BH q values remain available in its source tables.
- In-phase and anti-phase windows are ±3 h around 0 h and 12 h.

The original non-GEO TMM matrix and `Merged/CHIRAL_dual_penalty_phi_12h.csv` reproduce the saved regression phases and peak-to-trough amplitudes to numerical precision, including CXCL12. All modalities use this verified phase input. Gene selection and phase remain authoritative; a separate TMM fit supplies the displayed CXCL12 curve, with agreement exported. Host-gene and target joins use ENSG IDs; curated symbol panels choose a representative by lowest q then ENSG.

## Inputs and environment

The default environment is the existing Docker Jupyter runtime. Packages: numpy, pandas, scipy, statsmodels, matplotlib, and IPython. S10D additionally needs `pysam` when raw BAM extraction is enabled. Set `CHRONO_PROJECT_ROOT` and/or `CHRONO_DATA_DIR` for a different directory layout.

Required inputs, relative to the data directory:

- `Merged/CHIRAL_dual_penalty_phi_12h.csv`
- `Merged/GTEx_TCGA_merged_TMM_unfiltered_expression.csv`
- `Merged/GTEx_TCGA_merged_meta.csv`
- `Merged/regression/{all_genes,sig_genes}_regression_periodic_12h_{TCGA-PAAD_Tumor,CPTAC-3_Tumor,GTEx_Normal}.txt`
- `Ref/t2g.txt`
- `polyAPA/{TCGA,CPTAC,GTEx}/results/*/*_result_All_Prediction_Results.txt`
- `miRNA/TCGA/*.tsv` and `miRNA/CPTAC-3_Tumor/*.tsv`, containing sample-named processed GDC quantifications
- `miRNA/GTEx/raw/miRNA_TPM_matrix_PORTAL_2025_03_17.txt`
- `miRNA/GTEx/raw/rnacentral_mirbase.tsv`
- `Merged/GTEx_TCGA_merged_meta.csv`
- `miRNA/targets/predicted_targets.tsv`, the annotated processed TargetScan reference
- `enrichment/human_info.txt`, `enrichment/go_name.txt`, and `enrichment/goa_human.gpa`
- `Ref/GRCh38.primary_assembly.genome.fa` and its `.fai`, for the prespecified distal-UTR sequence check

The notebook accepts the existing cohort APA quantifications as processed inputs. Input manifests, sample mappings, source tables, and package versions are exported. Data and reference files are not embedded in the public notebook; distribute them separately or document access routes.

## S10D: TMEM183A BAM coverage

S10D is enabled and completes Supplementary Figure S10A–D. The first setup cell exposes `TMEM183A_BAM_DIR`, `BAM_PATHS`, and the GENCODE v49 GTF setting; local configured BAM paths take precedence. Each selected TCGA BAM must already have a readable `.bai` or `.csi` index. The notebook never creates an index.

`TMEM183A_USE_SOURCE_CACHE=True` is the public-reproduction default. It reads the bundled, hash-validated record of per-base depth extracted from the real indexed BAMs, together with its BAM/index provenance manifest and GENCODE v49 gene model. This permits a public rerun without redistributing the large controlled BAMs. Set it to `False` to force a fresh raw-BAM extraction after configuring the paths. The notebook verifies that the current, newly selected four samples agree with the manifest before it accepts the cached observations.

S10D writes `S10D_TMEM183A_coverage_manifest.tsv`, `S10D_TMEM183A_per_base_coverage.tsv`, and `S10D_TMEM183A_gene_model.tsv`. Coverage is observed, unsmoothed per-base depth across the 2,000-bp TMEM183A locus on chr1; exported coordinates are 1-based (`203022501`–`203024500`) while the `pysam` fetch interval is the equivalent 0-based half-open range. The gene model identifies TMEM183A transcript ENST00000965092 from GENCODE v49. No simulated or PDUI-derived coverage is used.

## Outputs

Figure 4D compares the TCGA ∪ CPTAC union with GTEx (a gene is rhythmic in the union if it passes the mRNA threshold in either tumor cohort, and is counted once) using the percentage of tested genes rhythmic in each APA-machinery category. Rhythmicity remains BH q < 0.05 and stored peak-to-trough amplitude > 0.1 from TMM log2 CPM. Each bracket uses a two-sided Fisher exact test of rhythmic/nonrhythmic gene counts in the tumor union versus GTEx. Raw-P annotations are * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001; nonsignificant comparisons have no bracket or annotation. The displayed P values are unadjusted. Source counts, percentages, and P values are in `Figure4D_union_gtex_apa_machinery_expression_source.tsv`.

- Notebook: `notebook_publication/04_figure4_apa_mirna_coordination.ipynb`
- Source tables and QC: `data/publication/fig04_apa_mirna_coordination/`
- Figures: `data/figures/publication/fig04_apa_mirna_coordination/`, divided into `figure4/` and `supplementary10/`–`supplementary13/`
- Validation logs: `data/logs/publication/fig04_apa_mirna_coordination/`

The notebook exports **39 PNG/PDF pairs**: 34 individual panels and five assemblies. All are displayed inline. PNGs omit panel letters; PDFs include them. Exported dimensions retain exact multiples of 52.5 mm, without tight-bounding-box resizing. Global font roles are 7–10 pt.

## Figure 4 assembly and panel conventions

Figure 4 uses a 4 × 6 grid of 52.5 mm modules (210 × 315 mm). Its rows are A–D (each 1 × 1), E and F (each 1 × 1) with G (2 × 1), H (4 × 2, retaining its current height), I–L (each 1 × 1), and M–P (each 1 × 1).

Panels B and C annotate their comparisons with raw P values; the corresponding source tables retain BH q values. Panels E and F use the same APA-machinery category colors, with one shared category legend beneath E/F in the assembly; the individual E export carries the key and F does not repeat it. Panel G uses a stronger blue for the normal-only group so it remains distinguishable at its reduced assembled size.

The panel subjects follow [the provenance map](figure4_s10_s13_provenance.md). The CXCL12/miR-135b example is prespecified, and the notebook reports whether its recomputed rhythmicity, anticorrelation, phase relationship, and predicted seed match satisfy the configured criteria. A predicted seed match does not establish experimental regulation.

## Changes in the threshold revision

- Restore nominal P < 0.05 for APA/miRNA, replacing the earlier universal-q revision.
- Replace newly refitted log2(TPM+1) genes with the exact prior-figure TMM gene memberships, phases and significance; retain their additional stored-amplitude gate.
- Join APA hosts and predicted miRNA targets by ENSG; resolve missing TargetScan ENSG values only by unambiguous exact authoritative gene-symbol matching. Export mapping and regression-reproduction QC.
- Retain the existing assay QC, corrected GO enrichment, median proximal PAS, layout and export conventions. Counts may therefore differ slightly from historical screenshots even with nominal-P selection.
- Keep the CXCL12/miR-135b example prespecified and report each eligibility criterion. Its predicted seed span is `chr10:44377322–44377327`, 105 nt distal to the median proximal PAS.
- Complete S10A–D with observed TMEM183A coverage, its BAM/index provenance manifest, and a GENCODE v49 gene-model source table. The cache exists only to distribute this BAM-derived evidence reproducibly when raw BAM access is unavailable.

Results and validation are refreshed after execution in the source tables and `validation_summary.json`.

## Validated threshold-revision results

| Cohort | APA events: P < 0.05 | miRNAs: P < 0.05 | Genes: q < 0.05 and peak-to-trough amplitude > 0.1 |
| --- | ---: | ---: | ---: |
| TCGA | 11,210 | 110 | 12,905 |
| CPTAC | 17,392 | 70 | 6,750 |
| GTEx | 6,605 | 134 | 10,530 |

All three gene memberships match the previous significant-gene lists exactly. The original TMM/non-GEO inputs reproduce stored gene phases and amplitudes to numerical precision. The prespecified CXCL12 triad passes all configured criteria: APA P = 0.0133, miRNA P = 7.81 × 10⁻¹⁸, gene q = 1.14 × 10⁻⁹, and miRNA–gene rho = −0.619 (correlation q = 8.48 × 10⁻¹⁸).

The target reference lacks ENSG values; 214,878 reference rows resolve by exact, unambiguous authoritative symbols and 18,702 remain unmapped. Those exclusions are exported explicitly. Sample overlap for target correlations is TCGA 178, CPTAC 146 and GTEx 346. Batched correlations were checked against 96 direct SciPy calculations, and all 51,663 adjusted correlation q-values were checked against statsmodels BH correction. Validation records are saved under the analysis log directory.

## Panel H selection and layout

Use the original whole-GOA annotated-protein background (19,716 UniProt accessions in the supplied resource), shared across cohorts. Foregrounds remain the current mapped rhythmic APA host genes. Keep the corrected inclusive hypergeometric upper tail P(X ≥ observed overlap), evaluated as sf(k − 1), and exclude NOT-qualified BP memberships. Select raw P < 0.01 and retain the top 20 terms per cohort, ranked by P then GO ID; retain BH q values in the source tables.

Order the common term axis globally by mean −log10(P) across the selected-cohort matrix, with absent entries scored as zero. Stronger terms run left-to-right as in the reference image. Cap display scores at 22, without changing exported P values. Restore the gray background, white grid, red dots with black edges and alpha 0.8, and proportional dot areas (uniformly scaled to fit the publication width). Wrap only long labels while retaining the global minimum font size of 7 pt. The panel remains 210 × 105 mm. Source tables include `apa_go_whole_proteome_background.tsv`, `Figure4H_go_selected_terms.tsv`, and `Figure4H_go_plot_data.tsv`, with explicit selection, rank, mean-score, and plotting-order columns.

## Panel I: tumor-union APA host genes

Compare TCGA ∪ CPTAC with GTEx. A splicing-machinery host gene is rhythmic if any tested APA event has P < 0.05. Form the union of tested tumor host genes for the denominator and the union of rhythmic tumor host genes for the numerator, counting each gene once. Use the analogous GTEx sets for normal percentages. Annotate raw two-sided Fisher P values with * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001; omit nonsignificant labels/brackets. Save counts, percentages, P values, and annotations in `Figure4I_union_gtex_splicing_host_apa_source.tsv`. Other supplementary comparisons retain their existing methods.


S10A–C are adjacent in the complete assembly and use one shared row-z-score color bar (−2 to +2). Individual heatmap exports retain their own color bar for standalone interpretation. Panel dimensions and underlying values are unchanged.


Statistical typography: all numeric significance annotations and scales in Figure 4 and S10–S13 use italic lowercase *p* for raw values and *q* for FDR-adjusted values. This is a formatting change; statistics, thresholds, and asterisk annotations are unchanged.

Figure S12C omits nonsignificant labels; any asterisk denotes BH-adjusted q < 0.05 from the existing two-sided Fisher exact comparisons. Statistical values and thresholds are unchanged.
