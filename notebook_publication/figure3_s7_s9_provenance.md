# Figure 3 and Supplementary Figures S7–S9: source map

Implementation follow-up: the publication notebook is now available as
[03_figure3_rhythmic_alternative_splicing.ipynb](03_figure3_rhythmic_alternative_splicing.ipynb).
See [its execution guide](03_figure3_readme.md) for inputs, validation, and the
real BAM-derived Figure 3E now included in the assembly. The audit below describes
the original provenance investigation, before this implementation.

Audited against the four supplied composites, legacy notebook code and stored outputs, saved figures, and analysis tables. This document records provenance only; no analysis or figures were regenerated. The screenshots are interpreted in the supplied order as Figure 3, S7, S8, and S9. Comments embedded in the screenshots are reference material, not adopted revision instructions.

Cell numbers are **one-based notebook positions, including Markdown cells**, not execution counts or embedded `# Cell` labels. Figure paths below are relative to `data/figures/`.

## Main source notebooks

- **GEO:** [2.1.3_CHIRAL_dual_penalty_pipeline_geo_tpm.ipynb](../2.1.3_CHIRAL_dual_penalty_pipeline_geo_tpm.ipynb).
- **AS:** [3.1.2_alternative_splicing_rhythmic_analysis.ipynb](../3.1.2_alternative_splicing_rhythmic_analysis.ipynb).

## Pathway enrichment panels

| Panel | Source | Saved output, available as PDF and PNG |
|---|---|---|
| Figure 3A | GEO cell **59**, GO-only heatmap block | `pathway_grouping_12h_geo/domain_crosscheck_heatmap_GO_only` |
| S7A | GEO cell **59**, selected tumor-pathway dotplot block | `pathway_grouping_12h_geo/highlighted_tumor_pathway_dotplot` |

GEO cells **53, 55, 57** load enrichment results, cluster pathway text, and assign highlighted domains. The actual input files recorded by cell 53 are:

- `data/enrichment/full_term_GO_Regression_CHIRAL_12h_geo.txt`
- `data/enrichment/full_term_Reactome_Regression_CHIRAL_12h_plus_geo.txt`
- `data/enrichment/full_term_KEGG_Regression_CHIRAL_12h_plus_geo.txt`

The GO-only summary uses nominal **p ≤ 0.05**; cell color is maximum −log10(p), and numbers count pathways. The current saved table exactly matches Figure 3A, in GTEx/CPTAC/TCGA order: splicing **9/13/5**, APA **1/8/2**, miRNA **4/5/3**, transcription **37/41/30**. DNA/chromatin methylation is explicitly omitted from this GO-only panel.

S7A uses **−log10(p) ≥ 2**, selects up to four pathways per domain/library/tumor cohort, then includes qualifying GTEx points for those selected pathways. Dot color represents −log10(p), size represents enrichment ratio, and shape represents library. Domain membership is a regex overlay after TF-IDF/LSA/K-means clustering; the domains do not define the unsupervised clusters.

Regenerating tables are under `data/enrichment/pathway_grouping_12h_geo/`, particularly `domain_crosscheck_summary.tsv`, `highlighted_domain_pathways_primary_threshold.tsv`, `pathway_grouping_input_all_sources.tsv`, and `pathway_text_clusters.tsv`. Upstream GO enrichment is called in GEO cell **23** through `hypergeometry.go_analysis`; KEGG/Reactome table concatenation is in cells **28/32**. The latter tables combine `_geo` and earlier non-`_geo` enrichment results, so their input versions need explicit treatment in a new notebook.

## Splicing machinery panels

| Panel | Source | Saved output, available as PDF and PNG |
|---|---|---|
| Figure 3G | AS cell **33**, TCGA subplot | `AS_splicing_machinery/splicing_machinery_polar_acrophase` |
| Figure 3H | AS cell **33**, GTEx subplot | Same file, GTEx subplot |
| Figure 3I | AS cell **34**, right subplot | `AS_splicing_machinery/splicing_machinery_consistency_tumor_normal` |
| S8A | AS cell **33**, CPTAC subplot | `AS_splicing_machinery/splicing_machinery_polar_acrophase` |
| S8B | AS cell **32**, significance and acrophase heatmaps | `AS_splicing_machinery/splicing_machinery_rhythmicity_heatmap` |

AS cell **31** defines 81 genes across four categories and reads `data/Merged/regression/all_genes_regression_periodic_12h_{cohort}.txt` plus `data/Ref/t2g.txt`. Current tables reproduce **44 CPTAC, 54 TCGA, and 50 GTEx** significant genes using **raw p < 0.05**, matching the polar-plot titles. The radius is regression amplitude. S8B color is −log10(p) × sign(amplitude), clipped to [−5,5]; its phase heatmap masks nonsignificant genes.

Figure 3I counts a gene as tumor rhythmic when significant in **at least one enabled tumor cohort**, rather than using TCGA alone. The current notebook enables CPTAC and TCGA tumors. Its saved plot and stored table give tumor/normal counts **30/20 of 34**, **9/9 of 14**, **11/9 of 13**, and **16/12 of 20** for core spliceosome, SR proteins, hnRNP proteins, and other regulators. Those bar heights differ from the supplied Figure 3I screenshot. The plotting source is identified, but the exact older data/configuration producing the screenshot's category bars has not been established.

## SE events, SLK example, and expression comparisons

| Panel | AS notebook cell | Existing output under `data/figures/` |
|---|---|---|
| Figure 3B | **17** defines density; **24/25** call it | `AS_rhythmicity/GTEx/GTEx_Normal/SE/phase_distribution.pdf/png` |
| Figure 3C | Same | `AS_rhythmicity/TCGA/TCGA_Tumor/SE/phase_distribution.pdf/png` |
| Figure 3D | **29** | `AS_rhythmicity/TCGA/TCGA_Tumor/SE/SLK_demo/scatter_SLK_representatives.pdf/png` |
| Figure 3E | **29** selects samples and writes real-sashimi settings/commands | Real coverage export in screenshot not located. Related local export is `AS_rhythmicity/TCGA/TCGA_Tumor/SE/SLK_demo/sashimi_pseudo_SLK.pdf/png`, which is a different visualization. |
| Figure 3F | **16** defines heatmap; **24/25** call it | `AS_rhythmicity/TCGA/TCGA_Tumor/SE/rhythmic_heatmap.pdf/png` |
| S7B | Same | `AS_rhythmicity/CPTAC/CPTAC_Tumor/SE/rhythmic_heatmap.pdf/png` |
| S7C | Same | `AS_rhythmicity/GTEx/GTEx_Normal/SE/rhythmic_heatmap.pdf/png` |
| Figure 3J/K | **37**, GTEx/TCGA subplots | `AS_rhythmicity/phase_distribution_comparison/overall_SE_as_vs_gene_phase_density.pdf/png` |
| Figure 3L | **59**, raw-amplitude subplot | `AS_rhythmicity/se_rhythm_stratified/gtex_tcga_panel_a.pdf/png` |
| Figure 3M | **59**, effect-size subplot | `AS_rhythmicity/se_rhythm_stratified/gtex_tcga_panel_c.pdf/png` |
| Figure 3N | **79** | `AS_clock_tf_enrichment/clock_tf_peak_fraction_SE_bar.pdf/png` |

The saved SE significant tables under `data/AS/{project}/rhythmic_results/{dataset}/SE/significant_pvalue_0.05.csv` contain **162 CPTAC, 387 TCGA, and 242 GTEx** events. The phase-density counts match Figure 3B/C. The saved heatmap function z-scores each event's PSI, bins samples into 24 phase bins, and sorts events by acrophase. S7B/C are rotated crops of the CPTAC/GTEx heatmaps; the CPTAC missing-bin bands distinguish its panel.

SLK is **SE event 150506**. `data/AS/TCGA/rhythmic_results/TCGA_Tumor/SE/SLK_demo/` contains `representative_samples.csv`, `SE_SLK.MATS.JC.txt`, and `sashimi_plot_settings.txt`. The plot-generation cell reports fitted acrophase **4.58 h** and selects four representatives near 4, 10, 16, and 22 h.

Figure 3J/K compare the overall marginal distributions, **242 SE events versus 10,530 rhythmic expression genes** in GTEx and **387 versus 12,905** in TCGA. Their saved statistics are `data/AS/phase_distribution_comparison/overall_SE_as_vs_gene_phase_stats.tsv`. They are not restricted to an event's own rhythmic host gene.

Figure 3L/M source and bootstrap results are `data/AS/se_rhythm_stratified/gtex_tcga_peak_to_trough_{source,effects}.tsv`. The full original three-panel figure is `AS_rhythmicity/se_rhythm_stratified/gtex_tcga_peak_to_trough_comparison.pdf/png`; its middle panel is reused as S9E.

## Supplementary Figure S9

| Panel | AS notebook cell | Existing output under `data/figures/` |
|---|---|---|
| A–D | **66**, four subplots; circular helpers in **64** | `AS_rhythmicity/phase_distribution_comparison/as_se_vs_gene_phase_metrics.pdf/png` |
| E | **59**, GTEx/TCGA normalized accumulation | `AS_rhythmicity/se_rhythm_stratified/gtex_tcga_panel_b.pdf/png` |
| F | **59**, GTEx/CPTAC normalized accumulation | `AS_rhythmicity/se_rhythm_stratified/gtex_cptac_panel_b.pdf/png` |
| G | **59**, GTEx/CPTAC effect size | `AS_rhythmicity/se_rhythm_stratified/gtex_cptac_panel_c.pdf/png` |
| H | **78** | `AS_clock_tf_enrichment/clock_tf_peak_fraction_bar.pdf/png` |
| I–K | **83**, bottom polar subplots in CPTAC/TCGA/GTEx order; statistics in **81/82** | `AS_clock_tf_enrichment/phase_diff_by_clock_peaks.pdf` |

S9A–D use SE events matched to rhythmic expression host genes, constructed through `all_phase_data` in cell **62**. Their event/unique-host counts are **41/36 CPTAC, 236/189 TCGA, and 108/91 GTEx**. `data/AS/phase_distribution_comparison/as_se_vs_gene_phase_tests.tsv` reproduces the screenshot's permutation p values **0.289, 0.0001, 0.608** and Watson–Wheeler p values **0.0368, 0.0473, 0.0369**. Additional bootstrap summaries are in `as_se_vs_gene_phase_metrics.tsv` in the same directory.

S9F/G sources are `data/AS/se_rhythm_stratified/gtex_cptac_peak_to_trough_{source,effects}.tsv`. S9I–K stratify the paired event–host phase differences by any clock-TF promoter peak; they use all included AS event types, rather than SE alone. The comparison table is `data/AS/phase_diff_clock_peak_comparison.csv`.

Figure 3N and S9H use `data/enrichment/chip_atlas_enrichment/peak_gene_matrix.txt`, loaded in cell **76**. Their saved counts match the screenshots:

- Figure 3N, SE only: **15/41 CPTAC, 73/236 TCGA, 60/108 GTEx**.
- S9H, all included AS types: **163/483 CPTAC, 544/1685 TCGA, 701/1084 GTEx**.

The fraction tables are `AS_clock_tf_enrichment/clock_tf_peak_fraction_SE.tsv` and `clock_tf_peak_fraction.tsv`. Despite the titles referring to host genes, these bars count **event rows**, so a host with multiple events is counted more than once. Both use the expression-matched `all_phase_data` population. Preserve this distinction when deciding the publication analysis unit.

## Processed-input provenance and migration boundaries

The legacy AS analysis starts with `data/AS/{CPTAC,TCGA,GTEx}/post/*.MATS.JC.txt`, positional sample lists in `data/AS/{project}/all_samples.txt`, `data/Merged/CHIRAL_dual_penalty_phi_12h_geo.csv`, and `data/Merged/GTEx_TCGA_GEO_merged_meta.csv`. Gene-expression comparisons also use the all-gene/significant-gene regression tables in `data/Merged/regression/`. SLK junction details use `individualCounts.SE.txt`.

The checked-in root `process.smk` sets `BATCH='geo_batch1'`. It declares kallisto abundance files, STAR sorted BAM/BAI and junction outputs, and `data/geo_batch1/featurecounts/gene_counts.tsv`. Its rMATS post rule declares a completion marker, while the JC and individual-count tables are tool-generated side effects. It does **not** directly declare the TCGA/CPTAC/GTEx `data/AS/` paths consumed by these figures. Those cohort artifacts therefore need their own explicit provenance rather than being described as verified direct products of this exact workflow.

The rMATS sample order is essential: legacy code indexes comma-separated PSI/read-count vectors using the BAM-list order. A public notebook should validate vector lengths and sample identities before analysis. The legacy AS significance configuration is nominal **p < 0.05**, despite also calculating BH q values; migration must not silently substitute a different threshold.

True SLK read-coverage rendering requires four BAM/BAI pairs referenced under `data/bam/tcga_genomic/`; they are absent locally. Saved plots/settings exist, and a count-based pseudo-sashimi implementation exists, but that implementation is not equivalent to the screenshot's real coverage tracks. External annotation inputs include the ChIP-Atlas BED files, their derived peak-gene matrix, and GENCODE v49 annotation; these are not produced by `process.smk`.

An older generator, [scripts/build_fig03_publication_notebooks.py](../scripts/build_fig03_publication_notebooks.py), targets `notebooks_publish/` and numbers related supplements **S5–S7**, corresponding by topic to the supplied **S7–S9**. It is a proposed reimplementation, not the primary provenance for these screenshots. It recomputes CHIRAL from TPM and uses different calling parameters; it must not be run unchanged to reproduce the legacy results. No current Figure 3 publication notebook has been created by this audit.

## Publication conventions

The persistent rules in [AGENTS.md](AGENTS.md) apply to subsequent implementation: self-contained notebooks in `notebook_publication/`; processed-data inputs with documented external resources; individual and assembled figures displayed inline; outputs under `data/figures/publication/`; names `Figure{number}{panel}_{content_related}.png/pdf`; 52.5 mm square layout modules; minimal padding; globally configured fonts of at least 7 pt; consistent styling; panel letters in PDFs only. Existing notebooks, data, and figures must be preserved.

## Figure 3I SR-count audit (2026-09-22)

The current category comparison uses nominal expression-regression P < 0.05.
The tumor bar is the union of genes passing in TCGA or CPTAC, not TCGA alone.
The saved executed output in `3.1.2_alternative_splicing_rhythmic_analysis.ipynb`,
cell index 33 (one-based cell 34, embedded heading “Cell 25: Cross-dataset
consistency”), already reports equal SR fractions: 9/14 for the tumor union and
9/14 for GTEx. The current publication `machinery.tsv` reproduces these counts.

| Category | TCGA | CPTAC | TCGA/CPTAC union | GTEx |
|---|---:|---:|---:|---:|
| SR proteins | 6/14 | 8/14 | 9/14 (64.3%) | 9/14 (64.3%) |
| Core spliceosome | 25/34 | 19/34 | 30/34 (88.2%) | 20/34 (58.8%) |

SR membership in the tumor union is SRSF1, SRSF2, SRSF3, SRSF6, SRSF7, SRSF10,
SRSF11, SRSF12, and TRA2A. GTEx membership is SRSF1, SRSF3, SRSF6, SRSF8, SRSF9,
SRSF10, SRSF11, SRSF12, and TRA2B. Equal counts therefore do not imply identical
memberships. The historical screenshot appears to show about 86% (12/14) for
tumor, but no inspected saved configuration or source table reproduced that bar.
Its exact cause remains unresolved; it is not an effect of the current layout
or label revision.

The notebook exports `Figure3I_machinery_category_counts.tsv` with all cohort
counts and two-sided Fisher P values, and `Figure3I_machinery_gene_membership.tsv`
with per-gene values and rhythmic calls. I labels the tumor union “TCGA/CPTAC”
and displays only a Core asterisk when its Fisher P < 0.05. G–I use the requested
“APA machinery” titles while retaining the original splicing-machinery gene set.


## Real panel E added from supplied BAMs (2026-09-24)

The user supplied the four original genomic BAMs and indexes on GRC under
`data/figure3E_bams/`. Docker samtools 1.21 quickcheck passes all four; read groups
match the selected TCGA samples and chr10 length is 133,797,422 (GRCh38).
Using primary, QC-pass, properly paired, NH=1, unclipped 48-nt reads reproduces
all 12 saved rMATS UTC/TDC/UDC counts exactly. The new individual panel plots
actual unsmoothed per-base read depth and those directly observed junctions.
The legacy phase and PSI labels are retained. Coverage is raw read depth with
independent y scales; it does not claim to reproduce the unknown normalization
of the historical screenshot. Full BAMs remain on GRC; the extracted source
bundle includes input metadata, header/index hashes, QC and table hashes.
The subsequent layout revision inserts this real panel into Figure 3, moves the
first-to-third exon skipping arc above coverage, and exposes hypothetical BAM
paths in the first notebook code cell.
