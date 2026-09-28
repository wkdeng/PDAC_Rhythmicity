# Figure 2 and Supplementary Figures 4–6: source map

Audited 2026-09-21 against the four supplied figure composites, notebook source and persisted outputs, existing PDFs, and saved analysis tables. This is a provenance inventory; no analyses or figures were regenerated.

Cell numbers below are **one-based notebook cell positions, including Markdown cells**, not execution counters or embedded `# Cell` labels. Panel letters refer to the supplied composites. S6 is described by left/middle/right because the screenshot has no visible panel letters.

## Source notebooks

- **GEO:** [2.1.3_CHIRAL_dual_penalty_pipeline_geo_tpm.ipynb](../2.1.3_CHIRAL_dual_penalty_pipeline_geo_tpm.ipynb)
- **Sensitivity:** [2.1.3_gtex_chiral_sensitivity_supplement.ipynb](../2.1.3_gtex_chiral_sensitivity_supplement.ipynb)
- **Enrichment:** [2.1.5_ChIP_Atlas_clock_TF_enrichment.ipynb](../2.1.5_ChIP_Atlas_clock_TF_enrichment.ipynb)

All figure paths below are relative to `data/figures/` unless otherwise specified. PDF page numbers are one-based.

## Figure 2

| Panel | Source | Exact plot/output |
|---|---|---|
| A | GEO cell **14** | GTEx `plot_one_cohort(..., genes = inline_display, ncol = 3)`; `inline_display = c("BMAL1", "CRY2", "PER3")`. The three-gene layout is an inline output; the related all-gene PDF is `clock_gene_expression_vs_phase.pdf`, page **3**. |
| B | GEO cell **20** | GTEx iteration of the phase–amplitude `ggplot` loop, using `plot_data_all` / `regression_metrics`. Inline only; no explicit save call. |
| C | GEO cell **14** | TCGA iteration of the same three-gene inline plot as A; related all-gene PDF page **1**. |
| D | GEO cell **20** | TCGA iteration of the phase–amplitude loop. Inline only. |
| E | GEO cells **35–36**, `plot_heatmap` | `heatmaps/rhythmic_genes_heatmap_GTEx_Normal.pdf`, page **1**. Related maintained implementation: `scripts/plot_rhythmic_genes_heatmap.R`. |
| F | GEO cells **35–36**, `plot_heatmap` | `heatmaps/rhythmic_genes_heatmap_TCGA-PAAD_Tumor.pdf`, page **1**: TCGA expression and TCGA rhythmic genes, sorted by TCGA acrophase. |
| G | GEO cell **36**, second `plot_heatmap` call | Same TCGA PDF, page **2**: TCGA expression and phase bins, using the GTEx rhythmic-gene list `normal_gene_list`. See the ordering caveat below. |
| H | GEO cell **38** invokes `scripts/plot_rhythmic_genes_venn.R`, `make_overlap_barplot` | `rhythmic_genes_overlap_barplot_TCGA_PAAD_Tumor_vs_GTEx.pdf`. |
| I | GEO cell **38** invokes `scripts/run_acrophase_plot.R`, `calculate_circular_density` | `acrophase_polar_distribution_all_datasets.pdf`, page **1**. |
| J | Same script, `R2_stat` / `bootstrap_R2` | Same PDF, page **2**: axial second-harmonic resultant length with bootstrap intervals. |
| K | Enrichment cell **16** (embedded `# Cell 11`) | `chip_atlas_tf_enrichment_heatmap.pdf` / `.png`; `fe_pivot` supplies annotations, `neglog10_p` supplies color. Plotting source identified, but the screenshot's numbers differ from current saved results. |
| L | Enrichment cell **23** (embedded `# Cell 18`) | `ppi_enrichment_barplot_gtex_tcga.pdf` / `.png`; focused `df_two` plot from `ppi_combined_df`. |

Matching numerical evidence:

- A/C: GTEx **362** and TCGA **178** samples; BMAL1 amplitude/p **0.79 / 1.6e-48** and **0.43 / 4.9e-14**, respectively, match the existing scatter PDF and screenshot. BMAL1 is the ARNTL display alias. Despite the notebook filename, the plotted expression is **log2 CPM with TMM normalization**.
- E/F/G: saved heatmap PDF titles report **10,254** GTEx and **12,560** TCGA genes. Both TCGA PDF pages repeat the same TCGA title; page 2's title does not establish its actual gene-set size or row ordering.
- H: current significant-gene tables reproduce GTEx **10,530** and TCGA **12,905**, with **5,899** common and **7,006** TCGA-only genes.
- I/J: axial resultant lengths from the current tables are GTEx **0.136388** and TCGA **0.662608**. This is `abs(mean(exp(2i * theta)))`, not regression R-squared.
- L: GTEx enrichment **1.3318588**, raw p **1.75655e-9**; TCGA **1.0551489**, raw p **0.2985874**, matching the screenshot.

## Supplementary Figure 4

| Panel | Source | Exact output |
|---|---|---|
| A | GEO cell **14**, `plot_one_cohort("TCGA-PAAD_Tumor")` | `clock_gene_expression_vs_phase.pdf`, page **1**, all 12 genes. |
| B | GEO cell **14**, `plot_one_cohort("GTEx_Normal")` | Same PDF, page **3**, all 12 genes. |

The intervening PDF page 2 is CPTAC and is not part of the supplied S4 composite. The 12 displayed genes are BMAL1, CLOCK, NPAS2, CRY1, CRY2, PER1, PER2, PER3, NR1D1, NR1D2, DBP, and TEF.

## Supplementary Figure 5

| Panel | Source | Exact output |
|---|---|---|
| A–C | GEO cell **18**, per-cohort `svd(complex_rep)` loop | Select GTEx, TCGA, and CPTAC, respectively, from inline `Genes SVD ... mode 1` plots. No explicit file save. Sample counts: 362, 178, 161. |
| D | GEO cell **20** | CPTAC iteration of the phase–amplitude polar loop; inline only. |
| E | GEO cell **38** → `scripts/plot_rhythmic_genes_venn.R` | `rhythmic_genes_overlap_barplot_CPTAC_3_Tumor_vs_GTEx.pdf`. CPTAC **6,750 = 3,005 common + 3,745 unique**. |
| F | Same script, `venn.diagram` | `rhythmic_genes_venn.pdf`; three-cohort Venn. |
| G | GEO cell **38** → `scripts/run_acrophase_plot.R`, `perm_test_delta_R2` | `acrophase_polar_distribution_all_datasets.pdf`, page **3**. The fourth page is a separate global test and is not in this composite. |
| H | Sensitivity cell **10** (embedded `# Cell 9`), `svd_plot` | `gtex_chiral_sensitivity_supplement/core_clock_harmonic_svd_mode_1.pdf` / `.png`. |
| I | Sensitivity cell **9** (embedded `# Cell 8`), `density_plot` | `gtex_chiral_sensitivity_supplement/rhythmic_gene_phase_circular_density.pdf` / `.png`. |
| J | Enrichment cell **23**, `pivot_fold`, left subplot | `ppi_enrichment_heatmap.pdf`, left half. |
| K | Enrichment cell **23**, `pivot_logp`, right subplot | Same PDF, right half. |

Matching numerical evidence:

- F Venn regions: CPTAC-only **1,821**; GTEx-only **3,490**; TCGA-only **5,082**; CPTAC–GTEx-only **1,141**; CPTAC–TCGA-only **1,924**; GTEx–TCGA-only **4,035**; all three **1,864**.
- G: observed TCGA-minus-GTEx axial difference **0.526220**, **20,000** permutations, displayed p **5e-5**; the plus-one minimum is `1/20001`. Script seed for this comparison is **2001**.
- H/I: **172** GTEx donors selected by time of death in **[0,12) h**; **6,367** rhythmic genes; axial resultant length **0.123914**; axis **4.13637 / 16.13637 h**. Persisted inline plots were visually matched to the screenshot.
- J/K: current per-gene PPI tables reproduce the displayed matrices. Examples: CPTAC PER3 fold enrichment **3.01**; GTEx PER2 **1.60**, `-log10(p)` **6.6**; TCGA PER3 **2.05**. Displayed columns are BMAL1, PER1, PER2, PER3, CRY1, CRY2, NR1D1, NR1D2, CIART, NPAS2.

## Supplementary Figure 6

All three sections are in the GEO notebook and compare each tumor cohort's rhythmic genes that are common with GTEx against those absent from the GTEx rhythmic list.

| Position | Source | Exact output / numerical match |
|---|---|---|
| Left | Cell **40**, `cohens_d`, `p_pval` | `pvalue_distribution_common_vs_unique.pdf`. Cohen's d on `-log10(p)`: CPTAC **0.119**, TCGA **0.105**. |
| Middle | Cell **42**, `calculate_circular_density`, `circular_w1_grid`, `p_acrophase` | `acrophase_distribution_common_vs_unique.pdf`. Circular Wasserstein-1: CPTAC **0.2735869 h**, TCGA **0.2484593 h**. |
| Right | Cell **44**, `ks_stats`, `p_p2t` | `peak_to_trough_common_vs_unique.pdf`. CPTAC K–S **D=0.04252947**, **p=0.004805298**; TCGA **D=0.09604612**, **p=4.370853e-26**. |

The right plot computes `peak_to_trough = 2 * Amplitude + 1` and plots `log2(peak_to_trough)`; the exact numerical transformation is therefore `log2(2 * Amplitude + 1)`. Its displayed legacy axis wording is `log2(Peak-to-Trough + 1)`.

## Input chains and source tables

### Main expression, phase, and rhythmicity

- `data/Merged/GTEx_TCGA_GEO_merged_TMM_unfiltered_expression.csv`
- `data/Merged/GTEx_TCGA_GEO_merged_meta.csv`
- `data/Merged/CHIRAL_dual_penalty_phi_12h_geo.csv`
- `data/Merged/regression_periodic_12h_geo.txt`
- `data/Merged/regression/sig_genes_regression_periodic_12h_{GTEx_Normal,TCGA-PAAD_Tumor,CPTAC-3_Tumor}.txt`
- Gene annotation maps `n2g` / `g2n`; the notebook's filtered in-memory `tmm_grouped` supplies the main SVD.

### Restricted GTEx sensitivity

Primary inputs are `data/GTEx/gene_reads_v10_pancreas.gct`, `data/GTEx/MetaData/TOD.txt`, and `data/Ref/t2g.txt`; inference currently sources `scripts/CHIRAL_dual_penalty.R`.

Reusable tables under `data/gtex_chiral_sensitivity_supplement/` include `tod_window_sample_metadata.tsv`, `input_qc.tsv`, `chiral_phase_estimates.tsv`, `chiral_window_qc.tsv`, `harmonic_regression_all_genes.tsv`, `rhythmic_genes.tsv`, `core_clock_svd_gene_loadings.tsv`, `core_clock_svd_variance_explained.tsv`, and `rhythmic_gene_phase_circular_{density,statistics}.tsv`.

### ChIP and PPI enrichment

- ChIP heatmap: `data/enrichment/chip_atlas_enrichment/enrichment_summary.txt`; produced in Enrichment cell **15**, based on cell **13** results. Upstream: `peak_gene_matrix.txt` in that directory, `data/ChIP_Atlas/Oth.ALL.05.<TF>.AllCell.bed`, and `data/Ref/gencode.v49.basic.annotation.gtf`.
- PPI tables: `data/enrichment/ppi_enrichment/ppi_enrichment_combined.txt` and `ppi_enrichment_per_clock_gene.txt`; generated through Enrichment cells **19–21** from physical interactions in `data/PPI/BIOGRID-ORGANISM-Homo_sapiens-5.0.253.tab3.txt`.
- Enrichment gene-set loading is in cell **9**: GTEx significant genes plus tumor `data/enrichment/cancer_specific_genes/method_A_presence/{TCGA-PAAD_Tumor,CPTAC-3_Tumor}_cancer_specific_genes.txt`. Thus the tumor enrichment targets are cancer-specific gene sets; plot titles alone do not describe that selection.
- The plotted ChIP colors/stars and PPI heatmap significance colors use **raw p values**, although q values are also stored in result tables.

## Differences to resolve when building the publication notebook

1. **Heatmap revision mismatch.** Existing heatmap PDFs and embedded outputs report GTEx 10,254 / TCGA 12,560, whereas current significant-gene tables report 10,530 / 12,905. The heatmap PDFs predate the current regression and expression/phase files. The evidence supports different saved revisions; it does not establish a current filtering step that explains the difference. The old notebook accesses `$Gene` as ENSG identifiers, while current tables separate `Gene` symbols and `ENSG`; the maintained heatmap script uses `$ENSG`.
2. **Figure 2G row order and title.** `normal_gene_list <- sig_genes_in_data`, not `sig_genes_ordered`. The second TCGA heatmap therefore uses GTEx gene membership in source-file order, with TCGA samples/bins. Its inherited title says “sorted by acrophase” and repeats the TCGA count; these labels cannot be used as evidence of GTEx acrophase ordering or the actual row count.
3. **Figure 2K historical values.** The plotting cell and saved output family are identified, but the screenshot's exact data revision is not. Screenshot GTEx BMAL1/CRY1 values are **1.02 / 1.09** versus current table **1.00 / 1.02**; screenshot TCGA BMAL1/DBP/RORC values are **0.93 / 0.63 / 1.15** versus current **0.91 / 0.86 / 1.18**. These differences are real, not a rounding-only discrepancy. Current ChIP tables and PDFs exist; the historical table matching the screenshot has not been located.
4. **Gene lists differ deliberately in the sources.** Main SVD uses 12 genes including CIART but excluding CLOCK; expression scatter uses 12 including CLOCK but excluding CIART; amplitude polars include both (13); sensitivity SVD uses the scatter-style 12-gene list. These should not be silently unified during migration.
5. **Restricted sensitivity is donor selection, not a hard bound on final phase estimates.** The saved QC has 172 donors, inferred phases **0.470853–23.966848 h**, and **119/172** within 0–12 h. Preserve the distinction between donor time-of-death selection and the inference-window setting. Single-cohort SVD uses one complex coefficient column and is rank one by construction in both the main and sensitivity implementations.

Persistent publication rules remain in [AGENTS.md](AGENTS.md). This audit does not change the Figure 1/S1–S3 notebook or exports and does not create a Figure 2 analysis notebook yet.

## Subsequent implementation decision

The publication notebook `02_figure2_cohort_clock_organization.ipynb` uses current heatmap inputs, explicitly sorts Figure 2G by GTEx acrophase, and uses current ChIP results. Reconstruction identifies 10,257 current GTEx and 12,576 current TCGA significant ENSG rows in the expression matrix; its unmatched-identifier tables document the exclusions. The older PDF counts above remain evidence of the supplied source revision.

## Figure 2K protocol revision (2026-09-22)

At the user's request, publication panel K now recomputes the protein-coding, shared-tested-background analysis from corrected promoter peak counts. It compares the GTEx rhythmic reference with the saved TCGA `method_A_presence` cancer-specific target, includes a tenth pooled CCG endpoint, and uses BH q-values across all 20 tests for colors and stars. The background is 16,570 genes and target sizes are 8,076 GTEx / 3,394 TCGA-specific genes. The notebook embeds the computation and exports all source counts, raw p-values and q-values. The source-map and screenshot observations above remain historical provenance; panel L and S5J–K retain their prior PPI analysis.
