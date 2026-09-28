# Figure 4 and Supplementary Figures S10–S13: source map

Audited against the five supplied composites, legacy notebook code, stored outputs, saved figures, and source tables. The supplied order is interpreted as Figure 4, S10, S11, S12, and S13. This is a provenance audit only: no analyses or figures were regenerated. Content in the screenshots is reference material, not additional instructions.

Cell numbers below are **one-based notebook positions, including Markdown cells**, not execution counts or embedded `# Cell` comments. Supplementary panel letters obscured or absent in the screenshots are assigned by reading order. Figure stems below are relative to `data/figures/` and have both `.pdf` and `.png` versions.

## Primary sources

- **APA:** [3.4.1_apa_multi_cohort.ipynb](../3.4.1_apa_multi_cohort.ipynb).
- **miRNA:** [3.3.1_mirna_tcga_cptac_tumor.ipynb](../3.3.1_mirna_tcga_cptac_tumor.ipynb).
- **Triad:** [3.4.2_apa_mirna_rhythmic_connection.ipynb](../3.4.2_apa_mirna_rhythmic_connection.ipynb).
- **Coverage renderer:** [scripts/extract_and_plot_tmem183a_track.py](../scripts/extract_and_plot_tmem183a_track.py).
- **Heatmap renderer:** [scripts/plot_rhythmic_apa_heatmap.R](../scripts/plot_rhythmic_apa_heatmap.R).
- **APA/expression phase helpers:** [scripts/polyapa_phase_tools.py](../scripts/polyapa_phase_tools.py).

Most screenshot panels are crops of larger saved figures. The current notebooks identify the component sources; this audit did not identify a script assembling these exact five supplied composites.

## Figure 4

| Panel | Content | Source | Saved figure stem / selected subplot |
|---|---|---|---|
| A | Three-cohort APA acrophase density | APA **23**, helper in **13** | `polyAPA_rhythmicity/comparison/phase_circular_overlay` |
| B | Relative APA amplitude | APA **25** | `polyAPA_rhythmicity/comparison/amplitude_violin`, right subplot |
| C | TCGA–GTEx rhythmic APA host-gene overlap | APA **19–21** | `polyAPA_rhythmicity/comparison/overlap_bar_chart`, TCGA/GTEx subplot |
| D | APA machinery gene-expression rhythmicity: pooled tumors vs normal | APA **43, 46** | `apa_machinery_rhythmicity/apa_machinery_consistency_tumor_vs_normal`, top-right subplot |
| E | GTEx APA machinery expression amplitude/acrophase | APA **45** | `apa_machinery_rhythmicity/apa_machinery_polar_acrophase`, GTEx subplot |
| F | TCGA APA machinery expression amplitude/acrophase | APA **45** | Same stem, TCGA subplot |
| G | Rhythmic APA vs rhythmic host-expression phases, TCGA/GTEx | APA **33, 34, 38, 40** | `polyAPA_rhythmicity/comparison/apa_expression_overlap_phase_scatter_tcga_vs_gtex` |
| H | GO enrichment of rhythmic APA host genes across cohorts | APA **27, 29** | `polyAPA_rhythmicity/comparison/go_enrichment_dotplot_compare`; the screenshot appears to rotate this tall source plot into a wide panel |
| I | Splicing-machinery genes carrying rhythmic APA, TCGA/GTEx | APA **31** | `apa_event_heatmap/splicing_machinery_tumor_vs_normal_barplot`, TCGA subplot |
| J | GTEx miRNA vs gene-expression phase density | miRNA **23** | `miRNA_rhythmicity/landscape/gtex_normal/gtex_normal_mirna_vs_gene_phase_density` |
| K | TCGA miRNA vs gene-expression phase density | miRNA **23** | `miRNA_rhythmicity/landscape/tcga_paad_tumor/tcga_paad_tumor_mirna_vs_gene_phase_density` |
| L | Predicted miRNA-target phase alignment across cohorts | miRNA **27, 43** | `miRNA_rhythmicity/landscape/mirna_target_functional_landscape`, target phase-alignment subplot |
| M | CXCL12 expression vs phase | Triad **19–20** | `apa_mirna_rhythmic_connection/selected_example_panel`, top-left |
| N | CXCL12 APA PDUI vs phase | Triad **19–20** | Same stem, top-right |
| O | hsa-mir-135b expression vs phase | Triad **19–20** | Same stem, middle-left |
| P | miRNA, expression, and APA acrophase timeline | Triad **19–20** | Same stem, middle-right |

**D is pooled:** tumor rhythmic means nominal p < 0.05 in **either TCGA or CPTAC**. Its seven categories use curated category genes represented in the expression tables as denominators, retaining cross-listed members for each category; this is not a TCGA-only comparison. **E/F are expression of APA machinery genes**, whereas **I is APA events hosted by splicing-machinery genes**. The latter uses cohort-specific tested host-gene denominators, not all genes in each curated category.

## Supplementary figures

| Panel | Content | Source | Saved figure stem / subplot |
|---|---|---|---|
| S10A | GTEx rhythmic APA PDUI heatmap | APA **49**, R heatmap renderer | `apa_event_heatmap/GTEx_apa_pdui_heatmap` |
| S10B | TCGA rhythmic APA PDUI heatmap | Same | `apa_event_heatmap/TCGA_apa_pdui_heatmap` |
| S10C | CPTAC rhythmic APA PDUI heatmap | Same | `apa_event_heatmap/CPTAC_apa_pdui_heatmap` |
| S10D | TMEM183A four-sample RNA-seq coverage | APA **53** selects samples; coverage renderer draws tracks | `apa_event_heatmap/TMEM183A_apa_track_4samples` |
| S11A | Absolute APA amplitude | APA **25** | `polyAPA_rhythmicity/comparison/amplitude_violin`, left subplot |
| S11B | CPTAC–GTEx rhythmic APA host-gene overlap | APA **19–21** | `polyAPA_rhythmicity/comparison/overlap_bar_chart`, CPTAC/GTEx subplot |
| S11C | CPTAC APA machinery expression amplitude/acrophase | APA **45** | `apa_machinery_rhythmicity/apa_machinery_polar_acrophase`, CPTAC subplot |
| S11D | APA machinery expression rhythmicity heatmap | APA **43–44** | `apa_machinery_rhythmicity/apa_machinery_heatmap`, significance component |
| S11E | Significant machinery-gene acrophases | APA **44** | Same stem, companion phase component |
| S12A | APA minus host-expression phase lag, TCGA/GTEx | APA **38, 40** | `polyAPA_rhythmicity/comparison/apa_expression_phase_lag_polar_density_by_group_tcga_vs_gtex` |
| S12B | APA minus host-expression phase lag, CPTAC/GTEx | APA **39, 41** | `polyAPA_rhythmicity/comparison/apa_expression_phase_lag_polar_density_by_group_cptac_vs_gtex` |
| S12C | Splicing-machinery genes carrying rhythmic APA, CPTAC/GTEx | APA **31** | `apa_event_heatmap/splicing_machinery_tumor_vs_normal_barplot`, CPTAC subplot |
| S13A | TCGA miRNA amplitude vs −log10(p) | miRNA **17, 19** | `miRNA_rhythmicity/landscape/tcga_paad_tumor/tcga_paad_tumor_mirna_rhythmicity_volcano` |
| S13B | GTEx predicted miRNA-target expression correlations | miRNA **37** | `miRNA_rhythmicity/landscape/gtex_normal/gtex_normal_predicted_mirna_target_expression_correlation` |
| S13C | TCGA predicted miRNA-target expression correlations | miRNA **37** | `miRNA_rhythmicity/landscape/tcga_paad_tumor/tcga_paad_tumor_predicted_mirna_target_expression_correlation` |
| S13D | CPTAC predicted miRNA-target expression correlations | miRNA **37** | `miRNA_rhythmicity/landscape/cptac3_tumor/cptac3_tumor_predicted_mirna_target_expression_correlation` |
| S13E | hsa-mir-135b vs CXCL12 expression | Triad **19–20** | `apa_mirna_rhythmic_connection/selected_example_panel`, bottom-left |
| S13F | Predicted seed match in CXCL12 distal aUTR | Triad **19–21** | Same stem, bottom-right |

S12 phase-lag plots are also saved under equivalent `...phase_lag_by_group...`, `...phase_lag_circular_density_by_group...`, and `...phase_lag_second_moment_by_group...` stems. They are aliases from the same plotting loop, not independent analyses.

## Verified screenshot fingerprints

| Feature | Values in screenshot and current saved sources |
|---|---|
| Figure 4C, TCGA-only / shared / GTEx-only genes | **1,475 / 656 / 1,256**, Fisher p = 9.0836e−8 |
| S11B, CPTAC-only / shared / GTEx-only genes | **2,792 / 918 / 994**, Fisher p = 3.2857e−47 |
| S10 heatmap event counts, GTEx / TCGA / CPTAC | **6,605 / 11,244 / 17,339** |
| Machinery polar significant genes, GTEx / TCGA / CPTAC | **16 / 17 / 23**; 30 of the 32 curated genes occur in the expression tables |
| Figure 4G event-expression pairs, TCGA / GTEx | **2,774 / 1,890** |
| S12A tumor-only / shared / normal-only event rows | **1,730 / 1,861 / 1,073**; R2 = **0.572 / 0.281 / 0.135** |
| S12B tumor-only / shared / normal-only event rows | **1,876 / 1,776 / 920**; R2 = **0.394 / 0.282 / 0.121** |
| CXCL12 expression | Phase **17.142 h**, p **3.5176e−11**, q **1.1370e−9**, n **178** |
| CXCL12 APA, ENST00000395794.2 | Phase **22.773 h**, p **0.0132635**, FDR **0.0383607**, n **156** |
| hsa-mir-135b | Phase **3.605 h**, n **178**; legacy p/FDR stored as zero from numerical underflow |
| miRNA–CXCL12 correlation | Spearman rho **−0.61895**, BH q **8.2076e−18**, n **178** |
| Triad timing | Target-minus-miRNA **−10.463 h**; APA-minus-expression **+5.631 h** |
| S13F site | Predicted **7mer-A1**, 111 nt into a 273-nt distal aUTR, minus-strand CXCL12 |

Source tables include `data/polyAPA/comparison/pairwise_overlap.csv`, `cohort_summary.csv`, `apa_expression_overlap_events_{tcga,cptac}_vs_gtex.tsv`, `apa_expression_phase_lag_second_moment_{tcga,cptac}_vs_gtex.tsv`; `data/apa_machinery_rhythmicity/machinery_expression_long.tsv`; and `data/apa_event_heatmap/splicing_machinery_tumor_vs_normal.tsv`.

The six CXCL12/miRNA panels are one original figure split between Figure 4M–P and S13E–F. Their source tables are under `data/apa_mirna_rhythmic_connection/`, including `ranked_example_candidates.tsv`, `selected_example_source_data.tsv`, `selected_example_mirna_site_apa_side_check.tsv`, and `selected_example_panel_legend.md`. These panels support a predicted regulatory hypothesis, not experimentally demonstrated miRNA regulation.

## Analysis settings and input lineage

### APA, machinery, and host-expression comparisons

APA cells **5–15** read per-sample `data/polyAPA/{TCGA,CPTAC,GTEx}/results/*/*_result_All_Prediction_Results.txt`, match CHIRAL sample phases, filter events, and fit intercept + cosine + sine models with a 24 h period. The local directories contain **183 TCGA, 213 CPTAC, and 362 GTEx** prediction files; the recorded analysis matches **178, 153, and 362** phased samples, respectively. Filtered/tested event counts are **23,789 / 80,031 / 50,856** for TCGA/CPTAC/GTEx.

- The configured phase input is **`data/Merged/CHIRAL_dual_penalty_phi_12h.csv`**, not `_12h_geo.csv`. The heatmap R script uses this same non-geo file, reads `sample` and `phi`, and expands comma-delimited sample aliases. The two phase files differ and must not be substituted silently.
- Expression-supported filtering requires both long- and short-isoform expression ≥20 in at least `floor(0.5 × sample count)` samples. If these columns lack usable expression values, the code instead requires nonmissing PDUI prevalence ≥50%. It then requires PDUI SD ≥0.05 and range ≥0.10, with at least 10 observations for fitting. The fallback is part of the existing cohort analysis and should be documented.
- APA figures primarily use **nominal p < 0.05**, despite also calculating BH FDR. Absolute amplitude is `sqrt(beta_cos² + beta_sin²)`; relative amplitude is that value divided by PDUI SD, **not by mean PDUI**. Violin comparisons are two-sided Mann–Whitney tests.
- Gene overlaps count unique mapped host ENSG IDs, not transcript events. Mapping uses `data/Ref/t2g.txt`. Their Fisher tests use `alternative='greater'` and the union of tested host genes as the universe.
- Machinery plots use expression regression tables `data/Merged/regression/all_genes_regression_periodic_12h_{cohort}.txt` and nominal p < 0.05. The significance heatmap uses −log10(p) × sign(amplitude), clipped to [−5,5], and the phase panel masks nonsignificant entries.
- Figure 4G/S12 instead match significant APA events to rhythmic expression genes from the `sig_genes_regression_periodic_12h_{cohort}.txt` files, additionally filtered to **q < 0.05 and amplitude >0.1** in `polyapa_phase_tools.load_expression_table`. Shared-group S12 rows combine both cohorts; these are event rows and can repeat host genes. Circular lag is centered on [−12,12), with von Mises bandwidth 0.3, 1,000 bootstrap samples, and 5,000 uniform-null simulations for axial R2.
- Figure 4A uses the older Gaussian-on-circular-distance density helper (bandwidth 0.3, 200 grid points); it is distinct from the later von Mises helper used for the lag plots. A unified implementation must make any kernel change explicit.
- S10 heatmaps use 24 one-hour phase bins and row z-scores of binned PDUI, ordered by event acrophase. Their saved binned means, z-scores, and event lists are under `data/apa_event_heatmap/`.

The non-geo phase filename has writers in the older `2.1.3_CHIRAL_dual_penalty_pipeline_{tpm_archive,cpm_backup}.ipynb` notebooks, while the current GEO notebook writes `_12h_geo.csv`. The cohort expression regression filenames are reused by multiple analyses, including GEO notebook cell **16**. Matching filenames alone does not establish a single consistent upstream execution; pin and verify phase/expression versions before migration.

### GO panel H: implementation details to preserve or explicitly revise

APA cell **27** computes a tested-gene background variable but does **not pass it** to `hypergeometry.count_numbers`. The helper therefore uses the GO-annotated proteome background from `data/enrichment/goa_human.gpa`, with `human_info.txt` and `go_name.txt` for mapping/labels. It calculates `1 - hypergeom.cdf(m, ...)`, i.e. the strict upper tail rather than the usual inclusive enrichment tail, and clips extremely small probabilities to 1e−22.

Cell **29** includes terms with −log10(p) >2 in **at least one cohort**, although comments/titles say two cohorts. The source dotplot uses cohort-specific top-term tables and significance for both dot size and color. It is **not** the later tumor-only/shared/normal-only host-set GO analysis in cells 34–35. Source tables are `data/polyAPA/{cohort}/enrichment/to_plot_GO_rhythmic_apa.txt` and `data/polyAPA/comparison/go_enrichment_comparison.csv`.

### miRNA landscape and selected triad

Both miRNA cell **5** and Triad cell **4** use the same **non-geo** `CHIRAL_dual_penalty_phi_12h.csv`. These are inferred phases of unordered samples, not observed collection times.

- TCGA and CPTAC starting inputs are processed, sample-named GDC quantification tables in `data/miRNA/TCGA/*.tsv` and `data/miRNA/CPTAC-3_Tumor/*.tsv`. The recorded runs load 178 and 146 valid tables. The loader reads root-level files with the expected miRNA quantification header, not nested original downloads. GDC sample sheets and the CPTAC rename manifest provide supporting provenance.
- GDC miRNA counts become CPM and then `log2(CPM+1)`. Filters require ≥50 reads and ≥1 CPM in ≥60% of samples, exclude miRNAs cross-mapped in any sample, and require transformed SD ≥0.15 and range ≥1.
- GTEx uses `data/miRNA/GTEx/raw/miRNA_TPM_matrix_PORTAL_2025_03_17.txt`, `rnacentral_mirbase.tsv`, and `data/Merged/GTEx_TCGA_merged_meta.csv`. It is a separate portal dataset, not a tumor-cohort proxy. Biological-sample-prefix matching yields **346 phased samples**. It uses `log2(TPM+1)`, ≥1 TPM prevalence ≥60%, SD ≥0.15, and range ≥1, without the GDC cross-mapping filter. The recorded ID audit maps 453 of 893 RNAcentral rows to human miRBase and retains 440 unmapped RNAcentral IDs. Saved row/sample mapping tables make this inspectable.
- Figure 4J/K miRNAs **and genes** use **nominal p <0.05**. This differs from the q/amplitude-filtered expression set in Figure 4G/S12. Gene phases come from `all_genes_regression_periodic_12h_{cohort}.txt`. Current phase-density inputs contain 134 miRNAs/12,689 genes for GTEx and 111 miRNAs/14,431 genes for TCGA.
- Predicted targets come from TargetScanHuman 8.0 `data/miRNA/targets/Summary_Counts.default_predictions.txt.zip` and `miR_Family_Info.txt.zip`, restricted to human, with `data/Ref/t2g.txt` for identifiers. The resolved `predicted_targets.tsv` contains 233,580 records. A separate miRDB reference exists but is not the plotted target tier; the validated-target tier is empty. Identifier normalization strips mature-arm suffixes and some locus suffixes, which must be recorded when joining GDC precursor names and GTEx mature miRNA names.
- Figure 4L uses pairs for which both miRNA and target have p <0.05. In-phase is within ±3 h; anti-phase is within 3 h of a 12 h offset. In TCGA/CPTAC/GTEx order, both-rhythmic pair counts are **48,241 / 20,097 / 30,559**; anti-phase counts **17,483 / 4,097 / 6,905**; in-phase counts **16,672 / 4,284 / 9,810**. Regenerating summaries are `data/miRNA/rhythmic_landscape/mirna_target_phase_relationship_summary.tsv` and the per-cohort phase-pair tables.
- S13B–D are **not all predicted targets or a random sample**: eligible pairs are ranked by `−log10(miRNA p) −log10(target p)`, deduplicated by normalized miRNA/target symbol, and capped at **20,000 per cohort/tier**. The primary eligible set requires both partners nominally rhythmic; the code falls back to rhythmic miRNAs alone only if that set is empty. Gene matrices `data/Merged/Group_TPM_{cohort}.csv` are collapsed by gene symbol and transformed to `log2(TPM+1)` for correlation. Spearman tests use common samples and at least 10 valid pairs; BH correction is over retained comparisons. Saved correlation tables and `mirna_target_expression_correlation_summary.tsv` record the actual population.
- The triad combines `all_resolved_mirna_target_pairs.tsv` with `data/polyAPA/comparison/apa_expression_phase_lag_{tcga,cptac}_vs_gtex.tsv`. Although its own nominal APA/host-expression filter is p <0.05, those overlap files already inherit the stricter expression selection described above. Candidate selection includes anticorrelation FDR <0.05, anti-phase timing, and a predicted distal functional seed of at least 7mer; the CXCL12 seed check and genomic annotations are saved with the selected-example tables. Plotting reloads filtered miRNA matrices, gene TPM, and filtered APA PDUI. This is a selected example, not an independent validation experiment.
- S13A clips −log10(p) at 15 and marks over-cap values with triangles. Numerical zero p values in the legacy F-test calculation must be treated as underflow, not mathematical zero.

### TMEM183A coverage

The manifest `data/apa_event_heatmap/TMEM183A_4samples_manifest.tsv` exactly matches screenshot samples/phases: **TCGA-2J-AAB6-01A (2.7 h), TCGA-YB-A89D-01A (8.7 h), TCGA-IB-AAUW-01A (14.1 h), and TCGA-HZ-7926-01A (21.0 h)**. APA cell 53 creates the selection; the separate renderer reads BAM coverage using `pysam` and GENCODE v49 annotation for ENST00000965092.1. The renderer's prose still refers to an older embedded cell label; use notebook position 53.

The original coverage PDF/PNG and manifest remain historical references. The four genomic BAMs and indexes are now supplied on GRC in `jptsvr`, under the project’s `data/figureS10D_bams/` directory. S10D is enabled in the publication notebook, with configurable paths in the first code cell. Coverage is newly extracted with the notebook’s `pysam.count_coverage` implementation and exported as 8,000 observed base-depth rows, together with BAM/index QC, a GENCODE v49 transcript model, and source checksums. Current PDUI-based selection independently reproduces the original four samples. Where the original BAMs are unavailable, the notebook can validate and plot these extracted source tables; it never substitutes a legacy image or a PDUI-derived depth schematic. The current median DaPars proximal marker is 203023501, replacing the historical fixed 203023481; GENCODE supplies the stop-codon start 203023038 and distal boundary 203023908.

### Pipeline boundary

The checked-in root `process.smk` currently sets `BATCH='geo_batch1'`, GENCODE v43 references, and the older `scripts/dapars/src/DaPars_main.py`, producing batch-specific `data/geo_batch1/...` paths. It does not directly declare the three cohort `data/polyAPA/` paths used here and does not generate miRNA quantifications or target annotations.

`pipeline/Snakefile_rmats.smk` and its Singularity counterpart contain the DaPars2 workflow and `aggregate_dapars_chromosomes` rule that declares `{sample}_result_All_Prediction_Results.txt`. `scripts/aggregate_dapars_chromosome_results.py` creates the cohort-compatible files. There are **213 CPTAC and 362 GTEx aggregation manifests** locally, but no such TCGA manifests. Example manifests record remote source paths with different directory names; configuration and manifests establish the relevant workflow family, not proof that the present default configuration generated every historical input. Preserve these cohort-specific provenance records and resolve the annotation/version differences before claiming reproduction from root `process.smk` alone.

## Earlier implementations to avoid confusing with these sources

[3.4.0_apa_tcga.ipynb](../3.4.0_apa_tcga.ipynb) is a TCGA-only predecessor, with **10,689** nominally rhythmic events in its stored run, versus **11,244** in the supplied figures/current multi-cohort artifacts.

[scripts/build_fig04_publication_notebooks.py](../scripts/build_fig04_publication_notebooks.py) is an earlier proposed reimplementation targeting `notebooks_publish/`. It calls the supplements **S8–S11**, corresponding by topic to the user's **S10–S13**. It recomputes phases and changes thresholds/selection, reduces some phase-lag analyses to genes, and substitutes an abundance schematic for the TMEM183A coverage track. It must not be run unchanged to reproduce these screenshots or the current publication conventions.

## Publication conventions for subsequent work

Follow [AGENTS.md](AGENTS.md): self-contained public notebooks in `notebook_publication/`; all analysis/plotting code included; processed-data inputs with external resources and upstream gaps explicitly documented; individual panels and assembled figures displayed inline; source tables under `data/publication/`; figures under analysis-specific `data/figures/publication/` subfolders; names `Figure{number}{panel}_{content_related}.png/pdf` with supplementary numbers S10–S13; 52.5 mm square modules and widths of 1–4 modules; minimal padding; globally configured fonts ≥7 pt and consistent styling; panel letters in PDFs only. Preserve all existing data, legacy notebooks, and figures. Carry these conventions into summaries and handoffs.

## Accepted publication threshold revision (2026-09-22)

The current Figure 4 notebook uses nominal P < 0.05 for APA events and miRNAs. All mRNA panels and cross-modality gene eligibility use the exact Figures 2/3 versionless ENSG sets: BH q < 0.05 and stored `Amplitude > 0.1`, from the TMM-normalized logCPM regression tables. Counts are TCGA 12,905, CPTAC 6,750 and GTEx 10,530. This supersedes the first generated notebook's universal-q rule and log2(TPM+1) gene refits. The historical source descriptions above remain unchanged as provenance, not current execution instructions. The original non-GEO TMM matrix and non-GEO CHIRAL phases numerically reproduce the stored gene regressions, including CXCL12; all modalities therefore use that verified phase input. The GEO TMM matrix does not reproduce these regression files and omits CXCL12. Gene selection remains exactly the same as Figures 2/3; input/fit concordance is exported.

## Panel D revision (2026-09-23)

The current publication notebook compares **TCGA ∪ CPTAC versus GTEx** (logical union of rhythmic genes, counted once), displaying the percentage of tested category genes passing the current q/amplitude mRNA rule. Each category uses a raw two-sided Fisher exact test of rhythmic/nonrhythmic counts, with asterisk thresholds at P < 0.05, 0.01, 0.001, and 0.0001, and no annotations for nonsignificant comparisons. The seven inherited category definitions are retained. ER calculations and paired-background randomization code have been removed from the active notebook. Superseded ER outputs are archived; historical descriptions above describe the original source, not the current panel.

## Panel H display revision

Panel H restores the original whole-GOA protein background (19,716 accessions), top 20 raw-P-ranked terms per cohort at P < 0.01, and global mean-significance ordering (absent selected entries = 0). It restores the gray plotting background, white grid, red scale, black dot edges, alpha 0.8, and proportional dot areas within the 4 × 2 module publication panel. Long labels are wrapped to retain ≥7 pt text. Display scores cap at 22, while exported P values remain uncapped. The corrected inclusive hypergeometric tail and current processed APA foregrounds are retained, so identical historical P values or term identities are not assumed. The active annotation parser continues to exclude NOT-qualified BP memberships.

## Panel I comparison revision

Panel I now compares the TCGA/CPTAC union with GTEx, instead of the historical TCGA-only comparison. Count unique tested and rhythmic APA host genes in the four splicing categories. Use nominal APA P < 0.05 for membership and raw two-sided Fisher P for asterisk annotations; omit nonsignificant labels/brackets. This changes the cohort contrast and its significance annotation, while keeping the inherited splicing category definitions.
