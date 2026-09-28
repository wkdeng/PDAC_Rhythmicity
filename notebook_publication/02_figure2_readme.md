# Figure 2 and Supplementary Figures 4–6

Open `02_figure2_cohort_clock_organization.ipynb` and run all cells with the R (`ir`) kernel. The notebook is self-contained: statistical functions, plotting functions, layout, and exports are embedded. It does not source legacy notebooks or scripts.

## Execution and inputs

The default environment is the existing Docker Jupyter runtime. `CHRONO_PROJECT_ROOT` selects a portable repository root; `CHRONO_DATA_DIR` optionally selects a separate processed-data directory. Required R packages are checked in the first cell: dplyr, tidyr, tibble, readr, data.table, ggplot2, patchwork, ragg, ggrepel, scales, digest, jsonlite, and stringr. IRkernel and IRdisplay provide notebook execution and inline PNGs.

The input boundary is the **current processed analysis products**, not raw sequencing or raw peak/interaction files. The notebook recomputes harmonic fits, SVD, heatmaps, gene overlaps, circular density, bootstrap/permutation statistics, common/unique comparisons, ChIP enrichment, and PPI enrichment. It retains saved CHIRAL phases, verifies current significant-gene membership, and rebuilds tumor-specific targets from current all-gene regressions. It does not rerun upstream inference, read quantification, peak-to-gene assignment, or raw BioGRID parsing. A direct `process.smk` provenance link is not assumed.

Files required relative to the data directory:

- `Merged/GTEx_TCGA_GEO_merged_meta.csv`
- `Merged/GTEx_TCGA_GEO_merged_TMM_unfiltered_expression.csv`
- `Merged/CHIRAL_dual_penalty_phi_12h_geo.csv`
- `Merged/regression_periodic_12h_geo.txt`
- `Merged/regression/sig_genes_regression_periodic_12h_<cohort>.txt` for `GTEx_Normal`, `TCGA-PAAD_Tumor`, and `CPTAC-3_Tumor`
- `Ref/genes.bed` for identifier checks
- `gtex_chiral_sensitivity_supplement/harmonic_regression_all_genes.tsv`
- `gtex_chiral_sensitivity_supplement/rhythmic_genes.tsv`
- `gtex_chiral_sensitivity_supplement/chiral_phase_estimates.tsv`
- `chip_background_sensitivity/assessed_promoter_peak_counts.tsv` (corrected ENSG mapping, including assessed zero-binding promoters)
- `chip_background_sensitivity/eligible_gene_membership.tsv` (shared tested/promoter-assessed gene membership)
- `reference/gencode_gene_coordinates.tsv` (GENCODE v49 gene biotypes)
- `Merged/regression/all_genes_regression_periodic_12h_<cohort>.txt` for GTEx and CPTAC normals and TCGA and CPTAC tumors
- `ppi_panel_l_background_comparison/background_symbols.tsv` (processed BioGRID physical-interaction universe)
- `enrichment/ppi_enrichment/clock_gene_ppi_partners.txt`

Include these files, or documented access routes, when publishing the repository. The generated input manifest stores file sizes and SHA-256 hashes.

## Outputs and conventions

- Figures: `data/figures/publication/fig02_cohort_clock_organization/`, divided into `figure2/` and `supplementary4/`–`supplementary6/`.
- Source tables and QC: `data/publication/fig02_cohort_clock_organization/`.
- Validation logs: `data/logs/publication/fig02_cohort_clock_organization/`.
- 32 PNG/PDF pairs: 12 Figure 2 panels plus assembly, two S4 panels plus assembly, 11 S5 panels plus assembly, three S6 panels plus assembly.
- All 32 PNGs appear inline. PNGs omit panel letters; PDFs include them. Text roles are globally configured at 7 pt or larger, except B/D gene annotations at the user-requested 5 pt.
- Panel widths are multiples of 52.5 mm. The export manifest records nominal sizes; the R PDF device rounds page dimensions to whole PostScript points (less than 0.36 mm difference).

Figure 2E/G use the same current GTEx gene membership and exactly the same GTEx acrophase ordering. Figure 2G displays TCGA expression. Figure 2F uses TCGA genes and TCGA ordering. Current matches are 10,257 GTEx rows and 12,576 TCGA rows; the unmatched-identifier table records why these differ from the 10,530/12,905 significant-table totals. Heatmap colors are gene-wise z-scores across 1 h bin means, not log fold changes.

The clock-gene mapping is checked against annotation before fitting or plotting. Main SVD retains the 12-gene set with CIART and without CLOCK; expression and restricted-GTEx SVD use the 12-gene set with CLOCK and without CIART; amplitude polars use all 13. A one-cohort SVD has one column and therefore explains 100% by construction. Restricted-GTEx phases remain cached inputs; the donor selection interval is not asserted to be a hard bound on final inferred phases.

Scatter annotations use harmonic half-amplitude; phase/amplitude polars use twice that value (peak-to-trough amplitude). S6 preserves the original numerical threshold transformation, `log2(2 * stored Amplitude + 1)`, and labels it explicitly. Figure 2K recomputes ChIP enrichment from corrected promoter peak counts using the shared tested, promoter-assessed protein-coding background (16,570 genes). GTEx targets are 8,076 protein-coding rhythmic genes; TCGA targets are 3,526 protein-coding genes from the current 4,953-gene tumor-specific ENSG set. All cohorts use BH q < 0.05 and peak-to-trough amplitude > 0.1 on the log2 CPM scale. A tumor-specific gene is rhythmic in that tumor, tested with finite q and amplitude in both GTEx and CPTAC normal cohorts, and does not pass the same rhythmicity rule in either normal. This differs from pairwise tumor-minus-GTEx membership and does not establish differential rhythmicity. Both tumor-specific sets are rebuilt from current regression tables; saved method_A_presence lists are no longer used. The CPTAC-specific set contains 2,081 ENSG identifiers. Gene-level membership, normal/tumor values, and mapping QC are exported.

ChIP binding is at least one ChIP-Atlas peak (score >=200) overlapping the GENCODE v49 GRCh38 promoter (TSS +/-2 kb). The pooled CCG column is the union across all nine TFs, counted once per gene. Upper-tail hypergeometric tests are independently checked against one-sided Fisher exact tests. BH correction covers all 20 tests (two cohorts x nine TFs plus pooled CCG). Numbers show fold enrichment; colors show -log10(BH q), and stars use q<0.05/0.01/0.001. Counts, raw p-values, BH q-values, target membership and QC are exported under `data/publication/fig02_cohort_clock_organization/`. The notebook embeds the full computation and does not require the standalone sensitivity scripts or their final result tables.

Panel K columns are ordered BMAL1, CLOCK, CRY1, NPAS2, NR1D1, RORA, RORC, DBP, TEF, then pooled CCG peaks. Panel L compares GTEx with TCGA-specific: 7,699 and 3,372 mapped genes respectively. S5J–K use these same current target definitions plus 1,362 mapped CPTAC-specific genes. All PPI counts and tests are recomputed from the processed partner sets and 22,695-symbol background. The existing ten-clock-gene pool (633 unique partners) and raw-P display convention remain. This broader legacy network universe has no protein-coding/shared-tested restriction; it includes symbols from physical interactions without filtering both interactors to human. It is retained for comparability, not asserted to be the preferred inferential universe. The notebook does not load saved enrichment statistics.

S5 retains its four-row modular arrangement. Polar labels are beside points without leader lines, at 7 pt; SVD axes identify mode-1 loading. Its overlap plot uses the same bar width and 0–15,000 range as Figure 2H. Restricted-GTEx rhythmicity is verified against the same thresholds and its density uses the same Gaussian circular-distance kernel (bandwidth 0.3 radians) as Figure 2I. The axial statistic and prespecified SVD gene sets are unchanged. PPI heatmaps explicitly label both tumor-specific sets.

## Cohort palette and compact Figure 2 layout

Cohort-coded main and supplementary plots use GTEx `#456BAB`, TCGA `#F2B770`, and CPTAC `#C46456`. Panel K retains the red BH-q scale by explicit request. S6 shared-with-GTEx genes are blue, TCGA-unique genes orange, and CPTAC-unique genes red; the target/statistical definitions remain unchanged.

Figure 2 retains the reference panel arrangement, with I/J each 52.5 × 52.5 mm beside K/L. K is 105 × 31.5 mm (2 × 0.6 modules), and L is 105 × 21 mm (2 × 0.4). These two heights are user-authorized exceptions to whole modules. Rows occupy 1, 1, 2 and 1 modules, making the assembly 210 × 262.5 mm. K uses angled gene labels and a right-side color bar labeled with lowercase italic q (BH-adjusted q-value); stars mean q<0.05, q<0.01 and q<0.001 for one, two and three stars. L retains the original raw-p stars. PNGs omit panel letters; PDFs include them.

All panel titles are centered. PDF letters overlay the upper-left margin without reserving a title row. B/D labels sit beside their dots with no leader lines. A/C subtitles include sample counts; heatmap subtitles combine acrophase order and gene count. H spans 0–15,000 genes; L spans 0–1.5-fold enrichment. I has a circadian-time axis and thicker density curves. J uses subscript 2 in its axial-polarization title and axis. Raw p-value scales use lowercase italic p.

Figure 2 A/C fit annotations use lowercase italic p. B/D use 5 pt gene labels,
a larger circular plotting region, a left-side amplitude axis title, and a bottom
phase axis title. H uses bar width 0.45 (previously 0.7). Panel dimensions and all
plotted values remain unchanged.

All displayed probability symbols in Figure 2 and S4–S6 use lowercase italic p (raw tests) or q (FDR-adjusted tests). This includes the S4 harmonic-fit and S6C KS annotations.
