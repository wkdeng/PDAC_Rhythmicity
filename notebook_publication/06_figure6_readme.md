# Figure 6 and Supplementary Figures S16–S17

Run `06_figure6_cellline_validation.ipynb` from top to bottom with the Docker Python 3 kernel. The notebook contains its own downstream analysis and plotting code; it does not import legacy project helpers or reuse legacy figure images.

## Current results and input boundary

Circadian panels use the current CircadianMultiOmics results. S16F uses the embedded, literature-audited driver reference annotations (23 September 2026). Figure 6E/F use Hs 766T; S16J uses PANC-1 and S16K uses Panc 10.05. The current rhythmic expression counts are 2,151, 821, and 3,529 respectively. Figure 6A–C show MIA PaCa-2, Hs 766T, and Panc 08.13 BMAL1/PER2 reporter traces from the supplied Kronos workbook.

Inputs are processed BioCycle gene/event rhythmicity calls, current validation metric/cohort-alignment/driver tables, harmonic comparison metrics, an RBP catalog, and three expression time-series matrices. Gene rhythmicity uses nominal *p* < 0.05. These are explicit upstream analysis inputs; the notebook does not refit external BioCycle predictions or rerun alignment, quantification, or ChIP/CLIP inference. Direct linkage to `ChronoTherapy/process.smk` has not been established. Inputs must be distributed or made accessible separately for a public release.

The notebook recomputes per-line support and ranking, harmonic/BioCycle comparisons, best-per-axis selections, circular summaries and von Mises densities, core-clock calls, module rhythmic fractions, and expression heatmaps. Expression heatmaps use numeric ZT, average same-ZT replicates and duplicate symbols, and apply per-gene z scoring after matching current rhythmic gene symbols. Existing polarization/enrichment p-values and cohort-axis annotations remain processed upstream results.

## Paths and packages

- `CHRONO_KRONOS_PATH`: optional override for `data/kronos/Kronos.xlsx` under the publication project.
- `CHRONO_PROJECT_ROOT`: optional publication repository root.
- `CHRONO_CELLLINE_DATA_DIR`: current CircadianMultiOmics `data/` directory. Default: sibling `CircadianMultiOmics/data/` next to the publication project.
- Python packages: NumPy, pandas, SciPy, Matplotlib, scikit-learn, IPython, openpyxl. Notebook execution also needs nbformat, nbclient and ipykernel.
- No personal host paths or credentials are embedded in the notebook source.

Required inputs within `CHRONO_CELLLINE_DATA_DIR`:

- `cellline_validation/1708_comprehensive_biocycle/perline_metric_matrix.tsv`
- `cellline_validation/1708_comprehensive_biocycle/overall_alignment_to_cohort.tsv`
- `cellline_validation/1708_comprehensive_biocycle/driver_panel_with_metrics.tsv`
- `cellline_validation/1708_comprehensive/perline_metric_matrix.tsv`
- `biocycle_analysis/biocycle_all_results.csv`
- `cellline_validation/biocycle_events/biocycle_se_all.csv`
- `cellline_validation/biocycle_events/biocycle_apa_<line>.csv` for each of the 10 retained lines.
- `ref/rbp_rna_processing_links/rbp_catalog.tsv`
- `biocycle_analysis/inputs/input_HS766T_Dexa.txt`
- `biocycle_analysis/inputs/input_PANC1_Dx.txt`
- `biocycle_analysis/inputs/input_PANC1005_Dexa.txt`

The generated input manifest lists each input root and relative file path, sizes and SHA-256 hashes. The retained lines and analytical sources are exported alongside the figures' data tables.

## Publication outputs

- Source tables: `data/publication/fig06_cellline_biocycle_validation/`
- PNG/PDF figures: `data/figures/publication/fig06_cellline_biocycle_validation/`
- Validation records: `data/logs/publication/fig06_cellline_biocycle_validation/`

There are 31 individual/assembled PNG/PDF pairs: Figure 6A–L plus assembly, S16A–K plus assembly, and S17A–E plus assembly. All 31 PNGs appear inline. Panel letters appear only in the PDFs. All text sizes are configured by role and are at least 7 pt. Fixed-size canvases preserve 52.5 mm module dimensions; exports are not cropped with `bbox_inches='tight'`.

S17's 40 polar plots are grouped into two consecutive five-cell-line blocks, each containing the four phase layers. Both blocks together remain S17A, with a single panel letter. This preserves the current rank order while giving each polar plot enough room for readable text.

The assembled layouts use the following module grid (one module = 52.5 mm):

| Figure | Panel layout | Assembly size |
| --- | --- | --- |
| 6 | A–D each 1×1; E/F at 3×2 and 1×2; G–J each 1×1; K/L each 2×1 | 4×5 modules |
| S16 | Row 1 A/B/C: 1/2/1; row 2 D/E/F: 1/1/2; row 3 G/H: 2/2; row 4 I: 4; row 5 J/K: 2/2. All rows one module high. | 4×5 modules |
| S17 | A at 4×8; B–E each 1×1 | 4×9 modules |

These are fixed publication canvases rather than A4 page exports. S17 uses additional height to preserve 7 pt labels across all 40 polar plots.

See `figure6_s16_s17_provenance.md` for the source audit and `AGENTS.md` for persistent publication conventions.

## Numerical validation

The downstream calculations were compared with the current source outputs: rank
matched exactly; phase summaries (40 line/layer combinations), 13,357 density
points, the original 140 core-clock rows, 60 module rows, eight method-comparison rows, and all
three expression heatmaps matched to floating-point precision (maximum absolute
difference below 2 × 10⁻¹⁵). BioCycle gene identifiers are unique within each
retained line. The comparison record is saved as
`data/logs/publication/fig06_cellline_biocycle_validation/numeric_validation.json`.

## September 2026 annotation revision

Figure 6F–J and S16G use subscript R₂. S16E also follows this notation for consistency. Figure 6I is labeled “Axial R₂”, matching the existing second-moment calculation; no circadian values were changed.

S16F is titled “Driver annotations”. Its 40 gene-by-cell-line entries are audited against primary papers and curated Cellosaurus/ATCC records. The notebook embeds the complete evidence snapshot and exports `fig6_s16_s17_driver_evidence.tsv` (original call, displayed call, detailed variant, source URL, PMID/DOI, and qualifications). `fig6_s16_s17_driver_annotations_upstream.tsv` preserves the inherited table. The original project also stored `06_figure6_driver_evidence.tsv`; this data sidecar is omitted from the code-only export because the evidence records are already embedded in the notebook source.

Figure labels distinguish reported WT, sequence variant (`mut`), gene/multiexon deletion (`del`), and unresolved (`nd`). These are reference annotations, not matched-stock genotyping or a pathogenicity assertion. The audit corrects both Capan-2 WT calls (TP53 splice variant and CDKN2A duplication) and adds missing SMAD4 annotations. Panc 03.27 CDKN2A is now displayed as a reported deletion based on explicit ATCC Table14 evidence, with historical copy-number disagreement retained in the notes. Panc 10.05 CDKN2A is **WT**, based on explicit p16/CDKN2A wild-type statements in the [2016 Panc10.05 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5097013/) and the [later Pa16C study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10010283/). [Cellosaurus](https://www.cellosaurus.org/CVCL_1639) lists Pa16C as a synonym. SMAD4 remains **nd (conflicting reports)**: the 2016 study describes a Panc10.05 deletion, whereas the later study describes Pa16C as WT. These are reference annotations; synonymous names do not establish that different laboratories used identical stocks.

The supplied CCLE report was reconciled in `06_figure6_ccle_report_comparison.md`: Hs 766T CDKN2A is now a supported splice variant (zygosity unverified), Panc 03.27/Panc 08.13 CDKN2A are ATCC-reported deletions, and Panc 08.13 TP53 is unresolved. The report comparison and source-file hash are embedded in the notebook.

## Kronos reporter panels A–C

A–C show BMAL1 (blue) and PER2 (orange) in MIA PaCa-2, Hs 766T, and Panc 08.13, respectively, read directly from `data/kronos/Kronos.xlsx`. Workbook identifiers are `MiaPaca2`, `Hs766t`, and `Paca0813`. This is a separate experimental workbook, not a `process.smk` output. Its time values are interpreted as hours; synchronization metadata, signal units, biological-replicate identifiers, and prior signal processing are not specified in the file.

The analysis reproduces the settings of `3.8.0_kronos.ipynb` with all required code embedded here. Retain 0–60 h for MIA PaCa-2 and Hs 766T, and 24–84 h for Panc 08.13 (inclusive bounds). Display time is relative to each window start; retain the original slightly offset BMAL1/PER2 timestamps and do not interpolate endpoints. These are prespecified legacy windows, not windows optimized in this notebook.

Fit each retained raw-signal trace by ordinary least squares to `M + bc*cos(2πt/24) + bs*sin(2πt/24)`. The 24-hour period is fixed, not estimated. No additional detrending, smoothing or outlier removal is applied. Standardize observations and fitted curves with the retained trace's mean and population SD (`ddof=0`). Dashed horizontal lines mark the standardized MESOR; dotted vertical lines mark the first fitted peak within the observed window. Export raw and standardized observations, 1,000-point fitted curves, parameters, and missingness/window QC.

The summary retains the legacy nominal zero-amplitude F test and BH adjustment across six traces for reproducibility. Dense serial measurements are not biological replicates, and the strong residual autocorrelation violates the independent-error assumption. These probabilities do not establish replicate-level rhythmicity; no significance stars are shown in A–C. The fixed-period fits describe waveform timing and amplitude within the selected windows.

The notebook exports six `fig6_kronos_*.tsv` source tables and records the workbook hash. Legacy notebooks, scripts, and saved Kronos results are not required at runtime.

Validation on 24 September 2026: all six Kronos source tables match the legacy results exactly (2,093 input rows, 1,444 retained measurements, six fit summaries, and 6,000 curve points). All 22 pre-existing publication TSVs remain byte-identical. The notebook executed without errors, displays 29 PNGs inline, and exports 29 PNG/PDF pairs with the required module sizes, minimum 7 pt text, PDF-only panel labels, and no out-of-page text. See `data/logs/publication/fig06_cellline_biocycle_validation/kronos_validation.json`.

## Statistical notation

**Statistical notation.** Use lowercase italic *p* for p-values and lowercase italic *q* for FDR-adjusted values. The current Figure 6/S16–S17 panels contain no numeric p- or q-value annotations. Figure 6G–J encode nominal *p* < 0.05 by point color. Numerical statistical results remain in the source tables.

## TCGA comparison and S16 relabeling — 28 September 2026

S16A adopts the transcriptome PCA and S16B adopts the core-clock abundance heatmap from `3.9.0_cellline_tcga_representativeness.ipynb` (formerly analysis panels A and F). The publication notebook embeds the necessary reconstruction from processed matrices and does not load or execute that notebook. It retains the 163-patient PDAC-compatible cohort, 9,467 shared expressed protein-coding genes, 2,000 patient-selected variable genes, first-cycle phase-balanced cell-line TPM profiles, TCGA-only PCA fitting, and twelve clock-gene percentiles. The PCA cell-line colors match B's column-label colors; B is reordered to the established S16 composite ranking.

Additional inputs are `data/Merged/Group_TPM_TCGA-PAAD_Tumor.csv`, `data/reference/TCGA/clinical.tsv`, `data/reference/gencode_gene_coordinates.tsv`, and all ten `biocycle_analysis/inputs/input_<line>.txt` files. New source/QC tables are prefixed `s16_tcga_`. No genetics API input is needed for these two panels. Cross-assay quantification and bulk-tumor composition remain interpretation limitations; clock abundance is not a rhythmicity test.

| Final panel | Content | Previous label | Width × height modules |
| --- | --- | --- | --- |
| A | TCGA-trained transcriptome PCA | 3.9.0 A | 1 × 1 |
| B | Core-clock TPM percentiles in TCGA | 3.9.0 F | 2 × 1 |
| C | Composite alignment rank | S16A | 1 × 1 |
| D | Harmonic/BioCycle method comparison | S16B | 1 × 1 |
| E | Best example per conclusion axis | S16C | 1 × 1 |
| F | Driver annotations | S16D | 2 × 1 |
| G | Scaled effect-size matrix | S16E | 2 × 1 |
| H | Phase-value availability | S16F | 2 × 1 |
| I | Core-clock BioCycle hits | S16G | 4 × 1 |
| J | PANC-1 rhythmic expression | S16H | 2 × 1 |
| K | Panc 10.05 rhythmic expression | S16I | 2 × 1 |

The core-clock hit grid is transposed (cell lines as rows, genes as columns) to fit one full-width row at ≥7 pt. Values and hit thresholds are unchanged. The assembled S16 canvas is 210 × 262.5 mm. Earlier S16 exports are preserved in the figure archive.

## BMAL1 naming revision

Displayed gene symbols and the embedded TCGA comparison now use BMAL1, accepting ARNTL as an input alias. The core-clock grid contains 13 distinct genes (130 line-by-gene rows): the 10 unmeasured ARNTL alias rows were removed, retaining the measured BMAL1 rows. PCA coordinates, explained variance, clock percentiles and module timing remain unchanged within numerical precision. The alias-aware checks are saved in `data/logs/publication/fig06_cellline_biocycle_validation/bmal1_numerical_validation.json`.
