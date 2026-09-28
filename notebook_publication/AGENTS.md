# Publication notebooks

These user-specified conventions take precedence over legacy notebook and figure defaults.

- Put new public-ready, self-contained notebooks in `notebook_publication/`.
- Include all analysis and plotting code needed to regenerate results from processed inputs. The intended starting point is `process.smk` outputs; explicitly document any external benchmark matrices and any unverified upstream provenance instead of claiming a pipeline connection.
- Do not require legacy notebooks or project scripts at notebook runtime. Parameterize input paths and avoid personal host paths.
- Save figures in analysis-specific subfolders under `data/figures/publication/`; save regenerating source tables under `data/publication/`.
- Display every required individual panel and assembled figure inline.
- Preserve previous styling where compatible, but make styling and layout consistent across the project.
- The panel module is a square with a 52.5 mm edge (one quarter of A4 width). Use widths of 1, 2, 3, or 4 modules and whole-square layout units unless explicitly instructed otherwise.
- Minimize padding without clipping text or legends.
- Export individual panels and an assembled figure in both PNG and PDF. Names follow `Figure{number}{panel}_{content_related}.png/pdf`; use `S1`, `S2`, etc. for supplementary numbers and omit the panel letter in assembled filenames, e.g. `Figure1_assembled_chiral_benchmark`.
- PNGs have no panel labels. PDFs include panel labels. Figure 1's analysis panels begin at D. Exception: Supplementary Figures S1 and S2 are whole figures with no panel letters in either format; display only each whole figure inline. Their auxiliary component exports also omit visible panel letters (filename identifiers are retained).
- Render raw p-values as lowercase italic *p* and FDR-adjusted values as lowercase italic *q* across all publication panels, including annotations, titles, legends, and color scales. This latest user convention supersedes earlier uppercase P/Q preferences.
- Use BMAL1 for human displayed gene symbols and Bmal1 for mouse displayed gene symbols in publication figures, notebook tables and legends. Figure 1 H/I use Bmal1 because the data are mouse. Accept ARNTL/Arntl and MOP3 as input aliases where needed; preserve source identifiers used for matching and avoid duplicate ARNTL/BMAL1 entries.
- Configure fonts globally. Equivalent elements (ticks, axis labels, legends, titles, annotations, etc.) use the same size across figures. No text below 7 pt unless explicitly authorized.
- Use numbered code-cell comments, deterministic seeds, input/QC checks, exported source tables, and runtime/package information.
- Preserve original data, legacy notebooks, and legacy figures. New outputs must not overwrite them.
- Preserve these rules in context summaries and handoffs.

## Figure 2 decisions

- Use current processed results for Figure 2 and S4–S6.
- Heatmaps must match current significant ENSG identifiers to the current expression matrix and report unmatched identifiers.
- Figure 2G displays TCGA expression in the same explicitly sorted GTEx-acrophase gene order as Figure 2E.
- Figure 2K uses the corrected protein-coding ChIP protocol: shared tested/promoter-assessed background; GTEx rhythmic reference versus current TCGA-specific ENSG identifiers rebuilt from current all-gene regression tables; nine TFs plus pooled CCG peaks; BH correction across all 20 tests. Recompute from processed promoter peak counts and GENCODE biotypes, and use BH q-values for colors and stars. Quantitative comparisons with historical screenshots belong in the conversation, not in the generated analysis notebook.

- Figure 2K display order: BMAL1, CLOCK, CRY1, NPAS2, NR1D1, RORA, RORC, DBP, TEF, pooled CCG peaks. Figure 2L explicitly labels its target as TCGA-specific; use the same current TCGA-specific ENSG set as K, intersected with the existing PPI background, and retain the ten-clock-gene pool. S5J/K use the same current GTEx and TCGA targets plus CPTAC-specific genes defined by the same rule.

- Figure 2 palette: GTEx `#456BAB`, TCGA `#F2B770`, CPTAC `#C46456` for cohort-coded main panels and related S4–S6 panels. User explicitly excludes panel K from recoloring: retain its red BH-q scale. Quantitative heatmap/polar significance scales are not cohort colors.
- Figure 2 layout: rows retain A/B (3:1), C/D (3:1), E–H, and I/J beside stacked K/L. I and J are 1 × 1 modules. K is 2 × 0.6 modules (105 × 31.5 mm); L is 2 × 0.4 (105 × 21 mm), an explicit height exception. Their combined block is 2 × 1. The assembly is 4 × 5 modules. Minimum text remains 7 pt except the explicit B/D gene-label exception below.

## Figure 3 decisions

- Figure 3 layout: row 1 A–D at 1 × 1 module each; row 2 E at 3 × 1.5 and F at 1 × 1.5; row 3 G–J at 1 × 1 each; row 4 K–N at 1 × 1 each. The assembled figure is 4 × 4.5 modules (210 × 236.25 mm). The 1.5-module height is explicitly requested.
- Include the real SLK BAM coverage panel E in both individual exports and the Figure 3 assembly at 3 × 1.5 modules. Draw the first-to-third exon skipping junction above coverage and the inclusion junctions below. Configure a hypothetical BAM directory and all four sample file paths in the first code cell. Retain extraction code, validated source bundle, and BAM input/QC provenance; no synthetic coverage.

- Figure 3 annotations: G/H/I titles are “TCGA - APA machinery”, “GTEx - APA machinery”, and “APA machinery rhythmicity” as requested; retain the existing splicing-machinery gene set. A/B/C/G/H P-value and n annotations use regular weight. I compares the TCGA/CPTAC union with GTEx, labels that union explicitly, and shows only a Core asterisk for two-sided Fisher P < 0.05. SR remains 9/14 versus 9/14; do not alter counts to match the unexplained historical screenshot.

- Figure 3N displays significance ranges as stars for the existing two-sided Fisher comparisons against GTEx: * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001 (use the most stringent matching threshold); ns P ≥ 0.05. State the ranges in the notebook caption.

- Figure 3A uses a fixed one-module layout: a wide heatmap, full domain labels, centered title with regular-weight subtitle, horizontal cohort labels, and a horizontal color scale below. Avoid automatic layout shrinking the heatmap into a narrow strip.

- Figure S9 layout: three rows, A–D then E–H then I–K. All panels are 1 × 1 module; I–K are centered at equal intervals across the full-width third row. Assembly dimensions are 210 × 157.5 mm.

## Figure 4 decisions

- Figure 4D compares the TCGA ∪ CPTAC union with GTEx using separate percentage-of-tested-genes-rhythmic bars. Union membership is the logical OR of the two tumor rhythmic gene sets, counting each gene once. Use current mRNA q/amplitude thresholds. Brackets show raw two-sided Fisher exact P values from each category’s rhythmic/nonrhythmic counts: * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001; omit nonsignificant annotations and brackets. No ER or background-normalization/randomization code remains in the active notebook. Preserve the inherited category membership unless separately requested.

- Use `04_figure4_apa_mirna_coordination.ipynb` for Figure 4 and S10–S13; retain their audited panel subjects.
- Use nominal P < 0.05 for rhythmic APA events and miRNAs. For mRNA, match the exact Figures 2/3 versionless ENSG lists: BH q < 0.05 and stored Amplitude > 0.1 from the TMM-normalized logCPM regression tables. Verify equality to the saved significant-gene sets. Use the original non-GEO TMM expression matrix and non-GEO CHIRAL phases, which numerically reproduce the authoritative regression tables. Canonicalize aliases and export regression-reproduction QC. Keep assay-specific input QC explicit.
- Accept the existing cohort APA quantifications as processed inputs; their pipeline origin is not a blocker.
- S10D is enabled using actual indexed TMEM183A RNA-seq BAMs, with configurable BAM directory, four sample paths, and GENCODE v49 GTF in the first code cell. Prefer direct BAM extraction; permit checksum-validated real BAM-derived source tables when BAMs are unavailable, checking agreement with the current sample selection and fit. Export observed per-base depth, sample/BAM QC, transcript features and checksums. Include S10D (4 × 3 modules) below adjacent S10A–C (each 1 × 3, one shared assembly color bar) in the complete 4 × 6 assembly. Do not fabricate coverage or use legacy images.
- Report methodological and numerical changes to the user after generation. Keep the CXCL12/miR-135b example prespecified and explicitly report its current eligibility.

## Figure 6 decisions

- Figure 6 and S16–S17: use lowercase italic p for p-values and lowercase italic q for FDR-adjusted values in all displayed statistical notation, overriding the general uppercase convention. Current panels have no numeric p/q annotations; G–J encode nominal p < 0.05 by color. Do not add tests or annotations solely to apply typography.

- Use current processed CircadianMultiOmics results consistently for Figure 6 and S16–S17.
- Figure 6A–C are Kronos BMAL1/PER2 reporter traces from `data/kronos/Kronos.xlsx`: A MIA PaCa-2, B Hs 766T, C Panc 08.13, each 1 × 1 module beside D. Preserve fixed 24 h cosinor fits, 0–60 h windows for A/B, 24–84 h for C displayed from zero, and per-retained-trace z scoring (ddof=0). Embed reader/fitting code; export source/QC tables; keep workbook provenance separate from process.smk. Use blue BMAL1/orange PER2, dashed MESOR and dotted first acrophase, with a shared legend in the assembly.
- Current S16J is PANC-1 and S16K is Panc 10.05.
- S16/S17 use the standard publication labeling rule (PDF letters, no PNG letters); the whole-figure no-letter exception is only for S1/S2.
- Use subscript R₂ for Figure 6F–J and S16 effect annotations; Figure 6I has axis title “Axial R₂” (the existing second-moment statistic is unchanged).
- S16F title is “Driver annotations”. Its genotype table must carry gene-by-cell-line literature/database evidence; distinguish sequence variants, deletions, supported WT calls, and unresolved calls. Do not infer WT from an omitted database annotation. Preserve upstream calls in the audit trail.
- After CCLE-report reconciliation: Hs 766T CDKN2A is a splice variant supported by NCI MutSpliceDB c.457+2T>C (zygosity unverified); Panc 03.27/Panc 08.13 CDKN2A are ATCC-reported homozygous deletions; Panc 08.13 TP53 is unresolved because a DepMap LoF subtype does not expose the exact alteration. Keep Panc 03.27 SMAD4 as ATCC-reported deletion and Panc 10.05 CDKN2A as reported WT (explicit Panc10.05/Pa16C statements in PMC5097013 and PMC10010283). Keep Panc 10.05 SMAD4 conflicting/unresolved: the former reports deletion and the latter Pa16C WT. Retain historical disagreements and do not claim current-release DepMap mutation/CN validation without the underlying records.

## Figure 1 decisions

- Figure 1 and S1–S3 only: render any p-value symbols as lowercase italic p and any FDR/q-value symbols as lowercase italic q. This scoped user instruction overrides the uppercase probability-symbol convention above. Their current figures contain no p-value or q-value annotations; do not add statistical tests or labels merely to apply this style.

- Figure 1 only: panel titles and facet titles are reduced by 1 pt (8 pt and 7 pt respectively); annotations/subtitles are black and reduced from 7 to 6 pt by explicit user request. Other text roles retain their configured sizes. Exception for panel J: heatmap numbers have transparent backgrounds and use black or white text based on the actual cell color for maximum contrast.
- Supplementary Figures S1 and S2 have no overall figure titles or subtitles; retain their existing subplot descriptions and supplementary font sizes.

- Figure 1 H–J height exception: use 1.33 modules (69.825 mm), with unchanged widths 1:1:2; D–G stay one module tall. The assembled Figure 1 is 4 × 2.33 modules.


- Figure 2 annotation refinements: center all panel titles; overlay PDF letters without a separate title row. B/D gene annotations alone use 5 pt (explicit user exception) and sit beside points without leader lines. Use lowercase italic p on raw significance scales and lowercase italic q (BH-adjusted) on K’s right-side color bar. A/C subtitles include sample n; E–G combine order and n on one line; G title identifies TCGA expression of GTEx rhythmic genes. H range is 0–15000; I has a title/time axis/thicker curves; J uses subscript R2; L range is 0–1.5.

- Figure 2 A/C fit annotations use lowercase italic p. B/D use left/bottom amplitude/phase titles, compact margins, and an enlarged circular plotting region. Figure 2H bar width is 0.45.

- Figure 2/S5 rhythmicity consistency: globally configure BH q < 0.05 and peak-to-trough amplitude > 0.1 (log2 CPM). Tumor-specific means rhythmic in the tumor, tested in both GTEx and CPTAC normals, and not passing that same rule in either normal. Rebuild from current regression tables; do not load stale cancer-specific lists or enrichment statistics. K retains its shared protein-coding promoter universe; L/S5J–K retain the processed BioGRID symbol universe and raw-P display. Export target membership and 2×2 counts. S5 uses Figure 2 cohort colors, direct polar labels (7 pt), Gaussian circular-distance phase density (bandwidth 0.3 radians), subscript axial R2, and explicit TCGA-specific/CPTAC-specific PPI row labels.

- Figure 4H uses the original whole-GOA annotated-protein background, raw enrichment P < 0.01, and top 20 terms per cohort (ties by GO ID). Keep the corrected inclusive hypergeometric tail. Order terms globally by mean capped −log10(P), assigning zero to absent selected-cohort entries; stronger terms run left-to-right. Restore gray background/white grid, red dots with black edges and alpha 0.8, and proportional sizes adapted to the 4 × 2 module panel. Display scores cap at 22; exported P values remain uncapped. Wrap long labels and retain ≥7 pt text. BH q values remain in source tables.

- Figure 4I compares TCGA ∪ CPTAC with GTEx using percentages of splicing-machinery genes hosting any rhythmic APA event (nominal P < 0.05). Use union tested host genes as the tumor denominator and union rhythmic host genes as numerator, counting genes once. Brackets use raw two-sided Fisher P thresholds (* .05, ** .01, *** .001, **** .0001); omit nonsignificant annotations. Keep this distinct from mRNA-expression rhythmicity and leave supplementary comparisons unchanged.


## Figure 5 decisions

- Figure S14 omits the “No term passes q < 0.05” plot annotation. Figure S15E omits the APA PDUI / SE PSI color legend; Figure 5J retains its legend. Apply these choices to individual exports, assemblies, and inline notebook outputs.

- B–D show regular-weight `n=48, R₂=0.34`-style annotations below the title, with subscript 2 and no axis-hour annotation. Apply the same presentation to the repeated CPTAC phase panel S15A.
- F/G use a regular-weight FDR threshold subtitle; keep the plot title bold. Apply the same threshold styling to the related S15B enrichment panel.
- Remove the explanatory “Prespecified evidence criteria; rows ordered by criteria met, then event support” axis text from Figure 5K and S15F, retaining the criterion names.

- Figure 4 layout: row 1 A–D, each 1 × 1 module; row 2 E/F each 1 × 1 and G 2 × 1; row 3 H 4 × 2; row 4 I–L each 1 × 1; row 5 M–P each 1 × 1. Assembly is 4 × 6 modules (210 × 315 mm). E/F use one shared machinery legend in the assembly and enlarged clock faces. Figure 4B/C display raw P values from their existing tests, retaining q values in source tables. G uses darker normal-only blue and more opaque points for readability.

- Figure 2 and S4–S6: use lowercase italic p for p-values and lowercase italic q for FDR-adjusted values everywhere, including fit annotations, permutation/KS labels, axes and color bars. This user instruction supersedes the general uppercase convention for these figures. Statistical values and star thresholds are unchanged.

- Figure 5 and S14–S15: use lowercase italic p for nominal p-values and lowercase italic q for FDR-adjusted values in every figure annotation, subtitle, axis, legend, and color scale. This scoped user instruction overrides the uppercase probability-symbol convention. Retain the existing numerical values and thresholds.

- Figure 4 and S10–S13: use lowercase italic p for raw p-values and lowercase italic q for FDR-adjusted values throughout annotations, titles, axes, legends, and color bars. This supersedes the uppercase convention for these figures. Preserve the existing raw/adjusted choice, numerical results, and star thresholds.


## Figure S16 TCGA integration (28 September 2026)

- S16A is the transcriptome PCA from root 3.9.0 analysis panel A; S16B is its core-clock TCGA-percentile heatmap (analysis F). Embed regeneration code in the Figure 6 publication notebook; do not require another notebook or its exported results at runtime. Retain the validated 163-patient primary cohort, feature selection, TCGA-trained PCA, first-cycle TPM aggregation and clock percentiles.
- Final S16 mapping: A PCA, B clock abundance, C composite rank, D method robustness, E best examples, F driver annotations, G effect sizes, H phase counts, I core-clock rhythmic hits, J PANC-1 heatmap, K Panc 10.05 heatmap. Earlier S16A–I shift to C–K.
- S16 five-row layout: A/B/C widths 1/2/1; D/E/F widths 1/1/2; G/H widths 2/2; I width 4; J/K widths 2/2. Every row is one 52.5-mm module high; assembly is 4×5 modules. Transpose I's grid to retain readable labels without changing values. All fonts remain ≥7 pt; individual and assembled PNG/PDF exports and inline displays remain required.
