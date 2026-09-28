# Figure 3 and Supplementary Figures S7–S9

Open [03_figure3_rhythmic_alternative_splicing.ipynb](03_figure3_rhythmic_alternative_splicing.ipynb)
in the project's Docker **Python 3** kernel and run all cells. The notebook contains
the analysis and plotting code; it does not execute legacy notebooks or project scripts.

## Inputs and execution

The default recomputes harmonic AS fits from processed rMATS JC tables. Set
`CHRONO_RECOMPUTE_AS=false` before launching the kernel to use the existing AS fit
tables and recompute downstream statistics and figures. The supplied executed
notebook uses this cached-fit mode. A separate direct TCGA SE refit reproduced
2,362 fitted events, 387 nominally significant events, and SLK's 4.5801-hour phase.
All other AS event-type/cohort refits have not been independently rerun during this
publication export.

Set `CHRONO_PROJECT_ROOT` for another repository location and `CHRONO_DATA_DIR` for
an external processed-data directory. Install the packages in
[03_figure3_requirements.txt](03_figure3_requirements.txt) in a Python 3 environment.
The recorded export environment uses Python 3.13.2. Notebook paths contain no
personal macOS directories.

Required inputs, relative to the data directory:

| Input | Location |
|---|---|
| Cohort metadata | `Merged/GTEx_TCGA_GEO_merged_meta.csv` |
| CHIRAL sample phases | `Merged/CHIRAL_dual_penalty_phi_12h_geo.csv` |
| rMATS sample order | `AS/{CPTAC,TCGA,GTEx}/all_samples.txt` |
| SE, A3SS, A5SS, MXE, RI quantifications | `AS/{project}/post/{event_type}.MATS.JC.txt` |
| Expression regression tables | `Merged/regression/{all,sig}_genes_regression_periodic_12h_{subcohort}.txt` |
| GO enrichment | `enrichment/full_term_GO_Regression_CHIRAL_12h_geo.txt` |
| Reactome enrichment | `enrichment/full_term_Reactome_Regression_CHIRAL_12h_plus_geo.txt` |
| KEGG enrichment | `enrichment/full_term_KEGG_Regression_CHIRAL_12h_plus_geo.txt` |
| Clock-TF peak–gene matrix | `enrichment/chip_atlas_enrichment/peak_gene_matrix.txt` |
| Optional cached AS fits | `AS/{project}/rhythmic_results/{cohort}/{event_type}/all_results.csv` |

Subcohorts are `CPTAC-3_Tumor`, `TCGA-PAAD_Tumor`, and `GTEx_Normal`; AS cache
cohorts are `CPTAC_Tumor`, `TCGA_Tumor`, and `GTEx_Normal`.

These are explicitly declared upstream processed inputs. The checked-in
`process.smk` is configured for a GEO batch and does not directly produce all of
these cohort files. CHIRAL inference, expression regression, enrichment testing,
and construction of the ChIP-Atlas matrix are upstream of this notebook. Supply
these inputs separately with a public repository; they are not embedded in the
notebook. Exact file sizes and SHA-256 hashes are exported in
`processed_input_manifest.tsv`.

## Analysis and figure details

- AS calls retain the legacy nominal `p < 0.05` definition and pre-cohort
  PSI-availability, variance, range, and junction-coverage filters. Parameters,
  sample order, phase overlap, and fit counts are exported for inspection.
- Figure 3J/K compare marginal SE and rhythmic-expression populations. S9A–D
  restrict to SE events with rhythmic hosts, deduplicating host genes for the
  host-phase distribution. S9I–K retain event-specific signed phase differences
  for all five AS types and show the complete −12 to +12-hour circle.
- Clock-TF fractions count matched **events**, so genes with several events
  contribute several rows. Intervals and tests retain that legacy analysis unit.
- Figure 3I uses significance in either tumor cohort. Its current counts differ
  from the supplied historical screenshot; the notebook uses the current data.
- Figure 3E is exported individually as real SLK BAM coverage and measured
  junction counts and is included in the Figure 3 assembly. The first-to-third
  exon skipping connection is above coverage; inclusion connections are below.

See [the panel provenance map](figure3_s7_s9_provenance.md) for legacy sources.
Comments embedded in the reference screenshots were not adopted as new analysis
or panel-rearrangement instructions.

## Outputs and validation

Figures are under `data/figures/publication/fig03_rhythmic_alternative_splicing/`,
in `figure3/`, `supplementary7/`, `supplementary8/`, and `supplementary9/`.
There are **34 PNG/PDF pairs**: 30 individual panels and four assemblies. All 34
PNGs are displayed inline. PNGs omit panel letters; PDFs include them. Global
fonts are at least 7 pt and panel dimensions use 52.5 mm modules. Tall pathway and
gene lists retain sufficient height for readable text.

Source tables, input hashes, export dimensions, thresholds, software versions,
and configuration are under `data/publication/fig03_rhythmic_alternative_splicing/`.
Validation checks cover notebook execution, expected exports, minimum font sizes,
circular-statistic identities, synthetic harmonic-fit recovery, Wilson intervals,
sample matching, and the direct TCGA SE refit. Legacy files are preserved.

### Figure 3 layout

| Row | Panels | Module sizes (width × height) |
|---|---|---|
| 1 | A–D | 1 × 1 each |
| 2 | E, F | E: 3 × 1.5, BAM coverage; F: 1 × 1.5 |
| 3 | G–J | 1 × 1 each |
| 4 | K–N | 1 × 1 each |

The Figure 3 assembly is **210 × 236.25 mm**. The 1.5-module height is an
explicit author-requested exception. The layout revision regenerates Figure 3
from the exported analysis tables and updates its inline images; it does not
refit the analyses or change Supplementary Figures S7–S9.

### Machinery annotations

Panel I labels the tumor gene union **TCGA/CPTAC** and displays a Core asterisk
for two-sided Fisher P < 0.05; other categories show no P-value annotation.
`Figure3I_machinery_category_counts.tsv` and
`Figure3I_machinery_gene_membership.tsv` retain the underlying statistics and
membership. SR is 9/14 (64.3%) in both the tumor union and GTEx, consistent with
the saved legacy notebook output; see the provenance map for the screenshot audit.
G/H/I use the author-requested “APA machinery” titles; the analyzed categories
remain the original splicing-machinery gene set. A/B/C/G/H annotations use regular
font weight, independently of the bold panel titles.

Figure 3N uses significance stars for the existing two-sided Fisher comparisons
against GTEx: * 0.01 ≤ P < 0.05; ** 0.001 ≤ P < 0.01;
*** 0.0001 ≤ P < 0.001; **** P < 0.0001; ns P ≥ 0.05.

Figure 3A reserves fixed space for domain labels and a wide heatmap within its
52.5 mm square; the color scale is horizontal below the cohort labels. Its
counts, enrichment values, and color normalization are unchanged.

### S9 layout and statistical typography

S9 is 210 × 157.5 mm: row 1 A–D, row 2 E–H, and row 3 I–K. Each panel
is 52.5 mm square; the third-row panel centers are equally spaced across the
full row. Individual exports use those same dimensions. Figure 3 and S7–S9
render raw p-values as lowercase italic *p* and FDR-adjusted values as lowercase italic *q* throughout. The displayed statistical values are raw p-values; expression selection uses BH-adjusted *q* < 0.05.

### SLK BAM panel E

The four indexed genomic BAMs supplied on GRC are read with pysam in its existing
Docker environment. All pass samtools quickcheck, sample read-group validation,
and the GRCh38 chr10 length check. The notebook contains the extraction code;
configure `SLK_BAM_DIR` and the four `SLK_BAM_PATHS` in the first code cell.
The hypothetical default directory is `/path/to/figure3E_bams`; it can also be
overridden with `CHRONO_SLK_BAM_DIR`. `SLK_USE_SOURCE_CACHE=False` requires BAMs.
When BAMs are unavailable, it checks SHA-256 hashes, parameters, representative
phases/PSI, and junction counts before loading the portable
`Figure3E_SLK_bam_{coverage,junctions,qc,exons,samples}.tsv` source bundle plus
`Figure3E_SLK_bam_provenance.json` from the analysis table directory.

Coverage uses primary, QC-pass, properly paired, uniquely mapped (NH=1),
unclipped 48-nt reads. It is raw per-base aligned-read depth with independent
track y scales, without smoothing or library-size normalization; overlapping
mates are counted separately. All 12 observed event junction counts match rMATS.
The BAM header/index hashes and file sizes/timestamps are recorded; full multi-GB
BAM hashes were not computed. These BAMs were supplied separately, and their
connection to the checked-in pipeline is not independently established.

The individual E PNG/PDF is 157.5 × 78.75 mm; its two exon models represent event
inclusion/skipping paths, not full transcripts. Introns are compressed fivefold.
The original BAMs stay on GRC. Only derived source tables are needed to render
this panel locally. The former blank individual export is retained in the figure
archive. The current Figure 3 assembly includes the real E panel.
