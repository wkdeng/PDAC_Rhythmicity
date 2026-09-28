# Figure 5 and Supplementary Figures S14–S15: source map

Audited 2026-09-22 against the three supplied screenshots, current source notebook, plotting script, saved figures, and source tables. This is a provenance audit only; no analyses or figures were regenerated. Screenshot text is reference material, not additional instructions.

The supplementary screenshots still carry the old numbers **S12** and **S13**. In the current requested numbering, these map to **S14** and **S15**, respectively. These RBP figures are distinct from the APA/miRNA S12–S13 mapped in `figure4_s10_s13_provenance.md`.

Cell numbers below are **one-based positions including Markdown cells**, not execution counts or embedded `# Cell` comments. Paths below are relative to the project root.

## Exact saved composites and renderer

| Requested figure | Matching saved PNG/PDF stem under `data/figures/rbp_rna_processing_links/` | Renderer |
|---|---|---|
| Figure 5A–N | `rbp_processing_tcga_primary_story` | `assemble_story(TCGA_CONFIG)` |
| Figure S14A–B (old S12) | `rbp_processing_figure_s12` | `assemble_reorganized_supplement`, Reactome specification |
| Figure S15A–H (old S13) | `rbp_processing_figure_s13` | `assemble_reorganized_supplement`, CPTAC specification |

All three saved images visually match the supplied content and panel arrangement. The common renderer is `scripts/assemble_rbp_processing_story_figures.py`: `story_panel_specs()` begins at line 928, main assembly at line 1104, and `REORGANIZED_SUPPLEMENTS` at line 1182. Existing separated panels are in `data/figures/rbp_rna_processing_links/separated_panels/`, named `<composite_stem>_panel_<lowercase_letter>.png/.pdf`.

## Source notebook locations

The direct source is `3.4.3_apa_se_rbp_rhythmic_connection.ipynb`. It orchestrates imported analysis modules; it does **not** embed all implementation code.

| Notebook position | Embedded comment | Role |
|---|---|---|
| 3 | Cell 1 | Paths, imports of the four analysis modules and assembly module, display helpers |
| 5 | Cell 2 | `generate_current_source_tables(..., n_lag_permutations=500, use_clip_cache=True)` |
| 7 | Cell 3 | Deep-dive permutation analysis (`n_permutations=2000`), site availability, and sample-level models (`max_models=48`, all eligible events) |
| 9 | Cell 4 | Source-table audit |
| 11 | Cell 5 | Builds TCGA/CPTAC stories, combined Reactome/CPTAC supplement, and reorganized old S12/S13; exports separated panels |
| 13 | Cell 6 | Displays the exact three requested composites and the two parent supplement composites |
| 15 | Cell 7 | Displays Figure 5A–N separated panels |
| 17 | Cell 8 | Displays old S12A–B and S13A–H separated panels |
| 19 | Cell 9 | Output validation |

## Figure 5 panel-to-source mapping

All TSV filenames in the following tables are under `data/rbp_rna_processing_links/`. Assembly creates cohort-specific aliases; their internal letters do **not** always equal the displayed panel letters.

| Panel | Content | Plot function | Canonical data / assembly alias |
|---|---|---|---|
| A | Rhythmic RBP counts across TCGA, GTEx, CPTAC | `plot_rhythmic_landscape` | `rbp_rhythmicity_summary.tsv`; aliases `rbp_{tcga_primary,cptac_parallel}_panel_a_rhythmicity.tsv` |
| B | TCGA tumor RBP phase density | `plot_rbp_phase_density` | `rbp_phase_circular_density.tsv`, `rbp_phase_polarization_summary.tsv`, `rbp_expression_rhythmicity.tsv`; `TCGA_Tumor` |
| C | GTEx normal RBP phase density | Same | Same three tables; `GTEx_Normal` |
| D | CPTAC tumor RBP phase density | Same | Same three tables; `CPTAC_Tumor` |
| E | Tumor–normal rhythmicity overlap classes | `plot_overlap` | `rbp_tumor_normal_rhythmic_overlap.tsv`; combines both `*_panel_b_overlap.tsv` aliases |
| F | TCGA APA CLIP-supported RBP enrichment | `plot_fdr_bars` | `rbp_apa_clip_enrichment.tsv` → `rbp_tcga_primary_panel_c_apa_enrichment.tsv` |
| G | TCGA SE CLIP-supported RBP enrichment | `plot_fdr_bars` | `rbp_se_clip_enrichment.tsv` → `rbp_tcga_primary_panel_d_se_enrichment.tsv` |
| H | TCGA target-minus-background rhythmicity | `plot_target_contrast` | `rbp_target_background_contrast.tsv` → `rbp_tcga_primary_panel_e_target_contrast.tsv` |
| I | TCGA phase-lag coherence | `plot_phase_coherence` | `rbp_phase_lag_coherence.tsv` → `rbp_tcga_primary_panel_f_phase_coherence.tsv` |
| J | TCGA sample-level PDUI/PSI usage models | `plot_sample_models` | `rbp_sample_level_usage_model_summary.tsv` → `rbp_tcga_primary_panel_j_sample_model.tsv` |
| K | TCGA module evidence matrix | `build_evidence_matrix` + `plot_evidence_matrix` | Candidate priority, site availability, and sample-model tables; writes `rbp_tcga_primary_panel_k_evidence_matrix.tsv` |
| L | TCGA RBP→event→mRNA phase modules | `plot_phase_modules` | `rbp_deep_candidate_priority.tsv` → `rbp_tcga_primary_panel_g_candidate_priority.tsv`, plus K's evidence-selected modules |
| M | TCGA SE positional CLIP support | `plot_se_region` | `rbp_se_region_bias.tsv` → `rbp_tcga_primary_panel_m_se_region.tsv`; loader also supports older j/h aliases |
| N | TCGA APA regulatory-site availability | `plot_site_availability` | `rbp_apa_regulatory_site_bias.tsv` → `rbp_tcga_primary_panel_n_site_availability.tsv`; loader also supports older k/i aliases |

Identity checks agree with the screenshot: **48/160 TCGA, 103/160 GTEx, 76/160 CPTAC** rhythmic RBPs; TCGA overlap categories are **17 tumor-only, 31 shared, 72 normal-only, 40 neither**. The phase summaries contain n=48/103/76 and axial R2 approximately 0.345/0.354/0.255.

K is an unweighted six-criterion evidence matrix, not a single fitted score. Its columns encode CLIP enrichment, target/background contrast, phase-lag coherence, cross-cohort support, APA/SE feature support, and sample-level model support. The renderer selects up to seven APA and five SE modules and ranks them by evidence and event support.

## Supplementary Figure S14 (old S12)

Both panels use `plot_reactome_overlap()` at line 786 and `rbp_overlap_reactome_enrichment_top.tsv`.

| New panel | Legacy panel | Data subset |
|---|---|---|
| S14A | S12A / parent Reactome-CPTAC supplement A | `comparison == "TCGA_vs_GTEx"` |
| S14B | S12B / parent Reactome-CPTAC supplement B | `comparison == "CPTAC_vs_GTEx"` |

The reorganized specification retains these two parent panels as a two-row figure. For each comparison, the plot selects the ten lowest nominal p-values (odds ratio breaks ties), displays −log10(nominal p), colors overlap classes, and scales point size by overlap count. This is **descriptive pathway context**: no term in the saved top table passes FDR < 0.05 (minimum FDR 0.512 for TCGA and 0.397 for CPTAC). Do not relabel the x-axis as adjusted significance.

## Supplementary Figure S15 (old S13)

This is a curated subset of `rbp_processing_cptac_parallel_supplement`, not the entire parent figure. The source mapping is explicitly encoded in `REORGANIZED_SUPPLEMENTS`.

| New panel | Old S13 panel | CPTAC parent panel | Content / data |
|---|---|---|---|
| S15A | A | D | CPTAC RBP phase; same cohort and tables as Figure 5D |
| S15B | B | F | APA CLIP enrichment; `rbp_cptac_parallel_panel_c_apa_enrichment.tsv` |
| S15C | C | H | Target/background contrast; `rbp_cptac_parallel_panel_e_target_contrast.tsv` |
| S15D | D | I | Phase-lag coherence; `rbp_cptac_parallel_panel_f_phase_coherence.tsv` |
| S15E | E | J | Sample-level usage models; `rbp_cptac_parallel_panel_j_sample_model.tsv` |
| S15F | F | K | Evidence matrix; `rbp_cptac_parallel_panel_k_evidence_matrix.tsv` |
| S15G | G | L | Phase modules; `rbp_cptac_parallel_panel_g_candidate_priority.tsv` plus evidence-selected modules |
| S15H | H | N | Site availability; `rbp_cptac_parallel_panel_n_site_availability.tsv` |

The parent-letter sequence is **D, F, H, I, J, K, L, N**. Figure 5D and S15A intentionally show the same CPTAC phase distribution. Despite the APA-centered figure title, the evidence matrix and phase-module panels retain some SE rows from their parent specifications.

## Other code found

`scripts/build_fig05_publication_notebooks.py` is a separate draft generator targeting old `notebooks_publish/fig05_rbp_processing_modules.ipynb`, `fig_s12_rbp_overlap_reactome_context.ipynb`, and `fig_s13_cptac_rbp_processing_validation.ipynb`. Those three target notebooks were absent at audit time. The generator uses a different reconstruction route and is not the renderer of these screenshots. New work should use the requested `notebook_publication/` location and S14/S15 numbering.

## Upstream analysis and portability boundary

The source-table chain is:

1. `scripts/generate_rbp_processing_source_tables.py::generate_current_source_tables` (line 933) orchestrates the RBP catalogue/expression, APA/SE inputs, CLIP links, phase summaries, overlap/Reactome, contrasts, lag coherence, and candidate tables. Its implementation relies on `scripts/rbp_rna_processing_links.py`.
2. `scripts/rbp_rna_processing_deep_dive.py::run_deep_dive` adds positional/rank bias, phase modules, and candidate evidence from those saved tables.
3. `scripts/rbp_rna_processing_site_modeling.py::run_site_availability` and `run_sample_level_models` add the regulatory-site proxy and usage-model results.
4. `scripts/assemble_rbp_processing_story_figures.py` refreshes cohort-specific aliases, constructs the evidence matrix, and renders the composites and components.

Important upstream inputs include:

| Input | Role |
|---|---|
| `data/Merged/regression_windowed/all_genes_regression_periodic_{TCGA-PAAD_Tumor,CPTAC-3_Tumor,GTEx_Normal}.txt` | Existing expression rhythmicity estimates used to classify/phase RBPs |
| `data/polyAPA/comparison/` | Existing APA event/host membership and phase-lag analysis results |
| `data/AS/<cohort>/filtered_events/` and associated rhythmicity outputs | Skipped-exon events, PSI, and rhythmicity |
| `data/Ref/clip_rbp/` | External CLIP references discovered locally, including POSTAR3/ENCODE eCLIP sources |
| `data/rbp_rna_processing_links/clip_supported_rbp_target_links_consensus.tsv` | Reusable consensus CLIP-link cache; notebook permits cache reuse |
| `data/Ref/hg38_extracted_3UTR.bed` | UTR intervals/rank/overlap context |
| `data/apa_mirna_rhythmic_connection/apa_mirna_gene_triads.tsv` | Existing TargetScan-derived miRNA site burden for the regulatory-site proxy |
| `data/Merged/Group_TPM_<cohort>.csv`, APA usage and SE usage matrices | Sample-level RBP expression versus PDUI/PSI models |
| `../SingleCellRhythmicity/data/Ref/Reactome_anno_HS.txt` | Reactome annotation **outside this repository**, explicitly loaded by `build_reactome_overlap_tables` at line 366 |

The site-availability axis represents a proxy based on miRNA burden, broad UTR CLIP burden, and annotated UTR rank, not direct measured site accessibility. Preserve this qualification in the later publication notebook.

`process.smk` is upstream of RNA/rMATS/DaPars processing but does not run this RBP/Reactome analysis or the figure assembly. A notebook that starts from its processed outputs will need the intervening rhythmicity/APA/SE analysis code and explicit external references. Reading only the existing panel TSVs would reproduce the plots but would not satisfy the requested complete analysis-from-processed-data contract.

`scripts/build_3_4_3_rbp_section_notebook.py` is a notebook template/wrapper. Its current source lacks the reorganized-supplement calls present in the actual notebook, so the actual notebook and assembly script are the authoritative sources for these screenshots.

Only this provenance Markdown was added for this request. Existing notebooks, data, figures, and unrelated working-tree files were left unchanged.
