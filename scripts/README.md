# Supporting code

The six publication notebooks embed their analysis functions; these helpers are provided for upstream preparation and standalone reuse. They are not an instruction to regenerate authoritative figure notebooks from old builders.

| Files | Purpose |
| --- | --- |
| `aggregate_dapars_chromosome_results.py` | Pipeline dependency: aggregate DaPars2 chromosome outputs with a manifest |
| `build_enst_to_gene_name.py` | Build a versionless transcript-to-gene lookup from the supplied GENCODE v49 GTF; run from the repository root |
| `CHIRAL_dual_penalty.R`, `CHIRAL_dual_penalty_benchmark_fixed.R` | Core CHIRAL implementations for unordered-sample phase inference and benchmark support; Figure 1 embeds its own authoritative engine |
| `kegg_enrichment_analysis.R`, `reactome_enrichment_analysis.R` | Standalone pathway-enrichment preparation from supplied significant-gene tables |
| `rbp_rna_processing_links.py`, `generate_rbp_processing_source_tables.py`, `rbp_rna_processing_deep_dive.py`, `rbp_rna_processing_site_modeling.py` | Figure 5 supporting RBP/APA/SE analysis and source-table modules, also embedded in the notebook |
| `run_chip_background_sensitivity.py`, `run_ppi_panel_l_background_comparison.py` | Prepare Figure 2 promoter-assessed ChIP and BioGRID background inputs; their historical target comparisons are not the final publication inference |
| `kronos_cosinor.py` | Kronos workbook parsing and cosinor support; Figure 6 embeds the required runtime code |
| `DaPars2/` | Third-party source and original license/README; test datasets and nested Git history excluded |
| `ecotyper/` | Third-party code, configuration and license/README; models, utility reference data, examples and outputs excluded |

Install the required external data/resources before using these helpers. Source code was copied from the current project working tree; this is not a new calibration of the scientific methods. Historical notebook builders, exploratory scripts, plotting assemblers, ad hoc audits and site-specific BAM transfer scripts are intentionally omitted. Direct BAM coverage extraction is already embedded in the publication notebooks.

The ChIP/PPI preparation helpers require the historical source inputs they declare, including saved `method_A_presence` gene lists. Use them to rebuild upstream support tables only after validating the intended inputs; notebook 02 recomputes the current tumor-specific targets and final statistics from the authoritative regression tables.
