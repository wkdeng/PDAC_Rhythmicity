# Publication notebooks

For Figure 4 and Supplementary Figures S10–S13, see [the Figure 4 execution guide](04_figure4_readme.md) and `04_figure4_apa_mirna_coordination.ipynb`. S10D is enabled in the current notebook and uses real indexed BAMs or its checksum-validated BAM-derived source bundle; supply these inputs separately. Historical provenance notes may describe earlier disabled-panel versions.

For Figure 3 and Supplementary Figures 7–9, see [the Figure 3 execution guide](03_figure3_readme.md) and [the self-contained notebook](03_figure3_rhythmic_alternative_splicing.ipynb).

For Figure 2 and Supplementary Figures 4–6, see [the Figure 2 execution guide](02_figure2_readme.md) and `02_figure2_cohort_clock_organization.ipynb`.

## Figure 1 and Supplementary Figures 1–3

Open `01_figure1_chiral_benchmark.ipynb` using the **R (`ir`) kernel** and run all cells.
The notebook contains the CHIRAL engine, expression preprocessing, scoring, benchmark
summaries, sampling-window diagnostics, and plotting code. It does not source another
notebook or project R script at runtime.

The default uses existing prediction caches and regenerates downstream analyses and
figures. Set `CHRONO_REFIT=true` before launching the kernel, or set `RUN_ANALYSIS <- TRUE`
in setup, to rerun the 1,500 primary fits and sampling-window perturbation fits from
processed count/TPM matrices. Full refitting is substantially more expensive than
regenerating figures from predictions.

### Environment and inputs

The intended environment is the project's Docker Jupyter environment. The setup cell
checks R packages. Core dependencies include dplyr, tidyr, purrr, tibble, readr, ggplot2,
patchwork, ragg, scales, edgeR, digest, jsonlite, foreach, and doParallel. IRkernel and
IRdisplay support notebook execution and inline PNGs. Package versions are recorded
in the executed notebook and `publication_session_info.txt`.

Paths are detected from the working directory or its parent. Set `CHRONO_PROJECT_ROOT`
for another repository location. `CHRONO_BENCHMARK_INPUT_DIR` can point to an alternate
benchmark input directory; `CHRONO_BENCHMARK_CACHE_DIR` can point to the legacy cache.
The notebook contains no personal macOS paths.

Required processed inputs under `data/benchmarking/`:

| Cohort | Expression | Metadata |
|---|---|---|
| Mouse tumor | `Testing_Tumor_Mouse_raw_counts.csv` | `Testing_Tumor_Mouse_sample_info.csv` |
| Mouse liver CR-spread | `mouse_liver/cr_spread_6mon.txt` | Parsed from matrix columns |
| Human tumor cell lines | `Testing_Tumor_Cell_raw_counts.csv` | `Testing_Tumor_Cell_sample_info.csv` |
| Human tumor PDO | `Testing_Tumor_PDO_raw_TPM.txt` | `Testing_Tumor_PDO_sample_info.csv` |
| Human normal whole blood | `Training_Normal_WB_raw_TPM.txt` | `Training_Normal_WB_sample_info.csv` |

These are existing processed benchmark inputs. Their derivation from `process.smk`
outputs has not been established. The notebook reports gene-symbol coverage and hashes
the input files; it does not invent an upstream sequencing or genome-build provenance.
The expression and cache files are not embedded in the notebook. Supply them separately
when distributing the public repository.

Cached mode additionally requires these files under `data/benchmarking/1.6.3_chiral/`:

- `chiral_1_6_3_predictions.csv`
- `chiral_1_6_3_scored_predictions.csv` (saved global alignment parameters)
- `point5_window_stress_test/point5_chiral_window_stress_test_predictions.csv`

The saved alignment is retained in cached mode because rounding phases to CSV can
change which effectively tied global rotation is selected. Errors are recomputed from
the raw predictions, and the stored rotation's MedAE is checked against fresh optimization.
Refit mode computes alignment from its new predictions.

### Outputs

- Notebook: `notebook_publication/01_figure1_chiral_benchmark.ipynb`
- Source tables and recomputed benchmark summaries: `data/publication/fig01_chiral_benchmark/`
- Figures: `data/figures/publication/fig01_chiral_benchmark/`, with `figure1/`,
  `supplementary1/`, `supplementary2/`, and `supplementary3/` subfolders.

There are 48 PNG/PDF pairs: seven Figure 1 panels plus assembly; sixteen panels plus
assembly for each of S1 and S2; and five panels plus assembly for S3. S1 and S2 are
presented as whole figures without panel letters in either PNG or PDF. Their auxiliary
component files retain A–P only as filename identifiers and are not displayed separately.
The notebook displays 16 PNGs inline. Other PNGs omit panel labels; other PDFs include them. Export dimensions and minimum text
sizes are recorded in `publication_export_manifest.csv`. Global dimensions, fonts,
and the save helper are defined once in setup. Nominal sizes use 52.5 mm modules;
the R PDF device rounds page dimensions to whole PostScript points (less than
0.36 mm difference).

Figure 1 uses two rows: D–G on the first; H–J on the second. H and I each occupy
52.5 × 69.825 mm, and J occupies 105 × 69.825 mm, with a 1:1:2
width ratio. The H–J row is 33% taller than one module by explicit request; D–G
remain 52.5 mm tall. The assembled figure is nominally 210 × 122.325 mm. Figure 1 titles are reduced by
1 pt (panel titles 8 pt, facet titles 7 pt), and annotations/subtitles are black
and 6 pt by explicit request. Panel J numbers have transparent backgrounds and
automatically use black or white text for contrast against each cell. Other text roles and supplementary fonts are unchanged.
S1 and S2 omit overall figure titles and subtitles.

### Source mapping

The legacy analysis is `1.6.3_benchmarking.ipynb`. Its cell ordinals 34, 35, and 36
(counting Markdown cells) contain the Figure 1, robustness, and rotation generators.
The publication notebook preserves Figure 1D–J; S1 is mouse tumor, S2 mouse liver,
and S3 the five-panel robustness figure. Existing notebooks, data, and figure exports
are preserved.

See `AGENTS.md` in this folder for the persistent publication layout and naming rules.
