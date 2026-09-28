# PDAC_Rhythmicity

Publication code for circadian organization and rhythmic RNA processing in pancreatic cancer across TCGA, GTEx, CPTAC, and experimental cell-line cohorts. The six notebooks in `notebook_publication/` are the authoritative figure entry points.

This repository is a code-only export from the current ChronoTherapy working tree. It contains no copied datasets, notebook execution outputs, rendered figures, legacy analysis notebooks, or original Git history. Curated gene sets, parameters, and literature-evidence records embedded in notebook source remain part of the analysis code.

## Layout

```text
notebook_publication/    Six self-contained figure notebooks, execution guides, provenance
pipeline/               RNA-seq and BAM-based rMATS/DaPars workflows and cluster launchers
scripts/                Supporting analysis/preparation helpers and third-party source
process.smk             Original quantification workflow retained for provenance
data/                   Directory scaffold only; see data/README.md for required inputs
PUBLICATION_EXPORT.json Source/export hashes and the scope of copying changes
```

| Notebook | Figures | Kernel |
| --- | --- | --- |
| `01_figure1_chiral_benchmark.ipynb` | Figure 1, S1–S3 | R (`ir`) |
| `02_figure2_cohort_clock_organization.ipynb` | Figure 2, S4–S6 | R (`ir`) |
| `03_figure3_rhythmic_alternative_splicing.ipynb` | Figure 3, S7–S9 | Python 3 |
| `04_figure4_apa_mirna_coordination.ipynb` | Figure 4, S10–S13 | Python 3 |
| `05_figure5_rbp_processing.ipynb` | Figure 5, S14–S15 | Python 3 |
| `06_figure6_cellline_validation.ipynb` | Figure 6, S16–S17 | Python 3 |

## Run the publication notebooks

1. Supply the inputs described in [data/README.md](data/README.md) and the notebook-specific execution guides. The empty directories alone are not sufficient to run the analyses.
2. Use the existing Docker-hosted Jupyter environment, or an environment with the packages below. From the repository root, set `export CHRONO_PROJECT_ROOT="$PWD"` before starting Jupyter. The historical `CHRONO_*` environment variable names remain supported. Project-name fallback paths now use `PDAC_Rhythmicity`.
3. Open a notebook with the indicated kernel and run cells in order. Figures 1–2 require R; Figures 3–6 require Python. They embed their runtime analysis functions and do not execute the omitted legacy notebooks.
4. Generated source tables, figures, and logs go to `data/publication/`, `data/figures/publication/`, and `data/logs/publication/`. These are ignored by Git.

Figure 1 uses supplied prediction caches by default; use `CHRONO_REFIT=true` to refit from the benchmark matrices. Figure 6 needs separately supplied CircadianMultiOmics processed inputs via `CHRONO_CELLLINE_DATA_DIR`, plus the local Kronos workbook and TCGA reference inputs. Figures 3/4 have real BAM coverage panels; supply indexed BAMs or the checksum-validated source bundles accepted by the notebooks.

Python dependencies include numpy, pandas, scipy, matplotlib, statsmodels, scikit-learn, openpyxl, IPython and Jupyter; direct BAM coverage extraction needs pysam. Companion `03_figure3_requirements.txt` and `06_figure6_requirements.txt` retain the original requirements. R dependencies include IRkernel, IRdisplay, dplyr, tidyr, purrr, tibble, readr, data.table, stringr, ggplot2, ggrepel, patchwork, ragg, scales, edgeR, digest, jsonlite, foreach, doParallel and tictoc. Standalone enrichment helpers additionally need their declared Bioconductor packages. These are dependency notes, not a newly validated environment lockfile.

## Generate upstream quantifications

Use [pipeline/README.md](pipeline/README.md). TACC/Slurm with explicit per-rule Apptainer execution is the default; the existing BSIC/OpenPBS/Singularity alternative is retained. Configure container paths, reference paths, allocation and cohort inputs, then dry-run before submitting any jobs.

The workflows produce sequencing quantifications, rMATS outputs and APA results. They do **not** produce every harmonized cohort matrix, CHIRAL result, external annotation or cell-line analysis used by the publication notebooks. Those processed inputs are an explicit distribution requirement; omitting legacy notebooks does not turn the repository into a complete raw-data-to-figure workflow.

## Export and validation

Only publication notebooks and relevant supporting code are included. Notebook source is preserved apart from renaming project path defaults; execution outputs, execution counts, transient execution metadata and attachments are removed. Source/provenance documents may name historical files that are intentionally absent. `PUBLICATION_EXPORT.json` records copied files and checksums; new repository guides are not source copies.

Pipeline configuration now uses a placeholder container location and allocation. The Slurm profile resolves its scripts relative to `pipeline/`, where the launchers run. Data scaffolding preserves inherited directory names; see the case-sensitive Linux path notes in the pipeline guide.

Syntax and export-integrity checks were performed; analysis execution requires the omitted data, and Snakemake dry runs require a configured workflow environment. No remote repository or data download is created by this export. Third-party licenses are preserved in their source folders; a license for the project-authored code has not been selected.
