# Upstream quantification pipelines

The default launchers run from this directory and target TACC/Slurm with Apptainer. The source workflows retain strict shell handling and explicit per-rule container execution.

| Entry point | Input | Main output |
| --- | --- | --- |
| `run.sh` / `Snakefile` / `config.yaml` | Paired `{sample}_1.fastq`, `{sample}_2.fastq` | FastQC/MultiQC, kallisto, STAR BAMs, featureCounts, rMATS and DaPars |
| `run_rmats.sh` / `Snakefile_rmats.smk` / `config_rmats.yaml` | Existing genomic BAMs | rMATS, bedGraph, DaPars2 chromosome/sample results, expression and EcoTyper recovery targets |
| `submit_singularity.sh` / `Snakefile_rmats_singularity.smk` | Existing genomic BAMs | Existing BSIC/OpenPBS/Singularity variant |
| `../process.smk` | Paired FASTQs | Original host-specific quantification workflow retained as a provenance reference |

## Configure before running

- Replace `/path/to/bioinformatics.sif` and `YOUR_ALLOCATION` in the YAML files. Configure the partition, cohort, raw/BAM folders, read length/type and tool locations inside the container.
- Supply GENCODE v43 references/indexes for the pipeline, the hg38 3′UTR BED and chromosome sizes. Publication coverage panels separately use GENCODE v49; do not assume the references are interchangeable.
- The inherited analysis uses `data/Ref` while pipelines use `data/ref`; similarly `GEO_Batch1` and `geo_batch1` occur. On case-sensitive Linux, provide the paths exactly as configured or explicitly update the configuration. On default macOS filesystems these may refer to the same directory. No reference data are bundled.
- `dapars_script` in `config.yaml` is an inherited container-internal **DaPars v1** path. Supply that software in the image and correct the path as necessary; the included `scripts/DaPars2/src` is a different implementation. rMATS configurations refer to DaPars2 inside the container.
- `scripts/aggregate_dapars_chromosome_results.py` combines DaPars2 chromosome outputs into the sample-level result contract. `hg38_chroms.txt` is a pipeline control file, not a dataset.
- EcoTyper source is included, but its carcinoma discovery models, signature matrices, utility reference tables, CIBERSORTx dependencies and datasets must be supplied separately. See `scripts/ecotyper/README.md` and `data/README.md`.
- Install the compatible Snakemake environment and cluster-generic executor used by the Slurm profile. The profile is invoked from `pipeline/`; its submit/status script paths are relative to that directory.

From the repository root:

```bash
cd pipeline
./run.sh -n
# Or, for pre-aligned BAMs:
./run_rmats.sh -n
```

Review detected sample counts: an empty data scaffold can yield zero samples and is not a successful validation of a cohort. Submit only after a meaningful dry run with populated inputs and configured software. This export did not submit any jobs. Snakemake was unavailable in the export environment, so the dry-run entry points could not be evaluated there.

`process.smk` retains historical host paths. It is included because the publication provenance refers to it; use the configurable `pipeline/` workflows for the supported cluster setup. The workflows do not generate all processed publication inputs; see `data/README.md` for that boundary.
