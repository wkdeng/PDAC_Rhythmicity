# =============================================================================
# ChronoTherapy rMATS/DaPars Pipeline (Pre-aligned BAMs)
# Designed for OpenPBS cluster with Singularity containers
#
# Workflow: BAM -> bedGraph -> rMATS (prep + post) + DaPars
# For cohorts with existing BAM files (e.g., TCGA genomic)
# =============================================================================
shell.prefix("set -eo pipefail; ")
import os
import glob

configfile: "config_rmats_singularity.yaml"

# --- Configuration ---
# Resolve all paths to absolute so they work inside Singularity containers
# and when PBS jobs run from different working directories
COHORT = config["cohort"]
BAM_DIR = os.path.abspath(config["bam_dir"])
GTF_PATH = os.path.abspath(config["gtf"])
UTR3_BED = os.path.abspath(config["utr3_bed"])
GENOME_SIZE_FILE = os.path.abspath(config["genome_size_file"])
SINGULARITY_SIF = config["singularity_sif"]  # already absolute

RMATS_PYTHON = config["rmats_python"]    # container-internal path
RMATS_PATH = config["rmats_path"]        # container-internal path
READ_LENGTH = config["read_length"]
READ_TYPE = config["read_type"]

DAPARS_SCRIPT = config["dapars_script"]  # container-internal path
DAPARS_COVERAGE = config["dapars_coverage"]
DAPARS_CHR_FILE = os.path.abspath(config.get("dapars_chr_file", "hg38_chroms.txt"))
DAPARS_ALLOW_MISSING_CHROMS = bool(config.get("dapars_allow_missing_chroms", False))

FEATURECOUNTS_PATH = config.get("featurecounts_path", "/root/subread-2.0.2-Linux-x86_64/bin/featureCounts")  # container-internal path

ECOTYPER_DIR = os.path.abspath(config.get("ecotyper_dir", "../scripts/ecotyper"))
ECOTYPER_DISCOVERY = config.get("ecotyper_discovery", "Carcinoma")
ECOTYPER_THREADS = config.get("ecotyper_threads", 10)

# Project root (parent of pipeline/) — bind this so container sees both
# pipeline/ and data/ directories
PROJECT_DIR = os.path.abspath("..")
DATA_DIR = os.path.abspath(f"../data/{COHORT}")
LOGS_DIR = os.path.abspath("../data/logs")
POLYAPA_DIR = os.path.abspath(config.get("polyapa_dir", f"../data/polyAPA/{COHORT}"))
DAPARS_CHROMOSOME_DIR = os.path.abspath(
    config.get("dapars_chromosome_results_dir", os.path.join(POLYAPA_DIR, "result_by_chromosome"))
)
DAPARS_RESULTS_DIR = os.path.abspath(
    config.get("dapars_results_dir", os.path.join(POLYAPA_DIR, "results"))
)
DAPARS_AGGREGATOR = os.path.abspath(
    config.get("dapars_aggregator", "../scripts/aggregate_dapars_chromosome_results.py")
)
DAPARS_ALLOW_MISSING_FLAG = "--allow-missing-chroms" if DAPARS_ALLOW_MISSING_CHROMS else ""

print(f"PROJECT_DIR: {PROJECT_DIR}")
print(f"DATA_DIR: {DATA_DIR}")
print(f"POLYAPA_DIR: {POLYAPA_DIR}")
print(f"DAPARS_CHROMOSOME_DIR: {DAPARS_CHROMOSOME_DIR}")
print(f"DAPARS_RESULTS_DIR: {DAPARS_RESULTS_DIR}")
print(f"BAM_DIR: {BAM_DIR}")

# --- Auto-detect samples from BAM files ---
# Known BAM file suffixes (longest first for greedy matching)
BAM_SUFFIXES = [
    ".Aligned.sortedByCoord.out.patched.md.bam",
    ".Aligned.sortedByCoord.out.bam",
    ".bam",
]


def get_samples_and_bam_paths():
    """Auto-detect samples from BAM files, handling various naming conventions."""
    sample_to_bam = {}
    for f in sorted(glob.glob(f"{BAM_DIR}/*.bam")):
        basename = os.path.basename(f)
        sample_name = basename
        for suffix in BAM_SUFFIXES:
            if basename.endswith(suffix):
                sample_name = basename[: -len(suffix)]
                break
        sample_to_bam[sample_name] = f
    return sample_to_bam


SAMPLE_TO_BAM = get_samples_and_bam_paths()
SAMPLES = sorted(SAMPLE_TO_BAM.keys())

print(f"Detected {len(SAMPLES)} samples from {BAM_DIR}")
if len(SAMPLES) > 0:
    print(f"First 5 samples: {SAMPLES[:5]}")


def get_bam(wildcards):
    """Look up the actual BAM file path for a sample."""
    return SAMPLE_TO_BAM[wildcards.sample]


# --- Local rules (run on head node, not submitted to PBS) ---
localrules: all, create_bam_list, create_combined_bam_list, merge_sequencing_depth, create_dapars_config, convert_counts_to_tpm


# =============================================================================
# Target rule
# =============================================================================

rule all:
    input:
        rmats_prep=expand(f"{DATA_DIR}/rmats/prep/{{sample}}/prep_completed.txt", sample=SAMPLES),
        bedgraph=expand(f"{DATA_DIR}/bedgraph/{{sample}}.bedgraph", sample=SAMPLES),
        rmats_post=f"{DATA_DIR}/rmats/post/post_completed.txt",
        dapars=expand(
            f"{DAPARS_RESULTS_DIR}/{{sample}}/{{sample}}_result_All_Prediction_Results.txt",
            sample=SAMPLES,
        ),
        ecotyper=f"{DATA_DIR}/ecotyper/ecotyper_completed.txt",


# =============================================================================
# BAM to BedGraph Conversion
# =============================================================================

rule convert_bedgraph:
    input:
        bam=get_bam,
    output:
        bedgraph=f"{DATA_DIR}/bedgraph/{{sample}}.bedgraph",
    log:
        f"{LOGS_DIR}/{COHORT}/bedgraph/{{sample}}.log",
    params:
        genome_size=GENOME_SIZE_FILE,
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        bam_dir=BAM_DIR,
    threads: 4
    resources:
        walltime="24:00:00",
        mem_mb=16384,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p $(dirname {output.bedgraph})
        mkdir -p $(dirname {log})

        echo "=== convert_bedgraph: {wildcards.sample} ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "Input BAM: {input.bam}"
        echo "Output bedgraph: {output.bedgraph}"
        echo "Genome size file: {params.genome_size}"
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir} and {params.bam_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            --bind {params.bam_dir}:{params.bam_dir} \
            {params.sif} \
            genomeCoverageBed -bg -ibam {input.bam} -g {params.genome_size} -split \
                > {output.bedgraph} 2> {log}

        echo "=== convert_bedgraph: {wildcards.sample} DONE ==="
        """


# =============================================================================
# rMATS - Alternative Splicing Analysis
# =============================================================================

rule create_bam_list:
    input:
        bam=get_bam,
    output:
        bam_list=f"{DATA_DIR}/rmats/bam_lists/per_sample/{{sample}}.txt",
    shell:
        """
        mkdir -p $(dirname {output.bam_list})
        echo "{input.bam}" > {output.bam_list}
        """


rule rmats_prep:
    input:
        bam_list=f"{DATA_DIR}/rmats/bam_lists/per_sample/{{sample}}.txt",
    output:
        completed=f"{DATA_DIR}/rmats/prep/{{sample}}/prep_completed.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/rmats/{{sample}}.log",
    params:
        gtf=GTF_PATH,
        rmats=RMATS_PATH,
        rmats_python=RMATS_PYTHON,
        read_length=READ_LENGTH,
        read_type=READ_TYPE,
        outdir=f"{DATA_DIR}/rmats/prep/{{sample}}",
        tmpdir=f"{DATA_DIR}/rmats/tmp/{{sample}}",
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        bam_dir=BAM_DIR,
    threads: 12
    resources:
        walltime="24:00:00",
        mem_mb=40960,
    shell:
        """
        set +u
        module load singularity
        set -eo pipefail

        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        rm -f {params.tmpdir}/*.rmats
        mkdir -p $(dirname {log})

        echo "=== rmats_prep: {wildcards.sample} ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "BAM list: {input.bam_list}"
        echo "BAM list content: $(cat {input.bam_list})"
        echo "GTF: {params.gtf}"
        echo "Output dir: {params.outdir}"
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir} and {params.bam_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            --bind {params.bam_dir}:{params.bam_dir} \
            {params.sif} \
            {params.rmats_python} {params.rmats} \
                --b1 {input.bam_list} \
                --gtf {params.gtf} \
                -t {params.read_type} \
                --readLength {params.read_length} \
                --nthread {threads} \
                --od {params.outdir} \
                --tmp {params.tmpdir} \
                --task prep \
                --variable-read-length \
                --individual-counts \
                > {log} 2>&1

        echo "=== rmats_prep: {wildcards.sample} DONE ==="
        echo "Completed at $(date)" > {output.completed}
        """


rule create_combined_bam_list:
    input:
        bams=ancient(list(SAMPLE_TO_BAM.values())),
    output:
        combined_list=f"{DATA_DIR}/rmats/bam_lists/all_samples.txt",
    shell:
        """
        mkdir -p $(dirname {output.combined_list})
        echo "{input.bams}" | tr ' ' ',' > {output.combined_list}
        """


rule rmats_post:
    input:
        combined_list=f"{DATA_DIR}/rmats/bam_lists/all_samples.txt",
        prep_completed=expand(f"{DATA_DIR}/rmats/prep/{{sample}}/prep_completed.txt", sample=SAMPLES),
    output:
        completed=f"{DATA_DIR}/rmats/post/post_completed.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/rmats/post.log",
    params:
        gtf=GTF_PATH,
        rmats=RMATS_PATH,
        rmats_python=RMATS_PYTHON,
        read_length=READ_LENGTH,
        read_type=READ_TYPE,
        outdir=f"{DATA_DIR}/rmats/post",
        tmpdir=f"{DATA_DIR}/rmats/tmp",
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        bam_dir=BAM_DIR,
    threads: 48
    resources:
        walltime="24:00:00",
        mem_mb=81920,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        mkdir -p $(dirname {log})

        echo "=== rmats_post ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "Combined BAM list: {input.combined_list}"
        echo "GTF: {params.gtf}"
        echo "Output dir: {params.outdir}"
        echo "Threads: {threads}"
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir} and {params.bam_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            --bind {params.bam_dir}:{params.bam_dir} \
            {params.sif} \
            {params.rmats_python} {params.rmats} \
                --b1 {input.combined_list} \
                --gtf {params.gtf} \
                -t {params.read_type} \
                --readLength {params.read_length} \
                --nthread {threads} \
                --od {params.outdir} \
                --tmp {params.tmpdir} \
                --task post \
                --statoff \
                --variable-read-length \
                --individual-counts \
                > {log} 2>&1

        echo "=== rmats_post DONE ==="
        echo "Completed at $(date)" > {output.completed}
        """


# =============================================================================
# DaPars - Alternative Polyadenylation Analysis
# =============================================================================

rule compute_sequencing_depth:
    input:
        bam=get_bam,
    output:
        depth=f"{DATA_DIR}/dapars/depth/{{sample}}_depth.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/dapars/depth/{{sample}}.log",
    params:
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        bam_dir=BAM_DIR,
    threads: 1
    resources:
        walltime="04:00:00",
        mem_mb=4096,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p $(dirname {output.depth})
        mkdir -p $(dirname {log})

        COUNT=$(singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            --bind {params.bam_dir}:{params.bam_dir} \
            {params.sif} \
            samtools view -c -F 260 {input.bam} 2> {log})

        printf "{wildcards.sample}\\t%s\\n" "$COUNT" > {output.depth}
        """


rule merge_sequencing_depth:
    input:
        depths=expand(f"{DATA_DIR}/dapars/depth/{{sample}}_depth.txt", sample=SAMPLES),
    output:
        merged=f"{DATA_DIR}/dapars/sequencing_depth.txt",
    shell:
        """
        cat {input.depths} > {output.merged}
        """


rule create_dapars_config:
    input:
        bedgraph=f"{DATA_DIR}/bedgraph/{{sample}}.bedgraph",
        utr3=UTR3_BED,
        sequencing_depth=f"{DATA_DIR}/dapars/sequencing_depth.txt",
    output:
        config=f"{DATA_DIR}/dapars/config/{{sample}}_config.txt",
    params:
        outdir=f"{DAPARS_CHROMOSOME_DIR}/{{sample}}",
        coverage=DAPARS_COVERAGE,
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.config})

        cat > {output.config} <<EOF
Annotated_3UTR={input.utr3}
Aligned_Wig_files={input.bedgraph}
Output_directory={params.outdir}
Output_result_file={wildcards.sample}_result
Coverage_threshold={params.coverage}
Num_Threads={threads}
sequencing_depth_file={input.sequencing_depth}
EOF
        """


rule run_dapars:
    input:
        config=f"{DATA_DIR}/dapars/config/{{sample}}_config.txt",
        bedgraph=f"{DATA_DIR}/bedgraph/{{sample}}.bedgraph",
        chr_file=DAPARS_CHR_FILE,
    output:
        completed=f"{DAPARS_CHROMOSOME_DIR}/{{sample}}/dapars_completed.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/dapars/{{sample}}.log",
    params:
        dapars_script=DAPARS_SCRIPT,
        outdir=f"{DAPARS_CHROMOSOME_DIR}/{{sample}}",
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
    threads: 4
    resources:
        walltime="24:00:00",
        mem_mb=16384,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        echo "=== run_dapars: {wildcards.sample} ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "Config file: {input.config}"
        echo "Config content:"
        cat {input.config}
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            {params.sif} \
            python {params.dapars_script} {input.config} {input.chr_file} >> {log} 2>&1

        echo "=== run_dapars: {wildcards.sample} DONE ==="
        echo "Completed at $(date)" > {output.completed}
        """


rule aggregate_dapars_chromosomes:
    input:
        completed=f"{DAPARS_CHROMOSOME_DIR}/{{sample}}/dapars_completed.txt",
        chr_file=DAPARS_CHR_FILE,
    output:
        result=f"{DAPARS_RESULTS_DIR}/{{sample}}/{{sample}}_result_All_Prediction_Results.txt",
        manifest=f"{DAPARS_RESULTS_DIR}/{{sample}}/{{sample}}_aggregation_manifest.tsv",
        completed=f"{DAPARS_RESULTS_DIR}/{{sample}}/dapars_completed.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/dapars/aggregate/{{sample}}.log",
    params:
        aggregator=DAPARS_AGGREGATOR,
        chrom_root=DAPARS_CHROMOSOME_DIR,
        output_root=DAPARS_RESULTS_DIR,
        allow_missing=DAPARS_ALLOW_MISSING_FLAG,
    threads: 1
    resources:
        walltime="02:00:00",
        mem_mb=4096,
    shell:
        """
        mkdir -p $(dirname {log})

        echo "=== aggregate_dapars_chromosomes: {wildcards.sample} ===" > {log}
        echo "Host: $(hostname)" >> {log}
        echo "Working dir: $(pwd)" >> {log}
        echo "Chromosome root: {params.chrom_root}" >> {log}
        echo "Output root: {params.output_root}" >> {log}
        echo "Aggregator: {params.aggregator}" >> {log}

        python {params.aggregator} \
            --chrom-root {params.chrom_root} \
            --output-root {params.output_root} \
            --chromosomes {input.chr_file} \
            --sample {wildcards.sample} \
            --force \
            {params.allow_missing} \
            >> {log} 2>&1

        echo "=== aggregate_dapars_chromosomes: {wildcards.sample} DONE ===" >> {log}
        """


# =============================================================================
# EcoTyper - Cell State and Ecotype Recovery
# =============================================================================

rule featurecounts:
    input:
        bam=ancient(list(SAMPLE_TO_BAM.values())),
    output:
        counts=f"{DATA_DIR}/featurecounts/gene_counts.tsv",
        summary=f"{DATA_DIR}/featurecounts/gene_counts.tsv.summary",
    log:
        f"{LOGS_DIR}/{COHORT}/featurecounts/gene_counts.log",
    params:
        gtf=GTF_PATH,
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        bam_dir=BAM_DIR,
        featurecounts=FEATURECOUNTS_PATH
    threads: 8
    resources:
        walltime="04:00:00",
        mem_mb=32768,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p $(dirname {output.counts})
        mkdir -p $(dirname {log})

        echo "=== featurecounts ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "GTF: {params.gtf}"
        echo "Output: {output.counts}"
        echo "Number of BAMs: $(echo {input.bam} | wc -w)"
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir} and {params.bam_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            --bind {params.bam_dir}:{params.bam_dir} \
            {params.sif} \
            {params.featurecounts} \
                -a {params.gtf} \
                -o {output.counts} \
                -T {threads} \
                -p --countReadPairs \
                {input.bam} \
                > {log} 2>&1

        echo "=== featurecounts DONE ==="
        """


rule convert_counts_to_tpm:
    input:
        counts=f"{DATA_DIR}/featurecounts/gene_counts.tsv",
    output:
        tpm=f"{DATA_DIR}/featurecounts/gene_tpm.tsv",
    run:
        import pandas as pd
        import numpy as np

        # featureCounts output: first line is a comment (#), then header row
        df = pd.read_csv(input.counts, sep="\t", comment="#")

        # Columns: Geneid, Chr, Start, End, Strand, Length, then sample columns
        gene_ids = df["Geneid"]
        lengths = df["Length"].values.astype(float)
        sample_cols = df.columns[6:]  # all columns after metadata

        counts = df[sample_cols].values.astype(float)

        # TPM: (counts / length) normalized per sample
        rate = counts / lengths[:, np.newaxis]
        tpm = rate / rate.sum(axis=0)[np.newaxis, :] * 1e6

        # Build output with "Gene" as first column, clean sample names
        out = pd.DataFrame(tpm, columns=sample_cols)
        out.insert(0, "Gene", gene_ids.values)

        # Clean sample names: strip path prefixes and BAM suffixes
        rename = {}
        for col in sample_cols:
            name = os.path.basename(col)
            for suffix in BAM_SUFFIXES:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            rename[col] = name
        out.rename(columns=rename, inplace=True)

        out.to_csv(output.tpm, sep="\t", index=False)


rule ecotyper_recovery:
    input:
        tpm=f"{DATA_DIR}/featurecounts/gene_tpm.tsv",
    output:
        completed=f"{DATA_DIR}/ecotyper/ecotyper_completed.txt",
    log:
        f"{LOGS_DIR}/{COHORT}/ecotyper/recovery.log",
    params:
        ecotyper_dir=ECOTYPER_DIR,
        discovery=ECOTYPER_DISCOVERY,
        outdir=f"{DATA_DIR}/ecotyper",
        sif=SINGULARITY_SIF,
        project_dir=PROJECT_DIR,
        threads=ECOTYPER_THREADS,
    threads: 10
    resources:
        walltime="24:00:00",
        mem_mb=65536,
    shell:
        """
        module load singularity
        set -eo pipefail

        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        echo "=== ecotyper_recovery ==="
        echo "Host: $(hostname)"
        echo "Working dir: $(pwd)"
        echo "EcoTyper dir: {params.ecotyper_dir}"
        echo "Discovery: {params.discovery}"
        echo "Input TPM: {input.tpm}"
        echo "Output dir: {params.outdir}"
        echo "Threads: {params.threads}"
        echo "Singularity SIF: {params.sif}"
        echo "Bind: {params.project_dir}"

        singularity exec \
            --bind {params.project_dir}:{params.project_dir} \
            {params.sif} \
            bash -c "cd {params.ecotyper_dir} && Rscript EcoTyper_recovery_bulk.R \
                -d {params.discovery} \
                -m {input.tpm} \
                -t {params.threads} \
                -o {params.outdir}" \
                > {log} 2>&1

        echo "=== ecotyper_recovery DONE ==="
        echo "Completed at $(date)" > {output.completed}
        """
