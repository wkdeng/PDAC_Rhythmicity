import os
import re
from pathlib import Path
import glob

###**--- CONFIG ---**###
KALLISTO_IDX = "data/ref/gencode.v43.transcripts.kallisto.idx"
STAR_IDX = "data/ref/STAR"
GTF_PATH = "data/ref/gencode.v43.basic.annotation.gtf"
BATCH = 'geo_batch1'
RAW_DIR = f"data/{BATCH}/raw"

# rMATS configuration
RMATS_PYTHON = "/data/brutus_data33/wdeng/workspace/software/anaconda/envs/rmats/bin/python"
RMATS_PATH = "/data/brutus_data33/wdeng/workspace/software/rmats-turbo/rmats.py"
READ_LENGTH = 50
READ_TYPE = "paired"

# DaPars configuration
DAPARS_SCRIPT = "scripts/dapars/src/DaPars_main.py"
UTR3_BED = "data/ref/hg38_extracted_3UTR.bed"
GENOME_SIZE_FILE = "data/ref/hg38_chr_size.txt"
DAPARS_COVERAGE = 30
DAPARS_FDR = 0.05
DAPARS_PDUI = 0.5
DAPARS_FOLD_CHANGE = 0.59

# Extract sample information from fastq files
# Pattern: {sample}_1.fastq and {sample}_2.fastq (e.g., SRR30694397_1.fastq)
def get_sample_info():
    """Extract sample names from fastq files"""
    samples = []
    fastq_files = glob.glob(f"{RAW_DIR}/*_1.fastq")

    for f in fastq_files:
        basename = os.path.basename(f)
        # Extract sample name (e.g., SRR30694397 from SRR30694397_1.fastq)
        sample_name = basename.replace("_1.fastq", "")
        samples.append(sample_name)

    return sorted(samples)

SAMPLES = get_sample_info()

# Print samples for verification
print(f"Found {len(SAMPLES)} samples")
if len(SAMPLES) > 0:
    print(f"First 5 samples: {SAMPLES[:5]}")

def get_r1(wildcards):
    """Get R1 fastq file for a sample"""
    return f"{RAW_DIR}/{wildcards.sample}_1.fastq"

def get_r2(wildcards):
    """Get R2 fastq file for a sample"""
    return f"{RAW_DIR}/{wildcards.sample}_2.fastq"

###**--- RULES ---**###

rule all:
    input:
        fastqc_r1=expand(f"data/{BATCH}/fastqc/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        fastqc_r2=expand(f"data/{BATCH}/fastqc/{{sample}}_R2_fastqc.html", sample=SAMPLES),
        kallisto=expand(f"data/{BATCH}/kallisto/{{sample}}/abundance.tsv", sample=SAMPLES),
        star_bam=expand(f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        star_log=expand(f"data/{BATCH}/star/{{sample}}/Log.final.out", sample=SAMPLES),
        featurecounts=f"data/{BATCH}/featurecounts/gene_counts.tsv",
        bedgraph=expand(f"data/{BATCH}/bedgraph/{{sample}}.bedgraph", sample=SAMPLES),
        rmats_prep=expand(f"data/{BATCH}/rmats/prep/{{sample}}/prep_completed.txt", sample=SAMPLES),
        rmats_post=f"data/{BATCH}/rmats/post/post_completed.txt",
        dapars=expand(f"data/{BATCH}/dapars/results/{{sample}}/dapars_completed.txt", sample=SAMPLES),
        multiqc=f"data/{BATCH}/multiqc/multiqc_report.html"

rule fastqc:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        r1_html=f"data/{BATCH}/fastqc/{{sample}}_R1_fastqc.html",
        r1_zip=f"data/{BATCH}/fastqc/{{sample}}_R1_fastqc.zip",
        r2_html=f"data/{BATCH}/fastqc/{{sample}}_R2_fastqc.html",
        r2_zip=f"data/{BATCH}/fastqc/{{sample}}_R2_fastqc.zip"
    log:
        f"data/{BATCH}/logs/fastqc/{{sample}}.log"
    params:
        out_dir=f"data/{BATCH}/fastqc",
        sample="{sample}"
    threads: 10
    shell:
        """
        mkdir -p {params.out_dir}
        # Run fastqc and rename outputs to match expected names
        fastqc -o {params.out_dir} -t {threads} {input.r1} {input.r2} > {log} 2>&1

        # Get the actual output filenames
        r1_base=$(basename {input.r1} .fastq)
        r2_base=$(basename {input.r2} .fastq)

        # Create symlinks with simplified names
        if [ -f {params.out_dir}/${{r1_base}}_fastqc.html ]; then
            ln -sf ${{r1_base}}_fastqc.html {output.r1_html}
            ln -sf ${{r1_base}}_fastqc.zip {output.r1_zip}
        fi
        if [ -f {params.out_dir}/${{r2_base}}_fastqc.html ]; then
            ln -sf ${{r2_base}}_fastqc.html {output.r2_html}
            ln -sf ${{r2_base}}_fastqc.zip {output.r2_zip}
        fi
        """

rule kallisto_quant:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        h5=f"data/{BATCH}/kallisto/{{sample}}/abundance.h5",
        tsv=f"data/{BATCH}/kallisto/{{sample}}/abundance.tsv",
        json=f"data/{BATCH}/kallisto/{{sample}}/run_info.json"
    log:
        f"data/{BATCH}/logs/kallisto/{{sample}}.log"
    params:
        index=KALLISTO_IDX,
        out_dir=f"data/{BATCH}/kallisto/{{sample}}"
    threads: 40
    shell:
        """
        mkdir -p {params.out_dir}
        kallisto quant -i {params.index} -o {params.out_dir} -t {threads} {input.r1} {input.r2} > {log} 2>&1
        """

rule star_align:
    input:
        r1=get_r1,
        r2=get_r2
    output:
        bam=f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam",
        sj=f"data/{BATCH}/star/{{sample}}/SJ.out.tab",
        log_final=f"data/{BATCH}/star/{{sample}}/Log.final.out",
        log_progress=f"data/{BATCH}/star/{{sample}}/Log.progress.out",
        log_out=f"data/{BATCH}/star/{{sample}}/Log.out",
        unmapped_1=f"data/{BATCH}/star/{{sample}}/Unmapped.out.mate1",
        unmapped_2=f"data/{BATCH}/star/{{sample}}/Unmapped.out.mate2"
    log:
        f"data/{BATCH}/logs/star/{{sample}}.log"
    params:
        index=STAR_IDX,
        gtf=GTF_PATH,
        prefix=f"data/{BATCH}/star/{{sample}}/"
    threads: 40
    shell:
        """
        mkdir -p {params.prefix}
        STAR --genomeDir {params.index} \
            --readFilesIn {input.r1} {input.r2} \
            --sjdbGTFfile {params.gtf} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --twopassMode Basic \
            --runThreadN {threads} \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outReadsUnmapped Fastx \
            > {log} 2>&1
        """

rule star_index_bam:
    input:
        f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam.bai"
    log:
        f"data/{BATCH}/logs/star/{{sample}}_index.log"
    threads: 1
    shell:
        """
        samtools index {input} > {log} 2>&1
        """

rule featurecounts:
    input:
        bam=expand(f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bai=expand(f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)
    output:
        counts=f"data/{BATCH}/featurecounts/gene_counts.tsv",
        summary=f"data/{BATCH}/featurecounts/gene_counts.tsv.summary"
    log:
        f"data/{BATCH}/logs/featurecounts/gene_counts.log"
    params:
        gtf=GTF_PATH,
        out_dir=f"data/{BATCH}/featurecounts"
    threads: 8
    shell:
        """
        mkdir -p {params.out_dir}
        featureCounts \
            -a {params.gtf} \
            -o {output.counts} \
            -T {threads} \
            {input.bam} \
            > {log} 2>&1
        """

rule multiqc:
    input:
        fastqc=expand(f"data/{BATCH}/fastqc/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        kallisto=expand(f"data/{BATCH}/kallisto/{{sample}}/abundance.tsv", sample=SAMPLES),
        star_log=expand(f"data/{BATCH}/star/{{sample}}/Log.final.out", sample=SAMPLES),
        featurecounts_summary=f"data/{BATCH}/featurecounts/gene_counts.tsv.summary"
    output:
        f"data/{BATCH}/multiqc/multiqc_report.html"
    log:
        f"data/{BATCH}/logs/multiqc.log"
    params:
        out_dir=f"data/{BATCH}/multiqc",
        search_dir=f"data/{BATCH}"
    shell:
        """
        mkdir -p {params.out_dir}
        multiqc {params.search_dir} -o {params.out_dir} -f > {log} 2>&1
        """

###**--- rMATS and DaPars RULES ---**###

# Create BAM list file for each sample (required for rMATS)
rule create_bam_list:
    input:
        bam=f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        bam_list=f"data/{BATCH}/rmats/bam_lists/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {output.bam_list})
        echo "{input.bam}" > {output.bam_list}
        """

# Convert BAM to bedGraph for DaPars
rule convert_bedgraph:
    input:
        bam=f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        bedgraph=f"data/{BATCH}/bedgraph/{{sample}}.bedgraph"
    params:
        genome_size=GENOME_SIZE_FILE
    threads: 1
    log:
        f"data/{BATCH}/logs/bedgraph/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {output.bedgraph})
        mkdir -p $(dirname {log})
        genomeCoverageBed -bg -ibam {input.bam} -g {params.genome_size} -split > {output.bedgraph} 2> {log}
        """

# Create DaPars configuration file per sample
rule create_dapars_config:
    input:
        bedgraph=f"data/{BATCH}/bedgraph/{{sample}}.bedgraph",
        utr3=UTR3_BED
    output:
        config=f"data/{BATCH}/dapars/config/{{sample}}_config.txt"
    params:
        outdir=f"data/{BATCH}/dapars/results/{{sample}}",
        coverage=DAPARS_COVERAGE,
        fdr=DAPARS_FDR,
        pdui=DAPARS_PDUI,
        fold_change=DAPARS_FOLD_CHANGE
    shell:
        """
        mkdir -p $(dirname {output.config})

        cat > {output.config} <<EOF
Annotated_3UTR={input.utr3}
Group1_Tophat_aligned_Wig={input.bedgraph}
Group2_Tophat_aligned_Wig={input.bedgraph}
Output_directory={params.outdir}
Output_result_file={wildcards.sample}_result
Num_least_in_group1=1
Num_least_in_group2=1
Coverage_cutoff={params.coverage}
FDR_cutoff={params.fdr}
PDUI_cutoff={params.pdui}
Fold_change_cutoff={params.fold_change}
EOF
        """

# Run DaPars alternative polyadenylation analysis per sample
rule run_dapars:
    input:
        config=f"data/{BATCH}/dapars/config/{{sample}}_config.txt",
        bedgraph=f"data/{BATCH}/bedgraph/{{sample}}.bedgraph"
    output:
        completed=f"data/{BATCH}/dapars/results/{{sample}}/dapars_completed.txt"
    params:
        dapars_script=DAPARS_SCRIPT,
        outdir=f"data/{BATCH}/dapars/results/{{sample}}"
    threads: 1
    log:
        f"data/{BATCH}/logs/dapars/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        # Run DaPars
        python {params.dapars_script} {input.config} 2>&1 | tee {log}

        # Create completion marker
        echo "Completed at $(date)" > {output.completed}
        """

# Run rMATS prep for each sample
rule rmats_prep:
    input:
        bam_list=f"data/{BATCH}/rmats/bam_lists/{{sample}}.txt",
        bam=f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        completed=f"data/{BATCH}/rmats/prep/{{sample}}/prep_completed.txt"
    params:
        gtf=GTF_PATH,
        rmats=RMATS_PATH,
        read_length=READ_LENGTH,
        read_type=READ_TYPE,
        outdir=f"data/{BATCH}/rmats/prep/{{sample}}",
        tmpdir=f"data/{BATCH}/rmats/tmp/{{sample}}",
        python_env=RMATS_PYTHON
    threads: 1
    log:
        f"data/{BATCH}/logs/rmats/{{sample}}.log"
    shell:
        """
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        mkdir -p $(dirname {log})

        # Run rMATS prep
        {params.python_env} {params.rmats} \
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
            2>&1 | tee {log}

        # Create completion marker
        echo "Completed at $(date)" > {output.completed}
        """

# Create combined BAM list for rMATS post step
rule create_combined_bam_list:
    input:
        bams=expand(f"data/{BATCH}/star/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        combined_list=f"data/{BATCH}/rmats/bam_lists/all_samples.txt"
    shell:
        """
        mkdir -p $(dirname {output.combined_list})
        echo "{input.bams}" | tr ' ' ',' > {output.combined_list}
        """

# Run rMATS post for all samples combined
rule rmats_post:
    input:
        combined_list=f"data/{BATCH}/rmats/bam_lists/all_samples.txt",
        prep_completed=expand(f"data/{BATCH}/rmats/prep/{{sample}}/prep_completed.txt", sample=SAMPLES)
    output:
        completed=f"data/{BATCH}/rmats/post/post_completed.txt"
    params:
        gtf=GTF_PATH,
        rmats=RMATS_PATH,
        read_length=READ_LENGTH,
        read_type=READ_TYPE,
        outdir=f"data/{BATCH}/rmats/post",
        tmpdir=f"data/{BATCH}/rmats/tmp",
        python_env=RMATS_PYTHON
    threads: 80
    log:
        f"data/{BATCH}/logs/rmats/post.log"
    shell:
        """
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        mkdir -p $(dirname {log})

        # Run rMATS post
        {params.python_env} {params.rmats} \
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
            2>&1 | tee {log}

        # Create completion marker
        echo "Completed at $(date)" > {output.completed}
        """
