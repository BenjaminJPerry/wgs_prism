# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz
# configfile: "config/config.yaml"


import os


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


FIDs, = glob_wildcards("fastq/{samples}_R1_001.fastq.gz")


rule all:
    input:
        "results/Sequence_Production_and_QC_Report.multiqc.html",


rule downsample_for_QC:
    input:
        read1 = "fastq/{samples}_R1_001.fastq.gz",
        read2 = "fastq/{samples}_R2_001.fastq.gz"
    output:
        downed1 = "results/00_downsample/{samples}.DS.R1.fastq.gz",
        downed2 = "results/00_downsample/{samples}.DS.R2.fastq.gz",
    conda:
        "seqtk"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 6 + ((attempt - 1) * 6),
        partition="compute",
    shell:
        "seqtk sample -s1953 {input.read1} 1000000 | gzip > {output.downed1} "
        "&& "
        "seqtk sample -s1953 {input.read2} 1000000 | gzip > {output.downed2} "


rule bbduk_read_trim:
    input:
        downed1 = "results/00_downsample/{samples}.DS.R1.fastq.gz",
        downed2 = "results/00_downsample/{samples}.DS.R2.fastq.gz",
        adapters = "resources/adapters.fa",
    output:
        bbdukRead1 = "results/01_readMasking/{samples}_R1_bbduk.fastq.gz",
        bbdukRead2 = "results/01_readMasking/{samples}_R2_bbduk.fastq.gz"
    log:
        "logs/bbduk/{samples}.bbduk.log"
    benchmark:
        "benchmarks/bbduk.{samples}.align.log"
    conda:
        "bbduk"
    threads: 8
    resources: #TODO Update
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 6),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
        partition='compute',
    shell:
        "bbduk.sh "
        "threads={threads} "
        "in1={input.downed1} "
        "in2={input.downed2} "
        "ref={input.adapters} "
        "ktrim=rl "
        "k=19 mink=9 hdist=1 "
        "rcomp=t " 
        "trimpolyg=2 "
        "trimpolya=10 "
        "qtrim=f "
        "minlen=10 "
        "out1={output.bbdukRead1} "
        "out2={output.bbdukRead2} "
        "2>&1 | tee {log} "


rule bowtie2_SILVA_alignment_read1:
    input:
        bbdukRead1 = "results/01_readMasking/{samples}_R1_bbduk.fastq.gz",
    output:
        bowtie2_R1 = "results/02_SILVA/{samples}.DS.R1.bowtie2.log",
        silva_R1 = "results/02_SILVA/{samples}_R1_bbduk_silva.fastq",
    benchmark:
        "benchmarks/bowtie2_SILVA_alignment_read1.{samples}.txt"
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 20),
        partition = "compute"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x /datasets/2024-silva-rrna/SILVA138.1 "
        "--un {output.silva_R1} "
        "-U {input.bbdukRead1} "
        "1> /dev/null "
        "2> {output.bowtie2_R1} "


rule bowtie2_SILVA_alignment_read2:
    input:
        bbdukRead2 = "results/01_readMasking/{samples}_R2_bbduk.fastq.gz"
    output:
        bowtie2_R2 = "results/02_SILVA/{samples}.DS.R2.bowtie2.log",
        silva_R2 = "results/02_SILVA/{samples}_R2_bbduk_silva.fastq",
    benchmark:
        "benchmarks/bowtie2_SILVA_alignment_read2.{samples}.txt"
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 20),
        partition = "compute"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x /datasets/2024-silva-rrna/SILVA138.1 "
        "--un {output.silva_R2} "
        "-U {input.bbdukRead2} "
        "1> /dev/null "
        "2> {output.bowtie2_R2} "


rule kraken2_read_composition_read1:
    input:
        filtered_read1 = "results/02_SILVA/{samples}_R1_bbduk_silva.fastq",
    output:
        k2Output = temp("results/03_kraken2/{samples}.DS.R1.nt.k2"),
        k2Report_R1 = "results/03_kraken2/{samples}.DS.R1.nt.report.kraken2",
    log:
        "logs/kraken2/kraken2_read_composition.{samples}.NCBI-nt.log",
    benchmark:
        "benchmarks/kraken2_read_composition.{samples}.txt"
    conda:
        "kraken2"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 714 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 45 + ((attempt - 1) * 40),
        partition = "hugemem"
    shell:
        "kraken2 "
        "--use-names "
        "--db /datasets/2024-kraken2-indices/k2_nt_20231129 " 
        "-t {threads} "
        "--report {output.k2Report_R1} "
        "--report-minimizer-data "
        "--output {output.k2Output} "
        "{input.filtered_read1} "
        "2>&1 | tee {log} "


rule kraken2_read_composition_read2:
    input:
        filtered_read2 = "results/02_SILVA/{samples}_R2_bbduk_silva.fastq",
    output:
        k2Output = temp("results/03_kraken2/{samples}.DS.R2.nt.k2"),
        k2Report_R2 = "results/03_kraken2/{samples}.DS.R2.nt.report.kraken2",
    log:
        "logs/kraken2/kraken2_read_composition_read2.{samples}.NCBI-nt.log",
    benchmark:
        "benchmarks/kraken2_read_composition_read2.{samples}.txt"
    conda:
        "kraken2"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 714 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 45 + ((attempt - 1) * 40),
        partition = "hugemem"
    shell:
        "kraken2 "
        "--use-names "
        "--db /datasets/2024-kraken2-indices/k2_nt_20231129 " 
        "-t {threads} "
        "--report {output.k2Report_R2} "
        "--report-minimizer-data "
        "--output {output.k2Output} "
        " {input.filtered_read2} "
        "2>&1 | tee {log} "


rule fastqc_read1:
    input:
        read1 = "fastq/{samples}_R1_001.fastq.gz",
    output:
        html = "results/00_QC/fastqc/{samples}_R1_001_fastqc.html",
        zip = "results/00_QC/fastqc/{samples}_R1_001_fastqc.zip"
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
        partition="compute",
    shell:
        "fastqc "
        "-o results/00_QC/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.read1}"


rule fastqc_read2:
    input:
        read2 = "fastq/{samples}_R2_001.fastq.gz",
    output:
        html = "results/00_QC/fastqc/{samples}_R2_001_fastqc.html",
        zip = "results/00_QC/fastqc/{samples}_R2_001_fastqc.zip"
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
        partition="compute",
    shell:
        "fastqc "
        "-o results/00_QC/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.read2}"


rule fastqc_filtered_read1: #TODO
    input:
        filtered_read1 = "results/02_SILVA/{samples}_R1_bbduk_silva.fastq",
    output:
        html = "results/00_QC/fastqc/{samples}_R1_bbduk_silva_fastqc.html",
        zip = "results/00_QC/fastqc/{samples}_R1_bbduk_silva_fastqc.zip"
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
        partition="compute",
    shell:
        "fastqc "
        "-o results/00_QC/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.filtered_read1}"


rule fastqc_filtered_read2: #TODO
    input:
        filtered_read2 = "results/02_SILVA/{samples}_R2_bbduk_silva.fastq",
    output:
        html = "results/00_QC/fastqc/{samples}_R2_bbduk_silva_fastqc.html",
        zip = "results/00_QC/fastqc/{samples}_R2_bbduk_silva_fastqc.zip"
    conda:
        "fastqc-0.12.1"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 6 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 15),
        partition="compute",
    shell:
        "fastqc "
        "-o results/00_QC/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.filtered_read2}"


rule genome_alignment_check_R1:
    input:
        silva_R1 = "results/02_SILVA/{samples}_R1_bbduk_silva.fastq",
    output:
        bowtie2_genome = "results/02_REF/{samples}.DS.genome_alignment.bowtie2.R1.log",
    benchmark:
        "benchmarks/genome_alignment_check.R1.{samples}.txt"
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 30),
        partition = "compute"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x resources/GRCh38 "
        "-U {input.silva_R1} "
        "1> /dev/null "
        "2> {output.bowtie2_genome} "


rule genome_alignment_check_R2:
    input:
        silva_R2 = "results/02_SILVA/{samples}_R2_bbduk_silva.fastq",
    output:
        bowtie2_genome = "results/02_REF/{samples}.DS.genome_alignment.bowtie2.R2.log",
    benchmark:
        "benchmarks/genome_alignment_check.R2.{samples}.txt"
    conda:
        "bowtie2-2.5.1"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 30),
        partition = "compute"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x resources/GRCh38 "
        "-U {input.silva_R2} "
        "1> /dev/null "
        "2> {output.bowtie2_genome} "


rule multiQC_report:
    input:
        bclconvert_runinfo = "results/00_bclconvert/RunInfo.xml",
        bclconvert_demux_stats = "results/00_bclconvert/Demultiplex_Stats.csv",
        bclconvert_qual_metrics = "results/00_bclconvert/Quality_Metrics.csv",
        bclconvert_adpt_metrics = "results/00_bclconvert/Adapter_Metrics.csv",
        bclconvert_top_unknown = "results/00_bclconvert/Top_Unknown_Barcodes.csv",
        fastqc_read1 = expand("results/00_QC/fastqc/{samples}_R1_001_fastqc.zip", samples = FIDs),
        fastqc_read2 = expand("results/00_QC/fastqc/{samples}_R2_001_fastqc.zip", samples = FIDs),
        bbduk_log = expand("logs/bbduk/{samples}.bbduk.log", samples = FIDs),
        bowtie2_R1 = expand("results/02_SILVA/{samples}.DS.R1.bowtie2.log", samples = FIDs),
        bowtie2_R2 = expand("results/02_SILVA/{samples}.DS.R2.bowtie2.log", samples = FIDs),
        kraken2_R1 = expand("results/03_kraken2/{samples}.DS.R1.nt.report.kraken2", samples = FIDs),
        kraken2_R2 = expand("results/03_kraken2/{samples}.DS.R2.nt.report.kraken2", samples = FIDs),
        bowtie2_genome_R1 = expand("results/02_REF/{samples}.DS.genome_alignment.bowtie2.R1.log", samples = FIDs),
        bowtie2_genome_R2 = expand("results/02_REF/{samples}.DS.genome_alignment.bowtie2.R2.log", samples = FIDs),
        fastqc_filtered_read1 = expand("results/00_QC/fastqc/{samples}_R1_bbduk_silva_fastqc.zip", samples = FIDs), 
        fastqc_filtered_read2 = expand("results/00_QC/fastqc/{samples}_R2_bbduk_silva_fastqc.zip", samples = FIDs), 
        multiQC_config = "resources/multiQC_config.yaml",
    output:
        multiQC ="results/Sequence_Production_and_QC_Report.multiqc.html"
    conda:
        "multiqc"
    log:
        "results/logs/multiQC/multiQC_report.log",
    benchmark:
        "results/benchmarks/multiQC_report_stats.txt",
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 15),
        partition="compute",
    shell:
        "multiqc "
        "-n results/Sequence_Production_and_QC_Report.multiqc "
        "-s "
        "-f "
        "-c {input.multiQC_config} "
        "--interactive "
        "{input.bclconvert_runinfo} "
        "{input.bclconvert_demux_stats} "
        "{input.bclconvert_qual_metrics} "
        "{input.bclconvert_adpt_metrics} "
        "{input.bclconvert_top_unknown} "
        "{input.fastqc_read1} "
        "{input.fastqc_read2} "
        "{input.bbduk_log} "
        "{input.fastqc_filtered_read1} "
        "{input.fastqc_filtered_read2} "
        "{input.bowtie2_R1} "
        "{input.bowtie2_R2} "
        "{input.kraken2_R1} "
        "{input.kraken2_R2} "
        "{input.bowtie2_genome_R1} "
        "{input.bowtie2_genome_R2} "
        "2>&1 | tee {log}"



