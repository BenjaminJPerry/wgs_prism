# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/pipeline_config.yaml'

import os
import pandas as pd


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
    print("Found: ")


# wildcard_constraints:
#     samples="\w+"

# Global variables


# FIDs, = glob_wildcards('results/01_cutadapt/{samples}.fastq.gz')


rule all:
    input:
        'results/00_QC/seqkit.report.KDTrim.txt',
     


rule get_samplesheet:


rule bclconvert:
    input:
        "results/01_cutadapt/{samples}.fastq.gz",
    output:
        temp("results/01_readMasking/{samples}.sana.fastq.gz"),
    log:
        "logs/sana/sana.{samples}.log"
    conda:
        "seqkit"
    benchmark:
        "benchmarks/sana.{samples}.log"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 8 + ((attempt - 1) * 10),
        partition='compute',
    shell:
        "seqkit sana "
        "-j {threads} "
        "{input} "
        "| gzip --fast -c > {output} 2> {log} "


checkpoint seqkitRaw:
    input:
        expand('results/01_readMasking/{samples}.sana.fastq.gz', samples = FIDs),
    output:
        'results/00_QC/seqkit.report.raw.txt'
    benchmark:
        'benchmarks/seqkitRaw.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0' 
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
        partition="compute"
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '


# STANDARD READ FILTERING AND QC RULES
rule bbduk:
    input:
        reads = 'results/01_readMasking/{samples}.sana.fastq.gz',
    output:
        bbdukReads = temp('results/01_readMasking/{samples}.bbduk.fastq.gz')
    log:
        'logs/bbduk/{samples}.bbduk.log'
    conda:
        'bbduk'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 8 + ((attempt - 1) * 10),
        partition='compute',
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.reads} '
        'entropy=0.3 '
        'entropywindow=50 '
        'trimpolygright=5 '
        'qtrim=r '
        'trimq=20 '
        'out={output.bbdukReads} '
        '2>&1 | tee {log}'



def get_seqkitMaskingBBDukReads_passing_samples(wildcards, minReads=min_reads):
    file = checkpoints.seqkitRaw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/01_readMasking/{samples}.bbduk.fastq.gz", samples = passed)


rule seqkitMaskingBBDukReads:
    input:
        bbdukReads = get_seqkitMaskingBBDukReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.bbduk.txt'
    benchmark:
        'benchmarks/seqkitMaskingBBDukReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 90 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.bbdukReads} > {output} '
