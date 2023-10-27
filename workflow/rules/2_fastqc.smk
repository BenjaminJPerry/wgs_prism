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


wildcard_constraints: sample = "(?!Undetermined).+"

# Global variables

# config dictionary values to be defined on running snakemake with --config flag
fastqc_in_root = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert")
fastqc_in_samples = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/{sample}.fastq.gz")

(SAMPLES,) = glob_wildcards( os.path.join(fastqc_in_root,"{sample,(?!Undetermined).*}.fastq.gz") ).sample

fastqc_out_root = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
fastqc_out_samples_zips = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.zip")
fastqc_out_samples_htmls = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.html")

fastqc_log = os.path.join(config["OUT_ROOT"], "logs/2.1_run_fastqc.{sample}.log")
fastqc_benchmark = os.path.join(config["OUT_ROOT"], "benchmarks/run_fastqc.{sample}.log")


rule targets:
    input:
        expand(fastqc_out_samples_zips, sample = SAMPLES),
        expand(fastqc_out_samples_htmls, sample = SAMPLES),


rule run_fastqc:
    input:
        fastq = fastqc_in_samples
    output:
        zip = fastqc_out_samples_zips,
        html = fastqc_out_samples_htmls
    log:
        fastqc_log
    singularity:
        'docker://biocontainers/fastqc:v0.11.9_cv8'
    benchmark:
        fastqc_benchmark
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda w: config["fastqc_walltime"],
    shell:
        """ 
        mkdir -p {fastqc_out_root}

        fastqc -t {threads} -o {fastqc_out_root} {input.fastq} > {log} 2>&1

        success_landmark={output.zip}

        if [ ! -f $success_landmark ]; then
        echo "fastqc  of {wildcards.sample} did not generate the expected output file {output.zip} "
        exit 1
        fi 

        """

