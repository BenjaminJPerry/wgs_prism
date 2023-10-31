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


# wildcard_constraints:
#     samples="\w+"

# Global variables
# config dictionary values to be defined on running snakemake with --config flag

# logs and reports for multiQC
# BCLConvert Reports
bclconvert_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/Reports")
# FastQC Reports
fastqc_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
# kmer_prism.py Reports
kmer_agg_plot_data_dir = os.path.join(config["OUT_ROOT"], "SampleSheet")

multiQC_config = config["multiqc_config"]



rule targets:
    input:
        fastq_complete_path,
        top_unknown_path


rule run_bclconvert:
    input:
        run_in = bclconvert_in_path,
        sample_sheet = sample_sheet_path,
    output:
        bclconvert_out = directory(bclconvert_out_path),
        fastq_complete = fastq_complete_path,
        top_unknown = top_unknown_path
    log:
        bclconvert_log
    singularity:
        "docker://nfcore/bclconvert:3.9.3"
    benchmark:
        bclconvert_benchmark
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
    shell:
        """
        
        
        
        """
        
