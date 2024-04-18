# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/pipeline_config.yaml'

import os
import pandas as pd


# onstart:
#     print(f"Working directory: {os.getcwd()}")
#     print("TOOLS: ")
#     os.system('echo "  bash: $(which bash)"')
#     os.system('echo "  PYTHON: $(which python)"')
#     os.system('echo "  CONDA: $(which conda)"')
#     os.system('echo "  SNAKEMAKE: $(which snakemake)"')
#     print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
#     os.system('echo "  PYTHON VERSION: $(python --version)"')
#     os.system('echo "  CONDA VERSION: $(conda --version)"')


# wildcard_constraints:
#     samples="\w+"

# Global variables
# config dictionary values to be defined on running snakemake with --config flag
run_name = config["RUN"]
# logs and reports for multiQC
# BCLConvert Reports
bclconvert_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/Reports")
# FastQC Reports
fastqc_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
# kmer_prism.py Reports
kmer_reports_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_run/kmer_analysis")

multiqc_report_file = run_name + ".multiqc.html"
multiqc_report_path = os.path.join(config["OUT_ROOT"], "SampleSheet", "multiqc", multiqc_report_file)

multiqc_data_dir = "multiqc_data"
multiqc_data_dir_path = os.path.join(config["OUT_ROOT"], "SampleSheet", "multiqc", multiqc_data_dir)


multiqc_log = "logs/3.0.0_run_multiqc.log"
multiqc_log_path = os.path.join(config["OUT_ROOT"], multiqc_log)

multiqc_benchmark = "benchmarks/run_multiqc.txt"
multiqc_benchmark_path = os.path.join(config["OUT_ROOT"], multiqc_benchmark)


rule targets:
    input:
        multiqc_report_path

rule run_multiqc:
    input:
        bclconvert_in = bclconvert_reports_dir,
        fastqc_in = fastqc_reports_dir,
        #kmer_in = kmer_reports_dir
    output:
        report = multiqc_report_path,
    log:
        multiqc_log_path
    conda:
        "multiqc"
    benchmark:
        multiqc_benchmark_path
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
    params:
        multiqc_config = config["multiqc_config"]
    shell:
        """

        multiqc --interactive --outdir {multiqc_data_dir_path} --filename {output.report} --force -c {params.multiqc_config} --data-dir --data-format tsv {input.bclconvert_in} {input.fastqc_in} > {log} 2>&1
        """

        # {input.kmer_in}
