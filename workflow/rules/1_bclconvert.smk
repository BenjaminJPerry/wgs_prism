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
# bclconvert_in_path = os.path.join(config["IN_ROOT"], config["RUN"])
# bclconvert_out_path = os.path.join(config["OUT_ROOT"], config["RUN"], "SampleSheet")
# sample_sheet_path = os.path.join(config["OUT_ROOT"], config["RUN"], "SampleSheet.csv")


rule targets:
    input:
        expand("{out_root}/{run}/SampleSheet/bclconvert/Reports/Top_Unknown_Barcodes.csv", out_root = config["OUT_ROOT"], run = config["RUN"]),
        expand("{out_root}/{run}/SampleSheet/bclconvert/Logs/FastqComplete.txt", out_root = config["OUT_ROOT"], run = config["RUN"]),


rule run_bclconvert:
    input:
        run_in = expand("{in_root}/{run}", in_root = config["IN_ROOT"], run = config["RUN"]),
        sample_sheet = expand("{out_root}/{run}/SampleSheet.csv", out_root = config["OUT_ROOT"], run = config["RUN"])
    output:
        bclconvert_out = directory(expand("{out_root}/{run}/SampleSheet/bclconvert", out_root = config["OUT_ROOT"], run = config["RUN"])),
        fastq_complete = expand("{out_root}/{run}/SampleSheet/bclconvert/Logs/FastqComplete.txt", out_root = config["OUT_ROOT"], run = config["RUN"]),
        top_unknown = expand("{out_root}/{run}/SampleSheet/bclconvert/Reports/Top_Unknown_Barcodes.csv", out_root = config["OUT_ROOT"], run = config["RUN"])
    log:
        expand("{out_root}/{run}/logs/1_run_bclconvert.log", out_root = config["OUT_ROOT"], run = config["RUN"])
    conda:
        "bclconvert" # Container OR module in the future!
    benchmark:
        expand("{out_root}/{run}/benchmarks/run_bclconvert.log", out_root = config["OUT_ROOT"], run = config["RUN"])
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
    shell:
        "bcl-convert --bcl-input-directory {input.run_in} --sample-sheet {input.sample_sheet} --output-directory {output.bclconvert_out} > {log} 2>&1"


        # """
        # # run bcl-convert
        # # report version 
        
        # # echo "bcl-convert version in use:"
        
        # bcl-convert -V 
        
        # bcl-convert --bcl-input-directory {input.run_in} --sample-sheet {input.sample_sheet} --output-directory {output.bclconvert_out} > {log} 2>&1
        
        # # if [ $? != 0 ]; then
        # # echo "error: bclconvert of {input.sample_sheet} - returned an error code."
        # # exit 1
        # # fi
        
        # """
        
