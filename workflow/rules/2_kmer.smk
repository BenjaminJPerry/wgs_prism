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


#wildcard_constraints: sample = "(?!Undetermined).+"

# Global variables

# config dictionary values to be defined on running snakemake with --config flag
fastqc_in_root = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert")
fastqc_in_samples = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/{sample}.fastq.gz")

SAMPLES = glob_wildcards(os.path.join(fastqc_in_root,"{sample, (?!Undetermined).*}.fastq.gz")).sample

fastqc_out_root = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc")
fastqc_out_samples_zips = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.zip")
fastqc_out_samples_htmls = os.path.join(config["OUT_ROOT"], "SampleSheet/fastqc_run/fastqc/{sample}_fastqc.html")

fastqc_log = os.path.join(config["OUT_ROOT"], "logs/2_run_fastqc.{sample}.log")
fastqc_benchmark = os.path.join(config["OUT_ROOT"], "benchmarks/run_fastqc.{sample}.log")


rule targets:
    input:
        expand(fastqc_out_samples_zips, sample = SAMPLES),
        expand(fastqc_out_samples_htmls, sample = SAMPLES),

rule fastqc:
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
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
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

sequencing_qc_prism.sh:
$OUT_ROOT/kmer_prism.sh 
    -W $MAX_WALL_TIME 
    -B $MEM_PER_CPU 
    -a fastq 
    -p \"-k 6 -A\" 
    -O $OUT_ROOT/kmer_analysis $OUT_ROOT/fastq_sample/*.fastq.gz  
    >  $OUT_ROOT/kmer_analysis/kmer_analysis.log 2>&1

kmer_prism.sh:
function run_prism() {
   # this distributes the kmer distribtion builds for each file across the cluster
   make -f kmer_prism.mk -d -k  --no-builtin-rules -j $NUM_THREADS `cat $OUT_DIR/kmer_targets.txt` > $OUT_DIR/kmer_prism.log 2>&1
   # this uses the pickled distributions to make the final spectra
   # (note that the -k 6 arg here is not actually used , as the distributions have already been done by the make step)
   rm -f $OUT_DIR/kmer_summary_plus.${parameters_moniker}.txt
   tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -t zipfian -o $OUT_DIR/kmer_summary_plus.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1

   rm -f $OUT_DIR/kmer_frequency_plus.${parameters_moniker}.txt
   tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -t frequency -o $OUT_DIR/kmer_frequency_plus.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1
   
   rm -f $OUT_DIR/kmer_summary.${parameters_moniker}.txt
   tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -a CGAT -t zipfian -o $OUT_DIR/kmer_summary.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1

   rm -f  $OUT_DIR/kmer_frequency.${parameters_moniker}.txt
   tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -a CGAT -t frequency -o $OUT_DIR/kmer_frequency.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1

   ### Plotting the outputs from above ###
   # first do plots including N's , then rename and do plots excluding N's 
   for version in "" "_plus" ; do
      rm -f $OUT_DIR/kmer_summary.txt
      cp -s $OUT_DIR/kmer_summary${version}.${parameters_moniker}.txt $OUT_DIR/kmer_summary.txt
      tardis.py 
            --hpctype $HPC_TYPE 
            -d $OUT_DIR 
            --shell-include-file configure_bioconductor_env.src 
            Rscript --vanilla 
                $OUT_DIR/kmer_plots.r datafolder=$OUT_DIR >> $OUT_DIR/kmer_prism.log 2>&1
      
      for output in kmer_entropy kmer_zipfian_comparisons kmer_zipfian zipfian_distances; do
         if [ -f $OUT_DIR/${output}.jpg ]; then
            mv $OUT_DIR/${output}.jpg $OUT_DIR/${output}${version}.${parameters_moniker}.jpg
         fi
      done
      
      for output in heatmap_sample_clusters  zipfian_distances_fit ; do
         if [ -f $OUT_DIR/${output}.txt ]; then
            mv $OUT_DIR/${output}.txt $OUT_DIR/${output}${version}.${parameters_moniker}.txt
         fi
      done

   done
}

Where,
${parameters_moniker} = -p \"-k 6 -A\"
$OUT_ROOT = SampleSheet/kmer_run
-O $OUT_ROOT/kmer_analysis


