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
kmer_in_root = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert")
kmer_in_samples = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/{sample}.fastq.gz")

(SAMPLES,) = glob_wildcards( os.path.join(kmer_in_root, "{sample,(?!Undetermined).*}.fastq.gz") )

kmer_out_root = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_analysis")



# Path and file name construction for rule downsample_fastq
sampling_rate = str(config["SAMPLE_RATE"])
downsample_out_samples_root = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_run/fastq_sample")
downsample_out_samples = "{sample}.fastq.gz" + "." + "s" + sampling_rate + "." + "fastq.gz"
downsample_out_samples_path = os.path.join(downsample_out_samples_root, downsample_out_samples)

downsample_log_files = "logs/2.2.1_downsample_fastq.{sample}.s" + sampling_rate + "." + "log"
downsample_logs_path = os.path.join(config["OUT_ROOT"], downsample_log_files)

downsample_benchmark_files = "benchmarks/downsample_fastq.{sample}.s" + sampling_rate + "." + "txt"
downsample_benchmark_path = os.path.join(config["OUT_ROOT"], downsample_benchmark_files)


# Path and file name construction for rule fastq_to_fasta
kmer_parameters = config["KMER_PARAMETERS"]
kmer_moniker = kmer_parameters.replace(" ", "").replace("-", "")

kmer_fastq_to_fasta_root = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_run/kmer_analysis")
kmer_fastq_to_fasta_out_samples = downsample_out_samples + "." + kmer_moniker + "." + "1"
kmer_fastq_to_fasta_out_samples_path = os.path.join(kmer_fastq_to_fasta_root, kmer_fastq_to_fasta_out_samples)

kmer_fastq_to_fasta_out_log_files = "logs/2.2.2_kmer_fastq_to_fasta.{sample}.s" + sampling_rate + "." + "log"
kmer_fastq_to_fasta_out_logs_path = os.path.join(config["OUT_ROOT"], kmer_fastq_to_fasta_out_log_files)

kmer_fastq_to_fasta_benchmark_files = "benchmarks/kmer_fastq_to_fasta.{sample}.s" + sampling_rate + "." + "txt"
kmer_fastq_to_fasta_benchmark_path = os.path.join(config["OUT_ROOT"], kmer_fastq_to_fasta_benchmark_files)

# Path and file name construction for rule run_kmer_prism
kmer_prism_root = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_run/kmer_analysis")
kmer_prism_out_samples_frequency = kmer_fastq_to_fasta_out_samples + "." + "kmer_prism" + "." + "frequency.txt"
kmer_prism_out_samples_pickle = kmer_fastq_to_fasta_out_samples + "." + "kmderdist" + "." + "pickle"

kmer_prism_out_samples_frequency_path = os.path.join(kmer_prism_root, kmer_prism_out_samples_frequency)
kmer_prism_out_samples_pickle_path = os.path.join(kmer_prism_root, kmer_prism_out_samples_pickle)

kmer_prism_out_log_files = "logs/2.2.3_run_kmer_prism.pickle.{sample}" + ".log"
kmer_prism_out_logs_path = os.path.join(config["OUT_ROOT"], kmer_prism_out_log_files)

kmer_prism_out_benchmark_files = "benchmarks/run_kmer_prism.pickle.{sample}" + ".txt"
kmer_prism_out_benchmark_path = os.path.join(config["OUT_ROOT"], kmer_fastq_to_fasta_benchmark_files)


# Beging Snakemake rule definitions

rule targets:
    input:
        expand(kmer_prism_out_samples_pickle_path, sample = SAMPLES),


rule downsample_fastq:
    input:
        kmer_in_samples,
    output:
        temp(downsample_out_samples_path)
    log:
        downsample_logs_path
    conda:
        "envs/seqkit.yaml"
    benchmark:
        downsample_benchmark_path
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
    shell:
        """ 
        seqkit sample -s 1953 -p .0002 --threads {threads} {input} -o {output} > {log} 2>&1

        """ 


rule fastq_to_fasta:
    input:
        downsample_out_samples_path
    output:
        temp(kmer_fastq_to_fasta_out_samples_path)
    log:
        kmer_fastq_to_fasta_out_logs_path
    conda:
        "envs/seqkit.yaml"
    benchmark:
        kmer_fastq_to_fasta_benchmark_path
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
    shell:
        """ 
        seqkit fq2fa {input} -o {output} > {log} 2>&1

        """    


rule run_kmer_prism:
    input:
        kmer_fastq_to_fasta_out_samples_path
    output:
        txt = temp(kmer_prism_out_samples_frequency_path),
        pickle = kmer_prism_out_samples_pickle_path
    log:
        kmer_prism_out_logs_path
    conda:
        "envs/biopython.yaml"
    benchmark:
        kmer_prism_out_benchmark_path
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
    shell:   
        """ 
        
        workflow/scripts/kmer_prism.py -f fasta -k 6 -A -b {kmer_prism_root} -o {output.txt} {input} > {log} 2>&1

        """


# rule aggregate_kmer_spectra:
#     input:
#         fastq = fastqc_in_samples
#     output:
#         zip = fastqc_out_samples_zips,
#         html = fastqc_out_samples_htmls
#     log:
#         fastqc_log
#     conda:
#         'envs/biopython.yaml'
#     benchmark:
#         fastqc_benchmark
#     threads: 12
#     resources:
#         mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
#         time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
#     shell:
#         """ 
#         # below uses the .kmerdist.pickle distribution files to make the final spectra
#         # (note that the -k 6 arg here is not actually used , as the distributions have already been done by the make step)
#         rm -f $OUT_DIR/kmer_summary_plus.${parameters_moniker}.txt
#         tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -t zipfian -o $OUT_DIR/kmer_summary_plus.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1

#         rm -f $OUT_DIR/kmer_frequency_plus.${parameters_moniker}.txt
#         tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -t frequency -o $OUT_DIR/kmer_frequency_plus.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1
        
#         rm -f $OUT_DIR/kmer_summary.${parameters_moniker}.txt
#         tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -a CGAT -t zipfian -o $OUT_DIR/kmer_summary.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1

#         rm -f  $OUT_DIR/kmer_frequency.${parameters_moniker}.txt
#         tardis.py --hpctype $HPC_TYPE -d $OUT_DIR --shell-include-file configure_biopython_env.src kmer_prism.py -k 6 -a CGAT -t frequency -o $OUT_DIR/kmer_frequency.${parameters_moniker}.txt -b $OUT_DIR $SUMMARY_TARGETS >> $OUT_DIR/kmer_prism.log 2>&1
        
#         """



# rule plot_kmer_prism:
#     input:
#         fastq = fastqc_in_samples
#     output:
#         zip = fastqc_out_samples_zips,
#         html = fastqc_out_samples_htmls
#     log:
#         fastqc_log
#     conda:
#         'envs/biopython.yaml'
#     benchmark:
#         fastqc_benchmark
#     threads: 12
#     resources:
#         mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
#         time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
#     shell:
#         """ 
#         for version in "" "_plus" ; do
#             rm -f $OUT_DIR/kmer_summary.txt
#             cp -s $OUT_DIR/kmer_summary${version}.${parameters_moniker}.txt $OUT_DIR/kmer_summary.txt
#             tardis.py 
#                     --hpctype $HPC_TYPE 
#                     -d $OUT_DIR 
#                     --shell-include-file configure_bioconductor_env.src 
#                     Rscript --vanilla 
#                         $OUT_DIR/kmer_plots.r datafolder=$OUT_DIR >> $OUT_DIR/kmer_prism.log 2>&1
            
#             for output in kmer_entropy kmer_zipfian_comparisons kmer_zipfian zipfian_distances; do
#                 if [ -f $OUT_DIR/${output}.jpg ]; then
#                     mv $OUT_DIR/${output}.jpg $OUT_DIR/${output}${version}.${parameters_moniker}.jpg
#                 fi
#             done
            
#             for output in heatmap_sample_clusters  zipfian_distances_fit ; do
#                 if [ -f $OUT_DIR/${output}.txt ]; then
#                     mv $OUT_DIR/${output}.txt $OUT_DIR/${output}${version}.${parameters_moniker}.txt
#                 fi
#             done

#         done


#         """


