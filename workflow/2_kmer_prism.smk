# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/pipeline_config.yaml'


import os
import pandas as pd


wildcard_constraints: sample = "(?!Undetermined).+"

### Global variables ###
# config dictionary values to be defined on running snakemake with --config flag
kmer_in_root = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert")
kmer_in_samples = os.path.join(config["OUT_ROOT"], "SampleSheet/bclconvert/{sample}.fastq.gz")

(SAMPLES,) = glob_wildcards( os.path.join(kmer_in_root, "{sample,(?!Undetermined).*}.fastq.gz") )
min_reads = config["MIN_READS"]

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
kmer_prism_out_samples_pickle = kmer_fastq_to_fasta_out_samples + "." + "kmerdist" + "." + "pickle"

kmer_prism_summary_metrics = os.path.join(kmer_prism_root + "seqkit.downsample.summary.txt")

kmer_prism_out_samples_frequency_path = os.path.join(kmer_prism_root, kmer_prism_out_samples_frequency)
kmer_prism_out_samples_pickle_path = os.path.join(kmer_prism_root, kmer_prism_out_samples_pickle)

kmer_prism_out_log_files = "logs/2.2.3_run_kmer_prism.pickle.{sample}" + ".log"
kmer_prism_out_logs_path = os.path.join(config["OUT_ROOT"], kmer_prism_out_log_files)

kmer_prism_out_benchmark_files = "benchmarks/run_kmer_prism.pickle.{sample}" + ".txt"
kmer_prism_out_benchmark_path = os.path.join(config["OUT_ROOT"], kmer_fastq_to_fasta_benchmark_files)


# Path and file name construction for rule aggregate_kmer_spectra
kmer_agg_summary_plus = "kmer_summary_plus" + "." + kmer_moniker + ".txt"
kmer_agg_summary_plus_path = os.path.join(kmer_prism_root, kmer_agg_summary_plus)

kmer_agg_frequency_plus = "kmer_frequency_plus" + "." + kmer_moniker + ".txt"
kmer_agg_frequency_plus_path  = os.path.join(kmer_prism_root, kmer_agg_frequency_plus)

kmer_agg_summary = "kmer_summary" + "." + kmer_moniker + ".txt"
kmer_agg_summary_path  = os.path.join(kmer_prism_root, kmer_agg_summary)

kmer_agg_frequency = "kmer_frequency" + "." + kmer_moniker + ".txt"
kmer_agg_frequency_path  = os.path.join(kmer_prism_root, kmer_agg_frequency)

kmer_agg_plot_data = "kmer_summary.txt"
kmer_agg_plot_data_path  = os.path.join(kmer_prism_root, kmer_agg_plot_data)

kmer_agg_out_log_files = "logs/2.2.4_run_aggregate_kmer_spectra.log"
kmer_agg_out_logs_path = os.path.join(config["OUT_ROOT"], kmer_agg_out_log_files)

kmer_agg_out_benchmark_files = "benchmarks/run_aggregate_kmer_spectra.txt"
kmer_agg_out_benchmark_path = os.path.join(config["OUT_ROOT"], kmer_agg_out_benchmark_files)


# Path and file name construction for rule aggregate_kmer_spectra
kmer_agg_plot_data_dir = os.path.join(config["OUT_ROOT"], "SampleSheet/kmer_run/kmer_analysis")
kmer_zipfian_plot = "kmer_zipfian.jpg"
kmer_zipfian_plot_path = os.path.join(kmer_agg_plot_data_dir, kmer_zipfian_plot)

kmer_entropy_plot = "kmer_entropy.jpg"
kmer_entropy_plot_path = os.path.join(kmer_agg_plot_data_dir, kmer_entropy_plot)

kmer_zipfian_comparison_plot = "kmer_zipfian_comparisons.jpg"
kmer_zipfian_comparison_plot_path =os.path.join(kmer_agg_plot_data_dir, kmer_zipfian_comparison_plot)

kmer_zipfian_distances = "zipfian_distances.jpg"
kmer_zipfian_distances_path = os.path.join(kmer_agg_plot_data_dir, kmer_zipfian_distances)

plot_kmer_spectra_log = "logs/2.2.4_run_aggregate_kmer_spectra.log"
plot_kmer_spectra_log_path = os.path.join(config["OUT_ROOT"], plot_kmer_spectra_log)

plot_kmer_spectra_benchmark = "benchmarks/run_plot_kmer_spectra.txt"
plot_kmer_spectra_benchmark_path = os.path.join(config["OUT_ROOT"], plot_kmer_spectra_benchmark)


# Beging Snakemake rule definitions
rule targets:
    input:
        kmer_entropy_plot_path,
        kmer_zipfian_comparison_plot_path,
        kmer_zipfian_distances_path


rule downsample_fastq:
    input:
        kmer_in_samples,
    output:
        downsample_out_samples_path
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
        seqkit sample -s 1953 -p {sampling_rate} --threads {threads} {input} -o {output} > {log} 2>&1

        """ 


rule fastq_to_fasta:
    input:
        downsample_out_samples_path
    output:
        kmer_fastq_to_fasta_out_samples_path
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


checkpoint summary_QC_fasta:
    input:
        expand(kmer_fastq_to_fasta_out_samples_path, sample = SAMPLES)
    output:
        kmer_prism_summary_metrics
    conda:
        "envs/seqkit.yaml"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 24),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
    shell:
        """
        
        seqkit stats -j {threads} -a {input} > {output}
        
        """


rule run_kmer_prism:
    input:
        kmer_fastq_to_fasta_out_samples_path,
    output:
        txt = kmer_prism_out_samples_frequency_path,
        pickle = kmer_prism_out_samples_pickle_path,
    log:
        kmer_prism_out_logs_path
    conda:
        "envs/biopython.yaml"
    benchmark:
        kmer_prism_out_benchmark_path
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 120),
    shell:   
        """ 
    
        workflow/scripts/kmer_prism.py -f fasta -p {threads} -k 6 -A -b {kmer_prism_root} -o {output.txt} {input} > {log} 2>&1

        success_landmark={output.pickle}

        if [ ! -f $success_landmark ]
        then
            echo "error: kmer_prism.py did not generate the expected output file {output.pickle} "
            exit 1
        else
            exit 0
        fi

        """


def get_summary_QC_fasta_passing_pickles(wildcards, minReads = min_reads):
    import pandas as pd
    file = checkpoints.summary_QC_fasta.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(kmer_prism_out_samples_pickle_path, sample = passed)


def get_summary_QC_fasta_passing_fasta(wildcards, minReads = min_reads):
    import pandas as pd
    file = checkpoints.summary_QC_fasta.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(kmer_fastq_to_fasta_out_samples_path, sample = passed)


rule aggregate_kmer_spectra:
    input:
        pickles = get_summary_QC_fasta_passing_pickles,
        fastas = get_summary_QC_fasta_passing_fasta
    output:
        summary_plus = kmer_agg_summary_plus_path,
        frequency_plus = kmer_agg_frequency_plus_path,
        summary = kmer_agg_summary_path,
        frequency = kmer_agg_frequency_path,
        plot_data = kmer_agg_plot_data_path
    log:
        kmer_agg_out_logs_path
    conda:
        'envs/biopython.yaml'
    benchmark:
        kmer_agg_out_benchmark_path
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
    shell:
        """ 

        # below uses the .kmerdist.pickle distribution files to make the final spectra
        # (note that the -k 6 arg here is not actually used , as the distributions have already been done by the make step)

        rm -f {output.summary_plus}
        workflow/scripts/kmer_prism.py -p {threads} -k 6 -t zipfian -o {output.summary_plus} -b {kmer_prism_root} {input.fastas} >> {log} 2>&1
        if [ -s {output.summary_plus}]
        then
            echo "error: kmer_prism.py did not aggregate spectra into {output.summary_plus} " | tee >> {log}
            exit 1
        fi

        rm -f {output.frequency_plus}
        workflow/scripts/kmer_prism.py -p {threads} -k 6 -t frequency -o {output.frequency_plus} -b {kmer_prism_root} {input.fastas} >> {log} 2>&1
        if [ -s {output.frequency_plus}]
        then
            echo "error: kmer_prism.py did not aggregate spectra into {output.frequency_plus} " | tee >> {log}
            exit 1
        fi

        rm -f {output.summary}
        workflow/scripts/kmer_prism.py -p {threads} -k 6 -a CGAT -t zipfian -o {output.summary} -b {kmer_prism_root} {input.fastas} >> {log} 2>&1
        if [ -s {output.summary}]
        then
            echo "error: kmer_prism.py did not aggregate spectra into {output.summary} " | tee >> {log}
            exit 1
        fi

        rm -f  {output.frequency}
        workflow/scripts/kmer_prism.py -p {threads} -k 6 -a CGAT -t frequency -o {output.frequency} -b {kmer_prism_root} {input.fastas} >> {log} 2>&1
        if [ -s {output.frequency}]
        then
            echo "error: kmer_prism.py did not aggregate spectra into {output.frequency} " | tee >> {log}
            exit 1
        fi

        cp -s {output.summary} {output.plot_data}
        """


rule plot_kmer_spectra:
    input: 
        plot_data = kmer_agg_plot_data_path,
    output:
        #kmer_zipfian_plot = kmer_zipfian_plot_path, 
        kmer_entropy_plot = kmer_entropy_plot_path, 
        kmer_zipfian_comparison_plot = kmer_zipfian_comparison_plot_path,
        kmer_zipfian_distances = kmer_zipfian_distances_path,
    log:
        plot_kmer_spectra_log_path
    conda:
        'bioconductor'
    benchmark:
        plot_kmer_spectra_benchmark_path
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
    shell:
        """

        Rscript --vanilla --verbose workflow/scripts/kmer_plots.r datafolder={kmer_prism_root} > {log} 2>&1


        """

