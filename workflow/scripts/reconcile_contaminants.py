#!/usr/bin/env python
from __future__ import print_function
#########################################################################
# top level guery interface 
#########################################################################
import argparse
import sys
import os
import re
import itertools 

sys.path.append('/dataset/gseq_processing/active/bin/gquery')  # hack until gquery is packaged and has an installer  

def get_options():
    description = """
    """
    long_description = """

examples :

cat /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/adapter_list.txt /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/contaminant_list.txt | ./reconcile_contaminants.py /dataset/gseq_processing/scratch/illumina/novaseq/210702_A01439_0005_AHCC27DRXY/kmer_analysis/*.k6A.log
cat /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/adapter_list.txt /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/contaminant_list.txt | ./reconcile_contaminants.py /dataset/gseq_processing/scratch/illumina/hiseq/210628_D00390_0628_ACD9AUANXX/SampleSheet/gbs/kmer_analysis/*.k6A.log > hiseq_adapters_example.txt

cat /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/adapter_list.txt /stash/miniconda3/envs/bifo-essential/opt/fastqc-0.11.8/Configuration/contaminant_list.txt /dataset/gseq_processing/active/bin/gquery/database/t_BarcodePlates.csv  | /dataset/gseq_processing/active/bin/gbs_prism/reconcile_contaminants.py /dataset/2023_illumina_sequencing_a/scratch/postprocessing/illumina/novaseq/20230504_FS10001778_51_BSB09412-3119/SampleSheet/kmer_run/kmer_analysis/*.k6A.log | head -30


"""

    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', type=str, nargs='*',help='space-separated list of items to process (e.g. names of samples, subjects, libraries , sequence files, taxnames, taxids etc.')
    
    parser.add_argument('-t', '--task' , dest='task', required=False, type=str,
                        choices=[], help="what you want to get / do")
    
    parser.add_argument('-l','--list_name', dest='list_name', type=str, default=None, help='list name (if requesting list processing)')
    parser.add_argument('-n','--dry_run', dest='dry_run', action='store_const', default = False, const=True, help='dry run only')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_const', default = False, const=True, help='request (context-sensitive) verbosity')

    args = vars(parser.parse_args())

    return args

def get_fastqc_seqs():
    seqs = dict([ re.split("\t+", record.strip()) for record in sys.stdin if len( re.split("\t+", record.strip()) ) == 2])

    # one seq has a comment after a #
    for seq in seqs:
        seqs[seq] = re.split("#", seqs[seq])[0]
        
    return seqs

def safe_print(content, end='\n', outfile=sys.stdout):
    """
    workaround for having to run this currently under python 2 , and
    the lack of a flush option in the python 2 print function. This
    means that piping the output of print() statements (to e.g. head ) causes an exception:

    Example:
       ::

          print("\t".join([str(item) for item in result_tuple]))
          IOError: [Errno 32] Broken pipe
          close failed in file object destructor:
          sys.excepthook is missing
          lost sys.stderr

    """
    try:
        print(content, end=end, file=outfile)
        outfile.flush()
        sys.stderr.flush()
    except IOError as i:
        if i.errno != 32:
            raise
        else:
            outfile.flush()
            sys.stderr.flush()


def get_most_common_seqs(options, fastqc_seqs):
    # parse the log files to get these records 
    #assembled_by_distinct   CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTCGT  contained_in    CCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTCGT  counts= (231, 49)
    #assembled_by_distinct   CCCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTC   contained_in    CCCGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTC   counts= (149, 48)
    #assembled_by_distinct   GCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTCG   contained_in    GCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAACTGAGCGATCTCG   counts= (14, 48)
    #into a dictionary indexed by sequence, with count as value (count is first element of tuple)

    # we wil then search each fastq sequence in turn against all entries, an accumulate the counts , and list the entire fastqc set of seqs

    assembled_generators = []
    for filename in options["files"]:
        with open(filename,"r") as instream:
            assembled = (re.split("\t", record.strip()) for record in instream)
            assembled  = ((record[1], eval(record[5])[0]) for record in assembled if record[0] == "assembled_by_distinct")
            assembled_generators += assembled

    assembled_list = list(assembled_generators)
    #print(assembled_list) 

    #for record in itertools.chain(assembled_generators):
    #    print(record)
    fastqc_seq_counts = {}
    for name in fastqc_seqs:
        fastqc_seq_counts[name] = [ fastqc_seqs[name], 0]
        for (seq, count) in assembled_list:
            if re.search(fastqc_seqs[name], seq) is not None:
                fastqc_seq_counts[name][1] += 1

    fastqc_names = fastqc_seq_counts.keys()

    fastqc_names_sorted = sorted( fastqc_names , lambda n1,n2 : cmp( fastqc_seq_counts[n1][1], fastqc_seq_counts[n2][1]), reverse = True)

    for name in fastqc_names_sorted:
        safe_print("%s\t\t\t\t%s\t%s"%(name, fastqc_seq_counts[name][1], fastqc_seq_counts[name][0]))
        
def main():    
    options = get_options()
    fastqc_seqs = get_fastqc_seqs()

    get_most_common_seqs(options, fastqc_seqs) 
            
if __name__=='__main__':
    sys.exit(main())    

