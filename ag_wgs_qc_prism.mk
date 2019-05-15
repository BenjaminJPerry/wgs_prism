# ag_gbs_qc_prism.mk prism main makefile
#***************************************************************************************
# references:
#***************************************************************************************
# make: 
#     http://www.gnu.org/software/make/manual/make.html
#

########## non-standard analysis - these not (currently) part of "all" as expensive
%.annotation:   %.blast_analysis
	$@.sh > $@.mk.log 2>&1
	date > $@

%.blast_analysis:   %.fasta_sample
	$@.sh > $@.mk.log 2>&1
	date > $@

########## standard analysis 
%.all:  %.kmer_analysis %.fastqc  
	date > $@

%.kmer_analysis:   %.fasta_sample
	$@.sh > $@.mk.log 2>&1
	date > $@

%.fastqc:   %.bcl2fastq
	$@.sh > $@.mk.log 2>&1
	date > $@

%.fasta_sample:   %.bcl2fastq
	$@.sh > $@.mk.log 2>&1
	date > $@

%.fastq_sample:   %.bcl2fastq
	$@.sh > $@.mk.log 2>&1
	date > $@

%.bcl2fastq:
	$@.sh > $@.mk.log 2>&1
	date > $@

.PHONY: %.clean
%.clean: 
	$@.sh > $@.mk.log 2>&1


##############################################
# specify the intermediate files to keep 
##############################################
.PRECIOUS: %.log %.ag_wgs_qc_prism %.blast_analysis %.kmer_analysis  %.bcl2fastq %.fasta_sample %.fastq_sample %.fastqc %.annotation
