#!/usr/bin/env python


#
# this script creates an html document giving a mash-up of
# wgs production 
#
import os
import re
import itertools
import string
import exceptions
import argparse


header1="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
   "httpd://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
<head>
<title>
Overview of %(run_name)s
</title>
</head>
<body>
<h1> Q/C Summaries for <a href="http://agbrdf.agresearch.co.nz/cgi-bin/fetch.py?obid=%(run_name)s&context=default">%(run_name)s</a> </h1>
"""

overview_section="""
<p/>
<table width=90%% align=center>
<tr id=bcl2fastq>
<td> bcl2fastq reports  </td>
<td> <a href=bcl2fastq/index.html> bcl2fastq reports </a>  </td>
</tr>
</table>
<hr/>
<p/>
"""


footer1="""
</body>
</html>
"""

BASEDIR="/dataset/gseq_processing/scratch/illumina/hiseq"



def get_samples(options):
    # samples  are idenitified as subfolders of the run folder 
    html_folder=os.path.join(BASEDIR, options["run_name"] , "html")
    #print "DEBUG : "+run_folder

    sample_folders=[ node for node in os.listdir(html_folder) if re.search("^tardis", node) is None and  \
                     re.search("bcl2fastq", node) is None and \
                     os.path.isdir(os.path.join(html_folder, node)) ]
    
    return sample_folders
    
    
def generate_run_plot(options):
    stats = {
        "found file count" : 0,
        "no file count" : 0,
        "no sample count" : 0
    }

    file_group_iter = ( ("fastqc", "link"),\
                       ("kmer summary (plots)", "image"), ("text kmer summary (clusters)", "link"),\
                       ("Preview common sequence", "in-line"), ("All common sequence", "link"), \
                       )
    file_iters = {
        #"KGD" : ['KGD/MAFHWdgm.05.png', 'KGD/SNPDepthHist.png', 'KGD/AlleleFreq.png', 'KGD/GHWdgm.05-diag.png', 'KGD/SNPDepth.png', 'KGD/finplot.png', 'KGD/Heatmap-G5HWdgm.05.png', 'KGD/SampDepth.png', 'KGD/G-diag.png', 'KGD/Gdiagdepth.png', 'KGD/LRT-hist.png', 'KGD/MAF.png', 'KGD/GcompareHWdgm.05.png', 'KGD/Gcompare.png', 'KGD/SampDepthHist.png', 'KGD/CallRate.png', 'KGD/GHWdgm.05diagdepth.png', 'KGD/Heatmap-G5.png', 'KGD/SampDepth-scored.png', 'KGD/HWdisMAFsig.png', 'KGD/LRT-QQ.png', 'KGD/SampDepthCR.png', 'KGD/PC1v2G5HWdgm.05.png'],
        #"KGD plots" : ['KGD/AlleleFreq.png', 'KGD/finplot.png', 'KGD/G-diag.png', 'KGD/HWdisMAFsig.png', 'KGD/MAF.png', 'KGD/SampDepth.png', 'KGD/SNPDepth.png',
        #        'KGD/CallRate.png', 'KGD/GcompareHWdgm.05.png', 'KGD/GHWdgm.05diagdepth.png', 'KGD/LRT-hist.png', 'KGD/PC1v2G5HWdgm.05.png', 'KGD/SampDepth-scored.png'
        #        'KGD/Co-call-HWdgm.05.png', 'KGD/Gcompare.png', 'KGD/GHWdgm.05-diag.png', 'KGD/LRT-QQ.png', 'KGD/SampDepthCR.png', 'KGD/SNPCallRate.png'
        #        'KGD/Co-call-.png', 'KGD/Gdiagdepth.png', 'KGD/Heatmap-G5HWdgm.05.png', 'KGD/MAFHWdgm.05.png', 'KGD/SampDepthHist.png', 'KGD/SNPDepthHist.png'],
        "fastqc" : ["fastqc"],
        "Preview common sequence" : [ 'common_sequence/preview_common_sequence.txt']            ,
        "All common sequence" : [ 'common_sequence/all_common_sequence.txt']            ,        
        "kmer summary (plots)" : [ 'kmer_analysis/kmer_entropy.k6A.jpg', 'kmer_analysis/kmer_entropy_plus.k6A.jpg', 'kmer_analysis/kmer_zipfian_comparisons.k6A.jpg','kmer_analysis/zipfian_distances.k6A.jpg']            ,
        "text kmer summary (clusters)" : [ 'kmer_analysis/heatmap_sample_clusters.k6A.txt', 'kmer_analysis/heatmap_sample_clusters_plus.k6A.txt','kmer_analysis/zipfian_distances_fit.k6A.txt']        
    }

    
    with open(options["output_filename"],"w") as out_stream:

        print >> out_stream, header1%options

        print >> out_stream, overview_section

        print >> out_stream, """
Contents:
</p>
<ul>
"""


        for (file_group, file_type)  in file_group_iter:
            print >> out_stream, "<li> <a href=#%s> %s </a>"%(file_group, file_group)
        print >> out_stream, """
</ul>
<hr/>
</p>
"""
        #print "DEBUG : calling get_samples"
        samples = get_samples(options)

        for (file_group, file_type)  in file_group_iter:
            # output overall header for file group
            #print >> out_stream, "<h2> %s </h2>\n"%file_group
            # output sample column headings
            print >> out_stream, "<table width=90%% align=left id=%s>\n"%file_group, \
                                 "<tr>\n", \
                                 "<td> <h4> Sample name </h4> </td>\n",\
                                 "\n".join(["<td><h4> %s </h4></td>"%sample for sample in samples]), \
                                 "\n</tr>\n"
            for file_name in file_iters[file_group]:

                print >> out_stream , "<tr><td>%s</td>\n"%file_name
                for sample in samples:
                    file_path = os.path.join(BASEDIR, options["run_name"], "html", sample, file_name)

                    if file_type == "image":
                        image_relpath=os.path.join(sample, file_name)

                        if os.path.exists(file_path):
                            print >> out_stream, "<td> <img src=%s title=%s height=%s width=%s/> </td>\n"%(image_relpath, file_path, options["image_height"], options["image_width"])
                        else:
                            print >> out_stream, "<td> %s unavailable </td>\n"%file_path
                    elif file_type == "link":
                        link_relpath=os.path.join(sample, file_name)

                        if os.path.exists(file_path):
                            print >> out_stream, "<td width=300 align=left> <a href=%s target=%s> %s </a></td>\n"%(link_relpath, file_name, link_relpath)
                        else:
                            print >> out_stream, "<td width=300> %s unavailable </td>\n"%file_path
                    elif file_type == "in-line":
                        text = "(unavailable)"
                        if os.path.exists(file_path):
                            with open(file_path,"r") as infile:
                                text="\n".join((record.strip() for record in infile))
                                
                        print >> out_stream, "<td id=\"%s\"> <font size=-2> <pre>%s</pre> </font> </td>"%(file_group,text)
                            
                print >> out_stream , "</tr>\n"
            print >> out_stream, "</table>\n"



        print >> out_stream, footer1



    print stats
                
                
def get_options():
    description = """
    """
    long_description = """
example :

./make_sample_pages.py -r 190510_D00390_0457_AHYFCVBCX2 -o /dataset/gseq_processing/scratch/illumina/hiseq/190510_D00390_0457_AHYFCVBCX2/html/peacock.html
    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--run_name' , dest='run_name', required=True, type=str, help="run name")
    parser.add_argument('-H', '--image_height' , dest='image_height', default=800, type=int, help="image height")
    parser.add_argument('-W', '--image_width' , dest='image_width', default=800, type=int, help="image width")
    parser.add_argument('-o', '--output_filename' , dest='output_filename', default="peacock.html", type=str, help="name of output file")

    
    args = vars(parser.parse_args())

    return args


def main():

    options = get_options()
    print options 
    generate_run_plot(options)

    
if __name__ == "__main__":
   main()



        

