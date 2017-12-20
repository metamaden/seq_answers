This is a short guide on how the analysis for "A comparison of
control samples for ChIP-seq of histone modifications" can be
reproduced using the attached .R scripts.

#######################################################
#
# Data pre-processing outside R:
#
#######################################################
- Align all reads. (bowtie2 for ChIP-seq and tophat for
RNA-seq).
- set the output directory for tophat to the aligned fastq file. the R
#script will look for files of the form originalFastqName.fastq/accepted_hits.bam.
- Transform the ChIP .bams to .bed files.
- It is important to keep the body of the original fastq file names
#for both ChIP and RNA data, or it may not import properly into R.
- Merge ChIP replicates and run MACS on all groups, and on H3K27me3 with
WCE and H3 as control. Use the following names (-n option in macs) for
#the output to match import functions in R:
	- WCE
	- H3
	- K27
	- K27overWCE (with WCE as control)
	- K27overH3 (with merged H3 as control)

#######################################################
#
# The rest of the analysis is done in R:
#
#######################################################
- Set the path to input directories and desired output directory in
  the first few lines of the script "input.R".
- Make sure that you have the following R packages installed:
       - rtracklayer
       - GenomicRanges
       - parallel (can be avoided by setting cpus=1 in the
       downsample(...) calls in input.R.) 
       - limma
       - edgeR
       - Rsubread
       - org.Mm.eg.db
       - lme4
- run the script "input.R" in R in a place where you can leave it
overnight. The script will report as it progresses.
- The script will produce the plots in the publication.
- If you wish to perform further analysis, the script leaves all the
needed information in the global environment. Some of the more useful
variables are:
	  - cov100, cov1k: matrices of ChIP-seq read counts in bins of 100 and
	  1000 bp.
	  - geneCounts, promCounts: matrices of ChIP-seq read counts
	  in gene body and promoter. 
	  - rpkm: expression of the genes in RPKM.
	  - genes: GRanges of genes.
	  - chipfitGenes, chipfitBins: limma-voom "fit" objects from the
	  differential analysis of counts over genes and bins.
	  - WCEpeaks, H3peaks, K27peaks: GRanges of the MACS peaks.
	  - K27overWCEpeaks, K27overH3peaks: GRanges of the MACS peaks
	  with control.
- There are some potentially useful analysis functions in "inputFunctions.R".
- Note that most of the preprocessing in R is saved to file, so if you want
#to rerun the script, it will be significantly faster.
