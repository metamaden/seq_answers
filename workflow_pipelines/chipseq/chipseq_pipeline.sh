
# ChIP-Seq data pipeline for narrow peaks, single reads

#================================================================================================================================
#
# The below example shows preprocessing, QC and analysis of a cohort of normal colon crypts from cancer patients
# In this example, we will process H3K27ac and DNA input control samples together
# FASTQ files are available through Sequence Read Archive, and can be obtained using fastq-dump command from the SRA toolkit
#
#================================================================================================================================

#
# Pipeline operations
#
# 0. Module load bioinformatics programs
# for Fred Hutch gizmo clusters, note the list of available programs (stored here: https://github.com/metamaden/seq_answers/blob/master/misc/gizmo_programs_list.txt)

ml bowtie2
ml samtools
ml bedtools
ml MACS2/2.1.0.20151222-foss-2015b-Python-2.7.9

#
# 1. Obtain FASTQ Files from SRA using SRA toolkit 
#

# Details: use fastq-dump from SRA toolkit (must have sra toolkit installed and in current dir or PATH). 
# see FASTQC and scripts in this repo for QC instructions

# test samples (colon crypts)
fastq-dump SRR3157776 # C28 H3K27ac ChIP Seq
fastq-dump SRR3157777 # C29 H3K27ac ChIP Seq;
fastq-dump SRR3157778 # C37 H3K27ac ChIP Seq
fastq-dump SRR3157779 # Crypt5 H3K27ac ChIP Seq

# control samples (colon crypts)
fastq-dump SRR3157849 # C28 input DNA
fastq-dump SRR3157850 # C29 input DNA
fastq-dump SRR3157851 # C37 input DNA
fastq-dump SRR3157852 # Crypt5 input DNA

#
# 2. Map FASTQ files using bowtie 2
#
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157776.fastq > SRR3157776.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157777.fastq > SRR3157777.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157778.fastq > SRR3157778.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157779.fastq > SRR3157779.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157849.fastq > SRR3157849.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157850.fastq > SRR3157850.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157851.fastq > SRR3157851.sam
bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157852.fastq > SRR3157852.sam

# example Terminal output from mapping:
#>  maden@gizmoh2:/fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017$ bowtie2 -p 8 -x /fh/fast/grady_w/ChIP_Seq/bowtie2-2.3.3.1/indexes/hg19 SRR3157776.fastq > SRR3157776.sam
#>  69387874 reads; of these:
#>    69387874 (100.00%) were unpaired; of these:
#>    4442882 (6.40%) aligned 0 times
#>  53606999 (77.26%) aligned exactly 1 time
#>  11337993 (16.34%) aligned >1 times
#>  93.60% overall alignment rate

#
# 3. File conversion (for C29 test and control samples)
#

samtools view -bS /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157777.sam > /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157777.bam
samtools view -bS /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157777.sam > /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157777.bam

#
# 4. Peak calling with MACS2 (for C29, using test and control samples)
#
# EDIT: macs2 -t /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157777.bed -c /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/SRR3157850.bed -f BED -g '2.7e9' -n /fh/fast/grady_w/ChIP_Seq/fastq/cohen_2017/C29_h3k27ac_macs2
