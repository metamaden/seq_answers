# Purpose: Make a standard preprocessing pipeline for ChIP-Seq data run on an arbitrary NGS sequencer 

#install.packages("fastqcr")
library(fastqcr)

# navigate to folder containing fastq files to anlayze
fastqc()
qca <- qc_aggregate("./FASTQC")
qcf <- qc_fails(qca)
