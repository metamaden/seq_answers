# seq_answers
Sequencing data preprocessing

# Purpose
This repo is a knowledge base of preprocessing information and lessons learned in handing next generation sequencing (NGS) data. Some of the files here are recreated directly from certain resoures, such as the hg19 build shell script from bowtie 2, and all rights are preserved as stipulated by the original authors of those files. 

### What is Sequencing Data?

# Accessing and Handling Sequencing Data

### Platforms for Analysis of Sequencing Data
Many viable strategies abound for handling NGS data. The best strategy will depend on factors like your experimental design, comfort level with coding, and logistics of data access and manipulation. 

### Sequencing Pipelines
Various tools and methods will be added periodically to describe pipelines for preparing and analyzing sequencing data of various kinds.

### Sequencing Data Repositories
Sequence data can be stored in many [formats](https://galaxyproject.org/learn/datatypes/) (eg. fastq, fasta, bed, bam, sam, wig, bigwig, etc.  etc.) and can be found in many places online. Appendix 1 includes a list of many of the most commonly used public repositories for storing published sequenceing data. Several large public repositories that specialize in sequence data storage exist (SRA and EBA). Other repositories are more general but will occasionally include sequence data, typically in a more processed form (eg. [bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59395) and [bw](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) files on GEO).

### The Entrez System
The Entrez system, hosted and maintained by NCBI, is an interconnected system of repositories for biological and clinical research. It includes many commonly used databases, such as GEO, SRA, dbGap, PubMed, etc. Studies and study data are often indexed in many ways using a sophisticated informatics semantics that enables detailed and granular data storage and reference. 

To avoid getting overwhelmed, it's best to think not in terms of redundancies (eg. "there are many id's referring to the same information and I must memorize them all!"), but in terms of convenience (eg. "there is a distinct id for a distinct aspect of the data that will be particularly handy for certain bioinformatics applications, and less handy for other applications!").  

### Remote Access and Cloud Platforms, for Compute Performance
Galaxy is a powerful online platform, and it simplifies the handling of large amounts of data in memory by allocating around 1Tb of cloud space to the user, for free. Various iterations of the Galaxy platform are specially available for specific institutions, such as Fredhutch.io Galaxy. 

### Interfacing with the Sequence Read Archive and uploading to Galaxy
Files can be remotely acessed through a selection of repo's, or uploaded directly using FTP. Galaxy provides many useful resources for accessing data remotely from commonly used resources like the EBA and UCSC tables. However, some files may not be available through these resources. In these cases, it may be necesary to orchestrate some variation of the following in order to make a FASTQ file available to use on Galaxy:

1. Download the [SRA toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/). This is a set of functions to conveneient and remotely interface with various Entrez system databases, and especially the SRA. 
2. Use the fastq-dump command to deposit a desired FASTQ file locally or on a remove drive. You will need to identify the sample SRA ID (format ). 
3. Open a file transfer protocol (FTP) connection to the desired Galaxy server. You can use a program such as [FileZilla](https://filezilla-project.org/) to open a connection and transfer the files. To connect to your account server on the public Galaxy server, use the address "usegalaxy.org" and enter your account info. For large files or files not accessible using their remote database links, an FTP transfer is required.
4. Transfer FASTQ files to your Galaxy server directory.
5. Upload the FASTQ files to your Galaxy session (remember: specify a genome build in this step if you can! Otherwise analysis and visualization options will be limited or made more complicated).

# Coding, Scripting, and Programming
This section is intended for users more comfortable with computer programming, scripting, and command line work, and will overview steps and options for different types of pipelines and use cases for sequencing data.

### Order of Operations
Analysis involving sequencing data usually starts with an experimental design. Determine what sample groups you are comparing, what data is available, and the feasibility of carrying out your desired analysis. Especially in cases where you are re-analyzing public published data from multiple cohorts of samples, it is vital to plan out the experimental design before attempting to construct a workflow. In some cases, batch (in this case, experiment) will be confounded by factors such as treatment (eg. all normal samples are from experiment 1, and all disease samples are from experiment 2) or platform (eg. all normal samples were analyzed on an old sequencer platform, while disease samples were analyzed using a newer generation of that platform) or both treatment and platform!

A considerable amount of time can be spent determining whether preprocessed data is available, and whether it will suffice in an experiment that requires harmonized data between two studies. For instance, you might consider starting with preprocessed, pre-called peak bedfiles or bigwig files from GEO, only to find that files of the same format from a second study seem to show platform-specific differences (eg. experiment 2 samples were run on a different platform, and have 2-10x as many total peaks called, or a difference that we wouldn't expect to be accounted for by biological differences between samples in the experiments). Clearly, this cohort or batch effect will need to be reconciled before any cross-comparison of peak data between studies is possible.  

### General Considerations
Sequence data can be generated from many types of experiments, and the specific experiment conditions will likely dictate certain measures be taken in your workflow while excluding other possible measures. One important thing to note is that you will work with many types of sequence file formats, each containing different types of information and intended for different steps in the workflow. FASTQ is typically the least processed or most raw data form, while downstream analysis can be run with a myriad of other data types derived from the read file, including bed, bam, sam, wig, and bigwig. It is also important to note that the greater your knowledge a priori of the biological system studied, the better the workflow you can design. For instance, knowing you will be analyzing narrow peaks in a ChIP-Seq experiment will change the types of peak calling and analysis tools you will consider. 

### Early quality control
Read files in FASTQ or FASTA format should initially be inspected and analyzed using a quality control software, such as FASTQC and multiqc (both available in GALAXY, or through CRAN as an R package). This step enables an approximation of quality by basepair location across all reads in the samples. If reads have, for instance, systematic low quality at the start and end positions, it may be necessary to initially trim all of the reads. These initial QC measures can reduce the memory footprint of the files you work with, as well as the compute time for mapping and analysis.  

### Mapping
BWA, bowtie, bowtie2. With bowtie and other programs, an index genome must be built before an alignment can be performed. This could involve running a pre-written script, as shown in the misc directory of this repository, or running a command. 

### Peak Calling
Peak calling strategies depend on many factors, one of the most important being the anticipated sizes of the peaks. Narrow or small peaks should be called using different software from large or broad peaks, as different software is optimized for different run conditions. MACS, MACS2. Broad vs. narrow peak calling. 

### ChIP-Seq Preprocessing and Analysis Pipelines
Chromatin Immunoprecipitation followed by sequencing, or ChIP-Seq, has become a viable and important method of measuring various types of post-translational histone chemical modifications and enzymatic affinities, for the purpose of elucidating epigenetic mechanisms. Sequencing data generated in ChIP-Seq experiments are subjected to many of the same data handling and processing measures as sequencing data for RNAseq or other expression assay, with several important exceptions. First, ChIP-Seq data can a priori be anticipated to yield small or narrow peaks of overlapping reads, versus wide or broad overlapping read peaks. This difference changes the type of math and statistics that should be applied in the ChIP-Seq workflow. Another consideration is the potential for bias in read coverage due to open versus closed chromatin states. This potential for bias makes it even more important for ChIP-Seq experiments to include DNA input controls for tested samples, which enables control peak signal to be subtracted from the sample peaks (which include signal and noise), to increase the signal to noise ratio of the dataset. 

### Assembling Your Own Workflow for Sequencing Data

# Appendices
## Apendix 1: List of Public Sequence Data Repositories and Systems
Databases and Repositories:
1. [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) - resource for cataloging study, sample, and data platform data. Occasionally contains sequence data, but usually in a preprocessed or pre-filtered format (eg. bedfiles, bw, etc.). Handy modules and tools for remote access are available. 
2. [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) - resource for cataloging sequence data samples and meta-data, with occasional browser reference to a FASTQ file (if a FASTQ is unavailable through the SRA browser, check ENA, described below, or use the fastq-dump utility in the Entrez toolkit, described above).
3. [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena) - resource for published un-preprocessed FASTQ files, with handy access through Galaxy browser. 
4. [UCSC Genome Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?GALAXY_URL=https%3A//usegalaxy.org/tool_runner&tool_id=ucsc_table_direct1&hgta_compressType=none&sendToGalaxy=1&hgta_outputType=bed) - resource for mainly bed, bam, and wig-type files of sequence data commonly used for visualization in ideogram track format. 

Systems: 
1. The Entrez System - engine powering cross-database queries for NCBI and NLM ([more info](https://en.wikipedia.org/wiki/Entrez))
2. [Galaxy browser](https://usegalaxy.org/) - excellent and free cloud-based platform with pre-loaded utilities for preprocessing, analyzing, and visualizing all manner of sequence and array data. Considerable amount of free cloud storage and memory space for the user, with reasonable speed and highly convenient integrated database cross-talk. ([info](https://galaxyproject.org/learn/), [videos](https://vimeo.com/galaxyproject), [news](https://galaxyproject.org/news/))
3. [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) - toolkit for conveniently and remotely interfacing with Entrez system databases ([wiki](https://github.com/ncbi/sra-tools/wiki)), and especially for remotely downloading FASTQ files using sample SRA IDs ([info about the fastq-dump command](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump)).
