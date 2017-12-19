# seq_answers
Sequencing data preprocessing

# Purpose
This repo is a knowledge base of preprocessing information and lessons learned in handing NGS sequencing data. Some of the files here are recreated directly from certain resoures, such as the hg19 build shell script from bowtie 2, and all rights are preserved as stipulated by the original authors of those files. 

# Sequencing Pipelines
Various tools and methods will be added periodically to describe pipelines for preparing and analyzing sequencing data of various kinds.

# Platforms
Many viable strategies abound for handling NGS data. The best strategy will depend on factors like your experimental design, comfort level with coding, and logistics of data access and manipulation. 

# Repositories
Sequence data can be stored in many [formats](https://galaxyproject.org/learn/datatypes/) (eg. fastq, fasta, bed, bam, sam, wig, bigwig, etc.  etc.) and can be found in many places online. Appendix 1 includes a list of many of the most commonly used public repositories for storing published sequenceing data. Several large public repositories that specialize in sequence data storage exist (SRA and EBA). Other repositories are more general but will occasionally include sequence data, typically in a more processed form (eg. [bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59395) and [bw](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) files on GEO).

## The Entrez System
The Entrez system, hosted and maintained by NCBI, is an interconnected system of repositories for biological and clinical research. It includes many commonly used databases, such as GEO, SRA, dbGap, PubMed, etc. Studies and study data are often indexed in many ways using a sophisticated informatics semantics that enables detailed and granular data storage and reference. To avoid getting overwhelmed, it's best to think not in terms of redundancies (eg. "there are many id's referring to the same information and I must memorize them all!"), but in terms of convenience (eg. "there is a distinct id for a distinct aspect of the data that will be particularly handy for certain bioinformatics applications, and less handy for other applications!").  

# Remote/Cloud 
Galaxy is a powerful online platform, and it simplifies the handling of large amounts of data in memory by allocating around 1Tb of cloud space to the user, for free. Various iterations of the Galaxy platform are specially available for specific institutions, such as Fredhutch.io Galaxy. 

## Interfacing with the Sequence Read Archive and uploading to Galaxy
Files can be remotely acessed through a selection of repo's, or uploaded directly using FTP. Galaxy provides many useful resources for accessing data remotely from commonly used resources like the EBA and UCSC tables. However, some files may not be available through these resources. In these cases, it may be necesary to orchestrate some variation of the following in order to make a FASTQ file available to use on Galaxy:

1. Download the [SRA toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/). This is a set of functions to conveneient and remotely interface with various Entrez system databases, and especially the SRA. 
2. Use the fastq-dump command to deposit a desired FASTQ file locally or on a remove drive. You will need to identify the sample SRA ID (format ). 
3. Open a file transfer protocol (FTP) connection to the desired Galaxy server. You can use a program such as [FileZilla](https://filezilla-project.org/) to open a connection and transfer the files. To connect to your account server on the public Galaxy server, use the address "usegalaxy.org" and enter your account info. For large files or files not accessible using their remote database links, an FTP transfer is required.
4. Transfer FASTQ files to your Galaxy server directory.
5. Upload the FASTQ files to your Galaxy session (remember: specify a genome build in this step if you can! Otherwise analysis and visualization options will be limited or made more complicated).

# Coding, Scripting, and Programming

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
