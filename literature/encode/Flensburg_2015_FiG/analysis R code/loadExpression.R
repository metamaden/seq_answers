#this file loads the RNA expression data, remaking the featureCounts call if needed.
#output is stored in global variable 'RNAcounts' with symbol gene names as rownames.
#also fC is a featureCounts output object, containing information about annotation etc.



#prepare the files
expressionCountFiles = c(
    '1_LSK-GFP-nons_BH0AABADXX_ATCACG_L001_R1.fastq/accepted_hits.bam.fC',
    '1_LSK-GFP-nons_BH0AABADXX_ATCACG_L002_R1.fastq/accepted_hits.bam.fC',
    '1_LSK-GFP_nons_BH0AABADXX_CGATGT_L001_R1.fastq/accepted_hits.bam.fC',
    '1_LSK-GFP_nons_BH0AABADXX_CGATGT_L002_R1.fastq/accepted_hits.bam.fC',
    '2_Nons_LSK_GFP__1_Nons_LSK_GFP_BH0AABADXX_CAGATC_L001_R1.fastq/accepted_hits.bam.fC',
    '2_Nons_LSK_GFP__1_Nons_LSK_GFP_BH0AABADXX_CAGATC_L002_R1.fastq/accepted_hits.bam.fC')


if ( !filesInPlace(expressionDirectory, expressionCountFiles)) {
  cat('Did not find featureCount files. Remaking count files.\n')
  require(Rsubread)

  bamFiles = paste0(expressionDirectory, gsub('.fC', '', expressionCountFiles))
  if ( !filesInPlace(mouseExpressionDirectory, bamFiles) ) {
    stop('Could not find featureCount files or bam files in ', expressionDirectory, '.\n')
  }
  for ( file in bamFiles ) {
    fC = featureCounts(file, annot.inbuilt='mm10', nthreads=5)
    save(fC, file=paste0(file, '.fC'))
  }
}

#load the LSK data
cat('Loading expression data...')
load(paste0(expressionDirectory, expressionCountFiles[1]))
counts1.1 = fC$counts
load(paste0(expressionDirectory, expressionCountFiles[2]))
counts1.2 = fC$counts
load(paste0(expressionDirectory, expressionCountFiles[3]))
counts2.1 = fC$counts
load(paste0(expressionDirectory, expressionCountFiles[4]))
counts2.2 = fC$counts
load(paste0(expressionDirectory, expressionCountFiles[5]))
counts3.1 = fC$counts
load(paste0(expressionDirectory, expressionCountFiles[6]))
counts3.2 = fC$counts

#featureCount gives entrez gene names, but we want them in symbol. This function translates.
translateMouseGeneNamesFromEntrez = function(genes) {
  if(require("org.Mm.eg.db")){
    cat("org.Mm.eg.db loaded...\n")
  } else {
    cat("installing org.Mm.eg.db...\n")
    source("http://bioconductor.org/biocLite.R")
    biocLite("org.Mm.eg.db")
    if(require(lme4)){
        cat("org.Mm.eg.db installed and loaded...\n")
    } else {
        stop("could not install org.Mm.eg.db...\n")
    }
  } 
  
  # Get the mapping from Entrez 2 SYMBOL
  eg2sym=data.frame(symbol=unlist(as.list(org.Mm.egSYMBOL)))
 
  #return the translated genes in vector format
  return(as.vector(eg2sym[genes,]))
}

#keep only genes located on chromosome 1-19, X, Y, M.
use = gsub('chr','',gsub(';.*', '', fC$annotation$Chr)) %in% names(mouseChrLengths())

#merge technical replicates and summarise in matrix with symbol gene names.
RNAcounts = cbind(counts1.1 + counts1.2, counts2.1 + counts2.2, counts3.1 + counts3.2)
RNAnames = c(1, 2, 3)
rownames(RNAcounts) = fC$annotation$GeneID
rownames(RNAcounts) = translateMouseGeneNamesFromEntrez(rownames(RNAcounts))
colnames(RNAcounts) = RNAnames
RNAcountsOri = RNAcounts
RNAcounts = RNAcounts[use,]

exonLength = fC$annotation$Length[use]
expression = (rowSums(RNAcounts)/sum(RNAcounts) + 1e-6)/exonLength
rpkm = 0.01 + (0.5+rowSums(RNAcounts)/sum(RNAcounts)*1e6)/(exonLength/1000)


cat('Done loading expression data.\n')
