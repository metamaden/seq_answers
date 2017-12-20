#runs on all .bed files in the directory.
#Takes long time to read in the .bed files, up to or above an hour for a 50Mread file.
#creates reads files, and coverage files binned on 100 and 1000 bp bins.
#all you need to do with the mouse chip data for downstream analysis.
preprocessMouseChip = function(dir) {
  bedToR(dir)
  readsToCoverage(dir, binSize=1000)
  fillUpCoverageToSameLength(dir, binSize=1000)
  readsToCoverage(dir, binSize=100)
  fillUpCoverageToSameLength(dir, binSize=100)
}

#imports the coverages of provided bin size in the directory to file.
getAllMouseCovs = function(dir, binSize = 1000) {
  pattern = paste('*.Cov', binSize,  '.Rdata$', sep='')
  CovFiles = list.files(dir, pattern = pattern)
  cat('found files:\n')
  for ( i in 1:length(CovFiles) )
    cat(CovFiles[[i]], '\n')

  cat('Loading data.\n')
  for ( i in 1:length(CovFiles) )
    load(paste(dir, CovFiles[[i]], sep=''))

  covs = list()
  
  for ( i in 1:length(CovFiles) ) {
    varName = gsub('.Rdata', '', CovFiles[[i]])
    cat('Adding ', varName, '.\n', sep='')
    covList = get(varName)
    sortedList = lapply(names(mouseChrLengths()), function(chr) covList[[chr]])
    covs[[i]] = unlist(sortedList)
  }

  cat('merging to matrix...')
  ret = do.call(cbind, covs)
  cat('done.\n')

  colnames(ret) = gsub(pattern, '', CovFiles)

  return(ret)
}

#returns the length of the mouse chromosomes
#this function is updated after publication. previously it was rounded up to closest 100kbp.
mouseChrLengths = function() {
  lengths = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993,
    122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299)

  names(lengths) = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10' ,'11', '12' ,'13', '14' ,'15', '16' ,'17', '18', '19', 'X', 'Y', 'M')

  return(lengths)
}

#This reads all bed files in a directory to R format, and saves the data to files.
bedToR = function( dir ) {
  require(rtracklayer)
  bedFiles = list.files(dir, pattern = '*.bed$')
  #filter out MACS peak files
  bedFiles = bedFiles[!grepl('_peaks.bed$', bedFiles)]
  bedFiles = bedFiles[!grepl('_summits.bed$', bedFiles)]
  cat('Found files:\n')
  lapply(1:length(bedFiles), function(i) cat(bedFiles[i],'\n') )
  
  readAndSave = function(i) {
    cat('collecting garbage.\n')
    gc()
    cat('starting file ', i, ' of ', length(bedFiles), '.\n')
    varName = gsub('.*?\\.Input','Input',
      gsub('.*?\\.H3', 'H3',
           gsub('.bed', '', bedFiles[i])))
    varName = paste(varName, '.reads', sep='')
    fileName = paste(dir, bedFiles[i], sep='')
    cat('Importing file ', fileName, '.\n', sep='')
    assign(varName, import(fileName))
    cat('Found ', length(get(varName)), ' reads. ', sep='')
    saveFile = paste(dir, varName, '.Rdata', sep='')
    cat('Saving to ', saveFile, '.\n', sep='')
    save(list=varName, file = saveFile)
  }
  
  lapply(1:length(bedFiles), readAndSave )
  return(0)
}

#help function to readsToCoverage below.
getHitsList = function(reads, binsize=1000, chrLengths = 0, minScore=20) {
  scores = score(reads)
  keep = scores >= minScore
  cat('Keeping', sum(keep), 'out of', length(keep), 'reads due to quality.\n')
  
  chr = seqnames(reads)[keep]
  mid = (start(reads)[keep] + end(reads)[keep])/2

  chrNames = levels(chr)
  midByChr = lapply(chrNames, function(chrName) mid[as.logical(chr == chrName)])

  nBins = lapply(midByChr, function(mid) ceiling(max(mid)/binsize))

  hitCounts = function(mid) {
    if ( length(mid) == 0 ) return(0)
    h = hist(mid, breaks = seq(from=0, to=max(mid)+binsize, by=binsize), plot=F)
    return(h$counts)
  }
  hitsList = lapply(midByChr, hitCounts)
  names(hitsList) = gsub('MT', 'M', gsub('chr', '', chrNames))
  return(hitsList)
}

#takes a firectory, finds all the files with saved reads in R format, and bins the read density.
#saves the result to file.
readsToCoverage = function(dir, binSize = 1000) {
  RFiles = list.files(dir, pattern = '*.reads.Rdata')
  if ( sum(grepl('Cov', RFiles)) > 0 )
    RFiles = RFiles[-grep('Cov', RFiles)]
  if ( sum(grepl('Ns', RFiles)) > 0 )
    RFiles = RFiles[-grep('Ns', RFiles)]
  if ( sum(grepl('RNA', RFiles)) > 0 )
    RFiles = RFiles[-grep('RNA', RFiles)]
  cat('found files:\n')
  for ( i in 1:length(RFiles) )
    cat(RFiles[[i]], '\n')

  readAndSave = function(i) {
    cat('collecting garbage.\n')
    gc()
    
    cat('starting file ', i, ' of ', length(RFiles), '.\n')
    varName = gsub('.Rdata', '', RFiles[i])
    fileName = paste(dir, RFiles[i], sep='')
    
    cat('Loading reads ', fileName, '.\n', sep='')
    load(fileName)
    newVarName = paste(varName, '.Cov', binSize, sep='')
    
    cat('Calculating hits list.\n')
    assign(newVarName, getHitsList(get(varName), binSize))
    saveFile = paste(dir, newVarName, '.Rdata', sep='')
    
    cat('Saving to ', saveFile, '.\n', sep='')
    save(list=newVarName, file = saveFile)
  }
  
  lapply(1:length(RFiles), readAndSave )
  return(0)
}

#goes through coverage files and fills up with empty bins to the end of each cromosome.
fillUpCoverageToSameLength = function(dir, binSize = 1000) {
  CovFiles = list.files(dir, pattern = paste0('*.Cov', binSize, '.Rdata'))
  cat('found files:\n')
  for ( i in 1:length(CovFiles) )
    cat(CovFiles[[i]], '\n')

  chrLengths = ceiling(mouseChrLengths()/binSize)
  
  cat('Loading data.\n')
  for ( i in 1:length(CovFiles) )
    load(paste(dir, CovFiles[[i]], sep=''))

  for ( i in 1:length(CovFiles) ) {
    varName = gsub('.Rdata', '', CovFiles[[i]])
    cat('Fixing ', varName, '.\n', sep='')
    cov = get(varName)
    for ( chr in names(chrLengths) ) {
      if ( chrLengths[[chr]] > length(cov[[chr]]) ) {
        cat('Extending chr ', chr, ' from length ', length(cov[[chr]]), sep='')
        cov[[chr]] = c(cov[[chr]], rep(0, chrLengths[[chr]] - length(cov[[chr]])))
        cat(' to ', length(cov[[chr]]), '.\n', sep='')
      }
      else if ( chrLengths[[chr]] < length(cov[[chr]]) ) {
        cat('Trunkating chr ', chr, ' from length ', length(cov[[chr]]), sep='')
        lost = sum(cov[[chr]][(chrLengths[[chr]]+1):length(cov[[chr]])])
        cov[[chr]] = cov[[chr]][1:chrLengths[[chr]]]
        cat(' to ', length(cov[[chr]]), ', losing ',  lost, ' reads.\n', sep='')
      }
    }
    saveFile = paste(dir, CovFiles[[i]], sep='')
    cat('Savefile is ', saveFile, '.\n', sep='')
    assign(varName, cov)
    save(list=varName, file = saveFile)
  }

  return(0)
}

#calculates the counts over the gene or promoter.
promoterCount = function(cov, genes, binSize, species='mouse', range=2000) {
  mid = ceiling(genesToStartX(genes, species=species)/binSize)
  start = mid - round(range/binSize)
  end = mid + round(range/binSize)
  count = sapply(1:length(start), function(gene) sum(cov[start[gene]:end[gene]]))
  return(count)
}
geneCount = function(cov, genes, binSize, species='mouse') {
  start = ceiling(genesToStartX(genes, species=species)/binSize)
  end = ceiling(genesToEndX(genes, species=species)/binSize)
  count = sapply(1:length(start), function(gene) sum(cov[start[gene]:end[gene]]))
  return(count)
}

#calculates the enrichment (nreads/base/libSize (aka RPKM)) over the gene or promoter.
promoterEnrichment = function(cov, genes, binSize, species='mouse') {
  mid = ceiling(genesToStartX(genes, species=species)/binSize)
  start = mid - round(2000/binSize)
  end = mid + round(2000/binSize)
  count = sapply(1:length(start), function(gene) sum(cov[start[gene]:end[gene]]))
  nRpBL = 1e9*(count/sum(cov))/4000
  return(nRpBL)
}
geneEnrichment = function(cov, genes, binSize, species='mouse') {
  start = ceiling(genesToStartX(genes, species=species)/binSize)
  end = ceiling(genesToEndX(genes, species=species)/binSize)
  count = sapply(1:length(start), function(gene) sum(cov[start[gene]:end[gene]]))
  nRpBL = 1e9*(count/sum(cov))/((1+abs(end-start))*binSize)
  return(nRpBL)
}

#changes base from chr and bp to a single number x count bp from start of chr 1 through to chr Y.
genesToStartX = function( genes, species='mouse' ) {
  prevChrL = c(0, cumsum(mouseChrLengths()))
  names(prevChrL) = c(names(mouseChrLengths()), 'outside')
  chr = as.character(gsub('chr','',seqnames(genes)))
  if ( species == 'Arabidopsis' ) {
    prevChrL = c(0, cumsum(flowerChrLengths()))
    names(prevChrL) = c(names(flowerChrLengths()), 'outside')
    chr = as.character(seqnames(genes))
  }
  start = start(genes)
  end =  end(genes)
  x = ifelse(strand(genes) == '+', start, end) + prevChrL[chr]
  return(x)
}
genesToEndX = function( genes, species='mouse' ) {
  prevChrL = c(0, cumsum(mouseChrLengths()))
  names(prevChrL) = c(names(mouseChrLengths()), 'outside')
  chr = as.character(gsub('chr','',seqnames(genes)))
  if ( species == 'Arabidopsis' ) {
    prevChrL = c(0, cumsum(flowerChrLengths()))
    names(prevChrL) = c(names(flowerChrLengths()), 'outside')
    chr = as.character(seqnames(genes))
  }
  start = start(genes)
  end =  end(genes)
  x = ifelse(strand(genes) == '+', end, start) + prevChrL[chr]
  return(x)
}


#Reads genes from file and converts to granges
getGenes = function(filename) {
  require(rtracklayer)
  #import file
  cat("importing genes from ",filename,"...")
  genes = import(filename)
  cat("done.\n")

  return(genes)
}

#check if files already exist
mouseChipFilesInPlace = function(directory) {
  requiredFiles = c(
    'ChIP1_LSK_WCE_C2KJYACXX_CGATGT_L006_R1.reads.Cov100.Rdata',
    'ChIP2_LSK_H3_1_C2KJYACXX_ACAGTG_L006_R1.reads.Cov100.Rdata',
    'ChIP3_LSK_H3_2_C2KJYACXX_GCCAAT_L006_R1.reads.Cov100.Rdata',
    'ChIP4_LSK_H3K27me3_1_C2KJYACXX_CAGATC_L006_R1.reads.Cov100.Rdata',
    'ChIP5_LSK_H3K27me3_2_C2KJYACXX_CTTGTA_L006_R1.reads.Cov100.Rdata',
    'ChIP6_LSK_H3K27me3_3_C2KJYACXX_GTGAAA_L006_R1.reads.Cov100.Rdata',
    'ChIP1_LSK_WCE_C2KJYACXX_CGATGT_L006_R1.reads.Cov1000.Rdata',
    'ChIP2_LSK_H3_1_C2KJYACXX_ACAGTG_L006_R1.reads.Cov1000.Rdata',
    'ChIP3_LSK_H3_2_C2KJYACXX_GCCAAT_L006_R1.reads.Cov1000.Rdata',
    'ChIP4_LSK_H3K27me3_1_C2KJYACXX_CAGATC_L006_R1.reads.Cov1000.Rdata',
    'ChIP5_LSK_H3K27me3_2_C2KJYACXX_CTTGTA_L006_R1.reads.Cov1000.Rdata',
    'ChIP6_LSK_H3K27me3_3_C2KJYACXX_GTGAAA_L006_R1.reads.Cov1000.Rdata'
    )
  files = list.files(directory, pattern='*.Rdata$')
  return( sum(!(requiredFiles %in% files)) == 0 )
}

#a function that takes a featureCounts annotation and turns it into a GRanges object, covering the entire gene body.
annotationToGRanges = function(annotation, names=NA) {
  require(GenomicRanges)
  chrLengths = mouseChrLengths()
  chr = gsub('chr','',gsub(';.*', '',annotation$Chr))
  use = chr %in% names(chrLengths)
  chr = chr[use]
  start = as.numeric(gsub(';.*', '',annotation$Start))[use]
  end =  as.numeric(unlist(lapply(lapply(strsplit(annotation$End, ';'), as.numeric), max)))[use]
  strand = gsub(';.*', '',annotation$Strand)[use]

  if ( length(names) == 1 && is.na(names) ) names = annotation$GeneID
  seqlengths = chrLengths
  ret = GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), seqlengths=seqlengths, strand=strand)
  names(ret) = names[use]
  return(ret)
}

#check if files already exist
filesInPlace = function(directory, requiredFiles) {
  f = list.files(directory, pattern='*', recursive=T)
  return( sum(!(requiredFiles %in% f)) == 0 )
}

require(rtracklayer)

#import macs2 bed files to GRanges
#The extra columns in the bed files are not supported by rtracklayer
importPeaks = function(file) {
  raw = read.table(file)
  peaks = GRanges(seqnames = gsub('MT', 'M', raw$V1), ranges = IRanges(start=raw$V2, end=raw$V3), score=raw$V5)
  return(peaks)
}
