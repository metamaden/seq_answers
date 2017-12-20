#loads the chip data from 'chipDirectory'

#import the coverage in 1000 and 100 bins.
cov1k = getAllMouseCovs(chipDirectory, binSize=1000)
cov100 = getAllMouseCovs(chipDirectory, binSize=100)
colnames(cov100) = colnames(cov1k) = c('WCE', 'H3.1', 'H3.2', 'H3K27me3.1', 'H3K27me3.2', 'H3K27me3.3')

#sum the coverage over the gene bodies and promoters.
#this gives rpkm over gene body and promoter (-2kbp to +2kbp).
cat('Summarising chip counts over genes and promoters...\n')
annGR = annotationToGRanges(fC$annotation)
names(annGR) = rownames(RNAcounts)

#get counts over promoter and gene body.
#The wide promoter count is used to exclude genes from the coverage profile analysis, that stretches out to 6 kbp from TSS.
promCounts = do.call(cbind, lapply(1:ncol(cov100), function(col) promoterCount(cov100[,col], annGR, 100)))
widePromCounts = do.call(cbind, lapply(1:ncol(cov100), function(col) promoterCount(cov100[,col], annGR, 100, range=6000)))
geneCounts = do.call(cbind, lapply(1:ncol(cov100), function(col) geneCount(cov100[,col], annGR, 100)))
colnames(promCounts) = colnames(geneCounts) = colnames(widePromCounts) = colnames(cov100)
rownames(promCounts) = rownames(geneCounts) = rownames(widePromCounts) = rownames(RNAcounts)

geneLength = abs(genesToStartX(annGR, species='mouse') - genesToEndX(annGR, species='mouse'))

#public chip data from Adli et al.
pubcov1k = getAllMouseCovs('~/data/publicchip/', binSize=1000)
pubcov100 = getAllMouseCovs('~/data/publicchip/', binSize=100)
colnames(pubcov100) = colnames(pubcov1k) = c('K27.1', 'K27.2', 'K36.1', 'K36.2', 'K4.1', 'K4.2')

pubpromCounts = do.call(cbind, lapply(1:ncol(pubcov100), function(col) promoterCount(pubcov100[,col], annGR, 100)))
pubgeneCounts = do.call(cbind, lapply(1:ncol(pubcov100), function(col) geneCount(pubcov100[,col], annGR, 100)))
colnames(pubpromCounts) = colnames(pubgeneCounts) = c('K27.1', 'K27.2', 'K36.1', 'K36.2', 'K4.1', 'K4.2')
rownames(pubpromCounts) = rownames(pubgeneCounts) = rownames(RNAcounts)
