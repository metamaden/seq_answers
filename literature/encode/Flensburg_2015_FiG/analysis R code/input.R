########################
#
# Set up input directories.
#

chipDirectory = '~/data/sarah/chipseq/'
expressionDirectory = '~/data/sarah/goodRNA/'
RscriptDirectory = '~/input/R/'

#
# Set up output directory where plots will be placed.

outputDirectory = '~/plots/input/'

#
#####################

source(paste0(RscriptDirectory, 'importChipFunctions.R'))
if ( !mouseChipFilesInPlace(chipDirectory) )
  preprocessMouseChip(chipDirectory)

#load the mouse expression data
source(paste0(RscriptDirectory, 'loadExpression.R'))
#load the mouse chip data
source(paste0(RscriptDirectory, 'loadChip.R'))

#load help function
source(paste0(RscriptDirectory, 'inputFunctions.R'))

#First downsample the 1k bin data for the control samples, to be used in the plot of distribution of counts
cat('Downsampling.\n')
libSize = min(colSums(cov1k[,c('WCE', 'H3.1', 'H3.2')]))
dsWCE = downSample(cov1k[,'WCE'], libSize, cpus=10)
dsH31 = downSample(cov1k[,'H3.1'], libSize, cpus=10)
dsH32 = downSample(cov1k[,'H3.2'], libSize, cpus=10)
ds1k = cbind('WCE'=dsWCE, 'H3.1'=dsH31, 'H3.2'=dsH32)

#take count distributions, and prepare plotting
cat('plotting count distribution.\n')
Nmax = max(ds1k)
breaks = c(-0.5, unique(round(exp((0:101)*log(Nmax)/100)) - 0.5))
mids = (breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)])/2
binW = breaks[2:length(breaks)] - breaks[1:(length(breaks)-1)]
binWlog1 = binW/(mids+1)
hWCE = hist(ds1k[,'WCE'], breaks, plot=F)
hH31 = hist(ds1k[,'H3.1'], breaks, plot=F)
hH32 = hist(ds1k[,'H3.2'], breaks, plot=F)

#create poisson distributed (flat mean) random counts
cat('Comparing to poisson distributions.\n')
nonZeroBins = rowSums(ds1k[, c('WCE', 'H3.1', 'H3.2')]) != 0
NnonZeroBins = sum(nonZeroBins)
randomCov = rep(0, nrow(ds1k))
randomCov[nonZeroBins] = rpois(NnonZeroBins, libSize/NnonZeroBins)
hrandom = hist(randomCov, breaks, plot=F)


#run limma-voom on the counts over genes, and over bins.
cat('Run differential enrichment over bins and genes with limma-voom.\n')
colnames(geneCounts)[4:6] = c('K27.1', 'K27.2', 'K27.3')
chipfitGenes = runDE(geneCounts, groups = c('WCE', 'H3', 'K27'))
colnames(geneCounts)[4:6] = c('H3K27me3.1', 'H3K27me3.2', 'H3K27me3.3')

colnames(cov1k)[4:6] = c('K27.1', 'K27.2', 'K27.3')
chipfitBins = runDE(cov1k, groups = c('WCE', 'H3', 'K27'))
chipfitBinsRandom = runDE(cbind(cov1k, randomCov), groups = c('WCE', 'H3', 'K27', 'randomCov'))
colnames(cov1k)[4:6] = c('H3K27me3.1', 'H3K27me3.2', 'H3K27me3.3')

#AM plot between bin counts
cat('Make MA plot of mouse bins.\n')
#markt he bins in Rn45s
chrL = mouseChrLengths()/1000
rn45s = cumsum(chrL)['16'] + (round(39842997/1000)+1):round(39848825/1000)
#mark bins in chromsome M
chrM = (cumsum(chrL)['Y']+1):cumsum(chrL)['M']
#mark bins not in rn45s or chr M, but still FDR below 0.05.
fdr = p.adjust(chipfitBins$p.value[,1], method='fdr')
allDE = chipfitBins$rows[which(fdr < 0.05)]

#tag the mitochondrial, Rn45s and remaining DE bins.
binstatus = rep('null', nrow(chipfitBins))
binstatus[which(chipfitBins$rows %in% allDE)] = 'significantly different bins'
binstatus[which(chipfitBins$rows %in% chrM)] = 'mithocondrial bins'
binstatus[which(chipfitBins$rows %in% rn45s)] = 'Rn45s bins'

#the random DE analysis  has different indices.
binstatusr = rep('null', nrow(chipfitBinsRandom))
binstatusr[which(chipfitBinsRandom$rows %in% allDE)] = 'significantly different bins'
binstatusr[which(chipfitBinsRandom$rows %in% chrM)] = 'mithocondrial bins'
binstatusr[which(chipfitBinsRandom$rows %in% rn45s)] = 'Rn45s bins'





##################################
#
# MACS analysis
#
##################################
cat('MACS analysis.\n')
#import peaks
WCEpeaks = importPeaks(paste0(chipDirectory, 'WCE_peaks.narrowPeak'))
H3peaks = importPeaks(paste0(chipDirectory, 'H3_peaks.narrowPeak'))
K27overWCEpeaks = importPeaks(paste0(chipDirectory, 'K27overWCE_peaks.narrowPeak'))
K27overH3peaks = importPeaks(paste0(chipDirectory, 'K27overH3_peaks.narrowPeak'))
K27peaks = importPeaks(paste0(chipDirectory, 'K27_peaks.narrowPeak'))

#find overlaps between peaks
OL = as.matrix(findOverlaps(WCEpeaks, H3peaks))
OLK27vsWCE = as.matrix(findOverlaps(K27peaks, WCEpeaks))
OLK27vsH3 = as.matrix(findOverlaps(K27peaks, H3peaks))
OLK27 = as.matrix(findOverlaps(K27overWCEpeaks, K27overH3peaks))
OLK27noC1 = as.matrix(findOverlaps(K27peaks, K27overWCEpeaks))
OLK27noC2 = as.matrix(findOverlaps(K27peaks, K27overH3peaks))

#ignore overlaps where one of the peaks overlap a higher scored peak.
#for WCE vs H3 comparison
WCEo = unique(OL[,1])
WCEbestOverlap = sapply(WCEo, function(WCEpeak) max(H3peaks$score[OL[which(OL[,1]==WCEpeak),2]]))
names(WCEbestOverlap) = as.character(WCEo)
OL = OL[H3peaks$score[OL[,2]] == WCEbestOverlap[as.character(OL[,1])],]
H3o = unique(OL[,2])
H3bestOverlap = sapply(H3o, function(H3peak) max(WCEpeaks$score[OL[which(OL[,2]==H3peak),1]]))
names(H3bestOverlap) = as.character(H3o)
OL = OL[WCEpeaks$score[OL[,1]] == H3bestOverlap[as.character(OL[,2])],]

#and for K27 over WCE vs K27 over H3 comparison
K27overWCEo = unique(OLK27[,1])
K27overWCEbestOverlap = sapply(K27overWCEo, function(K27overWCEpeak) max(K27overH3peaks$score[OLK27[which(OLK27[,1]==K27overWCEpeak),2]]))
names(K27overWCEbestOverlap) = as.character(K27overWCEo)
OLK27 = OLK27[K27overH3peaks$score[OLK27[,2]] == K27overWCEbestOverlap[as.character(OLK27[,1])],]
K27overH3o = unique(OLK27[,2])
K27overH3bestOverlap = sapply(K27overH3o, function(K27overH3peak) max(K27overWCEpeaks$score[OLK27[which(OLK27[,2]==K27overH3peak),1]]))
names(K27overH3bestOverlap) = as.character(K27overH3o)
OLK27 = OLK27[K27overWCEpeaks$score[OLK27[,1]] == K27overH3bestOverlap[as.character(OLK27[,2])],]


#tag the peaks that do not overlap the other sample
noHitWCE = which(!(1:length(WCEpeaks) %in% OL[,1]))
noHitH3 = which(!(1:length(H3peaks) %in% OL[,2]))
noHitWCEc = which(!(1:length(K27overWCEpeaks) %in% OLK27[,1]))
noHitH3c = which(!(1:length(K27overH3peaks) %in% OLK27[,2]))
noHitK271 = which(!(1:length(K27peaks) %in% OLK27noC1[,1]))
noHitK272 = which(!(1:length(K27peaks) %in% OLK27noC2[,1]))

#identify the vanishing and new K27 peaks when introducing a control
K27goneWCE = K27peaks[noHitK271]
K27goneH3 = K27peaks[noHitK272]
K27newWCE = K27overWCEpeaks[-OLK27noC1[,2]]
K27newH3 = K27overH3peaks[-OLK27noC2[,2]]

#select cut for "high score" to be top 10% amongst the H3K27me3 peaks with control.
highScoreCut = quantile(c(K27overWCEpeaks$score,K27overH3peaks$score), probs=0.75)

#identify the genes with the new and vanishing peaks in the gene body
OLnewGeneWCE = as.matrix(findOverlaps(K27newWCE, genes))
OLnewGeneH3 = as.matrix(findOverlaps(K27newH3, genes))
OLgoneGeneWCE = as.matrix(findOverlaps(K27goneWCE, genes))
OLgoneGeneH3 = as.matrix(findOverlaps(K27goneH3, genes))
OLnewGeneWCEhigh = as.matrix(findOverlaps(K27newWCE[K27newWCE$score > highScoreCut], genes))
OLnewGeneH3high = as.matrix(findOverlaps(K27newH3[K27newH3$score > highScoreCut], genes))
OLgoneGeneWCEhigh = as.matrix(findOverlaps(K27goneWCE[K27goneWCE$score > highScoreCut], genes))
OLgoneGeneH3high = as.matrix(findOverlaps(K27goneH3[K27goneH3$score > highScoreCut], genes))

#identify the genes with the new and vanishing peaks in the promoters
proms = genesToProms(genes)
OLnewPromWCE = as.matrix(findOverlaps(K27newWCE, proms))
OLnewPromH3 = as.matrix(findOverlaps(K27newH3, proms))
OLgonePromWCE = as.matrix(findOverlaps(K27goneWCE, proms))
OLgonePromH3 = as.matrix(findOverlaps(K27goneH3, proms))
OLnewPromWCEhigh = as.matrix(findOverlaps(K27newWCE[K27newWCE$score > highScoreCut], proms))
OLnewPromH3high = as.matrix(findOverlaps(K27newH3[K27newH3$score > highScoreCut], proms))
OLgonePromWCEhigh = as.matrix(findOverlaps(K27goneWCE[K27goneWCE$score > highScoreCut], proms))
OLgonePromH3high = as.matrix(findOverlaps(K27goneH3[K27goneH3$score > highScoreCut], proms))

#print some statistics on what happens with the peak calling with and without control
{
  cat('Number of new peaks in K27, when using WCE as control:', length(K27newWCE))
  cat('\nNumber of new peaks overlapping gene or promoter:', length(unique(c(OLnewGeneWCE[,1], OLnewPromWCE[,1]))),
      ' (',length(unique(c(OLnewGeneWCE[,1], OLnewPromWCE[,1])))/length(K27newWCE), ')', sep='')
  cat('\nNumber of new peaks overlapping promoter:', length(unique(OLnewPromWCE[,1])),
            ' (',length(unique(OLnewPromWCE[,1]))/length(K27newWCE), ')', sep='')
  cat('\nNumber of new peaks in K27 with high score, when using WCE as control:', sum(K27newWCE$score > highScoreCut))
  cat('\nNumber of new high peaks overlapping gene or promoter: ', length(unique(c(OLnewGeneWCEhigh[,1], OLnewPromWCEhigh[,1]))),
            ' (',length(unique(c(OLnewGeneWCEhigh[,1], OLnewPromWCEhigh[,1])))/sum(K27newWCE$score > highScoreCut), ')', sep='')
  cat('\nNumber of new high peaks overlapping promoter: ', length(unique(OLnewPromWCEhigh[,1])),
                  ' (',length(unique(OLnewPromWCEhigh[,1]))/sum(K27newWCE$score > highScoreCut), ')', sep='')
  cat('\nNumber of promoters getting a new high score peak: ', length(unique(OLnewPromWCEhigh[,2])), sep='')
  cat('\n\n')
  cat('Number of new peaks in K27, when using H3 as control:', length(K27newH3))
  cat('\nNumber of new peaks overlapping gene or promoter:', length(unique(c(OLnewGeneH3[,1], OLnewPromH3[,1]))),
      ' (',length(unique(c(OLnewGeneH3[,1], OLnewPromH3[,1])))/length(K27newH3), ')', sep='')
  cat('\nNumber of new peaks overlapping promoter:', length(unique(OLnewPromH3[,1])),
            ' (',length(unique(OLnewPromH3[,1]))/length(K27newH3), ')', sep='')
  cat('\nNumber of new peaks in K27 with high score, when using H3 as control:', sum(K27newH3$score > highScoreCut))
  cat('\nNumber of new high peaks overlapping gene or promoter: ', length(unique(c(OLnewGeneH3high[,1], OLnewPromH3high[,1]))),
            ' (',length(unique(c(OLnewGeneH3high[,1], OLnewPromH3high[,1])))/sum(K27newH3$score > highScoreCut), ')', sep='')
  cat('\nNumber of new high peaks overlapping promoter: ', length(unique(OLnewPromH3high[,1])),
                  ' (',length(unique(OLnewPromH3high[,1]))/sum(K27newH3$score > highScoreCut), ')', sep='')
  cat('\nNumber of promoters getting a new high score peak: ', length(unique(OLnewPromH3high[,2])), sep='')
  cat('\n\n')
  cat('Number of gone peaks in K27, when using WCE as control:', length(K27goneWCE))
  cat('\nNumber of gone peaks overlapping gene or promoter:', length(unique(c(OLgoneGeneWCE[,1], OLgonePromWCE[,1]))),
      ' (',length(unique(c(OLgoneGeneWCE[,1], OLgonePromWCE[,1])))/length(K27goneWCE), ')', sep='')
  cat('\nNumber of gone peaks overlapping promoter:', length(unique(OLgonePromWCE[,1])),
            ' (',length(unique(OLgonePromWCE[,1]))/length(K27goneWCE), ')', sep='')
  cat('\nNumber of gone peaks in K27 with high score, when using WCE as control:', sum(K27goneWCE$score > highScoreCut))
  cat('\nNumber of gone high peaks overlapping gene or promoter: ', length(unique(c(OLgoneGeneWCEhigh[,1], OLgonePromWCEhigh[,1]))),
            ' (',length(unique(c(OLgoneGeneWCEhigh[,1], OLgonePromWCEhigh[,1])))/sum(K27goneWCE$score > highScoreCut), ')', sep='')
  cat('\nNumber of gone high peaks overlapping promoter: ', length(unique(OLgonePromWCEhigh[,1])),
                  ' (',length(unique(OLgonePromWCEhigh[,1]))/sum(K27goneWCE$score > highScoreCut), ')', sep='')
  cat('\nNumber of promoters getting a gone high score peak: ', length(unique(OLgonePromWCEhigh[,2])), sep='')
  cat('\n\n')
  cat('Number of gone peaks in K27, when using H3 as control:', length(K27goneH3))
  cat('\nNumber of gone peaks overlapping gene or promoter:', length(unique(c(OLgoneGeneH3[,1], OLgonePromH3[,1]))),
      ' (',length(unique(c(OLgoneGeneH3[,1], OLgonePromH3[,1])))/length(K27goneH3), ')', sep='')
  cat('\nNumber of gone peaks overlapping promoter:', length(unique(OLgonePromH3[,1])),
            ' (',length(unique(OLgonePromH3[,1]))/length(K27goneH3), ')', sep='')
  cat('\nNumber of gone peaks in K27 with high score, when using H3 as control:', sum(K27goneH3$score > highScoreCut))
  cat('\nNumber of gone high peaks overlapping gene or promoter: ', length(unique(c(OLgoneGeneH3high[,1], OLgonePromH3high[,1]))),
            ' (',length(unique(c(OLgoneGeneH3high[,1], OLgonePromH3high[,1])))/sum(K27goneH3$score > highScoreCut), ')', sep='')
  cat('\nNumber of gone high peaks overlapping promoter: ', length(unique(OLgonePromH3high[,1])),
                  ' (',length(unique(OLgonePromH3high[,1]))/sum(K27goneH3$score > highScoreCut), ')', sep='')
  cat('\nNumber of promoters getting a gone high score peak: ', length(unique(OLgonePromH3high[,2])), sep='')
  cat('\n')
}

#extract the genes that get new peaks in the promoter/gene body with control
breaks = ((-20):40)/10
hExp = hist(log10(rpkm), breaks=breaks, plot=F)
newExpWCE = rpkm[unique(OLnewPromWCEhigh[,2])]
hExpNewWCE = hist(log10(newExpWCE), breaks=hExp$breaks)
newExpH3 = rpkm[unique(OLnewPromH3high[,2])]
hExpNewH3 = hist(log10(newExpH3), breaks=hExp$breaks)

newExpWCEg = rpkm[unique(OLnewGeneWCEhigh[,2])]
hExpNewWCEg = hist(log10(newExpWCEg), breaks=hExp$breaks)
newExpH3g = rpkm[unique(OLnewGeneH3high[,2])]
hExpNewH3g = hist(log10(newExpH3g), breaks=hExp$breaks)



##################################
#
# prepare MA plots over genes
#
##################################
cat('Analyse Enrichment shape over genes and promoters.\n')
rn45sGene = which(rownames(geneCounts) == 'Rn45s')
#mark bins in chromsome M
chrMGene = which(rownames(geneCounts) %in% c("Rnr1", "Rnr2", "ND1", "ND2", "COX1", "ND4", "ND5", "ND6", "CYTB"))
#mark bins not in rn45s or chr M, but still FDR below 0.05.
derGene = chipfitGenes$rows[which(p.adjust(chipfitGenes$p.value[,1], method='fdr') < 0.05)]

#transform to a status vector that can be fed into edgeR's plotMA().
genestatus = rep('null', nrow(chipfitGenes))
genestatus[which(chipfitGenes$rows %in% derGene)] = 'significantly different genes'
genestatus[which(chipfitGenes$rows %in% chrMGene)] = 'mithocondrial genes'
genestatus[which(chipfitGenes$rows %in% rn45sGene)] = 'Rn45s'


##################################
#
# expression correlation
#
##################################
#First get the RPKMs. The added cpm and rpkm are to avoid noise from the logarithm
#of very low number (such as genes without reads) in the Pearson correlations.
cat('Correlations with expression.\n')
rpkm = (rowSums(RNAcounts)/sum(RNAcounts)*1e6 + 1)/(exonLength/1000)
K27 = 0.5 + rowSums(geneCounts[,4:6])*1e6/sum(geneCounts[,4:6])/(geneLength/1000)
pubK27 = 0.5 + rowSums(pubgeneCounts[,1:2])*1e6/sum(pubgeneCounts[,4:6])/(geneLength/1000)
pubK36 = 0.5 + rowSums(pubgeneCounts[,3:4])*1e6/sum(pubgeneCounts[,3:4])/(geneLength/1000)
WCE = 0.5 + geneCounts[,1]*1e6/sum(geneCounts[,1])/(geneLength/1000)
H3 = 0.5 + rowSums(geneCounts[,2:3])*1e6/sum(geneCounts[,2:3])/(geneLength/1000)

#RPKM over promoters
K27p = 0.5 + rowSums(promCounts[,4:6])*1e6/sum(promCounts[,4:6])/4
pubK27p = 0.5 + rowSums(pubpromCounts[,1:2])*1e6/sum(pubpromCounts[,4:6])/4
pubK36p = 0.5 + rowSums(pubpromCounts[,3:4])*1e6/sum(pubpromCounts[,3:4])/4
WCEp = 0.5 + promCounts[,1]*1e6/sum(promCounts[,1])/4
H3p = 0.5 + rowSums(promCounts[,2:3])*1e6/sum(promCounts[,2:3])/4

#scatter plot log(RPKMs) against each other.
#This will output the correlations (in the table in the paper) to the terminal as well.
plotColourScatter(log10(WCE), log10(rpkm), xlab = 'log(K27)', 'log(expression)')
plotColourScatter(log10(H3), log10(rpkm), xlab = 'log(K27)', 'log(expression)')
plotColourScatter(log10(WCEp), log10(rpkm), xlab = 'log(K27)', 'log(expression)')
plotColourScatter(log10(H3p), log10(rpkm), xlab = 'log(K27)', 'log(expression)')

plotColourScatter(log10(K27), log10(rpkm), xlab = 'log(K27)', 'log(expression)')
plotColourScatter(log10(K27/WCE), log10(rpkm), xlab = 'log(K27/WCE)', 'log(expression)')
plotColourScatter(log10(K27/H3), log10(rpkm), xlab = 'log(K27/H3)', 'log(expression)')

plotColourScatter(log10(pubK27), log10(rpkm), xlab = 'log(adli K27)', 'log(expression)')
plotColourScatter(log10(pubK27/WCE), log10(rpkm), xlab = 'log(adli K27/WCE)', 'log(expression)')
plotColourScatter(log10(pubK27/H3), log10(rpkm), xlab = 'log(adli K27/H3)', 'log(expression)')

plotColourScatter(log10(pubK36), log10(rpkm), xlab = 'log(adli K36)', 'log(expression)')
plotColourScatter(log10(pubK36/WCE), log10(rpkm), xlab = 'log(adli K36/WCE)', 'log(expression)')
plotColourScatter(log10(pubK36/H3), log10(rpkm), xlab = 'log(adli K36/H3)', 'log(expression)')

plotColourScatter(log10(K27p), log10(rpkm), xlab = 'log(K27p)', 'log(expression)')
plotColourScatter(log10(K27p/WCEp), log10(rpkm), xlab = 'log(K27/WCE)', 'log(expression)')
plotColourScatter(log10(K27p/H3p), log10(rpkm), xlab = 'log(K27/H3)', 'log(expression)')

plotColourScatter(log10(pubK27p), log10(rpkm), xlab = 'log(adli K27)', 'log(expression)')
plotColourScatter(log10(pubK27p/WCEp), log10(rpkm), xlab = 'log(adli K27/WCE)', 'log(expression)')
plotColourScatter(log10(pubK27p/H3p), log10(rpkm), xlab = 'log(adli K27/H3)', 'log(expression)')

plotColourScatter(log10(pubK36p), log10(rpkm), xlab = 'log(adli K36)', 'log(expression)')
plotColourScatter(log10(pubK36p/WCEp), log10(rpkm), xlab = 'log(adli K36/WCE)', 'log(expression)')
plotColourScatter(log10(pubK36p/H3p), log10(rpkm), xlab = 'log(adli K36/H3)', 'log(expression)')


#remove very high rpkms, to be able to do non-logged averages without being dominated by a few genes.
#these are all ribosomal and mitochondrial in our data
WCEuse = 1e9*geneCounts[,1]/width(annGR)/sum(geneCounts[,1]) < 100
H3use = 1e9*rowSums(geneCounts[,2:3])/width(annGR)/sum(geneCounts[,2:3]) < 100
K27use = 1e9*rowSums(geneCounts[,4:6])/width(annGR)/sum(geneCounts[,4:6]) < 100

WCEpuse = 1e9*widePromCounts[,1]/4000/sum(widePromCounts[,1]) < 100
H3puse = 1e9*rowSums(widePromCounts[,2:3])/4000/sum(widePromCounts[,2:3]) < 100
K27puse = 1e9*rowSums(widePromCounts[,4:6])/4000/sum(widePromCounts[,4:6]) < 100

#divide genes in four quantiles based on expression.
qs = quantile(rpkm, probs=c(0.25, 0.5, 0.75))
q1 = rpkm < qs[1]
q2 = rpkm > qs[1] & rpkm < qs[2]
q3 = rpkm > qs[2] & rpkm < qs[3]
q4 = rpkm > qs[3]

#calculate average gene coverage for each quantile for each ChIP-sample.
WCElist = list(q1 & WCEuse, q2 & WCEuse, q3 & WCEuse, q4 & WCEuse)
H3list = list(q1 & H3use, q2 & H3use, q3 & H3use, q4 & H3use)
K27list = list(q1 & K27use, q2 & K27use, q3 & K27use, q4 & K27use)
WCEgeneCovs = sapply(WCElist, function(subset) colMeans(geneCoverage(cov100[,1], annGR[subset], 100)))
H3geneCovs = sapply(H3list, function(subset) colMeans(geneCoverage(rowSums(cov100[,2:3]), annGR[subset], 100)))
K27geneCovs = sapply(K27list, function(subset) colMeans(geneCoverage(rowSums(cov100[,4:6]), annGR[subset], 100)))

#calculate average coverage around TSS for each quantile for each ChIP-sample.
WCEplist = list(q1 & WCEpuse, q2 & WCEpuse, q3 & WCEpuse, q4 & WCEpuse)
H3lpist = list(q1 & H3puse, q2 & H3puse, q3 & H3puse, q4 & H3puse)
K27plist = list(q1 & K27puse, q2 & K27puse, q3 & K27puse, q4 & K27puse)
WCEpromCovs = sapply(WCEplist, function(subset) colMeans(promCoverage(cov100[,1], annGR[subset], 100)))
H3promCovs = sapply(H3lpist, function(subset) colMeans(promCoverage(rowSums(cov100[,2:3]), annGR[subset], 100)))
K27promCovs = sapply(K27list, function(subset) colMeans(promCoverage(rowSums(cov100[,4:6]), annGR[subset], 100)))


######################################################################################
#
# Export counts to a .csv file.
#
######################################################################################

export = data.frame('entrezID'=names(genes), 'symbolID'=rownames(geneCounts),
  'WCE.gene'=geneCounts[,'WCE'],
  'H3.rep1.gene'=geneCounts[,'H3.1'], 'H3.rep2.gene'=geneCounts[,'H3.2'],
  'H3K27me3.rep1.gene'=geneCounts[,'K27.1'], 'H3K27me3.rep2.gene'=geneCounts[,'K27.2'], 'H3K27me3.rep3.gene'=geneCounts[,'K27.3'],
  'WCE.promoter'=promCounts[,'WCE'],
  'H3.rep1.promoter'=promCounts[,'H3.1'], 'H3.rep2.promoter'=promCounts[,'H3.2'],
  'H3K27me3.rep1.promoter'=promCounts[,'K27.1'], 'H3K27me3.rep2.promoter'=promCounts[,'K27.2'], 'H3K27me3.rep3.promoter'=promCounts[,'K27.3'],
  'RNA.rep1.exons'=RNAcounts[,1], 'RNA.rep2.exons'=RNAcounts[,2], 'RNA.rep3.exons'=RNAcounts[,3])

write.csv(export, paste0(outputDirectory, 'counts.csv'))


######################################################################################
#
# Submitted figures
#
######################################################################################



pdf(paste0(outputDirectory, 'fig1.pdf'), width=7, height=7)
x = (250:3300)/100
s = smooth.spline(hWCE$mids, hWCE$counts/binWlog1/1e6, spar=0.1)
pred = predict(s$fit, x)
plot(pred$x, pred$y, xlim=c(3,27), ylim=c(0, 3.5),
     type='l', ylab='density(millions of bins)', xlab='counts per bin', log='x',
     lwd=3, col=mcri('blue'))
s = smooth.spline(hH31$mids, hH31$counts/binWlog1/1e6, spar=0.1)
pred = predict(s$fit, x)
lines(pred$x, pred$y, col=mcri('red'), lwd=3)
s = smooth.spline(hH32$mids, hH32$counts/binWlog1/1e6, spar=0.1)
pred = predict(s$fit, x)
lines(pred$x, pred$y, col=mcri('red'), lwd=3)
s = smooth.spline(hrandom$mids, hrandom$counts/binWlog1/1e6, spar=0.1)
pred = predict(s$fit, x)
lines(pred$x, pred$y, col=mcri('darkblue'), lwd=2, lty=2)
legend('topleft', c('WCE', 'H3', 'poisson'), col=c(mcri('blue'), mcri('red'), mcri('darkblue')), lwd=c(5, 5, 2), lty=c(1,1,2))
dev.off()

pdf(paste0(outputDirectory, 'fig2.pdf'), width=7, height=10)
layout(matrix(c(1,2), nrow=2))
plotMA(chipfitBins[,1], status=binstatus, values = c('mithocondrial bins', 'Rn45s bins', 'significantly different bins'), col=mcri(c('red', 'blue', 'green')), main='(A) bin read counts in WCE vs H3', xlab=expression('log'[2]*'(average read count)'), ylab=expression('log'[2]*'(WCE/H3)'))

plotMA(chipfitGenes[,1], status=genestatus, values = c('mithocondrial genes', 'Rn45s', 'significantly different genes'), col=mcri(c('red', 'blue', 'green')), main='(B) gene read counts in WCE vs H3', xlab=expression('log'[2]*'(average read count)'), ylab=expression('log'[2]*'(WCE/H3)'))
dev.off()

jpeg(paste0(outputDirectory, 'fig3.jpg'), width=7, height=7, units = 'in', res=900)
layout(matrix(c(1,3,2,4), nrow=2))
boxplot(log10(WCEpeaks$score[-noHitWCE]), log10(WCEpeaks$score[noHitWCE]),
        log10(H3peaks$score[-noHitH3]), log10(H3peaks$score[noHitH3]),
        col=c(mcri('blue', 1), mcri('blue', 0.3), mcri('red'), mcri('red', 0.3)),
        names=c('WCE o', 'WCE n', 'H3 o', 'H3 n'),
        ylab='peak score', main='(A) peak scores: controls', pch=16)

plotColourScatter(log10(WCEpeaks$score[OL[,1]]), log10(H3peaks$score[OL[,2]]), xlab=expression('WCE log'[10]*'(score)'), ylab=expression('H3 log'[10]*'(score)'), main='(B) scores of overlapping peaks')

boxplot(log10(K27overWCEpeaks$score[-noHitWCEc]), log10(K27overWCEpeaks$score[noHitWCEc]),
        log10(K27overH3peaks$score[-noHitH3c]), log10(K27overH3peaks$score[noHitH3c]),
        col=c(mcri('blue', 1), mcri('blue', 0.3), mcri('red'), mcri('red', 0.3)),
        names=c('WCE o', 'WCE n', 'H3 o', 'H3 n'),
        ylab='peak score', main='(C) peak scores: H3K27me3 with control', pch=16)

plotColourScatter(log10(K27overWCEpeaks$score[OLK27[,1]]), log10(K27overH3peaks$score[OLK27[,2]]), xlab=expression('H3K27me3, WCE control log'[10]*'(score)'), ylab=expression('H3K27me3, H3 control log'[10]*'(score)'), main='(D) scores of overlapping peaks')
dev.off()

pdf(paste0(outputDirectory, 'fig4.pdf'), width=7, height=10)
layout(matrix(c(1,3,5,7,2,4,6,8), nrow=4))
par(mar = c(3, 4, 2, 2))
ymax = 1.3
lineW = 0.7
blur = 2
xGene = (-50:149 + 0.5)/100
plot(xGene, WCEgeneCovs[,1], type='n', lwd=lineW, xlim=c(-0.3,1.3), ylim=c(0, ymax),
     col=mcri('blue'), xlab='', ylab='average rpkm', main='(A) Gene body H3K27me3')
segments(c(0, 1), c(0, 0), c(0, 1), c(ymax*2, ymax*2), lwd=3, col='grey')
lines(xGene, multiBlur(K27geneCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('green'))
lines(xGene, multiBlur(K27geneCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('green'))
lines(xGene, multiBlur(K27geneCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('green'))
lines(xGene, multiBlur(K27geneCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('green'))

ymax = 1.6
blur = 1
xProm = (-100:99 + 0.5)*100
plot(xProm, WCEpromCovs[,1], type='n', lwd=lineW, xlim=c(-6000, 6000), ylim=c(0, ymax),
     col=mcri('blue'), xlab='bp from TSS', ylab='average rpkm', main='(B) Promoter H3K27me3')
segments(0, 0, 0, ymax*2, lwd=3, col='grey')
lines(xProm, multiBlur(K27promCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('green'))
lines(xProm, multiBlur(K27promCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('green'))
lines(xProm, multiBlur(K27promCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('green'))
lines(xProm, multiBlur(K27promCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('green'))

ymax = 1.3
lineW = 0.7
blur = 2
xGene = (-50:149 + 0.5)/100
plot(xGene, WCEgeneCovs[,1], type='n', lwd=lineW, xlim=c(-0.3,1.3), ylim=c(0, ymax),
     col=mcri('blue'), xlab='position in gene', ylab='average rpkm', main='(C) Gene body WCE')
segments(c(0, 1), c(0, 0), c(0, 1), c(ymax*2, ymax*2), lwd=3, col='grey')
lines(xGene, multiBlur(WCEgeneCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('blue'))
lines(xGene, multiBlur(WCEgeneCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('blue'))
lines(xGene, multiBlur(WCEgeneCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('blue'))
lines(xGene, multiBlur(WCEgeneCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('blue'))

ymax = 1.6
blur = 1
xProm = (-100:99 + 0.5)*100
plot(xProm, WCEpromCovs[,1], type='n', lwd=lineW, xlim=c(-6000, 6000), ylim=c(0, ymax),
     col=mcri('blue'), xlab='bp from TSS', ylab='average rpkm', main='(D) Promoter WCE')
segments(0, 0, 0, ymax*2, lwd=3, col='grey')
lines(xProm, multiBlur(WCEpromCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('blue'))
lines(xProm, multiBlur(WCEpromCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('blue'))
lines(xProm, multiBlur(WCEpromCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('blue'))
lines(xProm, multiBlur(WCEpromCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('blue'))

ymax = 1.3
lineW = 0.7
blur = 2
xGene = (-50:149 + 0.5)/100
plot(xGene, WCEgeneCovs[,1], type='n', lwd=lineW, xlim=c(-0.3,1.3), ylim=c(0, ymax),
     col=mcri('blue'), xlab='position in gene', ylab='average rpkm', main='(E) Gene body H3')
segments(c(0, 1), c(0, 0), c(0, 1), c(ymax*2, ymax*2), lwd=3, col='grey')
lines(xGene, multiBlur(H3geneCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('red'))
lines(xGene, multiBlur(H3geneCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('red'))
lines(xGene, multiBlur(H3geneCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('red'))
lines(xGene, multiBlur(H3geneCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('red'))

ymax = 1.6
blur = 1
xProm = (-100:99 + 0.5)*100
plot(xProm, WCEpromCovs[,1], type='n', lwd=lineW, xlim=c(-6000, 6000), ylim=c(0, ymax),
     col=mcri('blue'), xlab='bp from TSS', ylab='average rpkm', main='(F) Promoter H3')
segments(0, 0, 0, ymax*2, lwd=3, col='grey')
lines(xProm, multiBlur(H3promCovs[,1], 1, blur), type='l', lwd=1*lineW, col=mcri('red'))
lines(xProm, multiBlur(H3promCovs[,2], 1, blur), type='l', lwd=2*lineW, col=mcri('red'))
lines(xProm, multiBlur(H3promCovs[,3], 1, blur), type='l', lwd=3*lineW, col=mcri('red'))
lines(xProm, multiBlur(H3promCovs[,4], 1, blur), type='l', lwd=4*lineW, col=mcri('red'))

par(mar = c(5, 4, 2, 2))
ymax = 4
lineW = 0.5
blur = 2
plot(xGene, WCEgeneCovs[,1], type='n', lwd=lineW, xlim=c(-0.3,1.3), ylim=c(0, ymax),
     col=mcri('blue'), xlab='position in gene', ylab='relative enrichment', main='(G) Gene body H3K27me3/control')
segments(c(0, 1), c(0, 0), c(0, 1), c(ymax*2, ymax*2), lwd=3, col='grey')
lines(xGene, multiBlur(K27geneCovs[,1], 1, blur)/multiBlur(WCEgeneCovs[,1], 1, blur),
      type='l', lwd=1*lineW, col=mcri('blue'))
lines(xGene, multiBlur(K27geneCovs[,2], 1, blur)/multiBlur(WCEgeneCovs[,2], 1, blur),
      type='l', lwd=2*lineW, col=mcri('blue'))
lines(xGene, multiBlur(K27geneCovs[,3], 1, blur)/multiBlur(WCEgeneCovs[,3], 1, blur),
      type='l', lwd=3*lineW, col=mcri('blue'))
lines(xGene, multiBlur(K27geneCovs[,4], 1, blur)/multiBlur(WCEgeneCovs[,4], 1, blur),
      type='l', lwd=4*lineW, col=mcri('blue'))
lines(xGene, multiBlur(K27geneCovs[,1], 1, blur)/multiBlur(H3geneCovs[,1], 1, blur),
      type='l', lwd=1*lineW, col=mcri('red'))
lines(xGene, multiBlur(K27geneCovs[,2], 1, blur)/multiBlur(H3geneCovs[,2], 1, blur),
      type='l', lwd=2*lineW, col=mcri('red'))
lines(xGene, multiBlur(K27geneCovs[,3], 1, blur)/multiBlur(H3geneCovs[,3], 1, blur),
      type='l', lwd=3*lineW, col=mcri('red'))
lines(xGene, multiBlur(K27geneCovs[,4], 1, blur)/multiBlur(H3geneCovs[,4], 1, blur),type='l', lwd=4*lineW, col=mcri('red'))
legend('topright', c('K27/WCE', 'K27/H3'), lwd=3, col=c(mcri('blue'), mcri('red')), bg='white')

plot(xProm, WCEpromCovs[,1], type='n', lwd=lineW, xlim=c(-6000, 6000), ylim=c(0, ymax),
     col=mcri('blue'), xlab='bp from TSS', ylab='relative enrichment', main='(H) Promoter H3K27me3/control')
segments(0, 0, 0, ymax*2, lwd=3, col='grey')
lines(xProm, multiBlur(K27promCovs[,1], 1, blur)/multiBlur(WCEpromCovs[,1], 1, blur),
      type='l', lwd=1*lineW, col=mcri('blue'))
lines(xProm, multiBlur(K27promCovs[,2], 1, blur)/multiBlur(WCEpromCovs[,2], 1, blur),
      type='l', lwd=2*lineW, col=mcri('blue'))
lines(xProm, multiBlur(K27promCovs[,3], 1, blur)/multiBlur(WCEpromCovs[,3], 1, blur),
      type='l', lwd=3*lineW, col=mcri('blue'))
lines(xProm, multiBlur(K27promCovs[,4], 1, blur)/multiBlur(WCEpromCovs[,4], 1, blur),
      type='l', lwd=4*lineW, col=mcri('blue'))
lines(xProm, multiBlur(K27promCovs[,1], 1, blur)/multiBlur(H3promCovs[,1], 1, blur),
      type='l', lwd=1*lineW, col=mcri('red'))
lines(xProm, multiBlur(K27promCovs[,2], 1, blur)/multiBlur(H3promCovs[,2], 1, blur),
      type='l', lwd=2*lineW, col=mcri('red'))
lines(xProm, multiBlur(K27promCovs[,3], 1, blur)/multiBlur(H3promCovs[,3], 1, blur),
      type='l', lwd=3*lineW, col=mcri('red'))
lines(xProm, multiBlur(K27promCovs[,4], 1, blur)/multiBlur(H3promCovs[,4], 1, blur),type='l', lwd=4*lineW, col=mcri('red'))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

pdf(paste0(outputDirectory, 'supfig2.pdf'), width=7, height=7)
hist(log10(rpkm), breaks=hExp$breaks, col=rgb(0.9, 0.9, 0.9), xlim=c(-1, 2.5), main = 'New high score peaks in H3K27me3 from control', xlab = expression('expression log'[10]*'(RPKM)'), ylab='gene frequency')
lines(hExpNewWCEg$mids, multiBlur(hExpNewWCEg$counts*10, 1,3), lwd=5, col=mcri('blue'))
lines(hExpNewH3g$mids, multiBlur(hExpNewH3g$counts*10, 1,3), lwd=3, col=mcri('red'))
legend('topright', c('all genes', 'new WCE genes (x10)', 'new H3 genes (x10)'), lwd=c(10, 5, 5), col=c(rgb(0.7, 0.7, 0.7), mcri('blue'), mcri('red')))
dev.off()

jpeg(paste0(outputDirectory, 'supfig2.jpg'), width=7, height=7, units = 'in', res=900)
hist(log10(rpkm), breaks=hExp$breaks, col=rgb(0.9, 0.9, 0.9), xlim=c(-1, 2.5), main = 'New high score peaks in H3K27me3 from control', xlab = expression('expression log'[10]*'(RPKM)'), ylab='gene frequency')
lines(hExpNewWCEg$mids, multiBlur(hExpNewWCEg$counts*10, 1,3), lwd=5, col=mcri('blue'))
lines(hExpNewH3g$mids, multiBlur(hExpNewH3g$counts*10, 1,3), lwd=3, col=mcri('red'))
legend('topright', c('all genes', 'new WCE genes (x10)', 'new H3 genes (x10)'), lwd=c(10, 5, 5), col=c(rgb(0.7, 0.7, 0.7), mcri('blue'), mcri('red')))
dev.off()


jpeg(paste0(outputDirectory, 'supfig1.jpg'), width=7, height=15, units = 'in', res=900)
layout(matrix(c(1,2, 3), nrow=3))
limma::plotMA(chipfitBinsRandom[,4], status=binstatusr, values = c('mithocondrial bins', 'Rn45s bins', 'significantly different bins'), col=mcri(c('red', 'blue', 'green')), main='(A) bin read counts in WCE vs random', xlab='Amean', ylab='logFC', ylim=c(-5, 11))
limma::plotMA(chipfitBinsRandom[,5], status=binstatusr, values = c('mithocondrial bins', 'Rn45s bins', 'significantly different bins'), col=mcri(c('red', 'blue', 'green')), main='(B) bin read counts in H3 vs random', xlab='Amean', ylab='logFC', ylim=c(-5, 11))
limma::plotMA(chipfitBins[,1], status=binstatus, values = c('mithocondrial bins', 'Rn45s bins', 'significantly different bins'), col=mcri(c('red', 'blue', 'green')), main='(C) bin read counts in WCE vs H3', xlab='Amean', ylab='logFC', ylim=c(-5, 11))
dev.off()



