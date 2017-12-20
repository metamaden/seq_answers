#help function for main analysis performed in input.R

#poissonian downsampling the vector x to a library size around N
#Supports parallel processes.
downSample = function(x, N, cpus=1) {
  if ( cpus > 1 ) require(parallel)
  if ( N > sum(x) ) return(x)
  keep = N/sum(x)
  dump = 1-keep
  if ( cpus == 1 )
    ret = sapply(x, function(n) {
      sum(sample(c(0,1), n, replace=T, prob=c(dump, keep)))
    })
 else
    ret = unlist(mclapply(x, function(n) {
      sum(sample(c(0,1), n, replace=T, prob=c(dump, keep)))
    }, mc.cores=cpus))
  cat('done.\n')
  return(ret)
}

#converts the colour(s) to the MCRI version of the colour, if available, otherwise returns the input.
#al sets the opaqueness.
mcri = function(col=0, al=1) {
  if ( col[1] == 0 ) {
    cat('Use: mcri(\'colour\'), returning an official MCRI colour.\nAvailable MCRI colours are:\n\ndarkblue\nblue\nlightblue\nazure\ngreen\norange\nviolet\ncyan\nred\nmagenta (aka rose).\n\nReturning default blue.\n')
    return(mcri('blue'))
  }
  if ( length(col) > 1 ) return(sapply(col, function(c) mcri(c, al)))
  if ( is.numeric(col) ) {
    if ( col == 1 ) col = 'blue'
    else if ( col == 2 ) col = 'orange'
    else if ( col == 3 ) col = 'green'
    else if ( col == 4 ) col = 'magenta'
    else if ( col == 5 ) col = 'cyan'
    else if ( col == 6 ) col = 'red'
    else if ( col == 7 ) col = 'violet'
    else if ( col == 8 ) col = 'darkblue'
    else if ( col == 9 ) col = 'azure'
    else if ( col == 10 ) col = 'lightblue'
    else col = 'black'
  }
  ret = 0
  if ( col == 'darkblue') ret = rgb(9/255, 47/255, 94/255, al)
  if ( col == 'blue') ret = rgb(0, 83/255, 161/255, al)
  if ( col == 'lightblue') ret = rgb(0, 165/255, 210/255, al)
  if ( col == 'azure') ret = rgb(0, 173/255, 239/255, al)
  if ( col == 'green') ret = rgb(141/255, 198/255, 63/255, al)
  if ( col == 'orange') ret = rgb(244/255, 121/255, 32/255, al)  
  if ( col == 'violet') ret = rgb(122/255, 82/255, 199/255, al)  
  if ( col == 'cyan') ret = rgb(0/255, 183/255, 198/255, al)  
  if ( col == 'red') ret = rgb(192/255, 80/255, 77/255, al)  
  if ( col == 'magenta' | col == 'rose') ret = rgb(236/255, 0/255, 140/255, al)
  if ( ret == 0 ) ret = do.call(rgb, as.list(c(col2rgb(col)/255, al)))
  return(ret)
}

noneg = function(x) (x+abs(x))/2

#takes a count matrix and runs a standard limma-voom De analysis on the provided groups/contrasts.
#has the option of removing genes that are highly variable within groups.
runDE = function(counts, cleanLimit=Inf, groups = c('jarid', 'suz', 'nons'), contrasts = NA) {
  require(limma)
  require(edgeR)
  messy = c()
  if ( cleanLimit < Inf ) {
    for ( group in groups ) {
      mx = counts[,grepl(group, colnames(counts), ignore.case=T)]
      flag = lapply(colnames(mx), function(col1)
        lapply(colnames(mx), function(col2)
               plotMA(mx[,col1], mx[,col2], libNorm=T, labelLimit=cleanLimit, dontPlot=T))
        )
      flag = unique(unlist(flag))
      cat('Flagging ', length(flag), ' genes for being too messy between replicates in ', group,'.\n', sep='')
      messy = union(messy, flag)
    }
    counts = counts[!(rownames(counts) %in% messy),]
  }

  
  cpms = cpm(counts)
  lowCounts = rowSums(cpms>1) < 2
  counts = counts[!lowCounts,]

  groupIs = lapply(groups, function(group) as.integer(grepl(group, colnames(counts), ignore.case=T)))

  design = do.call(cbind, groupIs)
  colnames(design) = groups
  if ( is.na(contrasts[1]) ) {
    contrastStrings = as.list(unlist(sapply(2:length(groups), function(i) sapply(1:(i-1), function(j) paste0(groups[j], '-', groups[i])))))
    contrastStrings[['levels']] = groups
    contrasts = do.call(makeContrasts, contrastStrings)
  }

  v = voom(counts, design, plot=T)
  fit = lmFit(v, design)
  fit = contrasts.fit(fit, contrasts)
  fit = eBayes(fit)
  fit$rows = which(!lowCounts)
  return(fit)
}

#takes an annotation from featureCounts and transformes it into bin indices
#for bin coverage matrices.
annotationToXRange = function( annotation, i, binsize=1000 ) {
  prevChrL = c(0, cumsum(mouseChrLengths()))
  names(prevChrL) = c(names(mouseChrLengths()), 'outside')
  chr =gsub('chr','',gsub(';.*', '', as.character(annotation$Chr[i])))
  start = as.numeric(gsub(';.*', '', annotation$Start[i]))
  end =  as.numeric(gsub('.*;', '', annotation$End[i]))
  x = round(start/binsize):round(end/binsize) + prevChrL[chr]/binsize
  return(x)
}

#Essentially just "plot", with different defaults.
#overplots to give heatmap effect for dense regions.
plotColourScatter = function(x, y, xlab='', ylab='', col=mcri('darkblue'), main='cor',
  add=F, ...) {
  cat('Correlation is ', cor(x,y), '.\n', sep='')
  if ( main == 'cor' ) main = paste('Correlation is', signif(cor(x,y), 2))
  if ( !add ) plot(x, y, cex=0.6, pch=16, xlab=xlab, ylab=ylab,
                   col=col, main=main, ...)
  else points(x, y, cex=0.6, pch=16, col=col, ...)
  points(x, y, cex=0.4, pch=16, col=mcri('blue', 0.4))
  points(x, y, cex=0.3, pch=16, col=mcri('azure', 0.1))
  points(x, y, cex=0.2, pch=16, col=mcri('green', 0.02))
}

#finds the coverage over the gene body or promoter for the provided genes.
#return format is a matrix with genes on the rows, and relative positon on the columns.
#Average coverage over all genes can conveniently be taken by colMeans of the output.
geneCoverage = function(cov, genes, binSize=100, species='mouse') {
  cat('Finding coverage over ', length(genes), ' genes.\n', sep='')
  direction = ifelse(strand(genes) == '+', 1, -1)
  starts = lapply((-50:149)/100, function(offset) ceiling((genesToStartX(genes, species=species)+direction*offset*width(genes))/binSize))
  ends = lapply((-49:150)/100, function(offset) ceiling((genesToStartX(genes, species=species)+direction*offset*width(genes))/binSize))
  countsPerBP = sapply(1:200, function(bin) sapply(1:length(starts[[bin]]), function(gene) mean(cov[starts[[bin]][gene]:ends[[bin]][gene]])/binSize))
  rpkm = 1e9*(countsPerBP/sum(cov))
  return(rpkm)
}
promCoverage = function(cov, genes, binSize=100, species='mouse') {
  cat('Finding coverage over ', length(genes), ' genes.\n', sep='')
  proms = genesToProms(genes)
  direction = ifelse(strand(proms) == '+', 1, -1)
  starts = lapply((-80:119)/40, function(offset) ceiling((genesToStartX(proms, species=species)+direction*offset*width(proms))/binSize))
  ends = lapply((-79:120)/40, function(offset) ceiling((genesToStartX(proms, species=species)+direction*offset*width(proms))/binSize))
  countsPerBP = sapply(1:200, function(bin) sapply(1:length(starts[[bin]]), function(gene) mean(cov[starts[[bin]][gene]:ends[[bin]][gene]])/binSize))
  rpkm = 1e9*(countsPerBP/sum(cov))
  return(rpkm)
}

#smoothens the vector 'x' over bins within 'range' distance
#repeats this smoothening 'repeats' times.
multiBlur = function(x, range, repeats) {
  if ( repeats == 0 ) return(x)
  for ( i in 1:repeats )
    x = blurN(x, range)

  return(x)
}
#returns a vector with the average of bins within distance M.
blurN = function(vec, M) {
  N = length(vec)
  ret = 1:N
  
  for(i in 1:M) {
    ret[i] = sum(vec[1:(i+M)])/(i+M)
  }
  for(i in (N-M):N) {
    ret[i] = sum(vec[(i-M):N])/(N-i+1+M)
  }
  for(i in (M+1):(N-M-1)) {
    ret[i] = sum(vec[(i-M):(i+M)])/(2*M+1)
  }
  return(ret)
}

#returns a GRanges object from -2kbp to +2kbp around the TSS.
genesToProms = function(genes) {
  start = ifelse(strand(genes) == '-', end(genes), start(genes))
  start(genes) = start - 2000
  end(genes) = start + 2000
  return(genes)
}
