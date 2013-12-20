rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

registerDoMC(5)

sum.marks = function(counts){
  new.counts = counts
  new.counts[, grep('left', colnames(new.counts))] = rowMedians(counts[, grep('left', colnames(new.counts))])
  new.counts[, grep('right', colnames(new.counts))] = rowMedians(counts[, grep('right', colnames(new.counts))])
  return(new.counts)
}
read.counts = function(region.file, indir, pref, marks, indivs, q){
  counts = NULL
  good.rows = NULL
  for(m in marks){
    for(p in c('left', 'right')){
      signal.files = file.path(indir, paste(pref, '_', p, '_AT_SNYDER_HG19_', indivs, '_', m, '.txt', sep = ''))
      counts.tmp = asinh(load.avg.sig.data(region.file, signal.files, indivs)$signal)
      colnames(counts.tmp) = paste(m, p, colnames(counts.tmp), sep = '_')
      counts = cbind(counts, counts.tmp)
      good.rows.tmp = apply(counts, 1, function(x) !any(is.na(x))) & rowMedians(counts) > quantile(rowMedians(counts), q, na.rm = T)
      
      if(is.null(good.rows)){good.rows = good.rows.tmp
      }else{good.rows = good.rows | good.rows.tmp}
    }
  }
  regions = read.bed(region.file)
  counts = counts[good.rows, ]
  regions = regions[good.rows, ]
  return(list(regions = regions, counts = counts))
}
read.counts.2 = function(region.file, indir, pref, marks, indivs, q){
  counts = NULL
  for(m in marks){
    signal.files = file.path(indir, paste(gsub('.bed', '', basename(region.file)), '_AT_SNYDER_HG19_', indivs, '_', m, '.txt', sep = ''))
    counts.dat = load.avg.sig.data(region.file, signal.files, indivs)
    counts.dip = counts.dat$signal
    regions = counts.dat$regions
    for(p in c('left', 'right')){
      signal.files = file.path(indir, paste(pref, '_', p, '_AT_SNYDER_HG19_', indivs, '_', m, '.txt', sep = ''))
      counts.tmp = asinh(load.avg.sig.data(region.file, signal.files, indivs)$signal)
      counts.tmp[rowMeans(counts.tmp) < quantile(rowMeans(counts.tmp), q, na.rm = T), ] = NaN
      counts.tmp = asinh(counts.dip) - counts.tmp
      colnames(counts.tmp) = paste(m, p, colnames(counts.tmp), sep = '_')
      counts = cbind(counts, counts.tmp)
    }
  }
  return(list(regions = regions, counts = counts))
}

outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/'
outpref = ''

# Load SNP information and genotypes
load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load('../../rawdata/variants/all/snps/allNonSan/genot_pca.RData')
snp.pos = snp.pos[pca.rows, ] # This selects rows that were used for PCA (so had some min variance in genotype)
indivs = colnames(genot)
nindivs = ncol(genot)

# Load regions 
# count.dat = read.counts('../../rawdata/dhs/alan/pritchard_dhs.bed',
#             '../../rawdata/dhs/alan/combrep/extractSignal/fc/avgSig/textFiles/', 
#             'pritchard_dhs_200bp', c('H3K27AC'), indivs, 0.4)
# regions = count.dat$regions
# counts = normalize.quantiles(count.dat$counts)
# colnames(counts) = colnames(count.dat$counts)

marks = c('H3K27AC', 'H3K4ME3')
count.dat = read.counts.2('../../rawdata/dhs/alan/pritchard_dhs.bed',
                        '../../rawdata/dhs/alan/combrep/extractSignal/fc/avgSig/textFiles/', 
                        'pritchard_dhs_200bp', marks, indivs, 0.4)
regions = count.dat$regions
counts = count.dat$counts

# Remove regions overlapping indels
indels = read.bed('../../rawdata/variants/all/masks/all.blacklist.bed')
good.rows = is.na(findOverlaps(regions.to.ranges(regions), regions.to.ranges(indels), select = 'first', ignore.strand = T))
counts = counts[good.rows, ]
regions = regions[good.rows, ]

# Find regions with variants
ov = findOverlaps(regions.to.ranges(regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov))[1:10000, ]
counts = counts[ov.mat[, 1], ]
regions = regions[ov.mat[, 1], ]
genot.tmp = genot[ov.mat[, 2], ]

prefs = paste(rep(marks, times = 2), rep(c('left', 'right'), length(marks)), sep = '_')
true.cor = array(NaN, dim = c(nrow(regions), length(prefs)))
rand.genot = genot.tmp[seq(nrow(genot.tmp), 1, -1), ]
rand.cor = array(NaN, dim = c(nrow(regions), length(prefs)))
for(i in 1:length(prefs)){
  sel.counts = counts[, grep(prefs[i], colnames(counts))]
  good.rows = apply(sel.counts, 1, function(x) !any(is.na(x)))
  #sel.counts = sel.counts[good.rows, & rowMedians(sel.counts) > quantile(rowMedians(counts), 0.4, na.rm = T)
  true.cor[good.rows, i] = cor.par(genot.tmp[good.rows, ], sel.counts[good.rows, ], nchunk = 5, method = 'spearman')
  rand.cor[good.rows, i] = cor.par(rand.genot[good.rows, ], sel.counts[good.rows, ], nchunk = 5, method = 'spearman')
}
sum(rowMaxs(abs(rand.cor)) > 0.8 )

fit.dat = array(NaN, dim = c(nrow(counts), 6)) 

fit.dat = foreach(i = 1:nrow(counts), .combine = 'rbind', .inorder = T) %dopar% {
  y = genot.tmp[i, ] 
  x = t(matrix(counts[i, ], ncol = nindivs, byrow = T))
  
  # Use leave-one-out CV to pick the optimal lambda
  cv.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = 0.5, standardize = F, grouped = F, nfolds = length(y))
  k = which(cv.model$lambda == cv.model$lambda.min)
  #coef.dat[i, ] = append(as.numeric(cv.model$glmnet.fit$a0[k]), cv.model$glmnet.fit$beta[, k])
  c(deviance(cv.model$glmnet.fit)[k], cv.model$glmnet.fit$dev.ratio[k], cv.model$lambda.min,
    as.numeric(cv.model$glmnet.fit$a0[k]), cv.model$glmnet.fit$beta[, k])
}

rand.genot = genot.tmp[seq(nrow(genot.tmp), 1, -1), ]
rand.fit.dat = foreach(i = 1:nrow(counts), .combine = 'rbind', .inorder = T) %dopar% {
  y = rand.genot[i, ] 
  x = t(rbind(rank(counts[1,1:14]),rank(counts[1,15:28]))) #t(matrix(counts[i, ], ncol = nindivs, byrow = T))
  
  # Use leave-one-out CV to pick the optimal lambda
  cv.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = 0.5, standardize = F, grouped = F, nfolds = length(y))
  k = which(cv.model$lambda == cv.model$lambda.min)
  #coef.dat[i, ] = append(as.numeric(cv.model$glmnet.fit$a0[k]), cv.model$glmnet.fit$beta[, k])
  c(deviance(cv.model$glmnet.fit)[k], cv.model$glmnet.fit$dev.ratio[k], cv.model$lambda.min,
    as.numeric(cv.model$glmnet.fit$a0[k]), cv.model$glmnet.fit$beta[, k])
}
outfile = file.path(outdir, paste(outpref, '_', mark, '.RData', sep = ''))
if(!file.exists(outdir)) dir.create(outdir)
colnames(fit.dat) = c('deviance', 'dev.ratio', 'lambda', 'intercept')
tss = middles
save(ac.regions, ac.counts, rna.counts, gene.meta, fit.dat, coef.dat, tss, file = outfile)  