rm(list=ls())
library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(glmnet)
library(randomForest)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

signal.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/'
indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/rdata/'
score.files = list.files(indir, pattern = '.*_vs_.*scores.RData', full.names = T)
nfiles = length(score.files)
insuf = '_scores.RData'
inpref = 'Gm12878_allTFBS.sorted.noPol.'
indivs_tmp = gsub(paste(inpref, '|', insuf, '|.RData', sep = ''), '', basename(score.files))
indivs_tmp = unlist(strsplit(indivs_tmp, '_vs_'))
file.indivs1 = indivs_tmp[seq(1, length(indivs_tmp), 2)]
file.indivs2 = indivs_tmp[seq(2, length(indivs_tmp), 2)]

# Only compute correlations and regression on the subset of regions that overlap Gm12878_<sel.tf>sorted.bed.
# If sel.tf is '', then consider all regions
sel.tf = ''
regions = read.bed(file.path(signal.dir, paste(inpref, 'regions_H3K27AC.txt', sep = '')))
if(sel.tf  != ''){
  sel.regions = read.bed(paste('../../rawdata/TFs/Gm12878_', sel.tf, 'sorted.bed', sep = ''))
  ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(sel.regions), select = 'first', ignore.strand = T)
  sel.rows = !is.na(ov)  
}else{
  sel.rows = array(T, dim = c(nrow(regions), 1))
}

X = NULL
y = c()
ridx = c() # region index
indivs1 = c()
indivs2 = c()
 
for(i in 1:nfiles){
  signal.file = file.path(signal.dir, gsub(insuf, '_H3K27AC_v2.txt', basename(score.files[i])))
  if(!file.exists(signal.file)) next
  cat(basename(score.files[i]), '\n')
  env = new.env()
  load(score.files[i], env)
  stopifnot(length(env$good.rows) == nrow(regions))
  tmp.sel = sel.rows[env$good.rows]
  y.tmp = read.table(signal.file)[, 1]
  stopifnot(length(y.tmp) == nrow(regions))
  y.tmp = y.tmp[env$good.rows]
  #sel.rows = 1:nrow(X_tmp) #sample(1:nrow(X_tmp), nsel)
  if(is.null(X)){
    X = env$X[tmp.sel, ]
    y = y.tmp[tmp.sel]
    ridx = which(sel.rows & env$good.rows)
    indivs1 = rep(file.indivs1[i], length(y.tmp))
    indivs2 = rep(file.indivs2[i], length(y.tmp))
  }else{
    X = rBind(X, env$X[tmp.sel, ])
    y = append(y, y.tmp[tmp.sel]) 
    ridx = append(ridx, which(sel.rows & env$good.rows))
    indivs1 = append(indivs1, rep(file.indivs1[i], length(y.tmp)))
    indivs2 = append(indivs2, rep(file.indivs2[i], length(y.tmp)))
  }
}

indivs1 = factor(indivs1)
indivs2 = factor(indivs2)

#good.col = !apply(X, 2, function(x) all(abs(x) < 0.001 | is.na(x)))
#good.row = !apply(X, 1, function(x) all(abs(x) < 0.001) | any(is.na(x)) | any(is.infinite(abs(x)))) & !is.na(y) & !is.infinite(abs(y))
#X = X[good.row, good.col]
#y = y[good.row]
save(X, y, indivs1, indivs2, ridx, regions, file = file.path(indir, paste(inpref, 'scores_v2.RData', sep = '')))

load(file.path(indir, paste(inpref, 'scores.RData', sep = '')))
good.y = which(abs(y) > 2)
sel.train = sample(good.y, 50000)
sel.test = sample(good.y, 1000)

rf = randomForest(as.matrix(X[sel.train, ]), y[sel.train], as.matrix(X[sel.test, ]), y[sel.test], ntree = 100, nodesize = 2, maxnodes = 10, importance = T, keep.forest = T)
