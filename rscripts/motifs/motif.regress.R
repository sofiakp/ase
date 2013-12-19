rm(list=ls())
library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(glmnet)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# Regresses signal log-scores on motif log-ratios. See get.motif.features.R for a newer version
# that uses random forests instead.

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/'
inpref = 'Gm12878_allTFBS.sorted.noPol.'
insuf = '.jaspar.scores'
score.files = list.files(indir, pattern = paste(inpref, '.*', insuf, sep = ''), full.names = T)
nfiles = length(score.files)
tf.file = '../../rawdata/TFs/jaspar/motif_names.txt'
tfs = read.table(tf.file, stringsAsFactors = F)[,1]
#nsel = 10000

indivs_tmp = gsub(paste(inpref, '|', insuf, '|.txt', sep = ''), '', basename(score.files))
indivs_tmp = unlist(strsplit(indivs_tmp, '_vs_'))
indivs1 = indivs_tmp[seq(1, length(indivs_tmp), 2)]
indivs2 = indivs_tmp[seq(2, length(indivs_tmp), 2)]

# Only compute correlations and regression on the subset of regions that overlap Gm12878_<sel.tf>sorted.bed.
# If sel.tf is '', then consider all regions
sel.tf = 'Pu1.'
regions = read.bed(file.path(indir, paste(inpref, 'regions_H3K27AC.txt', sep = '')))
if(sel.tf  != ''){
  sel.regions = read.bed(paste('../../rawdata/TFs/Gm12878_', sel.tf, 'sorted.bed', sep = ''))
  ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(sel.regions), select = 'first', ignore.strand = T)
  sel.rows = which(!is.na(ov))  
}else{
  sel.rows = 1:length(regions)
}

X = NULL
y = c()
i2 = c()

for(i in 1:(nfiles-1)){
  signal.file = gsub(insuf, '_H3K27AC', score.files[i])
  if(!file.exists(signal.file)) next
  cat(basename(score.files[i]), '\n')
  X_tmp = as.matrix(read.table(score.files[i]))
  stopifnot(nrow(X_tmp) == )
  y_tmp = read.table(signal.file)[, 1]
  #sel.rows = 1:nrow(X_tmp) #sample(1:nrow(X_tmp), nsel)
  X = rbind(X, X_tmp[sel.rows, ])
  y = append(y, y_tmp[sel.rows])
  i2 = append(i2, rep(indivs2[i], length(sel.rows)))
}

i2 = factor(i2)
colnames(X) = tfs

good.col = !apply(X, 2, function(x) all(abs(x) < 0.001 | is.na(x)))
good.row = !apply(X, 1, function(x) all(abs(x) < 0.001) | any(is.na(x)) | any(is.infinite(abs(x)))) & !is.na(y) & !is.infinite(abs(y))
X = X[good.row, good.col]
y = y[good.row]
i2 = i2[good.row]

spear.cor = array(0, dim = c(ncol(X), 1))
pear.cor = array(0, dim = c(ncol(X), 1))
for(i in 1:ncol(X)){
  spear.cor[i] = cor(X[, i], y, method = 'spearman')
  pear.cor[i] = cor(X[, i], y, method = 'pearson')
}
cor.dat = data.frame(names = colnames(X), spear = spear.cor, pear = pear.cor)
cor.dat = cor.dat[order(cor.dat$spear), ]
cor.dat$names = ordered
p1 = ggplot(cor.dat) + geom_bar(aes(x = names, y = spear))

model = model.matrix(~ X)
m = glm(y ~ model)
cv = cv.glmnet(y = y, x = model, family = 'gaussian', alpha = 0.5, standardize = T, grouped = F, nfolds = 10)
coef = sort(cv$glmnet.fit$beta[, which(cv$lambda == cv$lambda.min)])
names(coef) = gsub('^X|^i2', '', names(coef))
coef.dat = data.frame(feat = ordered(factor(names(coef), levels = names(coef))), coef = coef)
p = ggplot(coef.dat) + geom_bar(aes(x = feat, y = coef), stat = 'identity') + 
  theme_bw() + theme(axis.text.x = element_text(size = 6, angle = -45, vjust = 1, hjust = 0), axis.text.y = element_text(size = 6))
ggsave(file.path(indir, 'gmlnet_test_coef.pdf'), p, width = 12, height = 9)
