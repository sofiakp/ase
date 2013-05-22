rm(list=ls())
library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(glmnet)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/'
inpref = 'Gm12878_allTFBS.sorted.noPol.'
insuf = '.jaspar.scores'
score.files = list.files(indir, pattern = paste(inpref, '.*12878.*', insuf, sep = ''), full.names = T)
nfiles = length(score.files)
tf.file = '../../rawdata/TFs/jaspar/motif_names.txt'
tfs = read.table(tf.file, stringsAsFactors = F)[,1]
nsel = 10000

indivs_tmp = gsub(paste(inpref, '|', insuf, '|.txt', sep = ''), '', basename(score.files))
indivs_tmp = unlist(strsplit(indivs_tmp, '_vs_'))
indivs1 = indivs_tmp[seq(1, length(indivs_tmp), 2)]
indivs2 = indivs_tmp[seq(2, length(indivs_tmp), 2)]

X = NULL
y = c()
i2 = c()

for(i in 1:nfiles){
  signal.file = gsub(insuf, '_H3K27AC', score.files[i])
  if(!file.exists(signal.file)) next
  cat(basename(score.files[i]), '\n')
  X_tmp = as.matrix(read.table(score.files[i]))
  y_tmp = read.table(signal.file)[, 1]
  sel.rows = sample(1:nrow(X_tmp), nsel)
  X = rbind(X, X_tmp[sel.rows, ])
  y = append(y, y_tmp[sel.rows])
  i2 = append(i2, rep(indivs2[i], length(sel.rows)))
}

i2 = factor(i2)
colnames(X) = tfs

good.col = !apply(X, 2, function(x) all(x == 0))
good.row = !apply(X, 1, function(x) all(x == 0) | any(is.na(x)) | any(is.infinite(x))) & !is.na(y) & !is.infinite(abs(y))
X = X[good.row, good.col]
y = y[good.row]
i2 = i2[good.row]

model = model.matrix(~ X + i2)
m = glm(y ~ model)
cv = cv.glmnet(y = y, x = model, family = 'gaussian', alpha = 0.5, standardize = T, grouped = F, nfolds = 10)
coef = sort(cv$glmnet.fit$beta[, which(cv$lambda == cv$lambda.min)])
names(coef) = gsub('^X', '', names(coef))
coef.dat = data.frame(feat = ordered(factor(names(coef), levels = names(coef))), coef = coef)
p = ggplot(coef.dat) + geom_bar(aes(x = feat, y = coef), stat = 'identity') + 
  theme_bw() + theme(axis.text.x = element_text(size = 12, angle = -45, vjust = 1, hjust = 0), axis.text.y = element_text(size = 12))
#X = X + matrix(runif(nrow(X) * ncol(X), -1, 1), nrow = nrow(X)) 
