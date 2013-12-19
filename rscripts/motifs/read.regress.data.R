library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(glmnet)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# For each pair of individuals, reads text files with log-ratios of motif scores
# and creates a matrix regions-by-motifs. This is then used in get.motif.features.R for
# correlating the log-ratios of chromatin signal to log-ratios of motif scores.

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/'
outdir = file.path(indir, 'rdata')
if(!file.exists(outdir)) dir.create(outdir)
inpref = 'Gm12878_allTFBS.sorted.noPol.'
insuf1 = '.jaspar.scores'
insuf2 = '.pouya.scores'
score.files = list.files(indir, pattern = paste(inpref, '.*', insuf1, sep = ''), full.names = T)
nfiles = length(score.files)
tf.file1 = '../../rawdata/TFs/jaspar/motif_names.txt'
tfs1 = read.table(tf.file1, stringsAsFactors = F)[,1]
tf.file2 = '../../rawdata/TFs/pouya/motif_names.txt'
tfs2 = read.table(tf.file2, stringsAsFactors = F)[,1]

indivs_tmp = gsub(paste(inpref, '|', insuf1, '|.txt', sep = ''), '', basename(score.files))
indivs_tmp = unlist(strsplit(indivs_tmp, '_vs_'))
indivs1 = indivs_tmp[seq(1, length(indivs_tmp), 2)]
indivs2 = indivs_tmp[seq(2, length(indivs_tmp), 2)]
overwrite = T

for(i in 1:nfiles){
  score.file2 = gsub(insuf1, insuf2, score.files[i])
  outfile = file.path(outdir, paste(inpref, indivs1[i], '_vs_', indivs2[i], '_scores.RData', sep = ''))
  if(!file.exists(score.file2) || (file.exists(outfile) && !overwrite)) next
  cat(basename(score.files[i]), '\n')
  X = as.matrix(read.table(score.files[i]))
  X = cbind(X, as.matrix(read.table(score.file2)))
  colnames(X) = append(tfs1, tfs2)
  good.rows = !apply(X, 1, function(x) all(abs(x) < 0.001) | any(is.na(x)) | any(is.infinite(abs(x))))
  X = Matrix(X[good.rows, ])
  save(X, good.rows, file = outfile)
}
