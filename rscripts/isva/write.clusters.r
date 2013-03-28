rm(list=ls())
library(ggplot2)
library(Matrix)
library(GenomicRanges)
source('utils/deseq.utils.r')

# Writes the regions in the population-specific clusters together with their signal

load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/SNYDER_HG19_all_reg_H3K27AC_qn_isvaNull.RData')
load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/SNYDER_HG19_all_reg_H3K27AC_qn_isvaNull_clust.RData')
outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/qn_isvaNull_fits_all_reg/'
outpref = 'SNYDER_HG19_all_reg_'
mark = 'H3K27AC'
k = max(kclusters$cluster)
regions = regions[isva.fit$deg, ]
counts = counts[isva.fit$deg, ]
colnames(counts)[colnames(counts) == 'SNYDER'] = 'MS1'
for(i in 1:k){
  sel = which(kclusters$cluster == kord[i])
  outfile = file.path(outdir, paste(outpref, 'isva_', mark, '_k',k, '_clust', kord[i], '_withSignal.txt', sep = ''))
  dat = data.frame(chr = regions$chr, start = regions$start - 1, end = regions$end)[sel, ]
  dat = cbind(dat, counts[sel, plot.cols])
  write.table(dat, file = outfile, quote = F, row.names = F, col.names = T, sep = "\t")
}
