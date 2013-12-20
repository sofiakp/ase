rm(list=ls())
library(reshape)
library(GenomicRanges)
library(ggplot2)
source('utils/deseq.utils.r')

var.files = c()
vardir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/anova/myfiles/'
outdir = vardir
var.files = append(var.files, list.files(vardir, pattern = 'withF.RData', full.names = T))
vardir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig_newNorm/rdata/anova/myfiles/'
var.files = append(var.files, list.files(vardir, pattern = 'withF.RData', full.names = T))
vardir = '../../rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig_newNorm/rdata/anova/myfiles/'
var.files = append(var.files, list.files(vardir, pattern = 'withF.RData', full.names = T))

enrich = array(0, dim = c(length(var.files), 2))
marks = array('', , dim = c(length(var.files), 1))

all.dat = NULL

for(i in 1:length(var.files)){
  base = basename(var.files[i])
  mark = gsub('SNYDER_HG19_|_qn_withF.RData', '', base)
  marks[i] = mark
  
  load(var.files[i])
  mean.sig = rowMeans(counts)
  var = !is.na(regions$logF) & regions$logF > 1
  non.var = !is.na(regions$logF) & regions$logF < 1
  var.sig = mean.sig[var]
  non.var.sig = mean.sig[non.var]
  
  enrich[i, 1] = wilcox.test(var.sig, non.var.sig)$p.value
  enrich[i, 2] = median(var.sig) / median(non.var.sig)
  
  dat = data.frame(c = append(var.sig, non.var.sig), t = rep(c('Variable', 'Non-variable'), c(sum(var), sum(non.var))), 
                   m = rep(mark, sum(var) + sum(non.var)))
  all.dat = rbind(all.dat, dat)
}

rownames(enrich) = marks
write.table(enrich, file.path(outdir, 'var_sig_enrich.txt'), sep = '\t', row.names = T, col.names = F, quote = F)
enrich = data.frame(enrich)
colnames(enrich) = c('wpval', 'wenrich')
enrich$mark = marks

p = ggplot(all.dat) + geom_boxplot(aes(x = t, y = c), size = 0.5, outlier.size = 1) + facet_wrap(~m, scales = 'free_y') + ylab('Average region signal') + xlab('') +
  theme_bw() + theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13), 
                     axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                     strip.text = element_text(size = 15))
ggsave(file.path(outdir, 'var_sig_enrich_box.pdf'), p, width = 8, height = 8)