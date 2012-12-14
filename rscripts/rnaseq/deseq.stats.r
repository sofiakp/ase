rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

# deseq.dir = '../../rawdata/geneCounts/rdata/repsComb/deseq/'
# plotdir = '../../rawdata/geneCounts/rdata/repsComb/plots'
# if(!file.exists(plotdir)) dir.create(plotdir)
# diff.counts = get.diff.count(list.files(deseq.dir, pattern = '*RZ_deseq.RData', full.names = T), 0.0001)
# 
# diff.sum = as.matrix(table(diff.counts)) # How many regions differ in X pairs of individuals? 
# n = as.numeric(rownames(diff.sum)) 
# diff.sum = cumsum(diff.sum[seq(length(diff.sum), 1, -1)])[seq(length(diff.sum), 1, -1)] # How many regions differ in >= X pairs of individuals?
# diff.sum = diff.sum[n > 0 & n < 11]
# diff.sum = data.frame(pairs = 1:10, count = diff.sum)
# p = ggplot(diff.sum) + geom_line(aes(x = pairs, y = count), size = 1) + geom_point(aes(x = pairs, y = count), size = 5) +
#   xlab('Number of pairwise comparisons') + ylab('Number of genes varying in at least X pairs') + theme_bw() + 
#   scale_x_continuous(breaks = 1:10) +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
#         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
# ggsave(file.path(plotdir, 'RZ_num_diff_pairs.png'), p, width = 6.5, height = 5.6)

deseq.dir = '../../rawdata/signal/combrep/countsAtPeaksBroad/repsComb/deseq/'
plotdir = '../../rawdata/signal/combrep/countsAtPeaksBroad/repsComb/plots/'
if(!file.exists(plotdir)) dir.create(plotdir)

marks = c('H3K27AC', 'H3K27ME3', 'H3K4ME3', 'H3K4ME1')
diff.all = NULL
diff.all.norm = NULL
for(m in marks){
  files = list.files(deseq.dir, pattern = paste(m, '*_deseq.RData', sep = ''), full.names = T)
  print(length(files))
  diff.counts = get.diff.count(files, 0.0001)
  
  diff.sum = as.matrix(table(diff.counts)) # How many regions differ in X pairs of individuals? 
  n = as.numeric(rownames(diff.sum)) 
  diff.sum = cumsum(diff.sum[seq(length(diff.sum), 1, -1)])[seq(length(diff.sum), 1, -1)] # How many regions differ in >= X pairs of individuals?
  diff.sum = diff.sum[n > 0 & n < 11]
  diff.sum = data.frame(pairs = 1:10, count = diff.sum)
  diff.sum$mark = rep(m, nrow(diff.sum))
  diff.all = rbind(diff.all, diff.sum)
  diff.sum$count = diff.sum$count / length(diff.counts)
  diff.all.norm = rbind(diff.all.norm, diff.sum)
}
p = ggplot(diff.all) + geom_line(aes(x = pairs, y = count, color = mark), size = 1) + geom_point(aes(x = pairs, y = count, color = mark), size = 5) + 
  xlab('Number of pairwise comparisons') + ylab('Number of regions varying in at least X pairs') + theme_bw() + 
  scale_x_continuous(breaks = 1:10) + scale_color_discrete('') + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        legend.position = c(.86, .82))
ggsave(file.path(plotdir, 'num_diff_pairs.png'), p, width = 6.5, height = 5.6)

p = ggplot(diff.all.norm) + geom_line(aes(x = pairs, y = count, color = mark), size = 1) + geom_point(aes(x = pairs, y = count, color = mark), size = 5) + 
  xlab('Number of pairwise comparisons') + ylab('Fraction of regions varying in at least X pairs') + theme_bw() + 
  scale_x_continuous(breaks = 1:10) + scale_color_discrete('') + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        legend.position = c(.86, .82))
ggsave(file.path(plotdir, 'frac_diff_pairs.png'), p, width = 6.5, height = 5.6)