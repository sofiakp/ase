rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
source('utils/sample.info.r')
source('utils/deseq.utils.r')
source('utils/binom.val.r')
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')

# This snippet selects and outputs differential genes
# load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
# deseq.dir = '../../rawdata/geneCounts/rdata/repsComb/deseq/'
# diff.counts = get.diff.count(list.files(deseq.dir, pattern = '*RZ_deseq.RData', full.names = T), 0.01)
# genes = unique(gene.meta$gene.name[diff.counts > 0])

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

files = c()
deseq.dir = '../../rawdata/signal/rep/countsAtPeaksBroad/merged_Mar13/repsComb/deseq2/'
files = list.files(deseq.dir, pattern = '*_deseq2.RData', full.names = T)
deseq.dir = '../../rawdata/genomeGrid/hg19_w10k/rep/counts_newNorm/repsComb/deseq2/'
files = append(files, list.files(deseq.dir, pattern = '*_deseq2.RData', full.names = T))
deseq.dir = '../../rawdata/transcriptomes/rep/counts_newNorm/repsComb/deseq2/'
files = append(files, list.files(deseq.dir, pattern = '*_deseq2.RData', full.names = T))
deseq.dir = '../../rawdata/geneCounts/rdata/repsComb/deseq2/'
files = append(files, list.files(deseq.dir, pattern = '*RZ_deseq2.RData', full.names = T))

files = files[!grepl('BUB|POL', files)]
plotdir = '../../rawdata/signal/rep/countsAtPeaksBroad/merged_Mar13/repsComb/plots/' #'../../rawdata/geneCounts/rdata/repsComb/deseq2/plots'
if(!file.exists(plotdir)) dir.create(plotdir)

# We will use the pop-clustering data to reorder individuals if possible
isva.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/'
isva.files = c() #list.files(isva.dir, pattern = 'SNYDER_HG19_all_reg_[a-zA-Z0-9]*_comp3_q0.01_K4_clust.RData', full.names = T)
isva.files = isva.files[!grepl('rand', isva.files)]

uniq.indivs = get.indivs()

marks = unlist(strsplit(basename(files), '_'))
marks = marks[seq(3, length(marks), 4)]
uniq.marks = unique(marks)
diff.all = NULL
diff.all.norm = NULL

for(m in uniq.marks){
  sel.files = files[grep(m, files)]
  print(length(sel.files))
  diff.dat = get.diff.count(sel.files, 1e-5, is.log = T, fold.cut = log2(3))
  diff.counts = diff.dat$diff.count
  
  diff.sum = as.matrix(table(diff.counts)) # How many regions differ in X pairs of individuals? 
  n = as.numeric(rownames(diff.sum)) 
  diff.sum = cumsum(diff.sum[seq(length(diff.sum), 1, -1)])[seq(length(diff.sum), 1, -1)] # How many regions differ in >= X pairs of individuals?
  diff.sum = diff.sum[n > 0 & n < 21]
  diff.sum = data.frame(pairs = 1:20, count = diff.sum)
  diff.sum$mark = rep(m, nrow(diff.sum))
  diff.all = rbind(diff.all, diff.sum)
  diff.sum$count = diff.sum$count / length(diff.counts)
  diff.all.norm = rbind(diff.all.norm, diff.sum)
  
  pair.diff = diff.dat$pair.diff
  colnames(pair.diff) = fix.indiv.names(colnames(pair.diff))
  rownames(pair.diff) = fix.indiv.names(rownames(pair.diff))  
  frac = pair.diff*100/length(diff.counts)
  sel.pop.file = isva.files[grep(m, isva.files)]
  if(length(sel.pop.file) > 0){
    load(sel.pop.file) 
    load(gsub('_comp.*clust', '_qn', sel.pop.file))
    out.indiv = colnames(counts)[plot.cols]
    indiv.idx = match(out.indiv, colnames(frac))
  }else{
    # No population biclustering. Just sort by population.
    out.indiv = uniq.indivs
    indiv.idx = match(out.indiv[out.indiv %in% colnames(frac)], colnames(frac))
  }
  frac = frac[indiv.idx, indiv.idx]
  plot.heatmap(frac, filt.thresh = NA, row.cluster = F, show.dendro = "none", col.cluster = F, row.title= '', col.title = '', 
               cellnote = matrix(sprintf('%.2f', frac), ncol = nrow(frac)), 
               ColSideColors = get.pop.col(get.pop(colnames(frac))), RowSideColors = get.pop.col(get.pop(colnames(frac))), 
               cex.row = 1.4, cex.col = 1.4, margins = c(9, 9), keysize = 1,
               dist.metric='spearman', clust.method = "average", break.type='linear', palette = brewer.pal(9, 'Reds')[2:7],
               to.file = file.path(plotdir, paste(m, 'num_diff.pdf', sep = '_')))
}
diff.all$mark = order.marks(diff.all$mark)
diff.all.norm$mark = order.marks(diff.all.norm$mark)
p = ggplot(diff.all, aes(group = mark)) + geom_line(aes(x = pairs, y = count, color = mark), size = 1) + 
  geom_point(aes(x = pairs, y = count, color = mark, shape = mark), size = 5) + 
  xlab('Number of pairwise comparisons') + ylab('Number of regions varying in at least X pairs') + theme_bw() + 
  scale_x_continuous(breaks = 1:20) + #scale_color_discrete('') + 
  scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17)) + 
  scale_color_manual(values =  mark.colors(unique(diff.all$mark))) + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        legend.position = c(.85, .75), legend.text = element_text(size = 14), legend.title = element_blank())
ggsave(file.path(plotdir, 'num_diff_pairs.pdf'), p, width = 6.5, height = 5.6)

q = ggplot(diff.all.norm, aes(group = mark)) + geom_line(aes(x = pairs, y = count, color = mark), size = 1) + 
  geom_point(aes(x = pairs, y = count, color = mark, shape = mark), size = 5) + 
  xlab('Number of pairwise comparisons') + ylab('Fraction of regions varying in at least X pairs') + theme_bw() + 
  scale_x_continuous(breaks = 1:20) + 
  scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17)) + 
  scale_color_manual(values =  mark.colors(unique(diff.all$mark))) + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        legend.position = c(.87, .75), legend.text = element_text(size = 14), legend.title = element_blank())
ggsave(file.path(plotdir, 'frac_diff_pairs.pdf'), q, width = 6.5, height = 5.6)