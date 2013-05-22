rm(list=ls())
library(ggplot2)
library(Matrix)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(foreach)
library(doMC)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# Gets the population regions of a mark and plots the signal at the overlapping genes.
plotdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
rdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/'
mark = 'H3K27AC'
pref = 'SNYDER_HG19_all_reg_'
comp = 3
qval = 0.01
K = 4
pref = paste(pref, mark, '_comp', comp, '_q', qval, sep = '')

load(file.path(rdir, paste(pref, '_qn_isvaNull.RData', sep = '')))
pop.regions = regions[isva.fit$deg, ]
clust.pref = paste(pref, '_K', K, sep = '')
load(file.path(rdir, paste(clust.pref, '_clust.RData', sep = '')))
pop.regions = pop.regions[plot.rows, ]

# Get gene annotations
c2 = new.env()
load('../../rawdata/geneCounts/rdata/repsComb/rdata/SNYDER_HG19_all_reg_RZ_qn.RData', c2)
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
# Remove these two rows if you want to plot the uncorrected signal.
gene.meta = gene.meta[c2$good.rows, ]
load('../../rawdata/geneCounts/rdata/repsComb/rdata/SNYDER_HG19_all_reg_RZ_comp1_q0.01_qn_isvaNull.RData', c2)

ov = findOverlaps(regions.to.ranges(pop.regions), regions.to.ranges(gene.meta), select = 'first', ignore.strand = T)
has.ov = !is.na(ov)
clust.tab = table(kclusters$cluster[plot.rows][has.ov])
ov = ov[!is.na(ov)]
matched.col = match(colnames(counts)[plot.cols], colnames(c2$counts))
counts2 = c2$counts[ov, matched.col[!is.na(matched.col)]]
counts2 = t(scale(t(counts2)))
h2 = plot.tile(counts2, x.ord.samples = colnames(counts2), xcolor = get.pop.col(get.pop(colnames(counts2))), 
               ysep = cumsum(clust.tab)[-K] + 1, ylabels = array('', dim = c(K-1,1)), lcex = 12, xcex = 13)
ggsave(file.path(plotdir, paste(clust.pref, '_vs_genes_biclust.pdf', sep = '')), h2, width = 7, height = 6.8)