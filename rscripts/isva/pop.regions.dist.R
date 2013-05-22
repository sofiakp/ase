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

# Tests whether pop-specific regions tend to cluster together or with pop-specific regions of other clusters

get.closest = function(dists, ov.mat){
  dat.tmp = data.frame(d = dists, region.idx1 = ov.mat[, 1], region.idx2 = ov.mat[, 2], k = ov.mat[, 3])
  cor.dat.tmp.2 = cast(dat.tmp, region.idx1~., function(x) min(x), value = 'd')
  colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
  # create matrix region, cor, snp_chr, snp_pos, reg_start, reg_end. 
  # This might have more than one SNP per region, if there were multiple SNPs in the same region
  # with equally high correlation
  cor.dat.true = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor')) 
  cor.dat.true = cor.dat.true[!duplicated(cor.dat.true[,1:2]), ]
  return(cor.dat.true)
}

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
no.pop.regions = regions[-isva.fit$deg, ]
all.regions = rbind(pop.regions, no.pop.regions)
indivs = colnames(counts)
nindivs = length(indivs)
clust.pref = paste(pref, '_K', K, sep = '')
load(file.path(rdir, paste(clust.pref, '_clust.RData', sep = '')))
all.clust = append(kclusters$cluster, array(0, dim = c(nrow(no.pop.regions), 1)))

# Get gene annotations
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
deseq.dir = '../../rawdata/geneCounts/rdata/repsComb/deseq2/'
files = list.files(deseq.dir, pattern = '*RZ_deseq2.RData', full.names = T)
diff.dat = get.diff.count(files, 1e-5, is.log = T, fold.cut = log2(3))

bino.p = array(NaN, dim = c(K + 1, K + 1))
counts = NULL
for(k in 0:K){
  sel.regions = all.regions[all.clust == k, ]
  comp.clust = 0:4
  minr = min(nrow(sel.regions), min(table(all.clust[all.clust %in% comp.clust])))
  sel.regions = sel.regions[sample(1:nrow(sel.regions), minr), ]
  other.regions = NULL
  other.clust = c()
  for(j in comp.clust){
    other.regions = rbind(other.regions, all.regions[all.clust == j, ][sample(1:sum(all.clust == j), minr), ])  
    other.clust = append(other.clust, array(j, dim = c(minr, 1)))
  }

  # Find all the regions in other.regions that are in a 1Mb window from each of the sel.regions.
  tmp.regions = sel.regions
  tmp.regions$start = pmax(sel.regions$start - 500000, 1)
  tmp.regions$end = sel.regions$end + 500000
  ov = findOverlaps(regions.to.ranges(tmp.regions), regions.to.ranges(other.regions), select = 'all', ignore.strand = T)
  ov.mat = cbind(queryHits(ov), subjectHits(ov), other.clust[subjectHits(ov)])
  
  # Get the distances between sel.regions and other.regions and then select the closest other.region for each sel.region.
  dists = pmax(sel.regions$start[ov.mat[, 1]], other.regions$start[ov.mat[, 2]]) - pmin(sel.regions$end[ov.mat[, 1]], other.regions$end[ov.mat[, 2]])
  self.ov = dists <= 0
  ov.mat = ov.mat[!self.ov, ]
  dists = dists[!self.ov]
  
  dat.tmp = data.frame(d = dists, region.idx1 = ov.mat[, 1], region.idx2 = ov.mat[, 2], k = ov.mat[, 3])
  dat.tmp.2 = cast(dat.tmp, region.idx1~., function(x) min(x), value = 'd')
  colnames(dat.tmp.2) = c('region.idx1', 'd')
  dat.out = merge(dat.tmp, dat.tmp.2, by = c('region.idx1', 'd'))
  
  counts = rbind(counts, table(dat.out$k))
  for(j in comp.clust){
    bino.p[k + 1, j + 1] = binom.test(counts[k + 1, j + 1], minr, 1/length(comp.clust), 'greater')$p.value
  }
}

# Find overlaps between regions and genes. 
ov = findOverlaps(regions.to.ranges(all.regions), regions.to.ranges(gene.meta), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov), all.clust[queryHits(ov)]) # region idx, gene idx, cluster idx of region
rand.genes = unique(ov.mat[ov.mat[, 3] == 0, 2]) # genes overlapping non population specific enhancers

# For each cluster i test if genes with regions from i are more likely to also have regions from cluster j than genes
# overlapping non-pop regions.
gene.enrich = array(NaN, dim = c(K, K)) 
gene.bino.p = array(NaN, dim = c(K, K)) 
gene.wilc.p = array(NaN, dim = c(K, K)) 
for(i in 1:K){
  for(j in 1:K){
    if(i == j) next
    sel.genes = unique(ov.mat[ov.mat[, 3] == i, 2])
    sel.genes2 = unique(ov.mat[ov.mat[, 3] == j, 2])
    common.genes = intersect(sel.genes, sel.genes2)
    #common.genes.ind = array(F, dim = c(nrow(gene.meta), 1))
    #common.genes.ind[common.genes] = T
    rand.common.genes = intersect(rand.genes, sel.genes2)
    gene.wilc.p[i, j] = wilcox.test(diff.dat$diff.count[rand.common.genes], diff.dat$diff.count[common.genes])$p.value   
    gene.enrich[i, j] = (length(common.genes)/length(sel.genes)) / (length(rand.common.genes)/length(rand.genes))
    gene.bino.p[i, j] = binom.test(length(common.genes), length(sel.genes), length(rand.common.genes)/length(rand.genes), 'greater')$p.value
  }
}

tab = rbind(c(length(common.genes), length(setdiff(sel.genes2, sel.genes))), c(length(setdiff(sel.genes, sel.genes2)), nrow(gene.meta) - length(union(sel.genes, sel.genes2))))

       # dat.tmp = data.frame(d = dists, region.idx1 = ov.mat[, 1], region.idx2 = ov.mat[, 2])
# dat.tmp.2 = cast(dat.tmp, region.idx1~., function(x) min(x), value = 'd')
# colnames(dat.tmp.2) = c('region.idx1', 'd')
# dat.out = merge(dat.tmp, dat.tmp.2, by = c('region.idx1', 'd')) 
# 
# sum(clust[dat.out$region.idx2] != 0)
# 
# sel.regions = pop.regions[kclusters$cluster == 1, ]
# other.regions = rbind(pop.regions[kclusters$cluster != 1, ], no.pop.regions)
# clust = array(0, dim = c(nrow(other.regions), 1))
# clust[1:sum(kclusters$cluster != 1)] = kclusters$cluster[kclusters$cluster != 1]
# tmp.regions = sel.regions
# tmp.regions$start = pmax(sel.regions$start - 500000, 1)
# tmp.regions$end = sel.regions$end + 500000
# ov = findOverlaps(regions.to.ranges(tmp.regions), regions.to.ranges(other.regions), select = 'all', ignore.strand = T)
# ov.mat = cbind(queryHits(ov), subjectHits(ov))
# dists = pmax(sel.regions$start[ov.mat[, 1]], other.regions$start[ov.mat[, 2]]) - pmin(sel.regions$end[ov.mat[, 1]], other.regions$end[ov.mat[, 2]])
# 
# dat.tmp = data.frame(d = dists, region.idx1 = ov.mat[, 1], region.idx2 = ov.mat[, 2])
# dat.tmp.2 = cast(dat.tmp, region.idx1~., function(x) min(x), value = 'd')
# colnames(dat.tmp.2) = c('region.idx1', 'd')
# dat.out = merge(dat.tmp, dat.tmp.2, by = c('region.idx1', 'd')) 
# 
# sum(clust[dat.out$region.idx2] != 0)