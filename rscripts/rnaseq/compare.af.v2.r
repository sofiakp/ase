rm(list=ls())
library(ggplot2)
library(Matrix)
library(reshape)
library(GenomicRanges)
library(matrixStats)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('utils/binom.val.r')

compute.dist = function(data){
  indivs = rownames(data)
  nindiv = length(indivs)
  pop = get.pop(indivs)
  uniq.pop = unique(pop)
  npop = length(uniq.pop)
  d = as.matrix(dist(data, method = 'euclidean'))
  pop.d = array(0, dim = c(npop, npop))
  for(i in 1:(nindiv-1)){
    for(j in (i+1):nindiv){
      p1 = match(pop[i], uniq.pop)
      p2 = match(pop[j], uniq.pop)
      pop.d[p1, p2] = (pop.d[p1, p2] + d[i, j])
    }
  }
  pop.d.tmp = triu(pop.d, 1) + tril(pop.d, -1)
  pop.d = pop.d + t(pop.d.tmp) # Add upper triangle to lower triangle without duplicating the diagonal
  colnames(pop.d) = uniq.pop
  rownames(pop.d) = uniq.pop
  norm.pop.d = pop.d
  for(i in 1:npop){
    for(j in 1:npop){
      norm.pop.d[i,j] = pop.d[i,j] / (pop.d[i,i] + pop.d[j,j])  
    }
  } 
  return(list(pop.d = as.matrix(pop.d), norm.pop.d = as.matrix(norm.pop.d)))
}

set.seed(1)
plotdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
rdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/'
mark = 'H3K27AC'
pref = 'SNYDER_HG19_all_reg_'
comp = 3
qval = 0.05
K = 6
nrand = 100
pref = paste(pref, mark, '_comp', comp, '_q', qval, sep = '')
clust.pref = paste(pref, '_K', K, sep = '')

load(file.path(rdir, paste(pref, '_qn_isvaNull.RData', sep = '')))
regions = regions[isva.fit$deg, ]
load(file.path(rdir, paste(clust.pref, '_clust.RData', sep = '')))

load('../../rawdata/variants/all_Mar13/genot.RData')
load('../../rawdata/variants/all_Mar13/snps.RData')
     
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')
genot = genot[, colnames(genot) != 'GM19193']
nindivs = ncol(genot)
#sel = rowSums(genot == 0) < nindivs - 2 & rowSums(genot == 1) < nindivs - 2 & rowSums(genot == 2) < nindivs - 2
#genot.var = genot[sel, ]
#snp.pos.var = snp.pos[sel, ]
ov = findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(gene.meta), select = 'first', ignore.strand = T) # Restrict to exons
snp.pos.var = snp.pos[!is.na(ov), ]
genot.var = genot[!is.na(ov), ]
ov = findOverlaps(snps.to.ranges(snp.pos.var), regions.to.ranges(regions), select = 'first', ignore.strand = T) 
snp.pos.ov = snp.pos.var[!is.na(ov), ]
genot.ov = genot.var[!is.na(ov), ]
genot.indivs = fix.indiv.names(colnames(genot))

#genot.norm = scale(genot[sample(1:nrow(genot), 1200), ])
#pca.fit = prcomp(t(genot.norm[, !(colnames(genot) %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
#p = plot.pcs(t(genot.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = colnames(genot), groups = get.pop(colnames(genot)), all = F, ndim = 4)
#pca.plot.all = p$p1
#bg.dist = compute.dist(t(genot.norm) %*% pca.fit$rotation[, 1:3])$norm.pop.d

pca.plots = list()
d.all = NULL
pops = unique(get.pop(colnames(counts)[plot.cols]))
for(k in 1:K){
  sel = kclusters$cluster == k
  # Get all the SNPs overlapping regions of the k-th cluster
  ov = findOverlaps(regions.to.ranges(regions[sel, ]), snps.to.ranges(snp.pos.ov),  select = 'all', ignore.strand = T)
  genot.norm = scale(genot.ov[subjectHits(ov), ])
  pca.fit = prcomp(t(genot.norm[, !(genot.indivs %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
  p = plot.pcs(t(genot.norm) %*% pca.fit$rotation, pca.fit$rotation, pca.fit$sdev, labels = genot.indivs, groups = get.pop(genot.indivs), all = T, ndim = 4)
  pca.plots[[k]] = p$p1
  ggsave(file.path(plotdir, paste(clust.pref, '_clust', k, '_genotPca.pdf', sep = '')), p$p1, width = 9, height = 6.8)
  d = compute.dist(t(genot.norm) %*% pca.fit$rotation[, 2:4])$pop.d
  
  avg.d = array(0, dim = dim(d))
  for(i in 1:nrand){
    genot.norm.bg = scale(genot.var[sample(1:nrow(genot.var), nrow(genot.norm), replace = T), ])
    pca.fit.bg = prcomp(t(genot.norm.bg[, !(genot.indivs %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
    avg.d = avg.d + compute.dist(t(genot.norm.bg) %*% pca.fit.bg$rotation[, 2:4])$pop.d
  }
  d = d / (avg.d / nrand)
  d = d[match(pops, rownames(d)), match(pops, colnames(d))]
  print(d)
  d.bin = as.matrix(triu(d, 1))
  d.bin[d.bin < max(rowMaxs(d.bin))] = 0
  d.bin[d.bin > 2] = 2
  p = plot.tile(d.bin, midpoint = 1, y.ord.samples = pops[seq(length(pops), 1, -1)], 
               ycex = 25, xcex = 25, ycolor = get.pop.col(pops[seq(length(pops), 1, -1)]), xcolor = get.pop.col(pops)) + 
   scale_fill_gradient2(low = 'white', high = 'grey20', mid = 'grey80', midpoint = 1, limits = c(1,2))
  ggsave(file.path(plotdir, paste(clust.pref, '_clust', k, '_dist.pdf', sep = '')), p, width = 7, height = 6)
  d.all = rbind(d.all, as.vector(tril(d)[tril(d, -1) > 0])) # Differences between all pairs of populations in a vector
}
rownames(d.all) = paste('cluster', 1:6)
cols = c()
for(i in 1:(ncol(d)-1)){
  for(j in (i+1):ncol(d)){
    cols = append(cols, paste(colnames(d)[i], colnames(d)[j]))
  }
}
for(i in 1:K){
  cat('Cluster', i, 'max diff: ', paste(cols[sort(d.all[i,], index.return=T, decreasing = T)$ix], collapse = ' , '), '\n')
}
colnames(d.all) = cols
d.all.norm = t(scale(t(d.all), center = T, scale = F))
p = plot.tile(d.all, x.ord.samples = colnames(d.all), y.ord.sample = rownames(d.all), midpoint = 1)