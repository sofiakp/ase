rm(list=ls())
library(ggplot2)
library(Matrix)
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
      pop.d[p1, p2] = pop.d[p1, p2] + d[i, j]
    }
  }
  pop.d = pop.d + t(pop.d)
  colnames(pop.d) = uniq.pop
  rownames(pop.d) = uniq.pop
  norm.pop.d = pop.d
  for(i in 1:npop){
    for(j in 1:npop){
      norm.pop.d[i,j] = pop.d[i,j] / (pop.d[i,i] + pop.d[j,j])  
    }
  } 
  return(list(pop.d = pop.d, norm.pop.d = norm.pop.d))
}

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
plotdir = indir
mark = 'H3K27AC'

load('../../rawdata/variants/all_Mar13/genot.RData')
load('../../rawdata/variants/all_Mar13/snps.RData')

load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/SNYDER_HG19_all_reg_H3K4ME1_qn_isvaNull.RData')
regions = regions[isva.fit$deg, ]
load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/SNYDER_HG19_all_reg_H3K4ME1_qn_isvaNull_clust6.RData')

load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')
#ov = findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(gene.meta), select = 'first', ignore.strand = T) # Restrict to exons
ov = findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(regions), select = 'first', ignore.strand = T) 
snp.pos = snp.pos[!is.na(ov), ]
genot = genot[!is.na(ov), ]
genot = genot[, colnames(genot) != 'GM19193']
genot.indivs = fix.indiv.names(colnames(genot))

#genot.norm = scale(genot[sample(1:nrow(genot), 1200), ])
#pca.fit = prcomp(t(genot.norm[, !(colnames(genot) %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
#p = plot.pcs(t(genot.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = colnames(genot), groups = get.pop(colnames(genot)), all = F, ndim = 4)
#pca.plot.all = p$p1
#bg.dist = compute.dist(t(genot.norm) %*% pca.fit$rotation[, 1:3])$norm.pop.d

pca.plots = list()
d.all = NULL
for(k in 1:6){
  sel = kclusters$cluster == k
  ov = findOverlaps(regions.to.ranges(regions[sel, ]), snps.to.ranges(snp.pos),  select = 'all', ignore.strand = T)
  genot.norm = scale(genot[subjectHits(ov), ])
  pca.fit = prcomp(t(genot.norm[, !(genot.indivs %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
  p = plot.pcs(t(genot.norm) %*% pca.fit$rotation, pca.fit$rotation, pca.fit$sdev, labels = genot.indivs, groups = get.pop(colnames(genot)), all = F, ndim = 4)
  pca.plots[[k]] = p$p1
  ggsave(file.path(indir, paste('test_H3K27AC_k', 6, '_clust', k, '_genotPca.pdf', sep = '')), p$p1, width = 9, height = 6.8)
  d = compute.dist(t(genot.norm) %*% pca.fit$rotation[, 1:2])$pop.d
  
  avg.d = array(0, dim = dim(d))
  for(i in 1:100){
    genot.norm.bg = scale(genot[sample(1:nrow(genot), nrow(genot.norm)), ])
    pca.fit.bg = prcomp(t(genot.norm.bg[, !(colnames(genot) %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
    avg.d = avg.d + compute.dist(t(genot.norm.bg) %*% pca.fit.bg$rotation[, 1:2])$pop.d
  }
  d = d / (avg.d / 100)
  print(d)
  d.all = rbind(d.all, as.vector(tril(d)[tril(d, -1) > 0]))
}
rownames(d.all) = paste('cluster', 1:6)
cols = c()
for(i in 1:(ncol(d)-1)){
  for(j in (i+1):ncol(d)){
    cols = append(cols, paste(colnames(d)[i], colnames(d)[j]))
  }
}
colnames(d.all) = cols
d.all.norm = scale(d.all, center = T, scale = F)
p = plot.tile(d.all, x.ord.samples = colnames(d.all), y.ord.sample = rownames(d.all), midpoint = 1)