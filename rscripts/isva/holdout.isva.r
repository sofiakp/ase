require(fastICA)
library(reshape)
library(ggplot2)
library(qvalue)
library(matrixStats)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/isvaFn.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/DoISVA.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/EstDimRMT.R'))

get.clust.ov = function(old.fit, new.fit, old.clust, new.clust){
  # Compute correlation between all pairs of ISVs
  ncomp = ncol(old.fit$isv)
  isv.cor = array(0, dim = c(ncomp, ncomp))
  for(i in 1:ncomp){
    for(j in 1:ncomp){
      isv.cor[i, j] = cor(new.fit$isv[, i], old.fit$isv[, j])
    }
  }
  
  k = length(old.clust$size)
  cluster.ov = array(0, dim = c(k, k)) # Overlap of clusters (fraction of minimum cluster that is in the intersection)
  centroid.cor = array(0, dim = c(k, k)) # Correlation between old and new centroids
  for(i in 1:k){
    for(j in i:k){
      cluster.ov[i, j] = length(intersect(old.fit$deg[old.clust$cluster == i], 
                                          new.fit$deg[new.clust$cluster == j])) / min(sum(old.clust$cluster == i), 
                                                                                       sum(new.clust$cluster == j))
      cluster.ov[j, i] = cluster.ov[i, j]
      centroid.cor[i, j] = cor(old.clust$centers[i, ], new.clust$centers[j, ])
      centroid.cor[j, i] = centroid.cor[i, j]
    }
  }
  abs.ov = length(intersect(old.fit$deg, new.fit$deg))
  ov = c(abs.ov / old.fit$ndeg, abs.ov / new.fit$ndeg, abs.ov / length(union(old.fit$deg, new.fit$deg)))
  return(list(isv.cor = isv.cor, cluster.ov = cluster.ov, centroid.cor = centroid.cor, ov = ov))
}

plot.stats = function(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats, row.names, outpref = NULL){
  rownames(isv.cor.stats) = row.names
  colnames(isv.cor.stats) = paste('ISV', 1:ncol(isv.cor.stats))
  p1 = plot.tile(isv.cor.stats, midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  rownames(centroid.cor.stats) = row.names
  colnames(centroid.cor.stats) = paste('cluster', c(1:k))
  p2 = plot.tile(centroid.cor.stats, y.ord.samples = rownames(centroid.cor.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  rownames(cluster.ov.stats) = row.names
  colnames(cluster.ov.stats) = paste('cluster', c(1:k))
  p3 = plot.tile(cluster.ov.stats, y.ord.samples = rownames(cluster.ov.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  rownames(ov.stats) = row.names
  colnames(ov.stats) = c('fraction old', 'fraction new', 'Jaccard')
  p4 = plot.tile(ov.stats, y.ord.samples = rownames(ov.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  
  if(!is.null(outpref)){
    ggsave(paste(outpref, 'isvs.pdf', sep = '_'), p1, width = 9, height = 8)
    ggsave(paste(outpref, 'centroidCor.pdf', sep = '_'), p2, width = 9, height = 8)
    ggsave(paste(outpref, 'clusterOv.pdf', sep = '_'), p3, width = 9, height = 8)
    ggsave(paste(outpref, 'ov.pdf', sep = '_'), p4, width = 9, height = 8)
  }
}

set.seed(2)
k = 4
comp = 3
qval = 0.01
mark = 'H3K27AC'
indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/' #'../../rawdata/geneCounts/rdata/repsComb/'
plot.dir = file.path(indir, 'plots', 'qn_isvaNull_fits_all_reg_v2', 'stability')
if(!file.exists(plot.dir)) dir.create(plot.dir)

pref = 'SNYDER_HG19_all_reg_' #'gencode.v13.annotation.noM.genes_all_reg_'
pref = paste(pref, mark, sep = '')
isv.pref = paste(pref, '_comp', comp, '_q', qval, sep = '')
clust.pref = paste(isv.pref, '_K', k, sep = '')

original.data = new.env()
load(file.path(indir, 'rdata', paste(pref, '_qn.RData', sep = '')), original.data)
original.fit = new.env()
load(file.path(indir, 'rdata', paste(isv.pref, '_qn_isvaNull.RData', sep = '')), original.fit)
original.clust = new.env()
load(file.path(indir, 'rdata', paste(clust.pref, '_clust.RData', sep = '')), original.clust)

run.repeats = T
run.single = T
run.pairs = T

indivs = colnames(original.data$counts)
nindiv = length(indivs)
counts = original.data$counts[original.data$good.rows, ]
counts.norm = t(scale(t(counts)))
regions = original.data$regions[original.data$good.rows, ]
pop = factor(get.pop(indivs))
stopifnot(ncol(original.fit$isva.fit$isv) == comp)

if(run.repeats){
  cluster.ov.stats = c()
  centroid.cor.stats = c()
  isv.cor.stats = c()
  ov.stats = c()
  nrepeats = 10
  for(r in 1:nrepeats){
    isva.fit = DoISVA(counts.norm, pop, th = qval, ncomp = comp, sel.col = 1:nindiv)
    
    new.counts = t(scale(t(normalize.quantiles(isva.fit$res.null))))
    colnames(new.counts) = colnames(counts)
    kclusters = kmeans(new.counts[isva.fit$deg, ], centers = k, iter.max = 1000, nstart = 10)
    
    cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
    # For each of the original ICs, get the maximum correlation with a new IC
    isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
    # For each of the original clusters, get the maximum overlap with a new cluster
    cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
    # For each of the original cluster centroids, get the maximum correlation with a new cluster centroid
    centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
    # Fraction of overlap between population specific regions overall
    ov.stats = rbind(ov.stats, cstats$ov)
  }
  p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats,
                 paste('repeat', 1:nrepeats), file.path(plot.dir, paste(clust.pref, 'repeatStats', sep = '_')))
  save(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats, file = file.path(plot.dir, paste(clust.pref, 'repeatStats.RData', sep = '_')))
}

if(run.single){
  cluster.ov.stats = c()
  centroid.cor.stats = c()
  isv.cor.stats = c()
  ov.stats = c()
  for(sel.out in 1:nindiv){
    sel.col = setdiff(1:nindiv, sel.out)
    isva.fit = DoISVA(counts.norm, pop, th = qval, ncomp = comp, sel.col = sel.col)
    
    new.counts = t(scale(t(normalize.quantiles(isva.fit$res.null))))
    colnames(new.counts) = colnames(counts)
    kclusters = kmeans(scale(new.counts[isva.fit$deg, ]), centers = k, iter.max = 1000, nstart = 10)
    
    cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
    isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
    cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
    centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
    ov.stats = rbind(ov.stats, cstats$ov)
  }
  
  p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats,
                 indivs, file.path(plot.dir, paste(clust.pref, 'hold1Out', sep = '_'))) 
  save(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats, file = file.path(plot.dir, paste(clust.pref, 'hold1Out.RData', sep = '_')))
}

if(run.pairs){
  cluster.ov.stats = c()
  centroid.cor.stats = c()
  isv.cor.stats = c()
  ov.stats = c()
  out.cols = c()
  all.pops = unique(get.pop(indivs))
  while(is.null(nrow(cluster.ov.stats)) || nrow(cluster.ov.stats) < 50){
    sel.out = sample(1:nindiv, 2)
    new.name = paste(indivs[sel.out], collapse = '_')
    if(length(unique(get.pop(indivs[-sel.out]))) < length(all.pops) || new.name %in% out.cols) next
    out.cols = append(out.cols, new.name)
    sel.col = setdiff(1:nindiv, sel.out)
    isva.fit = DoISVA(counts.norm, pop, th = qval, ncomp = comp, sel.col = sel.col)
    
    new.counts = t(scale(t(normalize.quantiles(isva.fit$res.null))))
    colnames(new.counts) = colnames(counts)
    kclusters = kmeans(scale(new.counts[isva.fit$deg, ]), centers = k, iter.max = 1000, nstart = 10)
    cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
    isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
    cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
    centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
    ov.stats = rbind(ov.stats, cstats$ov)
  }
  
  p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats,
                 out.cols, file.path(plot.dir, paste(clust.pref, 'hold2Out', sep = '_'))) 
  save(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, ov.stats, file = file.path(plot.dir, paste(clust.pref, 'hold2Out.RData', sep = '_')))
}
