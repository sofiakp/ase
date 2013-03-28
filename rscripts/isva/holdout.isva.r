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
  ncomp = ncol(old.fit$isv)
  isv.cor = array(0, dim = c(ncomp, ncomp))
  for(i in 1:ncomp){
    for(j in 1:ncomp){
      isv.cor[i, j] = cor(new.fit$isv[, i], old.fit$isv[, j])
    }
  }
  
  k = length(old.clust$size)
  cluster.ov = array(0, dim = c(k, k))
  centroid.cor = array(0, dim = c(k, k))
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
  return(list(isv.cor = isv.cor, cluster.ov = cluster.ov, centroid.cor = centroid.cor))
}

plot.stats = function(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, row.names, outpref = NULL){
  rownames(isv.cor.stats) = row.names
  colnames(isv.cor.stats) = paste('ISV', 1:ncomp)
  p1 = plot.tile(isv.cor.stats, y.ord.samples = rownames(isv.cor.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  rownames(centroid.cor.stats) = row.names
  colnames(centroid.cor.stats) = paste('cluster', c(1:k))
  p2 = plot.tile(centroid.cor.stats, y.ord.samples = rownames(centroid.cor.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  rownames(cluster.ov.stats) = row.names
  colnames(cluster.ov.stats) = paste('cluster', c(1:k))
  p3 = plot.tile(cluster.ov.stats, y.ord.samples = rownames(cluster.ov.stats), midpoint = 0.5, draw.y.line = F , xcex = 12, ycex = 12, lcex = 12)
  
  if(!is.null(outpref)){
    ggsave(paste(outpref, 'isvs.pdf'), p1, width = 9, height = 8)
    ggsave(paste(outpref, 'centroidCor.pdf'), p2, width = 9, height = 8)
    ggsave(paste(outpref, 'clusterOv.pdf'), p3, width = 9, height = 8)
  }
}

set.seed(1)

indir = '../../rawdata/geneCounts/rdata/repsComb/'
plot.dir = file.path(indir, 'plots', 'qn_isvaNull_fits_all_reg_v2')
mark = 'RZ'
pref = 'gencode.v13.annotation.noM.genes_all_reg'
original.data = new.env()
load(file.path(indir, 'rdata', paste(pref, mark, 'qn.RData', sep = '_')), original.data)
original.fit = new.env()
load(file.path(indir, 'rdata', paste(pref, mark, 'qn_isvaNull.RData', sep = '_')), original.fit)
original.clust = new.env()
load(file.path(indir, 'rdata', paste(pref, mark, 'qn_isvaNull_clust6.RData', sep = '_')), original.clust)
run.repeats = T

indivs = colnames(original.data$counts)
nindiv = length(indivs)
counts = original.data$counts[original.data$good.rows, ]
counts.norm = scale(counts)
regions = original.data$regions[original.data$good.rows, ]
pop = factor(get.pop(indivs))
ncomp = ncol(original.fit$isva.fit$isv)
k = max(original.clust$kclusters$cluster)

if(run.repeats){
  cluster.ov.stats = c()
  centroid.cor.stats = c()
  isv.cor.stats = c()
  nrepeats = 50
  for(r in 1:nrepeats){
    isva.fit = DoISVA(counts.norm, pop, th = 0.05, ncomp = ncomp, sel.col = 1:nindiv)
    
    new.counts = normalize.quantiles(isva.fit$res.null)
    colnames(new.counts) = colnames(counts)
    kclusters = kmeans(scale(new.counts[isva.fit$deg, ]), centers = k, iter.max = 1000, nstart = 10)
    
    cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
    isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
    cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
    centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
  }
  p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, paste('repeat', 1:nrepeats), file.path(plot.dir, paste(pref, 'repeatStats', mark, sep = '_')))
}

cluster.ov.stats = c()
centroid.cor.stats = c()
isv.cor.stats = c()
for(sel.out in 1:nindiv){
  sel.col = setdiff(1:nindiv, sel.out)
  isva.fit = DoISVA(counts.norm, pop, th = 0.05, ncomp = ncomp, sel.col = sel.col)
 
  new.counts = normalize.quantiles(isva.fit$res.null)
  colnames(new.counts) = colnames(counts)
  kclusters = kmeans(scale(new.counts[isva.fit$deg, ]), centers = k, iter.max = 1000, nstart = 10)
  
  cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
  isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
  cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
  centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
}

p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, indivs, file.path(plot.dir, paste(pref, 'hold1Out', mark, sep = '_')))

cluster.ov.stats = c()
centroid.cor.stats = c()
isv.cor.stats = c()
out.cols = c()
all.pops = unique(get.pop(indivs))
while(is.null(nrow(cluster.ov.stats)) || nrow(cluster.ov.stats) < 50){
  sel.out = sample(1:nindiv, 2)
  new.name = paste(indivs[sel.out], collapse = '_')
  if(length(unique(get.pop(indivs[-sel.out]))) < length(all.pops) || new.name %in% out.cols) next
  out.cols = append(out.cols, new.name)
  sel.col = setdiff(1:nindiv, sel.out)
  isva.fit = DoISVA(counts.norm, pop, th = 0.05, ncomp = ncomp, sel.col = sel.col)
 
  new.counts = normalize.quantiles(isva.fit$res.null)
  colnames(new.counts) = colnames(counts)
  kclusters = kmeans(scale(new.counts[isva.fit$deg, ]), centers = k, iter.max = 1000, nstart = 10)
  cstats = get.clust.ov(original.fit$isva.fit, isva.fit, original.clust$kclusters, kclusters)
  isv.cor.stats = rbind(isv.cor.stats, rowMaxs(abs(cstats$isv.cor)))
  cluster.ov.stats = rbind(cluster.ov.stats, rowMaxs(cstats$cluster.ov))
  centroid.cor.stats = rbind(centroid.cor.stats, rowMaxs(cstats$centroid.cor))
}

p = plot.stats(isv.cor.stats, centroid.cor.stats, cluster.ov.stats, out.cols, file.path(plot.dir, paste(pref, 'hold2Out', mark, sep = '_')))
