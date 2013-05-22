rm(list=ls())
library(DESeq)
library(matrixStats)
library(reshape)
library(ggplot2)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/isvaFn.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/DoISVA.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/EstDimRMT.R'))

set.seed(1)

#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/segSignal/14indiv/extractSignal/fc/avgSig/') 
#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig/') 
#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/extractSignal/fc/avgSig/')
counts.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/'
rdir = file.path(counts.dir, 'rdata')
#counts.dir= file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/repsComb/')
#deseq.dir = file.path(counts.dir, 'deseq/')

plotdir = file.path(counts.dir, 'plots', 'qn_isvaNull_fits_all_reg_v2', 'perm')
outdir = file.path(counts.dir, 'rdata', 'perm')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)
if(!file.exists(outdir)) dir.create(outdir)

nperm = 10
mark = 'H3K4ME1'
pref = 'SNYDER_HG19_all_reg_'
comp = 3
qval = 0.01
K = 4
nrand = 100
pref = paste(pref, mark, sep = '')

load(file.path(rdir, paste(pref, '_qn.RData', sep = '')))
pref = paste(pref, '_comp', comp, '_q', qval, sep = '')
counts = counts[good.rows, ]
counts.orig = t(scale(t(counts)))
regions = regions[good.rows, ]

indivs = colnames(counts)
nindivs = length(indivs)

pop = factor(get.pop(indivs))

for(n in 1:nperm){
  sam = sample(1:nindivs, nindivs)
  print(pop[sam])
  isva.fit = DoISVA(counts.orig[, sam], pop, th = qval, ncomp = comp, sel.col = 1:nindivs)
  counts = normalize.quantiles(isva.fit$res.null) # Get residuals after removing ISVs and renormalize
  colnames(counts) = indivs
  counts.norm = t(scale(t(counts)))
  colnames(counts.norm) = indivs[sam]
  indivs.tmp = indivs[sam]
  perm.pref = paste(pref, '_rand_pop', n, sep = '')
  if(!is.null(outdir)) save(regions, counts, isva.fit, sam, file = file.path(outdir, paste(perm.pref, '_qn_isvaNull.RData', sep = ''))) 
  
  perm.pref = paste(perm.pref, '_K', K, sep = '')
  if(isva.fit$ndeg > 10){
    kclusters = kmeans(counts.norm[isva.fit$deg, ], centers = K, iter.max = 1000, nstart = 10)
    kord = heatmap.2(t(scale(t(kclusters$centers))))$rowInd
    new.clusters = kclusters
    new.clusters$size = kclusters$size[kord]
    new.clusters$centers = kclusters$centers[kord, ]
    new.clusters$withinss = kclusters$withinss[kord]
    
    sel.rows = c()
    row.sep = c()
    for(i in 1:K){
      sel = which(kclusters$cluster == kord[i])
      new.clusters$cluster[sel] = i
      cluster.frac = length(sel) / length(isva.fit$deg) # fraction of regions belonging in the cluster
      sel.rows = append(sel.rows, sel[sample(1:length(sel), min(length(sel), round(5000 * cluster.frac)))]) # Sample the cluster proportionally to its size
      row.sep = append(row.sep, length(sel.rows))
    }
    kclusters = new.clusters
    plot.rows = sel.rows
    
    nj.tree = nj(dist(t(counts.norm[isva.fit$deg, ])))
    leaves.tmp = node.leaves(nj.tree, nindivs + 1)
    root = which(get.pop(leaves.tmp) == 'San')[1]
    leaves.tmp = append(leaves.tmp[root:nindivs], leaves.tmp[1:(root-1)])
    lpops = unique(get.pop(leaves.tmp))
    leaves = c()
    for(i in 1:length(lpops)){
      leaves = append(leaves, leaves.tmp[get.pop(leaves.tmp) == lpops[i]])
    }
    plot.cols = match(leaves, indivs)
    
    h2 = plot.tile(counts.norm[isva.fit$deg[sel.rows], plot.cols], x.ord.samples = indivs.tmp[plot.cols], xcolor = get.pop.col(get.pop(indivs.tmp[plot.cols])), 
                   ysep = row.sep[-K] + 1, ylabels = array('', dim = c(K-1,1)), lcex = 12, xcex = 13)
    ggsave(file.path(plotdir, paste(perm.pref, '_biclust.pdf', sep = '')), h2, width = 7, height = 6.8)
    save(kclusters, kord, plot.rows, plot.cols, file = file.path(outdir, paste(perm.pref, '_clust.RData', sep = '')))
  }
}