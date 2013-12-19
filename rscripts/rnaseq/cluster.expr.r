rm(list=ls())
library(DESeq)
library(matrixStats)
library(reshape)
library(ggplot2)
library(preprocessCore)
library(ape)
library(geiger)
source('utils/sample.info.r')
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('isva/isvaFn.R')
source('isva/DoISVA.R')
source('isva/EstDimRMT.R')

set.seed(1)

# Clustering and visualization of signal in regions (exons, peak regions etc)

#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/segSignal/14indiv/extractSignal/fc/avgSig/') 
#counts.dir = '../../rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig_newNorm/'
counts.dir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig_newNorm/'
#counts.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13'
#counts.dir = '../../rawdata/signal/combrep/extractSignal/rand/fc/avgSig/merged_Mar13/'
#counts.dir = '../../rawdata/dhs/alan/combrep/extractSignal/fc/avgSig/'
#counts.dir = '../../rawdata/geneCounts/rdata/repsComb/'
########### CHANGE THIS !!!!!!!!!!!!!!!!!!!!!
outpref = 'SNYDER_HG19_all_reg_' 
#outpref = 'gencode.v13.annotation.noM.genes_all_reg_'
#outpref = 'txStates_10_11_12_'
#outpref = 'hg19_w10k_all_reg_'
#outpref = 'all_reg_'
#outpref = 'pritchard_dhs_200bp_left_'

plotdir = file.path(counts.dir, 'plots', 'qn_isvaNull_fits_all_reg_v2')
outdir = file.path(counts.dir, 'rdata') # Set to null if you don't want to write the merged data to a file
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)
if(!file.exists(outdir)) dir.create(outdir)

k = 4
comp = 3
qval = 0.01
is.genes = F # T for RZ data
plot.only = F
quant = 0.4 # 0.4 for peak regions and transcriptomes
mark = 'H3K36ME3'
outpref = paste(outpref, mark, sep = '')
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

if(!plot.only){
  if(is.genes){
    indivs = unique(as.character(sample.info(list.files(counts.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
    nindivs = length(indivs)
    counts.dat = avg.counts(counts.dir, indivs, paste(mark, '_0.RData', sep = ''), meta = gene.meta, len.norm = F)
    regions = counts.dat$regions
    counts = asinh(counts.dat$counts) 
  }else{
    ############ Do this if you're using average signal (eg FC) in regions
    # region.file: BED file with regions to read. 
    # signal.files: should be txt files with just one column of values with the signal in each of the regions in region.file
    counts.dir = file.path(counts.dir, 'textFiles')
    region.file = file.path('../../rawdata/signal/combrep/peakFiles/merged_Mar13/', paste('SNYDER_HG19', mark, 'merged.bed.gz', sep = '_'))
    #region.file = paste('../../rawdata/signal/combrep/peakFiles/merged_Mar13/rand/SNYDER_HG19', mark, 'merged_rand.bed.gz', sep = '_')
    #region.file = '../../rawdata/genomeGrid/hg19_w10k.bed'
    region.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.bed'
    #region.file = '../../rawdata/segSignal/14indiv/txStates_10_11_12.bed'
    #region.file = '../../rawdata/dhs/alan/pritchard_dhs_200bp_left.bed'
    signal.files = list.files(counts.dir, pattern = paste(gsub('.bed|.bed.gz', '', basename(region.file)), '_AT_SNYDER_HG19_.*', mark, '.*.txt', sep = ''), full.names = T)
    indivs = unique(gsub(paste('.*_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
    sel.indivs = indivs != 'GM12890' #& get.pop(indivs) != 'San'
    indivs = indivs[sel.indivs]
    nindivs = length(indivs)
    counts.dat = load.avg.sig.data(region.file, signal.files[sel.indivs], indivs) 
    regions = counts.dat$regions
    counts = asinh(counts.dat$signal) 
    bad.ranges = regions.to.ranges(read.bed('../../rawdata/genomes_local/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed'))
    bad = countOverlaps(regions.to.ranges(regions), bad.ranges, ignore.strand = T) > 0
    counts = counts[!bad, ]
    regions = regions[!bad, ]
#     if(basename(region.file) == 'gencode.v13.annotation.noM.genes.bed'){
#       tmp = read.table(region.file, header = F, stringsAsFactors = T, sep = '\t')
#       load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
#       match.id = match(tmp[, 4], rownames(gene.meta)) # remove genes not in the annotations
#       regions = regions[!is.na(match.id), ]
#       counts = counts[!is.na(match.id), ]
#       regions$gene.name = gene.meta$gene.name[match.id[!is.na(match.id)]]
#       regions = regions[order(match.id[!is.na(match.id)]), ] # Put in the same order as the annotation
#       counts = counts[order(match.id[!is.na(match.id)]), ]
#     }
  }
  indivs = fix.indiv.names(indivs)
  
  ############# Quantile normalization
  counts = normalize.quantiles(counts)
  colnames(counts) = indivs 
  
  # Remove rows with low variance or NaNs
  good.rows = apply(counts, 1, function(x) !any(is.na(x)))
  counts = counts[good.rows, ]
  regions = regions[good.rows, ]
  row.means = rowMeans(counts)
  row.sds = rowSds(counts)
  cvs = row.sds / row.means
  ################## Increase this minimum cutoff if you're using non-enriched regions
  good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, quant, na.rm = T)
  if(!is.null(regions$chr)){
    good.rows = good.rows & !grepl('chr[XY]', regions$chr)
  }else{
    good.rows = good.rows & !grepl('chr[XY]', gene.meta$chr)
  }
  if(!is.null(outdir)) save(regions, counts, good.rows, file = file.path(outdir, paste(outpref, '_qn.RData', sep = '')))
  counts = counts[good.rows, ]
  counts.norm = t(scale(t(counts))) #apply(counts, 2, function(x) (x - row.means[good.rows])/row.sds[good.rows])
  regions = regions[good.rows, ]
  
  ############### Do PCA on the "un-corrected" data and plot eigenvalues
  counts.no.child = counts.norm[, !(indivs %in% c('GM12878', 'GM19240'))]
  pca.fit = prcomp(t(counts.no.child), center = F, scale = F)
  p.dat = plot.pcs(t(counts.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = array('', dim=c(ncol(counts.norm),1)), groups = get.pop(indivs), all = F, ndim = 2)
  ggsave(file.path(plotdir, paste(outpref, '_pca_preIsva_small.pdf', sep = '')), p.dat$p1, width = 4, height = 3)
  ggsave(file.path(plotdir, paste(outpref, '_eigen.pdf', sep = '')), p.dat$p2, width = 6.5, height = 5.6)
  
  ############### ISVA correction to remove batch effects
  pop = factor(get.pop(indivs))
  isva.fit = DoISVA(counts.norm, pop, th = qval, ncomp = comp, sel.col = 1:nindivs) # th = 0.05, ncomp = 2) for peaks and transcriptomes
  
  outpref2 = paste(outpref, '_comp', comp, '_q', qval, sep = '')
  counts = normalize.quantiles(isva.fit$res.null) # Get residuals after removing ISVs and renormalize
  colnames(counts) = indivs
  counts.norm = t(scale(t(counts)))
  if(!is.null(outdir)) save(regions, counts, isva.fit, file = file.path(outdir, paste(outpref2, '_qn_isvaNull.RData', sep = ''))) 
  
  sel = isva.fit$deg # Regions significantly correlated with population
  # Write the significant regions
  oldsc = options(scipen = 100) # prevent scientific notation in output
  outfile = file.path(plotdir, paste(outpref2, '_sign.txt', sep = ''))
  if(!is.null(regions$chr)){
    names = data.frame(chr = regions$chr, start = regions$start - 1, end = regions$end)[sel, ]
  }else{
    names = gene.meta$gene.name[match(rownames(regions)[sel], rownames(gene.meta))]
  }
  write.table(names, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")
  options(oldsc)
  
  # Plot ISVs
  isv.dat = data.frame(isva.fit$isv)
  colnames(isv.dat) = paste('ISV', 1:ncol(isv.dat))
  isv.dat$indiv = factor(indivs)
  isv.dat = melt(isv.dat, id.vars=c('indiv'))
  p = ggplot(isv.dat) + geom_bar(aes(x = indiv, y = value), position = 'dodge', stat = 'identity') + facet_wrap(~variable) +
    xlab('') + ylab('ISV value') + 
    theme(strip.text.x = element_text(size = 16), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 11, angle = -65, vjust = 1, hjust = 0))
  ggsave(file.path(plotdir, paste(outpref2, '_isvs.pdf', sep = '')), p, width = 9, height = 6.8)
}else{
  outpref2 = paste(outpref, '_comp', comp, '_q', qval, sep = '')
  load(file.path(outdir, paste(outpref2, '_qn_isvaNull.RData', sep = '')))
  indivs = colnames(counts)
  nindivs = length(indivs)
  counts.norm = t(scale(t(counts)))
}

############### FINAL PLOTTING - Can do it by reloading the SVA corrected data
no.child.cols = !(indivs %in% c('GM12878', 'GM19240'))
sel = isva.fit$deg

# pca.indiv = prcomp(counts, center = F, scale = F)
# pcs = c(1)
# pc.counts = counts
# for(i in 1:dim(counts)[1]){
#   d = data.frame(y = counts[i, ])
#   model = 'y ~ '
#   for(j in pcs){
#     d[[paste('x', j, sep = '')]] = pca.indiv$rotation[, j]
#     if(j == pcs[1]){ model = paste(model, ' ', 'x', j, sep = '')
#     }else{model = paste(model, ' + ', 'x', j, sep = '')}
#   }
#   pc.counts[i, ] = lm(model, data = d)$residuals
# }
# counts = normalize.quantiles(pc.counts)
# colnames(counts) = indivs

############### Pairwise correlations between individuals
corr.mat = array(1, dim = c(nindivs, nindivs))
for(i in 1:(nindivs - 1)){
  for(j in (i + 1):nindivs){
    corr.mat[i, j] = cor(counts.norm[, i], counts.norm[,j], use = 'pairwise.complete.obs', method = 'pearson')
    corr.mat[j, i] = corr.mat[i, j]
  }
}
plot.heatmap(corr.mat, filt.thresh = NA, symm.cluster = T, lab.row = indivs, lab.col = indivs, row.title= '', col.title = '', cex.row = 1.5, cex.col = 1.5,
             dist.metric='spearman', clust.method = "average", break.type='linear', palette = brewer.pal(9, 'Reds'), ColSideColors = get.pop.col(get.pop(indivs)),
             RowSideColors = get.pop.col(get.pop(indivs)), keysize = 1,
             to.file = file.path(plotdir, paste(outpref2, '_corrMat.pdf', sep = '')))

############### PCA
# Rows are observations (cells), columns are variables
# Remove daughters before PCA, but then project them too on the new dimensions
counts.no.child = counts.norm[, no.child.cols]
pca.fit = prcomp(t(counts.no.child), center = F, scale = F)
p=plot.pcs(t(counts.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = array('', dim=c(ncol(counts.norm),1)), groups = get.pop(indivs), all = F, ndim = 2)
ggsave(file.path(plotdir, paste(outpref2, '_pca_small.pdf', sep = '')), p$p1, width = 4, height = 3)

# Write the genes or regions with the largest loadings for enrichment analysis
# sel = order(pca.fit$rotation[, 1])[1:1000]
# outfile = file.path(plotdir, paste('pca_', mark, '.txt', sep = ''))
# if(is.null(meta)){
#   names = data.frame(chr = regions$chr[sel], start = regions$start[sel] - 1, end = regions$end[sel])
# }else{
#   names = meta$gene.name[sel]
# }
# write.table(names, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

############### Kmeans
outpref3 = paste(outpref2, '_K', k, sep = '')
if(isva.fit$ndeg > 10){
  kclusters = kmeans(counts.norm[isva.fit$deg, ], centers = k, iter.max = 1000, nstart = 10)
  kord = heatmap.2(t(scale(t(kclusters$centers))))$rowInd
  new.clusters = kclusters
  new.clusters$size = kclusters$size[kord]
  new.clusters$centers = kclusters$centers[kord, ]
  new.clusters$withinss = kclusters$withinss[kord]
  
  oldsc = options(scipen = 100)
  sel.rows = c()
  row.sep = c()
  for(i in 1:k){
    sel = which(kclusters$cluster == kord[i])
    new.clusters$cluster[sel] = i
    outfile = file.path(plotdir, paste(outpref3, '_clust', i, '.txt', sep = ''))
    if(!is.null(regions$chr)){
      names = data.frame(chr = regions$chr, start = regions$start - 1, end = regions$end)[isva.fit$deg[sel], ]
    }else{
      names = gene.meta$gene.name[match(rownames(regions)[isva.fit$deg][sel], rownames(gene.meta))]
    }
    write.table(names, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")
    cluster.frac = length(sel) / length(isva.fit$deg) # fraction of regions belonging in the cluster
    sel.rows = append(sel.rows, sel[sample(1:length(sel), min(length(sel), round(5000 * cluster.frac)))]) # Sample the cluster proportionally to its size
    row.sep = append(row.sep, length(sel.rows))
    
    nj.tree = nj(dist(t(counts[isva.fit$deg[sel], ])))
    pdf(file.path(plotdir, paste(outpref3, '_clust', i, '_dendro_nj.pdf', sep = '')))
    plot(nj.tree, 'u', cex = 1, edge.width = 0.5, no.margin = T, lab4ut='axial', label.offset = 0.5)
    dev.off()
  }
  options(oldsc)
  kclusters = new.clusters
  plot.rows = sel.rows
  
  #idx = sort(kclusters$cluster[sel.rows], index.return = T) # Sort by cluster idx
  #sel.rows = isva.fit$deg[sel.rows][idx$ix]
  #h = plot.heatmap(counts.norm[isva.fit$deg[sel.rows], ], row.cluster = F, col.cluster = T, show.dendro = "col", row.title= '', col.title = '', lab.row = NA, 
  #                 dist.metric = "euclidean", clust.method = "ward", 
  #             break.type='quantile', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], ColSideColors = get.pop.col(get.pop(indivs)), 
  #             RowSideColors = rep('white', length(sel.rows)), cex.col = 2, row.sep = row.sep, keysize = 1)
               #to.file = file.path(plotdir, paste(outpref, 'biclust_', mark, '_k', k, '.pdf', sep = '')))
  #plot.cols = h$colInd
  #pdf(file.path(plotdir, paste(outpref3, '_dendro.pdf', sep = '')))
  #plot(h$colDendrogram)
  #dev.off()
  
  nj.tree = nj(dist(t(counts[isva.fit$deg, ])))
  edges = nj.tree$edge
  edge.len = as.integer(nj.tree$edge.length * 100 / max(nj.tree$edge.length))
  sel.edges = edges[, 1] > nindivs & edges[, 2] > nindivs & edge.len > 10
  edge.lab = array('', dim = c(nrow(edges), 1))
  edge.lab[sel.edges] = edge.len[sel.edges]
  pdf(file.path(plotdir, paste(outpref2, '_dendro_nj.pdf', sep = '')))
  plot(nj.tree, 'u', cex = 1, edge.width = 0.5, no.margin = T, lab4ut='axial', label.offset = 0.5, tip.col = get.pop.col(get.pop(indivs)))
  edgelabels(edge.lab, frame = 'none', adj = c(1, 0.5), cex = 0.9)
  dev.off()
  
  leaves.tmp = node.leaves(nj.tree, nindivs + 1)
  root = which(get.pop(leaves.tmp) == 'San')[1]
  leaves.tmp = append(leaves.tmp[root:nindivs], leaves.tmp[1:(root-1)])
  lpops = unique(get.pop(leaves.tmp))
  leaves = c()
  for(i in 1:length(lpops)){
    leaves = append(leaves, leaves.tmp[get.pop(leaves.tmp) == lpops[i]])
  }
  plot.cols = match(leaves, indivs)
  
  h2 = plot.tile(counts.norm[isva.fit$deg[sel.rows], plot.cols], x.ord.samples = indivs[plot.cols], xcolor = get.pop.col(get.pop(indivs[plot.cols])), 
                 ysep = row.sep[-k] + 1, ylabels = array('', dim = c(k-1,1)), lcex = 12, xcex = 13)
  ggsave(file.path(plotdir, paste(outpref3, '_biclust.pdf', sep = '')), h2, width = 7, height = 6.8)
  save(kclusters, kord, plot.rows, plot.cols, file = file.path(outdir, paste(outpref3, '_clust.RData', sep = '')))
  
  orig = new.env()
  load(file.path(outdir, paste(outpref, '_qn.RData', sep = '')), orig)

  # First select the rows that participated in ISVA. Then, get the same rows that were used in the above clustering.
  plot.counts = orig$counts[orig$good.rows, plot.cols][isva.fit$deg[sel.rows], ]
  h3 = plot.tile(t(scale(t(plot.counts))), x.ord.samples = indivs[plot.cols], xcolor = get.pop.col(get.pop(indivs[plot.cols])), 
                 ysep = row.sep[-k] + 1, ylabels = array('', dim = c(k-1,1)), lcex = 12, xcex = 13)
  ggsave(file.path(plotdir, paste(outpref3, '_biclust_preIsva.pdf', sep = '')), h3, width = 7, height = 6.8)
  #plot.heatmap(t(scale(t(plot.counts))), row.cluster = F, col.cluster = F, show.dendro = "none", row.title= '', col.title = '', lab.row = NA, lab.col = indivs[plot.cols], dist.metric = "euclidean", clust.method = "ward", 
  #             break.type='quantile', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], ColSideColors = get.pop.col(get.pop(indivs[plot.cols])), 
  #             RowSideColors = rep('white', length(sel.rows)), cex.col = 2, row.sep = row.sep, keysize = 1)
#                to.file = file.path(plotdir, paste(outpref, 'biclust_preIsva_', mark, '_k', k, '.pdf', sep = '')))
}

########### Do this if you're using gene/peakcounts
# diff = get.diff.count(list.files(deseq.dir, pattern = paste('.*', mark, '_deseq.RData', sep = ''),  full.names = T), 1e-10) > 10
# indivs = unique(as.character(sample.info(list.files(counts.dir, pattern = paste('SNYDER_HG19_.*', mark, '_0.RData', sep = '')), '.RData')$indiv))
# counts.dat = avg.counts(counts.dir, indivs, paste(mark, '_0.RData', sep = ''), meta = meta) 
# regions = counts.dat$regions[diff, ]
# counts = counts.dat$counts[diff, ]
# meta = meta[diff, ]

############ Do this if you're using signal at dips
# dip.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/dips/llr/bed/SNYDER_HG19_H3K27AC_merged_dips.bed')
# dip.count.files = list.files(counts.dir, pattern = 'SNYDER_HG19_H3K27AC_merged_dips_AT_SNYDER_HG19_.*txt', full.names = T)
# indivs = unique(as.character(sample.info(gsub('SNYDER_HG19_H3K27AC_merged_dips_AT_', '', dip.count.files), '.txt')$indiv))
# nindivs = length(indivs)
# counts.dat = load.dip.data(dip.file, dip.count.files, indivs)
# regions = counts.dat$regions
# counts = (counts.dat$left.sig + counts.dat$right.sig) / 2
# sel = apply(counts, 1, function(x) !any(is.na(x)))
# regions = regions[sel, ]
# counts = counts[sel, ]

############ MAnorm correction
# Bed files with individual specific peaks. Fits between different individuals are performed on the common peaks only.
# bed.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles')
# bed.files = array('', dim = c(nindivs, 1))
# missing = array(F, dim = c(nindivs, 1))
# for(i in 1:nindivs){
#   tmp.files = list.files(bed.dir, pattern = paste('SNYDER_HG19', indivs[i], mark, 'VS_.*', sep = '_'), full.names = T)
#   if(length(tmp.files) != 1){
#     #missing[i] = T
#     tmp.files = list.files(bed.dir, pattern = paste('SNYDER_HG19', indivs[i], mark, '.*dedup_VS_.*', sep = '_'), full.names = T)
#   }
#   stopifnot(length(tmp.files) == 1)
#   bed.files[i] = tmp.files
# }
# 
# bed.files = bed.files[!missing]
# indivs = indivs[!missing]
# nindivs = length(indivs)
# counts = counts[, !missing]
# fits = fit.norm.counts(counts, bed.files, ref.idx = which(indivs == 'GM12878'), take.log = F)
# counts = fits$lr + rep(fits$avg[, which(indivs == 'GM12878')], nindivs)
