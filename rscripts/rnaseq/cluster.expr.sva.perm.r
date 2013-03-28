rm(list=ls())
library(DESeq)
library(matrixStats)
library(reshape)
library(ggplot2)
library(gplots)
library(preprocessCore)
library(sva)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/isvaFn.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/DoISVA.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/EstDimRMT.R'))

set.seed(1)

# Clustering, visualization, and SVA correction of signal in regions (exons, peak regions etc).
# Similar to cluster.expr.r, except that is uses SVA instead of ISVA.

run.sva = function(data, pheno, cf = NULL, n.sv = 1, th = 0.05){
  mod = model.matrix(~ pheno)
  if(is.null(n.sv)){
    n.sv = num.sv(data,mod,method="leek")
  }
  
  if(is.null(cf)){
    mod0 = model.matrix(~ 1, data = data.frame(t(data)))
  }else{
    mod0 = model.matrix(~ cf)
  }
  svobj = sva(data, mod, mod0, n.sv=n.sv, method = "two-step")
  
  modSv = cbind(mod,svobj$sv)
  mod0Sv = cbind(mod0,svobj$sv)
  
  n = dim(data)[2]
  m = dim(data)[1]
  df1 = dim(mod)[2]
  df0 = dim(mod0)[2]
  p = rep(0, m)
  Id = diag(n)
  
  resid <- data %*% (Id - modSv %*% solve(t(modSv) %*% modSv) %*% t(modSv))
  rss1 <- rowSums(resid*resid)  
  resid0 <- data %*% (Id - mod0Sv %*% solve(t(mod0Sv) %*% mod0Sv) %*% t(mod0Sv))
  rss0 <- rowSums(resid0*resid0)
  
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  
  q = p.adjust(p, method="BH")
  cat('Number of significant regions', sum(q < th), '\n')
  # Return: 
  # ndeg, deg: number and indices of selected regions (names compatible with ISVA)
  # res.null: "Corrected" signal, after removing the selected SVs
  # sv: SVs
  # sv.mem: which regions participated in SVs
  return(list(ndeg = sum(q < th), deg = which(q < th), res.null = resid0, sv = svobj$sv, sv.mem = svobj$pprob.gam))
}

#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/segSignal/14indiv/extractSignal/fc/avgSig/') 
#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig/') 
#counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/extractSignal/fc/avgSig/')
counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/')
#counts.dir= file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/repsComb/')
#deseq.dir = file.path(counts.dir, 'deseq/')
########### CHANGE THIS !!!!!!!!!!!!!!!!!!!!!
outpref = 'SNYDER_HG19_' 
#outpref = 'gencode.v13.annotation.noM.genes_all_reg_'
#outpref = 'txStates_10_11_12_'
#outpref = 'hg19_w10k_'
#outpref = 'all_reg_'

plotdir = file.path(counts.dir, 'plots', 'qn_svaNull_fits')
outdir = file.path(counts.dir, 'rdata')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)
if(!file.exists(outdir)) dir.create(outdir)
counts.dir = file.path(counts.dir, 'textFiles')
k = 5
nperm = 10
is.genes = F # T for RZ data
quant = 0.4 # 0.4 for peak regions and transcriptomes
mark = 'H3K27AC'

if(is.genes){
  load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
  indivs = unique(as.character(sample.info(list.files(counts.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
  nindivs = length(indivs)
  counts.dat = avg.counts(counts.dir, indivs, paste(mark, '_0.RData', sep = ''), meta = gene.meta, len.norm = F)
  regions = counts.dat$regions
  regions$gene.name = gene.meta$gene.name
  counts = asinh(counts.dat$counts) 
}else{
  ############ Do this if you're using average signal (eg FC) in regions
  # region.file: BED file with regions to read. 
  # signal.files: should be txt files with just one column of values with the signal in each of the regions in region.file
  
  region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/', paste('SNYDER_HG19', mark, 'merged.bed.gz', sep = '_'))
  #egion.file = paste('../../rawdata/signal/combrep/peakFiles/merged/rand/SNYDER_HG19', mark, 'merged_rand.bed.gz', sep = '_')
  #region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomeGrid/hg19_w10k.bed')
  #region.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.bed'
  #region.file = '../../rawdata/segSignal/14indiv/txStates_10_11_12.bed'
  signal.files = list.files(counts.dir, pattern = paste(gsub('.bed|.bed.gz', '', basename(region.file)), '_AT_SNYDER_HG19_.*', mark, '.*.txt', sep = ''), full.names = T)
  indivs = unique(gsub(paste('.*_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
  nindivs = length(indivs)
  counts.dat = load.avg.sig.data(region.file, signal.files, indivs) 
  regions = counts.dat$regions
  counts = asinh(counts.dat$signal) 
  if(basename(region.file) == 'gencode.v13.annotation.noM.genes.bed'){
    tmp = read.table(region.file, header = F, stringsAsFactors = T, sep = '\t')
    load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
    match.id = match(tmp[, 4], rownames(gene.meta)) # remove genes not in the annotations
    regions = regions[!is.na(match.id), ]
    counts = counts[!is.na(match.id), ]
    regions$gene.name = gene.meta$gene.name[match.id[!is.na(match.id)]]
    regions = regions[order(match.id[!is.na(match.id)]), ] # Put in the same order as the annotation
    counts = counts[order(match.id[!is.na(match.id)]), ]
  }
}

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
good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, quant, na.rm = T)
#if(!is.null(outdir)) save(regions, counts, good.rows, file = file.path(outdir, paste(outpref, mark, '_qn.RData', sep = '')))
orig.counts = counts[good.rows, ]
regions = regions[good.rows, ]

############### SVA correction to remove batch effects
pop = factor(get.pop(indivs))

for(n in 1:nperm){
  sva.fit = run.sva(scale(orig.counts[, sample(1:nindivs, nindivs)]), pop, th = 0.05, n.sv = 2)
  
  counts = normalize.quantiles(sva.fit$res.null) # Get residuals after removing SVs and renormalize
  colnames(counts) = indivs
  counts.norm = scale(counts)
  if(!is.null(outdir)) save(regions, counts, sva.fit, file = file.path(outdir, paste(outpref, 'rand_pop', n, '_', mark, '_qn_svaNull.RData', sep = ''))) 
  
  sel = sva.fit$deg # Regions significantly correlated with population
  
  # Kmeans clustering
  if(sva.fit$ndeg > 50){
    kclusters = kmeans(counts.norm[sva.fit$deg, ], centers = k, iter.max = 500, nstart = 10)
    kord = heatmap.2(scale(kclusters$centers))$rowInd
    sel.rows = c()
    row.sep = c()
    for(i in 1:k){
      sel = which(kclusters$cluster == kord[i])
      cluster.frac = length(sel) / length(sva.fit$deg) # fraction of regions belonging in the cluster
      sel.rows = append(sel.rows, sel[sample(1:length(sel), min(length(sel), round(5000 * cluster.frac)))]) # Sample the cluster proportionally to its size
      row.sep = append(row.sep, length(sel.rows))
    }
    
    # Plot the clusters
    h = plot.heatmap(counts.norm[sva.fit$deg[sel.rows], ], row.cluster = F, col.cluster = T, show.dendro = "none", row.title= '', col.title = '', lab.row = NA, dist.metric = "euclidean", clust.method = "ward", 
                     break.type='quantile', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], ColSideColors = get.pop.col(get.pop(indivs)), 
                     RowSideColors = rep('white', length(sel.rows)), cex.col = 2, row.sep = row.sep,
                     to.file = file.path(plotdir, paste(outpref, 'rand_pop', n, '_', 'biclust_', mark, '.png', sep = '')))
  }
}
