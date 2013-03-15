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

# Compares population-specific regions for two different marks

isva.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/') 
########### CHANGE THIS !!!
outpref = 'SNYDER_HG19_' #'hg19_w10k_' #gencode.v13.annotation.noM.genes_' #
suf = '_qn_svaNull'
plotdir = file.path(isva.dir, '..', 'plots', 'qn_svaNull_fits')
is.sva = !grepl('isva', suf) # ISVA or SVA? Use the suffix to figure out

if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)
k = 7
in.mark = 'H3K4ME1'
pair.mark = 'H3K27AC'

# Reference clusters
in.env = new.env()
load(file.path(isva.dir, paste(outpref, in.mark, suf, '.RData', sep = '')), in.env)
indivs = colnames(in.env$counts)
load(file.path(isva.dir, paste(outpref, in.mark, suf, '_clust.RData', sep = '')))
in.kord = kord # make sure this is not overwritten
if(is.sva){
  in.deg = in.env$sva.fit$deg
}else{
  in.deg = in.env$isva.fit$deg
}
load(file.path(isva.dir, paste(outpref, pair.mark, suf, '.RData', sep = '')))
if(is.sva){
  pair.deg = sva.fit$deg
}else{
  pair.deg = isva.fit$deg
}
stopifnot(all(colnames(counts) == indivs)) # make sure columns are in the same order

# Compute overlap statistics
match.idx = findOverlaps(regions.to.ranges(in.env$regions[in.deg, ]), regions.to.ranges(regions[pair.deg, ]), select = 'first', ignore.strand = T)
cat('Fraction of pop-specific regions of', in.mark, 'that are pop-specific for', pair.mark, sum(!is.na(match.idx)) / length(in.deg), '\n')

# Take reference regions that are differentially expressed. Then take the subset of them that was used for plotting.
# Finally, for each of these regions, find an overlapping region in the second set. This will be used in the new plot.
match.idx = findOverlaps(regions.to.ranges(in.env$regions[in.deg[plot.rows], ]), regions.to.ranges(regions), select = 'first', ignore.strand = T)
row.sep = which(!duplicated(in.env$kclusters$cluster[plot.rows[!is.na(match.idx)]])) # Replace the cluster separators

plot.counts = scale(counts[match.idx[!is.na(match.idx)], plot.cols])
plot.heatmap(plot.counts, row.cluster = F, col.cluster = F, show.dendro = "none", row.title= '', col.title = '', lab.row = NA, 
             dist.metric = "euclidean", clust.method = "ward", 
             break.type='quantile', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], 
             ColSideColors = get.pop.col(get.pop(indivs[plot.cols])), 
             RowSideColors = rep('white', nrow(plot.counts)), cex.col = 2, row.sep = row.sep,
             to.file = file.path(plotdir, paste('biclust_', in.mark, '_vs_', pair.mark, '.png', sep = '')))