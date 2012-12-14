rm(list=ls())
library('DESeq')
library(matrixStats)
library(fastICA)
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')

# Tests different ways of clustering variable regions

compute.corr = function(mat, nsamples, cols, outfile, corr = T){
  cols = as.ordered(cols)
  if(corr){
    corr.mat = array(1, dim = c(nsamples, nsamples))
    for(i in 1:(nsamples - 1)){
      for(j in (i + 1):nsamples){
        corr.mat[i, j] = cor(mat[, i], mat[, j])
        corr.mat[j, i] = corr.mat[i, j]
      }
    }
    header = 'Correlation'
  }else{
    corr.mat = as.matrix(dist(t(mat)))
    header = 'Distance'
  }
  corr.mat = data.frame(corr.mat, row.names = cols)
  colnames(corr.mat) = cols
  corr.mat$sample = rownames(corr.mat)
  corr.dat = melt(corr.mat, id.vars = c('sample'))
  #corr.dat = corr.dat[order(corr.dat$variable, corr.dat$sample), ]
  p = ggplot(corr.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
    scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient(header, low = "beige", high = "tomato") +
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey20', size = 10),
          axis.text.y = element_text(hjust = 1, colour = 'grey20', size = 10)) + 
            ggtitle(paste('Correlation across differential regions for', mark))
  ggsave(outfile)
  #return(corr.dat)
}

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/countsAtPeaks/repsComb')
plotdir = file.path(indir, 'plots_deseq_v3')
if(!file.exists(plotdir)) dir.create(plotdir)

mark = 'H3K4ME3'
# Filenames to cluster
filenames <- list.files(indir, pattern = paste('SNYDER_HG19_.*', mark, '_0.RData', sep = ''), full.name = T, recursive = F, include.dirs = F)
# Cluster on differential regions from DESeq
deseq.filenames = list.files(file.path(indir, 'deseq_v3'), pattern = paste(mark, '_deseq.RData', sep = ''), full.name = T, recursive = F, include.dirs = F)
load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

qcut = 1e-10
for(i in 1:length(deseq.filenames)){
  load(deseq.filenames[i])
  diff.count = !is.na(regions$qval) & regions$qval < qcut
}

samples = sample.info(filenames, '.RData')
nfiles = length(filenames)
counts.all = NULL

for(i in 1:nfiles){
  print(basename(filenames[i]))
  load(filenames[i])
  if(dim(counts)[2] > 1){counts = t(apply(counts, 1, function(x) x / size.factors))
  }else{counts = counts / size.factors}
  
  if(is.null(counts.all)){
    counts.all = counts
  }else{
    counts.all = cbind(counts.all, counts)
  }
}

widths = regions$end - regions$start # gene.meta$len #
if(length(deseq.filenames) < 1){
  sel = rowSums(counts.all) > quantile(rowSums(counts.all), 0.5) & !grepl('chr[XY]', regions$chr)
  cvs = rowSds(counts.all) / rowMeans(counts.all)
  sel = sel & cvs > quantile(cvs, 0.8)
}else{sel = diff.count}
counts.all = apply(counts.all[sel, ], 2, function(x) x / widths[sel])
nsamples = dim(counts.all)[2]

########################### Use the variance stabilized counts
vst.counts = vst.counts[sel, ]
compute.corr(vst.counts, nsamples, colnames(vst.counts), file.path(plotdir, paste('corrMat_vst_', mark, '.png', sep = '')))
compute.corr(vst.counts, nsamples, colnames(vst.counts), file.path(plotdir, paste('corrMat_vst_eucl_', mark, '.png', sep = '')), corr = F)

########################## Use the fits 
fits = fit.coef[diff.count, 4:dim(fit.coef)[2]]
colnames(fits) = gsub('condition', '', colnames(fits))
compute.corr(fits, dim(fits)[2], colnames(fits), file.path(plotdir, paste('corrMat_fits_', mark, '.png', sep = '')))
plot.heatmap(fits, row.cluster = T, col.cluster = T, show.dendro = "none", dist.metric = "euclidean", clust.method = "ward", 
             break.type = "quantile", filt.thresh = 0, replace.na = F, to.file = file.path(plotdir, paste('biclust_fits_', mark, '.png', sep = '')))

########################## Use the quantile normalized counts
counts.quant = normalize.quantiles(counts.all)
row.means = rowMeans(counts.quant)
counts.quant = apply(counts.quant, 2, function(x) (x - row.means))
colnames(counts.quant) = colnames(counts.all)
pca.fit = prcomp(t(counts.quant), center = T, scale = F)
png(filename = file.path(plotdir, 'pca2d_qn.png'), width = 800, height = 800)
plot(pca.fit$x[,1], pca.fit$x[,2], xlab = 'PC1', ylab = 'PC2', main = '')
text(pca.fit$x[,1] + 0.2, pca.fit$x[,2] - 0.25, colnames(counts.all), cex = 0.8)
dev.off()
plot.heatmap(counts.quant, row.cluster = T, col.cluster = T, show.dendro = "none", dist.metric = "euclidean", clust.method = "ward", 
             break.type = "quantile", filt.thresh = min(counts.quant), replace.na = F, to.file = file.path(plotdir, paste('biclust_qn_', mark, '.png', sep = '')))

########################## Use the counts (normalized by region length and size factors)
pca.fit = prcomp(t(counts.all), center = T, scale = F)
png(filename = file.path(plotdir, 'pca2d.png'), width = 800, height = 800)
plot(pca.fit$x[,1], pca.fit$x[,2], xlab = 'PC1', ylab = 'PC2', main = '')
text(pca.fit$x[,1] + 0.2, pca.fit$x[,2] - 0.25, colnames(counts.all), cex = 0.8)
dev.off()
png(filename = file.path(plotdir, 'pca2d_2.png'), width = 800, height = 800)
plot(pca.fit$x[,2], pca.fit$x[,3], xlab = 'PC2', ylab = 'PC3', main = '')
text(pca.fit$x[,2] + 0.2, pca.fit$x[,3] - 0.25, colnames(counts.all), cex = 0.8)
dev.off()

compute.corr(asinh(counts.all), nsamples, colnames(counts.all), file.path(plotdir, paste('corrMat_', mark, '.png', sep = '')))
plot.heatmap(asinh(counts.all), row.cluster = T, col.cluster = T, show.dendro = "none", dist.metric = "euclidean", clust.method = "ward", 
             break.type = "quantile", filt.thresh = 0, replace.na = F, to.file = file.path(plotdir, paste('biclust_', mark, '.png', sep = '')))
plot.heatmap(pca.fit$x, row.cluster = T, col.cluster = T, show.dendro = "none", dist.metric = "euclidean", clust.method = "ward", 
             break.type = "quantile", filt.thresh = min(pca.fit$x), replace.na = F, to.file = file.path(plotdir, paste('biclust_pc_', mark, '.png', sep = '')))

ica.fit = fastICA(pca.fit$x, n.comp = dim(counts.all)[2], method = 'C', verbose = T)
png(filename = file.path(plotdir, 'ica2d.png'), width = 800, height = 800)
plot(ica.fit$S[,1], ica.fit$S[,2], xlab = 'IC1', ylab = 'IC2', main = '')
text(ica.fit$S[,1] + 0.2, ica.fit$S[,2] - 0.25, colnames(counts.all), cex = 0.8)
dev.off()
########################## Use mean normalized counts
row.means = rowMeans(counts.all)
row.sds = rowSds(counts.all)
counts.norm = apply(counts.all, 2, function(x) (x - row.means)/row.sds)
compute.corr(counts.norm, nsamples, colnames(counts.norm), file.path(plotdir, paste('corrMat_norm_eucl_', mark, '.png', sep = '')), corr = F)
plot.heatmap(counts.norm, row.cluster = T, col.cluster = T, show.dendro = "none", dist.metric = "euclidean", clust.method = "ward", 
             break.type = "quantile", filt.thresh = min(counts.norm), replace.na = F, to.file = file.path(plotdir, paste('biclust_norm_', mark, '.png', sep = '')))

