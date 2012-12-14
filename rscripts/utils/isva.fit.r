rm(list=ls())
library('DESeq')
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

# Clustering and visualization of signal in regions (exons, peak regions etc)

counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomeGrid/hg19_w1k/fc/avgSig/') #'rawdata/signal/combrep/extractSignal/fc/avgSig/')
#rawdata/geneCounts/rdata/repsComb/') #'rawdata/signal/combrep/countsAtPeaksBroad/repsComb')
outdir = file.path(counts.dir, 'rdata') # Set to null if you don't want to write the merged data to a file
plotdir = file.path(counts.dir, 'plots', 'qn_isvaNull_fits')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
counts.dir = file.path(counts.dir, 'textFiles')
outpref = 'hg19_w1k_'  #'SNYDER_HG19_'
mark = 'H3K4ME3'

# region.file: BED file with regions to read. 
# signal.files: should be txt files with just one column of values with the signal in each of the regions in region.file

#region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/', 
#                        paste('SNYDER_HG19', mark, 'merged.bed.gz', sep = '_'))
#signal.files = list.files(counts.dir, pattern = paste('SNYDER_HG19', mark, 'merged_AT_.*txt', sep = '_'), full.names=T)
#indivs = unique(gsub(paste('SNYDER_HG19_', mark, '_merged_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomeGrid/hg19_w1k.bed')
signal.files = list.files(counts.dir, pattern = paste('hg19_w1k_AT_SNYDER_HG19_.*', mark, '.txt', sep = ''), full.names = T)
indivs = unique(gsub(paste('hg19_w1k_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
nindivs = length(indivs)
counts.dat = load.avg.sig.data(region.file, signal.files, indivs)
regions = counts.dat$regions
counts = asinh(counts.dat$signal) 
counts[is.na(counts)] = 0

############# Quantile normalization
counts = normalize.quantiles(counts)
colnames(counts) = indivs
counts.norm = scale(counts)

############### ISVA correction to remove batch effects
pop = factor(get.pop(indivs))
isva.fit = DoISVA(counts.norm, pop, pvth = 0.05, th = 0.1)
counts = normalize.quantiles(isva.fit$res.null) # Get residuals after removing ISVs and renormalize
colnames(counts) = indivs
save(regions, counts, isva.fit, file = file.path(outdir, paste(outpref, mark, '_qn_isvaNull_all.RData', sep = ''))) 

isv.dat = data.frame(isva.fit$isv)
colnames(isv.dat) = paste('ISV', 1:ncol(isv.dat))
isv.dat$indiv = factor(indivs)
isv.dat = melt(isv.dat, id.vars=c('indiv'))
p = ggplot(isv.dat) + geom_bar(aes(x = indiv, y = value), position = 'dodge', stat = 'identity') + facet_wrap(~variable) +
  xlab('') + ylab('ISV value') + 
  theme(strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = -45, vjust = 1, hjust = 0))
ggsave(file.path(plotdir, paste('isvs_', mark, '_all.png', sep = '')), p, width = 13.6, height = 11.8)
