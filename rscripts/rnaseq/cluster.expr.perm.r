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
counts.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig_withSan/')
#counts.dir= file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/repsComb/')
#deseq.dir = file.path(counts.dir, 'deseq/')
########### CHANGE THIS !!!!!!!!!!!!!!!!!!!!!
outpref = 'SNYDER_HG19_all_reg_' 
#outpref = 'gencode.v13.annotation.noM.genes_all_reg_'
#outpref = 'txStates_10_11_12_'
#outpref = 'hg19_w10k_'
#outpref = 'all_reg_'

plotdir = file.path(counts.dir, 'plots', 'qn_isvaNull_fits_all_reg_withSan')
outdir = file.path(counts.dir, 'rdata')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)
if(!file.exists(outdir)) dir.create(outdir)
counts.dir = file.path(counts.dir, 'textFiles')
nperm = 10
is.genes = F # T for RZ data
quant = 0.4 # 0.4 for peak regions and transcriptomes
mark = 'H3K4ME1'

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
  
  region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged_withSan/', paste('SNYDER_HG19', mark, 'merged.bed.gz', sep = '_'))
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

indivs[indivs == 'SNYDER'] = 'MS1'
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

############### ISVA correction to remove batch effects
pop = factor(get.pop(indivs))
print(pop)

for(n in 1:nperm){
  sam = sample(1:nindivs, nindivs)
  print(pop[sam])
  isva.fit = DoISVA(scale(orig.counts[, sam]), pop, pvth = 0.05, th = 0.05, ncomp = 3)
  counts = normalize.quantiles(isva.fit$res.null) # Get residuals after removing SVs and renormalize
  colnames(counts) = indivs
  counts.norm = scale(counts)
  if(!is.null(outdir)) save(regions, counts, isva.fit, sam, file = file.path(outdir, paste(outpref, 'rand_pop', n, '_', mark, '_qn_isvaNull.RData', sep = ''))) 
}