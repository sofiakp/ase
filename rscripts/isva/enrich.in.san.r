rm(list=ls())
library(matrixStats)
library(reshape)
library(ggplot2)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# Compares population specific clusters to segmentation states

isva.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig_withSan/rdata/') 
outpref = 'SNYDER_HG19_all_reg_' #'hg19_w10k_' #gencode.v13.annotation.noM.genes_' #
suf = '_qn_isvaNull'
plotdir = file.path(isva.dir, '..', 'plots', 'qn_isvaNull_fits_all_reg_withSan')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)

# Read San specific SNPs
san.spec = read.table('../../rawdata/mapped/Illumina/round2/snps/sanConsensus/sanSpecific.txt', header = T, sep = '\t')
colnames(san.spec) = c('chr', 'pos')

mark = 'H3K27AC'

# Reference clusters
load(file.path(isva.dir, paste(outpref, mark, suf, '.RData', sep = '')))
indivs = colnames(counts)
nindiv = length(indivs)
load(file.path(isva.dir, paste(outpref, mark, suf, '_clust.RData', sep = '')))

indivs[indivs == 'SNYDER'] = 'MS1'

all.san.ov = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), snps.to.ranges(san.spec), select = 'first', ignore.strand = T)
have.hit = sum(!is.na(all.san.ov))

enrich = array(1, dim = c(max(kord), 2))
for(k in 1:max(kord)){
  mem = kclusters$cluster == kord[k]
  cregions = regions[isva.fit$deg[mem], ]
  deg = findOverlaps(regions.to.ranges(cregions), snps.to.ranges(san.spec), ignore.strand = T, select = 'first')  
  ndeg = sum(!is.na(deg))
  print(ndeg / sum(mem))
  enrich[k, 1] = (ndeg / sum(mem)) / (have.hit / length(all.san.ov))
  enrich[k, 2] = binom.val(ndeg, sum(mem), have.hit / length(all.san.ov))
}