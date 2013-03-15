rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

nchunks = 5
registerDoMC(nchunks)

peak.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/rdata')
deseq.dir = '../../rawdata/signal/combrep/countsAtPeaksBroad/repsComb/deseq/'
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')
plotdir = file.path(peak.dir, '../plots')
if(!file.exists(plotdir)) dir.create(plotdir)

# Load SNPs
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

# Load signal
files = list.files(peak.dir, pattern = 'SNYDER_HG19_[A-Z0-9]*_qn.RData', full.names = T)
marks = gsub('SNYDER_HG19_|_qn.RData', '', basename(files))
div.all = NULL
for(m in 1:length(files)){
  load(files[m])
  # sel.col = !(colnames(counts) %in% c('GM12878', 'GM19240'))
  # counts = counts[, sel.col]
  indivs = colnames(counts)
  nindivs = length(indivs)
  good = apply(counts, 1, function(x) !any(is.na(x)))
  counts = counts[good, ]
  regions = regions[good, ]
  # Get a random subset of regions
  sel = sample(1:nrow(counts), min(nrow(counts), 5000))
  counts = counts[sel, ]
  regions = regions[sel, ]
  widths = regions$end - regions$start + 1
  count.var = rowVars(counts)
  
  # Load DESeq results
  # files = list.files(deseq.dir, pattern = paste(mark, '*_deseq.RData', sep = ''), full.names = T)
  # print(length(files))
  # diff.counts = get.diff.count(files, 0.0001)
  # diff.counts = diff.counts[good][sel]
  
  # Overlap the SNPs and the regions
  ranges = regions.to.ranges(regions)
  ov = findOverlaps(ranges, snps.to.ranges(snp.pos), ignore.strand = T, select = 'all')
  ov = cbind(queryHits(ov), subjectHits(ov)) # Convert to matrix for fast access. Col 1 - region idx, Col 2 - SNP idx
  nov = dim(ov)[1]
  
  # Load the genotypes for the SNPs that do have overlap
  genot = array(F, dim = c(dim(ov)[1], nindivs * 2))
  snp.env = new.env()
  for(i in 1:nindivs){
    cat('Loading genotypes for', indivs[i], '\n')
    load(file.path(geno.dir, paste(indivs[i], 'snps.RData', sep = '.')), snp.env)
    genot[, (i - 1) * 2 + 1] = as.vector(snp.env$geno.info$mat[ov[,2]])
    genot[, i * 2] = as.vector(snp.env$geno.info$pat[ov[,2]])
  }
  # genot = array(F, dim = c(dim(ov)[1], nindivs))
  # snp.env = new.env()
  # for(i in 1:nindivs){
  #   cat('Loading genotypes for', indivs[i], '\n')
  #   load(file.path(geno.dir, paste(indivs[i], 'snps.RData', sep = '.')), snp.env)
  #   genot[, i] = as.vector(snp.env$geno.info$mat[ov[,2]]) + as.vector(snp.env$geno.info$pat[ov[,2]])
  # }
  
  # snp.counts.tmp = table(ov[, 1])
  # snp.counts = array(0, dim = c(length(widths), 1))
  # snp.counts[as.numeric(names(snp.counts.tmp))] = snp.counts.tmp
  # var.dat = data.frame(snps = nhaps, var = diff.counts)
  # var.dat = var.dat[var.dat$snps < 6, ]
  # var.dat$snps = factor(var.dat$snps)
  # p = ggplot(var.dat) + geom_density(aes(color = snps, x = var)) + xlim(c(0, 10))
  # 
  # pind = unique(ov[, 1])
  # nhap = array(0, dim = c(nrow(counts), 1))
  # for(i in 1:length(pind)){ # For each region...
  #   sel.pos = findInterval(c(pind[i] - 1, pind[i]), ov[, 1]) # Indices of SNPs overlapping the pind[i]-th region
  #   snp.ind = (sel.pos[1] + 1):sel.pos[2] # pind[i] appears in ov, so this interval makes sense
  #   nsel = length(snp.ind) # Number of snps overlapping the region
  #   if(nsel > 1){
  #     nhap[pind[i]] = nrow(unique(genot[snp.ind, ]))
  #   }else{
  #     nhap[pind[i]] = length(unique(genot[snp.ind, ]))
  #   }
  # }
  
  div = get.div.par(genot, ov[, 1], widths, nchunks)
  #div.dat = data.frame(div = div, var = count.var)
  #p = ggplot(div.dat) + 
  div.dat = NULL
  q = c(0, .25, .5, .75, 1)
  for(i in 2:length(q)){
    sel.var = count.var[div > quantile(div, q[i - 1]) & div <= quantile(div, q[i])]
    div.dat = rbind(div.dat, data.frame(quant = rep(i - 1, length(sel.var)), var = sel.var))
    if(i > 2 && length(div.dat$var[div.dat$quant == i - 2]) > 0) print(wilcox.test(sel.var, div.dat$var[div.dat$quant == i - 2]))
  }
  div.dat = div.dat[div.dat$quant == 1 | div.dat$quant == 4, ]
  div.dat$quant = droplevels(factor(div.dat$quant, levels = c(1,4), labels = c('ND < 25th percentile', 'ND > 75th percentile')))
  if(length(levels(div.dat$quant)) > 1){
    div.dat$mark = rep(marks[m], nrow(div.dat))
    p = ggplot(div.dat) + geom_density(aes(x = var, color = quant), adjust = 2) + 
      xlab('Signal variance') + ylab('density') + theme_bw() + scale_color_discrete('') + xlim(c(0, 1)) + ggtitle(marks[m]) + 
      theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), 
            axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
            legend.position = c(.8, .82), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
    ggsave(file.path(plotdir, paste(marks[m], '_var_vs_nd.png', sep = '')), p, width = 6.5, height = 5.6)
    div.all = rbind(div.all, div.dat) 
  }
}
p = ggplot(div.all) + geom_density(aes(x = var, color = quant), size = 1, adjust = 2) + facet_wrap(~mark) +
  xlab('Signal variance') + ylab('density') + theme_bw() + scale_color_discrete('') + xlim(c(0, 1))  +
  theme(axis.text.x = element_text(size = 14, angle = -45, hjust = 0, vjust = 1), axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        strip.text.x = element_text(size = 16),
        legend.position = c(.82, .2), legend.text = element_text(size = 14), legend.title = element_text(size = 14))
ggsave(file.path(plotdir, paste('all', '_var_vs_nd.png', sep = '')), p, width = 9.75, height = 8.4)

