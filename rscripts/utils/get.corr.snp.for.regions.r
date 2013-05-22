rm(list=ls())
library(reshape)
library(GenomicRanges)
library(ggplot2)
library(foreach)
library(doMC)
library(preprocessCore)
library(matrixStats)
source('utils/deseq.utils.r')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

get.cor.snp = function(cor, ov.mat, regions, snp.pos){
  # For each region, get the SNP with the highest correlation.
  cor.dat.tmp = data.frame(cor = cor, region.idx = ov.mat[, 1], chr = snp.pos$chr[ov.mat[, 2]], pos = snp.pos$pos[ov.mat[, 2]],
                           start = regions$start[ov.mat[, 1]], end = regions$end[ov.mat[, 1]])
  cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x)) * sign(x[which.max(abs(x))]), value = 'cor') # get max cor per region
  colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
  # create matrix region, cor, snp_chr, snp_pos, reg_start, reg_end. 
  # This might have more than one SNP per region, if there were multiple SNPs in the same region
  # with equally high correlation
  cor.dat.true = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor')) 
  cor.dat.true = cor.dat.true[!duplicated(cor.dat.true[,1:2]), ]
  return(cor.dat.true)
}
# For a set of regions and a corresponding null set, finds the SNP most correlated with the region's signal.

set.seed(1)

nchunks = 3
registerDoMC(nchunks)

# Load population specific regions and corresponding "null" regions
# indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata'
# outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/qn_isvaNull_fits_all_reg/'
# if(!file.exists(outdir)) dir.create(outdir)
# mark = 'H3K27AC'
# load(file.path(indir, paste('SNYDER_HG19_all_reg', mark, 'qn_isvaNull.RData', sep = '_')))
# pos.pref = paste('SNYDER_HG19_all_reg', mark, 'qn_isvaNull_sign', sep = '_')
# neg.pref = paste('SNYDER_HG19_all_reg', mark, 'qn_isvaNull_non_sign', sep = '_')
# pos.regions = regions[isva.fit$deg, ] # Population specific regions
# neg.regions = regions[-isva.fit$deg, ] # Negative set

# pos.counts = counts[isva.fit$deg, ]
# neg.counts = counts[-isva.fit$deg, ]
# pos.name = 'Population-specific'
# neg.name = 'Other'
# indivs = colnames(counts)
# nindivs = length(indivs)

# Load variable and non variable regions
indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata'
outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/anova/'
counts.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/textFiles/'
if(!file.exists(outdir)) dir.create(outdir)
mark = 'H3K27AC'
get.max.cor = T
shuf.regions = T
ext = 1000
half.len = 1000

load(file.path(indir, 'anova', paste(mark, 'regionsFvalues.rda', sep = '.')))
pos.pref = paste(mark, 'sig0.75_up0.15F', sep = '_')
neg.pref = paste(mark, 'sig0.75_low0.15F', sep = '_')
out.pref = paste(mark, 'sig0.75_up_vs_low_0.15F_shufRegions', sep = '_')
region.file = paste('../../rawdata/signal/combrep/peakFiles/merged_Mar13/SNYDER_HG19', mark, 'merged.bed.gz', sep = '_')
signal.files = list.files(counts.dir, pattern = paste(gsub('.bed|.bed.gz', '', basename(region.file)), '_AT_SNYDER_HG19_.*', mark, '.*.txt', sep = ''), full.names = T)
indivs = unique(gsub(paste('.*_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
#sel.indivs = get.pop(indivs) != 'San'
#indivs = indivs[sel.indivs]
#signal.files = signal.files[sel.indivs]
nindivs = length(indivs)
counts.dat = load.avg.sig.data(region.file, signal.files, indivs) 
regions = counts.dat$regions

counts = asinh(counts.dat$signal) #sqrt(counts.dat$signal) + sqrt(counts.dat$signal + 1)
counts = normalize.quantiles(counts)
colnames(counts) = indivs 
# large.dev = array(0, dim = c(nrow(regions), 1))
# for(i in 1:nrow(regions)){
#   tmp.sig = as.numeric(counts[i, ])
#   tmp.sig[is.na(tmp.sig)] = 0
#   tmp.sort = sort(tmp.sig + 1)
#   large.dev[i] = tmp.sort[nindivs - 2] / tmp.sort[3]
# }
colnames(regionsFvalues)[2] = 'start'
colnames(regionsFvalues)[3] = 'end'
ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(regionsFvalues), select = 'all', ignore.strand = T)
stopifnot(length(queryHits(ov)) == min(nrow(regions), nrow(regionsFvalues)))
regions = regions[queryHits(ov), ]
counts = counts[queryHits(ov), ]
regions$logF = regionsFvalues$logF[subjectHits(ov)]

row.means = rowMeans(counts)
print(quantile(row.means, na.rm = T))
good.rows = apply(counts, 1, function(x) !any(is.na(x))) & row.means > quantile(row.means, 0.75, na.rm = T)
counts = counts[good.rows, ]
regions = regions[good.rows, ]
row.sds = rowSds(counts)
row.means = row.means[good.rows]
cvs = row.sds / row.means
#good.rows = !is.na(row.means) & !is.na(cvs) #& cvs > quantile(cvs, 0.4, na.rm = T)
#row.means = row.means[good.rows]
#counts = counts[good.rows, ]
#regions = regions[good.rows, ]

#regions$dev = large.dev[good.rows]
#regions$sd = rowSds(counts)
# write.table(regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(out.pref, 'regionsWithF.txt', sep = '_')))
# 
# # For each non-variable region, pick a variable region with similar signal.
neg.idx = which(!is.na(regions$logF) & regions$logF < quantile(regions$logF, 0.15, na.rm = T)) #which(cvs < quantile(cvs, 0.15, na.rm = T)) #
# # Shuffle neg.idx. If good.idx has much fewer elements than neg.idx, and you don't shuffle, then the 
# # resulting lists will be restricted towards the beginning of the genome.
# neg.idx = neg.idx[sample(1:length(neg.idx), length(neg.idx), replace = F)] 
good.idx = !is.na(regions$logF) & regions$logF > quantile(regions$logF, 0.85, na.rm = T) # Variable regions to choose from #cvs > quantile(cvs, 0.85, na.rm = T)
cat('# neg regions: ', length(neg.idx), '\n')
cat('# candidate pos regions: ', sum(good.idx), '\n')
# 
# pos.idx = array(0, dim = c(length(neg.idx), 1))
# for (i in 1:length(neg.idx)){
#   if(!any(good.idx)) break
#   min.diff = min(abs(row.means[neg.idx[i]] - row.means[good.idx]))
#   if (min.diff < 0.001){
#     sel.idx = which.min(abs(row.means[neg.idx[i]] - row.means[good.idx]))
#     pos.idx[i] = which(good.idx)[sel.idx]
#     good.idx[which(good.idx)[sel.idx]] = F # Make the selected region non-available for future selection
#   }
# }
# neg.idx = neg.idx[pos.idx > 0] # Remove the non-variable regions for which you didn't find a good match
# pos.idx = pos.idx[pos.idx > 0]
# cat('# matched neg regions: ', length(neg.idx), '\n')
# cat('# matched pos regions: ', length(pos.idx), '\n')
# #pos.idx = !is.na(regionsFvalues$logF) & regionsFvalues$logF >= 1 #quantile(regionsFvalues$logF, 0.8, na.rm = T)
# 
pos.idx = good.idx ####

regions.mid = round((regions$start + regions$end) / 2)
regions$start = pmax(regions.mid - half.len, 1)
regions$end = regions.mid + half.len

pos.regions = regions[pos.idx, ]
ext.pos.regions = pos.regions
ext.pos.regions$start = pmax(ext.pos.regions$start - ext, 1)
ext.pos.regions$end = ext.pos.regions$end + ext
neg.regions = regions[neg.idx, ]
pos.counts = counts[pos.idx, ]
neg.counts = counts[neg.idx, ]
pos.name = 'variable'
neg.name = 'non-variable'
cat('Signal diff p-value: ', wilcox.test(rowMeans(pos.counts), rowMeans(neg.counts))$p.value, '\n')
# cat('Signal diff p-value (paired): ', wilcox.test(rowMeans(pos.counts), rowMeans(neg.counts), paired = T)$p.value, '\n')
# write.table(pos.regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(pos.pref, 'regionsWithF.txt', sep = '_')))
# write.table(neg.regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(neg.pref, 'regionsWithF.txt', sep = '_')))

# Read SNP information
snp.pos.file = '../../rawdata/variants/all_Mar13/snps.RData'
load(snp.pos.file)

# Load genotypes
load('../../rawdata/variants/all_Mar13/genot.RData')
genot = genot[, match(indivs, colnames(genot))]
sel.genot = rowSums(genot == 0) < nindivs - 1 & rowSums(genot == 1) < nindivs - 1 & rowSums(genot == 2) < nindivs - 1 # Select positions where individuals differ
genot = genot[sel.genot, ]
snp.pos = snp.pos[sel.genot, ]
shuf.genot = genot[sample(1:nrow(genot), nrow(genot)), ]
shuf.pos.counts = pos.counts[sample(1:nrow(pos.counts), nrow(pos.counts)), ]
shuf.neg.counts = neg.counts[sample(1:nrow(neg.counts), nrow(neg.counts)), ]

# Find all overlaps between regions and SNPs
ov = findOverlaps(regions.to.ranges(pos.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov))
ext.ov = findOverlaps(regions.to.ranges(ext.pos.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ext.ov.mat = cbind(queryHits(ext.ov), subjectHits(ext.ov))
ov.neg = findOverlaps(regions.to.ranges(neg.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat.neg = cbind(queryHits(ov.neg), subjectHits(ov.neg))

# Correlation between signal at regions and genotype of all overlapping SNPs
cor = cor.par(genot[ov.mat[, 2], ], pos.counts[ov.mat[, 1], ], nchunks, 'spearman')
ext.cor = cor.par(genot[ext.ov.mat[, 2], ], pos.counts[ext.ov.mat[, 1], ], nchunks, 'spearman')
neg.cor = cor.par(genot[ov.mat.neg[, 2], ], neg.counts[ov.mat.neg[, 1], ], nchunks, 'spearman')
if(shuf.regions){
  shuf.cor = cor.par(genot[ov.mat[, 2], ], shuf.pos.counts[ov.mat[, 1], ], nchunks, 'spearman')  
  shuf.ext.cor = cor.par(genot[ext.ov.mat[, 2], ], shuf.pos.counts[ext.ov.mat[, 1], ], nchunks, 'spearman')  
  shuf.neg.cor = cor.par(genot[ov.mat.neg[, 2], ], shuf.neg.counts[ov.mat.neg[, 1], ], nchunks, 'spearman')
}else{
  shuf.cor = cor.par(shuf.genot[ov.mat[, 2], ], pos.counts[ov.mat[, 1], ], nchunks, 'spearman')
  shuf.neg.cor = cor.par(shuf.genot[ov.mat.neg[, 2], ], neg.counts[ov.mat.neg[, 1], ], nchunks, 'spearman')
}

if(get.max.cor){
  cor.dat.true = get.cor.snp(cor, ov.mat, pos.regions, snp.pos)
  cor.dat.ext = get.cor.snp(ext.cor, ext.ov.mat, ext.pos.regions, snp.pos)
  cor.dat.shuf = get.cor.snp(shuf.cor, ov.mat, pos.regions, snp.pos)
  cor.dat.shuf.ext = get.cor.snp(shuf.ext.cor, ext.ov.mat, ext.pos.regions, snp.pos)
}else{
  cor.dat.true = data.frame(region.idx = ov.mat[, 1], cor = abs(cor))
  cor.dat.shuf = data.frame(region.idx = ov.mat[, 1], cor = abs(shuf.cor))
}
cor.dat = rbind(data.frame(cor = abs(cor.dat.true[, 2]), type = rep(pos.name, nrow(cor.dat.true))), 
                data.frame(cor = abs(cor.dat.shuf[, 2]), type = rep(paste('shuffled', pos.name), nrow(cor.dat.shuf))))

wilc.p = wilcox.test(abs(cor.dat.true$cor), abs(cor.dat.shuf$cor))$p.value
cat('Pos set, true vs shuffled wilcoxon and medians', wilc.p, median(abs(cor.dat.true$cor)), median(abs(cor.dat.shuf$cor)), '\n')
cat('Pos set, true vs extended wilcoxon and medians', wilcox.test(abs(cor.dat.true$cor), abs(cor.dat.ext$cor))$p.value, 
    median(abs(cor.dat.true$cor)), median(abs(cor.dat.ext$cor)), '\n')
cat('Pos set, shuffled true vs shuffled extended wilcoxon and medians', wilcox.test(abs(cor.dat.shuf$cor), abs(cor.dat.shuf.ext$cor))$p.value, 
    median(abs(cor.dat.shuf$cor)), median(abs(cor.dat.shuf.ext$cor)), '\n')

# p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type), size = 1, adjust = 1) +
#   annotate('text', 0.75, 1.5, label = sprintf('Wilcoxon P = %g', wilc.p), size = 5) +
#   xlab('Correlation between genotype and signal') + ylab('Density') + theme_bw() + 
#   theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
#         legend.position = c(.75, .85), legend.title = element_blank(), legend.text = element_text(size = 16))

if(get.max.cor){
  cor.dat.neg = get.cor.snp(neg.cor, ov.mat.neg, neg.regions, snp.pos)
  cor.dat.neg.shuf = get.cor.snp(shuf.neg.cor, ov.mat.neg, neg.regions, snp.pos)
}else{
  cor.dat.neg = data.frame(region.idx = ov.mat.neg[, 1], cor = abs(neg.cor))
  cor.dat.neg.shuf = data.frame(region.idx = ov.mat.neg[, 1], cor = abs(shuf.neg.cor))
}
cor.dat = rbind(cor.dat, data.frame(cor = abs(cor.dat.neg[, 2]), type = rep(neg.name, nrow(cor.dat.neg))), 
                data.frame(cor = abs(cor.dat.neg.shuf[, 2]), type = rep(paste('shuffled', neg.name), nrow(cor.dat.neg.shuf))))

wilc.p = wilcox.test(abs(cor.dat.neg$cor), abs(cor.dat.neg.shuf$cor))$p.value
cat('Neg set, true vs shuffled wilcoxon and medians', wilc.p, median(abs(cor.dat.neg$cor)), median(abs(cor.dat.neg.shuf$cor)), '\n')
cat('Pos vs neg set wilcoxon and medians', wilcox.test(abs(cor.dat.true$cor), abs(cor.dat.neg$cor))$p.value, median(abs(cor.dat.true$cor)), median(abs(cor.dat.neg$cor)), '\n')

p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type), size = 1, adjust = 2) +
  #annotate('text', 0.75, 1.5, label = sprintf('Wilcoxon P = %g', wilc.p), size = 5) +
  xlab('Correlation between genotype and signal') + ylab('Density') + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 16))

#write.table(cor.dat.true[, c(3,5,6,4,2)], row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(pos.pref, '_bestSNP.txt', sep = '')))
ggsave(file = file.path(outdir, paste(out.pref, '_snpGeneCorr.pdf', sep = '')), p, width = 8, height = 5)

p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type), size = 1, stat = 'ecdf', adjust = 2) +
  #annotate('text', 0.75, 1.5, label = sprintf('Wilcoxon P = %g', wilc.p), size = 5) +
  xlab('Correlation between genotype and signal') + ylab('Density') + theme_bw() + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 16))

#write.table(cor.dat.true[, c(3,5,6,4,2)], row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(pos.pref, '_bestSNP.txt', sep = '')))
ggsave(file = file.path(outdir, paste(out.pref, '_snpGeneCorr_cdf.pdf', sep = '')), p, width = 8, height = 5)