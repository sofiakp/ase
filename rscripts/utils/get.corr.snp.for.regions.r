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

# For a set of regions and a corresponding null set, finds the SNP most correlated with the region's signal.

set.seed(1)

nchunks = 6
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
indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata'
outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/anova/'
counts.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/textFiles/'
if(!file.exists(outdir)) dir.create(outdir)
mark = 'H3K4ME3'
load(file.path(indir, 'anova', paste(mark, 'regionsFvalues.rda', sep = '.')))
pos.pref = paste(mark, 'sqrt_sig4_Fgt1_SDgt1', sep = '_')
neg.pref = paste(mark, 'sqrt_sig4_Flt1_SDlt1', sep = '_')
out.pref = paste(mark, 'sqrt_sig4', sep = '_')
region.file = paste('../../rawdata/signal/combrep/peakFiles/merged/SNYDER_HG19', mark, 'merged.bed.gz', sep = '_')
signal.files = list.files(counts.dir, pattern = paste(gsub('.bed|.bed.gz', '', basename(region.file)), '_AT_SNYDER_HG19_.*', mark, '.*.txt', sep = ''), full.names = T)
indivs = unique(gsub(paste('.*_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
nindivs = length(indivs)
counts.dat = load.avg.sig.data(region.file, signal.files, indivs) 
regions = counts.dat$regions

counts = sqrt(counts.dat$signal) + sqrt(counts.dat$signal + 1)
counts = normalize.quantiles(counts)
colnames(counts) = indivs 
row.means = rowMeans(counts)
print(quantile(row.means, na.rm = T))
large.dev = array(0, dim = c(nrow(regions), 1))
for(i in 1:nrow(regions)){
  tmp.sig = as.numeric(counts[i, ])
  tmp.sig[is.na(tmp.sig)] = 0
  tmp.sort = sort(tmp.sig + 1)
  large.dev[i] = tmp.sort[nindivs - 2] / tmp.sort[3]
}
good.rows = apply(counts, 1, function(x) !any(is.na(x))) & row.means > 4 #quantile(row.means, 0.75, na.rm = T)

#counts = counts[good.rows, ]
#regions = regions[good.rows, ]
#regionsFvalues = regionsFvalues[good.rows, ]
#row.sds = rowSds(counts)
# cvs = row.sds / row.means
# good.rows = !is.na(row.means) & 
#!is.na(cvs) #asinh(0.2) & cvs > quantile(cvs, 0.4, na.rm = T)

row.means = row.means[good.rows]
counts = counts[good.rows, ]
regions = regions[good.rows, ]
regionsFvalues = regionsFvalues[good.rows, ]
regions$logF = regionsFvalues$logF
regions$dev = large.dev[good.rows]
regions$sd = rowSds(counts)
# write.table(regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(out.pref, 'regionsWithF.txt', sep = '_')))
# 
# # For each non-variable region, pick a variable region with similar signal.
neg.idx = which(!is.na(regions$logF) & (regions$logF < 0.5))# & regions$sd < 0.5)) #quantile(regionsFvalues$logF, 0.2, na.rm = T)
# # Shuffle neg.idx. If good.idx has much fewer elements than neg.idx, and you don't shuffle, then the 
# # resulting lists will be restricted towards the beginning of the genome.
# neg.idx = neg.idx[sample(1:length(neg.idx), length(neg.idx), replace = F)] 
good.idx = !is.na(regions$logF) & regions$logF > 1.5 #& regions$sd > 1 # Variable regions to choose from
# cat('# neg regions: ', length(neg.idx), '\n')
# cat('# candidate pos regions: ', sum(good.idx), '\n')
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
pos.regions = regions[pos.idx, ]
neg.regions = regions[neg.idx, ]
pos.counts = counts[pos.idx, ]
neg.counts = counts[neg.idx, ]
pos.name = 'Variable'
neg.name = 'Non-variable'
cat('Signal diff p-value: ', wilcox.test(rowMeans(pos.counts), rowMeans(neg.counts))$p.value, '\n')
# cat('Signal diff p-value (paired): ', wilcox.test(rowMeans(pos.counts), rowMeans(neg.counts), paired = T)$p.value, '\n')
# write.table(pos.regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(pos.pref, 'regionsWithF.txt', sep = '_')))
# write.table(neg.regions, row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(neg.pref, 'regionsWithF.txt', sep = '_')))

# # Read SNP information
# snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
# load(snp.pos.file)
# 
# # Load genotypes
# load('../../rawdata/variants/all/snps/allNonSan/all_genot.RData')
# genot = genot[, match(indivs, colnames(genot))]
# sel.genot = rowSums(genot == 0) < nindivs & rowSums(genot == 1) < nindivs & rowSums(genot == 2) < nindivs # Select positions where individuals differ
# genot = genot[sel.genot, ]
# snp.pos = snp.pos[sel.genot, ]
# 
# # Find all overlaps between regions and SNPs
# ov = findOverlaps(regions.to.ranges(pos.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
# ov.mat = cbind(queryHits(ov), subjectHits(ov))
# ov.neg = findOverlaps(regions.to.ranges(neg.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
# ov.mat.neg = cbind(queryHits(ov.neg), subjectHits(ov.neg))
# 
# # Correlation between signal at regions and genotype of all overlapping SNPs
# cor = cor.par(genot[ov.mat[, 2], ], pos.counts[ov.mat[, 1], ], nchunks, 'spearman')
# neg.cor = cor.par(genot[ov.mat.neg[, 2], ], neg.counts[ov.mat.neg[, 1], ], nchunks, 'spearman')
# 
# # For each region, get the SNP with the highest correlation. Then, test if the SNPs in the true regions have higher correlations than
# # the SNPs in the control regions.
# cor.dat.tmp = data.frame(cor = cor, region.idx = ov.mat[, 1], chr = snp.pos$chr[ov.mat[, 2]], pos = snp.pos$pos[ov.mat[, 2]],
#                          start = pos.regions$start[ov.mat[, 1]], end = pos.regions$end[ov.mat[, 1]])
# cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x)) * sign(x[which.max(abs(x))]), value = 'cor') # get max cor per region
# colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
# # create matrix region, cor, snp_chr, snp_pos, reg_start, reg_end. 
# # This might have more than one SNP per region, if there were multiple SNPs in the same region
# # with equally high correlation
# cor.dat.true = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor')) 
# cor.dat.true = cor.dat.true[!duplicated(cor.dat.true[,1:2]), ]
# #write.table(cor.dat.true[, c(3,5,6,4,2)], row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(pos.pref, '_bestSNP.txt', sep = '')))
# 
# cor.dat.tmp = data.frame(cor = neg.cor, region.idx = ov.mat.neg[, 1], chr = snp.pos$chr[ov.mat.neg[, 2]], 
#                          pos = snp.pos$pos[ov.mat.neg[, 2]], 
#                          start = neg.regions$start[ov.mat.neg[, 1]], end = neg.regions$end[ov.mat.neg[, 1]])
# cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x)) * sign(x[which.max(abs(x))]), value = 'cor')
# colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
# cor.dat.perm = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor'))
# cor.dat.perm = cor.dat.perm[!duplicated(cor.dat.perm[,1:2]), ]
# #write.table(cor.dat.perm[, c(3,5,6,4,2)], row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(outdir, paste(neg.pref, '_bestSNP.txt', sep = '')))
# 
# cor.dat = rbind(data.frame(cor = abs(cor.dat.true[, 2]), type = rep(pos.name, nrow(cor.dat.true))), 
#                 data.frame(cor = abs(cor.dat.perm[, 2]), type = rep(neg.name, nrow(cor.dat.perm))))
# 
# wilc.p = wilcox.test(abs(cor.dat.true$cor), abs(cor.dat.perm$cor))$p.value
# p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type), size = 1, adjust = 2) +
#   annotate('text', 0.75, 1.5, label = sprintf('Wilcoxon P = %g', wilc.p), size = 5) +
#   xlab('Correlation between genotype and signal') + ylab('Density') + theme_bw() + 
#   theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
#         legend.position = c(.75, .85), legend.title = element_blank(), legend.text = element_text(size = 16))
# #ggsave(file = file.path(outdir, paste(pos.pref, '_snpGeneCorr.pdf', sep = '')), p, width= 6.5, height = 5.6)