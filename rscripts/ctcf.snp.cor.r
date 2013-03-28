rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/'
outpref = 'CTCF_cor_'

# Load SNP information and genotypes
load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load('../../rawdata/variants/all/snps/allNonSan/genot_pca.RData')
snp.pos = snp.pos[pca.rows, ] # This selects rows that were used for PCA (so had some min variance in genotype)
nindivs = ncol(genot)
#sel = rowSums(genot == 0) < nindivs - 3 & rowSums(genot == 1) < nindivs - 3 & rowSums(genot == 2) < nindivs - 3 # Make min variance cutoff stricter
#genot = genot[sel, ]
#snp.pos = snp.pos[sel, ]

# Load normalized signal at peaks
load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/SNYDER_HG19_CTCF_qn.RData')
regions = regions[good.rows, ] # Remove rows with low signal or variance
counts = counts[good.rows, ]
peak.middles = round((regions$start + regions$end) / 2)
stopifnot(all(colnames(counts) == colnames(genot)))

# Read motif locations
mot = read.table('../../rawdata/motifs/personal_genome_motifs_CTCF.txt', header = F, sep = '\t')
mot = mot[, c(1:4, 6)]
colnames(mot) = c('chr', 'start', 'end', 'pos', 'strand') # motif position and position of SNP inside it
mot.middles = round((mot$start + mot$end) / 2)

# Get all the motifs overapping each peak
ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(mot), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov)) # Col 1 - region idx, Col 2 - motif idx
# For each peak, select the motif closest to the peak middle
dist = abs(peak.middles[ov.mat[, 1]] - mot.middles[ov.mat[, 2]])
sel = array(F, dim = c(length(dist), 1))
non.dup = which(!duplicated(ov.mat[, 1])) # rows of ov.mat where region index changes
for(i in 1:length(non.dup)){
  if(i < length(non.dup)){min.pos = which.min(dist[non.dup[i]:(non.dup[i + 1] - 1)])
  }else{min.pos = which.min(dist[non.dup[i]:length(dist)])}
  sel[non.dup[i] + min.pos - 1] = T
}
ov.mat = ov.mat[sel, ]
colnames(ov.mat) = c('region.idx', 'mot.idx')

# Get all the SNPs overlapping the motifs overlapping the peaks
ov2 = findOverlaps(regions.to.ranges(mot[ov.mat[, 2], ]), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat2 = cbind(ov.mat[queryHits(ov2), 2], subjectHits(ov2))
colnames(ov.mat2) = c('mot.idx', 'snp.idx')

# Col 1 - region idx, Col 2 - motif idx, Col 3 - snp idx
ov.idx = merge(ov.mat, ov.mat2, by = 'mot.idx')[, c(2,1,3)]
sel.mot = mot[ov.idx[, 2], ]

dist = sel.mot$pos - sel.mot$start # Position of the SNP in the motif
neg = sel.mot$strand == '-'
dist[neg] = sel.mot[neg, 3] - sel.mot[neg, 4] + 1 # For motifs on the negative strand: end - pos + 1
dist.dat = data.frame(table(dist))
p1 = ggplot(dist.dat) + geom_bar(aes(x = dist, y = Freq)) + 
  xlab('Motif position') + ylab('Number of SNPs') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
ggsave(file.path(outdir, paste(outpref, 'freq_bar.png', sep = '')), p1, width = 6.5, height = 5.6)

# Correlation between peak signal and genotype
cor = cor.par(genot[ov.idx[, 3], ], counts[ov.idx[, 1], ], nchunks = 1, method = 'spearman')
cor.dat = cbind(regions[ov.idx[, 1], ], mot[ov.idx[, 2], c(2:5)], dist, cor)
colnames(cor.dat)[4:5] = c('mot_start', 'mot_end')
cor.dat = cor.dat[order(abs(cor.dat$cor), decreasing = T), ]
write.table(cor.dat, file = file.path(outdir, paste(outpref, 'snps.txt', sep = '')), quote = F, sep = '\t', col.names = T, row.names = F)

cor.dat2 = data.frame(cor = abs(cor), dist = factor(dist))
p2 = ggplot(cor.dat2) + geom_boxplot(aes(x = dist, y = cor)) + 
  xlab('Motif position') + ylab('abs(correlaiton)') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
ggsave(file.path(outdir, paste(outpref, 'cor_box.png', sep = '')), p2, width = 6.5, height = 5.6)

# Overlap with segmentation states
pc = read.bed('../../rawdata/segSignal/14indiv/ctcfPc.bed') # CTCF with K27me3
pc.ov = !is.na(findOverlaps(regions.to.ranges(regions[ov.idx[, 1], ]), regions.to.ranges(pc), select = 'first', ignore.strand = T))
ct = read.bed('../../rawdata/segSignal/14indiv/ctcf.bed') # CTCF state
ct.ov = !is.na(findOverlaps(regions.to.ranges(regions[ov.idx[, 1], ]), regions.to.ranges(ct), select = 'first', ignore.strand = T))

# Divide into 3 types:
# overlapping only CTCF-K27me3 states, overlapping only CTCF states and other
cor.dat3 = rbind(cor.dat2[!pc.ov & ct.ov, ], cor.dat2[pc.ov & !ct.ov, ], cor.dat2[(pc.ov & ct.ov) | (!pc.ov & !ct.ov), ])
cor.dat3$type = factor(rep(c('CTCF', 'CTCF rep', 'mixed'), c(sum(!pc.ov & ct.ov), sum(pc.ov & !ct.ov), sum((pc.ov & ct.ov) | (!pc.ov & !ct.ov)))))
dist.by.type = cast(cor.dat3, dist~type, function(x) length(x), value = 'cor')
# for(i in 1:max(dist)){
#   cat(wilcox.test(cor.dat3$cor[cor.dat3$type == 'CTCF' & as.numeric(cor.dat3$dist) == i], 
#                   cor.dat3$cor[cor.dat3$type == 'CTCF rep' & as.numeric(cor.dat3$dist) == i])$p.value, '\n')
# }
p3 = ggplot(cor.dat3) + geom_boxplot(aes(x = dist, y = cor, fill = type), outlier.shape = NA) + 
  xlab('Motif position') + ylab('abs(correlaiton)') + scale_fill_discrete('') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.position = c(0.12, 0.88), legend.text = element_text(size = 14))
ggsave(file.path(outdir, paste(outpref, 'cor_box_by_type.png', sep = '')), p3, width = 6.5, height = 5.6)