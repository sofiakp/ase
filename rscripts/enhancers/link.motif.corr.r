rm(list=ls())
library(GenomicRanges)
library(ggplot2)
library(foreach)
library(doMC)
source('utils/deseq.utils.r')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

nchunks = 5
registerDoMC(nchunks)

# Load links
link.file = '../../rawdata/enhancers/rdata/enhancer_coef_elastic0.5_100kb_asinh0.2_cv0.2_H3K27AC_links.RData'
load(link.file)
plotdir = '../../rawdata/enhancers/plots/'
pref = gsub('.RData', '', basename(link.file))
indivs = colnames(rna.counts)
nindivs = length(indivs)
# Extend peaks to left and right to account for missing dips
ac.regions$start = ac.regions$start - 200
ac.regions$end = ac.regions$end + 200
nregions = length(unique(links$region.idx)) # unique regions involved in associations

# Read SNP information
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/') # Where to read genotype information from
big.meta = new.env()
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file, big.meta)

ov = findOverlaps(regions.to.ranges(ac.regions[links$region.idx, ]), snps.to.ranges(big.meta$snp.pos), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov)) # Col 1 - link idx, Col 2 - snp idx
rna.counts.ov = rna.counts[links$gene.idx[ov.mat[, 1]], ]
# ac.counts.ov = ac.counts[ov.mat[, 1], ]

nsnps = nrow(ov.mat)
genot = array(0, dim = c(nsnps, nindivs))
for(i in 1:nindivs){
  load(file.path(geno.dir, paste(indivs[i], '.snps.RData', sep = '')))
  genot[, i] = as.integer(as.vector(geno.info$mat[ov.mat[, 2]])) + as.integer(as.vector(geno.info$pat[ov.mat[, 2]]))
}
colnames(genot) = indivs
sel.genot = rowSums(genot == 0) < nindivs & rowSums(genot == 1) < nindivs & rowSums(genot == 2) < nindivs # Select positions where individuals differ
ov.mat = ov.mat[sel.genot, ]
genot = genot[sel.genot, ]
rna.counts.ov = rna.counts.ov[sel.genot, ]
cat('Number of distal elements', nregions, '\n')
hits = length(unique(links$region.idx[ov.mat[, 1]]))
cat('Distal elements with variants', hits, '( ', hits * 100 / nregions, '% )\n')

motif.hits = read.table('../../rawdata/motifs/personal_genome_motifs.txt')
motif.hits = motif.hits[, c(1,4,11)]
colnames(motif.hits) = c('chr', 'pos', 'name')
# Get the SNPs that overlapped the distal elements and overlap them with motifs
mov = findOverlaps(snps.to.ranges(big.meta$snp.pos[ov.mat[, 2], ]), snps.to.ranges(motif.hits), select = 'first', ignore.strand = T)
hits = length(unique(links$region.idx[ov.mat[!is.na(mov), 1]]))
cat('Distal elements with variants in motifs', hits, '( ', hits  * 100 / nregions, '% )\n')

snp.cor = cor.par(genot, rna.counts.ov, 5, 'spearman')
# cor.dat.tmp = data.frame(cor = snp.cor, region.idx = links$region.idx[ov.mat[, 1]], gene.idx = links$gene.idx[ov.mat[, 1]])
# cor.dat.tmp = cast(cor.dat.tmp, region.idx+gene.idx~., function(x) max(abs(x)) * sign(which.max(abs(x))), value = 'cor')
# cor.dat = data.frame(cor = cor.dat.tmp[, 3], type = rep('true', nrow(cor.dat.tmp)))
# rand.cor = cor.par(genot, rna.counts[sample(1:nrow(rna.counts), nrow(genot)), ], nchunks, 'spearman')
# cor.dat.tmp = data.frame(cor = rand.cor, region.idx = links$region.idx[ov.mat[, 1]], gene.idx = links$gene.idx[ov.mat[, 1]])
# cor.dat.tmp = cast(cor.dat.tmp, region.idx+gene.idx~., function(x) max(abs(x)) * sign(which.max(abs(x))), value = 'cor')
# cor.dat = rbind(cor.dat, data.frame(cor = cor.dat.tmp[, 3], type = rep('random', nrow(cor.dat.tmp))))
# 
# p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type)) +
#   xlab('Correlation') + ylab('density') + theme_bw() + 
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
#         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
#         legend.position = c(.1, .8))
#ggsave(file = file.path(plotdir, paste(pref, '_snpGeneCorr.png', sep = '')), p, width= 6.5, height = 5.6)

# Write correlated SNPs and genes to text file
output.dat = ac.regions[links$region.idx[ov.mat[, 1]], ]
output.dat$gene.name = gene.meta$gene.name[links$gene.idx[ov.mat[, 1]]]
output.dat$pos = big.meta$snp.pos$pos[ov.mat[, 2]]
output.dat$ref = big.meta$snp.pos$ref[ov.mat[, 2]]
output.dat = cbind(output.dat, genot)
mot = array('-', dim = c(length(mov), 1))
mot[!is.na(mov)] = as.character(motif.hits$name[mov[!is.na(mov)]])
output.dat$mot = mot
output.dat$cor = snp.cor
output.dat = output.dat[order(abs(output.dat$cor), decreasing = T), ]
write.table(output.dat, quote = F, col.names = T, row.names = F, sep = '\t', file = file.path(plotdir, paste(pref, '_corr.txt', sep = '')))