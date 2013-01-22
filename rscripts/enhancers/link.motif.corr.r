rm(list=ls())
library(reshape)
library(GenomicRanges)
library(ggplot2)
library(foreach)
library(doMC)
source('utils/deseq.utils.r')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

# Computes correlations between signal and genotype at linked regions and enrichments in SNPs and motifs

set.seed(1)

nchunks = 5
registerDoMC(nchunks)

# Load links
link.file = '../../rawdata/enhancers/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData'
load(link.file)
plotdir = '../../rawdata/enhancers/plots/'
pref = gsub('.RData', '', basename(link.file))
indivs = colnames(rna.counts)
nindivs = length(indivs)
# Extend peaks to left and right to account for missing dips
ac.regions$start = ac.regions$start - 200
ac.regions$end = ac.regions$end + 200
nregions = length(unique(links$region.idx)) # unique regions involved in associations

# Get some random non-linked regions for enrichments
rand.idx = sample(setdiff(1:nrow(ac.regions), unique(links$region.idx)), 10 * nregions)
rand.regions = ac.regions[rand.idx, ]
rand.counts = ac.counts[rand.idx, ]

ac.regions = ac.regions[links$region.idx, ]
ac.counts = ac.counts[links$region.idx, ]
rna.counts = rna.counts[links$gene.idx, ]

# Read SNP information
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

# Read indel information
indels = read.bed('../../rawdata/variants/all/masks/all.blacklist.bed')

# Load genotypes
load('../../rawdata/variants/all/snps/allNonSan/all_genot.RData')
genot = genot[, match(indivs, colnames(genot))]
sel.genot = rowSums(genot == 0) < nindivs & rowSums(genot == 1) < nindivs & rowSums(genot == 2) < nindivs # Select positions where individuals differ
genot = genot[sel.genot, ]
snp.pos = snp.pos[sel.genot, ]

# Read motif information
motif.hits = read.table('../../rawdata/motifs/personal_genome_motifs.txt')
motif.hits = motif.hits[, c(1,4,11)]
colnames(motif.hits) = c('chr', 'pos', 'name')

cat('Number of distal elements', nregions, '\n')

# Are regions enriched in indels?
ov = findOverlaps(regions.to.ranges(ac.regions), regions.to.ranges(indels), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov)) # Col 1 - linked region idx, Col 2 - snp idx
ov.rand = findOverlaps(regions.to.ranges(rand.regions), regions.to.ranges(indels), select = 'all', ignore.strand = T)
hits = length(unique(links$region.idx[ov.mat[, 1]])) # Number of unique regions with indels
rand.hits = length(unique(queryHits(ov.rand)))
indel.enrich = c((hits / nregions) / (rand.hits / nrow(rand.regions)), 
                 binom.val(hits, nregions, rand.hits / nrow(rand.regions)))
cat('Distal elements with indels', hits, '( ', hits * 100 / nregions, '% )\n')
print(indel.enrich)

# Are regions enriched in SNPs?
# Find linked regions that overlap SNPs and check for enrichment
ov = findOverlaps(regions.to.ranges(ac.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov)) # Col 1 - linked region idx, Col 2 - snp idx
ov.rand = findOverlaps(regions.to.ranges(rand.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)

reg.hits = length(unique(links$region.idx[ov.mat[, 1]])) # Number of unique regions with SNPs
reg.rand.hits = length(unique(queryHits(ov.rand)))
snp.enrich = c((reg.hits / nregions) / (reg.rand.hits / nrow(rand.regions)), 
               binom.val(reg.hits, nregions, reg.rand.hits / nrow(rand.regions)))
cat('Distal elements with variants', reg.hits, '( ', reg.hits * 100 / nregions, '% )\n')
print(snp.enrich)

# Are regions enriched in SNPs in motifs?
# Get the SNPs that overlapped the distal elements and overlap them with motifs
mov = findOverlaps(snps.to.ranges(snp.pos[ov.mat[, 2], ]), snps.to.ranges(motif.hits), select = 'first', ignore.strand = T)
mov.rand = findOverlaps(snps.to.ranges(snp.pos[subjectHits(ov.rand), ]), snps.to.ranges(motif.hits), select = 'first', ignore.strand = T)
hits = length(unique(links$region.idx[ov.mat[!is.na(mov), 1]])) # Number of unique regions with a SNP in a motif
rand.hits = length(unique(queryHits(ov.rand)[!is.na(mov.rand)]))
motif.enrich = c((hits / reg.hits) / (rand.hits / reg.rand.hits), 
                 binom.val(hits, reg.hits, rand.hits / reg.rand.hits))
cat('Distal elements with variants in motifs', hits, '( ', hits  * 100 / nregions, '% )\n')
print(motif.enrich)

enrich.dat = data.frame(p = c(indel.enrich[2], snp.enrich[2], motif.enrich[2]), e = c(indel.enrich[1], snp.enrich[1], motif.enrich[1]), annot = factor(c('Indels', 'SNPs', 'Motifs')))
p1 = ggplot(enrich.dat) + geom_bar(aes(x = annot, y = e)) + xlab('') + 
  annotate('text', x = enrich.dat$annot, y = enrich.dat$e + 0.05, label = sprintf('%.4f', enrich.dat$p), size = 6) + 
  ylab('Enrichments of linked regions') + theme_bw() +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = -68, hjust = 0, vjust = 1))
ggsave(file = file.path(plotdir, paste(pref, '_snpMotEnrich.png', sep = '')), p1, width= 6.5, height = 5.6)

# Correlation between signal at regions and genotype of all overlapping SNPs
snp.cor = cor.par(genot[ov.mat[, 2], ], ac.counts[ov.mat[, 1], ], 5, 'spearman')
rand.cor = cor.par(genot[subjectHits(ov.rand), ], rand.counts[queryHits(ov.rand), ], nchunks, 'spearman')

# Write correlated SNPs and genes to text file
output.dat = ac.regions[ov.mat[, 1], ]
output.dat$gene.name = gene.meta$gene.name[links$gene.idx[ov.mat[, 1]]]
output.dat$pos = snp.pos$pos[ov.mat[, 2]]
output.dat$ref = snp.pos$ref[ov.mat[, 2]]
output.dat = cbind(output.dat, as.matrix(genot[ov.mat[, 2], ]))
mot = array('-', dim = c(length(mov), 1))
mot[!is.na(mov)] = as.character(motif.hits$name[mov[!is.na(mov)]])
output.dat$mot = mot
output.dat$cor = snp.cor
output.dat = output.dat[order(abs(output.dat$cor), decreasing = T), ]
#write.table(output.dat, quote = F, col.names = T, row.names = F, sep = '\t', file = file.path(plotdir, paste(pref, '_signal_corr.txt', sep = '')))

# For each region, get the SNP with the highest correlation. Then, test if the SNPs in the true regions have higher correlations than
# the SNPs in the control regions.
cor.dat.tmp = data.frame(cor = snp.cor, region.idx = links$region.idx[ov.mat[, 1]],
                         chr = snp.pos$chr[ov.mat[, 2]], pos = snp.pos$pos[ov.mat[, 2]])
cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x)) * sign(x[which.max(abs(x))]), value = 'cor')
colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
cor.dat.true = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor'))
write.table(cor.dat.true[, 3:4], row.names = F, col.names = F, quote = F, sep = '\t', file = file.path(plotdir, paste(pref, '_sign_bestSNP.txt', sep = '')))

cor.dat.tmp = data.frame(cor = rand.cor, region.idx = rand.idx[queryHits(ov.rand)],
                         chr = snp.pos$chr[subjectHits(ov.rand)], pos = snp.pos$pos[subjectHits(ov.rand)])
cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x)) * sign(x[which.max(abs(x))]), value = 'cor')
colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
cor.dat.perm = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor'))
cor.dat = rbind(data.frame(cor = abs(cor.dat.true[, 2]), type = rep('true', nrow(cor.dat.true))), 
                data.frame(cor = abs(cor.dat.perm[, 2]), type = rep('random', nrow(cor.dat.perm))))

wilc.p = wilcox.test(abs(cor.dat.true$cor), abs(cor.dat.perm$cor))$p.value
p = ggplot(cor.dat) + geom_density(aes(x = cor, color = type), adjust = 1.5) +
  annotate('text', 0.75, 2, label = sprintf('Wilcoxon P = %.4f', wilc.p)) +
  xlab('Correlation between genotype and signal') + ylab('Density') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.position = c(.15, .8), legend.title = element_blank(), legend.text = element_text(size = 14))
ggsave(file = file.path(plotdir, paste(pref, '_snpGeneCorr.png', sep = '')), p, width= 6.5, height = 5.6)

# Overlap the "best" motif in each region with motifs
mov = findOverlaps(snps.to.ranges(cor.dat.true), snps.to.ranges(motif.hits), select = 'first', ignore.strand = T)
mov.rand = findOverlaps(snps.to.ranges(cor.dat.perm), snps.to.ranges(motif.hits), select = 'first', ignore.strand = T)
uniq.mot = unique(motif.hits$name)
mot.enrich = array(NaN, dim = c(length(uniq.mot), 2))
for(m in 1:length(uniq.mot)){
  hits = sum(motif.hits$name[mov[!is.na(mov)]] == uniq.mot[m])
  if(hits > 10){
    rand.hits = sum(motif.hits$name[mov.rand[!is.na(mov.rand)]] == uniq.mot[m])
    mot.enrich[m, ] = c((hits / nrow(cor.dat.true)) / (rand.hits / nrow(cor.dat.perm)), 
                        binom.val(hits, nrow(cor.dat.true), rand.hits / nrow(cor.dat.perm))) 
  }
}
rownames(mot.enrich) = uniq.mot

# Correlation with expression of associated genes
# Write correlated SNPs and genes to text file
# output.dat = ac.regions[links$region.idx[ov.mat[, 1]], ]
# output.dat$gene.name = gene.meta$gene.name[links$gene.idx[ov.mat[, 1]]]
# output.dat$pos = big.meta$snp.pos$pos[ov.mat[, 2]]
# output.dat$ref = big.meta$snp.pos$ref[ov.mat[, 2]]
# output.dat = cbind(output.dat, genot)
# mot = array('-', dim = c(length(mov), 1))
# mot[!is.na(mov)] = as.character(motif.hits$name[mov[!is.na(mov)]])
# output.dat$mot = mot
# output.dat$cor = snp.cor
# output.dat = output.dat[order(abs(output.dat$cor), decreasing = T), ]
# write.table(output.dat, quote = F, col.names = T, row.names = F, sep = '\t', file = file.path(plotdir, paste(pref, '_corr.txt', sep = '')))