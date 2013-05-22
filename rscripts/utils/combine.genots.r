rm(list=ls())
library(Matrix)
library(ggplot2)
library(ape)
library(geiger)
source('utils/deseq.utils.r')

non.san = new.env()
load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData', non.san)
non.san.genot = new.env()
load('../../rawdata/variants/all/snps/allNonSan/all_genot.RData', non.san.genot)
san = new.env()
load('../../rawdata/variants/sanConsensus/snps/san.snps.RData', san)
san.genot = new.env()
load('../../rawdata/variants/sanConsensus/snps/all_genot.RData', san.genot)

# Overlaps between non-san and san SNPs
ov = findOverlaps(snps.to.ranges(non.san$snp.pos), snps.to.ranges(san$snp.pos), select = 'first', ignore.strand = T)

genot = non.san.genot$genot
genot.cols = colnames(genot)
# For each San individual, append it's non-San-specific SNPs to genot
for(i in 1:ncol(san.genot$genot)){
  cat(colnames(san.genot$genot)[i], '\n')
  new.genot = array(0, dim = c(length(ov), 1))
  new.genot[!is.na(ov)] = san.genot$genot[ov[!is.na(ov)], i]
  genot = cBind(genot, new.genot)
  genot.cols = append(genot.cols, colnames(san.genot$genot)[i])
}

# San-specific SNPs that haven't been added so far
ov = findOverlaps(snps.to.ranges(san$snp.pos), snps.to.ranges(non.san$snp.pos), select = 'first', ignore.strand = T)
non.ov = is.na(ov)
snp.pos = non.san$snp.pos
snp.pos = rbind(snp.pos, san$snp.pos[non.ov, ])

genot2 = san.genot$genot[non.ov, ]
genot2.cols = colnames(genot2)
for(i in 1:ncol(non.san.genot$genot)){
  cat(colnames(non.san.genot$genot)[i], '\n')
  new.genot = array(0, dim = c(sum(non.ov), 1)) # These are not present in the non-San, so just zeros
  genot2 = cBind(genot2, new.genot)
  genot2.cols = append(genot2.cols, colnames(non.san.genot$genot)[i])
}

genot = rBind(genot, genot2[, match(genot.cols, genot2.cols)])

# Finally, add GM19193
gm = new.env()
load('../../rawdata/variants/novelCalls/filtered/snps/gm19193.snps.RData', gm)
load('../../rawdata/variants/novelCalls/filtered/snps/GM19193.snps.RData')
gm.genot = geno.info$mat + geno.info$pat

ov = findOverlaps(snps.to.ranges(gm$snp.pos), snps.to.ranges(snp.pos), select = 'first', ignore.strand = T)
new.gm.genot = array(0, dim = c(nrow(snp.pos), 1))
new.gm.genot[ov[!is.na(ov)]] = gm.genot[!is.na(ov)]
genot = cBind(genot, new.gm.genot)
genot.cols = append(genot.cols, 'GM19193')
colnames(genot) = genot.cols

snp.pos = rbind(snp.pos, gm$snp.pos[is.na(ov), ])
tmp.genot = array(0, dim = c(sum(is.na(ov)), ncol(genot)))
tmp.genot[, ncol(genot)] = gm.genot[is.na(ov)]
genot = rBind(genot, tmp.genot)

save(genot, file = '../../rawdata/variants/all_Mar13/genot.RData')
save(snp.pos, file = '../../rawdata/variants/all_Mar13/snps.RData')

load('../../rawdata/variants/all_Mar13/genot.RData')
load('../../rawdata/variants/all_Mar13/snps.RData')
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')

# SNPs on exons
ov = findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(gene.meta), select = 'first', ignore.strand = T) 
sel = !is.na(ov)
genot.sample = genot[sel, ][sample(1:sum(sel), 100000), colnames(genot) != 'GM19193' & colnames(genot) != 'GM12890'] 
colnames(genot.sample) = fix.indiv.names(colnames(genot.sample))
nj.tree = nj(dist(t(genot.sample), method = 'manhattan'))
edges = nj.tree$edge
edge.len = as.integer(nj.tree$edge.length * 100 / max(nj.tree$edge.length))
sel.edges = edges[, 1] > ncol(genot) & edges[, 2] > ncol(genot) & edge.len > 10
edge.lab = array('', dim = c(nrow(edges), 1))
edge.lab[sel.edges] = edge.len[sel.edges]

pdf('../../rawdata/variants/all_Mar13/genot_pca_noGM19193_noGM12890_exons_nj.pdf')
plot(nj.tree, 'u', cex = 1, edge.width = 0.5, no.margin = T, lab4ut='axial', label.offset = 0.5, tip.col = get.pop.col(get.pop(colnames(genot.sample))))
edgelabels(edge.lab, frame = 'none', adj = c(1, 0.5), cex = 0.9)
dev.off()

# Select variable sites
sel = rowSums(genot.sample == 0) < ncol(genot.sample) - 3 & rowSums(genot.sample == 1) < ncol(genot.sample) - 3 & rowSums(genot.sample == 2) < ncol(genot.sample) - 3
#genot.sample[genot.sample == 2] = 1 # Trios have much fewer homozygous alternative calls. This causes biases
genot.norm = scale(genot.sample[sel, ])
colnames(genot.norm) = colnames(genot.sample)
pca.fit = prcomp(t(genot.norm[, !(colnames(genot.sample) %in% c('GM12878', 'GM19240'))]), center = F, scale = F)
p=plot.pcs(t(genot.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = array('', dim=c(ncol(genot.sample),1)), groups = get.pop(colnames(genot.sample)), all = F, ndim = 2)
ggsave('../../rawdata/variants/all_Mar13/genot_pca_noGM19193_noGM12890_exonsVariable_small.pdf', p$p1, width = 4, height = 3)
save(genot.norm, pca.fit, file = '../../rawdata/variants/all_Mar13/genot_pca_noGM19193_exonsVariable_pca.RData')
ggsave('../../rawdata/variants/all_Mar13/genot_pca_noGM19193_exonsVariable.pdf', p$p1, width = 13.6, height = 11.8)
ggsave('../../rawdata/variants/all_Mar13/genot_eigen_noGM19193_exonsVariable.pdf', p$p2, width = 6.5, height = 5.6)

