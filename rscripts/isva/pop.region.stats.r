rm(list=ls())
library(ggplot2)
library(Matrix)
library(GenomicRanges)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

get.enrich.corr = function(regions, isva.fit, dat, ov.deg, ov.non.deg, is.snp = T){
  enrich = array(NaN, dim = c(1, 2))
  if(is.snp){
    deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), snps.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), snps.to.ranges(dat), ignore.strand = T, select = 'first')
  }else{
    deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  }
  deg.sum = sum(!is.na(ov.deg) & !is.na(deg)) # Regions with variants (in our set of samples) AND high Fst SNPs
  non.deg.sum = sum(!is.na(ov.non.deg) & !is.na(non.deg))
  enrich[1, 1] = (deg.sum / ov.deg.sum) / (non.deg.sum / ov.non.deg.sum)
  enrich[1, 2] = binom.val(deg.sum, ov.deg.sum, non.deg.sum / ov.non.deg.sum)
  return(enrich)
}

get.enrich = function(regions, isva.fit, dat, is.snp = T){
  enrich = array(NaN, dim = c(1, 2))
  if(is.snp){
    deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), snps.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), snps.to.ranges(dat), ignore.strand = T, select = 'first')
  }else{
    deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  }
  deg.sum = sum(!is.na(deg))
  non.deg.sum = sum(!is.na(non.deg))
  enrich[1, 1] = (deg.sum / isva.fit$ndeg) / (non.deg.sum / (nrow(regions) - isva.fit$ndeg))
  enrich[1, 2] = binom.val(deg.sum, isva.fit$ndeg, non.deg.sum / (nrow(regions) - isva.fit$ndeg))
  return(enrich)
}
get.wilc.enrich = function(regions, isva.fit, dat){
  deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
  non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  values.deg = dat[deg[!is.na(deg)], 4]
  values.non.deg = dat[non.deg[!is.na(non.deg)], 4]
  enrich = c(mean(values.deg) / mean(values.non.deg), wilcox.test(values.deg, values.non.deg)$p.value)
  return(enrich)
}
isva.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig_withSan/rdata/' #'../../rawdata/segSignal/14indiv/extractSignal/fc/avgSig/rdata'
plotdir = file.path(isva.dir, '..', 'plots')
files = list.files(isva.dir, pattern = 'all_reg_.*H3K4ME1.*_qn_isvaNull.RData', full.names = T)
isva.dir = '../../rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig/rdata/'
#files = append(files, list.files(isva.dir, pattern = 'hg19_w10k_all_reg.*_H3K27ME3.*_qn_isvaNull.RData', full.names = T))
isva.dir = '../../rawdata/geneCounts/rdata/repsComb/rdata/'
#files = append(files, list.files(isva.dir, pattern = 'all_reg.*_RZ.*_qn_isvaNull.RData', full.names = T))
isva.dir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig/rdata/'
#files = append(files, list.files(isva.dir, pattern = 'gencode.v13.annotation.noM.genes_all_reg.*_H3K36ME3.*_qn_isvaNull.RData', full.names = T))
isva.dir = '../../rawdata/signal/combrep/extractSignal/rand/fc/avgSig_withSan/rdata/'
files = append(files, list.files(isva.dir, pattern = 'all_reg.*H3K4ME1.*_qn_isvaNull.RData', full.names = T))
outpref = 'all_reg_H3K4ME1_rand_'
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')

#files = files[!grepl('rand|th0', files)]

# Load SNPs and genotype PCA
load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
#gp = new.env()
#load('../../rawdata/variants/all/snps/allNonSan/genot_pca.RData', gp)
#snp.pos = snp.pos[gp$pca.rows, 1:2] # SNPs that are variant

# Load gene metadata
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

ihs.dat = read.table('../../rawdata/selection/iHS_hg19.bed', header = F, sep = '\t')
colnames(ihs.dat) = c('chr', 'start', 'end', 'value')
#ihs.dat$pos = ihs.dat$end

# Read SNPs with high Fst
fst.dat = read.table('../../rawdata/variants/fst/FstSNP-HapMap3_selPop.txt')[, 1:3]
colnames(fst.dat) = c('chr', 'pos', 'fst')

# Read eQTLs
eqtl.dat = read.table('../../rawdata/QTL/eQTLs_uniq.txt')[, 1:2]
colnames(eqtl.dat) = c('chr', 'pos')

marks = gsub('SNYDER_HG19_|hg19_w.*k_|allEnhStates_|.*all_reg_|_qn_.*svaNull.RData', '', basename(files))
nmarks = length(marks)
tot.len = array(0, dim = c(nmarks, 1))
npeaks = array(0, dim = c(nmarks, 2))
snp.enrich = array(1, dim = c(nmarks, 2))
snp.enrich[, 1] = 0
fenrich = array(1, dim = c(nmarks, 2))
fenrich[, 1] = 0
qtl.enrich = array(1, dim = c(nmarks, 2))
qtl.enrich[, 1] = 0
ihs.enrich = array(1, dim = c(nmarks, 2))
ihs.enrich[, 1] = 0
signal.wilc = array(1, dim = c(nmarks, 2))
signal.wilc[, 1] = 0

penrich = NULL
all.deg.snps = NULL
all.non.deg.snps = NULL
for(i in 1:length(files)){
  load(gsub('_isvaNull|_rand_pop[0-9]*', '', files[i]))
  orig.counts = counts[good.rows, ]
  npeaks[i, 2] = length(good.rows) # Total number of regions, including those that were removed before ISVA
  load(files[i]) # Load the file with the ISVA results. This only has the regions used for ISVA.
  #isva.fit = sva.fit
  npeaks[i, 1] = isva.fit$ndeg 
  if(grepl('RZ|RNA', files[i])){
    regions = gene.meta[good.rows, 2:4]
  }
  tot.len[i] = sum((regions$end - regions$start + 1)[isva.fit$deg])
  
  if(isva.fit$ndeg > 1){
    #signal.wilc[i, 1] = median(rowMeans(orig.counts[isva.fit$deg, ])) / median(rowMeans(orig.counts[-isva.fit$deg, ]))
    #signal.wilc[i, 2] = wilcox.test(rowMeans(orig.counts[isva.fit$deg, ]), rowMeans(orig.counts[-isva.fit$deg, ]))$p.value
    
    ######## Are selected regions enriched in variants?
    # Count how many regions overlap at least one variant
    ov.deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), snps.to.ranges(snp.pos), select = 'first', ignore.strand = T)
    ov.deg.sum = sum(!is.na(ov.deg))
    ov.non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), snps.to.ranges(snp.pos), select = 'first', ignore.strand = T)
    ov.non.deg.sum = sum(!is.na(ov.non.deg))
    snp.enrich[i, 1] = (ov.deg.sum / isva.fit$ndeg) / (ov.non.deg.sum / (nrow(regions) - isva.fit$ndeg))
    snp.enrich[i, 2] = binom.val(ov.deg.sum, isva.fit$ndeg, ov.non.deg.sum / (nrow(regions) - isva.fit$ndeg))
    
    #ihs.enrich[i, ] = get.wilc.enrich(regions, isva.fit, ihs.dat) #get.enrich(regions, isva.fit, ihs.dat, is.snp = F)
    
    ######## Are the selected regions enriched in high-Fst variants?
    fenrich[i, ] = get.enrich.corr(regions, isva.fit, fst.dat, ov.deg, ov.non.deg, is.snp = T)
    
    ######## Are the selected regions enriched in eQTLs?    
    qtl.enrich[i, ] = get.enrich.corr(regions, isva.fit, eqtl.dat, ov.deg, ov.non.deg, is.snp = T)
    
    ######## Do the selected SNPs tend to have larger absolute loadings in the genotype PCA?
#     deg.snps = unique(ov.deg[!is.na(ov.deg)])
#     non.deg.snps = unique(ov.non.deg[!is.na(ov.non.deg)])
#     all.deg.snps = unique(append(all.deg.snps, deg.snps))
#     all.non.deg.snps = unique(append(all.non.deg.snps, non.deg.snps))
#     for(p in 1:12){
#       tmp.dat = tmp.dat = data.frame(mark = marks[i], pc = p,
#                                      enrich = median(abs(gp$pca.fit$rotation[deg.snps, p])) / median(abs(gp$pca.fit$rotation[non.deg.snps, p])),
#                                      pval = binom.test(sum(abs(gp$pca.fit$rotation[deg.snps, p]) > 0.0005), length(deg.snps),
#                                                        sum(abs(gp$pca.fit$rotation[non.deg.snps, p]) > 0.0005) / length(non.deg.snps))$p.value)
#                                      #pval = wilcox.test(abs(gp$pca.fit$rotation[deg.snps, p]), abs(gp$pca.fit$rotation[non.deg.snps, p]))$p.value)
#       penrich = rbind(penrich, tmp.dat)
#     }
  }
}

# for(p in 1:12){
#   tmp.dat = tmp.dat = data.frame(mark = 'all', pc = p,
#                                  enrich = median(abs(gp$pca.fit$rotation[all.deg.snps, p])) / median(abs(gp$pca.fit$rotation[all.non.deg.snps, p])),
#                                  pval = binom.val(sum(abs(gp$pca.fit$rotation[all.deg.snps, p]) > 0.0005), length(all.deg.snps),
#                                                    sum(abs(gp$pca.fit$rotation[all.non.deg.snps, p]) > 0.0005) / length(all.non.deg.snps)))
#   #pval = wilcox.test(abs(gp$pca.fit$rotation[deg.snps, p]), abs(gp$pca.fit$rotation[non.deg.snps, p]))$p.value)
#   penrich = rbind(penrich, tmp.dat)
# }

# genot.norm = scale(gp$genot)
# indivs = colnames(gp$genot)
# pca.fit = prcomp(t(genot.norm[all.deg.snps, !(colnames(gp$genot) %in% c('GM12878', 'SNYDER'))]), center = F, scale = F)
# p=plot.pcs(t(genot.norm[all.deg.snps, ]) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = indivs, groups = get.pop(indivs), all = F, ndim = 10)
# 
# pca.fit.n = prcomp(t(genot.norm[all.non.deg.snps, !(colnames(gp$genot) %in% c('GM12878', 'SNYDER'))]), center = F, scale = F)
# p=plot.pcs(t(genot.norm[all.non.deg.snps, ]) %*% pca.fit.n$rotation,  pca.fit.n$rotation, pca.fit.n$sdev, labels = indivs, groups = get.pop(indivs), all = F, ndim = 10)

#frac.peaks = npeaks[, 1] / npeaks[, 2]
sel.marks = npeaks[, 2] > -1 & npeaks[, 1] > -1
marks = gsub('_H3K4ME1', '', marks[sel.marks])
#marks = order.marks(marks)
npeaks = npeaks[sel.marks, ]
tot.len = tot.len[sel.marks]
fenrich = fenrich[sel.marks, ]
snp.enrich = snp.enrich[sel.marks, ]
qtl.enrich = qtl.enrich[sel.marks, ]
len.dat = data.frame(mark = marks, len = tot.len / 1e6, tot = npeaks[,1], npeaks = npeaks[, 1]/npeaks[, 2], f = fenrich[, 1], fp = fenrich[, 2], 
                     p = snp.enrich[, 1], pp = snp.enrich[, 2], q = qtl.enrich[, 1], qp = qtl.enrich[, 2])
#save(len.dat, file = file.path(plotdir, paste(outpref, 'enrich.RData', sep = '')))
write.table(len.dat, file = file.path(plotdir, paste(outpref, 'enrich.txt', sep = '')), col.names = T, row.names = F, quote = F, sep = '\t')
p1 = ggplot(len.dat) + geom_bar(aes(x = mark, y = len, fill = mark)) + xlab('') + 
  ylab('Length of population specific regions (Mb)') + theme_bw() +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = -68, hjust = 0, vjust = 1))
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_len.pdf', sep = '')), p1, width = 6.5, height = 5.6)

p2 = ggplot(len.dat) + geom_bar(aes(x = mark, y = npeaks, fill = mark)) + 
  #annotate('text', x = len.dat$mark, y = len.dat$npeaks + 0.01, label = len.dat$tot, size = 5) + 
  xlab('') + ylab('Fraction of population specific regions') + theme_bw() +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = -68, hjust = 0, vjust = 1))
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_frac.pdf', sep = '')), p2, width = 6.5, height = 5.6)

annot = array('', dim = c(nrow(len.dat), 1))
annot[fenrich[, 2] < 0.1 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$f
#annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p3 = ggplot(len.dat) + geom_bar(aes(x = mark, y = f, fill = mark)) + xlab('') + annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in SNPs with high Fst') + theme_bw() + coord_cartesian(ylim = c(0.5, max(annot.y) + 0.25)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = -68, hjust = 0, vjust = 1))
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_fst.pdf', sep = '')), p3, width = 6.5, height = 5.6)

annot = array('', dim = c(nrow(len.dat), 1))
annot[snp.enrich[, 2] < 0.1 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$p
#annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p4 = ggplot(len.dat) + geom_bar(aes(x = mark, y = p, fill = mark)) + xlab('') + annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in SNPs') + theme_bw() + ylim(c(0, 1.2)) + coord_cartesian(ylim = c(0.75, max(annot.y) + 0.05)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = -68, hjust = 0, vjust = 1))
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_snp_enrich.pdf', sep = '')), p4, width = 6.5, height = 5.6)

annot = array('', dim = c(nrow(len.dat), 1))
annot[qtl.enrich[, 2] < 0.1 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$q
#annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p5 = ggplot(len.dat) + geom_bar(aes(x = mark, y = q, fill = mark)) + xlab('') + annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in eQTLs') + theme_bw() +  coord_cartesian(ylim = c(0.5, max(annot.y) + 0.25)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) +
  theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = -68, hjust = 0, vjust = 1))
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_eQTL_enrich.pdf', sep = '')), p5, width = 6.5, height = 5.6)
