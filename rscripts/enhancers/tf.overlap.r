rm(list=ls())
library(GenomicRanges)
library(matrixStats)
source('utils/binom.val.r')
source('utils/deseq.utils.r')

# Statistics for overlaps between TF binding sites and sets of associations (or any set of regions)
tf.ov.stats = function(tf.dat, true.regions, rand.regions, is.pos = F){  
  # Overlaps between positive set and TFBSs
  if(is.pos){
    true.ov = findOverlaps(regions.to.ranges(true.regions), snps.to.ranges(tf.dat), select = 'all', ignore.strand = T)
    rand.ov = findOverlaps(regions.to.ranges(rand.regions), snps.to.ranges(tf.dat), select = 'all', ignore.strand = T)
  }else{
    true.ov = findOverlaps(regions.to.ranges(true.regions), regions.to.ranges(tf.dat), select = 'all', ignore.strand = T)
    rand.ov = findOverlaps(regions.to.ranges(rand.regions), regions.to.ranges(tf.dat), select = 'all', ignore.strand = T)
  }  
  true.hits = table(tf.dat$tf[subjectHits(true.ov)]) # Count overlaps with each TF
  rand.hits = table(tf.dat$tf[subjectHits(rand.ov)])
  
  tfs = names(true.hits) # TFs appearing in the positive set
  ratios = array(1, dim = c(length(tfs), 1))  
  pvals = array(1, dim = c(length(tfs), 1))
  for(i in 1:length(tfs)){
    if(true.hits[i] > 20){
      idx = names(rand.hits) == tfs[i]
      if(any(idx)){
        r = rand.hits[idx]
      }else{r = 0}
      pvals[i] = binom.val(true.hits[i], nrow(true.regions), r / nrow(rand.regions), alternative = 'greater')
      ratios[i] = (true.hits[i] / nrow(true.regions)) / (r / nrow(rand.regions))
    }
  }
  return(data.frame(tf = tfs, ratio = ratios, pval = pvals))
}

link.file = '../../rawdata/enhancers/merged_Mar13/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData'
load(link.file)
outpref = gsub('.RData', '', basename(link.file))
plotdir = '../../rawdata/enhancers/merged_Mar13/plots/'

indiv.idx = match('GM12878', colnames(ars.score.norm))
bg.indiv.idx = match(c('GM12878', 'HG2255', 'HG2588', 'HG2610', 'HG2630'), colnames(ars.score.norm))
#sel.pos = unique(links$region.idx[links$max.line == indiv.idx])
sel.pos = unique(coef.dat$region.idx[ars.score.norm[, indiv.idx] > 0.5 & coef.dat$ars.max > min(links$coef) / 4])
#sel.pos = unique(coef.dat$region.idx[coef.dat$max.line == 2][1:1000])
true.regions = ac.regions[sel.pos, ] # Take 1000 regions with highest ARS
true.regions$start = true.regions$start - 200
true.regions$end = true.regions$end + 200

sel.neg = unique(coef.dat$region.idx[ars.score.norm[, indiv.idx] < 0.3 & coef.dat$ars.max > min(links$coef)]) #[1:min(1000, sum(ars.score.norm[, indiv.idx] < 0.3))])
#sel.neg = unique(coef.dat$region.idx[rowMaxs(ars.score.norm[, bg.indiv.idx]) < 0.3 & coef.dat$ars.max > min(links$coef) / 4])#[1:min(1000, sum(ars.score.norm[, indiv.idx] < 0.3))]
#sel.neg = unique(coef.dat$region.idx[coef.dat$max.line == 6][1:1000])
rand.regions = ac.regions[sel.neg, ]
rand.regions$start = rand.regions$start - 200
rand.regions$end = rand.regions$end + 200
print(wilcox.test(rowMaxs(ac.counts[sel.pos, ]), rowMaxs(ac.counts[sel.neg, ]))) # Is there a significant difference in maximum signal?

tf.dat = read.table('../../rawdata/TFs/Gm12878_allTFBS.sorted.noPol.bed', header = F, sep = '\t')
tf.dat[, 4] = gsub('_Rank_[0-9]+|Iggmus.*|Std.*|Pcr.*|Iggrab.*|ucla|Haib|Sydh|V04*', '', tf.dat[, 4])
colnames(tf.dat) = c('chr', 'start', 'end', 'tf')

stats = tf.ov.stats(tf.dat, true.regions, rand.regions)
write.table(stats, quote = F, sep = '\t', col.names = T, row.names = F, file = file.path(plotdir, paste(outpref, '_tfEnrich.txt', sep = '')))
stats.cp = stats[order(stats$ratio, decreasing = T), ]
stats.cp$tf = ordered(stats.cp$tf, levels = stats.cp$tf)
stats.cp$ratio = log2(stats.cp$ratio)
#stats.cp = stats[stats$ratio > 1, ]

annot = array('', dim = c(nrow(stats.cp), 1))
annot[stats.cp$pval < 0.001] = '*'
annot.y = stats.cp$ratio
p = ggplot(stats.cp) + geom_bar(aes(x = tf, y = ratio), fill = 'darkorchid4', stat = 'identity', position = 'dodge') +
  annotate('text', x = stats.cp$tf, y = annot.y, label = annot, size = 5) +
  xlab('Transcription factor') + ylab('Enrichment') + theme_bw() + 
  theme(axis.text.x = element_text(size = 6, angle = -90, hjust = 0, vjust = 0.5), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file.path(plotdir, paste(outpref, '_tfEnrich_sortRatio.pdf', sep = '')), p, width = 6.5, height = 5.6)

# indel.dir = '../../rawdata/variants/all/masks/'
# indel.stats = NULL
# for(i in 1:ncol(ars.score.norm)){
#   sel.pos = unique(coef.dat$region.idx[ars.score.norm[, i] > 0.7][1:1000])
#   true.regions = ac.regions[sel.pos, ]
#   true.regions$start = true.regions$start - 200
#   true.regions$end = true.regions$end + 200
#   
#   sel.neg = unique(coef.dat$region.idx[ars.score.norm[, i] < 0.3][1:1000])
#   rand.regions = ac.regions[sel.neg, ]
#   rand.regions$start = rand.regions$start - 200
#   rand.regions$end = rand.regions$end + 200
#   
#   indels = read.bed(file.path(indel.dir, paste(colnames(ars.score.norm)[i], '.blacklist.bed', sep = '')))
#   indels$tf = rep('Indel', nrow(indels))
#   tmp.stats = tf.ov.stats(indels, true.regions, rand.regions, is.pos = F) 
#   tmp.stats$indiv = colnames(ars.score.norm)[i]
#   indel.stats = rbind(indel.stats, tmp.stats)
# }

# load('../../rawdata/variants/all/snps/allNonSan/all_genot.RData')
# load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
# snp.pos$tf = factor(rep('SNP', nrow(snp.pos)))
# 
# snp.stats = NULL
# for(i in 1:ncol(genot)){
#   indiv.idx = match(colnames(genot)[i], colnames(ars.score.norm))
#   sel.pos = unique(coef.dat$region.idx[ars.score.norm[, indiv.idx] > 0.7][1:1000])
#   true.regions = ac.regions[sel.pos, ]
#   true.regions$start = true.regions$start - 200
#   true.regions$end = true.regions$end + 200
#   
#   sel.neg = unique(coef.dat$region.idx[ars.score.norm[, indiv.idx] < 0.3][1:1000])
#   rand.regions = ac.regions[sel.neg, ]
#   rand.regions$start = rand.regions$start - 200
#   rand.regions$end = rand.regions$end + 200
#   
#   sel = (genot[, i] == 2 & rowSums(genot == 2) == 1) | ((genot[, i] == 0 & rowSums(genot == 0) == 1)) 
#   tmp.stats = tf.ov.stats(snp.pos[sel, ], true.regions, rand.regions, is.pos = T) 
#   tmp.stats$indiv = colnames(genot)[i]
#   snp.stats = rbind(snp.stats, tmp.stats)
# }