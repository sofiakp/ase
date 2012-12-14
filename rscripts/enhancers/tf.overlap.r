rm(list=ls())
library(GenomicRanges)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# Statistics for overlaps between TF binding sites and sets of associations (or any set of regions)
tf.ov.stats = function(tf.dat, true.regions, rand.regions){  
  # Overlaps between positive set and TFBSs
  true.ov = findOverlaps(regions.to.ranges(true.regions), regions.to.ranges(tf.dat), select = 'all', ignore.strand = T)
  true.hits = table(tf.dat$tf[subjectHits(true.ov)]) # Count overlaps with each TF
  rand.ov = findOverlaps(regions.to.ranges(rand.regions), regions.to.ranges(tf.dat), select = 'all', ignore.strand = T)
  rand.hits = table(tf.dat$tf[subjectHits(rand.ov)])
  
  tfs = names(true.hits)[true.hits > 20] # TFs appearing in the positive set
  ratios = array(1, dim = c(length(tfs), 1))  
  pvals = array(1, dim = c(length(tfs), 1))
  for(i in 1:length(tfs)){
    if(true.hits[i] > 20){
      idx = names(rand.hits) == tfs[i]
      if(any(idx)){
        r = rand.hits[idx]
      }else{r = 0}
      pvals[i] = binom.test(true.hits[i], nrow(true.regions), r / nrow(rand.regions))$p.value
      ratios[i] = (true.hits[i] / nrow(true.regions)) / (r / nrow(rand.regions))
    }
  }
  return(data.frame(tf = tfs, ratio = ratios, pval = pvals))
}

tf.dat = read.table('../../rawdata/Gm12878_allTFBS.sorted.noPol.bed', header = F, sep = '\t')
tf.dat[, 4] = gsub('_Rank_[0-9]+', '', tf.dat[, 4])
colnames(tf.dat) = c('chr', 'start', 'end', 'tf')

load('../../rawdata/enhancers/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData')
sel.pos = unique(coef.dat$region.idx[ars.score.norm[, 6] > 0.7][1:1000])
#sel.pos = unique(coef.dat$region.idx[coef.dat$max.line == 2][1:1000])
true.regions = ac.regions[sel.pos, ] # Take 1000 regions with highest ARS
true.regions$start = true.regions$start - 200
true.regions$end = true.regions$end + 200
sel.neg = unique(coef.dat$region.idx[ars.score.norm[, 6] < 0.3][1:1000])
#sel.neg = unique(coef.dat$region.idx[coef.dat$max.line == 6][1:1000])
rand.regions = ac.regions[sel.neg, ]
rand.regions$start = rand.regions$start - 200
rand.regions$end = rand.regions$end + 200
print(wilcox.test(rowMaxs(ac.counts[sel.pos, ]), rowMaxs(ac.counts[sel.neg, ]))) # Is there a significant difference in maximum signal?

stats = tf.ov.stats(tf.dat, true.regions, rand.regions)
print(stats[stats$pval < 0.001, ])
