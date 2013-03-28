rm(list=ls())
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(matrixStats)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('utils/binom.val.r')

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
plotdir = indir
mark = 'H3K27AC'
overwrite = F

bg.af.file = '../../rawdata/variants/af/1000gp_af_1in100.txt'
bg.af = as.matrix(read.table(bg.af.file, sep = '\t', header = F)[, 4:6])

af.file = '../../rawdata/variants/af/1000gp_af.txt'
regions.file = file.path(indir, paste('SNYDER_HG19_all_reg_isva_sign_', mark, '.txt', sep = ''))
tmp.file = file.path(indir, paste('SNYDER_HG19_all_reg_isva_sign_', mark, '_af.txt', sep = ''))
if(!file.exists(tmp.file) || overwrite) system(paste("awk 'BEGIN{OFS=\"\t\"}{print $1,$2-1,$2,\".\",\"0\",\"+\",$3,$4,$5,$6}'", af.file, '| intersectBed -a stdin -b', regions.file, ' -wa >', tmp.file))
af = read.table(tmp.file)
af.pos = af[, c(1,3)]
colnames(af.pos) = c('chr', 'pos')
af = af[, 8:10]
pop = c('YRI', 'CEU', 'Asian')
# Read allele frequencies at SNPs overlapping positive and negative regions
#af = read.table(file.path(indir, 'alleleFreq', paste('SNYDER_HG19_all_reg_isva_sign_', mark, '_AF.txt', sep = '')), header= F, sep='\t')
#af.pos = af[, 1:2]
#colnames(af.pos) = c('chr', 'pos')
#af = matrix(as.numeric(unlist(strsplit(as.character(af[,3]), ','))), ncol = 3, byrow = T)[, c(3,2,1)] # Reorder the populations        
#af.neg = read.table(file.path(indir, 'alleleFreq', paste('SNYDER_HG19_all_reg_', mark, '_qn_isvaNull_non_sign.bed', sep = '')), header= F, sep='\t')
#af.neg.pos = af.neg[, 1:2]
#colnames(af.neg.pos) = c('chr', 'pos')
#af.neg = matrix(as.numeric(unlist(strsplit(as.character(af.neg[,3]), ','))), ncol = 3, byrow = T)[, c(3,2,1)]
#pop = c('CEU', 'YRI', 'Asian')
cols = c()
K = 6

# Pairwise differences between all pairs of populations
#af.diff = NULL
#af.neg.diff = NULL
for(i in 1:2){
  for(j in (i + 1):3){
    cols = append(cols, paste(pop[i], pop[j], sep = '-'))
    #af.diff = cbind(af.diff, bg.af[, i] - bg.af[, j])
#    #af.neg.diff = cbind(af.neg.diff, af.neg[, i] - af.neg[, j])
  }
}
#af.diff.abs = rowMaxs(abs(af.diff), na.rm = T)
#af.neg.diff.abs = rowMaxs(abs(af.neg.diff), na.rm = T)

# pvals = c()
# for(i in 1:2){
#   for(j in (i + 1):3){
#     sel = !is.na(bg.af[, i]) & !is.na(bg.af[, j]) & bg.af[, i] > 0 & bg.af[, j] > 0
#     pvals = append(pvals, median(abs(bg.af[sel, i] - bg.af[sel, j]))) #append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T)$p.value)
#   }
# }

# For each cluster find the SNPs overlapping the regions of the cluster, and compute the allele frequencies between 
# all pairs of individuals.
pop.diff.p = NULL
pop.af = NULL
san.freq = array(0, dim = c(K, 2))
for(k in 1:K){
  cregions = read.bed(file.path(indir, paste('SNYDER_HG19_all_reg_isva_', mark, '_k', K, '_clust', k, '.txt', sep = '')))
  print(dim(cregions))
  ov = findOverlaps(snps.to.ranges(af.pos), regions.to.ranges(cregions), select = 'first', ignore.strand = T)
  #ov = ov[!is.na(ov)]
  pvals = c()
  for(i in 1:2){
    for(j in (i + 1):3){
      #sel = !is.na(af[ov, i]) & !is.na(af[ov, j]) & af[ov, i] > 0 & af[ov, j] > 0
      sel = !is.na(ov) & !is.na(af[, i]) & !is.na(af[, j]) & (af[, i] > 0 & af[, j] > 0)
      bg.sel = bg.af[, i] > 0 & bg.af[, j] > 0
      pvals = append(pvals, wilcox.test(af[sel, i] - af[sel, j], 
                                        bg.af[bg.sel, i] - bg.af[bg.sel, j], paired = F, alternative = 'less')$p.value) #append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T)$p.value)
      #pvals = append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T, alternative = 'two')$p.value) #append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T)$p.value)
    }
  }
  pop.diff.p = rbind(pop.diff.p, pvals)
  
  med.af = c()
  for(i in 1:3){
    sel = !is.na(ov) & af[, i] > 0
    bg.sel = bg.af[, i] > 0
    med.af = append(med.af, wilcox.test(af[sel, i], bg.af[bg.sel, i], alternative = 'less')$p.value)
  }
  pop.af = rbind(pop.af, med.af)
}

# pop.diff.dat = data.frame(-log10(pop.diff.p), row.names = 1:K)
# colnames(pop.diff.dat) = cols
# plot.heatmap(pop.diff.dat, filt.thresh = NA, row.cluster = F, show.dendro = "none", cellnote = matrix(as.character(as.integer(-log10(pop.diff.p))), ncol = 3), 
#              col.cluster = F, row.title= '', col.title = '', cex.row = 4, cex.col = 4, margins = c(18, 18), keysize = 1,
#              dist.metric='spearman', clust.method = "average", break.type='quantile', palette = brewer.pal(9, 'Reds')[2:7],
#              to.file = file.path(indir, paste('SNYDER_HG19_all_reg_isva_', mark, '_k6_af.pdf', sep = '')))