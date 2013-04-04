rm(list=ls())
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(matrixStats)
library(reshape)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('utils/binom.val.r')

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
rdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/'
plotdir = indir
overwrite = F
mark = 'H3K27AC'
pref = 'SNYDER_HG19_all_reg_'
comp = 3
qval = 0.01
K = 4
nrand = 100
pref = paste(pref, mark, '_comp', comp, '_q', qval, sep = '')
clust.pref = paste(pref, '_K', K, sep = '')

bg.af.file = '../../rawdata/variants/af/1000gp_af_1in100.txt'
bg.af = as.matrix(read.table(bg.af.file, sep = '\t', header = F)[, 4:6])
#sel = bg.af[, 1] > 0 & bg.af[, 2] > 0 & bg.af[, 3] > 0 & bg.af[, 1] < 1 & bg.af[, 2] < 1 & bg.af[, 3] < 1
#bg.af = bg.af[sel, ]
bg.af.norm = t(scale(t(bg.af), scale = F))

af.file = '../../rawdata/variants/af/1000gp_af.txt'
regions.file = file.path(indir, paste(pref, '_sign.txt', sep = ''))
tmp.file = file.path(indir, paste(pref, '_sign_af.txt', sep = ''))
#if(!file.exists(tmp.file) || overwrite) system(paste('intersectBed -a', af.file, '-b', regions.file, ' -wa >', tmp.file))
if(!file.exists(tmp.file) || overwrite) system(paste("awk 'BEGIN{OFS=\"\t\"}{print $1,$2-1,$2,\".\",\"0\",\"+\",$3,$4,$5,$6}'", af.file, '| intersectBed -a stdin -b', regions.file, ' -wa >', tmp.file))
af = read.table(tmp.file)
af.pos = af[, c(1,3)]
colnames(af.pos) = c('chr', 'pos')
af = as.matrix(af[, 8:10])
#sel = af[, 1] > 0 & af[, 2] > 0 & af[, 3] > 0 & af[, 1] < 1 & af[, 2] < 1 & af[, 3] < 1
#af = af[sel, ]
#af.pos = af.pos[sel, ]
colnames(af) = c('YRI', 'CEU', 'Asian')
af.norm = t(scale(t(af), scale = F))
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

load(file.path(rdir, paste(pref, '_qn_isvaNull.RData', sep = '')))
regions = regions[isva.fit$deg, ]
load(file.path(rdir, paste(clust.pref, '_clust.RData', sep = '')))
pops = unique(get.pop(colnames(counts)[plot.cols]))
npop = length(pops)

# Pairwise differences between all pairs of populations
#af.diff = NULL
#af.neg.diff = NULL
for(i in 1:(npop-1)){
  for(j in (i + 1):npop){
    cols = append(cols, paste(pops[i], pops[j], sep = '-'))
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
pop.diff.cor = NULL
pop.af = NULL
for(k in 1:K){
  sel = kclusters$cluster == k
  cregions = regions[sel, ]
  print(dim(cregions))
  ov = findOverlaps(snps.to.ranges(af.pos), regions.to.ranges(cregions), select = 'first', ignore.strand = T)
  
  sel = !is.na(ov) & af[, 1] > 0 & af[, 2] > 0 & af[, 3] > 0 & af[, 1] < 1 & af[, 2] < 1 & af[, 3] < 1
  nj.tree = nj(dist(t(af.norm[sel, ]), method = 'manhattan'))
  print(dist(t(af[sel, ]), method = 'manhattan') / sum(sel))
  #pdf(file.path(plotdir, paste(clust.pref, '_clust', k, '_genot_dendro.pdf', sep = '')))
  plot(nj.tree, 'u', cex = 1, edge.width = 0.5, no.margin = T, lab4ut='axial', label.offset = 0.5)
  #dev.off()
  
  #pca.fit = prcomp(t(af.norm[sel, ]), center = F, scale = F)
  #p = plot.pcs(t(af.norm) %*% pca.fit$rotation, pca.fit$rotation, pca.fit$sdev, labels = pops, all = T, ndim = 4)
  #bg.sel = sample(1:nrow(bg.af), sum(sel))
  #pca.fit.bg = prcomp(t(bg.af.norm[bg.sel, ]), center = F, scale = F)
  #print(as.matrix(dist(t(af.norm[sel, ]) %*% pca.fit$rotation)) %/% as.matrix(dist(t(bg.af.norm[bg.sel, ]) %*% pca.fit.bg$rotation)))
  
  #ov = ov[!is.na(ov)]
  pvals = array(NA, dim = c(npop, npop)) 
  rs = array(NA, dim = c(npop, npop)) 
  for(i in 1:2){
    for(j in (i + 1):3){
      #sel = !is.na(af[ov, i]) & !is.na(af[ov, j]) & af[ov, i] > 0 & af[ov, j] > 0
      sel = !is.na(ov) & !is.na(af[, i]) & !is.na(af[, j]) & af[, i] > 0 & af[, j] > 0 & af[, i] < 1 & af[, j] < 1
      bg.sel = bg.af[, i] > 0 & bg.af[, j] > 0 & bg.af[, i] < 1 & bg.af[, j] < 1
      print(sum(bg.sel))
      row.idx = match(colnames(af)[i], pops)
      col.idx = match(colnames(af)[j], pops)
      
      rs[row.idx, col.idx] = mean(abs(af[sel, i] - af[sel, j]) / rowMins(af[sel, c(i,j)])) / 
        mean(abs(bg.af[bg.sel, i] - bg.af[bg.sel, j]) / rowMins(bg.af[bg.sel, c(i,j)]))
      rs[col.idx, row.idx] = rs[row.idx, col.idx]
      
      pvals[row.idx, col.idx] = wilcox.test(af[sel, i], af[sel, j], paired = T, alternative = 'two')$p.value#wilcox.test(abs(af[sel, i] - af[sel, j]) / rowMins(af[sel, c(i,j)]), 
                                #abs(bg.af[bg.sel, i] - bg.af[bg.sel, j]) / rowMins(bg.af[bg.sel, c(i,j)]), 
                                #        paired = F, alternative = 'greater')$p.value #append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T)$p.value)
      pvals[col.idx, row.idx] = pvals[row.idx, col.idx]
      #pvals = append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T, alternative = 'two')$p.value) #append(pvals, wilcox.test(af[sel, i], af[sel, j], paired = T)$p.value)
    }
  }
  rownames(rs) = pops
  colnames(rs) = pops
  rownames(pvals) = pops
  colnames(pvals) = pops
  print(rs)
  print(pvals)
  log.pvals = -log10(pvals)
  log.pvals[log.pvals > 16] = 16
  p = plot.tile(log.pvals, midpoint = 1, y.ord.samples = pops[seq(length(pops), 1, -1)], 
                ycex = 25, xcex = 25, ycolor = get.pop.col(pops[seq(length(pops), 1, -1)]), xcolor = get.pop.col(pops)) + 
                  scale_fill_gradient2(na.value = 'grey60', midpoint = -log10(0.05), limits = c(1,6))
  ggsave(file.path(plotdir, paste(clust.pref, '_clust', k, '_wilc_af.pdf', sep = '')), p, width = 7, height = 6)
  pop.diff.p = rbind(pop.diff.p, as.vector(tril(pvals)[tril(pvals, -1) > 0]))
  pop.diff.cor = rbind(pop.diff.cor, as.vector(tril(rs)[tril(rs, -1) > 0]))
  
  #med.af = c()
  #for(i in 1:3){
  #  sel = !is.na(ov) & af[, i] > 0
  #  bg.sel = bg.af[, i] > 0
  #  med.af = append(med.af, wilcox.test(af[sel, i], bg.af[bg.sel, i], alternative = 'less')$p.value)
  #}
  #pop.af = rbind(pop.af, med.af)
}

# pop.diff.dat = data.frame(-log10(pop.diff.p), row.names = paste('cluster', 1:K, sep = '_'))
# colnames(pop.diff.dat) = cols
# pop.diff.cor = data.frame(pop.diff.cor, row.names = paste('cluster', 1:K, sep = '_'))
# colnames(pop.diff.cor) = cols
# plot.heatmap(pop.diff.dat, filt.thresh = NA, row.cluster = F, show.dendro = "none", cellnote = matrix(as.character(as.integer(-log10(pop.diff.p))), ncol = 3), 
#              col.cluster = F, row.title= '', col.title = '', cex.row = 4, cex.col = 4, margins = c(18, 18), keysize = 1,
#              palette = brewer.pal(9, 'Reds')[2:7])#,
#              #to.file = file.path(indir, paste(clust.pref, '_af.pdf', sep = '')))
# plot.heatmap(pop.diff.cor, filt.thresh = NA, row.cluster = F, show.dendro = "none", 
#              col.cluster = F, row.title= '', col.title = '', cex.row = 4, cex.col = 4, margins = c(18, 18), keysize = 1,
#              dist.metric='spearman', clust.method = "average", break.type='quantile', palette = brewer.pal(9, 'Reds')[2:7])#,
#              #to.file = file.path(indir, paste(clust.pref, '_af_ratio2Bg.pdf', sep = '')))