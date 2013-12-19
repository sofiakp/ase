rm(list=ls())
library(GenomicRanges)
library(preprocessCore)
library(matrixStats)
library(reshape)
library(ggplot2)
library(glmnet)
library(foreach)
library(doMC)
source('utils/sample.info.r')
source('utils/deseq.utils.r')

# Uses the ARS statistic to link enhancers to genes. You must run enhancer.regress.r first to get the mappings
# between genes and promoter regions.

# Computes the ARS statistic. x and y are matrices NxM where N is the number of associations (pairs) 
# for which the statistic will be computed, and M is the number of cell lines.
# You should divide each row of x and y by the row max and median center before calling this.
get.ars = function(x, y){
  d = sqrt(x^2 + y^2)
  r = d / rep(rowMedians(d), ncol(d))
  theta = 180 * acos(x / d) / pi # The angle between (x,y) and the x-axis is acos(x/d). acos returns radians
  theta = pmin(abs(45 - theta), abs(135 - theta))
  theta[is.na(theta)] = 0
  return(list(r = r, theta = theta))
}

# Checks whether all the cell lines that have ARS > 0.5 for a gene agree in the direction of influence. 
get.ars.direction = function(ac.counts.norm, rna.counts.norm, ars.score.norm){
  # The first product will be > 0 if ac and rna agree in direction and < 0 if the disagree. 
  # Non-significant cell lines are zeroed out. 
  mult = ac.counts.norm * rna.counts.norm * (ars.score.norm > 0.5)
  pos.effect.ind = apply(mult, 1, function(x) all(x >= 0)) # All significant lines suggest a positive interaction between ac and rna.
  neg.effect.ind = apply(mult, 1, function(x) all(x <= 0)) # All significant lines suggest a negative interaction between ac and rna.
  return(list(pos.effect.ind = pos.effect.ind, neg.effect.ind = neg.effect.ind))
}

#registerDoMC(5)
reg = new.env()
infile = '../../rawdata/enhancers/merged_Mar13/rdata/enhancer_coef_elastic0.5_100kb_asinh0.2_cv0.2_H3K27AC.RData' 
load(infile, reg)
outdir = '../../rawdata/enhancers/merged_Mar13/rdata'
plotdir = '../../rawdata/enhancers/merged_Mar13/plots'
outpref = 'enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01'
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)

ac.counts = reg$ac.counts[reg$coef.dat$region.idx, ]
rna.counts = reg$rna.counts[reg$coef.dat$gene.idx, ]
good.rows = !grepl('chr[XY]', reg$ac.regions$chr[reg$coef.dat$region.idx]) & rowMaxs(rna.counts) > asinh(50) # Remove sex chromosomes
ac.counts = ac.counts[good.rows, ]
rna.counts = rna.counts[good.rows, ]

row.max = rowMaxs(ac.counts) # There shouldn't be any NaNs
ac.counts = apply(ac.counts, 2, function(x) x / row.max) # Divide each row by the maximum of the row
row.median = rowMedians(ac.counts)
ac.counts.norm = apply(ac.counts, 2, function(x) x - row.median) # Subtract the row median
row.max = rowMaxs(rna.counts) 
rna.counts = apply(rna.counts, 2, function(x) x / row.max)
row.median = rowMedians(rna.counts)
rna.counts.norm = apply(rna.counts, 2, function(x) x - row.median)
nindiv = ncol(rna.counts)
ncand = nrow(rna.counts)
indivs = colnames(rna.counts)
fdr = 0.01 # FDR cutoff for choosing significant associations
nperm = 100
perm = 'gene_pairs' # Type of permutations
outpref = paste(outpref, '_perm_', perm, sep = '')
      
ars = get.ars(ac.counts.norm, rna.counts.norm)
ars.score = ars$r * exp(-0.01 * ars$theta)
ars.max = rowMaxs(ars.score, na.rm = T) # Max ARS across cell lines
ars.score.norm = apply(ars.score, 2, function(x) x / ars.max) # Normalized ARS scores
max.line = apply(ars.score, 1, function(x) which.max(x)) # Get the cell line with the max ARS for each link
ars.dir = get.ars.direction(ac.counts.norm, rna.counts.norm, ars.score.norm)

# Do permutations to compute null distribution of max ARS
#ars.score.perm = array(NaN, dim = c(nperm * ncand, 1))
ars.score.perm = foreach(n = 1:nperm, .combine = 'append') %do% {
  # if(n %% 10 == 0) cat('Permutation', n, '\n')
  if(perm == 'gene_pairs'){
    ac.counts.perm = ac.counts # Start with the max normalized (but not median centered) data.
    rna.counts.perm = rna.counts
    for(i in 1:nindiv){
      # Shuffle the rows within each column. Shuffle genes and regions in the same way.
      shuf = sample(1:ncand, ncand)
      ac.counts.perm[, i] = ac.counts[shuf, i]
      rna.counts.perm[, i] = rna.counts[shuf, i]
    }  
    row.median = rowMedians(ac.counts.perm)
    ac.counts.perm = apply(ac.counts.perm, 2, function(x) x - row.median)
  }else if(perm == 'gene'){
    shuf = sample(1:ncand, ncand)
    rna.counts.perm = rna.counts[shuf, ]
  }
  row.median = rowMedians(rna.counts.perm)
  rna.counts.perm = apply(rna.counts.perm, 2, function(x) x - row.median)
  ars.perm = get.ars(ac.counts.perm, rna.counts.perm)
  # ars.score.perm[((n - 1) * ncand + 1):(n * ncand)] = rowMaxs(ars.perm$r * exp(-0.01 * ars.perm$theta), na.rm = T)
  rowMaxs(ars.perm$r * exp(-0.01 * ars.perm$theta), na.rm = T)
}

test.values = median(ars.max):max(ars.max) # Possible cutoffs for significance
fdrs = array(NaN, dim = c(length(test.values), 1))
for(i in 1:length(test.values)){
  fdrs[i] = sum(ars.score.perm > test.values[i]) / (nperm * sum(ars.max > test.values[i]))
}
sel.cut = which.min(abs(fdrs - fdr)) # Find the cutoff that yields an FDR as close as possible to the desired FDR
effect.fdr = fdrs[sel.cut]

# Order the links by decreasing ARS (BE CAREFUL WITH THE REORDERING)
ord = order(ars.max, decreasing = T)
ars.score = ars.score[ord, ]
colnames(ars.score) = indivs
ars.score.norm = ars.score.norm[ord, ]
colnames(ars.score.norm) = indivs
coef.dat = reg$coef.dat[good.rows, ]
coef.dat$max.line = max.line
coef.dat$ars.max = ars.max
coef.dat = coef.dat[ord, ]
links = data.frame(gene.idx = coef.dat$gene.idx, region.idx = coef.dat$region.idx, 
                   coef = coef.dat$ars.max, max.line = coef.dat$max.line, pos.effect = ars.dir$pos.effect.ind[ord])

# Select the associations with large ARS that are also consistent across cell lines.
sel.assoc = coef.dat$ars.max > test.values[sel.cut] & (ars.dir$pos.effect.ind | ars.dir$neg.effect.ind)[ord]
links = links[sel.assoc, ]

# Make some plots
max.counts = data.frame(as.matrix(table(links$max.line))) # Count how many times each individual is the maximum
colnames(max.counts) = c('hits')
max.counts$indiv = indivs[as.numeric(rownames(max.counts))]
p = ggplot(max.counts) + geom_bar(aes(x = indiv, y = hits), stat = 'identity') + xlab('') + ylab('# times individual has max ARS') +
  theme(axis.text.x = element_text(size = 12, angle = -45, vjust = 1, hjust = 0), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 15))
ggsave(file.path(plotdir, paste(outpref, '_hitsPerIndiv.pdf', sep = '')), p, width = 6.5, height = 5.6)

sum.sign = data.frame(n = table(rowSums(ars.score.norm[sel.assoc, ] > 0.5))) # For each significant association how many individuals are "outliers"?
colnames(sum.sign) = c('n', 'hits')
p = ggplot(sum.sign) + geom_bar(aes(x = n, y = hits), stat = 'identity') + 
  xlab('# individuals with ARS > max(ARS)/2') + ylab('# links') + 
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 15))
ggsave(file.path(plotdir, paste(outpref, '_indivPerHit.pdf', sep = '')), p, width = 6.5, height = 5.6)

# Density plots for normalized ARS scores for significant and non-significant associations
score.dat = melt(data.frame(ars.score.norm[sel.assoc, ]))
score.dat$type = rep('significant', nrow(score.dat))
score.dat.tmp = data.frame(ars.score.norm[!sel.assoc, ])
score.dat.tmp = melt(score.dat.tmp[sample(1:nrow(score.dat.tmp), min(nrow(score.dat.tmp), sum(sel.assoc))), ])
score.dat.tmp$type = rep('non-significant', nrow(score.dat.tmp))
score.dat = rbind(score.dat, score.dat.tmp)
score.dat$type = factor(score.dat$type)
p = ggplot(score.dat) + geom_density(aes(x= value, color = variable, linetype = type)) + 
  xlab('ARS score normalized by max in region') + ylab('density') + scale_color_discrete('') + scale_linetype('') +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 15), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12))
ggsave(file.path(plotdir, paste(outpref, '_normARSdens.pdf', sep = '')), p, width = 9.75, height = 8.4)

# Save results
ac.counts = reg$ac.counts
ac.regions = reg$ac.regions
rna.counts = reg$rna.counts
gene.meta = reg$gene.meta
fdr = effect.fdr
tss = reg$tss
save(ac.regions, ac.counts, rna.counts, gene.meta, coef.dat, ars.score, ars.score.norm, fdr, links, sel.assoc, tss, 
     file = file.path(outdir, paste(outpref, '.RData', sep = '')))
