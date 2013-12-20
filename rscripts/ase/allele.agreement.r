rm(list=ls())
library('reshape')
require('ggplot2')
source('utils/sample.info.r')
source('utils/deseq.utils.r')

# Computes correlation of allelic biases and overlaps of AS SNPs between selected pairs of files
# (cell lines or marks).

indir = '../../rawdata/alleleCounts/allNonSan/rdata/reps'
countdir = file.path(indir, 'qvals')
plotdir = file.path(indir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = 'active_' # eg. foo_

# Read overlap matrices from here instead of regenerating them. Set to NULL to generate from scratch.
read.from = NULL #'rawdata/alleleCounts/allNonSan/rdata/reps/plots/mergedRep_overlaps.RData')
filter.regions = read.bed('../../rawdata/segSignal/14indiv/allEnhStates.bed')
load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
good.snps =  !is.na(findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(filter.regions), select = 'first', ignore.strand = T))

# GM(19238|19239|19240|18505|19099|18486) all YRI
# (GM(12878|12891|12892|12890|10847)|SNYDER) all CEU
# activating (H3K4ME1|H3K4ME3|H3K27AC|H3K9AC)
filenames = list.files(countdir, pattern = 'SNYDER_HG19.*(H3K27AC).*rep\\.RData', full.name = T)
samples = sample.info(filenames, '\\.RData$')

epsilon = 1
nsamples = length(filenames)
q.cut = -2 # in log10

if(is.null(read.from)){
  ov.mat = array(0, dim = c(nsamples, nsamples)) # What fraction of the AS SNPs for one dataset overlap the AS SNPs of another?
  corr.mat = array(NaN, dim = c(nsamples, nsamples)) # Correlation of allellic biases from different datasets
  dir.mat = array(NaN, dim = c(nsamples, nsamples)) # Fraction of AS SNPs that are in the same direction between datasets
  
  for(i in 1:nsamples){
    print(basename(filenames[i]))
    load(filenames[i])
    hits1 = snp.info$qval < q.cut # Significant hits
    pass1 = snp.info$pass # heterozygous, unmasked SNPs
    hits1 = hits1 & pass1
    ov.mat[i, i] = sum(hits1)
    dir.mat[i, i] = 1
    corr.mat[i, i] = 1
    ratios1 = log((counts[, 1] + 1) / (counts[, 2] + 1), base = 2)
    #hits1 = hits1 & ratios1 > log2(1.2)
    #if(i < nsamples){
    for(j in 1:nsamples){
      if(i == j) next
      load(filenames[j])
      hits2 = snp.info$qval < q.cut & snp.info$pass
      ratios2 = log((counts[, 1] + 1) / (counts[, 2] + 1), base = 2)
      #hits2 = hits2 & ratios2 > log2(1.2)
      ov.mat[i, j] = sum(hits1 & hits2) # SNPs that were significantly AS in both datasets
      # ov.mat[j, i] = ov.mat[i, j]
      # For the SNPs that are AS in at least one of the two datasets, and heterozygous (and unmasked)
      # in the other dataset, compute the correlation of allelic biases.
      hit.union = as.vector(((hits1 | hits2) & good.snps)) # | (hits2 & pass1))
      if(sum(hit.union) > 10){
        corr.mat[i, j] = cor(ratios1[hit.union], ratios2[hit.union]) 
        dir.mat[i, j] = sum(ratios1[hit.union] * ratios2[hit.union] > 0) / sum(hit.union)
        #corr.mat[j, i] = corr.mat[i, j]
        #dir.mat[j, i] = dir.mat[i, j]
        if(i < j){
          hit.dat = data.frame(cbind(x = ratios1[hit.union], y = ratios2[hit.union]))
          q0 = ggplot(hit.dat) + geom_point(aes(x = x, y = y)) + xlab(samples$mark[i]) + ylab(samples$mark[j]) + 
            theme_bw() + ggtitle(paste('Correlation of allelic biases', sprintf('%.4f', corr.mat[i,j]))) +
            theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14))
          ggsave(file.path(plotdir, paste('mergedRep_', outpref, samples$indiv[i], '_', samples$mark[i], '_', samples$indiv[i], '_', samples$mark[j], '_corr.pdf', sep = '')), q0, width = 6.5, height = 5.6)                    
        }
      }
    } 
  }
  d = t(diag(ov.mat))
  bad = d < 10
  ov.mat = ov.mat[!bad, !bad]
  corr.mat = corr.mat[!bad, !bad]
  dir.mat = dir.mat[!bad, !bad]
  d = d[!bad]
  samples = samples[!bad, ]
  names = paste(samples[, 1], samples[, 2], sep = '_')
  ov.dat = data.frame(apply(ov.mat, 2, function(x) x / d))
  colnames(ov.dat) = names
  ov.dat$sample = names
  ov.dat = melt(ov.dat)
  
  corr.dat = data.frame(corr.mat)
  colnames(corr.dat) = names
  corr.dat$sample = names
  corr.dat = melt(corr.dat)
  
  dir.dat = data.frame(dir.mat)
  colnames(dir.dat) = names
  dir.dat$sample = names
  dir.dat = melt(dir.dat)
  save(dir.dat, corr.dat, ov.dat, file = file.path(plotdir, paste('mergedRep_', outpref, 'overlaps.RData', sep = '')))
}else{
  load(read.from)
  sel = dir.dat$sample %in% paste(samples$indiv, samples$mark, sep = '_') & dir.dat$variable %in% paste(samples$indiv, samples$mark, sep = '_')
  dir.dat = dir.dat[sel, ]
  corr.dat = corr.dat[sel, ]
  ov.dat = ov.dat[sel, ]
}

out.size = 14
if(nsamples > 50){out.size = 5}
if(nsamples > 100){out.size = 4}

#pop = sort(get.pop(as.character(samples$indiv)), index.return = T)

# Plot[i,j]: what fraction of AS SNPs of dataset i overlap with dataset j
q1 = ggplot(ov.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
  scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient('Fraction', limits = c(0, 1), low = "beige", high = "tomato") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30', size = out.size),
        axis.text.y = element_text(hjust = 1, colour = 'grey30', size = out.size), plot.title = element_text(size = 14)) + 
          ggtitle('Fraction of overlapping AS SNPs')
ggsave(file.path(plotdir, paste('mergedRep_', outpref, 'overlaps.pdf', sep = '')), q1, width = 6.5, height = 5.6)
q2 = ggplot(corr.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
  scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient2('Fraction', limits = c(-0.5, 0.5), low = "blue", mid = "white", high = "tomato") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30', size = out.size),
        axis.text.y = element_text(hjust = 1, colour = 'grey30', size = out.size), plot.title = element_text(size = 14)) +
          ggtitle('Correlation of biases in heterozygous SNPs')
ggsave(file.path(plotdir, paste('mergedRep_', outpref, 'overlap_corr.pdf', sep = '')), q2, width = 6.5, height = 5.6)
q3 = ggplot(dir.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) + ggtitle('Fraction of overlapping AS SNPs with the same direction') +
  scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient2('Fraction', limits = c(0, 1), midpoint=0.5,low = "blue", mid = "white", high = "tomato") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30', size = out.size),
       axis.text.y = element_text(hjust = 1, colour = 'grey30', size = out.size), plot.title = element_text(size = 14))
ggsave(file.path(plotdir, paste('mergedRep_', outpref, 'overlap_agreement.pdf', sep = '')), q3, width = 6.5, height = 5.6)


# filenames <- list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.hits\\.RData', full.name = T)
# samples = sample.info(filenames, '\\.hits\\.RData$')    
# epsilon = 1
# nsamples = length(filenames)
# q.cut = 0.01
# 
# indivs = unique(as.character(samples$indiv))
# all.hits = list()
# ov.mat = array(0, dim = c(nsamples, nsamples))
# corr.mat = array(NaN, dim = c(nsamples, nsamples))
# for(i in 1:nsamples){
#   load(filenames[i])
#   ov.mat[i, i] = dim(hits)[1]
#   hits1 = sort(paste(hits$chr, hits$pos, sep = ':'), index.return = T)
#   indiv = as.character(samples$indiv[i])
#   if(is.null(all.hits[[indiv]])){
#     all.hits[[indiv]] = hits1$x
#   }else{
#     all.hits[[indiv]] = union(all.hits[[indiv]], hits1$x)
#   }
#   ratios1 = log((hits$ref.counts + 1) / (hits$alt.counts + 1), base = 2)[hits1$ix]
#   if(i < nsamples){
#     for(j in (i + 1):nsamples){
#       load(filenames[j])
#       hits2 = sort(paste(hits$chr, hits$pos, sep = ':'), index.return = T)
#       ratios2 = log((hits$ref.counts + 1) / (hits$alt.counts + 1), base = 2)[hits2$ix]
#       ov1 = hits1$x %in% hits2$x
#       ov.mat[i, j] = sum(ov1)
#       ov.mat[j, i] = ov.mat[i, j]
#       
#       ov2 = hits2$x %in% hits1$x
#       if(sum(ov1) > 10){
#         # Since hits are sorted, overlaps will be sorted, so we know they are in the same order in the two samples
#         corr.mat[i, j] = sum(ratios1[ov1] * ratios2[ov2] > 0) / sum(ov1) #cor(ratios[[i]][ov1], ratios[[j]][ov2])
#         corr.mat[j, i] = corr.mat[i, j]
#       }
#     } 
#   }
# }
# out.size = 10
# if(nsamples > 50){
#   out.size = 5
# }
# 
# len.dat = data.frame(len = diag(ov.mat), indiv = samples$indiv, mark = samples$mark)
# for(d in indivs){
#   len.dat = rbind(len.dat, data.frame(len = length(all.hits[[d]]), indiv = factor(d), mark = 'union'))
# }
# p <- ggplot(len.dat) + geom_bar(aes(x = indiv, y = len), stat = "identity") + 
#   facet_wrap(~mark, scales = "free_y") + scale_y_continuous('# AS SNPs') + scale_x_discrete('') +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
# ggsave(file.path(plotdir, 'nhits.png'))
# 
# d = t(diag(ov.mat))
# bad = d < 10
# ov.mat = ov.mat[!bad, !bad]
# corr.mat = corr.mat[!bad, !bad]
# d = d[!bad]
# samples = samples[!bad, ]
# names = paste(samples[, 1], samples[, 2], sep = '_')
# ov.dat = data.frame(apply(ov.mat, 2, function(x) x / d))
# colnames(ov.dat) = names
# ov.dat$sample = names
# ov.dat = melt(ov.dat)
# # Plot[i,j]: what fraction of AS SNPs of dataset i overlap with dataset j
# p = ggplot(ov.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
#   scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient('Fraction', limits = c(0, 1), low = "beige", high = "tomato") +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30', size = out.size),
#         axis.text.y = element_text(hjust = 1, colour = 'grey30', size = out.size)) + 
#           ggtitle('Fraction of AS SNPs in each dataset overlapping every other dataset')
# ggsave(file.path(plotdir, paste('mergedRep', outpref, 'overlaps.png', sep = '_')))
# 
# corr.dat = data.frame(corr.mat)
# colnames(corr.dat) = names
# corr.dat$sample = names
# corr.dat = melt(corr.dat)
# p = ggplot(corr.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
#   scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient('Fraction', limits = c(0, 1), low = "beige", high = "tomato") +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30', size = out.size),
#         axis.text.y = element_text(hjust = 1, colour = 'grey30', size = out.size)) +
#           ggtitle('Fraction of overlapping AS SNPs with the same direction')
# ggsave(file.path(plotdir, paste('mergedRep', outpref, 'overlap_agreement.png', sep = '_')))
