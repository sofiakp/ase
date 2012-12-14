rm(list=ls())
library('reshape')
require('ggplot2')
#library('fBasics')
#library('VennDiagram')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps')
countdir = file.path(indir, 'qvals')
plotdir = file.path(indir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = ''

filenames <- list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.genecounts\\.RData', full.name = T)
samples = sample.info(filenames, '\\.genecounts\\.RData$')    
epsilon = 1
nsamples = length(filenames)
q.cut = 0.01

indivs = unique(as.character(samples$indiv))
all.hits = list()
ov.mat = array(0, dim = c(nsamples, nsamples))
corr.mat = array(NaN, dim = c(nsamples, nsamples))
for(i in 1:nsamples){
  load(filenames[i])
  pass1 = !is.nan(gene.info$qval) & gene.info$qval < q.cut
  ov.mat[i, i] = sum(pass1)
  genes = rownames(gene.info)
  #hits1 = sort(paste(gene.info$chr[pass], gene.info$pos[pass], sep = ':'), index.return = T)
  indiv = as.character(samples$indiv[i])
  if(is.null(all.hits[[indiv]])){
    all.hits[[indiv]] = genes[pass1]
  }else{
    all.hits[[indiv]] = union(all.hits[[indiv]], genes[pass1])
  }
  if(i < nsamples){
    for(j in (i + 1):nsamples){
      load(filenames[j])
      stopifnot(all(genes == rownames(gene.info)))
      pass2 = !is.nan(gene.info$qval) & gene.info$qval < q.cut
      #hits2 = sort(paste(hits$chr, hits$pos, sep = ':'), index.return = T)
      #ratios2 = log((hits$ref.counts + 1) / (hits$alt.counts + 1), base = 2)[hits2$ix]
      #ov1 = hits1$x %in% hits2$x
      ov.mat[i, j] = sum(pass1 & pass2)
      ov.mat[j, i] = ov.mat[i, j]
    } 
  }
}
out.size = 11
if(nsamples > 50){
  out.size = 5
}

len.dat = data.frame(len = diag(ov.mat), indiv = samples$indiv, mark = samples$mark)
for(d in indivs){
  len.dat = rbind(len.dat, data.frame(len = length(all.hits[[d]]), indiv = factor(d), mark = 'union'))
}
p <- ggplot(len.dat) + geom_bar(aes(x = indiv, y = len), stat = "identity") + 
  facet_wrap(~mark, scales = "free_y") + scale_y_continuous('# AS SNPs') + scale_x_discrete('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 12), axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 17), strip.text.x = element_text(size = 17))
ggsave(file.path(plotdir, 'gene_nhits.png'))

d = t(diag(ov.mat))
bad = d < 10
ov.mat = ov.mat[!bad, !bad]
d = d[!bad]
samples = samples[!bad, ]
names = paste(samples[, 1], samples[, 2], sep = '_')
ov.dat = data.frame(apply(ov.mat, 2, function(x) x / d))
colnames(ov.dat) = names
ov.dat$sample = names
ov.dat = melt(ov.dat)
# Plot[i,j]: what fraction of AS SNPs of dataset i overlap with dataset j
p = ggplot(ov.dat) + geom_tile(aes(x = variable, y = sample, fill = value)) +
  scale_x_discrete('') + scale_y_discrete('') + scale_fill_gradient('Fraction', limits = c(0, 1), low = "beige", high = "tomato") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey20', size = out.size),
        axis.text.y = element_text(hjust = 1, colour = 'grey20', size = out.size), 
        legend.title = element_blank(), legend.text = element_text(size = 13)) + 
          ggtitle('Fraction of AS SNPs in each dataset overlapping every other dataset')
ggsave(file.path(plotdir, paste('mergedRep_', outpref, 'gene_as_overlaps.png', sep = '')))
