rm(list=ls())
library(GenomicRanges)
library(preprocessCore)
library(matrixStats)
library(reshape)
library(ggplot2)
library(glmnet)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

indir = '../../rawdata/enhancers/rdata/' 
plotdir = '../../rawdata/enhancers/plots/'
pref = 'enhancer_coef_elastic0.5_100kb_asinh0.2_cv0.2'
mark = 'H3K27AC'
nperm = 100
dev.dat = NULL
for(n in 1:nperm){
  load(file.path(indir, 'perm_gene', paste(pref, '_perm_gene', n, '_', mark, '.RData', sep = '')))
  dev.dat.tmp = data.frame(fit.dat)
  dev.dat.tmp$type = rep(paste('perm.gene', n, sep = ''), nrow(fit.dat))
  dev.dat.tmp$gene.idx = 1:nrow(fit.dat)
  dev.dat = rbind(dev.dat, dev.dat.tmp)
}
load(file.path(indir, paste(pref, '_', mark, '.RData', sep = '')))
dev.dat.tmp = data.frame(fit.dat)
dev.dat.tmp$type = rep(paste('true', sep = ''), nrow(fit.dat))
dev.dat.tmp$gene.idx = 1:nrow(fit.dat)
dev.dat = rbind(dev.dat, dev.dat.tmp)
dev.dat2 = dev.dat[, c(2,5,6)]
dev.dat2 = data.frame(cast(dev.dat2, gene.idx~type, value = 'dev.ratio'))

dev.dat$type = factor(gsub('perm.gene[0-9]+', 'perm.gene', as.character(dev.dat$type)))
for(i in 1:4){
  dev.dat.tmp = dev.dat[, c(i, 5, 6)]
  colnames(dev.dat.tmp) = c('x', 'type', 'gene.idx')
  p = ggplot(dev.dat.tmp) + geom_density(aes(x = x, color = type)) + 
    xlab(colnames(dev.dat)[i]) + ylab('density') + scale_color_discrete('') +
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 15), 
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 14))
  ggsave(file = file.path(plotdir, paste(pref, '_', mark, '_', colnames(dev.dat)[i], '.png')), p, width = 6.5, height = 5.6)
}

max.dev = apply(as.matrix(dev.dat2[, grep('perm', colnames(dev.dat2))]), 1, function(x) quantile(x, 0.99, na.rm = T))
#rowMaxs(as.matrix(dev.dat2[, grep('perm', colnames(dev.dat2))]))
sel.genes = !is.na(dev.dat2$true) & dev.dat2$true > max.dev & dev.dat2$true > 0.8

orig.coef = coef.dat
coef.dat = coef.dat[coef.dat$gene.idx %in% which(sel.genes) & abs(coef.dat$coef) > 0, ]
links = NULL
for(i in 1:nrow(coef.dat)){
  model.dat = data.frame(y = rna.counts[coef.dat$gene.idx[i], ], x = ac.counts[coef.dat$region.idx[i], ])
  simple.model = glm(y~x, family = 'gaussian', model.dat)
  model.sum = summary(simple.model)
  pval = model.sum$coefficients[2, 4]
  if(pval < 0.01){
    if(is.null(links)){links = cbind(coef.dat$gene.idx[i], coef.dat$region.idx[i], pval)
    }else{links = rbind(links, cbind(coef.dat$gene.idx[i], coef.dat$region.idx[i], pval))} 
  }
}
links = data.frame(links)
colnames(links) = c('gene.idx', 'region.idx', 'pval')
links = links[order(links$pval), ]
coef.dat = orig.coef
save(links, ac.counts, rna.counts, gene.meta, ac.regions, coef.dat, tss,
     file = file.path(indir, paste(pref, '_', mark, '_links.RData', sep = '')))
# load(paste(pref, '_perm_indiv_', mark, '.RData', sep = ''))
# coef.dat$type = 'perm_indiv'
# coef.all = coef.dat
# dev.dat = data.frame(gene.idx = 1:length(dev), dev = dev, type = rep('perm_indiv', length(dev)))
# load(paste(pref, '_perm_gene_', mark, '.RData', sep = ''))
# coef.dat$type = 'perm_gene'
# coef.all = rbind(coef.all, coef.dat)
# dev.dat = rbind(dev.dat, data.frame(gene.idx = 1:length(dev), dev = dev, type = rep('perm_gene', length(dev))))
# load(paste(pref, '_', mark, '.RData', sep = ''))
# coef.dat$type = 'true'
# coef.all = rbind(coef.all, coef.dat)
# coef.all$coef = asinh(coef.all$coef)
# dev.dat = rbind(dev.dat, data.frame(gene.idx = 1:length(dev), dev = dev, type = rep('true', length(dev))))
#dev.dat = 
#num.coef = cast(coef.all, gene.idx~type, function(x) sum(abs(x) > 0), value = 'coef')

# num.coef = array(NaN, dim = c(nrow(rna.counts), 3))
# num.coef[num.coef.tmp$gene.idx, ] = as.matrix(num.coef.tmp)
# num.coef = data.frame(num.coef)
# colnames(num.coef) = colnames(num.coef.tmp)[2:4]
# num.coef$gene.idx = 1:nrow(rna.counts)
# num.coef = num.coef[, c(4,1,2,3)]
