rm(list=ls())
library('GenomicRanges')
library('reshape')
require('ggplot2')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

marks = c('H3K27AC', 'H3K4ME1', 'H3K4ME3')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)
mark.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'alleleCounts', 'reps_15Sep12/merged')
rna.files = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'reps/qvals/'),
                       pattern = paste('SNYDER_HG19_.*_H3K36ME3_rep.genecounts.RData', sep  = ''), full.name = T)
indivs = as.character(sample.info(rna.files, '.genecounts.RData')$indiv)
nindivs = length(indivs)

gene.counts = NULL
for(i in 1:nindivs){
  load(rna.files[i])
  widths = gene.info$end - gene.info$start
  if(is.null(gene.counts)){
    gene.counts = rowSums(counts) / widths
  }else{
    gene.counts = cbind(gene.counts, rowSums(counts) / widths)
  }
}
good.rows = rowMeans(gene.counts) > quantile(rowMeans(gene.counts), 0.4)
gene.counts = gene.counts[good.rows, ]
gene.info = gene.info[good.rows, ]
rand.gene.counts = gene.counts
for(i in 1:dim(rand.gene.counts)[2]){
  rand.gene.counts[, i] = sample(gene.counts[, i])  
}

for(j in 1:length(marks)){
  merged = new.env()
  load(file.path(mark.dir, paste('SNYDER_HG19_', marks[j], '_counts.RData', sep  = '')), merged)
  #merged.counts = merged$counts
  #colnames(merged.counts) = merged$samples$indiv
  tmp.signal = merged$counts[, merged$samples$indiv %in% indivs][, order(indivs)]
  if(j == 1){
    mark.signal = array(0, dim = c(dim(merged$counts)[1], nindivs, length(marks)))
    good.rows = rowSums(tmp.signal) > 4 * nindivs
    snp.info = merged$snp.info
  }else{
    good.rows = good.rows & rowSums(tmp.signal) > 4 * nindivs
    stopifnot(snp.info$chr == merged$snp.info$chr, snp.info$pos == merged$snp.info$pos)
  }
  mark.signal[,,j] = tmp.signal
}
snp.info = snp.info[good.rows, ]
mark.signal = mark.signal[good.rows,,]

snp.cor = NULL
rand.snp.cor = NULL
for(i in 1:dim(snp.info)[1]){
  snp.chr = as.character(snp.info$chr[i])
  dists = abs(snp.info$pos[i] - gene.info$start)
  sel.genes = as.character(gene.info$chr) == snp.chr & dists > 5000 & dists < 25000
  if(any(sel.genes)){
    if(sum(sel.genes) == 1){
      co = t(cor(mark.signal[i,,], gene.counts[sel.genes, ]))
      randc = t(cor(mark.signal[i,,], rand.gene.counts[sel.genes, ]))
    }else{
      co = t(cor(mark.signal[i,,], t(gene.counts[sel.genes, ])))
      randc = t(cor(mark.signal[i,,], t(rand.gene.counts[sel.genes, ])))
    }
    if(is.null(snp.cor)){
      snp.cor = as.matrix(co)
      rand.snp.cor = randc
    }else{
      snp.cor = rbind(snp.cor, co)
      rand.snp.cor = rbind(rand.snp.cor, randc)
    }
  }
}
cor.dat = data.frame(rbind(snp.cor, rand.snp.cor))
colnames(cor.dat) = marks
cor.dat$labels = factor(append(rep(1, dim(snp.cor)[1]), rep(0, dim(rand.snp.cor)[1])), levels = c(0, 1), labels = c('random', 'neighbor'))
log.fit = glm(paste('labels ~', paste(marks, collapse = '+')), data = cor.dat, family = "binomial", y = T)

fit.dat = data.frame(label = cor.dat$labels, fit = predict(log.fit))
p = ggplot(fit.dat) + geom_density(aes(x = fit, y = ..density.., color = label), adjust = 1) +
  scale_x_continuous('Fitted log(p/1-p)') + scale_y_continuous('Density') + ggtitle(paste('Logistic regression using ', paste(marks, collapse = ' '))) +
  scale_color_discrete('') + theme_bw()
ggsave(file.path(plotdir, paste('gene_fit_at_snps_', paste(marks, collapse = '_'), '.png', sep = '')), width = 13.6, height = 11.8)