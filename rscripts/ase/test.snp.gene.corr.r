rm(list=ls())
library('GenomicRanges')
library('reshape')
require('ggplot2')

mark = 'H3K27AC'
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

merged = new.env()
load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/merged/',
               paste('SNYDER_HG19', mark, 'counts.RData', sep = '_')), merged)
samples = merged$samples
good.rows = rowSums(merged$counts) > 4 * dim(samples)[1]
snp.counts = merged$counts[good.rows, ]
snp.info = merged$snp.info[good.rows, ]
sel.samples = array(F, dim = c(dim(snp.counts)[2], 1))

rna.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'reps/qvals/')
gene.counts = NULL
for(i in 1:dim(snp.counts)[2]){
  rna.file = file.path(rna.dir, paste('SNYDER_HG19_', samples$indiv[i], '_RNA_rep.genecounts.RData', sep  = ''))
  if(file.exists(rna.file)){
    sel.samples[i] = T
    load(rna.file)
    widths = gene.info$end - gene.info$start
    if(is.null(gene.counts)){
      gene.counts = rowSums(counts) / widths
    }else{
      gene.counts = cbind(gene.counts, rowSums(counts) / widths)
    }
  }
}
snp.counts = snp.counts[, sel.samples]
samples = samples[sel.samples, ]
good.rows = rowMeans(gene.counts) > quantile(rowMeans(gene.counts), 0.4)
gene.counts = gene.counts[good.rows, ]
gene.info = gene.info[good.rows, ]

rand.gene.counts = gene.counts
for(i in 1:dim(rand.gene.counts)[2]){
  rand.gene.counts[, i] = sample(gene.counts[, i])  
}

snp.cor = NULL
for(i in 1:dim(snp.counts)[1]){
  snp.chr = as.character(snp.info$chr[i])
  dists = abs(snp.info$pos[i] - gene.info$start)
  sel.genes = as.character(gene.info$chr) == snp.chr & dists > 5000 & dists < 25000
  if(any(sel.genes)){
    if(sum(sel.genes) == 1){
      co = t(cor(snp.counts[i, ], gene.counts[sel.genes, ]))
      randc = t(cor(snp.counts[i, ], rand.gene.counts[sel.genes, ]))
    }else{
      co = t(cor(snp.counts[i, ], t(gene.counts[sel.genes, ])))
      randc = t(cor(snp.counts[i, ], t(rand.gene.counts[sel.genes, ])))
    }
    if(is.null(snp.cor)){
      snp.cor = cbind(co, randc)
    }else{snp.cor = rbind(snp.cor, cbind(co, randc))}
  }
}
snp.cor = data.frame(neighbor = snp.cor[, 1], random = snp.cor[, 2])
melt.cor = melt(snp.cor)
p = ggplot(melt.cor) + geom_density(aes(x = value, y = ..density.., color = variable), adjust = 1) +
  scale_x_continuous('Correlation') + scale_y_continuous('Density') + ggtitle(paste('Correlation between', mark, 'at SNPs and gene expression')) +
  scale_color_discrete('') + theme_bw()
ggsave(file.path(plotdir, paste(mark, '_gene_cor_at_snps.png', sep = '')), width = 13.6, height = 11.8)