rm(list=ls())
library('GenomicRanges')
library('reshape')
require('ggplot2')

mark = 'H3K27AC'
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'rdata/reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

peak.tab = read.table('rawdata/spp.optimal.wgEncodeBroadHistoneGm12878CtcfStdAlnRep0_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.narrowPeak')
peak.ranges = GRanges(seqnames = Rle(peak.tab[,1]), ranges = IRanges(start = peak.tab[, 2], end = peak.tab[, 3]), strand = Rle(rep('+', dim(peak.tab)[1])))

load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

merged = new.env()
load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/merged/',
               paste('SNYDER_HG19', mark, 'counts.RData', sep = '_')), merged)
samples = merged$samples
good.rows = rowSums(merged$counts) > 4 * dim(samples)[1]
snp.counts = merged$counts[good.rows, ]
snp.info = merged$snp.info[good.rows, ]
sel.samples = array(F, dim = c(dim(snp.counts)[2], 1))

rna.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'rdata/reps/qvals/')
gene.counts = NULL
for(i in 1:dim(snp.counts)[2]){
  rna.file = file.path(rna.dir, paste('SNYDER_HG19_', samples$indiv[i], '_RZ_rep.genecounts.RData', sep  = ''))
  if(file.exists(rna.file)){
    sel.samples[i] = T
    load(rna.file)
    widths = gene.meta$len
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
gene.meta = gene.meta[good.rows, ]
rand.gene.counts = gene.counts
for(i in 1:dim(rand.gene.counts)[2]){
  rand.gene.counts[, i] = sample(gene.counts[, i])  
}

snp.cor = NULL
for(i in 1:dim(snp.counts)[1]){
  snp.chr = as.character(snp.info$chr[i])
  dists = abs(snp.info$pos[i] - gene.meta$start)
  sel.genes = as.character(gene.meta$chr) == snp.chr & dists > 5000 & dists < 25000
  if(any(sel.genes)){
    ov.region = GRanges(seqnames = snp.chr, ranges = IRanges(start = min(snp.info$pos[i], gene.meta$start[sel.genes]), 
                                                             end = max(snp.info$pos[i], gene.meta$start[sel.genes])), 
                        strand = Rle(rep('+', sum(sel.genes)))) 
    if(sum(sel.genes) == 1){
      co = t(cor(snp.counts[i, ], gene.counts[sel.genes, ]))
      randc = t(cor(snp.counts[i, ], rand.gene.counts[sel.genes, ]))
    }else{
      co = t(cor(snp.counts[i, ], t(gene.counts[sel.genes, ])))
      randc = t(cor(snp.counts[i, ], t(rand.gene.counts[sel.genes, ])))
    }
    if(is.null(snp.cor)){
      snp.cor = cbind(co, randc)
      has.ctcf = ov.region %in% peak.ranges
    }else{
      snp.cor = rbind(snp.cor, cbind(co, randc))
      has.ctcf = append(has.ctcf, ov.region %in% peak.ranges)
    }
  }
}
#snp.cor = data.frame(neighbor = snp.cor[, 1], random = snp.cor[, 2])
melt.cor = data.frame(variable = factor(append(append(rep('no CTCF', sum(!has.ctcf)), rep('CTCF', sum(has.ctcf))), rep('random', length(has.ctcf)))),
                       value = append(append(snp.cor[!has.ctcf, 1], snp.cor[has.ctcf, 1]), snp.cor[, 2]))
p = ggplot(melt.cor) + geom_density(aes(x = value, y = ..density.., color = variable), size = 2, adjust = 1) +
  scale_x_continuous('Correlation') + scale_y_continuous('Density') + ggtitle(paste('Correlation between', mark, 'at SNPs and gene expression')) +
  scale_color_discrete('') + theme_bw() +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), title = element_text(size = 23),
        axis.title.x = element_text(size = 23), axis.title.y = element_text(size = 23), legend.text = element_text(size = 20))
ggsave(file.path(plotdir, paste(mark, '_gene_cor_at_snps.png', sep = '')), p, width = 13.6, height = 11.8)
q = ggplot(melt.cor) + geom_boxplot(aes(x = variable, y = value, fill = variable)) +
  #scale_x_continuous('Correlation') + scale_y_continuous('Density') + ggtitle(paste('Correlation between', mark, 'at SNPs and gene expression')) +
  scale_fill_discrete('') + theme_bw()
ggsave(file.path(plotdir, paste(mark, '_gene_cor_at_snps_boxplot.png', sep = '')), q, width = 13.6, height = 11.8)