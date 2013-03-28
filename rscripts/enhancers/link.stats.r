rm(list=ls())
library(GenomicRanges)
library(ggplot2)
source('utils/deseq.utils.r')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/enhancers/get.as.ov.r'))

# Computes statistics about the discovered links

plotdir = '../../rawdata/enhancers/plots/'

# Make plots for ARS hits first
link.file = '../../rawdata/enhancers/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData'
load(link.file)
indivs = colnames(ac.counts)
indivs[indivs == 'SNYDER'] = 'MS1'
outpref = gsub('.RData', '', basename(link.file))

assoc = ac.regions[links$region.idx, ]
assoc = cbind(gene.meta$gene.name[links$gene.idx], assoc)
colnames(assoc) = c('gene.name', 'chr', 'start', 'end')

dist.dat = data.frame(distance = round((assoc$start + assoc$end) / 2 - tss[links$gene.idx]),
                      type = rep('ARS', nrow(links)))
is.neg = gene.meta$strand[links$gene.idx] == '-'
dist.dat$distance[is.neg] = -dist.dat$distance[is.neg]

max.counts = data.frame(as.matrix(table(links$max.line))) # Count how many times each individual is the maximum
colnames(max.counts) = c('hits')
max.counts$indiv = indivs[as.numeric(rownames(max.counts))]
p1 = ggplot(max.counts) + geom_bar(aes(x = indiv, y = hits, stat = 'identity')) + 
  xlab('') + ylab('# times individual has max ARS') + theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 1, hjust = 0), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
ggsave(file.path(plotdir, paste(outpref, '_hitsPerIndiv.png', sep = '')), p1, width = 6.5, height = 5.6)

sum.sign = data.frame(n = table(rowSums(ars.score.norm[sel.assoc, ] > 0.5))) # For each significant association how many individuals are "outliers"?
colnames(sum.sign) = c('n', 'hits')
p2 = ggplot(sum.sign) + geom_bar(aes(x = n, y = hits), stat = 'identity') + 
  xlab('# individuals with ARS > max(ARS)/2') + ylab('# links') + theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
ggsave(file.path(plotdir, paste(outpref, '_indivPerHit.png', sep = '')), p2, width = 6.5, height = 5.6)

# Make common plots
link.file = '../../rawdata/enhancers/rdata/enhancer_coef_elastic0.5_100kb_asinh0.2_cv0.2_H3K27AC_links.RData'
load(link.file)
outpref = 'enhancer_coef_elastic0.5_vs_ars_100kb_asinh0.2_cv0.2_H3K27AC'

assoc.tmp = ac.regions[links$region.idx, ]
assoc.tmp = cbind(gene.meta$gene.name[links$gene.idx], assoc.tmp)
colnames(assoc.tmp) = c('gene.name', 'chr', 'start', 'end')

dist.dat.tmp = data.frame(distance = round((assoc.tmp$start + assoc.tmp$end) / 2) - tss[links$gene.idx],
                      type = rep('regression', nrow(links)))
is.neg = gene.meta$strand[links$gene.idx] == '-'
dist.dat.tmp$distance[is.neg] = -dist.dat.tmp$distance[is.neg]
dist.dat = rbind(dist.dat, dist.dat.tmp)

dist.dat$distance = dist.dat$distance / 1000
p4 = ggplot(dist.dat) + geom_density(aes(x = distance, linetype = type)) +
  xlab('Distance enhancer - TSS (Kb)') + ylab('Density') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = c(.2, .8))
ggsave(file = file.path(plotdir, paste(outpref, '_distToTSS.png', sep = '')), p4, width = 6.5, height = 5.6)

ov = findOverlaps(regions.to.ranges(assoc), regions.to.ranges(assoc.tmp), select = 'all', ignore.strand = T)
common = sum(assoc$gene.name[queryHits(ov)] == assoc.tmp$gene.name[subjectHits(ov)])
nhits = data.frame(hits = c(nrow(assoc), nrow(assoc.tmp), common), 
                   type = ordered(factor(c('ARS', 'regression', 'common'), levels = c('ARS', 'regression', 'common'))))
p3 = ggplot(nhits) + geom_bar(aes(x = type, y = hits), stat = 'identity') + 
  xlab('') + ylab('Number of links') + theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
ggsave(file.path(plotdir, paste(outpref, '_nhits.png', sep = '')), p3, width = 6.5, height = 5.6)
