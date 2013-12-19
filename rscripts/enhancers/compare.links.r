rm(list=ls())
library(GenomicRanges)
library(ggplot2)
source('utils/deseq.utils.r')
source('utils/sample.info.r')
#source('enhancers/get.as.ov.r')

# Computes overlaps between different sets of enhancer-gene links.

# Stam links
stam.assoc = read.table('../../rawdata/enhancers/external/stam/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8', header = F, sep = '\t', na.strings = '?')[, 4:8]
colnames(stam.assoc) = c('gene.name', 'chr', 'start', 'end', 'corr')
# Jason GM12878 links
j.assoc = read.table('../../rawdata/enhancers/external/jason/GM12878_links_hg19_merged.bed', header = F, sep = '\t', na.strings = '?')[, c(4,1,2,3,5)]
colnames(j.assoc) = c('gene.name', 'chr', 'start', 'end', 'corr')
# Jason all links 
j.all.assoc = read.table('../../rawdata/enhancers/external/jason/all_links_hg19_merged.bed', header = F, sep = '\t', na.strings = '?')[, c(4,1,2,3,5)]
colnames(j.all.assoc) = c('gene.name', 'chr', 'start', 'end', 'corr')

ext.assoc = rbind(stam.assoc, j.all.assoc)
ext.assoc$type = factor(append(rep('Thurman et al.', nrow(stam.assoc)), 
                               #rep('Ernst GM12878', nrow(j.assoc)), 
                                rep('Ernst et al.', nrow(j.all.assoc))))
                  
# Discovered links
link.files = c('../../rawdata/enhancers/merged_Mar13/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData')
types = c('ARS')

plotdir = '../../rawdata/enhancers/merged_Mar13/plots/'
outpref = 'enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_distToTSS' #gsub('.RData', '', basename(link.file))
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData') # gene metadata

dist.dat.all = NULL
ov.summary.all = NULL
for(i in 1:length(link.files)){
  ov.summary = NULL # Summary of overlaps between discovered set and other annotations
  all.ov = NULL
  
  links = new.env()
  load(link.files[i], links)
  assoc = links$ac.regions[links$links$region.idx, ]
  assoc = cbind(links$gene.meta$gene.name[links$links$gene.idx], assoc)
  colnames(assoc) = c('gene.name', 'chr', 'start', 'end')
  
  # Remove associations with non-annotated genes
  match.gene.idx = match(as.character(assoc$gene.name), gene.meta$gene.name)
  assoc = assoc[!is.na(match.gene.idx), ]
  gene.meta.tmp = gene.meta[match.gene.idx[!is.na(match.gene.idx)], ] # gene info for our associations, one line per association
  cat('# associations', length(match.gene.idx), ' # unique regions', length(unique(links$links$region.idx)), ' # unique genes', length(unique(links$links$gene.idx)), '\n')
  cat('# associations in annotations', nrow(assoc), '\n')
  cat('Avg length of distal elements', mean(assoc$end - assoc$start + 1), '\n')
  assoc$start = assoc$start - 200
  assoc$end = assoc$end + 200
  
  dist.dat = data.frame(distance = round((assoc$start + assoc$end) / 2 - links$tss[links$links$gene.idx]))
  is.neg = links$gene.meta$strand[links$links$gene.idx] == '-'
  dist.dat$distance[is.neg] = -dist.dat$distance[is.neg]
  dist.dat$type = rep(types[i], nrow(dist.dat))
  dist.dat.all = rbind(dist.dat.all, dist.dat)
  
  ##### Overlaps between my links and other sets of links
  ov = findOverlaps(regions.to.ranges(assoc), regions.to.ranges(ext.assoc), select = 'all', ignore.strand = T)
  ov.mat = data.frame(cbind(queryHits(ov), subjectHits(ov)))
  ov.mat$type = factor(ext.assoc$type[subjectHits(ov)])
  ov.mat = ov.mat[as.character(assoc$gene.name[ov.mat[, 1]]) == ext.assoc$gene.name[ov.mat[, 2]], ]
  for(t in levels(ov.mat$type)){
    ov.ind = unique(ov.mat[ov.mat$type == t, 1]) # unique regions overlapping external associations of type t
    cat('Overlapping', t, length(ov.ind), '(', length(ov.ind) * 100 / nrow(assoc), '% )\n')
    ov.summary = rbind(ov.summary, data.frame(hits = length(ov.ind), annotation = t))
    all.ov.tmp = array(F, dim = c(nrow(assoc), 1))
    all.ov.tmp[ov.ind] = T
    all.ov = cbind(all.ov, all.ov.tmp)
  }
  ov.ind = unique(ov.mat[,1])
  cat('Overlapping any link set', length(ov.ind), '(', length(ov.ind) * 100 / nrow(assoc), '% )\n')
  
  ##### dsQTLs
  ds = read.table('../../rawdata/enhancers/dsQTL/dsQtl_hg19_merged.bed', header = T, sep = '\t')
  colnames(ds) = c('chr', 'start', 'end', 'pos')
  # Extend the DHS sites a little bit. We want to find cases where the dsQTL SNP overlaps our predicted "enhancer"
  # and the dsQTL DHS is near the associated predicted gene.
  ds$start = ds$start - 1000 
  ds$end = ds$end + 1000
  
  # Overlaps between our predicted genes and the dsQTL DHSs
  ds.ov.regions = findOverlaps(regions.to.ranges(gene.meta.tmp), regions.to.ranges(ds), select = 'all', ignore.strand = T)
  ds.ov.regions = cbind(queryHits(ds.ov.regions), subjectHits(ds.ov.regions))
  cat('Overlapping a dsQTL DHS', length(unique(ds.ov.regions[, 1])), '(', length(unique(ds.ov.regions[, 1])) * 100 / nrow(assoc), '% )\n')
  
  # Overlaps between our predicted enhancers and the dsQTLs SNPs
  ds.ov.snps = findOverlaps(regions.to.ranges(assoc), snps.to.ranges(ds), select = 'all', ignore.strand = T)
  ds.ov.snps = cbind(queryHits(ds.ov.snps), subjectHits(ds.ov.snps))
  cat('Overlapping a dsQTL SNP', length(unique(ds.ov.snps[, 1])), '(', length(unique(ds.ov.snps[, 1])) * 100 / nrow(assoc), '% )\n')
  
  # Find common overlaps
  common = merge(ds.ov.regions, ds.ov.snps)
  cat('Overlapping dsQTLs', length(unique(common[, 1])), '(', length(unique(common[, 1])) * 100 / nrow(assoc), '% )\n')
  ov.summary = rbind(ov.summary, data.frame(hits = length(unique(common[, 1])), annotation = 'dsQTL'))
  
  ov.ind = unique(append(ov.ind, common[, 1]))
  all.ov.tmp = array(F, dim = c(nrow(assoc), 1))
  all.ov.tmp[common[, 1]] = T
  all.ov = cbind(all.ov, all.ov.tmp)
  cat('Overlapping any link set or dsQTL', length(ov.ind), '(', length(ov.ind) * 100 / nrow(assoc), '% )\n')
  
  ##### eQTLs
  eqtl.ov = unique(read.table('../../rawdata/signal/combrep/peakFiles/merged/LD/merged_0.9/SNYDER_HG19_H3K27AC_merged_LD0.9.txt', header = F, sep = '\t', na.strings = '?')[, c(1:4, 6)])
  colnames(eqtl.ov) = c('chr', 'start', 'end', 'pos', 'gene.name')
  gene.idx = match(as.character(eqtl.ov$gene.name), gene.meta$gene.name) 
  eqtl.ov = eqtl.ov[!is.na(gene.idx), ]
  gene.idx = gene.idx[!is.na(gene.idx)]
  
  # Find SNPs overlapping my predicted enhancers (or more precisely, SNPs in LD with these regions)
  ov = findOverlaps(regions.to.ranges(assoc), regions.to.ranges(eqtl.ov), select = 'all', ignore.strand = T)
  ov = data.frame(cbind(queryHits(ov), subjectHits(ov)))
  # Check if the corresponding gene is the same as the eQTL gene
  ov.sel = ov[as.character(assoc$gene.name[ov[, 1]]) == eqtl.ov$gene.name[ov[, 2]], ]
  cat('Overlapping eQTL links (exact gene)', length(unique(ov.sel[, 1])), '(', length(unique(ov.sel[, 1])) * 100 / nrow(assoc), '% )\n')
  
  gene.ov = findOverlaps(regions.to.ranges(gene.meta.tmp), regions.to.ranges(gene.meta[gene.idx, ]), select = 'all', ignore.strand = T)
  gene.ov = data.frame(cbind(queryHits(gene.ov), subjectHits(gene.ov)))
  common.ov = merge(ov, gene.ov)
  ov.ind = unique(append(ov.ind, common.ov[, 1]))
  ov.summary = rbind(ov.summary, data.frame(hits = length(unique(common.ov[, 1])), annotation = 'eQTL'))
  all.ov.tmp = array(F, dim = c(nrow(assoc), 1))
  all.ov.tmp[unique(common.ov[, 1])] = T
  all.ov = cbind(all.ov, all.ov.tmp)
  
  cat('Overlapping eQTL links (overlapping gene)', length(unique(common.ov[, 1])), '(', length(unique(common.ov[, 1])) * 100 / nrow(assoc), '% )\n')
  cat('Overlapping any link set or QTL', length(ov.ind), '(', length(ov.ind) * 100 / nrow(assoc), '% )\n')
  cat('Overlapping any link set or QTL or AS SNP', length(ov.ind), '(', length(ov.ind) * 100 / nrow(assoc), '% )\n')
  
  ov.summary = rbind(ov.summary, data.frame(hits = length(ov.ind), annotation = 'any'))
  ov.summary$type = rep(types[i], nrow(ov.summary))
  ov.summary.all = rbind(ov.summary.all, ov.summary)
}

col = c('darkblue', 'deepskyblue4')
names(col) = types
dist.dat.all$distance = dist.dat.all$distance / 1000
p = ggplot(dist.dat.all, group = type) + geom_density(aes(x = distance), size = 1.5) + 
  xlab('Distance enhancer - TSS (Kb)') + ylab('Density') + theme_bw() + 
  #scale_color_manual(values = col) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 15), legend.position = c(.2, .8))
ggsave(file = file.path(plotdir, paste(outpref, '_distToTSS.pdf', sep = '')), p, width = 5.2, height = 4.5)

q = ggplot(ov.summary.all) + geom_bar(aes(x = annotation, y = hits), position = 'dodge') +
  xlab('Overlapping annotation') + ylab('# associations') + scale_fill_manual(values = col) + theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = -60, hjust = 0, vjust = 1), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 15), legend.title = element_blank())
ggsave(file = file.path(plotdir, paste(outpref, '_annotOv.pdf', sep = '')), q, width = 6.5, height = 5.6)


# ### Associations supported by AS SNPs for H3K27AC and AS genes
# snp.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
# as.files = list.files(snp.dir, pattern = 'H3K27AC_rep.RData$', full.names = T)
# indivs = as.character(sample.info(as.files, '.RData')$indiv)
# gene.files = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps/qvals', 
#                        paste('SNYDER_HG19', indivs, 'RZ_rep.RData', sep = '_'))
# good.files = file.exists(gene.files)
# 
# as.ov = get.as.ov(assoc, as.files[good.files], gene.files[good.files], gene.compare = T)
# as.ov.ind = apply(as.ov$same.dir, 1, function(x) any(x)) | apply(as.ov$diff.dir, 1, function(x) any(x))
# ov.summary = rbind(ov.summary, data.frame(hits = sum(as.ov.ind), annotation = 'AS-H3K27AC-gene'))
# ov.ind = unique(append(ov.ind, which(as.ov.ind)))
# all.ov.tmp = array(F, dim = c(nrow(assoc), 1))
# all.ov.tmp[as.ov.ind] = T
# all.ov = cbind(all.ov, all.ov.tmp)
# 
# ###### Associations supported by AS SNPs for H3K4ME1 and AS genes
# as.files = list.files(snp.dir, pattern = 'H3K4ME1_rep.RData$', full.names = T)
# indivs = as.character(sample.info(as.files, '.RData')$indiv)
# gene.files = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps/qvals', 
#                        paste('SNYDER_HG19', indivs, 'RZ_rep.RData', sep = '_'))
# good.files = file.exists(gene.files)
# 
# as.ov = get.as.ov(assoc, as.files[good.files], gene.files[good.files], gene.compare = T)
# as.ov.ind = apply(as.ov$same.dir, 1, function(x) any(x)) | apply(as.ov$diff.dir, 1, function(x) any(x))
# ov.summary = rbind(ov.summary, data.frame(hits = sum(as.ov.ind), annotation = 'AS-H3K4ME1-gene'))
# ov.ind = unique(append(ov.ind, which(as.ov.ind)))
# all.ov.tmp = array(F, dim = c(nrow(assoc), 1))
# all.ov.tmp[as.ov.ind] = T
# all.ov = cbind(all.ov, all.ov.tmp)


