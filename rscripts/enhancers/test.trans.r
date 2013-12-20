rm(list=ls())
library(GenomicRanges)
library(preprocessCore)
library(foreach)
library(doMC)
library(matrixStats)
library(ggplot2)
source('utils/deseq.utils.r')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

make.cor.plot = function(cor.dat, title){
  return(ggplot(cor.dat) + geom_density(aes(x = cor, color = is.prox), size = 1.5) + 
    theme_bw() + xlab('Correlation') + ylab('Density') + scale_color_discrete('') + ggtitle(title) +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 15), plot.title = element_text(size = 16), legend.position = c(0.2, 0.8)))
}
get.cor = function(rna.counts, ac.counts, accoc){
  cor.rand = cor.par(ac.counts[sample(1:nrow(ac.counts), sample.size),], rna.counts[sample(1:nrow(rna.counts), sample.size), ], method = 'pearson')
  cor.prox = cor.par(ac.counts[assoc$ac.idx[is.prox],], rna.counts[assoc$gene.idx[is.prox], ], nchunks = 1, method = 'pearson')
  cor.dist = cor.par(ac.counts[assoc$ac.idx[!is.prox],], rna.counts[assoc$gene.idx[!is.prox], ], nchunks = 1, method = 'pearson')
  cor.all = cor.par(ac.counts[assoc$ac.idx,], rna.counts[assoc$gene.idx, ], nchunks = 1, method = 'pearson')
  cat('proximal correlation, mean =', mean(cor.prox, na.rm = T), 'median =', median(cor.prox, na.rm = T), '\n')
  cat('distal correlation, mean =', mean(cor.dist, na.rm = T), 'median =', median(cor.dist, na.rm = T), '\n')
  cat('average correlation, mean =', mean(cor.dist, na.rm = T), 'median =', median(cor.all, na.rm = T), '\n')
  cat(sprintf('Diff between proximal and distal P = %g\n', wilcox.test(cor.prox, cor.dist)$p.value))
  cat(sprintf('Diff between proximal and random P = %g\n', wilcox.test(cor.prox, cor.rand)$p.value))
  cat(sprintf('Diff between distal and random P = %g\n', wilcox.test(cor.rand, cor.dist)$p.value))
  cat(sprintf('Diff between all and random P = %g\n', wilcox.test(cor.rand, cor.all)$p.value))
  cor.dat = data.frame(cor = append(cor.rand, append(cor.prox, cor.dist)), 
                       is.prox = append(rep('random', sample.size), append(rep('Proximal', sum(is.prox)), rep('Distal', sum(!is.prox)))))
  p1 = make.cor.plot(cor.dat, 'Correlation with expression at TSS')
  cor.dat = data.frame(cor = append(cor.rand, cor.all), 
                       is.prox = append(rep('Random', sample.size), rep('True', length(cor.all))))
  p2 = make.cor.plot(cor.dat, 'Correlation with expression at TSS')
  return(list(p1 = p1, p2 = p2))
}

registerDoMC(2)

# Dir with counts for RNA. Should contain files SNYDER_HG19_<indiv>_<mark>_0.RData with read counts for genes.
# Each count matrix should have one column per replicate. The file should also have a vector size.factors with
# the size factor (library size) for each column.
rna.dir = '../../rawdata/geneCounts/rdata/repsComb'
peak.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata' 
outdir = file.path('../../rawdata/enhancers/external/jason/plots')
if(!file.exists(outdir)) dir.create(outdir)
mark = 'H3K27AC' # signal in distal elements to correlate with expression
outpref = 'all_jason_GM12878_'

# Load TSSs. gene.meta has gene information and trans has positions of all TSSs as well as the index of 
# the corresponding gene in gene.meta.
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.trans.RData')
trans$pos = trans$start

# Get a matrix of average RPKMS for each gene (rows) and individual (cols)
indivs = unique(as.character(sample.info(list.files(rna.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
indivs = indivs[!(indivs %in% c('GM12878', 'GM19240'))] # REMOVE daughters
nindiv = length(indivs)
rna.counts.dat = avg.counts(rna.dir, indivs, 'RZ_0.RData', meta = gene.meta, len.norm = F)
rna.counts = asinh(rna.counts.dat$counts)
rna.means = rowMeans(rna.counts, na.rm = T)
rna.sds = rowSds(rna.counts, na.rm = T)
rna.cv = rna.sds / rna.means
expr.genes = which(rna.means > asinh(0.2) & !is.na(rna.cv) & rna.cv > quantile(rna.cv, 0.2, na.rm = T))
#rna.counts = normalize.quantiles(scale(rna.counts, center = T, scale = F))
rna.counts = normalize.quantiles(rna.counts)
trans = trans[trans$gene.idx %in% expr.genes, ] # Remove TSSs for lowly-expressed genes. Leave gene.meta unchanged for now.
colnames(rna.counts) = indivs
cat('# genes:', nrow(rna.counts), ', selected', length(expr.genes), '\n')

# Get the H3K4ME3 peak overlapping each differential TSS. These peaks give the "promoter domains".
load(file.path(peak.dir, 'SNYDER_HG19_H3K4ME3_qn.RData'))
row.means = rowMeans(counts)
row.sds = rowSds(counts)
cvs = row.sds / row.means
good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.2, na.rm = T)
me.regions = regions[good.rows, ]
me.counts = counts[good.rows, colnames(counts) %in% colnames(rna.counts)]
stopifnot(all(colnames(me.counts) == colnames(rna.counts)))

# region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/', 
#                         paste('SNYDER_HG19_H3K4ME3_merged.bed.gz', sep = '_'))
# signal.files = list.files(peak.dir, pattern = paste('SNYDER_HG19', 'H3K4ME3', 'merged_AT_.*txt', sep = '_'), full.names=T)
# indivs = unique(gsub(paste('SNYDER_HG19_', 'H3K4ME3', '_merged_AT_SNYDER_HG19_|_', 'H3K4ME3', '.*.txt', sep = ''), '', basename(signal.files)))
# counts.dat = load.avg.sig.data(region.file, signal.files, indivs)
# counts = asinh(counts.dat$signal)
# regions = counts.dat$regions
# row.means = rowMeans(counts)
# row.sds = rowSds(counts)
# cvs = row.sds / row.means
# good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.4, na.rm = T)
# me.regions = regions[good.rows, ]
# me.counts = counts[good.rows, colnames(counts) %in% colnames(rna.counts)]
# oldcol = colnames(me.counts)
# me.counts = normalize.quantiles(scale(me.counts))
# colnames(me.counts) = oldcol
# stopifnot(all(colnames(me.counts) == colnames(rna.counts))) # Make sure the columns are in the same order

me.ov = findOverlaps(snps.to.ranges(trans), regions.to.ranges(me.regions), select = 'all', ignore.strand = T)
me.ov.mat = unique(cbind(trans$gene.idx[queryHits(me.ov)], subjectHits(me.ov))) # Col 1 - gene idx, Col 2 - me3 peak idx
nov = table(me.ov.mat[,1]) # number of me3 peaks overlapping each gene's TSSs
# Remove genes whose TSSs overlap more than one H3K4me3 peak
sel.genes = as.numeric(names(nov[nov == 1]))
me.ov.mat = me.ov.mat[me.ov.mat[,1] %in% sel.genes, ]

# sel.genes = unique(me.ov.mat[,1])
stopifnot(all(sel.genes %in% expr.genes))
cat('# unique active genes (overlapping H3K4me3 peaks):', length(nov), '\n')
cat('# genes with exactly one TSS-overlapping peak:', sum(nov == 1), '\n')

me.match = list() # For each gene, the list of all H3K4me3 peaks overlapping the gene's TSSs
for(i in 1:length(sel.genes)) me.match[[i]] = me.ov.mat[me.ov.mat[,1] == sel.genes[i], 2]
all.gene.meta = gene.meta # Keep the original annotations
gene.meta = gene.meta[sel.genes, ] # Select genes that do have an overlapping peak.
rna.counts = rna.counts[sel.genes, ]
me.match.counts = me.counts[unlist(lapply(me.match, function(x) x[1])), ] # signal of one of the gene's me3 peaks

# Load H3K27ac signal and peaks
load(file.path(peak.dir, paste('SNYDER_HG19', mark, 'qn.RData', sep = '_')))
row.means = rowMeans(counts)
row.sds = rowSds(counts)
cvs = row.sds / row.means
good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.2, na.rm = T)
ac.regions = regions[good.rows, ]
ac.counts = counts[good.rows, colnames(counts) %in% colnames(rna.counts)]
stopifnot(all(colnames(ac.counts) == colnames(rna.counts)))

# region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/', 
#                         paste('SNYDER_HG19', mark, 'merged.bed.gz', sep = '_'))
# signal.files = list.files(peak.dir, pattern = paste('SNYDER_HG19', mark, 'merged_AT_.*txt', sep = '_'), full.names=T)
# indivs = unique(gsub(paste('SNYDER_HG19_', mark, '_merged_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
# counts.dat = load.avg.sig.data(region.file, signal.files, indivs)
# counts = asinh(counts.dat$signal)
# regions = counts.dat$regions
# row.means = rowMeans(counts)
# row.sds = rowSds(counts)
# cvs = row.sds / row.means
# good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.4, na.rm = T)
# ac.regions = regions[good.rows, ]
# ac.counts = counts[good.rows, colnames(counts) %in% colnames(rna.counts)]
# oldcol = colnames(ac.counts)
# ac.counts = normalize.quantiles(scale(ac.counts))
# colnames(ac.counts) = oldcol
# stopifnot(all(colnames(ac.counts) == colnames(rna.counts)))

# Load enhancer-gene associations
# Stam links
# assoc = read.table('../../rawdata/enhancers/external/stam/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8', header = F, sep = '\t')[, 4:8]
# Jason links
assoc = read.table('../../rawdata/enhancers/external/jason/all_links_hg19_merged.bed', header = F, sep = '\t')[, c(4,1,2,3,5)]
# Discovered links
# links = new.env()
# load('../../rawdata/enhancers/rdata/enhancer_pred_l1_H3K27AC.RData', links)
# assoc = links$dist.regions[links$links$region.idx, ]
# assoc = cbind(links$gene.meta$gene.name[links$links$gene.idx], assoc)
# assoc = cbind(assoc, links$links$pval)
colnames(assoc) = c('gene.name', 'chr', 'start', 'end', 'corr')
# For each distal element, try to find an overlapping H3K27ac peak
assoc.ov = findOverlaps(regions.to.ranges(assoc), regions.to.ranges(ac.regions), select = 'first', ignore.strand = T)
gene.ov = match(as.character(assoc$gene.name), gene.meta$gene.name)
sel = !is.na(assoc.ov) & !is.na(gene.ov) # get associations for which the gene is active and the "distal" element is also active (overlapping H3K27ac peak)
cat('# associations:', nrow(assoc), '\n')
cat('overlapping', mark, 'peak:', sum(!is.na(assoc.ov)), '(', sum(!is.na(assoc.ov)) * 100 / nrow(assoc),'% )\n')
in.annot = sum(as.character(assoc$gene.name) %in% all.gene.meta$gene.name)
cat('overlapping annotated gene:', in.annot, '(', in.annot * 100 / nrow(assoc),'% )\n')
cat('overlapping an active gene:', sum(!is.na(gene.ov)), '(', sum(!is.na(gene.ov)) * 100 / nrow(assoc),'% )\n')
cat('overlapping', mark, 'peak and an active gene:', sum(sel), '(', sum(sel) * 100 / nrow(assoc),'% )\n')
assoc = assoc[sel, ]
assoc$ac.idx = assoc.ov[sel] # index of associated region in ac.counts
assoc$gene.idx = gene.ov[sel] # index of associated gene in gene.meta and rna.counts
assoc = assoc[!duplicated(cbind(assoc$ac.idx, assoc$gene.idx)), ]

# Find the "distal" elements that overlap a proximal H3K4me3 peak
ov = findOverlaps(regions.to.ranges(assoc), regions.to.ranges(me.regions[me.ov.mat[, 2], ]), select = 'first', ignore.strand = T)
is.prox = !is.na(ov)
cat('# proximal associations:', sum(is.prox), '(', sum(is.prox) * 100 / nrow(assoc),'% )\n')
sample.size = min(c(5000, nrow(ac.counts), nrow(rna.counts)))

cat('Correlations with RNA\n')
p = get.cor(rna.counts, ac.counts, assoc)
ggsave(file.path(outdir, paste(outpref, mark, '_cor_with_RZ.pdf', sep = '')), p$p1, width = 6.5, height = 5.6)
ggsave(file.path(outdir, paste(outpref, mark, '_cor_with_RZ_merged.pdf', sep = '')), p$p2, width = 6.5, height = 5.6)

cat('Correlations with H3K4me3 at TSSs\n')
q = get.cor(me.match.counts, ac.counts, assoc)
ggsave(file.path(outdir, paste(outpref, mark, '_cor_with_H3K4ME3.pdf', sep = '')), q$p1, width = 6.5, height = 5.6)
ggsave(file.path(outdir, paste(outpref, mark, '_cor_with_H3K4ME3_merged.pdf', sep = '')), q$p2, width = 6.5, height = 5.6)

# cor.prox = cor.par(ac.counts[assoc$ac.idx[is.prox],], me.match.counts[assoc$gene.idx[is.prox], ], nchunks = 1, method = 'pearson')
# cor.dist = cor.par(ac.counts[assoc$ac.idx[!is.prox],], me.match.counts[assoc$gene.idx[!is.prox], ], nchunks = 1, method = 'pearson')
# cat('proximal correlation, mean =', mean(cor.prox, na.rm = T), 'median =', median(cor.prox, na.rm = T), '\n')
# cat('distal correlation, mean =', mean(cor.dist, na.rm = T), 'median =', median(cor.dist, na.rm = T), '\n')
# cat(sprintf('Diff between proximal and distal P = %g\n', wilcox.test(cor.prox, cor.dist)$p.value))
# cat(sprintf('Diff between proximal and random P = %g\n', wilcox.test(cor.prox, cor.rand)$p.value))
# cat(sprintf('Diff between distal and random P = %g\n', wilcox.test(cor.rand, cor.dist)$p.value))
# cor.dat = data.frame(cor = append(cor.rand, append(cor.prox, cor.dist)), 
#                      is.prox = append(rep('random', sample.size), append(rep('proximal', sum(is.prox)), rep('distal', sum(!is.prox)))))
# p2 = make.cor.plot(cor.dat, 'Correlation with H3K4me3 at TSS')
# ggsave(file.path(outdir, paste(outpref, mark, '_cor_with_H3K4ME3.pdf', sep = '')), p2, width = 6.5, height = 5.6)

# is.prox = array(F, dim = c(nrow(assoc), 1))
# for(i in 1:length(is.prox)){
#   ov = findOverlaps(regions.to.ranges(assoc[i, ]), regions.to.ranges(me.regions[me.match[[assoc$gene.idx[i]]], ]), select = 'first', ignore.strand = T)
#   is.prox[i] = !is.na(ov[1])
# }
