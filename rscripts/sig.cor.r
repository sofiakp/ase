rm(list=ls())
library(GenomicRanges)
library(preprocessCore)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
library(DESeq)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

registerDoMC(10)

# Dir with counts for RNA. Should contain files SNYDER_HG19_<indiv>_<mark>_0.RData with read counts for genes.
# Each count matrix should have one column per replicate. The file should also have a vector size.factors with
# the size factor (library size) for each column.
rna.dir = '../../rawdata/geneCounts/rdata/repsComb'
peak.dir = '../../rawdata/signal/mergedInputs/combrep/extractSignal/fc/avgSig/rdata'
peak.dir = '../../rawdata/signal/rep/countsAtPeaksBroad/repsComb'

load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.trans.RData')
trans$pos = trans$start

# Get a matrix of RPKMS for each gene (rows) and individual (cols)
indivs = unique(as.character(sample.info(list.files(rna.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
nindiv = length(indivs)
rna.counts.dat = avg.counts(rna.dir, indivs, 'RZ_0.RData', meta = gene.meta, len.norm = F)
rna.counts = asinh(rna.counts.dat$counts)
rna.means = rowMeans(rna.counts, na.rm = T)
rna.sds = rowSds(rna.counts, na.rm = T)
rna.cv = rna.sds / rna.means
expr.genes = which(rna.means > asinh(1) & !is.na(rna.cv) & rna.cv < quantile(rna.cv, 0.5, na.rm = T))
rna.counts = normalize.quantiles(rna.counts)
trans = trans[trans$gene.idx %in% expr.genes, ] # Remove TSSs for lowly-expressed genes. Leave gene.meta unchanged for now.
#indivs[indivs=='SNYDER'] = 'MS1'
colnames(rna.counts) = indivs
cat('# genes:', nrow(rna.counts), '\n')

# Get the histone mark signal at variable peaks
# load(file.path(peak.dir, 'SNYDER_HG19_all_reg_H3K27AC_qn.RData'))
# sig.counts = counts[, colnames(counts) %in% colnames(rna.counts)]

# Get raw counts
regions = read.bed('../../rawdata/signal/combrep/peakFiles/merged/SNYDER_HG19_H3K27AC_merged.bed.gz')
sig.counts = NULL
for(i in 1:length(indivs)){
  load(file.path(peak.dir, paste('SNYDER_HG19', indivs[i], 'H3K27AC_0.RData', sep = '_')))
  sig.counts = cbind(sig.counts, rowSums(counts))  
}
#sf = estimateSizeFactorsForMatrix(sig.counts)
#sig.counts = t(apply(sig.counts, 1, function(x) x / sf))
sig.counts = normalize.quantiles(sig.counts)
colnames(sig.counts) = indivs

row.means = rowMeans(sig.counts)
row.sds = rowSds(sig.counts)
cvs = row.sds / row.means
good.rows = !is.na(cvs) #& row.means > asinh(0.2) & cvs > quantile(cvs, 0.2, na.rm = T)
regions$start = regions$start - 200
regions$end = regions$end + 200
sig.regions = regions[good.rows, ]
sig.counts = sig.counts[good.rows, ]
stopifnot(all(colnames(sig.counts) == colnames(rna.counts)))

# Get the peak overlapping each TSS
ov = findOverlaps(snps.to.ranges(trans), regions.to.ranges(sig.regions), select = 'all', ignore.strand = T)
# Col 1 - gene idx, Col 2 - me3 peak idx
ov.mat = unique(cbind(trans$gene.idx[queryHits(ov)], subjectHits(ov))) 
nov = table(ov.mat[,1]) # number of me3 peaks overlapping each gene's TSSs
# Remove genes whose TSSs overlap more than one H3K4me3 peak
sel.genes = as.numeric(names(nov[nov == 1]))
ov.mat = ov.mat[ov.mat[,1] %in% sel.genes, ]

cat('# unique genes with overlap', length(nov), '\n')
cat('# genes with exactly one TSS-overlapping peak:', sum(nov == 1), '\n')

sig.cor = cor.par(rna.counts[ov.mat[,1], ], sig.counts[ov.mat[,2], ], nchunks = 10, method = 'pearson')
cat('Mean,median abs cor:', mean(abs(sig.cor)), ' ', median(abs(sig.cor)), '\n')