rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# Regresses gene's expression on H3K4me3 and H3K27ac.

load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData') # Gene metadata
rna.deseq.dir = 'rawdata/geneCounts/rdata/repsComb/deseq/'
rna.dir = 'rawdata/geneCounts/rdata/repsComb'
peak.dir = 'rawdata/signal/combrep/countsAtPeaksBroad/repsComb/'

# Load the DESeq results for RNA and get a list of genes that are differential between several pairs of individuals.
rna.diff = get.diff.count(list.files(rna.deseq.dir, pattern = '*RZ_deseq.RData', full.names = T), 0.00001) > 10

# Get a matrix of average RPKMS for each gene (rows) and individual (cols)
# TODO: the gene metafile has only the longest transcript per gene. We need to consider all of them and compute the corresponding transcript lengths.
indivs = unique(as.character(sample.info(list.files(rna.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
nindiv = length(indivs)
rna.counts.dat = avg.counts(rna.dir, indivs, 'RZ_0.RData', meta = gene.meta) # This will normalize counts by the transcript length
rna.counts = rna.counts.dat$counts[rna.diff, ]
gene.meta = gene.meta[rna.diff, ]

# TSSs. We will overlap these with the H3K4me3 peaks. 
tss.meta = gene.meta 
pos = tss.meta$strand == '+'
tss.meta$end[pos] = tss.meta$start[pos]
tss.meta$start[!pos] = tss.meta$end[!pos]
tss.ranges = regions.to.ranges(tss.meta)

# Get the H3K4ME3 peak overlapping each differential TSS. These peaks give the "promoter domains".
load(list.files(peak.dir, 'SNYDER_HG19_GM12878_H3K4ME3_0.RData')) # Load any H3K4ME3 file. We just need the regions.
me3.regions = regions 

# For each TSS, find the overlapping H3K4ME3 peak (sincepeaks don't overlap, there should be at most one overlapping peak for each TSS)
me3.ov = findOverlaps(tss.ranges, regions.to.ranges(me3.regions), ignore.strand = T, select = 'first')
gene.meta = gene.meta[!is.na(me3.ov), ] # Select genes that do have an overlapping peak.
tss.meta = tss.meta[!is.na(me3.ov), ]
tss.ranges = regions.to.ranges(tss.meta)
rna.counts = rna.counts[!is.na(me3.ov), ]
me3.regions = me3.regions[me3.ov[!is.na(me3.ov)], ]
me3.ov = me3.ov[!is.na(me3.ov)]

co = array(NaN, dim = c(length(me3.ov), 1))
for(i in 1:length(me3.ov)){
  if(!is.na(me3.ov[i])){co[i] = cor(t(rna.counts[i, ]), t(me3.counts[i, ]))}
}

ac.diff = get.diff.count(list.files('rawdata/signal/combrep/countsAtPeaksBroad/repsComb/deseq/', pattern = '*H3K27AC_deseq.RData', full.names = T), 0.01) > 10
ac.counts.dat = avg.counts(peak.dir, indivs, 'H3K27AC_0.RData', meta = NULL)
ac.regions = ac.counts.dat$regions[ac.diff, ]
ac.counts = ac.counts.dat$counts[ac.diff, ]
ac.ranges = regions.to.ranges(ac.regions)

# For each gene, get the "proximal" and "distal" K27AC peaks.
ov.ranges = flank(tss.ranges, width = 500000, both = T)
tss.ov = findOverlaps(me3.ranges, ac.ranges, ignore.strand = T, select = 'all')
all.ov = findOverlaps(ov.ranges, ac.ranges, ignore.strand = T, select = 'all')

pval = array(1, dim=c(length(me3.ov), 1))
for(i in 1:length(pval)){
  pred = data.frame(rna = t(rna.counts[i, ]))
  colnames(pred) = c('rna')
  sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i]
  if(length(sel.prox) > 0){
    for(j in 1:length(sel.prox)){
      name = rownames(ac.counts)[sel.prox[j]]
      pred[[name]] = t(ac.counts[sel.prox[j], ])
      if(j == 1){null = paste('rna', name, sep = '~')
      }else{null = paste(null, name, sep = '+')}
    }
  }else{null = 'rna ~ 1'}
  sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
  if(length(sel.dist) > 0){
    alt = null
    ac.co = cor(t(rna.counts[i, ]), t(ac.counts[sel.dist, ]))
    sel.dist = sel.dist[order(ac.co)[1:min(nindiv - 3, length(sel.dist))]]
    for(j in 1:length(sel.dist)){
      name = rownames(ac.counts)[sel.dist[j]]
      pred[[name]] = t(ac.counts[sel.dist[j], ])
      alt = paste(alt, name, sep = '+')
    }
    fit.null = glm(null, family='gaussian', data = pred)
    fit.alt = glm(alt, family='gaussian', data = pred)
    pval[i] = 1-pchisq(fit.null$deviance-fit.alt$deviance,2)
    if(pval[i] < 0.1){
      print(pval[i])
      print(alt)
    }
  }
}
