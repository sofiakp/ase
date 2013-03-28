rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# File with bins around TSSs created with load.diff.regress.data.r
bin.file = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig/81bins_diff/rdata/gencode.v13.annotation.noM.genes.81bins_diff.RData'
# Directory with signal at peaks. Will be used to compute extra features with the overall signal at peaks 
# in the vicinity of each gene.
signal.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/'
rna.dir = '../../rawdata/geneCounts/rdata/repsComb' # RNA counts
outpref = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig/81bins_diff/rdata/gencode.v13.annotation.noM.genes.81bins_diff'
  
load(bin.file)
uniq.indivs = setdiff(unique(indivs), 'GM18486') # Remove individuals that don't have expression data
nindivs = length(uniq.indivs) 
sel = indivs %in% uniq.indivs
marks = marks[sel]
indivs = indivs[sel]
counts = counts[, sel]
cols = cols[sel]
nmarks = length(unique(marks))
uniq.bins = as.character(sort(unique(bins$bin))) # Unique bins
nbins = length(uniq.bins) # Number of bins for each gene
stopifnot(nrow(bins) %% nbins == 0) # Size of matrix is multiple of bin number
ngenes = length(unique(bins$gene.name)) # Number of unique genes

# The bins are in coordinate order, so bins for genes on the - strand are reversed. 
# Reverse them, so that the bins are in the same order for all bins
# for(i in seq(1, nrow(bins), nbins)){
#   start = i * nbins + 1
#   end = (i + 1) * nbins
#   stopifnot(all(bins$gene.name[start:end] == bins$gene.name[start]))
#   if(bins$strand[start] == '-'){
#     new.bins = bins[seq(end, start, -1), ]
#     bins[start:end, ] = new.bins
#     counts[start:end, ] = counts[seq(end, start, -1), ]
#   }
# }

# Input count matrix has one column per individual per mark and one row per bin per gene.
# We want one row per gene and one column per bin per mark for all pairs of individuals,
# so a row will be 'indiv1_bin1_mark1, indiv1_bin1_mark2,...,indiv1_binN_markM, indiv2..."
npairs = nindivs * (nindivs - 1) / 2 # Number of pairs of individuals
new.counts = NULL # array(0, dim = c(ngenes * npairs, new.ncol * 2))
cat('Reading pairs\n', file = stderr())
for(i in 1:(nindivs - 1)){
  cat(uniq.indivs[i], '\n', file = stderr())
  for(j in (i + 1):nindivs){
    sel.cols = append(which(indivs == uniq.indivs[i]), which(indivs == uniq.indivs[j]))
    pair.counts = NULL
    for(b in 1:length(uniq.bins)){
      tmp.counts = counts[bins$bin == uniq.bins[b], sel.cols] # all rows for this type of bin, so one row per gene
      colnames(tmp.counts) = paste(rep(c('ind1', 'ind2'), times = c(nmarks, nmarks)), marks[sel.cols], uniq.bins[b], sep = '_')
      if(is.null(pair.counts)){pair.counts = tmp.counts
      }else{pair.counts = cBind(pair.counts, tmp.counts)}
    }
    if(is.null(new.counts)){new.counts = pair.counts
    }else{new.counts = rBind(new.counts, pair.counts)}
  }
}

# for(i in 1:1){ #(nindivs - 1)){
#     col1 = grep(indivs[i], colnames(counts))
#     sig1 = matrix(t(counts[, col1]), ncol = new.ncol, byrow = T)
#     #sig1[sig1 < 1] = 0
#     sig1 = Matrix(sig1)
#   for(j in (i + 1):nindivs){
#         col2 = grep(indivs[j], colnames(counts))
#         sig2 = matrix(t(counts[, col2]), ncol = new.ncol, byrow = T)
#         #sig2[sig2 < 1] = 0
#         sig2 = Matrix(sig2)
#         if(is.null(new.counts)){new.counts = cBind(sig1, sig2)
#         }else{new.counts = rBind(new.counts, cBind(sig1, sig2))}
#   }
# }

# Get the vicinity of each gene (i.e. union of all bins for the gene)
tss.idx = bins$bin == 'TSS_1'
tts.idx = bins$bin == 'TTS_1'
starts = bins$start[tss.idx]
strands = bins$strand[tss.idx]
chroms = bins$chr[tss.idx]
genes = bins$gene.name[tss.idx]
neg = strands == '-'
starts[neg] = bins$start[tts.idx][neg]
ends = bins$end[tts.idx]
ends[neg] = bins$end[tss.idx][neg]

gene.regions = data.frame(chr = chroms, start = starts, end = ends)

# For each gene, we will compute the aggregate signal of peaks in its vicinity
marks = c('H3K27AC')#, 'H3K4ME3', 'H3K4ME1', 'CTCF')
aggr.sig = NULL
cat('Reading aggregate signal\n', file = stderr())
for(m in marks){
  cat(m, '\n', file = stderr())
  # Load signal at merged peaks for the mark
  peaks = new.env()
  load(file.path(signal.dir, paste('SNYDER_HG19', m, 'qn.RData', sep = '_')), peaks)
  
  ov = findOverlaps(regions.to.ranges(gene.regions), regions.to.ranges(peaks$regions), select = 'all', ignore.strand = T)
  ov.mat = cbind(queryHits(ov), subjectHits(ov))
  merged.sig = array(0, dim = c(nrow(gene.regions), ncol(peaks$counts)))
  ov.peaks = unique(ov.mat[, 1])
  for(i in 1:length(ov.peaks)){
    idx = ov.mat[, 1] == ov.peaks[i]
    merged.sig[ov.peaks[i], ] = colSums(peaks$counts[idx, ])
  }
  pair.sig = NULL
  for(i in 1:(nindivs - 1)){
    for(j in (i + 1):nindivs){
      tmp.sig = cbind(merged.sig[, colnames(peaks$counts) == uniq.indivs[i]],
              merged.sig[, colnames(peaks$counts) == uniq.indivs[j]])
      pair.sig = rbind(pair.sig, tmp.sig)
    }
  }
  colnames(pair.sig) = paste(c('ind1', 'ind2'), m, c('aggr_peaks', 'aggr_peaks'), sep = '_')
  if(is.null(aggr.sig)){aggr.sig = pair.sig  
  }else{aggr.sig = cbind(aggr.sig, pair.sig)}
}
aggr.sig[aggr.sig < 1] = 0
aggr.sig[is.na(aggr.sig)] = 0
aggr.sig = Matrix(aggr.sig)
pair.counts = cBind(aggr.sig, new.counts)
save(pair.counts, bins, genes, file = paste(outpref, '_designMat.RData', sep = ''))

# Read RNA counts and compute log-fold-changes
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
rna.counts = avg.counts(rna.dir, uniq.indivs, 'RZ_0.RData', meta = gene.meta, len.norm = F)
rna.counts = rna.counts + 1
rna.counts = normalize.quantiles(rna.counts$counts)
rna.counts = log(rna.counts)
match.idx = match(genes, gene.meta$gene.name)
stopifnot(all(!is.na(match.idx)))
expr.fc = NULL
for(i in 1:(nindivs - 1)){
  for(j in (i + 1):nindivs){
    tmp.fc = rna.counts[match.idx, i] - rna.counts[match.idx, j]
    expr.fc = rBind(expr.fc, tmp.fc)
  }
}
save(expr.fc, genes, file = paste(outpref, '_diffExpr.RData', sep = ''))