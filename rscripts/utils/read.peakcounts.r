rm(list=ls())
library(DESeq)
library(GenomicRanges)
source('utils/sample.info.r')
source('utils/deseq.utils.r')

# Reads files with counts in regions (eg. peaks). Replicates are assumed to be in different files, eg.
# SNYDER_HG19_GM12878_H3K4ME3_[123]_.*
# All replicate files are combined into one matrix with one column per replicate and written
# in SNYDER_HG19_indiv_mark_0.RData.
# Optionally, size factors for each replicate are computed with DESeq. The size factors 
# are based on ALL datasets for the same mark, eg. SNYDER_HG19_.*_H3K4ME3.*

indir = '../../rawdata/genomeGrid/hg19_w10k/rep/counts_newNorm/' #'../../rawdata/transcriptomes/rep/counts_newNorm/' #'../../rawdata/signal/rep/countsAtPeaksBroad/merged_Mar13/'
outdir = file.path(indir, 'repsComb')
if(!file.exists(outdir)) dir.create(outdir)

get.sf = F # Get size factors with DESeq using all datasets for the same mark.
is.gene = F # T if you're reading counts at gene bodies. In this case, the 4th column is the gene name

filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*.bed', sep = ''), full.name = F)
fsplit = strsplit(gsub('.*_AT_', '', filenames), '_')
marks = array(0, dim = c(length(filenames), 1))
indivs = array(0, dim = c(length(filenames), 1))
samples = array(0, dim = c(length(filenames), 1))
for(i in 1:length(filenames)){
  marks[i] = fsplit[[i]][4]
  indivs[i] = fsplit[[i]][3]
  samples[i] = paste(fsplit[[i]][3:5], collapse = '_') # indiv_mark_rep - replicate ids, will also be used as headers
}
uniq.marks = unique(marks)

# ENCODE blacklisted regions: These will be completely removed from the output files
bad.ranges = regions.to.ranges(read.bed('../../rawdata/genomes_local/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed'))

for(i in 1:length(uniq.marks)){
  uniq.indivs = unique(indivs[marks == uniq.marks[i]])
  if(length(uniq.indivs) < 1) next
  
  # For each individual with that mark
  for(k in 1:length(uniq.indivs)){
    sel = grep(paste('_AT_SNYDER_HG19', uniq.indivs[k], uniq.marks[i], sep = '_'), filenames) # Select all replicates for that individual and mark
    rep.filenames = filenames[sel]
    print(rep.filenames)
    nfiles = length(rep.filenames)
    outfile = file.path(outdir, paste('SNYDER_HG19', uniq.indivs[k], uniq.marks[i], '0.RData', sep = '_'))
    for(j in 1:nfiles){
      tab = read.table(file.path(indir, rep.filenames[j]), header = F)
      if(j == 1){
        peak.ranges = GRanges(seqnames = Rle(tab[,1]), 
                              ranges = IRanges(start = tab[, 2] + 1, end = tab[, 3]),
                              strand = Rle(rep('+', dim(tab)[1])))
        bad = peak.ranges %in% bad.ranges
        regions = data.frame(chr = tab[!bad, 1], start = tab[!bad, 2] + 1, end = tab[!bad, 3])
        if(is.gene){rownames(regions) = tab[!bad, 4]
        }else{rownames(regions) = paste(regions$chr, regions$start, sep = '_')}
        counts = array(0, dim = c(sum(!bad), nfiles))
      }
      stopifnot(all(tab[!bad, 1] == regions$chr), all(tab[!bad, 2] + 1 == regions$start), all(tab[!bad, 3] == regions$end))
      counts[, j] = tab[!bad, ncol(tab)]
    }
    counts = data.frame(counts, row.names = rownames(regions))
    colnames(counts) = samples[sel]
    # If we want size factors, then we need to concatenate all datasets for the same mark.
    if(k == 1 && get.sf){counts.all = counts
    }else if(get.sf){counts.all = cbind(counts.all, counts)}
    save(regions, counts, file = outfile)    
  }
  
  # Compute size factors if necessary
  if(get.sf){
    all.sf = estimateSizeFactorsForMatrix(counts.all[!grepl('chr[XY]', regions$chr), ])
    for(k in 1:length(uniq.indivs)){
      outfile = file.path(outdir, paste('SNYDER_HG19', uniq.indivs[k], uniq.marks[i], '0.RData', sep = '_'))
      load(outfile)
      sel.cols = colnames(counts.all) %in% colnames(counts)
      size.factors = all.sf[sel.cols]
      stopifnot(all(colnames(counts.all)[sel.cols] == colnames(counts))) # The order of the columns should be the same
      save(regions, size.factors, counts, file = outfile)
    }
  }
}
