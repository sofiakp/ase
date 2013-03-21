rm(list=ls())
library('foreach')
library('doMC')
library(DESeq)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

concat.rna.reps = function(filenames, samples, outdir, meta = NULL){
  # Combines replicates, summing up read groups. Size factors are computed using DESeq, but 
  # counts are kept unnormalized. Size factors are computed from all replicates of all individuals
  # provided, but the output is split per individual.
  
  stopifnot(length(unique(samples$mark)) == 1) # Doesn't make sense to compute size factors if this is true
  indivs = unique(as.character(samples$indiv))
  mark = as.character(samples$mark[1])
  
  en = new.env()
  nfiles = length(filenames)
  for(i in 1:nfiles){
    load(filenames[i], en)
    if(i == 1) agg.counts = array(0, c(dim(en$counts)[1], nfiles))
    agg.counts[, i] = rowSums(en$counts)
  } # COMPUTING SIZE FACTORS ON ALL DATA MIGHT TAKE FOREVER
  #print(head(agg.counts))
  #if(!is.null(meta)){
  #  pass = !grepl('chr[XY]', meta$chr) # Do not use sex chromosomes in the computation of size factors
  #}else{pass = array(T, dim = c(dim(agg.counts)[1], 1))}
  #all.size.factors = estimateSizeFactorsForMatrix(agg.counts[pass, ])
  
  for(i in 1:length(indivs)){
    outfile = file.path(outdir, paste('SNYDER_HG19', indivs[i], mark, '0.RData', sep = '_'))
    sel = samples$indiv == indivs[i]
    load(filenames[which(sel)[1]], en) # Reload one file for the individual to get the regions info.
    regions = data.frame(mask = en$regions$mask, het = en$regions$het, phased = en$regions$phased, row.names = rownames(en$regions))
    counts = data.frame(agg.counts[, sel], row.names = rownames(en$counts))
    colnames(counts) = paste(indivs[i], mark, as.character(samples$rep[sel]), sep = '_')
    #size.factors = all.size.factors[sel]
    save(counts, regions, file = outfile)
  }
}

combine.rna.reps = function(filenames, outfile, is.male = F, norm = F, meta = NULL){
  # Sums up replicates, keeping read groups separately and computes p-values for mat vs pat

  nfiles = length(filenames)
  if(norm){
    # Adjust replicates by their size factors before summing them up
    for(i in 1:nfiles){
      load(filenames[i])
      if(i == 1) agg.counts = array(0, c(dim(counts)[1], nfiles))
      agg.counts[, i] = rowSums(counts)
    }
    size.factors = estimateSizeFactorsForMatrix(agg.counts)    
  }
  else{size.factors = array(1, dim = c(nfiles, 1))}

  for(i in 1:nfiles){
    load(filenames[i])
    if(i == 1){
      all.counts = counts / size.factors[i]
      all.sample = sample
      all.info = regions
    }else{
      stopifnot(all.sample$indiv == sample$indiv)
      stopifnot(all.sample$mark == sample$mark)
      all.sample = rbind(all.sample, sample)
      
      stopifnot(all(colnames(all.counts) == colnames(counts)))
      stopifnot(all(rownames(all.counts) == rownames(counts)))
      all.counts = all.counts + counts / size.factors[i]
      
      stopifnot(all(colnames(all.info) == colnames(regions)))
      stopifnot(all(rownames(all.info) == rownames(regions)))
    }
  }
  if(is.male){
    if(!is.null(meta)){
      # Do not compute p-values on chrX for males. These genes might have very small p-values which will affect
      # the p-value correction.
      pass = !grepl('chr[XY]', meta$chr)
    }else{pass = !grepl('chr[XY]', all.info$chr)}
  }else{ pass = array(T, dim = c(dim(all.counts)[1], 1))} #!all.info$mask
  counts = round(all.counts) # Rounding will only have an effect if we computed size factors
  regions = all.info
  sample = all.sample
  mat.idx = grep('mat', colnames(counts))
  pat.idx = grep('pat', colnames(counts))
  pval = array(NaN, dim = c(length(pass), 1))
  pval[pass] = binom.val.par(counts[pass, mat.idx], rowSums(counts[pass, c(mat.idx, pat.idx)]))
  regions$pval = pval
  qval = array(NaN, dim = c(length(pval), 1))
  qval[pass] = p.adjust(regions$pval[pass], method = 'BH')
  print(sum(qval < 0.01, na.rm = T))
  regions$qval = qval
  save(counts, regions, sample, file = outfile)
}

# Reads files with counts in genes or exons. Replicates are assumed to be in different files, eg.
# SNYDER_HG19_GM12878_RNA_[123]_(exon|gene)counts.RData. Each file should contain a counts matrix
# and a regions data frame.
# 
# If rg is False then all replicate files are combined into one matrix with one column per replicate and written
# in SNYDER_HG19_indiv_mark_0.RData.
# Size factors for each replicate are computed with DESeq. The size factors 
# are based on ALL datasets for the same mark, eg. SNYDER_HG19_.*RNA.*
#
# If rg is True then replicate files are summed up and the counts matrix has 
# one column per read group with the sum of counts for the read group over all replicates. 
# p-values (and corrected p-values) for maternal vs paternal are 
# computed and appended to the regions data frame.

registerDoMC(2)
rg = F # T: sup up replicates, keep read groups separate, F: one file per indiv and mark with one column per replicate and RGs combined
is.gene = T # Read genecounts vs exoncounts
if(is.gene){
  subdir = 'geneCounts'
  suf = 'genecounts'
  load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData') # meta-data file
}else{
  subdir = 'exonCounts'
  suf = 'exoncounts'
  load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')
}
indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', subdir, 'rdata')
filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*(2255|2588|2610|2630|19193).*\\.', suf, '\\.RData', sep = ''), full.name = T, recursive = F, include.dirs = F)
if(rg){
  outdir = file.path(indir, 'reps/qvals')
}else{outdir = file.path(indir, 'repsComb')}
if(!file.exists(outdir)) dir.create(outdir, recursive = T)

all.info <- sample.info(filenames, paste('\\.', suf, '\\.RData$', sep = ''))
overwrite = F
males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486', 'GM2255', 'GM2588', 'GM2610', 'GM2630')

if(rg){
  types = unique(all.info[, 1:2])
  for(t in 1:dim(types)[1]){
    reps <- which(all.info$indiv == types$indiv[t] & all.info$mark == types$mark[t]) # Get all datasets for the same mark and individual
    if(length(reps) > 0){ # set to 0 to create rep files even for cases when there are no replicates
      print(basename(filenames[reps]))
      #sel.reps = reps[select.reps(all.info$rep[reps])]
      #print(basename(filenames[sel.reps]))
      # Write all the replicates (including technical)
      #if(length(sel.reps) != length(reps)){
      #  outfile <- file.path(outdir, gsub('[0-9\\.]+\\.counts.r$', 'reptr.counts.r', basename(filenames[reps[1]])))
      #  if(!file.exists(outfile) || overwrite) combine.reps(filenames[reps], outfile, get.q = F, get.p = F)
      #}
      outfile = file.path(outdir, gsub('INPUT', 'INPUT_rep', gsub(paste('[A-B0-9\\.]+\\.', suf, '.RData$', sep = ''), 
                                                                  paste('rep.RData', sep = ''), basename(filenames[reps[1]]))))
      if(!file.exists(outfile) || overwrite) combine.rna.reps(filenames[reps], outfile, is.male = types[t, 1] %in% males, norm = F, meta = gene.meta)
    }
  }
}else{
  marks = unique(as.character(all.info$mark)) 
  for(m in 1:length(marks)){
    reps <- which(all.info$mark == marks[m]) # Get all datasets for that mark
    print(basename(filenames[reps]))
    concat.rna.reps(filenames[reps], droplevels(all.info[reps, ]), outdir, meta = gene.meta)
  }
} 
