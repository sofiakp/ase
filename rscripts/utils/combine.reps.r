rm(list=ls())
library('foreach')
library('doMC')
library(Matrix)
library(DESeq)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

concat.reps = function(filenames, samples, outdir, meta = NULL){
  # Combines information about replicates. For each individual, one file is written for all replicates. 
  # The file has 3 matrices, tot, ref, and alt. Each matrix has one column per replicate.
  
  stopifnot(length(unique(samples$mark)) == 1) # Doesn't make sense to compute size factors if this is true
  indivs = unique(as.character(samples$indiv))
  mark = as.character(samples$mark[1])
  
  # Read all files to get total counts and compute size factors.
#   en = new.env()
#   nfiles = length(filenames)
#   for(i in 1:nfiles){
#     load(filenames[i], en)
#     if(i == 1){agg.counts = en$counts[, 3] #array(0, c(dim(en$counts)[1], nfiles))
#     }else{agg.counts = cbind(agg.counts, en$counts[, 3]) # 3rd column is the total counts
#   }
#   print(head(agg.counts))
#   if(!is.null(meta)){
#     pass = !grepl('chr[XY]', meta$chr) # Do not use sex chromosomes in the computation of size factors
#   }else{pass = array(T, dim = c(dim(agg.counts)[1], 1))}
#   all.size.factors = estimateSizeFactorsForMatrix(agg.counts[pass, ])

  # For each individual, read all files corresponding to that individual.
  for(i in 1:2){ #length(indivs)){
    print(indivs[i])
    outfile = file.path(outdir, paste('SNYDER_HG19', indivs[i], mark, '0.RData', sep = '_'))
    sel = which(samples$indiv == indivs[i])
    indiv.filenames = filenames[sel]
    for(j in 1:length(indiv.filenames)){
      load(filenames[j])
      if(j == 1){
        all.counts = counts
        bad = snp.info$bad
        print(head(all.counts))
      }else{
        all.counts = cBind(all.counts, counts)
        bad = bad | snp.info$bad
      } 
    }
    snp.info$bad = Matrix(bad)
    print(head(all.counts))
    ref = all.counts[, seq(1, dim(all.counts)[2], 3)]
    colnames(ref) = paste(indivs[i], mark, as.character(samples$rep[sel]), sep = '_')
    alt = all.counts[, seq(2, dim(all.counts)[2], 3)]
    colnames(alt) = paste(indivs[i], mark, as.character(samples$rep[sel]), sep = '_')
    tot = all.counts[, seq(3, dim(all.counts)[2], 3)]
    colnames(tot) = paste(indivs[i], mark, as.character(samples$rep[sel]), sep = '_')
    #size.factors = all.size.factors[sel]
    save(tot, ref, alt, snp.info, file = outfile)
  }
}

combine.reps <- function(filenames, samples, outfile, norm = F, meta = NULL){
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
      bad = snp.info$bad
    }else{
      all.counts = all.counts + counts / size.factors[i]
      bad = bad | snp.info$bad
    }
  }
  pass = as.vector(meta$mat != meta$pat & !bad)
  print(dim(all.counts))
  snp.info = list()
  snp.info$bad = Matrix(bad)
  counts = round(all.counts) # Rounding will only have an effect if we computed size factors
  pval = array(1, dim = c(length(pass), 1))
  pval[pass] = binom.val.par(counts[pass, 1], rowSums(counts[pass, 1:2]))
  snp.info$pval = Matrix(log(pval, base = 10))
  save(counts, snp.info, samples, file = outfile)
}

select.reps = function(reps){
  good.reps = array(T, dim = c(length(reps), 1))
  for(i in 1:5){
    ci = as.character(i)
    if(ci %in% reps){
      good.reps[grep(paste(ci, '\\.', sep = ''), reps)] = F
    }else if(paste(ci, '1', sep = '.') %in% reps){
      good.reps[grep(paste(ci, '\\.[2-9]', sep = ''), reps)] = F
    }
  }
  return(good.reps)
}

registerDoMC(3)
rg = T # T: sup up replicates, keep read groups separate, F: one file per indiv and mark

indir = '../../rawdata/alleleCounts/san/rdata/'
filenames = list.files(indir, pattern = 'SNYDER_HG19_.*RData$', full.name = T, recursive = F, include.dirs = F)
if(rg){
  outdir = file.path(indir, 'reps')
}else{outdir = file.path(indir, 'repsComb')}
if(!file.exists(outdir)) dir.create(outdir, recursive = T)

geno.dir = '../../rawdata/variants/sanConsensus/snps/' # genotype RData files should be here

all.info <- sample.info(filenames, '\\.RData$')
uniq.indivs = unique(as.character(all.info$indiv))
types <- unique(all.info[, 1:2])
overwrite = F

if(rg){
  for(i in uniq.indivs){
    uniq.marks = unique(as.character(all.info$mark[all.info$indiv == i]))
    # Load genotype data for the individual
    geno.file = file.path(geno.dir, paste(i, 'snps.RData', sep = '.'))
    load(geno.file)
    
    for(m in uniq.marks){
      reps = which(all.info$indiv ==  i & all.info$mark == m)
      outfile = file.path(outdir, paste('SNYDER_HG19', i, m, 'rep.RData', sep = '_'))
      if(!file.exists(outfile) || overwrite){
        print(basename(filenames[reps]))
        combine.reps(filenames[reps], droplevels(all.info[reps, ]), outfile, meta = geno.info)
      }else if(file.exists(outfile)){
        # Load the file and check which replicates it contains
        load(outfile)
        if(nrow(samples) != length(reps) ||
          !all(as.character(sort(samples$rep)) == as.character(sort(all.info[reps, 3])))){
          print(basename(filenames[reps]))
          combine.reps(filenames[reps], droplevels(all.info[reps, ]), outfile, meta = geno.info)
        }
      } 
    }
  }
}else{
  geno.file = file.path(geno.dir, 'allNonSan.snps.RData')
  load(geno.file)
  marks = unique(as.character(all.info$mark)) 
  for(m in 1:length(marks)){
    reps <- which(all.info$mark == marks[m]) # Get all datasets for that mark
    print(basename(filenames[reps]))
    concat.reps(filenames[reps], droplevels(all.info[reps, ]), outdir, meta = snp.pos)
  }
}
