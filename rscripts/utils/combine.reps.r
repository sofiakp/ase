rm(list=ls())
library('foreach')
library('doMC')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

combine.reps <- function(filenames, outfile, get.q = T, get.p = T){
  nfiles <- length(filenames)
  for(i in 1:nfiles){
    load(filenames[i])
    if(i == 1){
      all.counts <- counts
      all.sample <- sample
      all.snp.info <- snp.info
    }else{
      stopifnot(all.sample$indiv == sample$indiv)
      stopifnot(all.sample$mark == sample$mark)
      all.sample <- rbind(all.sample, sample)
      
      stopifnot(all(colnames(all.counts) == colnames(counts)))
      all.counts <- all.counts + counts
      
      stopifnot(all(all.snp.info$chr == snp.info$chr), all(all.snp.info$pos == snp.info$pos),
                all(all.snp.info$ref == snp.info$ref), all(all.snp.info$alt == snp.info$alt),
                all(all.snp.info$pass == snp.info$pass))
    }
  }
  counts <- all.counts
  snp.info <- all.snp.info
  sample <- all.sample
  if(length(filenames) > 1){
    mat.idx <- grep('mat', colnames(counts))
    pat.idx <- grep('pat', colnames(counts))
    fwd.idx <- grep('fwd', colnames(counts))
    sb <- binom.val(rowSums(counts[, fwd.idx]), rowSums(counts))
    pval <- array(1, dim = c(dim(snp.info)[1], 1))
    if(get.p) pval <- binom.val(rowSums(counts[, mat.idx]), rowSums(counts[, c(mat.idx, pat.idx)]))
    qval <- array(1, dim = c(dim(snp.info)[1], 1))
    if(get.q) qval[snp.info$pass] = p.adjust(pval[snp.info$pass], method = 'BH')
    snp.info$qval = qval
    snp.info$pval = pval
    snp.info$sb = sb
  }  
  save(counts, snp.info, sample, file = outfile)
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

registerDoMC(2)
filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/'), 
                        pattern = '[0-9]+\\.counts\\.r', full.name = T, recursive = T, include.dirs = F)
outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/', 'reps_6Aug12')
if(!file.exists(outdir)) dir.create(outdir)

all.info <- sample.info(filenames, '\\.counts\\.r$')
types <- unique(all.info[, 1:2])
overwrite = F
foreach(t = 1:dim(types)[1]) %dopar%{
  reps <- which(all.info$indiv == types$indiv[t] & all.info$mark == types$mark[t])
  if(length(reps) > 0){ # set to 0 to create rep files even for cases when there are no replicates
    print(basename(filenames[reps]))
    sel.reps = reps[select.reps(all.info$rep[reps])]
    # Write all the replicates (including technical)
    #if(length(sel.reps) != length(reps)){
    #  outfile <- file.path(outdir, gsub('[0-9\\.]+\\.counts.r$', 'reptr.counts.r', basename(filenames[reps[1]])))
    #  if(!file.exists(outfile) || overwrite) combine.reps(filenames[reps], outfile, get.q = F, get.p = F)
    #}
    outfile <- file.path(outdir, gsub('[0-9\\.]+\\.counts.r$', 'rep.counts.r', basename(filenames[reps[1]])))
    if(!file.exists(outfile) || overwrite) combine.reps(filenames[sel.reps], outfile, get.q = T)
  }
}