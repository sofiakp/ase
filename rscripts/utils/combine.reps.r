rm(list=ls())
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

combine.reps <- function(filenames, outfile, get.q = T){
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
  mat.idx <- grep('mat', colnames(counts))
  pat.idx <- grep('pat', colnames(counts))
  fwd.idx <- grep('fwd', colnames(counts))
  sb <- binom.val(rowSums(counts[, fwd.idx]), rowSums(counts))
  pval <- binom.val(rowSums(counts[, mat.idx]), rowSums(counts[, c(mat.idx, pat.idx)]))
  qval <- array(1, dim = dim(pval))
  if(get.q) qval[snp.info$pass] = p.adjust(pval[snp.info$pass], method = 'BH')
  snp.info$qval = qval
  snp.info$pval = pval
  snp.info$sb = sb
  save(counts, snp.info, sample, file = outfile)
}

filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/'), 
                        pattern = '[0-9]+\\.counts\\.r', full.name = T)
all.info <- sample.info(filenames, '\\.counts\\.r$')
types <- unique(all.info[, 1:2])
for(t in 1:dim(types)[1]){
  reps <- which(all.info$indiv == types$indiv[t] & all.info$mark == types$mark[t])
  if(length(reps) > 1){
    outfile <- gsub('[0-9\\.]+\\.counts.r$', 'rep.counts.r', filenames[reps[1]])
    combine.reps(filenames[reps], outfile, get.q = T)
  }
}