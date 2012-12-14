rm(list=ls())
library('GenomicRanges')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

mark = 'CTCF'
filenames = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'alleleCounts', 'reps_15Sep12/qvals/'),
                       pattern = paste('SNYDER_HG19_.*_', mark, '_rep.counts.RData', sep  = ''), full.name = T)
outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'alleleCounts', 'reps_15Sep12/merged')
if(!file.exists(outdir)) dir.create(outdir)
outfile = file.path(outdir, paste('SNYDER_HG19', mark, 'counts.RData', sep = '_'))

samples = sample.info(filenames, '.counts.RData')

nfiles = length(filenames)
for(i in 1:nfiles){
  load(filenames[i])
  print(basename(filenames[i]))
  new.ranges = GRanges(seqnames = Rle(snp.info$chr), ranges = IRanges(start = snp.info$pos, end = snp.info$pos), strand = Rle(rep('+', dim(snp.info)[1])))
  if(i == 1){
    joint.snp.ranges = new.ranges
    snp.names = paste(snp.info$chr, snp.info$pos, collapse = ':')
  }else{
    hits1 = joint.snp.ranges %in% new.ranges
    joint.snp.ranges = joint.snp.ranges[hits1, ]
  }
}

joint.counts = array(0, dim = c(length(seqnames(joint.snp.ranges)), nfiles))
for(i in 1:nfiles){
  load(filenames[i])
  new.ranges = GRanges(seqnames = Rle(snp.info$chr), ranges = IRanges(start = snp.info$pos, end = snp.info$pos), strand = Rle(rep('+', dim(snp.info)[1])))
  ov = findOverlaps(joint.snp.ranges, new.ranges)
  joint.counts[, i] = rowSums(counts[subjectHits(ov), ])
}
counts = joint.counts
snp.info = data.frame(chr = as.factor(seqnames(joint.snp.ranges)), pos = start(joint.snp.ranges))
save(counts, snp.info, samples, file = outfile)