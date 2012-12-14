rm(list=ls())
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

infile = 'rawdata/alleleCounts/reps_6Aug12/SNYDER_HG19_GM12878_H3K27AC_rep.counts.r'

load(infile)
het.cov = rowSums(counts[snp.info$pass, ])
sign.cov = rowSums(counts[snp.info$qval < 0.01, ])