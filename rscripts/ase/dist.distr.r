rm(list=ls())

countdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts', 'reps_6Aug12')
filenames <- list.files(countdir, pattern = 'GM12878.*rep\\.counts\\.r', full.name = T)
#samples = sample.info(filenames, '\\.counts\\.r$')    
nsamples = length(filenames)

qcut = 0.01
sb.cut = 0.01

pvals = array(NaN, dim = c(nsamples, 1))
for(i in 1:nsamples){
  load(filenames[i])
  uniq.chr = unique(snp.info$chr)
  dist = array(dim = c(0, 0))
  dist.rand = array(dim = c(0, 0))
  for(c in 1:length(uniq.chr)){
    cidx = snp.info$chr == uniq.chr[c]    
    good = snp.info$pass & cidx
    signif = good & snp.info$qval < qcut & snp.info$sb > sb.cut
    sam = sort(sample(1:sum(good), sum(signif))) # random set of SNPs (among those that qualify for qval computation)
    dist = append(dist, diff(snp.info$pos[signif])) # distance between consecutive significant SNPs
    dist.rand = append(dist.rand, diff(snp.info$pos[which(good)[sam]]))
  }
  print(median(dist)/median(dist.rand))
  if(length(dist) > 0) pvals[i] = wilcox.test(dist, dist.rand, paired = F)$p.value
}