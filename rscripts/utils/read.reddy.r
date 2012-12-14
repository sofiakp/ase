rm(list=ls())
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/reddy/chip/liftover/')
infiles = list.files(indir, pattern = 'snpcov$', full.name = F, recursive = F, include.dirs = F)

nfiles = length(infiles)
hit.list = list()
ratios = list()
tfs = gsub('.uniq.filter.snpcov', '', infiles)
for(i in 1:nfiles){
  print(infiles[i])
  tab = read.table(file.path(indir, infiles[i]), header = F, sep = '\t')
  tf = tfs[i]
  pat = tab[, 3]
  mat = tab[, 4]
  pvals = binom.val(mat, mat + pat)
  qvals = p.adjust(pvals, method = 'BH')
  snp.info = data.frame(chr = tab[, 1], pos = tab[, 2], pval = pvals, qval = qvals, mat = mat, pat = pat)
  save(tf, snp.info, file = file.path(indir, gsub('.snpcov', '.RData', infiles[i])))
  pass = snp.info$qval < 0.05
  hit.list[[tfs[i]]] = paste(snp.info$chr[pass], snp.info$pos[pass], sep = ':')
  ratios[[tfs[i]]] = (mat[pass] + 1) / (pat[pass] + 1)
}
save(hit.list, ratios, file = file.path(indir, 'all_tfs.RData'))