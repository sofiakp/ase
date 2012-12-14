rm(list=ls())
#library('fBasics')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

# DEPRECATED, USE compute.qvals.r

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals/')
filenames <- list.files(indir, pattern = 'SNYDER_HG19_GM128.*\\.counts\\.RData', full.name = T, recursive = F, include.dirs = F)
outdir = file.path(indir, 'hitLists')
if(!file.exists(outdir)) dir.create(outdir)
samples = sample.info(filenames, '\\.counts\\.RData$')    

q.cut = 0.01
nfiles = length(filenames)

for(t in 1:nfiles){
  print(basename(filenames[t]))
  hits.outfile = file.path(outdir, gsub('counts', 'hits', basename(filenames[t])))
  het.outfile = file.path(outdir, gsub('counts', 'het', basename(filenames[t])))
  if(!file.exists(hits.outfile) || !file.exists(het.outfile)){
    load(filenames[t])
    pass = !is.nan(snp.info$qval)
    idx1 = grep('ref', colnames(counts))
    idx2 = grep('alt', colnames(counts))
    ref.counts = rowSums(counts[, idx1])
    alt.counts = rowSums(counts[, idx2])
    tmp.counts = data.frame(chr = snp.info$chr, pos = snp.info$pos, ref.counts = ref.counts, alt.counts = alt.counts, 
                            mat = snp.info$mat, pat = snp.info$pat, phased = snp.info$phased)
    het.counts = tmp.counts[pass, ]
    signif = snp.info$qval[pass] < q.cut
    hits = het.counts[signif, ]
    het.counts$signif = signif
    print(sum(signif))
    save(het.counts, file = het.outfile)
    save(hits, file = hits.outfile)
  }
}
