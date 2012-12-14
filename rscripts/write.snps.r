rm(list=ls())

snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/')
filenames = list.files(indir, pattern = 'SNYDER_HG19_.*\\.hitInd\\.RData', full.name = T, recursive = F, include.dirs = F)
outdir = file.path(indir, 'text')
if(!file.exists(outdir)) dir.create(outdir)
nfiles = length(filenames)
overwrite = T

for(t in 1:nfiles){
  print(basename(filenames[t]))
  outfile = file.path(outdir, gsub('.RData', '.txt', basename(filenames[t])))
  if(!file.exists(outfile) || overwrite ){
    load(filenames[t])
    write.table(snp.pos[as.vector(snp.info$hits), 1:2], file = outfile, sep = '\t', quote = F, row.names = F, col.names = F) 
  }
}
