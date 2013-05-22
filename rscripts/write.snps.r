rm(list=ls())
library(Matrix)
source('utils/sample.info.r')

snp.pos.file = '../../rawdata/variants/sanConsensus/snps/san.snps.RData'
load(snp.pos.file)

indir = '../../rawdata/alleleCounts/san/rdata/reps/qvals/'
filenames = list.files(indir, pattern = 'SNYDER_HG19_.*\\_rep.*\\.RData', full.name = T, recursive = F, include.dirs = F)
outdir = file.path(indir, 'hitLists', 'text_withCounts')
if(!file.exists(outdir)) dir.create(outdir)
geno.dir = '../../rawdata/variants/sanConsensus/snps/'
overwrite = T
indivs = as.character(sample.info(filenames, '.RData')$indiv)
uniq.indiv = unique(indivs)

for(i in 1:length(uniq.indiv)){
  sel.files = filenames[grep(paste('SNYDER_HG19_', uniq.indiv[i], sep = ''), filenames)]
  nfiles = length(sel.files)
  load(file.path(geno.dir, paste(uniq.indiv[i], '.snps.RData', sep = '')))
  for(t in 1:nfiles){
    print(basename(sel.files[t]))
    outfile = file.path(outdir, gsub('.RData', '.hitInd.txt', basename(sel.files[t])))
    if(!file.exists(outfile) || overwrite ){
      load(sel.files[t])
      sel = as.vector(!is.na(snp.info$qval) & snp.info$qval < -2)
      out.dat = snp.pos[sel, ]
      out.dat$alt = geno.info$alt[sel]
      out.dat$mat = geno.info$mat[sel]
      out.dat$pat = geno.info$pat[sel]
      out.dat$unphased = geno.info$unphased[sel]
      out.dat = cbind(out.dat, as.matrix(counts[sel, ]))
      write.table(out.dat, file = outfile, sep = '\t', quote = F, row.names = F, col.names = T) 
      #write.table(snp.pos[as.vector(snp.info$hits), 1:2], file = outfile, sep = '\t', quote = F, row.names = F, col.names = F) 
    }
  } 
}
