rm(list=ls())

indir = file.path('../../rawdata/geneCounts/rdata/repsComb/')
filenames <- list.files(indir, pattern = 'RZ.*\\.RData', full.name = T, recursive = F, include.dirs = F)
type = 'rpkm'
outdir = file.path(indir, type)
if(!file.exists(outdir)) dir.create(outdir)
nfiles = length(filenames)

load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
overwrite = T

for(t in 1:nfiles){
  print(basename(filenames[t]))
  outfile = file.path(outdir, gsub('.RData', '.txt', basename(filenames[t])))
  if(overwrite || !file.exists(outfile)){
    load(filenames[t])
    if(type == 'rpkm'){
      if(dim(counts)[2] > 1){
        #counts = t(apply(counts, 1, function(x) x / size.factors))
        rpkm = rowSums(counts)
        #reg = colSums(counts)
        #counts = t(apply(counts, 1, function(x) x * 10^6 / reg))
      }else{
        rpkm = counts
        #counts = counts / size.factors
        #counts = counts * 10^6 / sum(counts)
      }
      widths = gene.meta$len
      rpkm = (rpkm * 10^6 / widths) / sum(rpkm)
      cat(sum(rpkm > 0.5), '\n')
      #counts = apply(counts, 2, function(x) x / widths)
      write.table(cbind(gene.meta, rpkm), file = outfile, sep = '\t', quote = F, row.names = F, col.names = T)
    }else{
      pass = !is.nan(gene.info$qval) & gene.info$qval < 0.01
      #write.table(unique(gsub('.[0-9]+$', '', gene.info$gene[pass])), file = outfile, sep = '\t', quote = F, row.names = F, col.names = F) 
    }
  }
}