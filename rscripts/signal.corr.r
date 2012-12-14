rm(list = ls())

filenames = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'geneCounts', 'reps/qvals/'),
                       pattern = paste('SNYDER_HG19_.*_RNA_rep.genecounts.RData', sep  = ''), full.name = T)

nfiles = length(filenames)
for(i in 1:nfiles){
  load(filenames[i])
  if(i == 1){
    all.counts = rowSums(counts)
  }else{
    all.counts = cbind(all.counts, rowSums(counts))
  }
}
good.rows = rowSums(all.counts) > 10
all.counts = all.counts[good.rows, ]
gene.info = gene.info[good.rows, ]
chr.counts = all.counts[gene.info$chr == 'chr1', ]