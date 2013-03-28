rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData') # Gene metadata
region.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.81bins_diff.bed'
count.dir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig/81bins_diff/textFiles/'
outdir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig/81bins_diff/rdata/'
if(!file.exists(outdir)) dir.create(outdir)
outpref = gsub('.bed', '', basename(region.file))
files = list.files(count.dir, pattern = outpref)

fsplit = strsplit(gsub('.*_AT_|.txt', '', files), '_')
marks = array(0, dim = c(length(files), 1))
indivs = array(0, dim = c(length(files), 1))
for(i in 1:length(files)){
  marks[i] = fsplit[[i]][4]
  indivs[i] = fsplit[[i]][3]
}
uniq.marks = unique(marks)
uniq.indivs = unique(indivs)

regions = read.table(region.file, header = F, sep = '\t')[, c(1,2,3,4,6,7)]
colnames(regions) = c('chr', 'start', 'end', 'gene.name', 'strand', 'bin')
good.regions = !duplicated(regions[, c(4,6)])

counts = NULL
cols = c()
for(m in uniq.marks){
  sel.files = grep(m, files)
  if(length(sel.files) < length(uniq.indivs)){
    cat('Skipping', m, '\n', file = stderr())
    indivs = indivs[!grepl(m, marks)]
    marks = marks[!grepl(m, marks)]
    next
  }
  for(i in uniq.indivs){
    cols = append(cols, paste(i, m, sep = '_'))
    sel.file = files[grep(paste(i, m, sep = '_'), files)]  
    print(sel.file)
    c = read.table(file.path(count.dir, sel.file), header = F, sep = '\t')
    stopifnot(nrow(c) == nrow(regions))
    c = c[good.regions, ]
    if(is.null(counts)){
      counts = c #Matrix(as.matrix(c))
    }else{
      counts = cbind(counts, c) #cBind(counts, Matrix(as.matrix(c))) 
    }
  }  
}
bins = regions[good.regions, ]
counts = normalize.quantiles(counts)
counts[counts < 1] = 0
counts[is.na(counts)] = 0
colnames(counts) = cols
counts = Matrix(counts)
save(bins, counts, indivs, marks, cols, file = file.path(outdir, paste(outpref, '.RData', sep = '')))

# for(i in 1:5){#length(files)){
#   print(files[i])
#   c = read.table(file.path(count.dir, files[i]), header = F, sep = '\t')
#   stopifnot(nrow(c) == nrow(regions))
#   tmp.counts = data.frame(gene = regions$gene.name, bin = regions$bin, c = c[, 1])[good.regions, ]
#   tmp.counts = cast(tmp.counts, gene~bin, function(x) mean(x, na.rm = T), value = 'c')
#   colnames(tmp.counts)[2:ncol(tmp.counts)] = paste(indivs[i], marks[i], colnames(tmp.counts)[2:ncol(tmp.counts)], sep = '_')
#   if(is.null(all.counts)){
#     sel.genes = tmp.counts$gene
#     tmp.counts = tmp.counts[, 2:ncol(tmp.counts)]
#     tmp.counts[tmp.counts < 0.5] = 0
#     all.counts = Matrix(data.matrix(tmp.counts))
#   }else{
#     stopifnot(all(sel.genes == tmp.counts$gene))
#     tmp.counts = tmp.counts[, 2:ncol(tmp.counts)]
#     tmp.counts[tmp.counts < 0.5] = 0
#     all.counts = cBind(all.counts, Matrix(data.matrix(tmp.counts)))
#   }
# }
