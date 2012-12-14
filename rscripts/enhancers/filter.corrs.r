rm(list=ls())

# Find prox-distal pairs that are consistent for all prox peaks within the same TSS domain

element.corrs = function(prox.pos, prox, distal.pos, distal, start.ind, end.ind){
  co.mat = NULL
  for(i in start.ind:end.ind){
    sel = prox.pos$start[i] - distal.pos$end < max.dist & distal.pos$start - prox.pos$end[i] < max.dist
    if(any(sel)){
      co.vals = as.vector(cor(t(prox[i, ]), t(distal[sel, ])))
      co.tmp = data.frame(prox = rep(i, sum(sel)), distal = which(sel), co = co.vals)
      if(is.null(co.mat)){
        co.mat = co.tmp
      }else{co.mat = rbind(co.mat, co.tmp)}
    }
  }
  return(co.mat)
}

indir = 'rawdata/enhancers/external/degner/normalized/hg19'
outdir = 'rawdata/enhancers/external/degner/normalized/hg19/corr'
if(!file.exists(outdir)) dir.create(outdir)

chrom = 'chr16'
prox.file = file.path(indir, paste('merged.hg19.prox', chrom, 'txt', sep = '.'))

# Read proximal regions
prox = read.table(prox.file, header = F, skip = 1)
headers = as.matrix(read.table(prox.file, nrows = 1, as.is = T)[1,])
prox.pos = data.frame(prox[, 1:4])
colnames(prox.pos) = c('chr', 'start', 'end', 'name')
prox = prox[, 5:dim(prox)[2]]
indivs = headers[5:length(headers)]
nindivs = length(indivs)

annot = read.table('rawdata/enhancers/external/degner/normalized/dhs_wins.hg19.prox.bed', header = F)
annot = annot[annot[, 1] == chrom, ]
gene.map = list()
start = proc.time()
for(i in 1:dim(annot)[1]){
  gene = as.character(annot[i, 5])
  if(is.null(gene.map[[gene]])){
    gene.map[[gene]] = c()
  }
  gene.map[[gene]] = append(gene.map[[gene]], i)
}
print(proc.time() - start)

nprox = dim(prox)[1]
co = c()
for(i in 1:nprox){
  if(i %% 500 == 0) cat(i, '\n')
  gene = as.character(annot[i, 5])
  if(!is.null(gene.map[[gene]])){
    sel.list = gene.map[[gene]]
    sel.list = sel.list[sel.list > i]
    if(length(sel.list) > 0){
      co.vals = as.vector(cor(t(prox[i, ]), t(prox[sel, ])))
      co = append(co, co.vals)
    }
  }
}

#nchunks = 4
#chunk.size = ceiling(nprox / nchunks)
#registerDoMC(nchunks)
# co.mat = foreach(i = 1:nchunks, .combine = 'rbind', .inorder = T) %dopar% {
#   chunk.start = (i - 1) * chunk.size + 1
#   chunk.end = min(nprox, i * chunk.size)
#   element.corrs(prox.pos, prox, distal.pos, distal, chunk.start, chunk.end)
# }

#load('rawdata/enhancers/external/degner/normalized/hg19/corr/merged.hg19.corr.chr18.RData')
#sel.co = co.mat[co.mat$co > 0.5, ]