rm(list=ls())
library(foreach)
library(doMC)

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

f1 = read.table('rawdata/enhancers/external/degner/normalized/NA19238.qnorm.bed', header= F)
f2 = read.table('rawdata/enhancers/external/degner/normalized/NA19099.qnorm.bed', header= F)

# chrom = 'chr16'
# prox.file = file.path(indir, paste('merged.hg19.prox', chrom, 'txt', sep = '.'))
# dist.file = file.path(indir, paste('merged.hg19.distal', chrom, 'txt', sep = '.'))
# outfile = file.path(outdir, paste('merged.hg19.corr', chrom, 'RData', sep = '.'))
# max.dist = 500000 # maximum distance between a proximal-distal pair to compute correlations
# 
# # Read proximal regions
# prox = read.table(prox.file, header = F, skip = 1)
# headers = as.matrix(read.table(prox.file, nrows = 1, as.is = T)[1,])
# prox.pos = data.frame(prox[, 1:4])
# colnames(prox.pos) = c('chr', 'start', 'end', 'name')
# prox = prox[, 5:dim(prox)[2]]
# indivs = headers[5:length(headers)]
# nindivs = length(indivs)
# 
# # Read distal regions
# distal = read.table(dist.file, header = F, skip = 1)
# stopifnot(all(read.table(dist.file, nrows = 1) == headers)) # make sure the individuals are in the same order
# distal.pos = data.frame(distal[, 1:4])
# colnames(distal.pos) = c('chr', 'start', 'end', 'name')
# distal = distal[, 5:dim(distal)[2]]
# 
# nprox = dim(prox)[1]
# co.mat = NULL
# nchunks = 4
# chunk.size = ceiling(nprox / nchunks)
# registerDoMC(nchunks)
# co.mat = foreach(i = 1:nchunks, .combine = 'rbind', .inorder = T) %dopar% {
#   chunk.start = (i - 1) * chunk.size + 1
#   chunk.end = min(nprox, i * chunk.size)
#   element.corrs(prox.pos, prox, distal.pos, distal, chunk.start, chunk.end)
# }
# # foreach(i = 1:nprox, .combine = 'rbind', .inorder = T) %dopar% {
# #   if(i %% 500 == 0) cat(i / nprox, '\n')
# #   sel = prox.pos$start[i] - distal.pos$end < max.dist & distal.pos$start - prox.pos$end[i] < max.dist
# #   if(any(sel)){
# #     co.vals = as.vector(cor(t(prox[i, ]), t(distal[sel, ])))
# #     co.tmp = data.frame(prox = rep(i, sum(sel)), distal = which(sel), co = co.vals)
# #     if(is.null(co.mat)){
# #       co.mat = co.tmp
# #     }else{co.mat = rbind(co.mat, co.tmp)}
# #   }
# # }
# save(co.mat, prox.pos, distal.pos, prox, distal, file = outfile)