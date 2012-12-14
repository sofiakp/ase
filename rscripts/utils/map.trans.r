rm(list=ls())

load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
tab = read.table('rawdata/transcriptomes/gencode.v13.annotation.noM.trans.txt', header = F, sep = '\t')
ntrans = dim(tab)[1]
trans = data.frame(chr = tab[,1], start = tab[,2] + 1)
neg = tab[, 4] == '-'
trans$start[neg] = tab[neg, 3]
names = as.character(gene.meta$gene.name)

gene.idx = array(0, dim = c(ntrans, 1))
for(i in 1:ntrans){
  idx = which(tab[i,5] == names)
  if(length(idx) == 1) gene.idx[i] = idx
}
trans$gene.idx = gene.idx
trans = trans[gene.idx != 0, ]
save(trans, gene.meta, file = 'rawdata/transcriptomes/gencode.v13.annotation.noM.trans.RData')