rm(list=ls())
library(Matrix)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# PCA on genotypes 

geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')
geno.files = list.files(geno.dir, pattern = '(GM|SN).*.snps.RData', full.names = T)
indivs = gsub('.snps.RData', '', basename(geno.files))
nindivs = length(indivs)

for(i in 1:nindivs){
  cat('Loading genotypes for', indivs[i], '\n')
  load(geno.files[i])
  if(i == 1){
    genot = geno.info$mat + geno.info$pat
  }else{
    genot = cBind(genot, geno.info$mat + geno.info$pat)
  }
}

# Select the SNPs that have enough variance in genotype
sel = rowSums(genot == 0) < nindivs - 2 & rowSums(genot == 1) < nindivs - 2 & rowSums(genot == 2) < nindivs - 2
genot = genot[sel, ]
genot.sample = genot[sample(1:nrow(genot), 100000), ]

# Remove daughters from PCA and then project them
genot.norm = scale(genot.sample)
pca.fit = prcomp(t(genot.norm[, !(indivs %in% c('GM12878', 'SNYDER'))]), center = F, scale = F)
p=plot.pcs(t(genot.norm) %*% pca.fit$rotation,  pca.fit$rotation, pca.fit$sdev, labels = indivs, groups = get.pop(indivs), all = T)
ggsave(file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/', 'genot_pca.png'), p$p1, width = 13.6, height = 11.8)
