rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
library(Matrix)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

# Computes overlaps between CAGT clusters of AS SNPs

geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/') # Where to read genotype information from
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

# Clustering of the signal of all individuals around the same set of AS SNPs 
clust.dat = read.table('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_cagt_clusters.bed')
clust = data.frame(chr = clust.dat[, 1], pos = clust.dat[, 3], clust = clust.dat[, 9])
nclust = max(clust$clust) # number of clusters
names = unlist(strsplit(as.character(clust.dat[, 4]), '_'))
clust$names = names[seq(1,length(names), 2)]
# Create a matrix SNP-by-individual, with the cluster index of the SNP signal for each individual.
clust = cast(clust, chr + pos ~ names, value = 'clust')
indivs = as.character(colnames(clust)[3:dim(clust)[2]])
nindivs = length(indivs)

# Select the SNPs from the genotype files that are in the cluster file
ov = findOverlaps(snps.to.ranges(snp.pos), snps.to.ranges(clust), select = 'first', ignore.strand = T)
sel.pos = as.vector(!is.na(ov))
snp.pos = snp.pos[sel.pos, ]
clust = clust[ov[sel.pos], ]
stopifnot(all(snp.pos$chr == clust$chr), all(snp.pos$pos == clust$pos)) # make sure they are in the same order

nsnps = dim(snp.pos)[1]
clust.hom.ref = array(0, dim = c(nsnps, nclust)) # [i,j]: Number of individuals that are homozygous ref for SNP i and their signal is in cluster j
clust.hom.alt = array(0, dim = c(nsnps, nclust))
clust.het = array(0, dim = c(nsnps, nclust))
genot = array(0, dim = c(nsnps, nindivs))

for(i in 1:nindivs){
  load(file.path(geno.dir, paste(indivs[i], '.snps.RData', sep = '')))
  hom.sites = as.vector(geno.info$mat[sel.pos] == geno.info$pat[sel.pos])
  hom.ref.sites = as.vector(hom.sites & !geno.info$mat[sel.pos])
  hom.alt.sites = as.vector(hom.sites & geno.info$mat[sel.pos])
  genot[, i] = as.integer(as.vector(geno.info$mat[sel.pos])) + as.integer(as.vector(geno.info$pat[sel.pos]))
  
  for(c in 1:nclust){
    in.clust = clust[, colnames(clust) == indivs[i]] == c
    clust.hom.ref[hom.ref.sites & in.clust, c] = clust.hom.ref[hom.ref.sites & in.clust, c] + 1
    clust.hom.alt[hom.alt.sites & in.clust, c] = clust.hom.alt[hom.alt.sites & in.clust, c] + 1
    clust.het[!hom.sites & in.clust, c] = clust.hom.alt[!hom.sites & in.clust, c] + 1
  }
}

have.cov = rowSums(clust.hom.ref) > 2 & rowSums(clust.hom.alt) > 2
are.excl = array(F, dim = c(nsnps, 1))

# Find SNPs such that the homozygous reference and homozygous alternative occupy mutually exlusive clusters
for(i in 1:nsnps){
  if(have.cov[i]){
    are.excl[i] = all(clust.hom.alt[i, clust.hom.ref[i, ] > 0] == 0) & all(clust.hom.ref[i, clust.hom.alt[i, ] > 0] == 0)    
  }
}