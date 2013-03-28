rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(Matrix)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

# Reads CAGT results with the clustering of the signal of all individuals around the same set of AS SNPs 
# and computes the "mobility" of SNPs - a measure of how much they switch between different clusters.

##### Read SNP information
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/') # Where to read genotype information from
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

##### Read clustering information
load('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_cagt_results.RData')
nclust = max(clust$clust)
names = unlist(strsplit(as.character(clust$name), '_'))
clust$name = names[seq(1,length(names), 2)]
# Create a matrix SNP-by-individual, with the cluster index of the SNP signal for each individual.
clust = cast(clust, chr + pos ~ name, value = 'clust')
indivs = as.character(colnames(clust)[3:dim(clust)[2]])
nindivs = length(indivs)

##### Select the SNPs from the genotype files that participated in clustering
ov = findOverlaps(snps.to.ranges(snp.pos), snps.to.ranges(clust), select = 'first', ignore.strand = T)
sel.pos = as.vector(!is.na(ov))
snp.pos = snp.pos[sel.pos, ]
clust = clust[ov[sel.pos], ]
stopifnot(all(snp.pos$chr == clust$chr), all(snp.pos$pos == clust$pos)) # make sure they are in the same order

nsnps = dim(snp.pos)[1]
genot = array(0, dim = c(nsnps, nindivs))

##### Load genotype data
for(i in 1:nindivs){
  load(file.path(geno.dir, paste(indivs[i], '.snps.RData', sep = '')))
  genot[, i] = as.integer(as.vector(geno.info$mat[sel.pos])) + as.integer(as.vector(geno.info$pat[sel.pos]))
}

clust.mat = as.matrix(clust[, 3:ncol(clust)]) + 1 # Add one to all clusters. Now low signal is cluster 1
# The distance from the low signal to any othe cluster will be 1-0=1 (low signal cluster has 0 correlation with any other cluster)
tmp.centroid.dist = rbind(1, cbind(1, hc.centroid.dist)) 
diag(tmp.centroid.dist) = 0
mobility = array(0, dim = c(nsnps, 1)) # sum of mobilities for all pairs of individuals
mobility.same = array(0, dim = c(nsnps, 3)) # sum for pairs of individuals with the same genotype (all, hom only, het only)
mobility.diff = array(0, dim = c(nsnps, 1)) # sum for pairs of individuals with different genotypes
sum.diff = array(0, dim = c(nsnps, 1)) # sum of genotype differences
pairs = which(upper.tri(matrix(1, nindivs, nindivs)), arr.ind=T) # all pairs of indices of a row of genot
npairs = nrow(pairs)

for(i in 1:nsnps){
  # Get linear indices inside the matrix of cluster distances for each pair of individuals
  clust.ind = (clust.mat[i, pairs[, 1]] - 1) * nclust + clust.mat[i, pairs[, 2]]
  tot = sum(tmp.centroid.dist[clust.ind])
  # Get the pairs of individuals that have the same genotype
  sel = genot[i, pairs[, 1]] == genot[i, pairs[, 2]]
  same.dist = c(sum(tmp.centroid.dist[clust.ind[sel]]), 
                sum(tmp.centroid.dist[clust.ind[sel & genot[i, pairs[, 1]] != 1]]),
                sum(tmp.centroid.dist[clust.ind[sel & genot[i, pairs[, 1]] == 1]]))
  diff.dist = tot - same.dist[1]
  mobility[i] = tot
  mobility.same[i, ] = same.dist / c(sum(sel), sum(sel & genot[i, pairs[, 1]] != 1), sum(sel & genot[i, pairs[, 1]] == 1))
  mobility.diff[i] = diff.dist / (npairs - sum(sel))
  sum.diff[i] = sum(abs(genot[i, pairs[, 1]] - genot[i, pairs[, 2]]))
}

sel = !is.na(mobility.diff) & !is.na(mobility.same[, 1])
plot.dat = data.frame(same = mobility.same[sel, 1], diff = mobility.diff[sel])
p = ggplot(plot.dat) + geom_point(aes(x = same, y = diff), size = 1) + 
  xlab('Cluster dissimilarity for individuals with the same genotype') + ylab('Cluster dissimilarity for individuals with diff genotypes') + 
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/plots/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_mobile.png', 
       p, width = 6.5, height = 5.6)