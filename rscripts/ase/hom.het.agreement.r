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

snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)
nsnps = dim(snp.pos)[1]

geno.dir = '../../rawdata/variants/all/snps/allNonSan/' # Directory with genotype data, should have a file <indiv>.snps.RData for each individual.
count.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

filenames = list.files(count.dir, pattern = paste('.*rep\\.RData', sep = ''), full.name = T)
samples = sample.info(filenames)
marks = as.character(samples$mark)
indivs = as.character(samples$indiv)
uniq.indivs = unique(indivs)
uniq.marks = unique(marks)

# Read genotypes
mat = NULL # One column per individual with the maternal haplotype
pat = NULL # One column per individual with the paternal haplotype
for(i in 1:length(uniq.indivs)){
  load(file.path(geno.dir, paste(uniq.indivs[i], '.snps.RData', sep = '')))
  if(is.null(mat)){
    mat = geno.info$mat
    pat = geno.info$pat
  }else{
    mat = cBind(mat, geno.info$mat)
    pat = cBind(pat, geno.info$pat)
  }
}

rdat = NULL
rdat2 = NULL
eps = 1
for(i in 1:length(uniq.marks)){
  file.idx = which(marks == uniq.marks[i]) # Get all the files for the mark
  pass = NULL # Indicators of AS SNPs for the mark in any individual
  good = NULL
  het.counts = array(0, dim = c(nsnps, 2)) # Normalized number of reads for the ref and alt alleles in heterozygotes
  hom.counts = array(0, dim = c(nsnps, 2)) # Normalized number of reads for homozygotes ref and alt
  num.hom = array(0, dim = c(nsnps, 2)) # Number of homozygous individuals for ref and alt
  for(j in file.idx){
    print(basename(filenames[j]))
    load(filenames[j])
    indiv.idx = which(uniq.indivs == indivs[j])
    if(is.null(pass)){
      pass = as.vector(snp.info$qval < -2)
      good = as.vector(!snp.info$bad)
    }else{
      pass = as.vector(pass | snp.info$qval < -2)
      good = as.vector(good & !snp.info$bad)
    }
    sum.counts = counts[, 3] # 3rd column is total
    norm.factor = sum(sum.counts) / 1e6
    het.sites = mat[, indiv.idx] != pat[, indiv.idx]
    het.counts[het.sites, 1] = het.counts[het.sites, 1] + asinh(counts[het.sites, 1] / norm.factor)
    het.counts[het.sites, 2] = het.counts[het.sites, 2] + asinh(counts[het.sites, 2] / norm.factor)
    
    ref.sites = as.vector(!het.sites & !mat[, indiv.idx])
    hom.counts[ref.sites, 1] = hom.counts[ref.sites, 1] + asinh(sum.counts[ref.sites] / norm.factor)
    num.hom[ref.sites, 1] = num.hom[ref.sites, 1] + 1
    alt.sites = as.vector(!het.sites & mat[, indiv.idx])
    hom.counts[alt.sites, 2] = hom.counts[alt.sites, 2] + asinh(sum.counts[alt.sites] / norm.factor)
    num.hom[alt.sites, 2] = num.hom[alt.sites, 2] + 1
  }
  pass = pass & good & num.hom[, 1] > 1 & num.hom[, 2] > 1 & rowSums(num.hom) < 13 # Make sure the SNP is AS and there are homozygotes of both types
  if(sum(pass) > 50){
    het.ratio = log((het.counts[pass, 1] + eps)/ (het.counts[pass, 2] + eps), base = 2)
    hom.ratio = log((hom.counts[pass, 1] + eps) / (hom.counts[pass, 2] + eps) * (num.hom[pass, 2] / num.hom[pass, 1]), base = 2)
    corr = cor(hom.ratio, het.ratio)
    agree.frac = sum(hom.ratio * het.ratio > 0) / sum(pass)
    sel = sample(1:length(het.ratio), min(1000, length(het.ratio)))
    rtmp = data.frame(hom = hom.ratio[sel], het = het.ratio[sel], mark = paste(uniq.marks[i], sprintf('%.4f', corr)))
    rtmp2 = data.frame(mark = uniq.marks[i], f = agree.frac)
    rdat = rbind(rdat, rtmp)
    rdat2 = rbind(rdat2, rtmp2)
  }
}

p = ggplot(rdat) + geom_point(aes(x = hom, y = het, color = mark), size = 2) + theme_bw() + 
  scale_x_continuous('homozygotes log(hom_ref/hom_alt)') + scale_y_continuous('heterozygotes log(ref/alt)') +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14))
ggsave(file = file.path(plotdir, paste('hom_het_cor.png', sep = '')), p, width = 6.5, height = 5.6)
q = ggplot(rdat2) + geom_bar(aes(x = mark, y = f), stat = "identity") + scale_x_discrete('') + theme_bw() + 
  scale_y_continuous('Fraction of AS SNPs with the same direction') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave(file = file.path(plotdir, paste('hom_het_agreement.png', sep = '')), q, width = 6.5, height = 5.6)

# dip.regions = read.table('../../rawdata/signal/combrep/dips/llr/bed/SNYDER_HG19_H3K27AC_merged_dips.bed', header = F, sep = '\t')
# dip.regions = dip.regions[, 1:3]
# colnames(dip.regions) = c('chr', 'start', 'end')
# dip.ranges = regions.to.ranges(dip.regions)
# dip.ov = findOverlaps(snps.to.ranges(snp.pos), dip.ranges, select = 'first', ignore.strand = T)
# 
# clust.dat = read.table('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_cagt_clusters.bed')
# clust = data.frame(chr = clust.dat[, 1], pos = clust.dat[, 3], clust = clust.dat[, 9])
# 
# snp.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
# geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')
# 
# nhom = array(0, dim = c(nsnps, 2))
# files = list.files(snp.dir, pattern = 'H3K27AC_rep.RData$', full.names = T)
# samples = sample.info(files, '.RData')
# sel.files = !(as.character(samples$indiv) %in% c('GM12878', 'GM19240'))
# samples = samples[sel.files, ]
# files = files[sel.files]
# nfiles = length(files)
# is.hit = array(F, dim = c(nsnps, 1))
# pass = array(T, dim = c(nsnps, 1))
# eps = 1
# 
# for(i in 1:nfiles){
#   print(basename(files[i]))
#   load(files[i])
#   load(file.path(geno.dir, paste(as.character(samples$indiv[i]), '.snps.RData', sep = '')))
#   is.hit = is.hit | snp.info$qval < -2
#   #if(i == 2) is.hit = snp.info$qval < -2
#   pass = pass & !snp.info$bad
#   norm.fact = sum(counts[, 3]) / 1000000
#   #norm.fact.par = sum(counts[, 1:2]) / 1000000
#   hom.sites = as.vector(geno.info$mat == geno.info$pat)
#   hom.ref.sites = as.vector(hom.sites & !geno.info$mat)
#   hom.alt.sites = as.vector(hom.sites & geno.info$mat)
#   
#   if(i==1){c = array(0, dim = c(nsnps, 1))}else{c = as.array(hom.ref.counts)}
#   c[hom.ref.sites, ] = c[hom.ref.sites, ] + asinh(counts[hom.ref.sites, 3] / norm.fact)
#   hom.ref.counts = Matrix(c)
#   
#   if(i==1){c = array(0, dim = c(nsnps, 1))}else{c = as.array(hom.alt.counts)}
#   c[hom.alt.sites, ] = c[hom.alt.sites, ] + asinh(counts[hom.alt.sites, 3] / norm.fact)
#   hom.alt.counts = Matrix(c)
#   
#   if(i==1){c = array(0, dim = c(nsnps, 2))}else{c = as.array(het.counts)}
#   c[!hom.sites, 1] = c[!hom.sites, 1] + asinh(counts[!hom.sites, 1] / norm.fact)
#   c[!hom.sites, 2] = c[!hom.sites, 2] + asinh(counts[!hom.sites, 2] / norm.fact)
#   het.counts = Matrix(c)
#   
#   nhom[hom.ref.sites, 1] = nhom[hom.ref.sites, 1] + 1
#   nhom[hom.alt.sites, 2] = nhom[hom.alt.sites, 2] + 1
# }
# 
# het.ratios = (het.counts[, 1] + eps) / (het.counts[, 2] + eps)
# hom.ratios = ((hom.ref.counts + eps) / (hom.alt.counts + eps)) * (nhom[, 2] / nhom[, 1])
# have.var = as.vector(nhom[, 1] > 1 & nhom[, 2] > 1 & nfiles - rowSums(nhom) > 2 & is.hit & pass)
# r1 = hom.ratios[have.var]
# r2 = het.ratios[have.var]
# co = cor(r1, r2)
# cat('all', co, sum(log(r1) * log(r2) > 0) / sum(have.var), '\n')
# 
# nclust = max(clust$clust)
# biasp = array(1, dim = c(nclust, nclust))
# for(i in 1:nclust){
#   clust.ov = countOverlaps(snps.to.ranges(snp.pos), snps.to.ranges(clust[clust$clust == i, ]), ignore.strand = T)
#   sel = have.var & clust.ov > 2
#   r1 = hom.ratios[sel]
#   r2 = het.ratios[sel]
#   co = cor(r1, r2)
#   cat(i, co, sum(log(r1) * log(r2) > 0) / sum(sel), median(abs(log(r1))), median(abs(log(r2))), '\n')
#   if(i < nclust){
#     for(j in (i + 1):nclust){
#       selj = have.var & countOverlaps(snps.to.ranges(snp.pos), snps.to.ranges(clust[clust$clust == j, ]), ignore.strand = T) > 2
#       biasp[i, j] = wilcox.test(abs(log(r2)), abs(log(het.ratios[selj])))$p.value
#     } 
#   }
# }

# load('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/SNYDER_HG19_GM12878_H3K27AC_rep.hitInd.RData')
# snp.pos = snp.pos[as.vector(snp.info$hits), ]
# 
# dip.regions = read.table('../../rawdata/signal/combrep/dips/llr/bed/SNYDER_HG19_H3K27AC_merged_dips.bed', header = F, sep = '\t')
# dip.regions = dip.regions[, 1:3]
# colnames(dip.regions) = c('chr', 'start', 'end')
# dip.ranges = regions.to.ranges(dip.regions)
# ext.regions = dip.regions
# ext.regions$start = ext.regions$start - 200
# ext.regions$end = ext.regions$end + 200
# ext.dip.ranges = regions.to.ranges(ext.regions)
# 
# ov = findOverlaps(snps.to.ranges(snp.pos), ext.dip.ranges, select = 'first', ignore.strand = T)