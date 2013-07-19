rm(list=ls())
library(ggplot2)
library(Matrix)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(foreach)
library(doMC)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# Tests how much variance is attributable to genotype vs population

get.cor.snp = function(cor, ov.mat, regions, snp.pos){
  # For each region, get the SNP with the highest correlation.
  cor.dat.tmp = data.frame(cor = cor, region.idx = ov.mat[, 1], chr = snp.pos$chr[ov.mat[, 2]], pos = snp.pos$pos[ov.mat[, 2]], pos.idx = ov.mat[, 2],
                           start = regions$start[ov.mat[, 1]], end = regions$end[ov.mat[, 1]])
  cor.dat.tmp.2 = cast(cor.dat.tmp, region.idx~., function(x) max(abs(x), na.rm = T) * sign(x[which.max(abs(x))]), value = 'cor') # get max cor per region
  colnames(cor.dat.tmp.2) = c('region.idx', 'cor')
  # create matrix region, cor, snp_chr, snp_pos, reg_start, reg_end. 
  # This might have more than one SNP per region, if there were multiple SNPs in the same region
  # with equally high correlation
  cor.dat.true = merge(cor.dat.tmp, cor.dat.tmp.2, by = c('region.idx', 'cor')) 
  cor.dat.true = cor.dat.true[!duplicated(cor.dat.true[,1:2]), ]
  return(cor.dat.true)
}

nchunks = 3
registerDoMC(nchunks)

plotdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots/qn_isvaNull_fits_all_reg_v2/'
rdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/'
mark = 'H3K4ME1'
pref = 'SNYDER_HG19_all_reg_'
comp = 3
qval = 0.01
K = 4
pref = paste(pref, mark, sep = '')
ext = 0
half.len = 1000 # Will adjust all regions so that their length is twice that

orig = new.env()
load(file.path(rdir, paste(pref, '_qn.RData', sep = '')), orig)
pref = paste(pref, '_comp', comp, '_q', qval, sep = '')
load(file.path(rdir, paste(pref, '_qn_isvaNull.RData', sep = '')))
regions.mid = round((regions$start + regions$end) / 2)
regions$start = pmax(regions.mid - half.len, 1)
regions$end = regions.mid + half.len
regions$start = pmax(regions$start - ext, 1)
regions$end = regions$end + ext
pop.regions = regions[isva.fit$deg, ]
no.pop.regions = regions[-isva.fit$deg, ]
counts = orig$counts[orig$good.rows, colnames(orig$counts) != 'GM12878' & colnames(orig$counts) != 'GM19240'][isva.fit$deg, ]
no.pop.counts = orig$counts[orig$good.rows, colnames(orig$counts) != 'GM12878' & colnames(orig$counts) != 'GM19240'][-isva.fit$deg, ]
indivs = colnames(counts)
nindivs = length(indivs)
sel.no.pop = sample(1:nrow(no.pop.regions), nrow(pop.regions))
no.pop.regions = no.pop.regions[sel.no.pop, ]
no.pop.counts = no.pop.counts[sel.no.pop, ]
regions = pop.regions
clust.pref = paste(pref, '_K', K, sep = '')
load(file.path(rdir, paste(clust.pref, '_clust.RData', sep = '')))

# Read SNP information
snp.pos.file = '../../rawdata/variants/all_Mar13/snps.RData'
load(snp.pos.file)

# Load genotypes
load('../../rawdata/variants/all_Mar13/genot.RData')
genot = genot[, match(indivs, fix.indiv.names(colnames(genot)))]
sel.genot = rowSums(genot == 0) < nindivs - 3 & rowSums(genot == 1) < nindivs - 3 & rowSums(genot == 2) < nindivs - 3 # Select positions where individuals differ
genot = genot[sel.genot, ]
snp.pos = snp.pos[sel.genot, ]

# Find all overlaps between regions and SNPs
ov = findOverlaps(regions.to.ranges(regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov))
ov.neg = findOverlaps(regions.to.ranges(no.pop.regions), snps.to.ranges(snp.pos), select = 'all', ignore.strand = T)
ov.mat.neg = cbind(queryHits(ov.neg), subjectHits(ov.neg))

# Correlation between signal at regions and genotype of all overlapping SNPs
pop.cor = cor.par(genot[ov.mat[, 2], ], counts[ov.mat[, 1], ], nchunks, 'spearman')
neg.cor = cor.par(genot[ov.mat.neg[, 2], ], no.pop.counts[ov.mat.neg[, 1], ], nchunks, 'spearman')

cor.dat.pop = get.cor.snp(pop.cor, ov.mat, regions, snp.pos)
cor.dat.neg = get.cor.snp(neg.cor, ov.mat.neg, no.pop.regions, snp.pos)
cor.dat = rbind(data.frame(cor = abs(cor.dat.pop[, 2]), type = rep('Haplotype-specific', nrow(cor.dat.pop))), 
                data.frame(cor = abs(cor.dat.neg[, 2]), type = rep('Other', nrow(cor.dat.neg))))

wilc.p = wilcox.test(abs(cor.dat.pop$cor), abs(cor.dat.neg$cor))$p.value
cat('Haplotype-specific vs other wilcoxon and medians', wilc.p, median(abs(cor.dat.pop$cor), na.rm = T), median(abs(cor.dat.neg$cor), na.rm = T), '\n')

pop = get.pop(indivs)
uniq.pop = c('San', 'YRI', 'Asian', 'CEU')
pop.genot = genot[cor.dat.pop$pos.idx, ] > 0
clust = kclusters$clust[cor.dat.pop$region.idx]
frac.div = array(0, dim = c(K, length(uniq.pop)))
for(i in 1:K){
  for(p in 1:length(uniq.pop)){
    sel.regions = clust == i
    sel.pop = pop == uniq.pop[p]
    g1 = rowSums(pop.genot[sel.regions, sel.pop] == 0) >= 0.75 * sum(sel.pop) & rowSums(pop.genot[sel.regions, !sel.pop] == 0) <= 0.25 * sum(!sel.pop)
    g2 = rowSums(pop.genot[sel.regions, sel.pop] == 1) >= 0.75 * sum(sel.pop) & rowSums(pop.genot[sel.regions, !sel.pop] == 1) <= 0.25 * sum(!sel.pop)
    frac.div[i, p] = sum(g1 | g2) / sum(sel.regions)
  }  
}
colnames(frac.div) = uniq.pop
rownames(frac.div) = paste('cluster', 1:K)
p = plot.tile(frac.div, x.ord.samples = colnames(frac.div), y.ord.sample = rownames(frac.div), midpoint = 0)
ggsave(file.path(plotdir, paste(clust.pref, '_genot_frac.pdf', sep = '')), p, width = 7, height = 6)

frac.div.n = frac.div / matrix(rep(colMaxs(frac.div), each = 4), nrow = K)
p = plot.tile(frac.div.n, x.ord.samples = colnames(frac.div), y.ord.sample = rownames(frac.div), midpoint = 0.8, ycex = 20, xcex = 20, lcex = 16)
ggsave(file.path(plotdir, paste(clust.pref, '_genot_frac_norm.pdf', sep = '')), p, width = 7, height = 6)

# sel.indivs = grep('CEU|YRI', pop) # Only consider these two populations
# pop = get.pop(indivs[sel.indivs])
# min.indiv = min(table(pop)) # Get the same number of individuals for each population
# sel.indivs = append(sample(sel.indivs[pop == 'CEU'], min.indiv), sample(sel.indivs[pop == 'YRI'], min.indiv))
# pop.genot = as.matrix(genot[cor.dat.pop$pos.idx, sel.indivs])
# pop.match.counts = counts[cor.dat.pop$region.idx, sel.indivs]
# no.pop.genot = as.matrix(genot[cor.dat.neg$pos.idx, sel.indivs])
# no.pop.match.counts = no.pop.counts[cor.dat.neg$region.idx, sel.indivs]
# pop = factor(get.pop(indivs[sel.indivs]))
# pop.coefs = array(NaN, dim = c(nrow(pop.genot), 2))
# no.pop.coefs = array(NaN, dim = c(nrow(no.pop.genot), 2))
# #sum(rowSums(pop.genot[, pop == 'CEU'] > 0) > 0 & rowSums(pop.genot[, pop == 'YRI'] > 0) > 0 & rowSums(pop.genot[, pop == 'CEU'] == 0) > 0 & rowSums(pop.genot[, pop == 'YRI'] == 0) > 0)
# 
# # Compute correlations signal-genot and signal-pop
# for(i in 1:nrow(no.pop.coefs)){
#   # Only take regions with 2 different alleles present in equal proportions
#   if(sum(no.pop.genot[i, ] == 0) < 5 && sum(no.pop.genot[i, ] == 1) < 5 && sum(no.pop.genot[i, ] == 2) < 5 && length(unique(no.pop.genot[i, ])) == 2){
#     no.pop.coefs[i, 1] = cor(no.pop.match.counts[i, ], no.pop.genot[i, ], method = 'spearman')
#     no.pop.coefs[i, 2] = cor(no.pop.match.counts[i, ], as.numeric(pop), method = 'spearman')
#     #go = glm(pop.match.counts[i, sel.pop] ~ (pop.genot[i, sel.pop] > 0) + pop[sel.pop])
#     #coefs[i, 1] = coefficients(go)[2]
#     #coefs[i, 2] = coefficients(go)[3]
#   }
# }
# 
# for(i in 1:nrow(pop.coefs)){
#   if(sum(pop.genot[i, ] == 0) < 5 && sum(pop.genot[i, ] == 1) < 5 && sum(pop.genot[i, ] == 2) < 5 && length(unique(pop.genot[i, ])) == 2){
#     pop.coefs[i, 1] = cor(pop.match.counts[i, ], pop.genot[i, ], method = 'spearman')
#     pop.coefs[i, 2] = cor(pop.match.counts[i, ], as.numeric(pop), method = 'spearman')
#     #go = glm(pop.match.counts[i, sel.pop] ~ (pop.genot[i, sel.pop] > 0) + pop[sel.pop])
#     #coefs[i, 1] = coefficients(go)[2]
#     #coefs[i, 2] = coefficients(go)[3]
#   }
# }
# sum(abs(no.pop.coefs[,1]) < abs(no.pop.coefs[,2]), na.rm=T) / sum(abs(no.pop.coefs[,1]) != abs(no.pop.coefs[,2]), na.rm=T)
# sum(abs(pop.coefs[,1]) < abs(pop.coefs[,2]), na.rm=T) / sum(abs(pop.coefs[,1]) != abs(pop.coefs[,2]), na.rm=T)

# for(i in 1:100){
#   if(sum(pop.genot[i, ] == 1) > 3){
#     sel.cols = pop.genot[i, ] == 1
#     coefs[i, 1] = sd(pop.match.counts[i, sel.cols]) / mean(pop.match.counts[i, sel.cols])  
#     sel.cols = pop == 'San'
#     coefs[i, 2] = sd(pop.match.counts[i, sel.cols]) / mean(pop.match.counts[i, sel.cols])  
#   }  
# }
# sel.pop = grepl('CEU|San', as.character(pop))
# coefs = array(NaN, dim = c(nrow(pop.genot), 2))
# for(i in 1:100){
#   if(length(unique(pop.genot[i, sel.pop])) > 1){
#     go = glm(pop.match.counts[i, sel.pop] ~ (pop.genot[i, sel.pop] > 0) + pop[sel.pop])
#     coefs[i, 1] = coefficients(go)[2]
#     coefs[i, 2] = coefficients(go)[3]
#   }
# }