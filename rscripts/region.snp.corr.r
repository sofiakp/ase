rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(foreach)
library(doMC)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source('/home/sofiakp/Documents/rpackages/myComp/isva.r')

# Computes correlation between signal at regions (counts, input normalized signal etc) and genotypes (0/1/2 alternative alleles)
# at overlapping SNPs
# Regions must be non-overlapping (so each SNP will overlap at most 1 region). 

wilc.par = function(x, y, nchunks = 10, mc = 3){
  stopifnot(all(dim(x) == dim(y)))
  nvals = dim(x)[1]
  chunk.size = ceiling(nvals / nchunks)
  pval = foreach(i = 1:nchunks, .combine = 'append', .inorder = T) %dopar%{
    chunk.start = (i - 1) * chunk.size + 1
    chunk.end = min(i * chunk.size, nvals)
    chunk.true.size = chunk.end - chunk.start + 1
    ctmp = array(0, dim = c(chunk.true.size, 1))
    for(j in 0:(chunk.true.size - 1)){
      g = x[chunk.start + j, ]
      if(sum(g > 0) > mc && sum(g == 0) > mc){
        p1 = -log(wilcox.test(as.numeric(y[chunk.start + j, g > 0]), as.numeric(y[chunk.start + j, g == 0]), alternative = 'two.sided', paired = F)$p.value, base = 10)
      }else{p1 = 0}
      if(sum(g == 1) > mc && sum(g == 2) > mc){
        p2 = -log(wilcox.test(as.numeric(y[chunk.start + j, g == 1]), as.numeric(y[chunk.start + j, g == 2]), alternative = 'two.sided', paired = F)$p.value, base = 10)
      }else{p2 = 0}
      ctmp[j + 1] = max(p1, p2)
    }
    ctmp
  }
  return(pval)
}

compute.hits = function(counts, genot, thresh.list, nperm = 20, method = 'kendall'){
  # cor.true = array(1, dim = c(dim(genot)[1], 1))
  # for(i in 1:dim(genot)[1]){
  #   dat = data.frame(y = as.numeric(counts[i, ]), x = genot[i, ])
  #   alt = glm('y ~ x', family = 'gaussian', data = dat)
  #   cor.true[i] = summary(alt)$coefficients[2, 4]
  # }
  cor.true = cor.par(genot, counts, nchunks, 'kendall') #wilc.par(genot, counts, nchunks)
  cor.perm = c()
  for(p in 1:nperm){
    perm = sample(1:nindivs)
    cor.perm.tmp = cor.par(genot, counts[, perm], nchunks, 'kendall') #wilc.par(genot, counts[, perm], nchunks)
    cor.perm = append(cor.perm, cor.perm.tmp)
  }
  
  fdrs = array(0, dim = c(length(thresh.list), 1))
  hits = array(0, dim = c(length(thresh.list), 1))
  for(i in 1:length(thresh.list)){
    hits[i] = sum(abs(cor.true) > thresh.list[i])
    fdrs[i] = sum(abs(cor.perm) > thresh.list[i]) / (nperm * hits[i])
  }
  return(list(cor.true = cor.true, sel.thresh, hits = hits))
}

nchunks = 5
registerDoMC(nchunks)

peak.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/') #'rawdata/signal/combrep/countsAtPeaksBroad/repsComb/')  #rawdata/signal/combrep/dips/llr/avgSig')
# bed.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles')
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')
plotdir = file.path(peak.dir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = 'allNonSan_vs_GM12878_allTFBS_noPol_H3K27AC_enrich'

# List of TF binding sites or DNase sites, or other locations. Only SNPs in these regions will be considered.
# Set to NULL if you want all regions.
tf.regions = read.bed(file.path(Sys.getenv('MAYAROOT'), 'rawdata/Gm12878_allTFBS.sorted.noPol.bed'))

get.all.hits = T # Compute the number of hits after each preprocessing step?
hits = c()

# Load SNPs
snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)

# p-value thresholds to try in order to find the threshold that gives FDR ~ 10% 
thresh.list = seq(0, 1, 0.01)
nperm = 20 # number of permutations for p-value computation

mark = 'H3K27AC'
region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/', 
                        paste('SNYDER_HG19', mark, 'merged.encodePeak', sep = '_'))
signal.files = list.files(peak.dir, pattern = paste('SNYDER_HG19', mark, 'merged_AT_.*txt', sep = '_'), full.names=T)
indivs = unique(gsub(paste('SNYDER_HG19_', mark, '_merged_AT_SNYDER_HG19_|_', mark, '.*.txt', sep = ''), '', basename(signal.files)))
sel.files = !(indivs %in% c('GM12878', 'GM19240')) # DO NOT consider daughters for computing correlations
indivs = indivs[sel.files]
nindivs = length(indivs)
counts.dat = load.avg.sig.data(region.file, signal.files[sel.files], indivs)
regions = counts.dat$regions
counts = asinh(counts.dat$signal)
good.rows = apply(counts, 1, function(x) !any(is.na(x)))
counts = counts[good.rows, ]
regions = regions[good.rows, ]
row.means = rowMeans(counts)
row.sds = rowSds(counts)
cvs = row.sds / row.means
good.rows = !is.na(cvs) & cvs > quantile(cvs, 0.4)
counts = counts[good.rows, ]
regions = regions[good.rows, ]

# Overlap the SNPs and the regions
# Each SNP should overlap at most one region, assuming the regions are non-overlapping.
ranges = regions.to.ranges(regions)
ov = findOverlaps(snps.to.ranges(snp.pos), ranges, ignore.strand = T, select = 'all')
ov = cbind(queryHits(ov), subjectHits(ov)) # Convert to matrix for fast access. Col 1 - SNP idx, Col 2 - region idx

# Now remove regions with too many overlapping SNPs. 
nhits = table(ov[, 2]) # For each region, count how many SNPs it overlaps
good.regions = as.numeric(names(nhits[nhits < 100])) # This will return the indices of regions (colnames of nhits) that don't have too many SNPs
ov = ov[ov[, 2] %in% good.regions, ]
has.ov = 1:dim(snp.pos)[1] %in% ov[, 1] # ov[, 1] is sorted so indexing by has.ov or ov[, 1] is the same thing

# Load the genotypes for the SNPs that do have overlap
genot = array(0, dim = c(sum(has.ov), nindivs))
snp.env = new.env()
for(i in 1:nindivs){
  cat('Loading genotypes for', indivs[i], '\n')
  load(file.path(geno.dir, paste(indivs[i], 'snps.RData', sep = '.')), snp.env)
  genot[, i] = as.integer(as.vector(snp.env$geno.info$mat[has.ov])) + as.integer(as.vector(snp.env$geno.info$pat[has.ov]))
}

# Select the SNPs that have enough variance in genotype
sel = rowSums(genot == 0) < nindivs - 3 & rowSums(genot == 1) < nindivs - 3 & rowSums(genot == 2) < nindivs - 3
sel = sel & !is.na(findOverlaps(snps.to.ranges(snp.pos[has.ov, ]), regions.to.ranges(tf.regions), select = 'first', ignore.strand = T))
genot = genot[sel, ]
regions = regions[ov[sel, 2], ]
# has.ov will be an indicator vector of which of the initial SNPs are in the genotype matrix
# so I turn to false all the positions that had overlaps but were removed in the step above
has.ov[which(has.ov)[!sel]] = F

if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh', sum(sel.f)))
}
counts = normalize.quantiles(counts)
if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh-qn', sum(sel.f))))
}
counts = scale(counts, center = T, scale = F)
if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh-qn-mean', sum(sel.f))))
}
isva.fit = isva(counts, maxq = 0.1)
for(i = 1:isva.fit$ncomp){
  mod = model.matrix( ~isva.fit$isv[, i])
  df = dim(mod)[2]
  Id = diag(ncol(counts))
  resid = counts %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  normalize.quantiles(resid)
  hit.dat = compute.hits(resid[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], 
                               method = rep(paste('asinh-qn-mean-isv', i, '-sign', sep = ''), sum(sel.f))))
}
isva.fit = isva(counts, maxq = 1)
for(i = 1:isva.fit$ncomp){
  mod = model.matrix( ~isva.fit$isv[, i])
  df = dim(mod)[2]
  Id = diag(ncol(counts))
  resid = counts %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  normalize.quantiles(resid)
  hit.dat = compute.hits(resid[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], 
                               method = rep(paste('asinh-qn-mean-isv', i, sep = ''), sum(sel.f))))
}
#colnames(genot) = indivs
#colnames(counts) = indivs

# idx = which.min(abs(fdrs - 0.1))
# sel.thresh = thresh.list[idx]
# fdr = fdrs[idx]
# hits = sum(abs(cor.true) > sel.thresh)
# 
# plot.sample = min(length(cor.true), 20000) 
# plot.dat = data.frame(cor = append(sample(cor.true, plot.sample), sample(cor.perm, plot.sample)),
#                       type = factor(append(rep('true', plot.sample), rep('permuted', plot.sample))))
# p = ggplot(plot.dat) + geom_density(aes(x = cor, color = type)) + xlab('Density') + ylab('Correlation') +
#   ggtitle('Correlation between levels of H3K27ac and genotypes at overlapping SNPs') +
#   theme(legend.title = element_blank(), legend.text = element_text(size = 12),
#         axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
#ggsave(file.path(plotdir, paste(outpref, '_pvals.png', sep = '')), p, width = 13.6, height = 11.8)

# save(counts, genot, has.ov, snp.pos.file, cor.true, sel.thresh, file = file.path(plotdir, paste(outpref, '.RData', sep = '')))

# Load signal in regions (dips, peak regions etc)
# dip.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/dips/llr/bed/SNYDER_HG19_H3K27AC_merged_dips.bed')
# dip.count.files = list.files(peak.dir, pattern = 'SNYDER_HG19_H3K27AC_merged_dips_AT_SNYDER_HG19_.*txt', full.names = T)
# indivs = unique(as.character(sample.info(gsub('SNYDER_HG19_H3K27AC_merged_dips_AT_', '', dip.count.files), '.txt')$indiv))
# sel.files = !(indivs %in% c('GM12878', 'GM19240')) # DO NOT consider daughters for computing correlations
# indivs = indivs[sel.files] 
# dip.count.files = dip.count.files[sel.files]
# nindivs = length(indivs)
# counts.dat = load.dip.data(dip.file, dip.count.files, indivs)
# counts.dat$dip.sig[counts.dat$dip.sig < 0] = 0
# counts.dat$left.sig[counts.dat$left.sig < 0] = 0
# counts.dat$right.sig[counts.dat$right.sig < 0] = 0
# 
# # Select dips that have evidence for being dips in at least a few individuals. This will remove some false positives.
# is.dip = counts.dat$dip.sig < counts.dat$left.sig & counts.dat$dip.sig < counts.dat$right.sig # & !is.na(findOverlaps(regions.to.ranges(counts.dat$regions), regions.to.ranges(tf.regions), select = 'first', ignore.strand = T))
# sel = rowSums(is.dip, na.rm = T) > 2
# regions = counts.dat$regions[sel, ]
# counts = (exp(counts.dat$left.sig[sel, ]) + exp(counts.dat$right.sig[sel, ])) / 2
# 
# sel = apply(counts, 1, function(x) !any(is.na(x)))
# regions = regions[sel, ]
# counts = counts[sel, ]

# indivs = unique(as.character(sample.info(list.files(peak.dir, pattern = 'SNYDER_HG19_.*H3K27AC_0.RData'), '.RData')$indiv))
# sel.files = !(indivs %in% c('GM12878', 'GM19240')) # DO NOT consider daughters for computing correlations
# indivs = indivs[sel.files]
# nindivs = length(indivs)
# bed.files = array('', dim = c(nindivs, 1))
# for(i in 1:nindivs){
#   tmp.files = list.files(bed.dir, pattern = paste('SNYDER_HG19_', indivs[i], '_H3K27AC_VS_.*', sep = ''), full.names = T)
#   stopifnot(length(tmp.files) == 1)
#   bed.files[i] = tmp.files
# }
# # ac.diff = get.diff.count(list.files(peak.deseq.dir, pattern = '*H3K27AC_deseq.RData', full.names = T), 0.01) > 5
# counts.dat = avg.counts(peak.dir, indivs, 'H3K27AC_0.RData', meta = NULL)
# regions = counts.dat$regions
# counts = asinh(counts.dat$counts)
# counts = fit.norm.counts(counts, bed.files, fit = 'rlm', ref.idx = which(indivs == 'GM12891'))
