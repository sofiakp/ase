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

# Evaluates different methods for normalizing/correcting the signal based on the number of significant associations between
# genotype and signal they lead to.

compute.hits = function(counts, genot, thresh.list, nperm = 20, method = 'kendall'){
  cor.true = cor.par(genot, counts, nchunks, method)
  cor.perm = c()
  for(p in 1:nperm){
    perm = sample(1:nindivs)
    cor.perm.tmp = cor.par(genot, counts[, perm], nchunks, method)
    cor.perm = append(cor.perm, cor.perm.tmp)
  }
  
  fdrs = array(0, dim = c(length(thresh.list), 1))
  hits = array(0, dim = c(length(thresh.list), 1))
  for(i in 1:length(thresh.list)){
    hits[i] = sum(abs(cor.true) > thresh.list[i])
    fdrs[i] = sum(abs(cor.perm) > thresh.list[i]) / (nperm * hits[i])
  }
  return(list(cor.true = cor.true, fdrs = fdrs, hits = hits))
}

nchunks = 1
registerDoMC(nchunks)

peak.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/')
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')
plotdir = file.path(peak.dir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = 'allNonSan_vs_GM12878_allTFBS_noPol_H3K27AC_enrich_eval3'

# Will only consider SNPs that overlap TF binding sites
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
counts = asinh(counts.dat$signal) # asinh normalize
good.rows = apply(counts, 1, function(x) !any(is.na(x)))
counts = counts[good.rows, ]
regions = regions[good.rows, ]
row.means = rowMeans(counts)
row.sds = rowSds(counts)
cvs = row.sds / row.means
good.rows = !is.na(cvs) & cvs > quantile(cvs, 0.5)
counts = counts[good.rows, ]
regions = regions[good.rows, ]

# Overlap the SNPs and the regions
# Each SNP should overlap at most one region, assuming the regions are non-overlapping.
ov = findOverlaps(snps.to.ranges(snp.pos), regions.to.ranges(regions), ignore.strand = T, select = 'all')
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

if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh', sum(sel.f)))
}
counts = normalize.quantiles(counts) ######### Quantile normalization
if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh-qn', sum(sel.f))))
}

counts = scale(counts, center = T, scale = T) ######### Mean centering
if(get.all.hits){
  hit.dat = compute.hits(counts[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], method = rep('asinh-qn-stand', sum(sel.f))))
}
isva.fit = isva(counts, maxq = 0.1) ######### ISVA considering only regions correlated with ISVs
Id = diag(ncol(counts))
for(i in 1:isva.fit$ncomp){
  mod = model.matrix( ~isva.fit$isv[, i] - 1)
  resid = counts %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid = normalize.quantiles(resid)
  hit.dat = compute.hits(resid[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], 
                               method = rep(paste('asinh-qn-stand-isv', i, '-sign', sep = ''), sum(sel.f))))
}
isva.fit = isva(counts, maxq = 1) ########## ISVA considering all regions
for(i in 1:isva.fit$ncomp){
  mod = model.matrix( ~ isva.fit$isv[, i] - 1)
  resid = counts %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid = normalize.quantiles(resid)
  hit.dat = compute.hits(resid[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], 
                               method = rep(paste('asinh-qn-stand-isv', i, sep = ''), sum(sel.f))))
}

pca.fit = prcomp(counts, center = F, scale = F) ########### PCA correction
nload = ncol(pca.fit$rotation)
for(i in 1:min(isva.fit$ncomp, nload)){
  mod = model.matrix( ~pca.fit$rotation[, 1:i] - 1)
  resid = counts %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid = normalize.quantiles(resid)
  hit.dat = compute.hits(resid[ov[sel, 2], ], genot, thresh.list, nperm)
  sel.f = !duplicated(hit.dat$fdrs)
  hit.stats = rbind(hit.stats, 
                    data.frame(hits = hit.dat$hits[sel.f], fdrs = hit.dat$fdrs[sel.f], 
                               method = rep(paste('asinh-qn-stand-pc', i, sep = ''), sum(sel.f))))
}
hit.stats = hit.stats[!is.na(hit.stats$fdrs) & !is.infinite(hit.stats$fdrs), ]
save(hit.stats, regions, file = file.path(plotdir, paste(outpref, '.RData', sep = '')))
hit.stats = hit.stats[hit.stats$fdrs < 0.25, ]
p = ggplot(hit.stats) + geom_line(aes(y = hits, x = fdrs, color = method, linetype = method)) + 
  ylab('# significant associations') + xlab('FDR') + 
  scale_linetype_manual(values = c(rep(c("solid", "dashed"), 7))) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10), legend.title = element_blank())
ggsave(file = file.path(plotdir, paste(outpref, '.png', sep = '')), p, width = 13.6, height = 11.8)