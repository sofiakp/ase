rm(list=ls())
library(MASS)
library(affy)
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

# Computes a robust fit in the mean vs log-ratio plot and plots the points before and after the fit
plot.fit = function(c1, c2, name1, name2, outfile){
  M = c2 - c1
  A = (c1 + c2) / 2
  b = rlm(M ~ A)$coefficients
  print(b)
  resc = (2-b[2]) * c2 /(2+b[2]) - 2*b[1]/(2+b[2])
  M.new = resc - c1
  A.new = (resc + c1) / 2;
  
  sample = sample(1:length(M), min(length(M), 5000))
  intercepts = array(NA, dim = c(2 * length(sample), 1))
  intercepts[1] = b[1]
  slopes = array(NA, dim = c(2 * length(sample), 1))
  slopes[1] = b[2]
  dat = data.frame(x = append(A[sample], A.new[sample]), y = append(M[sample], M.new[sample]),
                   i = intercepts, s = slopes,
                   type = factor(append(rep('before', length(sample)), rep('after', length(sample)))))
  line.dat = data.frame(s = c(b[2], NA), i = c(b[1], NA), type = factor('before', 'after'))
  p = ggplot(dat) + geom_point(aes(x = x, y = y)) + geom_abline(aes(intercept = i, slope = s), color = 'green') + facet_wrap(~type) +
    xlab('log(a*b)/2') + ylab('log(b/a)')
  ggsave(outfile, p)
}

plot.comp = function(c1a, c2a, c1b, c2b, name1, name2, outfile){
  M = c2a - c1a
  A = (c1a + c2a)/2
  M.new = c2b - c1b
  A.new = (c1b + c2b)/2
  
  sample = sample(1:length(M), min(length(M), 5000))
  dat = data.frame(x = append(A[sample], A.new[sample]), y = append(M[sample], M.new[sample]),
                   type = factor(append(rep('before', length(sample)), rep('after', length(sample)))))
  p = ggplot(dat) + geom_point(aes(x = x, y = y)) + facet_wrap(~type) + xlab('log(a*b)/2') + ylab('log(b/a)')
  ggsave(outfile, p)  
}

signal.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/')
plotdir = file.path(signal.dir, 'plots', 'qn_fits')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)

region.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/merged/SNYDER_HG19_H3K27AC_merged.encodePeak.gz')
signal.files = list.files(signal.dir, pattern = 'SNYDER_HG19_H3K27AC_merged_AT_.*txt', full.names=T)
indivs = unique(gsub('_H3K27AC|.txt', '', gsub('SNYDER_HG19_H3K27AC_merged_AT_SNYDER_HG19_', '', basename(signal.files))))
nindivs = length(indivs)
counts.dat = load.avg.sig.data(region.file, signal.files, indivs)
regions = counts.dat$regions
counts = log2(counts.dat$signal + 1)

# Bed files with individual specific peaks. Fits between different individuals are performed on the common peaks only.
bed.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles')

ref.indiv = 2 # Reference individual. The counts for every other individual will be fit against this one.
# Find the file that has this individual's peaks.
peak.file = list.files(bed.dir, pattern = paste('SNYDER_HG19_', indivs[ref.indiv], '_H3K27AC_VS_.*', sep = ''), full.names = T)
stopifnot(length(peak.file) == 1)
ref.peaks = regions.to.ranges(read.bed(peak.file))

for(i in 1:nindivs){
  if(i == ref.indiv) next
  # Find the file that has the peaks for the specific individual
  peak.file = list.files(bed.dir, pattern = paste('SNYDER_HG19_', indivs[i], '_H3K27AC_VS_', sep = ''), full.names = T)
  stopifnot(length(peak.file) == 1)
  tmp.peaks = regions.to.ranges(read.bed(peak.file))
  sel = !is.na(findOverlaps(ref.peaks, tmp.peaks, ignore.strand = T, select = 'first'))
  plot.comp(counts[sel, ref.indiv], counts[sel, i], indivs[ref.indiv], indivs[i],
           file.path(plotdir, paste(indivs[ref.indiv], 'vs', indivs[i], 'manorm_fit.png', sep = '_')))
}

# norm.counts = normalize.quantiles(counts)
# for(i in 1:nindivs){
#   if(i == ref.indiv) next
#   plot.comp(counts[, ref.indiv], counts[, i], norm.counts[, ref.indiv], norm.counts[, i], indivs[ref.indiv], indivs[i],
#            file.path(plotdir, paste(indivs[ref.indiv], 'vs', indivs[i], 'qn_fit.png', sep = '_')))
# }
# # Reads counts at peaks and plots the scatter plots for pairs of individuals. It also computes lowess and robust regression fits
# # on the counts.
# 
# # Computes LOWESS and robust regression fits of two sets of counts
# plot.fit = function(c1, c2, name1, name2, outfile){
#   b.rlm = rlm(c2~c1)$coefficients
#   lo = lowess(c2~c1)
#   x.lo = seq(min(c1), max(c1), length.out = 1000)
#   y.lo = approx(lo$x, lo$y, x.lo)
#   png(outfile)
#   plot(c1, c2, cex = 0.3, xlab = name1, ylab = name2, 
#        main = paste('Normalized asinh counts'))
#   lines(x.lo, y.lo$y, col = 'red', lwd = 1)
#   abline(b.rlm, col = 'green')
#   legend("topleft", c('lowess', 'robust'), cex=1, lty=1, col = c('red', 'green'), bty = 'n');
#   dev.off()
# }
# 
# # Merged files with all the replicates for each individual. Each file in this directory has a matrix of counts
# # (regions x replicates). It also contains size factors for each replicate. The counts are fragment counts in a common set
# # of peaks for each mark.
# peak.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/countsAtPeaksBroad/repsComb/')
# plotdir = file.path(peak.dir, 'plots', 'linfits')
# if(!file.exists(plotdir)) dir.create(plotdir)
# indivs = unique(as.character(sample.info(list.files(peak.dir, pattern = 'SNYDER_HG19_.*H3K27AC_0.RData'), '.RData')$indiv))
# nindivs = length(indivs)
# 
# # Bed files with individual specific peaks. Fits between different individuals are performed on the common peaks only.
# bed.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles')
# 
# ref.indiv = 2 # Reference individual. The counts for every other individual will be fit against this one.
# # Find the file that has this individual's peaks.
# peak.file = list.files(bed.dir, pattern = paste('SNYDER_HG19_', indivs[ref.indiv], '_H3K27AC_VS_.*', sep = ''), full.names = T)
# stopifnot(length(peak.file) == 1)
# ref.peaks = regions.to.ranges(read.bed(peak.file))
# # For the reference individual, plot all fits between all pair of replicates.
# count.file = list.files(peak.dir, pattern = paste('SNYDER_HG19', indivs[ref.indiv], 'H3K27AC_0.RData', sep = '_'), full.names = T)
# stopifnot(length(count.file) == 1)
# load(count.file)
# nreps = dim(counts)[2]
# counts = t(apply(counts, 1, function(x) x / size.factors))
# #counts = apply(counts, 2, function(x) x / (regions$end - regions$start + 1))
# counts = asinh(counts) # Variance stabilize data
# 
# for(j in 1:(nreps - 1)){
#   for(k in (j + 1):nreps){
#     plot.fit(counts[, j], counts[, k], colnames(counts)[j], colnames(counts)[k], 
#              file.path(plotdir, paste(colnames(counts)[j], 'vs', colnames(counts)[k], 'linfit.png', sep = '_')))
#   }
# }
# 
# # Now read the counts for all individuals. Counts will be averaged across replicates. 
# counts.dat = avg.counts(peak.dir, indivs, 'H3K27AC_0.RData', len.norm = F, meta = NULL)
# regions = counts.dat$regions
# counts = counts.dat$counts
# counts = asinh(counts)
# 
# for(i in 1:nindivs){
#   if(i == ref.indiv) next
#   # Find the file that has the peaks for the specific individual
#   peak.file = list.files(bed.dir, pattern = paste('SNYDER_HG19_', indivs[i], '_H3K27AC_VS_', sep = ''), full.names = T)
#   stopifnot(length(peak.file) == 1)
#   tmp.peaks = regions.to.ranges(read.bed(peak.file))
#   sel = !is.na(findOverlaps(ref.peaks, tmp.peaks, ignore.strand = T, select = 'first'))
#   plot.fit(counts[sel, ref.indiv], counts[sel, i], indivs[ref.indiv], indivs[i],
#            file.path(plotdir, paste(indivs[ref.indiv], 'vs', indivs[i], 'linfit.png', sep = '_')))
# }