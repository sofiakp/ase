rm(list=ls())
library(reshape)
library(ggplot2)
library(Matrix)
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')

clust.file = '../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_cagt_clusters.bed'
sig.file = '../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/extractSignal/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC.txt'
outdir = '../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/plots'
eq.size = F
orient.shapes = T

if(eq.size){
  outfile = file.path(outdir, paste(gsub('.bed', '', basename(clust.file)), '_heat_eq.png', sep = ''))
}else{
  outfile = file.path(outdir, paste(gsub('.bed', '', basename(clust.file)), '_heat.png', sep = ''))
}
clust.dat = read.table(clust.file, header = F, sep = '\t')
clust = data.frame(chr = clust.dat[, 1], pos = clust.dat[, 3], flip = clust.dat[,7] == 1, clust = clust.dat[, 9])
nclust = max(clust$clust) # number of clusters

signal = read.table(sig.file, header = F, sep = ',')

stopifnot(nrow(signal) == nrow(clust))
sig.len = ncol(signal)
signal[clust$flip, ] = signal[clust$flip, seq(sig.len, 1, -1)]
signal = scale(signal[clust$clust > 0, ])
clust = clust[clust$clust > 0, ]
nsig = nrow(clust)
middle = ceiling(sig.len / 2)

if(orient.shapes){
  for(n in 1:nclust){
    sel = clust$clust == n
    left = median(rowMaxs(signal[sel, 1:(middle - 1)], na.rm = T))
    right = median(rowMaxs(signal[sel, (middle + 1):sig.len], na.rm = T))
    if(left > right) signal[sel, ] = signal[sel, seq(sig.len, 1, -1)]
  }
}
if(eq.size){
  counts = table(clust$clust)
  msize = min(min(counts), 100)
  for(n in 1:nclust){
    sel = sample(1:counts[n], msize)
    if(n == 1){ sort.sig = signal[which(clust$clust == n)[sel], ]
    }else{sort.sig = rbind(sort.sig, signal[which(clust$clust == n)[sel], ])}
  }
  row.sep = seq(msize + 1, nrow(sort.sig), by = msize)
}else{
  sel = sample(1:nsig, min(nsig, 1000))
  sidx = sort(clust$clust[sel], index.return = T)
  sort.sig = signal[sel[sidx$ix], ]
  row.sep = which(!duplicated(sidx$x))
}

col.sep = c(middle)
png(filename = outfile, width = 8, height = 11, units="in", res=600)
plot.heatmap(sort.sig, row.cluster = F, col.cluster = F, show.dendro = "none", col.title = '', row.title = '', key = T,
             lab.row = NA, lab.col = NA, palette = brewer.pal(11,  "RdYlBu")[seq(9, 2, -1)], row.sep = row.sep, col.sep = col.sep, col.col = 'black',
             break.type='quantile', filt.thresh = NA, replace.na = F, margins = c(1,1)) #, to.file = file.path(plotdir, paste('biclust_', mark, '.png', sep = '')))
dev.off()

