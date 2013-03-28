rm(list=ls())

# Reads CAGT results and writes them as an RData file.

clust.pref = '../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/SNYDER_HG19_all_H3K27AC_AS_AT_H3K27AC_cagt_'
clust = read.table(paste(clust.pref, 'clusters.bed', sep = ''), header = F, sep = '\t')
clust = data.frame(chr = clust[, 1], pos = clust[, 3], name = clust[, 4],
                       flip = clust[,7] == 1, overseg = clust[, 8], clust = clust[, 9])

hc.centroid.dist = as.matrix(read.table(paste(clust.pref, 'clustDist.txt', sep = ''), header = F, sep = '\t'))
save(clust, hc.centroid.dist, file = paste(clust.pref, 'results.RData', sep = ''))