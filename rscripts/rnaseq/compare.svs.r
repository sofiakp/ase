# Computes overlaps between sets of population specific regions identified with different methods

marks = c('CTCF', 'H3K27AC', 'H3K4ME1', 'H3K4ME3')
stats = array(NaN, dim = c(length(marks), 4))

for(i in 1:length(marks)){
  load(file.path('../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/', paste('SNYDER_HG19_all_reg', marks[i], 'qn_isvaNull.RData', sep = '_')))
  regions1 = regions[isva.fit$deg, ]
  
  load(file.path('../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/', paste('SNYDER_HG19', marks[i], 'qn_svaNull.RData', sep = '_')))
  regions2 = regions[sva.fit$deg, ]
  
  ov = findOverlaps(regions.to.ranges(regions1), regions.to.ranges(regions2), select = 'all', ignore.strand = T)
  ov.mat = cbind(queryHits(ov), subjectHits(ov))
  
  stats[i, ] = c(nrow(regions1), nrow(regions2), length(unique(ov.mat[, 1])) / nrow(regions1), length(unique(ov.mat[, 2])) / nrow(regions2))
  
  cat(marks[i], '\n')
  cat('# regions in set 1', nrow(regions1), '\n')
  cat('# regions in set 2', nrow(regions2), '\n')
  cat('Fraction of set 1 that overlap', length(unique(ov.mat[, 1])) / nrow(regions1), '\n')
  cat('Fraction of set 2 that overlap', length(unique(ov.mat[, 2])) / nrow(regions2), '\n')  
}

dat = data.frame(stats)
colnames(dat) = c('n1', 'n2', 'ov1', 'ov2')
dat$mark = marks
