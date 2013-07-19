require(GenomicRanges)
require(MASS)

get.diff.count = function(filenames, qcut, fold.cut = NULL, is.log = T, sel.row = NULL){
  # Reads a set of files produced by doing pairwise comparisons with DESeq
  # and counts the number of times each region was found to be differential.
  # Each file must have a counts and a regions data.frame.
  # qcut is a qvalue cutoff for calling a region differential. 
  # fold.cut is a cutoff on fold change (if !is.log) or log-fold-change (if is.log)
  # for calling a region differential.
  
  if(is.null(fold.cut)){
    if(is.log){fold.cut = 0}else{fold.cut = 1}
  }
  tmp = unlist(strsplit(basename(filenames), '_'))
  indivs1 = tmp[seq(1, length(tmp), 4)]
  indivs2 = tmp[seq(2, length(tmp), 4)]
  indivs = unique(append(indivs1, indivs2))
  nindiv = length(indivs)
  pair.diff = array(0, dim = c(nindiv, nindiv))
  
  for(i in 1:length(filenames)){
    load(filenames[i])
    if(i == 1){
      diff.count = array(0, dim = c(nrow(regions), 1))
      names = rownames(counts)
    }else{
      stopifnot(all(rownames(counts) == names))
    }
    sel = !is.nan(regions$qval) & regions$qval < qcut
    if(is.log){
      sel = sel & (regions$log2Fold > fold.cut | regions$log2Fold < -fold.cut) 
    }else{
      sel = sel & (regions$fold > fold.cut | regions$fold < 1/fold.cut)
    }
    if(!is.null(sel.row)) sel = sel & sel.row
    diff.count = diff.count + as.numeric(sel)
    pair.diff[which(indivs == indivs1[i]), which(indivs == indivs2[i])] = sum(sel)
  }
  if(!is.null(sel.row)) diff.count = diff.count[sel.row]
  pair.diff = pair.diff + t(pair.diff)
  colnames(pair.diff) = indivs
  rownames(pair.diff) = indivs
  return(list(diff.count = diff.count, pair.diff = pair.diff))
}

avg.counts = function(indir, indivs, suf, len.norm = T, meta = NULL){
  # Concatenates count files produced by concat.reps. Concat reps concatenates replicates. 
  # This function can be used to concatenate the counts for all individuals for the same mark.
  # Counts from the same file (replicates) are AVERAGED after size factor normalization.
  # Counts are then normalized by the regions' lengths.
  # Filenames are assumed to be <indir>/SNYDER_HG19_<indiv[i]>_<suf>
  # Each file must have a counts data.frame and a size.factors vector with as many elements
  # as columns in the counts.
  
  counts.all = NULL
  for(i in 1:length(indivs)){
    filename = file.path(indir, paste('SNYDER_HG19', indivs[i], suf, sep = '_'))
    if(!file.exists(filename)){
      cat('File not found ', filename, '\n', file = stderr())
      next
    }
    print(filename)
    load(filename)
    if(!is.null(meta)){
      widths = meta$len
    }else{
      widths = regions$end - regions$start + 1
    }
    if(dim(counts)[2] > 1){
      counts.norm = t(apply(counts, 1, function(x) x / size.factors))
      counts.norm = rowMeans(counts.norm) # Take the mean over replicates after normalizing with the size factors.
    }else{counts.norm = counts / size.factors}
   if(is.null(counts.all)){
      counts.all = array(NaN, dim = c(dim(counts)[1], length(indivs)))
    }
    if(len.norm){
      counts.all[, i] = counts.norm / widths
    }else{
      counts.all[, i] = counts.norm
    }
  }
  #counts.all = data.frame(counts.all, row.names = rownames(counts))
  #colnames(counts.all) = indivs
  return(list(counts = counts.all, regions = regions))
}

regions.to.ranges = function(regions){
  ranges = GRanges(seqnames = Rle(regions$chr), 
                   ranges = IRanges(start = regions$start, end = regions$end),
                   strand = Rle(rep('+', dim(regions)[1])))
  return(ranges)
}

snps.to.ranges = function(regions){
  ranges = GRanges(seqnames = Rle(regions$chr), 
                   ranges = IRanges(start = regions$pos, end = regions$pos),
                   strand = Rle(rep('+', dim(regions)[1])))
  return(ranges)
}

load.avg.sig.data = function(region.file, sig.files, labels){
  regions = read.bed(region.file)
  nregions = dim(regions)[1]
  nfiles = length(sig.files)
  
  sig = array(0, dim = c(nregions, nfiles))
  print(dim(sig))
  for(i in 1:nfiles){
    print(basename(sig.files[i]))
    tmp.dat = read.table(sig.files[i], header = F)
    sig[, i] = tmp.dat[, 1]
  }
  colnames(sig) = labels
  return(list(regions = regions, signal = sig))
}

load.dip.data = function(region.file, sig.files, labels, reg = 'max'){
  regions = read.table(region.file, header = F, sep = '\t')
  regions = regions[, 1:3]
  colnames(regions) = c('chr', 'start', 'end')
  nregions = dim(regions)[1]
  nfiles = length(sig.files)
  
  dip.sig = array(0, dim = c(nregions, nfiles))
  left.sig = array(0, dim = c(nregions, nfiles))
  right.sig = array(0, dim = c(nregions, nfiles))
  for(i in 1:nfiles){
    print(basename(sig.files[i]))
    counts = read.table(sig.files[i], header = F, sep = '\t')
    dip.sig[, i] = counts[, 1]
    left.sig[, i] = counts[, 2]
    right.sig[, i] = counts[, 3]
  }
  if(reg == 'max'){
    m = colMaxs(dip.sig, na.rm = T)
    print(m)
    dip.sig = t(apply(dip.sig, 1, function(x) x / m))
    m = colMaxs(left.sig, na.rm = T)
    print(m)
    left.sig = t(apply(left.sig, 1, function(x) x / m))
    m = colMaxs(right.sig, na.rm = T)
    print(m)
    right.sig = t(apply(right.sig, 1, function(x) x / m))
  }
  dip.sig = data.frame(dip.sig)
  colnames(dip.sig) = labels
  left.sig = data.frame(left.sig)
  colnames(left.sig) = labels
  right.sig = data.frame(right.sig)
  colnames(right.sig) = labels
  return(list(regions = regions, dip.sig = dip.sig, left.sig = left.sig, right.sig = right.sig))
}

plot.pcs = function(data, eigen, dev, labels, groups = labels, all = F, ndim = 4){
  ###### Plot PC1 vs PC2 etc
  nplots = min(ndim - 1, dim(data)[2] - 1)
  nrows = dim(data)[1]
  
  if(all){
    for(i in 1:nplots){
      if(i == nplots + 1) next
      for(j in (i + 1):(nplots + 1)){
        tmp = data.frame(x = data[, i], y = data[, j], type = rep(paste('C', i, ' vs C', j, sep = ''), nrows), labels = labels, groups = groups)
        if(i == 1 && j == 2){
          d = tmp
        }else{
          d = rbind(d, tmp)
        }   
      }
    }
  }else{
    for(i in 1:nplots){
      tmp = data.frame(x = data[, i], y = data[, i + 1], type = rep(paste('C', i, ' vs C', i + 1, sep = ''), nrows), labels = labels, groups = groups)
      if(i == 1){
        d = tmp
      }else{
        d = rbind(d, tmp)
      }
    }
  }
  p1 = ggplot(d) + geom_point(aes(x = x, y = y, color = groups, shape = groups), size = 4) + theme_bw() + 
    geom_text(aes(x = x, y = y, label = labels, color = groups), show_guide = F, size = 4, vjust = 2) + 
    facet_wrap(~type, scales = 'free') + xlab('') + ylab('') + scale_color_manual('', values = get.pop.col(unique(groups))) + scale_shape_discrete('') +
    
    theme(strip.text.x = element_text(size = 12), legend.text = element_text(size = 12),
          #axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  ###### Plot eigenvalues 
#   pc.dat = data.frame(pca.fit$rotation[, 1:isva.fit$nsv])
#   pc.dat$indiv = rownames(pc.dat)
#   pc.dat = melt(pc.dat, id.vars = c('indiv'))
#   p3 = ggplot(pc.dat) + geom_bar(aes(x = indiv, y = value), position = 'dodge', stat = 'identity') + facet_wrap(~variable) +
#     xlab('') + ylab('PC value') + 
#     theme(strip.text.x = element_text(size = 12),
#           axis.text.x = element_text(size = 10, angle = -45, vjust = 1, hjust = 0))
  
  ###### Plot variances
  pc.labels = paste('PC', 1:ncol(data), sep = '')
  dev.dat = data.frame(dev = dev / sum(dev), pc = ordered(pc.labels, levels = pc.labels))
  p2 = ggplot(dev.dat) + geom_bar(aes(x = pc, y = dev), position = 'dodge', stat = 'identity') + xlab('') + ylab('') + 
    theme(axis.text.x = element_text(size = 12, angle = -45, vjust = 1, hjust = 0), axis.text.y = element_text(size = 12))
  
  return(list(p1 = p1, p2 = p2))
}

read.bed = function(bedfile){
  regions = read.table(bedfile, header = F, sep = '\t')
  regions = regions[, 1:3]
  regions[, 2] = regions[, 2] + 1
  colnames(regions) = c('chr', 'start', 'end')
  return(regions)
}

read.narrowPeak = function(peakfile){
  regions = read.table(peakfile, header = F, sep = '\t')
  regions = regions[, c(1, 2, 3, 7, 10)]
  colnames(regions) = c('chr', 'start', 'end', 'signal', 'summit')
  regions$summit = regions$start + regions$summit 
}

# Normalizes a set of counts using a linear fit between them (either robust regression or lowess regression)
# Input values:
# counts: counts (or other values) to normalize, one row per region, one column per sample
# bed.files: array of paths to peaks for each sample (column of counts). Only peaks that are common between two samples
# are used for the fit.
# Return value:
# A list with the following elements
# lr: log-ratio after the fit
# avg: log(a+b)/2 after the fit
# fits: coefficients for the fitted model, a ncol-by-2 array, where ncol is the number of columns of the input counts. 
fit.norm.counts = function(counts, bed.files, ref.idx = 1, take.log = T, fit = 'rlm'){
  ncol = dim(counts)[2]
  if(take.log) counts = log2(counts + 1)
  
  new.folds = counts
  new.folds[, ref.idx] = 0 # Log-ratio of 0 for the reference individual
  new.avg = counts
  npeaks = array(0, dim = c(ncol, 1))
  fits = array(0, dim = c(ncol, 2))
  
  # Find the file that has the peaks for the reference individual.
  ref.peaks = regions.to.ranges(read.bed(bed.files[ref.idx]))

  for(i in 1:ncol){
    if(i == ref.idx) next
    tmp.peaks = regions.to.ranges(read.bed(bed.files[i]))
    sel = !is.na(findOverlaps(ref.peaks, tmp.peaks, ignore.strand = T, select = 'first'))
    npeaks[i] = sum(sel)
    if(fit == 'rlm'){
      # Fit robust regression using the common peaks
      M = counts[sel, i] - counts[sel, ref.idx]
      A = (counts[sel, i] + counts[sel, ref.idx]) / 2
      b = rlm(M~A)$coefficients
      # Now adjust all peaks based on the fit
      resc = (2-b[2]) * counts[, i] /(2+b[2]) - 2*b[1]/(2+b[2]);
      new.folds[, i] = resc - counts[, ref.idx]
      new.avg[, i] = (resc + counts[, ref.idx]) / 2;
      fits[i, ] = b
    }else{
      # NOT IMPLEMENTED
      lo = lowess(c2~c1)
      x.lo = seq(min(c1), max(c1), length.out = 1000)
      y.lo = approx(lo$x, lo$y, x.lo)
    }
  }
  return(list(lr = new.folds, avg = new.avg, fits = fits, npeaks = npeaks))
}

# Computes the nucleotide diversity in a set of regions.
# The nucleotide diversity of the region is sum(f_i * f_j * x_ij)
# where the sum is taken over all pairs of haplotypes, f_i is the frequency of haplotype i and x_ij is the fraction of nucleotides
# in the region that differ between haplotypes i and j.
#
# Inputs:
# genot: boolean matrix SNP-by-haplotype. [i, j] is SNP i in haplotype j
# ov: overlaps between regions and SNPs. If ov[i] = j, then the i-th SNP in the genot matrix overlaps the j-th region
# (which corresponds to the j-th index in the widths array).
# eg. if ov is c(2,2,2,5,5), then the first three SNPs in genot overlap region 2 and the next 2 overlap region 5.
# Some regions might not have SNPs, in which case their indices won't appear in ov. ov MUST BE SORTED in non-decreasing order.
# widths: region widths
#
# Returning value:
# an array of the same length as widths with the nucleotide diversity for each region. If a region had no SNPs, its diversity will be 0.
get.div = function(genot, ov, widths){
  nregions = length(widths) # The output will be the diversity for all regions, even those without SNPs.
  div = array(0, dim = c(nregions, 1))
  pind = unique(ov) # Only consider regions that have some SNPs
  for(i in 1:length(pind)){ # For each region...
    sel.pos = findInterval(c(pind[i] - 1, pind[i]), ov) # Indices of SNPs overlapping the pind[i]-th region
    snp.ind = (sel.pos[1] + 1):sel.pos[2] # pind[i] appears in ov, so this interval makes sense
    nsel = length(snp.ind) # Number of snps overlapping the region
    # Get the distinct haplotypes (i.e. distinct columns of genot[snp.ind, ]) and their occurence counts.
    # hap.counts will have one row per unique haplotype. The last column will be the number of occurrences of the haplotype.
    if(nsel > 1){
      hap.counts = count(data.frame(t(genot[snp.ind, ]))) # the count function will count the distinct rows, but we want columns
    }else{
      hap.counts = count(data.frame(genot[snp.ind, ])) # If there's only one SNP, then just count Ts and Fs.
    }
    nhaps = dim(hap.counts)[1] # Number of distinct haplotypes
    tot.haps = sum(hap.counts$freq) # Number of haplotypes
    stopifnot(tot.haps == dim(genot)[2])
    hap.counts$freq = hap.counts$freq / tot.haps
    for(j in 1:nhaps){
      if(j == nhaps) next
      for(k in (j + 1):nhaps){
        div[pind[i]] = div[pind[i]] + hap.counts$freq[j] * hap.counts$freq[k] * 
          sum(xor(hap.counts[j, 1:nsel], hap.counts[k, 1:nsel])) / widths[pind[i]]
      }
    }
  }
  return(2 * div)
}

get.div.par = function(genot, ov, widths, nchunks = 10){
  vals = unique(ov)
  nvals = length(vals)
  chunk.size = ceiling(nvals / nchunks)
  div = foreach(i = 1:nchunks, .combine = 'cbind', .inorder = F) %dopar%{
    chunk.start = (i - 1) * chunk.size + 1
    chunk.end = min(i * chunk.size, nvals)
    interval = findInterval(c(vals[chunk.start] - 1, vals[chunk.end]), ov)
    #interval = ov >= vals[chunk.start] & ov <= vals[chunk.end]
    #if(any(interval)){
    if(interval[1] < interval[2]){
      interval = (interval[1] + 1):interval[2]
      #print(dim(ov[ov[,1] >= vals[chunk.start] & ov[, 1] <= vals[chunk.end], ]))
      get.div(genot[interval,], ov[interval], widths)
    }
 }
 return(rowSums(div, na.rm = T))
}

get.pop = function(indivs){
  pop = indivs
  for(i in 1:length(indivs)){
    pop[i] = switch(as.character(indivs[i]),
                    'GM12878' = 'CEU',
                    'GM12891' = 'CEU',
                    'GM12892' = 'CEU',
                    'GM12890' = 'CEU',
                    'GM10847' = 'CEU',
                    'GM18505' = 'YRI',
                    'GM18526' = 'Asian',
                    'GM18951' = 'Asian',
                    'GM19099' = 'YRI',
                    'GM19193' = 'YRI',
                    'GM19238' = 'YRI',
                    'GM19239' = 'YRI',
                    'GM19240' = 'YRI',
                    'SNYDER' = 'CEU',
                    'MS1' = 'CEU',
                    'GM18486' = 'YRI',
                    'GM2255' = 'San',
                    'GM2588' = 'San',
                    'GM2610' = 'San',
                    'GM2630' = 'San',
                    'HG2255' = 'San',
                    'HG2588' = 'San',
                    'HG2610' = 'San',
                    'HG2630' = 'San')
  }
  return(pop)
}

get.sex = function(indivs){
  pop = indivs
  for(i in 1:length(indivs)){
    pop[i] = switch(as.character(indivs[i]),
                    'GM12878' = 'female',
                    'GM12891' = 'male',
                    'GM12892' = 'female',
                    'GM12890' = 'female',
                    'GM10847' = 'female',
                    'GM18505' = 'female',
                    'GM18526' = 'female',
                    'GM18951' = 'female',
                    'GM19099' = 'female',
                    'GM19193' = 'female',
                    'GM19238' = 'male',
                    'GM19239' = 'female',
                    'GM19240' = 'female',
                    'SNYDER' = 'male',
                    'MS1' = 'male',
                    'GM18486' = 'male',
                    'GM2255' = 'male',
                    'GM2588' = 'male',
                    'GM2610' = 'male',
                    'GM2630' = 'male',
                    'HG2255' = 'male',
                    'HG2588' = 'male',
                    'HG2610' = 'male',
                    'HG2630' = 'male')
  }
  return(pop)
}

get.pop.col = function(pop){
  colors = as.character(pop)
  for(i in 1:length(pop)){
    colors[i] = switch(as.character(pop[i]),
                       'CEU' = '#E69F00',
                       'YRI' = '#009E73',
                       'Asian' = '#0072B2',
                       'San' = '#CC79A7')
  }
  names(colors) = pop
  return(colors)
}

# Given overlaps between regions and AS SNPs, find regions where all SNPs are in the same direction.
# ov.mat: Col 1 - region index, Col 2 - direction (eg log(ref/alt))
# Return value:
# list with elements:
# assoc: selected subset of Col 1 of ov.mat
# ratios: signs of SNPs for selected subset
get.consistent.ov = function(ov.mat){
  colnames(ov.mat) = c('region.idx', 'ratio')
  # Col 1 - region idx, Col 2 - T if all SNPs overlapping the region have ref > alt
  sel.pos = cast(ov.mat, region.idx~., function(x) all(x > 0), value = 'ratio')
  # Col 1 - region idx, Col 2 - T if all SNPs overlapping the region have ref < alt
  sel.neg = cast(ov.mat, region.idx~., function(x) all(x < 0), value = 'ratio')
  # Indices of regions such that all overlapping SNPs have the same direction
  sel.assoc = append(sel.pos[sel.pos[,2], 1], sel.neg[sel.neg[,2], 1])
  ratios = append(rep(1, sum(sel.pos[, 2])), rep(-1, sum(sel.neg[, 2])))
  return(list(assoc = sel.assoc, ratios = ratios))
}

order.marks = function(marks, sub.rna = T){
  marks = as.character(marks)
  if(sub.rna){
    marks[marks == 'RNA'] = 'polyA-RNA'
    marks[marks == 'RZ'] = 'RNA' 
  }
  levels = c('CTCF', 'PU1', 'SA1', 'H3K27AC', 'H3K4ME1', 'H3K4ME3', 'H2AZ', 'H3K9AC', 
             'RNA', 'polyA-RNA', 'POL4H8', 'EPOL', 'H3K36ME3', 'BUB', 'H3K27ME3','H3K9ME3', 'INPUT')
  if(!all(marks %in% levels)){
    return(marks)
  }
  return(ordered(factor(marks, levels = levels)))
}

get.hex = function(cols){
  hex.cols = as.character(cols)
  for(i in 1:length(cols)){
    h = as.numeric(unlist(strsplit(as.character(cols[i]), ',')))
    hex.cols[i] = rgb(h[1], h[2], h[3], maxColorValue = 255)
  }
  return(hex.cols)
}

mark.colors = function(marks){
  marks = as.character(marks)
  cols = marks
  for(i in 1:length(marks)){
    cols[i] = switch(as.character(marks[i]),
                    'CTCF' = rgb(0,0,0, maxColorValue = 255),
                    'PU1' = rgb(0,0,0, maxColorValue = 255),
                    'SA1' = rgb(104, 34, 139, maxColorValue = 255), # darkorchid4
                    'H3K27AC' = rgb(255,165,0, maxColorValue = 255),
                    'H3K4ME1' = rgb(255,215,0, maxColorValue = 255),
                    'H3K4ME3' = rgb(255,0,0, maxColorValue = 255),
                    'H2AZ' = rgb(255,127,80, maxColorValue = 255), # coral
                    'H3K9AC' = rgb(220,20,60, maxColorValue = 255),
                    'RNA' = rgb(0,205,205, maxColorValue = 255), # darkcyan
                    'polyA-RNA' = rgb(0,139,139, maxColorValue = 255),
                    'POL4H8' = rgb(0,200,0, maxColorValue = 255),
                    'EPOL' = rgb(0,200,0, maxColorValue = 255),
                    'H3K36ME3' = rgb(0,150,0, maxColorValue = 255),
                    'BUB' = rgb(0,150,0, maxColorValue = 255),
                    'H3K27ME3' = rgb(105,105,105, maxColorValue = 255),
                    'H3K9ME3' = rgb(169,169,169, maxColorValue = 255), 
                     rgb(0,0,0, maxColorValue = 255))
  }
  names(cols) = marks
  return(cols)
}

get.indivs = function(){
  return(c("GM10847", "GM12878", "GM12891", "GM12892", "GM12890", "MS1", 
           'GM18486', "GM18505", "GM19099", 'GM19193', "GM19238", "GM19239", "GM19240", "GM18526", "GM18951",
           'HG2255', 'HG2588', 'HG2610', 'HG2630'))
}
fix.indiv.names = function(indivs){
  indivs[indivs == 'SNYDER'] = 'MS1'
  indivs[indivs == 'GM2588'] = 'HG2588'
  indivs[indivs == 'GM2255'] = 'HG2255'
  indivs[indivs == 'GM2610'] = 'HG2610'
  indivs[indivs == 'GM2630'] = 'HG2630'
  return(indivs)
}
plot.tile = function(mat, x.ord.samples = NULL, y.ord.samples = NULL, xsep = NULL, ysep = NULL, ylabels = NULL, 
                     low = 'blue', high = 'red', mid = 'beige',
                     midpoint = 1, ycolor = 'black', xcolor = 'black', draw.y.line = T,
                     xcex = 12, ycex = 12, lcex = 14, xtitle = '', ytitle = ''){
  dat = data.frame(melt(as.matrix(mat)))
  if(!is.factor(dat$X1)){
    dat$X1 = -dat$X1
  }else{
    if(is.null(y.ord.samples)){
      dat$X1 = ordered(factor(dat$X1, levels = rownames(mat)[seq(length(rownames(mat)), 1, -1)]))
    }else{
      dat$X1 = droplevels(ordered(factor(dat$X1, levels = y.ord.samples[seq(length(y.ord.samples), 1, -1)])))
    }
    stopifnot(is.character(ycolor) || length(ycolor) == length(levels(dat$X1)))
    if(length(ycolor) == length(levels(dat$X1))) ycolor = ycolor[seq(length(levels(dat$X1)), 1, -1)]
  }
  if(!is.factor(dat$X2)){
    dat$X2 = -dat$X2
  }else{
    if(is.null(x.ord.samples)){
      dat$X2 = ordered(factor(dat$X2, levels = colnames(mat)))
    }else{
      dat$X2 = droplevels(ordered(factor(dat$X2, levels = x.ord.samples)))
    }
  }
  p = ggplot(dat) + geom_raster(aes(x = X2, y = X1, fill = value)) +
    scale_fill_gradient2(low = low, high = high, mid = mid, midpoint = midpoint)
  if(is.factor(dat$X2)){
    p = p + scale_x_discrete(expand = c(0, 0))
  }else{
    p = p + scale_x_continuous(expand = c(0, 0))
  }
  if(is.factor(dat$X1)){
    p = p + scale_y_discrete(expand = c(0, 0))
  }else{
    if(!is.null(ysep)){
      p = p + scale_y_continuous(breaks = ysep, labels = ylabels, expand = c(0, 0)) 
      if(draw.y.line) p = p + geom_hline(yintercept = -ysep, size = 0.6)
    }else{p = p + scale_y_continuous(expand = c(0, 0))}
  }
  p = p + theme_bw() + xlab(xtitle) + ylab(ytitle) + 
    theme(axis.text.x = element_text(size = xcex, angle = 90, vjust = 0.5, hjust = 1, color = xcolor),
          axis.title.x = element_text(size = lcex),
          axis.text.y = element_text(size = ycex, color = ycolor), axis.title.y = element_text(size = lcex),
          legend.text = element_text(size = lcex), legend.title = element_blank())
  return(p)
}
