require(GenomicRanges)
require(MASS)

get.diff.count = function(filenames, qcut){
  # Reads a set of files produced by doing pairwise comparisons with DESeq
  # and counts the number of times each region was found to be differential.
  # Each file must have a counts and a regions data.frame. 
  
  for(i in 1:length(filenames)){
    load(filenames[i])
    if(i == 1){
      diff.count = array(0, dim = c(dim(regions)[1], 1))
      names = rownames(counts)
    }else{
      stopifnot(all(rownames(counts) == names))
    }
    diff.count = diff.count + as.numeric(!is.nan(regions$qval) & regions$qval < qcut)
  }
  return(diff.count)
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

plot.pcs = function(data, eigen, dev, labels, groups = labels, all = F){
  ###### Plot PC1 vs PC2 etc
  nplots = min(3, dim(data)[2] - 1)
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
  p1 = ggplot(d) + geom_point(aes(x = x, y = y, color = groups, shape = groups), size = 5) + 
    geom_text(aes(x = x, y = y, label = labels, color = groups), show_guide = F, size = 5, vjust = 2) + 
    facet_wrap(~type, scales = 'free') + xlab('') + ylab('') + scale_color_manual('', values = get.pop.col(unique(groups))) + scale_shape_discrete('') +
    theme(strip.text.x = element_text(size = 14), legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  
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
  dev.dat = data.frame(dev = pca.fit$sdev / sum(dev), pc = ordered(pc.labels, levels = pc.labels))
  p2 = ggplot(dev.dat) + geom_bar(aes(x = pc, y = dev), position = 'dodge', stat = 'identity') + xlab('') + ylab('') + 
    theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 1, hjust = 0), axis.text.y = element_text(size = 14))
  
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
                    'GM18486' = 'YRI')
  }
  return(pop)
}

get.pop.col = function(pop){
  colors = pop
  for(i in 1:length(pop)){
    colors[i] = switch(as.character(pop[i]),
                       'CEU' = '#E69F00',
                       'YRI' = '#009E73',
                       'Asian' = '#0072B2')
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