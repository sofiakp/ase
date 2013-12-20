rm(list=ls())
#library('foreach')
#library('doMC')
library(GenomicRanges)
source('utils/binom.val.r')
source('utils/sample.info.r')
source('utils/deseq.utils.r')

############# Parameter setting 
# Files with counts and pvalues. Each file should have a snp.info list with the following elements:
# a (sparse) vector "bad" of problematic SNPs (see combine.reps) that were not included in the p-value computation.
# a (sparse) vector of log10 p-values
indir = '../../rawdata/alleleCounts/san/rdata/reps/'
outdir = file.path(indir, 'qvals') 
# Indicators of significant SNPs will go here.
if(!file.exists(file.path(outdir, 'hitLists'))) dir.create(file.path(outdir, 'hitLists'), recursive = T)

# INPUT files will be used to filter out SNPs with (even moderate) p-values in the control ChIP.
input.files = list.files(indir, pattern = 'SNYDER_HG19_.*_INPUT_rep.RData', full.names = T)
ninput = length(input.files)

# Peak files will be used to mark (but not filter before q-value computation) SNPs that are not within 
# enriched regions. These might result from high levels of background noise.
peak.files = list.files('../../rawdata/signal/combrep/peakFiles/', pattern = 'SNYDER_HG19_.*.encodePeak.gz', full.names = T)

# SNP positions. All count files read MUST correspond to these positions.
snp.pos.file = '../../rawdata/variants/sanConsensus/snps/san.snps.RData'
load(snp.pos.file)
snp.ranges = snps.to.ranges(snp.pos)

# For each individual, there must be a file <geno.dir>/<indiv>.snps.RData, with genotype information.
geno.dir = '../../rawdata/variants/sanConsensus/snps/'

q.cut = -2 # Log-10 pvalue cutoff
overwrite = F
males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486')
############# End parameter setting

for(n in 1:ninput){
  # Load INPUT and genotype file for the individual
  indiv = as.character(sample.info(input.files[n], '.RData')$indiv)
  maskfile = file.path(geno.dir, paste(indiv, '.snps.RData', sep = ''))
  
  load(maskfile)
  het = as.vector(geno.info$mat != geno.info$pat)
  mask = as.vector(geno.info$mask)
   
  load(input.files[n])
  stopifnot(length(mask) == dim(counts)[1])
  inputhit = as.vector(snp.info$pval < -2)
  
  # CHANGE THE PATTERN HERE IF YOU WANT TO ONLY GET A SUBSET OF MARKS
  filenames = list.files(indir, pattern = paste('SNYDER_HG19_', indiv, '.*rep.*\\.RData', sep = ''), full.name = F)
  marks = as.character(sample.info(filenames, '.RData')$mark)
  nfiles = length(filenames)
  
  if(nfiles < 1) next
  
  # Read each file for the individual whose input you just read.
  for(i in 1:nfiles){
    outfile = file.path(outdir, filenames[i])
    hits.outfile = file.path(outdir, 'hitLists', gsub('.RData', '.hitInd.RData', basename(filenames[i])))
    
    if(!grepl('INPUT', filenames[i]) && (overwrite || !file.exists(outfile))){
      print(filenames[i]) 
      peak.file.idx = grep(paste('SNYDER_HG19_', indiv, '_', marks[i], '.*\\.encodePeak.gz', sep = ''), peak.files)
      if(length(peak.file.idx) == 0){
        cat('Missing BED file with enriched regions.\n')
        next
      }else{
        if(length(peak.file.idx) > 1){
          peak.file.idx = peak.file.idx[!grepl('reconcile.dedup_VS_', peak.files[peak.file.idx])]
        }
        cat('Reading peaks from ', peak.files[peak.file.idx], '\n')
        peak.tab = read.table(peak.files[peak.file.idx], header = F)
        peaks = GRanges(seqnames = Rle(peak.tab[,1]), 
                        ranges = IRanges(start = peak.tab[, 2] + 1, end = peak.tab[, 3]),
                        strand = Rle(rep('+', dim(peak.tab)[1])))  
      }      
      load(file.path(indir, filenames[i]))
      stopifnot(length(mask) == dim(counts)[1])
      # Sites to use for q-value correction: heterozygous, non-masked, not "weird" (bad), without AS evidence in the INPUT
      pass = het & !mask & !inputhit & as.vector(!snp.info$bad)
      qval = array(0, dim = c(length(pass), 1))
      qval[pass] = log(p.adjust(10^(as.vector(snp.info$pval)[pass]), method = 'BH'), base = 10) # Convert the p-values to full for efficiency
      
      snp.info$pass = Matrix(pass)
      snp.info$qval = Matrix(qval)
      snp.info$enriched = Matrix(snp.ranges %in% peaks)
      save(counts, snp.info, file = outfile)  
      
      snp.info$hits = Matrix(snp.info$qval < q.cut)
      
      save(snp.info, file = hits.outfile)
    }
  }
}

# for(n in 1:ninput){
#   indiv = sample.info(input.files[n], '.RData')$indiv 
#   maskfile = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/', paste(indiv, '.snps.RData', sep = ''))
#   
#   load(maskfile)
#   ref = geno.info$ref
#   alt = geno.info$alt
#   mat = geno.info[[paste(indiv, 'mat', sep = '.')]]
#   pat = geno.info[[paste(indiv, 'pat', sep = '.')]]
#   mask = geno.info[[paste(indiv, 'mask', sep = '.')]]
#   phased = geno.info[[paste(indiv, 'phased', sep = '.')]]
#   het = mat != pat
#   snp.ranges = GRanges(seqnames = Rle(geno.info$chr), 
#                        ranges = IRanges(start = geno.info$pos, end = geno.info$pos),
#                        strand = Rle(rep('+', dim(geno.info)[1]))) 
#   rm(geno.info)
#   
#   load(input.files[n])
#   stopifnot(length(mask) == dim(snp.info)[1])
#   inputhit = !is.nan(snp.info$pval) & snp.info$pval < 0.01
#   
#   filenames = list.files(indir, pattern = paste('SNYDER_HG19_', indiv, '.*rep.*\\.counts\\.RData', sep = ''), full.name = F)
#   marks = as.character(sample.info(filenames, '.count.RData')$mark)
#   nfiles = length(filenames)
#   if(nfiles < 1) next
#   for(i in 1:nfiles){
#     print(filenames[i])
#     peak.file.idx = grep(paste('SNYDER_HG19_', indiv, '_', marks[i], '.*\\.bed', sep = ''), peak.files)
#     if(length(peak.file.idx) != 1){
#       cat('Missing BED file with enriched regions.\n')
#       next
#     }else{
#       peak.tab = read.table(peak.files[peak.file.idx], header = F)
#       peaks = GRanges(seqnames = Rle(peak.tab[,1]), 
#                       ranges = IRanges(start = peak.tab[, 2] + 1, end = peak.tab[, 3]),
#                       strand = Rle(rep('+', dim(peak.tab)[1])))  
#     }
#     outfile = file.path(outdir, filenames[i])
#     hits.outfile = file.path(outdir, 'hitLists', gsub('counts', 'hits', basename(filenames[i])))
#     het.outfile = file.path(outdir, 'hitLists', gsub('counts', 'het', basename(filenames[i])))
#     if(!grepl('INPUT', filenames[i]) && (overwrite || !file.exists(outfile))){
#       load(file.path(indir, filenames[i]))
#       stopifnot(length(mask) == dim(snp.info)[1])
#       is.male = indiv %in% males
#       pass = het & !mask & !inputhit
#       qval = array(NaN, dim = c(length(pass), 1))
#       qval[pass] = p.adjust(snp.info$pval[pass], method = 'BH')
#       
#       # Remove SNPs where there is evidence in the same allele from both paternal and maternal reads
#       # This might happen when there are several SNPs close to each other
#       # These will not be masked out, but I won't report them as AS
#       bad = rowSums(counts[, grep('pat.*ref', colnames(counts))]) > 0 & rowSums(counts[, grep('mat.*ref', colnames(counts))]) > 0
#       bad = bad | (rowSums(counts[, grep('pat.*alt', colnames(counts))]) > 0 & rowSums(counts[, grep('mat.*alt', colnames(counts))]) > 0)
#       #bad = bad | rowSums(counts) < quantile(rowSums(counts), 0.4)
#       qval[bad] = NaN
#       pass = pass & !bad
#       snp.info$qval = qval
#       snp.info$mat = mat
#       snp.info$pat = pat
#       snp.info$mask = mask
#       snp.info$phased = phased
#       snp.info$enriched = snp.ranges %in% peaks
#       save(counts, snp.info, sample, file = outfile)  
#       
#       # Now write hit lists
#       ref.counts = rowSums(counts[, grep('ref', colnames(counts))])
#       alt.counts = rowSums(counts[, grep('alt', colnames(counts))])
#       
#       tmp.counts = data.frame(chr = snp.info$chr, pos = snp.info$pos, ref.counts = ref.counts, alt.counts = alt.counts, 
#                               mat = snp.info$mat, pat = snp.info$pat, phased = snp.info$phased, enriched = snp.info$enriched)
#       het.counts = tmp.counts[pass, ]
#       signif = snp.info$qval[pass] < q.cut
#       hits = het.counts[signif, ]
#       het.counts$signif = signif
#       print(sum(signif))
#       save(het.counts, file = het.outfile)
#       save(hits, file = hits.outfile)
#     }
#   }
# }