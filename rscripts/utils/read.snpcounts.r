rm(list=ls())
library(GenomicRanges)
library(Matrix)
source('utils/sample.info.r')

# Reads SNP counts produced by AseQuantMultiRG.cpp.
# All count files read are assumed to have the SAME SNPs in the SAME order.
# There must exist an RData file (snp.pos.file) that has the SNPs that were used to generate the 
# counts files (see read.genotypes.r) as well as 
# a vector "black" that gives the SNPs that should be removed from the counts files.

indir = '../../rawdata/alleleCounts/san/' # count files should be here
filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*.counts.gz$', sep = ''), full.name = T)
outdir = file.path(indir, 'rdata')
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir))

snp.pos.file = '../../rawdata/variants/sanConsensus/snps/san.snps.RData'
load(snp.pos.file)
nrows = length(black) # All counts files read must have this number of rows.

samples = sample.info(filenames, '.counts.gz')
uniq.indivs = unique(as.character(samples$indiv))
overwrite = F

# Counts MUST be in this order in all the counts files
pat.idx = 1
mat.idx = 2
amb.idx = 3
counts.colnames = c('ref', 'alt', 'tot')

for(i in 1:length(uniq.indivs)){

  sel.files = filenames[samples$indiv == uniq.indivs[i]]  
  for(f in sel.files){
    outfile = file.path(outdir, gsub('.counts.*', '.RData', basename(f)))
    if(file.exists(outfile) && !overwrite){
      cat(paste('File', outfile, 'exists. Skipping sample.\n'), file = stderr())
      next
    }
    cat(f, '\n', file = stderr())
    start = proc.time()
    tmp.file = paste(f, 'tmp', sep = '.')
    # Separate the counts into individual columns
    system(paste('python', '../python/parseCounts.py', f, '>', tmp.file))
    # Read all the counts into one big vector
    tab = scan(tmp.file, what = integer(), sep = '\t') #header = T, sep = '\t', comment.char = '')
    cat('Finished reading counts:', proc.time()[3] - start[3], 'sec\n', file = stderr())
    start = proc.time()
    len = length(tab)
    
    # Read counts
    # 1st dimension is SNP, 2nd is read group, 3rd is fwd vs rev strand
    ref.tmp = array(0, dim = c(nrows, 3, 2))
    alt.tmp = array(0, dim = c(nrows, 3, 2))
    oth.tmp = array(0, dim = c(nrows, 3, 2))
    
    for(i in 0:2){
      for(j in 0:1){
        ref.tmp[, i + 1, j + 1] = tab[seq(6 * i + 3 * j + 1, len, 18)]
        alt.tmp[, i + 1, j + 1] = tab[seq(6 * i + 3 * j + 2, len, 18)]
        oth.tmp[, i + 1, j + 1] = tab[seq(6 * i + 3 * j + 3, len, 18)]
      }
    }
    cat('Finished building allele count matrix', proc.time()[3] - start[3], 'sec\n', file = stderr())
    
#     nrows = dim(tab)[1]
#     start.col = 4 # Counts start after this column
#     mat.idx = which(colnames(tab) == 'maternal') - start.col
#     pat.idx = which(colnames(tab) == 'paternal') - start.col
#     amb.idx = which(colnames(tab) == 'ambiguous') - start.col
#     
#     tmp.ranges = GRanges(seqnames = Rle(tab$X.CHROM), ranges = IRanges(start = tab$POS, end = tab$POS, width = rep(1, nrows)),
#                          strand = Rle(rep('+', nrows)))
#     cat('Building allele count matrix\n', file = stderr())
#     for(i in 1:3){
#       fields = unlist(strsplit(as.character(tab[, start.col + i]), ','))
#       for(j in 0:1){
#         ref.tmp[, i, j + 1] <- as.numeric(fields[seq(3 * j + 1, length(fields), 6)])
#         alt.tmp[, i, j + 1] <- as.numeric(fields[seq(3 * j + 2, length(fields), 6)])
#         oth.tmp[, i, j + 1] <- as.numeric(fields[seq(3 * j + 3, length(fields), 6)])
#       }
#     }
    
    ref.tmp = ref.tmp[!black, , ]
    alt.tmp = alt.tmp[!black, , ]
    oth.tmp = oth.tmp[!black, , ]
    # Make sure that the resulting SNPs are EXACTLY the snps in the genotype file.
    # stopifnot(all(tab[good, 2] == geno.info$pos))
    
    # Remove SNPs where there is evidence in the same allele from both paternal and maternal reads
    # This might happen when there are several SNPs close to each other
    bad = rowSums(ref.tmp[, mat.idx, ]) > 0 & rowSums(ref.tmp[, pat.idx, ]) > 0
    bad = bad | rowSums(alt.tmp[, mat.idx, ]) > 0 & rowSums(alt.tmp[, pat.idx, ]) > 0
    
    ref = rowSums(ref.tmp[, c(mat.idx, pat.idx), ]) # all non-ambiguous ref (this will be the same as the ref maternal or ref paternal for non-bad SNPs)
    alt = rowSums(alt.tmp[, c(mat.idx, pat.idx), ])
    tot = rowSums(ref.tmp + alt.tmp + oth.tmp)
    counts = Matrix(data = cbind(ref, cbind(alt, tot)))
    
    # Create snp.info with bad pos
    snp.info = list()
    snp.info$bad = Matrix(bad)
    save(counts, snp.info, file = outfile)
  }
}