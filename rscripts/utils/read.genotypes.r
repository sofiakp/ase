rm(list=ls())
library('GenomicRanges')

# VCF file with genotypes for the individuals in the population
# You'd better remove all the metainfo (## headers)!! Keep the header with the sample names though!
infile <- file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/supporting/trio/snps/test.vcf')
# Directory with blacklisted SNPs for each individual (<indiv>.blacklist.txt)
blackdir <- file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/trio/snps')
outfile <- file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/trio/snps/CEU.trio.2010_09.genotypes.hg19.r')
startcol = 10

geno.tab <- read.table(infile, header = T, sep = '\t', comment.char = '')
colnames(geno.tab) = gsub('^NA', 'GM', colnames(geno.tab))
geno.info <- data.frame(chr = geno.tab$X.CHROM, pos = geno.tab$POS, id = geno.tab$ID, ref = geno.tab$REF, alt = geno.tab$ALT)
nrows = dim(geno.info)[1]

# Make a GRanges object for easy overlapping
geno.ranges = GRanges(seqnames = Rle(geno.info$chr), ranges = IRanges(start = geno.info$pos, end = geno.info$pos, width = rep(1, nrows)),
                      strand = Rle(rep('+', nrows)))

for(i in 10:dim(geno.tab)[2]){
  indiv = colnames(geno.tab)[i]
  tmp.list = unlist(strsplit(as.character(geno.tab[, i]), ':'))
  genos = tmp.list[seq(1, length(tmp.list), 3)]
  geno.info[[paste(indiv, 'phased', sep = '.')]] = grepl('\\|', genos)
  genos = gsub('/', '\\|', genos)
  tmp.list <- unlist(strsplit(genos, '\\|'))
  geno.info[[paste(indiv, 'pat', sep = '.')]] = as.numeric(tmp.list[seq(1, length(tmp.list), 2)]) == 1
  geno.info[[paste(indiv, 'mat', sep = '.')]] = as.numeric(tmp.list[seq(2, length(tmp.list), 2)]) == 1
  maskfile = file.path(blackdir, paste(indiv, '.blacklist.txt', sep = ''))
  if(file.exists(maskfile)){
    cat(paste('Reading mask file for', indiv, '\n'), file = stderr())
    mask <- read.table(maskfile, header = F, sep = '\t', comment.char = '#')
    mask.ranges = GRanges(seqnames = Rle(mask[,1]), 
                          ranges = IRanges(start = mask[, 2], end = mask[, 2], width = rep(1, dim(mask)[1])),
                          strand = Rle(rep('+', dim(mask)[1])))    
    geno.info[[paste(indiv, 'mask', sep = '.')]] = geno.ranges %in% mask.ranges
  }
}
save(geno.info, file = outfile)
