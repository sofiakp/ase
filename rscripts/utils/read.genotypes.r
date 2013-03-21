rm(list=ls())
library('GenomicRanges')
library(Matrix)

# Genotype files, one per individual. The format should be: 
# CHROM, POS, ID, REF, ALT, GENOTYPE
# The name of the files should be <indiv>suf.
# The output will be written in <indiv>.snps.RData.
suf = '.san.txt'
snp.dir = '../../rawdata/variants/sanConsensus/snps/'
snp.files = list.files(snp.dir, pattern = paste(suf, sep = ''), full.names = T)
# Bed files with regions to mask. Should be named <indiv>.blacklist.bed.
mask.dir = '../../rawdata/variants/sanConsensus/masks/'
overwrite = F

# SNP positions will be written (or read from) here. The genotype files for all individuals are 
# assumed to have the SAME SNPs in the SAME ORDER. However, the alternative allele for each individual
# might differ. Therefore, we store the SNP positions, IDs, and reference alleles at one place, 
# and for each individual only store the alternative allele, genotypes, and individual-specific masks.
snp.pos.file = file.path(snp.dir, 'san.snps.RData')
snp.pos = NULL

# These will be completely removed from the output files
bad.tab = read.table(file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomes_local/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed'), header = F, sep = '\t')
bad.ranges = GRanges(seqnames = Rle(bad.tab[,1]), 
                     ranges = IRanges(start = bad.tab[, 2] + 1, end = bad.tab[, 3]),
                     strand = Rle(rep('+', dim(bad.tab)[1])))

for(f in snp.files){
  outfile = gsub(suf, '.snps.RData', f)
  if(file.exists(outfile) && !overwrite){
    cat(paste('File', outfile, 'exists. Skipping sample.'), file = stderr())
    next
  }
  # Read the mask (blacklisted regions) for the corresponding individual.
  indiv = gsub(suf, '', basename(f))
  mask.file = file.path(mask.dir, paste(indiv, '.blacklist.bed', sep = ''))
  if(!file.exists(mask.file)){
    cat(paste('Missing file:', mask.file), file = stderr())
    next
  }
  mask = read.table(mask.file, header = F, sep = '\t', comment.char = '#')
  mask.ranges = GRanges(seqnames = Rle(mask[,1]), 
                        ranges = IRanges(start = mask[, 2] + 1, end = mask[, 3]),
                        strand = Rle(rep('+', dim(mask)[1])))
  
  cat(paste('Reading data for', indiv, '\n'), file = stderr())
  
  # Read the genotypes
  geno.tab = read.table(f, header = T, sep = '\t', comment.char = '')
  
  # Load the pre-stored snp positions if they exist. Otherwise, create this file from the geno.tab you just read.
  if(is.null(snp.pos)){
    if(file.exists(snp.pos.file)){
      load(snp.pos.file)
      nrows = dim(snp.pos)[1]
      
      # Make a GRanges object for easy overlapping
      geno.ranges = GRanges(seqnames = Rle(snp.pos$chr), ranges = IRanges(start = snp.pos$pos, end = snp.pos$pos, width = rep(1, nrows)),
                            strand = Rle(rep('+', nrows)))
    }else{
      snp.pos = data.frame(chr = geno.tab$X.CHROM, pos = geno.tab$POS, id = geno.tab$ID, ref = geno.tab$REF)
      nrows = dim(snp.pos)[1]
      
      # Make a GRanges object for easy overlapping
      geno.ranges = GRanges(seqnames = Rle(snp.pos$chr), ranges = IRanges(start = snp.pos$pos, end = snp.pos$pos, width = rep(1, nrows)),
                            strand = Rle(rep('+', nrows)))
      black = countOverlaps(geno.ranges, bad.ranges, ignore.strand = T) > 0
      snp.pos = snp.pos[!black, ]
      geno.ranges = geno.ranges[!black, ]
      save(snp.pos, black, file = snp.pos.file)
    }
  }
  geno.tab = geno.tab[!black, ]
  stopifnot(all(snp.pos$chr == as.character(geno.tab$X.CHROM)), all(snp.pos$pos == geno.tab$POS), 
            all(snp.pos$id == as.character(geno.tab$ID)), all(snp.pos$ref == as.character(geno.tab$REF)))
  # Store as list, so you can store sparse matrices.
  geno.info = list()
  geno.info$alt = factor(geno.tab$ALT)
   
  genos = geno.tab[, 6]
  if(any(!grepl('^[01][\\|/][01]$', genos))){
    cat(paste('Haploid or incorrectly formatted lines for', indiv), file = stderr())
    next
  }
  geno.info$unphased = Matrix(!grepl('\\|', genos)) # Since most SNPs are phased, this should be very sparse.
  genos = gsub('/', '\\|', genos)
  tmp.list = unlist(strsplit(genos, '\\|'))
  geno.info$pat = Matrix(as.numeric(tmp.list[seq(1, length(tmp.list), 2)]) == 1)
  geno.info$mat = Matrix(as.numeric(tmp.list[seq(2, length(tmp.list), 2)]) == 1)
  geno.info$mask = Matrix(countOverlaps(geno.ranges, mask.ranges, ignore.strand = T) > 0)
  
  #geno.info = geno.info[!black, ]
  save(geno.info, file = outfile)
}
