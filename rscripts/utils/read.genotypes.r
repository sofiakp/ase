rm(list=ls())
library('GenomicRanges')

# Input directory should have VCFs WITHOUT metadata rows (##), but with the main header.
snp.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps')
# Bed files with regions to mask. Should be named <indiv>.blacklist.bed.
mask.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/masks')
overwrite = F

snp.files = list.files(snp.dir, recursive = F, full.names = F, pattern = '.snps.vcf')
for(f in snp.files[14]){
  geno.tab <- read.table(file.path(snp.dir, f), header = T, sep = '\t', comment.char = '')
  #colnames(geno.tab) = gsub('^NA', 'GM', colnames(geno.tab))
  main.geno.info <- data.frame(chr = geno.tab$X.CHROM, pos = geno.tab$POS, id = geno.tab$ID, ref = geno.tab$REF, alt = geno.tab$ALT)
  nrows = dim(main.geno.info)[1]
  
  # Make a GRanges object for easy overlapping
  geno.ranges = GRanges(seqnames = Rle(main.geno.info$chr), ranges = IRanges(start = main.geno.info$pos, end = main.geno.info$pos, width = rep(1, nrows)),
                        strand = Rle(rep('+', nrows)))  
  for(i in 10:dim(geno.tab)[2]){
    indiv = colnames(geno.tab)[i]
    mask.file = file.path(mask.dir, paste(indiv, '.blacklist.bed', sep = ''))
    outfile = file.path(snp.dir, paste(indiv, '.snps.r', sep = ''))
    if(file.exists(outfile) && !overwrite){
      cat(paste('File', outfile, 'exists. Skipping sample.'), file = stderr())
      next
    }
    if(!file.exists(mask.file)){
      cat(paste('Missing file:', mask.file), file = stderr())
      next
    }
    cat(paste('Reading data for', indiv, '\n'), file = stderr())
    
    geno.info = main.geno.info
    #tmp.list = unlist(strsplit(as.character(geno.tab[, i]), ':'))
    genos = geno.tab[, i] #tmp.list[seq(1, length(tmp.list), 3)]
    if(any(!grepl('^[01][\\|]|[01]$', geno.tab[, i]))){
      cat(paste('Haploid lines for', indiv), file = stderr())
      next
    }
    geno.info[[paste(indiv, 'phased', sep = '.')]] = grepl('\\|', genos)
    genos = gsub('/', '\\|', genos)
    tmp.list <- unlist(strsplit(genos, '\\|'))
    geno.info[[paste(indiv, 'pat', sep = '.')]] = as.numeric(tmp.list[seq(1, length(tmp.list), 2)]) == 1
    geno.info[[paste(indiv, 'mat', sep = '.')]] = as.numeric(tmp.list[seq(2, length(tmp.list), 2)]) == 1
    
    mask = read.table(mask.file, header = F, sep = '\t', comment.char = '#')
    mask.ranges = GRanges(seqnames = Rle(mask[,1]), 
                          ranges = IRanges(start = mask[, 2] + 1, end = mask[, 3]),
                          strand = Rle(rep('+', dim(mask)[1])))
    geno.info[[paste(indiv, 'mask', sep = '.')]] = geno.ranges %in% mask.ranges
    save(geno.info, file = outfile)
  }
}