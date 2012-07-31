rm(list=ls())
datafiles <- list.files('/media/fusion10/work/chromatinVariation/src/rscripts/utils/', pattern = '\\.[rR]', full.name = T)
for(f in seq(length(datafiles))) source(datafiles[f])

usage <- function(){
  cat('Rscript createAseData.r [options]\n', file = stderr())
  cat('Reads allele counts, computes p-values, and writes to R data file.\n', file = stderr())
  cat('Options:\n', file =stderr())
  cat('-h\tPrint this message and exit\n', file = stderr())
  cat('-i=FILE\tInput file.\n', file = stderr())
  cat('-m=FILE\tMask file. SNPs in this file will not be included in q-value computation.\n', file = stderr())
  cat('-o=FILE\tOutput file.\n', file = stderr())
  cat('-q\tCompute q-values\n', file  = stderr())
}

args <- commandArgs(trailingOnly = T)
infile <- ''
maskfile <- ''
outfile <- ''
get.q <- F
for(arg in args){
  if(grepl('^-m=', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2){usage(); stop(paste('Invalid option', arg))}
    maskfile <- arg.split[2]
  }else if(grepl('^-i=', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2){usage(); stop(paste('Invalid option', arg))}
    infile <- arg.split[2]
  }else if(grepl('^-o=', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2){usage(); stop(paste('Invalid option', arg))}
    outfile <- arg.split[2]
  }else if('-h' == arg){
    usage()
    stop('', call = F)
  }else if(arg == '-q'){
    get.q <- T
  }else if(grepl('^-', arg)){
    stop(paste('Invalid option', arg))
  }else break
}

if(infile == '') stop('Missing input file')
if(maskfile == '') stop('Missing mask file')
if(outfile == '') stop('Missing output file')

# Read filtered regions
cat('Reading mask file\n', file = stderr())
pass.tab <- read.table(maskfile, header = F, sep = '\t', comment.char = '#')
name.map <- new.env(hash = T)
names <- paste(pass.tab[,1], pass.tab[,2], sep = ':')
for(i in seq(length(names))){
  name.map[[names[i]]] <- T
}

# Get sample info
sample <- sample.info(infile)

cat('Reading allele counts\n', file = stderr())
tab <- read.table(infile, header = T, sep = '\t', comment.char = '')
nrows <- dim(tab)[1]
mat.idx = which(colnames(tab) == 'maternal') - 4
pat.idx = which(colnames(tab) == 'paternal') - 4
amb.idx = which(colnames(tab) == 'ambiguous') - 4

# Read counts
# 1st dimension is SNP, 2nd is read group, 3rd is fwd vs rev strand
ref = array(0, dim = c(nrows, 3, 2))
alt = array(0, dim = c(nrows, 3, 2))
oth = array(0, dim = c(nrows, 3, 2))

cat('Building allele count matrices\n', file = stderr())
for(i in 1:3){
  fields <- unlist(strsplit(as.character(tab[, 4 + i]), ','))
  for(j in 0:1){
    ref[, i, j + 1] <- as.numeric(fields[seq(3 * j + 1, length(fields), 6)])
    alt[, i, j + 1] <- as.numeric(fields[seq(3 * j + 2, length(fields), 6)])
    oth[, i, j + 1] <- as.numeric(fields[seq(3 * j + 3, length(fields), 6)])
  }
}

tot <- rowSums(ref + alt + oth) # Total read count
fwd.tot <- rowSums(ref[,,1] + alt[,,1] + oth[,,1])
cat('Computing p-values\n', file = stderr())
sb <- binom.val(fwd.tot, tot)
tot.unamb <- rowSums(ref[,c(mat.idx, pat.idx),] + alt[,c(mat.idx, pat.idx),] + oth[,c(mat.idx, pat.idx),]) # Total reads in mat or pat
pval <- binom.val(rowSums(ref[,mat.idx,] + alt[,mat.idx,] + oth[,mat.idx,]), tot.unamb)

# Indicators of which SNPs are in the filter file
pass <- array(F, dim = c(nrows, 1))
for(i in 1:nrows) pass[i] <- !is.null(name.map[[paste(tab$X.CHROM[i], tab$POS[i], sep = ':')]])

qval <- array(1, dim = dim(pass))
if(get.q){
  cat('Computing q-values\n', file = stderr())
  qval[pass] <- p.adjust(pval[pass], method = 'BH')
}

headers <- array('', dim = c(1, 6))
headers[mat.idx + c(0, 3)] <- 'mat'  
headers[pat.idx + c(0, 3)] <- 'pat'  
headers[amb.idx + c(0, 3)] <- 'amb'  
headers <- paste(headers, cbind(rep('fwd',3), rep('rev',3)), sep = '.')
counts <- data.frame(cbind(ref[,,1], ref[,,2], alt[,,1], alt[,,2], oth[,,1], oth[,,2]))
colnames(counts) <- append(append(paste(headers, 'ref', sep = '.'), paste(headers, 'alt', sep = '.')), paste(headers, 'oth', sep = '.'))

snp.info <- data.frame(chr = tab$X.CHROM, pos = tab$POS, ref = tab$REF, alt = tab$ALT, pass = pass, sb = sb, pval = pval, qval = qval)

# Save counts, snp.info, sample
save(counts, snp.info, sample, file = outfile)