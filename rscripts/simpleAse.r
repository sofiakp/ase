rm(list=ls())

binom.val <- function(x, y, p = 0.5){
  if(y > 0){
    pval = binom.test(x, y, p, alternative = 'two.sided')
    return(pval$p.value)
  }
  return(1.0)
}

args <- commandArgs(trailingOnly = T)
max.nlines <- Inf
infile <- ''
for(arg in args){
  if(grepl('^-n', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2) stop(paste('Invalid option', arg))
    max.nlines <- int(arg.split[2])
  }else if(grepl('^-', arg)){
    stop(paste('Invalid option', arg))
  }else{
    infile <- arg
  }
}

if(infile != ''){
  infile <- file(infile, open = 'r')
}else{
  infile = file('stdin', open = 'r')
}

ref.idx = c(1, 4) 
alt.idx = c(2, 5)
oth.idx = c(3, 6)
fwd.idx = c(1, 2, 3)
rev.idx = c(4, 5, 6)
# If the input file has more lines, this will be really slow
# all_pvals <- array(0, dim = c(5000000, 4))
nlines = 0

start.time <- proc.time()
while(length(line <- readLines(infile, n = 1, warn = F)) > 0){
  if(grepl('^#', line)){
    write(paste(line, 'SB', 'pvalue', 'total', 'other_allele', 'allele_overlap', sep = '\t'), stdout())
    next
  }
  fields = unlist(strsplit(line, '\t'))
  nlines = nlines + 1
  counts = array(0, dim = c(3, 6))   
  for(i in seq(1, 3)){
    counts[i, ] = as.numeric(unlist(strsplit(fields[4 + i], ',')))
  }
  flags <- ''
  # Overall, what fraction of reads comes from the fwd strand?
  sb <- binom.val(sum(counts[, fwd.idx]), sum(counts))
  tot = sum(counts[c(1,2), ])
  pval = binom.val(sum(counts[1, ]), tot) # How many come from one parent vs total?
  oth <- ifelse(tot > 0, sum(counts[, oth.idx]) / tot, 0) # Fraction in other
  # Maximum overlap of alleles between maternal and paternal
  diff = sum(apply(counts[c(1,2), ], 2, function(x) min(x)))
  #flags <- ifelse(sum(diff), '', 'ALLELES_OV')
  #if(flags == '') flags <- 'PASS'
  
  write(paste(line, sb, pval, tot, oth, diff, sep = '\t'), stdout())
  if(nlines %% 100000 == 0){
    end.time <- proc.time()
    write(paste('Read', nlines, 'lines', end.time[3] - start.time[3]), stderr())
  }
  if(nlines > max.nlines) break
}
close(infile)

  #if(nlines > dim(all.pvals)[1]){
  #  all.pvals <- rbind(all.pvals, t(c(sb, pval, pass, 1)))
  #}else{
  #  all.pvals[nlines, ] <- t(c(sb, pval, pass, 1))
  #}

#all_pvals <- all_pvals[seq(1,nlines), ]
#all_pvals[all_pvals[, 3] == 1, 4] <- p.adjust(all_pvals[all_pvals[, 3] == 1, 2], method = 'BH')

#infile <- file(args[1], 'r')
#line <- readLines(infile, n = 1)
#nlines = 0
#while(length(line <- readLines(infile, n = 1)) > 0){
#    nlines = nlines + 1
#    lineout <- sprintf('%s', line)	
#    write(paste(line, paste(all_pvals[nlines,], collapse = '\t'), sep = '\t'), stdout())
#    #if(nlines > 100000) break
#}    
#close(infile)
#tab <- read.table(infile, skip = 1, sep = '\t')
#print(tab) #unlist(strsplit(tab[4], ','))
