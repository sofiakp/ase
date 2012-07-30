rm(list=ls())

binom.val <- function(x, y, p = 0.5){
  if(y > 0){
    pval = binom.test(x, y, p, alternative = 'two.sided')
    return(pval$p.value)
  }
  return(1.0)
}

usage <- function(){
  cat('Rscript simpleAse.r [options]\n', file = stderr())
  cat('Computes some basic statistics on allele counts.\n', file = stderr())
  cat('Options:\n', file =stderr())
  cat('-h\tPrint this message and exit\n', file = stderr())
  cat('-i STR\tInput file. If not given, it will read from STDIN\n', file = stderr())
  cat('-q\tCompute q-values\n', file  = stderr())
  cat('-n INT\tNumber of lines to process\n', file = stderr())
}
args <- commandArgs(trailingOnly = T)
max.nlines <- Inf
filename <- ''
get.q <- F
for(arg in args){
  if(grepl('^-n=', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2){usage(); stop(paste('Invalid option', arg))}
    max.nlines <- as.numeric(arg.split[2])
  }else if(grepl('^-i=', arg)){
    arg.split <- unlist(strsplit(arg, '='))
    if(length(arg.split) != 2){usage(); stop(paste('Invalid option', arg))}
    filename <- arg.split[2]
  }else if('-h' == arg){
    usage()
    stop('', call = F)
  }else if(arg == '-q'){
    cat('Will compute q-values\n', file = stderr())
    get.q <- T
  }else if(grepl('^-', arg)){
    stop(paste('Invalid option', arg))
  }else break
}

if(filename == ''){
  filename = 'stdin'
  if(get.q) stop('Cannot read from STDIN and compute q-values. Specify a file with -i=FILE.\n')
}
infile = file(filename, open = 'r')

ref.idx = c(1, 4) 
alt.idx = c(2, 5)
oth.idx = c(3, 6)
fwd.idx = c(1, 2, 3)
rev.idx = c(4, 5, 6)
if(get.q){
  # If the input file has more lines, this will be really slow
  all.pvals <- array(0, dim = c(5000000, 1))
  #input.lines <- array('', dim = c(5000000, 1))
}
nlines = 0

start.time <- proc.time()
while(length(line <- readLines(infile, n = 1, warn = F)) > 0){
  if(grepl('^#', line)){
    if(!get.q) write(paste(line, 'SB', 'pvalue', 'total', 'other_allele', 'allele_overlap', sep = '\t'), stdout())
    next
  }
  fields = unlist(strsplit(line, '\t'))
  nlines = nlines + 1
  if(length(fields) > 11){
    pval <- as.numeric(fields[9])
  }
  else{
    counts = array(0, dim = c(3, 6))   
    for(i in seq(1, 3)){
      counts[i, ] = as.numeric(unlist(strsplit(fields[4 + i], ',')))
    }
    flags <- ''
    tot = sum(counts[c(1,2), ])
    pval = binom.val(sum(counts[1, ]), tot) # How many come from one parent vs total?
    tot = sum(counts)
    # Overall, what fraction of reads comes from the fwd strand?
    sb <- binom.val(sum(counts[, fwd.idx]), tot)
    oth <- sum(counts[, oth.idx]) # Reads in other allele
    # Maximum overlap of alleles between maternal and paternal
    diff = sum(apply(counts[c(1,2), cbind(t(ref.idx), t(alt.idx))], 2, function(x) min(x)))
    #flags <- ifelse(sum(diff), '', 'ALLELES_OV')
    #if(flags == '') flags <- 'PASS'   
  }
  if(get.q){
    if(nlines > dim(all.pvals)[1]){
      all.pvals <- rbind(all.pvals, pval)
      #input.lines <- rbind(input.lines, line)
    }else{
      all.pvals[nlines] <- pval
      #input.lines[nlines] <- line
    }
  }else write(paste(line, sb, pval, tot, oth, diff, sep = '\t'), stdout())
  if(nlines %% 100000 == 0){
    end.time <- proc.time()
    write(paste('Read', nlines, 'lines', end.time[3] - start.time[3]), stderr())
  }
  if(nlines > max.nlines) break
}

if(get.q){
  cat('Computing q-values\n', file = stderr())
  all.pvals <- all.pvals[1:nlines]
  all.pvals <- p.adjust(all.pvals, method = 'BH')
  
  #infile = file(filename, open = 'r')
  seek(infile, where = 0, origin = "start")
  nlines = 0
  
  while(length(line <- readLines(infile, n = 1, warn = F)) > 0){
    if(grepl('^#', line)){
      write(paste(line, 'qvalue', sep = '\t'), stdout())
      next
    }
    nlines = nlines + 1
    write(paste(line, all.pvals[nlines], sep = '\t'), stdout())
    if(nlines > max.nlines) break
  }    
  #close(infile)
}
close(infile)