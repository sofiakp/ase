rm(list=ls())

binom.val <- function(x, y, p = 0.5){
    if(y > 0){
        pval = binom.test(x, y, p)
	return(pval$p.value)
    }
    return(1.0)
}
args <- commandArgs(trailingOnly = T)
ref_idx = c(1, 4) 
alt_idx = c(2, 5)
oth_idx = c(3, 6)
fwd_idx = c(1, 2, 3)
rev_idx = c(4, 5, 6)

infile <- file(args[1], open = 'r')
line <- readLines(infile, n = 1)
all_pvals <- array(0, dim = c(5000000, 4))
nlines = 0

while(length(line <- readLines(infile, n = 1, warn = F)) > 0){
    start_time <- proc.time()
    fields = unlist(strsplit(line, '\t'))
    nlines = nlines + 1
    counts = array(0, dim = c(3, 6))   
    for(i in seq(1, 3)){
        counts[i, ] = as.numeric(unlist(strsplit(fields[4 + i], ',')))
    }
    # Overall, what fraction of reads comes from the fwd strand?
    sb <- binom.val(sum(counts[, fwd_idx]), sum(counts))
    tot = sum(counts[c(1,2), ])
    pval = binom.val(sum(counts[1, ]), tot)
    # This should be an input argument
    # Mark as non passing if any of the parents has evidence for more than one of REF and ALT
    if(sum(counts[, oth_idx]) == 0){
    	# && xor(sum(counts[1, ref_idx]), sum(counts[1, alt_idx])) && xor(sum(counts[2, ref_idx]), sum(counts[2, alt_idx]))){
        pass = 1	
    }else{pass = 0}
    if(nlines > dim(all_pvals)[1]){
        all_pvals <- rbind(all_pvals, t(c(sb, pval, pass, 1)))
    }else{
	all_pvals[nlines, ] <- t(c(sb, pval, pass, 1))
    }
    end_time <- proc.time()
    write(paste(line, paste(all_pvals[nlines,], collapse = '\t'), end_time[3] - start_time[3], sep = '\t'), stdout())
    if(nlines %% 100000 == 0) write(paste('Read', nlines, 'lines'), stderr())
    #if(nlines > 100000) break
}
close(infile)

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