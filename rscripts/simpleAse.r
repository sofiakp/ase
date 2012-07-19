args <- commandArgs(trailingOnly = T)
ref_idx = c(1, 4)
alt_idx = c(2, 5)
infile <- file(args[1], 'r')
line <- readLines(infile, n = 1)

count = 0
while(length(line <- readLines(infile, n = 1)) > 0){
    count = count + 1
    fields = unlist(strsplit(line, '\t'))
    counts1 = as.numeric(unlist(strsplit(fields[5], ',')))
    counts2 = as.numeric(unlist(strsplit(fields[6], ',')))
    ref_tot = sum(counts1[ref_idx] + counts2[ref_idx])
    pval = binom.test(counts1[ref_idx], ref_tot)
    lineout <- sprintf('%s', line)
    write(paste(line, pval$p.value, sep = '\t'), stdout())
}

#tab <- read.table(infile, skip = 1, sep = '\t')
#print(tab) #unlist(strsplit(tab[4], ','))