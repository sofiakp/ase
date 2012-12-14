rm(list=ls())
require(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

indiv = 'GM18505'
filenames = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/'),
                       pattern = paste('SNYDER_HG19_', indiv, '.*\\.counts\\.RData', sep = ''), 
                       full.name = T, recursive = F, include.dirs = F)
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/plots/')
if(!file.exists(plotdir)) dir.create(plotdir)
rep.info = sample.info(filenames, '\\.counts\\.RData$')
nfiles = length(filenames)
group1 = 'mat'
group2 = 'pat'

ratio.dat = NULL

for(i in 1:nfiles){
  load(filenames[i])
  print(filenames[i])
  idx1 <- grep(group1, colnames(counts))
  idx2 <- grep(group2, colnames(counts))
  
  signif.tmp = !is.nan(snp.info$pval) # unmasked and heterozygous
  signif = signif.tmp & !grepl('chr[XY]', snp.info$chr)
  counts1 = rowSums(counts[signif, idx1])
  counts2 = rowSums(counts[signif, idx2])
  good = counts1 + counts2 > 5 
  tmp.dat = data.frame(r = counts1[good] / (counts1 + counts2)[good], mark = rep(rep.info$mark[i], sum(good)), 
                       indiv = rep(rep.info$indiv[i], sum(good)), rep = rep(rep.info$rep[i], sum(good)))
  if(is.null(ratio.dat)){
    ratio.dat = tmp.dat
  }else{
    ratio.dat = rbind(ratio.dat, tmp.dat)
  }
}

p <- ggplot(ratio.dat) + geom_density(aes(x = r, y = ..density.., color = rep)) +
  facet_wrap(~mark) + scale_x_continuous('') + labs(title = 'Ratio of maternal reads to total reads at het SNPs')
ggsave(file.path(plotdir, paste(indiv, '_count_ratio_density.png', sep = '')))
# indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/')
# outdir = file.path(indir, 'reps_6Aug12/plots')
# if(!file.exists(outdir)) file.create(outdir)
# filenames <- list.files(indir, pattern = '\\.counts\\.r', full.name = T, recursive = T)
# 
# ceu = new.env()
# load('/media/fusion10/work/chromatinVariation/rawdata/variants/trio/snps/CEU.trio.2010_09.genotypes.hg19.r', ceu)
# yri = new.env()
# load('/media/fusion10/work/chromatinVariation/rawdata/variants/trio/snps/YRI.trio.2010_09.genotypes.hg19.r', yri)
# 
# qcut = 0.01
# sb.cut = 0.01
# group1 = 'ref'
# group2 = 'alt'
# nfiles = length(filenames)
# # First col is autosomes second is chrX
# pvals = array(1, dim = c(nfiles, 2))
# ratios = array(1, dim = c(nfiles, 2))
# indivs = array('', dim = c(nfiles, 1))
# marks = array('', dim = c(nfiles, 1))
# reps = array('', dim = c(nfiles, 1))
# 
# for(i in 1:nfiles){
#   load(filenames[i])
#   reps[i] = ifelse(dim(sample)[1] > 1, 'rep', as.character(sample$rep[1]))
#   indivs[i] = as.character(sample$indiv[1])
#   marks[i] = as.character(sample$mark[1])
#   idx1 <- grep(group1, colnames(counts))
#   idx2 <- grep(group2, colnames(counts))
#   if(!is.null(ceu$geno.info[[paste(indivs[i], 'phased', sep = '.')]])){
#     phased = ceu$geno.info[[paste(indivs[i], 'phased', sep = '.')]]
#   }else{
#     phased = yri$geno.info[[paste(indivs[i], 'phased', sep = '.')]]
#   }
#   signif.tmp = phased & snp.info$pass & snp.info$sb > sb.cut & snp.info$qval < qcut
#   signif.X = signif.tmp & snp.info$chr == 'chrX'
#   signif = signif.tmp & snp.info$chr != "chrX"
#   counts1 = rowSums(counts[signif, idx1])
#   counts2 = rowSums(counts[signif, idx2])
#   ratios[i, 1] = sum(counts1 > counts2) / sum(counts1 != counts2)
#   pvals[i, 1] = binom.val(sum(counts1 > counts2), sum(counts1 != counts2))
#   counts1 = rowSums(counts[signif.X, idx1])
#   counts2 = rowSums(counts[signif.X, idx2])
#   ratios[i, 2] = sum(counts1 > counts2) / sum(counts1 != counts2)
#   pvals[i, 2] = binom.val(sum(counts1 > counts2), sum(counts1 != counts2))
# }
# 
# for(i in 1:2){
#   dat <- data.frame(indiv = indivs, mark = marks, rep = reps, pval = pvals[, i], ratio = ratios[, i])
#   is.signif = array('', dim = c(nfiles, 1))
#   is.signif[pvals[, i] < 0.001] = '*'
#   dat$signif = factor(is.signif)
#   title = ifelse(i == 1, 'autosomes', 'chrX')
#   p <- ggplot(dat) + geom_bar(aes(x = indiv, y = ratio, fill = rep), stat = "identity", position = "dodge") +
#     facet_wrap(~mark) + #scale_fill_discrete("p-value", labels = c('> 0.001', '< 0.001')) + scale_x_discrete("") + 
#     scale_y_continuous(paste("Fraction of AS events in", group1)) +
#     opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0), 
#          title = title)
#   ggsave(file.path(outdir, paste('fracInRef_reps', title, 'png', sep = '.')))
# }