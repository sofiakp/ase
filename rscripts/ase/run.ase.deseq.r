rm(list=ls())
library(DESeq)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

get.counts = function(ref, alt, tot){  
  sel = !duplicated(gsub('\\.[0-9].*', '\\.1', colnames(ref))) # replace all technical replicate identifiers with 1.1. Then only keep one for each set of technical reps.
  counts.all = cBind(ref[, sel], alt[, sel])
  size.factors = estimateSizeFactorsForMatrix(as.matrix(tot[sample(dim(tot)[1], 10000), sel]))
  size.factors = rep(size.factors, 2)
  conds = data.frame(condition = factor(append(rep('ref', sum(sel)), rep('alt', sum(sel)))), rep = factor(colnames(counts.all)))
  rownames(conds) = paste(conds$rep, conds$condition, sep = '_') 
  
  colnames(counts.all) = rownames(conds)
  return(list(counts = counts.all, conds = conds, sf = size.factors))
}

# Plots mean vs dispersion
plotDispEsts <- function(cds){
  plot(rowMeans(counts(cds, normalized = TRUE)), fitInfo(cds)$perGeneDispEsts,
       pch = '.', log = "xy", xlab = 'Normalized mean counts', ylab = 'Estimated dispersion and fit')
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

# Plots mean vs fold-change
plotDE <- function(res){
  plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.3,
       col = ifelse(res$padj < .1, "red", "black" ), xlab = 'Normalized mean counts', ylab = 'Log-fold-change' )
}

load('../../rawdata/variants/all/snps/allNonSan/GM12878.snps.RData')
het = as.vector(geno.info$mat != geno.info$pat)
mask = as.vector(geno.info$mask)

load('../../rawdata/alleleCounts/allNonSan/rdata/repsComb/SNYDER_HG19_GM12878_H3K27AC_0.RData')
count.dat = get.counts(ref, alt, tot)
sums = rowSums(count.dat$counts)
pass = het & !mask & sums > 10

cds = newCountDataSet(count.dat$counts[pass, ], count.dat$conds)
sizeFactors(cds) = count.dat$sf
cds = estimateDispersions(cds, method = "pooled-CR", modelFormula = count ~ rep + condition) # method = "blind", sharingMode="fit-only"
#       #png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
plotDispEsts(cds)
#       #dev.off()
fit = fitNbinomGLMs(cds, count ~ rep + condition )
fit.null = fitNbinomGLMs(cds, count ~ rep)

pval = array(NaN, dim = c(length(pass), 1))
pval[pass] = nbinomGLMTest(fit, fit.null)
#qval = array(NaN, dim = c(length(pass), 1))
#qval[pass] = p.adjust(pval[pass], method="BH" )
#       fit.coef = array(NaN, dim = c(length(pass), dim(fit)[2] - 2))
#       fit.coef[pass, ] = as.matrix(fit[, 1:(dim(fit)[2] - 2)])
#       fit.coef = data.frame(fit.coef, row.names = rownames(counts.all))
#       colnames(fit.coef) = colnames(fit)[1:(dim(fit)[2] - 2)]
#       vst.counts.tmp = getVarianceStabilizedData(cds)
#       vst.counts = array(NaN, dim = dim(counts))
#       vst.counts[pass, ] = vst.counts.tmp
#       vst.counts = data.frame(vst.counts, row.names = rownames(counts.all))
#       colnames(vst.counts) = colnames(counts.all)


# indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/repsComb/')
# sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*_', mark, '_0.RData', sep = ''))
# outdir = file.path(indir, 'deseq')
# plotdir = file.path(outdir, 'plots')
# if(!file.exists(outdir)) dir.create(outdir)
# if(!file.exists(plotdir)) dir.create(plotdir)
# overwrite = F
# 
# # INPUT files will be used to filter out SNPs with (even moderate) p-values in the control ChIP.
# input.files = list.files(indir, pattern = 'SNYDER_HG19_.*_INPUT_rep.RData', full.names = T)
# ninput = length(input.files)
# 
# # Peak files will be used to mark (but not filter before q-value computation) SNPs that are not within 
# # enriched regions. These might result from high levels of background noise.
# peak.files = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/peakFiles/bed/'), pattern = 'SNYDER_HG19_.*.bed', full.names = T)
# 
# # SNP positions. All count files read MUST correspond to these positions.
# snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
# load(snp.pos.file)
# snp.ranges = GRanges(seqnames = Rle(snp.pos$chr), 
#                      ranges = IRanges(start = snp.pos$pos, end = snp.pos$pos),
#                      strand = Rle(rep('+', dim(snp.pos)[1])))
# 
# # For each individual, there must be a file <geno.dir>/<indiv>.snps.RData, with genotype information.
# geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan')
# 
# q.cut = -2 # Log-10 pvalue cutoff
# overwrite = F
# males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486')
# ############# End parameter setting
# 
# for(n in 1:ninput){
#   # Load INPUT and genotype file for the individual
#   indiv = as.character(sample.info(input.files[n], '.RData')$indiv)
#   maskfile = file.path(geno.dir, paste(indiv, '.snps.RData', sep = ''))
#   
#   load(maskfile)
#   het = as.vector(geno.info$mat != geno.info$pat)
#   mask = as.vector(geno.info$mask)
#   
#   load(input.files[n])
#   stopifnot(length(mask) == dim(counts)[1])
#   inputhit = as.vector(snp.info$pval < -2)
#   
#   # CHANGE THE PATTERN HERE IF YOU WANT TO ONLY GET A SUBSET OF MARKS
#   filenames = list.files(indir, pattern = paste('SNYDER_HG19_', indiv, '.*rep.*\\.RData', sep = ''), full.name = F)
#   marks = as.character(sample.info(filenames, '.RData')$mark)
#   nfiles = length(filenames)
#   
#   if(nfiles < 1) next
#   
#   # Read each file for the individual whose input you just read.
#   for(i in 1:nfiles){
#     outfile = file.path(outdir, filenames[i])
#     hits.outfile = file.path(outdir, 'hitLists', gsub('.RData', '.hitInd.RData', basename(filenames[i])))
#     
#     if(!grepl('INPUT', filenames[i]) && (overwrite || !file.exists(outfile))){
#       print(filenames[i]) 
#       peak.file.idx = grep(paste('SNYDER_HG19_', indiv, '_', marks[i], '.*\\.bed', sep = ''), peak.files)
#       if(length(peak.file.idx) == 0){
#         cat('Missing BED file with enriched regions.\n')
#         next
#       }else{
#         if(length(peak.file.idx) > 1){
#           peak.file.idx = peak.file.idx[!grepl('reconcile.dedup_VS_', peak.files[peak.file.idx])]
#         }
#         cat('Reading peaks from ', peak.files[peak.file.idx], '\n')
#         peak.tab = read.table(peak.files[peak.file.idx], header = F)
#         peaks = GRanges(seqnames = Rle(peak.tab[,1]), 
#                         ranges = IRanges(start = peak.tab[, 2] + 1, end = peak.tab[, 3]),
#                         strand = Rle(rep('+', dim(peak.tab)[1])))  
#       }      
#       load(file.path(indir, filenames[i]))
#       stopifnot(length(mask) == dim(counts)[1])
#       
#       # Sites to use for q-value correction: heterozygous, non-masked, not "weird" (bad), without AS evidence in the INPUT
#       pass = het & !mask & !inputhit & as.vector(!snp.info$bad)
#       
#       ######## Change this
#       cds = newCountDataSet(counts[pass, ], conds)
#       sizeFactors(cds) = sf
#       cds = estimateDispersions(cds, method = "pooled") #, modelFormula = count ~ batch + condition)
#       #png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
#       plotDispEsts(cds)
#       #dev.off()
#       fit = fitNbinomGLMs(cds, count ~ batch + condition )
#       fit.null = fitNbinomGLMs(cds, count ~ batch)
#       
#       # #png(file.path(plotdir, paste(outpref, '_diff_exp', '.png', sep = '')))
#       # plotDE(res)
#       # #dev.off()
#       pval = array(NaN, dim = c(length(pass), 1))
#       pval[pass] = nbinomGLMTest(fit, fit.null)
#       qval = array(NaN, dim = c(length(pass), 1))
#       qval[pass] = p.adjust(pval[pass], method="BH" )
#       fit.coef = array(NaN, dim = c(length(pass), dim(fit)[2] - 2))
#       fit.coef[pass, ] = as.matrix(fit[, 1:(dim(fit)[2] - 2)])
#       fit.coef = data.frame(fit.coef, row.names = rownames(counts.all))
#       colnames(fit.coef) = colnames(fit)[1:(dim(fit)[2] - 2)]
#       vst.counts.tmp = getVarianceStabilizedData(cds)
#       vst.counts = array(NaN, dim = dim(counts))
#       vst.counts[pass, ] = vst.counts.tmp
#       vst.counts = data.frame(vst.counts, row.names = rownames(counts.all))
#       colnames(vst.counts) = colnames(counts.all)
#       
#       #####################
#       
#       qval = array(0, dim = c(length(pass), 1))
#       qval[pass] = log(p.adjust(10^(as.vector(snp.info$pval)[pass]), method = 'BH'), base = 10) # Convert the p-values to full for efficiency
#       
#       snp.info$pass = Matrix(pass)
#       snp.info$qval = Matrix(qval)
#       snp.info$enriched = Matrix(snp.ranges %in% peaks)
#       save(counts, snp.info, file = outfile)  
#       
#       snp.info$hits = Matrix(snp.info$qval < q.cut)
#       
#       save(snp.info, file = hits.outfile)
#     }
#   }
# }
# 
# name.expr = "rownames(regions)" #paste(as.character(regions$chr), regions$start, sep = '_')" # Expression for getting row names for the count matrix
# sel.cols = NULL # c() Set to NULL to select all columns, set to c() to select no columns ################## CHANGE!!!!!!!!!!!!!
# 
# indivs = sapply(unique(sample.info(sel.filenames, '.RData')$indiv), as.character)
# males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486')
# 
# # One DESeq run for every pair of individuals.
# for(i in 1:(length(indivs) - 1)){
#   for(j in (i + 1):length(indivs)){
#     cat(indivs[i], ' vs ', indivs[j], '\n')
#     outpref = paste(indivs[i], indivs[j], mark, 'deseq', sep = '_')
#     outfile = file.path(outdir, paste(outpref, '.RData', sep = ''))
#     if(overwrite || !file.exists(outfile)){
#       filenames = file.path(indir, sel.filenames[grep(paste('SNYDER_HG19_(', indivs[i], '|', indivs[j], ')', sep = ''), sel.filenames)])
#       c = concat.counts(filenames)
#       counts = c$counts
#       regions = c$regions
#       pass = rowSums(counts) > quantile(rowSums(counts), 0.6) # 0.4 for expression 0.6 for marks with many peaks
#       
#       if(indivs[i] %in% males || indivs[j] %in% males){
#         # Don't consider sex chromosomes, they are expected to have high p-values in comparisons between
#         # males and females and will mess up the q-value computation.
#         if(is.null(meta)){
#           pass = pass & !grepl('chr[XY]', regions$chr)
#         }else{pass = pass & !grepl('chr[XY]', meta$chr)}
#       }
#       
#       c.dat = counts[pass, ]
#       cds = newCountDataSet(c.dat, c$conds)
#       # cds = estimateSizeFactors(cds)
#       sizeFactors(cds) = c$size.factors
#       cds = estimateDispersions(cds, method = "pooled")
#       png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
#       plotDispEsts(cds)
#       dev.off()
#       res = nbinomTest(cds, indivs[i], indivs[j]) 
#       png(file.path(plotdir, paste(outpref, '_diff_exp', '.png', sep = '')))
#       plotDE(res)
#       dev.off()
#       pval = array(NaN, dim = c(length(pass), 1))
#       pval[pass] = res$pval
#       qval = array(NaN, dim = c(length(pass), 1))
#       qval[pass] = res$padj
#       fold = array(NaN, dim = c(length(pass), 1))
#       fold[pass] = res$foldChange
#       if(!is.null(sel.cols)){
#         regions = regions[, colnames(regions) %in% sel.cols]
#       }
#       regions$pval = pval
#       regions$qval = qval
#       regions$fold = fold
#       size.factors = sizeFactors(cds)
#       save(size.factors, counts, regions, file = outfile)
#     }
#   }
# }