rm(list=ls())
library(DESeq)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

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

concat.counts = function(filenames){
  for(f in 1:length(filenames)){
    load(filenames[f])
    if(f == 1){
      counts.all = counts
      region.names = eval(parse(text = name.expr))
      sf = size.factors
      # conds = c(as.character(info$indiv[f]))
      # headers = c(paste(sapply(info[f, ], as.character), collapse = '_'))
    }else{
      stopifnot(all(region.names == eval(parse(text = name.expr))))
      counts.all = cbind(counts.all, counts)
      sf = append(sf, size.factors)
      # conds = append(conds, as.character(info$indiv[f]))
      # headers = append(headers, paste(sapply(info[f, ], as.character), collapse = '_'))
    } 
  }
  rownames(counts.all) = region.names
  # Remove technical replicates for correct dispersion estimation.
  sel = !duplicated(gsub('\\.[0-9].*', '\\.1', sample.info(colnames(counts.all), ''))) # replace all technical replicate identifiers with 1.1. Then only keep one for each set of technical reps.
  counts.all = counts.all[, sel]
  #counts = data.frame(counts.all)
  #colnames(counts) = headers
  return(list(counts = counts.all, conds = sample.info(colnames(counts.all), '')$indiv, regions = regions, size.factors = sf[sel]))
}

# Will run DESeq for pairs of individuals. For each individual, there should be a file
# SNYDER_HG19_indiv_mark_0.RData with one replicate per column.
# The file should contain a counts matrix with columns of the form indiv_mark_rep, a 
# regions data.frame and a size.factors array with size factors corresponding to the 
# columns of counts.
mark = 'H3K4ME3'
indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/countsAtPeaks/repsComb/')
load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData'))
meta = NULL # gene.meta ############### CHANGE!!!!!!!!!!!!
sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*_', mark, '_0.RData', sep = ''))
outdir = file.path(indir, 'deseq_v3')
plotdir = file.path(outdir, 'plots')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
overwrite = F

name.expr = "rownames(regions)" #paste(as.character(regions$chr), regions$start, sep = '_')" # Expression for getting row names for the count matrix
sel.cols = NULL # c() Set to NULL to select all columns, set to c() to select no columns ################## CHANGE!!!!!!!!!!!!!

uniq.indivs = c('GM10847', 'GM18505', 'GM18526', 'GM18951', 'GM19099', 'GM12891', 'GM12892', 'GM12878',
                'GM19238', 'GM19239', 'GM19240', 'GM18486', 'GM12890', 'SNYDER')
batch = factor(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3))

counts.all = NULL

for(i in 1:length(uniq.indivs)){
  cat(uniq.indivs[i], '\n')
  filenames = file.path(indir, sel.filenames[grep(paste('SNYDER_HG19_', uniq.indivs[i], sep = ''), sel.filenames)])
  for(j in 1:length(filenames)){
    load(filenames[j])
    cond.tmp = data.frame(condition = rep(uniq.indivs[i], dim(counts)[2]), batch = rep(batch[i], dim(counts)[2]),
                          row.names = colnames(counts))
    if(is.null(counts.all)){
      counts.all = counts
      region.names = eval(parse(text = name.expr))
      sf = size.factors
      conds = cond.tmp   
    }else{
      stopifnot(all(region.names == eval(parse(text = name.expr))))
      counts.all = cbind(counts.all, counts)
      sf = append(sf, size.factors)
      conds = rbind(conds, cond.tmp)
    } 
  }
}

# rownames(counts.all) = region.names
# Remove technical replicates for correct dispersion estimation.
sel = !duplicated(gsub('\\.[0-9].*', '\\.1', sample.info(colnames(counts.all), ''))) # replace all technical replicate identifiers with 1.1. Then only keep one for each set of technical reps.
counts = counts.all[, sel]
conds = conds[sel, ]
pass = rowSums(counts) > quantile(rowSums(counts), 0.6)
if(is.null(meta)){
  pass = pass & !grepl('chr[XY]', regions$chr)
}else{pass = pass & !grepl('chr[XY]', meta$chr)}


cds = newCountDataSet(counts[pass, ], conds)
sizeFactors(cds) = sf
cds = estimateDispersions(cds, method = "pooled") #, modelFormula = count ~ batch + condition)
#png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
plotDispEsts(cds)
#dev.off()
fit = fitNbinomGLMs(cds, count ~ batch + condition )
fit.null = fitNbinomGLMs(cds, count ~ batch)

# #png(file.path(plotdir, paste(outpref, '_diff_exp', '.png', sep = '')))
# plotDE(res)
# #dev.off()
pval = array(NaN, dim = c(length(pass), 1))
pval[pass] = nbinomGLMTest(fit, fit.null)
qval = array(NaN, dim = c(length(pass), 1))
qval[pass] = p.adjust(pval[pass], method="BH" )
fit.coef = array(NaN, dim = c(length(pass), dim(fit)[2] - 2))
fit.coef[pass, ] = as.matrix(fit[, 1:(dim(fit)[2] - 2)])
fit.coef = data.frame(fit.coef, row.names = rownames(counts.all))
colnames(fit.coef) = colnames(fit)[1:(dim(fit)[2] - 2)]
vst.counts.tmp = getVarianceStabilizedData(cds)
vst.counts = array(NaN, dim = dim(counts))
vst.counts[pass, ] = vst.counts.tmp
vst.counts = data.frame(vst.counts, row.names = rownames(counts.all))
colnames(vst.counts) = colnames(counts.all)
# fold = array(NaN, dim = c(length(pass), 1))
# fold[pass] = res$foldChange
# if(!is.null(sel.cols)){
#   regions = regions[, colnames(regions) %in% sel.cols]
# }
regions$pval = pval
regions$qval = qval
# regions$fold = fold
size.factors = sizeFactors(cds)
outfile = file.path(outdir, 'SNYDER_HG19_H3K4ME3_deseq.RData')
save(size.factors, counts, vst.counts, regions, fit.coef, file = outfile)
