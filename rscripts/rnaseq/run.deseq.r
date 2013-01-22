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

concat.counts = function(filenames, name.expr = "rownames(regions)"){
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
  sel = !duplicated(gsub('\\.[0-9].*', '', colnames(counts.all))) # replace all technical replicate identifiers with 1.1. Then only keep one for each set of technical reps.
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
mark = 'POL4H8'
method = 'blind' # Set to pooled unless you don't have replicates
if(method == 'blind'){
  sharingMode = 'fit-only'
}else{
  sharingMode = 'maximum'
}
indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/rep/counts/repsComb')# 'rawdata/genomeGrid/hg19_w10k/rep/counts/repsComb') #rawdata/signal/rep/countsAtPeaksBroad/repsComb/')
#load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData'))
meta = NULL # gene.meta ############### CHANGE!!!!!!!!!!!!
sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*_', mark, '_0.RData', sep = ''))
outdir = file.path(indir, 'deseq')
plotdir = file.path(outdir, 'plots')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
overwrite = T

#name.expr = "rownames(regions)" #paste(as.character(regions$chr), regions$start, sep = '_')" # Expression for getting row names for the count matrix
sel.cols = NULL # c() Set to NULL to select all columns, set to c() to select no columns ################## CHANGE!!!!!!!!!!!!!

indivs = sapply(unique(sample.info(sel.filenames, '.RData')$indiv), as.character)
males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486')

# One DESeq run for every pair of individuals.
for(i in 1:(length(indivs) - 1)){
  for(j in (i + 1):length(indivs)){
    cat(indivs[i], ' vs ', indivs[j], '\n')
    outpref = paste(indivs[i], indivs[j], mark, 'deseq', sep = '_')
    outfile = file.path(outdir, paste(outpref, '.RData', sep = ''))
    if(overwrite || !file.exists(outfile)){
      filenames = file.path(indir, sel.filenames[grep(paste('SNYDER_HG19_(', indivs[i], '|', indivs[j], ')', sep = ''), sel.filenames)])
      c = concat.counts(filenames)
      counts = c$counts
      regions = c$regions
      pass = rowSums(counts) > quantile(rowSums(counts), 0.4) # 0.4 for expression 0.4 for marks with many peaks
      
      if(indivs[i] %in% males || indivs[j] %in% males){
        # Don't consider sex chromosomes, they are expected to have high p-values in comparisons between
        # males and females and will mess up the q-value computation.
        if(is.null(meta)){
          pass = pass & !grepl('chr[XY]', regions$chr)
        }else{pass = pass & !grepl('chr[XY]', meta$chr)}
      }
      
      c.dat = counts[pass, ]
      cds = newCountDataSet(c.dat, c$conds)
      # cds = estimateSizeFactors(cds)
      sizeFactors(cds) = c$size.factors
      res = try({cds = estimateDispersions(cds, method = method, sharingMode = sharingMode)}, silent = TRUE)
      if(class(res) == 'try-error'){
        cat('Using local fit\n', file = stderr())
        cds = estimateDispersions(cds, method = method, sharingMode = sharingMode, fitType = 'local')
      }#else{
      #  cds = estimateDispersions(cds, method = "pooled")
      #}
      png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
      plotDispEsts(cds)
      dev.off()
      res = nbinomTest(cds, indivs[i], indivs[j]) 
      png(file.path(plotdir, paste(outpref, '_diff_exp', '.png', sep = '')))
      plotDE(res)
      dev.off()
      pval = array(NaN, dim = c(length(pass), 1))
      pval[pass] = res$pval
      qval = array(NaN, dim = c(length(pass), 1))
      qval[pass] = res$padj
      fold = array(NaN, dim = c(length(pass), 1))
      fold[pass] = res$foldChange
      if(!is.null(sel.cols)){
        regions = regions[, colnames(regions) %in% sel.cols]
      }
      regions$pval = pval
      regions$qval = qval
      regions$fold = fold
      size.factors = sizeFactors(cds)
      save(size.factors, counts, regions, file = outfile)
    }
  }
}
