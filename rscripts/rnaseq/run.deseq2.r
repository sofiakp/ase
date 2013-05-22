rm(list=ls())
library(DESeq2)
source('utils/sample.info.r')

concat.counts = function(filenames, name.expr = "rownames(regions)"){
  for(f in 1:length(filenames)){
    load(filenames[f])
    if(f == 1){
      counts.all = counts
      region.names = eval(parse(text = name.expr))
      # sf = size.factors
      # conds = c(as.character(info$indiv[f]))
      # headers = c(paste(sapply(info[f, ], as.character), collapse = '_'))
    }else{
      stopifnot(all(region.names == eval(parse(text = name.expr))))
      counts.all = cbind(counts.all, counts)
      # sf = append(sf, size.factors)
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
  return(list(counts = counts.all, conds = data.frame(cond = sample.info(colnames(counts.all), '')$indiv), 
              regions = regions))
}

# Will run DESeq2 for pairs of individuals. For each individual, there should be a file
# SNYDER_HG19_indiv_mark_0.RData with one replicate per column.
# The file should contain a counts matrix with columns of the form indiv_mark_rep, a 
# regions data.frame and a size.factors array with size factors corresponding to the 
# columns of counts.
mark = 'RZ'
fit.type = 'parametric' 
beta.prior = T
indir = '../../rawdata/geneCounts/rdata/repsComb/' #'../../rawdata/signal/rep/countsAtPeaksBroad/merged_Mar13/repsComb/' #'../../rawdata/genomeGrid/hg19_w10k/rep/counts_newNorm/repsComb/' #../../rawdata/transcriptomes/rep/counts_newNorm/repsComb/'
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
meta = gene.meta # SET TO NULL UNLESS YOU'RE USING GENECOUNTS
sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*_', mark, '_0.RData', sep = ''))
outdir = file.path(indir, 'deseq2_May13')
plotdir = file.path(outdir, 'plots')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
overwrite = F

#name.expr = "rownames(regions)" #paste(as.character(regions$chr), regions$start, sep = '_')" # Expression for getting row names for the count matrix
sel.cols = c() # c() Set to NULL to select all columns, set to c() to select no columns ################## CHANGE!!!!!!!!!!!!!

indivs = sapply(unique(sample.info(sel.filenames, '.RData')$indiv), as.character)
males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486', 'GM2255', 'GM2588', 'GM2610', 'GM2630')

# One DESeq run for every pair of individuals.
for(i in 1:(length(indivs) - 1)){
  for(j in (i + 1):length(indivs)){
    cat(indivs[i], ' vs ', indivs[j], '\n')
    outpref = paste(indivs[i], indivs[j], mark, 'deseq2', sep = '_')
    outfile = file.path(outdir, paste(outpref, '.RData', sep = ''))
    if(overwrite || !file.exists(outfile)){
      filenames = file.path(indir, sel.filenames[grep(paste('SNYDER_HG19_(', indivs[i], '|', indivs[j], ')', sep = ''), sel.filenames)])
      c = concat.counts(filenames)
      counts = c$counts
      regions = c$regions
      pass = rowSums(counts) > 10 #quantile(rowSums(counts), 0.4) # 0.4 for expression 0.4 for marks with many peaks
      if(!is.null(gene.meta)) pass = rowMeans(counts) / gene.meta$len > 0.1
      
      if(indivs[i] %in% males || indivs[j] %in% males){
        # Don't consider sex chromosomes, they are expected to have high p-values in comparisons between
        # males and females and will mess up the q-value computation.
        if(is.null(meta)){
          pass = pass & !grepl('chr[XY]', regions$chr)
        }else{pass = pass & !grepl('chr[XY]', meta$chr)}
      }
      
      c.dat = counts[pass, ]
      c$conds$cond = factor(c$conds$cond, levels = c(indivs[i], indivs[j]))
      cds = DESeqSummarizedExperimentFromMatrix(countData = c.dat, colData = c$conds, design = ~cond)
      dse = DESeq(cds, fitType = fit.type, betaPrior = beta.prior, pAdjustMethod = 'BH')
      pdf(file.path(plotdir, paste(outpref, '_dispersion', '.pdf', sep = '')))
      plotDispEsts(dse, cex = 0.55, cex.lab = 1.2)
      dev.off()
      pdf(file.path(plotdir, paste(outpref, '_diff_exp', '.pdf', sep = '')))
      plotMA(dse, cex = 0.55, cex.lab = 1.2)
      dev.off()
      res = results(dse)
      pval = array(NaN, dim = c(length(pass), 1))
      pval[pass] = res$pvalue
      qval = array(NaN, dim = c(length(pass), 1))
      qval[pass] = res$FDR
      fold = array(NaN, dim = c(length(pass), 1))
      fold[pass] = res$log2FoldChange
      if(!is.null(sel.cols)){
        regions = regions[, colnames(regions) %in% sel.cols]
      }
      regions$pval = pval
      regions$qval = qval
      regions$log2Fold = fold
      size.factors = sizeFactors(dse)
      save(size.factors, counts, regions, file = outfile)
    }
  }
}
