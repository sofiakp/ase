rm(list=ls())
library(DESeq2)
source('utils/sample.info.r')
source('utils/deseq.utils.r')

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

marks = c('H3K4ME1') #, 'H3K4ME3', 'CTCF', 'SA1')
fit.type = 'parametric' 
beta.prior = T
indir = '../../rawdata/signal/rep/countsAtPeaksBroad/merged_Mar13/repsComb/' #'../../rawdata/geneCounts/rdata/repsComb/' #'../../rawdata/genomeGrid/hg19_w10k/rep/counts_newNorm/repsComb/' #../../rawdata/transcriptomes/rep/counts_newNorm/repsComb/'
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
meta = NULL # SET TO NULL UNLESS YOU'RE USING GENECOUNTS
outdir = file.path(indir, 'deseq_pop_May13')
plotdir = file.path(outdir, 'plots')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
overwrite = F

for(m in marks){
  sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*_', m, '_0.RData', sep = ''), full.names = T)
  indivs = sapply(unique(sample.info(sel.filenames, '.RData')$indiv), as.character)
  
  c = concat.counts(sel.filenames)
  counts = c$counts
  regions = c$regions
  pass = rowSums(counts) > quantile(rowSums(counts), 0.4) # for expression 0.4 for marks with many peaks
  if(!is.null(meta)) pass = rowMeans(counts) / gene.meta$len > 0.1
  if(is.null(meta)){
    pass = pass & !grepl('chr[XY]', regions$chr)
  }else{pass = pass & !grepl('chr[XY]', meta$chr)}
  
  regions = regions[pass, ]
  c.dat = counts[pass, ]
  conds = c$conds
  conds$cond = factor(conds$cond, levels = indivs)
  conds$pop = factor(get.pop(as.character(conds$cond)))
  cds = DESeqSummarizedExperimentFromMatrix(countData = c.dat, colData = conds, design = ~ pop)
  dse = DESeq(cds, fitType = fit.type, betaPrior = beta.prior, pAdjustMethod = 'BH')
  
  cds = newCountDataSet(c.dat, conds)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  fit1 = fitNbinomGLMs(cds, count ~ pop + cond)
  fit0 = fitNbinomGLMs(cds, count ~ pop)
  pvalsGLM = nbinomGLMTest(fit1, fit0)
  qvalsGLM = -log10(p.adjust(pvalsGLM, method = 'BH'))
  
  regions$pval = pvalsGLM
  regions$qval = qvalsGLM
  outfile = file.path(outdir, paste('pop_deseq_', m, '.RData', sep = ''))
  save(regions, file = outfile)
}
q.cuts = seq(10, 0, -1)
q.num = array(0, dim = c(length(q.cuts), 1))
for(i in 1:length(q.num)){
  q.num[i] = sum(qvalsGLM > q.cuts[i], na.rm = T) / nrow(c.dat)
}
q.dat = data.frame(n = q.num, c = q.cuts)
p = ggplot(q.dat, aes(x = c, y = n)) + geom_point() + geom_line() + theme_bw() + 
  xlab('q-value cutoff') + ylab('Fraction of regions')
