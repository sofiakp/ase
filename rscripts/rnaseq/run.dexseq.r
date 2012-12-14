rm(list=ls())
library('DEXSeq')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

plot.fit = function(cdat){
  pass = fData(cdat)$testable
  meanvalues <- rowMeans(counts(cdat)[pass, ])
  plot(meanvalues, fData(cdat)$dispBeforeSharing[pass], log="xy", main="mean vs CR dispersion")
  x <- 0.01:max(meanvalues)
  y <- cdat@dispFitCoefs[1] + cdat@dispFitCoefs[2] / x
  lines(x, y, col="red")
}

concat.counts = function(filenames, info){
  for(f in 1:length(filenames)){
    load(filenames[f])
    if(f == 1){
      counts.all = rowSums(counts)
      conds = c(as.character(info$indiv[f]))
      headers = c(paste(sapply(info[f, ], as.character), collapse = '_'))
    }else{
      stopifnot(all(rownames(counts.all) == rownames(counts)))
      counts.all = cbind(counts.all, rowSums(counts))
      conds = append(conds, as.character(info$indiv[f]))
      headers = append(headers, paste(sapply(info[f, ], as.character), collapse = '_'))
    } 
  }
  conds = factor(conds)
  return(list(counts = counts.all, conds = conds, headers = headers, gene.info = gene.info))
}

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/exonCounts/')
# List of files that will be used for differential analyses. Change the individuals and mark.
# For example to do differential tests between the CEU trio for RNA, specify the pattern
# SNYDER_HG19_GM128(78|91|92).*_RNA_.*.exoncounts.RData
# This will enumerate all files that will be used later.
mark = 'RNA'
sel.filenames = list.files(indir, pattern = paste('SNYDER_HG19_GM1289[12].*_', mark, '_.*.exoncounts.RData', sep = ''))
outdir = file.path(indir, 'dexseq')
plotdir = file.path(indir, 'dexseq', 'plots')
htmldir = file.path(indir, 'dexseq', 'html')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)
if(!file.exists(htmldir)) dir.create(htmldir)
overwrite = T

load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData'))
exon.intervals = data.frame(chr = gene.meta$chr, start = gene.meta$start, end = gene.meta$end, strand = gene.meta$strand)

indivs = sapply(unique(sample.info(sel.filenames, '.exoncounts.RData')$indiv), as.character)
males = c('SNYDER', 'GM12891', 'GM19239', 'GM18486')

for(i in 1:(length(indivs) - 1)){
  print(indivs[i])
  for(j in (i + 1):length(indivs)){
    filenames = file.path(indir, sel.filenames[grep(paste('SNYDER_HG19_(', indivs[i], '|', indivs[j], ')', sep = ''), sel.filenames)])
    #list.files(indir, pattern = paste('SNYDER_HG19_(', indivs[i], '|', indivs[j], ')_RNA_.*counts.RData', sep = ''), full.names = T)
    outpref = paste(indivs[i], indivs[j], mark, 'dexseq', sep = '_')
    outfile = file.path(outdir, paste(outpref, '.RData', sep = ''))
    if(overwrite || !file.exists(outfile)){
      samples = sample.info(filenames, '.exoncounts.RData')
      c = concat.counts(filenames, samples)
      counts = c$counts
      gene.info = c$gene.info
      pass = rowSums(counts) > quantile(rowSums(counts), 0.4)
      
      if(indivs[i] %in% males || indivs[j] %in% males){
        # Don't consider sex chromosomes, they are expected to have high p-values in comparisons between
        # males and females and will mess up the q-value computation.
        pass = pass & !grepl('chr[XY]', gene.meta$chr)
      }
      
      c.dat = data.frame(counts[pass, ])
      colnames(c.dat) = c$headers
      exon.set = newExonCountSet(c.dat, c$conds, gene.meta$gene[pass], gene.meta$exon[pass], exon.intervals[pass, ], 
                                 as.character(gene.meta$transcripts[pass]))
      exon.set = estimateSizeFactors(exon.set)
      exon.set = estimateDispersions(exon.set, nCores = 10, minCount = 10, maxExon = 50)
      exon.set = fitDispersionFunction(exon.set)
      png(file.path(plotdir, paste(outpref, '_dispersion', '.png', sep = '')))
      plot.fit(exon.set)
      dev.off()
      
      exon.set = testForDEU(exon.set, nCores = 10)
      exon.set = estimatelog2FoldChanges(exon.set, nCores=10)
      res = DEUresultTable(exon.set) # this has one row per exon
      
      testable = array(F, dim = c(length(pass), 1))
      testable[pass] = fData(exon.set)$testable
      pass = pass & testable
      pval = array(NaN, dim = c(length(pass), 1))
      pval[pass] = res$pvalue[fData(exon.set)$testable]
      qval = array(NaN, dim = c(length(pass), 1))
      qval[pass] = res$padjust[fData(exon.set)$testable]
      fold = array(NaN, dim = c(length(pass), 1))
      fold[pass] = res[fData(exon.set)$testable, 7]
      gene.info$pval = pval
      gene.info$qval = qval
      gene.info$fold = fold
      size.factors = sizeFactors(exon.set)
      save(size.factors, counts, samples, gene.info, exon.set, file = outfile)
      DEXSeqHTML(exon.set, FDR = 0.01, path = file.path(htmldir, outpref), file = paste(outpref, '_results.html', sep = ''), color=c("#FF000080", "#0000FF80"))
    }
  }
}