rm(list = ls())
library('DEXSeq')
library('multicore')

plot.fit = function(cdat){
  pass = fData(cdat)$testable
  meanvalues <- rowMeans(counts(cdat)[pass, ])
  plot(meanvalues, fData(cdat)$dispBeforeSharing[pass], log="xy", main="mean vs CR dispersion")
  x <- 0.01:max(meanvalues)
  y <- cdat@dispFitCoefs[1] + cdat@dispFitCoefs[2] / x
  lines(x, y, col="red")
}
infiles = c('rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_1_maternal.exoncounts',
            'rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_2_maternal.exoncounts',
            'rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_3_maternal.exoncounts',
            'rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_1_paternal.exoncounts',
            'rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_2_paternal.exoncounts',
            'rawdata/exonCounts/SNYDER_HG19_GM12878_RNA_3_paternal.exoncounts')
exon.dat = read.HTSeqCounts(infiles, data.frame(condition = factor(c('mat','mat','mat','pat','pat','pat'))),
                            flattenedfile = 'rawdata/transcriptomes/gencode.v13.annotation.noM.flat.gff')

exon.dat = estimateSizeFactors(exon.dat)
exon.dat = estimateDispersions(exon.dat, nCores = 10, minCount = 10, maxExon = 20)
exon.dat = fitDispersionFunction(exon.dat)
plot.fit(exon.dat)

exon.dat = testForDEU(exon.dat, nCores = 10)
exon.dat = estimatelog2FoldChanges(exon.dat, nCores=10)
res = DEUresultTable(exon.dat) # this has one row per exon