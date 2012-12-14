rm(list=ls())
library('GenomicRanges')
library('reshape')
require('ggplot2')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

# Find overlaps between AS SNPs and AS genes and plot the fraction of AS SNPs that 
# have the same bias (i.e. in the same direction) as the gene they overlap.
is.gene = T # genecounts vs exoncounts
if(is.gene){
  subdir = 'geneCounts'
  suf = 'genecounts'
  load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
}else{
  subdir = 'exonCounts'
  suf = 'exoncounts'  
  load('rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')
}
outsuf = 'RNA'
rna.files = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata', subdir, 'rdata/reps/qvals/'),
                       pattern = paste('SNYDER_HG19_.*', outsuf, '_rep.', suf, '.RData', sep  = ''), full.name = T)
hitdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals/hitLists/')
outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', subdir, 'rdata/reps', 'plots')
if(!file.exists(outdir)) dir.create(outdir)

rna.samples = sample.info(rna.files, paste('.', suf, '.RData', sep = ''))
agree.dat = NULL

for(n in 1:length(rna.files)){
  load(rna.files[n]) 
  print(basename(rna.files[n]))
  pass = !is.nan(gene.info$qval) & gene.info$qval < 0.01
  gene.ranges = GRanges(seqnames = Rle(gene.meta$chr), 
                          ranges = IRanges(start = gene.meta$start, end = gene.meta$end),
                          strand = Rle(rep('+', dim(gene.meta)[1])))
  gene.ranges = gene.ranges[pass, ]
  counts = counts[pass, ]
  gene.ratios = log((counts$mat + 1) / (counts$pat + 1), base = 2)
  
  # Get all lists of AS SNPs for the individual
  filenames = list.files(hitdir, pattern = paste('SNYDER_HG19_', rna.samples$indiv[n], '.*\\.hits\\.RData', sep = ''), 
                         full.name = T, recursive = F, include.dirs = F)
  
  nfiles = length(filenames)
  for(i in 1:nfiles){
    load(filenames[i])
    mark = sample.info(filenames[i], '.hits.RData')$mark
    #hits = het.counts
    hits = hits[hits$phased, ] # We want phased SNPs since we're looking at mat/pat ratios 
    snp.ranges = GRanges(seqnames = Rle(hits$chr), 
                         ranges = IRanges(start = hits$pos, end = hits$pos),
                         strand = Rle(rep('+', dim(hits)[1])))
    snp.ratios = (hits$ref.counts + 1) / (hits$alt.counts + 1)
    snp.ratios[hits$mat] = 1 / snp.ratios[hits$mat]
    snp.ratios = log(snp.ratios, base = 2)
    
    # 2 columns: 
    # queryHits: index in gene.ranges
    # subjectHits: index of overlapping SNP in snp.ranges
    ov.dat = findOverlaps(gene.ranges, snp.ranges)
    avg.ratios = array(NaN, dim = c(length(gene.ratios), 1)) 
    min.dist = array(-100000000, dim = c(length(gene.ratios), 1))
    
    for(g in 1:length(gene.ratios)){
      # Distance between the gene start and the closest phased AS SNP
      d = start(snp.ranges)[as.character(seqnames(snp.ranges)) == as.character(seqnames(gene.ranges)[g])] - start(gene.ranges)[g]
      if(length(d) > 0){
        min.dist[g] = sign(d[which.min(abs(d))]) * min(abs(d))
        sel.ratios = snp.ratios[subjectHits(ov.dat)[queryHits(ov.dat) == g]]
      }
      # If there are overlapping SNPs and the SNPs agree in direction...
      if(length(sel.ratios) > 0 && (all(sel.ratios >= 0)|| all(sel.ratios <= 0))) avg.ratios[g] = mean(sel.ratios)
    }
    agree = sum(gene.ratios[!is.nan(avg.ratios)] * avg.ratios[!is.nan(avg.ratios)] > 0)
    disagree = sum(gene.ratios[!is.nan(avg.ratios)] * avg.ratios[!is.nan(avg.ratios)] < 0)
    tmp.dat = data.frame(indiv = rna.samples$indiv[n], mark = mark,
                         agree = agree, dis = disagree)
    tmp.dist = data.frame(indiv = rep(rna.samples$indiv[n], length(min.dist)), mark = rep(mark, length(min.dist)), dist = min.dist)
    if(is.null(agree.dat)){
      agree.dat = tmp.dat
      snp.dist = tmp.dist
    }else{
      agree.dat = rbind(agree.dat, tmp.dat)
      snp.dist = rbind(snp.dist, tmp.dist)
    }
  }  
}

agree.dat = melt(agree.dat, id.vars = c('indiv', 'mark'))
q = ggplot(agree.dat) + geom_bar(aes(x = indiv, y = value, fill = variable), stat = "identity", position = "stack") + scale_x_discrete('') + 
  facet_wrap(~mark) + scale_y_continuous('# of AS genes overlapping a AS SNPs with the same or different bias') +
  scale_fill_discrete('', labels = c('agree', 'disagree')) +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 11), axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 17), legend.text = element_text(size = 14))
ggsave(file.path(outdir, paste(outsuf, '_vs_as_snps.png', sep = '')))
p = ggplot(snp.dist) + geom_density(aes(x = dist, y = ..density.., color = indiv), adjust = 1) + facet_wrap(~mark) + 
  scale_x_continuous(lim = c(-10^5, 10^5))
ggsave(file.path(outdir, paste(outsuf, '_to_as_snps_dist.png', sep = '')))