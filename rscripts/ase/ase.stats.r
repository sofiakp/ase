rm(list=ls())
library(GenomicRanges)
library(matrixStats)
library(reshape)
library(ggplot2)
library(Matrix)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

snp.pos.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')
load(snp.pos.file)
nsnps = dim(snp.pos)[1]

geno.dir = '../../rawdata/variants/all/snps/allNonSan/' # Directory with genotype data, should have a file <indiv>.snps.RData for each individual.
count.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

filenames = list.files(count.dir, pattern = paste('.*rep\\.RData', sep = ''), full.name = T)
samples = sample.info(filenames)
marks = as.character(samples$mark)
indivs = as.character(samples$indiv)
uniq.indivs = unique(indivs)
uniq.marks = unique(marks)

# motifs
features = read.table('../../rawdata/motifs/personal_genome_motifs.txt', header = F, sep = '\t')[, c(1:3)]
colnames(features) = c('chr', 'start', 'end')
feat.name = 'motif'

dist.dat = NULL
for(i in 1:length(uniq.marks)){
  file.idx = which(marks == uniq.marks[i]) # Get all the files for the mark
  pass = NULL # Indicators of AS SNPs for the mark in any individual
  for(j in file.idx){
    print(basename(filenames[j]))
    load(filenames[j])
    if(is.null(pass)){
      pass = as.vector(!snp.info$bad & snp.info$qval < -2)
    }else{
      pass = as.vector(pass | (!snp.info$bad & snp.info$qval < -2))
    }
  }
  if(sum(pass) < 100) next
  snp.pos.tmp = snp.pos[pass, ]
  dist = array(NaN, dim = c(sum(pass), 1))
  uniq.chr = unique(as.character(snp.pos.tmp$chr))
  for(c in uniq.chr){
    sel = snp.pos.tmp$chr == c
    dist.tmp = distanceToNearest(ranges(snps.to.ranges(snp.pos.tmp[sel, ])), 
                                 ranges(regions.to.ranges(features[features$chr == c, ])))
    dist[which(sel)[dist.tmp$queryHits]] = dist.tmp$distance
  }
  rtmp = data.frame(dist = dist, mark = rep(uniq.marks[i], sum(pass)))
  dist.dat = rbind(dist.dat, rtmp)
}

p = ggplot(dist.dat) + geom_density(aes(x = log(dist, base = 10), y = ..density.., color = mark)) + 
  scale_x_continuous(paste('Distance to the nearest', feat.name)) + theme_bw() + scale_y_continuous('Density') +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14))
ggsave(file = file.path(plotdir, paste('as_snp_', feat.name, '_dist.png', sep = '')), p, width = 6.5, height = 5.6)