rm(list=ls())
library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
library(glmnet)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# Computes overlaps between the motif score matrix and the GM12878 TF peaks

#indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/rdata/'
score.file = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/rdata/Gm12878_allTFBS.sorted.noPol.scores.RData'
load(score.file)

tf.regions = read.table('../../rawdata/TFs/Gm12878_allTFBS.sorted.noPol.bed', sep = '\t')
colnames(tf.regions) = c('chr', 'start', 'end', 'name')
tf.peak.names = unlist(strsplit(as.character(tf.regions[, 4]), '_'))
tf.peak.names = tf.peak.names[seq(1, length(tf.peak.names), 3)]
tf.names = gsub('^Haib|^Sydh|^Broad|Std$|^UtaChip|^UtahChip|Iggmus$|Iggrab$|Pcr1x|Pcr2x$|V0416101$|iknucla|^Uw', '', tf.peak.names)
uniq.tfs = unique(tf.names)
ntfs = length(uniq.tfs)

tf.ov = NULL
for(i in 1:ntfs){
  sel.regions = tf.regions[tf.names == uniq.tfs[i], ] # Get peaks of the TF
  ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(sel.regions), select = 'first', ignore.strand = T)
  ov = which(!is.na(ov)) # incides of regions with overlaps
  if(is.null(tf.ov)){
    # ridx is an array of indices in "regions" 
    tf.ov = Matrix(ridx %in% ov) 
  }else{
    tf.ov = cBind(tf.ov, Matrix(ridx %in% ov))
  }
}
colnames(tf.ov) = uniiq.tfs
save(tf.ov, file = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_H3K27AC/rdata/Gm12878_allTFBS.sorted.noPol.tfOv.RData')