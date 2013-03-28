rm(list=ls())
library(matrixStats)
library(reshape)
library(ggplot2)
library(preprocessCore)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source('/media/fusion10/work/sofiakp/scott/rfiles/plot.heatmap.R')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/isvaFn.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/DoISVA.R'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/isva/EstDimRMT.R'))

# Compares population specific clusters to segmentation states

isva.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/signal/combrep/extractSignal/fc/avgSig/rdata/') 
outpref = 'SNYDER_HG19_all_reg_' #'hg19_w10k_' #gencode.v13.annotation.noM.genes_' #
suf = '_qn_isvaNull'
plotdir = file.path(isva.dir, '..', 'plots', 'qn_isvaNull_fits_all_reg')
if(!file.exists(plotdir)) dir.create(plotdir, recursive=T)

seg.dir = '../../rawdata/segmentations/chmmResults/14indivCore/final/'
seg.files = list.files(seg.dir, pattern = 'mnemonics.bed.gz', full.names = T)
state.ids = read.table('../../rawdata/segmentations/chmmResults/14indivCore/final/states.txt', header = F, sep = '\t')
colnames(state.ids) = c('id', 'name', 'color')
state.ids$full.name = paste(state.ids[,1], state.ids[,2], sep = '_')
nstates = nrow(state.ids)

mark = 'H3K27AC'

# Reference clusters
load(file.path(isva.dir, paste(outpref, mark, suf, '.RData', sep = '')))
indivs = colnames(counts)
nindiv = length(indivs)
load(file.path(isva.dir, paste(outpref, mark, suf, '_clust.RData', sep = '')))

# State overlaping each population specific region in each individual
seg.ov = array(0, dim = c(isva.fit$ndeg, nindiv))
for(i in 1:nindiv){
  sel.file = seg.files[grep(gsub('MS1', 'SNYDER', indivs[i]), seg.files)]
  segments = read.table(sel.file, header = F, sep = '\t')
  colnames(segments) = c('chr', 'start', 'end', 'seg')
  segments$start = segments$start + 1 # convert to 1-based 
  ov = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(segments), select = 'first', ignore.strand = T)
  states = unlist(strsplit(as.character(segments$seg[ov[!is.na(ov)]]), '_'))
  seg.ov[!is.na(ov), i] = as.numeric(states[seq(1, length(states), 2)])
}

indivs[indivs == 'SNYDER'] = 'MS1'

# For each cluster, plot a heatmap with the number of overlaps for each individual and state
stats.all = NULL
for(k in 1:max(kord)){
  mem = kclusters$cluster == kord[k]
  stats = array(0, dim = c(nindiv, nstates))
  for(i in 1:nindiv){
    t = table(seg.ov[mem, i])
    stats[i, as.numeric(names(t))] = t
  }
  rownames(stats) = indivs
  colnames(stats) = state.ids$full.name
  if(is.null(stats.all)){
    stats.all = stats / sum(mem)
  }else{
    stats.all = cbind(stats.all, stats / sum(mem))
  }
  plot.heatmap(t(stats[plot.cols, ]), row.cluster = F, col.cluster = F, show.dendro = "none", row.title= '', col.title = '',
               dist.metric = "euclidean", clust.method = "ward", margins = c(15, 15), keysize = 1,
               break.type = 'linear', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], 
               ColSideColors = get.pop.col(get.pop(indivs[plot.cols])), cex.row = 3, cex.col = 3,
               to.file = file.path(plotdir, paste(outpref, mark, suf, '_k6_clust', kord[k], '_statesHeat.pdf', sep = '')))
}
for(s in 1:nstates){
  sel.col = grep(paste('^', s, '_', sep = ''), colnames(stats.all))
  plot.heatmap(t(stats.all[plot.cols, sel.col]), row.cluster = F, col.cluster = F, show.dendro = "none", row.title= '', col.title = '', 
               dist.metric = "euclidean", clust.method = "ward", lab.row = paste('cluster', 1:max(kord)), title.name = state.ids$full.name[s],
               break.type = 'linear', filt.thresh = NA, replace.na = F, palette = brewer.pal(9,  "RdYlBu")[seq(9,1,-1)], 
               ColSideColors = get.pop.col(get.pop(indivs[plot.cols])), cex.row = 3, cex.col = 3, margins = c(15, 15), keysize = 1,
               to.file = file.path(plotdir, paste(outpref, mark, suf, '_k6_state', s, '_clustHeat.pdf', sep = '')))
}

seg.2 = data.frame(seg.ov[plot.rows,])
rownames(seg.2) = 1:nrow(seg.2)
colnames(seg.2) = indivs
seg.2$r = rownames(seg.2)
seg.flat = melt(seg.2, id.vars = 'r')
seg.flat$variable = ordered(factor(seg.flat$variable, levels = indivs[plot.cols]))
seg.flat$r = ordered(factor(seg.flat$r, levels = seq(nrow(seg.2), 1, -1)))
seg.flat$value = factor(seg.flat$value, levels = 1:15, labels = state.ids$full.name)
colors = get.hex(state.ids$color)
names(colors) = state.ids$full.name
p = ggplot(seg.flat) + geom_raster(aes(x = variable, y = r, fill = value)) + xlab('') + ylab('') +
  scale_fill_manual(values = colors, name = 'State') + theme_bw() + 
  theme(axis.text.x = element_text(size = 16, angle = -90, vjust = 0, hjust = 0), axis.title.x = element_text(size = 15), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 16), legend.title = element_blank())
ggsave(file.path(plotdir, paste(outpref, mark, suf, '_k6_states.pdf', sep = '')), p, width = 6.5, height = 5.6)