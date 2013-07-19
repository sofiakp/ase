rm(list=ls())
library(ggplot2)
library(reshape)
library(GenomicRanges)
library(matrixStats)
source('utils/deseq.utils.r')
source('utils/binom.val.r')

# Reads a set of TF regions, overlaps them with a set of signal regions (eg. H3K27AC peaks) and outputs the log-odds
# signal between each pair of individuals in all signal regions with overlaps.

bed.file = '../../rawdata/TFs/Gm12878_allTFBS.sorted.noPol.merged.bed'
#good.idx.file = '../../rawdata/TFs/jaspar/Gm12878_allTFBS.sorted.noPol.diffRows.txt'
inpref = 'Gm12878_allTFBS.sorted.noPol.'
#insuf = '.jaspar.scores'
#indir = '../../rawdata/TFs/jaspar/'
#tf.file = '../../rawdata/TFs/jaspar/motif_names.txt'
signal.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/'
mark = 'H3K27AC'
outdir = paste('../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/tfRegress_at_', mark, sep = '')
if(!file.exists(outdir)) dir.create(outdir)

bed.regions = read.bed(bed.file)
#good.idx = read.table(good.idx.file)
#bed.regions = bed.regions[good.idx[, 1] == 1, ]
#hit.files = list.files(indir, pattern = inpref, full.names = T)
#input.indivs = fix.indiv.names(gsub(paste(inpref, '[mp]aternal', insuf, '.txt', sep = ''), '', basename(hit.files)))
#nindivs = length(input.indivs)

#tfs = read.table(tf.file)
#ntfs = length(tfs)

load(file.path(signal.dir, 'rdata', paste('SNYDER_HG19_all_reg_', mark, '_qn.RData', sep = '')))
regions$start = pmax(regions$start - 200, 1)
regions$end = regions$end + 200
regions = regions[good.rows, ]
counts = counts[good.rows, ]
indivs = fix.indiv.names(colnames(counts))
nindivs = length(indivs)

# We only care about the signal at the signal regions with overlaps with TFs 
ov = findOverlaps(regions.to.ranges(regions), regions.to.ranges(bed.regions), select = 'all', ignore.strand = T)
ov.mat = cbind(queryHits(ov), subjectHits(ov))
stopifnot(all(ov.mat[,1]==sort(ov.mat[,1]))) # Make sure the first column is already sorted.
sel.regions = sort(unique(ov.mat[, 1]))
cat('Percentage of regions with overlap ', length(sel.regions) * 100 / nrow(regions), '\n')
cat('Average number of overlaps per region ', nrow(ov.mat) / length(sel.regions), '\n')

region.max = rowMaxs(counts[sel.regions, ])
for(i in 1:(nindivs - 1)){  
  for(j in (i + 1):nindivs){
    log.scores = counts[sel.regions, i] - counts[sel.regions, j]
    weights = 2 / (1 + exp(region.max - pmax(counts[sel.regions, i], counts[sel.regions, j])))
    log.scores = weights * log.scores
    outfile = file.path(outdir, paste(inpref, indivs[i], '_vs_', indivs[j], '_', mark, '_v2.txt', sep = ''))
    write.table(log.scores, outfile, sep = '\t', row.names = F, col.names = F, quote = F)
  }
}
    
idx.map = array(0, dim = c(nrow(regions), 1))
idx.map[1:nrow(regions) %in% sel.regions] = 1:length(sel.regions)
ov.mat[, 1] = idx.map[ov.mat[, 1]]
write.table(ov.mat, file.path(outdir, paste(inpref, 'idx_', mark, '.txt', sep = '')), 
            sep = '\t', row.names = F, col.names = F, quote = F)
out.regions = regions[sel.regions, ]
write.table(out.regions, file.path(outdir, paste(inpref, 'regions_', mark, '.txt', sep = '')), 
            sep = '\t', row.names = F, col.names = F, quote = F)
# Finally, output the overlaps between signal regions and TF regions, so we can create the log-odd motif scores
# in matlab.

#scores1 = (read.table(file.path(indir, paste(inpref, input.indivs[i], 'maternal', insuf, '.txt', sep = '\t')), sep = '\t') + 
#  read.table(file.path(indir, paste(inpref, input.indivs[i], 'paternal', insuf, '.txt')), sep = '\t')) / 2