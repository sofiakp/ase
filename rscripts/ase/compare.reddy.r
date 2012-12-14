rm(list=ls())
require(ggplot2)
library(reshape)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/', 'reps_15Sep12', 'plots')
if(!file.exists(outdir)) dir.create(outdir)
reddydat = file.path(Sys.getenv('MAYAROOT'), 'rawdata/reddy/chip/liftover/', 'all_tfs.RData')
load(reddydat)
tfs = names(hit.list)
filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals/hitLists/'), 
                        pattern = 'GM12878.*\\.hits\\.RData', full.name = T, recursive = F, include.dirs = F)
nfiles = length(filenames)
ntfs = length(tfs)
samples = sample.info(filenames, '.hits.RData')
overlaps = array(0, dim = c(nfiles, ntfs))
agree = array(0, dim = c(nfiles, ntfs))

uniq.hits = list()
for(j in 1:ntfs) uniq.hits[[j]] = array(F, dim = c(length(hit.list[[j]]), 1))
for(i in 1:nfiles){
  load(filenames[i])
  print(basename(filenames[i]))
  hits1 = sort(paste(hits$chr, hits$pos, sep = ':'), index.return = T)
  ratios1 = (hits$ref.counts + 1) / (hits$alt.counts + 1)
  ratios1[hits$mat] = 1 / ratios1[hits$mat]
  ratios1 = log(ratios1, base = 2)[hits1$ix]
  
  #pass = !is.nan(snp.info$qval) & snp.info$qval < 0.01
  #mat.idx = grep('mat', colnames(counts))
  #pat.idx = grep('pat', colnames(counts))
  #hist.ratio = log((rowSums(counts[, mat.idx]) + 1) / (rowSums(counts[, pat.idx]) + 1))[pass]
  #hist.names = paste(snp.info$chr, snp.info$pos, sep = ':')[pass]
  
  for(j in 1:ntfs){
    hits2 = sort(hit.list[[j]], index.return = T)
    ratios2 = log(ratios[[j]], base = 2)[hits2$ix]
    sel1 = hits1$x %in% hits2$x
    sel2 = hits2$x %in% hits1$x
    uniq.hits[[j]] = uniq.hits[[j]] | sel2
    overlaps[i, j] = sum(sel1)
    agree[i, j] = sum(ratios1[sel1] * ratios2[sel2] > 0)
    #sel1 = which(hi %in% hit.list[[j]])
    #sel2 = array(0, dim = c(length(sel1), 1))
    #ind2 = array(F, dim = c(length(hit.list[[j]]), 1))
    #for(k in 1:length(hit.list[[j]])){
    #  idx = which(hit.list[[j]][k] == hist.names[sel1])
    #  if(length(idx) > 0){
    #    sel2[idx] = k
    #    ind2[k] = T
    #  }
    #}
    #uniq.hits[[j]] = uniq.hits[[j]] | ind2
    #ratios1 = hist.ratio[sel1]
    #ratios2 = log(ratios[[j]][sel2])
    #overlaps[i, j] = length(sel1)
    #agree[i, j] = sum(ratios1 * ratios2 > 0)
  }
}

ov.dat = data.frame(overlaps)
colnames(ov.dat) = tfs
ov.dat$mark = samples$mark
ov.dat = melt(ov.dat)
p = ggplot(ov.dat) + geom_tile(aes(x = mark, y = variable, fill = value)) +
  scale_x_discrete('') + scale_y_discrete('') + scale_fill_continuous('Overlaps') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30'),
        axis.text.y = element_text(hjust = 1, colour = 'grey30'),
        title = 'Number of overlapping SNPs with allelic bias')
ggsave(file.path(outdir, 'reddy_overlaps.png'))
nov.tfs = data.frame(num = unlist(lapply(hit.list, length)), tf = tfs, type = factor(rep('total', ntfs)))
nov.tfs = rbind(nov.tfs, data.frame(num = unlist(lapply(uniq.hits, sum)), tf = tfs, type = factor(rep('overlapping', ntfs))))
p = ggplot(nov.tfs) + geom_bar(aes(x = tf, y = num, fill = type), stat = "identity", position = "dodge") + 
  scale_y_continuous('# AS SNPs') + scale_x_discrete('') + scale_fill_discrete('') + 
  theme(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file.path(outdir, 'reddy_overlaps_bar.png'))