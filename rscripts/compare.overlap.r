rm(list=ls())
require('ggplot2')
library('foreach')
library('doMC')
#source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

registerDoMC(6)
# filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/'), 
#                         pattern = '[0-9]+\\.counts\\.r', full.name = T, recursive = F, include.dirs = F)
# outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/', 'reps_6Aug12/plots/')
# if(!file.exists(outdir)) dir.create(outdir)
# 
# all.info <- sample.info(filenames, '\\.counts\\.r$')
# types <- unique(all.info[, 1:2])
# overwrite = F
# # cuts = seq(0, 0.1, by = 0.01)
# hit.dat = NULL
# for(t in 1:dim(types)[1]){
#   reps <- which(all.info$indiv == types$indiv[t] & all.info$mark == types$mark[t])
#   if(length(reps) > 1){
#     for(i in 1:(length(reps) - 1)){
#       rep1 = new.env()
#       rep1.name = all.info$rep[reps[i]]
#       print(rep1.name)
#       load(filenames[reps[i]], rep1)
#       good1 = rep1$snp.info$pass & rep1$snp.info$qval < 0.5 & rep1$snp.info$sb > 0.01
#       for(j in (i + 1):length(reps)){
#         rep2 = new.env()
#         load(filenames[reps[j]], rep2)
#         good2 = rep2$snp.info$pass & rep2$snp.info$qval < 0.5 & rep2$snp.info$sb > 0.01
#         sidx1 = sort(rep1$snp.info$qval[good1 | good2], index.return = T)
#         sidx2 = sort(rep2$snp.info$qval[good1 | good2], index.return = T)
#         rep2.name = all.info$rep[reps[j]]
#         print(rep2.name)
#         
#         ridx1 = sidx1$ix
#         ridx2 = sidx2$ix
#         for(i in 1:length(ridx1)){
#           ridx1[sidx1$ix[i]] = i
#           ridx2[sidx2$ix[i]] = i
#         }
#         cuts = seq(1, length(sidx1$ix), 10)
#         hits = array(0, dim = c(length(cuts), 1))
#         hits = foreach(i = 1:length(cuts), .combine = c) %dopar% {
#           sum(ridx1 <= cuts[i] & ridx2 <= cuts[i])
#         }
#         
# #         hits = array(0, dim = c(length(qcuts), 1))
# #         for(q in 1:length(qcuts)){
# #           hits1 = which(rep1$snp.info$qval < qcuts[q])
# #           hits2 = which(rep2$snp.info$qval < qcuts[q])
# #           hits[q] = sum(hits1 %in% hits2) / min(length(hits1), length(hits2))
# #         }
#         tmp.dat = data.frame(mark = rep(types$mark[t], length(hits)),
#                              indiv = rep(types$indiv[t], length(hits)),
#                              sets = rep(paste(rep1.name, rep2.name, sep = '_vs_'), length(hits)),
#                              hits = hits,
#                              qvals = cuts)
#         if(is.null(hit.dat)){
#           hit.dat = tmp.dat
#         }else{
#           hit.dat = rbind(hit.dat, tmp.dat)
#         }
#       }
#     }
#   }
# }
# 
# lt = factor(gsub('_vs_([0-9]\\.?)+', '', hit.dat$sets))
# col = factor(gsub('([0-9]\\.?)+_vs_', '', hit.dat$sets))
# hit.dat$lt = lt
# hit.dat$col = col
# indivs = unique(all.info$indiv)
# for(i in 1:length(indivs)){
#   tmp.dat = hit.dat[hit.dat$indiv == indivs[i], ]
#   p <- ggplot(tmp.dat) + geom_line(aes(x = qvals, y = hits, color = sets)) +
#     facet_wrap(~mark) + scale_color_discrete('Replicates') + scale_x_continuous("# Elements") + scale_y_continuous("# Overlap") +
#     opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0),
#          title = indivs[i])
#   ggsave(file.path(outdir, paste(indivs[i], '_rep_cum_ovelaps.png', sep = '')))
# }

rep1 = new.env()
load('rawdata/alleleCounts/SNYDER_HG19_GM12878_H3K27AC_1.counts.r', rep1)
good1 = rep1$snp.info$pass & rep1$snp.info$qval < 0.05 & rep1$snp.info$sb > 0.01
rep2 = new.env()
load('rawdata/alleleCounts/SNYDER_HG19_GM12878_H3K27AC_2.counts.r', rep2)
good2 = rep2$snp.info$pass & rep2$snp.info$qval < 0.05 & rep2$snp.info$sb > 0.01

pvals1 = -log(rep1$snp.info$pval[good1 | good2])
pvals2 = -log(rep2$snp.info$pval[good1 | good2])

rank1 = rank(-log(rep1$snp.info$pval[good1 | good2]))
rank2 = rank(-log(rep2$snp.info$pval[good1 | good2]))

c = get.correspondence(rank1, rank2, seq(0.01, 0.99, by=1/28))
sidx1 = sort(rep1$snp.info$pval[good1 | good2], index.return = T)
sidx2 = sort(rep2$snp.info$pval[good1 | good2], index.return = T)
# 
# ridx1 = sidx1$ix
# ridx2 = sidx2$ix
# for(i in 1:length(ridx1)){
#   ridx1[sidx1$ix[i]] = i
#   ridx2[sidx2$ix[i]] = i
# }
# cuts = seq(1, length(sidx1$ix), 10)
# hits = array(0, dim = c(length(cuts), 1))
# hits = foreach(i = 1:length(cuts), .combine = c) %dopar% {
#   sum(ridx1 <= cuts[i] & ridx2 <= cuts[i])
# }