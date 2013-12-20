# Basic QC plots (duplication rate, # reads, NSC...)

rm(list=ls())
library(ggplot2)
library(reshape)
source('utils/sample.info.r')
source('utils/binom.val.r')
source('utils/deseq.utils.r')

outpref = '../../rawdata/mapped/bam/personal/stats/qcStats_Jun13'
qc.dat = read.table('../../rawdata/mapped/bam/personal/stats/qc.stats', header = F, sep = '\t')[, c(1,3,9)]
qc.dat = rbind(qc.dat, read.table('../../rawdata/mapped/bam/reference/stats/qc.stats', header = F, sep = '\t')[, c(1,3,9)])
qc.dat[, 1] = gsub('_reconcile.dedup.bam|_dedup.bam', '', qc.dat[, 1])
qc.dat[, 2] = as.numeric(gsub(',.*', '', qc.dat[, 2]))
colnames(qc.dat) = c('dataset', 'fraglen', 'nsc')
qc.dat = qc.dat[order(qc.dat$dataset), ]
qc.dat = qc.dat[!duplicated(qc.dat$dataset), ]

read.counts = read.table('../../rawdata/mapped/bam/personal/stats/readCounts_q30_14Mar13.txt', header = F, sep = '\t')
#read.counts = rbind(read.counts, read.table('../../rawdata/mapped/bam/reference/stats/readCounts.txt', header = F, sep = '\t'))
read.counts[, 1] = gsub('_reconcile.dedup.bed|_dedup.bed', '', read.counts[, 1])
colnames(read.counts) = c('dataset', 'reads')
qc.dat$reads = read.counts[match(qc.dat[,1], read.counts[,1]), 2]

rna.counts = read.table('../../rawdata/mapped/bam/personal/stats/readCounts_rna_Apr13.txt', header = F, sep = '\t')
rna.counts[, 1] = gsub('_reconcile.dedup', '', rna.counts[, 1])
colnames(rna.counts) = c('dataset', 'reads')
rna.counts$fraglen = rep(NA, nrow(rna.counts))
rna.counts$nsc = rep(NA, nrow(rna.counts))
qc.dat = rbind(qc.dat, rna.counts)

qc.dat$mark = order.marks(sample.info(qc.dat[,1], '')$mark)
qc.dat$indiv = as.character(sample.info(qc.dat[,1], '')$indiv)
qc.dat$indiv = fix.indiv.names(qc.dat$indiv)
qc.dat$rep = sample.info(qc.dat[,1], '')$rep
qc.dat = qc.dat[!(qc.dat$mark %in% c('BUB', 'EPOL', 'POL4H8', 'H2AZ', 'RZ', 'H3K9AC', 'H3K9ME3', 'polyA-RNA', 'PU1')), ]
write.table(qc.dat[is.na(qc.dat$fraglen), -c(2,3)], file = paste(outpref, '_rna.txt', sep = ''), col.names = T, row.names = F, quote = F, sep = '\t')

qc.dat2 = cast(qc.dat, mark+indiv~., function(x) sum(x), value = 'reads')
colnames(qc.dat2)[3] = 'reads'
q4 = ggplot(qc.dat2) + geom_bar(aes(x = indiv, y = reads/1e6), stat = "identity", position = "dodge") + theme_bw() +
  facet_wrap(~mark, scales = 'free_y') + scale_x_discrete("") + scale_y_continuous("# Q30 reads after duplicate removal (millions)") +
  theme(axis.text.x = element_text(angle = -65, vjust = 1, hjust = 0, size = 14, color = get.pop.col(get.pop(unique(qc.dat$indiv)))), 
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20), strip.text = element_text(size = 20))
ggsave(paste(outpref, '_qreads.pdf', sep = ''), width = 13.6, height = 11.8)

qc.dat = qc.dat[!is.na(qc.dat$fraglen), ]
q3 = ggplot(qc.dat) + geom_boxplot(aes(x = mark, y = nsc), outlier.size = 1) + 
  geom_jitter(aes(x = mark, y = nsc, color = indiv), size = 1) + scale_y_continuous('NSC') + scale_color_discrete('') +
  scale_x_discrete('') + theme_bw() + 
  theme(axis.text.x = element_text(size = 14, angle = -45, vjust = 1, hjust = 0),axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_blank(), legend.text = element_text(size = 14))
ggsave(paste(outpref, '_qc.pdf', sep = ''), width = 6.5, height = 5.6)
write.table(qc.dat, file = paste(outpref, '.txt', sep = ''), col.names = T, row.names = F, quote = F, sep = '\t')

# indir = '/media/fusion10/work/chromatinVariation/rawdata/mapped/bam/personal/stats/'
# dedup.file = file.path(indir, 'dedupStats_5Oct12.txt')
# reconcile.file = file.path(indir, 'reconcileStats_5Oct12.txt')
# counts.file = file.path(indir, 'readCounts_q30_5Oct12.txt')
# qc.file = file.path(indir, 'qc.stats')
# outpref = file.path(indir, 'stats_5Oct12')
# 
# # Duplication rates
# tab = read.table(dedup.file, header = F, skip = 1)
# tab[, 1] = gsub('_reconcile.dedup', '', tab[, 1])
# tab = tab[order(tab[, 1]), ]
# expt = tab[, 1]
# samples = sample.info(as.character(expt), '')
# pat.reads = tab[, 2] * (1 - tab[, 3]) # paternal reads after duplication removal # TODO THIS IS NOT CORRECT
# mat.reads = tab[, 4] * (1 - tab[, 5])
# amb.reads = tab[, 6] * (1 - tab[, 7])
# dup <- (tab[,2] * tab[,3] + tab[,4] * tab[,5] + tab[,6] * tab[,7]) / (tab[,2] + tab[,4] + tab[,6]) # avg duplication rate for the three read groups
# ndup.reads <- (pat.reads + mat.reads + amb.reads)
# nreads = tab[, 2] + tab[, 4] + tab[, 6]
# 
# # Reconciliation stats
# tab = read.table(reconcile.file, header = F, skip = 1)
# tab[, 1] = gsub('_reconcile', '', tab[, 1])
# tab = tab[order(tab[, 1]), ]
# stopifnot(all(expt == tab[, 1]))
# #nreads = tab[, 2]
# pat.ratio = as.numeric(gsub('%', '', tab[, 4]))
# mat.ratio = as.numeric(gsub('%', '', tab[, 6]))
# 
# counts = data.frame(indiv = samples$indiv, mark = samples$mark, rep = samples$rep, nreads = nreads, dup = dup, 
#                     pat.ratio = pat.ratio, mat.ratio = mat.ratio, ndup.reads = ndup.reads)
# counts2 = melt(data.frame(indiv = samples$indiv, mark = samples$mark, rep = samples$rep, 
#                           pat.ratio = pat.ratio, mat.ratio = mat.ratio), id.vars = c('indiv', 'mark', 'rep'))
# # Plot the duplication rate per replicate
# q1 = ggplot(counts) + geom_bar(aes(x = indiv, y = dup, fill = rep), stat = "identity", position = "dodge") +
#   facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Duplication rate") +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
# ggsave(paste(outpref, '_dup.png', sep = ''), width = 13.6, height = 11.8)
# 
# # Plot the fraction of total reads that were unambiguously assigned to each parent
# q2 = ggplot(counts2) + geom_bar(aes(x = indiv, y = value, color = rep, fill = variable), stat = "identity", position = "dodge") +
#   facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Read group percentage") +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
# ggsave(paste(outpref, '_unambig_ratios.png', sep = ''), width = 13.6, height = 11.8)
# 
# if(file.exists(qc.file)){
#   # Plot the normalized strand cross-correlation per repilcate
#   tab = read.table(qc.file, header = F, skip = 0)
#   tab[, 1] = gsub('_reconcile.dedup.bam', '', tab[, 1])
#   tab = tab[order(tab[, 1]), ]
#   stopifnot(all(expt == tab[, 1]))
#   counts$nsc = tab[, 9]
#   q3 = ggplot(counts) + geom_boxplot(aes(x = mark, y = nsc)) + 
#     geom_jitter(aes(x = mark, y = nsc, color = indiv)) + scale_y_continuous('NSC') + scale_color_discrete('') +
#     scale_x_discrete('') + theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
#   ggsave(paste(outpref, '_qc.png', sep = ''), width = 13.6, height = 11.8)
# }
# if(file.exists(counts.file)){
#   tab = read.table(counts.file, header = F, skip = 0)
#   tab[, 1] = gsub('_reconcile.dedup', '', tab[, 1])
#   sel = tab[, 1] %in% expt
#   tab = tab[sel, ] # Remove maternal and paternal read counts if present
#   tab = tab[order(tab[, 1]), ]
#   stopifnot(all(expt == tab[, 1]))
#   counts$nqreads = tab[, 2]
#   counts$nqfrac = counts$nqreads / counts$ndup.reads
#   
#   # Plot the number and fraction of Q30 reads
# #   q4 = ggplot(counts) + geom_bar(aes(x = indiv, y = nqfrac, fill = rep), stat = "identity", position = "dodge") +
# #     facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Fraction of Q30 reads (out of non-duplicate reads)") +
# #     theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
# #   ggsave(paste(outpref, '_qreads_rate.png', sep = ''), width = 13.6, height = 11.8)
#   q4 = ggplot(counts) + geom_bar(aes(x = indiv, y = nqreads , fill = rep), stat = "identity", position = "dodge") +
#     facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("# Q30 reads after duplicate removal") +
#     theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 11), axis.text.y = element_text(size = 12),
#           axis.title.y = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 17))
#   ggsave(paste(outpref, '_qreads.png', sep = ''), width = 13.6, height = 11.8)
# }