# Basic QC plots (duplication rate, # reads, NSC...)

rm(list=ls())
library(ggplot2)
library(reshape)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

indir = '/media/fusion10/work/chromatinVariation/rawdata/mapped/bam/personal/stats/'
dedup.file = file.path(indir, 'dedupStats_5Oct12.txt')
reconcile.file = file.path(indir, 'reconcileStats_5Oct12.txt')
counts.file = file.path(indir, 'readCounts_q30_5Oct12.txt')
qc.file = file.path(indir, 'qc.stats')
outpref = file.path(indir, 'stats_5Oct12')

# Duplication rates
tab = read.table(dedup.file, header = F, skip = 1)
tab[, 1] = gsub('_reconcile.dedup', '', tab[, 1])
tab = tab[order(tab[, 1]), ]
expt = tab[, 1]
samples = sample.info(as.character(expt), '')
pat.reads = tab[, 2] * (1 - tab[, 3]) # paternal reads after duplication removal # TODO THIS IS NOT CORRECT
mat.reads = tab[, 4] * (1 - tab[, 5])
amb.reads = tab[, 6] * (1 - tab[, 7])
dup <- (tab[,2] * tab[,3] + tab[,4] * tab[,5] + tab[,6] * tab[,7]) / (tab[,2] + tab[,4] + tab[,6]) # avg duplication rate for the three read groups
ndup.reads <- (pat.reads + mat.reads + amb.reads)
nreads = tab[, 2] + tab[, 4] + tab[, 6]

# Reconciliation stats
tab = read.table(reconcile.file, header = F, skip = 1)
tab[, 1] = gsub('_reconcile', '', tab[, 1])
tab = tab[order(tab[, 1]), ]
stopifnot(all(expt == tab[, 1]))
#nreads = tab[, 2]
pat.ratio = as.numeric(gsub('%', '', tab[, 4]))
mat.ratio = as.numeric(gsub('%', '', tab[, 6]))

counts = data.frame(indiv = samples$indiv, mark = samples$mark, rep = samples$rep, nreads = nreads, dup = dup, 
                    pat.ratio = pat.ratio, mat.ratio = mat.ratio, ndup.reads = ndup.reads)
counts2 = melt(data.frame(indiv = samples$indiv, mark = samples$mark, rep = samples$rep, 
                          pat.ratio = pat.ratio, mat.ratio = mat.ratio), id.vars = c('indiv', 'mark', 'rep'))
# Plot the duplication rate per replicate
q1 = ggplot(counts) + geom_bar(aes(x = indiv, y = dup, fill = rep), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Duplication rate") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave(paste(outpref, '_dup.png', sep = ''), width = 13.6, height = 11.8)

# Plot the fraction of total reads that were unambiguously assigned to each parent
q2 = ggplot(counts2) + geom_bar(aes(x = indiv, y = value, color = rep, fill = variable), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Read group percentage") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave(paste(outpref, '_unambig_ratios.png', sep = ''), width = 13.6, height = 11.8)

if(file.exists(qc.file)){
  # Plot the normalized strand cross-correlation per repilcate
  tab = read.table(qc.file, header = F, skip = 0)
  tab[, 1] = gsub('_reconcile.dedup.bam', '', tab[, 1])
  tab = tab[order(tab[, 1]), ]
  stopifnot(all(expt == tab[, 1]))
  counts$nsc = tab[, 9]
  q3 = ggplot(counts) + geom_boxplot(aes(x = mark, y = nsc)) + 
    geom_jitter(aes(x = mark, y = nsc, color = indiv)) + scale_y_continuous('NSC') + scale_color_discrete('') +
    scale_x_discrete('') + theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
  ggsave(paste(outpref, '_qc.png', sep = ''), width = 13.6, height = 11.8)
}
if(file.exists(counts.file)){
  tab = read.table(counts.file, header = F, skip = 0)
  tab[, 1] = gsub('_reconcile.dedup', '', tab[, 1])
  sel = tab[, 1] %in% expt
  tab = tab[sel, ] # Remove maternal and paternal read counts if present
  tab = tab[order(tab[, 1]), ]
  stopifnot(all(expt == tab[, 1]))
  counts$nqreads = tab[, 2]
  counts$nqfrac = counts$nqreads / counts$ndup.reads
  
  # Plot the number and fraction of Q30 reads
#   q4 = ggplot(counts) + geom_bar(aes(x = indiv, y = nqfrac, fill = rep), stat = "identity", position = "dodge") +
#     facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Fraction of Q30 reads (out of non-duplicate reads)") +
#     theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
#   ggsave(paste(outpref, '_qreads_rate.png', sep = ''), width = 13.6, height = 11.8)
  q4 = ggplot(counts) + geom_bar(aes(x = indiv, y = nqreads , fill = rep), stat = "identity", position = "dodge") +
    facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("# Q30 reads after duplicate removal") +
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 11), axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 17))
  ggsave(paste(outpref, '_qreads.png', sep = ''), width = 13.6, height = 11.8)
}