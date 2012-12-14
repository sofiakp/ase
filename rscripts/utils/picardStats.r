rm(list=ls())
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/mapped/bam/personal/picardStats/'), 
                        pattern = '_qual.txt', full.name = T, recursive = F, include.dirs = F)
samples = sample.info(filenames, '_reconcile_qual.txt')
nsamples = length(filenames)

mean.qual = array(0, dim = c(nsamples, 1))
q20 = array(0, dim = c(nsamples, 1))

for(i in 1:nsamples){
  tab = read.table(filenames[i], header = T)
  q = as.numeric(tab[, 1])
  count = as.numeric(tab[, 2])
  nbases = sum(count)
  mean.qual[i] = sum(q * count) / nbases
  q20[i] = sum(count[q >= 20]) / nbases
}

counts = data.frame(mark = samples$mark, indiv = samples$indiv, rep = samples$rep, 
                    q20 = q20, m = mean.qual)
q <- ggplot(counts) + geom_bar(aes(x = indiv, y = q20, fill = rep), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Fraction of q20 bases") +
  opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/mapped/bam/personal/picardStats/q20_bases.png'))
p <- ggplot(counts) + geom_bar(aes(x = indiv, y = m, fill = rep), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Avg base qual") +
  opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/mapped/bam/personal/picardStats/mean_base_qual.png'))

tab = read.table(file.path(Sys.getenv('MAYAROOT'), 'rawdata/mapped/bam/personal/picardStats/all_insSize.txt'))
samples = sample.info(as.character(tab[,1]), '_reconcile_insSize.txt')
counts = data.frame(mark = samples$mark, indiv = samples$indiv, rep = samples$rep, 
                    m = tab[, 2], sd = tab[, 3])
p <- ggplot(counts) + geom_bar(aes(x = indiv, y = m, fill = rep), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + scale_y_continuous("Mean fragment size") +
  opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/mapped/bam/personal/picardStats/mean_ins_size.png'))