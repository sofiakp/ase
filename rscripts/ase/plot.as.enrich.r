rm(list=ls())
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

eq = read.table('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/enrich/hit_eQTL_enrich_LD.tsv', header = T, sep = '\t')[, c(1,6)]
eq = eq[!is.na(eq[, 2]) & !(grepl('Merged', eq[, 1])), ]
tmp = unlist(strsplit(as.character(eq[, 1]), '_'))
indiv = tmp[seq(1,length(tmp),2)]
mark = tmp[seq(2,length(tmp),2)]

pvals = data.frame(indiv = indiv, mark = mark, p = -log10(eq[,2]))

# eq = read.table('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/enrich/hit_dsQTL_enrich_LD.tsv', header = T, sep = '\t')[, c(1,6)]
# eq = eq[!is.na(eq[, 2]) & !(grepl('Merged', eq[, 1])), ]
# tmp = unlist(strsplit(as.character(eq[, 1]), '_'))
# indiv = tmp[seq(1,length(tmp),2)]
# mark = tmp[seq(2,length(tmp),2)]
# 
# pvals = rbind(pvals, data.frame(indiv = indiv, mark = mark, p = -log10(eq[,2]), type = rep('dsQTLs', nrow(eq))))

pvals$mark = order.marks(pvals$mark)
q = ggplot(pvals) + geom_boxplot(aes(x = mark, y = p, fill = mark)) + 
  scale_y_continuous('-log10(p-value)') + scale_color_discrete('') + theme_bw() +
  scale_fill_manual(values =  mark.colors(unique(pvals$mark)), guide = F) +
  scale_x_discrete('') + 
  theme(axis.text.x = element_text(size = 15, angle = -65, vjust = 1, hjust = 0), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file.path('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/enrich/', 'eQTL_pval.pdf'), width = 6.5, height = 5.6)