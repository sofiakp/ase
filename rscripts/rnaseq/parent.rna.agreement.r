rm(list=ls())
require(ggplot2)
library('DESeq')
source('utils/binom.val.r')
source('utils/sample.info.r')
source('utils/deseq.utils.r')

countdir = '../../rawdata/geneCounts/rdata/reps/qvals'
plotdir = '../../rawdata/geneCounts/rdata/reps/plots'
if(!file.exists(plotdir)) dir.create(plotdir)
pop = 'YRI'

if(pop == 'CEU'){
  child = 'GM12878'
  mom = 'GM12892'
  dad = 'GM12891'
}else{
  child = 'GM19240'
  mom = 'GM19238'
  dad = 'GM19239'
}

load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

filenames <- list.files(countdir, pattern = paste(child, '.*rep\\.RData', sep = ''), full.name = T)
filenames = filenames[!grepl(x = filenames, pattern = '_RNA_rep')]
epsilon <- 1
nfiles = length(filenames)
rdat = NULL
for(i in 1:nfiles){
  mat.file = gsub(child, mom, filenames[i])
  pat.file = gsub(child, dad, filenames[i])
  
  if(!file.exists(mat.file) || !file.exists(pat.file)) next
  
  c = new.env()
  load(filenames[i], envir = c)
  mat.idx <- grep('mat', colnames(c$counts))
  pat.idx <- grep('pat', colnames(c$counts))
  child.mat <- c$counts[, mat.idx]
  child.pat <- c$counts[, pat.idx]
  child.ratio <- log((child.mat + epsilon) / (child.pat + epsilon), base = 2)
  
  m = new.env()
  load(mat.file, m)
  mat.counts <- rowSums(m$counts)
  p = new.env()
  load(pat.file, p)
  pat.counts <- rowSums(p$counts)
  sf = estimateSizeFactorsForMatrix(cbind(mat.counts, pat.counts))
  par.ratio <- log((mat.counts / sf[1] + epsilon) / (pat.counts / sf[2] + epsilon), base = 2)
  
  pass.c =  !is.nan(c$regions$qval) & c$regions$qval < 0.01 & c$regions$phased &
    !(as.character(gene.meta$chr) %in% c('chrX', 'chrY')) & c$regions$het #& !m$regions$het
  print(sum(pass.c))
  if(sum(pass.c) < 50) next
  
  child.ratio = child.ratio[pass.c]
  par.ratio = par.ratio[pass.c]
  corr = cor(child.ratio, par.ratio)
  agree.frac = sum(child.ratio * par.ratio > 0) / length(child.ratio)
  rtmp = data.frame(child = child.ratio, par = par.ratio, mark1 = sample.info(filenames[i])$mark[1], mark = paste(sample.info(filenames[i])$mark[1], sprintf('%.4f', corr)))
  rtmp2 = data.frame(mark = sample.info(filenames[i])$mark, f = agree.frac, c = corr)
  if(is.null(rdat)){ 
    rdat = rtmp
    rdat2 = rtmp2
  }else{ 
    rdat = rbind(rdat, rtmp)
    rdat2 = rbind(rdat2, rtmp2)
  }
}

rdat$mark1 = factor(gsub('RNA', 'polyA-RNA', as.character(rdat$mark1)))
rdat$mark1 = factor(gsub('RZ', 'RNA', as.character(rdat$mark1)))
rdat$mark = factor(gsub('RNA', 'polyA-RNA', as.character(rdat$mark)))
rdat$mark = factor(gsub('RZ', 'RNA', as.character(rdat$mark)))
rdat2$mark = factor(gsub('RNA', 'polyA-RNA', as.character(rdat2$mark)))
rdat2$mark = factor(gsub('RZ', 'RNA', as.character(rdat2$mark)))

uniq.marks = order.marks(levels(rdat$mark1), sub.rna = F)
uniq.marks.corr = levels(rdat$mark)
tmp = unlist(strsplit(uniq.marks.corr, ' ')) 
uniq.marks.corr = uniq.marks.corr[match(uniq.marks, tmp[seq(1, length(tmp), 2)])]

other = new.env()
load(paste('../../rawdata/alleleCounts/allNonSan/rdata/reps/plots/', pop, '.child_par_cor_plot.RData', sep = ''), other)
other$rdat = droplevels(other$rdat[!(as.character(other$rdat$mark1) %in% c('H2AZ', 'H3K9AC', 'POL4H8')), ])
other$rdat2 = droplevels(other$rdat2[!(as.character(other$rdat2$mark) %in% c('H2AZ', 'H3K9AC', 'POL4H8')), ])
merged = rbind(other$rdat, rdat)
merged2 = rbind(other$rdat2, rdat2)
uniq.marks = order.marks(levels(merged$mark1), sub.rna = F)
uniq.marks.corr = levels(merged$mark)
tmp = unlist(strsplit(uniq.marks.corr, ' ')) 
uniq.marks.corr = uniq.marks.corr[match(uniq.marks, tmp[seq(1, length(tmp), 2)])]

q1 = ggplot(merged) + geom_point(aes(x = child, y = par, color = mark1), shape = 1, alpha = 0.8, size = 2) + scale_x_continuous('log(maternal/paternal)') + 
  scale_y_continuous('child:log(maternal/paternal)') + theme_bw() + 
  scale_color_manual(values =  mark.colors(levels(merged$mark1)), 
                     breaks = uniq.marks, labels = uniq.marks.corr) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14), panel.grid.major = element_blank())
for(i in 1:nrow(merged2)){
  co = coef(lm(merged[merged$mark1 == as.character(merged2$mark[i]), 2] ~ merged[merged$mark1 == as.character(merged2$mark[i]), 1]))
  q1 = q1 + geom_abline(intercept = co[1], slope = co[2], color = mark.colors(merged2$mark[i]), size = 0.5)
  #q1 = q1 + geom_abline(intercept = co[1], slope = co[2], color = mark.colors(merged2$mark[i]), size = 0.3)
}
ggsave(file = file.path(plotdir, paste(pop, '.child_par_gene_cor.pdf', sep = '')), q1, width = 8, height = 5.6)
q2 = ggplot(rdat2) + geom_bar(aes(x = mark, y = f, fill = mark), stat = "identity") + scale_x_discrete('') + 
  scale_y_continuous('Fraction of AS genes with same direction') + theme_bw() + 
  scale_fill_manual(values =  mark.colors(levels(rdat2$mark)), guide = F) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 14), axis.text.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14))
ggsave(file = file.path(plotdir, paste(pop, '.child_par_gene_agreement.pdf', sep = '')), q2, width = 6.5, height = 5.6)