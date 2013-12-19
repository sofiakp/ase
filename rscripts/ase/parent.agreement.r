rm(list=ls())
require(ggplot2)
library(DESeq)
source('utils/binom.val.r')
source('utils/sample.info.r')
source('utils/deseq.utils.r')

geno.dir = '../../rawdata/variants/all/snps/allNonSan/' # Directory with genotype data, should have a file <indiv>.snps.RData for each individual.
count.dir = '../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals'
plotdir = '../../rawdata/alleleCounts/allNonSan/rdata/reps/plots'
if(!file.exists(plotdir)) dir.create(plotdir)
pop = 'CEU'
if(pop == 'CEU'){
  child = 'GM12878'
  mom = 'GM12892'
  dad = 'GM12891'
}else{
  child = 'GM19240'
  mom = 'GM19238'
  dad = 'GM19239'
}

load('../../rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData') # Load SNP positions
# Select heterozygous positions in the child, phased and not masked
load(file.path(geno.dir, paste(child, '.snps.RData', sep = '')))
child.het = as.vector(geno.info$mat != geno.info$pat & !geno.info$mask & !geno.info$unphased)
child.mat = as.vector(geno.info$mat)
load(file.path(geno.dir, paste(mom, '.snps.RData', sep = '')))
mat.hom = as.vector(geno.info$mat == geno.info$pat & !geno.info$mask)
load(file.path(geno.dir, paste(dad, '.snps.RData', sep = '')))
pat.hom = as.vector(geno.info$mat == geno.info$pat & !geno.info$mask)
pass = child.het & mat.hom & pat.hom & !(as.character(snp.pos$chr) %in% c('chrX', 'chrY'))

filenames = list.files(count.dir, pattern = paste(child, '.*rep\\.RData', sep = ''), full.name = T)
epsilon = 1
nfiles = length(filenames)
rdat = NULL
rdat2 = NULL
for(i in 1:nfiles){
  mat.file = gsub(child, mom, filenames[i])
  pat.file = gsub(child, dad, filenames[i])
  
  if(!file.exists(mat.file) || !file.exists(pat.file)) next
  
  print(basename(filenames[i]))
  load(filenames[i])
  child.ratio = (counts[, 1] + 1) / (counts[, 2] + 1) # ref / alt
  child.ratio[child.mat] = 1 / child.ratio[child.mat] # mat / pat
  child.ratio = log(child.ratio, base = 2)
  pass.c =  as.vector(pass & snp.info$qval < -2) 
  print(sum(pass.c))
  if(sum(pass.c) < 50) next
  
  load(mat.file)
  mat.counts = counts[, 3]
  
  load(pat.file)
  pat.counts = counts[, 3]
  sf = estimateSizeFactorsForMatrix(cbind(mat.counts, pat.counts))
  par.ratio = log((mat.counts / sf[1] + 1) / (pat.counts / sf[2] + 1), base = 2)
  
  child.ratio = child.ratio[pass.c]
  par.ratio = par.ratio[pass.c]
  corr = cor(child.ratio, par.ratio, use = 'na.or.complete')
  agree.frac = sum(child.ratio * par.ratio > 0) / length(child.ratio)
  rtmp = data.frame(child = child.ratio, par = par.ratio, mark1 = sample.info(filenames[i])$mark[1],
                    mark = paste(sample.info(filenames[i])$mark[1], sprintf('%.4f', corr)))
  rtmp2 = data.frame(mark = sample.info(filenames[i])$mark, f = agree.frac, c = corr)
  rdat = rbind(rdat, rtmp)
  rdat2 = rbind(rdat2, rtmp2)
}

rdat$mark1 = order.marks(rdat$mark1)
rdat2$mark = order.marks(rdat2$mark)
uniq.marks = order.marks(levels(rdat$mark1))
uniq.marks.corr = levels(rdat$mark)
tmp = unlist(strsplit(uniq.marks.corr, ' ')) 
uniq.marks.corr = uniq.marks.corr[match(uniq.marks, tmp[seq(1, length(tmp), 2)])]

p = ggplot(rdat) + geom_point(aes(x = child, y = par, color = mark1), shape = 1, size = 1, alpha = 0.6) + theme_bw() + 
  scale_x_continuous('log(maternal/paternal)') + scale_y_continuous('child:log(maternal/paternal)') +
  scale_color_manual(values =  mark.colors(levels(rdat$mark1)), 
                     breaks = uniq.marks, labels = uniq.marks.corr) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 14), panel.grid.major = element_blank())
for(i in 1:nrow(rdat2)){
  p = p + geom_abline(intercept = 0, slope = rdat2$c, color = mark.colors(rdat2$mark), size = 0.3)
}
save(rdat, rdat2, p, file = file.path(plotdir, paste(pop, '.child_par_cor_plot.RData', sep = '')))

rdat = droplevels(rdat[!(as.character(rdat$mark1) %in% c('H2AZ', 'H3K9AC', 'POL4H8')), ])
rdat2 = droplevels(rdat2[!(as.character(rdat2$mark) %in% c('H2AZ', 'H3K9AC', 'POL4H8')), ])
ggsave(file = file.path(plotdir, paste(pop, '.child_par_cor.pdf', sep = '')), p, width =8, height = 5.6)
q = ggplot(rdat2) + geom_bar(aes(x = mark, y = f, fill = mark, width = 0.3), color = 'black', stat = "identity", position = 'identity') + scale_x_discrete('') + theme_bw() + 
  scale_y_continuous('Fraction of AS SNPs') +
  scale_fill_manual(values =  mark.colors(levels(rdat2$mark)), guide = F) +
  theme(axis.text.x = element_text(angle = -65, vjust = 1, hjust = 0, size = 14), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = file.path(plotdir, paste(pop, '.child_par_agreement_new.pdf', sep = '')), q, width = 4, height = 5)

# countdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals')
# plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/plots')
# if(!file.exists(plotdir)) dir.create(plotdir)
# pop = 'YRI'
# load(paste('/media/fusion10/work/chromatinVariation/rawdata/variants/trio/snps/', pop, '.trio.2010_09.genotypes.hg19.r', sep = ''))
# if(pop == 'CEU'){
#   child = 'GM12878'
#   mom = 'GM12892'
#   dad = 'GM12891'
# }else{
#   child = 'GM19240'
#   mom = 'GM19238'
#   dad = 'GM19239'
# }
# 
# child.het = geno.info[[paste(child, '.phased', sep = '')]] & 
#   geno.info[[paste(child, '.mat', sep = '')]] != geno.info[[paste(child, '.pat', sep = '')]] & 
#   !geno.info[[paste(child, '.mask', sep = '')]]
# mat.hom = geno.info[[paste(mom, '.mat', sep = '')]] == geno.info[[paste(mom, '.pat', sep = '')]] & 
#   !geno.info[[paste(mom, '.mask', sep = '')]]
# pat.hom = geno.info[[paste(dad, '.mat', sep = '')]] == geno.info[[paste(dad, '.pat', sep = '')]] & 
#   !geno.info[[paste(dad, '.mask', sep = '')]]
# pass = child.het & mat.hom & pat.hom
# # make sure they are homozygous for different alleles
# stopifnot(all(geno.info[[paste(mom, '.mat', sep = '')]][pass] != geno.info[[paste(dad, '.mat', sep = '')]][pass]))
# 
# filenames <- list.files(countdir, pattern = paste(child, '.*\\.counts\\.RData', sep = ''), full.name = T)
# epsilon <- 1
# nfiles = length(filenames)
# rdat = NULL
# for(i in 1:nfiles){
#   #mat.file = file.path('rawdata/alleleCounts/reps_27Aug12/qvals', gsub('_no1.2_no2.2', '', basename(gsub(child, mom, filenames[i]))))
#   #pat.file = file.path('rawdata/alleleCounts/reps_27Aug12/qvals', gsub('_no1.2_no2.2', '', basename(gsub(child, dad, filenames[i]))))
#   mat.file = gsub(child, mom, filenames[i])
#   pat.file = gsub(child, dad, filenames[i])
#   
#   if(!file.exists(mat.file) || !file.exists(pat.file)) next
#   
#   c = new.env()
#   load(filenames[i], envir = c)
#   print(filenames[i])
#   mat.idx <- grep('mat', colnames(c$counts))
#   pat.idx <- grep('pat', colnames(c$counts))
#   child.mat <- rowSums(c$counts[, mat.idx])
#   child.pat <- rowSums(c$counts[, pat.idx])
#   child.ratio <- log((child.mat + epsilon) / (child.pat + epsilon), base = 2)
#   
#   m = new.env()
#   load(mat.file, m)
#   mat.counts <- rowSums(m$counts) / mean(rowSums(m$counts))
#   p = new.env()
#   load(pat.file, p)
#   pat.counts <- rowSums(p$counts) / mean(rowSums(p$counts))
#   par.ratio <- log((mat.counts + epsilon) / (pat.counts + epsilon), base = 2)
# 
#   pass.c =  pass & !is.nan(c$snp.info$qval) & geno.info$chr != "chrX" & rowSums(c$counts) > 1 & c$snp.info$qval < 0.01 
#   print(sum(pass.c))
#   print(cor(mat.counts[pass.c], child.mat[pass.c]) / cor(pat.counts[pass.c], child.pat[pass.c]))
#   if(sum(pass.c) < 10) next
#   
#   child.ratio = child.ratio[pass.c]
#   par.ratio = par.ratio[pass.c]
#   corr = cor(child.ratio, par.ratio)
#   agree.frac = sum(child.ratio * par.ratio > 0) / length(child.ratio)
#   rtmp = data.frame(child = child.ratio, par = par.ratio, mark = paste(sample.info(filenames[i])$mark[1], sprintf('%.4f', corr)))
#   rtmp2 = data.frame(mark = sample.info(filenames[i])$mark, f = agree.frac)
#   if(is.null(rdat)){ 
#     rdat = rtmp
#     rdat2 = rtmp2
#   }else{ 
#     rdat = rbind(rdat, rtmp)
#     rdat2 = rbind(rdat2, rtmp2)
#   }
# }
# 
# p = ggplot(rdat) + geom_point(aes(x = child, y = par, color = mark), size = 5) + 
#   scale_x_continuous('log(maternal/paternal)') + scale_y_continuous('child:log(maternal/paternal)') +
#   theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
#         axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
#         legend.title=element_blank(), legend.text = element_text(size = 14))
# ggsave(file = file.path(plotdir, paste(pop, '.child_par_cor.png', sep = '')), p)
# q = ggplot(rdat2) + geom_bar(aes(x = mark, y = f), stat = "identity") + scale_x_discrete('') + 
#   scale_y_continuous('Fraction of AS SNPs with same direction in child and parents') +
#   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 10))
# #ggsave(file = file.path(plotdir, paste(pop, '.child_par_agreement.png', sep = ''), q))