rm(list=ls())
require(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

countdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/plots')
if(!file.exists(plotdir)) dir.create(plotdir)
pop = 'YRI'
load(paste('/media/fusion10/work/chromatinVariation/rawdata/variants/trio/snps/', pop, '.trio.2010_09.genotypes.hg19.r', sep = ''))
if(pop == 'CEU'){
  child = 'GM12878'
  mom = 'GM12892'
  dad = 'GM12891'
}else{
  child = 'GM19240'
  mom = 'GM19238'
  dad = 'GM19239'
}

child.het = geno.info[[paste(child, '.phased', sep = '')]] & 
  geno.info[[paste(child, '.mat', sep = '')]] != geno.info[[paste(child, '.pat', sep = '')]] & 
  !geno.info[[paste(child, '.mask', sep = '')]]
mat.hom = geno.info[[paste(mom, '.mat', sep = '')]] == geno.info[[paste(mom, '.pat', sep = '')]] & 
  !geno.info[[paste(mom, '.mask', sep = '')]]
pat.hom = geno.info[[paste(dad, '.mat', sep = '')]] == geno.info[[paste(dad, '.pat', sep = '')]] & 
  !geno.info[[paste(dad, '.mask', sep = '')]]
pass = child.het & mat.hom & pat.hom
# make sure they are homozygous for different alleles
stopifnot(all(geno.info[[paste(mom, '.mat', sep = '')]][pass] != geno.info[[paste(dad, '.mat', sep = '')]][pass]))

filenames <- list.files(countdir, pattern = paste(child, '.*\\.counts\\.RData', sep = ''), full.name = T)
epsilon <- 1
nfiles = length(filenames)
rdat = NULL
for(i in 1:nfiles){
  #mat.file = file.path('rawdata/alleleCounts/reps_27Aug12/qvals', gsub('_no1.2_no2.2', '', basename(gsub(child, mom, filenames[i]))))
  #pat.file = file.path('rawdata/alleleCounts/reps_27Aug12/qvals', gsub('_no1.2_no2.2', '', basename(gsub(child, dad, filenames[i]))))
  mat.file = gsub(child, mom, filenames[i])
  pat.file = gsub(child, dad, filenames[i])
  
  if(!file.exists(mat.file) || !file.exists(pat.file)) next
  
  c = new.env()
  load(filenames[i], envir = c)
  print(filenames[i])
  mat.idx <- grep('mat', colnames(c$counts))
  pat.idx <- grep('pat', colnames(c$counts))
  child.mat <- rowSums(c$counts[, mat.idx])
  child.pat <- rowSums(c$counts[, pat.idx])
  child.ratio <- log((child.mat + epsilon) / (child.pat + epsilon), base = 2)
  
  m = new.env()
  load(mat.file, m)
  mat.counts <- rowSums(m$counts) / mean(rowSums(m$counts))
  p = new.env()
  load(pat.file, p)
  pat.counts <- rowSums(p$counts) / mean(rowSums(p$counts))
  par.ratio <- log((mat.counts + epsilon) / (pat.counts + epsilon), base = 2)

  pass.c =  pass & !is.nan(c$snp.info$qval) & geno.info$chr != "chrX" & rowSums(c$counts) > 1 & c$snp.info$qval < 0.01 
  print(sum(pass.c))
  print(cor(mat.counts[pass.c], child.mat[pass.c]) / cor(pat.counts[pass.c], child.pat[pass.c]))
  if(sum(pass.c) < 10) next
  
  child.ratio = child.ratio[pass.c]
  par.ratio = par.ratio[pass.c]
  corr = cor(child.ratio, par.ratio)
  agree.frac = sum(child.ratio * par.ratio > 0) / length(child.ratio)
  rtmp = data.frame(child = child.ratio, par = par.ratio, mark = paste(sample.info(filenames[i])$mark[1], sprintf('%.4f', corr)))
  rtmp2 = data.frame(mark = sample.info(filenames[i])$mark, f = agree.frac)
  if(is.null(rdat)){ 
    rdat = rtmp
    rdat2 = rtmp2
  }else{ 
    rdat = rbind(rdat, rtmp)
    rdat2 = rbind(rdat2, rtmp2)
  }
}

p = ggplot(rdat) + geom_point(aes(x = child, y = par, color = mark), size = 5) + 
  scale_x_continuous('log(maternal/paternal)') + scale_y_continuous('child:log(maternal/paternal)') +
  theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title=element_blank(), legend.text = element_text(size = 14))
ggsave(file = file.path(plotdir, paste(pop, '.child_par_cor.png', sep = '')), p)
q = ggplot(rdat2) + geom_bar(aes(x = mark, y = f), stat = "identity") + scale_x_discrete('') + 
  scale_y_continuous('Fraction of AS SNPs with same direction in child and parents') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 10))
#ggsave(file = file.path(plotdir, paste(pop, '.child_par_agreement.png', sep = ''), q))