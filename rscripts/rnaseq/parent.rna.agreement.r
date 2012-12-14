rm(list=ls())
require(ggplot2)
library('DESeq')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

countdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps/qvals')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps/plots')
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

filenames <- list.files(countdir, pattern = paste(child, '.*\\.genecounts\\.RData', sep = ''), full.name = T)
epsilon <- 1
nfiles = length(filenames)
rdat = NULL
for(i in 1:nfiles){
  mat.file = gsub(child, mom, filenames[i])
  pat.file = gsub(child, dad, filenames[i])
  
  if(!file.exists(mat.file) || !file.exists(pat.file)) next
  
  c = new.env()
  load(filenames[i], envir = c)
  mat.idx <- grep('mat', colnames(c$norm.counts))
  pat.idx <- grep('pat', colnames(c$norm.counts))
  good = rowSums(c$norm.counts) > 10
  child.genes = c$gene.info[good, ]  
  child.mat <- c$norm.counts[good, mat.idx]
  child.pat <- c$norm.counts[good, pat.idx]
  child.ratio <- log((child.mat + epsilon) / (child.pat + epsilon), base = 2)
  
  m = new.env()
  load(mat.file, m)
  mat.counts <- rowSums(m$norm.counts)[good]
  p = new.env()
  load(pat.file, p)
  pat.counts <- rowSums(p$norm.counts)[good]
  #conds = c('mat', 'pat')
  #cds = newCountDataSet(cbind(mat.counts, pat.counts), conds)
  #cds = estimateSizeFactors(cds)
  par.ratio <- log((mat.counts + epsilon) / (pat.counts + epsilon), base = 2)
  
  pass.c =  !is.nan(child.genes$qval) & child.genes$qval < 0.01 & child.genes$chr != "chrX" & child.genes$het #& !m$gene.info$het
  print(sum(pass.c))
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

q1 = ggplot(rdat) + geom_point(aes(x = child, y = par, color = mark)) + scale_x_continuous('log(maternal/paternal)') + 
  scale_y_continuous('child:log(maternal/paternal)') + coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5))
#ggsave(file = file.path(plotdir, paste(pop, '.child_par_cor.png', sep = '')))
q2 = ggplot(rdat2) + geom_bar(aes(x = mark, y = f), stat = "identity") + scale_x_discrete('') + 
  scale_y_continuous('Fraction of AS SNPs with same direction in child and parents') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
#ggsave(file = file.path(plotdir, paste(pop, '.child_par_agreement.png', sep = '')))