rm(list=ls())
library('ggplot2')
library('reshape')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

filenames = list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/reps_15Sep12/qvals/'),
                       pattern = 'SNYDER_HG19_GM12878_.*_rep.counts.RData', full.name = T)
samples = sample.info(filenames, '.counts.RData')
nfiles = length(filenames)
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata', 'alleleCounts', 'reps_15Sep12/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

#marks = as.character(unique(samples$mark))
depth.dat = NULL
wilcox = array(NaN, dim = c(nfiles, 1))
ratios = array(NaN, dim = c(nfiles, 1))
for(f in 1:nfiles){
  load(filenames[f])  
  depths = rowSums(counts)
  signif = which(!is.nan(snp.info$pval) & snp.info$qval < 0.01)
  if(length(signif) < 10) next
  pass =  !is.nan(snp.info$pval) & depths >= min(depths[signif])
  true.pvals = snp.info$pval[signif]
  pass.depths = depths[pass]
  min.pvals = binom.val.par(rep(0, sum(pass)), pass.depths) # minimum p-values that the passing sites could possibly give
  
  match.depths = array(NaN, dim = c(length(signif), 1))
  for(i in 1:length(signif)){
    cand = which(min.pvals <= true.pvals[i])
    if(length(cand) > 1) match.depths[i] = pass.depths[cand[sample(1:length(cand), 1)]]
  }
  #good = !is.nan(match.depths)
  #match.depths = match.depths[good]
  #sdepths = depths[signif][good]
  wilcox[f] = wilcox.test(match.depths,depths[signif])$p.value
  ratios[f] = mean(depths[signif], na.rm = T) / mean(match.depths, na.rm = T)
  tmp.dat = data.frame(indiv = rep(samples$indiv[f], length(signif)), mark = rep(samples$mark[f], length(signif)), 
                       AS = depths[signif], matched = match.depths, pval = rep(wilcox[f], length(signif)), ratios = rep(ratios[f], length(signif)))
  if(is.null(depth.dat)){
    depth.dat = tmp.dat
  }else{depth.dat = rbind(depth.dat, tmp.dat)}
}
depth.dat = melt(depth.dat[,1:4], id.vars = c('indiv', 'mark'))
p = ggplot(depth.dat) + geom_density(aes(x = log(value), y = ..density.., color = variable)) + facet_wrap(~mark) + 
  scale_x_continuous('log(depth) at SNP') + scale_color_discrete('')
ggsave(file.path(plotdir, 'depth_distr_at_ASsnps_GM12878.png'), width = 13.6, height = 11.8)