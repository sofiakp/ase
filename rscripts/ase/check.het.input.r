rm(list=ls())
require('ggplot2')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

is.genecounts = T
if(is.genecounts){suf = 'genecounts'} else{suf = 'counts'}
indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rna/reps_15Sep12')
repfiles = list.files(file.path(indir, 'qvals'), pattern = 'SNYDER_HG19_.*rep.*\\.RData', full.name = T, recursive = F, include.dirs = F)
rep.info = sample.info(repfiles, paste('\\.', suf, '\\.RData$', sep = ''))
filenames <- list.files(file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rna'), 
                        pattern = paste(suf, '\\.RData', sep = ''), full.name = T, recursive = F, include.dirs = F)

all.info <- sample.info(filenames, paste('\\.', suf, '\\.RData$', sep = ''))

outdir = file.path(indir, 'qvals', 'repDisag')
plotdir = file.path(indir, 'plots')
if(!file.exists(outdir)) dir.create(outdir)
if(!file.exists(plotdir)) dir.create(plotdir)

types <- unique(rep.info[, 1:2])
rdat = NULL

# For each combination of individual and mark...
for(t in 1:dim(types)[1]){
  reps <- which(as.character(all.info$indiv) == as.character(types$indiv[t]) & 
    as.character(all.info$mark) == as.character(types$mark[t])) 
  samples = all.info[reps, ]
  outfile = file.path(outdir, paste('SNYDER_HG19', types$indiv[t], types$mark[t], 'rep_disag.RData', sep = '_'))
  sel.rep = which(rep.info$indiv == types$indiv[t] & rep.info$mark == types$mark[t])
  
  # Find all the datasets for this individual and mark...
  if(length(reps) > 1 && length(sel.rep) == 1){
    if(file.exists(outfile)){
      load(outfile)
    }else{
      
      print(repfiles[sel.rep])
      load(repfiles[sel.rep])
      if(is.genecounts){
        snp.info = gene.info
        rep.ratio = log((counts$mat + 1) / (counts$pat + 1))
      }else{
        rep.ratio = log((rowSums(counts[, grep('mat', colnames(counts))]) + 1) / (rowSums(counts[, grep('pat', colnames(counts))]) + 1))
      }
      good = !is.nan(snp.info$qval) & snp.info$chr != "chrX" & snp.info$chr != "chrY" # not masked, het
      pass = good & snp.info$qval < 0.01 # not masked, het, and AS
      npass = sum(pass)
      ngood = sum(good)
      
      # Fraction of SNPs in each replicate that have different AS direction from merged
      # and the number of SNPs considered in this comparison.
      frac.diff = array(0, dim = c(length(reps), 2)) 
      #ratios.good = array(0, dim = c(ngood, length(reps))) 
      #ratios.pass = array(0, dim = c(npass, length(reps)))
      for(i in 1:length(reps)){
        print(filenames[reps[i]])
        load(filenames[reps[i]])
        if(is.genecounts){
          r = log((counts$mat + 1) / (counts$pat + 1))
        }else{
          r = log((rowSums(counts[, grep('mat', colnames(counts))]) + 1) / (rowSums(counts[, grep('pat', colnames(counts))]) + 1))
        }
        considered = pass & rowSums(counts) >= 10
        frac.diff[i, 1] = sum(rep.ratio * r <= 0 & considered) / sum(considered) 
        frac.diff[i, 2] = sum(considered)
        #ratios.good[, i] = r[good]
        #ratios.pass[, i] = r[pass]
      }
      
      if(!is.genecounts) save(frac.diff, samples, file = outfile)
      #     npairs = 0
      #     rep.cor = array(0, dim = c(length(reps) * (length(reps) - 1) / 2, 2)) # AS correlation between replicates
      #     rep.names = array('', dim = c(length(reps) * (length(reps) - 1) / 2, 1))
      #     for(i in 1:(length(reps) - 1)){
      #       for(j in (i + 1):length(reps)){
      #         npairs = npairs + 1
      #         rep.cor[npairs, 1] = cor(ratios.good[, i], ratios.good[, j])
      #         rep.cor[npairs, 2] = cor(ratios.pass[, i], ratios.pass[, j])
      #         rep.names[npairs] = paste(all.info$rep[reps[i]], all.info$rep[reps[j]], sep = '_vs_')
      #       }
      #     }
    }
    tmp.dat = data.frame(indiv = all.info$indiv[reps], mark = all.info$mark[reps], rep = all.info$rep[reps], f = frac.diff[, 1], n = frac.diff[, 2])
    #tmp.dat2 = data.frame(indiv = rep(types$indiv[t], npairs), mark = rep(types$mark[t], npairs), r.good = rep.cor[, 1], r.pass = rep.cor[, 2], name = rep.names)
    if(is.null(rdat)){
      rdat = tmp.dat
      #rdat2 = tmp.dat2
    }else{
      rdat = rbind(rdat, tmp.dat)
      #rdat2 = rbind(rdat2, tmp.dat2)
    }
  }
}

p <- ggplot(rdat) + geom_bar(aes(x = indiv, y = 100*f, fill = rep), stat = "identity", position = "dodge") +
  facet_wrap(~mark) + scale_fill_discrete('Replicate') + scale_x_discrete("") + 
  scale_y_continuous("% of AS SNPs that disagree with merged replicates") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file.path(plotdir, 'rep_disagreement.png'), width = 14.8, height = 10.6)

# p <- ggplot(rdat2) + geom_bar(aes(x = name, y = r.pass), stat = "identity", position = "dodge") +
#     facet_grid(indiv~mark, scales = "free_x") + scale_fill_discrete('Replicate') + scale_x_discrete("") + 
#     scale_y_continuous("Correlation of log(mat/pat) between replicates") +
#     opts(axis.text.x = theme_text(angle = -70, vjust = 1, hjust = 0, size = 6))
# ggsave(file.path(outdir, 'rep_signif_corr.png'), width = 14.8, height = 10.6)
# save(rdat, rdat2, file = file.path(outdir, 'rep_corr.RData'))
# rep1 = new.env()
# load('rawdata/alleleCounts/SNYDER_HG19_GM12878_H3K27AC_1.counts.r', rep1)
# rep2 = new.env()
# load('rawdata/alleleCounts/SNYDER_HG19_GM12878_H3K27AC_2.counts.r', rep2)
# rep3 = new.env()
# load('rawdata/alleleCounts/SNYDER_HG19_GM12878_H3K27AC_3.counts.r', rep3)
# rep0= new.env()
# load('rawdata/alleleCounts/reps_6Aug12/SNYDER_HG19_GM12878_H3K27AC_rep.counts.r', rep0)
# mat.idx = grep('mat', colnames(rep0$counts))
# pat.idx = grep('pat', colnames(rep0$counts))
# pass = rep0$snp.info$qval < 0.01 & rep0$snp.info$sb > 0.1
# mat = rowSums(rep0$counts[, mat.idx]) > rowSums(rep0$counts[, pat.idx])
