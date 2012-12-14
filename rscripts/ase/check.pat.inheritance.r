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
mark = 'H3K27ME3'
c = new.env()
eps = 0.00001
load(paste('rawdata/geneCounts/reps/norm_qvals/SNYDER_HG19', child, mark, 'rep.genecounts.RData', sep = '_'), c)
c.mat = c$norm.counts$mat + eps
c.pat = c$norm.counts$pat + eps
c.all = rowSums(c$norm.counts)
load(paste('rawdata/geneCounts/reps/norm_qvals/SNYDER_HG19', mom, mark, 'rep.genecounts.RData', sep = '_'))
mom.counts = rowSums(norm.counts)
load(paste('rawdata/geneCounts/reps/norm_qvals/SNYDER_HG19', dad, mark, 'rep.genecounts.RData', sep = '_'))
dad.counts = rowSums(norm.counts)
good = c$gene.info$qval < 0.1 & c$gene.info$chr != "chrX" & c$gene.info$chr != "chrY" #& c$gene.info$qval < 0.01
cor(c.all[good], mom.counts[good])
cor(c.all[good], dad.counts[good])
#good = good & mom.counts > 10
cor(c.mat[good], mom.counts[good])
#good = good & dad.counts > 10
cor(c.pat[good], dad.counts[good])