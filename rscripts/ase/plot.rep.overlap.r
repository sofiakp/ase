rm(list=ls())
library('reshape')
library('fBasics')
library('VennDiagram')

infile = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts', 'reps_6Aug12', 'all_signif_pos.r')
outdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/', 'reps_6Aug12', 'plots', 'overlaps')
if(!file.exists(outdir)) dir.create(outdir)
load(infile)
types <- unique(samples[samples$rep == 'rep', 1:2])

# Overlap between replicates
for(t in 1:dim(types)[1]){
  reps <- which(samples$indiv == types$indiv[t] & samples$mark == types$mark[t])
  if(length(reps) > 2){
    if(length(reps) > 4) reps = reps[1:4]
    print(samples[reps, ])
    nfiles = length(reps)
    cols = divPalette(nfiles)
    pass.list = list()
    for(f in 1:nfiles){
      pass.list[[paste(samples$mark[reps[f]], samples$rep[reps[f]], sep = '_')]] = hit.list[[reps[f]]]
    }
    outfile <- file.path(outdir, paste(paste(as.matrix(samples[reps[1], 1:2]), collapse = '_'), '_overlaps.tiff', sep = ''))
    venn.diagram(pass.list, outfile, scaled = T, fill = cols, cat.cex = 0.9, margin = 0.2, 
              main = paste(as.matrix(samples[reps[1], 1:2]), collapse = '_'))
  }
}

# Overlap between all "rep" datasets
rep.idx = which(samples$rep == 'rep')
samples = samples[rep.idx, ]
nsamples = dim(samples)[1]
names = paste(samples[, 1], samples[, 2], sep = '_')
ov.mat = array(0, dim = c(nsamples, nsamples))
tot = array(0, dim = c(nsamples, 1))
for(i in 1:(nsamples - 1)){
  tot[i] = length(hit.list[[rep.idx[i]]])
  for(j in (i + 1):nsamples){
    ov.mat[i, j] = sum(hit.list[[rep.idx[i]]] %in% hit.list[[rep.idx[j]]])
  }
}
tot[nsamples] = length(hit.list[[rep.idx[nsamples]]])
ov.mat = ov.mat + t(ov.mat)
ov.mat = apply(ov.mat, 2, function(x) x / tot)
ov.dat = data.frame(ov.mat)
colnames(ov.dat) = names
ov.dat$sample = names
ov.dat = melt(ov.dat)
p = ggplot(ov.dat) + geom_tile(aes(x = sample, y = variable, fill = value)) +
  scale_x_discrete('') + scale_y_discrete('') +
  opts(axis.text.x = theme_text(angle = -45, vjust = 1, hjust = 0, colour = 'grey30'),
       axis.text.y = theme_text(colour = 'grey30')) 

