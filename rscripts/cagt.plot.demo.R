rm(list=ls())
source('/home/sofiakp/Documents/MATLAB/cagt/matlab/scripts/plot.cagt.multivariate.figs.R')

input.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt')
output.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals/hitLists/cagt/plots')
if(!file.exists(output.dir)) dir.create(output.dir)
all.cagt.table.filenames = list.files(path=input.dir, pattern = "SNYDER_HG19_GM12878_H3K27AC_AS_AT_SNYDER_HG19_GM12878_.*multi.txt", full.names = T) # Get names of cagt plot names

for (i in all.cagt.table.filenames) {
  cat(i,"\n")
  plot.cagt.multivariate.figs(cagt.table.filename = i,
                              output.dir = output.dir,
                              output.file = NULL, # Let it name the figures
                              output.format = "pdf",
                              support.thresh = 0.0, # remove shapes whose support is < .02 * fraction of peaks in high signal component
                              orient.shapes = T, # If required invert shapes to make them Low (left) to High (right)                                        
                              replace.flag = T, 
                              partner.marks = 'SIGNAL', # Only plot the signal used for clustering without partner marks
                              zscore = F, # Use the original scale and not the standardized signal.
                              common.scale = F, # All plots on a common scale
                              plot.over = T, # plot them all on the same plot or one below the other
                              whiskers = T, # Plot percentile whiskers around the cluster signal.
                              mag.shape = "shape" # Plot both the clusters and the averaged aggregation plots. 
                              )
#   plot.cagt.multivariate.figs(cagt.table.filename = i,
#                               output.dir = output.dir,
#                               output.format = "png",
#                               output.file = NULL,
#                               support.thresh = 0.0,
#                               orient.shapes = T, 
#                               replace.flag = T,
#                               partner.marks = 'all', 
#                               zscore = F,
#                               common.scale = F, 
#                               plot.over = F,
#                               whiskers = F,
#                               mag.shape = "all"
#                               )
}