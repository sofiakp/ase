rm(list=ls())

indir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/rdata'
outdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/plots/qn_isvaNull_fits_all_reg/'
if(!file.exists(outdir)) dir.create(outdir)

files = list.files(indir, pattern = 'SNYDER_HG19_all_reg_.*_qn_isvaNull_clust.RData', full.names = T)

oldsc = options(scipen = 100)
for(f in files){
  load(gsub('_clust.RData', '.RData', f))
  load(f)
  # Write NON-population specific regions
  reg = regions[-isva.fit$deg, 1:3]
  reg$start = reg$start - 1 # convert back to 0-based
  outfile = file.path(outdir, gsub('_clust.RData', '_non_sign.bed', basename(f)))
  write.table(reg, quote = F, row.names = F, col.names = F, sep = '\t', file = outfile)
  # Write population specific regions
  reg = regions[isva.fit$deg, 1:3]
  reg$start = reg$start - 1 # convert back to 0-based
  outfile = file.path(outdir, gsub('_clust.RData', '_sign.bed', basename(f)))
  write.table(reg, quote = F, row.names = F, col.names = F, sep = '\t', file = outfile)
  # Write the regions that were used for plotting
  reg = reg[plot.rows, ]
  outfile = file.path(outdir, gsub('_clust.RData', '_sign_clust.bed', basename(f)))
  write.table(reg, quote = F, row.names = F, col.names = F, sep = '\t', file = outfile)
}
options(oldsc)