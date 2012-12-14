rm(list=ls())
library(reshape)
library(ggplot2)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata', 'reps')
countdir = file.path(indir, 'qvals')
plotdir = file.path(indir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = ''

geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan') # genotype RData files should be here

read.lists = T # T: input is assumed to be indicators of significant positions, F: input is assumed to be data.frame with signif positions
if(read.lists){
  filenames = list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.RData', full.name = T)
  samples = sample.info(filenames, '\\.RData$')  
}else{
  filenames = list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.het\\.RData', full.name = T)
  samples = sample.info(filenames, '\\.het\\.RData$')  
}
nsamples = length(filenames)

indivs = unique(as.character(samples$indiv))
all.hits = list()
nhits = array(NaN, dim = c(nsamples, 1))
ncand = array(NaN, dim = c(nsamples, 1))

if(read.lists){
  for(i in 1:nsamples){
    load(filenames[i])
    pass = as.vector(snp.info$pass & snp.info$enriched & counts[, 3] > 7)
    hits = snp.info$qval < -2
    nhits[i] = sum(hits)
    ncand[i] = sum(hits & pass) * 100 / sum(pass)
    indiv = as.character(samples$indiv[i])
    if(is.null(all.hits[[indiv]])){
      all.hits[[indiv]] = hits
    }else{
      all.hits[[indiv]] = all.hits[[indiv]] | hits
    }
  }  
}else{
  for(i in 1:nsamples){
    load(filenames[i])
    tot.counts = het.counts$ref.counts + het.counts$alt.counts
    pass = het.counts$enriched & tot.counts > 7 # het.counts are already restricted to heterozygous and unmasked
    nhits[i] = sum(het.counts$signif)
    ncand[i] = sum(het.counts$signif & pass) * 100 / sum(pass)
    hits1 = sort(paste(het.counts$chr[het.counts$signif], het.counts$pos[het.counts$signif], sep = ':'), index.return = T)
    indiv = as.character(samples$indiv[i])
    if(is.null(all.hits[[indiv]])){
      all.hits[[indiv]] = hits1$x
    }else{
      all.hits[[indiv]] = union(all.hits[[indiv]], hits1$x)
    }
  }  
}

out.size = 10
if(nsamples > 50){
  out.size = 5
}

hit.dat = data.frame(indiv = samples$indiv, mark = samples$mark, hits = nhits, norm.hits = ncand)
p = ggplot(hit.dat) + geom_bar(aes(x = indiv, y = norm.hits), stat = "identity") + 
  facet_wrap(~mark, scales = "free_y") + scale_y_continuous('AS-SNPs / het SNPs (%)') + scale_x_discrete('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file.path(plotdir, 'num_as_snps_norm.png'))

if(read.lists){
  for(d in indivs){
    hit.dat = rbind(hit.dat, data.frame(indiv = factor(d), mark = 'union', hits = sum(all.hits[[d]]), norm.hits = 1))
  }
}else{
  for(d in indivs){
    hit.dat = rbind(hit.dat, data.frame(indiv = factor(d), mark = 'union', hits = length(all.hits[[d]]), norm.hits = 1))
  }  
}

q = ggplot(hit.dat) + geom_bar(aes(x = indiv, y = hits), stat = "identity") + 
  facet_wrap(~mark, scales = "free_y") + scale_y_continuous('# AS-SNPs') + scale_x_discrete('') +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave(file.path(plotdir, 'num_as_snps.png'))
