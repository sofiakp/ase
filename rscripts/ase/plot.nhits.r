rm(list=ls())
library(reshape)
library(ggplot2)
source('utils/sample.info.r')
source('utils/deseq.utils.r')

indir = '../../rawdata/alleleCounts/san/rdata/reps/'
countdir = file.path(indir, 'qvals')
plotdir = file.path(indir, 'plots')
if(!file.exists(plotdir)) dir.create(plotdir)
outpref = ''

#geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan') # genotype RData files should be here

read.lists = T # T: input is assumed to be indicators of significant positions, F: input is assumed to be data.frame with signif positions
if(read.lists){
  filenames = list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.RData', full.name = T)
  filenames = append(filenames, list.files('../../rawdata/alleleCounts/allNonSan/rdata/reps/qvals/', pattern = 'SNYDER_HG19.*rep\\.RData', full.name = T))
  samples = sample.info(filenames, '\\.RData$')  
  filenames = filenames[!(samples$mark %in% c('BUB', 'EPOL', 'POL4H8', 'H2AZ', 'RZ', 'H3K9AC', 'H3K9ME3', 'polyA-RNA', 'PU1'))]
  samples = sample.info(filenames, '\\.RData$')  
}else{
  filenames = list.files(countdir, pattern = 'SNYDER_HG19.*rep\\.het\\.RData', full.name = T)
  
  filenames = filenames[!grepl('BUB', filenames)]
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

new.indiv = as.character(fix.indiv.names(samples$indiv))
hit.dat = data.frame(indiv = new.indiv, mark = order.marks(samples$mark), hits = nhits, norm.hits = ncand)
#write.table(hit.dat, file = file.path(plotdir, 'num_as_snps_stats.txt'), row.names = F, col.names = T, sep = '\t', quote = F)

p1 = ggplot(hit.dat) + geom_bar(aes(x = indiv, y = norm.hits, fill = mark), stat = "identity") + 
  facet_wrap(~mark, scales = "free_y") + scale_y_continuous('AS-SNPs / het SNPs (%)') + scale_x_discrete('') + theme_bw() +
  scale_fill_manual(values =  mark.colors(hit.dat$mark), guide = F) +
  theme(axis.text.x = element_text(size = 15, angle = 50, vjust = 1, hjust = 1), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
#ggsave(file.path(plotdir, 'num_as_snps_norm.pdf'), p1, width = 13.6, height = 11.8)
p2 = ggplot(hit.dat) + geom_boxplot(aes(x = mark, y = norm.hits, fill = mark)) + 
  xlab('') + ylab('AS-SNPs / het SNPs (%)') + theme_bw() + 
  scale_fill_manual(values =  mark.colors(hit.dat$mark), guide = F) +
  theme(axis.text.x = element_text(size = 18, angle = 50, vjust = 1, hjust = 1), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, 'num_as_snps_norm_box_new.pdf'), p2, width = 6.5, height = 5.6)

# if(read.lists){
#   for(d in new.indiv){
#     hit.dat = rbind(hit.dat, data.frame(indiv = factor(d), mark = 'union', hits = sum(all.hits[[d]]), norm.hits = 1))
#   }
# }else{
#   for(d in new.indiv){
#     hit.dat = rbind(hit.dat, data.frame(indiv = factor(d), mark = 'union', hits = length(all.hits[[d]]), norm.hits = 1))
#   }  
# }
# 
# lev = levels(hit.dat$mark)
# q1 = ggplot(hit.dat) + geom_bar(aes(x = indiv, y = hits, fill = mark), stat = "identity") + 
#   facet_wrap(~mark, scales = "free_y") + scale_y_continuous('# AS-SNPs') + scale_x_discrete('') + theme_bw() +
#   scale_fill_manual(values =  mark.colors(hit.dat$mark), guide = F) +
#   theme(axis.text.x = element_text(size = 15, angle = -65, vjust = 1, hjust = 0), axis.title.x = element_text(size = 16),
#         axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 16),
#         strip.text.x = element_text(size = 16))
# #ggsave(file.path(plotdir, 'num_as_snps.pdf'), q1, width = 13.6, height = 11.8)
# q2 = ggplot(hit.dat) + geom_boxplot(aes(x = mark, y = hits, fill = mark)) + 
#   xlab('') + ylab('# AS-SNPs') + theme_bw() + 
#   scale_fill_manual(values = mark.colors(hit.dat$mark), guide = F) +
#   theme(axis.text.x = element_text(size = 15, angle = -65, vjust = 1, hjust = 0), axis.title.x = element_text(size = 16),
#         axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 16),
#         strip.text.x = element_text(size = 16))
# #ggsave(file.path(plotdir, 'num_as_snps_box.pdf'), q2, width = 6.5, height = 5.6)