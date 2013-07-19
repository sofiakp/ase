rm(list=ls())
library(ggplot2)
library(Matrix)
library(GenomicRanges)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))

get.enrich.corr = function(pos.regions, neg.regions, dat, ov.deg, ov.non.deg, is.snp = T){
  enrich = array(NaN, dim = c(1, 2))
  if(is.snp){
    deg = findOverlaps(regions.to.ranges(pos.regions), snps.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(neg.regions), snps.to.ranges(dat), ignore.strand = T, select = 'first')
  }else{
    deg = findOverlaps(regions.to.ranges(pos.regions), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(neg.regions), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  }
  deg.sum = sum(!is.na(ov.deg) & !is.na(deg)) # Regions with variants (in our set of samples) AND high Fst SNPs
  non.deg.sum = sum(!is.na(ov.non.deg) & !is.na(non.deg))
  enrich[1, 1] = (deg.sum / ov.deg.sum) / (non.deg.sum / ov.non.deg.sum)
  enrich[1, 2] = binom.val(deg.sum, ov.deg.sum, non.deg.sum / ov.non.deg.sum, alternative = 'greater')
  return(enrich)
}

get.enrich = function(pos.regions, neg.regions, dat, is.snp = T){
  enrich = array(NaN, dim = c(1, 2))
  if(is.snp){
    deg = findOverlaps(regions.to.ranges(pos.regions), snps.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(neg.regions), snps.to.ranges(dat), ignore.strand = T, select = 'first')
  }else{
    deg = findOverlaps(regions.to.ranges(pos.regions), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
    non.deg = findOverlaps(regions.to.ranges(neg.regions), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  }
  deg.sum = sum(!is.na(deg))
  non.deg.sum = sum(!is.na(non.deg))
  enrich[1, 1] = (deg.sum / nrow(pos.regions)) / (non.deg.sum / nrow(neg.regions))
  enrich[1, 2] = binom.val(deg.sum, nrow(pos.regions), non.deg.sum / nrow(neg.regions), alternative = 'greater')
  return(enrich)
}
get.wilc.enrich = function(regions, isva.fit, dat){
  deg = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')    
  non.deg = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), regions.to.ranges(dat), ignore.strand = T, select = 'first')
  values.deg = dat[deg[!is.na(deg)], 4]
  values.non.deg = dat[non.deg[!is.na(non.deg)], 4]
  enrich = c(mean(values.deg) / mean(values.non.deg), wilcox.test(values.deg, values.non.deg)$p.value)
  return(enrich)
}

is.rand = F
mark = '' # only used in rand
files = c()
isva.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata' #'../../rawdata/segSignal/14indiv/extractSignal/fc/avgSig/rdata'
plotdir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/plots' # CHANGE THIS!
files = append(files, list.files(isva.dir, pattern = paste('all_reg_', mark, '.*comp3_q0.01.*_qn_isvaNull.RData', sep = ''), full.names = T))
files = files[!grepl('daughters', files)]
outpref = 'SNYDER_HG19_all_reg_'
if(is.rand){
  isva.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/perm/'
  files = append(files, list.files(isva.dir, pattern = paste('all_reg_', mark, '.*comp3_q0.01.*_qn_isvaNull.RData', sep = ''), full.names = T))
  isva.dir = '../../rawdata/signal/combrep/extractSignal/rand/fc/avgSig/merged_Mar13/rdata/'
  files = append(files, list.files(isva.dir, pattern = paste('all_reg_rand_', mark, '.*comp3_q0.01.*_qn_isvaNull.RData', sep = ''), full.names = T))
  outpref = paste(outpref, mark, '_rand_', sep = '')
}else{
  isva.dir = '../../rawdata/genomeGrid/hg19_w10k/combrep/fc/avgSig_newNorm/rdata/'
  files = append(files, list.files(isva.dir, pattern = 'all_reg.*_H3K27ME3.*comp3_q0.01.*_qn_isvaNull.RData', full.names = T))
  isva.dir = '../../rawdata/geneCounts/rdata/repsComb/rdata/'
  files = append(files, list.files(isva.dir, pattern = '.*_RZ.*comp3_q0.01.*_qn_isvaNull.RData', full.names = T))
  isva.dir = '../../rawdata/transcriptomes/combrep/extractSignal/fc/avgSig_newNorm/rdata/'
  files = append(files, list.files(isva.dir, pattern = 'all_reg.*_(H3K36ME3).*comp3_q0.01.*_qn_isvaNull.RData', full.names = T))
}

#geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/')

if(!is.rand) files = files[!grepl('rand|th0', files)]

# Load SNPs and genotype PCA
load('../../rawdata/variants/all_Mar13/snps.RData')
#gp = new.env()
#load('../../rawdata/variants/all/snps/allNonSan/genot_pca.RData', gp)
#snp.pos = snp.pos[gp$pca.rows, 1:2] # SNPs that are variant

# Load gene metadata
load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')

#ihs.dat = read.table('../../rawdata/selection/iHS_hg19.bed', header = F, sep = '\t')
#colnames(ihs.dat) = c('chr', 'start', 'end', 'value')
#ihs.dat$pos = ihs.dat$end

# Read SNPs with high Fst
fst.dat = read.table('../../rawdata/variants/fst/FstSNP-HapMap3_selPop.txt')[, 1:3]
colnames(fst.dat) = c('chr', 'pos', 'fst')

# Read eQTLs
eqtl.dat = read.table('../../rawdata/QTL/eQTLs_uniq.txt')[, 1:2]
colnames(eqtl.dat) = c('chr', 'pos')

# 
good.10k.bins = read.bed('../../rawdata/signal/combrep/peakFiles/merged_Mar13/H3K27me3_10kbwin_overlap_SNYDER_HG19_H3K27ME3_merged.bed.gz')

marks = gsub('SNYDER_HG19_|hg19_w.*k_|allEnhStates_|.*all_reg_|_comp[0-9]*_q0.[0-9]*|_qn_.*svaNull.RData', '', basename(files))
nmarks = length(marks)
tot.len = array(0, dim = c(nmarks, 1))
npeaks = array(0, dim = c(nmarks, 2))
snp.enrich = array(1, dim = c(nmarks, 2))
snp.enrich[, 1] = 0
fenrich = array(1, dim = c(nmarks, 2))
fenrich[, 1] = 0
qtl.enrich = array(1, dim = c(nmarks, 2))
qtl.enrich[, 1] = 0
ihs.enrich = array(1, dim = c(nmarks, 2))
ihs.enrich[, 1] = 0
signal.wilc = array(1, dim = c(nmarks, 2))
signal.wilc[, 1] = 0

penrich = NULL
all.deg.snps = NULL
all.non.deg.snps = NULL
for(i in 1:length(files)){
  qn.file = gsub('_comp.*', '_qn.RData', files[i])
  if(is.rand & !file.exists(qn.file)){
    qn.file = gsub('_comp.*|_isvaNull|_rand_pop[0-9]*', '_qn.RData', files[1])
  }
  print(basename(files[i]))
  print(basename(qn.file))
  load(qn.file)
  #orig.counts = counts[good.rows, ]
  
  if(grepl('H3K27ME3', files[i])){
    npeaks[i, 2] = nrow(good.10k.bins)
  }else{
    npeaks[i, 2] = length(good.rows) # Total number of regions, including those that were removed before ISVA 
  }
  load(files[i]) # Load the file with the ISVA results. This only has the regions used for ISVA.
  if(is.null(isva.fit$deg)){
    next
  }
  if(grepl('H3K27ME3', files[i])){
    ov = findOverlaps(regions.to.ranges(regions[isva.fit$deg, ]), regions.to.ranges(good.10k.bins), select = 'first', ignore.strand = T)
    pos.regions = regions[isva.fit$deg, ][!is.na(ov), ]
    ov = findOverlaps(regions.to.ranges(regions[-isva.fit$deg, ]), regions.to.ranges(good.10k.bins), select = 'first', ignore.strand = T)
    neg.regions = regions[-isva.fit$deg, ][!is.na(ov), ]
  }else{
    if(grepl('RZ|RNA', files[i])){
      regions = gene.meta[good.rows, 2:4]
    }
    pos.regions = regions[isva.fit$deg, ]
    neg.regions = regions[-isva.fit$deg, ]
  }
  npos = nrow(pos.regions)
  nneg = nrow(neg.regions)
  npeaks[i, 1] = npos
  tot.len[i] = sum(pos.regions$end - pos.regions$start + 1)
  
  if(npos > 1){
    #signal.wilc[i, 1] = median(rowMeans(orig.counts[isva.fit$deg, ])) / median(rowMeans(orig.counts[-isva.fit$deg, ]))
    #signal.wilc[i, 2] = wilcox.test(rowMeans(orig.counts[isva.fit$deg, ]), rowMeans(orig.counts[-isva.fit$deg, ]))$p.value
    
    ######## Are selected regions enriched in variants?
    # Count how many regions overlap at least one variant
    ov.deg = findOverlaps(regions.to.ranges(pos.regions), snps.to.ranges(snp.pos), select = 'first', ignore.strand = T)
    ov.deg.sum = sum(!is.na(ov.deg))
    ov.non.deg = findOverlaps(regions.to.ranges(neg.regions), snps.to.ranges(snp.pos), select = 'first', ignore.strand = T)
    ov.non.deg.sum = sum(!is.na(ov.non.deg))
    snp.enrich[i, 1] = (ov.deg.sum / npos) / (ov.non.deg.sum / nneg)
    snp.enrich[i, 2] = binom.val(ov.deg.sum, npos, ov.non.deg.sum / nneg, alternative = 'greater')
    
    #ihs.enrich[i, ] = get.wilc.enrich(regions, isva.fit, ihs.dat) #get.enrich(regions, isva.fit, ihs.dat, is.snp = F)
    
    ######## Are the selected regions enriched in high-Fst variants?
    fenrich[i, ] = get.enrich.corr(pos.regions, neg.regions, fst.dat, ov.deg, ov.non.deg, is.snp = T)
    
    ######## Are the selected regions enriched in eQTLs?    
    qtl.enrich[i, ] = get.enrich.corr(pos.regions, neg.regions, eqtl.dat, ov.deg, ov.non.deg, is.snp = T)
  }
}

#frac.peaks = npeaks[, 1] / npeaks[, 2]
sel.marks = npeaks[, 2] > -1 & npeaks[, 1] > -1
marks = marks[sel.marks] #gsub('_H3K4ME1', '', marks[sel.marks])
marks = order.marks(marks)
npeaks = npeaks[sel.marks, ]
tot.len = tot.len[sel.marks]
fenrich = fenrich[sel.marks, ]
snp.enrich = snp.enrich[sel.marks, ]
qtl.enrich = qtl.enrich[sel.marks, ]
len.dat = data.frame(mark = marks, len = tot.len / 1e6, tot = npeaks[,1], npeaks = npeaks[, 1]/npeaks[, 2], f = fenrich[, 1], fp = fenrich[, 2], 
                     p = snp.enrich[, 1], pp = snp.enrich[, 2], q = qtl.enrich[, 1], qp = qtl.enrich[, 2])
save(len.dat, file = file.path(plotdir, paste(outpref, 'enrich.RData', sep = '')))
write.table(len.dat, file = file.path(plotdir, paste(outpref, 'enrich.txt', sep = '')), col.names = T, row.names = F, quote = F, sep = '\t')
p1 = ggplot(len.dat) + geom_bar(aes(x = mark, y = len, fill = mark), color = '#000000', stat = 'identity') + xlab('') + 
  ylab('Length of population specific regions (Mb)') + theme_bw() +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_color_discrete(guide = F) + 
  theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_len.pdf', sep = '')), p1, width = 6.5, height = 5.6)

p2 = ggplot(len.dat) + geom_bar(aes(x = mark, y = npeaks, fill = mark), color = '#000000', stat = 'identity') + 
  annotate('text', x = len.dat$mark, y = len.dat$npeaks + 0.005, label = len.dat$tot, size = 5) + 
  xlab('') + ylab('Fraction of population specific regions') + theme_bw() +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_color_discrete(guide = F) +
  theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_frac.pdf', sep = '')), p2, width = 6.5, height = 5.6)


# annot = array('', dim = c(nrow(len.dat), 1))
# annot[fenrich[, 2] < 0.05 & npeaks[, 1] > 0] = '*'
# annot.y = len.dat$f
# tmp.len.dat = len.dat[, c(1,5,6)]
# colnames(tmp.len.dat) = c('mark', 'enrich', 'pval')
# tmp.len.dat$type= rep('Fst', nrow(len.dat))
# tmp.len.dat.2 = len.dat[, c(1,9,10)]
# colnames(tmp.len.dat.2) = c('mark', 'enrich', 'pval')
# tmp.len.dat.2$type= rep('eQTL', nrow(len.dat))
# tmp.len.dat = rbind(tmp.len.dat, tmp.len.dat.2)
# 
# annot = array('', dim = c(nrow(tmp.len.dat), 1))
# annot[tmp.len.dat$pval < 0.05] = '*'
# tmp.len.dat$annot = annot
# annot.y = tmp.len.dat$enrich
# annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
# tmp.len.dat$annot.y = annot.y
# p3 = ggplot(tmp.len.dat) + geom_bar(aes(x = mark, y = enrich, fill = mark, alpha = type), stat = 'identity', position = 'dodge') +
#   geom_text(aes(x = mark, y = annot.y), label = annot, position = position_dodge(width = 2)) + 
#   scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_alpha_manual(values = c(1, 0.6), guide = F)

annot = array('', dim = c(nrow(len.dat), 1))
annot[len.dat$fp < 0.05 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$f
#tmp.len.dat$annot.p
#annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p3 = ggplot(len.dat) + geom_bar(aes(x = mark, y = f, fill = mark), color = '#000000', stat = 'identity') + xlab('') + 
  annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in SNPs with high Fst') + theme_bw() + coord_cartesian(ylim = c(0.5, max(annot.y) + 0.25)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_color_discrete(guide = F) +
  theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_fst.pdf', sep = '')), p3, width = 6.5, height = 5.6)

annot = array('', dim = c(nrow(len.dat), 1))
annot[len.dat$pp < 0.05 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$p
annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p4 = ggplot(len.dat) + geom_bar(aes(x = mark, y = p, fill = mark), color = '#000000', stat = 'identity') + xlab('') + annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in SNPs') + theme_bw() + ylim(c(0, 1.2)) + coord_cartesian(ylim = c(0.75, max(annot.y) + 0.05)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_color_discrete(guide = F) +
  theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_snp_enrich.pdf', sep = '')), p4, width = 6.5, height = 5.6)

annot = array('', dim = c(nrow(len.dat), 1))
annot[len.dat$qp < 0.05 & npeaks[, 1] > 0] = '*'
annot.y = len.dat$q
#annot.y[annot.y < 0] = annot.y[annot.y < 0] - 0.01
p5 = ggplot(len.dat) + geom_bar(aes(x = mark, y = q, fill = mark), color = '#000000', stat = 'identity') + xlab('') + annotate('text', x = len.dat$mark, y = annot.y, label = annot, size = 10) + 
  ylab('Enrichment in eQTLs') + theme_bw() +  coord_cartesian(ylim = c(0.5, max(annot.y) + 0.25)) +
  scale_fill_manual(values =  mark.colors(len.dat$mark), guide = F) + scale_color_discrete(guide = F) +
  theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 50, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(plotdir, paste(outpref, 'pop_spec_eQTL_enrich.pdf', sep = '')), p5, width = 6.5, height = 5.6)
