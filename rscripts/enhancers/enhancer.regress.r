rm(list=ls())
library(GenomicRanges)
library(preprocessCore)
library(matrixStats)
library(reshape)
library(ggplot2)
library(glmnet)
library(foreach)
library(doMC)
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/deseq.utils.r'))

registerDoMC(10)

# Runs elastic net regression to link enhancers to genes.

# Dir with counts for RNA. Should contain files SNYDER_HG19_<indiv>_<mark>_0.RData with read counts for genes.
# Each count matrix should have one column per replicate. The file should also have a vector size.factors with
# the size factor (library size) for each column.
rna.dir = '../../rawdata/geneCounts/rdata/repsComb'
peak.dir = '../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata' 
outdir = file.path('../../rawdata/enhancers/merged_Mar13/rdata')
if(!file.exists(outdir)) dir.create(outdir, recursive = T)
mark = 'H3K27AC' # signal in distal elements to correlate with expression
win.len = 100000 # will consider associations in a window [TSS-win.len, TSS+win.len]
perm = ''
nperm = 100
alpha = 0.5
outpref = 'enhancer_coef_elastic0.5_100kb_asinh0.2_cv0.2'

if(perm == ''){
  ######## Create datasets for regression
  
  # Load TSSs. gene.meta has gene information and trans has positions of all TSSs as well as the index of 
  # the corresponding gene in gene.meta.
  load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.trans.RData')
  trans$pos = trans$start
  
  # Get a matrix of average RPKMS for each gene (rows) and individual (cols)
  indivs = unique(as.character(sample.info(list.files(rna.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
  nindiv = length(indivs)
  rna.counts.dat = avg.counts(rna.dir, indivs, 'RZ_0.RData', meta = gene.meta, len.norm = F)
  rna.counts = asinh(rna.counts.dat$counts)
  rna.means = rowMeans(rna.counts, na.rm = T)
  rna.sds = rowSds(rna.counts, na.rm = T)
  rna.cv = rna.sds / rna.means
  expr.genes = which(rna.means > asinh(0.2) & !is.na(rna.cv) & rna.cv > quantile(rna.cv, 0.2, na.rm = T))
  #rna.counts = normalize.quantiles(scale(rna.counts, center = T, scale = F))
  rna.counts = normalize.quantiles(rna.counts)
  trans = trans[trans$gene.idx %in% expr.genes, ] # Remove TSSs for lowly-expressed genes. Leave gene.meta unchanged for now.
  colnames(rna.counts) = fix.indiv.names(indivs)
  indivs = fix.indiv.names(indivs)
  cat('# genes:', nrow(rna.counts), '\n')
  
  # Get the H3K4ME3 peak overlapping each differential TSS. These peaks give the "promoter domains".
  load(file.path(peak.dir, 'SNYDER_HG19_all_reg_H3K4ME3_qn.RData'))
  row.means = rowMeans(counts)
  row.sds = rowSds(counts)
  cvs = row.sds / row.means
  good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.2, na.rm = T)
  me.regions = regions[good.rows, ]
  rna.counts = rna.counts[, colnames(rna.counts) %in% colnames(counts)]
  me.counts = counts[good.rows, match(colnames(rna.counts), colnames(counts))]
  stopifnot(all(colnames(me.counts) == colnames(rna.counts)))
  
  me.ov = findOverlaps(snps.to.ranges(trans), regions.to.ranges(me.regions), select = 'all', ignore.strand = T)
  # Col 1 - gene idx, Col 2 - me3 peak idx
  me.ov.mat = unique(cbind(trans$gene.idx[queryHits(me.ov)], subjectHits(me.ov))) 
  nov = table(me.ov.mat[,1]) # number of me3 peaks overlapping each gene's TSSs
  # Remove genes whose TSSs overlap more than one H3K4me3 peak
  sel.genes = as.numeric(names(nov[nov == 1]))
  me.ov.mat = me.ov.mat[me.ov.mat[,1] %in% sel.genes, ]
  
  cat('# unique active genes (overlapping H3K4me3 peaks):', length(nov), '\n')
  cat('# genes with exactly one TSS-overlapping peak:', sum(nov == 1), '\n')
  
  gene.meta = gene.meta[me.ov.mat[,1], ] # Select genes that do have an overlapping peak.
  rna.counts = rna.counts[me.ov.mat[,1], ]
  me.regions = me.regions[me.ov.mat[, 2], ]
  
  # Load H3K27ac signal and peaks
  load(file.path(peak.dir, paste('SNYDER_HG19_all_reg', mark, 'qn.RData', sep = '_')))
  row.means = rowMeans(counts)
  row.sds = rowSds(counts)
  cvs = row.sds / row.means
  good.rows = !is.na(cvs) & row.means > asinh(0.2) & cvs > quantile(cvs, 0.2, na.rm = T)
  ac.regions = regions[good.rows, ]
  ac.counts = counts[good.rows, match(colnames(rna.counts), colnames(counts))]
  stopifnot(all(colnames(ac.counts) == colnames(rna.counts)))
  
  # For each gene, get the "proximal" K27ac peaks (those within the TSS associated H3K4me3 peak), 
  # as well as all K27ac peaks in a large window around the gene.
  ov.regions = me.regions
  middles = round((me.regions$end + me.regions$start) / 2)
  ov.regions$start = middles - win.len
  ov.regions$end = middles + win.len
  tss.ov = findOverlaps(regions.to.ranges(me.regions), regions.to.ranges(ac.regions), ignore.strand = T, select = 'all')
  all.ov = findOverlaps(regions.to.ranges(ov.regions), regions.to.ranges(ac.regions), ignore.strand = T, select = 'all')
  
  coef.dat = NULL
  # Col 1 - Sum of squared errors, deviance(model) = sum((model$a0[i] + model$beta[,i]%*%t(x) - y)^2)
  # Col 2 - dev.ratio (see glmnet doc) 1 - deviance(model) / nulldev, where nulldev is the deviance of the intercept model so sum((y - mean(y))^2)
  # Col 3 - lambda
  # Col 4 - intercept
  fit.dat = array(NaN, dim = c(nrow(rna.counts), 4)) 
  for(i in 1:nrow(rna.counts)){
    y = rna.counts[i, ] 
    sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
    sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
    
    if(length(sel.dist) > 0){
      if(length(sel.dist) > 1){x = as.matrix(t(ac.counts[sel.dist, ]))
      }else{x = as.matrix(ac.counts[sel.dist, ])}
      # en.model = glmnet(y = y, x = x, family = 'gaussian', alpha = alpha, standardize = T)
      # deviance(en.model)[i] is the sum of squared errors: sum((en.model$a0[i] + en.model$beta[,i]%*%t(x) - y)^2)
      # where i indexes the values of lambda. So this is what the model minimizes
      # 0.5 * deviance(en.model) / en.model$nobs + en.model$lambda * (0.5 * (1 - alpha) * colSums(en.model$beta^2) + alpha * colSums(abs(en.model$beta)))
      
      # Use leave-one-out CV to pick the optimal lambda
      cv.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = alpha, standardize = T, grouped = F, nfolds = length(y))
      k = which(cv.model$lambda == cv.model$lambda.min)
      coef.dat = rbind(coef.dat, data.frame(gene.idx = rep(i, length(sel.dist)),
                                            region.idx = sel.dist, 
                                            coef = cv.model$glmnet.fit$beta[, k])) # Coefficients for the optimal lambda
      fit.dat[i, ] = c(deviance(cv.model$glmnet.fit)[k], cv.model$glmnet.fit$dev.ratio[k], cv.model$lambda.min, as.numeric(cv.model$glmnet.fit$a0[k]))
    }
  }
  outfile = file.path(outdir, paste(outpref, '_', mark, '.RData', sep = ''))
  if(!file.exists(outdir)) dir.create(outdir)
  colnames(fit.dat) = c('deviance', 'dev.ratio', 'lambda', 'intercept')
  tss = middles
  save(ac.regions, ac.counts, rna.counts, gene.meta, fit.dat, coef.dat, tss, file = outfile)  
}else{
  # Load data. You should have run with perm = '' first.
  load(file.path(outdir, paste(outpref, '_', mark, '.RData', sep = '')))
  outdir = file.path(outdir, perm)
  if(!file.exists(outdir)) dir.create(outdir)
  
  foreach(n = 1:nperm) %dopar% {
    
    cat('Round', n, '\n')
    if(perm == 'perm_gene'){
      # Shuffle the selected genes.
      shuf = sample(1:nrow(rna.counts), nrow(rna.counts))
      rna.counts = rna.counts[shuf, ]
    }else if(perm == 'perm_indiv'){
      # For each gene shuffle the individuals
      for(i in 1:nrow(rna.counts)){
        shuf = sample(1:nindiv, nindiv)
        rna.counts[i, ] = rna.counts[i, shuf]
      }
    }
    fit.dat = array(NaN, dim = c(nrow(rna.counts), 4))
    scan.genes = unique(coef.dat$gene.idx) # Not all genes were fitted, so not all of them are in coef.dat.
    for(i in 1:length(scan.genes)){
      y = rna.counts[scan.genes[i], ]
      sel.ind = findInterval(c(scan.genes[i] - 1, scan.genes[i]), coef.dat$gene.idx) # Indices of coefficients for the i-th gene in coef.dat
      sel.ind = (sel.ind[1] + 1):sel.ind[2]
      sel.dist = coef.dat$region.idx[sel.ind]
      if(length(sel.dist) > 1){x = as.matrix(t(ac.counts[sel.dist, ]))
      }else{x = as.matrix(ac.counts[sel.dist, ])}
      cv.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = alpha, standardize = T, grouped = F, nfolds = length(y))
      k = which(cv.model$lambda == cv.model$lambda.min)
      fit.dat[scan.genes[i], ] = c(deviance(cv.model$glmnet.fit)[k], cv.model$glmnet.fit$dev.ratio[k], cv.model$lambda.min, as.numeric(cv.model$glmnet.fit$a0[k]))
    }
    outfile = file.path(outdir, paste(outpref, '_', perm, n, '_', mark, '.RData', sep = ''))
    colnames(fit.dat) = c('deviance', 'dev.ratio', 'lambda', 'intercept')
    save(fit.dat, file = outfile)
  }
}


# for(n in 1:nperm){
#   cat('Round', n, '\n')
#   if(perm == 'perm_gene'){
#     # Shuffle the selected genes.
#     shuf = sample(1:nrow(rna.counts), nrow(rna.counts))
#     rna.counts = rna.counts[shuf, ]
#   }else if(perm == 'perm_indiv'){
#     # For each gene shuffle the individuals
#     for(i in 1:nrow(rna.counts)){
#       shuf = sample(1:nindiv, nindiv)
#       rna.counts[i, ] = rna.counts[i, shuf]
#     }
#   }else if(perm != '') stop('Invalid perm value')
#   
#   coef.dat = NULL
#   dev = array(NaN, dim = c(nrow(rna.counts), 1))
#   for(i in 1:nrow(rna.counts)){
#     y = rna.counts[i, ] 
#     sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
#     sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
#     
#     if(length(sel.dist) > 0){
#       if(length(sel.dist) > 1){x = as.matrix(t(ac.counts[sel.dist, ]))
#       }else{x = as.matrix(ac.counts[sel.dist, ])}
#       en.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = 0.5, standardize = T, grouped = F)
#       k = which(en.model$lambda == en.model$lambda.min)
#       #cat('#variables', length(sel.dist), '#selected', en.model$glmnet.fit$df[k], '\n')
#       coef.dat = rbind(coef.dat, data.frame(gene.idx = rep(i, length(sel.dist)),
#                                             region.idx = sel.dist,
#                                             coef = en.model$glmnet.fit$beta[, k]))
#       dev[i] = en.model$glmnet.fit$dev.ratio[k]
#     }
#   }
#   if(perm == ''){
#     outfile = file.path(outdir, paste(outpref, '_', mark, '.RData', sep = ''))
#   }else{
#     outfile = file.path(outdir, paste(outpref, n, '_', mark, '.RData', sep = ''))
#   }
#   save(ac.regions, ac.counts, rna.counts, gene.meta, dev, coef.dat, file = outfile)
# }

# links = NULL
# for(i in 1:nrow(rna.counts)){
#   y = rna.counts[i, ] 
#   sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
#   sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
#   
#   if(length(sel.dist) > 0){
#     if(length(sel.dist) > 1){x = as.matrix(t(ac.counts[sel.dist, ]))
#     }else{x = as.matrix(ac.counts[sel.dist, ])}
#     l1.model = cv.glmnet(y = y, x = x, family = 'gaussian', alpha = 1, standardize = T, grouped = F)
#     k = which(l1.model$lambda == l1.model$lambda.min)
#     coef = sort(abs(l1.model$glmnet.fit$beta[, k]), decreasing = T, index.return = T)
#     sel = coef$ix[coef$x > 1e-4]
#     if(length(sel) > 0){
#       sel = sel[1:min(length(sel), 10)]
#       glm.mat = data.frame(x[, sel])
#       sel.dist = sel.dist[sel]
#       glm.mat$rna = y  
#       glm.model = glm(paste('rna~', paste(colnames(glm.mat[1:length(sel.dist)]), collapse ='+'), sep=''), family = 'gaussian', data = glm.mat)
#       glm.sum = summary(glm.model)
#       pvals = glm.sum$coefficients[2:nrow(glm.sum$coefficients), 4]
#       selp = pvals < 1e-2
#       if(any(selp)){
#         if(is.null(links)){links = cbind(rep(i, sum(selp)), sel.dist[selp], pvals[selp])
#         }else{links = rbind(links, cbind(rep(i, sum(selp)), sel.dist[selp], pvals[selp]))} 
#       }
#     }
#   }
# }
# links = data.frame(links)
# colnames(links) = c('gene.idx', 'region.idx', 'pval')
# dist.regions = ac.regions
# save(dist.regions, rna.counts, gene.meta, links, file = file.path(outdir, paste('enhancer_pred_l1_', mark, '.RData', sep = '')))

#### OLD version
# load('../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData') # Gene metadata
# # Dir with DESeq results for RNA (see run.deseq.r)
# rna.deseq.dir = '../../rawdata/geneCounts/rdata/repsComb/deseq/'
# # Dir with counts for RNA. Should contain files SNYDER_HG19_<indiv>_<mark>_0.RData with read counts for genes.
# # Each count matrix should have one column per replicate. The file should also have a vector size.factors with
# # the size factor (library size) for each column.
# rna.dir = '../../rawdata/geneCounts/rdata/repsComb'
# peak.deseq.dir = '../../rawdata/signal/combrep/countsAtPeaksBroad/repsComb/deseq'
# peak.dir = '../../rawdata/signal/combrep/countsAtPeaksBroad/repsComb/' 
# outdir = file.path(rna.dir, 'plots')
# 
# # Load the DESeq results for RNA and get a list of genes that are differential between several pairs of individuals.
# rna.diff = get.diff.count(list.files(rna.deseq.dir, pattern = '*RZ_deseq.RData', full.names = T), 0.01) > 0
# 
# # Get a matrix of average RPKMS for each gene (rows) and individual (cols)
# # TODO: the gene metafile has only the longest transcript per gene. We need to consider all of them and compute the corresponding transcript lengths.
# indivs = unique(as.character(sample.info(list.files(rna.dir, pattern = 'SNYDER_HG19_.*RZ_0.RData'), '.RData')$indiv))
# indivs = indivs[!(indivs %in% c('GM12878', 'GM19240'))] # REMOVE daughters
# nindiv = length(indivs)
# rna.counts.dat = avg.counts(rna.dir, indivs, 'RZ_0.RData', meta = gene.meta) # This will normalize counts by the transcript length
# rna.counts = rna.counts.dat$counts[rna.diff, ]
# gene.meta = gene.meta[rna.diff, ]
# rna.counts = get.medoid(as.matrix(rna.counts))
# 
# #row.means = rowMeans(rna.counts) # Center data by subtracting row means
# #pca.fit = prcomp(t(rna.counts), center = F, scale = F)
# #rna.counts = rna.counts - rep(rowSums(pca.fit$rotation[, 1:2]), dim(rna.counts)[2])
#   
# # TSSs. We will overlap these with the H3K4me3 peaks. 
# tss.meta = gene.meta 
# pos = tss.meta$strand == '+'
# tss.meta$end[pos] = tss.meta$start[pos]
# tss.meta$start[!pos] = tss.meta$end[!pos]
# tss.ranges = regions.to.ranges(tss.meta)
# 
# # Get the H3K4ME3 peak overlapping each differential TSS. These peaks give the "promoter domains".
# load(file.path(peak.dir, 'SNYDER_HG19_GM12878_H3K4ME3_0.RData')) # Load any H3K4ME3 file. We just need the regions.
# me3.regions = regions 
# 
# # For each TSS, find the overlapping H3K4ME3 peak (sincepeaks don't overlap, there should be at most one overlapping peak for each TSS)
# me3.ov = findOverlaps(tss.ranges, regions.to.ranges(me3.regions), ignore.strand = T, select = 'first')
# gene.meta = gene.meta[!is.na(me3.ov), ] # Select genes that do have an overlapping peak.
# tss.meta = tss.meta[!is.na(me3.ov), ]
# tss.ranges = regions.to.ranges(tss.meta)
# rna.counts = rna.counts[!is.na(me3.ov), ]
# me3.regions = me3.regions[me3.ov[!is.na(me3.ov)], ] 
# me3.ov = me3.ov[!is.na(me3.ov)]
# 
# # Get the differential H3K27ac peaks and their corresponding counts. Use a permissive p-value cutoff.
# ac.diff = get.diff.count(list.files(peak.deseq.dir, pattern = '*H3K27AC_deseq.RData', full.names = T), 0.01) > 0
# ac.counts.dat = avg.counts(peak.dir, indivs, 'H3K27AC_0.RData', meta = NULL)
# ac.regions = ac.counts.dat$regions[ac.diff, ]
# ac.counts = ac.counts.dat$counts[ac.diff, ]
# ac.ranges = regions.to.ranges(ac.regions)
# ac.counts = get.medoid(as.matrix(ac.counts))
# 
# #row.means = rowMeans(ac.counts) # Center data by subtracting row means
# #pca.fit = prcomp(t(ac.counts), center = F, scale = F)
# #ac.counts = ac.counts - rep(rowSums(pca.fit$rotation[, 1:3]), dim(ac.counts)[2])
# 
# # For each gene, get the "proximal" K27ac peaks (those within the TSS associated H3K4me3 peak), 
# # as well as all K27ac peaks in a large window around the gene.
# ov.ranges = flank(tss.ranges, width = 500000, both = T)
# tss.ov = findOverlaps(regions.to.ranges(me3.regions), ac.ranges, ignore.strand = T, select = 'all')
# all.ov = findOverlaps(ov.ranges, ac.ranges, ignore.strand = T, select = 'all')
# 
# score.true = array(NaN, dim=c(length(me3.ov), 1))
# score.perm = array(NaN, dim=c(10 * length(me3.ov), 1))
# best.cell = array(0, dim=c(length(me3.ov), 1))
# best.reg = array(0, dim=c(length(me3.ov), 1))
# 
# for(i in 1:length(score.true)){
#   sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
#   sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
# 
#   if(length(sel.dist) > 0){
#     x = as.matrix(rna.counts[i, ])
#     x = t(x[, rep(1, length(sel.dist))])
#     y = ac.counts[sel.dist, ]
# 
#     score.tmp = ars(x,y)
#     score.true[i] = max(score.tmp, na.rm = T)
#     ind = which(score.tmp == score.true[i], arr.ind = T)
#     if(length(sel.dist) > 1){
#       best.cell[i] = ind[1, 2]
#       best.reg[i] = sel.dist[ind[1, 1]]
#     }else{
#       best.cell[i] = ind[1]
#       best.reg[i] = sel.dist[1]
#     }
# #     for(p in 1:10){
# #       perm = sample(1:nindiv)
# #       if(length(sel.dist) > 1){
# #         score.perm[10 * (i - 1) + p] = max(ars(x, y[, perm]), na.rm = T)
# #       }else{score.perm[10 * (i - 1) + p] = max(ars(x, y[perm]), na.rm = T)}
# #     }
#   }  
# }
# 
# ord = order(score.true, decreasing = T, na.last = T)[1:500]
# names = gene.meta$gene.name[ord]
# outfile = file.path(outdir, paste('asr_genes.txt', sep = ''))
# write.table(names, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")
# outreg = ac.regions[best.reg[ord],]
# outreg$gene = names
# outfile = file.path(outdir, paste('asr_assoc.txt', sep = ''))
# write.table(outreg, file = outfile, quote = F, row.names = F, col.names = F, sep = "\t")

# Solution 2: Use an L1-regression to get the non-proximal peaks that affect the gene's expression
# pval = array(1, dim=c(length(me3.ov), 1))
# for(i in 1:length(pval)){
#   pred = data.frame(rna = t(rna.counts[i, ]))
#   colnames(pred) = c('rna')
#   sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
#   sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
#   if(length(sel.dist) > 0){
#     for(j in 1:length(sel.dist)){
#       name = rownames(ac.counts)[sel.dist[j]]
#       pred[[name]] = t(ac.counts[sel.dist[j], ])
#     }
#     fit.alt = glmnet(as.matrix(pred[, 2:length(sel.dist)]), pred[, 1], family='gaussian', alpha = 1, standardize = F)
#     tmp = array(1, dim = c(length(fit.alt$df), 1))
#     for(j in 1:length(tmp)) tmp[j] = 1 - pchisq(2 * fit.alt$dev.ratio[j] * fit.alt$nulldev, fit.alt$df[j])
#     pval[i] = min(tmp)
#     if(pval[i] < 0.1){
#       cat(i, pval[i], '\n')
#       #print(fit.alt)
#     }
#   }
# }

# Solution 1: For each gene, test if the non-proximal enhancer regions significantly improve the fit of the gene's expression
# compared to the non-proximal ones.
# pval = array(1, dim=c(length(me3.ov), 1))
# for(i in 1:length(pval)){
#   pred = data.frame(rna = t(rna.counts[i, ]))
#   colnames(pred) = c('rna')
#   sel.prox = subjectHits(tss.ov)[queryHits(tss.ov) == i] # Get the proximal peaks for gene i.
#   if(length(sel.prox) > 0){
#     for(j in 1:length(sel.prox)){
#       name = rownames(ac.counts)[sel.prox[j]]
#       pred[[name]] = t(ac.counts[sel.prox[j], ])
#       if(j == 1){null = paste('rna', name, sep = '~')
#       }else{null = paste(null, name, sep = '+')}
#     }
#   }else{null = 'rna ~ 1'}
#   sel.dist = setdiff(subjectHits(all.ov)[queryHits(all.ov) == i], sel.prox)
#   if(length(sel.dist) > 0){
#     alt = null
#     # Select the enhancers most strongly correlated with the gene's expression.
#     if(length(sel.dist) > 1){
#       ac.co = cor(t(rna.counts[i, ]), t(ac.counts[sel.dist, ]))      
#     }else{
#       ac.co = cor(rna.counts[i, ], ac.counts[sel.dist, ])
#     }
#     sel.dist = sel.dist[order(ac.co, decreasing = T)[1:min(nindiv - 1 - length(sel.prox), length(sel.dist))]]
#     for(j in 1:length(sel.dist)){
#       name = rownames(ac.counts)[sel.dist[j]]
#       pred[[name]] = t(ac.counts[sel.dist[j], ])
#       alt = paste(alt, name, sep = '+')
#     }
#     fit.null = glm(null, family='gaussian', data = pred)
#     fit.alt = glm(alt, family='gaussian', data = pred)
#     pval[i] = 1 - pchisq(fit.null$deviance - fit.alt$deviance, fit.null$df.residual - fit.alt$df.residual)
#     if(pval[i] < 0.1){
#       cat(i, pval[i], '\n')
#       print(fit.alt)
#     }
#   }
# }
