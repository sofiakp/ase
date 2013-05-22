rm(list=ls())
library(reshape)
library(GenomicRanges)
source('utils/sample.info.r')
source('utils/deseq.utils.r')

# Reads a set of region-gene associations (eg enhancer-gene links) and finds the pairs that are supported by allelic biases.
# These are pairs such as the region has AS SNPs and the gene has a bias in the same direction as the SNPs.
# as.files: files with AS SNPs to overlap with associated regions
# gene.files: either files with allelic biases in genes to overlap with genes in associations (if gene.compare is T) 
#             or files with AS SNPs to overlap with genes in associations (if gene.compare is F)
get.as.ov = function(assoc, as.files, gene.files, gene.compare){
    
  load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan/allNonSan.snps.RData')) # SNP metadata
  load(file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')) # gene metadata
  
  # Remove associations where the gene is not in the gene metadata (this can happen with external associations)
  good = assoc$gene.name %in% gene.meta$gene.name
  assoc = assoc[good, ]
  cat('# associations', nrow(assoc), '\n')
  cat('# associations with genes in metadata', sum(good), '\n')
  
  geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/allNonSan')
  indivs = as.character(sample.info(as.files, '.RData')$indiv)
  nindiv = length(indivs)
  
  same.dir.tmp = array(F, dim = c(nrow(assoc), nindiv)) # associations that agree with allelic biases
  diff.dir.tmp = array(F, dim = c(nrow(assoc), nindiv)) # associations that disagree with allelic biases
  
  for(i in 1:nindiv){
    load(file.path(geno.dir, paste(indivs[i], '.snps.RData', sep = '')))
    load(as.files[i])
    cat(indivs[i], '\n')
    sel = as.vector(snp.info$qval < -2) # Get AS SNPs 
    # Find overlaps between enhancer regions and AS SNPs
    ov = findOverlaps(regions.to.ranges(assoc), snps.to.ranges(snp.pos[sel, ]), select = 'all', ignore.strand = T)
    if(length(subjectHits(ov)) == 0){
      cat('No overlaps with AS SNPs\n')
      next
    }
    # Allelic biases for overlapping AS SNPs
    ratios = (counts[sel, 1][subjectHits(ov)] + 1) / (counts[sel, 2][subjectHits(ov)] + 1) # ref / alt 
    mat = as.vector(geno.info$mat[sel])[subjectHits(ov)]
    ratios[mat] = 1 / ratios[mat] # Convert to mat / pat
    ratios = log(ratios)
    
    ov.mat = data.frame(cbind(queryHits(ov), ratios)) # Col 1 - region idx, Col 2 - ratio of overlapping SNP
    sel.assoc = get.consistent.ov(ov.mat) # Get regions where all SNPs agree in direction
    ratios = sel.assoc$ratios
    sel.assoc = sel.assoc$assoc
    
    #cat('# regions with overlapping SNPs:', length(unique(ov.mat[,1])), '\n')
    #cat('# ovelapping SNPs/region:', nrow(ov.mat)/length(unique(ov.mat[,1])), '\n')
    #cat('# regions with all ovelapping SNPs in the same direction:', length(sel.assoc), '\n')
    
    # Get the genes associated with the selected regions and find their indices in the expression counts
    gene.idx = match(assoc$gene.name[sel.assoc], gene.meta$gene.name)
    load(gene.files[i])
    
    if(gene.compare){
      gene.ratios = log((counts$mat[gene.idx] + 1) / (counts$pat[gene.idx] + 1))
      indiv.common = sel.assoc[ratios * gene.ratios > 0 & regions$qval[gene.idx] < 0.01]
      indiv.diff = sel.assoc[ratios * gene.ratios < 0 & regions$qval[gene.idx] < 0.01]
    }else{
      sel = as.vector(snp.info$qval < -2)
      ov = findOverlaps(regions.to.ranges(gene.meta[gene.idx, ]), snps.to.ranges(snp.pos[sel, ]), select = 'all', ignore.strand = T)
      if(length(subjectHits(ov)) == 0){
        indiv.common = c()
        indiv.diff = c()
      }else{
        gene.ratios = (counts[sel, 1][subjectHits(ov)] + 1) / (counts[sel, 2][subjectHits(ov)] + 1) # ref / alt 
        mat = as.vector(geno.info$mat[sel][subjectHits(ov)])
        gene.ratios[mat] = 1 / gene.ratios[mat] # Convert to mat / pat
        gene.ratios = log(gene.ratios)
        
        ov.mat = data.frame(cbind(queryHits(ov), gene.ratios)) # Col 1 - region idx, Col 2 - ratio of overlapping SNP
        gene.assoc = get.consistent.ov(ov.mat) # Get regions where all SNPs agree in direction
        gene.ratios = gene.assoc$ratios
        gene.assoc = gene.assoc$assoc
        # Subset the original regions
        ratios = ratios[gene.assoc]
        sel.assoc = sel.assoc[gene.assoc]
        indiv.common = sel.assoc[ratios * gene.ratios > 0] 
        indiv.diff = sel.assoc[ratios * gene.ratios < 0]
      }  
    }      
    same.dir.tmp[indiv.common, i] = T
    diff.dir.tmp[indiv.diff, i] = T
    #cat('# regions with agreeing allelic biases:', length(indiv.common), '(', length(indiv.common) * 100 / nrow(assoc), '% )\n')
    #cat('# regions with different allelic biases:', length(indiv.diff), '(', length(indiv.diff) * 100 / nrow(assoc), '% )\n')
  }
  tot = sum(apply(same.dir.tmp, 1, function(x) any(x)))
  cat('# associations supported by allelic biases in any individual (same dir):', tot, '(', tot * 100 / nrow(assoc), '% )\n')
  tot = sum(apply(diff.dir.tmp, 1, function(x) any(x)))
  cat('# associations supported by allelic biases in any individual (diff dir):', tot, '(', tot * 100 / nrow(assoc), '% )\n')
  
  # Fix the indices, since some of the associations might have been removed
  same.dir = array(F, dim = c(length(good), nindiv))
  same.dir[good, ] = same.dir.tmp
  diff.dir = array(F, dim = c(length(good), nindiv))
  diff.dir[good, ] = diff.dir.tmp
  return(list(same.dir = same.dir, diff.dir = diff.dir))
}

# Load associations
# load('../../rawdata/enhancers/rdata/enhancer_pred_l1_H3K27AC.RData')
# assoc = dist.regions[links$region.idx, ]
# assoc = cbind(gene.meta$gene.name[links$gene.idx], assoc)
# colnames(assoc) = c('gene.name', 'chr', 'start', 'end')
# assoc$start = assoc$start - 200
# assoc$end = assoc$end + 200
# 
# snp.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
# as.files = list.files(snp.dir, pattern = 'H3K4ME1_rep.RData$', full.names = T)
# indivs = as.character(sample.info(as.files, '.RData')$indiv)
# gene.files = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/rdata/reps/qvals', 
#                        paste('SNYDER_HG19', indivs, 'RZ_rep.RData', sep = '_'))
# good.files = file.exists(gene.files)
# 
# ov = get.as.ov(assoc, as.files[good.files], gene.files[good.files], gene.compare = T)

##### OLD VERSION 
# # Some statistics on overlaps between Jason's links for GM12878 and AS SNPs for H3K27ac / AS genes.
# 
# ench.tab = read.table('rawdata/enhancers/external/GM12878_links_hg19_merged.bed', header = F, sep = '\t')
# ench.ranges = GRanges(seqnames = Rle(ench.tab[,1]), 
#                       ranges = IRanges(start = ench.tab[, 2] + 1, end = ench.tab[, 3]),
#                       strand = Rle(rep('+', dim(ench.tab)[1]))) 
# 
# # AS SNPs
# load('rawdata/alleleCounts/reps_15Sep12/qvals/hitLists/SNYDER_HG19_GM12878_H3K27AC_rep.hits.RData')
# hit.ranges = GRanges(seqnames = Rle(hits$chr), 
#                      ranges = IRanges(start = hits$pos, end = hits$pos),
#                      strand = Rle(rep('+', dim(hits)[1]))) 
# hit.ench.ov = hit.ranges %in% ench.ranges
# # mat / pat for all AS SNPs
# snp.ratios = (hits$ref.counts + 1) / (hits$alt.counts + 1)
# snp.ratios[hits$mat] = 1 / snp.ratios[hits$mat] # convert to mat/pat
# 
# ench.hit.ov.dat = findOverlaps(ench.ranges, hit.ranges)
# avg.ratios = array(NaN, dim = c(dim(ench.tab)[1], 1)) 
# # For each enhancer region, get the overlapping AS SNPs
# for(g in 1:dim(ench.tab)[1]){
#   sel.ratios = snp.ratios[subjectHits(ench.hit.ov.dat)[queryHits(ench.hit.ov.dat) == g]]
#   # If there are overlapping SNPs and the SNPs agree in direction...
#   if(length(sel.ratios) > 0 && (all(sel.ratios >= 0)|| all(sel.ratios <= 0))) avg.ratios[g] = mean(sel.ratios)
# }
# ench.hit.ov = !is.nan(avg.ratios) # Indicators of enhancer regions that overlap AS SNPs
# ench.genes = as.character(ench.tab[ench.hit.ov, 4])
# avg.ratios = avg.ratios[ench.hit.ov]
# cat('Non-distict regions with AS SNPs: ', length(ench.genes), '\n')
# cat('Distict regions with AS SNPs: ', dim(unique(ench.tab[ench.hit.ov, 1:3]))[1], '\n')
# 
# #dim(unique(ench.tab[ench.hit.ov, 1:3]))[1] / dim(unique(ench.tab[, 1:3]))[1] # Fraction of enchancer regions with AS SNPs
# #sum(hit.ench.ov) / dim(hits)[1] # Fraction of AS SNPs in enhancer regions
# 
# load('rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData')
# load('rawdata/geneCounts/rdata/reps/qvals/SNYDER_HG19_GM12878_RNA_rep.genecounts.RData')
# pass.genes = !is.nan(gene.info$qval) & gene.info$qval < 0.01
# pass.gene.names = gene.meta$gene.name[pass.genes] 
# pass.gene.ratios = ((counts$mat + 1) / (counts$pat + 1))[pass.genes]
# 
# # For each enhancer region with AS SNPs find the mat/pat ratio of the associated gene (if it's also AS)
# gene.ratios = array(NaN, dim = c(length(ench.genes), 1))
# for(g in 1:length(ench.genes)){
#   sel = which(ench.genes[g] == pass.gene.names)
#   if(length(sel) == 1) gene.ratios[g] = pass.gene.ratios[sel]
# }
# sum(!is.nan(gene.ratios) & log(gene.ratios) * log(avg.ratios) > 0)
# ovs = ench.hit.ov & ench.tab[, 4] %in% pass.gene.names
