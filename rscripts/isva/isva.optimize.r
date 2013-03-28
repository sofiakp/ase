rm(list=ls())
source('utils/deseq.utils.r')

genot.pca = new.env()
load('../../rawdata/variants/all_Mar13/genot_pca_noGM19193_exonsVariable_pca.RData', genot.pca)
genot.dims = t(genot.pca$genot.norm) %*% genot.pca$pca.fit$rotation
rownames(genot.dims) = fix.indiv.names(rownames(genot.dims))

load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/SNYDER_HG19_all_reg_H3K27AC_qn.RData')
counts = counts[good.rows, ]
counts.norm = scale(counts) #apply(counts, 2, function(x) (x - row.means[good.rows])/row.sds[good.rows])
regions = regions[good.rows, ]
indivs = colnames(counts)
pop = factor(get.pop(indivs))

isv = isvaFn(counts.norm, pop, ncol(counts), 1:ncol(counts))

for(i in 1:ncol(isv)){
  
}
selisv.m = matrix(isv[, selisv.idx], ncol = length(selisv.idx))
print("Running final multivariate regressions with selected ISVs")
mod.sel = model.matrix( ~ pheno[sel.col] + selisv.m[sel.col, ])
modNULL.sel = model.matrix( ~ selisv.m[sel.col,])

df1 = dim(mod.sel)[2]
df0 = dim(modNULL.sel)[2]

theta = data.mat[, sel.col] %*% mod.sel %*% solve(t(mod.sel) %*% mod.sel)
resid = data.mat[, sel.col] - theta %*% t(mod.sel)
rss1 = rowSums(resid * resid)
thetaNULL = data.mat[, sel.col] %*% modNULL.sel %*% solve(t(modNULL.sel) %*% modNULL.sel)
residNULL = data.mat[, sel.col] - thetaNULL %*% t(modNULL.sel)
rssNULL = rowSums(residNULL * residNULL)

# Get F-statistics for comparing the null and full model
fstats = ((rssNULL - rss1)/(df1 - df0)) / (rss1/(length(sel.col) - df1))
pv.v = 1 - pf(fstats, df1 = (df1 - df0), df2 = (length(sel.col) - df1))
pv.s = sort(pv.v, decreasing = FALSE, index.return = TRUE)
qv.v = p.adjust(pv.s$x, method = 'BH') #qvalue(pv.s$x)$qvalue # TODO: change this to p.adjust
ntop = sum(qv.v < th)
cat('Number of DEGs after ISV adjustment =', ntop, '\n')

rm(list=ls())
source('utils/deseq.utils.r')
load('../../rawdata/signal/combrep/extractSignal/fc/avgSig/merged_Mar13/rdata/SNYDER_HG19_all_reg_H3K27AC_qn_isvaNull.RData')
counts = counts[isva.fit$deg, ]
indivs = colnames(counts)
counts.norm = scale(counts)
counts.no.child = counts.norm[, !(indivs %in% c('GM12878', 'GM19240'))]
pca.fit = prcomp(t(counts.no.child), center = F, scale = F)
pca.dims = t(counts.no.child) %*% pca.fit$rotation
  
genot.pca = new.env()
load('../../rawdata/variants/all_Mar13/genot_pca_noGM19193_exonsVariable_pca.RData', genot.pca)
genot.dims = t(genot.pca$genot.norm) %*% genot.pca$pca.fit$rotation
rownames(genot.dims) = fix.indiv.names(rownames(genot.dims))

common.names = intersect(rownames(pca.dims), rownames(genot.dims))
for(i in 1:ncol(pca.dims)){
  c = cor(pca.dims[match(common.names, rownames(pca.dims)), i], genot.dims[match(common.names, rownames(genot.dims)), 2], method = 'spearman')
  cat(c, '\n')
}