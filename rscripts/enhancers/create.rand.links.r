rm(list=ls())

infile = '../../rawdata/enhancers/rdata/enhancer_coef_ars_100kb_asinh0.2_cv0.2_H3K27AC_links_fdr0.01_perm_gene_pairs.RData'
rand.type = 'rand'
load(infile)

if(rand.type == 'nonlinks'){
  neg.links = coef.dat[coef.dat$gene.idx %in% links$gene.idx & !(coef.dat$region.idx %in% links$region.idx), 1:2]
  neg.links = neg.links[sample(1:nrow(neg.links), min(nrow(neg.links), nrow(links))), ]
  links = neg.links
}else if(rand.type == 'rand'){
  links = coef.dat[sample(1:nrow(coef.dat), nrow(links)), ]
}
save(ac.regions, ac.counts, rna.counts, gene.meta, links, file = gsub('.RData', paste('_', rand.type, '.RData', sep = ''), infile))
