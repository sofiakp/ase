rm(list=ls())
library('GenomicRanges')
#library(foreach)
#library(doMC)
source('utils/binom.val.r')
source('utils/sample.info.r')
source('utils/deseq.utils.r')

#registerDoMC(10)
indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/geneCounts/')
filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*19193.*_ambiguous.genecounts', sep = ''), full.name = T)
outdir = file.path(indir, 'rdata')
if(!file.exists(outdir)) dir.create(outdir, recursive = T)
# Bed files with regions to mask. Should be named <indiv>.blacklist.bed.
mask.dir = '../../rawdata/variants/novelCalls/filtered//masks'
# Files with phased SNPs - we want to mark genes that contain unphased SNPs
geno.dir = '../../rawdata/variants/novelCalls/filtered/snps/'
# File with SNP metadata
load('../../rawdata/variants/sanConsensus/snps/san.snps.RData')
# Gene metadata: If the RData file does not exist, it will read the gtf and create the RData. Otherwise, it will read straight from the RData.
gtf.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.gtf'
gene.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.genes.RData'
# Exon metadata. MUST EXIST (See read.exoncounts)
exon.file = '../../rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData'

samples = sample.info(filenames, '_(maternal|paternal|ambiguous).genecounts')
overwrite = F

# These will be completely removed from the output files
bad.tab = read.table('../../rawdata/genomes_local/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed', header = F, sep = '\t')
bad.ranges = GRanges(seqnames = Rle(bad.tab[,1]), 
                     ranges = IRanges(start = bad.tab[, 2] + 1, end = bad.tab[, 3]),
                     strand = Rle(rep('+', dim(bad.tab)[1])))  

last.rows = c('ambiguous', 'no_feature', 'too_low_aQual', 'not_aligned', 'alignment_not_unique')
nfiles = length(filenames)
filesufs =   c('paternal', 'maternal', 'ambiguous')
varsufs = c('pat', 'mat', 'amb')

for(i in 1:nfiles){
  print(filenames[i])
  outfile = file.path(outdir, gsub('_ambiguous.genecounts', '.genecounts.RData', basename(filenames[i])))
  
  if(overwrite || !file.exists(outfile)){
    # Read metadata if necessary. Otherwise, just read the counts (which is much faster)
    if(!file.exists(gene.file)){
      gff = read.table(gtf.file, header = F, sep = '\t')
      gene.chr = gff[, 1]
      gene.starts = gff[, 4] + 1
      gene.ends = gff[, 5]
      gene.strands = gff[, 7]
      gene.ids = gsub(';.*$', '', gsub('gene_id ', '', gff[, 9]))
      gene.names = gsub(';.*$', '', gsub('.*gene_name ', '', gff[, 9]))
      gene.meta = data.frame(gene.name = gene.names, chr = gene.chr, start = gene.starts, end = gene.ends, strand = gene.strands, row.names = gene.ids)
      gene.ranges = GRanges(seqnames = Rle(gene.chr), 
                            ranges = IRanges(start = gene.starts, end = gene.ends),
                            strand = Rle(gene.strands))
      exons = new.env()
      load(exon.file, exons)
      all.genes = as.character(exons$gene.meta$gene)
      len = array(0, dim = c(length(gene.ids), 1))
      for(j in 1:length(all.genes)){
        if(j %% 1000 == 0) print(j)
        sel = which(gene.ids %in% unlist(strsplit(all.genes[j], '\\+')))
        stopifnot(length(sel) >= 1)
        len[sel] = len[sel] + exons$gene.meta$end[j] - exons$gene.meta$start[j] + 1 
      }
      rem = len == 0
      len[rem] = gene.ends[rem] - gene.starts[rem] + 1 
      gene.meta$len = len
      bad = countOverlaps(gene.ranges, bad.ranges, ignore.strand = T) > 0
      gene.meta = droplevels(gene.meta[!bad, ])
      gene.ranges = gene.ranges[!bad, ]
      the.order = order(rownames(gene.meta))
      gene.meta = gene.meta[the.order, ]
      gene.ranges = gene.ranges[the.order, ]
      save(gene.meta, gene.ranges, file = gene.file)
    }else{
      load(gene.file)
    }    
    counts = array(0, dim = c(dim(gene.meta)[1], length(varsufs)))
    for(g in 1:length(filesufs)){
      count.tab = read.table(gsub('ambiguous', filesufs[g], filenames[i]), header = F)
      #count.tab = count.tab[!(count.tab[, 1] %in% last.rows), ]
      count.tab = count.tab[order(count.tab[, 1]), ]
      count.tab = count.tab[count.tab[, 1] %in% rownames(gene.meta), ]
      stopifnot(all(rownames(gene.meta) == count.tab[, 1]))
      counts[, g] = count.tab[, 2]
      
      # Make sure the genes are in the same order as in the GFF
      #genes.tmp = as.character(count.tab[, 1])
      #rank.tab = sort(genes.tmp, index.return = T)
      #correct.order = rank.tab$ix[gene.rank]
      #stopifnot(all(genes.tmp[correct.order] == gene.ids))
      #counts[, g] = count.tab[correct.order, 2]
      
      if(g == 1){
        indiv = samples$indiv[i]
        mask.file = file.path(mask.dir, paste(indiv, '.blacklist.bed', sep = ''))
        mask.ranges = regions.to.ranges(read.bed(mask.file))  
        mask = countOverlaps(gene.ranges, mask.ranges, ignore.strand = T) > 0
        
        geno.file = file.path(geno.dir, paste(indiv, '.snps.RData', sep = ''))
        load(geno.file)
        het = as.vector(geno.info$mat != geno.info$pat)
        unphased = as.vector(geno.info$unphased)
        het.ov = countOverlaps(gene.ranges, snps.to.ranges(snp.pos[het, ]), ignore.strand = T) > 0
        unphased.ov = countOverlaps(gene.ranges, snps.to.ranges(snp.pos[unphased, ]), ignore.strand = T) > 0
      }
    }
    counts = data.frame(counts, row.names = rownames(gene.meta))
    colnames(counts) = varsufs
    pval = binom.val(counts$mat, counts$mat + counts$pat)
    regions = data.frame(pval = pval, mask = mask, het = het.ov, phased = !unphased.ov, row.names = rownames(gene.meta))
    sample = droplevels(samples[i, ])
    save(counts, regions, sample, file = outfile)
  }
}
