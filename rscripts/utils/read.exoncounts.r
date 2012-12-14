rm(list=ls())
library('GenomicRanges')
library('DEXSeq')
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/binom.val.r'))
source(file.path(Sys.getenv('MAYAROOT'), 'src/rscripts/utils/sample.info.r'))

indir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/exonCounts/')
filenames = list.files(indir, pattern = paste('SNYDER_HG19_.*H3K36ME3_.*_ambiguous.exoncounts', sep = ''), full.name = T)
outdir = indir
# Bed files with regions to mask. Should be named <indiv>.blacklist.bed.
mask.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/masks')
# Files with phased SNPs - we want to mark genes that contain unphased SNPs
geno.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/variants/all/snps/')
# Exon metadata: If the RData file does not exist, it will read the gff and create the RData. Otherwise, it will read straight from the RData.
gff.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.flat.gff')
exon.file = file.path(Sys.getenv('MAYAROOT'), 'rawdata/transcriptomes/gencode.v13.annotation.noM.flat.RData')

samples = sample.info(filenames, '_(maternal|paternal|ambiguous).exoncounts')
overwrite = T

# These will be completely removed from the output files
bad.tab = read.table(file.path(Sys.getenv('MAYAROOT'), 'rawdata/genomes_local/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed'), header = F, sep = '\t')
bad.ranges = GRanges(seqnames = Rle(bad.tab[,1]), 
                     ranges = IRanges(start = bad.tab[, 2] + 1, end = bad.tab[, 3]),
                     strand = Rle(rep('+', dim(bad.tab)[1]))) 

nfiles = length(filenames)
filesufs =   c('paternal', 'maternal', 'ambiguous')
varsufs = c('pat', 'mat', 'amb')

for(i in 1:nfiles){
  print(filenames[i])
  outfile = file.path(outdir, gsub('_ambiguous.*counts', '.exoncounts.RData', basename(filenames[i])))
  
  if(overwrite || !file.exists(outfile)){
    # Read metadata if necessary. Otherwise, just read the counts (which is much faster)
    if(!file.exists(exon.file)){
      exon.dat = read.HTSeqCounts(c(filenames[i], gsub('ambiguous', 'maternal', filenames[i])), data.frame(condition = c('amb', 'other')), flattenedfile = gff.file)      
      gene.meta = data.frame(gene = featureData(exon.dat)$geneID, exon = featureData(exon.dat)$exonID, 
                             chr = featureData(exon.dat)$chr, start = featureData(exon.dat)$start,
                             end = featureData(exon.dat)$end, strand = featureData(exon.dat)$strand,
                             transcripts = featureData(exon.dat)$transcripts)
      rownames(gene.meta) = rownames(fData(exon.dat))
      exon.ranges = GRanges(seqnames = Rle(gene.meta$chr), 
                            ranges = IRanges(start = gene.meta$start, end = gene.meta$end),
                            strand = gene.meta$strand) #Rle(rep('+', dim(gene.meta)[1])))
      bad = countOverlaps(exon.ranges, bad.ranges, ignore.strand = T) > 0
      gene.meta = gene.meta[!bad, ]
      gene.meta = droplevels(gene.meta)
      exon.ranges = exon.ranges[!bad, ]
      save(gene.meta, exon.ranges, file = exon.file)
    }else{
      load(exon.file)
    }
    
    counts = array(0, dim = c(dim(gene.meta)[1], length(varsufs)))
    for(g in 1:length(filesufs)){
      count.tab = read.table(gsub('ambiguous', filesufs[g], filenames[i]), header = F)
      ids = count.tab[, 1]
      good = ids %in% rownames(gene.meta)
      stopifnot(all(rownames(gene.meta) == ids[good])) # Make sure the order is the same
      counts[, g] = count.tab[good, 2]
      
      if(g == 1){
        # Check if exon overlaps indels
        indiv = samples$indiv[i]
        mask.file = file.path(mask.dir, paste(indiv, '.blacklist.bed', sep = ''))
        mask.tab = read.table(mask.file, header = F, sep = '\t', comment.char = '#')
        mask.ranges = GRanges(seqnames = Rle(mask.tab[,1]), 
                              ranges = IRanges(start = mask.tab[, 2] + 1, end = mask.tab[, 3]),
                              strand = Rle(rep('+', dim(mask.tab)[1])))  
        mask = countOverlaps(exon.ranges, mask.ranges, ignore.strand = T) > 0
        
        # Check if exons overlap het or unphased SNPs
        geno.file = file.path(geno.dir, paste(indiv, '.snps.r', sep = ''))
        load(geno.file)
        het = geno.info[[paste(indiv, 'mat', sep = '.')]] != geno.info[[paste(indiv, 'pat', sep = '.')]]
        unphased = !geno.info[[paste(indiv, 'phased', sep = '.')]]
        snp.ranges = GRanges(seqnames = Rle(geno.info$chr), 
                             ranges = IRanges(start = geno.info$pos, end = geno.info$pos),
                             strand = Rle(rep('+', length(het))))
        het.ov = countOverlaps(exon.ranges, snp.ranges[het, ], ignore.strand = T) > 0
        unphased.ov = countOverlaps(exon.ranges, snp.ranges[unphased, ], ignore.strand = T) > 0 
      }
    }
    counts = data.frame(counts, row.names = rownames(gene.meta))
    colnames(counts) = varsufs
    pval = binom.val(counts$mat, counts$mat + counts$pat)
    gene.info = data.frame(pval = pval, mask = mask, het = het.ov, phased = !unphased.ov, row.names = rownames(gene.meta))
    sample = samples[i, ]
    save(counts, gene.info, sample, file = outfile)
  }
}
