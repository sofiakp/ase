geno.dir = '../../rawdata/variants/all/snps/allNonSan/' # Directory with genotype data, should have a file <indiv>.snps.RData for each individual.
count.dir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/qvals')
plotdir = file.path(Sys.getenv('MAYAROOT'), 'rawdata/alleleCounts/allNonSan/rdata/reps/plots')
if(!file.exists(plotdir)) dir.create(plotdir)

motifs = read.table('../../rawdata/motifs/personal_genome_motifs.txt', header = F, sep = '\t')
neg = as.character(motifs[, 6]) == '-'
motifs = motifs[, c(1,4,7:10)]

filenames = list.files(count.dir, pattern = paste('.*H3K27AC.*rep\\.RData', sep = ''), full.name = T)

for(i in 1:)