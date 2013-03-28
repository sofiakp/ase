import sys
import re
import argparse

# Creates metafiles for peak calling, signal track generation, and segmentations.

parser = argparse.ArgumentParser()
parser.add_argument('metafile')
parser.add_argument('genderfile') # genders.txt
parser.add_argument('qcfile')
parser.add_argument('outpref')
args = parser.parse_args()

# Map each individual to its ID for segmentations
segmap = {}
with open(args.genderfile, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        indiv = fields[0]
        if indiv == 'Snyder':
            indiv = indiv.upper()
        else:
            indiv = 'GM' + indiv
        segmap[indiv] = fields[4]

# Read the metadata file for the datasets and get a list of the filenames for all the paired-end datasets.
# The single-end ones will be ignored.
paired = {}
with open(args.metafile, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        if fields[8] != 'NA':
            name = fields[0].upper() + '_reconcile.dedup.bam' #'_reconcile.dedup.bam'
            if fields[3] != 'Snyder':
                name = 'GM' + name
            name = 'SNYDER_HG19_' + name
            paired[name] = True

combfile = open(args.outpref + '_combrep.tab', 'w') # Metafile for signal track generation and peak calling - all reps in one line, columns 1. rep1;rep2;... 2. input file 3. short pref (eg. SNYDER_HG19_GM12878_CTCF)
#repfile = open(args.outpref + '_rep.tab', 'w') # Metafile for signal track generation and peak calling - reps in different lines
#segfile = open(args.outpref + '_seg_combrep.tab', 'w') # Metafile for segmentation calls. For this, replicates are assumed to be pre-merged into a single file 
namefile = open(args.outpref + '_combrep_names.tab', 'w') # Mappings from signal track and peak filenames to cell line and mark. Signal files (created using combfile) will have names with prefixes col1_vs_col2 where col1,col2 are the first 2 columns of namefile
#namefile2 = open(args.outpref + '_rep_names.tab', 'w') # Same as above, one line per replicate

inputmap = {}
repmap = {}
# Read the qc file to get the fragment lengths. Also keep track of mappings from dataset types to all the 
# corresponding replicates and the corresponding input.
with open(args.qcfile, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        filename = fields[0]
        if not filename in paired:
            print >> sys.stderr, 'Sample', filename, 'is not paired. Skipping.'
            continue
        namefields = filename.split('_')[0:4]
        inputfile = '_'.join(['_'.join(namefields[0:3]), 'INPUT', 'reconcile.dedup.bam']) #'reconcile.dedup.bam'])
        filepref = '_'.join(namefields[0:4]) # Type of dataset, eg SNYDER_HG19_GM12878_H3K27AC
        if namefields[3] != 'INPUT':
            if not filepref in inputmap:
                inputmap[filepref] = inputfile # Map the dataset to its input
                repmap[filepref] = []
            repmap[filepref].append(filename)

all_inputs = ';'.join(list(set(inputmap.values())))

for pref, fs in repmap.iteritems():
    indiv = pref.split('_')[2]
    mark = pref.split('_')[3]
    if not indiv in segmap:
        print >> sys.stderr, indiv, 'is missing from the genders file.'
        continue
    combfile.write('\t'.join([';'.join(fs), all_inputs, pref]) + '\n')
    # If there's more than one replicate, then the segmentations will use files with merged replicates (rep number 0, eg SNYDER_HG19_GM12878_H3K27AC_0_suf)
    if len(fs) > 1:
        #segfile.write('\t'.join([segmap[indiv] + '.' + indiv, mark, pref + '_0_reconcile.dedup.bed.gz', re.sub('.bam', '.bed.gz', inputfile)]) + '\n')
        #segfile.write('\t'.join([segmap[indiv] + '.' + indiv, mark, pref + '_0_dedup.bed.gz', re.sub('.bam', '.bed.gz', inputfile)]) + '\n')
        nameStart = pref
    else:
        #segfile.write('\t'.join([segmap[indiv] + '.' + indiv, mark, re.sub('.bam', '.bed.gz', fs[0]), re.sub('.bam', '.bed.gz', inputfile)]) + '\n')
        nameStart = re.sub('.bam', '', fs[0])
    namefile.write('\t'.join([nameStart, 'SNYDER_HG19_all_INPUT', indiv, mark]) + '\n')
    #for f in fs:
    #    repfile.write('\t'.join([f, inputfile, pref]) + '\n')
    #    namefile2.write('\t'.join([re.sub('.bam', '', f), re.sub('.bam', '', inputfile), indiv, mark]) + '\n')

#fragfile.close()
combfile.close()
#repfile.close()
#segfile.close()
namefile.close()
#namefile2.close()
