import sys
import os
import os.path
import re

# Removes deprecated files from combrep. Eg. if a dataset initially had one replicate only,
# its name would include the replicate number as in hg19_w1k_AT_SNYDER_HG19_GM12892_BUB_1_reconcile.dedup.txt.
# When a second replicate is added, the new file will be hg19_w1k_AT_SNYDER_HG19_GM12892_BUB.txt (because these
# are the naming conventions in generateMacs2NormSignalTracks.sh), so the old file won't be overwritten.
# Run like this: python removeRedundantFiles.py ../metadata/chromatinVariation_combrep_names.tab .

meta = sys.argv[1] # Metafile. First column has the dataset names to be used. Files with other names are considered to be deprecated and removed.
indir = sys.argv[2]

filepref = {}
with open(meta, 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        filepref[fields[0]] = True

subdirs = ['fc/avgSig/textFiles', 'fc/avgSig/matFiles']
sufs = ['.txt', '.mat']

for i, s in enumerate(subdirs):
    for c in os.listdir(os.path.join(indir, s)):
        if os.path.isfile(os.path.join(indir, s, c)):
            pref = os.path.basename(c).split('_AT_')[1]
            pref = re.sub(sufs[i], '', pref)
            if not pref in filepref:
                print 'Will remove', os.path.join(indir, s, c)
                #os.remove(os.path.join(indir, s, c))
