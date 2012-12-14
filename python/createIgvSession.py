import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('segMeta')
parser.add_argument('bwMeta')
args = parser.parse_args()

ftpdir = ''
outdir = 'http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/users/anshul/temp/chromatinVariation/segmentations/chmmResults/14indiv/final'
suf = '_25_14indiv_dense.bb'

seg_names = []
with open(args.seg_meta, 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        if fields[0] != 'CELLTYPE':
            if fields[0] != 'SNYDER':
                fields[0] = 'GM' + fields[0]
            name = fields[5] + '.' + fields[0]
            if os.path.exists(os.path.join(ftpdir, name + '_25_14indiv_dense.bb')):
                seg_names.append(name)

print '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
print '<Session genome="hg19" locus="chr20:45598936-45623455" version="4">'
print '\t<Resources>'

for s in seg_names:
    print '\t\tResource path = "' + outdir + '/' +  s + suf + '"/>'

print '\t</Resources>'

print '\t<Panel height="3050" name="DataPanel" width="1516">'

for s in seg_names:
    print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="15" height="25" id="' + outdir + '/' + s + suf + ' name="' + re.sub('^C[0-9]*.', '', name) + '" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>'

print '\t</Panel>'
print '\t<Panel height="372" name="FeaturePanel" width="1516">'
print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="138" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>'
print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="10" height="139" id="http://www.broadinstitute.org/igvdata/annotations/hg19/EnsemblGenes.ensGene" name="Ensemble Genes" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>'
print '\t</Panel>'
print '\t<PanelLayout dividerFractions="0.898021308980213"/>'
print '</Session>'
