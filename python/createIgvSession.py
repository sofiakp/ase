import sys
import argparse
import os.path
import re
import fileinput

parser = argparse.ArgumentParser(description = 'Creates IGV session file. Reads a list of signal tracks from STDIN and a list of individuals from a file. Makes a session with the segmentations for these individuals and all the signal tracks read')
parser.add_argument('segMeta', help = 'Metafile with ids for segmentations (like genders.txt)')
parser.add_argument('segdir', help = 'Segmentation files should be in this location')
parser.add_argument('bwdir', help = 'Signal tracks should be in this location')
args = parser.parse_args()

segdir = args.segdir
bwdir = args.bwdir
suf = '_25_14indiv_dense.bb'

colors = {'H3K27AC':'255,204,0', 'H3K4ME3':'255,102,102', 'H3K4ME1':'204,0,51'} # Signal track colors - blue is default

# Read the names of the cell lines used for segmentation
seg_names = []
with open(args.segMeta, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        if fields[0] != 'CELLTYPE':
            if fields[0] != 'Snyder':
                name = fields[4] + '.GM' + fields[0]
            else:
                name = fields[4] + '.SNYDER'
            seg_names.append(name)

print '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
print '<Session genome="hg19" locus="chr20:45598936-45623455" version="4">'
print '\t<Resources>'
print '\t\t<Resource name="Ensemble Genes" path="http://www.broadinstitute.org/igvdata/annotations/hg19/EnsemblGenes.ensGene"/>'

# Add list of segmentation files to the "Resources"
for s in seg_names:
    print '\t\t<Resource path="' + segdir + '/' +  s + suf + '"/>'

# Read the list of signal tracks and their corresponding cell lines and marks and add them to the "Resources"
lines = {}
marks = {}
for line in fileinput.input([]):
    fields = line.strip().split()
    filename = fields[0] + '_VS_' + fields[1] + '.fc.signal.bw'
    lines[filename] = fields[2]
    marks[filename] = fields[3]
    print '\t\t<Resource path="' + bwdir + '/' + filename + '"/>'

print '\t</Resources>'

# Create track lines
print '\t<Panel height="2050" name="DataPanel" width="1516">'
print '\t\t<Track altColor="0,0,178" color="0,0,178" colorScale="ContinuousColorScale;0.0;236.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="25" id="hg19_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count">\n\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="236.0" minimum="0.0" type="LINEAR"/>\n\t\t</Track>'

for s in seg_names:
    print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="35" id="' + segdir + '/' + s + suf + '" name="' + re.sub('^C[0-9]*.', '', s) + '" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>'

for f, m in marks.iteritems():
    if m in colors:
        c = colors[m]
    else:
        c = '0,0,255' # Default color is blue
    print '\t\t<Track altColor="0,0,178" autoscale="false" color="' + c + '" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="35" id="' + bwdir + '/' + f + '" name="' + lines[f] + '_' + marks[f] + '" renderer="BAR_CHART" showDataRange="true" visible="true" windowFunction="mean">\n\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>\n\t\t</Track>'

print '\t</Panel>'
# Gene tracks
print '\t<Panel height="372" name="FeaturePanel" width="1516">'
print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="138" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>'
print '\t\t<Track altColor="0,0,178" color="0,0,178" displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="12" height="60" id="http://www.broadinstitute.org/igvdata/annotations/hg19/EnsemblGenes.ensGene" name="Ensemble Genes" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>'
print '\t</Panel>'
print '\t<PanelLayout dividerFractions="0.8"/>'
print '</Session>'
