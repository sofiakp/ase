import sys

inf = open("gwascatalog")
inf.readline()

chr_rs = {}
for line in inf:
    if line.split[22][:2] == 'rs':
        chr_rs.add(line.split[11] + '\t' + line.split[22])
inf.close()

outf = open("filteredGWASSet",'w')
for item in chr_rs:
    outs = item + '\n'
    outf.write(outs)
outf.close()