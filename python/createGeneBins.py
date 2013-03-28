import sys

nbins = 40
bin_len = 100
names = ['TSS_' + str(n+1) for n in range(nbins)]
names.append('gene_body')
names.extend(['TTS_' + str(n) for n in range(nbins, 0, -1)])

with open(sys.argv[1], 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        chrom = fields[0]
        gene_start = int(fields[1])
        gene_end = int(fields[2])
        gene_len = gene_end - gene_start
        strand = fields[5]
        if gene_len < 4100:
            continue
        bins = []
        bin_start = gene_start - int(nbins * bin_len / 2)
        for i in range(nbins):
            start = bin_start + i * bin_len
            end = bin_start + (i + 1) * bin_len 
            bins.append('\t'.join([str(start), str(end)]))
        # The +1 is to ensure that the TSSs and TTSs of positive and negative genes always lie on the same bins.
        bin_start = gene_end - int(nbins * bin_len / 2) + 1
        bins.append('\t'.join([str(end), str(bin_start)]))
        for i in range(nbins):
            start = bin_start + i * bin_len
            end = bin_start + (i + 1) * bin_len 
            bins.append('\t'.join([str(start), str(end)]))
        if strand == '-':
            out_names = names[::-1]
        else:
            out_names = names
        for i,b in enumerate(bins):
            print '\t'.join([chrom, b, fields[6], '.', strand, out_names[i]])
