import sys
import argparse

parser = argparse.ArgumentParser(description = 'Outputs the sequence for the regions in a bed file, one FASTA entry per region. All letters are converted to uppercase. The input sequence can be either a single FASTA file with all the chromosomes, or one file per chromosome, named chr*.fa. In the latter case, the input should be sorted by chromosome.')
parser.add_argument('fapref', help = 'Prefix of fasta files. If you are providing a single FASTA, this should be the path to the fasta without  the .fa extension')
parser.add_argument('bed', help = 'Input bed file')
parser.add_argument('-m', '--mergedFa', action = 'store_true', default = False, help = 'The input is a merged FASTA.')
args = parser.parse_args()
mergedFa = args.mergedFa
fapref = args.fapref
bed = args.bed

# Read all the chromosomes (this requires a lot of memory).
if mergedFa:
    seqMap = {}
    name = ''
    seq = ''
    with open(fapref + '.fa', 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if name != '':
                    # Save the previous chromosome and reset the sequence
                    seqMap[name] = seq
                    seq = ''
                name = line.strip().strip('>')
                print >> sys.stderr, 'Reading', name
            else:
                seq = seq + line.strip().upper()

# Read the bed file
lastSeq = ''
with open(args.bed, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        print '>' + '_'.join([chrom, str(start), str(end)])
        if mergedFa: 
            print seqMap[chrom][start:end]
        else:
            if lastSeq != chrom:
                seq = ''
                lastSeq = chrom
                print >> sys.stderr, 'Reading', chrom
                # Read the chromosome you just encountered
                with open(fapref + chrom + '.fa', 'r') as fasta:
                    for line in fasta: 
                        if line.startswith('>'):
                            continue
                        else:
                            seq = seq + line.strip().upper()
            print seq[start:end]
