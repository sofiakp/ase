import sys
import scipy
from scipy.stats import binom
from numpy import *

infile = sys.argv[1]

ref_idx = [0, 3]
alt_idx = [1, 4]
for line in infile:
    if line[0] == '#':
        continue
    fields = line.strip().split()
    counts1 = array([int(n) for n in fields[4].split(',')])
    counts2 = array([int(n) for n in fields[5].split(',')])
    ref_tot = sum(counts1[ref_idx] + counts2[ref_idx])
    pval = binom.cdf(counts1[ref_idx], ref_tot, 0.5)
    print line.strip() + '\t' + str(pval)
