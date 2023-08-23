import argparse
import time
from collections import Counter
import numpy as np
import sys

start = time.time()

parser = argparse.ArgumentParser(description='Calculates the horizontal and vertical coverage from a paf alignment file')
parser.add_argument('paf', help='Path to the .paf file.')
parser.add_argument('-q', type=int, default = 0, help='Quality score cutoff. [0]')
#parser.add_argument('-c', '--circular_ref', type = int, default = 0, help = "extra basepairs that have been pasted at the end and beginning of the sequence to account for the circularity of the DNA")
#parser.add_argument('-o','--out_dir', type=str, default='/nfs/home/mibohl/project/scratch', help='path the output directory, default: ./scratch')
#parser.add_argument('-p','--percent', action='store_true', default=False, help='return the vertical coverage of all covered bases in % (True or False)')
args = parser.parse_args()

coverages = {}
lengths = {}

with open(args.paf, 'r') as fi:
    for line in fi:
        fields = line.strip().split('\t')
        if int(fields[11]) >= args.q:
            if fields[5] not in coverages.keys():
                coverages[fields[5]] = Counter()
                lengths[fields[5]] = int(fields[6])
            coverages[fields[5]][int(fields[7])] += 1
            coverages[fields[5]][int(fields[8])] -= 1

results = {}
for seq in coverages.keys():
    d_cov = [coverages[seq][i] for i in range(1,lengths[seq]+1)]
    cov = np.cumsum(d_cov)
    h = np.count_nonzero(cov)/lengths[seq]
    v = np.sum(cov)/lengths[seq]
    results[seq] = (h, v)

for seq in results.keys():
    sys.stdout.write(f'{seq}\t{results[seq][0]}\t{results[seq][1]}\n')

sys.stderr.write(f'Completed in {time.time()-start}s\n')
