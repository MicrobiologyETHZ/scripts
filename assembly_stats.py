# Produces basic stats on assemblies with optional graph output

import os,sys
import argparse
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Calculate basic assembly statistics.')
parser.add_argument('files',nargs='+',metavar='contigs_files',help='Contigs or scaffolds file(s) in fasta format')
parser.add_argument('-c','--cutoffs',nargs='+',metavar='cutoffs',type=int,help='Size cutoff(s) below which contigs are discarded')
parser.add_argument('-p','--plot',metavar='plot_file_prefix',help='Output various plots')
parser.add_argument('-l','--log',action='store_true',help='Set the y-axis to log scale')

args = parser.parse_args()

# Import contigs
allContigs = dict()
for file in args.files:
	allContigs[file] = list(SeqIO.parse(file,'fasta'))

def getStats(contigs):
	lengths = [len(x.seq) for x in contigs]
	lengths = np.sort(lengths)
	cumsum = np.cumsum(lengths)
	num = len(lengths)
	tot = sum(lengths)
	n50 = [lengths[i] for i,x in enumerate(cumsum) if x>(tot*0.5)][0]
	n90 = [lengths[i] for i,x in enumerate(cumsum) if x>(tot*0.9)][0]
	max = np.max(lengths)

	stats = [lengths,cumsum,num,tot,n50,n90,max]
	return(stats)

# Get stats for each set of contigs for each cutoff
allStats = dict()
for file in args.files:
	allStats[file+str(0)] = getStats(allContigs[file])
	if args.cutoffs is not None:
		for cutoff in args.cutoffs:
			subset = [contig for contig in allContigs[file] if len(contig.seq)>cutoff]
			allStats[file+str(cutoff)] = getStats(subset)

# Text output
print("File\tCutoff\tCount\tLength\tN50\tN90\tMax")
for file in args.files:
	stats = allStats[file+str(0)]
	print(file+'\t'+str(0)+'\t'+'\t'.join(str(x) for x in stats[2:]))
if args.cutoffs is not None:
	for cutoff in args.cutoffs:
		for file in args.files:
			stats = allStats[file+str(cutoff)]
			print(file+'\t'+str(cutoff)+'\t'+'\t'.join(str(x) for x in stats[2:]))

# Graphical output

if args.plot:
	# Contig size vs. cumulative size
	xmax = np.max([x[6] for x in allStats.values()])
	ymax = np.max([x[3] for x in allStats.values()])
	for file in args.files:
		stats = allStats[file+str(0)]
		plt.xlim(10,xmax*1.05)
		plt.ylim(0,ymax*1.05)
		if args.log:
			plt.loglog(stats[0],stats[1],'-')
		else:
			plt.semilogx(stats[0],stats[1],'-')
	if args.cutoffs is not None:
		for cutoff in args.cutoffs:
			plt.plot([cutoff,cutoff],plt.ylim(),'k--')
	plt.xlabel("Contig size (bp)")
	plt.ylabel("Cumulative metagenome size (bp)")
	plt.legend(args.files,prop={'size':6})
	plt.savefig(args.plot+"_cum.png",dpi=200)
	plt.close()

	# Contig size histogram
	ymax = np.max([x[2] for x in allStats.values()])
	xmin = np.min([x[0][0] for x in allStats.values()])
	for file in args.files:
		stats = allStats[file+str(0)]
		plt.xlim(xmin*0.95,xmax*1.05)
		plt.ylim(0.1,ymax)
		counts,breaks = np.histogram(stats[0],bins=np.logspace(np.log10(xmin),np.log10(xmax),21))
		centers = [np.sqrt(breaks[x]*breaks[x+1]) for x in range(0,len(breaks)-1)]
		plt.loglog(centers,counts,"-")
	if args.cutoffs is not None:
                for cutoff in args.cutoffs:
                        plt.plot([cutoff,cutoff],plt.ylim(),'k--')
	plt.xlabel("Contig size (bp)")
	plt.ylabel("Frequency")
	plt.legend(args.files,prop={'size':6})
	plt.savefig(args.plot+"_hist.png",dpi=200)
	plt.close()
