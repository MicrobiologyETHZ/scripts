import os, sys
import numpy as np
import argparse
import itertools as it
import gzip
from Bio import Seq, SeqIO

def __main__():
    # Demultiplex Illumina reads, dual-indexed, given index and read files
    parser = argparse.ArgumentParser(description='A demultiplexer for Illumina reads')
    parser.add_argument('-r1', metavar='fwd_reads', required=True, help='Forward or single-end reads in either fastq or fastq.gz format')
    parser.add_argument('-r2', metavar='rev_reads', help='Reverse reads from paired-end sequencing in either fastq or fastq.gz format')
    parser.add_argument('-i1', metavar='fwd_index', help='Forward or single index reads in either fastq or fastq.gz format')
    parser.add_argument('-i2', metavar='rev_index', help='Reverse index from dual-indexed sequencing in either fastq or fastq.gz format')
    parser.add_argument('-ref', metavar='barcode_ref', required=True, help='Barcode reference file which should be tab-separated columnar, with sample name first, then forward index, then reverse index')
    parser.add_argument('-swap', action='store_true', help='Barcodes in the reference file are swapped, reverse index then forward index')
    parser.add_argument('-rc1', action='store_true', help='Use the reverse complement of the forward index reference sequences')
    parser.add_argument('-rc2', action='store_true', help='Use the reverse complement of the reverse index reference sequences')
    parser.add_argument('-m1', type=int, default=0, metavar='fwd_mismatch', help='Mismatches allowed in the forward index')
    parser.add_argument('-m2', type=int, default=0, metavar='rev_mismatch', help='Mismatches allowed in the reverse index')
    parser.add_argument('-nallow', action='store_true', help='Ns are allowed as mismatches, else indexes containing an N are discarded')
    parser.add_argument('-out', metavar='output_dir', default="demux", help='Directory for output')
    parser.add_argument('-v', action='store_true', help='Verbose output')
    parser.add_argument('-embed', action='store_true', help='Barcodes embedded in read names')

    args = parser.parse_args()

    if args.i2 is not None and args.r2 is None and not args.embed:
        sys.exit("No read file has been provided for dual-indexing.")

    # Import indexes and reads
    fileFormat = os.path.splitext(args.r1)[1]
    if fileFormat==".gz":
        if args.embed:
            i1gz = gzip.open(args.r1, 'rt')
            index1 = SeqIO.parse(i1gz, 'fastq')
        else:
            i1gz = gzip.open(args.i1, 'rt')
            index1 = SeqIO.parse(i1gz, 'fastq')
        if args.i2 is not None:
            i2gz = gzip.open(args.i2, 'rt')
            index2 = SeqIO.parse(i2gz, 'fastq')
        r1gz = gzip.open(args.r1, 'rt')
        reads1 = SeqIO.parse(r1gz, 'fastq')
        if args.r2 is not None:
            if args.embed:
                i2gz = gzip.open(args.r2, 'rt')
                index2 = SeqIO.parse(i2gz, 'fastq')
            r2gz = gzip.open(args.r2, 'rt')
            reads2 = SeqIO.parse(r2gz, 'fastq')
    elif fileFormat=='.fastq':
        if args.embed:
            index1 = SeqIO.parse(args.r1, 'fastq')
        else:
            index1 = SeqIO.parse(args.i1, 'fastq')
        if args.i2 is not None:
            index2 = SeqIO.parse(args.i2, 'fastq')
        reads1 = SeqIO.parse(args.r1, 'fastq')
        if args.r2 is not None:
            if args.embed:
                index2 = SeqIO.parse(args.r2, 'fastq')
            reads2 = SeqIO.parse(args.r2, 'fastq')
    else:
        sys.exit("Unrecognised sequence format.")

    # Import barcode reference file
    indexList = np.genfromtxt(args.ref, dtype=str)
    if args.swap:
        indexList = indexList[:, [0, 2, 1]]
    len1 = len(indexList[0, 1])
    len2 = len(indexList[0, 2])

    if args.rc1:
        barcodes1 = [str(Seq.Seq(bc).reverse_complement()).upper() for bc in indexList[:, 1]]
    else:
        barcodes1 = [bc.upper() for bc in indexList[:, 1]]
    if args.r2 is not None:
        if args.rc2:
            barcodes2 = [str(Seq.Seq(bc).reverse_complement()).upper() for bc in indexList[:, 2]]
        else:
            barcodes2 = [bc.upper() for bc in indexList[:, 2]]
    else:
        barcodes2 = ["N"*len2 for bc in indexList[:, 1]]

    barcodes = dict(zip(zip(barcodes1, barcodes2), indexList[:, 0]))
    if len(barcodes)<len(indexList):
        sys.exit("Two or more barcodes are identical.")

    # Generate a list of sequences n mutations from the original
    def mutate(seq, n):
        for poss in it.combinations(range(len(seq)), n):
            new = [[x] for x in seq]
            for pos in poss:
                o = seq[pos]
                if args.nallow:
                    new[pos] = [b for b in "ACGTN" if b!=o]
                else:
                    new[pos] = [b for b in "ACGT" if b!=o]
            for out in it.product(*new):
                yield ''.join(out)

    # Generate a dictionary of valid barcodes
    barcodeDict = dict()
    ambiguous_keys = list()
    for bc1, bc2 in barcodes.keys():
        for n1 in range(0, args.m1+1):
            for n2 in range(0, args.m2+1):
                for v1, v2 in it.product(mutate(bc1, n1), mutate(bc2, n2)):
                    if (v1, v2) in barcodeDict.keys():
                        ambiguous_keys.append((v1, v2))
                    else:
                        barcodeDict[(v1, v2)] = barcodes[(bc1, bc2)]
    ambiguous_keys = list(set(ambiguous_keys))
    sys.stderr.write("Allowing for errors, {} possible barcodes are ambiguous and will not be assigned to a sample.\n".format(len(ambiguous_keys)))
    for key in ambiguous_keys:
        del barcodeDict[key]
    sys.stderr.write("{} valid barcodes found.\n".format(len(barcodeDict)))
    sys.stderr.flush()

    sys.stderr.write("Assigning reads..\n")
    sys.stderr.flush()

    # Create output filestreams
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    files1 = dict()
    if args.r2 is not None:
        files2 = dict()
    if fileFormat==".gz":
        for id in barcodes.values():
            files1[id] = gzip.open(os.path.join(args.out, str(id)+"_R1.fastq.gz"), 'wt')
            if args.r2 is not None:
                files2[id] = gzip.open(os.path.join(args.out, str(id)+"_R2.fastq.gz"), 'wt')
        reject1 = gzip.open(os.path.join(args.out, "undetermined_R1.fastq.gz"), 'wt')
        if args.r2 is not None:
            reject2 = gzip.open(os.path.join(args.out, "undetermined_R2.fastq.gz"), 'wt')
    else:
        for id in barcodes.values():
            files1[id] = open(os.path.join(args.out, str(id)+"_R1.fastq"), 'w')
            if args.r2 is not None:
                files2[id] = open(os.path.join(args.out, str(id)+"_R2.fastq"), 'w')
        reject1 = open(os.path.join(args.out, "undetermined_R1.fastq"), 'w')
        if args.r2 is not None:
            reject2 = open(os.path.join(args.out, "undetermined_R2.fastq"), 'w')

    # Output reads to the correct files where possible, otherwise they go to undetermined
    rejected = 0
    total = 0
    if args.r2 is not None:
        if args.i2 is not None or args.embed:
            ir = it.zip_longest(index1, reads1, index2, reads2)
        else:
            ir = it.zip_longest(index1, reads1, [SeqIO.SeqRecord("N"*len2)], reads2, fillvalue=SeqIO.SeqRecord("N"*len2))
    else:
        ir = it.zip_longest(index1, reads1, [SeqIO.SeqRecord("N"*len2)], [SeqIO.SeqRecord("N"*len2)], fillvalue=SeqIO.SeqRecord("N"*len2))

    for i1, r1, i2, r2 in ir:
        total+=1
        try:
            if args.embed:
                bfwd = (i1.description.split(":")[-1]).split("+")[0]
                brev = (i2.description.split(":")[-1]).split("+")[1]
            else:
                bfwd = str(i1.seq)
                brev = str(i2.seq)
            files1[barcodeDict[(bfwd[0:len1], brev[0:len2])]].write(r1.format('fastq'))
            if args.r2 is not None:
                files2[barcodeDict[(bfwd[0:len1], brev[0:len2])]].write(r2.format('fastq'))
        except(KeyError):
            reject1.write(r1.format('fastq'))
            if args.r2 is not None:
                reject2.write(r2.format('fastq'))
            rejected+=1
        if args.v and not total%1000:
            sys.stderr.write(str(total-rejected)+"/"+str(total)+"\n")

    for id in barcodes.values():
        files1[id].close()
        if args.r2 is not None:
            files2[id].close()
    reject1.close()
    if args.r2 is not None:
        reject2.close()

    sys.stderr.write("\n"+str(total-rejected)+"/"+str(total)+" barcodes identified.\n")
    sys.stderr.flush()

if __name__ == "__main__":
    __main__()
