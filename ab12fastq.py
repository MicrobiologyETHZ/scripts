import argparse, os, sys
from Bio import SeqIO

def __main__():
    parser = argparse.ArgumentParser(description='Convert ab1 Sanger sequencing files to fastq')
    parser.add_argument('files', nargs='+', help='One or more ab1 format files to convert')
    parser.add_argument('--separate', action='store_true', help='Output to individual files instead of stdout')
    args = parser.parse_args()

    for file in args.files:
        record = SeqIO.parse(file, "abi")
        if args.separate:
            path_noext = os.path.splitext(file)[0]
            SeqIO.write(record, f'{path_noext}.fastq', "fastq")
        SeqIO.write(record, sys.stdout, "fastq")

if __name__ == "__main__":
    __main__()

