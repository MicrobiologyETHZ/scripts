import argparse, os, subprocess
from Bio import Seq, SeqIO

if True:
#def __main__():
    parser = argparse.ArgumentParser(description='Extract nucleotide and/or protein sequences from GlimmerHMM output')
    parser.add_argument('contigs',metavar='contigs_file',help='File containing contigs in fasta format')
    parser.add_argument('gff',metavar='gff_file',help='GFF file with gene predictions')
    parser.add_argument('-p','--prefix',metavar='prefix',help='Output file prefix')

    args = parser.parse_args()
    if not args.prefix:
        args.prefix = os.path.splitext(args.contigs)[0]

    # Import sequences
    sequences = SeqIO.to_dict(SeqIO.parse(args.contigs,'fasta'))

    # Manually parse GFF file to construct sequences
    mrnas = {}
    name = ""
    strand = "+"
    subseq = Seq.Seq("")
    with open(args.gff,'r') as gff:
        for record in gff.readlines():
            if record[0] != "#":
                fields = record.split("\t")
                if fields[2] == "mRNA":
                    if strand == "+":
                        mrnas[name] = subseq 
                    else:
                        mrnas[name] = subseq.reverse_complement()
                    name = "{}_{}:{}".format(fields[0],fields[3],fields[4])
                    strand = fields[6]
                    subseq = Seq.Seq("") 
                else:
                    subseq = subseq + sequences[fields[0]].seq[(int(fields[3])-1):int(fields[4])]

    # Output
    with open(args.prefix+".fna",'w') as fo:
        for k,v in mrnas.items():
            fo.write(">{}\n{}\n".format(k,v))
    with open(args.prefix+".faa",'w') as fo:
        for k,v in mrnas.items():
            fo.write(">{}\n{}\n".format(k,v.translate()))

#if __name__ == "__main__":
#    __main__()

