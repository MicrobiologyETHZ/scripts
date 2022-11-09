# Useful Scripts

| Script | Description |
| ------ | ----------- |
| ab12fastq.py | Converts ab1 Sanger sequencing files to fastq |
| assembly_stats.py | Summarises (meta)genome assembly statistics |
| demultiplex.py | Demultiplex illumina reads |

## ab12fastq.py

```
usage: ab12fastq.py [-h] [--separate] files [files ...]

Convert ab1 Sanger sequencing files to fastq

positional arguments:
  files       One or more ab1 format files to convert

optional arguments:
  -h, --help  show this help message and exit
  --separate  Output to individual files instead of stdout
```

## assembly_stats.py

```
usage: assembly_stats.py [-h] [-c cutoffs [cutoffs ...]] [-p plot_file_prefix] [-l] contigs_files [contigs_files ...]

Calculate basic assembly statistics.

positional arguments:
  contigs_files         Contigs or scaffolds file(s) in fasta format

optional arguments:
  -h, --help            show this help message and exit
  -c cutoffs [cutoffs ...], --cutoffs cutoffs [cutoffs ...]
                        Size cutoff(s) below which contigs are discarded
  -p plot_file_prefix, --plot plot_file_prefix
                        Output various plots
  -l, --log             Set the y-axis to log scale
```

## demultiplex.py

```
usage: demultiplex.py [-h] -r1 fwd_reads [-r2 rev_reads] [-i1 fwd_index] [-i2 rev_index] -ref barcode_ref [-swap] [-rc1] [-rc2] [-m1 fwd_mismatch] [-m2 rev_mismatch] [-nallow] [-out output_dir] [-v]
                      [-embed]

A demultiplexer for Illumina reads

optional arguments:
  -h, --help        show this help message and exit
  -r1 fwd_reads     Forward or single-end reads in either fastq or fastq.gz format
  -r2 rev_reads     Reverse reads from paired-end sequencing in either fastq or fastq.gz format
  -i1 fwd_index     Forward or single index reads in either fastq or fastq.gz format
  -i2 rev_index     Reverse index from dual-indexed sequencing in either fastq or fastq.gz format
  -ref barcode_ref  Barcode reference file which should be tab-separated columnar, with sample name first, then forward index, then reverse index
  -swap             Barcodes in the reference file are swapped, reverse index then forward index
  -rc1              Use the reverse complement of the forward index reference sequences
  -rc2              Use the reverse complement of the reverse index reference sequences
  -m1 fwd_mismatch  Mismatches allowed in the forward index
  -m2 rev_mismatch  Mismatches allowed in the reverse index
  -nallow           Ns are allowed as mismatches, else indexes containing an N are discarded
  -out output_dir   Directory for output
  -v                Verbose output
  -embed            Barcodes embedded in read names
```
