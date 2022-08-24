#!/usr/bin/env python3

# yupf05@gmail.com
# modified by Andy Tock (ajt200@cam.ac.uk) for use with fastq.gz files
### AT's changes based on https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
### and https://bioinformatics.stackexchange.com/questions/892/how-do-you-write-a-gz-fastq-file-with-biopython/894

# Usage:
# deduplicate_gz_python3.py ./data/libName_R1.fastq.gz ./data/libName_R2.fastq.gz ./data/dedup/libName_R1_dedup.fastq.gz ./data/dedup/libName_R2_dedup.fastq.gz

from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
import os,sys,argparse
import gzip

def ParseArg():
    p=argparse.ArgumentParser(description = 'Remove duplicated reads which have same sequences for both forward and reverse reads. Choose the one appears first.', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('input1', type = str, metavar = 'reads1', help = 'forward input fastq/fasta file')
    p.add_argument('input2', type = str, metavar = 'reads2', help = 'reverse input fastq/fasta file')
    p.add_argument('output1', type = str, metavar = 'reads1', help = 'deduplicated forward output fastq/fasta file')
    p.add_argument('output2', type = str, metavar = 'reads2', help = 'deduplicated reverse output fastq/fasta file')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def Main():

    unique_seqs = set()
    args = ParseArg()
    fileName1 = os.path.basename(args.input1)
    fileName2 = os.path.basename(args.input2)
    indexOfDot1 = fileName1.index(".")
    indexOfDot2 = fileName2.index(".")
    filePrefix1 = fileName1[:indexOfDot1]
    filePrefix2 = fileName2[:indexOfDot2]
    dataPath = os.path.dirname(os.path.realpath(args.input1))
    dedupPath = dataPath+"/dedup/"
    #try:
    #    os.mkdir(dedupPath)
    #except OSError:
    #    print("Creation of the directory %s failed" % dedupPath)
    #else:
    #    print("Successfully created the directory %s " % dedupPath) 
    outfile1 = gzip.open(args.output1, "wt")
    outfile2 = gzip.open(args.output2, "wt")
    #outfile1 = gzip.open(dedupPath+filePrefix1+"_dedup.fastq.gz", "wt")
    #outfile2 = gzip.open(dedupPath+filePrefix2+"_dedup.fastq.gz", "wt")
    fastq_iter1 = SeqIO.parse(gzip.open(args.input1, "rt"), "fastq")
    fastq_iter2 = SeqIO.parse(gzip.open(args.input2, "rt"), "fastq")
    for rec1, rec2 in zip(fastq_iter1, fastq_iter2):
        if str((rec1+rec2).seq) not in unique_seqs:
                outfile1.write(rec1.format("fastq"))
                outfile2.write(rec2.format("fastq"))
                unique_seqs.add(str((rec1+rec2).seq))
    outfile1.close()
    outfile2.close()

if __name__ == '__main__':
    Main()
