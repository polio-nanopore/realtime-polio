from Bio import SeqIO
from Bio import AlignIO
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Checking snps and generating report.')
    parser.add_argument("-i", action="store", type=str, dest="i")
    parser.add_argument("-o", action="store", type=str, dest="o")
    parser.add_argument("--sample", action="store", type=str, dest="sample")


    return parser.parse_args()

def get_snp_locs(aln, outfile):
    alignment = AlignIO.read(aln, "fasta")
    outfile.write(f"Closest reference in the database is:\t{alignment[0].id}\n")
    print(f"Reading in {aln}.")
    outfile.write(f"Length of alignment is:\t{len(alignment[0])}\n\n")
    
    snps = 0
    for i in range(len(alignment[0])):
        col = alignment[:,i]
        if len(set(col))>1:
            
            fw.write(f"- Position {i+1}:\tReference:\t{col[0]}\tConsensus:\t{col[1]}\n")

            if 'N' not in col and 'n' not in col and '-' not in col:
                snps+=1

    fw.write(f"\nTotal number of snps from {alignment[0].id}:\t{snps}\n\n")

if __name__ == '__main__':

    args = parse_args()
    with open(str(args.o), "w") as fw:
        fw.write(f"## Report for {args.sample}\n")
        get_snp_locs(args.i, fw)
        
