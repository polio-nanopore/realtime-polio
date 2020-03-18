from Bio import SeqIO
from Bio import AlignIO
import sys
import argparse
import parasail

def parse_args():
    parser = argparse.ArgumentParser(description='Checking snps and generating report.')
    parser.add_argument("-i", action="store", type=str, dest="i")
    parser.add_argument("-r", action="store", type=str, dest="r")
    parser.add_argument("-o", action="store", type=str, dest="o")
    parser.add_argument("--sample", action="store", type=str, dest="sample")

    return parser.parse_args()

def get_snp_locs(traceback, outfile):

    outfile.write(f"Length of alignment is:\t{len(traceback.ref)}\n\n")
    
    snps = 0
    for i in range(len(traceback.ref)):
        col = [traceback.ref[i], traceback.query[i]]
        if len(set(col))>1:
            
            outfile.write(f"- Position {i+1}:\tReference:\t{col[0]}\tConsensus:\t{col[1]}\n")

            if 'N' not in col and 'n' not in col and '-' not in col:
                snps+=1
    return snps

def align_sequences(query,reference,o):

    nuc_matrix = parasail.matrix_create("ACGT", 2, -1)

    result_trace = parasail.nw_trace_striped_sat(query, reference, 3, 2, nuc_matrix)
    traceback = result_trace.get_traceback('|', '.', ' ')
    print(traceback.ref)
    print(traceback.comp)
    print(traceback.query)
    snps = get_snp_locs(traceback, o)
    alignment_string = ""
    for i in range(0,len(traceback.ref),60):

        alignment_string += f"Ref\t{i}\t{traceback.ref[i:i+60]}\n"

        gap = ''
        for j in range(len(str(i))):
            gap += ' ' 

        alignment_string += f"   \t{gap}\t{traceback.comp[i:i+60]}\n"
        alignment_string += f"CNS\t{i}\t{traceback.query[i:i+60]}\n"
        alignment_string += "\n"
    o.write("\n*****\n\n### Alignment with best reference\n\n")
    o.write(alignment_string + "\n*****\n")

if __name__ == '__main__':

    args = parse_args()
    cns_name = ""
    cns_seq = ""
    ref_seq = ""
    
    with open(str(args.o), "w") as fw:
        fw.write(f"## Report for {args.sample}\n\n")
        # get_snp_locs(args.i, fw)

        for record in SeqIO.parse(str(args.r),"fasta"):
            ref_seq = str(record.seq)
            fw.write(f"Closest reference in the database is:\t{record.id}\n")

        for record in SeqIO.parse(str(args.i),"fasta"):
            cns_seq = str(record.seq)
            cns_name = str(record.description)

        align_sequences(cns_seq,ref_seq, fw)
    
        
