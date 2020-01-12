import os
from Bio import SeqIO
import collections
import sys
import argparse
import datetime
import subprocess

parser = argparse.ArgumentParser(description='Generate blast report.')
parser.add_argument("--consensus", action="store", type=str, dest="consensus")
parser.add_argument("--blast_db", action="store", type=str, dest="blast_db")
parser.add_argument("--detailed_blast_db", action="store", type=str, dest="detailed_blast_db")
parser.add_argument("--blast_file", action="store", type=str, dest="blast_file")
parser.add_argument("--detailed_blast_file", action="store", type=str, dest="detailed_blast_file")
parser.add_argument("--output_report", action="store", type=str, dest="output_report")
parser.add_argument("--output_seqs", action="store", type=str, dest="output_seqs")
parser.add_argument("--sample", action="store", type=str, dest="sample")
args = parser.parse_args()

def parse_blast_for_top_hit(blast_csv):
    hits = []
    with open(blast_csv, "r") as f:
        for l in f:
            tokens= l.rstrip('\n').split(',')
            query = tokens[0]
            hit = tokens[1]
            score= tokens[-1] 
            hits.append((hit,score,l.rstrip('\n')))

    top = sorted(hits, key = lambda x : float(x[1]), reverse = True)[0]
    return(top[2].split(','))

top_hit = parse_blast_for_top_hit(args.blast_file)
detailed_top_hit = parse_blast_for_top_hit(args.detailed_blast_file)

query,subject,pid,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore= top_hit
query,subject_detailed,pid_detailed,length_detailed,mismatch_detailed,gapopen_detailed,qstart_detailed,qend_detailed,sstart_detailed,send_detailed,evalue_detailed,bitscore_detailed= detailed_top_hit
seq_dict={
"Sabin1_vacc":"GGGTTAGGTCAGATGCTTGAAAGCATGATTGACAACACAGTCCGTGAAACGGTGGGGGCGGCAACGTCTAGAGACGCTCTCCCAAACACTGAAGCCAGTGGACCAGCACACTCCAAGGAAATTCCGGCACTCACCGCAGTGGAAACTGGGGCCACAAATCCACTAGTCCCTTCTGATACAGTGCAAACCAGACATGTTGTACAACATAGGTCAAGGTCAGAGTCTAGCATAGAGTCTTTCTTCGCGCGGGGTGCATGCGTGGCCATTATAACCGTGGATAACTCAGCTTCCACCAAGAATAAGGATAAGCTATTTACAGTGTGGAAGATCACTTATAAAGATACTGTCCAGTTACGGAGGAAATTGGAGTTCTTCACCTATTCTAGATTTGATATGGAATTTACCTTTGTGGTTACTGCAAATTTCACTGAGACTAACAATGGGCATGCCTTAAATCAAGTGTACCAAATTATGTACGTACCACCAGGCGCTCCAGTGCCCGAGAAATGGGACGACTACACATGGCAAACCTCATCAAATCCATCAATCTTTTACACCTACGGAACAGCTCCAGCCCGGATCTCGGTACCGTATGTTGGTATTTCGAACGCCTATTCACACTTTTACGACGGTTTTTCCAAAGTACCACTGAAGGACCAGTCGGCAGCACTAGGTGACTCCCTCTATGGTGCAGCATCTCTAAATGACTTCGGTATTTTGGCTGTTAGAGTAGTCAATGATCACAACCCGACCAAGGTCACCTCCAAAATCAGAGTGTATCTAAAACCCAAACACATCAGAGTCTGGTGCCCGCGTCCACCGAGGGCAGTGGCGTACTACGGCCCTGGAGTGGATTACAAGGATGGTACGCTTACACCCCTCTCCACCAAGGATCTGACCACATAT",
"Sabin2_vacc":"GGAATTGGTGACATGATTGAGGGGGCCGTTGAAGGGATTACTAAAAATGCATTGGTTCCCCCGACTTCCACCAATAGCCTGCCTGACACAAAGCCGAGCGGTCCAGCCCACTCCAAGGAGATACCTGCATTGACAGCCGTGGAGACAGGGGCTACCAATCCGTTGGTGCCTTCGGACACCGTGCAAACGCGCCATGTCATCCAGAGACGAACGCGATCAGAGTCCACGGTTGAGTCATTCTTTGCAAGAGGGGCTTGCGTGGCTATCATTGAGGTGGACAATGATGCACCGACAAAGCGCGCCAGCAGATTGTTTTCGGTTTGGAAAATAACTTACAAAGATACTGTTCAACTGAGACGCAAACTGGAATTTTTCACATATTCGAGATTTGACATGGAGTTCACTTTTGTGGTCACCTCAAACTACATTGATGCAAATAACGGACATGCATTGAACCAAGTTTATCAGATAATGTATATACCACCCGGAGCACCTATCCCTGGTAAATGGAATGACTATACGTGGCAGACGTCCTCTAACCCGTCGGTGTTTTACACCTATGGGGCGCCCCCAGCAAGAATATCAGTGCCCTACGTGGGAATTGCTAATGCGTATTCCCACTTTTATGATGGGTTTGCAAAAGTACCACTAGCGGGTCAAGCCTCAACTGAAGGCGATTCGTTGTACGGTGCTGCCTCACTGAATGATTTTGGATCACTGGCTGTTCGCGTGGTAAATGATCACAACCCCACGCGGCTCACCTCCAAGATCAGAGTGTACATGAAGCCAAAGCATGTCAGAGTCTGGTGCCCACGACCTCCACGAGCAGTCCCATACTTCGGACCAGGTGTTGATTATAAAGATGGGCTCACCCCACTACCAGAAAAGGGATTAACGACTTAT",
"Sabin3_vacc":"ATTGAAGATTTGACTTCTGAAGTTGCACAGGGCGCCCTAACTTTGTCACTCCCGAAGCAACAGGATAGCTTACCTGATACTAAGGCCAGTGGCCCGGCGCATTCCAAGGAGGTACCTGCACTCACTGCAGTCGAGACTGGAGCCACCAATCCTCTGGCACCATCCGACACAGTTCAAACGCGCCACGTAGTCCAACGACGCAGCAGGTCAGAGTCCACAATAGAATCATTCTTCGCACGCGGGGCGTGCGTCGCTATTATTGAGGTGGACAATGAACAACCAACCACCCGGGCACAGAAACTATTTGCCATGTGGCGCATTACATACAAAGATACAGTGCAGTTGCGCCGTAAGTTGGAGTTTTTCACATACTCTCGTTTTGACATGGAATTCACCTTCGTGGTAACCGCCAACTTCACCAACGCTAATAATGGGCATGCACTCAACCAGGTGTACCAGATAATGTACATCCCCCCAGGGGCACCCACACCAAAGTCATGGGACGACTACACTTGGCAAACATCTTCCAACCCGTCCATATTTTACACCTATGGGGCTGCCCCGGCGCGAATCTCAGTGCCATACGTGGGGTTAGCCAATGCTTACTCGCACTTTTACGACGGCTTCGCCAAGGTGCCATTGAAGACAGATGCCAATGACCAGATTGGTGATTCCTTGTACAGCGCCATGACAGTTGATGACTTTGGTGTATTGGCAGTTCGTGTTGTCAATGATCACAACCCCACTAAAGTAACCTCCAAAGTCCGCATTTACATGAAACCCAAACACGTACGTGTCTGGTGCCCTAGACCGCCGCGCGCGGTACCTTATTATGGACCAGGGGTGGACTATAGGAACAACTTGGACCCCTTATCTGAGAAAGGTTTGACCACATAT"
}
for record in SeqIO.parse(args.blast_db,"fasta"):
    if record.id==subject:
        seq_dict[record.id]=record.seq
if subject!=subject_detailed:
    for record in SeqIO.parse(args.detailed_blast_db,"fasta"):
        if record.id==subject_detailed:
            seq_dict[record.id]=record.seq
for record in SeqIO.parse(args.consensus,"fasta"):
    seq_dict[record.id]=record.seq

fw= open(args.output_seqs,"w")
for i in seq_dict:
    fw.write(">{}\n{}\n".format(i, seq_dict[i]))
fw.close()

fw=open(args.output_report,"w")
fw.write("## Report for sample: {}\n\n".format(args.sample))

now=datetime.datetime.now()
fw.write("{}/{}/{}\n\n".format(now.year,now.month,now.day))

fw.write("The type of virus that sample {} sequence is most similar to is {}.\n\n".format(args.sample, subject.split('_')[0]))
fw.write("Percentage ID: {}\n\nLength: {}\n\nNo. Mismatch: {}\n\nNo. Gap Opening: {}\n\nE-value: {}\n\n".format(pid,length,mismatch,gapopen,evalue))

fw.write("The closest hit in the detailed database to sample {} is {}.\n\n".format(args.sample, subject_detailed))
fw.write("Percentage ID: {}\n\nLength: {}\n\nNo. Mismatch: {}\n\nNo. Gap Opening: {}\n\nE-value: {}\n\n".format(pid_detailed,length_detailed,mismatch_detailed,gapopen_detailed,evalue_detailed))

fw.close()