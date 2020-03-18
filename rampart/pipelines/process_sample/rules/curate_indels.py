from Bio import SeqIO
from Bio import AlignIO
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Clean up coding consensus.')
    parser.add_argument("-a", action="store", type=str, dest="alignment")
    parser.add_argument("-o", action="store", type=str, dest="output_seq")
    parser.add_argument("-r", action="store", type=str, dest="round")
    parser.add_argument("--name", action="store", type=str, dest="name")

    return parser.parse_args()

def trim_trailing_gaps(alignment):
    start_position = 0
    end_position = 0
    for col in range(alignment.get_alignment_length()):
        if not "-" in alignment[:, col]:
            start_position = col
            break

    for col in range(alignment.get_alignment_length()):
        end_index = col+1
        if not "-" in alignment[:, -end_index]:
            end_position = col
            break

    print(f"\nTrimming trailing gaps in alignment.\nAlignment now from {start_position} to {alignment.get_alignment_length()-end_position}.\n")   

    if end_position == 0:
        return alignment[:,start_position:]
    else:
        return alignment[:,start_position:-end_position]

def remove_gaps(col):    
    if len(set(col))>1 and '-' in col:
        if col[0] == '-':
            return ""
        else:
            return 'N'


def check_frame_in_window(window):
    frame = 0
    last_position = 0

    for i in range(len(window[0])):
        
        col = window[:,i]
        
        if len(set(col))>1 and '-' in col:

            if '-' == col[0]:
                frame +=1

            elif col[1] == '-':
                frame -=1

            last_position = i
    if frame % 3 == 0: 
        return True
    else:
        return False

def find_gaps(aln):

    trimmed = trim_trailing_gaps(aln)

    print(f"Reading in {aln}.\n\nGaps found:")
    indels_to_remove = []
    for i in range(len(trimmed[0])):
        col = trimmed[:,i]
        if len(set(col))>1:
            print(i, col)
        if len(set(col))>1 and '-' in col:

            print(f"Position {i+1}:\tReference:\t{col[0]}\tConsensus:\t{col[1]}")
            
            window = trimmed[:,i-12:i+12]
            if not check_frame_in_window(window):
                indels_to_remove.append(i)
                print("Will remove.")
            else:
                print("Frame maintained.")

    cns_string = ""
 
    for i in indels_to_remove:
        col = trimmed[:,i]
        remove_gaps(col)

    cns_string = ""
    for i in range(len(trimmed[0])):
        col = trimmed[:,i]
        if i in indels_to_remove:
            print(f"Removing position {i+1}:\tReference:\t{col[0]}\tConsensus:\t{col[1]}")
            cns_string += remove_gaps(col)
        else:
            cns_string+= col[1]

    return trimmed[1].id, cns_string.replace("-","")

if __name__ == '__main__':

    args = parse_args()

    round_name = ''
    if args.round:
        round_name = f" round_name={args.round}"
    #the rule is to replace a gap in the query with 'N' and to force delete a base that causes a gap in the reference
    
    with open(args.output_seq, "w") as fw:

        
        alignment = AlignIO.read(str(args.alignment), "fasta")

        cns_id, new_consensus = find_gaps(alignment)
        
        fw.write(f">{args.name} accession={cns_id.split(':')[0]}{round_name} length={len(new_consensus)}\n{new_consensus.upper()}\n")
