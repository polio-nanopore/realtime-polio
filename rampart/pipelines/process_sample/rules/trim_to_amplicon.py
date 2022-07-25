0import csv
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import parasail
# Need to adapt to trimming more specifically


def parse_args():
    parser = argparse.ArgumentParser(description='Trim primers.')
    parser.add_argument("--reads", action="store", type=str, dest="reads")

    parser.add_argument("--output_reads", action="store", type=str, dest="output_reads")
    parser.add_argument("--substitution_matrix", action="store", type=str, dest="substitution_matrix")

#    parser.add_argument("--forward", action="store",type=str, dest="forward", default="TTTAACCTGGGTTTGTGTCAGCCTGTAATGA")
#    parser.add_argument("--reverse", action="store",type=str, dest="reverse",default="TACACCTTRTCTCTGGAGAATCCAATTT")

    parser.add_argument("--forward", action="store",type=str, dest="forward", default="TGGCGGAACCGACTACTTTGGGTG")
    parser.add_argument("--reverse", action="store",type=str, dest="reverse",default="TCAATACGGTGTTTGCTCTTGAACTG")

    return parser.parse_args()

def get_best_reference(query, ref_dict, matrix):
    best_reference_alignment = {
        "reference": "None_NA",
        "identity": 0,
        "coverage": 0
        }
    for ref in ref_dict:
        
        new_reference_alignment = align_read(query, ref, ref_dict[ref], matrix)
        if new_reference_alignment["coverage"] > 0.7:
            if best_reference_alignment["identity"] < new_reference_alignment["identity"]: 
                best_reference_alignment = new_reference_alignment
            else:
                pass
    return best_reference_alignment

def align_read(read, primers, matrix, gap_open=3, gap_extension=2):
    result_trace = None
    traceback = None
    result_trace = parasail.sw_trace_striped_sat(query, primers[0], gap_open, gap_extension, matrix)
    profile = parasail.ssw_init("asdf", parasail.blosum62, score_size)
    result = parasail.ssw_profile(profile, "asdf", 10, 1)
    traceback = result_trace.get_traceback('|', '.', ' ')
    # print(traceback.ref)
    # print(traceback.comp)
    # print(traceback.query)
    query_start = result_trace.cigar.beg_query
    reference_start = result_trace.cigar.beg_ref
    # cigar = result.cigar.decode.decode("UTF-8")
    result_stats = None
    result_stats = parasail.sw_stats_striped_sat(query, reference, gap_open, gap_extension, matrix)
    alignment_covers = int(result_stats.length) / len(primer[0])
    # print(ref_id, len(reference), alignment_covers, result_stats.matches, len(reference))
    return {
            "reference":ref_id,
            "query_start": query_start,
            "reference_start": reference_start,
            "matches": result_stats.matches,
            "coverage": alignment_covers,
            "aln_len": result_stats.length,
            "len": result_stats.len_ref,
            "identity": result_stats.matches / result_stats.len_ref,
            "ref": traceback.ref,
            "comp": traceback.comp,
            "query": traceback.query
        }


def process_file(reads,primers,nuc_matrix):

    counts = Counter()
    # error = collections.defaultdict(list)

    for record in SeqIO.parse(reads, "fastq"):
        read_f = str(record.seq)
        read_rc = str(record.seq.reverse_complement())
        record_seqs = {"read_forward": read_f, "read_reverse": read_rc}

        # print("*****")
        # print("the record id is",record.id)
        stats = best_alignment(primers, record_seqs, nuc_matrix)

        best_ref,direction = stats["reference"].split("_")

        background_error_rate = get_background_error_rate(stats)

        # output_file.write('{},{},{},{},{},{},{}\n'.format(sample,best_ref,direction,stats["identity"],stats["query_start"],stats["reference_start"],background_error_rate))
        alignment_covers = int(stats["aln_len"]) / int(stats["len"]) # doesn't account for gaps so can be > 1

        if stats["identity"] > 0.75:

            read_seq = ''
            if direction == "forward":
                read_seq = str(record.seq)
            else:
                read_seq = str(record.seq.reverse_complement())
            ref_seq = references[best_ref + "_forward"]
            
            alignment = align_read(read_seq, best_ref, ref_seq, nuc_matrix)
            
            # print('\n', alignment["reference"], '\n', alignment["query"],'\n', alignment["comp"],'\n', alignment["ref"])
            # print(best_ref,'\t', direction, '\t',alignment_covers, '\t', alignment["len"], '\t',alignment["aln_len"], '\t',alignment["identity"],'\t', alignment["query_start"],'\t', alignment["reference_start"])



def get_read_list(read):

    read_list = []
    for record in SeqIO.parse(ref_file, "fasta"):
         references[record.id + "_forward"] = str(record.seq)
         references[record.id + "_reverse"] = str(record.seq.reverse_complement())
    return references


def file_writer(reads,output_reads,trim):
    trim = int(trim)
    with open(output_reads, "w") as handle:
        for title, seq, qual in FastqGeneralIterator(open(reads)):
            handle.write(f"@{title}\n{seq[trim:-trim]}\n+\n{qual[trim:-trim]}\n")


if __name__ == '__main__':

    args = parse_args()

    primers = (args.forward, args.reverse)


    fw = open(str(args.output_reads),"w")

    nuc_matrix = parasail.Matrix(str(args.substitution_matrix))

    process_file(str(args.reads), primers, nuc_matrix)

    file_writer(args.reads,args.output_reads,58)
