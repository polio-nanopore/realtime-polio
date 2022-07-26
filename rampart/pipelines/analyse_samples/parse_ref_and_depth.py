import argparse
from Bio import SeqIO
from collections import OrderedDict
from collections import Counter
from collections import defaultdict
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
import sys
import csv

def parse_args():
    parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

    parser.add_argument("--csv", action="store", type=str, dest="csv")
    parser.add_argument("--reads", action="store", type=str, dest="reads")
    parser.add_argument("--references", action="store", type=str, dest="references")
    parser.add_argument("--out_counts", action="store", type=str, dest="out_counts")
    parser.add_argument("--sample", action="store", type=str, dest="sample")
    parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
    parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

    parser.add_argument("--output_path", action="store", type=str, dest="output_path")

    return parser.parse_args()

def make_ref_dict(references):
    refs = {}
    for record in SeqIO.parse(references,"fasta"):
        refs[record.id]=record.seq
    return refs

def count_and_return_analysis_dict(report,csv_out,sample):

    counts = OrderedDict()
    counts["Enterovirus_D111"]=0
    counts["Coxsackievirus_A1"]=0
    counts["Coxsackievirus_A10"]=0
    counts["Coxsackievirus_A11"]=0
    counts["Coxsackievirus_A12"]=0
    counts["Coxsackievirus_A13"]=0
    counts["Coxsackievirus_A14"]=0
    counts["Coxsackievirus_A15"]=0
    counts["Coxsackievirus_A16"]=0
    counts["Coxsackievirus_A17"]=0
    counts["Coxsackievirus_A18"]=0
    counts["Coxsackievirus_A19"]=0
    counts["Coxsackievirus_A2"]=0
    counts["Coxsackievirus_A20"]=0
    counts["Coxsackievirus_A21"]=0
    counts["Coxsackievirus_A22"]=0
    counts["Coxsackievirus_A24"]=0
    counts["Coxsackievirus_A3"]=0
    counts["Coxsackievirus_A4"]=0
    counts["Coxsackievirus_A5"]=0
    counts["Coxsackievirus_A6"]=0
    counts["Coxsackievirus_A7"]=0
    counts["Coxsackievirus_A8"]=0
    counts["Coxsackievirus_A9"]=0
    counts["Coxsackievirus_B1"]=0
    counts["Coxsackievirus_B2"]=0
    counts["Coxsackievirus_B3"]=0
    counts["Coxsackievirus_B4"]=0
    counts["Coxsackievirus_B5"]=0
    counts["Coxsackievirus_B6"]=0
    counts["Echovirus_E1"]=0
    counts["Echovirus_E11"]=0
    counts["Echovirus_E12"]=0
    counts["Echovirus_E13"]=0
    counts["Echovirus_E14"]=0
    counts["Echovirus_E15"]=0
    counts["Echovirus_E16"]=0
    counts["Echovirus_E17"]=0
    counts["Echovirus_E18"]=0
    counts["Echovirus_E19"]=0
    counts["Echovirus_E2"]=0
    counts["Echovirus_E20"]=0
    counts["Echovirus_E21"]=0
    counts["Echovirus_E24"]=0
    counts["Echovirus_E25"]=0
    counts["Echovirus_E26"]=0
    counts["Echovirus_E27"]=0
    counts["Echovirus_E29"]=0
    counts["Echovirus_E3"]=0
    counts["Echovirus_E30"]=0
    counts["Echovirus_E31"]=0
    counts["Echovirus_E32"]=0
    counts["Echovirus_E33"]=0
    counts["Echovirus_E4"]=0
    counts["Echovirus_E5"]=0
    counts["Echovirus_E6"]=0
    counts["Echovirus_E7"]=0
    counts["Echovirus_E9"]=0
    counts["Enterovirus_A114"]=0
    counts["Enterovirus_A119"]=0
    counts["Enterovirus_A120"]=0
    counts["Enterovirus_A121"]=0
    counts["Enterovirus_A71"]=0
    counts["Enterovirus_A76"]=0
    counts["Enterovirus_A89"]=0
    counts["Enterovirus_A90"]=0
    counts["Enterovirus_A91"]=0
    counts["Enterovirus_A92"]=0
    counts["Enterovirus_B100"]=0
    counts["Enterovirus_B101"]=0
    counts["Enterovirus_B106"]=0
    counts["Enterovirus_B107"]=0
    counts["Enterovirus_B111"]=0
    counts["Enterovirus_B112"]=0
    counts["Enterovirus_B69"]=0
    counts["Enterovirus_B73"]=0
    counts["Enterovirus_B74"]=0
    counts["Enterovirus_B75"]=0
    counts["Enterovirus_B77"]=0
    counts["Enterovirus_B78"]=0
    counts["Enterovirus_B79"]=0
    counts["Enterovirus_B80"]=0
    counts["Enterovirus_B81"]=0
    counts["Enterovirus_B82"]=0
    counts["Enterovirus_B83"]=0
    counts["Enterovirus_B84"]=0
    counts["Enterovirus_B85"]=0
    counts["Enterovirus_B86"]=0
    counts["Enterovirus_B87"]=0
    counts["Enterovirus_B88"]=0
    counts["Enterovirus_B93"]=0
    counts["Enterovirus_B97"]=0
    counts["Enterovirus_B98"]=0
    counts["Enterovirus_C102"]=0
    counts["Enterovirus_C104"]=0
    counts["Enterovirus_C105"]=0
    counts["Enterovirus_C109"]=0
    counts["Enterovirus_C113"]=0
    counts["Enterovirus_C116"]=0
    counts["Enterovirus_C117"]=0
    counts["Enterovirus_C118"]=0
    counts["Enterovirus_C96"]=0
    counts["Enterovirus_C99"]=0
    counts["Enterovirus_D68"]=0
    counts["Enterovirus_D70"]=0
    counts["Enterovirus_D94"]=0
    counts["Poliovirus_Sabinlike1"]=0
    counts["Poliovirus_Sabinlike2"]=0
    counts["Poliovirus_Sabinlike3"]=0
    counts["Poliovirus_wt1"]=0
    counts["Poliovirus_wt2"]=0
    counts["Poliovirus_wt3"]=0
    counts["*"]=0
    counts["?"]=0

    detail_dict= {
    "Enterovirus_D111": Counter(),
    "Coxsackievirus_A1": Counter(),
    "Coxsackievirus_A10": Counter(),
    "Coxsackievirus_A11": Counter(),
    "Coxsackievirus_A12": Counter(),
    "Coxsackievirus_A13": Counter(),
    "Coxsackievirus_A14": Counter(),
    "Coxsackievirus_A15": Counter(),
    "Coxsackievirus_A16": Counter(),
    "Coxsackievirus_A17": Counter(),
    "Coxsackievirus_A18": Counter(),
    "Coxsackievirus_A19": Counter(),
    "Coxsackievirus_A2": Counter(),
    "Coxsackievirus_A20": Counter(),
    "Coxsackievirus_A21": Counter(),
    "Coxsackievirus_A22": Counter(),
    "Coxsackievirus_A24": Counter(),
    "Coxsackievirus_A3": Counter(),
    "Coxsackievirus_A4": Counter(),
    "Coxsackievirus_A5": Counter(),
    "Coxsackievirus_A6": Counter(),
    "Coxsackievirus_A7": Counter(),
    "Coxsackievirus_A8": Counter(),
    "Coxsackievirus_A9": Counter(),
    "Coxsackievirus_B1": Counter(),
    "Coxsackievirus_B2": Counter(),
    "Coxsackievirus_B3": Counter(),
    "Coxsackievirus_B4": Counter(),
    "Coxsackievirus_B5": Counter(),
    "Coxsackievirus_B6": Counter(),
    "Echovirus_E1": Counter(),
    "Echovirus_E11": Counter(),
    "Echovirus_E12": Counter(),
    "Echovirus_E13": Counter(),
    "Echovirus_E14": Counter(),
    "Echovirus_E15": Counter(),
    "Echovirus_E16": Counter(),
    "Echovirus_E17": Counter(),
    "Echovirus_E18": Counter(),
    "Echovirus_E19": Counter(),
    "Echovirus_E2": Counter(),
    "Echovirus_E20": Counter(),
    "Echovirus_E21": Counter(),
    "Echovirus_E24": Counter(),
    "Echovirus_E25": Counter(),
    "Echovirus_E26": Counter(),
    "Echovirus_E27": Counter(),
    "Echovirus_E29": Counter(),
    "Echovirus_E3": Counter(),
    "Echovirus_E30": Counter(),
    "Echovirus_E31": Counter(),
    "Echovirus_E32": Counter(),
    "Echovirus_E33": Counter(),
    "Echovirus_E4": Counter(),
    "Echovirus_E5": Counter(),
    "Echovirus_E6": Counter(),
    "Echovirus_E7": Counter(),
    "Echovirus_E9": Counter(),
    "Enterovirus_A114": Counter(),
    "Enterovirus_A119": Counter(),
    "Enterovirus_A120": Counter(),
    "Enterovirus_A121": Counter(),
    "Enterovirus_A71": Counter(),
    "Enterovirus_A76": Counter(),
    "Enterovirus_A89": Counter(),
    "Enterovirus_A90": Counter(),
    "Enterovirus_A91": Counter(),
    "Enterovirus_A92": Counter(),
    "Enterovirus_B100": Counter(),
    "Enterovirus_B101": Counter(),
    "Enterovirus_B106": Counter(),
    "Enterovirus_B107": Counter(),
    "Enterovirus_B111": Counter(),
    "Enterovirus_B112": Counter(),
    "Enterovirus_B69": Counter(),
    "Enterovirus_B73": Counter(),
    "Enterovirus_B74": Counter(),
    "Enterovirus_B75": Counter(),
    "Enterovirus_B77": Counter(),
    "Enterovirus_B78": Counter(),
    "Enterovirus_B79": Counter(),
    "Enterovirus_B80": Counter(),
    "Enterovirus_B81": Counter(),
    "Enterovirus_B82": Counter(),
    "Enterovirus_B83": Counter(),
    "Enterovirus_B84": Counter(),
    "Enterovirus_B85": Counter(),
    "Enterovirus_B86": Counter(),
    "Enterovirus_B87": Counter(),
    "Enterovirus_B88": Counter(),
    "Enterovirus_B93": Counter(),
    "Enterovirus_B97": Counter(),
    "Enterovirus_B98": Counter(),
    "Enterovirus_C102": Counter(),
    "Enterovirus_C104": Counter(),
    "Enterovirus_C105": Counter(),
    "Enterovirus_C109": Counter(),
    "Enterovirus_C113": Counter(),
    "Enterovirus_C116": Counter(),
    "Enterovirus_C117": Counter(),
    "Enterovirus_C118": Counter(),
    "Enterovirus_C96": Counter(),
    "Enterovirus_C99": Counter(),
    "Enterovirus_D68": Counter(),
    "Enterovirus_D70": Counter(),
    "Enterovirus_D94": Counter(),
    "Poliovirus_Sabinlike1": Counter(),
    "Poliovirus_Sabinlike2": Counter(),
    "Poliovirus_Sabinlike3": Counter(),
    "Poliovirus_wt1": Counter(),
    "Poliovirus_wt2": Counter(),
    "Poliovirus_wt3": Counter(),
    "*": Counter(),
    "?":     Counter()
    }

    total = 0

    with open(str(report),"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total +=1
            try:
                if int(row["mapping_len"]) > 1000:
                    detail_dict[row["display_name"]][row["best_reference"]]+=1
                    counts[row["display_name"]]+=1
            except:
                if int(row["aln_block_len"]) > 1000:
                    detail_dict[row["display_name"]][row["best_reference"]]+=1
                    counts[row["display_name"]]+=1

    count_str = f"{sample},"
    for i in counts:
        count_str += f"{counts[i]},"
    count_str = count_str.rstrip(',') + '\n'
    csv_out.write(count_str)

    analysis_dict = {}
    for key in counts:
        if key not in ['*',"?"]:
            if counts[key] > args.min_reads and 100*(counts[key]/total)> args.min_pcent:
            
                top = detail_dict[key].most_common()[0]
                analysis_dict[key] = (top[0],f">{key} best_reference={top[0]} num_reads={counts[key]}")

    read_dict = defaultdict(list)

    with open(str(report),"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["display_name"] in analysis_dict:
                read_dict[row["display_name"]].append(row["read_name"])
    print(','.join(analysis_dict.keys()))

    return analysis_dict, read_dict
    

if __name__ == '__main__':

    args = parse_args()

    ref_dict = make_ref_dict(str(args.references))

    csv_report = open(str(args.out_counts), "w")

    analysis_dict,read_dict = count_and_return_analysis_dict(args.csv, csv_report, args.sample)

    for ref in analysis_dict:
        ref_file = args.output_path + "/" + ref + ".fasta"
        best_ref = analysis_dict[ref][0]
        header = analysis_dict[ref][1]
        with open(ref_file,"w") as fw:
            fw.write(f"{header}\n{ref_dict[best_ref]}\n")
        
        read_file = args.output_path + "/" + ref + ".fastq"
        with open(read_file,"w") as fw:
            records = []
            for record in SeqIO.parse(args.reads,"fastq"):
                if record.id in read_dict[ref]:
                    records.append(record)

            SeqIO.write(records, fw, "fastq")

    csv_report.close()

