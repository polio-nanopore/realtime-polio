import argparse
from Bio import SeqIO
import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import csv

parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

parser.add_argument("--csv", action="store", type=str, dest="report")
parser.add_argument("--reads", action="store", type=str, dest="reads")
parser.add_argument("--references", action="store", type=str, dest="references")

parser.add_argument("--sample", action="store", type=str, dest="sample")
parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

parser.add_argument("--output_path", action="store", type=str, dest="output_path")


args = parser.parse_args()

def make_ref_dict(references):
    refs = {}
    for record in SeqIO.parse(references,"fasta"):
        refs[record.id]=record.seq
    return refs

ref_dict = make_ref_dict(str(args.references))

detail_dict= {}
counts = {
    "Sabin1":0,
    "Sabin2":0,
    "Sabin3":0,
    "Poliovirus1-wt":0,
    "Poliovirus2-wt":0,
    "Poliovirus3-wt":0,
    "NonPolioEV":0,
    "*":0,
    "?":0
}
total = 0
with open(str(args.report),"w") as f:
    reader = csv.DictReader(f)
    for row in reader:
        total +=1
        detail_dict[row["best_reference"]] = row["display_name"]
        counts[row["display_name"]]+=1

to_analyse = []
for key in counts:
    if counts[key] > args.min_reads:
        if 100*(counts[key]/total)> args.min_pcent:
            if i not in ['*',"?","NonPolioEV"]:
                to_analyse.append(key)


fw = open(args.output_path + "/reference_count_{}.csv".format(args.sample),"w")
header = "sample,"
count_str = f"{args.sample},"
for i in sorted(counts):
    header += f"{i},"
    count_str += f"{counts[i]},"
header = header.rstrip(',') + '\n'
count_str = count_str.rstrip(',') + '\n'
fw.write(header)
fw.write(count_str)
fw.close()

report = pd.read_csv(args.report)
report["ref_stem"]=report["best_reference"].split("_")[0]
detailed_ref_count = report['best_reference'].value_counts()
ref_count  = report['ref_stem'].value_counts()

# if len(ref_count) > 1:
#     fig, ax = plt.subplots(figsize=(15,11))
#     sns.barplot(ref_count.index, ref_count.values, alpha=0.8)
#     plt.title('Reference profile of sample {}'.format(args.sample))
#     plt.ylabel('Read Count', fontsize=12)
#     plt.xlabel('Reference', fontsize=12)
#     plt.xticks(rotation=20)
#     fig.savefig(args.output_path + "/reference_count_{}.pdf".format(args.sample))

detail_dict= {}
for i,x in zip(list(detailed_ref_count.index), list(detailed_ref_count.values)):
    stem = i.split("_")[0]
    if stem not in detail_dict:
        detail_dict[stem] = i

total = len(report)
refs = []
for i,x in zip(list(ref_count.index), list(ref_count.values)):
    pcent = 100*(x/total)
    if x>args.min_reads and pcent > args.min_pcent:
        if i.startswith("Sabin") or i.startswith("Poliovirus"):

            refs.append(i.split('_')[0])

print(",".join(refs))

for ref in refs:
    with open(args.output_path + "/" + ref+ ".fasta","w") as fw:
        
        fw.write(">{} detail={}\n{}\n".format(ref, detail_dict[ref], ref_dict[detail_dict[ref]]))

    filtered_df = report.loc[(report["ref_stem"]==ref)]
    read_names = list(filtered_df["read_name"].values)
    new_file = args.output_path + "/" + ref + ".fastq"
    with open(new_file,"w") as fw:
        records = []
        for record in SeqIO.parse(args.reads,"fastq"):
            if record.id in read_names:
                records.append(record)

        SeqIO.write(records, fw, "fastq")

none_df = report.loc[(report["best_reference"]=='*')]
read_names = list(none_df["read_name"].values)
new_file = args.output_path + "/no_hit.fastq"
with open(new_file,"w") as fw:
    records = []
    for record in SeqIO.parse(args.reads,"fastq"):
        if record.id in read_names:
            records.append(record)

    SeqIO.write(records, fw, "fastq")