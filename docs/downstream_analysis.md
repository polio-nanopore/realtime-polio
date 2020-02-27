
## Downstream analysis

### Quick usage

Recommended: all samples can be analysed in parallel by editing the following command to give the path to realtime-polio and then typing it into the command line:

```
postbox -p path/to/realtime-polio
```

```
usage: postbox [-h] -p PROTOCOL [-q PIPELINE] [-d RUN_DIRECTORY]
               [-r RUN_CONFIGURATION] [-c CSV] [-t THREADS]
```

Alternatively, for each sample, the downstream analysis can be performed within the RAMPART GUI by clicking on the button to 'Analyse to consensus'.

## Pipeline description

The bioinformatic pipeline was developed using [snakemake](https://snakemake.readthedocs.io/en/stable/).

1. The server process of ``RAMPART`` watches the directory where the reads will be produced.
2. This snakemake takes each file produced in real-time and identifies the barcodes using [``porechop``](https://github.com/rambaut/Porechop).
3. Reads are mapped against a panel of references using [``minimap2``](https://github.com/lh3/minimap2). 
4. This information is collected into a csv file corresponding to each read file and the information is visualised in a web-browser, with depth of coverage and composition for each sample shown.
5. Once sufficient depth is achieved, the anaysis pipeline can be started for one sample at a time by clicking in the web browser or, to run analysis for all samples, type ``postbox -p path/to/realtime-polio`` on the command line.
6. The downstream analysis pipeline runs the following steps:
    - [``binlorry``](https://github.com/rambaut/binlorry) parses through the fastq files with barcode labels, pulling out the relevant reads and binning them into a single fastq file for each sample. It also applies a read-length filter (pre-set in the config file to only include full length amplicons).
    - The number of reads mapping to distinct viruses is assessed with a custom python script (``parse_ref_and_depth.py``) and reports whether multiple types of viruses are present in the sample and the number of corresponding reads. 
    - An iterative neural-net based polishing cycle is performed per virus type to provide a consensus sequence in ``.fasta`` format.  [``racon``](https://github.com/isovic/racon) and [``minimap2``](https://github.com/lh3/minimap2) are run iteratively four times, with gap removal in each round, against the fastq reads and then a final polishing consensus-generation step is performed using [``medaka consensus``](https://github.com/nanoporetech/medaka). 
    - The pipeline produces summary of the composition of each sample is provided and a report for each virus found, including distance to Sabin if vaccine-related.

### Output

By default the downstream analysis output will be put in a directory called ``analysis``. 

Within that directory will be:
- a ``reports`` directory with a file corresponding to each sample, detailing the virus present in the sample and how many SNPS from the closest reference in the database they are.
- a ``consensus_sequences`` directory with ``.fasta`` files for each sample. If the sample contained a mixture of viruses, all viruse sequences present at high enough levels in the sample will be in that file.
- ``sample_composition_summary.csv`` is a summary file that gives the raw read counts for each sample that have mapped to particular virus sequences. 

These are the main output files with summary information and the consensus sequences can be taken for further analysis at this point (i.e. alignments and phylogenetic trees). This directory also contains detailed output of the different steps performed in the analysis.

- ``binned_sample.csv`` and ``binned_sample.fastq`` are present for each sample. These are the output of ``BinLorry``. The csv file contains the mapping information from ``RAMPART``
- Within each ``binned_sample`` directory are many of the intermediate files produced during the analysis, including outputs of the rounds of racon polishing and medaka consensus generation. 
