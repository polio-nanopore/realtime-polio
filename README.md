# real-time polio
This pipeline complements [``RAMPART``](https://github.com/artic-network/rampart) and continues downstream analysis to consensus level.

<img src="https://github.com/aineniamh/realtime_polio/blob/master/rampart/figures/rampart_polio.png">



## Table of contents


  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Setting up your run](#setting-up-your-run)
  * [Running RAMPART](#running-rampart)
  * [RAMPART command line options](#rampart-command-line-options)
  * [Downstream analysis](#downstream-analysis)
     * [Quick usage](#quick-usage)
     * [Pipeline description](#pipeline-description)
     * [Reference FASTA](#reference-fasta)
     * [Output](#output)
  * [Troubleshooting](#troubleshooting)
  * [License](#license)

## Requirements
This pipeline will run on MacOS and Linux. An install of Miniconda will make the setup of this pipeline on your local machine much more streamlined. To install Miniconda, visit here https://conda.io/docs/user-guide/install/ in a browser, select your type of machine (mac or linux) and follow the link to the download instructions. We recommend to install the 64-bit Python 3.6 version of Miniconda. 

## Installation
Clone this repository:

```
git clone https://github.com/aineniamh/realtime-polio.git
```

1. Create the conda environment.
This may take some time, but will only need to be done once. It allows the pipeline to access all the software it needs, including RAMPART.

```
cd realtime-polio
conda env create -f environment.yml
```

2. Activate the conda environment.

```
conda activate realtime-polio
```

## Setting up your run


## Running RAMPART

Create run folder:

```
mkdir [run_name]
cd [run_name]
```

Where `[run_name]` is whatever you are calling todays run (as specified in MinKNOW).

If you have a run_configuration.json file and a barcodes.csv file, you can run RAMPART with very few command line options.

The run_configuration.json can specify the path to your basecalled reads or alternatively you can input that information on the command line. `basecalledPath` should be set to wherever MinKNOW/guppy is going to write its basecalled files.
You may also supply a barcodes.csv file that gives information about which sample corresponds to which barcode.

With this setup, to run RAMPART:

```
rampart --protocol ~/Documents/realtime-polio/rampart/ --referencesLabel display_name
```

Open a web browser to view [http://localhost:3000](http://localhost:3000)

More information about RAMPART can be found [here](https://github.com/artic-network/rampart).

## RAMPART command line options

```
usage: rampart [-h] [-v] [--verbose] [--ports PORTS PORTS]
               [--protocol PROTOCOL] [--title TITLE]
               [--basecalledPath BASECALLEDPATH]
               [--annotatedPath ANNOTATEDPATH]
               [--referencesPath REFERENCESPATH]
               [--referencesLabel REFERENCESLABEL]
               [--barcodeNames BARCODENAMES [BARCODENAMES ...]]
               [--annotationOptions ANNOTATIONOPTIONS [ANNOTATIONOPTIONS ...]]
               [--clearAnnotated] [--simulateRealTime SIMULATEREALTIME]
               [--devClient] [--mockFailures]
```

## Downstream analysis


Can be performed within the RAMPART GUI for each sample by clicking on the button to 'Analyse to consensus'. The bioinformatic pipeline was developed using [snakemake](https://snakemake.readthedocs.io/en/stable/). 


Alternatively, all samples can be analysed in parallel by running the following:

```
postbox -p path/to/realtime_polio -t 3
```

-t specifices the number of parallel processes to run and will speed up your analysis if running multiple samples. 

## Pipeline description

1. The server process of ``RAMPART`` watches the directory where the reads will be produced.
2. This snakemake takes each file produced in real-time and identifies the barcodes using [``porechop``](https://github.com/rambaut/Porechop).
3. Reads are mapped against a panel of references using [``minimap2``](https://github.com/lh3/minimap2). 
4. This information is collected into a csv file corresponding to each read file and the information is visualised in a web-browser, with depth of coverage and composition for each sample shown.
5. Once sufficient depth is achieved, the anaysis pipeline can be started for one sample at a time by clicking in the web browser or, to run analysis for all samples, type ``artic-rampart --pipeline analyse_all`` on the command line.
6. The downstream analysis pipeline runs the following steps:
    - [``binlorry``](https://github.com/rambaut/binlorry) parses through the fastq files with barcode labels, pulling out the relevant reads and binning them into a single fastq file for each sample. It also applies a read-length filter (pre-set in the config file to only include full length amplicons).
    - The number of reads mapping to distinct viruses is assessed with a custom python script (``parse_ref_and_depth.py``) and reports whether multiple types of viruses are present in the sample and the number of corresponding reads. 
    - An iterative neural-net based polishing cycle is performed per virus type to provide a consensus sequence in ``.fasta`` format.  [``racon``](https://github.com/isovic/racon) and [``minimap2``](https://github.com/lh3/minimap2) are run iteratively four times, with gap removal in each round, against the fastq reads and then a final polishing consensus-generation step is performed using [``medaka consensus``](https://github.com/nanoporetech/medaka). 
    - The pipeline produces summary of the composition of each sample is provided and a report for each virus found, including distance to Sabin if vaccine-related.

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)