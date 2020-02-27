
## Setting up your run

You can familiarise yourself with some basic bash commands with the document found in [resources](https://github.com/polio-nanopore/realtime-polio/blob/master/resources/basic_command_line.md).


### Step 1:

Create run folder:

```
mkdir [run_name]
cd [run_name]
```

Where `[run_name]` is whatever you are calling todays run (as specified in MinKNOW).

### Step 2:

Create a ``run_configuration.json`` and a ``barcodes.csv`` file (examples are shown in [resources](https://github.com/polio-nanopore/realtime-polio/tree/master/resources/examples)).

If you have a ``run_configuration.json`` file and a ``barcodes.csv`` file, you can run RAMPART with very few command line options. A template of the configuration files needed to run both RAMPART and the downstream analysis pipeline is provided in the examples directory.

The run_configuration.json can specify the path to your basecalled reads or alternatively you can input that information on the command line. `basecalledPath` should be set to wherever MinKNOW/guppy is going to write its basecalled files. If you want alter where the annotations files from RAMPART or the analysis files from the downstream pipeline are put, you can add the optional ``"annotatedPath"`` and ``"outputPath"`` options. By default the annotations are written to a directory called ``annotations`` and the analysis output is written to a directory called ``analysis``.

```
run_configuration.json

{
  "title": "MinION_run_example",
  "basecalledPath": "fastq_pass"
}
```

Optional for RAMPART, but required for the downstream analysis pipeline, the ``barcodes.csv`` file describes which barcode corresponds to which sample. 

```
barcodes.csv

sample,barcode
sample1,BC01
sample2,BC02
sample3,BC03
sample4,BC04
```

> *Note*: the barcode name needs to be accurate and start BC. 

> *Note*: the sample names should be unique, or they will try to merge and may cause the pipeline to exit. 

## Checklist

- The conda environment ``realtime-polio`` is active.
- ``barcodes.csv`` file with sample to barcode mapping either in the current directory or the path to it will need to be provided.
- ``annotations`` directory with csv files from RAMPART
- The path to basecalled ``.fastq`` files is provided either in the ``run_configuration.json`` or it will need to be specified on the command line.

Note that these file names and locations are the default setup. If you have files in other locations or with different names, this is fine, but you'll need to give RAMPART that information using the options detailed in [RAMPART command line options](#rampart-command-line-options).
