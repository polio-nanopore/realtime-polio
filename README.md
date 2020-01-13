# real-time polio
This pipeline complements [``RAMPART``](https://github.com/artic-network/rampart) and continues downstream analysis to consensus level.

## Installing
Clone this repository:

```
git clone https://github.com/aineniamh/realtime-polio.git
```

Create conda environment and activate it:

```
cd realtime-polio
conda env create -f environment.yml
conda activate realtime-polio
```

To deactivate the environment, simply type:

```
conda deactivate 
```

## Running RAMPART

Create run folder:

```
mkdir [run_name]
cd [run_name]
```

Where `[run_name]` is whatever you are calling todays run (as specified in MinKNOW).

Run RAMPART:

```
rampart --protocol ../rampart-polio --basecalledPath ~/MinKNOW/data/reads/[run_name]/pass \
--annotationOptions barcode_set=[native | rapid | pcr | all]
```

`basecalledPath` should be set to whereever MinKNOW/guppy is going to write its basecalled files.

Open a web browser to view [http://localhost:3000](http://localhost:3000)

## Downstream analysis

Can be performed within the RAMPART GUI for each sample by clicking on the button to 'analyse to consensus'.


(Will change by Wednesday 13/Jan/20)
Alternatively, all samples can be analysed by running the following (giving the path to the realtime_polio directory as demonstrated):
```
snakemake --snakefile path/to/realtime_polio/rampart/pipelines/analyse_all/Snakefile --configfile snake_configuration.yaml --cores 10
```
