# real-time polio
This pipeline complements [``RAMPART``](https://github.com/artic-network/rampart) and continues downstream analysis to consensus level.

## Installing
Clone this repository:

```
git clone https://github.com/aineniamh/rampart-polio.git
```

Create conda environment and activate it:

```
cd rampart-polio
conda env create -f environment.yml
conda activate rampart-polio
```

## Running

Create run folder:

```
mkdir [run_name]
cd [run_name]
```

Where `[run_name]` is whatever you are calling todays run (as specified in MinKNOW).

Run RAMPART:

```
rampart --protocol ../rampart-polio --basecalledPath ~/MinKNOW/data/reads/[run_name]/pass --annotationOptions barcode_set=[native | rapid | pcr | all]
```

`basecalledPath` should be set to whereever MinKNOW/guppy is going to write its basecalled files.

Open a web browser to view [http://localhost:3000](http://localhost:3000)

## Downstream analysis

Can be performed within the RAMPART GUI by clicking on the button to 'analyse to consensus' or by running the analysis in parallel for all samples from the command line. 
