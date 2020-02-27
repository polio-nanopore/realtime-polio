# realtime-polio
This pipeline complements [``RAMPART``](https://github.com/artic-network/rampart) and continues downstream analysis to consensus level.

<img src="https://github.com/aineniamh/realtime-polio/blob/master/rampart/figures/rampart_polio.png">

## Documentation

  * [Requirements](#requirements)
  * [Installation](docs/installation.md)
  * [Check the install worked](docs/installation.md)
  * [Setting up your run](docs/setup.md)
  * [Checklist](docs/setup.md)
  * [Running RAMPART](docs/running_rampart.md)
  * [RAMPART command line options](docs/running_rampart.md)
  * [Stopping RAMPART](docs/running_rampart.md)
  * [Downstream analysis](docs/downstream_analysis.md)
     * [Quick usage](docs/downstream_analysis.md)
     * [Pipeline description](docs/downstream_analysis.md)
     * [Output](docs/downstream_analysis.md)
  * [Reference FASTA](#reference-fasta)
  * [Troubleshooting](#troubleshooting)
  * [Updating your repository](#updating-your-repository)
  * [License](#license)


## Requirements
This pipeline will run on MacOSX or Linux. An install of Miniconda will make the setup of this pipeline on your local machine much more streamlined. To install Miniconda, visit [here](https://docs.conda.io/en/latest/miniconda.html) in a browser, scroll to see the install link specific for your operating system (MacOSX or Linux) and follow the link to the download instructions. We recommend to install the 64-bit Python 3.7 version of Miniconda. 

## Reference FASTA

The ``references.fasta`` file is a detailed database containing VP1 sequences from an array of poliovirus sequences, including a representative sequence for each of the Sabin vaccine strains, and also other enterovirus VP1 sequences.

## Troubleshooting

If you're having issues running the pipeline, and you're certain the checklist has been tickeds, check out some common [issues on the wiki](https://github.com/aineniamh/realtime-polio/wiki)

## Updating your repository

If it's been a while since you cloned the github repository to your local machine and you want to check if there has been any changes, you can simply update to the latest version of the pipeline by typing in a terminal in the realtime-polio directory:

```
git pull
```

This will sync any changes that have been made to your local computer.

If you see any changes in the ``environment.yml`` file, make sure you're in the realtime-polio environment and then, you may need to type:

```
conda env update --file environment.yml 
```

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)