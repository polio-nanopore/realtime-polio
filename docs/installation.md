
## Requirements
This pipeline will run on MacOSX or Linux. An install of Miniconda will make the setup of this pipeline on your local machine much more streamlined. To install Miniconda, visit [here](https://docs.conda.io/en/latest/miniconda.html) in a browser, scroll to see the install link specific for your operating system (MacOSX or Linux) and follow the link to the download instructions. We recommend to install the 64-bit Python 3.7 version of Miniconda. 

To clone the repository you will need to use a ``git`` command, many computers will have this available to use already, however for first time terminal users on a MacOSX system, you may need to install Xcode Command Line Developer Tools. You will know if this is the case when you try to run the installation instructions below and you see this screen pop up:

<img src="https://github.com/aineniamh/realtime-polio/blob/master/rampart/figures/xcode_popup.png" width="500">

Follow the instructions to install Xcode and then start the installation below again. 


## Installation
Clone this repository:

```
git clone https://github.com/polio-nanopore/realtime-polio
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

More information about conda environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/index.html) and a conda 'cheatsheet' can be found [here](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

## Check the install worked 

To check you have everything you need to run RAMPART, type:

```
rampart --help
```

If you see the menu shown in [RAMPART command line options](#rampart-command-line-options), you have successfully installed RAMPART. 

To check you have successfully installed all the necessary software for the downstream analysis, type:

```
postbox -h
```

If the menu shown [here](#quick-usage) appears, you know you've successfully installed all the necessary software. If you get an error message at this stage, the [wiki](https://github.com/polio-nanopore/realtime-polio/wiki/Cryptic-numpy-error) contains some common troubleshooting steps. 

Once you confirm everything has installed correctly, you're ready to set up your first run! 