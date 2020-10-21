# metis_star_utils
Utilities to plot star maps visible by Metis, calculate PSF, find stars from Metis images.

## Installation
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have Python installed.

1. Donwload the source code of the script from this site: use the green *Code* button or clone the repository with `git clone...`

1. Open a conda console and go to the directory where you downloaded the source code.

1. Optional, if you want to avoid conflicts with your other Python projects: after downloading the files, create a virtual environment, eg:  
`conda create --name metis_star_utils pip`

Activate the environment:  
`conda activate metis_star_utils`

1. Download library dependencies:  
`pip install -r requirements.txt`

1. Download Solar Orbiter SPICE Kernel dataset here: [Solar-Orbiter SPICE Kernel](https://repos.cosmos.esa.int/socci/rest/api/latest/projects/SPICE_KERNELS/repos/solar-orbiter/archive?format=zip) and unzip it somewhere. Make note of the directory, and locate the file `kernels/solar-orbiter/kernels/mk/solo_ANC_soc-flown-mk.tm`. You will need to provide the location to the script.

More information on SO SPICE Kernel [here](https://www.cosmos.esa.int/web/spice/solar-orbiter)

1. The installation is now ready. Read on for usage instructions.

## starfield.py
Plots a star field at specific UTC time.  
If you used a virtual environment to install the program, remember to activate it before launching the script.
The script needs to know the location of the SPICE metakernel to load. Make sure the SPICE kernels used contain the specific date and are up-to-date.  

To plot using flown data, use:  `solo_ANC_soc-flown-mk.tm`  

To plot field with dates in the future, use predicted kernels: `solo_ANC_soc-pred-mk.mk`

Example:  
`python starfield.py "2020 MAY 15 07:03:33.7649" -k ../kernels/solar-orbiter/kernels/mk/solo_ANC_soc-flown-mk.tm`

Type `python starfield.py -h` for help.

### Notes:
1. UV does not use latest calibration data yet
