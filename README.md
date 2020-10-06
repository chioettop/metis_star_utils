# metis_star_utils
Utils to plot star maps visible by Metis, calculate PSF, find stars from Metis images.

## Installation
Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have Python installed. Open a Miniconda console.

After downloading the files, create a virtual environment, eg:
`conda create --name metis_star_utils`

Activate the environment:
`conda activate metis_star_utils`

Download dependencies:
`pip install -r requirements.txt`

Download Solar Orbiter SPICE Kernel dataset here: [Solar-Orbiter SPICE Kernel](https://repos.cosmos.esa.int/socci/rest/api/latest/projects/SPICE_KERNELS/repos/solar-orbiter/archive?format=zip)

More information on SO SPICE Kernel [here](https://www.cosmos.esa.int/web/spice/solar-orbiter)

## starfield.py
Plots a star field at specific UTC time. Example:  
`python starfield.py "2020 MAY 15 07:03:33.7649" -k ../kernels/solar-orbiter/kernels/mk/solo_ANC_soc-flown-mk.tm`

Type `python starfield.py -h` for help.
