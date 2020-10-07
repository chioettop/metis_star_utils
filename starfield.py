# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:59:07 2020

@author: P. Chioetto 
"""
import matplotlib.pyplot as plt
import argparse
from os import path
from spicelib import MetisSpice
import starlib as sl

# default Spice meta kernel to load
META_KERNEL = '../kernels/solar-orbiter/kernels/mk/solo_ANC_soc-flown-mk.tm'

# command line parser
parser = argparse.ArgumentParser(description='Metis VL sensor view of stars from Simbad catalog.\n \
                                 Donwload Metis SPICE kernels here: ftp://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels')
parser.add_argument('time', metavar='time', 
                    help='UTC time of observation, eg. "2020 MAY 15 07:03:33.7649".')
parser.add_argument('-k', metavar='META_KERNEL', default=META_KERNEL,
                    help='Full path of Metis SPICE meta kernel to load. Use solo_ANC_soc-flown-mk.tm or solo_ANC_soc-pred-mk')
parser.add_argument('-m', metavar='magnitude', type=float, default=6,
                    help='Magnitude limit for star catalog search. Default is 6.')
parser.add_argument('-s', '--sky', action="store_true",
                    help='Set view orientation as Sky.')
parser.add_argument('--UV', action="store_true",
                    help='Metis UV sensor.')
args = parser.parse_args()

# load SPICE Kernels
kernel_path, meta_kernel = path.split(args.k)
spice = MetisSpice(meta_kernel, kernel_path)

# et timestamp from UTC time
et = spice.utc2et(args.time)

ra, dec, roll = spice.boresight(et, UV=args.UV)

# build Worls Coordinate System object
wcs = sl.wcs_from_boresight(ra, dec, roll, UV=args.UV)

# query SIMBAD catalog
stars = sl.simbad_search(ra, dec, args.m)

ax = sl.simbad_plot(stars, wcs, color='black', scatter=True)

# x axis scale is negative to invert axis (maybe should use inver_axis?)
scale = -wcs.wcs.cdelt[0]
sl.plot_fov_circle(wcs.wcs.crpix, sl.METIS_fov/scale, 'r', ax)
sl.plot_fov_circle(wcs.wcs.crpix, sl.METIS_fov_min/scale, 'r', ax)

fig = ax.figure
fig.suptitle("UTC Time " + spice.et2utc(et))
ax.set_title(f"RA {ra}, dec {dec}")
plt.show()

# add sensor coordinates to star list and print
x, y = wcs.wcs_world2pix(stars['ra'], stars['dec'], 0)
stars['xsensor'] = x
stars['ysensor'] = y
stars.pprint_all()
