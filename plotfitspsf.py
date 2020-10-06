# -*- coding: utf-8 -*-
"""
Created on 

@author: P. Chioetto 
"""

import sys
import numpy as np
from spicelib import MetisSpice
import starlib as sl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from astropy.io import fits
    
if len(sys.argv) < 2:
    raise SystemExit('Please provide the name of the METIS FITS file to use.')

file = sys.argv[1]
"""
# Examples of files to use
file = "http://metisarchive.oato.inaf.it/marc//20200515_01_PFM_IT-6B1-IOM-Adjustment/OUTPUT/solo_l0_metis-vl-image_0642850450_v01.fits"
file = '../PFM_RSICW\\solo_l0_metis-vl-image_0645545693_v01.fits'
file = "G:\Il mio Drive\METIS\Stellar fields identification\PFM_RSICW\solo_l0_metis-vl-image_0645742792_v01.fits"
"""
# Image from Metis UV instrument?
UV = "uv" in file

hdul = fits.open(file)

META_KERNEL = 'solo_ANC_soc-flown-mk.tm'
KERNEL_PATH = '../kernels/solar-orbiter/kernels/mk'

spice = MetisSpice(META_KERNEL, KERNEL_PATH)

obt = sl.obt_from_fits_hdr(hdul)
et = spice.scs2et(obt)

print("UTC time: ", spice.et2utc(et))

ra, dec, roll = spice.boresight(et, UV=UV)

wcs = sl.wcs_from_boresight(ra, dec, roll, UV)

catalog_stars = sl.simbad_search(ra, dec, 8) # find stars through catalog search

image_data = np.rot90(hdul[0].data, 3)
#image_data = remove_dark(image_data, hdul[1].header['DIT']/1000)    

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

sl.plot_fits(file, ax=ax, coor=(ra, dec), rotate=True, dark=False)

"""    
# file with WCS calculated from star positions
hdulist = fits.open('WCS for solo_l0_metis-vl-image_0642850455_v01.fits')
# Parse the WCS keywords in the primary HDU
wcs = astropy.wcs.WCS(hdulist[0].header) 
"""

sl.simbad_plot(catalog_stars, wcs, ax=ax)
scale = -wcs.wcs.cdelt[0]
sl.plot_fov_circle(wcs.wcs.crpix, sl.METIS_fov/scale, 'r', ax)
sl.plot_fov_circle(wcs.wcs.crpix, sl.METIS_fov_min/scale, 'r', ax)

vis_stars = sl.find_stars(image_data, fwhm=3.5, ax=ax, plot=True)

def clbk_fit(event):
    sl.plot_gaussian_fit(image_data, *ax.viewLim.bounds)
    
def clbk_hdr(event):
    print(hdul[1].header)

def handle_close(evt):
    hdul.close()

fig.canvas.mpl_connect('close_event', handle_close)    # when the plot window closes, close the FITS file

axfit = plt.axes([0.7, 0.05, 0.1, 0.075])
axhdr = plt.axes([0.81, 0.05, 0.1, 0.075])
fitbtn = Button(axfit, 'Fit')
fitbtn.on_clicked(clbk_fit)
hdrbtn = Button(axhdr, 'FITS hdr')
hdrbtn.on_clicked(clbk_hdr)

plt.show()
#hdul.close()
