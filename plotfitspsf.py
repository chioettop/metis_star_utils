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
    
'''
if len(sys.argv) < 2:
    raise SystemExit('Please provide the name of the METIS FITS file to use.')

file = sys.argv[1]
'''
"""
# Examples of files to use
file = "http://metisarchive.oato.inaf.it/marc//20200515_01_PFM_IT-6B1-IOM-Adjustment/OUTPUT/solo_l0_metis-vl-image_0642850450_v01.fits"
file = '../PFM_RSICW\\solo_l0_metis-vl-image_0645545693_v01.fits'
file = "G:\Il mio Drive\METIS\Stellar fields identification\PFM_RSICW\solo_l0_metis-vl-image_0645742792_v01.fits"
file = "http://metisarchive.oato.inaf.it/marc//20200515_01_PFM_IT-6B1-IOM-Adjustment/OUTPUT/solo_l0_metis-vl-image_0642850450_v01.fits"
# file UV che Andretta dice mostrare problemi in OBT
file = "http://metisarchive.oato.inaf.it/marc//20200618_01_PFM_RSICW/L0/solo_l0_metis-uv-image_0645937059_v01.fits"
"""
file = 'G:\Il mio Drive\METIS\Stellar fields identification\PFM_RSICW\solo_l0_metis-uv-image_0645545692_v01.fits'

# Image from Metis UV instrument?
UV = "uv" in file
ref_frame = 'detector'
dark = True

hdul = fits.open(file)

META_KERNEL = 'solo_ANC_soc-flown-mk.tm'
KERNEL_PATH = '../kernels/solar-orbiter/kernels/mk'

spice = MetisSpice(META_KERNEL, KERNEL_PATH)

obt = sl.obt_from_fits_hdr(hdul)
et = spice.scs2et(obt)

print("UTC time: ", spice.et2utc(et))

ra, dec, roll = spice.boresight(et, UV=UV)

wcs = sl.wcs_from_boresight(ra, dec, roll, UV)

catalog_stars = sl.simbad_search(ra, dec, max_mag=8) # find stars through catalog search
# add sensor coordinates to catalog star list
x, y = wcs.wcs_world2pix(catalog_stars['ra'], catalog_stars['dec'], 0)
catalog_stars['xsensor'] = x
catalog_stars['ysensor'] = y

if ref_frame=='sky':
    image_data = np.rot90(hdul[0].data, 3)
if ref_frame=='detector':
    image_data = hdul[0].data    

fig, ax = plt.subplots(figsize=(10,8))
plt.subplots_adjust(bottom=0.2)

sl.plot_fits(file, ax=ax, coor=(ra, dec), utc=spice.et2utc(et), ref_frame=ref_frame, dark=dark)

"""    
# file with WCS calculated from star positions
hdulist = fits.open('WCS for solo_l0_metis-vl-image_0642850455_v01.fits')
# Parse the WCS keywords in the primary HDU
wcs = astropy.wcs.WCS(hdulist[0].header) 
"""

sl.simbad_plot(catalog_stars, wcs, ax=ax)
scale = -wcs.wcs.cdelt[0]
fov = (sl.METIS_fov_min/scale, sl.METIS_fov/scale)
center = wcs.wcs.crpix
sl.plot_fov_circle(center, fov[0], 'r', ax)
sl.plot_fov_circle(center, fov[1], 'r', ax)
ax.scatter(center[0], center[1], s=40, c='r', marker='+')

vis_stars = sl.find_stars(image_data, fwhm=2, threshold=130, ax=ax, plot=True, fov=fov)

# distance of each catalog star with the closest SPF on sensor
catalog_stars['closest_SPF'] = [
    np.argmin([np.hypot(cs['xsensor']-xv, cs['ysensor']-yv) for xv, yv in vis_stars['xcentroid','ycentroid']])
    for cs in catalog_stars]
i = catalog_stars['closest_SPF']
catalog_stars['dist_to_SPF'] = np.hypot( 
    catalog_stars['xsensor'] - vis_stars[i]['xcentroid'], 
    catalog_stars['ysensor'] - vis_stars[i]['ycentroid'])

# display buttons below plot and handle click events
def clbk_fit(event):
    sl.plot_gaussian_fit(image_data, *ax.viewLim.bounds)
    plt.show()
    
def clbk_hdr(event):
    print(hdul[1].header)
    
def clbk_dist(event):
    fig = plt.figure()
    plt.suptitle('Star distance from PSF in pixels')
    # stars with more than 50 pixels of distance from SPF are in fact not found
    cat = catalog_stars[catalog_stars['dist_to_SPF']<50]
    plt.title(f"Average = {np.round(np.mean(cat['dist_to_SPF']), 1)}")
    plt.axis('square')
    size = wcs.pixel_shape[0]
    xlim = ylim = (0, size)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.plot(cat['xsensor'], cat['ysensor'], 'o')
    for c in cat:
        plt.annotate(np.round(c['dist_to_SPF'], 1), (c['xsensor']+4, c['ysensor']+4))
    plt.show()

def handle_close(evt):
    hdul.close()

fig.canvas.mpl_connect('close_event', handle_close)    # when the plot window closes, close the FITS file

buttons = {
    'Fit PSF': clbk_fit,
    'FITS hdr': clbk_hdr,
    'SPF dist': clbk_dist}

btn = []
for n, b in enumerate(buttons):
    bax = plt.axes([0.2+n/5, 0.05, 0.1, 0.075])
    btn.append(Button(bax, b))
    btn[n].on_clicked(buttons[b])
    
    
plt.show()