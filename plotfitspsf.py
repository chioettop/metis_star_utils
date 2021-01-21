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
from astropy.table import Table
from astropy.table import hstack, vstack

import pandas as pd
import os.path

    
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
file = 'http://metisarchive.oato.inaf.it/marc//20200515_01_PFM_IT-6B1_IOM-coarseAdj/L0/solo_l0_metis-vl-image_0642857088_v01.fits'
# file UV che Andretta dice mostrare problemi in OBT
file = "http://metisarchive.oato.inaf.it/marc//20200618_01_PFM_RSICW/L0/solo_l0_metis-vl-image_0645742792_v01.fits"
"""
file = "G:\Il mio Drive\METIS\Stellar fields identification\PFM_RSICW\solo_l0_metis-vl-image_0646069288_v01.fits"

# Image from Metis UV instrument?
UV = "uv" in file
ref_frame = 'detector'
dark = False

hdul = fits.open(file)

# load raw image
if ref_frame=='sky':
    image_data = np.rot90(hdul[0].data, 3)
if ref_frame=='detector':
    image_data = hdul[0].data    
    
if dark:
    t_exp = sl.texp_from_fits_hdr(hdul)
    image_data = sl.remove_dark(image_data, t_exp, UV)

META_KERNEL = 'solo_ANC_soc-flown-mk.tm'
KERNEL_PATH = '../kernels/solar-orbiter/kernels/mk'

spice = MetisSpice(META_KERNEL, KERNEL_PATH)

obt = sl.obt_from_fits_hdr(hdul)

et = spice.scs2et(obt)

print("UTC time: ", spice.et2utc(et))

ra, dec, roll = spice.boresight(et, UV=UV)

# plot fits image
fig, ax = plt.subplots(figsize=(10,8))
plt.subplots_adjust(bottom=0.2)

sl.plot_fits(file, ax=ax, coor=(ra, dec), utc=spice.et2utc(et), ref_frame=ref_frame, dark=dark)

# build wcs from nominal Metis pointing
wcs = sl.wcs_from_boresight(ra, dec, roll, UV)

# adjust boresight with nominal wcs
# need to check if proper motion should be included in the
# calculation: https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_cali/ssec_cu3ast_cali_frame.html
ra_adj, dec_adj = wcs.wcs_pix2world(np.array([wcs.pixel_shape])/2, 0).flatten()

# find stars through catalog search
catalog_stars = sl.simbad_search(ra, dec, max_mag=7) 

# add sensor coordinates to catalog star list
x, y = wcs.wcs_world2pix(catalog_stars['ra'], catalog_stars['dec'], 0)
catalog_stars['xsensor'] = x
catalog_stars['ysensor'] = y

# remove stars that fall outside of the sensor
catalog_stars.remove_rows((x < 0) | (y < 0) | 
                          (x > wcs.pixel_shape[0]) | (y > wcs.pixel_shape[0]))

"""    
# file with WCS calculated from star positions
hdulist = fits.open('WCS for solo_l0_metis-vl-image_0642850455_v01.fits')
# Parse the WCS keywords in the primary HDU
wcs = astropy.wcs.WCS(hdulist[0].header) 
"""

sl.simbad_plot(catalog_stars, wcs, ax=ax)

# plot FOV circles and center of FOV
scale = -wcs.wcs.cdelt[0]
fov = (sl.METIS_fov_min/scale, sl.METIS_fov/scale)
center = wcs.wcs.crpix
sl.plot_fov_circle(center, fov[0], 'r', ax)
sl.plot_fov_circle(center, fov[1], 'r', ax)
ax.scatter(center[0], center[1], s=40, c='r', marker='+')

# find peaks close to catalog stars
vis_stars = sl.find_nearby_stars(image_data, 
                                 catalog_stars['xsensor'], 
                                 catalog_stars['ysensor'], 
                                 size=5)

catalog_stars = hstack([catalog_stars, vis_stars])  

# subset of stars with matched SPF
matched_stars = catalog_stars[~np.isnan(catalog_stars['peak_value'])]  
                  
ax.plot(vis_stars['x_peak'], vis_stars['y_peak'], 'b*')

sl.plot_fitted_psf(catalog_stars, image_data, size=5) 
   
# distance of each catalog star with the closest SPF on sensor
# closest_SPF is the index to the vis_stars table, 0 based (not the 'id', which starts from 1)
#sl.comp_star_distance(catalog_stars, vis_stars)

# display buttons below plot and handle click events
def clbk_fit(event):
    sl.plot_gaussian_fit(image_data, *ax.viewLim.bounds)
    plt.show()
    
def clbk_hdr(event):
    print(hdul[1].header)
    
def clbk_dist(event):
    sl.plot_star_distance(
        (matched_stars['xsensor'], matched_stars['ysensor']),
        (matched_stars['x_peak'], matched_stars['y_peak']), 
        wcs.pixel_shape[0]
        )
    
def clbk_fit_wcs(event):
    from astropy.coordinates import SkyCoord
    
    # build list of sky coordinates
    world_coords = SkyCoord(matched_stars['ra'], matched_stars['dec'], 
                            frame="icrs", unit="deg")

    xy = (matched_stars['x_peak'], matched_stars['y_peak'])

    proj_point = SkyCoord(ra, dec, frame="icrs", unit="deg")

    
    fit_wcs, fit_res = sl.wcs_fit(
        xy=xy,
        world_coords=world_coords,
        proj_point=proj_point,
        projection='TAN',
        UV=False
    )
    
    # update catalog stars positions on detector
    x, y = fit_wcs.wcs_world2pix(catalog_stars['ra'], catalog_stars['dec'], 0)
    catalog_stars['xsensor'] = x
    catalog_stars['ysensor'] = y
    x, y = fit_wcs.wcs_world2pix(matched_stars['ra'], matched_stars['dec'], 0)
    matched_stars['xsensor'] = x
    matched_stars['ysensor'] = y

    
    fig = plt.figure()
    ax = fig.gca()
    sl.plot_image(image_data, ax=ax)
    sl.simbad_plot(catalog_stars, fit_wcs, ax=ax)
    
    print("Fitted WCS")
    print("\tOriginal pixel coords\n", xy)
    print("\tComputed\n", fit_wcs.world_to_pixel(world_coords))
    world_coords_new=fit_wcs.pixel_to_world(*xy)
    print("\tDifference\n", world_coords.separation(world_coords_new))
    # roll angle computed reversing eq. in wcs_from_boresight()
    print("Roll angle: ", fit_res.x[0] -0.077 + 90)
    print(fit_wcs)
    print()
    
    print("\tOriginal WCS")
    print("\tComputed pixel coords\n", wcs.world_to_pixel(world_coords))
    world_coords_new=wcs.pixel_to_world(*xy)
    print("\tDifference\n", world_coords.separation(world_coords_new))
    
def clbk_daofind(event):
    global vis_stars
    found = sl.find_stars_sub(image_data, *ax.viewLim.bounds, len(vis_stars), 
                              fwhm=2, ax=ax, plot=True)
    if found:
        vis_stars = vstack([vis_stars, found])
        # recomputes distances and closest star
        sl.comp_star_distance(catalog_stars, vis_stars)

def handle_close(evt):
    hdul.close()

fig.canvas.mpl_connect('close_event', handle_close)    # when the plot window closes, close the FITS file

buttons = {
    'Fit PSF': clbk_fit,
    'FITS hdr': clbk_hdr,
    'SPF dist': clbk_dist,
    'Fit WCS':  clbk_fit_wcs,
    'DAOFind':  clbk_daofind
    }

btn = []
for n, b in enumerate(buttons):
    bax = plt.axes([n/5, 0.05, 0.1, 0.075])
    btn.append(Button(bax, b))
    btn[n].on_clicked(buttons[b])
    
    
plt.show()