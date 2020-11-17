# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:59:07 2020

@author: P. Chioetto 
"""
import starlib as sl
from spicelib import MetisSpice
from astropy.io import fits
from astropy.table import Table
from astropy.table import hstack, vstack

import os.path
import glob

import numpy as np
import pandas as pd

import warnings

import plotly.express as px
import plotly.io as pio


META_KERNEL = 'solo_ANC_soc-flown-mk.tm'
KERNEL_PATH = '../kernels/solar-orbiter/kernels/mk'

def file_range(fits_range, base_dir):
    # find files within a range
    start, stop = fits_range
    
    if start > stop:
        start, stop = stop, start
        
    files = glob.iglob(base_dir+'*.fits')
    
    found = []
    for file in files:
        beg = file.find('image_') + len('image_')
        end = file.find('_v', beg)
        seq = int(file[beg:end])
        if (seq <= stop) and (seq >= start):
            found.append(file)
            
    return found

base_path = "../PFM_RSICW/"
files = glob.glob(base_path + "*image_*.fits")
#fits_range = (645545693, 645755992) # alfa Leo
#fits_range = (645742792, 645755992) # rho Leo
#files = file_range(fits_range, base_dir="../PFM_RSICW/") # range of vl of 22Mb with a Leo visible

assert files, "No files found."

# initialize Dataframes
# that will contain info on stars and images
stars_df = pd.DataFrame()
img_info_df = pd.DataFrame()

spice = MetisSpice(META_KERNEL, KERNEL_PATH)

for n, file in enumerate(files):
    
    print("Processing ", n, "/", len(files))
    
    UV = 'uv' in file
    hdul = fits.open(file)
    image_data = hdul[0].data

    obt = sl.obt_from_fits_hdr(hdul)
    et = spice.scs2et(obt)

    #print("UTC time: ", spice.et2utc(et))

    ra, dec, roll = spice.boresight(et, UV=UV)

    # build wcs from nominal Metis pointing
    wcs = sl.wcs_from_boresight(ra, dec, roll, UV)

    # adjust boresight with nominal wcs
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

    # find peaks close to catalog stars
    vis_stars = sl.find_nearby_stars(image_data, 
                                 catalog_stars['xsensor'], 
                                 catalog_stars['ysensor'], 
                                 size=5)

    catalog_stars = hstack([catalog_stars, vis_stars])                     

    # create and save dataframe with star info, dropping rows with nan
    temp_stars_df = catalog_stars.to_pandas().dropna(subset=['peak_value'])

    _, file_name = os.path.split(file)
    
    # write filename
    temp_stars_df['source_file'] = file_name
    temp_stars_df['uv'] = 'uv' in file_name
    temp_stars_df['cutout'] = [
        sl.image_cutout(image_data, x, y, size=11) \
            for x, y in temp_stars_df[['x_peak', 'y_peak']].values
            ]
        
    stars_df = stars_df.append(temp_stars_df, ignore_index=True)
    
    # create dict with both fits headers, with extra info
    fits_hdrs = {**dict(hdul[0].header), **dict(hdul[1].header),
                 'UTC_time': spice.et2utc(et), 
                 'channel': 'vl' if 'vl' in file else 'uv',
                 'ra': ra,     # nominal RA from SPICE
                 'dec': dec,   # nominal dec from SPICE
                 'roll_angle': roll}
    
    # convert COMMENT and HISTORY to str, otherwise it won't convert to df
    fits_hdrs['COMMENT'] = str(fits_hdrs['COMMENT'])
    fits_hdrs['HISTORY'] = str(fits_hdrs['HISTORY'])
    
    img_info_df = img_info_df.append(fits_hdrs, ignore_index=True)
    
    hdul.close()
        
stars_df.to_pickle(base_path + 'spf_db.gz')

img_info_df.set_index('FILENAME', inplace=True)
img_info_df.to_pickle(base_path + 'fits_info_db.gz')

pio.renderers.default = "browser"

# plot stars in vl images
fig = px.scatter(stars_df, x='x_peak', y='y_peak', color='MAIN_ID',
                 facet_col='uv',
                 hover_name='source_file', 
                 hover_data=['MAIN_ID', 'mag'])

# decouple axes so I can change ranges separately
fig.update_xaxes(matches=None)
fig.update_yaxes(matches=None)

fig.update_xaxes(range=(0, 1024), row=1,col=1)
fig.update_yaxes(range=(0, 1024), row=1,col=1)
fig.update_xaxes(range=(0, 2048), row=1,col=2)
fig.update_yaxes(range=(0, 2048), row=1,col=2)

fig.show()
#fig.write_html("RSICW_all_stars_plot.html")

'''
# plot ra, dec of single star
stars_df = pd.read_pickle('../PFM_RSICW/spf_db.gz')
img_info_df = pd.read_pickle('../PFM_RSICW/PFM_RSICW_info.gz')
star = '* c Leo'
cleo=img_info_df.loc[stars_df[stars_df.MAIN_ID == star].source_file]
plt.scatter(cleo.ra, cleo.dec, c=cleo.roll_angle)
plt.colorbar(label='roll angle (deg)')
plt.ylabel('Dec (deg)')
plt.title('Metis boresight for ' + star +' in RSICW')
plt.xlabel('RA (deg)')

fig, axes = plt.subplots(5,5)

for s, ax in zip(stars_df[stars_df.MAIN_ID == star].iloc[0:25].iterrows(), axes.ravel()):
    s = s[1]
    extent = (s.xcentroid-11, s.xcentroid+11, s.ycentroid-11, s.ycentroid+11)
    ax.imshow(s.cutout, origin='lower', extent=extent)
    ax.plot(s.xcentroid, s.ycentroid, 'b+')
fig.show()    
'''