# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:59:07 2020

@author: A. Slemer, P. Chioetto 
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.wcs
from scipy.stats import scoreatpercentile
from astroquery.simbad import Simbad
from photutils.centroids import fit_2dgaussian
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats


# common parameters
METIS_fov = 3.4 # METIS Field of View in deg, for astrocatalog searches
METIS_fov_min = 1.6 
SUN_RADIUS = 6.957e5   # Nominal Solar radius (km)

def obt_from_fits_hdr(hdul):    
    # returns ET time from Metis fits image header
    
    obt_beg = hdul[0].header['OBT_BEG'] # observation start time
    obt_end = hdul[0].header['OBT_END'] # observation end time 
    
    obt_avg = (obt_beg + obt_end) / 2
    
    # convert to SCLK string
    frac, whole = math.modf(obt_avg)   
    frac *= 65536.
    obt_str = str(np.int(whole))+':'+str(np.int(frac))
    
    return obt_str

def plot_fov_circle(center, radius, color, ax):
    # plot FOV circle on an existing ax in pixels
    
    fov = plt.Circle(center, radius, fill=None, color=color)
    ax.add_artist(fov)


def simbad_search(ra, dec, max_mag=6):
    
    # ref https://astroquery.readthedocs.io/en/latest/simbad/simbad.html    
    cs = Simbad()  # custom query fields
    cs.add_votable_fields('ra(d;A;ICRS;J2000;2000)', 'dec(d;D;ICRS;J2000;2000)') # ICRS decimal coords
    cs.add_votable_fields('flux(V)') # Visible flux
    cs.add_votable_fields('id(HD)')
    cs.remove_votable_fields('coordinates')
    # cs.get_votable_fields()
    # Simbad.get_field_description('id')
    
    # query SIMBAD database
    res = cs.query_criteria(f'region(circle, ICRS, {ra} +{dec}, {METIS_fov}d) & Vmag < {max_mag}')
    
    # rename columns for plotting routines
    res.rename_columns(['RA_d_A_ICRS_J2000_2000', 'DEC_d_D_ICRS_J2000_2000', 'FLUX_V'], ['ra', 'dec', 'mag'])
    
    # remove stars in the inner Metis FOV
    res = res[np.sqrt((res['ra'] - ra)**2 + (res['dec'] - dec)**2) > METIS_fov_min]    

    return res


def simbad_plot(st, wcs, ax=None, color='green', scatter=False, ref_frame='sky'):
               
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal', 'box')    
    
    # sets plot limits (sensor size)
    size = wcs.pixel_shape[0]
    xlim = ylim = (0, size)
    plt.xlim(xlim)
    plt.ylim(ylim)

    x, y = wcs.wcs_world2pix(st['ra'], st['dec'], 0)
    max_mag = np.max(st['mag'])
    sizes = (max_mag-st['mag'])**2 + 10
    
    if scatter:
        # draw scatter circles that do not rescale when zooming
        scatter = plt.scatter(x, y, c=color, s=sizes)
        # there is a bug in legend_elements that prevents showing reverse order scales
        '''handles, labels = scatter.legend_elements(prop="sizes", 
                                                  num=8, 
                                                  func=lambda s: max_mag-np.sqrt(s))
        ax.legend(handles, labels, title="Mag", bbox_to_anchor=(1.05, 1.0), loc='upper left')
        '''
    else:
        [ax.add_patch(plt.Circle((xi, yi), radius=10, color=color, fill=False)) \
                   for xi, yi in zip(x, y)]
    
    plt.xlabel('detector x')
    plt.ylabel('detector y')
    
    l_off = 5 # pixels
    
    for xp, yp, label, mag in zip(x, y, st['MAIN_ID'], st['mag']):
       if label:
          plt.annotate(f"{label.decode('UTF-8')}, {mag:.1f}", (xp+l_off, yp+l_off), color=color) # label is a byte string
          
    #plt.gca().invert_xaxis()  # in case this function is used alone
    
    return ax

def wcs_from_boresight(ra, dec, roll, UV=False):
    '''
    Returns WCS coordinates transformation object from RA, Dec and Roll of Metis.
    The WCS converts from celestial and sensor coordinates (and vice versa).
    VL sensor coordinates have inverted x axis and are rotated 90 deg clockwise
    wrt. celestial.
    https://fits.gsfc.nasa.gov/fits_wcs.html

    '''
    
    if UV:
        scx = 20/3600      # platescale in deg
        scy = 20/3600       
        bx, by = 512, 512  # boresight position in px
        det_angle = 0

    else:    
        scx = -10.137/3600  # platescale in deg (negative to invert axis)
        scy = 10.137/3600       
        bx = 966.075  # boresight position in px
        by = 2048-1049.130 
        det_angle = -90     # deg, detector base rotation wrt. sky 
    
    roll_offset = 0.077
    
    # build a WCS between sensor and sky coords
    # trasformation order pix->sky: matrix (rotation) -> translation -> scale
    w = astropy.wcs.WCS(naxis=2)
    w.pixel_shape = (1024, 1024) if UV else (2048, 2048)
    w.wcs.crpix = [bx, by]     # for boresight translation (center in pixels)         
    w.wcs.cdelt = [scx, scy]   # for scaling
    w.wcs.crval = [ra, dec]    # boresight in celestial coords
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]  # projection type from plane to spherical (TAN is gnomonic)
    w.wcs.cunit = ["deg", "deg"]            # unit for crdelt and crval
    t = -np.deg2rad(roll+roll_offset+det_angle)  # negative wrt the one returned by the rot matrix. Why?
    c, s = np.cos(t), np.sin(t)
    w.wcs.pc = np.array([[c, -s], [s, c]]) # rotation matrix accounting for roll angle
    
    return w      

def remove_dark(image_data, t_exp, UV=False):
    if UV:
        return image_data
    
    dark = fits.getdata('dark_vlda_it2.fits')
    bias = fits.getdata('bias_it6b1_test.fits')
    return image_data - bias - dark*t_exp

def plot_image(image_data, ax=None):
    low, high = scoreatpercentile(image_data, (10, 95))
    
    if ax is None:
        fig, ax = plt.subplots()
        
    im = ax.imshow(
        image_data,
        vmin=low, vmax=high, 
        cmap='gray', origin='lower')  
     
    plt.colorbar(im, ax=ax)
    
    return ax

def plot_fits(file, ax=None, coor=None, dark=True, utc=None, ref_frame='sky'):
    # plt.style.use(astropy_mpl_style)
    hdul = fits.open(file)
    image_data = hdul[0].data
    
    if ref_frame=='sky':
        image_data = np.rot90(image_data, 3)
    
    if coor:
        ra, dec = coor
        coor_str = f"RA {ra:.4f}, dec {dec:.4f}, "
    else:
        coor_str = ''
        
    if dark:
        image_data = remove_dark(image_data, hdul[1].header['DIT']/1000, "uv" in file)
        
    ax = plot_image(image_data, ax)
    title = hdul[0].header['FILENAME']
    if type(utc) is str:
        title += "\n" + utc
    ax.get_figure().suptitle(title)
    ax.set_title(coor_str + f"DIT {hdul[1].header['DIT']};")
    
    hdul.close()
    
    return ax
      
def find_stars(image_data, fwhm=3., threshold=None, ax=None, plot=True, fov=None):
  
    image_data = np.copy(image_data)
    mask = (image_data == 0)
    
    if fov:
        w, h = image_data.shape
        i, j = np.indices(image_data.shape, sparse=True)
        r = np.sqrt((i - w/2)**2 + (j - h/2)**2)
        mask = mask | (r < fov[0]) | (r > fov[1])
    
    image_data[mask] = np.nan   # the algorithm gets fooled by zeros
        
    mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, mask_value=0.0)
    
    if not threshold:
        threshold = 3 * std
    
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(image_data - median)  
    
    if plot: 
        if sources and (len(sources) < 25):
            if ax is None:
                fig, ax = plt.subplots()
                
            [ax.add_patch(plt.Circle((xi, yi), radius=12, color='blue', fill=False)) \
                   for xi, yi in zip(sources['xcentroid'], sources['ycentroid'])]

            for s in sources:
                ax.annotate(s['id'], (s['xcentroid'] + 7, s['ycentroid'] + 7), color='blue')
    
    return sources

def plot_gaussian_fit(image_data, xstart, ystart, xdelta, ydelta):
        xend = int(xstart + xdelta)
        yend = int(ystart + ydelta)
        xstart, ystart = int(xstart), int(ystart)
        
        # plot area of interest
        gridsize = (2, 3)
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid(gridsize, (0, 0))
        ax2 = plt.subplot2grid(gridsize, (0, 1))
        ax3 = plt.subplot2grid(gridsize, (0, 2))
        ax4 = plt.subplot2grid(gridsize, (1, 1))
        
        sub = image_data[ystart:yend+1, xstart:xend+1]
        vmin, vmax = np.min(sub), np.max(sub)
        im = ax1.imshow(sub, origin='lower', cmap='gray', vmin=vmin, vmax=vmax) # data
        ax1.set_title('Data')
        plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
        
        # fit 2D gaussian and plot
        f = fit_2dgaussian(sub)
        y, x = np.indices(sub.shape)
        im = ax2.imshow(f(x,y), origin='lower', cmap='gray', vmin=vmin, vmax=vmax) 
        ax2.set_title('Model')
        plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
        
        # residuals
        im = ax3.imshow(sub - f(x,y), origin='lower', cmap='gray')
        ax3.set_title('Residuals')
        plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
        
        #plt.tight_layout()
        #fig.colorbar(im, ax=axes.flat)      
        
        # print fitting params
        plt.setp(ax4, frame_on=False, xticks=(), yticks=())
        rowlabels = ["Centroid"] + list(f.param_names) + ["FWHM x", "FWHM y"]
        cell_text = [[f"{xstart + f.x_mean.value:.2f}, {ystart + f.y_mean.value:.2f}"]]
        cell_text += [[f"{val:.2f}"] for val in f.parameters]
        cell_text.append([f"{f.x_stddev .value*2.35482:.2f}"])
        cell_text.append([f"{f.y_stddev.value*2.35482:.2f}"])
        the_table = plt.table(cellText=cell_text, rowLabels=rowlabels, loc='center')
        
        plt.tight_layout()
        
        return f, fig     

def gaussiand2D_fit(data):
    # experimental. Using fit_2dgaussian from photutils instead
    # following this recipe: https://stackoverflow.com/questions/50522464/astropy-model-2dgaussian-issue
    
    from astropy.modeling import models, fitting
    
    data_m = data - np.mean(data) # subtract mean to shift the gaussian to the floor. May not work if the gaussian is large compared to the image
    
    y0, x0 = np.unravel_index(np.argmax(data), data.shape) # find coordinates of highest pixel
    sigma = np.std(data_m) 
    amp = np.max(data_m)

    p_init = models.Gaussian2D(amp, x0, y0, 1, 1) # using 1 for sigma seems to work better for the optimizer
    fit_p = fitting.LevMarLSQFitter()
    
    y, x = np.indices(data_m.shape)
    p = fit_p(p_init, x, y, data_m)
    plt.imshow(p(x, y), cmap='gray', origin='lower')
    

    
    
