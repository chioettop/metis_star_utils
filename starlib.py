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
from photutils import centroids
from photutils import DAOStarFinder
from photutils import find_peaks
from astropy.table import Table

from astropy.stats import sigma_clipped_stats


# common parameters
METIS_fov = 3.4 # METIS Field of View in deg, for astrocatalog searches
METIS_fov_min = 1.6 
SUN_RADIUS = 6.957e5   # Nominal Solar radius (km)
FWHM = 2.35482         # To convert from std to fwhm

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

def texp_from_fits_hdr(hdul):    
    # returns exposition time from Metis fits image header
    return hdul[1].header['DIT']/1000

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
    
    # remove stars outside Metis FOV
    sdist = (res['ra'] - ra)**2 + (res['dec'] - dec)**2
    res = res[sdist > METIS_fov_min**2]    

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
        # there is a bug in legend_elements that prevents showing reverse 
        # order scales
        '''handles, labels = scatter.legend_elements(prop="sizes", 
                                                  num=8, 
                                                  func=lambda s: max_mag-np.sqrt(s))
        ax.legend(handles, labels, title="Mag", bbox_to_anchor=(1.05, 1.0), loc='upper left')
        '''
    else:
        ids = [f"{ID}, {mag:.1f}" for ID, mag in zip(st['MAIN_ID'], st['mag'])]
        plot_stars(ids, x, y, ax=ax, color=color, radius=12)
    
    plt.xlabel('detector x')
    plt.ylabel('detector y')
        
    return ax

def wcs_from_boresight(ra, dec, roll, UV=False):
    '''
    Returns WCS coordinates transformation object from RA, Dec and Roll of Metis.
    The WCS converts from celestial and sensor coordinates (and vice versa).
    VL sensor coordinates have inverted x axis and are rotated 90 deg clockwise
    wrt. celestial.
    https://fits.gsfc.nasa.gov/fits_wcs.html

    Parameters
    ----------
    ra : float
        boresight right ascension in degrees.
    dec : float
        boresight declination in degrees.
    roll : float
        Metis roll angle in degrees.
    UV : boolean, optional
        Set to True for UV channel. The default is False.

    Returns
    -------
    w : WCS object
        WCS for transforming celestial coordinates to sensor coordinates.

    '''
    
    if UV:
        scx = -20.401/3600      # platescale in deg
        scy = 20.401/3600       
        bx, by = 512+2.6, 512-4.2  # boresight position in px
        det_angle = 0
        flip_xy = True          # UV detector appears to be flipped

    else:    
        scx = -10.137/3600  # platescale in deg (negative to invert axis)
        scy = 10.137/3600       
        bx = 966.075        # boresight position in px
        by = 2048-1049.130 
        det_angle = -90     # deg, detector base rotation wrt. sky
        flip_xy = False
    
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
    if flip_xy:
        w.wcs.pc = w.wcs.pc @ np.array([[0,-1],[-1,0]])
    
    return w      

def _linear_wcs_fit(params, lon, lat, x, y, w_obj):
    """
    Objective function for fitting linear terms.
    Parameters
    ----------
    params : array
        5 element array. First element is the roll angle (deg), 
        then CDELT, then CRPIX.
    lon, lat: array
        Sky coordinates.
    x, y: array
        Pixel coordinates
    w_obj: `~astropy.wcs.WCS`
        WCS object
        """
    roll = params[0]
    cdelt = params[1:3]
    crpix = params[3:5]

    t = -np.deg2rad(roll)  # following the sign convention in wcs_from_boresight()
    c, s = np.cos(t), np.sin(t)
    w_obj.wcs.pc = np.array([[c, -s], [s, c]]) # rotation matrix accounting for roll angle
    w_obj.wcs.cdelt = cdelt
    w_obj.wcs.crpix = crpix
    lon2, lat2 = w_obj.wcs_pix2world(x, y, 0)

    lat_resids = lat - lat2
    lon_resids = lon - lon2
    # In case the longitude has wrapped around
    lon_resids = np.mod(lon_resids - 180.0, 360.0) - 180.0

    resids = np.concatenate((lon_resids * np.cos(np.radians(lat)), lat_resids))

    return resids

def wcs_fit(xy, world_coords, proj_point,
                        projection='TAN', sip_degree=None,
                        UV=False):
    """
    From astropy.utils, modified to have CDELT and rot.matrix separated.
    
    Given two matching sets of coordinates on detector and sky,
    compute the WCS.
    Fits a WCS object to matched set of input detector and sky coordinates.
    Optionally, a SIP can be fit to account for geometric
    distortion. Returns an `~astropy.wcs.WCS` object with the best fit
    parameters for mapping between input pixel and sky coordinates.
    The projection type (default 'TAN') can passed in as a string, one of
    the valid three-letter projection codes - or as a WCS object with
    projection keywords already set. Note that if an input WCS has any
    non-polynomial distortion, this will be applied and reflected in the
    fit terms and coefficients. Passing in a WCS object in this way essentially
    allows it to be refit based on the matched input coordinates and projection
    point, but take care when using this option as non-projection related
    keywords in the input might cause unexpected behavior.
    Notes
    ------
    - The fiducial point for the spherical projection can be set to 'center'
      to use the mean position of input sky coordinates, or as an
      `~astropy.coordinates.SkyCoord` object.
    - Units in all output WCS objects will always be in degrees.
    - If the coordinate frame differs between `~astropy.coordinates.SkyCoord`
      objects passed in for ``world_coords`` and ``proj_point``, the frame for
      ``world_coords``  will override as the frame for the output WCS.
    - If a WCS object is passed in to ``projection`` the CD/PC matrix will
      be used as an initial guess for the fit. If this is known to be
      significantly off and may throw off the fit, set to the identity matrix
      (for example, by doing wcs.wcs.pc = [(1., 0.,), (0., 1.)])
    Parameters
    ----------
    xy : tuple of two `numpy.ndarray`
        x & y pixel coordinates.
    world_coords : `~astropy.coordinates.SkyCoord`
        Skycoord object with world coordinates.
    proj_point : ~astropy.coordinates.SkyCoord`
        A Skycoord object with a coordinate pair For consistency, the units 
        and frame of these coordinates will be transformed to match 
        ``world_coords`` if they don't.
    projection : str or `~astropy.wcs.WCS`
        Three letter projection code, of any of standard projections defined
        in the FITS WCS standard. Optionally, a WCS object with projection
        keywords set may be passed in.
    sip_degree : None or int
        If set to a non-zero integer value, will fit SIP of degree
        ``sip_degree`` to model geometric distortion. Defaults to None, meaning
        no distortion corrections will be fit.
    UV : boolean
        Used to set pixel_shape for Metis sensor size
    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        The best-fit WCS to the points given.
    """

    from astropy.coordinates import SkyCoord # here to avoid circular import
    import astropy.units as u
    from astropy.wcs import Sip
    from astropy.wcs.utils import celestial_frame_to_wcs
    from scipy.optimize import least_squares

    xp, yp = xy
    try:
        lon, lat = world_coords.data.lon.deg, world_coords.data.lat.deg
    except AttributeError:
        unit_sph = world_coords.unit_spherical
        lon, lat = unit_sph.lon.deg, unit_sph.lat.deg

    # verify input
    if type(proj_point) != type(world_coords):
        raise ValueError("proj_point must be set to an" +
                         "`~astropy.coordinates.SkyCoord` object with " +
                         "a pair of points.")
        assert proj_point.size == 1

    proj_codes = [
        'AZP', 'SZP', 'TAN', 'STG', 'SIN', 'ARC', 'ZEA', 'AIR', 'CYP',
        'CEA', 'CAR', 'MER', 'SFL', 'PAR', 'MOL', 'AIT', 'COP', 'COE',
        'COD', 'COO', 'BON', 'PCO', 'TSC', 'CSC', 'QSC', 'HPX', 'XPH'
    ]
    if type(projection) == str:
        if projection not in proj_codes:
            raise ValueError("Must specify valid projection code from list of "
                             + "supported types: ", ', '.join(proj_codes))
        # empty wcs to fill in with fit values
        wcs = celestial_frame_to_wcs(frame=world_coords.frame,
                                     projection=projection)
    else: #if projection is not string, should be wcs object. use as template.
        wcs = copy.deepcopy(projection)
        wcs.cdelt = (1., 1.) # make sure cdelt is 1
        wcs.sip = None
        
    # I am in fact using PC and CDELT. Not sure if this check should be reversed
    # Change PC to CD, since cdelt will be set to 1
    # if wcs.wcs.has_pc():
    #     wcs.wcs.cd = wcs.wcs.pc
    #     wcs.wcs.__delattr__('pc')

    if (type(sip_degree) != type(None)) and (type(sip_degree) != int):
        raise ValueError("sip_degree must be None, or integer.")

    # set pixel_shape to span of input points
    wcs.pixel_shape = (1024, 1024) if UV else (2048, 2048)

    # determine CRVAL from input
    if proj_point is not None:  # convert units, initial guess for crpix
        proj_point.transform_to(world_coords)
        wcs.wcs.crval = (proj_point.data.lon.deg, proj_point.data.lat.deg)
        wcs.wcs.crpix = (512, 512) if UV else (1024, 1024)

    # fit linear terms, assign to wcs
    # use (1, 0, 0, 1) as initial guess, in case input wcs was passed in
    # and cd terms are way off.
    # Use bounds to require that the fit center pixel is on the input image
    xpmin, xpmax, ypmin, ypmax = xp.min(), xp.max(), yp.min(), yp.max()
    if xpmin == xpmax:
        xpmin, xpmax = xpmin - 0.5, xpmax + 0.5
    if ypmin == ypmax:
        ypmin, ypmax = ypmin - 0.5, ypmax + 0.5
    
    # instead of a cd matrix, I use CDELT and a rotation matrix, with angle alpha
    roll = 0
    p0 = np.concatenate([[roll], wcs.wcs.cdelt.flatten(), wcs.wcs.crpix.flatten()])
    fit = least_squares(_linear_wcs_fit, p0,
                        args=(lon, lat, xp, yp, wcs))#,
#                        bounds=[[-180, -np.inf, -np.inf, xpmin, ypmin],
#                                [ 180, np.inf, np.inf, xpmax, ypmax]])

    t = -np.deg2rad(fit.x[0])  # following the sign convention in wcs_from_boresight()
    c, s = np.cos(t), np.sin(t)
    wcs.wcs.pc = np.array([[c, -s], [s, c]]) # rotation matrix accounting for roll angle
    wcs.wcs.cdelt = np.array(fit.x[1:3])
    wcs.wcs.crpix = np.array(fit.x[3:5])

    # fit SIP, if specified. Only fit forward coefficients
    if sip_degree:
        degree = sip_degree
        if '-SIP' not in wcs.wcs.ctype[0]:
            wcs.wcs.ctype = [x + '-SIP' for x in wcs.wcs.ctype]

        coef_names = [f'{i}_{j}' for i in range(degree+1)
                      for j in range(degree+1) if (i+j) < (degree+1) and
                      (i+j) > 1]
        p0 = np.concatenate((np.array(wcs.wcs.crpix), wcs.wcs.cd.flatten(),
                             np.zeros(2*len(coef_names))))

        fit = least_squares(_sip_fit, p0,
                            args=(lon, lat, xp, yp, wcs, degree, coef_names),
                            bounds=[[xpmin, ypmin] + [-np.inf]*(4 + 2*len(coef_names)),
                                    [xpmax, ypmax] + [np.inf]*(4 + 2*len(coef_names))])
        coef_fit = (list(fit.x[6:6+len(coef_names)]),
                    list(fit.x[6+len(coef_names):]))

        # put fit values in wcs
        wcs.wcs.cd = fit.x[2:6].reshape((2, 2))
        wcs.wcs.crpix = fit.x[0:2]

        a_vals = np.zeros((degree+1, degree+1))
        b_vals = np.zeros((degree+1, degree+1))

        for coef_name in coef_names:
            a_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[0].pop(0)
            b_vals[int(coef_name[0])][int(coef_name[2])] = coef_fit[1].pop(0)

        wcs.sip = Sip(a_vals, b_vals, np.zeros((degree+1, degree+1)),
                      np.zeros((degree+1, degree+1)), wcs.wcs.crpix)

    return wcs, fit

def remove_dark(image_data, t_exp, UV=False):
    if UV:
        return image_data
    
    dark = fits.getdata('dark_vlda_it2.fits')
    bias = fits.getdata('bias_it6b1_test.fits')
    return image_data - bias - dark*t_exp

def plot_image(image_data, ax=None):
    # clip bottom 10% and top 1%, excluding mask (0 pixels)
    low, high = scoreatpercentile(image_data, per=(10, 99), limit=(1,np.inf))
    
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
        image_data = remove_dark(image_data, texp_from_fits_hdr(hdul), "uv" in file)
        
    ax = plot_image(image_data, ax)
    title = hdul[0].header['FILENAME']
    if type(utc) is str:
        title += "\n" + utc
    ax.get_figure().suptitle(title)
    ax.set_title(coor_str + f"DIT {hdul[1].header['DIT']};")
    
    hdul.close()
    
    return ax
      
def plot_stars(ids, xs, ys, ax=None, color='blue', radius=10):
    
    if ax is None:
        fig, ax = plt.subplots()
        
    for id, x, y in zip(ids, xs, ys):
        ax.add_patch(plt.Circle((x, y), radius=radius, color=color, fill=False))
        ax.annotate(id, (x + 7, y + 7), color=color)


def find_stars(image_data, fwhm=3., threshold=None, ax=None, plot=False, fov=None):
  
    image_data = np.copy(image_data)
    mask = (image_data == 0)
    
    if fov:
        w, h = image_data.shape
        i, j = np.indices(image_data.shape, sparse=True)
        r = np.sqrt((i - w/2)**2 + (j - h/2)**2)
        mask = mask | (r < fov[0]) | (r > fov[1])
    
    image_data[mask] = np.nan   # the algorithm gets fooled by zeros
        
    mean, median, std = sigma_clipped_stats(image_data[~mask], sigma=3.0, mask_value=0.0)
    
    if not threshold:
        threshold = 3 * std
    
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(image_data - median)  
    
    if plot: 
        plot_stars(sources['id'], sources['x_centroid'], sources['y_centroid'], ax=ax)
        
    return sources

def image_cutout(image_data, x, y, size):
    """
    Returns a cutout at pixel x, y of size x size pixels

    Parameters
    ----------
    image_data : ndarray
        array with the image.
    x, y : float
        center of the cutout.
    size : int
        size of the cutout.

    Returns
    -------
    ndarray
        cutout of the image.

    """
    if np.isnan(x) | np.isnan(y):
        return None

    x, y = int(x), int(y)
    xmax, ymax = image_data.shape
    x0, y0 = x-size, y-size
    x1, y1 = x+size+1, y+size+1

    if (x0 < 0) | (y0 < 0) | (x1 > xmax) | (y1 > ymax):
        # cutout outside of boundaries
        return None
    
    return image_data[y0:y1,x0:x1]


def find_nearby_stars(image_data, xs, ys, size=10, 
                      centroid_func=centroids.centroid_2dg):
    """
    Finds a star nearby a given location in pixels, using a centroiding
    function if provided (eg. the ones in photutils.centroids)


    """
    
    assert xs.shape == ys.shape
    
    names = ['x_peak', 'y_peak', 'peak_value', 'x_centroid', 'y_centroid']
    stars = Table(names=names)
    
    for x, y in zip(xs, ys):
        # get a cutout of the probable star
        sub = image_cutout(image_data, x, y, size)
    
        if sub is None:
            stars.add_row([np.nan]*len(names))
        else:
    
            mean, median, std = sigma_clipped_stats(sub, sigma=3.0)
            f = find_peaks(sub-median, threshold=3*std, box_size=size*2, npeaks=1,
                           centroid_func=centroid_func)
            if f:
                # translate to global sensor coordinates
                f['x_peak'] += int(x) - size 
                f['y_peak'] += int(y) - size
                
                # centroids outside of the sub are clearly a bad fit
                if (f['x_centroid'] < 0)      | \
                   (f['y_centroid'] < 0)      | \
                   (f['x_centroid'] > size*2) | \
                   (f['y_centroid'] > size*2):
                       f['x_centroid'] = np.nan 
                       f['y_centroid'] = np.nan
                else:
                    f['x_centroid'] += int(x) - size 
                    f['y_centroid'] += int(y) - size 
                
                stars.add_row(dict(f))
            else:
                stars.add_row([np.nan]*len(names))
        
    return stars



def find_stars_sub(image_data, xstart, ystart, xdelta, ydelta, pos, fwhm=3., threshold=None, ax=None, plot=False):
    xend = int(xstart + xdelta)
    yend = int(ystart + ydelta)
    xstart, ystart = int(xstart), int(ystart)
           
    sub = image_data[ystart:yend+1, xstart:xend+1]

    image_data = np.copy(sub)
    mask = (image_data == 0)
    
    image_data[mask] = np.nan   # the algorithm gets fooled by zeros
        
    mean, median, std = sigma_clipped_stats(image_data[~mask], sigma=3.0, mask_value=0.0)
    
    if not threshold:
        threshold = 3 * std
    
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    # median to be added back to results?
    if sources:
        sources = daofind(image_data - median)
        sources['x_centroid'] += xstart
        sources['y_centroid'] += ystart
        sources['id'] += pos
    
        if plot: 
            plot_stars(sources['id'], sources['x_centroid'], sources['y_centroid'], ax=ax)
        
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
        cell_text.append([f"{f.x_stddev .value*FWHM:.2f}"])
        cell_text.append([f"{f.y_stddev.value*FWHM:.2f}"])
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
        
def plot_star_distance(s1, s2, side):
    x1, y1 = s1
    x2, y2 = s2
    fig = plt.figure()
    plt.suptitle('Star distance from PSF in pixels')

    distance = np.hypot(
        x2 - x1,
        y2 - y1
        )
    
    plt.title(f"Average = {np.round(np.mean(distance), 2)}")
    plt.axis('square')
    
    xlim = ylim = (0, side)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.quiver(x1, y1,
               x2 - x1,
               y2 - y1)
    for s, d in zip(s1, distance):
        plt.annotate(np.round(d, 2), (s[0], s[1]))
    plt.show()

def plot_fitted_psf(stars, image_data, size):
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.io as pio

    pio.renderers.default = "browser"

    # set up a 5 by x grid of subplots    
    cols = 6
    rows = int(np.ceil(len(stars)/cols))
    rowscols = [(x+1, y+1) for x in range(rows) for y in range(cols)]
    fig = make_subplots(
             rows=rows, cols=cols,
             subplot_titles=[f"{ID}, mag {mag:.1f}" \
                             for ID, mag in stars[['MAIN_ID', 'mag']] ]
             )

    for star, idx in zip(stars, rowscols):
        xs, ys = star['xsensor'], star['ysensor']
        x = int(xs)
        y = int(ys)
        sub = image_cutout(image_data, xs, ys, size)
        fig.add_trace(
            go.Heatmap(z=sub, x0=x-size, y0=y-size, 
                       colorscale='gray', showscale=False),
            row=idx[0], col=idx[1]
        )
           
        fig.add_trace(
            go.Scatter(x=[xs], y=[ys],
                       mode='markers', 
                       marker_color='blue', marker_symbol='star' 
                      ),
            row=idx[0], col=idx[1]
            )
        
        if ~np.isnan(star['peak_value']):
            xc, yc = star['x_centroid'], star['y_centroid']
            # xf, yf = np.abs(star['x_stddev'])*FWHM, np.abs(star['y_stddev'])*FWHM
            # fig.add_shape(
            #     dict(type="circle", line_color='orange',
            #          x0=xc-xf, y0=yc-yf, x1=xc+xf, y1=yc+yf),
            #     row=idx[0], col=idx[1]
            #     )
            fig.add_trace(
                go.Scatter(x=[xc], y=[yc],
                       mode='markers', 
                       marker_color='orange', marker_symbol='cross'
                       ),
                row=idx[0], col=idx[1]
                )
            fig.add_trace(
                go.Scatter(x=[star['x_peak']], y=[star['y_peak']],
                       mode='markers', 
                       marker_color='red', marker_symbol='cross'
                       ),
                row=idx[0], col=idx[1]
                )
    # make square pixels by setting scaleanchor for each y axis                       
    fig.for_each_yaxis(lambda yaxis: yaxis.update(dict(scaleanchor=yaxis.anchor)))
    
    fig.update_layout(title_text="Side By Side Subplots", showlegend=False)
    fig.show()
    
    
