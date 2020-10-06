# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:59:07 2020

@author: P. Chioetto 
"""
from starlib import *
from matplotlib.backends.backend_pdf import PdfPages
import glob
import pandas as pd
import warnings

def file_range(fits_range, base_dir, UV=False):
    # find files within a range
    start, stop = fits_range
    
    if start > stop:
        start, stop = stop, start
        
    files = glob.iglob(base_dir+'*' + ('uv' if UV else 'vl' + '*.fits'))
    
    found = []
    for file in files:
        beg = file.find('image_') + len('image_')
        end = file.find('_v', beg)
        seq = int(file[beg:end])
        if (seq <= stop) and (seq >= start):
            found.append(file)
            
    return found

#files = glob.glob("../PFM_RSICW/*.fits")
#fits_range = (645545693, 645755992) # alfa Leo
fits_range = (645742792, 645755992) # rho Leo
files = file_range(fits_range, base_dir="../PFM_RSICW/") # range of vl of 22Mb with a Leo visible

assert files, "No files found."

t = 0
no_stars = 0
too_many = 0
fit_warning = 0
data = []

with PdfPages(f'Metis PSF {fits_range}.pdf') as pdf:
    for file in files:
        hdul = fits.open(file)
        image_data = hdul[0].data
        
    
        # first page with whole image and stars found
        fig, ax = plt.subplots(figsize=(10,9))#, dpi=300)
        plot_image(image_data, ax=ax)
        with warnings.catch_warnings(record=True) as w:
            sources = find_stars(image_data, fwhm=1.5, plot=True, ax=ax)
            warning = ('\n'+str(w[0].message) if w else '')
            fig.suptitle(str(t) + ' - ' + hdul[0].header['FILENAME'] + warning)
            
        ax.set_title(f"DIT {hdul[1].header['DIT']}, VLFPFILT {hdul[1].header['VLFPFILT']}")
        pdf.savefig(fig)
        plt.close(fig)
       
        print(hdul[0].header['FILENAME'], ", DIT ", hdul[1].header['DIT'])
        
        if sources:
            print(f"\t Num sources = {len(sources)}")    
            
            for star in sources: 
                x, y = star['xcentroid'], star['ycentroid']
                # avoid stars too close to the border
                if (x>10) and (y>10):
                    with warnings.catch_warnings(record=True) as w:  # catch fitting warnings
                        f, fig = plot_gaussian_fit(image_data, x-5, y-5, 10, 10)
    
                        if w:
                            assert (str(w[0].message) == "The fit may be unsuccessful; check fit_info['message'] for more information.")
                            fit_warning += 1
                            data.append([file]+list(star)+f.parameters.tolist()+[True])
                        else:
                            data.append([file]+list(star)+f.parameters.tolist()+[False])
                        print(f"\t", ''.join(f"{param} = {val:.2f}, " for param, val in zip(f.param_names, f.parameters)))
                        
                        fig.suptitle(f"PSF {star['id']} - Fit {'warning' if w else 'successful'}")
                        pdf.savefig(fig)
                        plt.close(fig)

                        t += 1
                        #input()
            else:
                too_many += 1
        else:
            no_stars += 1    
                
        hdul.close()
        
print("Total images processed: ", t)

res = pd.DataFrame(data, columns=['file']+star.colnames+list(f.param_names)+['fit_warn'])
res.to_csv(f'Metis PSF {fits_range}.csv')

#%% 
fig, ax = plt.subplots()
fres = res[res.xcentroid < 1486]
y = (fres.ycentroid-1024)*10.2/3600
x_fwhm = np.abs(fres.x_stddev*2.355)
y_fwhm = np.abs(fres.y_stddev*2.355)
plt.plot(y[~res.fit_warn], x_fwhm[~res.fit_warn], '.', label='FWHM_x')
plt.plot(y[~res.fit_warn], y_fwhm[~res.fit_warn], '.', label='FWHM_y')
plt.plot(y[res.fit_warn], x_fwhm[res.fit_warn], 'x', label='FWHM_x (fit warning)')
plt.plot(y[res.fit_warn], y_fwhm[res.fit_warn], 'x', label='FWHM_y (fit warning)')

plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, borderaxespad=0.)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')

# guardare questo per i grafici mediati
#fres.groupby(pd.qcut(fres.ycentroid, 6)).y_stddev.mean()
#%%