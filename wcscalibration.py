# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:19:03 2020

@author: chioettop
"""
#%%
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

matched_stars = catalog_stars[catalog_stars['dist_to_SPF']<10]
cols = 4
rows = int(np.ceil(len(matched_stars)/cols))
rowscols = [(x+1, y+1) for x in range(rows) for y in range(cols)]
fig = make_subplots(rows=rows, cols=cols,
                    subplot_titles=matched_stars['MAIN_ID'])

for star, idx in zip(matched_stars, rowscols):
    vs = vis_stars[star['closest_SPF']]
    dist = np.hypot(vs['xcentroid'] - star['xsensor'], 
                    vs['ycentroid'] - star['ysensor']) 
    x = int(np.round(vs['xcentroid']))
    y = int(np.round(vs['ycentroid']))
    sub = image_data[y-5:y+6,x-5:x+6]
    fig.add_trace(
        go.Heatmap(z=sub, x0=x-5, y0=y-5, 
                   colorscale='gray', showscale=False),
        row=idx[0], col=idx[1]
        )
    
    # fig.update_yaxes(
    #     scaleanchor = "x",         # to work, should be x1, x2 etc.
    #     row=idx[0], col=idx[1]
    #   )

    fig.add_trace(
        go.Scatter(x=[vs['xcentroid']], y=[vs['ycentroid']],
                   mode='markers+text', marker_color='blue', marker_symbol='cross',
                   text=f"{dist:.1f}", 
                   textposition="bottom center",
                   textfont_color='blue'
 
                  ),
        row=idx[0], col=idx[1]
        )

    fig.add_trace(
        go.Scatter(x=[star['xsensor']], y=[star['ysensor']],
                   mode='markers+text', marker_color='orange', marker_symbol='star',
                   text=f"{star['mag']:.1f}", 
                   textposition="bottom center",
                   textfont_color='orange'
                  ),
        row=idx[0], col=idx[1]
        )


fig.update_layout(title_text="Side By Side Subplots", showlegend=False)
fig.show()
#%%
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
 
fit_st = matched_stars[[1, 2, 3, 4]]
world_coords = SkyCoord(fit_st['ra'], fit_st['dec'], frame="icrs", unit="deg")

fit_vis = vis_stars[fit_st['closest_SPF']]
xy = (fit_vis['xcentroid'], fit_vis['ycentroid'])

proj_point = SkyCoord(ra, dec, frame="icrs", unit="deg")

fit_wcs_from_points(
    xy=xy, 
    world_coords=world_coords, 
    proj_point=proj_point, 
    projection='TAN'
)

fit_wcs, fit_res = metis_fit_wcs(
    xy=xy,
    world_coords=world_coords,
    proj_point=proj_point,
    projection='TAN',
    UV=False
)

print("Fitted WCS")
print("Original pixel coords\n", xy)
print("Computed\n", fit_wcs.world_to_pixel(world_coords))
world_coords_new=fit_wcs.pixel_to_world(*xy)
print("Difference\n", world_coords.separation(world_coords_new))

print("Original WCS")
print("Computed pixel coords\n", wcs.world_to_pixel(world_coords))
world_coords_new=wcs.pixel_to_world(*xy)
print("Difference\n", world_coords.separation(world_coords_new))

# Singular Values Decomposition
#u, s, vh = np.linalg.svd(fit_wcs.wcs.cd)
#print('Scale: ', s[0], s[1])
