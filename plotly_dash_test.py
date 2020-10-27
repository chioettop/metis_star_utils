# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 18:18:51 2020

@author: chioettop
"""

import numpy as np
from spicelib import MetisSpice
import starlib as sl
from astropy.io import fits
import plotly.express as px

file = 'G:\Il mio Drive\METIS\Stellar fields identification\PFM_RSICW\solo_l0_metis-uv-image_0645545692_v01.fits'
hdul = fits.open(file)
image_data = hdul[0].data

fig = px.imshow(image_data, origin='lower')
vis_stars = sl.find_stars(image_data, fwhm=4)

for star in vis_stars:
    x, y = star['xcentroid'], star['ycentroid']
    fig.add_shape(
        type='circle',
        x0=x-5, x1=x+5, 
        y0=y-5, y1=y+5,
        xref='x', yref='y',
        line_color='cyan'
    )

import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()
app.layout = html.Div([
    dcc.Graph(figure=fig)
])

app.run_server(debug=True, use_reloader=False)  # Turn off reloader if inside Jupyter