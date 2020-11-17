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
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from functools import lru_cache

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.Div([
        dcc.Markdown(
            """
            # Metis Star Finder
            """
        )
    ]),
    
    html.Div([
        html.Label('Input image '),
        dcc.Input(
            id='image_file',
            type='url',
            placeholder="enter url of image to analyze",    
            debounce=True
        )]),
    
    html.Div(id='plot_controls')
        
    ])

@lru_cache(maxsize=32)
def load_image_data(image_file):
    hdul = fits.open(image_file)
    return hdul[0].data

@app.callback(
    Output('plot_controls', 'children'),
    [Input('image_file', 'value')]
)
def load_image(image_file):
    if image_file is not None:
        image_data = load_image_data(image_file)
        min_th = int(np.min(image_data))
        max_th = int(np.max(image_data))
        
        return [
            html.Label('FWHM'),
            dcc.Slider(
                id='fwhm',
                min=1,
                max=5,
                marks={i: str(i) for i in range(1, 6)},
                value=3,
                ),
            html.Label('Threshold'),
            dcc.Slider(
                id='threshold',
                min=min_th,
                max=max_th,
                marks={min_th:str(min_th), max_th:str(max_th)},
                value=3 * np.std(image_data),
                ),
            dcc.Graph(id='star_finder')
        ]

@app.callback(
    Output('star_finder', 'figure'),
    [Input('image_file', 'value'), Input('fwhm', 'value'), Input('threshold', 'value')]
)
def update_graph(image_file, fwhm, threshold):
    image_data = load_image_data(image_file)
    fig = px.imshow(image_data, origin='lower')
    vis_stars = sl.find_stars(image_data, fwhm=fwhm, threshold=threshold)
    
    if vis_stars:
        for star in vis_stars:
            x, y = star['xcentroid'], star['ycentroid']
            fig.add_shape(
                type='circle',
                x0=x-5, x1=x+5, 
                y0=y-5, y1=y+5,
                xref='x', yref='y',
                line_color='cyan'
            )
        
    return fig
 
app.run_server(debug=True)#, use_reloader=False)  # Turn off reloader if inside Jupyter