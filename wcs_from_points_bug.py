# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 12:20:05 2020

@author: chioettop
"""
import numpy as np
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord

xy = np.array([[1766.88276168,  662.96432257,  171.50212526,  120.70924648],
               [1706.69832901, 1788.85480559, 1216.98949653, 1307.41843381]])

world_coords = SkyCoord([(66.3542367 , 22.20000162), (67.15416174, 19.18042906),
                         (65.73375432, 17.54251555), (66.02400512, 17.44413253)],
                        frame="icrs", unit="deg")

proj_point = SkyCoord(64.67514918, 19.63389538,
                      frame="icrs", unit="deg")

fit_wcs_from_points(
    xy=xy, 
    world_coords=world_coords, 
    proj_point=proj_point, 
    projection='TAN'
)