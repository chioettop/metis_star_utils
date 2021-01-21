import spiceypy
import numpy as np
import contextlib
import os
import sys

@contextlib.contextmanager
def chdir(path):
    CWD = os.getcwd()
    os.chdir(path)
    try:
        yield
    except:
        print('Exception caught: ', sys.exc_info()[0])
    finally:
        os.chdir(CWD)

class MetisSpice():
 
    # SOLAR ORBITER naif identifier
    solar_orbiter_naif_id = -144
 
    def __init__(self, meta_kernel, path='.'):
        self.path = path
        with chdir(self.path):
            spiceypy.furnsh(meta_kernel)
            
    def utc2et(self, utc_string):
        with chdir(self.path):
            return spiceypy.utc2et(utc_string)
       
    def scs2et(self, obt_string):
        with chdir(self.path):
            return spiceypy.scs2e(self.solar_orbiter_naif_id, obt_string)
 
    def obt2utc(self, obt_string):
        # Obt to Ephemeris time (seconds past J2000)
        with chdir(self.path):
            ephemeris_time = spiceypy.scs2e(self.solar_orbiter_naif_id, obt_string)
        # Ephemeris time to Utc
        # Format of output epoch: ISOC (ISO Calendar format, UTC)
        # Digits of precision in fractional seconds: 3
        return spiceypy.et2utc(ephemeris_time, "ISOC", 3)
 
    def utc2obt(self, utc_string):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc_string)
        # Ephemeris time to Obt
        with chdir(self.path):
            return spiceypy.sce2s(self.solar_orbiter_naif_id, ephemeris_time)
    
    def boresight(self, et, UV=False):
    # Returns Metis boresight and roll from ET time
     
        channel = 'SOLO_METIS_EUV_ILS' if UV else 'SOLO_METIS_VIS_ILS'
        with chdir(self.path):
            # rotation matrix from Metis local frame to J2000
            rot = spiceypy.pxform(channel, 'J2000', et) 
            # transform Metis boresight to J2000
            res = np.dot(rot, [-1, 0, 0])
            # returns [distance, RA, DEC] of vector
            radec = spiceypy.recrad(res)               
            # ra, dec in degrees
            ra, dec = np.rad2deg(radec[1:3])            
            # convert rot matrix to Euler angles 
            # the one corresponding to x axis is Metis roll
            _, _, roll = np.rad2deg(spiceypy.m2eul(rot, 3, 2, 1)) 
        
        return ra, dec, roll

    def sun_distance(self, et):
    # Returns Metis-Sun distance in Km from et time
        channel = 'SOLO_METIS_VIS_ILS'  # works also for the UV channel
        # sun state vector in Metis reference frame
        with chdir(self.path):
            sv = spiceypy.spkezr('SUN', et, channel, 'LT+S', 'SOLAR ORBITER')
        # sun position vector 
        pv = sv[0][0:3]   
        return np.linalg.norm(pv)
    
    def __del__(self):
        spiceypy.kclear()      
    
    def et2utc(self, et):
        return spiceypy.et2utc(et, 'C', 4)



