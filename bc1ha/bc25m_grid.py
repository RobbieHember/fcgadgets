
#%% Import modules

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point,Polygon
from shapely import geometry
import cv2
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcgadgets.bc1ha.bc25m_utilities as u25m

#%% Import paths and look-up-tables

# Initialize
meta=u25m.Init()

# Import look up tables
meta=u1ha.ImportLUTs(meta)

# Import graphics settings
gp=gu.SetGraphics('Manuscript')

#%% Define reference grid

flg=0
if flg==1:
    # Create reference grid
    # Manually copy BC1ha LandMask to bc25m database folder, then resample below
    gis.ResampleRaster(meta['Paths']['bc25m Ref Grid'],4)
zRef=gis.OpenGeoTiff(meta['Paths']['bc25m Ref Grid'])

#%% Reproject CEC Land Use Map 2010

fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
fout=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010c.tif'
gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

#z=gis.OpenGeoTiff(fin)

fin=r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
fout=meta['Paths']['bc25m'] + '\\LandCoverUse\\LandCover_CEC_2010.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc25m Ref Grid'],meta['Geos']['crs'])

z=gis.OpenGeoTiff(fout)
z=gis.ClipToRaster(z,zRef)

ind=np.where(zRef['Data']==0)
z['Data'][ind]=0

gis.SaveGeoTiff(z,)

#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\CEC_LandCover2010.tif')
#plt.matshow(z['Data'][0::5,0::5])

#z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif')

