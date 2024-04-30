'''
GAIA COMMAND
'''
#%% Import modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import rasterio
from rasterio.transform import from_origin
from scipy.interpolate import griddata
import cv2
import copy
import time
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.gaia.gaia_util as gaia
import na1k.na1k_util as u1k
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Create reference grid for bc5k
gaia.CreateRefGrid_BC5k(meta)

#%% Prepare nitrogen deposition
gaia.Process_Ndep_ISIMIP(meta)

#%% Import CMIP6
gaia.Process_CMIP6(meta)

#%%
tv=np.arange(1850,2151,1)
d=gaia.ImportEnvironment(meta,tv)

plt.plot(tv,d['ndep'],'b-')
plt.plot(tv,d['co2'],'b-')

A=100*np.ones(d['co2'].size)
fA=gaia.fA(A,140,0.032,4.4)

fC=gaia.fCO2(d['co2'],[1,0.5,400])
plt.plot(tv,fC,'b-')
fC=gaia.fCO2(d['co2'],[1,0.3,400])
plt.plot(d['co2'],fC,'b-')
fC=gaia.fCO2(d['co2'],[1,0.15,400])
plt.plot(d['co2'],fC,'r--')


