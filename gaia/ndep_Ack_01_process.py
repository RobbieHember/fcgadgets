
"""
PROCESS N DEPOSITION
"""

#%% IMPORT MODULES

import sys
import numpy as np
import gdal
import osr
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
from rasterio import features
#import cv2
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#%% Import data

# Ackerman et al. (2018) (https://conservancy.umn.edu/handle/11299/197613)
data=pd.read_csv(r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\inorganic_N_deposition.csv')

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

#%%

#nam='tot_1985'
nam='tot_2016'

lat=np.unique(data['latitude'].to_numpy())
lon=np.unique(data['longitude'].to_numpy())
z0=data[nam].to_numpy()/100
m=lat.size
n=lon.size
z0=np.flip(z0.reshape(n,m).T,0)

z=zTSA.copy()
z['Data']=z0
z['Data']=z['Data'].astype('float32')

srs=gis.ImportSRSs()
#z['Projection']=srs['Proj']['Geographic']
z['Proj4_String']=srs['String']['Geographic']

sr=osr.SpatialReference()
sr.ImportFromEPSG(4326)
z['Projection']=sr.ExportToWkt()

z['Y']=np.tile(np.reshape(lat,(-1,1)),(1,n))
z['X']=np.tile(np.reshape(lon,(1,-1)),(m,1))
z['m']=m
z['n']=n

# extent labels (left, right, bottom, top)
arra=np.array([np.min(lon),np.max(lat),np.min(lat),np.max(lat)])
z['Extent']=tuple(arra)

# Transform
arra=np.array([np.min(lon),2.5,0.0,np.max(lat),0.0,-2.0])
z['gt']=tuple(arra)

# Transform
z['Transform']=from_origin(np.min(lon),np.max(lat),2.5,2.0)

fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '.tif'
gis.SaveGeoTiff(z,fout)

#%% Clip to NA

fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '.tif'
z=gis.OpenGeoTiff(fout)

#plt.close('all')
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
#im=ax.matshow(z['Data'],clim=(0,12),extent=z['Extent'])

zC=gis.ClipRaster(z,[-140,-40],[0,90])

fout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip.tif'
gis.SaveGeoTiff(zC,fout)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
im=ax.matshow(zC['Data'],clim=(0,12),extent=zC['Extent'])


#%% Reproject

pthin=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip.tif'
pthout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip_projected.tif'
gis.ReprojectGeoTiff(pthin,pthout,srs['String']['BC1ha'])
#gis.ReprojectGeoTiff(pthin,pthout,srs['String']['NACID'])

z=gis.OpenGeoTiff(pthout)

z['X']=z['X'][:,0:-1]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
im=ax.matshow(z['Data'],clim=(0,12),extent=z['Extent'])

#%%

zND=zTSA.copy()

ivl=10
zND['X']=zTSA['X'][0::ivl,0::ivl]
zND['Y']=zTSA['Y'][0::ivl,0::ivl]
zND['Data']=griddata( (z['X'].flatten(),z['Y'].flatten()) ,z['Data'].flatten(),(zND['X'],zND['Y']),method='cubic')

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
im=ax.matshow(zND['Data'],clim=(0,5),extent=z['Extent'])




#%% Plot

# See bc1ha_map_full for plots
