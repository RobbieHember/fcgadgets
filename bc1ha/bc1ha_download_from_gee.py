'''
DOWNLOAD FROM GOOGLE EARTH ENGINE
'''

#%% Import modules
import ee
import numpy as np
import time
import matplotlib.pyplot as plt
import pandas as pd
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha

#%% Import canopy height
ee.Authenticate(force=True)
ee.Initialize(project='cdf-zone1')

#%% Canopy height (Lang, N., Jetz, W., Schindler, K., & Wegner, J. D. 2023). A high-resolution canopy height model of the Earth. Nature Ecology & Evolution, 1-12.
ch0=ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
proj0=ch0.projection().getInfo()

# Get canopy height at 100 m
geomPoly=ee.Geometry.BBox(-139.061,48.3,-114.05,60.01) # Province wide (East,south,west,north)
ch1=ch0.clip(geomPoly)

# Reproject to standard BC projection and scale at 100m
ch1=ch1.reproject(crs='EPSG:3005',scale=100)

# Start task
task=ee.batch.Export.image.toDrive(image=ch1,description='h_Lang22',scale=100,maxPixels=1e13,region=geomPoly,fileFormat="GeoTIFF") #
task.start()
# task.status()

# Manually move to bc1ha database, then reproject
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
fin=r'C:\Data\BC1ha\CanopyHeight\h_Lang22_raw.tif'
fout=r'C:\Data\BC1ha\CanopyHeight\h_Lang22.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

#%% Canopy height (Tolan et al. 2023)
ch0=ee.ImageCollection('projects/meta-forest-monitoring-okw37/assets/CanopyHeight').mosaic() # Tolan 2023
proj0=ch0.projection().getInfo()

# Get canopy height at 100 m
geomPoly=ee.Geometry.BBox(-139.061,48.3,-114.05,60.01) # Province wide (East,south,west,north)
ch1=ch0.clip(geomPoly)

# Reproject to standard BC projection and scale at 100m
ch1=ch1.reproject(crs='EPSG:3005',scale=100)

# Start task
task=ee.batch.Export.image.toDrive(image=ch1,description='h_Tolan23_raw',scale=100,maxPixels=1e13,region=geomPoly,fileFormat="GeoTIFF") #
task.start()
# task.status()

# Manually move to bc1ha database, then reproject
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
fin=r'C:\Data\BC1ha\CanopyHeight\h_Tolan23_raw.tif'
fout=r'C:\Data\BC1ha\CanopyHeight\h_Tolan23.tif'
gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])

#%%