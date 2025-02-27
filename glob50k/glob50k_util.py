#%% Import modules
import os
import numpy as np
import gc
import ee
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import rasterio
from rasterio.transform import from_origin
import fiona
import urllib
import time
import copy
import cv2
import scipy.io as spio
import tarfile
from scipy.interpolate import NearestNDInterpolator
import gzip
import shutil
import string
from scipy.interpolate import griddata
import scipy.ndimage
import netCDF4 as nc
from cartopy import crs as ccrs
from scipy.interpolate import griddata
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.gaia.gaia_util as gaia

#%% Initialize project
def Init(*argv):
	meta={}
	meta['Paths']={}
	meta['Paths']['glob50k']=r'C:\Data\Global50k'
	meta['Paths']['glob50k Ref Grid']=meta['Paths']['glob50k'] + '\\LandMask.tif'
	#meta['Paths']['glob50k Ref Grid Buf']=meta['Paths']['glob50k'] + '\\LandCoverUse\LandMaskBuf.tif'
	#meta['Paths']['glob50k Ref Grid Buf Rev']=meta['Paths']['glob50k'] + '\\LandCoverUse\LandMaskBufRev.tif'

	# Defaults assume province-wide map
	meta['Graphics']={}
	meta['Graphics']['Map']={}
	meta['Graphics']['Map']['RGSF']=1
	meta['Graphics']['Map']['Fig Width']=16
	meta['Graphics']['Map']['Side Space']=0
	meta['Graphics']['Map']['yxrat']=0.5
	meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.0,1]
	meta['Graphics']['Map']['Map Axis Vis']='off'
	meta['Graphics']['Map']['Map Grid Vis']=False
	meta['Graphics']['Map']['Legend X']=0.01
	meta['Graphics']['Map']['Legend Width']=0.03
	meta['Graphics']['Map']['Legend Font Size']=7
	meta['Graphics']['Map']['Legend Text Space']=0.035
	meta['Graphics']['Map']['Show Bound Land Mask']='On'
	meta['Graphics']['Map']['Show Bound Within']='Off'
	meta['Graphics']['Map']['Show Lakes']='Off'
	meta['Graphics']['Map']['Show Rivers']='Off'
	meta['Graphics']['Map']['Show Roads']='Off'
	meta['Graphics']['Map']['Show Rail']='Off'
	meta['Graphics']['Map']['Show Cities']='Off'
	meta['Graphics']['Map']['Show TPFs']='Off'
	meta['Graphics']['Map']['Show Symbol Labels']='Off'

	# Import vector data
	robinson=ccrs.Robinson().proj4_init
	gdf=gpd.read_file(r'C:\Data\Geodatabases\Global\World_Countries_Generalized.geojson')
	gdf=gdf.to_crs(robinson)

	return meta,gdf

#%%
def CreateReferenceGrid(meta):
	# Save a blank file with the desired extent, resolution and CRS
	Cellsize=25000
	X=np.arange(-12890000,15760000+Cellsize,Cellsize)
	Y=np.arange(-5960000,7796000+Cellsize,Cellsize)
	robinson=ccrs.Robinson().proj4_init
	crs=ccrs.Robinson().to_string
	z=np.random.randint(5,size=(Y.size,X.size),dtype='int16')
	with rasterio.open(
		meta['Paths']['glob50k Ref Grid'],
		mode="w",
		driver="GTiff",
		height=Y.size,
		width=X.size,
		count=1,
		dtype=z.dtype,
		crs=ccrs.Robinson(),
		transform=from_origin(X[0],Y[-1],Cellsize,Cellsize), # Inputs (west,north,xsize,ysize)
		) as new_dataset:
		new_dataset.write(z,1)

	# Reproject land cover
	fin=meta['Paths']['glob50k'] + '\\Downloads\\lcc_copern.tif'
	fout=meta['Paths']['glob50k'] + '\\lcc1_tmp.tif'
	gis.ReprojectGeoTiff(fin,fout,ccrs.Robinson())
	fout2=meta['Paths']['glob50k'] + '\\lcc1.tif'
	gis.ReprojectRasterAndClipToRaster(fout,fout2,meta['Paths']['glob50k Ref Grid'],ccrs.Robinson())

	z1=gis.OpenGeoTiff(fout2)
	z1=gis.ClipRasterByXYLimits(z1,[-12900000,15800000],[-6000000,7800000])
	gis.SaveGeoTiff(z1,fout2)

# 	zRef=gis.OpenGeoTiff(meta['Paths']['glob50k Ref Grid'])
# 	shapes=((geom,value) for geom, value in zip(meta['pb'].geometry,meta['pb']['ID']))
# 	z0=np.zeros(zRef['Data'].shape,dtype=float)
# 	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
# 	
# 	# Mask
# 	z1=copy.deepcopy(zRef)
# 	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
# 	ind=np.where(burned>0)
# 	z1['Data'][ind]=1
# 	plt.close('all');plt.matshow(z1['Data'])
# 	gis.SaveGeoTiff(z1,meta['Paths']['glob50k Ref Grid'])

	return

#%%
def Download_LandCover(meta):
	# Import canopy height
	ee.Authenticate(force=True)
	ee.Initialize(project='cdf-zone1')
	y0=ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global").select('discrete_classification').first()
	#proj0=y0.projection().getInfo()
	scale=10000
	# set geographic projection
	#crs='EPSG:3857'
	#y1=y0.reproject(crs=crs,scale=scale)
	task=ee.batch.Export.image.toDrive(image=y0,description='lcc_copern',scale=scale,maxPixels=1e13,fileFormat="GeoTIFF")
	task.start()
	# task.status()
	return