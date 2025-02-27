'''
GENERATE CUSTOM MASKS
'''

#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
from scipy.interpolate import griddata
import copy
from shapely.geometry import Polygon,Point,box,shape
import pyproj
import rasterio
from rasterio import features
from scipy.interpolate import NearestNDInterpolator
import fiona
import time
import cv2
import gzip
from datetime import datetime
import netCDF4 as nc
from scipy.interpolate import griddata
import scipy.ndimage
from perlin_noise import PerlinNoise
from numpy.random import default_rng
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.gaia.gaia_util as gaia
import fcexplore.field_plots.Processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def Mask_CDF(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z=gis.OpenGeoTiff(r'C:\Data\BC1ha\BEC_BIOGEOCLIMATIC_POLY\ZONE_GapFilled.tif')['Data']
	sz=gis.OpenGeoTiff(r'C:\Data\BC1ha\BEC_BIOGEOCLIMATIC_POLY\SUBZONE.tif')['Data']
	
	Mask=np.zeros(z.shape,dtype='int8')
	ind=np.where(zRef['Data']==1); Mask[ind]=1
	ind=np.where(z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']); Mask[ind]=2
	ind=np.where( (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (sz==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['SUBZONE']['xm'])  ); Mask[ind]=3
	plt.matshow(Mask)

	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(z.shape,dtype='int8')
	ind=np.where( (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']) | (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (sz==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['SUBZONE']['xm']))
	z1['Data'][ind]=1
	gis.SaveGeoTiff(z1,r'C:\Data\BC1ha\Masks\CDF\Mask.tif')
	return

#%%
def Mask_IDF(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z=gis.OpenGeoTiff(r'C:\Data\BC1ha\BEC_BIOGEOCLIMATIC_POLY\ZONE_GapFilled.tif')['Data']
	sz=gis.OpenGeoTiff(r'C:\Data\BC1ha\BEC_BIOGEOCLIMATIC_POLY\SUBZONE.tif')['Data']
	
	Mask=np.zeros(z.shape,dtype='int8')
	ind=np.where(zRef['Data']==1); Mask[ind]=1
	ind=np.where(z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']); Mask[ind]=2
	ind=np.where( (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (sz==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['SUBZONE']['xm'])  ); Mask[ind]=3
	plt.matshow(Mask)
	
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(z.shape,dtype='int8')
	ind=np.where( (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']) | (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BG']) | (z==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['PP']) )
	z1['Data'][ind]=1
	gis.SaveGeoTiff(z1,r'C:\Data\BC1ha\Masks\IDF\Mask.tif')
	return

#%%
def Mask_LowerThompson(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import watersheds
	pthin=meta['Paths']['GDB']['GDB'] + '\\WatershedsNamed.gdb'
	df=gpd.read_file(pthin,layer='FWA_NAMED_WATERSHEDS_POLY')

	# Extract Thompson
	df0=df[ (np.isin(df['STREAM_ORDER'],9)==True) & (df['GNIS_NAME']=='Thompson River')]
	df0=df0.reset_index()
	df0['ID']=1
	shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	plt.matshow(z0)

	# Remove upper watersheds
	df0=df[ (np.isin(df['STREAM_ORDER'],8)==True) & (df['GNIS_NAME']=='North Thompson River') | \
		(np.isin(df['STREAM_ORDER'],8)==True) & (df['GNIS_NAME']=='South Thompson River') | \
		(np.isin(df['STREAM_ORDER'],8)==True) & (df['GNIS_NAME']=='Nicola River') | \
		(np.isin(df['STREAM_ORDER'],8)==True) & (df['GNIS_NAME']=='Clearwater River')]
	df0=df0.reset_index()
	df0['ID']=1
	shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['ID']))
	z1=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z1,transform=zRef['Transform'])
	plt.matshow(z1)

	z2=copy.deepcopy(zRef)
	z2['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	z2['Data'][(z0==1) & (z1==0)]=1
	plt.matshow(z2['Data'])

	gis.SaveGeoTiff(z2,r'C:\Data\BC1ha\Masks\Lower Thompson\Mask.tif')
	return