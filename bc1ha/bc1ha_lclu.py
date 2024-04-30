#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from scipy.interpolate import griddata
import copy
from shapely.geometry import Polygon,Point,box,shape
import rasterio
from rasterio import features
from scipy.interpolate import NearestNDInterpolator
import fiona
import time
import cv2
import gzip
import netCDF4 as nc
from scipy.interpolate import griddata
import scipy.ndimage
from perlin_noise import PerlinNoise
from numpy.random import default_rng
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.gaia.gaia_util as gaia
import fcexplore.field_plots.Processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def DeriveLandCoverComp1(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['lc_vri_l3','lc_vri_l4','lc_cec_2020','lc_ntems_2019','harv_yr_cc','harv_yr_ntems']
	z=Import_Raster(meta,[],vList,'Extract Grid')
	lut_comp1=meta['LUT']['Derived']['lc_comp1']
	lut_ntems=meta['LUT']['Derived']['lc_ntems_2019']	
	lut_lc_vri_l3=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']
	lut_lc_vri_l4=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']
	#--------------------------------------------------------------------------
	# 2019 Land Cover
	# In the past, lc_vri_l2 treed area has differed from that of lc_vri_l4 (treatment of TFLs and private land), but that appears to be fixed in 2023
	#--------------------------------------------------------------------------	
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Coniferous'],lut_ntems['Broadleaf'],lut_ntems['Mixedwood'],lut_ntems['Wetland-treed']])==True) )
	z1['Data'][ind]=lut_comp1['Forest']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Shrubs']])==True) )
	z1['Data'][ind]=lut_comp1['Shrub']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Herbs']])==True) )
	z1['Data'][ind]=lut_comp1['Herb']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Wetland']])==True) )
	z1['Data'][ind]=lut_comp1['Wetland']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Bryoids']])==True) )
	z1['Data'][ind]=lut_comp1['Bryoid']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Rock and rubble'],lut_ntems['Exposed barren land']])==True) )
	z1['Data'][ind]=lut_comp1['Earth and Rock']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Snow and ice']])==True) )
	z1['Data'][ind]=lut_comp1['Snow and Ice']
	ind=np.where( (np.isin(z['lc_ntems_2019'],[lut_ntems['Water']])==True) )
	z1['Data'][ind]=lut_comp1['Water']
	
	# If NTEM says harvest and CEC says not settlement, assume forest
	ind=np.where( (z['harv_yr_ntems']>0) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Urban']) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Cropland']) )
	z1['Data'][ind]=lut_comp1['Forest']
	
	# Override with forest from VRI lc_vri_l4
	# *** Stopped doing this because it leads to unrealistically high forest area ***
	#ind=np.where( (np.isin(z['lc_vri_l4'],[lut_lc_vri_l4['TC'],lut_lc_vri_l4['TM'],lut_lc_vri_l4['TB']])==True) )
	#z1['Data'][ind]=lut_comp1['Forest']
	
	# Convert areas  where CEC map says urban
	ind=np.where( (z['lc_cec_2020']==meta['LUT']['Derived']['lc_cec_c']['Urban']) )
	z1['Data'][ind]=lut_comp1['Built-up']
	
	# CEC is more reliable indicator of cropland -> convert to shrubland
	ind=np.where( (z1['Data']==lut_comp1['Herb']) & (z['lc_cec_2020']!=meta['LUT']['Derived']['lc_cec_c']['Cropland']) )
	z1['Data'][ind]=lut_comp1['Shrub']
	
	# Check for any cells left unclassified
	# Very small number, likely water, reclassify as water
	ind=np.where( (zRef['Data']==1) & (z1['Data']==0) )
	print(ind[0].size)
	z1['Data'][ind]=lut_comp1['Water']

	# Check forest area
	ind=np.where( (z1['Data']==lut_comp1['Forest']) )
	A['For 2019']=ind[0].size/1e6

	# Patches get left as earth and rock, possibly associated with wildfire and
	# fire guards. Fix by converting to forest outside alpine areas
	zE=Import_Raster(meta,[],['elev'],'Extract Grid')
	
	MaskBuf=np.zeros(zRef['Data'].shape)
	ind=np.where(z['lc_vri_l3']==lut_lc_vri_l3['A'])
	MaskBuf[ind]=1
	MaskBuf=gis.BufferRasterMask(MaskBuf,5)
	ind=np.where( (z1['Data']==lut_comp1['Earth and Rock']) & (MaskBuf==0) & (zE['elev']<1450) )
	A_ToAdd=ind[0].size/1e6
	print(A_ToAdd)
	# a=np.zeros(zRef['Data'].shape); a[ind]=1
	# plt.matshow(a)
	z1['Data'][ind]=lut_comp1['Forest']

	# Reclassify outside land mask as water
	ind=np.where( (zRef['Data']==0) )
	z1['Data'][ind]=lut_comp1['Water']
	
	# Save
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')

	# Look at result
	flg=0
	if flg==1:
		plt.close('all')
		plt.matshow(z1['Data'])
		
		a=zRef['Data'].copy()
		ind=np.where(z1['Data']==lut_comp1['Herb'])
		a[ind]=2
		plt.matshow(a)	
	
		a=zRef['Data'].copy()
		ind=np.where(z['lc_cec_2020']['Data']==meta['LUT']['Derived']['lc_cec_c']['Cropland'])
		a[ind]=2
		plt.matshow(a)
	
	return

#%%
def DeriveLandUseComp1(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	lut_comp1=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	lut_ntems=meta['LUT']['Derived']['lc_ntems_2019']
	lut_lc_vri_l3=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']
	lut_cec=meta['LUT']['Derived']['lc_cec_c']
	
	vList=['lc_comp1_2019','lc_vri_l3','munic'] 
	z=Import_Raster(meta,[],vList,'Extract Grid')
	
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	
	#--------------------------------------------------------------------------
	# Conservation
	#--------------------------------------------------------------------------
	vList=['protected','park','ogma','ogdef','parknat','harv_yr_comp2','own']
	zC=Import_Raster(meta,[],vList,'Extract Grid')
	ind=np.where( (zC['harv_yr_comp2']==0) & (zC['protected']>0) | 
				 (zC['harv_yr_comp2']==0) & (zC['park']>0) |
				 (zC['harv_yr_comp2']==0) & (zC['parknat']>0) |
				 (zC['harv_yr_comp2']==0) & (zC['ogma']>0) |
				 (zC['harv_yr_comp2']==0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Conservancy Area, Ecological Reserve, Protected Area, Provincial Park']) |
				 (zC['harv_yr_comp2']==0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Local/Regional Park']) )
	z1['Data'][ind]=lut_lu['Conservation Natural']
	
	ind=np.where( (zC['harv_yr_comp2']>0) & (zC['protected']>0) | 
				 (zC['harv_yr_comp2']>0) & (zC['park']>0) |
				 (zC['harv_yr_comp2']>0) & (zC['parknat']>0) |
				 (zC['harv_yr_comp2']>0) & (zC['ogma']>0) |
				 (zC['harv_yr_comp2']>0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Conservancy Area, Ecological Reserve, Protected Area, Provincial Park']) |
				 (zC['harv_yr_comp2']>0) & (zC['own']==meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION']['Crown - Local/Regional Park']) )
	z1['Data'][ind]=lut_lu['Conservation Consistent']
	del zC
	
	#--------------------------------------------------------------------------
	# Cropland
	#--------------------------------------------------------------------------
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='OATS_ALR_POLYS')
	gdf['geometry']=gdf.geometry.buffer(10000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])	
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Herb']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Cropland']
	
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\crop_per.shp')
	#gdf['geometry']=gdf.geometry.buffer(10000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])	
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Herb'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Cropland']
	
	#--------------------------------------------------------------------------
	# Pasture
	#--------------------------------------------------------------------------
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='RANGE_PASTURE_POLY_SVW')
	gdf=gdf[(gdf['geometry'].area<1e8)]
	gdf['geometry']=gdf.geometry.buffer(5000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])	
	ind=np.where( (z1['Data']==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Pasture']
	
	#--------------------------------------------------------------------------
	# Residential and commercial
	#--------------------------------------------------------------------------
	zRC=np.zeros(zRef['Data'].shape,dtype='int8') # Keep track of all residential and commercial search areas
	z['munic'][0:3500,:]=0 # Remove rediculously big municipality in Fort Nelson
	#ind=np.where(z['munic']>0); z['munic'][ind]=1
	#z['munic']=gis.BufferRasterMask(z['munic'],1)
	#plt.matshow(z['munic'])
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up']])==True) & (z['munic']>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[z['munic']>0]=1
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(2000)
	#plt.close('all'); gdf.plot()
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])	
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[zBuf>0]=1
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POPULATED_PLACES_1M_SP')
	gdf['geometry']=gdf.geometry.buffer(2000)
	#plt.close('all'); gdf.plot()
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[zBuf>0]=1
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
	gdf['geometry']=gdf.geometry.buffer(500)
	#plt.close('all'); gdf.plot()
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform'])	
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[zBuf>0]=1
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='REC_TENURE_ALPINE_SKI_AREAS_SP')
	gdf['geometry']=gdf.geometry.buffer(50)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[zBuf>0]=1
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='CLAB_INDIAN_RESERVES')
	gdf['geometry']=gdf.geometry.buffer(1)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	zRC[zBuf>0]=1
	
	# Golf courses
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_GOLF_COURSES_SV')
	gdf['geometry']=gdf.geometry.buffer(500)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Residential and Commercial']
	
	#--------------------------------------------------------------------------
	# Energy and mines
	#--------------------------------------------------------------------------
	
	# Mines
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='HSP_MJR_MINES_PERMTTD_AREAS_SP')
	gdf['geometry']=gdf.geometry.buffer(500)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Energy and Mines']
	
	# Transmission lines
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\GBA_TRANSMISSION_LINES_SP.geojson')
	gdf['geometry']=gdf.geometry.buffer(200)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Energy and Mines']
	
	# BC Hydro facilities
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='DRP_BCHYDRO_FACILITY_PNTS_SP')
	gdf['geometry']=gdf.geometry.buffer(500)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Energy and Mines']
	
	# Oil and gas facilities
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\ogf.geojson')
	gdf['geometry']=gdf.geometry.buffer(200)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) &  (zBuf>0) )
	z1['Data'][ind]=lut_lu['Energy and Mines']
	
	#--------------------------------------------------------------------------		
	# Transportation
	#--------------------------------------------------------------------------   
	buffer_width=55 # metres
	
	# Major roads
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\road.geojson')
	gdf['geometry']=gdf.geometry.buffer(buffer_width)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	# Forestry roads
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\FTEN_ROAD_SEGMENT_LINES_SVW.geojson')
	gdf['geometry']=gdf.geometry.buffer(buffer_width)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	# Oil and gas roads
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\OG_ROAD_SEGMENT_PERMIT_SP.geojson')
	gdf['geometry']=gdf.geometry.buffer(buffer_width)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	# Oil and gas pipelines
	gdf=gpd.read_file(meta['Paths']['GDB']['GDB'] + '\\DRP_OIL_GAS_PIPELINES_BC_SP.geojson')
	gdf['geometry']=gdf.geometry.buffer(buffer_width)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']	
	
	# Airports
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_AIRPORTS_SVW')
	gdf['geometry']=gdf.geometry.buffer(500)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	# Ports
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_PORTS_TERMINALS_SVW')
	gdf['geometry']=gdf.geometry.buffer(500)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	#--------------------------------------------------------------------------
	# Water Management
	#--------------------------------------------------------------------------	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandCover'],layer='FWA_MANMADE_WATERBODIES_POLY')
	gdf['geometry']=gdf.geometry.buffer(50)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zBuf=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zBuf,transform=zRef['Transform']) 
	ind=np.where( (zRC==0) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Built-up'],lut_comp1['Earth and Rock'],lut_comp1['Bryoid'],lut_comp1['Shrub'],lut_comp1['Water']])==True) & (zBuf>0) )
	z1['Data'][ind]=lut_lu['Water Management']
	ind=np.where( (np.isin(z['lc_comp1_2019'],[lut_comp1['Water']])==True) )
	z1['Data'][ind]=lut_lu['Water Management']
	
	#--------------------------------------------------------------------------
	# Timber production
	#--------------------------------------------------------------------------
	ind=np.where( (z['lc_comp1_2019']==lut_comp1['Forest']) & (z1['Data']==0) )
	z1['Data'][ind]=lut_lu['Timber Production Passive']
	
	#--------------------------------------------------------------------------
	# Transporation 2
	# There's a lot of logging roads that don't get picked up. We can take advanage
	# of the change detected by CEC between 2010 and 2020, applying to areas
	# where no other designation is found above.
	#--------------------------------------------------------------------------
	vList=['lc_cec_2010','lc_cec_2020','lc_vri_l3'] 
	z2=Import_Raster(meta,[],vList,'Extract Grid')
	
	MaskA=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where(z2['lc_vri_l3']==lut_lc_vri_l3['A'])
	MaskA[ind]=1
	#plt.matshow(MaskA)
	MaskA=gis.BufferRasterMask(MaskA,5)
	#plt.matshow(MaskA2)
	#ind=np.where( (MaskA==0) & (MaskA2==1) )
	#MaskA[ind]=2
	
	ind=np.where( (z2['lc_cec_2010']==lut_cec['Forest']) & (z2['lc_cec_2020']!=lut_cec['Forest']) & (np.isin(z['lc_comp1_2019'],[lut_comp1['Shrub'],lut_comp1['Wetland'],lut_comp1['Herb']])==False) & (MaskA==0) & (z1['Data']==0) )
	z1['Data'][ind]=lut_lu['Transportation']
	
	#--------------------------------------------------------------------------
	# No designation
	#--------------------------------------------------------------------------
	ind=np.where( (z1['Data']==0) & (zRef['Data']==1) )
	z1['Data'][ind]=lut_lu['No Designation']
	
	#--------------------------------------------------------------------------
	# Non-land
	#--------------------------------------------------------------------------
	ind=np.where( (zRef['Data']==0) )
	z1['Data'][ind]=lut_lu['No Designation']+1
	
	#--------------------------------------------------------------------------
	# Save
	#--------------------------------------------------------------------------
	# plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])  
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2019.tif')
	
	return

#%% This needs to be run after the LU comp!
def DeriveLandCoverComp1_1800(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['lc_comp1_2019','lu_comp1_2019']
	z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	lc_comp1=meta['LUT']['Derived']['lc_comp1']
	lu_comp1=meta['LUT']['Derived']['lu_comp1']
	
	z1=copy.deepcopy(zRef)
	z1['Data']=z['lc_comp1_2019'].copy()
	ind=np.where( (z['lc_comp1_2019']!=lc_comp1['Forest']) & (np.isin(z['lu_comp1_2019'],[lu_comp1['Residential and Commercial'],lu_comp1['Energy and Mines'],lu_comp1['Transportation']])==True) )
	z1['Data'][ind]=lc_comp1['Forest']
	
	ind=np.where( (z['lc_comp1_2019']!=lc_comp1['Forest']) & (np.isin(z['lu_comp1_2019'],[lu_comp1['Cropland'],lu_comp1['Pasture']])==True) )
	rn=np.random.random(ind[0].size)
	indF=np.where(rn<0.7)[0]
	z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Forest']
	indF=np.where( (rn>=0.7) & (rn<0.9) )[0]
	z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Shrub']
	indF=np.where( (rn>=0.9) )[0]
	z1['Data'][ ind[0][indF],ind[1][indF] ]=lc_comp1['Herb']
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_1800.tif')
	return

#%% Derive land use compilation 2049 scenarios
def DeriveLandCoverLandUseComp1_2020to2049_Scenarios(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])		
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']		
	lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']
	vList=['lc_comp1_2019','lu_comp1_2019','harv_yr_comp2','ogdef','d2road','d2fac','d2set'] 
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	tv=np.arange(2020,2050,1)
	
	# Annual prob deforestation
	P_Def0=0.00065
	
	# Track area
	A={}
	A['For 2019']=np.where(z0['lc_comp1_2019']==lut_lc['Forest'])[0].size/1e6
	ind=np.where( (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Natural']) | (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) )
	A['Tot Con 2019']=ind[0].size/1e6
	ind=np.where( (z0['lc_comp1_2019']==lut_lc['Forest']) & (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Natural']) | (z0['lc_comp1_2019']==lut_lc['Forest']) & (zRef['Data']==1) & (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) )
	A['For Con 2019']=ind[0].size/1e6
	print('Current fraction of forest conserved:' + str(A['For Con 2019']/A['For 2019']))
	#A['Con 2049']=2.0*A['Con 2019']
	A['For Con 2049']=(0.3*A['For 2019'])
	A['OGD']=np.where( (z0['ogdef']>0) )[0].size/1e6
		
	#==========================================================================
	# LCLU-Comp1 2020to2049 Scenario 1: LUC-BAU / C-Hist / WM-None / FM-None
	#==========================================================================
	
	#ind=np.where( (z0['lu_comp1_2019']==lut_lu['Conservation Consistent']) & (z1['LU 2049']==lut_lu['Conservation Consistent']) )
	
	z1={}	
	z1['LC 2049']=z0['lc_comp1_2019'].copy()
	z1['LU 2049']=z0['lu_comp1_2019'].copy()
	z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
	z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
	for iT in range(tv.size):
		print(tv[iT])
		#--------------------------------------------------------------------------
		# Add LUC
		# Average NIR/PIR annual forest conversion = 7955ha/yr (0.01%/yr)
		# Modelling: x=np.arange(0,100,1); y=(1/(1+np.exp(0.14*(x-25)))); plt.close('all'); plt.plot(x,y,'r-')
		#--------------------------------------------------------------------------					
		rn0=np.random.random(zRef['Data'].shape)
		P_Def=P_Def0*(1/(1+np.exp(0.14*(z0['d2road']-25))))*(1/(1+np.exp(0.14*(z0['d2set']-25))))
		#plt.close('all'); plt.matshow(P_Def); plt.colorbar()	
		#ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (rn0<P_Def) )
		#A_Def=ind[0].size
		#print(A_Def/A['For 2019']*100)
		#rn=-1*np.ones(zRef['Data'].shape)
		rn1=np.random.random(zRef['Data'].shape)
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0) & (rn1<0.15) )
		z1['LC 2049'][ind]=lut_lc['Herb']
		z1['LU 2049'][ind]=lut_lu['Cropland']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-CL']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.15) & (rn1<0.45) )
		z1['LC 2049'][ind]=lut_lc['Shrub']
		z1['LU 2049'][ind]=lut_lu['Pasture']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-PA']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.45) & (rn1<0.51) )
		z1['LC 2049'][ind]=lut_lc['Built-up']
		z1['LU 2049'][ind]=lut_lu['Residential and Commercial']		
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-RC']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.51) & (rn1<0.85) )
		z1['LC 2049'][ind]=lut_lc['Shrub']
		z1['LU 2049'][ind]=lut_lu['Transportation']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-TR']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & (rn0<P_Def) & (rn1>=0.85) )
		z1['LC 2049'][ind]=lut_lc['Built-up']
		z1['LU 2049'][ind]=lut_lu['Energy and Mines']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-EM']
	
	#plt.close('all'); plt.matshow(z1['LU 2049']); plt.colorbar()	
	z2=copy.deepcopy(zRef); z2['Data']=z1['LC 2049']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn1.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['LU 2049']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn1.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['Type Change']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Type.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['Year Change']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn1_Year.tif')
	
	#==========================================================================
	# LCLU-Comp1 2020to2049 Scenario 2: LUC-BAU / C-CAP / FM-Hist
	#==========================================================================
	
	z1={}	
	z1['LC 2049']=z0['lc_comp1_2019'].copy()
	z1['LU 2049']=z0['lu_comp1_2019'].copy()
	z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
	z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
	
	# Add old growth deferrals
	ind1=np.where( (z0['harv_yr_comp2']==0) & (z0['ogdef']>0) )
	z1['LU 2049'][ind1]=lut_lu['Conservation Natural']
	z1['Year Change'][ind1]=2020
	ind2=np.where( (z0['harv_yr_comp2']>0) & (z0['ogdef']>0) )
	z1['LU 2049'][ind2]=lut_lu['Conservation Consistent']
	z1['Year Change'][ind2]=2020
	   
	for iT in range(tv.size):
		print(tv[iT])
		#--------------------------------------------------------------------------
		# Add LUC
		# Average NIR/PIR annual forest conversion = 7955ha/yr (0.01%/yr)
		# Modelling: x=np.arange(0,100,1); y=(1/(1+np.exp(0.14*(x-25)))); plt.close('all'); plt.plot(x,y,'r-')
		#--------------------------------------------------------------------------					
		rn0=np.random.random(zRef['Data'].shape)
		P_Def=P_Def0*(1/(1+np.exp(0.14*(z0['d2road']-25))))*(1/(1+np.exp(0.14*(z0['d2set']-25))))
		#plt.close('all'); plt.matshow(P_Def); plt.colorbar()	
		#ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & (rn0<P_Def) )
		#A_Def=ind[0].size
		#print(A_Def/A['For 2019']*100)
		#rn=-1*np.ones(zRef['Data'].shape)
		rn1=np.random.random(zRef['Data'].shape)
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
					 (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
					 (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
					 (rn0<P_Def) & (rn1>=0) & (rn1<0.15) )
		z1['LC 2049'][ind]=lut_lc['Herb']
		z1['LU 2049'][ind]=lut_lu['Cropland']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-CL']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
					 (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
					 (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
					 (rn0<P_Def) & (rn1>=0.15) & (rn1<0.45) )
		z1['LC 2049'][ind]=lut_lc['Shrub']
		z1['LU 2049'][ind]=lut_lu['Pasture']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-PA']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
					 (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
					 (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
					 (rn0<P_Def) & (rn1>=0.45) & (rn1<0.51) )
		z1['LC 2049'][ind]=lut_lc['Built-up']
		z1['LU 2049'][ind]=lut_lu['Residential and Commercial']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-RC']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
					 (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
					 (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
					 (rn0<P_Def) & (rn1>=0.51) & (rn1<0.85) )
		z1['LC 2049'][ind]=lut_lc['Shrub']
		z1['LU 2049'][ind]=lut_lu['Transportation']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-TR']
		
		ind=np.where( (z1['LC 2049']==lut_lc['Forest']) & \
					 (z0['lu_comp1_2019']!=lut_lu['Conservation Natural']) & (z0['lu_comp1_2019']!=lut_lu['Conservation Consistent']) & \
					 (z1['LU 2049']!=lut_lu['Conservation Natural']) & (z1['LU 2049']!=lut_lu['Conservation Consistent']) & \
					 (rn0<P_Def) & (rn1>=0.85) )
		z1['LC 2049'][ind]=lut_lc['Built-up']
		z1['LU 2049'][ind]=lut_lu['Energy and Mines']
		z1['Year Change'][ind]=tv[iT]
		z1['Type Change'][ind]=lut_chng['FL-EM']
	
	#--------------------------------------------------------------------------
	# Add Conservation
	#--------------------------------------------------------------------------

	# Add small areas
	seedval=10
	rng=default_rng(seed=seedval)
	sig=21
	noiseS=rng.integers(0,1000,zRef['Data'].shape,'int16',True)
	blurS=cv2.GaussianBlur(noiseS,(0,0),sigmaX=sig,sigmaY=sig,borderType=cv2.BORDER_DEFAULT)
	# plt.close('all'); plt.matshow(blurS)
	u=np.unique(blurS)
	#print(u)
	
	MaskS=np.zeros(blurS.shape,dtype='int8')
	ind=np.where(blurS<=u[15]); MaskS[ind]=1
	plt.close('all'); plt.matshow(MaskS)
	ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z1['LU 2049'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) ); #Mask2[ind]=1
	A_PropS=ind[0].size/1e6
	print(A_PropS)
	print( (A['For Con 2019']+A['OGD']+A_PropS)/A['For Con 2049'] )
	
	# Add new areas
	ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z0['lu_comp1_2019'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) & (z0['harv_yr_comp2']==0) )	
	z1['LU 2049'][ind]=lut_lu['Conservation Natural']
	ind=np.where( (MaskS==1) & (zRef['Data']==1) & (np.isin(z0['lu_comp1_2019'],[lut_lu['Timber Production Passive'],lut_lu['Timber Production Intense']])==True) & (z0['harv_yr_comp2']>0) )
	z1['LU 2049'][ind]=lut_lu['Conservation Consistent']
	
	# Save
	z2=copy.deepcopy(zRef); z2['Data']=z1['LC 2049']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['LU 2049']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['Type Change']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn2_Type.tif')
	z2=copy.deepcopy(zRef); z2['Data']=z1['Year Change']
	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2020to2049_Scn2_Year.tif')
   
	#==========================================================================
	# LCLU-Comp1 2020to2049 Scenario 3: LUC-BAU / C-CAP / FM-CAP1
	#==========================================================================	
	# Import required data
	zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
	zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
	zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
	#z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
	#z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
	
	# Create Perlin noise
	noise=PerlinNoise(octaves=10,seed=3)
	xpix,ypix=1000,1000
	zNoise0=np.zeros((xpix,ypix))
	for i in range(xpix):
		for j in range(ypix):
			zNoise0[i,j]=noise([i/xpix, j/ypix])			
	zNoise2=np.zeros(zNoise0.shape); 
	ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
	zNoise2[ind]=1
	zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
	zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
	# plt.close('all'); plt.matshow(zNoise2)

	#--------------------------------------------------------------------------
	# Add Bands of intense forestry zone
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(7*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zInner=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(55*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zOuter=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(56*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zOuter2=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zOuter2,transform=zRef['Transform'])
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(6*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zInner2=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zInner2,transform=zRef['Transform'])
	
	# zCb=np.zeros(zRef['Data'].shape,dtype='int8')
	# ind=np.where( (zF==1) & (zOuter==1) )
	# zCb[ind]=1
	# plt.close('all'); plt.matshow(zCb)
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
	gdf['geometry']=gdf.geometry.buffer(75*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zF=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
	
	zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
	zMaskIFM[ind]=1
	#plt.close('all'); plt.matshow(zMaskIFM)
		
	# Add intensee forest management
	ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
	zLU['Data'][ind]=lut_lu['Timber Production Intense']
	
	# Add wildfire management zones
	ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zInner==1) )
	zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
	
	# Save	
	gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn3.tif')
   
	#==========================================================================
	# LCLU-Comp1 2020to2049 Scenario 4: LUC-BAU / C-CAP / FM-CAP2
	#==========================================================================	
	# Import required data
	zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
	zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
	zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
	
	#--------------------------------------------------------------------------
	# Add Bands of intense forestry zone
	noise=PerlinNoise(octaves=10,seed=3)
	xpix,ypix=1000,1000
	zNoise0=np.zeros((xpix,ypix))
	for i in range(xpix):
		for j in range(ypix):
			zNoise0[i,j]=noise([i/xpix, j/ypix])			
	zNoise2=np.zeros(zNoise0.shape); 
	ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
	zNoise2[ind]=1
	zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
	zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
	# plt.close('all'); plt.matshow(zNoise2)

	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(7*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zInner=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(55*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zOuter=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
	gdf['geometry']=gdf.geometry.buffer(75*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zF=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
	
	zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
	zMaskIFM[ind]=1
	#plt.close('all'); plt.matshow(zMaskIFM)
		
	# Add intensee forest management
	ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
	zLU['Data'][ind]=lut_lu['Timber Production Intense']
	
	# Add wildfire management zones
	ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zInner==1) )
	zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
	
	#--------------------------------------------------------------------------
	# Add bioenergy production
	
	# Fibrous pattern for bioenergy 
	noise=PerlinNoise(octaves=35,seed=3)
	xpix,ypix=1000,1000
	zNoise0=np.zeros((xpix,ypix))
	for i in range(xpix):
		for j in range(ypix):
			zNoise0[i,j]=noise([i/xpix, j/ypix])
	zNoise1=np.zeros(zNoise0.shape); 
	ind=np.where(np.abs(zNoise0)<0.05); 
	zNoise1[ind]=1	
	zNoise1=np.tile(zNoise1,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
	zNoise1=zNoise1[0:zRef['m'],0:zRef['n']]  
	#plt.close('all'); plt.matshow(zNoise1)	
	
	gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	gdf['geometry']=gdf.geometry.buffer(60*1000)
	gdf=gdf[gdf.geometry!=None]
	gdf=gdf.reset_index()
	gdf['ID']=np.ones(len(gdf))
	shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	zOuterBE=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=zOuterBE,transform=zRef['Transform'])
	
	Mask1=np.zeros(zRef['Data'].shape,dtype='int8')
	ListLU=[lut_lu['Timber Production Intense'],lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<20) & (zNoise1==1) & (zOuterBE==1) )
	Mask1[ind]=1
	plt.close('all'); plt.matshow(Mask1)
	print(ind[0].size/1e6)
	zLU['Data'][ind]=lut_lu['Bioenergy Production']
	
	# Save	
	gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn4.tif')
	
	
	
	
	# # Import required data
	# zLC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2049_Scn2.tif')
	# zLU=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn2.tif')
	# zSlope=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
	# #z1['Year Change']=np.zeros(zRef['Data'].shape,dtype='int16')
	# #z1['Type Change']=np.zeros(zRef['Data'].shape,dtype='int8')
	
	# # Create Perlin noise
	# noise=PerlinNoise(octaves=35,seed=3)
	# xpix,ypix=1000,1000
	# zNoise0=np.zeros((xpix,ypix))
	# for i in range(xpix):
	#	 for j in range(ypix):
	#		 zNoise0[i,j]=noise([i/xpix, j/ypix])
	
	# # Fibrous pattern for bioenergy 
	# zNoise1=np.zeros(zNoise0.shape); 
	# ind=np.where(np.abs(zNoise0)<0.05); 
	# zNoise1[ind]=1	
	# zNoise1=np.tile(zNoise1,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
	# zNoise1=zNoise1[0:zRef['m'],0:zRef['n']]  
	# #plt.close('all'); plt.matshow(zNoise1)
	
	# # Exterior
	# noise=PerlinNoise(octaves=10,seed=3)
	# xpix,ypix=1000,1000
	# zNoise0=np.zeros((xpix,ypix))
	# for i in range(xpix):
	#	 for j in range(ypix):
	#		 zNoise0[i,j]=noise([i/xpix, j/ypix])			
	# zNoise2=np.zeros(zNoise0.shape); 
	# ind=np.where( (zNoise0<-0.12) | (zNoise0>0.12) ); 
	# zNoise2[ind]=1
	# zNoise2=np.tile(zNoise2,(int(np.ceil(zRef['m']/1000)),int(np.ceil(zRef['n']/1000))))
	# zNoise2=zNoise2[0:zRef['m'],0:zRef['n']]  
	# # plt.close('all'); plt.matshow(zNoise2)
	
	# #--------------------------------------------------------------------------
	# # Add bioenergy production	
	# Mask1=np.zeros(zRef['Data'].shape,dtype='int8')
	# ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	# ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<20) & (zNoise1==1) )
	# Mask1[ind]=1
	# #plt.close('all'); plt.matshow(Mask1)
	# print(ind[0].size/1e6)
	# zLU['Data'][ind]=lut_lu['Bioenergy Production']

	# # #--------------------------------------------------------------------------
	# # # Wildfire management
	# # zMaskWM=np.zeros(zRef['Data'].shape,dtype=float)
	# # # Add FN communities	
	# # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='FN_COMMUNITY_LOCATIONS_SP')
	# # gdf['geometry']=gdf.geometry.buffer(2*1000)	
	# # gdf=gdf[gdf.geometry!=None]
	# # gdf=gdf.reset_index()
	# # gdf['ID']=np.ones(len(gdf))
	# # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))	
	# # burned=features.rasterize(shapes=shapes,fill=0,out=zMaskWM,transform=zRef['Transform'])		
	# # ind=np.where( (zMaskWM>0) )		
	# # zMaskWM[ind]=1	
	# # # Add major cities
	# # gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# # gdf['geometry']=gdf.geometry.buffer(2*1000)
	# # gdf=gdf[gdf.geometry!=None]
	# # gdf=gdf.reset_index()
	# # gdf['ID']=np.ones(len(gdf))
	# # shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))	
	# # burned=features.rasterize(shapes=shapes,fill=0,out=zMaskWM,transform=zRef['Transform'])		
	# # # Add to LU map
	# # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	# # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zSlope<30) & (zMaskWM==1) & (zNoise2==1) )
	# # zLU['Data'][ind]=lut_lu['Wildfire Risk Management']

	# #--------------------------------------------------------------------------
	# # Add Bands of intense forestry zone
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(7*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zInner=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zInner,transform=zRef['Transform'])
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(55*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zOuter=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zOuter,transform=zRef['Transform'])
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(56*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zOuter2=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zOuter2,transform=zRef['Transform'])
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	# gdf['geometry']=gdf.geometry.buffer(6*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zInner2=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zInner2,transform=zRef['Transform'])
	
	# zCb=np.zeros(zRef['Data'].shape,dtype='int8')
	# ind=np.where( (zF==1) & (zOuter==1) )
	# zCb[ind]=1
	# plt.close('all'); plt.matshow(zCb)
	
	# gdf=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='GSR_TMBR_PRCSSING_FAC_SV')
	# gdf['geometry']=gdf.geometry.buffer(75*1000)
	# gdf=gdf[gdf.geometry!=None]
	# gdf=gdf.reset_index()
	# gdf['ID']=np.ones(len(gdf))
	# shapes=((geom,value) for geom, value in zip(gdf['geometry'],gdf['ID']))
	# zF=np.zeros(zRef['Data'].shape,dtype=float)
	# burned=features.rasterize(shapes=shapes,fill=0,out=zF,transform=zRef['Transform'])
	
	# zMaskIFM=np.zeros(zRef['Data'].shape,dtype='int8')
	# ind=np.where( (zInner==0) & (zOuter==1) & (zF==1) )
	# zMaskIFM[ind]=1
	# #plt.close('all'); plt.matshow(zMaskIFM)
	
	# zMaskWMRing=np.zeros(zRef['Data'].shape,dtype='int8')
	# ind=np.where( (zInner==0) & (zOuter==0) & (zOuter2==1) )
	# zMaskWMRing[ind]=1
	# ind=np.where( (zInner2==0) & (zInner==1) )
	# zMaskWMRing[ind]=1
	# #plt.close('all'); plt.matshow(zMaskWMRing)
		
	# # Add intensee forest management
	# ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	# ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskIFM==1) & (zNoise2==1) )
	# zLU['Data'][ind]=lut_lu['Timber Production Intense']
	
	# # # Add outter protective WM ring
	# # ListLU=[lut_lu['Timber Production Passive'],lut_lu['No Designation']]
	# # ind=np.where( (zLC['Data']==lut_lc['Forest']) & (np.isin(zLU['Data'],ListLU)==True) & (zMaskWMRing==1) )
	# # zLU['Data'][ind]=lut_lu['Wildfire Risk Management']
	
	# # Save	
	# gis.SaveGeoTiff(zLU,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandUse_Comp1_2049_Scn3.tif') 
   
	return

#%% Derive deforestation mask from land use compilation 1
def DeriveLandCoverLandUseComp1_DeforstationMask(meta):	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])		
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	vList=['lc_comp1_1800','lc_comp1_2019','lc_comp1_2049s1','lc_comp1_2049s2']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	z1=copy.deepcopy(zRef); 
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where( (z0['lc_comp1_2049s1']!=lut_lc['Forest']) & (z0['lc_comp1_1800']==lut_lc['Forest']) | (z0['lc_comp1_2049s2']!=lut_lc['Forest']) & (z0['lc_comp1_1800']==lut_lc['Forest']) )
	z1['Data'][ind]=1
	# plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandCover_Comp1_DeforestationMaskAll.tif')
	return

#%%
def ExportSummaryLCLU(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	
	vList=['lc_comp1_1800','lc_comp1_2019','lu_comp1_2019'] 
	z=Import_Raster(meta,[],vList)
	for k in z.keys(): 
		z[k]=z[k]['Data']
	
	bc=np.array(list(lut_lc.values()))
	bu=np.array(list(lut_lu.values()))
	
	# Matrix of 1800 land cover / 2019 land cover
	M=np.zeros((bc.size,bc.size) )
	for i in range(bc.size):
		for j in range(bc.size):
			ind=np.where( (zRef['Data']==1) & (z['lc_comp1_1800']==bc[i]) & (z['lc_comp1_2019']==bu[j]) )		
			if ind[0].size>0:
				M[i,j]=ind[0].size	
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lc.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandCover2019.xlsx')
	
	# Matrix of 2019 land cover / 2019 land use
	M=np.zeros((bu.size,bc.size) )
	for i in range(bu.size):
		for j in range(bc.size):
			ind=np.where( (z['lc_comp1_2019']==bc[j]) & (z['lu_comp1_2019']==bu[i]) )
			if ind[0].size>0:
				M[i,j]=ind[0].size	
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lu.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover2019_LandUse2019.xlsx')
	
	# Matrix of 1800 land cover / 2019 land use
	M=np.zeros( (bu.size,bc.size) )
	for i in range(bu.size):
		for j in range(bc.size):
			ind=np.where( (z['lc_comp1_1800']==bc[j]) & (z['lu_comp1_2019']==bu[i]) )		
			if ind[0].size>0:
				M[i,j]=ind[0].size	
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lu.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandUse2019.xlsx')	
	
	return

#%%
def ExportSummaryLCLU(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	
	vList=['lc_comp1_1800','lc_comp1_2019','lu_comp1_2019','fire_2023']
	z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

	bc=np.array(list(lut_lc.values()))
	bu=np.array(list(lut_lu.values()))
	
	# Matrix of 1800 land cover / 2019 land cover
	M=np.zeros((bc.size,bc.size) )
	for i in range(bc.size):
		for j in range(bc.size):
			ind=np.where( (zRef['Data']==1) & (z['lc_comp1_1800']==bc[i]) & (z['lc_comp1_2019']==bu[j]) )		
			if ind[0].size>0:
				M[i,j]=ind[0].size	
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lc.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandCover2019.xlsx')
	
	# Matrix of 2019 land cover / 2019 land use (*** FOR 2023 WILDFIRE AREA***)
	M=np.zeros((bu.size,bc.size) )
	for i in range(bu.size):
		for j in range(bc.size):
			ind=np.where( (z['lc_comp1_2019']==bc[j]) & (z['lu_comp1_2019']==bu[i]) & (z['fire_2023']>0) )
			if ind[0].size>0:
				M[i,j]=ind[0].size
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lu.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover2019_LandUse2019_Wildfire2023.xlsx')
	
	# Matrix of 1800 land cover / 2019 land use
	M=np.zeros( (bu.size,bc.size) )
	for i in range(bu.size):
		for j in range(bc.size):
			ind=np.where( (z['lc_comp1_1800']==bc[j]) & (z['lu_comp1_2019']==bu[i]) )
			if ind[0].size>0:
				M[i,j]=ind[0].size	
	df=pd.DataFrame(M,columns=lut_lc.keys())
	df['LC']=lut_lu.keys()
	df.to_excel(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover1800_LandUse2019.xlsx')	
	
	return

#%% Plot land use change time series
def Plot_LUC_TS(meta):
	tv=np.arange(1851,2050,1)
	vList=['luc1_hist_yr','luc1_hist_type','luc1_fut_s1_yr','luc1_fut_s1_type','luc1_1019_type','luc1_1019_yr'] # 
	z0=Import_Raster(meta,[],vList,'Extract Grid')
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']
	cL1=['Cropland','Pasture','Residential and Commercial','Transportation','Energy and Mines','Water Management','No Designation']
	#cL2=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND']
	
	tsHist={}
	for c in lut_chng.keys():
		tsHist[c]=np.zeros(tv.size)
		ind0=np.where( (z0['luc1_hist_yr']>0) & (z0['luc1_hist_type']==lut_chng[c]) )
		if ind0[0].size==0:
			continue
		idx=gu.IndicesFromUniqueArrayValues(z0['luc1_hist_yr'][ind0])
		for i in idx.keys():
			iT=np.where(tv==i)[0]
			tsHist[c][iT]=idx[i].size
		
	tsFut={}
	for c in lut_chng.keys():
		tsFut[c]=np.zeros(tv.size)
		ind0=np.where( (z0['luc1_fut_s1_yr']>0) & (z0['luc1_fut_s1_type']==lut_chng[c]) )
		if ind0[0].size==0:
			continue
		idx=gu.IndicesFromUniqueArrayValues(z0['luc1_fut_s1_yr'][ind0])
		for i in idx.keys():
			iT=np.where(tv==i)[0]
			tsFut[c][iT]=idx[i].size
	   
	ts1019={}
	for c in lut_chng.keys():		
		ts1019[c]=np.zeros(tv.size)
		ind0=np.where( (z0['luc1_1019_yr']>0) & (z0['luc1_1019_type']==lut_chng[c]) )
		if ind0[0].size==0:
			continue
		idx=gu.IndicesFromUniqueArrayValues(z0['luc1_1019_yr'][ind0])
		for i in idx.keys():
			iT=np.where(tv==i)[0]
			ts1019[c][iT]=idx[i].size
	
	cl=np.array([[0.95,0.9,0.7],[0.75,0.9,0.55],[0.85,0.4,0.4],[0,0.3,0.7],[0.65,0.1,0],[0.4,0.6,0.8],[0.75,0.75,0.75]])
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,6))
	yC=np.zeros(tv.size); cnt=0
	for c in lut_chng.keys():  
		ax.bar(tv,tsHist[c]/1000,1,bottom=yC,facecolor=cl[cnt,:],label=c)
		yC=yC+tsHist[c]/1000
		cnt=cnt+1
	cnt=0
	for c in lut_chng.keys():  
		ax.bar(tv,tsFut[c]/1000,1,bottom=yC,facecolor=cl[cnt,:])
		yC=yC+tsFut[c]/1000
		cnt=cnt+1
	# cnt=0
	# for c in lut_chng.keys():  
	#	 ax.bar(tv,ts1019[c]/1000,1,bottom=yC,facecolor=cl[cnt,:])
	#	 yC=yC+ts1019[c]/1000
	#	 cnt=cnt+1
	ax.set(yscale='linear',ylabel='Deforestation (Kha/year)',xlabel='Time, years',xticks=np.arange(1800,2140,10),ylim=[0,80],xlim=[1860,2050])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
	ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	plt.tight_layout()
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LUC_Comp1_TS','png',900)

	A_Pasture2019_All=np.where( (z0['lu_comp1_2019']==lut_lu['Pasture']) )[0].size/1e3
	A_Pasture2019_For=np.where( (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Pasture']) )[0].size/1e3
	
	A_NetDelta=(tsHist['FL-PA']+ts1019['FL-PA'])/1000
	A_GrossClearing=1.3*A_NetDelta.copy(); iTaGC=np.where(tv>=1945)[0]; A_GrossClearing[iTaGC]=20.5
	A_WoodyEncroachment=np.zeros(tv.size)
	A_NetDeltaSim=np.zeros(tv.size)
	A_Pasture=np.zeros(tv.size)
	for i in range(1,tv.size):			
		A_WoodyEncroachment[i]=0.028*A_Pasture[i-1]
		A_Pasture[i]=A_Pasture[i-1]+A_GrossClearing[i]-A_WoodyEncroachment[i]
		A_NetDeltaSim[i]=A_GrossClearing[i]-A_WoodyEncroachment[i]		
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,9))
	ax.bar(tv,np.cumsum(A_NetDelta),1,facecolor=[0.8,0.8,0.8],label='Reconstruction from current observations')
	ax.plot(2019,A_Pasture2019_All,'k^',ms=6,mfc=[1,0.5,0],mec=[0,0,0],mew=0.25,label='Pature area total (observed)')
	ax.plot(2019,A_Pasture2019_For,'ks',ms=6,mfc=[0.5,0,1],mec=[0,0,0],mew=0.25,label='Pature area, originally forest (observed)')		
	#ax.plot(tv,np.cumsum(A_NetDeltaSim),'b--',label='Net conversion, simulated')
	ax.set(yscale='linear',ylabel='Pasture area (Kha)',xlabel='Time, years',xticks=np.arange(1800,2140,10),xlim=[1900,2023])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
	ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25)		
	ax2=ax.twinx(); ax2.plot([1800,2100],[0,0],'k-',lw=0.5)		
	iT1=np.where( (tv>=2010) & (tv<=2018) )[0]; ax2.plot(tv[iT1],A_NetDelta[iT1],'k-',lw=2,label='Conversion (CEC Maps)')
	iT1=np.where( (tv<2010) )[0]; 
	ax2.plot(tv[iT1],A_NetDelta[iT1],'k--',label='Net conversion from reconstruction')
	ax2.plot(tv,A_NetDeltaSim,'b--',label='Net conversion, simulated')
	ax2.plot(tv,A_GrossClearing,'g--',color=[0.6,0.9,0],label='Gross clearing, simulated')
	ax2.plot(tv,-A_WoodyEncroachment,'r-',label='Woody encroachment')		
	ax2.set(ylabel='Annual change (Kha/year)',ylim=[-30,40])
	ax2.legend(loc='lower left',frameon=False,facecolor=[1,1,1])
	plt.tight_layout()
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LCLU_PastureForestDynamics','png',900)
	return

#%%
def CompareLandCover(meta):

	vList=['lc_comp1_2019','lc_comp1_vri','refyear']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

	ind=np.where( (z0['lc_comp1_2019']==1) | (z0['lc_comp1_vri']==1) )
	x=z0['refyear'][ind]
	a=z0['lc_comp1_2019'][ind]
	b=z0['lc_comp1_vri'][ind]

	y=np.zeros(ind[0].size)
	ind2=np.where( (a==1) & (b==1) )[0]
	y[ind2]=1
	ind2=np.where( (a==1) & (b!=1) )[0]
	y[ind2]=2
	ind2=np.where( (a!=1) & (b==1) )[0]
	y[ind2]=0

	bw=1; bin=np.arange(1970,2023,bw)
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,10));
	ax.plot([1950,2040],[1,1],'k-',lw=2.5,color=[0.8,0.8,0.8])
	ax.plot(bin,mu,'ko')
	ax.set(yticks=np.arange(0,2,0.1),ylabel='Agreement index',xticks=np.arange(1500,2200,5),xlabel='Time, years',ylim=[0.5,1.5],xlim=[1967,2025])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	ax.text(2000,1.4,'NTEMS = forest,\nVRI = non-forest',va='center',fontsize=10)
	ax.text(2000,0.6,'NTEMS = non-forest,\nVRI = forest',va='center',fontsize=10)

	return

#%%
def DeriveLandUseChangeYear_1800to2019(meta):
	
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['lc_comp1_1800','lc_comp1_2019','lu_comp1_2019','lc_cec_2010','lc_cec_2020']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	lut_lc=meta['LUT']['Derived']['lc_comp1']
	lut_lu=meta['LUT']['Derived']['lu_comp1']
	lut_cec=meta['LUT']['Derived']['lc_cec_c']
	lut_chng=meta['LUT']['Derived']['lclu_chng_comp1']
	
	cL1=['Cropland','Pasture','Residential and Commercial','Transportation','Energy and Mines','Water Management','No Designation']
	cL2=['FL-CL','FL-PA','FL-RC','FL-TR','FL-EM','FL-WM','FL-ND']
	
	A_For_1800=np.where(z0['lc_comp1_1800']==lut_lc['Forest'])[0].size/1e6
	A_For_2019=np.where(z0['lc_comp1_2019']==lut_lc['Forest'])[0].size/1e6
	A_For_Delta=A_For_1800-A_For_2019
	
	# Create rings from populated locations
	df=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='BC_MAJOR_CITIES_POINTS_500M')
	df=df[df['CITY_TYPE']!='DM']
	df=df.reset_index()
	df['Seed Year']=1865*np.ones(len(df['BCMJ_TAG']))
	dSet=gu.ReadExcel(r'C:\Data\Cities\BC_MAJOR_CITIES_POINTS_500M.xlsx')	
	for i in range(len(dSet['BCMJ_TAG'])):
		ind=np.where(df['BCMJ_TAG']==dSet['BCMJ_TAG'][i])[0]; 
		df['Seed Year'][ind]=dSet['Year Settled'][i]
	df=df[df['Seed Year']>0]
	df=df[df.geometry!=None]
	df=df.reset_index()
	
	tv=np.arange(1851,2020,1)
	binA=9000*(1/(1+(np.exp(-0.16*(tv-1945)))))
	binGA=np.diff(binA)
	plt.close('all'); plt.plot(tv,binA,'bo')
	plt.close('all'); plt.plot(tv[1:],binGA,'bo')
	binR=np.sqrt(binA/np.pi) #A=pi*r^2
	binRG=np.diff(binR)
	plt.close('all'); plt.plot(tv,binR,'bo')
	#plt.close('all'); plt.plot(binRG,'ko')   
	
	zR=np.zeros(zRef['Data'].shape,dtype='float')
	zY=np.zeros(zRef['Data'].shape,dtype='int16')
	for iT in range(tv.size):
		print(binR[iT])
		ind=np.where(df['Seed Year']<=tv[iT])[0]
		df0=df.loc[ind].copy()
		df0['geometry']=df0.geometry.buffer(1000*binR[iT])
		z2=np.zeros(zRef['Data'].shape,dtype=float)
		shapes=((geom,value) for geom, value in zip(df0['geometry'],df0['BCMJ_TAG']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z2,transform=zRef['Transform'])
		ind2=np.where( (burned>0) & (zR==0) )
		zR[ind2]=binR[iT]
		zY[ind2]=tv[iT]
	# indNZ=np.where(zY>0)
	# idx=gu.IndicesFromUniqueArrayValues(zY[indNZ])
	
	# #plt.close('all'); plt.matshow(zY,clim=[1850,2010])	   
	# A=np.zeros(tv.size)
	# for i in idx.keys():
	#	 iT=np.where(tv==i)[0]
	#	 A[iT]=idx[i].size
	# plt.close('all'); plt.plot(tv,A,'bo')
	
	# LUC over 1800-2019	
	zYear=copy.deepcopy(zRef); zYear['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zType=copy.deepcopy(zRef); zType['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iT in range(tv.size):		
		print(tv[iT])
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Cropland']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-CL']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Pasture']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-PA']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Residential and Commercial']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-RC']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Transportation']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-TR']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Energy and Mines']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-EM']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['Water Management']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-WM']
		ind=np.where( (zY==tv[iT]) & (z0['lc_comp1_1800']==lut_lc['Forest']) & (z0['lu_comp1_2019']==lut_lu['No Designation']) )
		zYear['Data'][ind]=tv[iT]; zType['Data'][ind]=lut_chng['FL-ND']
	A_For_Delta_Prime=np.where(zYear['Data']>0)[0].size/1e6
	zYearBK=copy.deepcopy(zYear)
	zTypeBK=copy.deepcopy(zType)
	#zYear=copy.deepcopy(zYearBK)
	#zType=copy.deepcopy(zTypeBK)	
	
	# Fill unclassified areas, but not areas with deforestaion over 2010-2019
	# Assume random uniform distribution over post-war period	
	for j in range(len(cL1)):
		ind=np.where( (zYear['Data']==0) &				  
				 (z0['lc_comp1_1800']==lut_lc['Forest']) & 
				 (z0['lc_comp1_2019']!=lut_lc['Forest']) &
				 (z0['lu_comp1_2019']==lut_lu[cL1[j]]) )
		zYear['Data'][ind]=np.random.randint(1955,2020,ind[0].size)
		zType['Data'][ind]=zType['Data'][ind]=lut_chng[cL2[j]]
	A_For_Delta_Prime=np.where(zYear['Data']>0)[0].size/1e6	
	#(zDef1019_Year['Data']==0) &
	
	idx=gu.IndicesFromUniqueArrayValues(zYear['Data'][zYear['Data']>0])
	tv=np.arange(1851,2020,1)
	A=np.zeros(tv.size)
	for i in idx.keys():
		iT=np.where(tv==i)[0]
		A[iT]=idx[i].size
	plt.close('all'); plt.plot(tv,A/1000,'-bo')
	
	# Deforestation indicated over 2010-2019
	zDef1019_Year=copy.deepcopy(zRef); zDef1019_Year['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zDef1019_Type=copy.deepcopy(zRef); zDef1019_Type['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where( (z0['lc_cec_2010']==meta['LUT']['Derived']['lc_cec_c']['Forest']) & 
				 (np.isin(z0['lu_comp1_2019'],[lut_lu['Cropland'],lut_lu['Pasture'],lut_lu['Residential and Commercial'],lut_lu['Energy and Mines'],lut_lu['Water Management'],lut_lu['Transportation'],lut_lu['No Designation']])==True) & 
				 (z0['lc_comp1_1800']==lut_lc['Forest']) & 
				 (z0['lc_comp1_2019']!=lut_lc['Forest']) )
	zDef1019_Year['Data'][ind]=np.random.randint(2010,2019,ind[0].size)
	zDef1019_Type['Data'][ind]=zType['Data'][ind]
	A_Def1019=np.where(zDef1019_Year['Data']>0)[0].size/1e6
	
	print(A_For_Delta)
	print(A_For_Delta_Prime)
	print(A_Def1019)	
	
	# Add year of transition to conservation
	zYear=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Year.tif')
	ind=np.where(z0['lu_comp1_2019']==lut_lu['Conservation Natural']); zYear['Data'][ind]=np.random.randint(1975,high=2019,size=ind[0].size)
	ind=np.where(z0['lu_comp1_2019']==lut_lu['Conservation Consistent']); zYear['Data'][ind]=np.random.randint(1975,high=2019,size=ind[0].size)
	
	# Save	
	gis.SaveGeoTiff(zYear,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Year.tif')
	gis.SaveGeoTiff(zType,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_1800to2019_Type.tif')
	gis.SaveGeoTiff(zDef1019_Year,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Year.tif')
	gis.SaveGeoTiff(zDef1019_Type,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Comp1_2010to2019_Type.tif')
		
	return

#%% Reclassify VRI land cover into LCC Comp 1 classes
def ReclassifyVRI_LandCover(meta):
	vList=['refg','lc_vri_l2','lc_vri_l3','lc_vri_l4']
	z=Import_Raster(meta,[],vList)	
	z1=z['refg'].copy()
	z1['Data']=np.zeros(z['refg']['Data'].shape,dtype='int8')
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Forest']
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SI'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['ST']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Shrub']
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HE'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HF'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['HG']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Herb']
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BM'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['BY']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Bryoid']	
	ind=np.where( (np.isin(z['lc_vri_l3']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_3']['W']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Wetland']	
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['EL'],meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['RO']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Earth and Rock']
	ind=np.where( (np.isin(z['lc_vri_l4']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['SI']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Snow and Ice']
	ind=np.where( (np.isin(z['lc_vri_l2']['Data'],[meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['W']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']
	ind=np.where( (z1['Data']==0) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']+1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2022_FromVRI2022.tif')
	return

#%% Prepare NTEMS data
def PrepareNTEMS(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	# Harvest year
	fin=r'C:\Data\Harvest\NTEMS\harv85to20.tif'
	a=gis.OpenGeoTiff(fin)
	fout=r'C:\Data\BC1ha\Disturbances\Harvest_NTEM_Year.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout)
	z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0)
	z['Data'][ind]=0
	gis.SaveGeoTiff(z,fout)
	# Land cover 
	fin=r'C:\Data\Land Cover\NTEMS\vlce2_2019.tif'
	fout=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout)
	z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0)
	z['Data'][ind]=0
	gis.SaveGeoTiff(z,fout)
	# Age
	fin=r'C:\Data\Age\NTEMS\age_ntem_c.tif'
	fout=meta['Paths']['bc1ha'] + '\\Age\\Age_NTEM_2019.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	return

#%% Reclassify NTEMS land cover into LCC Comp 1 classes
def ReclassifyNTEMS_LandCover(meta):

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	vList=['lc_ntems_2019']
	z0=Import_Raster(meta,[],vList)['lc_ntems_2019']['Data'] 

	z1=copy.deepcopy(zRef)
	z1['Data']=(meta['LUT']['Derived']['lc_comp1']['Water']+1)*np.ones(zRef['Data'].shape,dtype='int8')

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Coniferous'],meta['LUT']['Derived']['lc_ntems_2019']['Broadleaf'],meta['LUT']['Derived']['lc_ntems_2019']['Mixedwood'],meta['LUT']['Derived']['lc_ntems_2019']['Wetland-treed']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Forest']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Shrubs']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Shrub']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Herbs']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Herb']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Wetland']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Wetland']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Bryoids']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Bryoid']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Rock and rubble'],meta['LUT']['Derived']['lc_ntems_2019']['Exposed barren land']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Earth and Rock']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Snow and ice']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Snow and Ice']

	ind=np.where( (np.isin(z0,[meta['LUT']['Derived']['lc_ntems_2019']['Water']])==True) )
	z1['Data'][ind]=meta['LUT']['Derived']['lc_comp1']['Water']
	
	plt.close('all'); plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_NTEMS_2019_ReclassAsComp1.tif')
	return

#%% NALCMS CEC 2020 - Compress categories to remove tropics
def NALCMS_Compress(lut_in,zRef,zCEC10,zCEC20):

	cg=['Forest','Shrubland','Grassland','Lichen-Moss','Wetland','Cropland','Barren Ground','Urban','Water','Snow and Ice']
	d={'Value':(len(cg)+1)*np.ones(len(cg),dtype='uint8'),'Name':np.array(['empty' for _ in range(len(cg))],dtype=object)}
	lut_out={}
	cnt=1
	for k in cg:
		lut_out[k]=cnt
		d['Value'][cnt-1]=cnt
		d['Name'][cnt-1]=k
		cnt=cnt+1

	zCEC10c=copy.deepcopy(zRef)
	zCEC20c=copy.deepcopy(zRef)
	zCEC10c['Data']=(len(cg)+1)*np.ones(zRef['Data'].shape,dtype='uint8')
	zCEC20c['Data']=(len(cg)+1)*np.ones(zRef['Data'].shape,dtype='uint8')

	# Forest
	ind=np.where( (zCEC10['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
				 (zCEC10['Data']==lut_in['Sub-polar taiga needleleaf forest']) |
				 (zCEC10['Data']==lut_in['Tropical or sub-tropical broadleaf evergreen forest']) |
				 (zCEC10['Data']==lut_in['Tropical or sub-tropical broadleaf deciduous forest']) |
				 (zCEC10['Data']==lut_in['Temperate or sub-polar broadleaf deciduous forest']) |
				 (zCEC10['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
				 (zCEC10['Data']==lut_in['Mixed Forest']) )
	zCEC10c['Data'][ind]=lut_out['Forest']
	ind=np.where( (zCEC20['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
				 (zCEC20['Data']==lut_in['Sub-polar taiga needleleaf forest']) |
				 (zCEC20['Data']==lut_in['Tropical or sub-tropical broadleaf evergreen forest']) |
				 (zCEC20['Data']==lut_in['Tropical or sub-tropical broadleaf deciduous forest']) |
				 (zCEC20['Data']==lut_in['Temperate or sub-polar broadleaf deciduous forest']) |
				 (zCEC20['Data']==lut_in['Temperate or sub-polar needleleaf forest']) |
				 (zCEC20['Data']==lut_in['Mixed Forest']) )
	zCEC20c['Data'][ind]=lut_out['Forest']

	# Shrubland
	ind=np.where( (zCEC10['Data']==lut_in['Tropical or sub-tropical shrubland']) |
				 (zCEC10['Data']==lut_in['Temperate or sub-polar shrubland']) )
	zCEC10c['Data'][ind]=lut_out['Shrubland']
	ind=np.where( (zCEC20['Data']==lut_in['Tropical or sub-tropical shrubland']) |
				 (zCEC20['Data']==lut_in['Temperate or sub-polar shrubland']) )
	zCEC20c['Data'][ind]=lut_out['Shrubland']

	# Grassland
	ind=np.where( (zCEC10['Data']==lut_in['Tropical or sub-tropical grassland']) |
				 (zCEC10['Data']==lut_in['Temperate or sub-polar grassland']) )
	zCEC10c['Data'][ind]=lut_out['Grassland']
	ind=np.where( (zCEC20['Data']==lut_in['Tropical or sub-tropical grassland']) |
				 (zCEC20['Data']==lut_in['Temperate or sub-polar grassland']) )
	zCEC20c['Data'][ind]=lut_out['Grassland']

	# Liche-moss
	ind=np.where( (zCEC10['Data']==lut_in['Sub-polar or polar shrubland-lichen-moss']) |
				 (zCEC10['Data']==lut_in['Sub-polar or polar grassland-lichen-moss']) |
				 (zCEC10['Data']==lut_in['Sub-polar or polar barren-lichen-moss']) )
	zCEC10c['Data'][ind]=lut_out['Lichen-Moss']
	ind=np.where( (zCEC20['Data']==lut_in['Sub-polar or polar shrubland-lichen-moss']) |
				 (zCEC20['Data']==lut_in['Sub-polar or polar grassland-lichen-moss']) |
				 (zCEC20['Data']==lut_in['Sub-polar or polar barren-lichen-moss']) )
	zCEC20c['Data'][ind]=lut_out['Lichen-Moss']

	# Wetland
	ind=np.where( (zCEC10['Data']==lut_in['Wetland']) )
	zCEC10c['Data'][ind]=lut_out['Wetland']
	ind=np.where( (zCEC20['Data']==lut_in['Wetland']) )
	zCEC20c['Data'][ind]=lut_out['Wetland']
	cnt=cnt+1

	# Cropland
	ind=np.where( (zCEC10['Data']==lut_in['Cropland']) )
	zCEC10c['Data'][ind]=lut_out['Cropland']
	ind=np.where( (zCEC20['Data']==lut_in['Cropland']) )
	zCEC20c['Data'][ind]=lut_out['Cropland']

	# Barren Land
	ind=np.where( (zCEC10['Data']==lut_in['Barren lands']) )
	zCEC10c['Data'][ind]=lut_out['Barren Ground']
	ind=np.where( (zCEC20['Data']==lut_in['Barren lands']) )
	zCEC20c['Data'][ind]=lut_out['Barren Ground']

	# Urban
	ind=np.where( (zCEC10['Data']==lut_in['Urban']) )
	zCEC10c['Data'][ind]=lut_out['Urban']
	ind=np.where( (zCEC20['Data']==lut_in['Urban']) )
	zCEC20c['Data'][ind]=lut_out['Urban']

	# Water
	ind=np.where( (zCEC10['Data']==lut_in['Water']) )
	zCEC10c['Data'][ind]=lut_out['Water']
	ind=np.where( (zCEC20['Data']==lut_in['Water']) )
	zCEC20c['Data'][ind]=lut_out['Water']
	cnt=cnt+1

	# Snow and Ice
	ind=np.where( (zCEC10['Data']==lut_in['Snow and Ice']) )
	zCEC10c['Data'][ind]=lut_out['Snow and Ice']
	ind=np.where( (zCEC20['Data']==lut_in['Snow and Ice']) )
	zCEC20c['Data'][ind]=lut_out['Snow and Ice']
	cnt=cnt+1

	return lut_out,zCEC10c,zCEC20c

#%%
def PrepareLandUseCEC(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	fin=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2010.tif'
	fout=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2010_bc1ha.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	zLCC_CEC10=gis.OpenGeoTiff(fout)
	zLCC_CEC10=gis.ClipToRaster(zLCC_CEC10,zRef)	
	ind=np.where(zRef['Data']==0)
	zLCC_CEC10['Data'][ind]=0 
	fin=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2020.tif'
	fout=r'C:\Data\Land Cover\NALCMS\NALCMS_lc_2020_bc1ha.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])	
	zLCC_CEC20=gis.OpenGeoTiff(fout)
	zLCC_CEC20=gis.ClipToRaster(zLCC_CEC20,zRef)	
	ind=np.where(zRef['Data']==0)
	zLCC_CEC20['Data'][ind]=0		 
	# Compress categories to remove tropics	
	lut,zLCC_CEC10,zLCC_CEC20=NALCMS_Compress(meta['LUT']['Derived']['lcc_cec'],zRef,zLCC_CEC10,zLCC_CEC20)
	# Manually saved LUT to excell
	gis.SaveGeoTiff(zLCC_CEC10,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2010_Compressed.tif')
	gis.SaveGeoTiff(zLCC_CEC20,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_CEC_2020_Compressed.tif')
	return
