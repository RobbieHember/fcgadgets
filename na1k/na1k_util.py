#%% Import modules
import os
import numpy as np
import gc
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
from scipy.interpolate import griddata
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.gaia.gaia_util as gaia

#%% Initialize project
def Init(*argv):
	meta={}
	# Set paths
	meta['Paths']={}
	#meta['Paths']['na1k']=r'C:\Data\na1k'
	meta['Paths']['na1k']=r'D:\Data\na1k'
	meta['Paths']['na1k Ref Grid']=meta['Paths']['na1k'] + '\\refgrid.tif'
	meta['Paths']['na1k Ref Grid Buf']=meta['Paths']['na1k'] + '\\refgrid_buf.tif'
	meta['Paths']['na1k Ref Grid Buf Rev']=meta['Paths']['na1k'] + '\\refgrid_buf_rev.tif'
	meta['Paths']['na1k DEM']=r'C:\Data\DEM\North America\namer_dem1.tif'
	meta['Paths']['ECCD']=r'C:\Data\Climate\ECCD'
	meta['Paths']['GSOD']=r'C:\Data\Climate\GSOD'

	meta['Paths']['DB']={}
	meta['Paths']['DB']['Harvest']=r'C:\Data\Harvest'
	meta['Paths']['DB']['Waste']=r'C:\Data\Waste Wood\FromWasteSystem'
	meta['Paths']['DB']['HBS']=r'C:\Data\Harvest\HBS'
	meta['Paths']['DB']['CO2']=r'C:\Data\Carbon Dioxide Concentration\Carbon Dioxide Concentration Mid Year.xlsx'
	meta['Paths']['DB']['NDEP']=r'C:\Data\Nitrogen Deposition'
	meta['Paths']['DB']['CMIP6']=r'C:\Data\Climate\CMIP6\Monthly'
	meta['Paths']['DB']['20thCR']=r'C:\Data\Climate\Reanlysis\20th Century'
	meta['Paths']['DB']['CRU']=r'C:\Data\Climate\CRU'

	meta['Graphics']={}
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\na1k'
	meta['Graphics']['Map']={}
	meta['Graphics']['Map']['Fig Width']=9.75
	meta['Graphics']['Map']['Side Space']=0
	meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
	meta['Graphics']['Map']['Map Axis Vis']='off'
	meta['Graphics']['Map']['Map Grid Vis']=False
	meta['Graphics']['Map']['Legend X']=0.72
	meta['Graphics']['Map']['Legend Width']=0.03
	meta['Graphics']['Map']['Legend Font Size']=7
	meta['Graphics']['Map']['Legend Text Space']=0.035
	meta['Graphics']['Map']['Show Bound Land Mask']='On'
	meta['Graphics']['Map']['Show Bound Within']='Off'
	meta['Graphics']['Map']['Show Lakes']='Off'
	meta['Graphics']['Map']['Show Rivers']='Off'
	meta['Graphics']['Map']['Show Roads']='Off'
	meta['Graphics']['Map']['Show TPF and Cities']='Off'
	meta['Graphics']['Map']['Show Symbol Labels']='Off'

	meta['Graphics']['gp']=gu.SetGraphics('Manuscript')
	meta['Graphics']['gp']['AxesLetterStyle']='Caps'
	meta['Graphics']['gp']['AxesLetterFontWeight']='Bold'

	# Initiate geospatial info
	meta['Geos']={}

	# Import coordinate reference system
	pb=ImportPoliticalBoundaries(meta)
	meta['Geos']['crs']=pb.crs

	# Import LUTs
	if 'LUT' not in meta:
		meta['LUT']={}
		meta['LUT']['EZ']={}
		d=gu.ReadExcel(meta['Paths']['na1k'] + '\\LUTs\\LUT_Ecozones_CEC_L1.xlsx',sheet_name='Sheet1',skiprows=0)
		for i in range(d['ID'].size):
			meta['LUT']['EZ'][d['Name'][i]]=d['ID'][i]
		meta['LUT']['PB']={}
		d=gu.ReadExcel(meta['Paths']['na1k'] + '\\LUTs\\LUT_PoliticalBoundaries.xlsx',sheet_name='Sheet1',skiprows=0)
		for i in range(d['ID'].size):
			meta['LUT']['PB'][d['Name'][i]]=d['ID'][i]

	# Time vectors
	meta['tvd']=gu.tvec('d',1950,2025)
	meta['tvm']=gu.tvec('m',1950,2025)
	meta['tva']=np.unique(meta['tvd'][:,0])

	# Scale factors
	meta['Climate']={}
	meta['Climate']['SF']={}
	meta['Climate']['SF']['cwd']=0.1 # Stored as mm/month
	meta['Climate']['SF']['cmi']=0.1 # Stored as mm/month
	meta['Climate']['SF']['ea']=0.01
	meta['Climate']['SF']['ea_fromtmin']=0.01
	meta['Climate']['SF']['ea_biasadj']=0.01
	meta['Climate']['SF']['es']=0.01
	meta['Climate']['SF']['eta']=0.1 # Stored as mm/month
	meta['Climate']['SF']['etp']=0.1 # Stored as mm/month
	meta['Climate']['SF']['prcp']=0.1 # Stored as mm/month
	meta['Climate']['SF']['rswd']=0.01
	meta['Climate']['SF']['tmean']=0.01
	meta['Climate']['SF']['tmin']=0.01
	meta['Climate']['SF']['tmax']=0.01
	meta['Climate']['SF']['vpd']=0.01
	meta['Climate']['SF']['ws']=0.1 # Stored as mm/month
	meta['Climate']['SF']['wsp']=0.1 # Stored as mm/month
	meta['Climate']['SF']['runoff']=0.1 # Stored as mm/month
	meta['Climate']['SF']['melt']=0.1 # Stored as mm/month

	meta['Missing Number']=-99

	# Parameters from Gaia
	meta['Parameters']=gaia.HydroMetCon()

	# Import political boundaries
	meta['pb']=ImportPoliticalBoundaries(meta)

	return meta

#%%
def ImportPoliticalBoundaries(meta):
	srs=gis.ImportSRSs()
	pb=gpd.read_file(r'C:\Data\Geodatabases\North America\bound_p.shp').to_crs(srs['String']['NACID'])
	u=pb['STATEABB'].unique()[1:]
	pb['ID']=np.zeros(len(pb['STATEABB']),dtype='int16')
	for iU in range(u.size):
		pb['ID'][(pb['STATEABB']==u[iU])]=iU+1
	return pb

#%%
def CreateReferenceGrid(meta):
	z=np.random.randint(5,size=(6500,6400),dtype='int16')
	with rasterio.open(
		meta['Paths']['na1k Ref Grid'],
		mode="w",
		driver="GTiff",
		height=z.shape[0],
		width=z.shape[1],
		count=1,
		dtype=z.dtype,
		crs=meta['pb'].crs,
		transform=from_origin(-3.05e6,4.96e6,1000,1000), # Inputs (west,north,xsize,ysize)
		) as new_dataset:
		new_dataset.write(z,1)
	
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	shapes=((geom,value) for geom, value in zip(meta['pb'].geometry,meta['pb']['ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	
	# Mask
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where(burned>0)
	z1['Data'][ind]=1
	plt.close('all');plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['na1k Ref Grid'])
	
	# Buffer
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	kernel=np.ones((5,5))
	mask=cv2.dilate(zRef['Data'].astype(np.uint8),kernel,iterations=1)
	z1=copy.deepcopy(zRef)
	z1['Data']=mask.astype('int8')
	plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['na1k Ref Grid Buf'])

	kernel=np.ones((2,2))
	mask=cv2.erode(zRef['Data'].astype(np.uint8),kernel,iterations=1)
	z1=copy.deepcopy(zRef)
	z1['Data']=mask.astype('int8')
	plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['na1k Ref Grid Buf Rev'])
	return

#%%
def RasterizeEcozonesCEC(meta):
	srs=gis.ImportSRSs()
	df=gpd.read_file(r'C:\Data\Ecozones\CEC North America\NA_Terrestrial_Ecoregions_v2_level1.shp').to_crs(srs['String']['NACID'])
	df['LEVEL1']=pd.to_numeric(df['LEVEL1'])
	
	df2=pd.DataFrame({'ID':df['LEVEL1'].values,'Name':df['NameL1_En']})
	df2=df2.drop_duplicates()
	df2=df2.sort_values(by=['ID'])
	df2.to_excel(meta['Paths']['na1k'] + '\\LUTs\\LUT_Ecozones_CEC_L1.xlsx',index=False)
	
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	shapes=((geom,value) for geom, value in zip(df.geometry,df['LEVEL1']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where(burned>0)
	z1['Data'][ind]=burned[ind]
	plt.close('all');plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Ecozones_CEC_L1.tif')
	return

#%%
def RasterizePoliticalBoundaries(meta):

	srs=gis.ImportSRSs()
	pb=gpd.read_file(r'C:\Data\Geodatabases\North America\bound_p.shp').to_crs(srs['String']['NACID'])
	u=pb['STATEABB'].unique()[1:]
	pb['ID']=np.zeros(len(pb['STATEABB']),dtype='int16')
	for iU in range(u.size):
		pb['ID'][(pb['STATEABB']==u[iU])]=iU+1
	
	df2=pd.DataFrame({'ID':pb['ID'],'Name':pb['STATEABB']})
	df2=df2.drop_duplicates()
	df2=df2.sort_values(by=['ID'])
	df2=df2[df2.ID!=0]
	df2.to_excel(meta['Paths']['na1k'] + '\\LUTs\\LUT_PoliticalBoundaries.xlsx',index=False)
	
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	shapes=((geom,value) for geom, value in zip(pb.geometry,pb['ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind=np.where(burned>0)
	z1['Data'][ind]=burned[ind]
	plt.close('all');plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\PoliticalBoundaries.tif')
	return

#%% Process DEM
def ProcessDEM(meta):
	# DEM from DataBasin. Comes in Arc Layer file so convert to geotiff in Arc.
	#dem=gis.OpenGeoTiff(z1,r'C:\Users\rhember\Documents\Data\DEM\North America\dem.tif')
	
	fin=r'C:\Users\rhember\Documents\Data\DEM\North America\namer_dem1.tif'
	fout=meta['Paths']['na1k'] + '\\dem.tif'
	fref=meta['Paths']['na1k Ref Grid']
	gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])
	
	zE=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\dem.tif')
	zE['Data'][zE['Data']>30000]=-1
	zE['Data']=zE['Data'].astype('int16')
	#plt.close('all');plt.matshow(zE['Data'],clim=[0,3000])
	gis.SaveGeoTiff(zE,meta['Paths']['na1k'] + '\\dem.tif')
	return

#%%
def ImportNormalsFromClimateNA(meta):
	# Instructions for running ClimateNA:
	# Convert DEM to ASCII format (.asc) in ArcGIS. (Trying DIY has never worked easily!)
	# Run ClimateNA with the .asc file in geographic coordinates
	# Use this script to import and reproject the monthly climatologies.
	# Notes: 2024 - temperatures not given to any decimal place, limited precision
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	zBuf=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid Buf'])
	zBufR=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid Buf Rev'])
	zE=gis.OpenGeoTiff(meta['Paths']['na1k DEM'])
	iOut=np.where(zE['Data']>6000)
	#pthin=r'C:\Data\Climate\ClimateNA\NA1k'
	pthin=r'D:\Backup Data\Climate\ClimateNA\NA1k'
	vL=['Tmin','Tmax','PPT']
	vL2=['tmin','tmax','prcp']
	#vL=['PPT'];	vL2=['prcp']
	for iV in range(len(vL)):
		for mo in range(1,13):
			if mo<=9:
				mot='0' + str(mo)
			else:
				mot=str(mo)

			# Convert to tif and save in geographic coordinates
			flg=0
			if flg==1:
				# Convert to geotiff
				src=rasterio.open(pthin + '\\' + vL[iV] + mot + '.asc')
				z=copy.deepcopy(zE)
				z['Data']=src.read(1).astype('float')
				src.close()
				z['Data']=z['Data']/meta['SF'][vL2[iV]]
				z['Data']=z['Data'].astype('int16')
				gis.SaveGeoTiff(z,pthin + '\\tmp2.tif')
				# Reproject and change extent
				fin=pthin + '\\tmp2.tif'
				fout=meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + vL2[iV] + '_norm_1971to2000_' + str(mo) + '.tif'
				fref=meta['Paths']['na1k Ref Grid']
				gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])
			else:
				# Received as geotiffs from TW
				z0=gis.OpenGeoTiff(pthin + '\\' + vL[iV] + mot + '.tif')
				z=copy.deepcopy(zE)
				z['Data']=z0['Data']
				z['Data'][iOut]=meta['Missing Number']
				z['Data']=z['Data']/meta['SF'][vL2[iV]]
				z['Data']=z['Data'].astype('int16')
				gis.SaveGeoTiff(z,pthin + '\\tmp2.tif')

				# Reproject and change extent
				fin=pthin + '\\tmp2.tif'
				fout=meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + vL2[iV] + '_norm_1971to2000_' + str(mo) + '.tif'
				fref=meta['Paths']['na1k Ref Grid']
				gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])

			# Import reprojected map
			z0=gis.OpenGeoTiff(fout)['Data'].astype('float')*meta['SF'][vL2[iV]]

			# Remove exterior
			z0[zRef['Data']==0]=meta['Missing Number'] #plt.matshow(z0)

			# Find missing values within mask
			iMissing=np.where( (zRef['Data']==1) & (z0==meta['Missing Number']) )
			#print(iMissing[0].size)

			# Gap-fill edges
			MaskGap=-1*(zRef['Data'].astype(float)-zBuf['Data'].astype(float)); plt.close('all'); #plt.matshow(MaskGap); plt.colorbar()
			indGap=np.where( (MaskGap==1) | (zRef['Data']==1) & (z0==meta['Missing Number']) )
			MaskCal=zRef['Data']-zBufR['Data'] #plt.matshow(MaskCal)
			indCal=np.where( (MaskCal==1) & (z0!=meta['Missing Number']) )
			#a=MaskGap.copy();a[MaskCal==1]=2;plt.close('all'); plt.matshow(a)
			x=zRef['X'][indCal]
			y=zRef['Y'][indCal]
			z=z0[indCal]
			interp=NearestNDInterpolator(list(zip(x,y)),z)
			iz=interp(zRef['X'][indGap],zRef['Y'][indGap])
			z0[indGap]=iz # plt.close('all');plt.matshow(z0)

			# Save
			z=copy.deepcopy(zRef)
			z['Data']=z0
			#z['Data'][iOut]=meta['Missing Number']
			z['Data']=z['Data']/meta['SF'][vL2[iV]]
			z['Data']=z['Data'].astype('int16')
			gis.SaveGeoTiff(z,fout)

	for mo in range(1,13):
		zMin=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmin_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['tmin']
		zMax=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmax_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['tmax']
		z0=(zMin+zMax)/2
		z=copy.deepcopy(zRef)
		z['Data']=z0
		#z['Data'][iOut]=meta['Missing Number']
		z['Data']=z['Data']/meta['SF']['tmean']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcSolarRadiationNormals(meta):
	# The Matlab code has not been translated into Python.
	# For the time being, continue using the old normals, reprojected and clipped
	# to the new na1k spatial reference system.
	for mo in range(12):
		fin=r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif'
		fout=meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_rswd_norm_1971to2000_' + str(mo+1) + '.tif'
		gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['na1k Ref Grid'],meta['Geos']['crs'])

		zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
		z1=gis.OpenGeoTiff(fout)
		z1['Data']=z1['Data'].astype('float')
		z1['Data']=z1['Data']/10
		z1['Data']=z1['Data']/meta['SF']['rswd']
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.close('all'); plt.matshow(z1['Data'],clim=[0,3000])
		gis.SaveGeoTiff(z1,fout)
	return

#%%
def CalcSaturationVapourPressureNormal(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	for mo in range(1,13):
		zT=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['tmean']
		z=copy.deepcopy(zRef)
		z['Data']=gaia.GetEstar(zT)
		z['Data'][zRef['Data']==0]=meta['Missing Number']
		z['Data']=z['Data']/meta['SF']['es']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_es_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcActualVapourPressureNormalFromTemps(meta):
	# Reference: https://www.agraria.unirc.it/documentazione/materiale_didattico/1462_2016_412_24509.pdf
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	for mo in range(1,13):
		zT=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmin_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['tmean']
		z=copy.deepcopy(zRef)
		z['Data']=10*(0.6108*np.exp((17.27*zT)/(zT+237.3)))
		z['Data'][zRef['Data']==0]=meta['Missing Number']
		z['Data']=z['Data']/meta['SF']['ea']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcActualVapourPressureNormalFromCRU(meta):
	# Given in numpy mask arrays so you need to use the filled method
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	tv=gu.tvec('m',1901,2023)
	srs=gis.ImportSRSs()

	# Import netcdf data
	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.vap.dat.nc.gz'
	ea0={}
	with gzip.open(fin) as gz:
		with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
			for k in ds.variables.keys():
				ea0[k]=ds.variables[k][:]
	ea0['vap']=ea0['vap'].filled()
	ind=np.where(ea0['vap']>999)
	ea0['vap'][ind]=meta['Missing Number']
	#ind=np.where(ea0['vap']>-999)
	#np.mean(ea0['vap'][ind])

	# Isolate BC area  (so that interpolation is faster)
	lon,lat=np.meshgrid(ea0['lon'].filled(),ea0['lat'].filled(),sparse=False)
	x1,y1=srs['Proj']['NACID'](lon,lat)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))

	for mo in range(1,13):
		iT=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo) )[0]
		z0=np.mean(ea0['vap'][iT,:,:],axis=0)[indS]
		z1=copy.deepcopy(zRef)
		z1['Data']=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')
		z1['Data'][zRef['Data']==0]=meta['Missing Number']
		z1['Data']=z1['Data']/meta['SF']['ea']
		z1['Data']=z1['Data'].astype('int16')
		#plt.close('all'); plt.matshow(z1['Data'],clim=[0,3000])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_fromcru_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcActualVapourPressureNormalBiasCorrected(meta):
	# Import GSOD
	dS=gu.ipickle(meta['Paths']['GSOD'] + '\\Processed\\gsod.pkl')
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	iGrdS=gis.GetGridIndexToPoints(zRef,dS['X'],dS['Y'])
	iOut=np.where(zRef['Data'][iGrdS]==0)[0]
	dS['Data']['ea'][:,iOut]=np.nan
	srs=gis.ImportSRSs()

	# Import CRU
	tvC=gu.tvec('m',1901,2023)
	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.vap.dat.nc.gz'
	vap0={}
	with gzip.open(fin) as gz:
		with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
			for k in ds.variables.keys():
				vap0[k]=ds.variables[k][:]
	ea_cru=vap0['vap'].filled()
	ea_cru[ea_cru>999]=meta['Missing Number']
	#ind=np.where(ea_cru>meta['Missing Number']); print(np.mean(ea_cru[ind]))
	for mo in range(12):
		#itmu=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo) )
		mu=np.mean(ea_cru[mo::12,:,:],axis=0)
		ea_cru[mo::12,:,:]=ea_cru[mo::12,:,:]-np.tile(mu,(int(ea_cru.shape[0]/12),1,1))
	lonC,latC=np.meshgrid(vap0['lon'].filled(),vap0['lat'].filled(),sparse=False)
	xC,yC=srs['Proj']['NACID'](lonC,latC)

	for mo in range(1,13):
		# Find stations with enough data
		#iT=np.where( (meta['tvm'][:,0]>=1971) & (meta['tvm'][:,0]<=2019) & (meta['tvm'][:,1]==mo) )[0]
		iT=np.where( (meta['tvm'][:,0]>=1971) & (meta['tvm'][:,0]<=2000) & (meta['tvm'][:,1]==mo) )[0]
		N=np.sum(~np.isnan(dS['Data']['ea'][iT,:]),axis=0)
		iS=np.where(N>15)[0]
		# plt.plot(dS['X'][iS],dS['Y'][iS],'k.')

		ea_mu0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['ea']
		ea_mu0[zRef['Data']==0]=np.nan
		iGrd=gis.GetGridIndexToPoints(zRef,dS['X'][iS],dS['Y'][iS])
		ea_mu0=ea_mu0[iGrd]

		iGrdC=gis.GetGridIndexToPoints({'X':xC,'Y':yC},dS['X'][iS],dS['Y'][iS])
		ea_hat=np.nan*np.ones((meta['tvm'].shape[0],iS.size))
		for yr in range(1970,2020):
			iT1=np.where( (tvC[:,0]==yr) & (tvC[:,1]==mo) )[0]
			iT2=np.where( (meta['tvm'][:,0]==yr) & (meta['tvm'][:,1]==mo) )[0]
			ea_cruT=np.squeeze(ea_cru[iT1,:,:])
			#ea_hat[iT2,:]=ea_mu0
			ea_hat[iT2,:]=ea_mu0+ea_cruT[iGrdC]

		# Remove missing data, remove large outliers and reduced negative outliers (that affect mapping)
		x=dS['Data']['ea'][iT,:][:,iS]
		y=ea_hat[iT,:]
		E_th=12
		ind=np.where( (np.isnan(x)==True) | (np.isnan(y)==True) | (np.abs(x-y)>E_th) )
		x[ind]=np.nan
		y[ind]=np.nan
		E=np.nanmean(y-x,axis=0)
		ind=np.where(E<-5)[0]
		E[ind]=-2

		# Plot map of error at stations
		bw=1; bin=np.arange(0,12,bw)
		ms=np.linspace(0.5,12,bin.size)
		plt.close('all');fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,14))
		for i in range(bin.size):
			ind=np.where( (E>0) & (np.abs(E-bin[i])<=bw/2) )[0]
			plt.plot(dS['X'][iS][ind],dS['Y'][iS][ind],'ro',ms=ms[i],mew=0.25,mfc='none')
			ind=np.where( (E<0) & (np.abs(E--bin[i])<=bw/2) )[0]
			plt.plot(dS['X'][iS][ind],dS['Y'][iS][ind],'bo',ms=ms[i],mew=0.25,mfc='none')

		# Extract
		ikp=np.where( (np.isnan(E)==False) )[0]
		e0=E[ikp]
		x0=dS['X'][iS][ikp]
		y0=dS['Y'][iS][ikp]

		# Add water areas
		iWat=np.where(zRef['Data'][0::150,0::150]==0)
		x0=np.append(x0,zRef['X'][0::150,0::150][iWat])
		y0=np.append(y0,zRef['Y'][0::150,0::150][iWat])
		e0=np.append(e0,np.zeros(iWat[0].size))

		# Linear interpolation map of error
		flg=1
		if flg==1:
			ivl=5
			iz=griddata(np.column_stack((x0,y0)),e0,(zRef['X'][0::ivl,0::ivl],zRef['Y'][0::ivl,0::ivl]),method='linear')
			iz=scipy.ndimage.zoom(iz,ivl,order=0)
			iz[zRef['Data']==0]=0
			plt.close('all'); plt.matshow(iz,clim=[-12,12])
			CorrFac=np.tile(iz[iGrd],(iT.size,1))

		else:
			# Polynomial surface interpolation
			beta=np.array([0.000038,0.75,0.025])
	
			mask=np.zeros(zRef['Data'].shape,dtype='int8')
			ind=np.where(zRef['Data']>0)
			mask[ind]=1
			kernel=np.ones((50,50))
			mask=cv2.dilate(mask.astype(np.uint8),kernel,iterations=1)
			
			itvl=50
			type='Hyperbolic'
			xI=zRef['X'][0::itvl,0::itvl]
			yI=zRef['Y'][0::itvl,0::itvl]
			maskI=mask[0::itvl,0::itvl]
			
			iz,iE,iN=gu.PolynomialSurfaceFit(x0/1000,y0/1000,e0,xI/1000,yI/1000,maskI,beta)
			plt.close('all'); plt.matshow(iz,clim=[-12,12])
			iz=scipy.ndimage.zoom(iz,itvl,order=0)
			CorrFac=np.tile(iz[iGrd],(iT.size,1))

		# Skill without correction
		x=dS['Data']['ea'][iT,:][:,iS].flatten()
		y=ea_hat[iT,:].flatten()
		ikp=np.where( (x>0) & (y>0) & (np.abs(x-y)<20) )[0]
		x=x[ikp]
		y=y[ikp]
		rs,txt=gu.GetRegStats(x,y)
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
		ax.plot(x,y,'o',ms=3,mec='w',mfc='k',mew=0.5)
		ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
		ax.text(33,2,rs['txt'],fontsize=7,ha='right')
		ax.text(33,33,'1:1',fontsize=7,ha='center')
		ax.set(xlabel='GSOD measurements (hPa)',ylabel='NA1K grid (hPa)',xticks=np.arange(0,40,5),yticks=np.arange(0,40,5),xlim=[0,35],ylim=[0,35])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		plt.tight_layout()
		#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Name','png',900)

		# Skill with correction
		x=dS['Data']['ea'][iT,:][:,iS].flatten()
		ea_hatCF=ea_hat[iT,:]-CorrFac
		y=ea_hatCF.flatten()
		ikp=np.where( (x>0) & (y>0) & (np.abs(x-y)<20) )[0]
		x=x[ikp]
		y=y[ikp]
		rs,txt=gu.GetRegStats(x,y)
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
		ax.plot(x,y,'o',ms=3,mec='w',mfc='k',mew=0.5)
		ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
		ax.text(33,2,rs['txt'],fontsize=7,ha='right')
		ax.text(33,33,'1:1',fontsize=7,ha='center')
		ax.set(xlabel='GSOD measurements (hPa)',ylabel='NA1K grid (hPa)',xticks=np.arange(0,40,5),yticks=np.arange(0,40,5),xlim=[0,35],ylim=[0,35])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		plt.tight_layout()

		# Save new actual vapour pressure
		z1=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')
		tmp=np.maximum(0,(z1['Data'].astype('float')*meta['SF']['ea'])-iz)
		tmp=tmp/meta['SF']['ea']
		z1['Data']=tmp.astype('int16')
		# plt.matshow(z1['Data'],clim=[0,3000])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_biasadj_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def QA_VapourPressureCRU_Vs_BiasAdjTempMod(meta):
	ivl=50
	z1=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_ea_biasadj_norm_1971to2000_mjjas.tif')['Data'][0::ivl,0::ivl].astype('float')*meta['SF']['ea']
	z2=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_ea_fromcru_norm_1971to2000_mjjas.tif')['Data'][0::ivl,0::ivl].astype('float')*meta['SF']['ea']
	iMask=np.where(zRef['Data'][0::ivl,0::ivl]==1)
	x=z2[iMask]
	y=z1[iMask]
	ikp=np.where( (x>0.5) & (y>0.5) & (np.abs(x-y)<30) )[0]
	x=x[ikp]
	y=y[ikp]
	rs,txt=gu.GetRegStats(x,y)
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
	ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot(x,y,'o',ms=3,mec='w',mfc='k',mew=0.5)
	ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
	ax.text(33,2,rs['txt'],fontsize=7,ha='right')
	ax.text(33,33,'1:1',fontsize=7,ha='center')
	ax.set(xlabel='CRU (hPa)',ylabel='Temp-based model (hPa)',xticks=np.arange(0,40,5),yticks=np.arange(0,40,5),xlim=[0,35],ylim=[0,35])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	return

#%%
def CalcVapourPressureDeficitNormal(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	for mo in range(1,13):
		zEa=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_biasadj_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['ea']
		zEs=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_es_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['es']
		vpd=np.maximum(0,zEs-zEa)

		z=copy.deepcopy(zRef)
		z['Data']=vpd
		z['Data'][zRef['Data']==0]=meta['Missing Number']
		z['Data']=z['Data']/meta['SF']['vpd']
		z['Data']=z['Data'].astype('int16')
		# plt.matshow(z['Data'],clim=[0,2000])
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_vpd_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
# def CalcVapourPressureDeficitNormalFromCRU(meta):
# 	# Given in numpy mask arrays so you need to use the filled method
# 	srs=gis.ImportSRSs()
# 	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
# 	tv=gu.tvec('m',1901,2023)

# 	# Import netcdf data
# 	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.vap.dat.nc.gz'
# 	vap0={}
# 	with gzip.open(fin) as gz:
# 		with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
# 			for k in ds.variables.keys():
# 				vap0[k]=ds.variables[k][:]
# 	vap0['vap']=vap0['vap'].filled()
# 	ind=np.where(vap0['vap']>999)
# 	vap0['vap'][ind]=meta['Missing Number']
# 	ind=np.where(vap0['vap']>-999)
# 	np.mean(vap0['vap'][ind])

# 	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.tmp.dat.nc.gz'
# 	tmp0={}
# 	with gzip.open(fin) as gz:
# 		with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
# 			for k in ds.variables.keys():
# 				if k!='tmp':
# 					continue
# 				tmp0[k]=ds.variables[k][:]
# 				#print(np.mean(ds.variables[k][:]))
# 	tmp0=tmp0['tmp'].filled()
# 	ind=np.where(tmp0>999)
# 	tmp0[ind]=meta['Missing Number']
# 	esat0=gaia.GetEstar(tmp0)
# 	ind=np.where(esat0>999999)
# 	esat0[ind]=meta['Missing Number']

# 	vpd0=np.maximum(0,esat0-vap0['vap'])
# 	ind=np.where(vpd0>999)
# 	vpd0[ind]=meta['Missing Number']

# 	# Isolate BC area  (so that interpolation is faster)
# 	lon,lat=np.meshgrid(vap0['lon'].filled(),vap0['lat'].filled(),sparse=False)
# 	x1,y1=srs['Proj']['NACID'](lon,lat)
# 	bw=100000
# 	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
# 	xy0=np.column_stack((x1[indS],y1[indS]))

# 	for mo in range(12):
# 		iT=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo) )[0]
# 		z0=np.mean(vpd0[iT,:,:],axis=0)[indS]
# 		z1=copy.deepcopy(zRef)
# 		z1['Data']=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')
# 		ind=np.where(zRef['Data']==0)
# 		z1['Data'][ind]=0
# 		z1['Data']=z1['Data']/meta['SF']['vpd']
# 		z1['Data']=z1['Data'].astype('int16')
# 		#plt.close('all'); plt.matshow(z1['Data'],clim=[0,3000])
# 		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_vpd_norm_1971to2000_' + str(mo+1) + '.tif')
# 	return

#%% Download daily ECCD data
# Station inventory: https://collaboration.cmc.ec.gc.ca/cmc/climate/Get_More_Data_Plus_de_donnees/
# timeframe = 1: for hourly data, 2: for daily data, 3 for monthly data 
def DownloadEnvironmentCanadaClimateData(meta):

	# Starting from scratch or just updating?
	#mode='Start from Scratch'
	mode='Update Single Year'
	YearUpdate=2023
	
	# List of stations
	dSL=gu.ReadExcel(meta['Paths']['ECCD'] + '\\Daily\\Station Inventory EN.xlsx')
	
	# Fix climate ID in dictionary -> excel messes it up with scientific notation
	dfSL=pd.read_csv(meta['Paths']['ECCD'] + '\\Daily\\Station Inventory EN.csv',skiprows=3)
	dSL['Climate ID']=dfSL['Climate ID'].values
	
	# List of variables collected by Environment Canada
# 	vL_All=['Data Quality','Max Temp (°C)','Max Temp Flag','Min Temp (°C)','Min Temp Flag','Mean Temp (°C)',
# 	'Mean Temp Flag','Heat Deg Days (°C)','Heat Deg Days Flag','Cool Deg Days (°C)','Cool Deg Days Flag','Total Rain (mm)',
# 	'Total Rain Flag','Total Snow (cm)','Total Snow Flag','Total Precip (mm)','Total Precip Flag','Snow on Grnd (cm)',
# 	'Snow on Grnd Flag','Dir of Max Gust (10s deg)','Dir of Max Gust Flag','Spd of Max Gust (km/h)','Spd of Max Gust Flag']
	
	vL=['Max Temp (°C)','Min Temp (°C)','Mean Temp (°C)','Mean Temp Flag','Total Precip (mm)','Total Precip Flag']
	uS=np.unique(dSL['Station ID'])
	ListMissing=[]
	for iS in range(uS.size):
		indS=np.where(dSL['Station ID']==uS[iS])[0]
		fout=meta['Paths']['ECCD'] + '\\Daily\Raw\\Combined\\' + str(dSL['Climate ID'][indS[0]]) + '.pkl'
		tv=gu.tvec('d',dSL['First Year'][indS],dSL['Last Year'][indS])
		uY=np.unique(tv[:,0])
	
		if mode=='Start from Scratch':
			# Initialize combined station file if does not exist
			if os.path.isfile(fout)==False:
				d={}
				for v in vL:
					if v[-4:]=='Flag':
						d[v]=np.array(['' for _ in range(tv.shape[0])],dtype=object)
					else:
						d[v]=-999*np.ones(tv.shape[0])
			else:
				print('Already downloaded, move on...')
				continue
		else:
			if YearUpdate in uY:
				uY=[YearUpdate]
			else:
				continue
				d=gu.ipickle(fout)
	
		#if iS<39:
		#continue
		#break
		
		# Download annual files
		for yr in uY:
			if (yr<dSL['First Year'][indS]) | (yr>dSL['Last Year'][indS]):
				continue
			print('Working on: ' + str(uS[iS]) + ' ' + str(yr))
			url='http://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID=' + str(uS[iS]) + '&Year=' + str(yr) + '&Month=1&Day=14&timeframe=2&submit=Download+Data'
	
			try:
				testfile=urllib.request.urlretrieve(url,meta['Paths']['ECCD'] + '\\Daily\Raw\\file.csv')
			except:
				print('Failed to find download file.')
				ListMissing.append(uS[iS])
				continue
	
		df=pd.read_csv(meta['Paths']['ECCD'] + '\\Daily\Raw\\file.csv')
		ind=np.where( (tv[:,0]==yr) )[0]
		for v in vL:
			y=df[v].values
			if y.dtype==object:
				ikp=np.where(df[v].notna()==True)[0]
			else:
				ikp=np.where(np.isnan(y)==False)[0]
		
			try:
				d[v][ind[ikp]]=df[v].values[ikp]
			except:
				print('Unknown characters!')
		
		#ind=np.where(tv[:,0]==2023)[0]
		#plt.plot(d['Mean Temp (°C)'][ind],'b-')
		
		# Save
		gu.opickle(fout,d)
	
	return

#%% Gap-fill reference stations
def GapFillDailyECCD(meta):
	# Import reference stations to gap fill
	dRS=gu.ReadExcel(meta['Paths']['ECCD'] + '\\Daily\\Reference Station List for Gap Filling.xlsx')
	# Fix ID of reference stations
	for i in range(dRS['ID'].size):
		dRS['ID'][i]=dRS['ID'][i][1:-1]

	# Import station list
	dSL=gu.ReadExcel(meta['Paths']['ECCD'] + '\\Daily\\Station Inventory EN.xlsx')
	# Fix climate ID in dictionary -> excel messes it up with scientific notation
	dfSL=pd.read_csv(meta['Paths']['ECCD'] + '\\Daily\\Station Inventory EN.csv',skiprows=3)
	dSL['Climate ID']=dfSL['Climate ID'].values
	
	# Add coordinates for reference stations
	dRS['Lat']=np.zeros(dRS['ID'].size)
	dRS['Lon']=np.zeros(dRS['ID'].size)
	for i in range(dRS['ID'].size):
		ind=np.where(dSL['Climate ID']==dRS['ID'][i])[0]
		dRS['Lat'][i]=dSL['Latitude (Decimal Degrees)'][ind]
		dRS['Lon'][i]=dSL['Longitude (Decimal Degrees)'][ind]
	srs=gis.ImportSRSs()
	dRS['X'],dRS['Y']=gis.ReprojectCoordinates(srs['Proj']['Geographic'],srs['Proj']['NACID'],dRS['Lon'],dRS['Lat'])

	dRS['Neighbour List']=[None]*dRS['ID'].size
	for i in range(dRS['ID'].size):
		dRS['Neighbour List'][i]=np.zeros(10)

	# Import locations of all data
	f=os.listdir(meta['Paths']['ECCD'] + '\\Daily\\Raw\\Combined')
	dAll={'Name':np.array(['' for _ in range(len(f))],dtype=object),'Lat':np.zeros(len(f)),'Lon':np.zeros(len(f))}
	for i in range(len(f)):
		ind=np.where(dSL['Climate ID']==f[i][:-4])[0]
		if ind.size==0:
			print('Missing ' + str(f[i][:-4]))
			continue
		dAll['Name'][i]=dSL['Climate ID'][ind[0]]
		dAll['Lat'][i]=dSL['Latitude (Decimal Degrees)'][ind]
		dAll['Lon'][i]=dSL['Longitude (Decimal Degrees)'][ind]
	ind=np.where( (dAll['Lat']>0) & (dAll['Lon']>-180) )
	for k in dAll.keys():
		dAll[k]=dAll[k][ind]
	dAll['X'],dAll['Y']=gis.ReprojectCoordinates(srs['Proj']['Geographic'],srs['Proj']['NACID'],dAll['Lon'],dAll['Lat'])
	#gu.opickle(meta['Paths']['ECCD'] + '\\Daily\\Data All.pkl',dAll)

	# Get Climate IDs for ten closest neighbours for each subset station
	uID=np.array([])
	for i in range(dRS['ID'].size):
		Dist=np.sqrt((dRS['X'][i]-dAll['X'])**2+(dRS['Y'][i]-dAll['Y'])**2)
		Ord=np.argsort(Dist)
		dRS['Neighbour List'][i]=dAll['Name'][Ord[0:10]]
		uID=np.append(uID,dAll['Name'][Ord[0:10]])
	uID=np.unique(uID)

	# Import monthly data for neighbouring stations
	th_GapDaily=3
	#vL=['tmean','tmin','tmax','prcp','rain','sdep']
	vL=['tmean','prcp']
	dN={}
	for v in vL:
		dN[v]=meta['Missing Number']*np.ones((meta['tvm'].shape[0],uID.size))
	for iU in range(uID.size):
		print(iU)
		indSL=np.where(dSL['Climate ID']==uID[iU])[0]
		tvd0=gu.tvec('d',dSL['First Year'][indSL],dSL['Last Year'][indSL])
		d0=gu.ipickle(meta['Paths']['ECCD'] + '\\Daily\\Raw\\Combined\\' + uID[iU] + '.pkl')
		for Year in meta['tva']:
			iT2=np.where( (tvd0[:,0]==Year) )[0]
			tvd1=tvd0[iT2,:]
			d1={}
			d1['Mean Temp (°C)']=d0['Mean Temp (°C)'][iT2]
			d1['Total Precip (mm)']=d0['Total Precip (mm)'][iT2]
			for Month in range(1,13):
				iT1=np.where( (meta['tvm'][:,0]==Year) & (meta['tvm'][:,1]==Month) )[0]
				iT3=np.where( (tvd1[:,1]==Month) )[0]
				if iT3.size==0:
					continue
				ikp=np.where(d1['Mean Temp (°C)'][iT3]>-65)[0]
				if iT3.size-ikp.size<=th_GapDaily:
					dN['tmean'][iT1,iU]=np.mean(d1['Mean Temp (°C)'][iT3[ikp]])
				ikp=np.where(d1['Total Precip (mm)'][iT3]>=0)[0]
				if iT3.size-ikp.size<=th_GapDaily:
					dN['prcp'][iT1,iU]=np.sum(d1['Total Precip (mm)'][iT3[ikp]])
	gu.opickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Neighbours.pkl',dN)

	# Import monthly data for reference stations
	dRS['Data']={}
	for v in vL:
		dRS['Data'][v]=meta['Missing Number']*np.ones((meta['tvm'].shape[0],dRS['ID'].size))
	for i in range(dRS['ID'].size):
		indSL=np.where(dSL['Climate ID']==dRS['ID'][i])[0]
		tvd0=gu.tvec('d',dSL['First Year'][indSL],dSL['Last Year'][indSL])
		d0=gu.ipickle(meta['Paths']['ECCD'] + '\\Daily\\Raw\\Combined\\' + dRS['ID'][i] + '.pkl')
		for Year in meta['tva']:
			iT2=np.where( (tvd0[:,0]==Year) )[0]
			tvd1=tvd0[iT2,:]
			d1={}
			d1['Mean Temp (°C)']=d0['Mean Temp (°C)'][iT2]
			d1['Total Precip (mm)']=d0['Total Precip (mm)'][iT2]
			for Month in range(1,13):
				iT1=np.where( (meta['tvm'][:,0]==Year) & (meta['tvm'][:,1]==Month) )[0]
				iT3=np.where( (tvd1[:,1]==Month) )[0]
				if iT3.size==0:
					continue
				ikp=np.where(d1['Mean Temp (°C)'][iT3]>-999)[0]
				if iT3.size-ikp.size<=th_GapDaily:
					dRS['Data']['tmean'][iT1,i]=np.mean(d1['Mean Temp (°C)'][iT3[ikp]])
				ikp=np.where(d1['Total Precip (mm)'][iT3]>-999)[0]
				if iT3.size-ikp.size<=th_GapDaily:
					dRS['Data']['prcp'][iT1,i]=np.sum(d1['Total Precip (mm)'][iT3[ikp]])
	gu.opickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Ref Stations.pkl',dRS)

	# Gap fill reference stations
	dRS['Data GF']={}
	for v in vL:
		nN=dN[v].shape[1]
		dRS['Data GF'][v]=dRS['Data'][v].copy()
		for i in range(dRS['ID'].size):
			y=dRS['Data'][v][:,i]
			N_Overlap=np.zeros(nN)
			mec=np.zeros(nN)
			for j in range(nN):
				x=dN[v][:,j]

				iCal=np.where( (x>=-65) & (y>=-65) )[0]
				if iCal.size<12:
					continue
				rs,txt=gu.GetRegStats(x[iCal],y[iCal])
				mec[j]=rs['MEC']
				N_Overlap[j]=iCal.size
			Score=mec
			ord=np.flip(np.argsort(Score))
			
			for j in range(1,300):
				#if mec[ord[j]]<0.5:
				#print('Warning - using MEC < 0.5')
				x=dN[v][:,ord[j]]
				iFil=np.where( (x>=-65) & (y<-65) & (dRS['Data GF'][v][:,i]<-65) )[0]
				dRS['Data GF'][v][iFil,i]=x[iFil]
				#N_mis0=np.sum(dRS['Data GF'][v][:,i]<-65)
				#print(str(N_mis0) + ' ' + str(mec[ord[j]]))
			N_mis0=np.sum(dRS['Data GF'][v][:,i]<-65)
			print( str(i) + ' ' + str(N_mis0) )

		gu.opickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Ref Stations Gap Filled.pkl',dRS)

	return

#%%
def ImportUSHCN(meta):
	# Recent data (long term data nolonger needed as the recent data goes back)
	
	# Manually download (https://www.ncei.noaa.gov/pub/data/ushcn/v2.5/)
	tar=tarfile.open(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\ushcn.prcp.latest.FLs.52j.tar.gz',"r:gz")
	tar.extractall(r'D:\Data\Climate\US\USHCN\Recent Serial\2023')
	tar.close()
	tar=tarfile.open(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\ushcn.tavg.latest.FLs.52j.tar.gz',"r:gz")
	tar.extractall(r'D:\Data\Climate\US\USHCN\Recent Serial\2023')
	tar.close()
	
	# Import station list
	dSL=gu.ReadExcel(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\ushcn-v2.5-stations.xlsx')
	srs=gis.ImportSRSs()
	dSL['X'],dSL['Y']=gis.ReprojectCoordinates(srs['Proj']['Geographic'],srs['Proj']['NACID'],dSL['Lon'],dSL['Lat'])
	
	# Import data
	dSL['tmean']=meta['Missing Number']*np.ones((meta['tvm'].shape[0],dSL['ID'].size))
	dSL['prcp']=meta['Missing Number']*np.ones((meta['tvm'].shape[0],dSL['ID'].size))
	f=os.listdir(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\ushcn.v2.5.5.20230910')
	w=[11,1,4,1,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4,5,4]
	for fn in f:
		df=pd.read_fwf(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\ushcn.v2.5.5.20230910' + '\\' + fn,widths=w,header=None)
		iS=np.where(dSL['ID']==df[0][0])[0]
		Year=df[2].values
		Data=np.ones((Year.size,12))
	cnt=0
	for i in range(4,28,2):
		Data[:,cnt]=df[i].values
		cnt=cnt+1
	for iY in range(Year.size):
		iT=np.where(meta['tvm'][:,0]==Year[iY])[0]
		if iT.size==0:
			continue
		if fn[-4:]=='prcp':
			dSL['prcp'][iT,iS]=Data[iY,:]/10
		else:
			dSL['tmean'][iT,iS]=Data[iY,:]/100
	
	ind=np.where(dSL['prcp']<0)
	dSL['prcp'][ind]=np.nan
	plt.close('all')
	plt.plot(np.nanmean(dSL['prcp'],axis=1))
	
	ind=np.where(dSL['tmean']<-25)
	dSL['tmean'][ind]=np.nan
	plt.close('all')
	plt.plot(np.nanmean(dSL['tmean'],axis=1))
	
	# Save
	gu.opickle(r'D:\Data\Climate\US\USHCN\Recent Serial\2023\Data.pkl',dSL)
	
	return

#%%
def ImportReanalysis(meta,vr,form):
	# https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface_gauss/
	# http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.surfaceflux.html

	# Import data
	if vr=='tmean':
		ds=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\air.2m.mon.mean.nc')
	elif vr=='prcp':
		ds=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\prate.sfc.mon.mean.nc')
	elif vr=='ea':
		ds=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\shum.2m.mon.mean.nc')
		dsT=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\air.2m.mon.mean.nc')
	elif vr=='vpd':
		ds=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\shum.2m.mon.mean.nc')
		dsT=nc.Dataset(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\air.2m.mon.mean.nc')
	else:
		return
	# ds.variables.keys()
	lat=np.array(list(ds.variables['lat'][:]))
	lon=np.array(list(ds.variables['lon'][:]))
	ind=np.where(lon>180)[0]
	lon[ind]=lon[ind]-360
	lon,lat=np.meshgrid(lon,lat)
	srs=gis.ImportSRSs()
	xR,yR=gis.ReprojectCoordinates(srs['Proj']['Geographic'],srs['Proj']['NACID'],lon,lat)
	tvR=gu.tvec('m',1948,2035)
	tvR=tvR[:ds.variables['time'][:].size,:]

	# Conversions
	if vr=='tmean':
		zR=ds.variables['air'][:]-273.15
	elif vr=='prcp':
		# Convert from kg m-2 s-1 to mm month-1
		zR=ds.variables['prate'][:]*86400
		for mo in range(12):
			zR[mo::12,:,:]=zR[mo::12,:,:]*meta['Parameters']['DIM'][mo]
	elif vr=='ea':
		# Convert air temperature from K to C
		tmeanR=dsT.variables['air'][:]-273.15
		rhoAir=gaia.SF_get_rho_air(tmeanR) # Density of air (kg m-3)
		estarR=gaia.GetEstar(tmeanR) # Saturation vapour pressure (hPa)

		# Convert specific humidity (given in g/kg) to
		shumR=ds.variables['shum'][:]*1000
		eaR=shumR/rhoAir*461.5*(tmeanR+273.15)
		vpdR=np.maximum(0,(estarR-eaR)/100) # Pa
		zR=eaR/100 # Vapour pressure (hPa)
	elif vr=='vpd':
		# Convert air temperature from K to C
		tmeanR=dsT.variables['air'][:]-273.15
		rhoAir=gaia.GetDensityAir(tmeanR) # Density of air (kg m-3)
		#rhoAir=meta['Parameters']['RhoAir'] # Density of air (kg m-3)
		estarR=gaia.GetEstar(tmeanR) # Saturation vapour pressure (hPa)
		# Convert specific humidity (given in g/kg) to kg/kg
		shumR=ds.variables['shum'][:]/1000
		eaR=(shumR/rhoAir*461.5*(tmeanR+273.15))/100 # Converted to hPa
		zR=np.maximum(0,estarR-eaR) # Vapour pressure (hPa)

	# Interpolate
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	buff=200000
	indR=np.where( (xR>=zRef['xmin']-buff) & (xR<=zRef['xmax']+buff) & (yR>=zRef['ymin']-buff) & (yR<=zRef['ymax']+buff) )
	spc=50
	xR2=zRef['X'][np.ix_(np.linspace(0,zRef['m']-1,spc).astype(int),np.linspace(0,zRef['n']-1,spc).astype(int))]
	yR2=zRef['Y'][np.ix_(np.linspace(0,zRef['m']-1,spc).astype(int),np.linspace(0,zRef['n']-1,spc).astype(int))]
	zR2=np.zeros(shape=(zR.shape[0],xR2.shape[0],xR2.shape[1]))
	xy0=np.column_stack((xR[indR],yR[indR]))
	for i in range(tvR.shape[0]):
		z0=zR[i,:,:][indR]
		zR2[i,:,:]=griddata(xy0,z0,(xR2,yR2),method='linear')

	# Convert to anomalies
	if form=='anom':
		for mo in range(12):
			iMo=np.where( (tvR[:,1]==mo+1) )[0]
			ind=np.where( (tvR[:,0]>=1971) & (tvR[:,0]<=2000) & (tvR[:,1]==mo+1) )[0]
			mu=np.mean(zR2[ind,:,:],axis=0)
			zR2[iMo,:,:]=zR2[iMo,:,:]-np.tile(mu,(iMo.size,1,1))

	return xR2,yR2,zR2,tvR

#%%
def CalcTemperatureAnomalies(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])

	# Import station data and combine
	dCA=gu.ipickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Ref Stations Gap Filled.pkl')
	dUS=gu.ipickle(r'C:\Data\Climate\USHCN\Recent Serial\2023\Data.pkl')

	# Combine
	dS={}
	dS['tmean']=np.column_stack((dCA['Data GF']['tmean'],dUS['tmean']))
	dS['X']=np.append(dCA['X'],dUS['X'])
	dS['Y']=np.append(dCA['Y'],dUS['Y'])

	ind=np.where(dS['tmean']<-55)
	dS['tmean'][ind]=np.nan

	# Convert to anomalies 1971-2000 base period
	for mo in range(12):
		iMo=np.where( (meta['tvm'][:,1]==mo+1) )[0]
		ind=np.where( (meta['tvm'][:,0]>=1971) & (meta['tvm'][:,0]<=2000) & (meta['tvm'][:,1]==mo+1) )[0]
		mu=np.nanmean(dS['tmean'][ind,:],axis=0)
		dS['tmean'][iMo,:]=dS['tmean'][iMo,:]-np.tile(mu,(iMo.size,1))

	# Import reanalysis
	xR,yR,zR,tvR=ImportReanalysis(meta,'tmean','anom')

	# Interpolation
	for iT in range(meta['tvm'].shape[0]):

		zS0=dS['tmean'][iT,:]
		ikp=np.where( (zS0>-45) & (zS0<45) )[0]
		zS0=zS0[ikp]
		xS0,yS0=dS['X'][ikp],dS['Y'][ikp]

		indR=np.where( (tvR[:,0]==meta['tvm'][iT,0]) & (tvR[:,1]==meta['tvm'][iT,1]) )[0]
		zR0=zR[indR[0],:,:]
		MinDist=np.zeros(xR.shape)
		for i in range(xR.shape[0]):
			for j in range(xR.shape[1]):
				d=np.sqrt((xS0-xR[i,j])**2+(yS0-yR[i,j])**2)
				MinDist[i,j]=np.min(d)
		iXYR=np.where( (MinDist>100000) & (np.isnan(zR0)==False) )
		xR0=xR[iXYR]
		yR0=yR[iXYR]
		zR0=zR0[iXYR]

		x0=np.append(xS0,xR0)
		y0=np.append(yS0,yR0)
		z0=np.append(zS0,zR0)

		flg=0
		if flg==1:
			plt.close('all')
			plt.plot(xS0,yS0,'b.')
			plt.plot(xR0,yR0,'g.')

		# Interpolation parameters for errorgram
		beta=np.array([0.000038,0.75,0.025])

		mask=np.zeros(zRef['Data'].shape,dtype='int8')
		ind=np.where(zRef['Data']>0)
		mask[ind]=1
		kernel=np.ones((50,50))
		mask=cv2.dilate(mask.astype(np.uint8),kernel,iterations=1)

		itvl=50
		type='Hyperbolic'
		xI=zRef['X'][0::itvl,0::itvl]
		yI=zRef['Y'][0::itvl,0::itvl]
		maskI=mask[0::itvl,0::itvl]

		iZ,iE,iN=gu.PolynomialSurfaceFit(x0/1000,y0/1000,z0,xI/1000,yI/1000,maskI,beta)
		#plt.close('all'); plt.matshow(iZ,clim=[-2,6])

		iZ2=griddata(np.column_stack((xI.flatten(),yI.flatten())),iZ.flatten(),(zRef['X'],zRef['Y']),method='linear')
		ind=np.where(zRef['Data']==0)
		iZ2[ind]=meta['Missing Number']
		#plt.close('all'); plt.matshow(iZ2,clim=[-6,6])
		z=copy.deepcopy(zRef)
		z['Data']=iZ2/meta['SF']['tmean']
		z['Data']=z['Data'].astype('int16')
		z['Data']=z['Data']*zRef['Data']
		#plt.close('all'); plt.matshow(z['Data'],clim=[-800,800])
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_tmean_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')

	return

#%%
def CalcPrecipAnomalies(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	zRefB=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid Buf'])

	# Import station data and combine
	dCA=gu.ipickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Ref Stations Gap Filled.pkl')
	dUS=gu.ipickle(r'C:\Data\Climate\USHCN\Recent Serial\2023\Data.pkl')

	# Combine
	dS={}
	dS['prcp']=np.column_stack((dCA['Data GF']['prcp'],dUS['prcp']))
	dS['X']=np.append(dCA['X'],dUS['X'])
	dS['Y']=np.append(dCA['Y'],dUS['Y'])
	ind=np.where(dS['prcp']<0)
	dS['prcp'][ind]=np.nan

	# Convert to anomalies 1971-2000 base period
	for mo in range(12):
			iMo=np.where( (meta['tvm'][:,1]==mo+1) )[0]
			ind=np.where( (meta['tvm'][:,0]>=1971) & (meta['tvm'][:,0]<=2000) & (meta['tvm'][:,1]==mo+1) )[0]
			mu=np.mean(dS['prcp'][ind,:],axis=0)
			dS['prcp'][iMo,:]=dS['prcp'][iMo,:]-np.tile(mu,(iMo.size,1))

	# Import reanalysis
	xR,yR,zR,tvR=ImportReanalysis(meta,'prcp','anom')

	# Interpolation
	for iT in range(meta['tvm'].shape[0]):
	#t=np.flip(np.arange(0,meta['tvm'].shape[0]))
	#for iT0 in range(t.size):
		#iT=t[iT0]
		#print(iT)
		zS0=dS['prcp'][iT,:]
		ikp=np.where( (zS0>-500) & (zS0<500) )[0]
		zS0=zS0[ikp]
		xS0,yS0=dS['X'][ikp],dS['Y'][ikp]

		indR=np.where( (tvR[:,0]==meta['tvm'][iT,0]) & (tvR[:,1]==meta['tvm'][iT,1]) )[0]
		if indR.size==0:
			continue

		zR0=zR[indR[0],:,:]
		MinDist=np.zeros(xR.shape)
		for i in range(xR.shape[0]):
			for j in range(xR.shape[1]):
				d=np.sqrt((xS0-xR[i,j])**2+(yS0-yR[i,j])**2)
				MinDist[i,j]=np.min(d)
		iXYR=np.where( (MinDist>100000) & (np.isnan(zR0)==False) )
		xR0=xR[iXYR]
		yR0=yR[iXYR]
		zR0=zR0[iXYR]

		x0=np.append(xS0,xR0)
		y0=np.append(yS0,yR0)
		z0=np.append(zS0,zR0)

		flg=0
		if flg==1:
			plt.close('all')
			plt.plot(xS0,yS0,'b.')
			plt.plot(xR0,yR0,'g.')

		# Interpolation parameters for errorgram
		beta=np.array([0.0008,17,1.1])

		itvl=50
		#type='Hyperbolic'
		xI=zRef['X'][0::itvl,0::itvl]
		yI=zRef['Y'][0::itvl,0::itvl]
		maskI=zRefB['Data'][0::itvl,0::itvl]

		iZ,iE,iN=gu.PolynomialSurfaceFit(x0/1000,y0/1000,z0,xI/1000,yI/1000,maskI,beta)
		#plt.close('all'); plt.matshow(iZ,clim=[-200,200])

		iZ2=griddata(np.column_stack((xI.flatten(),yI.flatten())),iZ.flatten(),(zRef['X'],zRef['Y']),method='linear')
		iZ2[zRefB['Data']==0]=meta['Missing Number']
		#plt.close('all'); plt.matshow(iZ2,clim=[-200,200])
		z=copy.deepcopy(zRef)
		z['Data']=iZ2
		z['Data']=z['Data']/meta['SF']['prcp']
		z['Data']=z['Data'].astype('int16')
		#plt.close('all'); plt.matshow(z['Data'],clim=[-2000,2000])
		gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_prcp_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')
	return

# for iT in range(meta['tvm'].shape[0]):
# 	z=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_prcp_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')
# 	z['Data']=z['Data'].astype('float')*meta['SF']['prcp']
# 	z['Data']=z['Data']*meta['Parameters']['DIM'][mo]

# 	z['Data']=z['Data']/meta['SF']['prcp']
# 	z['Data']=z['Data'].astype('int16')
# 	gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_prcp_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')

#%% Import CRU Time Series dataset
def CalcTempAnomaliesFromCRU(metaNA):
	# Given in numpy mask arrays so you need to use the filled method
	srs=gis.ImportSRSs()
	zRef=gis.OpenGeoTiff(metaNA['Paths']['na1k Ref Grid'])
	tv=gu.tvec('m',1901,2023)
	
	# Import netcdf data
	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.tmp.dat.nc.gz'
	tmp0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					tmp0[k]=ds.variables[k][:]
					#print(np.mean(ds.variables[k][:]))
	tmp1=tmp0['tmp'].filled()
	ind=np.where(tmp1>999)
	tmp1[ind]=metaNA['Missing Number']

	# Convert to anomalies
	for mo in range(12):
		#itmu=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo) )
		mu=np.mean(tmp1[mo::12,:,:],axis=0)
		tmp1[mo::12,:,:]=tmp1[mo::12,:,:]-np.tile(mu,(int(tmp1.shape[0]/12),1,1))
	
	# Isolate BC area  (so that interpolation is faster)
	lon,lat=np.meshgrid(tmp0['lon'].filled(),tmp0['lat'].filled(),sparse=False)
	x1,y1=srs['Proj']['NACID'](lon,lat)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))
	
	for iT in range(tv.shape[0]):
		z0=tmp1[iT,:,:][indS]
		z1=copy.deepcopy(zRef)
		z1['Data']=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-6,6])
		z1['Data']=z1['Data']/metaNA['Climate']['SF']['vpd']
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-600,600])
		gis.SaveGeoTiff(z1,metaNA['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_tmean_anom_cru_' + str(tv[iT,0]) + '_' + str(tv[iT,1]) + '.tif')
	return

#%% Import CRU Time Series dataset
def CalcVapourPressureDefictAnomaliesFromCRU(meta):
	# Given in numpy mask arrays so you need to use the filled method
	srs=gis.ImportSRSs()
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	tv=gu.tvec('m',1901,2023)
	
	# Import netcdf data
	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.vap.dat.nc.gz'
	vap0={}
	with gzip.open(fin) as gz:
		with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
			for k in ds.variables.keys():
				vap0[k]=ds.variables[k][:]
	vap0['vap']=vap0['vap'].filled()
	ind=np.where(vap0['vap']>999)
	vap0['vap'][ind]=meta['Missing Number']
	ind=np.where(vap0['vap']>meta['Missing Number'])
	np.mean(vap0['vap'][ind])
	
	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.tmp.dat.nc.gz'
	tmp0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					if k!='tmp':
						continue
					tmp0[k]=ds.variables[k][:]
					#print(np.mean(ds.variables[k][:]))
	tmp0=tmp0['tmp'].filled()
	ind=np.where(tmp0>999)
	tmp0[ind]=meta['Missing Number']
	esat0=gaia.GetEstar(tmp0)
	ind=np.where(esat0>999999)
	esat0[ind]=meta['Missing Number']
	
	vpd0=np.maximum(0,esat0-vap0['vap'])
	ind=np.where(vpd0>999)
	vpd0[ind]=meta['Missing Number']
	
	# Convert to anomalies
	for mo in range(12):
		#itmu=np.where( (tv[:,0]>=1971) & (tv[:,0]<=2000) & (tv[:,1]==mo) )
		mu=np.mean(vpd0[mo::12,:,:],axis=0)
		vpd0[mo::12,:,:]=vpd0[mo::12,:,:]-np.tile(mu,(int(vpd0.shape[0]/12),1,1))
	
	# Isolate BC area  (so that interpolation is faster)
	lon,lat=np.meshgrid(vap0['lon'].filled(),vap0['lat'].filled(),sparse=False)
	x1,y1=srs['Proj']['NACID'](lon,lat)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))
	
	for iT in range(tv.shape[0]):
		if tv[iT,0]<1950:
			continue
		z0=vpd0[iT,:,:][indS]
		z1=copy.deepcopy(zRef)
		z1['Data']=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-6,6])
		z1['Data']=z1['Data']/meta['SF']['vpd']
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-600,600])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_vpd_anom_' + str(tv[iT,0]) + '_' + str(tv[iT,1]) + '.tif')
	return

#%%
def CalcVapourPressureDefictAnomaliesFromReanalysis(meta,year):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	#year=2023
	xR,yR,zR,tvR=ImportReanalysis(meta,'vpd','anom')
	for iT in range(meta['tvm'].shape[0]):
		if meta['tvm'][iT,0]!=year:
			continue
		indT=np.where( (tvR[:,0]==meta['tvm'][iT,0]) & (tvR[:,1]==meta['tvm'][iT,1]) )[0]
		zR0=zR[indT[0],:,:]
		z1=copy.deepcopy(zRef)
		z1['Data']=griddata(np.column_stack((xR.flatten(),yR.flatten())),zR0.flatten(),(zRef['X'],zRef['Y']),method='linear')
		ind=np.where(zRef['Data']==0)
		z1['Data'][ind]=0
		z1['Data']=z1['Data']/meta['SF']['vpd']
		z1['Data']=z1['Data'].astype('int16')
		#plt.matshow(z1['Data'],clim=[-20,20])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_vpd_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')
	return

#%%
def CalcSurfaceWaterBalance(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	iMask=np.where(zRef['Data']>0)

	# Parameters
	par={}
	par['Method']='Grid'
	par['Ws_max']=200.0 # mm
	par['Tmin']=-3.0
	par['Tmax']=3.0
	par['Ei_FracMax']=0.15
	par['Ei_ALMax']=5.0
	par['ETp Method']='Penman-Monteith'
	par['Daily_Interval']=5
	par['Include Rainfall Fraction']='No'

	# Import normals
	nrms={}
	nrms['tmean']=[None]*12
	nrms['rswn']=[None]*12
	nrms['prcp']=[None]*12
	nrms['vpd']=[None]*12
	for mo in range(12):
		nrms['tmean'][mo]=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype('float')*meta['SF']['tmean']
		nrms['rswn'][mo]=(1-meta['Parameters']['Albedo']['Forest Coniferous'])*gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_rswd_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype('float')*meta['SF']['rswd']
		nrms['prcp'][mo]=np.maximum(0,gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_prcp_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype('float')*meta['SF']['prcp'])
		nrms['vpd'][mo]=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_vpd_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype('float')*meta['SF']['vpd']

	#plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(14,14)) # See if it is working
	#cnt=0
	vi={}
	for iT in range(meta['tvm'].shape[0]):
		if meta['tvm'][iT,0]<1950:
			continue
		print('Year:' + str(meta['tvm'][iT,0]) + ', Month:' + str(meta['tvm'][iT,1]) )
		mo=meta['tvm'][iT,1]
		iMo=mo-1
		vi['Month']=mo
		vi['LAI']=5.0
		vi['Gs']=0.010
		vi['Ga']=0.058
		vi['tmean']=np.maximum(-50,nrms['tmean'][iMo]+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_tmean_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')['Data'][iMask].astype('float')*meta['SF']['tmean'])
		vi['prcp']=np.maximum(0,nrms['prcp'][iMo]+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_prcp_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')['Data'][iMask].astype('float')*meta['SF']['prcp'])
		vi['vpd']=np.maximum(0.01,nrms['vpd'][iMo]+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_vpd_anom_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')['Data'][iMask].astype('float')*meta['SF']['vpd'])
		vi['rswn']=nrms['rswn'][mo-1]

		vo=gaia.WBM_SparseGrid(par,vi)

		# Populate inputs for next time step
		vi['ws']=vo['ws']
		vi['wsp']=vo['wsp']

		# Add climatic water deficit
		vo['cwd']=vo['etp']-vo['eta']

		# Add climatic moisture index
		vo['cmi']=vi['prcp']-vo['etp']

		#ax[0,0].plot(cnt,np.mean(vo['etp']),'b.')
		#ax[0,0].plot(cnt,np.mean(vo['eta']),'go')
		#ax[0,1].plot(cnt,np.mean(vo['cwd']),'b.')
		#ax[1,0].plot(cnt,np.mean(vo['wsp']),'b.')
		#ax[1,1].plot(cnt,np.mean(vo['ws']),'b.')
		#cnt=cnt+1

		# Save variables to file
		for v in vo.keys():

			if v=='etp':
				# Scale factor cannot handle values above 325/month
				vo[v]=np.minimum(325,vo[v])

			z1=copy.deepcopy(zRef)
			z1['Data']=meta['Missing Number']*np.ones(zRef['Data'].shape,dtype='float')
			z1['Data'][iMask]=vo[v]
			z1['Data'][iMask]=z1['Data'][iMask]/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			# plt.close('all'); plt.matshow(z1['Data'],clim=[0,6])
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Monthly\\Actual\\na1k_' + v + '_actual_' + str(meta['tvm'][iT,0]) + '_' + str(meta['tvm'][iT,1]) + '.tif')

	return

#%% Calculate Thornthwaite's Climate Classification System
def CalcThornthwaitesClassification(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	
	pet=np.zeros(zRef['Data'].shape,dtype='float')
	prcp=np.zeros(zRef['Data'].shape,dtype='float')
	defi=np.zeros(zRef['Data'].shape,dtype='float')
	sur=np.zeros(zRef['Data'].shape,dtype='float')
	for mo in range(12):
	  # Version 2:
	  pet0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.tif')
	  pet0=pet0['Data'].astype('float32')/100*30
	  pet=pet+pet0
	
	  prcp0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
	  prcp0=prcp0['Data'].astype('float32')
	  prcp=prcp+prcp0
	
	  ind=np.where(pet0>prcp0)
	  defi[ind]=defi[ind]+pet0[ind]-prcp0[ind]
	
	  ind=np.where(pet0<prcp0)
	  sur[ind]=sur[ind]+prcp0[ind]-pet0[ind]
	
	# Classification
	
	z0=((sur-defi)/pet)*100
	
	lab=['Hyper humid','Humid','Subhumid wet','Subhumid dry','Semi arid','Arid','Outside Boundary']
	
	z1=7*np.ones(z0.shape,dtype='int8')
	ind=np.where(z0>=100); z1[ind]=1;
	ind=np.where( (z0>=20) & (z0<100) ); z1[ind]=2;
	ind=np.where( (z0>=0) & (z0<20) ); z1[ind]=3;
	ind=np.where( (z0>=-20) & (z0<0) ); z1[ind]=4;
	ind=np.where( (z0>=-40) & (z0<-20) ); z1[ind]=5;
	ind=np.where(z0<-40); z1[ind]=6;
	ind=np.where(zRef['Data']==0); z1[ind]=7;
	
	# plt.matshow(z1,clim=[1,7])
	
	# Save
	z2=zRef.copy()
	z2['Data']=z1
	gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ThornthwaiteClimateClass_norm_1971to2000.tif')
	
	# Classification (condensed)
	
	z0=((sur-defi)/pet)*100
	
	lab=['Hyper humid','Humid','Subhumid','Semi arid','Arid','Outside Boundary']
	
	z1=6*np.ones(z0.shape,dtype='int8')
	ind=np.where(z0>=100); z1[ind]=1;
	ind=np.where( (z0>=20) & (z0<100) ); z1[ind]=2;
	ind=np.where( (z0>=-20) & (z0<20) ); z1[ind]=3;
	ind=np.where( (z0>=-40) & (z0<-20) ); z1[ind]=4;
	ind=np.where(z0<-40); z1[ind]=5;
	ind=np.where(zRef['Data']==0); z1[ind]=6;
	
	# plt.matshow(z1,clim=[1,7])
	
	# Save
	z2=zRef.copy()
	z2['Data']=z1
	gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ThornthwaiteClimateClassCondensed_norm_1971to2000.tif')
	return

#%%
def CalcSummariesBasicVariables(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	#iOut=np.where(zRef['Data']==0)

	# Water Year
	mo=np.arange(1,13,1)
	vL=['tmean']
	for v in vL:
		# Import Normals
		zN=np.zeros(zRef['Data'].shape)
		for iM in range(mo.size):
			zN=zN+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
		zN=zN/mo.size

		for iT in range(meta['tva'].size):
			if (meta['tva'][iT]<2023):
				# Cannot calculate the first winter (missing previous december)
				continue
			if (meta['tva'][iT]>2023):
				# Cannot calculate the first winter (missing previous december)
				continue
			if meta['tva'][iT]<1951:
				continue
			# Import anomalies
			mo=np.array([10,11,12,1,2,3,4,5,6,7,8,9])
			zA=np.zeros(zRef['Data'].shape)
			for iM in range(mo.size):
				if iM<=2:
					zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]-1) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
				else:
					zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zA=zA/mo.size
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_wyr_' + str(meta['tva'][iT]) + '.tif')
	
			# Actuals
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zN+zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_wyr_' + str(meta['tva'][iT]) + '.tif')

	# Warm season
	mo=np.arange(5,10,1)
	vL=['tmean','prcp','rswd','vpd','ea_fromtmin','ea_biasadj']
	vL=['rswd','vpd','ea_fromtmin','ea_biasadj']
	vL=['prcp']
	for v in vL:
		# Import Normals
		zN=np.zeros(zRef['Data'].shape)
		for iM in range(mo.size):
			zN=zN+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
		zN=zN/mo.size

		# Save warm season normals
		z1=copy.deepcopy(zRef)
		z1['Data']=zN
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-10,12])
		z1['Data']=z1['Data']/meta['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_mjjas.tif')

		# Save warm season time series
		if np.isin(v,['rswd','vpd','ea_fromtmin','ea_biasadj'])==True:
			# Don't calculate time series for these variables
			continue

		for iT in range(meta['tva'].size):
			if meta['tva'][iT]>2023:
				continue
			if meta['tva'][iT]<1951:
				continue

			# Import anomalies
			zA=np.zeros(zRef['Data'].shape)
			for iM in range(mo.size):
				zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zA=zA/mo.size
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_mjjas_' + str(meta['tva'][iT]) + '.tif')

			# Calculate actuals from normals and anomalies
			zR=zN+zA
			if v=='prcp':
				# Scale factor cannot handle values above 225
				zR=np.minimum(325,zR)
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zR
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(meta['tva'][iT]) + '.tif')

	# Annual
	mo=np.arange(1,13,1)
	vL=['tmean','prcp']
	for v in vL:
		# Import Normals
		zN=np.zeros(zRef['Data'].shape)
		for iM in range(mo.size):
			zN=zN+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
		zN=zN/mo.size

		# Save normals
		z1=copy.deepcopy(zRef)
		z1['Data']=zN
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-10,12])
		z1['Data']=z1['Data']/meta['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_ann.tif')

		# 		if np.isin(v,['prcp'])==True:
		# 			# Don't calculate time series for these variables
		# 			continue

		for iT in range(meta['tva'].size):
			if (meta['tva'][iT]>2022):
				# Cannot calculate the first winter (missing previous december)
				continue
			if meta['tva'][iT]<1951:
				continue
			# Import anomalies
			zA=np.zeros(zRef['Data'].shape)
			for iM in range(mo.size):
				zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zA=zA/mo.size
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_ann_' + str(meta['tva'][iT]) + '.tif')
	
			# Actuals
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zN+zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_ann_' + str(meta['tva'][iT]) + '.tif')

	# Winter
	mo=np.array([12,1,2])
	vL=['tmean']
	for v in vL:
		# Import normals
		zN=np.zeros(zRef['Data'].shape)
		for iM in range(mo.size):
			zN=zN+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
		zN=zN/mo.size

		# Save normals
		z1=copy.deepcopy(zRef)
		z1['Data']=zN
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-10,12])
		z1['Data']=z1['Data']/meta['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_djf.tif')
	
		for iT in range(meta['tva'].size):
			if (meta['tva'][iT]>2023) | (iT==0):
				# Cannot calculate the first winter (missing previous december)
				continue
	
			# Import anomalies
			zA=np.zeros(zRef['Data'].shape)
			for iM in range(mo.size):
				if mo[iM]==12:
					zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]-1) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
				else:
					zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zA=zA/mo.size
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_djf_' + str(meta['tva'][iT]) + '.tif')
	
			# Actuals
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zN+zA
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_djf_' + str(meta['tva'][iT]) + '.tif')

	return

#%%
def CalcSummariesWB(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	iMask=np.where(zRef['Data']==1)
	iOut=np.where(zRef['Data']==0)
	tvN=np.arange(1971,2001,1)

	# Warm season (May-Sept)
	mo=np.arange(5,10,1)
	vL=['etp','ws','cwd','cmi']
	for v in vL:
		print(v)
		for iT in range(meta['tva'].size):
			if meta['tva'][iT]>2023:
				continue

			# Import actual monthly values directly from output of water balance model
			zR=np.zeros(zRef['Data'].shape)
			for iM in range(mo.size):
				zR=zR+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Actual\\na1k_' + v + '_actual_' + str(meta['tva'][iT]) + '_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zR=zR/mo.size
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=zR
			#plt.close('all'); plt.matshow(z1['Data'],clim=[0,12])
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(meta['tva'][iT]) + '.tif')

		# Save warm-season normals
		print(v)
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
		for iT in range(tvN.size):
 			z1['Data']=z1['Data']+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(tvN[iT]) + '.tif')['Data'].astype(float)*meta['SF'][v]
		z1['Data']=z1['Data']/tvN.size
		z1['Data']=z1['Data']/meta['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		#plt.matshow(z1['Data'],clim=[-40,40])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_mjjas.tif')

		# Save warm-season anomalies
		print('Anomalies:')
		zN=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_mjjas.tif')['Data'].astype(float)*meta['SF'][v]
		for iT in range(meta['tva'].size):
			if meta['tva'][iT]>2023:
				continue
			z0=(gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(meta['tva'][iT]) + '.tif')['Data'].astype(float)*meta['SF'][v])-zN
			z1=copy.deepcopy(zRef)
			z1['Data']=z0/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			#plt.matshow(z1['Data'],clim=[-40,40])
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_mjjas_' + str(meta['tva'][iT]) + '.tif')

		# Save warm-season standard deviation
		print('S.D.:')
		z0=np.zeros((tvN.size,iMask[0].size),dtype='float')
		for iT in range(tvN.size):
			z0[iT,:]=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(tvN[iT]) + '.tif')['Data'][iMask].astype(float)*meta['SF'][v]
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
		z1['Data'][iMask]=np.std(z0,axis=0)
		z1['Data']=z1['Data']/meta['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.matshow(z1['Data'],clim=[-40,40])
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_stdev_1971to2000_mjjas.tif')

		# Save warm-season z-scores
		print('Z-scores:')
		zN=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_mjjas.tif')['Data'].astype(float)*meta['SF'][v]
		zSD=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_stdev_1971to2000_mjjas.tif')['Data'].astype(float)*meta['SF'][v]
		for iT in range(meta['tva'].size):
			if meta['tva'][iT]>2023:
				continue
			z0=((gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Actual\\na1k_' + v + '_actual_mjjas_' + str(meta['tva'][iT]) + '.tif')['Data'].astype(float)*meta['SF'][v])-zN)/zSD
			z1=copy.deepcopy(zRef)
			z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
			z1['Data']=z0
			z1['Data']=z1['Data']/meta['SF'][v]
			z1['Data']=z1['Data'].astype('int16')
			#plt.matshow(z1['Data'],clim=[-40,40])
			gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Zscores\\na1k_' + v + '_zscore_mjjas_' + str(meta['tva'][iT]) + '.tif')

	# Save annual normal
	mo=np.arange(1,13,1)
	vL=['runoff']
	vL=['etp']
	for v in vL:
		z0=np.zeros(zRef['Data'].shape,dtype='float')
		for iM in range(mo.size):
			z0=z0+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['Climate']['SF'][v]
		z0=z0/mo.size
		z1=copy.deepcopy(zRef)
		z1['Data']=z0
		#plt.close('all'); plt.matshow(z1['Data'],clim=[-10,12])
		z1['Data']=z1['Data']/meta['Climate']['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_ann.tif')

	# Save winter normal
	mL=[12,1,2]
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	for i in range(len(mL)):
		z1['Data']=z1['Data']+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_wsp_norm_1971to2000_' + str(mL[i]) + '.tif')['Data'].astype(float)*meta['SF']['wsp']
	z1['Data']=z1['Data']/len(mL)
	z1['Data'][iOut]=meta['Missing Number']
	# plt.matshow(z1['Data'],clim=[0,500])
	z1['Data']=z1['Data']/meta['SF']['wsp']
	z1['Data']=z1['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_wsp_norm_1971to2000_djf.tif')

	# Save winter-spring normal
	mL=[1,2,3,4,5]
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	for i in range(len(mL)):
		z1['Data']=z1['Data']+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_melt_norm_1971to2000_' + str(mL[i]) + '.tif')['Data'].astype(float)*meta['SF']['melt']
	z1['Data']=z1['Data']/len(mL)
	z1['Data'][iOut]=meta['Missing Number']
	# plt.matshow(z1['Data'],clim=[0,5])
	z1['Data']=z1['Data']/meta['SF']['melt']
	z1['Data']=z1['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_melt_norm_1971to2000_jfmam.tif')

	return

#%%
def CalcTimeSeriesByEZ(meta):
	zPB=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\PoliticalBoundaries.tif')['Data']
	zEZ=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Ecozones_CEC_L1.tif')['Data']
	indPB=np.where( (zPB==meta['LUT']['PB']['CA-BC']) )
	zEZf=zEZ[indPB]
	uEZ=np.unique(zEZf)
	uEZ=uEZ[uEZ>0]
	d={}

	# Warm season
	vL=['tmean','prcp','cwd','cmi','ws']
	for i in range(uEZ.size):
		cd=lut_n2s(meta['LUT']['EZ'],uEZ[i])[0]
		d[cd]={'mjjas':{}}
		for v in vL:
			d[cd]['mjjas'][v]=np.nan*np.ones((meta['tva'].size))
	for v in vL:
		print(v)
		for iT in range(meta['tva'].size):
			try:
				z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_mjjas_' + str(meta['tva'][iT]) + '.tif')['Data'][indPB]
				if np.isin(v,['prcp','cwd','cmi'])==True:
					z0=z0/30.5
			except:
				pass
			for iU in range(uEZ.size):
				cd=lut_n2s(meta['LUT']['EZ'],uEZ[iU])[0]
				ind=np.where( (zEZf==uEZ[iU]) )
				d[cd]['mjjas'][v][iT]=np.mean(z0[ind].astype('float')*meta['SF'][v])

	# Annual
	vL=['tmean']
	for i in range(uEZ.size):
		cd=lut_n2s(meta['LUT']['EZ'],uEZ[i])[0]
		d[cd]['ann']={}
		for v in vL:
			d[cd]['ann'][v]=np.nan*np.ones((meta['tva'].size))
	for v in vL:
		print(v)
		for iT in range(meta['tva'].size):
			try:
				z0=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_' + v + '_anom_ann_' + str(meta['tva'][iT]) + '.tif')['Data'][indPB]
				if np.isin(v,['prcp','cwd','cmi'])==True:
					z0=z0/30.5
			except:
				pass
			for iU in range(uEZ.size):
				cd=lut_n2s(meta['LUT']['EZ'],uEZ[iU])[0]
				ind=np.where( (zEZf==uEZ[iU]) )
				d[cd]['ann'][v][iT]=np.mean(z0[ind].astype('float')*meta['SF'][v])

	# Save summary time series by ecozone
	gu.opickle(meta['Paths']['na1k'] + '\\Stats\\SummaryTimeSeries_ByEZ.pkl',d)

	# Time trend stats
	sts={'Region':np.array(list(d.keys()),dtype='object'),
	  'tmean_ann':np.zeros(uEZ.size),
	  'tmean_mjjas':np.zeros(uEZ.size),
	  'prcp_mjjas':np.zeros(uEZ.size),
	  'cmi_mjjas':np.zeros(uEZ.size),
	  'cwd_mjjas':np.zeros(uEZ.size),
	  'ws_mjjas':np.zeros(uEZ.size)}
	sea='mjjas'
	cnt=0
	for r in d.keys():
		for v in d[r][sea].keys():
			iT=np.where( (meta['tva']>=1951) & (meta['tva']<=2023) & (np.isnan(d[r][sea][v])==False) )[0]
			x=meta['tva'][iT]
			y=d[r][sea][v][iT]
			rs,txt=gu.GetRegStats(x,y)
			sts[v + '_' + sea][cnt]=np.round(rs['B'][1]*rs['N'],decimals=1)
		cnt=cnt+1
	
	sea='ann'
	cnt=0
	for r in d.keys():
		for v in d[r][sea].keys():
			iT=np.where( (meta['tva']>=1951) & (meta['tva']<=2023) & (np.isnan(d[r][sea][v])==False) )[0]
			x=meta['tva'][iT]
			y=d[r][sea][v][iT]
			rs,txt=gu.GetRegStats(x,y)
			sts[v + '_' + sea][cnt]=np.round(rs['B'][1]*rs['N'],decimals=1)
		cnt=cnt+1

	df=pd.DataFrame(sts)
	df.to_excel(meta['Paths']['na1k'] + '\\Stats\\SummaryTrends_ByEZ.xlsx')

	#plt.close('all'); plt.plot(meta['tva'],d['prcp'][:,2],'-bo')
	return

#%%
def lut_n2s(dc,numb):
	if numb!=-999:
		vals=np.fromiter(dc.values(),dtype=float)
		keys=np.fromiter(dc.keys(),dtype='<U70')
		ind=np.where(vals==numb)[0]
		s=keys[ind]
	else:
		s=np.array(['Unidentified'],ndmin=1)
	return s

#%%
def QA_ClimateNA(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	mo=7

	# Station data
	dS=gu.ipickle(meta['Paths']['ECCD'] + '\\Monthly\\GapFilled\\Data Ref Stations Gap Filled.pkl')
	idx=gis.GetGridIndexToPoints(zRef,dS['X'],dS['Y'])
	x0=dS['Data']['tmean']
	x0[x0<-30]=np.nan
	iT=np.where( (meta['tvm'][:,0]>=1971) & (meta['tvm'][:,0]<=2000) & (meta['tvm'][:,1]==mo) )[0]
	x=np.nanmean(x0[iT,:],axis=0)

	# ClimateNA
	y=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo) + '.tif')['Data'][idx].astype('float')*meta['SF']['tmean']
	#z=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['SF']['tmean']
	#plt.matshow(z,clim=[-25,20])

	ind=np.where(y<-60)[0]
	plt.close('all')
	plt.plot(dS['X'],dS['Y'],'k.')
	plt.plot(zRef['X'][idx][ind],zRef['Y'][idx][ind],'go')

	ikp=np.where( (zRef['Data'][idx]==1) & (np.isnan(x)==False) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
	ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot(x[ikp],y[ikp],'o',ms=3,mec='w',mfc='k',mew=0.5)
	rs,txt=gu.GetRegStats(x[ikp],y[ikp])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
	ax.text(18.5,-27,rs['txt'],fontsize=7,ha='right')
	ax.text(17,17,'1:1',fontsize=7,ha='center')
	ax.set(xlabel='Station observations',ylabel='ClimateNA',xticks=np.arange(-30,30,5),yticks=np.arange(-30,30,5),xlim=[-30,30],ylim=[-30,30])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Name','png',900)
	return

#%% Import GSOD
def ImportGSOD(meta):
	n=50000
	cnt=0
	dS={}
	dS['ID']=np.array(['' for _ in range(n)],dtype=object)
	dS['Name']=np.array(['' for _ in range(n)],dtype=object)
	dS['Lat']=np.zeros(n)
	dS['Lon']=np.zeros(n)

	dD={}
	dD['tmean']=np.nan*np.zeros((meta['tvm'].shape[0],n))
	dD['tmin']=np.nan*np.zeros((meta['tvm'].shape[0],n))
	dD['tmax']=np.nan*np.zeros((meta['tvm'].shape[0],n))
	dD['dpt']=np.nan*np.zeros((meta['tvm'].shape[0],n))

	for yr in range(2017,2020):
		print(yr)
		file=tarfile.open(meta['Paths']['GSOD'] + '\\Given\\Raw\\' + str(yr) + '.tar.gz')
		file.extractall(meta['Paths']['GSOD'] + '\\Given\\Raw\\temp')
		file.close()

		fL=os.listdir(meta['Paths']['GSOD'] + '\\Given\\Raw\\Temp')
		for f in fL:
			#break
			d=gu.ReadCSV(meta['Paths']['GSOD'] + '\\Given\\Raw\\Temp\\' + f)

			if (d['LATITUDE'][0]<24) | (d['LONGITUDE'][0]<-168) | (d['LONGITUDE'][0]>-50):
				os.remove(meta['Paths']['GSOD'] + '\\Given\\Raw\\Temp\\' + f)
				continue

			if np.isin(d['STATION'][0],dS['ID'])==False:
				dS['ID'][cnt]=d['STATION'][0]
				dS['Name'][cnt]=d['NAME'][0]
				dS['Lat'][cnt]=d['LATITUDE'][0]
				dS['Lon'][cnt]=d['LONGITUDE'][0]
				cnt=cnt+1

			indS=np.where(dS['ID']==d['STATION'][0])[0][0]
			if indS.size==0:
				os.remove(meta['Paths']['GSOD'] + '\\Given\\Raw\\Temp\\' + f)
				continue
			#print(indS)

			Month=np.zeros(d['DATE'].size)
			for i in range(d['DATE'].size):
				Month[i]=d['DATE'][i][5:7]

			for mo in range(12):
				iT0=np.where( (meta['tvm'][:,0]==yr) & (meta['tvm'][:,1]==mo+1) )[0]
				iT1=np.where( (Month==mo+1) )[0]
				if iT1.size<28:
					continue
				dD['tmean'][iT0,indS]=(np.mean(d['TEMP'][iT1])-32)*0.5555
				dD['dpt'][iT0,indS]=(np.mean(d['DEWP'][iT1])-32)*0.5555
				dD['tmax'][iT0,indS]=(np.mean(d['MAX'][iT1])-32)*0.5555
				dD['tmin'][iT0,indS]=(np.mean(d['MIN'][iT1])-32)*0.5555
			os.remove(meta['Paths']['GSOD'] + '\\Given\\Raw\\Temp\\' + f)

	# Truncate
	for k in dS.keys():
		dS[k]=dS[k][0:cnt]
	for k in dD.keys():
		dD[k]=dD[k][:,0:cnt]

	# Import spatial reference systems
	srs=gis.ImportSRSs()
	dS['X'],dS['Y']=srs['Proj']['NACID'](dS['Lon'],dS['Lat'])

	dD['dpt'][dD['dpt']>65]=np.nan

	dD['es']=gaia.GetEstar(dD['tmean']) #hPa
	dD['ea']=gaia.GetEstar(dD['dpt']) #hPa
	dD['vpd']=np.maximum(0,dD['es']-dD['ea'])

	# Save
	dS['Data']=dD
	gu.opickle(meta['Paths']['GSOD'] + '\\Processed\\gsod.pkl',dS)
	return

#%%
def CalcAnnualAnomalies(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	vL=['tmean']
	for v in vL:
		#zN=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_' + v + '_norm_1971to2000_ann.tif')['Data'].astype('float')*meta['SF'][v]
		for iT in range(meta['tva'].size):
			if meta['tva'][iT]>2023:
				continue
			zA=np.zeros(zRef['Data'].shape)
			for iM in range(12):
				zA=zA+gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\na1k_' + v + '_anom_' + str(meta['tva'][iT]) + '_' + str(iM+1) + '.tif')['Data'].astype('float')*meta['SF'][v]
			zA=zA/12
			z=copy.deepcopy(zRef)
			z['Data']=zA/meta['SF'][v]
			z['Data']=z['Data'].astype('int16')
			gis.SaveGeoTiff(z,meta['Paths']['na1k'] + '\\Summaries\\Anomalies\\na1k_' + v + '_anom_ann_' + str(meta['tva'][iT]) + '.tif')
	return

#%%
def HoldridgeLifeZones(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['na1k Ref Grid'])
	iMask=np.where( (zRef['Data']>0) )

	zMAT=np.zeros(iMask[0].size)
	for mo in range(12):
		a=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_tmean_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype(float)*meta['Climate']['SF']['tmean']
		ind=np.where(a<0)
		a[ind]=0
		zMAT=zMAT+a
	zMAT=zMAT/12

	zMAP=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_prcp_norm_1971to2000_ann.tif')['Data'][iMask].astype(float)
	zPET=gis.OpenGeoTiff(meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_etp_norm_1971to2000_ann.tif')['Data'][iMask].astype(float)
	zAI=gu.Clamp(zPET/zMAP,0,1000)
	lut=gu.ReadExcel(r'G:\My Drive\Code_Python\fcgadgets\cbrunner\Parameters\LUT_HoldgridgeLifeZone.xlsx',sheet_name='Sheet1',skiprows=0)

	zT,muT,sigT=gu.zscore(zMAT)
	zP,muP,sigP=gu.zscore(zMAP)
	zA,muA,sigA=gu.zscore(zAI)

	lut['zT']=(lut['MAT']-muT)/sigT
	lut['zP']=(lut['MAP']-muP)/sigP
	lut['zA']=(lut['AI']-muA)/sigA

	z1=np.zeros(iMask[0].size)
	e1=100000*np.ones(iMask[0].size)
	for i in range(lut['Name'].size):
		print(lut['Name'][i])
		#if lut['Name'][i]=='Wet Tundra':
		#	continue
		#eMAT=np.abs(zMAT-lut['MAT'][i])/np.mean(zMAT)
		#eMAP=np.abs(zMAP-lut['MAP'][i])/np.mean(zMAP)
		#eAI=np.abs(zAI-lut['AI'][i])/np.mean(zAI)
		eT=np.abs(zT-lut['zT'][i])
		eP=np.abs(zP-lut['zP'][i])
		eA=np.abs(zA-lut['zA'][i])
		eTot=eT+eP+eA
		ind=np.where(eTot<e1)
		z1[ind]=lut['ID'][i]
		e1[ind]=eTot[ind]

	z2=copy.deepcopy(zRef)
	z2['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	z2['Data'][iMask]=z1
	d=gu.CountByCategories(z2['Data'][iMask].flatten(),'Percent')
	# plt.close('all');plt.matshow(z2['Data'],clim=[0,15])
	gis.SaveGeoTiff(z2,meta['Paths']['na1k'] + '\\Summaries\\Normals\\na1k_HoldridgeLifeZones_1971to2000.tif')

	return

#%%
