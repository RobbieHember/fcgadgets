#%% Import modules
import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
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
import fcgadgets.na1k.na1k_util as u1k

#%% Create reference grid
def CreateRefGrid_BC5k(meta):
	# Manually copy the Land Mask from the bc1ha DB and then rescale it to 5k
	fin=meta['Paths']['bc5k Ref Grid']
	gis.ResampleRaster(fin,0.02)
	#Check that it worked: z=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	return

#%% Import CMIP6 data
def ImportCMIP6(meta):
	scnL=['ssp245','ssp585']
	for scn in scnL:
		t0=time.time()
		zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
		iLand=np.where(zRef['Data']==1)
		srs=gis.ImportSRSs()
		lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
		tvm=gu.tvec('m',1950,2101)
		tvH0=gu.tvec('m',1950,2014)
		tvF0=gu.tvec('m',2015,2100)
		press=1013.25
		con=gaia.HydroMetCon()
		vL=['tmean','ea','prcp','vpd','rswd']
		mdL=['CESM2','CNRM-CM6-1','GFDL-ESM4','HadGEM3-GC31-LL','MIROC-ES2L','MPI-ESM1-2-LR']

		# Import monthly data for each model
		dC={}
		for md in mdL:
			print(md)
			flL=os.listdir(meta['Paths']['DB']['CMIP6'] + '\\' + md)
			def GetFileName(flL,v,period):
				fl1=[]
				for fl in flL:
					if fl.startswith(v):
						if period in fl:
							fl1=fl
				return fl1
	
			# Initialize
			dC[md]={}
			dC[md]={}
			for v in vL:
				dC[md][v]=np.zeros((tvm.shape[0],zRef['m'],zRef['n']))
	
			# Temperature
			fl=GetFileName(flL,'tas','historical')
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['tas']
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
				dC[md]['tmean'][indT[0],:,:]=z1-273.15
			fl=GetFileName(flL,'tas',scn)
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['tas']
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
				dC[md]['tmean'][indT[0],:,:]=z1-273.15
			#plt.plot(dC[md]['tmean'][:,10,10])
	
			# Humidity
			fl=GetFileName(flL,'huss','historical')
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['huss']*press/(0.378*d0['huss']+0.622) # Convert to hPa
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
				dC[md]['ea'][indT[0],:,:]=z1
			fl=GetFileName(flL,'huss',scn)
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['huss']*press/(0.378*d0['huss']+0.622) # Convert to hPa
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
				dC[md]['ea'][indT[0],:,:]=z1
			#plt.plot(dC[md]['ea'][:,10,10])

			# Radiation
			fl=GetFileName(flL,'rsds','historical')
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['rsds']*86400/1e6 # Converted from W/m2 to MJ/m2/d-1
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
				dC[md]['rswd'][indT[0],:,:]=z1
			fl=GetFileName(flL,'rsds',scn)
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			z0=d0['rsds']*86400/1e6 # Converted from W/m2 to MJ/m2/d-1
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
				dC[md]['rswd'][indT[0],:,:]=z1
			#plt.plot(dC[md]['rswd'][:,10,10])

			# Precipitation
			fl=GetFileName(flL,'pr','historical')
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			# Convert from kg m-2 s-1 to mm/mo
			z0=d0['pr']
			for mo in range(12):
				z0[mo::12,:,:]=86400*z0[mo::12,:,:]*con['DIM'][mo]
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvH0[iT,0]) & (tvm[:,1]==tvH0[iT,1]) )[0]
				dC[md]['prcp'][indT[0],:,:]=z1
			fl=GetFileName(flL,'pr',scn)
			d0=gu.ReadNC(meta['Paths']['DB']['CMIP6'] + '\\' + md + '\\' + fl)
			lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
			ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
			indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
			lat1=lat0[indS]; lon1=lon0[indS]
			# Convert from kg m-2 s-1 to mm/mo
			z0=d0['pr']
			for mo in range(12):
				z0[mo::12,:,:]=86400*z0[mo::12,:,:]*con['DIM'][mo]
			for iT in range(z0.shape[0]):
				z1=z0[iT,:,:].T
				z1=z1[indS]
				z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
				indT=np.where( (tvm[:,0]==tvF0[iT,0]) & (tvm[:,1]==tvF0[iT,1]) )[0]
				dC[md]['prcp'][indT[0],:,:]=z1
			#plt.plot(dC[md]['rswd'][:,10,10])

		# Add vpd
		for md in mdL:
			dC[md]['vpd']=np.maximum(0,gaia.GetEstar(dC[md]['tmean'])-dC[md]['ea'])

		# Add water balance
		for md in mdL:
			vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
			for v in vwbL:
				dC[md][v]=np.zeros((tvm.shape[0],zRef['m'],zRef['n']))
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
			vi={}
			for iT in range(tvm.shape[0]):
				print('Year:' + str(tvm[iT,0]) + ', Month:' + str(tvm[iT,1]) )
				mo=tvm[iT,1]
				vi['Month']=mo
				vi['LAI']=5.0
				vi['Gs']=0.010
				vi['Ga']=0.058
				vi['tmean']=dC[md]['tmean'][iT,:,:][iLand]
				vi['prcp']=dC[md]['prcp'][iT,:,:][iLand]
				vi['vpd']=dC[md]['vpd'][iT,:,:][iLand]
				vi['rswn']=dC[md]['rswd'][iT,:,:][iLand]
				vo=gaia.WBM_SparseGrid(par,vi)
				vi['ws']=vo['ws']
				vi['wsp']=vo['wsp']
				vo['cwd']=vo['etp']-vo['eta']
				for v in vwbL:
					dC[md][v][iT,:,:][iLand]=vo[v]
		#plt.plot(dC[md]['ws'][:,150,150])

		# Calculate seasonal indicators
		tva=np.arange(1950,2101,1)
		iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
		dCS={}
		for md in mdL:
			dCS[md]={}
			dCS[md]={}
			dCS[md]['tmean']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),
					 'wyr':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			dCS[md]['ws']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			dCS[md]['cwd']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			dCS[md]['cmi']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			dCS[md]['prcp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			dCS[md]['etp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
			for iY in range(tva.size):
				iT=np.where( (tvm[:,0]==tva[iY]) )[0]
				dCS[md]['tmean']['ann'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)

				try:
					iT=np.where( (tvm[:,0]==tva[iY-1]) & (tvm[:,1]>=10) & (tvm[:,0]==tva[iY]) & (tvm[:,1]<=9) )[0]
				except:
					# It will crash in year 0, so just use full year, should not affect results that much
					iT=np.where( (tvm[:,0]==tva[iY]) )[0]
				dCS[md]['tmean']['wyr'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)

				iT=np.where( (tvm[:,0]==tva[iY]) & (tvm[:,1]>=5) & (tvm[:,1]<=9) )[0]
				dCS[md]['ws']['mjjas'][iY]=np.mean(dC[md]['ws'][iT,:,:],axis=0)
				dCS[md]['cwd']['mjjas'][iY]=np.mean(dC[md]['cwd'][iT,:,:],axis=0)
				dCS[md]['cmi']['mjjas'][iY]=np.mean(dC[md]['prcp'][iT,:,:]-dC[md]['etp'][iT,:,:],axis=0)
				dCS[md]['prcp']['mjjas'][iY]=np.mean(dC[md]['prcp'][iT,:,:],axis=0)
				dCS[md]['etp']['mjjas'][iY]=np.mean(dC[md]['etp'][iT,:,:],axis=0)

		# Convert to anomalies
		dCSa=copy.deepcopy(dCS)
		for md in mdL:
			for k in dCSa[md].keys():
				for j in dCSa[md][k].keys():
					mu=np.mean(dCSa[md][k][j][iTmu,:,:],axis=0)
					for iY in range(tva.size):
						dCSa[md][k][j][iY,:,:]=dCSa[md][k][j][iY,:,:]-mu
	
		# Compress
		for md in mdL:
			for k in dCS[md].keys():
				for j in dCS[md][k].keys():
					dCS[md][k][j]=dCS[md][k][j]/meta['Climate']['SF'][k]
					dCS[md][k][j]=dCS[md][k][j].astype('int16')
					dCSa[md][k][j]=dCSa[md][k][j]/meta['Climate']['SF'][k]
					dCSa[md][k][j]=dCSa[md][k][j].astype('int16')

		# Save
		gu.opickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_actuals_' + scn + '.pkl',dCS)
		gu.opickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_anomalies_' + scn + '.pkl',dCSa)
		print((time.time()-t0)/60)
	return

#%% Import 20thC reanlysis data
def Import20thCenturyReanalysis(meta):
	t0=time.time()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iLand=np.where(zRef['Data']==1)

	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	press=1013.25
	con=gaia.HydroMetCon()
	vL=['tmean','ea','prcp','vpd','rswd']

	dC={}
	dC['tv']=gu.tvec('m',1871,2012)
	for v in vL:
		dC[v]=np.zeros( (dC['tv'].shape[0],zRef['m'],zRef['n']) )

	# Temperature
	d0=gu.ReadNC(meta['Paths']['DB']['20thCR'] + '\\air.2m.mon.mean.nc')
	lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
	ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
	indS=np.where( (lat0>40) & (lat0<65) & (lon0>-140) & (lon0<-100) )
	lat1=lat0[indS]; lon1=lon0[indS]
	for iT in range(dC['tv'].shape[0]):
		z1=d0['air'][iT,:,:].T
		z1=z1[indS]
		z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		dC['tmean'][iT,:,:]=z1-273.15
	#plt.close('all'); plt.matshow(np.mean(dC['tmean'],axis=0),clim=[0,15]); plt.colorbar()

	# Humidity
	d0=gu.ReadNC(meta['Paths']['DB']['20thCR'] + '\\shum.2m.mon.mean.nc')
	for iT in range(dC['tv'].shape[0]):
		z1=d0['shum'][iT,:,:].T
		z1=z1[indS]
		z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		#dC['ea'][iT,:,:]=z1*1000/con['RhoAir']*461.5*(dC['tmean'][iT,:,:]+273.15)
		#eaR=shumR/rhoAir*461.5*(tmeanR+273.15)
		dC['ea'][iT,:,:]=z1*press/(0.378*z1+0.622) # Convert to hPa
	#plt.close('all'); plt.matshow(np.mean(dC['ea'],axis=0),clim=[0,15]); plt.colorbar()

	# Radiation
	d0=gu.ReadNC(meta['Paths']['DB']['20thCR'] + '\\dswrf.sfc.mon.mean.nc')
	for iT in range(dC['tv'].shape[0]):
		z1=d0['dswrf'][iT,:,:].T
		z1=z1[indS]
		z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		dC['rswd'][iT,:,:]=z1*86400/1e6 # Converted from W/m2 to MJ/m2/d-1
	#plt.close('all'); plt.matshow(np.mean(dC['rswd'],axis=0),clim=[0,15]); plt.colorbar()

	# Precipitation
	d0=gu.ReadNC(meta['Paths']['DB']['20thCR'] + '\\prate.mon.mean.nc')
	for iT in range(dC['tv'].shape[0]):
		z1=d0['prate'][iT,:,:].T
		z1=z1[indS]
		z1=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		dim=con['DIM'][dC['tv'][iT,1]-1]
		dC['prcp'][iT,:,:]=z1*86400*dim # Convert from kg m-2 s-1 to mm/mo
	#plt.close('all'); plt.matshow(np.mean(dC['prcp'],axis=0),clim=[0,150]); plt.colorbar()

	# Add vpd
	dC['vpd']=np.maximum(0,gaia.GetEstar(dC['tmean'])-dC['ea'])

	# Add water balance
	vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
	for v in vwbL:
		dC[v]=np.zeros( (dC['tv'].shape[0],zRef['m'],zRef['n']) )
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
	vi={}
	for iT in range(dC['tv'].shape[0]):
		print('Year:' + str(dC['tv'][iT,0]) + ', Month:' + str(dC['tv'][iT,1]) )
		mo=dC['tv'][iT,1]
		vi['Month']=mo
		vi['LAI']=5.0
		vi['Gs']=0.010
		vi['Ga']=0.058
		vi['tmean']=dC['tmean'][iT,:,:][iLand]
		vi['prcp']=dC['prcp'][iT,:,:][iLand]
		vi['vpd']=dC['vpd'][iT,:,:][iLand]
		vi['rswn']=dC['rswd'][iT,:,:][iLand]
		vo=gaia.WBM_SparseGrid(par,vi)
		vi['ws']=vo['ws']
		vi['wsp']=vo['wsp']
		vo['cwd']=vo['etp']-vo['eta']
		for v in vwbL:
			dC[v][iT,:,:][iLand]=vo[v]
	#plt.plot(dC[md]['ws'][:,150,150])

	# Calculate seasonal indicators
	tva=dC['tv'][0::12,0]
	iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
	dCS={}

	dCS={}
	dCS={}
	dCS['tmean']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),'wyr':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['ws']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['cwd']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['cmi']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['prcp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['etp']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	for iY in range(tva.size):
		iT=np.where( (dC['tv'][:,0]==tva[iY]) )[0]
		dCS['tmean']['ann'][iY]=np.mean(dC['tmean'][iT,:,:],axis=0)
		iT=np.where( (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]>=5) & (dC['tv'][:,1]<=9) )[0]
		dCS['ws']['mjjas'][iY]=np.mean(dC['ws'][iT,:,:],axis=0)
		dCS['cwd']['mjjas'][iY]=np.mean(dC['cwd'][iT,:,:],axis=0)
		dCS['cmi']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:]-dC['etp'][iT,:,:],axis=0)
		dCS['prcp']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:],axis=0)
		dCS['etp']['mjjas'][iY]=np.mean(dC['etp'][iT,:,:],axis=0)
		iT=np.where( (dC['tv'][:,0]==tva[iY]-1) & (dC['tv'][:,1]>=10) | (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]<=9) )[0]
		dCS['tmean']['wyr'][iY]=np.mean(dC['tmean'][iT,:,:],axis=0)

	# Convert to anomalies
	dCSa=copy.deepcopy(dCS)
	for k in dCSa.keys():
		for j in dCSa[k].keys():
			mu=np.mean(dCSa[k][j][iTmu,:,:],axis=0)
			for iY in range(tva.size):
				dCSa[k][j][iY,:,:]=dCSa[k][j][iY,:,:]-mu

	# Compress
	for k in dCS.keys():
		for j in dCS[k].keys():
			dCS[k][j]=dCS[k][j]/meta['Climate']['SF'][k]
			dCS[k][j]=dCS[k][j].astype('int16')
			dCSa[k][j]=dCSa[k][j]/meta['Climate']['SF'][k]
			dCSa[k][j]=dCSa[k][j].astype('int16')

	# Save
	gu.opickle(meta['Paths']['DB']['20thCR'] + '\\20thCR_actuals.pkl',dCS)
	gu.opickle(meta['Paths']['DB']['20thCR'] + '\\20thCR_anomalies.pkl',dCSa)
	print((time.time()-t0)/60)
	return

#%% Import CRU data
def ImportCRU(meta):
	t0=time.time()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iLand=np.where(zRef['Data']==1)

	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	press=1013.25
	con=gaia.HydroMetCon()

	# Import CRU data
	tvCRU=gu.tvec('m',1901,2022)
	tvCRUa=np.arange(1901,2023)

	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.tmp.dat.nc.gz'
	tmean0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					tmean0[k]=ds.variables[k][:]
	tmean1=tmean0['tmp'].filled() # (given in numpy mask arrays so you need to use the filled method)
	tmean1[tmean1>50]=np.nan

	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.pre.dat.nc.gz'
	prcp0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					prcp0[k]=ds.variables[k][:]
	prcp1=prcp0['pre'].filled()
	prcp1[prcp1>1500]=np.nan

	fin=r'C:\Data\Climate\CRU\cru_ts4.07.1901.2022.vap.dat.nc.gz'
	ea0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					ea0[k]=ds.variables[k][:]
	ea1=ea0['vap'].filled()
	ea1[ea1>100]=np.nan

	# Saturation vapour pressure
	esat1=gaia.GetEstar(tmean1)

	# Vapour pressure deficit
	vpd1=np.maximum(0,esat1-ea1)

	# Isolate BC area  (so that interpolation is faster)
	lon,lat=np.meshgrid(tmean0['lon'].filled(),tmean0['lat'].filled(),sparse=False)
	x1,y1=srs['Proj']['BC1ha'](lon,lat)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))

	# Interpolate to 5 km grid
	vL=['tmean','prcp','vpd','rswd']
	dC={}
	dC['tv']=gu.tvec('m',1851,2022)
	for v in vL:
		dC[v]=np.nan*np.ones( (dC['tv'].shape[0],zRef['m'],zRef['n']) )

	for iT in range(tvCRU.shape[0]):
		ind=np.where( (dC['tv'][:,0]==tvCRU[iT,0]) & (dC['tv'][:,1]==tvCRU[iT,1]) )
		if ind[0].size==0: continue
		dC['tmean'][ind,:,:]=griddata(xy0,tmean1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')
		dC['prcp'][ind,:,:]=griddata(xy0,prcp1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')
		dC['vpd'][ind,:,:]=griddata(xy0,vpd1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')

	# Radiation
	d0=gu.ReadNC(meta['Paths']['DB']['20thCR'] + '\\dswrf.sfc.mon.mean.nc')
	tvR=gu.tvec('m',1871,2012)
	lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
	ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
	x1,y1=srs['Proj']['BC1ha'](lon0,lat0)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))
	for iT in range(tvR.shape[0]):
		z1=d0['dswrf'][iT,:,:].T
		z1=griddata(xy0,z1[indS],(zRef['X'],zRef['Y']),method='linear')
		ind=np.where( (dC['tv'][:,0]==tvR[iT,0]) & (dC['tv'][:,1]==tvR[iT,1]) )[0]
		dC['rswd'][ind,:,:]=z1*86400/1e6 # Converted from W/m2 to MJ/m2/d-1

	d0=gu.ReadNC(r'C:\Data\Reanalysis\NCEPGlobal\Monthly\dswrf.sfc.mon.mean.nc')
	tvR=gu.tvec('m',1948,2024)
	lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
	ind=np.where(lon0>180)[0]; lon0[ind]=lon0[ind]-360
	x1,y1=srs['Proj']['BC1ha'](lon0,lat0)
	bw=100000
	indS=np.where( (x1>=zRef['xmin']-bw) & (x1<=zRef['xmax']+bw) & (y1>=zRef['ymin']-bw) & (y1<=zRef['ymax']+bw) )
	xy0=np.column_stack((x1[indS],y1[indS]))
	for iT in range(tvR.shape[0]):
		try:
			z1=d0['dswrf'][iT,:,:].T
		except:
			continue
		z1=griddata(xy0,z1[indS],(zRef['X'],zRef['Y']),method='linear')
		ind=np.where( (dC['tv'][:,0]==tvR[iT,0]) & (dC['tv'][:,1]==tvR[iT,1]) )[0]
		if ind.size==0: continue
		dC['rswd'][ind,:,:]=z1*86400/1e6 # Converted from W/m2 to MJ/m2/d-1

	#plt.plot(np.nanmean(dC['rswd'],axis=(1,2)))

	# Populate early record with the lowest temperature percentile data
	mu=np.nanmean(tmean1,axis=(1,2))
	dCtva=np.arange(1851,2023)
	for mo in range(12):
		p=np.nanpercentile(mu[mo::12],15)
		ind=np.where(mu[mo::12]<p)
		yr0=tvCRUa[ind]
		yr0=np.random.choice(yr0,50)
		for iY in range(50):
			ind1=np.where( (dC['tv'][:,0]==dCtva[iY]) & (dC['tv'][:,1]==mo+1) )[0]
			ind2=np.where( (dC['tv'][:,0]==yr0[iY]) & (dC['tv'][:,1]==mo+1) )[0]
			dC['tmean'][ind1,:,:]=dC['tmean'][ind2,:,:]
			dC['prcp'][ind1,:,:]=dC['prcp'][ind2,:,:]
			dC['vpd'][ind1,:,:]=dC['vpd'][ind2,:,:]
			dC['rswd'][ind1,:,:]=dC['rswd'][ind2,:,:]

	# Add water balance
	vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
	for v in vwbL:
		dC[v]=np.zeros( (dC['tv'].shape[0],zRef['m'],zRef['n']) )
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
	vi={}
	for iT in range(dC['tv'].shape[0]):
		print('Year:' + str(dC['tv'][iT,0]) + ', Month:' + str(dC['tv'][iT,1]) )
		vi['Month']=dC['tv'][iT,1]
		vi['LAI']=5.0
		vi['Gs']=0.010
		vi['Ga']=0.058
		vi['tmean']=np.nan_to_num(dC['tmean'][iT,:,:][iLand])
		vi['prcp']=np.nan_to_num(dC['prcp'][iT,:,:][iLand])
		vi['vpd']=np.nan_to_num(dC['vpd'][iT,:,:][iLand])
		vi['rswn']=np.nan_to_num(dC['rswd'][iT,:,:][iLand])
		vo=gaia.WBM_SparseGrid(par,vi)
		vi['ws']=vo['ws']
		vi['wsp']=vo['wsp']
		vo['cwd']=vo['etp']-vo['eta']
		for v in vwbL:
			dC[v][iT,:,:][iLand]=vo[v]

	flg=0
	if flg==1:
		plt.plot(np.nanmean(dC['tmean'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['prcp'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['etp'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['eta'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['rswd'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['vpd'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['ws'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['wsp'],axis=(1,2)))
		plt.plot(np.nanmean(dC['runoff'],axis=(1,2)))
		plt.plot(np.nanmean(dC['cwd'],axis=(1,2)))

	# Calculate seasonal indicators
	tva=dC['tv'][0::12,0]
	iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
	dCS={}
	dCS['tmean']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),
		'wyr':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['prcp']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),
		'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['eta']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),
		'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['etp']={'ann':np.zeros((tva.shape[0],zRef['m'],zRef['n'])),
		'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['cmi']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['ws']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	dCS['cwd']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	for iY in range(tva.size):
		iT=np.where( (dC['tv'][:,0]==tva[iY]) )[0]
		for v in ['tmean','prcp','etp','eta']:
			dCS[v]['ann'][iY]=np.mean(dC[v][iT,:,:],axis=0)

		iT=np.where( (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]>=5) & (dC['tv'][:,1]<=9) )[0]
		dCS['ws']['mjjas'][iY]=np.mean(dC['ws'][iT,:,:],axis=0)
		dCS['cwd']['mjjas'][iY]=np.mean(dC['cwd'][iT,:,:],axis=0)
		dCS['cmi']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:]-dC['etp'][iT,:,:],axis=0)
		dCS['prcp']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:],axis=0)
		dCS['eta']['mjjas'][iY]=np.mean(dC['eta'][iT,:,:],axis=0)
		dCS['etp']['mjjas'][iY]=np.mean(dC['etp'][iT,:,:],axis=0)

		iT=np.where( (dC['tv'][:,0]==tva[iY]-1) & (dC['tv'][:,1]>=10) | (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]<=9) )[0]
		dCS['tmean']['wyr'][iY]=np.mean(dC['tmean'][iT,:,:],axis=0)

	# Convert to anomalies
	dCSa=copy.deepcopy(dCS)
	for k in dCSa.keys():
		for j in dCSa[k].keys():
			mu=np.mean(dCSa[k][j][iTmu,:,:],axis=0)
			for iY in range(tva.size):
				dCSa[k][j][iY,:,:]=dCSa[k][j][iY,:,:]-mu

	# Compress
	for k in dCS.keys():
		for j in dCS[k].keys():
			dCS[k][j]=dCS[k][j]/meta['Climate']['SF'][k]
			dCS[k][j]=dCS[k][j].astype('int16')
			dCSa[k][j]=dCSa[k][j]/meta['Climate']['SF'][k]
			dCSa[k][j]=dCSa[k][j].astype('int16')

	# Save
	gu.opickle(meta['Paths']['DB']['CRU'] + '\\cru_actuals.pkl',dCS)
	gu.opickle(meta['Paths']['DB']['CRU'] + '\\cru_anomalies.pkl',dCSa)
	print((time.time()-t0)/60)
	return

#%% Import QFED wildfire emissions
# units: kg s-1 m-2
def Process_WildfireE_QFED(meta):

	# Conversion factor kg/s/m2 to kg/s/5km
	cf1=5*1000000
	# Conversion factor kg/s/5km to kg/5km/yr
	cf2=60*60*24*365
	# Conversion factor kg/5km/yr to t/5km/yr
	cf3=1/1e9
	# GWP for CH4
	gwp_ch4=28

	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iMaskOut=np.where(zRef['Data']!=1)
	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])

	# CO2 climatology
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	z2=copy.deepcopy(zRef)
	z2['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	ym=np.zeros(12)
	for i in range(1,13):
		if i<10:
			mo='0' + str(i)
		else:
			mo=str(i)
		d0=gu.ReadNC(r'C:\Data\Wildfire\QFED\qfed2.emis_co2.061.x3600_y1800.2003-2023' + mo + 'clm.nc4')
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		lat0=lat0.T
		lon0=lon0.T

		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]

		z0=d0['biomass'][0,:,:] # transposed to match lat and lon
		z0=z0[indS]
		z0=griddata( (lon1,lat1),z0,(lon,lat),method='linear')
		z0[iMaskOut]=0
		z0=z0*cf1*cf2*cf3
		z1['Data']=z1['Data']+z0
		if (i>=10):
			z2['Data']=z2['Data']+z0
		ym[i-1]=np.sum(z0)

	# CH4 climatology
	#z1=copy.deepcopy(zRef)
	#z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	#ym=np.zeros(12)
	for i in range(1,13):
		if i<10:
			mo='0' + str(i)
		else:
			mo=str(i)
		d0=gu.ReadNC(r'C:\Data\Wildfire\QFED\qfed2.emis_ch4.061.x3600_y1800.2003-2023' + mo + 'clm.nc4')
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		lat0=lat0.T
		lon0=lon0.T

		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]

		z0=d0['biomass'][0,:,:] # transposed to match lat and lon
		z0=z0[indS]
		z0=griddata( (lon1,lat1),z0,(lon,lat),method='linear')
		z0[iMaskOut]=0
		z0=z0*cf1*cf2*cf3*gwp_ch4
		z1['Data']=z1['Data']+z0
		if (i>=10):
			z2['Data']=z2['Data']+z0
		ym[i-1]=ym[i-1]+np.sum(z0)

	plt.close('all'); plt.matshow(z1['Data'],clim=[0,0.02])
	plt.close('all'); plt.matshow(z2['Data'],clim=[0,0.02])
	plt.close('all'); plt.bar(np.arange(1,13),ym,0.75)

	# CO2 2017
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	ym=np.zeros(12)
	for i in range(1,13):
		if i<10:
			mo='0' + str(i)
		else:
			mo=str(i)
		d0=gu.ReadNC(r'C:\Data\Wildfire\QFED\qfed2.emis_co2.061.x3600_y1800.2017' + mo + 'mm.nc4')
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		lat0=lat0.T
		lon0=lon0.T

		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]

		z0=d0['biomass'][0,:,:] # transposed to match lat and lon
		z0=z0[indS]
		z0=griddata( (lon1,lat1),z0,(lon,lat),method='linear')
		z0[iMaskOut]=0
		z0=z0*cf1*cf2*cf3
		z1['Data']=z1['Data']+z0
		ym[i-1]=np.sum(z0)

	# CH4 2017
	#z1=copy.deepcopy(zRef)
	#z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
	#ym=np.zeros(12)
	for i in range(1,13):
		if i<10:
			mo='0' + str(i)
		else:
			mo=str(i)
		d0=gu.ReadNC(r'C:\Data\Wildfire\QFED\qfed2.emis_ch4.061.x3600_y1800.2017' + mo + 'mm.nc4')
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		lat0=lat0.T
		lon0=lon0.T

		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]

		z0=d0['biomass'][0,:,:] # transposed to match lat and lon
		z0=z0[indS]
		z0=griddata( (lon1,lat1),z0,(lon,lat),method='linear')
		z0[iMaskOut]=0
		z0=z0*cf1*cf2*cf3*gwp_ch4
		z1['Data']=z1['Data']+z0
		ym[i-1]=ym[i-1]+np.sum(z0)

	print(np.sum(z1['Data']))
	plt.close('all'); plt.matshow(z1['Data'],clim=[0,0.02])
	plt.close('all'); plt.bar(np.arange(1,13),ym,0.75)

	#gu.opickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_' + rcp + '.pkl',d)
	return

#%% Import nitrogen deposition from ISIMIP
def Process_Ndep_ISIMIP(meta):
	tv=np.arange(1850,2101,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	rcpL=['26','60','85']
	for rcp in rcpL:
		d={'ndep':{}}
		d['ndep']['tv']=tv
		d['ndep']['Data']=np.zeros((tv.size,zRef['m'],zRef['n']),dtype='float')
		d['ndep']['TS Mean']=np.zeros(tv.size)
		#meta['Climate']['SF']['ndep']=0.01
	
		tv0=np.arange(1861,2006,1)
		d0=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_nhx_histsoc_annual_1861_2005.nc4')
		#print(np.mean(d0['nhx'][140,:,:]))
		d1=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_noy_histsoc_annual_1861_2005.nc4')
		#print(np.mean(d1['noy'][140,:,:]))
		d0['ntot']=10*(d0['nhx']+d1['noy'])
		lat0,lon0=np.meshgrid(d0['lat'],d0['lon'])
		#print(np.mean(d0['ntot'][140,:,:]))
		#plt.close('all'); plt.matshow(d0['ntot'][140,:,:],clim=[0,8]); plt.colorbar()
	
		# Isolate region
		indS=np.where( (lat0>42) & (lat0<62) & (lon0>-140) & (lon0<-100) )
		lat1=lat0[indS]
		lon1=lon0[indS]
		iLand=np.where(zRef['Data']==1)
	
		indT=np.where(tv<tv0[0])[0]
		z0=d0['ntot'][0,:,:].T # transposed to match lat and lon
		z1=z0[indS]
		z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		for j in indT:
			d['ndep']['Data'][j,:,:]=z2
			d['ndep']['TS Mean'][j]=np.mean(z2[iLand])
	
		for iT in range(tv0.size):
			indT=np.where(tv==tv0[iT])[0]
			z0=d0['ntot'][iT,:,:].T # transposed to match lat and lon
			z1=z0[indS]
			z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			d['ndep']['Data'][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['ndep']['TS Mean'][indT]=np.mean(z2[iLand])
			#print(np.mean(z2[ind]))

		# Future
		tv0=np.arange(2006,2100,1)
		d0=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_nhx_rcp' + rcp + 'soc_monthly_2006_2099.nc4')
		d0['nhx']=gu.BlockSum(d0['nhx'],12)
		d1=gu.ReadNC(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_noy_rcp' + rcp + 'soc_monthly_2006_2099.nc4')
		d1['noy']=gu.BlockSum(d1['noy'],12)
		d0['ntot']=10*(d0['nhx']+d1['noy'])
		for iT in range(tv0.size):
			indT=np.where(tv==tv0[iT])[0]
			z0=d0['ntot'][iT,:,:].T # transposed to match lat and lon
			z1=z0[indS]
			z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			d['ndep']['Data'][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['ndep']['TS Mean'][indT]=np.mean(z2[iLand])
	
		indT=np.where(tv>tv0[-1])[0]
		z0=d0['ntot'][-1,:,:].T # transposed to match lat and lon
		z1=z0[indS]
		z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		for j in indT:
			d['ndep']['Data'][j,:,:]=z2
			d['ndep']['TS Mean'][j]=np.mean(z2[iLand])
		plt.plot(tv,np.mean(np.mean(d['ndep']['Data'],axis=2),axis=1),'k-')
		#plt.close('all'); plt.plot(tv,d['ndep']['TS Mean'],'g-.')
		#plt.matshow(np.mean(d['ntot'],axis=0))
		d['ndep']['Data']=d['ndep']['Data']/meta['Climate']['SF']['ndep']
		d['ndep']['Data']=d['ndep']['Data'].astype('int16')
		gu.opickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_' + rcp + '.pkl',d)
	return

#%%
def Process_Ndep_FromAckerman2018(meta):
	# Ackerman et al. (2018) (https://conservancy.umn.edu/handle/11299/197613)
	data=pd.read_csv(r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\inorganic_N_deposition.csv')
	
	zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

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
	
	# Clip to NA
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

	# Reproject
	pthin=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip.tif'
	pthout=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_' + nam + '_clip_projected.tif'
	gis.ReprojectGeoTiff(pthin,pthout,srs['String']['BC1ha'])
	#gis.ReprojectGeoTiff(pthin,pthout,srs['String']['NACID'])
	
	z=gis.OpenGeoTiff(pthout)
	z['X']=z['X'][:,0:-1]
	
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	im=ax.matshow(z['Data'],clim=(0,12),extent=z['Extent'])
	
	#
	zND=zTSA.copy()
	ivl=10
	zND['X']=zTSA['X'][0::ivl,0::ivl]
	zND['Y']=zTSA['Y'][0::ivl,0::ivl]
	zND['Data']=griddata( (z['X'].flatten(),z['Y'].flatten()) ,z['Data'].flatten(),(zND['X'],zND['Y']),method='cubic')
	
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,16))
	im=ax.matshow(zND['Data'],clim=(0,5),extent=z['Extent'])
	
	# Plot
	# See bc1ha_map_full for plots
	return

#%%
def BiasCorrectCMIP6_FromCRU(meta):

	# Calibration period
	t0=1951
	t1=2022
	mat=10

	# Import CRU
	dH=gu.ipickle(meta['Paths']['DB']['CRU'] + '\\cru_anomalies.pkl')
	tvH=np.arange(1901,2023)
	itH=np.where( (tvH>=t0) & (tvH<=t1) )[0]

	# Import CMIP6 anomalies
	scnL=['ssp245','ssp585']
	tvM=np.arange(1950,2101,1)
	itM=np.where( (tvM>=t0) & (tvM<=t1) )[0]
	for scn in scnL:
		dM=gu.ipickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_anomalies_' + scn + '.pkl')

		# Temperature
		muH=np.mean(dH['tmean']['wyr'][itH,:,:],axis=0)
		for m in dM.keys():
			muM=np.mean(dM[m]['tmean']['wyr'][itM,:,:],axis=0)
			D=np.squeeze(muM-muH)
			for i in range(tvM.size):
				dM[m]['tmean']['wyr'][i,:,:]=dM[m]['tmean']['wyr'][i,:,:]-D
		muH=np.mean(dH['tmean']['ann'][itH,:,:],axis=0)
		for m in dM.keys():
			muM=np.mean(dM[m]['tmean']['ann'][itM,:,:],axis=0)
			D=np.squeeze(muM-muH)
			for i in range(tvM.size):
				dM[m]['tmean']['ann'][i,:,:]=dM[m]['tmean']['ann'][i,:,:]-D

		vL=['prcp', 'etp', 'cmi', 'ws', 'cwd']
		for v in vL:
			muH=np.mean(dH[v]['mjjas'][itH,:,:],axis=0)
			for m in dM.keys():
				muM=np.mean(dM[m][v]['mjjas'][itM,:,:],axis=0)
				D=np.squeeze(muM-muH)
				for i in range(tvM.size):
					dM[m][v]['mjjas'][i,:,:]=dM[m][v]['mjjas'][i,:,:]-D

		# Save bias-corrected values
		gu.opickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_anomalies_biascorrected_' + scn + '.pkl',dM)

	return