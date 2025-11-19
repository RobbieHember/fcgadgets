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
import fcexplore.field_plots.processing.fp_util as ufp
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.na1k.na1k_util as u1k

#%% Create reference grid
def CreateRefGrid_BC5k(meta):
	# Manually copy the Land Mask from the bc1ha DB and then rescale it to 5k
	fin=meta['Paths']['bc5k Ref Grid']
	gis.ResampleRaster(fin,0.02)
	#Check that it worked: z=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	return

#%%
def ProcessCMIP6(meta):
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
		esmL=['CESM2','CNRM-CM6-1','GFDL-ESM4','HadGEM3-GC31-LL','MIROC-ES2L','MPI-ESM1-2-LR']

		# Import monthly data for each model
		dC={}
		for md in esmL:
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
		for md in esmL:
			dC[md]['vpd']=np.maximum(0,gaia.GetEstar(dC[md]['tmean'])-dC[md]['ea'])

		# Add water balance
		for md in esmL:
			#vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
			vwbL=['cwd','etp','ws']
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

		for md in esmL:
			dC[md]['cmi']=dC[md]['prcp']-dC[md]['etp']

		# Calculate seasonal indicators
		tva=np.arange(1950,2101,1)
		iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
		dCS={}
		for md in esmL:

			dCS[md]={}
			seaL=['ann','mjjas','djf','wyr']
			for v in dC[md].keys():
				dCS[md][v]={}
				for s in seaL:
					dCS[md][v][s]=np.zeros( (tva.shape[0],zRef['m'],zRef['n']) )

			for iY in range(tva.size):

				iT=np.where( (tvm[:,0]==tva[iY]) )[0]
				dCS[md]['tmean']['ann'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)
				dCS[md]['prcp']['ann'][iY]=np.sum(dC[md]['prcp'][iT,:,:],axis=0)

				try:
					iT=np.where( (tvm[:,0]==tva[iY-1]) & (tvm[:,1]>=10) | (tvm[:,0]==tva[iY]) & (tvm[:,1]<=9) )[0]
				except:
					# It will crash in year 0, so just use full year, should not affect results that much
					iT=np.where( (tvm[:,0]==tva[iY]) )[0]
				dCS[md]['tmean']['wyr'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)

				try:
					iT=np.where( (tvm[:,0]==tva[iY-1]) & (tvm[:,1]==12) | (tvm[:,0]==tva[iY]) & (tvm[:,1]<=2) )[0]
				except:
					# It will crash in year 0, so just use full year, should not affect results that much
					iT=np.where( (tvm[:,0]==tva[iY]) & (tvm[:,1]<=2) )[0]
				dCS[md]['tmean']['djf'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)

				iT=np.where( (tvm[:,0]==tva[iY]) & (tvm[:,1]>=5) & (tvm[:,1]<=9) )[0]
				dCS[md]['ws']['mjjas'][iY]=np.mean(dC[md]['ws'][iT,:,:],axis=0)
				dCS[md]['cwd']['mjjas'][iY]=np.mean(dC[md]['cwd'][iT,:,:],axis=0)
				dCS[md]['cmi']['mjjas'][iY]=np.mean(dC[md]['cmi'][iT,:,:],axis=0)
				dCS[md]['prcp']['mjjas'][iY]=np.mean(dC[md]['prcp'][iT,:,:],axis=0)
				dCS[md]['etp']['mjjas'][iY]=np.mean(dC[md]['etp'][iT,:,:],axis=0)
				dCS[md]['tmean']['mjjas'][iY]=np.mean(dC[md]['tmean'][iT,:,:],axis=0)

		# Convert to anomalies
		dCSa=copy.deepcopy(dCS)
		for md in esmL:
			for k in dCSa[md].keys():
				for j in dCSa[md][k].keys():
					mu=np.mean(dCSa[md][k][j][iTmu,:,:],axis=0)
					for iY in range(tva.size):
						dCSa[md][k][j][iY,:,:]=dCSa[md][k][j][iY,:,:]-mu
	
		# Compress
		for md in esmL:
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

#%%
def Process20thCenturyReanalysis(meta):
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
	tva=gu.BlockMean(dC['tv'][:,0],12)
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
	dCS['eta']={'mjjas':np.zeros((tva.shape[0],zRef['m'],zRef['n']))}
	for iY in range(tva.size):
		iT=np.where( (dC['tv'][:,0]==tva[iY]) )[0]
		dCS['tmean']['ann'][iY]=np.mean(dC['tmean'][iT,:,:],axis=0)
		iT=np.where( (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]>=5) & (dC['tv'][:,1]<=9) )[0]
		dCS['ws']['mjjas'][iY]=np.mean(dC['ws'][iT,:,:],axis=0)
		dCS['cwd']['mjjas'][iY]=np.mean(dC['cwd'][iT,:,:],axis=0)
		dCS['cmi']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:]-dC['etp'][iT,:,:],axis=0)
		dCS['prcp']['mjjas'][iY]=np.mean(dC['prcp'][iT,:,:],axis=0)
		dCS['etp']['mjjas'][iY]=np.mean(dC['etp'][iT,:,:],axis=0)
		dCS['eta']['mjjas'][iY]=np.mean(dC['eta'][iT,:,:],axis=0)
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
	a={}
	a['Year']=tva
	a['Data']=dCS
	gu.opickle(meta['Paths']['bc5k'] + '\\20thCR_actuals.pkl',a)
	a={}
	a['Year']=tva
	a['Data']=dCSa
	gu.opickle(meta['Paths']['bc5k'] + '\\20thCR_anomalies.pkl',a)
	print((time.time()-t0)/60)
	return

#%%
def ProcessCRU(meta):
	t0=time.time()

	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	iLand=np.where(zRef['Data']==1)

	metaNA=u1k.Init()

	srs=gis.ImportSRSs()
	lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'],zRef['Y'])
	press=1013.25
	con=gaia.HydroMetCon()

	# Import CRU data

	# Import netcdf data
	# C:\Data\Climate\CRU\cru_ts4.09.1901.2024.tmp.dat.gz
	tvCRU=gu.tvec('m',1901,metaNA['YearLast'])
	tvCRUa=np.arange(1901,metaNA['YearLast']+1)

	fin=metaNA['Paths']['CRU'] + '\\' + metaNA['Current Versions']['CRU'] + '.tmp.dat.nc.gz'
	tmean0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					tmean0[k]=ds.variables[k][:]
	tmean1=tmean0['tmp'].filled() # (given in numpy mask arrays so you need to use the filled method)
	tmean1[tmean1>50]=np.nan
	#plt.matshow(tmean1[1000,:,:])

	fin=metaNA['Paths']['CRU'] + '\\' + metaNA['Current Versions']['CRU'] + '.pre.dat.nc.gz'
	prcp0={}
	with gzip.open(fin) as gz:
			with nc.Dataset('dummy',mode='r',memory=gz.read()) as ds:
				for k in ds.variables.keys():
					prcp0[k]=ds.variables[k][:]
	prcp1=prcp0['pre'].filled()
	prcp1[prcp1>1500]=np.nan

	fin=metaNA['Paths']['CRU'] + '\\' + metaNA['Current Versions']['CRU'] + '.vap.dat.nc.gz'
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
	dC['Data']={}
	for v in vL:
		dC['Data'][v]=np.nan*np.ones( (dC['tv'].shape[0],zRef['m'],zRef['n']) )

	for iT in range(tvCRU.shape[0]):
		ind=np.where( (dC['tv'][:,0]==tvCRU[iT,0]) & (dC['tv'][:,1]==tvCRU[iT,1]) )
		if ind[0].size==0: continue
		dC['Data']['tmean'][ind,:,:]=griddata(xy0,tmean1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')
		dC['Data']['prcp'][ind,:,:]=griddata(xy0,prcp1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')
		dC['Data']['vpd'][ind,:,:]=griddata(xy0,vpd1[iT,:,:][indS],(zRef['X'],zRef['Y']),method='linear')

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
		dC['Data']['rswd'][ind,:,:]=z1*86400/1e6 # Converted from W/m2 to MJ/m2/d-1

	d0=gu.ReadNC(r'C:\Data\Climate\Reanlysis\NCEP Global\Monthly\dswrf.sfc.mon.mean.nc')
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
		dC['Data']['rswd'][ind,:,:]=z1*86400/1e6 # Converted from W/m2 to MJ/m2/d-1

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
			dC['Data']['tmean'][ind1,:,:]=dC['Data']['tmean'][ind2,:,:]
			dC['Data']['prcp'][ind1,:,:]=dC['Data']['prcp'][ind2,:,:]
			dC['Data']['vpd'][ind1,:,:]=dC['Data']['vpd'][ind2,:,:]
			dC['Data']['rswd'][ind1,:,:]=dC['Data']['rswd'][ind2,:,:]

	# Add water balance
	vwbL=['cwd','eta','etp','melt','runoff','ws','wsp']
	for v in vwbL:
		dC['Data'][v]=np.zeros( (dC['tv'].shape[0],zRef['m'],zRef['n']) )
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
		vi['tmean']=np.nan_to_num(dC['Data']['tmean'][iT,:,:][iLand])
		vi['prcp']=np.nan_to_num(dC['Data']['prcp'][iT,:,:][iLand])
		vi['vpd']=np.nan_to_num(dC['Data']['vpd'][iT,:,:][iLand])
		vi['rswn']=np.nan_to_num(dC['Data']['rswd'][iT,:,:][iLand])
		vo=gaia.WBM_SparseGrid(par,vi)
		vi['ws']=vo['ws']
		vi['wsp']=vo['wsp']
		vo['cwd']=vo['etp']-vo['eta']

		for v in vwbL:
			dC['Data'][v][iT,:,:][iLand]=vo[v]

	dC['Data']['cmi']=dC['Data']['prcp']-dC['Data']['etp']

	flg=0
	if flg==1:
		plt.plot(np.nanmean(dC['Data']['tmean'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['prcp'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['etp'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['eta'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['rswd'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['vpd'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['ws'],axis=(1,2)),'-b.')
		plt.plot(np.nanmean(dC['Data']['wsp'],axis=(1,2)))
		plt.plot(np.nanmean(dC['Data']['runoff'],axis=(1,2)))
		plt.plot(np.nanmean(dC['Data']['cwd'],axis=(1,2)))

	# Calculate seasonal indicators
	tva=gu.BlockMean(dC['tv'][:,0],12)
	iTmu=np.where( (tva>=1971) & (tva<=2000) )[0]
	seaL=['ann','mjjas','djf','wyr']
	dCS={}
	dCS['Data']={}
	for v in dC['Data'].keys():
		dCS['Data'][v]={}
		for s in seaL:
			dCS['Data'][v][s]=np.zeros( (tva.shape[0],zRef['m'],zRef['n']) )

	for iY in range(tva.size):

		iT=np.where( (dC['tv'][:,0]==tva[iY]) )[0]
		for v in dC['Data'].keys():
			dCS['Data'][v]['ann'][iY]=np.mean(dC['Data'][v][iT,:,:],axis=0)

		iT=np.where( (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]>=5) & (dC['tv'][:,1]<=9) )[0]
		for v in dC['Data'].keys():
			dCS['Data'][v]['mjjas'][iY]=np.mean(dC['Data'][v][iT,:,:],axis=0)

		iT=np.where( (dC['tv'][:,0]==tva[iY]-1) & (dC['tv'][:,1]>=10) | (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]<=9) )[0]
		for v in dC['Data'].keys():
			dCS['Data'][v]['wyr'][iY]=np.mean(dC['Data'][v][iT,:,:],axis=0)

		iT=np.where( (dC['tv'][:,0]==tva[iY]-1) & (dC['tv'][:,1]==12) | (dC['tv'][:,0]==tva[iY]) & (dC['tv'][:,1]<=2) )[0]
		for v in dC['Data'].keys():
			dCS['Data'][v]['djf'][iY]=np.mean(dC['Data'][v][iT,:,:],axis=0)

	# Convert to anomalies
	dCSa=copy.deepcopy(dCS)
	for k in dCSa['Data'].keys():
		for j in dCSa['Data'][k].keys():
			mu=np.mean(dCSa['Data'][k][j][iTmu,:,:],axis=0)
			for iY in range(tva.size):
				dCSa['Data'][k][j][iY,:,:]=dCSa['Data'][k][j][iY,:,:]-mu

	# Compress
	for k in dCS['Data'].keys():
		for j in dCS['Data'][k].keys():
			dCS['Data'][k][j]=dCS['Data'][k][j]/meta['Climate']['SF'][k]
			dCS['Data'][k][j]=dCS['Data'][k][j].astype('int16')
			dCSa['Data'][k][j]=dCSa['Data'][k][j]/meta['Climate']['SF'][k]
			dCSa['Data'][k][j]=dCSa['Data'][k][j].astype('int16')

	# Save
	dCS['Year']=tva
	dCSa['Year']=tva
	gu.opickle(meta['Paths']['bc5k'] + '\\cru_actuals.pkl',dCS)
	gu.opickle(meta['Paths']['bc5k'] + '\\cru_anomalies.pkl',dCSa)
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
	d={}
	d['tva']=tv
	d['Data']={}
	d['TS Mean']=np.zeros(tv.size)
	for rcp in rcpL:
		d['Data'][rcp]=np.zeros((tv.size,zRef['m'],zRef['n']),dtype='float')
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
			d['Data'][rcp][j,:,:]=z2
			d['TS Mean'][j]=np.mean(z2[iLand])
	
		for iT in range(tv0.size):
			indT=np.where(tv==tv0[iT])[0]
			z0=d0['ntot'][iT,:,:].T # transposed to match lat and lon
			z1=z0[indS]
			z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
			d['Data'][rcp][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['TS Mean'][indT]=np.mean(z2[iLand])
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
			d['Data'][rcp][indT,:,:]=z2
			#plt.matshow(z2,clim=[0,8]); plt.colorbar()
			d['TS Mean'][indT]=np.mean(z2[iLand])
	
		indT=np.where(tv>tv0[-1])[0]
		z0=d0['ntot'][-1,:,:].T # transposed to match lat and lon
		z1=z0[indS]
		z2=griddata( (lon1,lat1),z1,(lon,lat),method='cubic')
		for j in indT:
			d['Data'][j,:,:]=z2
			d['TS Mean'][j]=np.mean(z2[iLand])
		#plt.plot(tv,np.mean(np.mean(d['Data'][rcp],axis=2),axis=1),'k-')
		#plt.close('all'); plt.plot(tv,d['ndep']['TS Mean'],'g-.')
		#plt.matshow(np.mean(d['ntot'],axis=0))
		d['Data'][rcp]=d['Data'][rcp]/meta['Climate']['SF']['ndep']
		d['Data'][rcp]=d['Data'][rcp].astype('int16')

		# Gap fill future missing after 2099
		iT=np.where(tv==2099)[0]
		d['Data'][rcp][-1,:,:]=d['Data'][rcp][iT,:,:]

	# Save
	gu.opickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_ISIMIP.pkl',d)
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
def ProcessNA1K(meta):
	srs=gis.ImportSRSs()
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])

	metaNA=u1k.Init()
	zRefNA=gis.OpenGeoTiff(metaNA['Paths']['na1k Ref Grid'])
	ivlNA=5
	x0,y0=pyproj.transform(srs['Proj']['NACID'],srs['Proj']['BC1ha'],zRefNA['X'],zRefNA['Y'])
	x0=x0[0::ivlNA,0::ivlNA]
	y0=y0[0::ivlNA,0::ivlNA]
	ikpNA=np.where( (zRefNA['Data'][0::ivlNA,0::ivlNA]==1) & (x0>=np.min(zRef['X'][0,:])) & (x0<=np.max(zRef['X'][0,:])) & (y0>=np.min(zRef['Y'][:,0])) & (y0<=np.max(zRef['Y'][:,0])) )
	xy0=np.column_stack((x0[ikpNA].flatten(),y0[ikpNA].flatten()))

	# Climate anomalies

	vL=['tmean','cwd','ws','prcp','etp','eta']
	typL=['wyr','mjjas']

	d={}
	d['Year']=metaNA['tva']
	d['Data']={}
	for v in vL:
		d['Data'][v]={}
		for typ in typL:
			d['Data'][v][typ]=np.zeros((meta['Env']['tva'].size,zRef['m'],zRef['n']),dtype='int16')

	for i,yr in enumerate(metaNA['tva']):
		if yr<1950:
			continue
		try:
			z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_tmean_anom_wyr_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['tmean']
		except:
			z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_tmean_anom_ann_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['tmean']
		zTA=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear') #plt.close('all'); plt.matshow(z1); ind=np.where(zRef['Data']==1); plt.hist(z1[ind].flatten())

		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_ws_anom_mjjas_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['ws']
		zWS=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')

		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_cwd_anom_mjjas_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['cwd']
		zCWD=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')

		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_prcp_anom_mjjas_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['cwd']
		zP=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')

		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_etp_anom_mjjas_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['cwd']
		zEp=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')

		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Summaries\\Time Series\\Anomalies\\na1k_eta_anom_mjjas_' + str(yr) + '.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]#.astype('float')*metaNA['Climate']['SF']['cwd']
		zEa=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear')

		d['Data']['tmean']['wyr'][i,:,:]=np.nan_to_num(zTA)
		d['Data']['etp']['mjjas'][i,:,:]=np.nan_to_num(zEp)
		d['Data']['eta']['mjjas'][i,:,:]=np.nan_to_num(zEa)
		d['Data']['cwd']['mjjas'][i,:,:]=np.nan_to_num(zCWD)
		d['Data']['ws']['mjjas'][i,:,:]=np.nan_to_num(zWS)
		d['Data']['prcp']['mjjas'][i,:,:]=np.nan_to_num(zP)

	del d['Data']['cwd']['wyr']
	del d['Data']['eta']['wyr']
	del d['Data']['etp']['wyr']
	del d['Data']['ws']['wyr']
	del d['Data']['prcp']['wyr']
	del d['Data']['tmean']['mjjas']

	gu.opickle(meta['Paths']['bc5k'] + '\\na1k_climate_anomalies.pkl',d)

	# Ndep (ECCC 2021)

	z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\ANN_2021_Accumulated_Deposition.tif')['Data'][0::ivlNA,0::ivlNA][ikpNA]
	zN=griddata(xy0,z0,(zRef['X'],zRef['Y']),method='linear') #plt.close('all'); plt.matshow(z1); ind=np.where(zRef['Data']==1); plt.hist(z1[ind].flatten())
	#zN=zN/meta['Climate']['SF']['ndep']
	zN=zN.astype('int16')
	gu.opickle(meta['Paths']['bc5k'] + '\\na1k_ndep_eccc_2021.pkl',zN)

	return

#%%
def CompileEnvironmentalData(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	metaNA=u1k.Init()

	#--------------------------------------------------------------------------
	# Import historical climate data
	#--------------------------------------------------------------------------
	dNA=gu.ipickle(meta['Paths']['bc5k'] + '\\na1k_climate_anomalies.pkl')
	dCRU=gu.ipickle(meta['Paths']['bc5k'] + '\\cru_anomalies.pkl')
	dR=gu.ipickle(meta['Paths']['bc5k'] + '\\20thCR_anomalies.pkl')

	typL=['CRU','20CR','NA1K','BC5K Comp1']
	vL=['tmean','cwd','ws','prcp','etp','eta']
	seaL=['wyr','mjjas']
	dHC={}
	for typ in typL:
		dHC[typ]={}
		for v in vL:
			dHC[typ][v]={}
			for sea in seaL:
				dHC[typ][v][sea]=-999*np.ones((meta['Env']['tva'].size,zRef['m'],zRef['n']),dtype='int16')

	for i,yr in enumerate(meta['Env']['tva']):
		typ='CRU'
		iT=np.where(dCRU['Year']==yr)[0]
		if iT.size>0:
			dHC[typ]['tmean']['wyr'][i,:,:]=dCRU['Data']['tmean']['wyr'][iT,:,:]#.astype('float')*meta['Climate']['SF']['tmean']
			dHC[typ]['cwd']['mjjas'][i,:,:]=dCRU['Data']['cwd']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['cwd']
			dHC[typ]['ws']['mjjas'][i,:,:]=dCRU['Data']['ws']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['prcp']['mjjas'][i,:,:]=dCRU['Data']['prcp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['etp']['mjjas'][i,:,:]=dCRU['Data']['etp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['eta']['mjjas'][i,:,:]=dCRU['Data']['eta']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
		typ='20CR'
		iT=np.where(dR['Year']==yr)[0]
		if iT.size>0:
			dHC[typ]['tmean']['wyr'][i,:,:]=dR['Data']['tmean']['wyr'][iT,:,:]#.astype('float')*meta['Climate']['SF']['tmean']
			dHC[typ]['cwd']['mjjas'][i,:,:]=dR['Data']['cwd']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['cwd']
			dHC[typ]['ws']['mjjas'][i,:,:]=dR['Data']['ws']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['prcp']['mjjas'][i,:,:]=dR['Data']['prcp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['etp']['mjjas'][i,:,:]=dR['Data']['etp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['eta']['mjjas'][i,:,:]=dR['Data']['eta']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
		typ='NA1K'
		iT=np.where(dNA['Year']==yr)[0]
		if iT.size>0:
			dHC[typ]['tmean']['wyr'][i,:,:]=dNA['Data']['tmean']['wyr'][iT,:,:]#.astype('float')*meta['Climate']['SF']['tmean']
			dHC[typ]['cwd']['mjjas'][i,:,:]=dNA['Data']['cwd']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['cwd']
			dHC[typ]['ws']['mjjas'][i,:,:]=dNA['Data']['ws']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['prcp']['mjjas'][i,:,:]=dNA['Data']['prcp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['etp']['mjjas'][i,:,:]=dNA['Data']['etp']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['eta']['mjjas'][i,:,:]=dNA['Data']['eta']['mjjas'][iT,:,:]#.astype('float')*meta['Climate']['SF']['ws']

		typ='BC5K Comp1'
		if (yr>=1880) & (yr<1950):
			v='tmean'; dHC[typ][v]['wyr'][i,:,:]=(dHC['CRU'][v]['wyr'][i,:,:]+dHC['20CR'][v]['wyr'][i,:,:])/2
			v='prcp'; dHC[typ][v]['mjjas'][i,:,:]=dHC['CRU'][v]['mjjas'][i,:,:]
			v='etp'; dHC[typ][v]['mjjas'][i,:,:]=dHC['CRU'][v]['mjjas'][i,:,:] # 20CR is unlikely
			#v='etp'; dHC[typ][v]['mjjas'][i,:,:]=(dHC['CRU'][v]['mjjas'][i,:,:]+dHC['20CR'][v]['mjjas'][i,:,:])/2
			v='eta'; dHC[typ][v]['mjjas'][i,:,:]=dHC['CRU'][v]['mjjas'][i,:,:]#+dHC['20CR'][v]['mjjas'][i,:,:])/2
			v='cwd'; dHC[typ][v]['mjjas'][i,:,:]=dHC['CRU'][v]['mjjas'][i,:,:]#+dHC['20CR'][v]['mjjas'][i,:,:])/2
			v='ws'; dHC[typ][v]['mjjas'][i,:,:]=dHC['CRU'][v]['mjjas'][i,:,:]#+dHC['20CR'][v]['mjjas'][i,:,:])/2
		if (yr>=1950) & (yr<=metaNA['YearLast']):
			dHC[typ]['tmean']['wyr'][i,:,:]=dNA['Data']['tmean']['wyr'][i,:,:]#.astype('float')*meta['Climate']['SF']['tmean']
			dHC[typ]['cwd']['mjjas'][i,:,:]=dNA['Data']['cwd']['mjjas'][i,:,:] #.astype('float')*meta['Climate']['SF']['cwd']
			dHC[typ]['ws']['mjjas'][i,:,:]=dNA['Data']['ws']['mjjas'][i,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['prcp']['mjjas'][i,:,:]=dNA['Data']['prcp']['mjjas'][i,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['etp']['mjjas'][i,:,:]=dNA['Data']['etp']['mjjas'][i,:,:]#.astype('float')*meta['Climate']['SF']['ws']
			dHC[typ]['eta']['mjjas'][i,:,:]=dNA['Data']['eta']['mjjas'][i,:,:]#.astype('float')*meta['Climate']['SF']['ws']

	# Delete empty variable summaries
	for k in dHC.keys():
		del dHC[k]['tmean']['mjjas']
		del dHC[k]['cwd']['wyr']
		del dHC[k]['ws']['wyr']
		del dHC[k]['prcp']['wyr']
		del dHC[k]['etp']['wyr']
		del dHC[k]['eta']['wyr']

	gu.opickle(meta['Paths']['bc5k'] + '\\bc5k_historical_climate_anomalies.pkl',dHC)

	#--------------------------------------------------------------------------
	# Add future climate data to historical data
	#--------------------------------------------------------------------------

	esmL=['HadGEM3-GC31-LL','GFDL-ESM4','MIROC-ES2L','MPI-ESM1-2-LR','CNRM-CM6-1','CESM2']
	scnL=['ssp245','ssp585']

	dEC={}
	dEC['Climate']={}
	for scn in scnL:
		dEC['Climate'][scn]={}
		for esm in esmL:
			dEC['Climate'][scn][esm]={}
			for v in vL:
				dEC['Climate'][scn][esm][v]={}
				for sea in seaL:
					dEC['Climate'][scn][esm][v][sea]=-999*np.ones((meta['Env']['tva'].size,zRef['m'],zRef['n']),dtype='int16')

	# Import future model data, apply bias correction
	typ='BC5K Comp1'
	bias_cor_period=[1951,metaNA['YearLast']]
	tvM=np.arange(1950,2101,1)
	iT0=np.where( (meta['Env']['tva']>=bias_cor_period[0]) & (meta['Env']['tva']<=bias_cor_period[1]) )[0]
	iT1=np.where( (tvM>=bias_cor_period[0]) & (tvM<=bias_cor_period[1]) )[0]
	iT2=np.where( (meta['Env']['tva']>=metaNA['YearLast']) & (meta['Env']['tva']<=tvM[-1]) )[0]
	iT3=np.where( (tvM>=metaNA['YearLast']) & (tvM<=meta['Env']['tva'][-1]) )[0]
	iT4=np.where( (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	for scn in scnL:
		dM=gu.ipickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_anomalies_' + scn + '.pkl')
		for esm in esmL:
			# Temp
			y=dM[esm]['tmean']['wyr']#.astype(float)*meta['Climate']['SF']['tmean']
			yH_mu=np.mean(dHC[typ]['tmean']['wyr'][iT0,:,:],axis=0)
			yM_mu=np.mean(y[iT1,:,:],axis=0)
			cf=np.tile(np.squeeze(yM_mu-yH_mu),(tvM.size,1,1))
			y=y-cf
			dEC['Climate'][scn][esm]['tmean']['wyr'][iT2,:,:]=y[iT3,:,:]
			dEC['Climate'][scn][esm]['tmean']['wyr'][iT4,:,:]=dHC[typ]['tmean']['wyr'][iT4,:,:]

			# CWD
			y=dM[esm]['cwd']['mjjas']#.astype(float)*meta['Climate']['SF']['cwd']
			yH_mu=np.mean(dHC[typ]['cwd']['mjjas'][iT0,:,:],axis=0)
			yM_mu=np.mean(y[iT1,:,:],axis=0)
			cf=np.tile(np.squeeze(yM_mu-yH_mu),(tvM.size,1,1))
			y=y-cf
			dEC['Climate'][scn][esm]['cwd']['mjjas'][iT2,:,:]=y[iT3,:,:]
			dEC['Climate'][scn][esm]['cwd']['mjjas'][iT4,:,:]=dHC[typ]['cwd']['mjjas'][iT4,:,:]

			# WS
			y=dM[esm]['ws']['mjjas']#.astype(float)*meta['Climate']['SF']['ws']
			yH_mu=np.mean(dHC[typ]['ws']['mjjas'][iT0,:,:],axis=0)
			yM_mu=np.mean(y[iT1,:,:],axis=0)
			cf=np.tile(np.squeeze(yM_mu-yH_mu),(tvM.size,1,1))
			y=y-cf
			dEC['Climate'][scn][esm]['ws']['mjjas'][iT2,:,:]=y[iT3,:,:]
			dEC['Climate'][scn][esm]['ws']['mjjas'][iT4,:,:]=dHC[typ]['ws']['mjjas'][iT4,:,:]

	# Delete empty variable summaries
	for k1 in dEC['Climate'].keys():
		for k2 in dEC['Climate'][k1].keys():
			del dEC['Climate'][k1][k2]['tmean']['mjjas']
			del dEC['Climate'][k1][k2]['cwd']['wyr']
			del dEC['Climate'][k1][k2]['ws']['wyr']

	# Import CO2 (ppm), source: https://www.pik-potsdam.de/~mmalte/rcps/ (ppm)
	dCO2=gu.ReadExcel(meta['Paths']['DB']['CO2'],sheet_name='Sheet1',skiprows=0)
	dEC['co2']={}
	dEC['co2']['ssp245']=-999*np.ones(meta['Env']['tva'].size,dtype='int16')
	dEC['co2']['ssp585']=-999*np.ones(meta['Env']['tva'].size,dtype='int16')
	for i,yr in enumerate(meta['Env']['tva']):
		iT=np.where(dCO2['Year']==yr)[0]
		dEC['co2']['ssp245'][i]=dCO2['CO2 RCP45'][iT]/meta['Climate']['SF']['co2']
		dEC['co2']['ssp585'][i]=dCO2['CO2 RCP85'][iT]/meta['Climate']['SF']['co2']

	# Import nitrogen deposition
	d0=gu.ipickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_ISIMIP.pkl')
	dEC['ndep']={}
	dEC['ndep']['ISIMIP']={}
	for s in ['26','60','85']:
		dEC['ndep']['ISIMIP'][s]=-999*np.ones((meta['Env']['tva'].size,zRef['Data'].shape[0],zRef['Data'].shape[1]),dtype='int16')
		for i,yr in enumerate(meta['Env']['tva']):
			iT=np.where(d0['tva']==yr)[0]
			if iT.size==0: continue
			dEC['ndep']['ISIMIP'][s][i,:,:]=d0['Data'][s][iT,:,:]#.astype('float')*meta['Climate']['SF']['ndep']

	dNEC=gu.ipickle(meta['Paths']['bc5k'] + '\\na1k_ndep_eccc_2021.pkl')#.astype('float')*meta['Climate']['SF']['ndep']
	dEC['ndep']['Comp 1']={}
	for s in ['26','60','85']:
		y0=dEC['ndep']['ISIMIP'][s]
		iT=np.where(meta['Env']['tva']==2021)[0]
		dY=y0[iT,:,:]-dNEC
		# print(np.nanmean(dY))
		dEC['ndep']['Comp 1'][s]=y0-np.tile(dY,(meta['Env']['tva'].size,1,1))

	# Save
	gu.opickle(meta['Paths']['bc5k'] + '\\bc5k_env_comp1.pkl',dEC)

	return

#%%
def LoadEnvironmentalComp1(meta):

	dEC=gu.ipickle(meta['Paths']['bc5k'] + '\\bc5k_env_comp1.pkl')

	for k1 in dEC['Climate'].keys():
		for k2 in dEC['Climate'][k1].keys():
			for k3 in dEC['Climate'][k1][k2].keys():
				for k4 in dEC['Climate'][k1][k2][k3].keys():
					# Convert to float
					dEC['Climate'][k1][k2][k3][k4]=dEC['Climate'][k1][k2][k3][k4].astype('float')

					# Convert missing number to NAN
					ind=np.where(dEC['Climate'][k1][k2][k3][k4]==-999)
					dEC['Climate'][k1][k2][k3][k4][ind]=np.nan

					# Apply scale factor
					dEC['Climate'][k1][k2][k3][k4]=dEC['Climate'][k1][k2][k3][k4]*meta['Climate']['SF'][k3]

	for k1 in dEC['co2'].keys():
		dEC['co2'][k1]=dEC['co2'][k1].astype('float')
		ind=np.where(dEC['co2'][k1]==-999)
		dEC['co2'][k1][ind]=np.nan
		dEC['co2'][k1]=dEC['co2'][k1]*meta['Climate']['SF']['co2']

	for k1 in dEC['ndep'].keys():
		for k2 in dEC['ndep'][k1].keys():
			dEC['ndep'][k1][k2]=dEC['ndep'][k1][k2].astype('float')
			ind=np.where(dEC['ndep'][k1][k2]==-999)
			dEC['ndep'][k1][k2][ind]=np.nan
			dEC['ndep'][k1][k2]=dEC['ndep'][k1][k2]*meta['Climate']['SF']['ndep']

	dHC=gu.ipickle(meta['Paths']['bc5k'] + '\\bc5k_historical_climate_anomalies.pkl')
	for k1 in dHC.keys():
		for k2 in dHC[k1].keys():
			for k3 in dHC[k1][k2].keys():
				# Convert to float
				dHC[k1][k2][k3]=dHC[k1][k2][k3].astype('float')
				ind=np.where(dHC[k1][k2][k3]==-999)
				dHC[k1][k2][k3][ind]=np.nan

				# Convert missing number to NAN
				ind=np.where(dHC[k1][k2][k3]==-999)
				dHC[k1][k2][k3][ind]=np.nan

				# Apply scale factor
				dHC[k1][k2][k3]=dHC[k1][k2][k3]*meta['Climate']['SF'][k2]

	return dHC,dEC

#%%
def Plot_HistoricalAssessment(meta):

	metaNA=u1k.Init()
	dHC,dEC=LoadEnvironmentalComp1(meta)
	dAHCCD=gu.ipickle(r'C:\Data\Climate\AHCCD\2025\ahccd.pkl')
	ma=10; xlim=[1871,2025]; lw=2;

	# Temperature

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['tmean']['wyr'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['tmean']['wyr'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')

	iBC=np.where( (dAHCCD['Stations Temp']['Prov']=='BC') & (dAHCCD['Stats']['tmean']['Percent Complete']>0) )[0]
	y=np.nanmean(dAHCCD['Data']['Seasonal']['wyr']['tmean_a'][:,iBC],axis=1)
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(dAHCCD['tva'],y_ma,'m-',lw=lw,color=[1,0.55,0.9],label='ECCC-AHCCD')
	
	y=np.nanmean(dHC['NA1K']['tmean']['wyr'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['tmean']['wyr'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')
	
	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.7,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' ($\circ$C)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average mean annual\ntemperature anomaly ($\circ$C)',xlabel='Time, years',yticks=np.arange(-3.2,3.2,0.4),
		xlim=xlim,ylim=[-1.6,1])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\tmean_HistoricalAssessment','png',900)

	# Precipitation

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['prcp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['prcp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')

	flg=1
	if flg==1:
		iBC=np.where( (dAHCCD['Stations Precip']['Prov']=='BC') & (dAHCCD['Stats']['prcp']['Percent Complete']>0) )[0]
		y=np.nanmean(dAHCCD['Data']['Seasonal']['mjjas']['prcp_a'][:,iBC],axis=1)/30.5
		y_ma=gu.movingave(y,ma,'historical')
		ax.plot(dAHCCD['tva'],y_ma,'m-',lw=lw,color=[1,0.55,0.9],label='ECCC-AHCCD')
	
	y=np.nanmean(dHC['NA1K']['prcp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['prcp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')

	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.4,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')

	ax.set(ylabel='Ten-year moving average\n MJJAS precipitation anomaly (mm d$^{-1}$)',xlabel='Time, years',yticks=np.arange(-30,30,0.2),
		xlim=xlim,ylim=[-0.5,0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\prcp_HistoricalAssessment','png',900)

	# Potential evapotranspiration

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')
	
	y=np.nanmean(dHC['NA1K']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')
	
	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.145,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average MJJAS potential\nevapotranspiration anomaly (mm d$^{-1}$)',xlabel='Time, years',yticks=np.arange(-3,3,0.05),
		xlim=xlim,ylim=[-0.17,0.17])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\etp_HistoricalAssessment','png',900)

	# Actual evapotranspiration

	y=np.nanmean(dHC['BC5K Comp1']['prcp']['mjjas'],axis=(1,2))/30.5
	plt.close('all'); plt.plot(meta['Env']['tva'],y,'.k-',lw=0.6)

	y1=np.nanmean(dHC['NA1K']['etp']['mjjas'],axis=(1,2))/30.5
	y2=np.nanmean(dHC['NA1K']['eta']['mjjas'],axis=(1,2))/30.5
	plt.plot(y1,y2,'bo')

	y1=np.nanmean(dHC['NA1K']['ws']['mjjas'],axis=(1,2))/30.5
	y2=np.nanmean(dHC['NA1K']['cwd']['mjjas'],axis=(1,2))/30.5
	plt.close('all'); plt.plot(y1,y2,'bo')

	y1=np.nanmean(dHC['NA1K']['prcp']['mjjas'],axis=(1,2))/30.5
	y2=np.nanmean(dHC['NA1K']['ws']['mjjas'],axis=(1,2))/30.5
	plt.close('all'); plt.plot(y1,y2,'bo')

	y1=np.nanmean(dHC['NA1K']['prcp']['mjjas'],axis=(1,2))/30.5
	y2=np.nanmean(dHC['NA1K']['cwd']['mjjas'],axis=(1,2))/30.5
	plt.close('all'); plt.plot(y1,y2,'bo')

	y1=np.nanmean(dHC['NA1K']['eta']['mjjas'],axis=(1,2))/30.5
	y2=np.nanmean(dHC['NA1K']['cwd']['mjjas'],axis=(1,2))/30.5
	plt.close('all'); plt.plot(y1,y2,'bo')

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['eta']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['eta']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')
	
	y=np.nanmean(dHC['NA1K']['eta']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['eta']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')
	
	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.145,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average MJJAS actual\nevapotranspiration anomaly (mm d$^{-1}$)',xlabel='Time, years',yticks=np.arange(-3,3,0.05),
		xlim=xlim,ylim=[-0.17,0.17])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\eta_HistoricalAssessment','png',900)

	# CMI

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['prcp']['mjjas']-dHC['20CR']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['prcp']['mjjas']-dHC['CRU']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')

	y=np.nanmean(dHC['NA1K']['prcp']['mjjas']-dHC['NA1K']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['prcp']['mjjas']-dHC['BC5K Comp1']['etp']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')

	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.22,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average MJJAS\nclimate moisture index anomaly (mm d$^{-1}$)',xlabel='Time, years',yticks=np.arange(-30,30,0.2),
		xlim=xlim,ylim=[-0.5,0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\cmi_HistoricalAssessment','png',900)

	# Climatic water deficit

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['cwd']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['cwd']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')
	
	y=np.nanmean(dHC['NA1K']['cwd']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['cwd']['mjjas'],axis=(1,2))/30.5
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')
	
	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,0.7,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average MJJAS climatic\nwater deficit anomaly (mm d$^{-1}$)',xlabel='Time, years',yticks=np.arange(-3.2,3.2,0.02),
		xlim=xlim,ylim=[-0.1,0.1])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\cwd_HistoricalAssessment','png',900)

	# Soil water

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
	ax.plot([1871,2025],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	
	y=np.nanmean(dHC['20CR']['ws']['mjjas'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'r-',lw=lw,color=[1,0.9,0.4],label='20thC reanalysis')
	
	y=np.nanmean(dHC['CRU']['ws']['mjjas'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'g-',lw=lw,color=[0.7,0.9,0.3],label='CRU-TS')
	
	y=np.nanmean(dHC['NA1K']['ws']['mjjas'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	iT1=np.where( (meta['Env']['tva']>1956) & (meta['Env']['tva']<=2022) )[0]
	ax.plot(meta['Env']['tva'][iT1],y_ma[iT1],'b-',lw=lw,color=[0.55,0.76,1],label='NA1K')
	
	y=np.nanmean(dHC['BC5K Comp1']['ws']['mjjas'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	ax.plot(meta['Env']['tva'],y_ma,'k--',lw=1.5,label='BC5K Comp 1')
	
	iT1=np.where( (meta['Env']['tva']>=1901) & (meta['Env']['tva']<=metaNA['YearLast']) )[0]
	rs,txt=gu.GetRegStats(meta['Env']['tva'][iT1],y[iT1])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',label='BC5K Comp 1 Best Fit Line')
	ax.text(1875,12.5,'Change in BC5K Comp1 over 1901-2024: ' + str(np.round(rs['B'][1]*iT1.size,decimals=1)) + ' (mm d$^{-1}$)',fontsize=8,ha='left')
	
	ax.set(ylabel='Ten-year moving average MJJAS soil\nwater content anomaly (mm)',xlabel='Time, years',yticks=np.arange(-500,500,5),
		xlim=xlim,ylim=[-15,15])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\ws_HistoricalAssessment','png',900)

	return

#%%
def Plot_EnvironmentalTimeSeries(meta):

	dHC,dEC=LoadEnvironmentalComp1(meta)

	metaNA=u1k.Init()
	tva=meta['Env']['tva']
	ma=10;
	xlim=[1880,2100];
	xt=np.arange(1900,2200+50,50)
	lw1=1.5; lw2=0.25
	cl=np.array([[0.2,1,1],[0.8,0,0]])
	scnN=['CMIP6 model mean SSP-4.5','CMIP6 model mean SSP-8.5']

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,10))

	# Temperature

	ax[0,0].plot([0,3000],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)

	y=np.nanmean(dHC['BC5K Comp1']['tmean']['wyr'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	iT=np.where(tva<=metaNA['YearLast'])[0]
	ax[0,0].plot(tva[iT],y_ma[iT],'k-',lw=lw1,label='Historical')

	iT=np.where(tva>metaNA['YearLast'])[0]
	for i,scn in enumerate(dEC['Climate'].keys()):
		y=np.zeros(meta['Env']['tva'].size)
		for esm in dEC['Climate'][scn].keys():
			y0=np.nanmean(dEC['Climate'][scn][esm]['tmean']['wyr'],axis=(1,2))
			y=y+y0
			y_ma=gu.movingave(y0,ma,'historical')
			ax[0,0].plot(tva[iT],y_ma[iT],'k-',lw=lw2,color=cl[i,:])
		y=y/len(dEC['Climate'][scn].keys())
		y_ma=gu.movingave(y,ma,'historical')
		ax[0,0].plot(tva[iT],y_ma[iT],'k-',lw=lw1,color=cl[i,:],label=scnN[i])

	ax[0,0].set(ylabel='Mean annual temperature\nanomaly ($^\circ$C)',xlabel='Time, years',xticks=xt,xlim=xlim)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=metaNA['Graphics']['gp']['tickl'])
	ax[0,0].legend(loc='center left',frameon=False,facecolor=None,edgecolor='w')

	# Hydrology

	ax[0,1].plot([0,3000],[0,0],'k-',color=[0.8,0.8,0.8],lw=1.5)
	v='cwd'; yl='Climatic water deficit\nMMJAS anomaly (mm d$^{-1}$)'
	#v='ws'; yl='Warm-season soil water\ncontent anomaly (mm)'
	y=np.nanmean(dHC['BC5K Comp1'][v]['mjjas'],axis=(1,2))
	y_ma=gu.movingave(y,ma,'historical')
	iT=np.where(tva<=metaNA['YearLast'])[0]
	ax[0,1].plot(tva[iT],y_ma[iT],'k-',lw=lw1,label='BC5K Comp 1')
	iT=np.where(tva>metaNA['YearLast'])[0]
	for i,scn in enumerate(dEC['Climate'].keys()):
		y=np.zeros(meta['Env']['tva'].size)
		for esm in dEC['Climate'][scn].keys():
			y0=np.nanmean(dEC['Climate'][scn][esm][v]['mjjas'],axis=(1,2))
			y=y+y0
			y_ma=gu.movingave(y0,ma,'historical')
			ax[0,1].plot(tva[iT],y_ma[iT],'k-',lw=lw2,color=cl[i,:])
		y=y/len(dEC['Climate'][scn].keys())
		y_ma=gu.movingave(y,ma,'historical')
		ax[0,1].plot(tva[iT],y_ma[iT],'k-',lw=lw1,color=cl[i,:],label=scnN[i])

	ax[0,1].set(ylabel=yl,xlabel='Time, years',xticks=xt,xlim=xlim)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=metaNA['Graphics']['gp']['tickl'])

	# CO2

	iTh=np.where(tva<=metaNA['YearLast'])[0]
	iTf=np.where(tva>metaNA['YearLast'])[0]

	ax[1,0].plot(tva[iTh],dEC['co2']['ssp245'][iTh],'k-',lw=lw1,label='Historical')
	ax[1,0].plot(tva[iTf],dEC['co2']['ssp245'][iTf],'c-',color=cl[0,:],lw=lw1,label='SSP245')
	ax[1,0].plot(tva[iTf],dEC['co2']['ssp585'][iTf],'r-',color=cl[1,:],lw=lw1,label='SSP585')
	ax[1,0].set(ylabel='Atmospheric carbon dioxide\nconcentration (ppm)',xlabel='Time, years',xticks=xt,xlim=xlim,ylim=[0,1000])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=metaNA['Graphics']['gp']['tickl'])

	# N deposition

	ax[1,1].plot(tva[iTh],np.nanmean(dEC['ndep']['Comp 1']['26'][iTh],axis=(1,2)),'k-',lw=lw1,label='Historical')
	ax[1,1].plot(tva[iTf],np.nanmean(dEC['ndep']['Comp 1']['26'][iTf],axis=(1,2)),'c-',color=cl[0,:],lw=lw1,label='SSP245')
	ax[1,1].plot(tva[iTf],np.nanmean(dEC['ndep']['Comp 1']['85'][iTf],axis=(1,2)),'r-',color=cl[1,:],lw=lw1,label='SSP585')
	ax[1,1].set(ylabel='Nitrogen deposition\n(kgN ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xticks=xt,xlim=xlim,ylim=[0,2.5])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=metaNA['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.03,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EnvironmentalData','png',900)

	return

#%%