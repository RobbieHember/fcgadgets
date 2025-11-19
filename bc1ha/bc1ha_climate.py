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

#%%
def ImportNormalsFromClimateNA(meta):
	# 2024 - these are tiffs from TW because the .asc temperatures were not giving decimals
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zBuf=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid Buf'])
	zBufR=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid Buf Rev'])
	zE=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation_geog5.tif')
	iOut=np.where(zE['Data']>6000)
	pthin=r'D:\Backup Data\Climate\ClimateNA\BC1ha'
	vL=['Tave','Tmin','Tmax','PPT']
	vL2=['tmean','tmin','tmax','prcp']
	#vL=['PPT']
	#vL2=['prcp']
	for iV in range(len(vL)):
		for mo in range(1,13):
			if mo<=9:
				mot='0' + str(mo)
			else:
				mot=str(mo)

			# Received as geotiffs from TW
			z0=gis.OpenGeoTiff(pthin + '\\' + vL[iV] + mot + '.tif')
			z=copy.deepcopy(zE)
			z['Data']=np.zeros(zE['Data'].shape,dtype='float')
			z['Data']=z0['Data'].astype('float')
			z['Data']=np.minimum(327,z['Data'])
			#plt.close('all'); plt.matshow(z['Data'],clim=[0,150])
			z['Data'][iOut]=meta['Climate']['Missing Number']
			z['Data']=z['Data']/meta['Climate']['SF'][vL2[iV]]
			z['Data']=z['Data'].astype('int16')
			gis.SaveGeoTiff(z,pthin + '\\tmp2.tif')

			# Reproject and change extent
			fin=pthin + '\\tmp2.tif'
			fout=meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_' + vL2[iV] + '_norm_1971to2000_' + str(mo) + '.tif'
			fref=meta['Paths']['bc1ha Ref Grid']
			gis.ReprojectRasterAndClipToRaster(fin,fout,fref,meta['Geos']['crs'])

			# Import reprojected map
			z0=gis.OpenGeoTiff(fout)['Data'].astype('float')*meta['Climate']['SF'][vL2[iV]]

			# Remove exterior
			z0[zRef['Data']==0]=meta['Climate']['Missing Number'] #plt.matshow(z0)

			# Find missing values within mask
			#iMissing=np.where( (zRef['Data']==1) & (z0==meta['Climate']['Missing Number']) )
			#print(iMissing[0].size)

			# Gap-fill edges
			MaskGap=-1*(zRef['Data'].astype(float)-zBuf['Data'].astype(float)); plt.close('all'); #plt.matshow(MaskGap); plt.colorbar()
			indGap=np.where( (MaskGap==1) | (zRef['Data']==1) & (z0==meta['Climate']['Missing Number']) )
			MaskCal=zRef['Data']-zBufR['Data'] #plt.matshow(MaskCal)
			indCal=np.where( (MaskCal==1) & (z0!=meta['Climate']['Missing Number']) )
			#a=MaskGap.copy();a[MaskCal==1]=2;plt.close('all'); plt.matshow(a)
			x=zRef['X'][indCal]
			y=zRef['Y'][indCal]
			z=z0[indCal]
			interp=NearestNDInterpolator(list(zip(x,y)),z)
			iz=interp(zRef['X'][indGap],zRef['Y'][indGap])
			z0[indGap]=iz # plt.close('all');plt.matshow(z0)
			plt.close('all'); plt.matshow(z0,clim=[0,350])

			# Save
			z=copy.deepcopy(zRef)
			z['Data']=z0
			z['Data']=z['Data']/meta['Climate']['SF'][vL2[iV]]
			z['Data']=z['Data'].astype('int16')
			gis.SaveGeoTiff(z,fout)
	return

#%%
def CalcSaturationVapourPressureNormal(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	for mo in range(1,13):
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_tmean_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['tmean']
		z=copy.deepcopy(zRef)
		z['Data']=gaia.GetEstar(zT)
		z['Data'][zRef['Data']==0]=meta['Climate']['Missing Number']
		z['Data']=z['Data']/meta['Climate']['SF']['es']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_es_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcActualVapourPressureNormalFromTemps(meta):
	# Reference: https://www.agraria.unirc.it/documentazione/materiale_didattico/1462_2016_412_24509.pdf
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	for mo in range(1,13):
		zT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_tmin_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['tmin']
		z=copy.deepcopy(zRef)
		z['Data']=10*(0.6108*np.exp((17.27*zT)/(zT+237.3)))
		z['Data'][zRef['Data']==0]=meta['Climate']['Missing Number']
		z['Data']=z['Data']/meta['Climate']['SF']['ea']
		z['Data']=z['Data'].astype('int16')
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcActualVapourPressureNormalFromNA1k(meta):
	metaNA=u1k.Init()
	fRef=meta['Paths']['bc1ha Ref Grid']
	xlim=[-2100000,-1000000]; ylim=[1100000,3000000]
	f1=meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\tmp.tif'
	for mo in range(1,13):
		z0=gis.OpenGeoTiff(metaNA['Paths']['na1k'] + '\\Monthly\\Normals\\na1k_ea_biasadj_norm_1971to2000_' + str(mo) + '.tif')
		z0=gis.ClipRasterByXYLimits(z0,xlim,ylim)
		gis.SaveGeoTiff(z0,f1)
		f2=meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_biasadj_from_na1k_norm_1971to2000_' + str(mo) + '.tif'
		gis.ReprojectRasterAndClipToRaster(f1,f2,fRef,meta['Geos']['crs'])
		os.remove(f1)
	return

#%%
def CalcActualVapourPressureNormalBiasCorrected(meta):
	meta['tvm']=gu.tvec('m',1950,2025)

	# Import GSOD
	meta['Paths']['GSOD']=r'C:\Data\Climate\GSOD'
	dS=gu.ipickle(meta['Paths']['GSOD'] + '\\Processed\\gsod.pkl')
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
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
	ea_cru[ea_cru>999]=meta['Climate']['Missing Number']
	#ind=np.where(ea_cru>meta['Climate']['Missing Number']); print(np.mean(ea_cru[ind]))
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

		ea_mu0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['ea']
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
			ivl=1
			iz=griddata(np.column_stack((x0,y0)),e0,(zRef['X'][0::ivl,0::ivl],zRef['Y'][0::ivl,0::ivl]),method='linear')
			#iz=scipy.ndimage.zoom(iz,ivl,order=0)
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
		ax.set(xlabel='GSOD measurements (hPa)',ylabel='bc1ha grid (hPa)',xticks=np.arange(0,40,5),yticks=np.arange(0,40,5),xlim=[0,35],ylim=[0,35])
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
		fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
		ax.plot(x,y,'o',ms=3,mec='w',mfc='k',mew=0.5)
		ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
		ax.text(33,2,rs['txt'],fontsize=7,ha='right')
		ax.text(33,33,'1:1',fontsize=7,ha='center')
		ax.set(xlabel='GSOD measurements (hPa)',ylabel='bc1ha grid (hPa)',xticks=np.arange(0,40,5),yticks=np.arange(0,40,5),xlim=[0,35],ylim=[0,35])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		plt.tight_layout()

		# Save new actual vapour pressure
		z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_fromtmin_norm_1971to2000_' + str(mo) + '.tif')
		tmp=np.maximum(0,(z1['Data'].astype('float')*meta['Climate']['SF']['ea'])-iz)
		tmp=tmp/meta['Climate']['SF']['ea']
		z1['Data']=tmp.astype('int16')
		# plt.matshow(z1['Data'],clim=[0,3000])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_biasadj_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def CalcVapourPressureDeficitNormal(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	for mo in range(1,13):
		#zEa=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_biasadj_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['ea']
		zEa=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_ea_biasadj_from_na1k_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['ea']

		zEs=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_es_norm_1971to2000_' + str(mo) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['es']
		vpd=np.maximum(0,zEs-zEa)

		z=copy.deepcopy(zRef)
		z['Data']=vpd
		z['Data'][zRef['Data']==0]=meta['Climate']['Missing Number']
		z['Data']=z['Data']/meta['Climate']['SF']['vpd']
		z['Data']=z['Data'].astype('int16')
		# plt.matshow(z['Data'],clim=[0,2000])
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_vpd_norm_1971to2000_' + str(mo) + '.tif')
	return

#%%
def ConvertBCElevationToGeographicForClimateNA(meta):
	# ClimateNA needs a .asc geog file, so export DEM to goeg, in ArcGIS create .asc, then run ClimateNA
	pthin=meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif'
	pthout=meta['Paths']['bc1ha'] + '\\Terrain\\elevation_geog5.tif'
	pb=gpd.read_file(r'C:\Data\Geodatabases\North America\bound_p.shp')#.to_crs(srs['String']['NACID'])
	gis.ReprojectGeoTiff(pthin,pthout,pb.crs)
	zE=gis.OpenGeoTiff(pthout);
	#plt.close('all'); plt.matshow(zE['Data'])
	zE=gis.ClipRasterByXYLimits(zE,[-139.15,-114],[48.25,60.05])
	gis.SaveGeoTiff(zE,pthout)
	return

#%%
def CalcSurfaceWaterBalanceNormals(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	ivl=1
	iMask=np.where(zRef['Data'][0::ivl,0::ivl]>0)

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

	con=gaia.HydroMetCon()

	# Import normals
	nrms={}
	nrms['tmean']=[None]*12
	nrms['rswn']=[None]*12
	nrms['prcp']=[None]*12
	nrms['vpd']=[None]*12
	for iM in range(12):
		nrms['tmean'][iM]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_tmean_norm_1971to2000_' + str(iM+1) + '.tif')['Data'][0::ivl,0::ivl][iMask].astype('float')*meta['Climate']['SF']['tmean']
		nrms['rswn'][iM]=(1-con['Albedo']['Forest Coniferous'])*gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_rswd_norm_1971to2000_' + str(iM+1) + '.tif')['Data'][iMask].astype('float')*meta['Climate']['SF']['rswd']
		nrms['prcp'][iM]=np.maximum(0,gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_prcp_norm_1971to2000_' + str(iM+1) + '.tif')['Data'][0::ivl,0::ivl][iMask].astype('float')*meta['Climate']['SF']['prcp'])
		nrms['vpd'][iM]=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_vpd_norm_1971to2000_' + str(iM+1) + '.tif')['Data'][0::ivl,0::ivl][iMask].astype('float')*meta['Climate']['SF']['vpd']

	# QA - compare BC1ha normals with NA1k normals:
	#iM=7
	#z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_vpd_norm_1971to2000_' + str(iM) + '.tif')['Data'].astype('float')*meta['Climate']['SF']['vpd']
	#iMask=np.where( (zRef['Data']>0) & (z>-90) )
	#np.mean(z[iMask])

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(14,14)) # See if it is working
	cnt=0
	N_run=3
	vi={}
	for iY in range(N_run):
		for iM in range(12):
			print(str(iY) + ' ' + str(iM))
			vi['Month']=iM+1
			vi['LAI']=5.0
			vi['Gs']=0.010
			vi['Ga']=0.058
			vi['tmean']=nrms['tmean'][iM]
			vi['prcp']=nrms['prcp'][iM]
			vi['vpd']=nrms['vpd'][iM]
			vi['rswn']=nrms['rswn'][iM]
			vo=gaia.WBM_SparseGrid(par,vi)

			# Populate inputs for next time step
			vi['ws']=vo['ws']
			vi['wsp']=vo['wsp']

			vo['cwd']=vo['etp']-vo['eta'] # Add climatic water deficit
			vo['cmi']=vi['prcp']-vo['etp'] # Add climatic moisture index

			ax[0,0].plot(cnt,np.mean(vo['etp']),'b.')
			ax[0,0].plot(cnt,np.mean(vo['eta']),'go')
			ax[0,1].plot(cnt,np.mean(vo['melt']),'b.')
			ax[1,0].plot(cnt,np.mean(vo['wsp']),'b.')
			ax[1,1].plot(cnt,np.mean(vo['ws']),'b.')
			cnt=cnt+1
			#print('Year:' + str(meta['tvm'][iT,0]) + ', Month:' + str(meta['tvm'][iT,1]) )

			if iY==N_run-1:
				for k in vo.keys():
					z1=copy.deepcopy(zRef)
					z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
					z1['Data'][iMask]=vo[k]
					#ind=np.where(z1['Data']>200)
					#z1['Data'][ind]=0
					z1['Data'][iMask]=z1['Data'][iMask]/meta['Climate']['SF'][k]
					z1['Data']=z1['Data'].astype('int16')
					z1['Data']=z1['Data']*zRef['Data']
					# plt.close('all'); plt.matshow(z1['Data'],clim=[0,6])
					gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_' + k + '_norm_1971to2000_' + str(iM+1) + '.tif')

	# Seasonal indices
	mo=np.arange(5,10,1)
	vL=['tmean','prcp','ws','etp','cwd','cmi']
	for v in vL:
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
		for iM in range(mo.size):
			z1['Data']=z1['Data']+gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['Climate']['SF'][v]
		z1['Data']=z1['Data']/mo.size
		z1['Data'][iMask]=z1['Data'][iMask]/meta['Climate']['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.matshow(z1['Data'],clim=[-20,20])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_' + v + '_mjjas_norm_1971to2000.tif')

	mo=np.arange(1,13,1)
	#vL=['tmean','prcp']
	vL=['tmean','prcp','etp','cwd','ws','wsp','cmi','runoff','melt']
	for v in vL:
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='float')
		for iM in range(mo.size):
			z1['Data']=z1['Data']+gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_' + v + '_norm_1971to2000_' + str(mo[iM]) + '.tif')['Data'].astype('float')*meta['Climate']['SF'][v]
		z1['Data']=z1['Data']/mo.size
		#plt.close('all'); plt.matshow(z1['Data'],clim=[0,500])
		z1['Data'][iMask]=z1['Data'][iMask]/meta['Climate']['SF'][v]
		z1['Data']=z1['Data'].astype('int16')
		z1['Data']=z1['Data']*zRef['Data']
		#plt.matshow(z1['Data'],clim=[-20,20])
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_' + v + '_ann_norm_1971to2000.tif')

	return

#%%
def SolarRadiationVsTerrain(meta):
	vList=['lc_comp1_2019','rswd_gs_n','aspect','slope','elev']
	z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
	bw=200; bin=np.arange(0,3000,bw)
	N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	return

#%% 
def ClimateStatsByBGCZone(meta):
	# *** Needs updating ***
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz.tif')
	zBGC['Data']=zBGC['Data'].flatten()
	lutBGC=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\VRI 2023\\becz_lut.xlsx')
	
	zMAT=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
	zMAT['Data']=zMAT['Data'].flatten().astype(float)/10
	
	zWS=gis.OpenGeoTiff(r'C:\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
	zWS['Data']=zWS['Data'].flatten().astype(float)
	
	zSI=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\si.tif')
	zSI['Data']=zSI['Data'].flatten().astype(float)
	
	zA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\proj_age_1.tif')
	zA['Data']=zA['Data'].flatten().astype(float)
	
	lutBGC['MAT']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['WS']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['SI']=np.zeros(lutBGC['VALUE'].size)
	lutBGC['Age']=np.zeros(lutBGC['VALUE'].size)
	for i in range(lutBGC['VALUE'].size):
		ind=np.where( (zBGC['Data']==lutBGC['VALUE'][i]) & (zMAT['Data']>=-50) & (zWS['Data']>=0) & (zWS['Data']<=200) & (zSI['Data']>0) & (zSI['Data']<100) & (zA['Data']>=0) & (zA['Data']<1000) )[0]
		lutBGC['MAT'][i]=np.mean(zMAT['Data'][ind])
		lutBGC['WS'][i]=np.mean(zWS['Data'][ind])
		lutBGC['SI'][i]=np.mean(zSI['Data'][ind])
		lutBGC['Age'][i]=np.mean(zA['Data'][ind])
	
	df=pd.DataFrame(lutBGC)
	df.to_excel(r'C:\Data\BC1ha\Climate\tmp.xlsx')
	return

#%% Plot climate space
def ClimateSpace(meta,vX,vY):
	vX='tmean_ann_n'
	vY='ws_mjjas_n'
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z0=u1ha.Import_Raster(meta,[],['lc_comp1_2019',vX,vY],'Extract Grid')
	pX=np.percentile(z0[vX],[0.25,99.75])
	pY=np.percentile(z0[vY],[0.25,99.75])
	return

#%%
def HoldridgeLifeZones(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	iMask=np.where( (zRef['Data']>0) )

	zMAT=np.zeros(iMask[0].size)
	for mo in range(12):
		a=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Monthly\\Normals\\bc1ha_tmean_norm_1971to2000_' + str(mo+1) + '.tif')['Data'][iMask].astype(float)*meta['Climate']['SF']['tmean']
		ind=np.where(a<0)
		a[ind]=0
		zMAT=zMAT+a
	zMAT=zMAT/12

	#zMAT=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_tmean_ann_norm_1971to2000.tif')['Data'][iMask].astype(float)*meta['Climate']['SF']['tmean']
	zMAP=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_prcp_ann_norm_1971to2000.tif')['Data'][iMask].astype(float)
	zPET=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_etp_ann_norm_1971to2000.tif')['Data'][iMask].astype(float)
	zAI=zPET/zMAP
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
	#plt.close('all');plt.matshow(z2['Data'],clim=[0,15])

	gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_HoldridgeLifeZones_1971to2000.tif')

	return
