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
from datetime import datetime
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
import fcgadgets.bc1ha.bc1ha_utils as u1ha

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 1.5 hours
# Raw look-up-table spreadsheets stored with model parameters, while the processed
# pickle files are are stored with bc1ha data

#fiona.listlayers(meta['Paths']['GDB']['LandCover'])
#fiona.listlayers(meta['Paths']['GDB']['LandUse'])
#fiona.listlayers(meta['Paths']['GDB']['VRI 2023'])
#fiona.listlayers(meta['Paths']['GDB']['Disturbance'])
#fiona.listlayers(meta['Paths']['GDB']['Results'])

def BuildLUTsFromSourceGDBs(meta):
	# Unique layers
	uL=np.unique(meta['Geos']['Variable Info']['Layer Name'])

	d={}
	for iL in range(uL.size):

		#if uL[iL]=='VEG_BURN_SEVERITY_SP': # 'BEC_NATURAL_DISTURBANCE_SV':
		#	break
		#if uL[iL]=='PEST_INFESTATION_POLY': #
		#	break
		#if uL[iL]=='VEG_COMP_LYR_R1_POLY': #
		#	break

		t_start=time.time()

		# Index to all variables from layer
		indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==uL[iL]) & (meta['Geos']['Variable Info']['LUT Required']=='Yes') )[0]
		#indL=np.where( (meta['Geos']['Variable Info']['Variable Name']=='FOR_MGMT_LAND_BASE_IND') & (meta['Geos']['Variable Info']['LUT Required']=='Yes') )[0]

		if indL.size==0:
			# Not categorical variables
			continue

		# Path to file
		pthL=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indL][0] ]

		# Initilize dictionary to store variable codes
		d[uL[iL]]={}
		for iV in range(indL.size):
			vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
			d[uL[iL]][vNam]=[]

		with fiona.open(pthL,layer=uL[iL]) as source:
			for feat in source:
				prop=dict(feat['properties'].items())

				for iV in range(indL.size):
					vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
					# Only continue if the variable is populated
					if prop[vNam]==None:
						continue
					# Add to code list
					d[uL[iL]][vNam].append(prop[vNam])

		t_ela=time.time()-t_start
		print(t_ela)

		# Create LUT for each variable with categorical data
		lut={}
		for iV in range(indL.size):
			vNam=meta['Geos']['Variable Info']['Variable Name'][indL[iV]]
			u=np.unique(d[uL[iL]][vNam])
			lut[vNam]={}
			for iU in range(u.size):
				lut[vNam][u[iU]]=iU+1#,dtype=meta['Geos']['Variable Info']['Precision'][indL[iV]])

		# Save
		gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_' + uL[iL] + '.pkl',lut)

	#--------------------------------------------------------------------------
	# Standardize species codes
	#--------------------------------------------------------------------------

	# Import data
	lut_vri=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_VEG_COMP_LYR_R1_POLY.pkl')
	lut_fci=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl')
	lut_pl=gu.ipickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_PLANTING_SVW.pkl')

	# Add everything to a list
	cd=list(lut_vri['SPECIES_CD_1'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_2'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_3'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_4'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_5'].keys())
	cd=cd+list(lut_vri['SPECIES_CD_6'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_1'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_2'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_3'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_4'].keys())
	cd=cd+list(lut_fci['I_SPECIES_CODE_5'].keys())
	cd=cd+list(lut_pl['SILV_TREE_SPECIES_CODE'].keys())

	# Get unique list
	uCode=np.unique(cd)

	lut_vri['SPECIES_CD_1']={}
	lut_vri['SPECIES_CD_2']={}
	lut_vri['SPECIES_CD_3']={}
	lut_vri['SPECIES_CD_4']={}
	lut_vri['SPECIES_CD_5']={}
	lut_vri['SPECIES_CD_6']={}
	for i in range(len(uCode)):
		lut_vri['SPECIES_CD_1'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_2'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_3'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_4'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_5'][uCode[i]]=i+1
		lut_vri['SPECIES_CD_6'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_VEG_COMP_LYR_R1_POLY.pkl',lut_vri)

	lut_fci['I_SPECIES_CODE_1']={}
	lut_fci['I_SPECIES_CODE_2']={}
	lut_fci['I_SPECIES_CODE_3']={}
	lut_fci['I_SPECIES_CODE_4']={}
	lut_fci['I_SPECIES_CODE_5']={}
	for i in range(len(uCode)):
		lut_fci['I_SPECIES_CODE_1'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_2'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_3'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_4'][uCode[i]]=i+1
		lut_fci['I_SPECIES_CODE_5'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_FOREST_COVER_INV_SVW.pkl',lut_fci)

	lut_pl['SILV_TREE_SPECIES_CODE']={}
	for i in range(len(uCode)):
		lut_pl['SILV_TREE_SPECIES_CODE'][uCode[i]]=i+1
	gu.opickle(meta['Paths']['bc1ha'] + '\\LUTs\\LUTs_RSLT_PLANTING_SVW.pkl',lut_pl)

	# Export excel versions of each LUT
	for k1 in meta['LUT'].keys():
		if (k1=='Derived') | (k1=='Raw') | (k1=='Event') | (k1=='Region') | (k1=='TIPSY'):
			continue
		for k2 in meta['LUT'][k1].keys():
			df=pd.DataFrame({'ID':meta['LUT'][k1][k2].values(),'Name':meta['LUT'][k1][k2].keys()})
			df.to_excel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_' + k1 + '_' + k2 + '.xlsx',index=False)

	return

#%% Rasterize VRI
#fiona.listlayers(meta['Paths']['GDB']['Disturbance'])#%% Rasterize VRI
# *** Takes 7.5 hours ***
def RasterizeVRI(meta,year):

	# Pick year
	year=2025
	#year=2015
	t0=time.time()
	lNam='VEG_COMP_LYR_R1_POLY'

	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Feature ID
	flg=0
	if flg==1:
		# *** Manually do Feature to Raster conversion in ArcGIS (rasterizing once
		# took > 12 hours, next time only took 1 hours) So if it is taking more
		# than 2 hours something is wrong. ***

		# Clip
		fin=meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\vri_feaid_orig.tif'
		fout=meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\vri_feaid.tif'
		gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])
		#z=gis.OpenGeoTiff(fout)

	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\vri_feaid.tif')
	fid=zFID['Data'].flatten()
	iu=gu.IndicesFromUniqueArrayValues(fid)

	# Index to variables
	indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]

	# Just specific variables
	#vL=['WHOLE_STEM_BIOMASS_PER_HA','BRANCH_BIOMASS_PER_HA','FOLIAGE_BIOMASS_PER_HA','BARK_BIOMASS_PER_HA']
	#indL=np.where( (np.isin(meta['Geos']['Variable Info']['Variable Name'],vL)==True) )[0]
	#indL=np.where( (meta['Geos']['Variable Info']['Variable Name']=='FOR_MGMT_LAND_BASE_IND') )[0]

	# Initialize variables
	d={}
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		Prec=meta['Geos']['Variable Info']['Precision'][i]
		Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]
		if Date=='Yes':
			d[vNam + '_Year']=np.zeros(fid.size,dtype=Prec)
			d[vNam + '_Month']=np.zeros(fid.size,dtype=Prec)
		d[vNam]=np.zeros(fid.size,dtype=Prec)

	c=fiona.open(meta['Paths']['GDB']['VRI ' + str(year)],layer=lNam)
	df=gpd.GeoDataFrame.from_features(c[10000:20000])

	# Keep track of instances where there is no crosswalk between GDB and rasterized feature ID
	N_Missing=0
	with fiona.open(meta['Paths']['GDB']['VRI ' + str(year)],layer=lNam) as source:
		for feat in source:
			prp=dict(feat['properties'].items())

			if prp['BCLCS_LEVEL_2']!='T':
				continue

			try:
				ind=iu[ prp['FEATURE_ID'] ]
			except:
				N_Missing=N_Missing+1
				continue

			if 'SPECIES_PCT_1' not in prp:
				print(prp['FEATURE_ID'])
				continue

			for i in indL:
				vNam=meta['Geos']['Variable Info']['Variable Name'][i]
				Cat=meta['Geos']['Variable Info']['LUT Required'][i]
				Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]

				if vNam not in prp.keys():
					# There was a change in 2025 and many SI values are Null now
					#print(vNam)
					continue

				if prp[vNam]==None:
					continue

				if Date=='Yes':
					d[vNam + '_Year'][ind]=int(prp[vNam][0:4])
					d[vNam + '_Month'][ind]=int(prp[vNam][5:7])
				else:
					if Cat=='Yes':
						# Categorical, requires LUT
						try:
							d[vNam][ind]=meta['LUT'][lNam][vNam][ prp[vNam] ]
						except:
							continue
					else:
						# Numerical
						d[vNam][ind]=prp[vNam]

	# Save
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		z=copy.deepcopy(zRef)
		z['Data']=np.reshape(d[vNam],zRef['Data'].shape)
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\VRI ' + str(year) + '\\' + vNam + '.tif')

	print((time.time()-t0)/60/60)

	return

#%% Rasterize VRI 2002
# *** Takes 7.5 hours ***
def RasterizeVRI_2002(meta):
	t0=time.time()
	lNam='VEG_COMP_LYR_R1_POLY'

	# Import reference grid
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import feature ID
	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2002\\featid.tif')
	fid=zFID['Data'].flatten()
	iu=gu.IndicesFromUniqueArrayValues(fid)

	# Index to variables
	indL=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) )[0]

	# Initialize variables
	d={}
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		Prec=meta['Geos']['Variable Info']['Precision'][i]
		Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]
		if Date=='Yes':
			d[vNam + '_Year']=np.zeros(fid.size,dtype=Prec)
			d[vNam + '_Month']=np.zeros(fid.size,dtype=Prec)
		d[vNam]=np.zeros(fid.size,dtype=Prec)

	# Keep track of instances where there is no crosswalk between GDB and rasterized feature ID
	cn=0

	with fiona.open(meta['Paths']['GDB']['VRI 2002'],layer='VEG_COMP_LYR_R1_POLY_FINALV4') as source:
		for feat in source:
			prp=dict(feat['properties'].items())
			
			# if prp['BCLCS_LEVEL_2']!='T':
			#	 continue

			try:
				ind=iu[ prp['FEATURE_ID'] ]
			except:
				cn=cn+1
				print(cn)

			for i in indL:
				vNam=meta['Geos']['Variable Info']['Variable Name'][i]
				Cat=meta['Geos']['Variable Info']['LUT Required'][i]
				Date=meta['Geos']['Variable Info']['Date Conversion Required'][i]

				if vNam not in prp.keys():
					continue

				if prp[vNam]==None:
					continue

				if Date=='Yes':
					d[vNam + '_Year'][ind]=int(prp[vNam][0:4])
					d[vNam + '_Month'][ind]=int(prp[vNam][5:7])
				else:
					if Cat=='Yes':
						# Categorical, requires LUT
						try:
							d[vNam][ind]=meta['LUT'][lNam][vNam][ prp[vNam] ]
						except:
							# Older 2002 elements missing from 2023 LUT elements
							pass
					else:
						# Numerical
						d[vNam][ind]=prp[vNam]

	# Save
	for i in indL:
		vNam=meta['Geos']['Variable Info']['Variable Name'][i]
		z=copy.deepcopy(zFID)
		z['Data']=np.reshape(d[vNam],zFID['Data'].shape)
		z=gis.ClipToRaster(z,zRef)
		gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\VRI 2002\\' + vNam + '.tif')
		
	print((time.time()-t0)/60/60)

	return

#%%
def RasterizeForestCoverInventory(meta):
	t0=time.time()

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Forest Cover ID
	flg=0
	if flg==1:
		# 1) Convert to raster in ArcGIS
		# 2) Make sure extent is corrected
		fin=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid_orig.tif'
		fout=r'C:\Data\BC1ha\RSLT_FOREST_COVER_INV_SVW\fcid.tif'
		gis.ClipToRaster_ByFile(fin,fout,meta['Paths']['bc1ha Ref Grid'])

	zFID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\fcid.tif')['Data'].flatten()
	uFID=gu.IndicesFromUniqueArrayValues(zFID)
	
	fc={}
	fc['Path']=meta['Paths']['GDB']['Results']
	fc['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(fc['Path'])
	fc['crs']=meta['Geos']['crs']
	fc['Keep Geom']='Off'
	fc['Select Openings']=np.array([])
	fc['SBC']=np.array([])
	fc['STC']=np.array([])
	fc['SMC']=np.array([])
	fc['FSC']=np.array([])
	fc['SOC1']=np.array([])
	fc['ROI']=[]
	fc['gdf']=qgdb.Query_Openings(fc,[])

	# Add categorical variables
	vL=['STOCKING_STATUS_CODE','STOCKING_TYPE_CODE','SILV_RESERVE_CODE','SILV_RESERVE_OBJECTIVE_CODE','TREE_COVER_PATTERN_CODE']
	for v in vL:
		fc['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_FOREST_COVER_INV_SVW',v,fc['gdf'])

	# Add continuous variables
	vL=['I_TOTAL_STEMS_PER_HA','I_TOTAL_WELL_SPACED_STEMS_HA','I_WELL_SPACED_STEMS_PER_HA','I_FREE_GROWING_STEMS_PER_HA','I_CROWN_CLOSURE_PERCENT',
		'REFERENCE_YEAR','OPENING_ID','ID_STOCKING_STATUS_CODE','ID_STOCKING_TYPE_CODE','ID_SILV_RESERVE_CODE','ID_SILV_RESERVE_OBJECTIVE_CODE',
		'ID_TREE_COVER_PATTERN_CODE']
	z={}
	for v in vL:
		z[v]=np.zeros(zFID.size,dtype='int32')
		
	for i in range(fc['gdf']['OPENING_ID'].size):
		print(i)
		try:
			fid=fc['gdf']['FOREST_COVER_ID'][i]
			for v in vL:
				z[v][ uFID[fid] ]=fc['gdf'][v][i]
			#print('Working!')
		except:
			#print('Missing!')
			pass
	
	for v in vL:
		z1=copy.deepcopy(zRef)
		z1['Data']=np.reshape(z[v],zRef['Data'].shape)
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\' + v + '.tif')
		
	return

#%%
def GapFillBGCZ(meta):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z=Import_Raster(meta,[],['lc_comp1_2019','bgcz'])
	ivl=5
	iGap=np.where( (zRef['Data']==1) & (z['bgcz']['Data']==0) )
	iCal=np.where( (zRef['Data'][0::ivl,0::ivl]==1) & (z['bgcz']['Data'][0::ivl,0::ivl]>0) )
	xy=np.column_stack([zRef['X'][0::ivl,0::ivl][iCal],zRef['Y'][0::ivl,0::ivl][iCal]])
	vals=z['bgcz']['Data'][0::ivl,0::ivl][iCal]
	zFill=griddata(xy,vals,(zRef['X'][iGap],zRef['Y'][iGap]),method='nearest')
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	z1['Data']=z['bgcz']['Data']
	z1['Data'][iGap]=zFill.astype('int8')
	z1['Data'][zRef['Data']==0]=0
	#plt.close('all'); plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
	return

#%%
def RasterizeWildfirePerimitersHistorical(meta,zRef,YearLast):

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Rasterize historical polygons
	lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
	vNam='FIRE_YEAR'
	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)
	ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]
	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()

	tv=np.arange(1917,YearLast+1,1)

	for iT in range(tv.size):

		try:
			df0=df[df[vNam]==tv[iT]].copy()
		except:
			continue

		# Month
		date0=df0['FIRE_DATE'].values
		date1=np.zeros(date0.size)
		for iDate in range(date0.size):
			try:
				date1[iDate]=int(str(date0[iDate])[5:7])
			except:
				date1[iDate]=8

		# Day
		date2=np.ones(date0.size)
		for iDate in range(date0.size):
			try:
				date2[iDate]=int(str(date0[iDate])[8:10])
			except:
				date2[iDate]=1

		# Year
		shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		if len(df0)>0:
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		z1=copy.deepcopy(zRef)
		z1['Data']=z0.astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '_Year.tif')

		# DOY
		DOY=np.zeros(date0.size)
		for iDate in range(date0.size):
			try:
				DOY[iDate]=datetime(int(tv[iT]),int(date1[iDate]),int(date2[iDate])).timetuple().tm_yday
			except:
				if date2[iDate]==29:
					DOY[iDate]=datetime(int(tv[iT]),int(date1[iDate]),int(28)).timetuple().tm_yday

		shapes=((geom,value) for geom, value in zip(df0.geometry,DOY))
		z0=np.zeros(zRef['Data'].shape,dtype=float)
		if len(df0)>0:
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
		z1=copy.deepcopy(zRef)
		z1['Data']=z0.astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '_DOY.tif')
	return

#%% Current-year wildfire (Added to historical files)
def RasterizeWildfirePerimitersCurrentYear(meta,YearCurrent):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	pthin=r'C:\Data\Wildfire\Wildfire.gdb'
	# fiona.listlayers(pthin)
	df=gpd.read_file(pthin,layer='PROT_CURRENT_FIRE_POLYS_SP')
	ind=np.where(df['FIRE_YEAR']==YearCurrent)[0]
	df=df.iloc[ind]
	df=df[df.geometry!=None]
	df=df.reset_index()

	# Year
	shapes=((geom,value) for geom, value in zip(df['geometry'],df['FIRE_YEAR']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	ind=np.where(burned>0); z1['Data'][ind]=1
	plt.close(); plt.matshow(z1['Data'])
	print(np.sum(z1['Data']))
	#gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\prot_current_fire_polys\prot_current_fire_polys.tif')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(YearCurrent) + '_Year.tif')

	# DOY
	# *** Current year has no date! ***
	z1=copy.deepcopy(zRef)
	ind=np.where(zRef['Data']>0); z1['Data'][ind]=212
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(YearCurrent) + '_DOY.tif')

	return

#%%
def DeriveWildfireComposite(meta,zRef,YearCurrent):
	# Mostly from perimiters dataset (historical and current year), plus a small
	# amount from RESULTS ATU layer (13056 ha in 2024)

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	tv=np.arange(1917,YearCurrent,1)

	# Get current perimeters year
	# *** No longer needed because it is added to historical files ***
	#zC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(YearCurrent) + '.tif')

	# Import fire from RESULTS ATU layer
	zR_Y=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_B_Year.tif')['Data']
	zR_DOY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_B_DOY.tif')['Data']

	# Initialize compacted rasters
	N_Year=6
	z={'Year':{},'DOY':{}}
	for iY in range(N_Year):
		z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')
		z['DOY'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	# Initialize summary masks
	zYearFirst=copy.deepcopy(zRef)
	zYearFirst['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zYearLast=copy.deepcopy(zRef)
	zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	zCount=copy.deepcopy(zRef)
	zCount['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):

		try:
			zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '_Year.tif')['Data']
			zDOY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '_DOY.tif')['Data']
		except:
			pass

		# Add any events from RESULTS to the perimiters
		ind=np.where(zR_Y==tv[iT])
		if ind[0].size>0:
			zY[ind]=tv[iT]
			zDOY[ind]=zR_DOY[ind]

		# Update summary variables
		ind=np.where( (zYearFirst['Data']==0) & (zY>0) )
		zYearFirst['Data'][ind]=tv[iT]

		ind=np.where( (zY>0) )
		zYearLast['Data'][ind]=tv[iT]
		zCount['Data'][ind]=zCount['Data'][ind]+1

		# Populate compacted
		ind=np.where((z['Year'][1]==0) & (zY!=0))
		z['Year'][1][ind]=tv[iT]
		z['DOY'][1][ind]=zDOY[ind]
		ind=np.where((z['Year'][1]!=0) & (z['Year'][1]!=tv[iT]) & (z['Year'][2]==0) & (zY!=0))
		z['Year'][2][ind]=tv[iT]
		z['DOY'][2][ind]=zDOY[ind]
		ind=np.where((z['Year'][2]!=0) & (z['Year'][2]!=tv[iT]) & (z['Year'][3]==0) & (zY!=0))
		z['Year'][3][ind]=tv[iT]
		z['DOY'][3][ind]=zDOY[ind]
		ind=np.where((z['Year'][3]!=0) & (z['Year'][3]!=tv[iT]) & (z['Year'][4]==0) & (zY!=0))
		z['Year'][4][ind]=tv[iT]
		z['DOY'][4][ind]=zDOY[ind]
		ind=np.where((z['Year'][4]!=0) & (z['Year'][4]!=tv[iT]) & (z['Year'][5]==0) & (zY!=0))
		z['Year'][5][ind]=tv[iT]
		z['DOY'][5][ind]=zDOY[ind]
		ind=np.where((z['Year'][5]!=0) & (z['Year'][5]!=tv[iT]) & (zY!=0))
		z['Year'][6][ind]=tv[iT]
		z['DOY'][6][ind]=zDOY[ind]

	# Save year of first occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=zYearFirst['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YearFirst.tif')

	# Save year of last occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=zYearLast['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YearLast.tif')

	# Save count
	z1=copy.deepcopy(zRef)
	z1['Data']=zCount['Data'].astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_Count.tif')

	# Save mask of occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	ind0=np.where(zYearLast['Data']>0)
	z1['Data'][ind0]=1
	# plt.matshow(z1['Data'])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_MaskAll.tif')

	# Save compacted layers
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=z['DOY'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_' + str(iY+1) + '_DOY.tif')

	# Plot time series to confirm it worked
	flg=0
	if flg==1:
		lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
		vNam='WildfireComp1'
		tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
		plt.plot(tv,N,'-bo')
	return

#%%
def RasterizeInsects(meta,zRef,YearLast):

	#YearLast=2023

	lNam='PEST_INFESTATION_POLY'
	vNam='PEST_SEVERITY_CODE'

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)

	indI=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]

	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indI[0]] ]

	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()

	df=u1ha.CreateIdForCategoricalVariable(meta,lNam,vNam,df)

	#pestL=['IBM']#,'IBS','IBB','IBD','IDW','IDL']
	pestL=['IBM','IBS','IBB','IBD','IDW','IDL']

	tv=np.arange(1951,YearLast+1,1)

	for pest in pestL:
		# pest='IBM'

		zYearLast=copy.deepcopy(zRef)
		zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

		zYearLastS=copy.deepcopy(zRef)
		zYearLastS['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

		for iT in range(tv.size):
			# iT=tv.size-12
			df0=df[ (df['PEST_SPECIES_CODE']==pest) & (df['CAPTURE_YEAR']==tv[iT]) ].copy()
			df0=df0[df0.geometry!=None]
			df0=df0.reset_index()

			z0=np.zeros(zRef['Data'].shape,dtype=float)
			if len(df0)>0:
				shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID_PEST_SEVERITY_CODE']))
				burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])

			z1=copy.deepcopy(zRef)
			z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][indI[0]])
			gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_' + str(tv[iT]) + '.tif')

			# Update by year grid
			zYearLast['Data'][ (z0==meta['LUT'][lNam][vNam]['M']) | (z0==meta['LUT'][lNam][vNam]['S']) | (z0==meta['LUT'][lNam][vNam]['V']) ]=tv[iT]
			zYearLastS['Data'][ (z0==meta['LUT'][lNam][vNam]['S']) | (z0==meta['LUT'][lNam][vNam]['V']) ]=tv[iT]

		# Year of last occurrence
		z1=copy.deepcopy(zRef)
		z1['Data']=zYearLast['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_YearLast.tif')

		z1=copy.deepcopy(zRef)
		z1['Data']=zYearLastS['Data'].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_YearLast_Severe.tif')

		# Mask of occurrence
		z1=copy.deepcopy(zRef)
		z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
		ind0=np.where(zYearLast['Data']>0)
		z1['Data'][ind0]=1
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + pest + '_MaskAll.tif')

	return

#%% Rasterize planting
def RasterizePlanting(meta,vNam,YearLast):
	
	tv=np.arange(1960,YearLast+1,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	# Start with planting with spatial from RESULTS (takes 15 min)
	flg=1
	if flg==1:
		t0=time.time()
		ats={}
		ats['Path']=meta['Paths']['GDB']['Results']
		ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
		ats['crs']=meta['Geos']['crs']
		ats['Keep Geom']='On'
		ats['Select Openings']=np.array([])
		ats['SBC']=np.array(['PL'])
		ats['STC']=np.array([])
		ats['SMC']=np.array([])
		ats['FSC']=np.array([])
		ats['SOC1']=np.array([])
		ats['ROI']=[]
		ats['gdf']=qgdb.Query_Openings(ats,[])
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()	
		ats['gdf']['Year']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Year'].size):
			ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		ats['gdf']['Month']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Month'].size):
			ats['gdf']['Month'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][5:7])
		ats['gdf']['Day']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Day'].size):
			ats['gdf']['Day'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][8:10])
		ats['gdf']['DOY']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Day'].size):
			ats['gdf']['DOY'][i]=datetime(int(ats['gdf']['Year'][i]),int(ats['gdf']['Month'][i]),int(ats['gdf']['Day'][i])).timetuple().tm_yday
		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
		NumTreesPlanted=ats['gdf']['ACTUAL_PLANTED_NUMBER']
		ats['gdf']['SPH_Planted']=NumTreesPlanted/AreaPlanted
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',ats['gdf'])
		ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_1',ats['gdf'])
		#ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_2',ats['gdf'])
		#ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_3',ats['gdf'])
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		#ats['gdf'].to_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson',driver="GeoJSON")
	else:
		ats={}
		ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats.geojson')

	# Add areas where FC is artificial
	flg=1
	if flg==1:
		at={}
		at['Path']=meta['Paths']['GDB']['Results']
		at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
		at['crs']=meta['Geos']['crs']
		at['Keep Geom']='Off'
		at['Select Openings']=np.array([])
		at['SBC']=np.array(['PL'])
		at['STC']=np.array([])
		at['SMC']=np.array([])
		at['FSC']=np.array([])
		at['SOC1']=np.array([])
		at['ROI']=[]
		at['gdf']=qgdb.Query_Openings(at,[])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',at['gdf'])
		at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_1',at['gdf'])
		#at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_2',at['gdf'])
		#at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_OBJECTIVE_CODE_3',at['gdf'])
		# Make sure to remove entries that we know did not occur (planned or layout)
		ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
		for k in at['gdf'].keys():
			at['gdf'][k]=at['gdf'][k][ikp]
		at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Year'].size):
			at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		at['gdf']['Month']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Month'].size):
			at['gdf']['Month'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][5:7])
		at['gdf']['Day']=np.ones(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(ats['gdf']['Day'].size):
			at['gdf']['Day'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][8:10])
		at['gdf']['DOY']=np.ones(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Day'].size):
			try:
				at['gdf']['DOY'][i]=datetime(int(at['gdf']['Year'][i]),int(at['gdf']['Month'][i]),int(at['gdf']['Day'][i])).timetuple().tm_yday
			except:
				pass
		#gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl',at)
	else:
		at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl')

	#ind=np.where(at['gdf']['Year']==2010)[0]
	#at['gdf']['OPENING_ID'][ind]

	flg=0
	if flg==1:
		A=np.zeros(tv.size)
		N=np.zeros(tv.size)
		for iT in range(tv.size):
			ind=np.where(at['gdf']['Year']==tv[iT])[0]
			A[iT]=np.nansum(at['gdf']['ACTUAL_TREATMENT_AREA'][ind])
			N[iT]=np.nansum(at['gdf']['ACTUAL_PLANTED_NUMBER'][ind])
		plt.close('all'); plt.plot(tv,A,'-bo')
		#plt.close('all'); plt.plot(tv,N/A,'-bo')

	# Import opening ID with spatial
	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']

	zFC_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\OPENING_ID.tif')['Data']
	zFC_STC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_TYPE_CODE.tif')['Data']

	zVRI_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\OPENING_ID.tif')['Data'] 

	# Reduce the size of rasters
	indOP1=np.where( (zOP1!=0) )
	zOP1s=zOP1[indOP1]
	indOP2=np.where( (zOP2!=0) )
	zOP2s=zOP2[indOP2]
	
	indFC=np.where( (zFC_OID!=0) )
	zFC_OIDs=zFC_OID[indFC]
	zFC_STCs=zFC_STC[indFC]
	
	indVRI=np.where( (zVRI_OID!=0) )
	zVRI_OIDs=zVRI_OID[indVRI]
	
	# Unique indices to Opening ID
	uOP1s=gu.IndicesFromUniqueArrayValues(zOP1s)
	uOP2s=gu.IndicesFromUniqueArrayValues(zOP2s)
	uFCs=gu.IndicesFromUniqueArrayValues(zFC_OIDs)
	uVRIs=gu.IndicesFromUniqueArrayValues(zVRI_OIDs)
	
	# Index to planting and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(4):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					'DOY':np.array([]),
					'Source FC':np.array([]),
					'SPH_Planted':np.array([]),
					'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
					'ID_SILV_TECHNIQUE_CODE':np.array([]),
					'ID_SILV_METHOD_CODE':np.array([]),
					'ID_SILV_OBJECTIVE_CODE_1':np.array([])}

	N_MissingArea=0 # 2 entries with missing area
	for iAT in range(at['gdf']['Year'].size):
		print(iAT)

		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		ID=at['gdf']['OPENING_ID'][iAT]
		DOY=at['gdf']['DOY'][iAT]
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]
		STC=at['gdf']['ID_SILV_TECHNIQUE_CODE'][iAT]
		SMC=at['gdf']['ID_SILV_METHOD_CODE'][iAT]
		OC1=at['gdf']['ID_SILV_OBJECTIVE_CODE_1'][iAT]
		#OC2=at['gdf']['ID_SILV_OBJECTIVE_CODE_2'][iAT]
		#OC3=at['gdf']['ID_SILV_OBJECTIVE_CODE_3'][iAT]
		A_Planted=at['gdf']['ACTUAL_TREATMENT_AREA'][iAT]
		NumTreesPlanted=at['gdf']['ACTUAL_PLANTED_NUMBER'][iAT]
		SPH_Planted=NumTreesPlanted/A_Planted
		if np.isnan(A_Planted)==True:
			N_MissingArea=N_MissingArea+1
			continue
	
		iS=0
		flg=1
		try:
			indArt=np.where(zFC_STCs[uFCs[ID]]==meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_TYPE_CODE']['ART'])[0]
		except:
			flg=0
		if flg==1:
			A_Art=indArt.size
			if A_Art>0:
				fA_fc=np.sum(A_Planted)/A_Art
				if (np.abs(fA_fc-1.0)<0.02):
					ind=uFCs[ID][indArt]
					dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
					dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
					dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],1*np.ones(ind.size))
					dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
					dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind.size))
					dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1'],OC1*np.ones(ind.size))
					#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2'],OC2*np.ones(ind.size))
					#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3'],OC3*np.ones(ind.size))
					#print('1')
					continue
		
		iS=1
		flg=1
		try:
			ind=uOP1s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],2*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1'],OC1*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2'],OC2*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3'],OC3*np.ones(ind2.size))
			#print('2')
			continue
		
		iS=2
		flg=1
		try:
			ind=uOP2s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],3*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1'],OC1*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2'],OC2*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3'],OC3*np.ones(ind2.size))
			#print('3')
			continue
		iS=3
		flg=1
		try:
			ind=uVRIs[ID]
			if ind.size==1:
				continue
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind2)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind2.size))
			dP[Year][iS]['Source FC']=np.append(dP[Year][iS]['Source FC'],4*np.ones(ind2.size))
			dP[Year][iS]['SPH_Planted']=np.append(dP[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dP[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_METHOD_CODE']=np.append(dP[Year][iS]['ID_SILV_METHOD_CODE'],STC*np.ones(ind2.size))
			dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_1'],OC1*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_2'],OC2*np.ones(ind2.size))
			#dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3']=np.append(dP[Year][iS]['ID_SILV_OBJECTIVE_CODE_3'],OC3*np.ones(ind2.size))
			print('From VRI')
	 
		#print('Missing')
	#gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl',dP)
	#dP=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dP.pkl')

	# Pack

	# Initialize rasters
	N_Year=6
	vL=['DOY','ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','SPH_Planted','ID_SILV_TECHNIQUE_CODE','ID_SILV_OBJECTIVE_CODE_1']
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zPac={'Year':{},
	   'DOY':{},
	   'ACTIVITY_TREATMENT_UNIT_ID':{},
	   'ID_SILV_FUND_SOURCE_CODE':{},
	   'ID_SILV_TECHNIQUE_CODE':{},
	   'ID_SILV_METHOD_CODE':{},
	   'ID_SILV_OBJECTIVE_CODE_1':{},
	   'SPH_Planted':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		# Add activities without spatial
		iS=0
		iA=indFC[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indFC[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		zYr['ID_SILV_OBJECTIVE_CODE_1'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_1']
		#zYr['ID_SILV_OBJECTIVE_CODE_2'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_2']
		#zYr['ID_SILV_OBJECTIVE_CODE_3'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_3']
		
		iS=1
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		zYr['ID_SILV_OBJECTIVE_CODE_1'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_1']
		#zYr['ID_SILV_OBJECTIVE_CODE_2'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_2']
		#zYr['ID_SILV_OBJECTIVE_CODE_3'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_3']
		
		iS=2
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dP[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		zYr['ID_SILV_METHOD_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_METHOD_CODE']
		zYr['ID_SILV_OBJECTIVE_CODE_1'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_1']
		#zYr['ID_SILV_OBJECTIVE_CODE_2'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_2']
		#zYr['ID_SILV_OBJECTIVE_CODE_3'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_OBJECTIVE_CODE_3']

		# Add activity layer with spatial
		ats0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
		ats0=ats0[ats0.geometry!=None]; #ats0=ats0.reset_index()
		if len(ats0)>0:
			for v in vL:
				shapes=((geom,value) for geom, value in zip(ats0['geometry'],ats0[v]))
				burned=features.rasterize(shapes=shapes,fill=0,out=zYr[v],transform=zRef['Transform'])

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) ) # (zCounter<=1) & 
		zPac['Year'][1][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][1][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][2][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][2][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][3][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][3][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][4][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][4][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][5][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][5][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][6][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][6][ind]=zYr[k][ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['DOY'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_DOY.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_METHOD_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_METHOD_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_OBJECTIVE_CODE_1'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_OBJECTIVE_CODE_1.tif')

	return

#%% Planting (last year)
def DerivePlantingYearLast(meta,zRef):
	zYearLast=copy.deepcopy(zRef)
	zYearLast['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iEY in range(6):
		zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iEY+1) + '_Year.tif')
		zYearLast['Data']=np.maximum(zYearLast['Data'],zYr['Data'])
	gis.SaveGeoTiff(zYearLast,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_YearLast.tif')
	return

#%% Rasterize direct seeding
def RasterizeDirectSeeding(meta):
	
	tv=np.arange(1960,2024,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	
	# Start with planting with spatial from RESULTS (takes 15 min)
	flg=0
	if flg==1:
		t0=time.time()
		ats={}
		ats['Path']=meta['Paths']['GDB']['Results']
		ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
		ats['crs']=meta['Geos']['crs']
		ats['Keep Geom']='On'
		ats['Select Openings']=np.array([])
		ats['SBC']=np.array(['DS'])
		ats['STC']=np.array([])
		ats['SMC']=np.array([])
		ats['FSC']=np.array([])
		ats['SOC1']=np.array([])
		ats['ROI']=[]
		ats['gdf']=qgdb.Query_Openings(ats,[])
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()	
		ats['gdf']['Year']=np.zeros(len(ats['gdf']))
		for i in range(ats['gdf']['Year'].size):
			ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		
		AreaPlanted=ats['gdf']['ACTUAL_TREATMENT_AREA']
		NumTreesPlanted=ats['gdf']['ACTUAL_PLANTED_NUMBER']
		ats['gdf']['SPH_Planted']=NumTreesPlanted/AreaPlanted
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
		ats['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])	
		ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
		ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
		ats['gdf']=ats['gdf'].reset_index()
		#ats['gdf'].to_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson',driver="GeoJSON")
	else:
		ats={}
		ats['gdf']=gpd.read_file(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ats_ds.geojson')

	# Add areas where FC is artificial
	flg=0
	if flg==1:
		at={}
		at['Path']=meta['Paths']['GDB']['Results']
		at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
		at['crs']=meta['Geos']['crs']
		at['Keep Geom']='Off'
		at['Select Openings']=np.array([])
		at['SBC']=np.array(['DS'])
		at['STC']=np.array([])
		at['SMC']=np.array([])
		at['FSC']=np.array([])
		at['SOC1']=np.array([])
		at['ROI']=[]
		at['gdf']=qgdb.Query_Openings(at,[])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
		at['gdf']=CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
		# Make sure to remove entries that we know did not occur (planned or layout)
		ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
		for k in at['gdf'].keys():
			at['gdf'][k]=at['gdf'][k][ikp]
		at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
		for i in range(at['gdf']['Year'].size):
			at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
		gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl',at)
	else:
		at=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at_ds.pkl')

	# Import opening ID with spatial   
	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']

	zFC_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\OPENING_ID.tif')['Data']
	zFC_STC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\ID_STOCKING_TYPE_CODE.tif')['Data']
   
	zVRI_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\OPENING_ID.tif')['Data'] 
   
	# Reduce the size of rasters
	indOP1=np.where( (zOP1!=0) )
	zOP1s=zOP1[indOP1]
	indOP2=np.where( (zOP2!=0) )
	zOP2s=zOP2[indOP2]
	
	indFC=np.where( (zFC_OID!=0) )
	zFC_OIDs=zFC_OID[indFC]
	zFC_STCs=zFC_STC[indFC]
	
	indVRI=np.where( (zVRI_OID!=0) )
	zVRI_OIDs=zVRI_OID[indVRI]	
	
	# Unique indices to Opening ID
	uOP1s=gu.IndicesFromUniqueArrayValues(zOP1s)
	uOP2s=gu.IndicesFromUniqueArrayValues(zOP2s)
	uFCs=gu.IndicesFromUniqueArrayValues(zFC_OIDs)
	uVRIs=gu.IndicesFromUniqueArrayValues(zVRI_OIDs)
	
	# Index to planting and year
	dDS={}
	for iT in range(tv.size):
		dDS[tv[iT]]={}
		for iS in range(4):
			dDS[tv[iT]][iS]={}
			dDS[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'Source FC':np.array([]),
					   'SPH_Planted':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([]),
					   'ID_SILV_TECHNIQUE_CODE':np.array([])}
	
	N_MissingArea=0 # 2 entries with missing area
	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		ID=at['gdf']['OPENING_ID'][iAT]		
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]
		STC=at['gdf']['ID_SILV_TECHNIQUE_CODE'][iAT]
		A_Planted=at['gdf']['ACTUAL_TREATMENT_AREA'][iAT]
		NumTreesPlanted=at['gdf']['ACTUAL_PLANTED_NUMBER'][iAT]
		SPH_Planted=NumTreesPlanted/A_Planted
		if np.isnan(A_Planted)==True:
			N_MissingArea=N_MissingArea+1
			continue
	
		iS=0
		flg=1
		try:
			indArt=np.where(zFC_STCs[uFCs[ID]]==meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_TYPE_CODE']['ART'])[0]
		except:
			flg=0
		if flg==1:
			A_Art=indArt.size
			if A_Art>0:
				fA_fc=np.sum(A_Planted)/A_Art
				if (np.abs(fA_fc-1.0)<0.02):
					ind=uFCs[ID][indArt]
					dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind)
					dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],1*np.ones(ind.size))
					dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind.size))
					dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
					dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
					dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind.size))
					#print('1')
					continue
		
		iS=1
		flg=1
		try:
			ind=uOP1s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],2*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			#print('2')
			continue
		
		iS=2
		flg=1
		try:			
			ind=uOP2s[ID]
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],3*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			#print('3')
			continue
		iS=3
		flg=1
		try:
			ind=uVRIs[ID]
			if ind.size==1:
				continue
		except:
			flg=0
		if flg==1:
			ind2=ind[0:int(np.minimum(ind.size,np.round(A_Planted)))]
			dDS[Year][iS]['IndexToGrid']=np.append(dDS[Year][iS]['IndexToGrid'],ind2)
			dDS[Year][iS]['Source FC']=np.append(dDS[Year][iS]['Source FC'],4*np.ones(ind2.size))
			dDS[Year][iS]['SPH_Planted']=np.append(dDS[Year][iS]['SPH_Planted'],SPH_Planted*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dDS[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind2.size))
			dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dDS[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind2.size))
			dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE']=np.append(dDS[Year][iS]['ID_SILV_TECHNIQUE_CODE'],STC*np.ones(ind2.size))
			print('From VRI')
	 
		#print('Missing')
	gu.opickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl',dDS)
	#dDS=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\dDS.pkl')

	# Pack

	# Initialize rasters
	N_Year=6
	vL=['ACTIVITY_TREATMENT_UNIT_ID','ID_SILV_FUND_SOURCE_CODE','SPH_Planted','ID_SILV_TECHNIQUE_CODE']
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])	
	zPac={'Year':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{},'SPH_Planted':{},'ID_SILV_TECHNIQUE_CODE':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	vNam='DS_All'
	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		# Add activities without spatial
		iS=0
		iA=indFC[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indFC[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=1
		iA=indOP1[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']
		
		iS=2
		iA=indOP2[0][dDS[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dDS[tv[iT]][iS]['IndexToGrid']]
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dDS[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']
		zYr['SPH_Planted'][ iA,iB  ]=dDS[tv[iT]][iS]['SPH_Planted']
		zYr['ID_SILV_TECHNIQUE_CODE'][ iA,iB  ]=dDS[tv[iT]][iS]['ID_SILV_TECHNIQUE_CODE']

		# Add activity layer with spatial
		ats0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
		ats0=ats0[ats0.geometry!=None]; #ats0=ats0.reset_index()
		if len(ats0)>0:			
			for v in vL:
				shapes=((geom,value) for geom, value in zip(ats0['geometry'],ats0[v]))
				burned=features.rasterize(shapes=shapes,fill=0,out=zYr[v],transform=zRef['Transform'])

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) ) # (zCounter<=1) & 
		zPac['Year'][1][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][1][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][2][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][2][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][3][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][3][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][4][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][4][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][5][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][5][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][6][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][6][ind]=zYr[k][ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['SPH_Planted'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SPH_Planted.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_TECHNIQUE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_TECHNIQUE_CODE.tif')
	
	return

#%%
def RasterizeBTKApplication(meta):

	YearLast=2023
	tv=np.arange(1987,YearLast+1,1)

	vNam='PCCA-IDW'

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Import historical record
	gdf_spray=gpd.read_file(r'C:\Data\Aerial Btk Spray\Processed\btk_spray_comp.geojson')

	# Import from RESULTS (started tracking in 2023)
	ats={}
	ats['Path']=meta['Paths']['GDB']['Results']
	ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
	ats['crs']=meta['Geos']['crs']
	ats['Keep Geom']='On'
	ats['Select Openings']=np.array([])
	ats['SBC']=np.array(['PC'])
	ats['STC']=np.array(['CA'])
	ats['SMC']=np.array([])
	ats['FSC']=np.array([])
	ats['SOC1']=np.array([])
	ats['ROI']=[]
	ats['gdf']=qgdb.Query_Openings(ats,[])
	ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
	ats['gdf']=ats['gdf'].reset_index()	
	ats['gdf']['Year']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Year'].size):
		ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
	ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',ats['gdf'])
	ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',ats['gdf'])
	ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',ats['gdf'])
	ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',ats['gdf'])
	ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
	ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
	ats['gdf']=ats['gdf'].reset_index()

	# Initialize packed rasters
	N_Year=7
	zPac={'Year':{},'ID_PEST':{},'ID_SILV_FUND_SOURCE_CODE':{}}
	for k in zPac.keys():
		for iY in range(N_Year):
			zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		
		zYr=np.zeros(zRef['Data'].shape,dtype='int16')
		zFSC=np.zeros(zRef['Data'].shape,dtype='int16')
		z0=np.zeros(zRef['Data'].shape,dtype=float)

		# Add from historical compilation
		df0=gdf_spray[ (gdf_spray['Year']==tv[iT]) ].copy()
		df0=df0[df0.geometry!=None]
		df0['Dummy']=1
		if len(df0)>0:
			shapes=((geom,value) for geom, value in zip(df0.geometry,df0['Dummy']))
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
			zYr=burned.astype('int16')
			zFSC[burned>0]=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE']['FIP']

		# Add from RESULTS
		df0=ats['gdf'][ (ats['gdf']['Year']==tv[iT]) ].copy()
		df0=df0[df0.geometry!=None]
		df0['Dummy']=1
		if len(df0)>0:
			shapes=((geom,value) for geom, value in zip(df0.geometry,df0['Dummy']))
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
			ind=np.where(burned>0)
			zYr[ind]=burned[ind].astype('int16')

			shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID_SILV_FUND_SOURCE_CODE']))
			burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
			zFSC[ind]=burned[ind].astype('int16')

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zFSC!=0) )
		zPac['Year'][1][ind]=tv[iT]
		zPac['ID_PEST'][1][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][1][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zFSC!=0) )
		zPac['Year'][2][ind]=tv[iT]
		zPac['ID_PEST'][2][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][2][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zFSC!=0) )
		zPac['Year'][3][ind]=tv[iT]
		zPac['ID_PEST'][3][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][3][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zFSC!=0) )
		zPac['Year'][4][ind]=tv[iT]
		zPac['ID_PEST'][4][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][4][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zFSC!=0) )
		zPac['Year'][5][ind]=tv[iT]
		zPac['ID_PEST'][5][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][5][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zFSC!=0) )
		zPac['Year'][6][ind]=tv[iT]
		zPac['ID_PEST'][6][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][6][ind]=zFSC[ind]

		ind=np.where( (zPac['Year'][6]!=0) & (zPac['Year'][6]!=tv[iT]) & (zFSC!=0) )
		zPac['Year'][7][ind]=tv[iT]
		zPac['ID_PEST'][7][ind]=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SPECIES_CODE']['IDW']
		zPac['ID_SILV_FUND_SOURCE_CODE'][7][ind]=zFSC[ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')

		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_PEST'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ID_PEST.tif')

		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')

	return

#%%
def RasterizeBurnSeverityCanada(meta):
	# 1) Resample in ArcGIS;
	# 2) Clip to BC in ArcGIS using a BC1ha grid for extent;
	# 3) export geotiff;
	# 4) reproject and clip to ensure BC1ha grid extent (below))
	# https://ftp.maps.canada.ca/pub/nrcan_rncan/Forest-fires_Incendie-de-foret/CanLaBS-Burned_Severity-Severite_des_feux/
	pthin=r'C:\Data\Wildfire\Burn Severity'
	fin=pthin + '\\bsn_yr1.tif'
	z=gis.OpenGeoTiff(fin); print(zRef['Data'].shape); print(z['Data'].shape)
	fout=pthin + '\\bsn_yr2.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout); print(zRef['Data'].shape); print(z['Data'].shape) #z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0); z['Data'][ind]=0
	fout2=meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_Year.tif'
	gis.SaveGeoTiff(z,fout2)
	
	fin=pthin + '\\bsn_dnbr1.tif'
	z=gis.OpenGeoTiff(fin); print(zRef['Data'].shape); print(z['Data'].shape)
	fout=pthin + '\\bsn_dnbr2.tif'
	gis.ReprojectRasterAndClipToRaster(fin,fout,meta['Paths']['bc1ha Ref Grid'],meta['Geos']['crs'])
	z=gis.OpenGeoTiff(fout); print(zRef['Data'].shape); print(z['Data'].shape) #z=gis.ClipToRaster(z,zRef)
	ind=np.where(zRef['Data']==0); z['Data'][ind]=0
	fout2=meta['Paths']['bc1ha'] + '\\Disturbances\\BurnSev_Can_dNBR.tif'
	gis.SaveGeoTiff(z,fout2)
	return

#%%
def RasterizeBurnSeverityCurrentYear(meta,yr):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	# fiona.listlayers(pthin)
	# Current year
	#yr=2023
	pthin=r'C:\Data\Wildfire\Wildfire.gdb'
	df=gpd.read_file(pthin,layer='VEG_BURN_SEVERITY_SAME_YR_SP')
	df['ID']=np.zeros(df['FIRE_YEAR'].size)
	for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
		ind=np.where( (df['BURN_SEVERITY_RATING']==k) )
		df['ID'][ind]=meta['LUT']['Derived']['burnsev_comp1'][k]
	ind=np.where(df['FIRE_YEAR']==yr)[0]
	df=df.iloc[ind]
	df=df[df.geometry!=None]
	df=df.reset_index()
	shapes=((geom,value) for geom, value in zip(df['geometry'],df['ID']))
	z0=np.zeros(zRef['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
	z1=copy.deepcopy(zRef)
	z1['Data']=burned.astype('int16')
	plt.close(); plt.matshow(z1['Data'],clim=[0,5])
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\VEG_BURN_SEVERITY_SP\\BURN_SEVERITY_RATING_Current_' + str(yr) + '.tif')
	return

#%% Planting layer (species and genetic worth)
def RasterizePlantingLayer(meta,zRef):
	# Notes:
	# Relies on ACTIVITY_TREATMENT_UNIT_ID as the crosswalk with the planting
	# info from the activity layer -> must run the RasterizePlanting function
	# before this.
	pl={}
	pl['Path']=meta['Paths']['GDB']['Results']
	pl['Layer']='RSLT_PLANTING_SVW'; # fiona.listlayers(at['Path'])
	pl['crs']=meta['Geos']['crs']
	pl['Keep Geom']='Off'
	pl['Select Openings']=np.array([])
	pl['SBC']=np.array([])
	pl['STC']=np.array([])
	pl['SMC']=np.array([])
	pl['FSC']=np.array([])
	pl['SOC1']=np.array([])
	pl['ROI']=[]
	pl['gdf']=qgdb.Query_Openings(pl,[])
	pl['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_PLANTING_SVW','SILV_TREE_SPECIES_CODE',pl['gdf'])
	
	dGW=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_Seedlot_GW.xlsx')
	
	def round_retain_sum(x):
		N=np.round(np.sum(x)).astype(int)
		y=x.astype(int)
		M=np.sum(y)
		K=N-M
		z=y-x
		if K!=0:
			idx=np.argpartition(z,K)[:K]
			y[idx]+=1
		return y
	
	# Get information by ATU ID
	uAT=np.unique(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID'])
	
	plAT={}
	plAT['ACTIVITY_TREATMENT_UNIT_ID']=uAT
	vL=['PL_SPECIES_CD','PL_SPECIES_PCT','PL_SPECIES_GW']
	for v in vL:
		for i in range(6):
			plAT[v + str(i+1)]=np.zeros(uAT.size,dtype='int16')
	
	# Get weighted average genetic worth by AT ID
	N_no_spc=0
	N_no_num_planted=0
	N_sn_not_found=0
	for iU in range(uAT.size):
		iAT=np.where(pl['gdf']['ACTIVITY_TREATMENT_UNIT_ID']==uAT[iU])[0]
	
		cd0=pl['gdf']['ID_SILV_TREE_SPECIES_CODE'][iAT]
		pct0=np.nan_to_num(pl['gdf']['NUMBER_PLANTED'][iAT]/np.sum(pl['gdf']['NUMBER_PLANTED'][iAT])*100)
		sn0=np.nan_to_num(pl['gdf']['SEEDLOT_NUMBER'][iAT])
	
		if np.sum(cd0)==0:
			N_no_spc=N_no_spc+1
	
		if np.sum(pct0)==0:
			N_no_num_planted=N_no_num_planted+1
			pct0=(100/pct0.size)*np.ones(pct0.size)
	
		# Genetic worth
		gw0=np.zeros(pct0.size)
		for j in range(pct0.size):
			ind=np.where(dGW['SEEDLOT_NUMBER']==sn0[j])[0]
			if ind.size!=0:
				gw0[j]=dGW['GENETIC_WORTH_RTNG'][ind[0]]
			else:
				N_sn_not_found=N_sn_not_found+1
	
		# Dissolve to unique species combinations and calculate weighted average genetic worth
		cd1=np.unique(cd0)
		pct1=np.zeros(cd1.size)
		gw1=np.zeros(cd1.size)
		for iCD in range(cd1.size):
			ind=np.where(cd0==cd1[iCD])[0]
			pct1[iCD]=np.sum(pct0[ind])
			gw1[iCD]=np.sum(pct0[ind]*gw0[ind])/np.sum(pct0[ind])

		# Divi
		gw1=np.nan_to_num(gw1)

		# Sort
		iSort=np.argsort(-pct1)
		cd1=cd1[iSort]
		pct1=pct1[iSort]
		gw1=gw1[iSort]

		# Percents will be stored as integers, ensure they will sum to 100
		pct1=np.round(pct1,decimals=0).astype('int16')
		if (np.sum(pct1)!=100) & (pct1.size>1):
			pct1=round_retain_sum(pct1)

		# Populate structure
		if cd1.size>0:
			plAT['PL_SPECIES_CD1'][iU]=cd1[0]
			plAT['PL_SPECIES_PCT1'][iU]=pct1[0]
			plAT['PL_SPECIES_GW1'][iU]=gw1[0]
		if cd1.size>1:
			plAT['PL_SPECIES_CD2'][iU]=cd1[1]
			plAT['PL_SPECIES_PCT2'][iU]=pct1[1]
			plAT['PL_SPECIES_GW2'][iU]=gw1[1]
		if cd1.size>2:
			plAT['PL_SPECIES_CD3'][iU]=cd1[2]
			plAT['PL_SPECIES_PCT3'][iU]=pct1[2]
			plAT['PL_SPECIES_GW3'][iU]=gw1[2]
		if cd1.size>3:
			plAT['PL_SPECIES_CD4'][iU]=cd1[3]
			plAT['PL_SPECIES_PCT4'][iU]=pct1[3]
			plAT['PL_SPECIES_GW4'][iU]=gw1[3]
		if cd1.size>4:
			plAT['PL_SPECIES_CD5'][iU]=cd1[4]
			plAT['PL_SPECIES_PCT5'][iU]=pct1[4]
			plAT['PL_SPECIES_GW5'][iU]=gw1[4]
		if cd1.size>5:
			plAT['PL_SPECIES_CD6'][iU]=cd1[5]
			plAT['PL_SPECIES_PCT6'][iU]=pct1[5]
			plAT['PL_SPECIES_GW6'][iU]=gw1[5]
	
	# Populate packed layers
	N_Year=6
	lNam='RSLT_ACTIVITY_TREATMENT_SVW'
	vNam='PL_All'
	
	for iY in range(N_Year):
		t0=time.time()
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		uAT1=gu.IndicesFromUniqueArrayValues(z['Data'].flatten())
	
		z0={}
		for iS in range(6):
			z0['PL_SPECIES_CD' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
			z0['PL_SPECIES_PCT' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
			z0['PL_SPECIES_GW' + str(iS+1)]=np.zeros(zRef['Data'].shape,dtype='int16').flatten()
	
		for k in uAT1.keys():
			ind=np.where(plAT['ACTIVITY_TREATMENT_UNIT_ID']==k)[0]
			if ind.size==0:
				continue
			iIndex=uAT1[k]
			for iS in range(6):
				z0['PL_SPECIES_CD' + str(iS+1)][iIndex]=plAT['PL_SPECIES_CD' + str(iS+1)][ind]
				z0['PL_SPECIES_PCT' + str(iS+1)][iIndex]=plAT['PL_SPECIES_PCT' + str(iS+1)][ind]
				z0['PL_SPECIES_GW' + str(iS+1)][iIndex]=plAT['PL_SPECIES_GW' + str(iS+1)][ind]
	
		for iS in range(6):
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_CD' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_CD' + str(iS+1) + '.tif')
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_PCT' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_PCT' + str(iS+1) + '.tif')
			z=copy.deepcopy(zRef)
			z['Data']=np.reshape(z0['PL_SPECIES_GW' + str(iS+1)],(zRef['Data'].shape))
			gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(iY+1) + '_PL_SPECIES_GW' + str(iS+1) + '.tif')
	
		print((time.time()-t0)/60)
	return

#%% Rasterize OPENING ID from OPENING LAYER (1 hour)
def RasterizeOpeningID2(meta):
	# Import opening ID with spatial
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	# Rasterize opening ID
	#RasterizeFromSource(meta,zRef,'RSLT_OPENING_SVW','OPENING_ID')
	
	# Open Opening ID 1
	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	iOP1=np.where(zOP1!=0)
	zOP2=zOP1[iOP1]
	
	# Import opening spatial (takes 26 min)
	t0=time.time()
	ops={}
	ops['Path']=meta['Paths']['GDB']['Results']
	ops['Layer']='RSLT_OPENING_SVW'; 
	ops['crs']=meta['Geos']['crs']
	ops['Keep Geom']='On'
	ops['Select Openings']=np.array([])
	ops['SBC']=np.array([])
	ops['STC']=np.array([])
	ops['SMC']=np.array([])
	ops['FSC']=np.array([])
	ops['SOC1']=np.array([])
	ops['ROI']=[]
	ops['gdf']=qgdb.Query_Openings(ops,[])
	ops['gdf']=ops['gdf'][ops['gdf'].geometry!=None]
	ops['gdf']=ops['gdf'].reset_index()   
	print((time.time()-t0)/60)
	
	flg=np.zeros(ops['gdf']['OPENING_ID'].size,dtype=int)
	for iOP in range(ops['gdf']['OPENING_ID'].size):
		ind=np.where(zOP2==ops['gdf']['OPENING_ID'][iOP])[0]
		if ind.size==0:
			flg[iOP]=1
	ikp=np.where(flg==1)[0]
	ops0=ops['gdf'].iloc[ikp].copy()
	ops0=ops0.reset_index()
	
	zOP0=np.zeros(zRef['Data'].shape,dtype=float)
	shapes=((geom,value) for geom, value in zip(ops0['geometry'],ops0['OPENING_ID']))
	burned=features.rasterize(shapes=shapes,fill=0,out=zOP0,transform=zRef['Transform'])
	
	ind=np.where(zOP0!=0)
	zOP2=copy.deepcopy(zRef)
	zOP2['Data']=np.zeros(zRef['Data'].shape,dtype='int32')
	zOP2['Data'][ind]=zOP0[ind]#.astype('int32')
	gis.SaveGeoTiff(zOP2,meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')
	return

#%%
def RasterizeInsectComp1(meta,YearLast):
	#YearLast=2024
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	#meta['LUT']['PEST_INFESTATION_POLY'].keys()
	
	# Import vector data
	lNam='PEST_INFESTATION_POLY'
	vNam='PEST_SEVERITY_CODE'
	if os.path.exists(meta['Paths']['bc1ha'] + '\\' + lNam)==False:
		os.mkdir(meta['Paths']['bc1ha'] + '\\' + lNam)
	indI=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
	pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][indI[0]] ]
	df=gpd.read_file(pthin,layer=lNam)
	df=df[df.geometry!=None]
	df=df.reset_index()
	df=u1ha.CreateIdForCategoricalVariable(meta,lNam,vNam,df)

	# Parameters
	pestL=['IBB','IBD','IBM','IBS','IDW','IDL']
	sevL=['M','S','V']
	tv=np.arange(1951,YearLast+1,1)

	# Initialize rasters
	N_Year=10
	z={'Year':{},'Type':{}}
	for iY in range(N_Year):
		z['Year'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')
		z['Type'][iY+1]=np.zeros(zRef['Data'].shape,dtype='int8')

	for iT in range(tv.size):
		print(tv[iT])
		for pest in pestL:
			for sev in sevL:
				ind=np.where( (meta['Param']['Raw']['InsectComp1']['Insect Code']==pest) & (meta['Param']['Raw']['InsectComp1']['Severity Code']==sev) )[0][0]
				id=meta['Param']['Raw']['InsectComp1']['ID'][ind]
				
				df0=df[ (df['PEST_SPECIES_CODE']==pest) & (df['PEST_SEVERITY_CODE']==sev) & (df['CAPTURE_YEAR']==tv[iT]) ].copy()
				df0=df0[df0.geometry!=None]
				df0=df0.reset_index()
				df0['ID']=np.ones(len(df0['PEST_SPECIES_CODE']))
				if len(df0)>0:
					shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID']))
					z0=np.zeros(zRef['Data'].shape,dtype=float)
					burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=zRef['Transform'])
				else:
					z0=np.zeros(zRef['Data'].shape,dtype=float)

				iY=1
				ind=np.where((z['Year'][iY]==0) & (z0!=0))
				z['Year'][iY][ind]=tv[iT]
				z['Type'][iY][ind]=id
				for iY in range(2,N_Year):
					ind=np.where((z['Year'][iY-1]!=0) & (z['Year'][iY-1]!=tv[iT]) & (z['Year'][iY]==0) & (z0!=0))
					z['Year'][iY][ind]=tv[iT]
					z['Type'][iY][ind]=id

	# Save
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\InsectComp1_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=z['Type'][iY+1].astype('int8')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\InsectComp1_' + str(iY+1) + '_Type.tif')

	return

#%% Digitize the boundary of TSAs (the original is organized at the sub-TSA level)
def DigitizeTSABoundaries(meta):
	zTSA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FADM_TSA\\TSA_NUMBER_DESCRIPTION.tif')
	u=np.unique(zTSA['Data'])
	u=u[u>0]
	gdf=gpd.GeoDataFrame(data=[],columns=['Value','Name','geometry'])
	cnt=0
	for iU in range(u.size):
		ind=np.where(zTSA['Data']==u[iU])
		#z=np.zeros((zTSA['m'],zTSA['n'],3),dtype=np.uint8)
		z=np.zeros(zTSA['Data'].shape,dtype=np.uint8)
		z[ind]=255
		z=np.dstack([z]*3)
		z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale
		cont=cv2.findContours(image=z,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)
		for j in range(len(cont[0])):
			cont_inner=cont[0][j].squeeze()
			if cont_inner.size==2:
				continue
			if cont_inner.shape[0]<3:
				continue
			pointList=[]
			for k in range(len(cont_inner)):
				c=cont_inner[k][0]
				r=cont_inner[k][1]
				x=int(zTSA['X'][0,c])
				y=int(zTSA['Y'][r,0])
				pointList.append(geometry.Point(x,y))
			gdf.loc[cnt,'Value']=int(u[iU])
			gdf.loc[cnt,'Name']=u1ha.lut_n2s(meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION'],u[iU])[0]
			gdf.loc[cnt,'geometry']=geometry.Polygon([[p.x,p.y] for p in pointList])
			cnt=cnt+1
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\tsa.geojson',driver='GeoJSON')
	return

#%%
def Query_RESULTS(meta,year):

	d={}
	d['Path']=meta['Paths']['GDB']['Results']
	d['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(atu['Path'])
	#d['crs']=roi['crs']
	d['Keep Geom']='Off'
	d['Select Openings']=np.array([])
	d['SBC']=np.array([])
	d['STC']=np.array([])
	d['SMC']=np.array([])
	d['FSC']=np.array([])
	d['ROI']=[]
	d['gdf']=qgdb.Query_Openings(d,[])

	# Without Surveys and denudations
	ind=np.where( (d['gdf']['Year']>=2005) & (d['gdf']['SILV_BASE_CODE']!='SU') & (d['gdf']['SILV_BASE_CODE']!='DN') )[0]
	d2={}
	for k in d['gdf'].keys():
		d2[k]=d['gdf'][k][ind]
	df=pd.DataFrame(d2)
	df.to_excel(r'C:\Data\RESULTS Queries\RESULTS_All_2005on_WithoutSurveys.xlsx')

	# By year
	#year=2017
	ind=np.where( (d['gdf']['Year']==year) & (d['gdf']['SILV_BASE_CODE']!='SU') & (d['gdf']['SILV_BASE_CODE']!='DN') )[0]
	dY={}
	for k in d['gdf'].keys():
		dY[k]=d['gdf'][k][ind]
	df=pd.DataFrame(dY)
	df.to_excel(r'C:\Data\RESULTS Queries\RESULTS_' + str(year) + '.xlsx')

	return

#%%
def RasterizePiling(meta,vNam,YearLast):
	# Differs from RasterizeSilviculture because the spatial is tied to cutblocks

	tv=np.arange(1960,YearLast+1,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	at={}
	at['Path']=meta['Paths']['GDB']['Results']
	at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
	at['crs']=meta['Geos']['crs']
	at['Keep Geom']='Off'
	at['Select Openings']=np.array([])
	at['SBC']=np.array(['SP'])
	at['STC']=np.array(['BU'])
	at['SMC']=np.array(['PILE','RPILE'])
	at['FSC']=np.array([])
	at['SOC1']=np.array([])
	at['ROI']=[]
	at['gdf']=qgdb.Query_Openings(at,[])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',at['gdf'])
	# Make sure to remove entries that we know did not occur (planned or layout)
	ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
	for k in at['gdf'].keys():
		at['gdf'][k]=at['gdf'][k][ikp]
	at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Year'].size):
		at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
	at['gdf']['Month']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Month'].size):
		at['gdf']['Month'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][5:7])
	at['gdf']['Day']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Day'].size):
		at['gdf']['Day'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][8:10])
	at['gdf']['DOY']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Day'].size):
		try:
			at['gdf']['DOY'][i]=datetime(int(at['gdf']['Year'][i]),int(at['gdf']['Month'][i]),int(at['gdf']['Day'][i])).timetuple().tm_yday
		except:
			pass

	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']
	zVRI_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\OPENING_ID.tif')['Data']

	# Reduce the size of rasters
	indOP1=np.where( (zOP1!=0) )
	zOP1s=zOP1[indOP1]
	indOP2=np.where( (zOP2!=0) )
	zOP2s=zOP2[indOP2]

	indVRI=np.where( (zVRI_OID!=0) )
	zVRI_OIDs=zVRI_OID[indVRI]
	
	# Unique indices to Opening ID
	uOP1s=gu.IndicesFromUniqueArrayValues(zOP1s)
	uOP2s=gu.IndicesFromUniqueArrayValues(zOP2s)
	uVRIs=gu.IndicesFromUniqueArrayValues(zVRI_OIDs)
	
	# Index to activity and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(3):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'DOY':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([])}

	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		DOY=at['gdf']['DOY'][iAT].astype(int)
		ID=at['gdf']['OPENING_ID'][iAT]
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]

		iS=0
		flg=1
		try:
			ind=uOP1s[ID]
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
			continue
		
		iS=1
		flg=1
		try:
			ind=uOP2s[ID]
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
			continue

		iS=2
		flg=1
		try:
			ind=uVRIs[ID]
			if ind.size==1:
				continue
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))

	#zH=[None]*6
	#for iY in range(6):
	#	zH[iY]=gis.OpenGeoTiff(r'C:\Data\BC1ha\VEG_CONSOLIDATED_CUT_BLOCKS_SP\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')['Data']

	#--------------------------------------------------------------------------
	# Pack
	#--------------------------------------------------------------------------

	# Initialize rasters
	N_Year=6

	zPac={'Year':{},'DOY':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		iS=0
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		iS=1
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		iS=2
		iA=indVRI[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indVRI[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) ) # (zCounter<=1) & 
		zPac['Year'][1][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][1][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][2][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][2][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][3][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][3][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][4][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][4][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][5][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][5][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][6][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][6][ind]=zYr[k][ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['DOY'][iY+1].astype('int8')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_DOY.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')

	# Mask of occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where(z['Data']>0)
		z1['Data'][ind]=1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_MaskAll.tif')

	# Last year
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where( (z['Data']>z1['Data']) )
		z1['Data'][ind]=z['Data'][ind]
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_YearLast.tif')

	return

#%%
def RasterizePileBurn(meta,vNam,YearLast):
	# Differs from RasterizeSilviculture because the spatial is tied to cutblocks

	tv=np.arange(1960,YearLast+1,1)
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	at={}
	at['Path']=meta['Paths']['GDB']['Results']
	at['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(at['Path'])
	at['crs']=meta['Geos']['crs']
	at['Keep Geom']='Off'
	at['Select Openings']=np.array([])
	at['SBC']=np.array(['SP'])
	at['STC']=np.array(['BU'])
	at['SMC']=np.array(['PBURN'])
	at['FSC']=np.array([])
	at['SOC1']=np.array([])
	at['ROI']=[]
	at['gdf']=qgdb.Query_Openings(at,[])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_FUND_SOURCE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE',at['gdf'])
	at['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_METHOD_CODE',at['gdf'])
	# Make sure to remove entries that we know did not occur (planned or layout)
	ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
	for k in at['gdf'].keys():
		at['gdf'][k]=at['gdf'][k][ikp]
	at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Year'].size):
		at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])
	at['gdf']['Month']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Month'].size):
		at['gdf']['Month'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][5:7])
	at['gdf']['Day']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Day'].size):
		at['gdf']['Day'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][8:10])
	at['gdf']['DOY']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
	for i in range(at['gdf']['Day'].size):
		try:
			at['gdf']['DOY'][i]=datetime(int(at['gdf']['Year'][i]),int(at['gdf']['Month'][i]),int(at['gdf']['Day'][i])).timetuple().tm_yday
		except:
			pass

	zOP1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')['Data']
	zOP2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID_2.tif')['Data']
	zVRI_OID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\OPENING_ID.tif')['Data']

	# Reduce the size of rasters
	indOP1=np.where( (zOP1!=0) )
	zOP1s=zOP1[indOP1]
	indOP2=np.where( (zOP2!=0) )
	zOP2s=zOP2[indOP2]

	indVRI=np.where( (zVRI_OID!=0) )
	zVRI_OIDs=zVRI_OID[indVRI]
	
	# Unique indices to Opening ID
	uOP1s=gu.IndicesFromUniqueArrayValues(zOP1s)
	uOP2s=gu.IndicesFromUniqueArrayValues(zOP2s)
	uVRIs=gu.IndicesFromUniqueArrayValues(zVRI_OIDs)
	
	# Index to activity and year
	dP={}
	for iT in range(tv.size):
		dP[tv[iT]]={}
		for iS in range(3):
			dP[tv[iT]][iS]={}
			dP[tv[iT]][iS]={'IndexToGrid':np.array([],dtype=int),
					   'DOY':np.array([]),
					   'ID_SILV_FUND_SOURCE_CODE':np.array([]),
					   'ACTIVITY_TREATMENT_UNIT_ID':np.array([])}

	for iAT in range(at['gdf']['Year'].size):
		print(iAT)
		Year=at['gdf']['Year'][iAT].astype(int)
		if (Year<tv[0]) | (Year>tv[-1]):
			continue
		DOY=at['gdf']['DOY'][iAT].astype(int)
		ID=at['gdf']['OPENING_ID'][iAT]
		FSC=at['gdf']['ID_SILV_FUND_SOURCE_CODE'][iAT]
		ATUID=at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'][iAT]

		iS=0
		flg=1
		try:
			ind=uOP1s[ID]
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
			continue
		
		iS=1
		flg=1
		try:
			ind=uOP2s[ID]
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))
			continue

		iS=2
		flg=1
		try:
			ind=uVRIs[ID]
			if ind.size==1:
				continue
		except:
			flg=0
		if flg==1:
			dP[Year][iS]['IndexToGrid']=np.append(dP[Year][iS]['IndexToGrid'],ind)
			dP[Year][iS]['DOY']=np.append(dP[Year][iS]['DOY'],DOY*np.ones(ind.size))
			dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE']=np.append(dP[Year][iS]['ID_SILV_FUND_SOURCE_CODE'],FSC*np.ones(ind.size))
			dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID']=np.append(dP[Year][iS]['ACTIVITY_TREATMENT_UNIT_ID'],ATUID*np.ones(ind.size))

	#zH=[None]*6
	#for iY in range(6):
	#	zH[iY]=gis.OpenGeoTiff(r'C:\Data\BC1ha\VEG_CONSOLIDATED_CUT_BLOCKS_SP\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')['Data']

	#--------------------------------------------------------------------------
	# Pack
	#--------------------------------------------------------------------------

	# Initialize rasters
	N_Year=6

	zPac={'Year':{},'DOY':{},'ACTIVITY_TREATMENT_UNIT_ID':{},'ID_SILV_FUND_SOURCE_CODE':{}}
	for iY in range(N_Year):
		for k in zPac.keys():
			if k=='ACTIVITY_TREATMENT_UNIT_ID':
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int32')
			else:
				zPac[k][iY+1]=np.zeros(zRef['Data'].shape,dtype='int16')

	for iT in range(tv.size):
		print(tv[iT])
		
		zYr={}
		for k in zPac.keys():
			zYr[k]=np.zeros(zRef['Data'].shape,dtype=float)

		iS=0
		iA=indOP1[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP1[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		iS=1
		iA=indOP2[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indOP2[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		iS=2
		iA=indVRI[0][dP[tv[iT]][iS]['IndexToGrid']]
		iB=indVRI[1][dP[tv[iT]][iS]['IndexToGrid']]
		zYr['DOY'][ iA,iB  ]=dP[tv[iT]][iS]['DOY']
		zYr['ACTIVITY_TREATMENT_UNIT_ID'][ iA,iB  ]=dP[tv[iT]][iS]['ACTIVITY_TREATMENT_UNIT_ID']
		zYr['ID_SILV_FUND_SOURCE_CODE'][ iA,iB  ]=dP[tv[iT]][iS]['ID_SILV_FUND_SOURCE_CODE']

		# Populate packed grids
		ind=np.where( (zPac['Year'][1]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) ) # (zCounter<=1) & 
		zPac['Year'][1][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][1][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][1]!=0) & (zPac['Year'][1]!=tv[iT]) & (zPac['Year'][2]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][2][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][2][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][2]!=0) & (zPac['Year'][2]!=tv[iT]) & (zPac['Year'][3]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][3][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][3][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][3]!=0) & (zPac['Year'][3]!=tv[iT]) & (zPac['Year'][4]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][4][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][4][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][4]!=0) & (zPac['Year'][4]!=tv[iT]) & (zPac['Year'][5]==0) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][5][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][5][ind]=zYr[k][ind]

		ind=np.where( (zPac['Year'][5]!=0) & (zPac['Year'][5]!=tv[iT]) & (zYr['ID_SILV_FUND_SOURCE_CODE']!=0) )
		zPac['Year'][6][ind]=tv[iT]
		for k in zYr.keys():
			if k=='Year':
				continue
			zPac[k][6][ind]=zYr[k][ind]

	# Save to file
	for iY in range(N_Year):
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['Year'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['DOY'][iY+1].astype('int8')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_DOY.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ACTIVITY_TREATMENT_UNIT_ID'][iY+1].astype('int32')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_ACTIVITY_TREATMENT_UNIT_ID.tif')
		z1=copy.deepcopy(zRef)
		z1['Data']=zPac['ID_SILV_FUND_SOURCE_CODE'][iY+1].astype('int16')
		gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_SILV_FUND_SOURCE_CODE.tif')

	# Mask of occurrence
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where(z['Data']>0)
		z1['Data'][ind]=1
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_MaskAll.tif')

	# Last year
	z1=copy.deepcopy(zRef)
	z1['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
	for iY in range(N_Year):
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_' + str(iY+1) + '_Year.tif')
		ind=np.where( (z['Data']>z1['Data']) )
		z1['Data'][ind]=z['Data'][ind]
	gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_YearLast.tif')

	return

#%% Create geojson files by activity for interactive mapping (used in NOSE modelling)
def VectorMapsByActivity_NO(meta):
	# Create geojson files by activity for interactive mapping (used in NOSE modelling)
	ats={}
	ats['Path']=meta['Paths']['GDB']['Results']
	ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
	ats['crs']=meta['Geos']['crs']
	ats['Keep Geom']='On'
	ats['Select Openings']=np.array([])
	ats['SBC']=np.array([])
	ats['STC']=np.array([])
	ats['SMC']=np.array([])
	ats['FSC']=meta['Param']['Raw']['FSC']['NO List Name']
	ats['SOC1']=np.array([])
	ats['ROI']=[]
	ats['gdf']=qgdb.Query_Openings(ats,[])
	
	ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
	ats['gdf']=ats['gdf'].reset_index()
	
	ats['gdf']['Year']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Year'].size):
		ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
	
	ats['gdf']=ats['gdf'][ (ats['gdf']['RESULTS_IND']=='Y') & (ats['gdf']['SILV_METHOD_CODE']!='LAYOT') ].copy()
	ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
	ats['gdf']=ats['gdf'].reset_index()
	
	vList=['OPENING_ID','Year','SILV_BASE_CODE','SILV_TECHNIQUE_CODE', 'SILV_METHOD_CODE',
		   'SILV_OBJECTIVE_CODE_1','SILV_OBJECTIVE_CODE_2', 'SILV_OBJECTIVE_CODE_3',
		   'SILV_FUND_SOURCE_CODE','FIA_PROJECT_ID','ACTUAL_TREATMENT_AREA','geometry']
	
	# Nutrient management (aerial)
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='FE') & (ats['gdf']['SILV_TECHNIQUE_CODE']=='CA')].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(20)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_FE_CA.geojson',driver='GeoJSON')
	
	# Planting
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='PL') & (np.isin(ats['gdf']['SILV_METHOD_CODE'],['LAYOT','CONE','FENCE','BAGS'])==False)].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_PL_PL.geojson',driver='GeoJSON')
	
	# Planting (road rehab)
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='PL') & (ats['gdf']['SILV_TECHNIQUE_CODE']=='RR')].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_PL_RR.geojson',driver='GeoJSON')
	
	# Direct seeding
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='DS')].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_DS.geojson',driver='GeoJSON')

	# Knockdown
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='SP') & (np.isin(ats['gdf']['SILV_METHOD_CODE'],['CABLE','GUARD','HARV','MDOWN','PUSH'])==True)].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_SP_KD.geojson',driver='GeoJSON')

	# Ripping
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='SP') & (ats['gdf']['SILV_TECHNIQUE_CODE']=='ME') & (np.isin(ats['gdf']['SILV_METHOD_CODE'],['DISC','EXCAV','MOUND','LRIP','RRIP','WING'])==True)].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_SP_Rip.geojson',driver='GeoJSON')

	# Prescribed burning
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='SP') & (ats['gdf']['SILV_TECHNIQUE_CODE']=='BU') & (ats['gdf']['SILV_METHOD_CODE']=='BROAD')].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_SP_BU_BROAD.geojson',driver='GeoJSON')
	
	# Pile burning
	gdf=ats['gdf'].loc[(ats['gdf']['SILV_BASE_CODE']=='SP') & (ats['gdf']['SILV_TECHNIQUE_CODE']=='BU') & (ats['gdf']['SILV_METHOD_CODE']=='PBURN')].drop(columns=['level_0','index'])
	gdf=gdf.explode()
	gdf=gdf.reset_index(drop=True)
	gdf.geometry=gdf.geometry.simplify(5)
	gdf=gdf.filter(vList)
	gdf.crs=meta['Geos']['crs']
	gdf.to_file(meta['Paths']['GDB']['GDB'] + '\\atu_SP_BU_PBURN.geojson',driver='GeoJSON')
	return

#%%
def RasterizeDenudationsFromATU(meta):
	meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['DISTURBANCE_CODE']
	
	YearLast=2024
	sbc=np.array(['DN']);stc=np.array([]); smc=np.array([]);soc1=np.array([]);vNam='DN'
	ats={}
	ats['Path']=meta['Paths']['GDB']['Results']
	ats['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(ats['Path'])
	ats['crs']=meta['Geos']['crs']
	ats['Keep Geom']='On'
	ats['Select Openings']=np.array([])
	ats['SBC']=sbc
	ats['STC']=stc
	ats['SMC']=smc
	ats['FSC']=np.array([])
	ats['SOC1']=soc1
	ats['ROI']=[]
	ats['gdf']=qgdb.Query_Openings(ats,[])
	ats['gdf']=ats['gdf'][ats['gdf'].geometry!=None]
	ats['gdf']=ats['gdf'].reset_index()
	ats['gdf']['Year']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Year'].size):
		ats['gdf']['Year'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][0:4])
	ats['gdf']['Month']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Month'].size):
		ats['gdf']['Month'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][5:7])
	ats['gdf']['Day']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Day'].size):
		ats['gdf']['Day'][i]=int(ats['gdf']['ATU_COMPLETION_DATE'][i][8:10])
	ats['gdf']['DOY']=np.zeros(len(ats['gdf']))
	for i in range(ats['gdf']['Day'].size):
		ats['gdf']['DOY'][i]=datetime(int(ats['gdf']['Year'][i]),int(ats['gdf']['Month'][i]),int(ats['gdf']['Day'][i])).timetuple().tm_yday
	
	ats['gdf']=u1ha.CreateIdForCategoricalVariable(meta,'RSLT_ACTIVITY_TREATMENT_SVW','DISTURBANCE_CODE',ats['gdf'])

	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

	d={}
	for k in meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['DISTURBANCE_CODE'].keys():
		# Year
		z1=np.zeros(zRef['Data'].shape,dtype='int16')
		a0=ats['gdf'].loc[ats['gdf']['DISTURBANCE_CODE']==k]
		shapes=((geom,value) for geom, value in zip(a0['geometry'],a0['Year']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z1,transform=zRef['Transform'])
		d[k]=np.sum(z1)
		z2=copy.deepcopy(zRef)
		z2['Data']=z1
		gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_' + k + '_Year.tif')

		# DOY
		z1=np.zeros(zRef['Data'].shape,dtype='int16')
		a0=ats['gdf'].loc[ats['gdf']['DISTURBANCE_CODE']==k]
		shapes=((geom,value) for geom, value in zip(a0['geometry'],a0['DOY']))
		burned=features.rasterize(shapes=shapes,fill=0,out=z1,transform=zRef['Transform'])
		d[k]=np.sum(z1)
		z2=copy.deepcopy(zRef)
		z2['Data']=z1
		gis.SaveGeoTiff(z2,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_' + k + '_DOY.tif')

	flg=0
	if flg==1:
		vL=['B','C','D','E','F','I','L','P','R','S','W']
		d={}
		for v in vL:
			zR=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DisturbanceFromATU_' + v + '_Year.tif')['Data']
			d[v]=np.where(zR>0)[0].size/1000
		d

	return