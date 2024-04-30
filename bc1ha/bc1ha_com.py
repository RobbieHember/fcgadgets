#%% Import modules
import os
import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import pandas as pd
import pyproj
import fiona
import geopandas as gpd
import scipy.io as spio
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Build look up tables (only do this once a year, takes 8 hours)
flg=0
if flg==1:
	u1ha.BuildLUTsFromSourceGDBs(meta)

#%% Major rivers
u1ha.GetRiversMajor(meta)

#%% Simplify geometry of provincial geodatabases (faster useage)
flg=0
if flg==1:
	u1ha.SimplifyProvincialGDBs(meta)

#%% Create land mask for BC
u1ha.GenerateLandMaskBC(meta)

#%% Digitize the boundary of TSAs
# The original is organized at the sub-TSA level so need to aggregate and then digitize
u1ha.DigitizeTSABoundaries(meta)

#%% Rasterize variables from source
# *** VRI and Forest Cover Inventory area too big - crash ***

# prp=u1ha.GetVariablesFromGDB(meta,'VEG_COMP_LYR_R1_POLY')
# prp=u1ha.GetVariablesFromGDB(meta,'OATS_ALR_POLYS')
# prp=u1ha.GetVariablesFromGDB(meta,'NRC_POPULATED_PLACES_1M_SP')
# prp=u1ha.GetVariablesFromGDB(meta,'F_OWN')
u1ha.RasterizeFromSource(meta,zRef,'F_OWN','OWNERSHIP_DESCRIPTION')
u1ha.RasterizeFromSource(meta,zRef,'NRC_POPULATED_PLACES_1M_SP','NAME')
u1ha.RasterizeFromSource(meta,zRef,'BC_MAJOR_CITIES_POINTS_500M','NAME')
u1ha.RasterizeFromSource(meta,zRef,'GBA_TRANSMISSION_LINES_SP','TRANSMISSION_LINE_ID')
u1ha.RasterizeFromSource(meta,zRef,'HSP_MJR_MINES_PERMTTD_AREAS_SP','STATUS_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'FWA_MANMADE_WATERBODIES_POLY','WATERBODY_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'TA_MUNICIPALITIES_SVW','MUNICIPALITY_NAME')
u1ha.RasterizeFromSource(meta,zRef,'OATS_ALR_POLYS','ALR_POLY_ID')
u1ha.RasterizeFromSource(meta,zRef,'DRA_DGTL_ROAD_ATLAS_MPAR_SP','DIGITAL_ROAD_ATLAS_LINE_ID')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_ROAD_SEGMENT_LINES_SVW','LIFE_CYCLE_STATUS_CODE')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_ROAD_SEGMENT_LINES_SVW','RETIREMENT_DATE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','ZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_BIOGEOCLIMATIC_POLY','SUBZONE')
u1ha.RasterizeFromSource(meta,zRef,'BEC_NATURAL_DISTURBANCE_SV','NATURAL_DISTURBANCE_TYPE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','TIMBER_MARK')
u1ha.RasterizeFromSource(meta,zRef,'FTEN_CUT_BLOCK_POLY_SVW','CUT_BLOCK_ID')
u1ha.RasterizeFromSource(meta,zRef,'FWA_WETLANDS_POLY','WATERBODY_TYPE')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','OGSR_TPDA_SYSID')
u1ha.RasterizeFromSource(meta,zRef,'OGSR_TAP_PRIORITY_DEF_AREA_SP','PRIORITY_DEFERRAL_ID')
u1ha.RasterizeFromSource(meta,zRef,'RMP_OGMA_LEGAL_ALL_SVW','LEGAL_OGMA_PROVID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_OPENING_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','OPENING_ID')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_BASE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_ACTIVITY_TREATMENT_SVW','SILV_TECHNIQUE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_RESERVE_SVW','SILV_RESERVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','SILV_RESERVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','SILV_RESERVE_OBJECTIVE_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','STOCKING_TYPE_CODE') # takes 100 min
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','TREE_COVER_PATTERN_CODE')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_TOTAL_STEMS_PER_HA')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_TOTAL_WELL_SPACED_STEMS_HA')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','I_CROWN_CLOSURE_PERCENT')
u1ha.RasterizeFromSource(meta,zRef,'RSLT_FOREST_COVER_INV_SVW','REFERENCE_YEAR')
u1ha.RasterizeFromSource(meta,zRef,'TA_REGIONAL_DISTRICTS_SVW','REGIONAL_DISTRICT_NAME')
u1ha.RasterizeFromSource(meta,zRef,'TA_PROTECTED_LANDS_SV','PROTECTED_LANDS_DESIGNATION')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','FIRE_YEAR')
u1ha.RasterizeFromSource(meta,zRef,'VEG_BURN_SEVERITY_SP','BURN_SEVERITY_RATING')

#%% Generate sparse inputs
# This should speed import by an order of magnitude
u1ha.GenerateSparseInputs(meta,100,'Province')
u1ha.GenerateSparseInputs(meta,200,'Province')
u1ha.GenerateSparseInputs(meta,50,'Province')
u1ha.GenerateSparseInputs(meta,25,'Province')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtPlots')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtCN')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalCoast')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_Eval')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_LUC')
u1ha.GenerateSparseInputs(meta,10,'NOSE')
u1ha.GenerateSparseInputs(meta,10,'BCFCS_NMC')
u1ha.GenerateSparseInputs(meta,1,'TSA_DawsonCreek')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_CWH')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_SBS')

#%% Gap-fill BGC Zone
u1ha.GapFillBGCZ(meta)

#%% Ecozones of Canada
u1ha.RasterizeEcozonesOfCanada(meta)

#%% Prepare CEC Land Use Map 2020
u1ha.PrepareLandUseCEC(meta)

#%% Prepare data from NTEMS (age, land cover, harvest year)
u1ha.PrepareNTEMS(meta)

#%% Derive land cover compilation 1 (year 2019)
u1ha.DeriveLandCoverComp1(meta)

#%% Derive land use compilation 1 (year 2019)
u1ha.DeriveLandUseComp1(meta)

#%% Derive land cover compilation 1 (year 1800)
u1ha.DeriveLandCoverComp1_1800(meta)

#%% Export summary of LC/LU areas (ha)
u1ha.ExportSummaryLCLU(meta)

#%% Derive land use change year for 1800-2019
u1ha.DeriveLandUseChangeYear_1800to2019(meta)

#%% Derive land cover / use compilation 1 (year 2029 scenarios)
u1ha.DeriveLandCoverLandUseComp1_2020to2049_Scenarios(meta)

#%% Derive deforestation mask
u1ha.DeriveLandCoverLandUseComp1_DeforstationMask(meta)

#%% Reclassify NTEMS and VRI according to LC Comp 1
u1ha.ReclassifyNTEMS_LandCover(meta)
u1ha.ReclassifyVRI_LandCover(meta)

#%% Upland-wetland Forest
u1ha.DeriveUplandWetlandForstMask(meta)
#z=u1ha.Import_Raster(meta,[],['upwetf_ntems'])
#gu.CountByCategories(z['upwetf_ntems']['Data'],'Percent')

#%% Rasterize OPENING ID FROM OPENING LAYER (1 hour)
u1ha.RasterizeOpeningID(meta)

#%% Rasterize forest cover attributes
# Use FC ID created in ArcGIS to connect with attributes
u1ha.RasterizeForestCoverInventory(meta)

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
# Need to add this to historical and then re-create the condensed layers
u1ha.RasterizeWildfireCurrentYear(meta)

#%% Rasterize hisorical wildfire
u1ha.RasterizeWildfire(meta,zRef)

#%% Reproject National Burn Severity data
# Step 1) Resample in ArcGIS; 2) Clip to BC in ArcGIS using a BC1ha grid for extent; 
# 3) export geotiff; 4) reproject and clip to ensure BC1ha grid extent (below))
u1ha.ImportBurnSeverityCanada(meta)

#%% Rasterize current-year burn severity rating
u1ha.RasterizeBurnSeverityCurrentYear(meta)

#%% Derive burn severity to pair with wildfire perimeters
u1ha.DeriveBurnSeverityCompilation(meta)

#%% Rasterize insects from AOS
u1ha.RasterizeInsects(meta,zRef)

#%% Rasterize insect compilation from AOS (5 hours)
u1ha.RasterizeInsectComp1(meta)
u1ha.CalcInsectComp1_TimeSeries(meta)

#%% Rasterize consolidated cutblocks database
u1ha.RasterizeHarvest_CC(meta,zRef)

#%% Harvest early reconstruction
u1ha.DeriveHarvestEarlyYears(meta)

#%% Harvest consolidated (with and w/o early reconstruction)
u1ha.DeriveHarvestCompilation(meta)

#%% Harvest silviculture system code
u1ha.DeriveHarvest_SILV_SYSTEM_CODE(meta)

#%% Rasterize planting
u1ha.RasterizePlanting(meta) # Compile location, SBC, STC from various sources
u1ha.DeriveASETComp1(meta) # Derive artificial stand establishment type (ASET)
u1ha.RasterizePlantingLayer(meta) # Get attributes from planting layer (to get species and genetic worth)
u1ha.DerivePlantingStatsByTime(meta) # AIL summary
u1ha.DeriveLastASET(meta) # Derive last instance of Artificial Stand Establishment Type (ASET)
u1ha.MaskPlantingNonOb(meta,zRef) # Derive non-obligation planting mask

#%%
def DerivePlantingStatsByTime(meta):
	# ASET is not included here because a summary time series is created during production
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	iMask=np.where(zRef['Data']==1)
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BEC_ZONE_CODE.tif')['Data'][iMask]

	ds={}
	ds['Year']=np.arange(1960,2023,1)
	ds['Global']={}
	ds['Global']['GW Mean']=np.zeros(ds['Year'].size)
	ds['Global']['SPH Mean']=np.zeros(ds['Year'].size)
	ds['ByFSC']={}
	ds['ByFSC']['FSC']=np.array(list(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'].keys()))
	ds['ByFSC']['Area']=np.zeros((ds['Year'].size,ds['ByFSC']['FSC'].size))

	ds['ByBGC']={}
	ds['ByBGC']['ZONE']=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys()))
	ds['ByBGC']['GW Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	ds['ByBGC']['SPH Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))

	# Need to track average over 6 compacted files and then calculate weighted average at the end
	trackG=[None]*6
	for i in range(6):
		trackG[i]={}
		trackG[i]['GW Mean']=np.zeros(ds['Year'].size)
		trackG[i]['SPH Mean']=np.zeros(ds['Year'].size)
		trackG[i]['N']=np.zeros(ds['Year'].size)
	trackBGC=[None]*6
	for i in range(6):
		trackBGC[i]={}
		trackBGC[i]['GW Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
		trackBGC[i]['SPH Mean']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
		trackBGC[i]['N']=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))

	for i in range(6):
		print(i)
		z={}
		z['yr']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')['Data'][iMask]
		ikp=np.where(z['yr']>0)
		z['yr']=z['yr'][ikp]
		z['BGC']=zBGC[ikp]
		z['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data'][iMask][ikp]
		z['SPH']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SPH_Planted.tif')['Data'][iMask][ikp]
		z['S1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD1.tif')['Data'][iMask][ikp]
		z['P1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT1.tif')['Data'][iMask][ikp]
		z['GW1']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW1.tif')['Data'][iMask][ikp]
		z['S2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD2.tif')['Data'][iMask][ikp]
		z['P2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT2.tif')['Data'][iMask][ikp]
		z['GW2']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW2.tif')['Data'][iMask][ikp]
		z['S3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD3.tif')['Data'][iMask][ikp]
		z['P3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT3.tif')['Data'][iMask][ikp]
		z['GW3']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW3.tif')['Data'][iMask][ikp]
		z['S4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_CD4.tif')['Data'][iMask][ikp]
		z['P4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_PCT4.tif')['Data'][iMask][ikp]
		z['GW4']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_PL_SPECIES_GW4.tif')['Data'][iMask][ikp]

		for iT in range(ds['Year'].size):
			indT=np.where( (z['yr']==ds['Year'][iT]) )[0]

			# By Funding Source
			fsc0=z['FSC'][indT]
			uFSC=np.unique(fsc0)
			uFSC=uFSC[uFSC>0]
			for iFSC in range(uFSC.size):
				ind1=np.where( (fsc0==uFSC[iFSC]) )[0]
				ind2=np.where(ds['ByFSC']['FSC']==u1ha.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],uFSC[iFSC])[0])[0]
				ds['ByFSC']['Area'][iT,ind2]=ds['ByFSC']['Area'][iT,ind2]+ind1.size

			# Global
			f1=z['P1'][indT].astype('float')/100
			f2=z['P2'][indT].astype('float')/100
			f3=z['P3'][indT].astype('float')/100
			f4=z['P4'][indT].astype('float')/100
			gw1=np.nan_to_num(z['GW1'][indT].astype('float'))
			gw2=np.nan_to_num(z['GW2'][indT].astype('float'))
			gw3=np.nan_to_num(z['GW3'][indT].astype('float'))
			gw4=np.nan_to_num(z['GW4'][indT].astype('float'))
			trackG[i]['GW Mean'][iT]=np.nanmean(( (f1*gw1)+(f2*gw2)+(f3*gw3)+(f4*gw4) )/(f1+f2+f3+f4))
			trackG[i]['SPH Mean'][iT]=np.nanmean(z['SPH'][indT].astype('float'))
			trackG[i]['N'][iT]=indT.size

			# By BGC zone
			bgc0=z['BGC'][indT]
			sph0=z['SPH'][indT]
			uBGC=np.unique(bgc0)
			uBGC=uBGC[uBGC>0]
			for iBGC in range(uBGC.size):
				ind1=np.where( (bgc0==uBGC[iBGC]) )[0]
				ind2=np.where(ds['ByBGC']['ZONE']==u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uBGC[iBGC])[0])[0]
				trackBGC[i]['GW Mean'][iT,ind2]=np.mean(( (f1[ind1]*gw1[ind1])+(f2[ind1]*gw2[ind1])+(f3[ind1]*gw3[ind1])+(f4[ind1]*gw4[ind1]) )/(f1[ind1]+f2[ind1]+f3[ind1]+f4[ind1]))
				trackBGC[i]['SPH Mean'][iT,ind2]=np.mean(sph0[ind1].astype('float'))
				trackBGC[i]['N'][iT,ind2]=ind1.size

	# Calculate global mean genetic worth
	gw_sum=np.zeros(ds['Year'].size)
	sph_sum=np.zeros(ds['Year'].size)
	n=np.zeros(ds['Year'].size)
	for i in range(6):
		gw_sum=gw_sum+trackG[i]['N']*trackG[i]['GW Mean']
		sph_sum=sph_sum+trackG[i]['N']*trackG[i]['SPH Mean']
		n=n+trackG[i]['N']
	ds['Global']['GW Mean']=gw_sum/n
	ds['Global']['SPH Mean']=sph_sum/n

	# Calculate mean genetic worth by BGC zone
	gw_sum=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	sph_sum=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	n=np.zeros((ds['Year'].size,ds['ByBGC']['ZONE'].size))
	for i in range(6):
		gw_sum=gw_sum+trackBGC[i]['N']*trackBGC[i]['GW Mean']
		sph_sum=sph_sum+trackBGC[i]['N']*trackBGC[i]['SPH Mean']
		n=n+trackBGC[i]['N']
	ds['ByBGC']['GW Mean']=gw_sum/n
	ds['ByBGC']['SPH Mean']=sph_sum/n

	flg=0
	if flg==1:
		iT=np.where( (ds['Year']>=1985) )[0]
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6.5)); ms=2; lw=0.5
		ax.plot(ds['Year'][iT],ds['Global']['GW Mean'][iT],'-ko',ms=ms,lw=lw,label='Province')
		ind=np.where(ds['ByBGC']['ZONE']=='CWH')[0]
		iGF=np.where(np.isnan(ds['ByBGC']['GW Mean'][iT,ind])==False)[0]
		y=np.interp(ds['Year'][iT],ds['Year'][iT[iGF]],ds['ByBGC']['GW Mean'][iT[iGF],ind])
		ax.plot(ds['Year'][iT],y,'-gs',ms=ms,lw=lw,label='CWH')
		ind=np.where(ds['ByBGC']['ZONE']=='SBS')[0]
		ax.plot(ds['Year'][iT],ds['ByBGC']['GW Mean'][iT,ind],'-c^',ms=ms,lw=lw,label='SBS')
		ax.set(xticks=np.arange(0,3000,5),yticks=np.arange(0,3000,2),ylabel='Average genetic worth (%)',xlabel='Time, years',xlim=[1985,2023],ylim=[0,18])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\GeneticWorth_ts','png',900)

		iT=np.where( (ds['Year']>=1985) )[0]
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6.5)); ms=2; lw=0.5
		ax.plot(ds['Year'][iT],ds['Global']['SPH Mean'][iT],'-ko',ms=ms,lw=lw,label='Province')
		ind=np.where(ds['ByBGC']['ZONE']=='CWH')[0]
		iGF=np.where(np.isnan(ds['ByBGC']['SPH Mean'][iT,ind])==False)[0]
		y=np.interp(ds['Year'][iT],ds['Year'][iT[iGF]],ds['ByBGC']['SPH Mean'][iT[iGF],ind])
		#np.mean(y[-10:])
		ax.plot(ds['Year'][iT],y,'-gs',ms=ms,lw=lw,label='CWH')
		ind=np.where(ds['ByBGC']['ZONE']=='SBS')[0]
		#np.mean(ds['ByBGC']['SPH Mean'][iT,ind][-10:])
		ax.plot(ds['Year'][iT],ds['ByBGC']['SPH Mean'][iT,ind],'-c^',ms=ms,lw=lw,label='SBS')
		ax.set(xticks=np.arange(0,3000,5),yticks=np.arange(0,3000,100),ylabel='Planting density (SPH)',xlabel='Time, years',xlim=[1985,2023],ylim=[0,1600])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(loc='lower right',facecolor=[1,1,1],frameon=False);
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\PlantingDensity_ts','png',900)

	# Remove FSCs with no activity
	ind=np.where(np.sum(ds['ByFSC']['Area'],axis=0)>0)[0]
	ds['Area Unique']=ds['ByFSC']['Area'][:,ind]
	ds['FSC Unique']=ds['ByFSC']['FSC'][ind]

	# Save
	gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PlantingSummaryByTime.pkl',ds)
	return

#%% Rasterize fertilization
u1ha.RasterizeSilviculture(meta,np.array(['FE']),np.array(['CA']),np.array([]),np.array([]),'FE-CA')

#%% Rasterize knockdown
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array([]),np.array(['CABLE','GUARD','HARV','MDOWN','PUSH']),np.array([]),'SP-KD')

#%% Rasterize mechanical site prep
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['ME']),np.array(['CABLE','GUARD','HARV','MDOWN','PUSH']),np.array([]),'SP-Ripping')

#%% Rasterize pile burning
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array([]),np.array([]),'SP-BU')

#%% Rasterize broadcast burning (237 min)
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array(['BROAD']),np.array([]),np.array(['SP-BU-BROAD']))

#%% Percent dead from cruises
u1ha.RasterizeCruisePercentDead(meta)

#%% Salvage harvest mask
u1ha.HarvestSalvageMaskFromCruise(meta)

#%% Create range tenure consolidated
u1ha.DeriveRangeTenureMask(meta)

#%% Create Crown land mask
u1ha.DeriveCrownLandMask(meta)

#%% Create BGC Zone / NDT Combination
u1ha.DeriveBGCZoneNDTCombo(meta)

#%% Consolidate Tree Density Class
# *** Ensure harvested areas are listed as dense ***
u1ha.DeriveTreeDensityClass(meta)
df=u1ha.CalculateStatsByTDC(meta) # Create a table of descriptive stats

#%% Rasterize early land use change year
u1ha.RasterizeEaryLandUseChangeYear(meta)

#%% Rasterize distance from roads
u1ha.DeriveDistanceFromRoads(meta)

#%% Rasterize distance from timber facilities
u1ha.DeriveDistanceFromFacility(meta)

#%% Derive Forest Cover Reserve compilation 1
u1ha.DeriveForestCoverReserveComp1(meta)

#%% Rasterize aerial spray treatment
u1ha.RasterizeBTKSpray(meta)

#%% Global Forest Change Loss Year (mosaic and reproject)
u1ha.RasterizeGFC_LossYear(meta)

#%% Filter Global Forest Change Loss Year (to remove known disturbances)
u1ha.FilterGFC_LossYear(meta)

#%% Climate data
u1ha.ImportNormalsFromClimateNA(meta)
u1ha.CalcSaturationVapourPressureNormal(meta)
u1ha.CalcActualVapourPressureNormalFromTemps(meta)
u1ha.CalcActualVapourPressureNormalBiasCorrected(meta)
u1ha.CalcVapourPressureDeficitNormal(meta)
u1ha.CalcSurfaceWaterBalanceNormals(meta)

#%% Extract mean climate data by BGC zone
u1ha.ClimateStatsByBGCZone(meta)

#%% Plot climate space
vX='tmean_ann_n'
vY='ws_mjjas_n'
def ClimateSpace(meta,vX,vY):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	z0=u1ha.Import_Raster(meta,[],['lc_comp1_2019',vX,vY],'Extract Grid')
	pX=np.percentile(z0[vX],[0.25,99.75])
	pY=np.percentile(z0[vY],[0.25,99.75])
	return

#%% Deciduous fraction
u1ha.DeriveBroadleafDeciduousFraction(meta)

#%% Access zones
u1ha.DeriveAccessZones(meta)

#%% Species groups
z0=u1ha.Import_Raster(meta,[],['refg','lc_comp1_2019','fire_yr','bsr_sc','spc1_vri23','geomorph'],'Extract Grid')

#%%

z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\FromBurtWilliam\\geomorphons_search50.tif')
plt.matshow(z['Data'],clim=[0,11])

#%% QA
flg=0
if flg==1:
	z=u1ha.Import_Raster(meta,[],['refg','lcc1_c','PROJ_AGE_1','harv_yr_con1','fire_yr'])
	yr=2023-z['PROJ_AGE_1']['Data']
	plt.close('all')
	plt.matshow(yr,clim=[1800,2000])

#%% QA
flg=0
if flg==1:
	z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','lu_comp1_2049s1','lu_comp1_2049s2','lu_comp1_2019','lc_comp1_2049s1','lc_comp1_2049s2'],'Extract Grid')

ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)
ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s1']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s1']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s1']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s1']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)
ind=np.where( (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s2']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s2']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) |
			  (z['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lc_comp1_2049s2']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z['lu_comp1_2049s2']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
print(ind[0].size/1e6)



#%%
p=np.zeros(10)
for i in range(10):
	ind1=np.where( (z0['geomorph']==i+1) & (z0['lc_comp1_2019']==1) )
	ind2=np.where( (z0['geomorph']==i+1) & (z0['lc_comp1_2019']==1) & (z0['fire_yr']>0) & (z0['bsr_sc']>=3) )
	p[i]=ind2[0].size/ind1[0].size*100
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p,0.8)

#%%
p=np.zeros((4,10))
for i in range(4):
	ind1=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) )
	tmp=z0['geomorph'][ind1]
	for j in range(10):
		ind2=np.where( (tmp==j+1) )
		p[i,j]=ind2[0].size/ind1[0].size
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p[3,:]/p[0,:],0.8)

fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(10),p[0,:],0.8)

#%%
p=np.zeros(4)
for i in range(4):
	ind1=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) )
	ind2=np.where( (z0['bsr_sc']==i+1) & (z0['fire_yr']>0) & (z0['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['DG']) )
	p[i]=ind2[0].size/ind1[0].size
plt.bar(np.arange(4),p,0.8)

#%%
p=np.zeros(10)
for i in range(1,11):
	ind1=np.where( (z0['geomorph']==i) & (z0['lc_comp1_2019']==1) )
	ind2=np.where( (z0['geomorph']==i) & (z0['lc_comp1_2019']==1) & (z0['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PLI']) )
	p[i-1]=ind2[0].size/ind1[0].size*100
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
plt.bar(np.arange(1,11),p,0.8)

#%% Balsam-Deciduous Compilation
vList=['refg','lc_comp1_2019','spc1_vri23','wshed','wshedgrp']
z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

#c,N=gu.CountByCategories(z0['geomorph'][(z0['geomorph']>0)])
#plt.bar(c,N)
plt.matshow(z0['wshedgrp'])

#%%
z1=np.zeros(zRef['Data'].shape)
ind=np.where( (z0['geomorph']==4) )
z1[ind]=1
plt.matshow(z1)

#%%

vList=['refg','refyear','lc_comp1_2019','fire_2023']
z=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

ind=np.where( (z['fire_2023']>0) & (z['refyear']>0) )
ind1=np.where( (z['fire_2023']>0) & (z['refyear']==0) )
ind2=np.where( (z['fire_2023']>0) )
ind3=np.where( (z['fire_2023']>0) & (z['refyear']>2019) )
ind3[0].size/ind2[0].size*100

u,N=gu.CountByCategories(z['lc_comp1_2019'][ind1],'Percent')

2023-np.mean(z['refyear'][ind])
plt.hist(z['refyear'][ind])