
#%% Import modules

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import copy
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables

meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

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

z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\FromBurtWilliam\\geomorphons_search50.tif')
plt.matshow(z['Data'],clim=[0,11])


#%% Build look up tables (only do this once a year, takes 8 hours)
flg=0
if flg==1:
    u1ha.BuildLUTsFromSourceGDBs(meta)

#%% Simplify geometry of provincial geodatabases (faster useage)
flg=0
if flg==1:
    u1ha.SimplifyProvincialGDBs(meta,gdf)

#%% Create land mask for BC
u1ha.GenerateLandMaskBC(meta)

#%% Digitize the boundary of TSAs
# The original is organized at the sub-TSA level so need to aggregate and then 
# digitize
u1ha.DigitizeTSABoundaries(meta)

#%% Rasterize variables from source
# *** VRI and Forest Cover Inventory area too big - crash ***
# prp=u1ha.GetVariablesFromGDB(meta,'BC_MAJOR_WATERSHEDS')
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
u1ha.RasterizeFromSource(meta,zRef,'BC_MAJOR_WATERSHEDS','MAJOR_WATERSHED_SYSTEM')
u1ha.RasterizeFromSource(meta,zRef,'FWA_WATERSHED_GROUPS_POLY','WATERSHED_GROUP_NAME')
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
u1ha.GenerateSparseInputs(meta,25,'Province')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtPlots')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalAtCN')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_EvalCoast')
u1ha.GenerateSparseInputs(meta,1,'BCFCS_Eval')
u1ha.GenerateSparseInputs(meta,50,'BCFCS_LUC')
u1ha.GenerateSparseInputs(meta,100,'NOSE')
u1ha.GenerateSparseInputs(meta,10,'BCFCS_NMC')

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

#%% Plot land use change time series

flg=0
if flg==1:        
    # Plot time series
    tv=np.arange(1851,2050,1)
    vList=['luc1_hist_yr','luc1_hist_type','luc1_fut_s1_yr','luc1_fut_s1_type','luc1_1019_type','luc1_1019_yr'] # 
    z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
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
    #     ax.bar(tv,ts1019[c]/1000,1,bottom=yC,facecolor=cl[cnt,:])
    #     yC=yC+ts1019[c]/1000
    #     cnt=cnt+1
    ax.set(yscale='linear',ylabel='Deforestation (Kha/year)',xlabel='Time, years',xticks=np.arange(1800,2140,10),ylim=[0,80],xlim=[1860,2050])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    plt.tight_layout()
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\LandCover\LUC_Comp1_TS','png',900)
    
    #        
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

#%% Reclassify NTEMS and VRI according to LC Comp 1
u1ha.ReclassifyNTEMS_LandCover(meta)
u1ha.ReclassifyVRI_LandCover(meta)

#%% Upland-wetland Forest
u1ha.DeriveUplandWetlandForstMask(meta)
#z=u1ha.Import_Raster(meta,[],['upwetf_ntems'])
#gu.CountByCategories(z['upwetf_ntems']['Data'],'Percent')

#%% Rasterize OPENING ID FROM OPENING LAYER (1 hour)
u1ha.RasterizeOpeningID(meta)

#%% Rasterize Forest Cover Attributes
# Use FC ID created in ArcGIS to connect with attributes
u1ha.RasterizeForestCoverInventory(meta)

#%% Rasterize hisorical wildfire
u1ha.RasterizeWildfire(meta)

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
# Need to add this to historical and then re-create the condensed layers
u1ha.RasterizeWildfireCurrentYear(meta)

#%% Reproject National Burn Severity data
# Step 1) Resample in ArcGIS; 2) Clip to BC in ArcGIS using a BC1ha grid for extent; 
# 3) export geotiff; 4) reproject and clip to ensure BC1ha grid extent (below))
u1ha.ImportBurnSeverityCanada(meta)

#%% Derive burn severity to pair with wildfire perimeters
u1ha.DeriveBurnSeverityCompilation(meta)

#%% Rasterize insects from AOS
u1ha.RasterizeInsects(meta,zRef)

#%% Rasterize insect compilation from AOS (5 hours)
u1ha.RasterizeInsectComp1(meta)

#%% Rasterize consolidated cutblocks database
u1ha.RasterizeHarvest_CC(meta,zRef)

#%% Harvest early reconstruction
u1ha.DeriveHarvestEarlyYears(meta)

#%% Harvest consolidated (with and w/o early reconstruction)
u1ha.DeriveHarvestCompilation(meta)

#%% Rasterize planting
u1ha.RasterizePlanting(meta)

#%% Derive non-obligation planting mask
u1ha.MaskPlantingNonOb(meta,zRef)

#%% Derive Regen Type (used for non-ob stand establishment)
u1ha.DeriveRegenTypeCompilation(meta)

#%% Planting layer (to get species and genetic worth)
u1ha.RasterizePlantingLayer(meta)

#%% Rasterize fertilization
u1ha.RasterizeSilviculture(meta,np.array(['FE']),np.array(['CA']),np.array([]),np.array([]),'FE-CA')

#%% Rasterize knockdown
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array([]),np.array(['CABLE','GUARD','HARV','MDOWN','PUSH']),np.array([]),'SP-KD')

#%% Rasterize pile burning
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array([]),np.array([]),'SP-BU')

#%% Rasterize broadcast burning (237 min)
u1ha.RasterizeSilviculture(meta,np.array(['SP']),np.array(['BU']),np.array(['BROAD']),np.array([]),np.array(['SP-BU-BROAD']))

#%% Percent dead from cruises
u1ha.RasterizeCruisePercentDead(meta)

#%% Salvage logging mask
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

#%% Rasterize early land use change year
u1ha.RasterizeEaryLandUseChangeYear(meta)

#%% Rasterize distance from roads
u1ha.DeriveDistanceFromRoads(meta)

#%% Rasterize distance from timber facilities
u1ha.DeriveDistanceFromFacility(meta)

#%% Derive harvest retention compilation 1
u1ha.DeriveHarvestRetentionComp1(meta)

#%% Harvest regeneration type
u1ha.DeriveHarvestRegenType(meta)

#%% Rasterize aerial spray treatment
u1ha.RasterizeBTKSpray(meta)

#%% Global Forest Change Loss Year (mosaic and reproject)
u1ha.RasterizeGFC_LossYear(meta)

#%% Filter Global Forest Change Loss Year (to remove known disturbances)
u1ha.FilterGFC_LossYear(meta)

#%% Extract mean climate data by BGC zone
u1ha.ClimateStatsByBGCZone(meta)

#%%

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\BC_MAJOR_WATERSHEDS\MAJOR_WATERSHED_SYSTEM.tif')
plt.matshow(z['Data'])

#%% Import climate from old BC1ha project in matlab

# Annual summary
z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\Terrain\BC_1ha_twi.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,meta['Paths']['bc1ha'] + '\Terrain\bc1ha_twi.tif')
del z
garc.collect()

# Annual summary
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1_c.tif')
z=gis.ClipRaster(z,[zS['xmin'],zS['xmax']],[zS['ymin'],zS['ymax']])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
del z; garc.collect()

#%% Convert old BC1ha .mat files to new BC1ha geotiffs with standardized extent

def ConvertClimateNormals():

    z2=zRef.copy()
    z2['Data']=np.zeros(zRef['Data'].shape,dtype='float')
    for mo in range(12):
        print(mo+1)
        fin=r'E:\Data\Climate\Canada\BC\Grids\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat'
        z=spio.loadmat(fin,squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat],axis=0) # #.astype(float)*z['z'][()][iSF]
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        z2['Data']=z2['Data']+z1['Data'].astype('float')
    z2['Data']=z2['Data']/12
    z2['Data']=z2['Data'].astype('int16')
    gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')

    z2=zRef.copy()
    z2['Data']=np.zeros(zRef['Data'].shape,dtype='int16')
    for mo in range(12):
        print(mo+1)
        fin=r'E:\Data\Climate\Canada\BC\Grids\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat'
        z=spio.loadmat(fin,squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        z2['Data']=z2['Data']+np.maximum(0,z1['Data'])
    gis.SaveGeoTiff(z2,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=z0
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z0=z0*100
        z1=zRef.copy()
        z1['Data']=z0.astype('int16')
        z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_etp_tmw_norm_1971to2000_comp_hist_v1\BC1ha_etp_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.tif')

    return

#%% Convert old NACID .mat files to new BC1ha geotiffs with standardized extent

def ConvertClimateNormals():
    zRef=gis.OpenGeoTiff(r'E:\Data\Climate\NACID\Geotiff\NACID\grid.tif')
    for mo in range(12):
        print(mo+1)
        z=spio.loadmat(r'E:\Data\Climate\NACID\Grids\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.mat',squeeze_me=True)
        idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
        iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
        z0=np.flip(z['z'][()][idat].astype(float)*z['z'][()][iSF],axis=0)
        z1=zRef.copy()
        z1['Data']=10*z0
        z1['Data']=z1['Data'].astype('int16')
        #z1=gis.ClipToRaster(z1,zRef)
        gis.SaveGeoTiff(z1,r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')

    return

z=gis.OpenGeoTiff(r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_rswd_mon_norm_1971to2000_si_hist_v1\NACID_rswd_mon_norm_1971to2000_si_hist_v1_1.tif')

