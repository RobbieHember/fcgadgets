
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
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables

meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Rasterize Forest Cover Attributes
# Use FC ID created in ArcGIS to connect with attributes
u1ha.RasterizeForestCover(meta)

#%%

z=u1ha.Import_Raster(meta,[],['refg','lcc1_c','PROJ_AGE_1','harv_yr_con1','fire_yr'])
yr=2023-z['PROJ_AGE_1']['Data']
plt.close('all')
plt.matshow(yr,clim=[1800,2000])

#%%

ikp=np.where( (z['fire_yr']['Data']>0) )
plt.close('all')
x=yr[ikp][0::200]
y=z['fire_yr']['Data'][ikp][0::200]
bw=10; bin=np.arange(1900,2030,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'bo')

#%%

ikp=np.where( (z['harv_yr_con1']['Data']>0) )
plt.close('all')
x=yr[ikp][0::200]
y=z['harv_yr_con1']['Data'][ikp][0::200]
bw=10; bin=np.arange(1900,2030,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'bo')

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

u1ha.RasterizeFromSource(meta,zRef,'BC_MAJOR_WATERSHEDS','MAJOR_WATERSHED_CODE')
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

#%% Rasterize OPENING ID FROM OPENING LAYER (1 hour)
u1ha.RasterizeOpeningID_2(meta)

#%% Generate sparse inputs
# This should speed import by an order of magnitude
rgsf=100
mask='All'
u1ha.GenerateSparseInputs(meta,rgsf,mask)

#%% Gap-fill BGC Zone
u1ha.GapFill_BGCZ(meta)

#%% Rasterize wildfire
u1ha.RasterizeWildfire(meta)

#%% Current-year wildfire
# Current year + preveous year (sometimes missing)
# Need to add this to historical and then re-create the condensed layers
u1ha.RasterizeWildfireCurrentYear(meta)

#%% Ecozones of Canada
u1ha.RasterizeEcozonesOfCanada(meta)

#%% Rasterize insects from AOS
u1ha.RasterizeInsects(meta,zRef)

#%% Rasterize planting from RESULTS
u1ha.RasterizePlanting(meta)

# Plot time series to confirm it worked
flg=0
if flg==1:
    lNam='RSLT_ACTIVITY_TREATMENT_SVW'
    vNam='PL_All'
    tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)
    plt.plot(tv,N,'-bo')

#%% Derive planting Type

import fcgadgets.cbrunner.cbrun_util as cbu
meta=cbu.Load_LUTs_Modelling(meta)

tv=np.arange(1960,2023,1)
ptNam=np.array(list(meta['LUT']['Derived']['RegenTypeNO'].keys()))
ptID=np.array(list(meta['LUT']['Derived']['RegenTypeNO'].values()))
zD=u1ha.Import_Raster(meta,[],['harv_yr_con1','kd_yr'])

zP={}
for i in range(6):
    zP[i]={}
    zP[i]['Type']=np.zeros(zRef['Data'].shape,dtype='int8')
    zP[i]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')['Data']
    zP[i]['STC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_TECHNIQUE_CODE.tif')['Data']
    zP[i]['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']

# vL=['PL','RP','FP','CG','RO','RR','SE','SL']
# stcN={}
# for v in vL:
#     stcN[v]=0
#     for iP in range(6):
#         ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
#         stcN[v]=stcN[v]+ind[0].size

cd=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']
stcInclude=[cd['PL'],cd['RP'],cd['FP'],cd['RR']]

A_Tot=np.zeros(tv.size)
A=np.zeros((tv.size,ptNam.size))
DL=np.zeros(zRef['Data'].shape,dtype='int8')
for iT in range(tv.size):
    print(tv[iT])
    zWF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tv[iT]-1) + '.tif')['Data']
    zIBM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(tv[iT]-1) + '.tif')['Data']
    ind=np.where(zIBM>=3); DL[ind]=meta['LUT']['Event']['IBM'] # Severity greater than light
    ind=np.where(zD['harv_yr_con1']['Data']==tv[iT]-1); DL[ind]=meta['LUT']['Event']['Harvest']; #H_Yr[ind]=zD['harv_yr_con1']['Data'][ind]
    ind=np.where(zD['kd_yr']['Data']==tv[iT]-1); DL[ind]=meta['LUT']['Event']['Knockdown']
    ind=np.where(zWF>0); DL[ind]=meta['LUT']['Event']['Wildfire']
    
    pl_oc=np.zeros(zRef['Data'].shape,dtype='int8')
    pl_stc=np.zeros(zRef['Data'].shape,dtype='int16')
    pl_fsc=np.zeros(zRef['Data'].shape,dtype='int16')
    for iP in range(6):
        ind=np.where( (zP[iP]['Year']==tv[iT]) & (np.isin(zP[iP]['STC'],stcInclude)==True) )
        #ind=np.where( (zP[iP]['Year']==tv[iT]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) & (np.isin(zP[iP]['STC'],stcInclude)==True) )
        pl_oc[ind]=1
        pl_stc[ind]=zP[iP]['STC'][ind]
        pl_fsc[ind]=zP[iP]['FSC'][ind]
    
    ind=np.where( (pl_oc==1) )
    A_Tot[iT]=ind[0].size
    
    z0=np.zeros(zRef['Data'].shape,dtype='int8')
    
    # Replanting
    ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Replanting']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Replanting']-1]=ind[0].size
    # Fill Planting
    ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Fill Planting']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Fill Planting']-1]=ind[0].size
    # Road rehab
    ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Road Rehabilitation']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Road Rehabilitation']-1]=ind[0].size
    # Back to back planting
    ind=np.where( (pl_oc==1) & (z0==0) & (DL==meta['LUT']['Event']['Planting']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Back-to-back Planting']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Back-to-back Planting']-1]=ind[0].size    
    # Salvage
    ind=np.where( (pl_oc==1) & (tv[iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (z0==0) & (DL==meta['LUT']['Event']['Harvest']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Salvage and Planting']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Salvage and Planting']-1]=ind[0].size
    # Knockdown
    ind=np.where( (pl_oc==1) & (z0==0) & (DL==meta['LUT']['Event']['Knockdown']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Knockdown and Planting']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Knockdown and Planting']-1]=ind[0].size
    # Straight fire
    ind=np.where( (pl_oc==1) & (z0==0) & (DL==meta['LUT']['Event']['Wildfire']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Straight-to-planting Post Wildfire']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Straight-to-planting Post Wildfire']-1]=ind[0].size
    # Straight insect
    ind=np.where( (pl_oc==1) & (z0==0) & (DL==meta['LUT']['Event']['IBM']) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Straight-to-planting Post Insect Outbreak']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Straight-to-planting Post Insect Outbreak']-1]=ind[0].size
    # NSR backlog
    ind=np.where( (pl_oc==1) & (DL==meta['LUT']['Event']['Harvest']) & (zD['harv_yr_con1']['Data']<=1987) & (z0==0) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['NSR Backlog']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['NSR Backlog']-1]=ind[0].size        
    # Unknown
    ind=np.where( (pl_oc==1) & (z0==0) )
    z0[ind]=meta['LUT']['Derived']['RegenTypeNO']['Unknown']
    A[iT,meta['LUT']['Derived']['RegenTypeNO']['Unknown']-1]=ind[0].size
    
    DL[(pl_oc==1)]=meta['LUT']['Event']['Planting']
    
    # Pack
    for i in range(6):
        ind=np.where( (zP[i]['Year']==tv[iT]) & (zP[i]['Type']!=0) )
        zP[i]['Type'][ind]=z0[ind]

gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\AreaSummary_All',A)

# Save packed regen type
for i in range(6):
    z1=zRef.copy()
    z1['Data']=zP[i]['Type']
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')

#%%

def Plot_AIL_TS_NOSE(meta,tv):
    
    iT=np.where( (tv>=1960) & (tv<=2050) )[0]
    cl=np.array([[0.6,0.75,1],[0,0,0.5],[0.6,1,0],[0.15,0.75,0],[1,1,0],[0.9,0.87,0.84],[1,0.75,0.25],[0.75,0,0],[0.7,0.65,0.9],[0.25,0.25,0.25]])
    
    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,9));
    A_cumu=np.zeros(tv.size)
    cnt=0
    for k in meta['LUT']['Derived']['RegenTypeNO'].keys():        
        plt.bar(tv,A[:,cnt]/1e3,0.8,bottom=A_cumu,facecolor=cl[cnt,:],label=k)
        A_cumu=A_cumu+A[:,cnt]/1e3; 
        cnt=cnt+1
    #plt.plot(tv,A_Tot/1e3,'ks',ms=2.5,mec='k',mfc='w',mew=0.5)
    ax.set(xticks=np.arange(1950,2225+1,10),ylabel='Implementation level (Kha yr$^{-1}$)',
           xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[tv[iT][0]-0.75,tv[iT][-1]+0+.75],ylim=[0,140]) #
    plt.legend(frameon=False,loc='upper right',facecolor=[1,1,1],labelspacing=0.25,ncol=2)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
    plt.tight_layout()
        
    flg=0
    if flg==0:
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
    else:
        at=gu.ipickle(r'C:\Users\rhember\Documents\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\at.pkl')
    
    #ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') & (np.isin(at['gdf']['SILV_FUND_SOURCE_CODE'],meta['Param']['BE']['FSC']['NO List Name'])==True) )[0]
    ikp=np.where(  (at['gdf']['RESULTS_IND']=='Y') & (at['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
    for k in at['gdf'].keys():
        at['gdf'][k]=at['gdf'][k][ikp]
    at['gdf']['Year']=np.zeros(at['gdf']['ACTIVITY_TREATMENT_UNIT_ID'].size)
    for i in range(at['gdf']['Year'].size):
        at['gdf']['Year'][i]=int(at['gdf']['ATU_COMPLETION_DATE'][i][0:4])    
    A_Tot2=np.zeros(tv.size)
    A_Tot2b=np.zeros(tv.size)
    for iT in range(tv.size):
        ind=np.where( (at['gdf']['Year']==tv[iT]) )[0]
        A_Tot2[iT]=np.sum(at['gdf']['ACTUAL_TREATMENT_AREA'][ind])
        A_Tot2b[iT]=np.sum(at['gdf']['GEOMETRY_Area'][ind]/1e4)
    plt.plot(tv,A_Tot2/1e3,'k^',ms=3,mec='k',mfc='w',mew=0.5)    
    plt.plot(tv,A_Tot2b/1e3,'ks',ms=3,mec='k',mfc='w',mew=0.5)
    #gu.PrintFig(meta['Paths'][pNamC]['Figures'] + '\\PL_NO_AIL_ByProjectType','png',900)
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS_NOSEC\PL_NO_AIL_ByProjectType','png',900)
    return

Plot_AIL_TS_NOSE(meta,mos,pNamC,tv)

#%% Planting non-obligation Mask
# Dabbled in doing it this way and then went back to doing it on the fly
u1ha.Planting_NonOb_Mask(meta,zRef)

#%% Generate sparse inputs for NOSE project
rgsf=10
mask='NOSE'
u1ha.GenerateSparseInputs(meta,rgsf,mask)

#%% Planting layer (species and genetic worth)
u1ha.RasterizePlantingLayer(meta)

#%% Rasterize fertilization and knockdown
u1ha.RasterizeActivities(meta)

#%% Rasterize consolidated cutblocks
u1ha.RasterizeConsolidatedCutblocks(meta,zRef)

#%% Harvest early reconstruction
u1ha.DeriveEarlyHarvest(meta)

#%% Harvest consolidated 1 (w/o early reconstruction)
u1ha.ConsolidateHarvesting(meta)

#%% Percent dead from cruises
u1ha.RasterizeCruisePctDead(meta)

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
u1ha.ConsolidateTreeDensityClass(meta)

#%% Consolidate protected lands
u1ha.ConsolidateProtectedLand(meta)

#%% Rasterize early land use change year
u1ha.RasterizeEaryLandUseChangeYear(meta)

#%% Rasterize distance from roads
u1ha.DeriveDistanceFromRoads(meta)

#%% Rasterize distance from timber facilities
u1ha.DeriveDistanceFromFacility(meta)

#%% Consolidated retention layer
u1ha.ConsolidateRetention(meta)

#%% Harvest regeneration type
u1ha.DeriveHarvestRegenType(meta)

#%% Rasterize aerial spray treatment
u1ha.RasterizeBTKSpray(meta)

#%% Reproject National Burn Severity data
# Step 1) Resample in ArcGIS; 2) Clip to BC in ArcGIS using a BC1ha grid for extent; 
# 3) export geotiff; 4) reproject and clip to ensure BC1ha grid extent (below))
u1ha.Import_BurnSeverityCanada(meta)

#%% Derive burn severity to pair with wildfire perimeters
u1ha.DeriveBurnSeverityConsolidated(meta)

#%% Global Forest Change Loss Year (mosaic and reproject)
u1ha.RasterizeGFC_LossYear(meta)

#%% Filter Global Forest Change Loss Year (to remove known disturbances)
u1ha.FilterGFC_LossYear(meta)

#%% Reproject CEC Land Use Map 2020
u1ha.ReprojectAndCompressCECLandUse(meta)

#%% Reproject data from NTEMS (age, land cover, harvest year)
u1ha.ReprojectDataFromNTEMS(meta)

#%% Extract mean climate data by BGC zone
u1ha.ClimateStatsByBGCZone(meta)

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

#%% Rasterize forest cover
fci={}
fci['Path']=meta['Paths']['GDB']['Results']
fci['Layer']='RSLT_FOREST_COVER_INV_SVW'; # fiona.listlayers(op['Path'])
fci['crs']=meta['Geos']['crs']
fci['Keep Geom']='Off'
fci['Select Openings']=np.array([])
fci['SBC']=np.array([])
fci['STC']=np.array([])
fci['SMC']=np.array([])
fci['FSC']=np.array([])
fci['ROI']=[]
fci['gdf']=qgdb.Query_Openings(fci,[])

zID=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_FOREST_COVER_INV_SVW\\fcid.tif')
iMask=np.where(zID['Data']!=0)
idx=gu.IndicesFromUniqueArrayValues(zID['Data'][iMask])
#zID2=zID['Data'][iMask]
#uID=np.unique(zID2)

d={}
d['SSC']=np.zeros(
for i in range(uID.size):
    ind=np.where(zID['Data']
    zSSC['Data']meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_STATUS_CODE'].keys()

zSSC=zRef.copy()
zSSC['Data']=np.zeros(zRef['Data'].shape,dtype='int16')

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

