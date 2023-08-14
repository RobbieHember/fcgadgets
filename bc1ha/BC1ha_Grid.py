
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
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables

meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Rasterize Forest Cover Attributes
# Use FC ID created in ArcGIS to connect with attributes
u1ha.RasterizeForestCoverInventory(meta)

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
u1ha.RasterizeOpeningID(meta)

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

#%%

def DeriveRegenTypeCompilation(meta):
    meta=cbu.Load_LUTs_Modelling(meta)
    
    d={}
    d['tv']=np.arange(1960,2023,1)
    ptNam=np.array(list(meta['LUT']['Derived']['RegenType'].keys()))
    ptID=np.array(list(meta['LUT']['Derived']['RegenType'].values()))
    zD=u1ha.Import_Raster(meta,[],['harv_yr_con1','kd_yr','pdead_cruise'])
    
    zP={}
    for iP in range(6):
        zP[iP]={}
        zP[iP]['Type']=np.zeros(zRef['Data'].shape,dtype='int8')
        zP[iP]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
        zP[iP]['STC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_TECHNIQUE_CODE.tif')['Data']
        zP[iP]['FSC']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
    
    zDS={}
    for iDS in range(3):
        zDS[iDS]={}
        zDS[iDS]['Year']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\DS_All_' + str(iP+1) + '_Year.tif')['Data']
    
    # Summarize frequency of each Stocking Type Code (All and NOSE)
    # d['STC Summary']={}
    # d['STC Summary']['stc']=['PL','RP','FP','CG','RO','RR','SE','SL']
    # d['STC Summary']['N']=np.zeros((len(d['STC Summary']['stc']),2))
    # for v in d['STC Summary']['stc']:
    #     for iP in range(6):
    #         ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) )
    #         d['STC Summary']['N'][cnt,0]=d['STC Summary']['N'][cnt,0]+ind[0].size
    #         ind=np.where( (zP[iP]['STC']==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE'][v]) & (np.isin(zP[iP]['FSC'],meta['Param']['BE']['FSC']['NO List ID'])==True) )
    #         d['STC Summary']['N'][cnt,1]=d['STC Summary']['N'][cnt,1]+ind[0].size
    
    cd=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']
    stcInclude=[cd['PL'],cd['RP'],cd['FP'],cd['RR']]
    
    d['A Tot']=np.zeros(d['tv'].size)
    d['A']=np.zeros((d['tv'].size,ptNam.size))
    d['A NOSE']=np.zeros((d['tv'].size,ptNam.size))
    DL=np.zeros(zRef['Data'].shape,dtype='int8')
    for iT in range(d['tv'].size):
        print(d['tv'][iT])
        zWF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(d['tv'][iT]-1) + '.tif')['Data']
        zIBM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(d['tv'][iT]-1) + '.tif')['Data']
        ind=np.where(zIBM>=3); DL[ind]=meta['LUT']['Event']['IBM'] # Severity greater than light
        ind=np.where(zD['harv_yr_con1']['Data']==d['tv'][iT]-1); DL[ind]=meta['LUT']['Event']['Harvest']; #H_Yr[ind]=zD['harv_yr_con1']['Data'][ind]
        ind=np.where(zD['kd_yr']['Data']==d['tv'][iT]-1); DL[ind]=meta['LUT']['Event']['Knockdown']
        ind=np.where(zWF>0); DL[ind]=meta['LUT']['Event']['Wildfire']
        
        pl_oc=np.zeros(zRef['Data'].shape,dtype='int8')
        pl_stc=np.zeros(zRef['Data'].shape,dtype='int16')
        pl_fsc=np.zeros(zRef['Data'].shape,dtype='int16')
        for iP in range(6):
            ind=np.where( (zP[iP]['Year']==d['tv'][iT]) & (np.isin(zP[iP]['STC'],stcInclude)==True) )            
            pl_oc[ind]=1
            pl_stc[ind]=zP[iP]['STC'][ind]
            pl_fsc[ind]=zP[iP]['FSC'][ind]
        
        ind=np.where( (pl_oc==1) )
        d['A Tot'][iT]=ind[0].size
        
        zType0=np.zeros(zRef['Data'].shape,dtype='int8')
        
        # Replanting        
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Replanting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RP']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Replanting']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Replanting']
        
        # Fill Planting        
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Fill Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['FP']) )        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Fill Planting']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Fill Planting']
        
        # Road rehab
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Road Rehabilitation']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (pl_stc==meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_TECHNIQUE_CODE']['RR']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Road Rehabilitation']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Road Rehabilitation']
        
        # Back to back planting
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Back-to-back Planting']-1]=ind[0].size    
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Planting']) )        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Back-to-back Planting']-1]=ind[0].size    
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Back-to-back Planting']
        
        # Salvage
        ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Salvage and Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (d['tv'][iT]>=2004) & (zD['harv_yr_con1']['Data']>1987) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) | \
                     (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['pdead_cruise']['Data']>=50) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Salvage and Planting']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Salvage and Planting']
        
        # Knockdown
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Knockdown and Planting']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Knockdown']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Knockdown and Planting']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Knockdown and Planting']
        
        # Straight fire
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Wildfire']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']        
        
        # Straight insect
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['IBM']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']-1]=ind[0].size
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Beetles']
        
        # Straight other
        
        # Harvest and Planting NSR backlog
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (zD['harv_yr_con1']['Data']<=1987) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']-1]=ind[0].size        
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Harvest and Planting NSR Backlog']
        
        # Harvest and Planting
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting']-1]=ind[0].size 
        ind=np.where( (pl_oc==1) & (zType0==0) & (DL==meta['LUT']['Event']['Harvest']) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Harvest and Planting']-1]=ind[0].size         
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Harvest and Planting']
        
        # Unknown
        ind=np.where( (pl_oc==1) & (zType0==0) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
        d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Unknown']-1]=ind[0].size
        ind=np.where( (pl_oc==1) & (zType0==0) )
        d['A'][iT,meta['LUT']['Derived']['RegenType']['Unknown']-1]=ind[0].size        
        zType0[ind]=meta['LUT']['Derived']['RegenType']['Unknown']
        
        # Update last disturbance if planting occurs
        DL[(pl_oc==1)]=meta['LUT']['Event']['Planting']
        
        # Direct seeding
        for iDS in range(3):
            ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) & (np.isin(pl_fsc,meta['Param']['BE']['FSC']['NO List ID'])==True) )
            d['A NOSE'][iT,meta['LUT']['Derived']['RegenType']['Direct Seeding']-1]=ind[0].size
            ind=np.where( (zDS[iDS]['Year']==d['tv'][iT]) )            
            d['A'][iT,meta['LUT']['Derived']['RegenType']['Direct Seeding']-1]=ind[0].size            
            zType0[ind]=meta['LUT']['Derived']['RegenType']['Direct Seeding']
        
        # Pack
        for iP in range(6):
            ind=np.where( (zP[iP]['Year']==d['tv'][iT]) )
            zP[iP]['Type'][ind]=zType0[ind]
    
    # Save summary
    gu.opickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\RegenTypeSummary_All',d)
    
    # Save packed regen type
    for i in range(6):
        z1=zRef.copy()
        z1['Data']=zP[i]['Type']
        gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_RegenType.tif')
    
    return

#%% Derive Regen Type (used for non-ob stand establishment)
u1ha.DeriveRegenTypeCompilation(meta)

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

