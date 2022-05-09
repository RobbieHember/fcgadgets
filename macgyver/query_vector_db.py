'''

QUERY RESULTS ACTIVITY TREATMENT LAYER

'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
from shapely.geometry import Polygon,Point
import time
import gc as garc
from fcgadgets.macgyver import utilities_general as gu
#from fcgadgets.macgyver import utilities_gis as gis

#%% Set path and layer

# List layers
flg=0
if flg==1:
    path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
    #path=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb'
    #path=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
    #path=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'
    fiona.listlayers(path)

with fiona.open(path,layer='RSLT_OPENING_SVW') as source:
    for feat in source:
        break

#%% Query wildfire perimiter
    
def GetWildfirePerimiter(meta,Year0,Year1):
    
    List=[]
    with fiona.open(meta['Paths']['Forest Inventory Disturbances'],layer='PROT_HISTORICAL_FIRE_POLYS_SP') as source:
        for feat in source:
            if feat['geometry']==None:
                continue  
            if (feat['properties']['FIRE_YEAR']>=Year0) & (feat['properties']['FIRE_YEAR']<=Year1):
                List.append(feat)
    gdf=gpd.GeoDataFrame.from_features(List,crs=meta['crs'])
    
    return gdf

#%% Query pest DB
    
def GetPestSev():
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
    
    d={}    
    with fiona.open(fin,layer='PEST_INFESTATION_POLY') as source:
        for feat in source:    
            a=feat['properties']['PEST_SEVERITY_CODE']            
            if a not in d:
                d[a]=1  
    return d
 
#d=GetPestSev()

#%%
    
def GetAnnualPestArea(Year0,Year1,sp_cd):
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'  
    
    SevName=np.array(['Trace','Light','Moderate','Severe','Very Severe','Grey Attack'])
    SevCD=np.array(['T','L','M','S','V','G'])
    
    tv=np.arange(Year0,Year1+1,1)
    
    d={}
    for k in sp_cd:
        d[k]={}    
        for k2 in SevCD:
            d[k][k2]=np.zeros(tv.size)
    
    with fiona.open(fin,layer='PEST_INFESTATION_POLY') as source:
        for feat in source:
            
            if feat['geometry']==None:
                continue  

            if (feat['properties']['CAPTURE_YEAR']<Year0) | (feat['properties']['CAPTURE_YEAR']>Year1):
                continue
            
            if feat['properties']['PEST_SPECIES_CODE'] not in sp_cd:
                continue
            
            nam=feat['properties']['PEST_SPECIES_CODE']
            
            iSev=np.where(SevCD==feat['properties']['PEST_SEVERITY_CODE'])[0]
            if iSev.size==0:
                print(feat['properties']['PEST_SEVERITY_CODE'])
                
            sev=SevCD[iSev][0]

            iT=np.where(tv==feat['properties']['CAPTURE_YEAR'])[0]
            
            d[nam][sev][iT]=d[nam][sev][iT]+feat['properties']['AREA_HA']
    
    return tv,d

#%% Query all openings
    
def GetOpeningInfo(ids):

    # Load dataset with CRS
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    # HARVESTING
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
    List=[]
    with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
        for feat in source:
            if feat['geometry']==None:
                continue  
            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
                List.append(feat)
                if len(ids)==1:
                    break
    gdf_cut=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    # RESULTS FOREST COVER INVENTORY
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[]
    with fiona.open(fin,layer='RSLT_FOREST_COVER_INV_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
                List.append(feat)
    gdf_fcinv=gpd.GeoDataFrame.from_features(List,crs=crs)
    
#    # RESULTS FOREST COVER INVENTORY
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
#    List=[]
#    with fiona.open(fin,layer='RSLT_FOREST_COVER_SILV_SVW') as source:
#        for feat in source:
#            if feat['geometry']==None:
#                continue
#            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
#                List.append(feat)
#    gdf_fcsilv=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    # RESULTS OPENING LAYER
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[]
    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue            
            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
                List.append(feat)   
    gdf_op=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    return gdf_cut,gdf_fcinv,gdf_op

#%% Query all openings
    
def GetOpeningsAll():
    
    # Load dataset with CRS
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue            
            List[cnt]=feat
            cnt=cnt+1
    List=List[0:cnt-1]
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    return gdf

#%% Query planting for specific years
    
def GetOpeningsWithPlanting(yr0,yr1):
    
    # Load dataset with CRS
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
    crs=gdf_bm.crs
    del gdf_bm
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    
    # Import array of openings with planting
    op_id_pl=np.zeros(50000)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            if feat['properties']['SILV_BASE_CODE']!='PL':
                continue
            if feat['properties']['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
            if (Year<yr0) | (Year>yr1):
                continue                
            op_id_pl[cnt]=feat['properties']['OPENING_ID']
            cnt=cnt+1
    
    op_id_pl=op_id_pl[0:cnt-1]
    
    # Get full list of opening ID
    op_id=np.zeros(1000000)
    cnt=0
    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
        for feat in source:             
            op_id[cnt]=feat['properties']['OPENING_ID']
            cnt=cnt+1
    op_id=op_id[0:cnt-1]
    
    c,ia,ib=np.intersect1d(op_id_pl,op_id,return_indices=True)
    
    List=[None]*int(1e5)
    cnt=0
    cnt0=0
    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
        for feat in source:
            if cnt0 not in ib:
                cnt0=cnt0+1
                continue
            if feat['geometry']==None:
                continue            
            List[cnt]=feat
            cnt=cnt+1
            cnt0=cnt0+1
    List=List[0:cnt-1]    
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=crs)
    
    return gdf


def GetSurveySpatial(Year1,Year2):
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if feat['properties']['SILV_BASE_CODE']!='SU':
                continue
            if feat['properties']['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
            if (Year<Year1) | (Year>Year2):
                continue
            if feat['properties']['GEOMETRY_Area']/10000>2000:
                continue
                
            List[cnt]=feat
            cnt=cnt+1   
    List=List[0:cnt-1]
    
    gdf=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)

    #gdf.plot()
    #gdf_fcinv2=gdf_fcinv.to_crs({'init':'epsg:4326'})
    # Save
    #gdf_fcinv.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')
    return gdf