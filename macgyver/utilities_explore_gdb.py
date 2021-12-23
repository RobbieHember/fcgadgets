'''

EXPLORE GDB FILE

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
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Set path and layer

# Path to forest cover archive geodatabase
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930\Results.gdb'
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb'
path=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930\Disturbances.gdb'

# List layers
fiona.listlayers(path)

# Set layer
lyr='RSLT_ACTIVITY_TREATMENT_SVW'
#lyr='RSLT_FOREST_COVER_INV_SVW'
#lyr='VEG_COMP_LYR_R1_POLY'

# Load dataset with CRS
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

#%% Loop through gdb

with fiona.open(path,layer=lyr) as source:
    for feat in source:
        break
        if feat['properties']['SILV_BASE_CODE']=='SU':
            break


with fiona.open(path,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
    for feat in source:
        break

A=np.array([0.0])
with fiona.open(path,layer='VEG_BURN_SEVERITY_SP') as source:
    for feat in source:
        if (feat['properties']['FIRE_YEAR']<2017) | (feat['properties']['FIRE_YEAR']>2018) | (feat['properties']['BURN_SEVERITY_RATING']=='Unburned'):
            continue
        A=A+feat['properties']['AREA_HA']
        
        
        break

N=0
N_all=0
with fiona.open(path,layer='RSLT_FOREST_COVER_INV_SVW') as source:
    for feat in source:
        break
        
    N_all=N_all+1
        if feat['geometry']==None:
            N=N+1
        
        #prp=feat['properties']
        #if prp['OPENING_ID']==15241:
        #    if feat['geometry']!=None:
        #        print('found')


#%% Get spatial of surveys

def GetSurveySpatial():
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
            if Year<2018:
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



#%%
    
def GetPlantingSpatial(op_id_pl):
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if np.isin(feat['properties']['OPENING_ID'],op_id_pl)==False:
                continue               
            List[cnt]=feat
            cnt=cnt+1
            print(cnt)
    List=List[0:cnt-1]    
    gdf=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)
    return gdf


#%% 
    
def GetWildfire():    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='PROT_HISTORICAL_FIRE_POLYS_SP') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if (feat['properties']['FIRE_YEAR']!=2017) & (feat['properties']['FIRE_YEAR']!=2018):
                continue                
            List[cnt]=feat
            cnt=cnt+1
    List=List[0:cnt-1]

    # Give it the BC spatial reference system
        
    gdf=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)

    #gdf.plot()
    #gdf_fcinv2=gdf_fcinv.to_crs({'init':'epsg:4326'})
    # Save
    #gdf_fcinv.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')
    return gdf



#%% 
    
gdf_su=GetSurveySpatial()

gdf_pl=GetPlanting()

gdf_wf=GetWildfire()

gdf_pl_wf=gpd.overlay(gdf_pl,gdf_wf,how='intersection')

gdf_su_wf=gpd.overlay(gdf_su,gdf_wf,how='intersection')

gdf_su_wf_notpl=gpd.overlay(gdf_su_wf,gdf_pl,how='difference')

gdf_su_wf_notpl.plot()
gdf_su_wf.plot()
gdf_pl.plot()

gdf_pl_wf.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\pl_wf.geojson',driver='GeoJSON')
gdf_su_wf_notpl.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\su_wf_notpl.geojson',driver='GeoJSON')
gdf_su_wf.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\su_wf.geojson',driver='GeoJSON')
gdf_pl.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\pl.geojson',driver='GeoJSON')
