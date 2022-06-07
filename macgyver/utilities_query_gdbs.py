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
path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb'
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'

# List layers
fiona.listlayers(path)

#%% Set layer

lyr='RSLT_ACTIVITY_TREATMENT_SVW'
#lyr='RSLT_OPENING_SVW'
#lyr='RSLT_FOREST_COVER_INV_SVW'
#lyr='VEG_COMP_LYR_R1_POLY'

#%% Have a quick look at layer to see what the fields are

with fiona.open(path,layer=lyr) as source:
    for feat in source:
        break
feat['properties']


#%%
path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
with fiona.open(path,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
    for feat in source:
        if feat['properties']['FIA_PROJECT_ID']=='FCI0000427':
            print('found')
            break

#%% Load dataset with CRS

gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
crs=gdf_bm.crs
del gdf_bm

#%% Annual harvest area from RESULTS
    
def GetAnnualHarvestAreaFromRESULTS():
           
    tv=np.arange(1950,2022,1)
    
    d={}
    d['Year']=tv
    d['Area Harvested']=np.zeros(tv.size)
    
    with fiona.open(path,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            
            #if feat['geometry']==None:
            #    continue  
            if feat['properties']['ATU_COMPLETION_DATE']==None:
                continue
            if feat['properties']['ACTUAL_TREATMENT_AREA']==None:
                continue
            
            if (feat['properties']['DISTURBANCE_CODE']=='L') | (feat['properties']['DISTURBANCE_CODE']=='S') | (feat['properties']['DISTURBANCE_CODE']=='R'):
                #break
            
                Year=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])
            
                iT=np.where(tv==Year)[0]
                if iT.size==0:
                    continue
                d['Area Harvested'][iT]=d['Area Harvested'][iT]+feat['properties']['ACTUAL_TREATMENT_AREA']
    
    # Save
    gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromRESULTS.pkl',d)
    
    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(9,6)); 
    plt.plot(tv,d['Area Harvested']/1e3,'-bo')
    ax.set(position=[0.085,0.125,0.88,0.84],xlim=[1949.5,2021.5],xticks=np.arange(1800,2120,5), \
           yticks=np.arange(0,275,25),ylabel='Area harvested (Thousand ha yr$^-$$^1$)',xlabel='Time, years')
    #ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area\AnnualAreaHarvested_FromRESULTS','png',900)
    
    return d

#%% Annual harvest area from consolidated cutblock DB
    
def GetAnnualHarvestAreaFromConCutblocksDB():
    
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'  
        
    tv=np.arange(1950,2022,1)
    
    d={}
    d['Year']=tv
    d['Area Harvested']=np.zeros(tv.size)
    
    with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
        for feat in source:
            
            if feat['geometry']==None:
                continue  
            #break
            iT=np.where(tv==feat['properties']['HARVEST_YEAR'])[0]
            if iT.size==0:
                continue
            d['Area Harvested'][iT]=d['Area Harvested'][iT]+feat['properties']['AREA_HA']
    
    # Save
    gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromConCutblocksDB.pkl',d)
    
    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(9,6)); 
    plt.plot(tv,d['Area Harvested']/1e3,'-bo')
    ax.set(position=[0.085,0.125,0.88,0.84],xlim=[1949.5,2021.5],xticks=np.arange(1800,2120,5), \
           yticks=np.arange(0,275,25),ylabel='Area harvested (Thousand ha yr$^-$$^1$)',xlabel='Time, years')
    #ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area\AnnualAreaHarvested_FromConsolidatedCutblocks','png',900)
    
    return d

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

#%% Annual pest-affected area
    
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



#%% Overlay analysis
    
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
