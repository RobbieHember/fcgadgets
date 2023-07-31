'''

QUERY GEODATABASES

'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
from shapely.geometry import Polygon,Point,box,shape
import time
import gc as garc
from fcgadgets.macgyver import util_general as gu
from fcgadgets.macgyver import util_gis as gis
from fcgadgets.macgyver import util_inventory as invu
from fcgadgets.cbrunner import cbrun_util as cbu
from fcgadgets.bc1ha import bc1ha_util as bc1hau

#%% Query openings

def Query_Openings(meta,roi):

    if type(roi)==dict:
        # Using ROI
        # Convert to shapely object
        roi_poly=roi['gdf']['bound'].iloc[0]['geometry']
    else:
        # Not using ROI
        pass

    if meta['Keep Geom']=='Off':

        #======================================================================
        # Only keep attribute data
        #======================================================================

        n=3000000
        d={}
        cnt=0
        with fiona.open(meta['Path'],layer=meta['Layer']) as source:
            for feat in source:
                #geom=dict(feat['geometry'].items())
                prop=dict(feat['properties'].items())

                #--------------------------------------------------------------
                # Filters
                #--------------------------------------------------------------

                #if prop==None:
                #    continue

                if meta['Select Openings'].size>0:
                    if np.isin(prop['OPENING_ID'],meta['Select Openings'])==False:
                        continue

                if type(meta['ROI'])==dict:
                    # Using ROI
                    s=shape(prop)
                    if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                        continue

                if meta['SBC'].size!=0:
                    if np.isin(prop['SILV_BASE_CODE'],meta['SBC'])==False:
                        continue
                    # Planting (isolate real planting)
                    #if meta['SBC']=='PL':
                    #    if (prop['SILV_TECHNIQUE_CODE']=='SE') | (prop['SILV_TECHNIQUE_CODE']=='CG') | (prop['SILV_METHOD_CODE']=='LAYOT'):
                    #        continue

                if meta['STC'].size!=0:
                    if np.isin(prop['SILV_TECHNIQUE_CODE'],meta['STC'])==False:
                        continue

                if meta['SMC'].size!=0:
                    if np.isin(prop['SILV_METHOD_CODE'],meta['SMC'])==False:
                        continue

                if meta['FSC'].size!=0:
                    if np.isin(prop['SILV_FUND_SOURCE_CODE'],meta['FSC'])==False:
                        continue

                if meta['Layer']=='RSLT_ACTIVITY_TREATMENT_SVW':

                    if (prop['RESULTS_IND']!='Y') | (prop['ATU_COMPLETION_DATE']==None):
                        continue

                    prop['Year']=int(prop['ATU_COMPLETION_DATE'][0:4])

                    # Planting (isolate real planting)
                    if meta['SBC']=='PL':
                        if (prop['SILV_TECHNIQUE_CODE']=='SE') | (prop['SILV_TECHNIQUE_CODE']=='CG') | (prop['SILV_METHOD_CODE']=='LAYOT'):
                            continue

                if cnt==0:

                    # Initialize
                    for k in prop.keys():
                        try:
                            d[k]=np.zeros(n)
                            d[k][cnt]=prop[k]
                        except:
                            d[k]=np.array(['' for _ in range(n)],dtype=object)
                            d[k][cnt]=prop[k]
                else:

                    for k in prop.keys():
                        try:
                            d[k][cnt]=prop[k]
                        except:
                            d[k]=np.array(['' for _ in range(n)],dtype=object)
                            d[k][cnt]=prop[k]
                cnt=cnt+1

        # Truncate
        for k in d.keys():
            d[k]=d[k][0:cnt]

    else:

        #======================================================================
        # Keep geometries
        #======================================================================

        List=[None]*int(1e6)
        cnt=0
        with fiona.open(meta['Path'],layer=meta['Layer']) as source:
            for feat in source:

                prop=dict(feat['properties'].items())

                if prop==None:
                    continue

                try:
                    geom=dict(feat['geometry'].items())
                except:
                    continue

                if meta['Layer']=='RSLT_ACTIVITY_TREATMENT_SVW':
                    if prop['SILV_BASE_CODE']=='SU':
                        continue

                if type(roi)==dict:
                    # Using ROI
                    s=shape(geom)
                    if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                        continue

                if meta['Select Openings'].size>0:
                    if np.isin(prop['OPENING_ID'],meta['Select Openings'])==False:
                        continue

                if meta['SBC'].size!=0:
                    if np.isin(prop['SILV_BASE_CODE'],meta['SBC'])==False:
                        continue

                    # Planting (isolate real planting)
                    if meta['SBC']=='PL':
                        if (prop['SILV_TECHNIQUE_CODE']=='SE') | (prop['SILV_TECHNIQUE_CODE']=='CG') | (prop['SILV_METHOD_CODE']=='LAYOT'):
                            continue

                if meta['STC'].size!=0:
                    if np.isin(prop['SILV_TECHNIQUE_CODE'],meta['STC'])==False:
                        continue

                if meta['SMC'].size!=0:
                    if np.isin(prop['SILV_METHOD_CODE'],meta['SMC'])==False:
                        continue

                if meta['FSC'].size!=0:
                    if np.isin(prop['SILV_FUND_SOURCE_CODE'],meta['FSC'])==False:
                        continue

                if meta['SOC1'].size!=0:
                    if np.isin(prop['SILV_OBJECTIVE_CODE_1'],meta['SOC1'])==False:
                        continue

                if 'Denudation' in meta.keys():
                    if (np.isin(prop['DENUDATION_1_DISTURBANCE_CODE'],meta['Denudation'])==False) & (np.isin(prop['DENUDATION_2_DISTURBANCE_CODE'],meta['Denudation'])==False):
                        continue
                    if np.isin(prop['DENUDATION_1_DISTURBANCE_CODE'],meta['Denudation'])==True:
                        if prop['DENUDATION_1_COMPLETION_DATE']!=None:
                            prop['Year']=int(prop['DENUDATION_1_COMPLETION_DATE'][0:4])
                            #prop['Month']=int(prop['DENUDATION_1_COMPLETION_DATE'][5:7])
                    if np.isin(prop['DENUDATION_2_DISTURBANCE_CODE'],meta['Denudation'])==True:
                        if prop['DENUDATION_2_COMPLETION_DATE']!=None:
                            prop['Year']=int(prop['DENUDATION_2_COMPLETION_DATE'][0:4])
                            #prop['Month']=int(prop['DENUDATION_2_COMPLETION_DATE'][5:7])

                if 'STOCKING_TYPE_CODE' in meta.keys():
                    if np.isin(prop['STOCKING_TYPE_CODE'],meta['STOCKING_TYPE_CODE'])==False:
                        continue

                if meta['Layer']=='RSLT_ACTIVITY_TREATMENT_SVW':
                    if (prop['RESULTS_IND']!='Y') | (prop['ATU_COMPLETION_DATE']==None):
                        continue
                    prop['Year']=int(prop['ATU_COMPLETION_DATE'][0:4])

                if 'Drop Props' in meta.keys():
                    a={'OPENING_ID':prop['OPENING_ID'],'GEOMETRY_Area':prop['GEOMETRY_Area']}
                    prop=a

                List[cnt]=feat
                cnt=cnt+1

        List=List[0:cnt-1]
        d=gpd.GeoDataFrame.from_features(List,crs=meta['crs'])

        if type(roi)==dict:
            d=bc1hau.ClipGDF_ByROI(d,roi)

    return d

#%% Query wildfire perimiter

def Query_Wildfire(meta,roi):

    List=[]
    with fiona.open(meta['Path'],'r',layer=meta['Layer']) as source:
        for feat in source:
            geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if geom==None:
                continue
            if (prop['FIRE_YEAR']>=meta['Year Start']) & (prop['FIRE_YEAR']<=meta['Year End']):
                List.append(feat)

    gdf=gpd.GeoDataFrame.from_features(List,crs=meta['crs'])

    if type(roi)==dict:
        gdf=bc1hau.ClipGDF_ByROI(gdf,roi)

    return gdf

#%% Query VRI
# This excludes polygons outside the ROI
# Because VRI is so big, this only works with relatively small ROIs

def Query_VRI(meta,roi):

    List=[None]*int(1e5)
    cnt=0
    with fiona.open(meta['Path'],layer=meta['Layer']) as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue

            shp=shape(prop)
            bnd=shp.bounds # (W, S, E, N)

            # Check for overlap
            if (bnd[2]<roi['grd']['Extent'][0]) | (bnd[0]>roi['grd']['Extent'][1]) | (bnd[1]>roi['grd']['Extent'][3]) | (bnd[3]<roi['grd']['Extent'][2]):
                continue

            #if (shp.overlaps(roi_poly)==False) & (shp.within(roi_poly)==False):
            #    continue

            List[cnt]=feat
            cnt=cnt+1
            print('working')

    List=List[0:cnt-1]

    gdf=gpd.GeoDataFrame.from_features(List,crs=meta['crs'])

    if type(roi)==dict:
        gdf=bc1hau.ClipGDF_ByROI(gdf,roi)

    return gdf

#%% Query cutblocks

def Query_ConsolidatedCutblocks(meta,roi):
    List=[]
    with fiona.open(meta['Path'],layer=meta['Layer']) as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue
            if (prop['HARVEST_YEAR']>=meta['Year Start']) & (prop['HARVEST_YEAR']<=meta['Year End']):
                List.append(feat)
    gdf=gpd.GeoDataFrame.from_features(List,crs=meta['crs'])
    if type(roi)==dict:
        gdf=bc1hau.ClipGDF_ByROI(gdf,roi)
    return gdf

#%% QA

def Check():
    #path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930'
    with fiona.open(path,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            if prop['FIA_PROJECT_ID']=='FCI0000427':
                print('found')
                break

#%% Query pest DB

def GetPestSev():

    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'

    d={}
    with fiona.open(fin,layer='PEST_INFESTATION_POLY') as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            a=prop['PEST_SEVERITY_CODE']
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
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue

            if (prop['CAPTURE_YEAR']<Year0) | (prop['CAPTURE_YEAR']>Year1):
                continue

            if prop['PEST_SPECIES_CODE'] not in sp_cd:
                continue

            nam=prop['PEST_SPECIES_CODE']

            iSev=np.where(SevCD==prop['PEST_SEVERITY_CODE'])[0]
            if iSev.size==0:
                print(prop['PEST_SEVERITY_CODE'])

            sev=SevCD[iSev][0]

            iT=np.where(tv==prop['CAPTURE_YEAR'])[0]

            d[nam][sev][iT]=d[nam][sev][iT]+prop['AREA_HA']

    return tv,d

#%% Query all openings

#def GetOpeningInfo(ids):
#
#    # HARVESTING
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
#    List=[]
#    with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
#        for feat in source:
#            if prop==None:
#                continue
#            if np.isin(prop['OPENING_ID'],np.array(ids))==True:
#                List.append(feat)
#                if len(ids)==1:
#                    break
#    gdf_cut=gpd.GeoDataFrame.from_features(List,crs=crs)
#
#    # RESULTS FOREST COVER INVENTORY
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
#    List=[]
#    with fiona.open(fin,layer='RSLT_FOREST_COVER_INV_SVW') as source:
#        for feat in source:
#            if prop==None:
#                continue
#            if np.isin(prop['OPENING_ID'],np.array(ids))==True:
#                List.append(feat)
#    gdf_fcinv=gpd.GeoDataFrame.from_features(List,crs=crs)
#
##    # RESULTS FOREST COVER INVENTORY
##    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
##    List=[]
##    with fiona.open(fin,layer='RSLT_FOREST_COVER_SILV_SVW') as source:
##        for feat in source:
##            if prop==None:
##                continue
##            if np.isin(prop['OPENING_ID'],np.array(ids))==True:
##                List.append(feat)
##    gdf_fcsilv=gpd.GeoDataFrame.from_features(List,crs=crs)
#
#    # RESULTS OPENING LAYER
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
#    List=[]
#    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
#        for feat in source:
#            if prop==None:
#                continue
#            if np.isin(prop['OPENING_ID'],np.array(ids))==True:
#                List.append(feat)
#    gdf_op=gpd.GeoDataFrame.from_features(List,crs=crs)
#
#    return gdf_cut,gdf_fcinv,gdf_op

def GetSurveySpatial(Year1,Year2):

    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue
            if prop['SILV_BASE_CODE']!='SU':
                continue
            if prop['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(prop['ATU_COMPLETION_DATE'][0:4])
            if (Year<Year1) | (Year>Year2):
                continue
            if prop['GEOMETRY_Area']/10000>2000:
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

def GetSurveySpatial(gdf_bm):
    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(fin,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue
            if prop['SILV_BASE_CODE']!='SU':
                continue
            if prop['ATU_COMPLETION_DATE']==None:
                continue
            Year=int(prop['ATU_COMPLETION_DATE'][0:4])
            if Year<2018:
                continue
            if prop['GEOMETRY_Area']/10000>2000:
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

#%% Calculate area of shrub cover

def GetShrubs():

    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404\VRI.gdb'
    fiona.listlayers(fin)

    n=6000000
    d={}
    d['LC4']=np.array(['' for _ in range(n)],dtype=object)
    d['Area']=np.zeros(n)
    cnt=0
    with fiona.open(fin,layer='VEG_COMP_LYR_R1_POLY') as source:
        for feat in source:
            #geom=dict(feat['geometry'].items())
            prop=dict(feat['properties'].items())
            if prop==None:
                continue
            #break
            d['LC4'][cnt]=prop['BCLCS_LEVEL_4']
            d['Area'][cnt]=prop['GEOMETRY_Area']
            cnt=cnt+1
    for k in d.keys():
        d[k]=d[k][0:cnt]

    ikp=np.where( (d['LC4']=='ST') | (d['LC4']=='SL') )[0]
    np.sum(d['Area'][ikp])/10000/1e6

    ikp=np.where( (d['LC4']=='TC') | (d['LC4']=='TB') | (d['LC4']=='TM') )[0]
    np.sum(d['Area'][ikp])/10000/1e6

    return

