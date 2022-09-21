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
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1hau

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

                #--------------------------------------------------------------
                # Filters
                #--------------------------------------------------------------

                #if feat['geometry']==None:
                #    continue

                if meta['Select Openings'].size>0:
                    if np.isin(feat['properties']['OPENING_ID'],meta['Select Openings'])==False:
                        continue

                if type(meta['ROI'])==dict:
                    # Using ROI
                    s=shape(feat['geometry'])
                    if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                        continue

                if meta['SBC'].size!=0:
                    if np.isin(feat['properties']['SILV_BASE_CODE'],meta['SBC'])==False:
                        continue
                    # Planting (isolate real planting)
                    if meta['SBC']=='PL':
                        if (feat['properties']['SILV_TECHNIQUE_CODE']=='SE') | (feat['properties']['SILV_TECHNIQUE_CODE']=='CG') | (feat['properties']['SILV_METHOD_CODE']=='LAYOT'):
                            continue

                if meta['FSC'].size!=0:
                    if np.isin(feat['properties']['SILV_FUND_SOURCE_CODE'],meta['FSC'])==False:
                        continue

                if meta['Layer']=='RSLT_ACTIVITY_TREATMENT_SVW':
                    if (feat['properties']['RESULTS_IND']!='Y') | (feat['properties']['ATU_COMPLETION_DATE']==None):
                        continue
                    feat['properties']['Year']=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])

                if cnt==0:

                    # Initialize
                    for k in feat['properties'].keys():
                        try:
                            d[k]=np.zeros(n)
                            d[k][cnt]=feat['properties'][k]
                        except:
                            d[k]=np.array(['' for _ in range(n)],dtype=object)
                            d[k][cnt]=feat['properties'][k]
                else:

                    for k in feat['properties'].keys():
                        try:
                            d[k][cnt]=feat['properties'][k]
                        except:
                            d[k]=np.array(['' for _ in range(n)],dtype=object)
                            d[k][cnt]=feat['properties'][k]
                cnt=cnt+1

        # Truncate
        for k in d.keys():
            d[k]=d[k][0:cnt]

    else:

        #======================================================================
        # Keep geometries
        #======================================================================

        List=[None]*int(1e5)
        cnt=0
        with fiona.open(meta['Path'],layer=meta['Layer']) as source:
            for feat in source:

                if feat['geometry']==None:
                    continue

                if meta['Select Openings'].size>0:
                    if np.isin(feat['properties']['OPENING_ID'],meta['Select Openings'])==False:
                        continue

                if type(roi)==dict:
                    # Using ROI
                    s=shape(feat['geometry'])
                    if (s.overlaps(roi_poly)==False) & (s.within(roi_poly)==False):
                        continue

                if meta['SBC'].size!=0:
                    if np.isin(feat['properties']['SILV_BASE_CODE'],meta['SBC'])==False:
                        continue
                    # Planting (isolate real planting)
                    if meta['SBC']=='PL':
                        if (feat['properties']['SILV_TECHNIQUE_CODE']=='SE') | (feat['properties']['SILV_TECHNIQUE_CODE']=='CG') | (feat['properties']['SILV_METHOD_CODE']=='LAYOT'):
                            continue

                if meta['FSC'].size!=0:
                    if np.isin(feat['properties']['SILV_FUND_SOURCE_CODE'],meta['FSC'])==False:
                        continue

                if meta['Layer']=='RSLT_ACTIVITY_TREATMENT_SVW':
                    if (feat['properties']['RESULTS_IND']!='Y') | (feat['properties']['ATU_COMPLETION_DATE']==None):
                        continue
                    feat['properties']['Year']=int(feat['properties']['ATU_COMPLETION_DATE'][0:4])

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
    with fiona.open(meta['Path'],layer=meta['Layer']) as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            if (feat['properties']['FIRE_YEAR']>=meta['Year Start']) & (feat['properties']['FIRE_YEAR']<=meta['Year End']):
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
            if feat['geometry']==None:
                continue

            shp=shape(feat['geometry'])
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
            if feat['geometry']==None:
                continue
            if (feat['properties']['HARVEST_YEAR']>=meta['Year Start']) & (feat['properties']['HARVEST_YEAR']<=meta['Year End']):
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
            if feat['properties']['FIA_PROJECT_ID']=='FCI0000427':
                print('found')
                break


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

#def GetOpeningInfo(ids):
#
#    # HARVESTING
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
#    List=[]
#    with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
#        for feat in source:
#            if feat['geometry']==None:
#                continue
#            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
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
#            if feat['geometry']==None:
#                continue
#            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
#                List.append(feat)
#    gdf_fcinv=gpd.GeoDataFrame.from_features(List,crs=crs)
#
##    # RESULTS FOREST COVER INVENTORY
##    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
##    List=[]
##    with fiona.open(fin,layer='RSLT_FOREST_COVER_SILV_SVW') as source:
##        for feat in source:
##            if feat['geometry']==None:
##                continue
##            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
##                List.append(feat)
##    gdf_fcsilv=gpd.GeoDataFrame.from_features(List,crs=crs)
#
#    # RESULTS OPENING LAYER
#    fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
#    List=[]
#    with fiona.open(fin,layer='RSLT_OPENING_SVW') as source:
#        for feat in source:
#            if feat['geometry']==None:
#                continue
#            if np.isin(feat['properties']['OPENING_ID'],np.array(ids))==True:
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

def GetSurveySpatial(gdf_bm):
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
            if feat['geometry']==None:
                continue
            #break
            d['LC4'][cnt]=feat['properties']['BCLCS_LEVEL_4']
            d['Area'][cnt]=feat['properties']['GEOMETRY_Area']
            cnt=cnt+1
    for k in d.keys():
        d[k]=d[k][0:cnt]

    ikp=np.where( (d['LC4']=='ST') | (d['LC4']=='SL') )[0]
    np.sum(d['Area'][ikp])/10000/1e6

    ikp=np.where( (d['LC4']=='TC') | (d['LC4']=='TB') | (d['LC4']=='TM') )[0]
    np.sum(d['Area'][ikp])/10000/1e6

    return





