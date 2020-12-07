'''
PREPARE INVENTORY FROM POLYGONS
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
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities

#%% Project name

#project_name='FCI_RollupFCI_Inv'
#project_name='FertSummary'
project_name='ReforestationSummary'

#%% Define paths

Paths={}
Paths['Project']=r'D:\Data\FCI_Projects' + '\\' + project_name
#Paths['Project']=r'D:\Data\FCI_Projects\FertilizationSummary'
#Paths['Project']=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv'
Paths['Geospatial']=Paths['Project'] + '\\Geospatial'
Paths['QA']=Paths['Project'] + '\\QA'
#Paths['Figures']=r'G:\My Drive\Figures\Fertilization'
Paths['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20201203'
Paths['VRI']=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20200430'
Paths['Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430'
Paths['LandUse']=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20200706'
Paths['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

# Save
gu.opickle(Paths['Project'] + '\\Paths.pkl',Paths)

#%% Import FCI Admin Table
#Notes:
#FES recipients sometimes use the funding source code, "FES", for FCI-funded
#projects. To include them in the query, import the FCI project list.

flg=0
if flg==1:
    Paths['FCI DB File']='Z:\!Workgrp\Forest Carbon\Forest Carbon Initiative\Program\RollupProjects\Live Run\FCI_RollupProjects_01_Admin.xlsx'
    df_FCI=pd.read_excel(Paths['FCI DB File'],sheet_name='Sheet1')

    # Unique PP numbers in the database
    uPP=df_FCI['PP Number'].unique()

#%% Define sparse grid based on geotiff from BC1ha database

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')

# Grid interval (m)
ivl=zTSA['gt'][1]

# Domain (spans all of BC)
x_bc_all=zTSA['X'][0,:]
y_bc_all=zTSA['Y'][:,0]

# Load basemap (for CRS)
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

#%% Compile inventory data at sparse grid locations

# Get layer names and variables
InvLyrInfo=invu.DefineInventoryLayersAndVariables()

#%% Import look-up tables for each layer

for iLyr in range(len(InvLyrInfo)):
    lut=gu.ipickle(InvLyrInfo[iLyr]['Path'] + '\\LUTs_' + InvLyrInfo[iLyr]['Layer Name'] +'.pkl')
    for key in lut.keys():
        InvLyrInfo[iLyr]['LUT'][key]=lut[key]
    del lut

#%% Open crosswalk between missing AT geometries and opening geometries
# If this doesn't work, you need to run the script that creates the crosswalk

atu_mis=gu.ipickle(Paths['Results'] + '\\atu_mis.pkl')
at_geo_from_op=gu.ipickle(Paths['Results'] + '\\at_geo_from_op.pkl')

#%% Get layer

# Get planting lyaer id
for iLyr in range(len(InvLyrInfo)):
    if InvLyrInfo[iLyr]['Layer Name']=='RSLT_ACTIVITY_TREATMENT_SVW':
        break

# Define path
path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']

# Define layer name
lyr_nam=InvLyrInfo[iLyr]['Layer Name']
print(lyr_nam)
    
#%% Initialize list for multipolygons
# Keep multipolygons (the only way to get accurate area due to mismatches between
# opening area and treatment area)

atu_multipolygons=[None]*5000000
cnt_atu_multipolygons=0

#%% Initialize a geodataframe that will store each polygon

cnt_atu_polygons=0
cn_atu_polygons=['ID_atu_polygons','ID_atu_multipolygons','OPENING_ID','FIA_PROJECT_ID','SILV_BASE_CODE','ACTUAL_TREATMENT_AREA','Year','geometry']
list_atu_polygons=[None]*5000000
#gdf_atu_polygons=gpd.GeoDataFrame(columns=cn_atu_polygons,crs=gdf_bm.crs)
    
#%% Initialize sparse grid coordinate dictionary

L=20000000
sxy={}
sxy['x']=-999*np.ones(L)
sxy['y']=-999*np.ones(L)
sxy['ID_atu_multipolygons']=-999*np.ones(L,dtype=int)
sxy['ID_atu_polygons']=-999*np.ones(L,dtype=int)
cnt_sxy=0

#%% Get polygons from file

# Loop through features in layer
t0=time.time()

# Track amount of missing spatial
N_Missing_Spatial=0

trip=0

cnt=0

with fiona.open(path,layer=lyr_nam) as source:
        
    for feat in source:
            
        # Extract attributes and geometry
        prp=feat['properties']
        geom=feat['geometry']
        
        # In rare instances, there are no dates - add dummy numbers so that it 
        # does not crash - it appears to be so rare that it should not be a big
        # problem
        if prp['ATU_COMPLETION_DATE']==None:
            prp['ATU_COMPLETION_DATE']='99990000'
        
        if project_name=='FertSummary':
            
            # Fert Summary query
            if (prp['SILV_BASE_CODE']!='FE') | (prp['SILV_TECHNIQUE_CODE']!='CA') | (prp['RESULTS_IND']!='Y'):
                continue
            
        elif project_name=='ReforestationSummary':
            
            # Reforestation Summary query
            Year=int(prp['ATU_COMPLETION_DATE'][0:4])
            flg=0
            if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') & (Year>=1990):
                flg=1
            elif (prp['SILV_BASE_CODE']=='DS') & (prp['SILV_METHOD_CODE']!='LAYOT') & (prp['RESULTS_IND']=='Y') & (Year>=1990):
                flg=1
            else:
                flg=0
            
            if flg==0:
                continue
        
        elif project_name=='FCI Rollup':
            
            # FCI query
            if prp['RESULTS_IND']=='N':
                continue
            if (prp['SILV_FUND_SOURCE_CODE']!='FCE') & (prp['SILV_FUND_SOURCE_CODE']!='FCM') & (np.isin(prp['FIA_PROJECT_ID'],uPP)==False):
                continue
            if prp['SILV_BASE_CODE']=='SU':
                continue
        
        cnt=cnt+1
        print(cnt/365884)
        
        # Populate missing ATU layer geometry with geometry from 
        # OPENING layer where possible.
        flg_geom_from_op=0
        if (geom==None):
            ind=np.where(atu_mis['OPENING_ID']==prp['OPENING_ID'])[0]
            if ind.size>0:
                for i_ind in range(len(ind)):
                    if at_geo_from_op[ind[i_ind]]!=None:
                        geom=at_geo_from_op[ind[i_ind]]
                flg_geom_from_op=1
            else:
                print('Checked for opening spatial, but no match')
            
        # Don't conitnue if no spatial data
        if (geom==None): 
            N_Missing_Spatial=N_Missing_Spatial+1
            continue
        
        prp['ID_atu_multipolygons']=cnt_atu_multipolygons
        prp['Year']=int(prp['ATU_COMPLETION_DATE'][0:4])
        prp['geometry']=geom
        prp['GeomFromOpLyr']=flg_geom_from_op
        
        # Extract vector geometries and store as a list of dictionaries.
        # The list contains each polygon ("block"). Within each block, a
        # dictionary stores the x and y coordinates and the area.
        polys=[]
        polys_inner=[]
        coords0=geom['coordinates']
        for iPoly in range(len(coords0)):
            coords1=coords0[iPoly]
            
            x=[]; 
            y=[];
            for k in range(len(coords1[0])):
                x.append(coords1[0][k][0])
                y.append(coords1[0][k][1])
            
            l=[]
            for k in range(len(x)):
                l.append([x[k],y[k]])
            lp=Polygon(l)             
                    
            A=gis.PolyArea(x,y)/10000
                
            dct={}
            dct['Hectares']=A.copy()
            dct['x']=np.array(x).copy()
            dct['y']=np.array(y).copy()
            dct['x_Centroid']=lp.centroid.x
            dct['y_Centroid']=lp.centroid.y   
            dct['ID_atu_polygons']=cnt_atu_polygons
            polys.append(dct.copy())
            
            # Add to geodataframe
            dp0={}
            dp0['ID_atu_polygons']=cnt_atu_polygons
            dp0['ID_atu_multipolygons']=cnt_atu_multipolygons            
            dp0['OPENING_ID']=prp['OPENING_ID']
            dp0['FIA_PROJECT_ID']=prp['FIA_PROJECT_ID']            
            dp0['SILV_BASE_CODE']=prp['SILV_BASE_CODE']
            dp0['ACTUAL_TREATMENT_AREA']=prp['ACTUAL_TREATMENT_AREA']            
            dp0['Year']=prp['Year']
#            gdf_atu_polygons.loc[cnt_atu_polygons,'ID_atu_polygons']=cnt_atu_polygons
#            gdf_atu_polygons.loc[cnt_atu_polygons,'ID_atu_multipolygons']=cnt_atu_multipolygons            
#            gdf_atu_polygons.loc[cnt_atu_polygons,'OPENING_ID']=prp['OPENING_ID']
#            gdf_atu_polygons.loc[cnt_atu_polygons,'FIA_PROJECT_ID']=prp['FIA_PROJECT_ID']            
#            gdf_atu_polygons.loc[cnt_atu_polygons,'SILV_BASE_CODE']=prp['SILV_BASE_CODE']
#            gdf_atu_polygons.loc[cnt_atu_polygons,'ACTUAL_TREATMENT_AREA']=prp['ACTUAL_TREATMENT_AREA']            
#            gdf_atu_polygons.loc[cnt_atu_polygons,'Year']=prp['Year']
            
            gdf_outer=gpd.GeoDataFrame(crs=gdf_bm.crs)
            gdf_outer.loc[0,'geometry']=Polygon(coords1[0])
            if len(coords1)>1:
                # Adjust polygons if they have an inner ring
                for j in range(1,len(coords1)):
                    gdf_inner=gpd.GeoDataFrame(crs=gdf_bm.crs)
                    gdf_inner.loc[0,'geometry']=Polygon(coords1[j])                
                    try:
                        # This very rarely receives "'NoneType' object has no attribute 'intersection'"
                        # Not sure why
                        gdf_outer=gpd.overlay(gdf_outer,gdf_inner,how='difference')
                    
                        x=[]; 
                        y=[];
                        for k in range(len(coords1[j])):
                            x.append(coords1[j][k][0])
                            y.append(coords1[j][k][1])
                        dct={}
                        dct['x']=np.array(x).copy()
                        dct['y']=np.array(y).copy()
                        polys_inner.append(dct.copy())
                    except:
                        pass
            
            dp0['geometry']=gdf_outer.loc[0,'geometry']
            list_atu_polygons[cnt_atu_polygons]=dp0
            cnt_atu_polygons=cnt_atu_polygons+1        
            
        # Add polygon info to atu query polygon list
        prp['polys']=polys
        
        # Extract grid cells that intersect the AT polygons
        for iPoly in range(len(polys)):
            
            if polys[iPoly]['Hectares']==None:
                print('Missing geometry')
                continue
            
            x_feat=polys[iPoly]['x']
            y_feat=polys[iPoly]['y']
        
            # Index to a bounding box around the feature geometry
            ix=np.where((x_bc_all>=np.min(x_feat)-2*ivl) & (x_bc_all<=np.max(x_feat)+2*ivl))[0]
            iy=np.where((y_bc_all>=np.min(y_feat)-2*ivl) & (y_bc_all<=np.max(y_feat)+2*ivl))[0]
        
            if (ix.size==0) | (iy.size==0):
                continue
    
            x_bb=np.reshape(x_bc_all[ix],(1,ix.shape[0]))
            y_bb=np.reshape(y_bc_all[iy],(iy.shape[0],1))
            x_bb=np.tile(x_bb,(y_bb.shape[0],1))
            y_bb=np.tile(y_bb,(1,x_bb.shape[1]))
            y_bb=np.flip(y_bb,0)
        
            # Indicator for whether a cell is within the geometry
            InPol=gis.InPolygon(x_bb,y_bb,x_feat,y_feat)
          
            # Index to cells within feature geometry
            ikp=np.where(InPol==1)
        
            # If grid cell(s) fall within the project area, use those grid cells
            if ikp[0].size!=0:        
                x_ikp=np.atleast_1d(x_bb[ikp])
                y_ikp=np.atleast_1d(y_bb[ikp])
                
                # Remove cells within inner rings
                for iPoly_inner in range(len(polys_inner)):                    
                    InPol=gis.InPolygon(x_ikp,y_ikp,polys_inner[iPoly_inner]['x'],polys_inner[iPoly_inner]['y'])
                    ind=np.where(InPol==0)[0]
                    if ind.size>0:
                        x_ikp=x_ikp[ind]
                        y_ikp=y_ikp[ind]
                
                for i in range(x_ikp.size):
                    sxy['x'][cnt_sxy]=x_ikp[i]
                    sxy['y'][cnt_sxy]=y_ikp[i]                                
                    sxy['ID_atu_multipolygons'][cnt_sxy]=cnt_atu_multipolygons
                    sxy['ID_atu_polygons'][cnt_sxy]=prp['polys'][iPoly]['ID_atu_polygons']
                    cnt_sxy=cnt_sxy+1
            else:
                # If no grid cells fall within the project area, don't use the 
                # polygon centriod (doesn't work) - make a higher res grid and
                # re-search for interior cells.
                # (This occurs in road plantings)
                # Go all the way down to 5 m, 10 m misses some small features
                bin=[10,5,1,0.5]
                for iBin in range(len(bin)):                
                    #print('Going down to higher res!')
                    x_bb2=np.arange(np.min(x_feat),np.max(x_feat),bin[iBin])
                    y_bb2=np.arange(np.min(y_feat),np.max(y_feat),bin[iBin])
                    x_bb2=np.matlib.repmat(np.reshape(x_bb2,(1,x_bb2.size)),y_bb2.size,1)
                    y_bb2=np.matlib.repmat(np.reshape(y_bb2,(y_bb2.size,1)),1,x_bb2.size)
                    InPol=gis.InPolygon(x_bb2,y_bb2,x_feat,y_feat)
                    ikp=np.where(InPol==1)
                    if ikp[0].size>0:
                        break
                
                r=np.random.permutation(ikp[0].size)
                x_ikp=x_bb2[ikp]
                y_ikp=y_bb2[ikp]
                for iR in range(np.minimum(ikp[0].size,5)):
                    x_ikp0=np.atleast_1d(x_ikp[r[iR]])
                    y_ikp0=np.atleast_1d(y_ikp[r[iR]])
                  
                    # Remove cells within inner rings
                    for iPoly_inner in range(len(polys_inner)):
                        InPol=gis.InPolygon(x_ikp0,y_ikp0,polys_inner[iPoly_inner]['x'],polys_inner[iPoly_inner]['y'])
                        ind=np.where(InPol==0)[0]
                        if ind.size>0:
                            x_ikp0=x_ikp0[ind]
                            y_ikp0=y_ikp0[ind]
                    
                    sxy['x'][cnt_sxy]=x_ikp0
                    sxy['y'][cnt_sxy]=y_ikp0
                    sxy['ID_atu_multipolygons'][cnt_sxy]=cnt_atu_multipolygons
                    sxy['ID_atu_polygons'][cnt_sxy]=prp['polys'][iPoly]['ID_atu_polygons']
                    cnt_sxy=cnt_sxy+1
                
        # Update counter for multipolygon
        atu_multipolygons[cnt_atu_multipolygons]=prp
        cnt_atu_multipolygons=cnt_atu_multipolygons+1

#%% Truncate data

# Vector geometries
atu_multipolygons=atu_multipolygons[0:cnt_atu_multipolygons]

list_atu_polygons=list_atu_polygons[0:cnt_atu_polygons]

# Sparse sample points
for k in sxy.keys():
    sxy[k]=sxy[k][0:cnt_sxy]

#%% Convert list of atu polygons to gdf

gdf_atu_polygons=gpd.GeoDataFrame(list_atu_polygons,columns=cn_atu_polygons,crs=gdf_bm.crs)

#%% Remove duplicate cells
# Notes: This means that some values of ID_atu_multipolygons may not exist in
# the sxy dictionaries - they were only retained for one of the IDs of the 
# overlapping multipolygons. 
              
# Unique cells (ind_xy is the index to the first instance of each unique row)
uxy,ind_xy,inv_xy=np.unique(np.column_stack((sxy['x'],sxy['y'])),return_index=True,return_inverse=True,axis=0)

for k in sxy.keys():
    sxy[k]=sxy[k][ind_xy]

#%% Check lengths and make note of (expected) discrepency in length between uniquq
# IDs from sxy dictionary and unique IDs from list of multipolygons

flg=0
if flg==1:
    atu_multipolygons=gu.ipickle(Paths['Geospatial'] + '\\atu_multipolygons.pkl')
    print(len(atu_multipolygons))
    a=[]
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]!=None:
            a.append(atu_multipolygons[i]['ID_atu_multipolygons'])
    u_mp=np.unique(np.array(a))
    print(u_mp.size)

    u_sxy=np.unique(sxy['ID_atu_multipolygons'])
    print(u_sxy.size)

    c=0
    for i in range(len(atu_multipolygons)):
        if atu_multipolygons[i]['OPENING_ID']==1702534:        
            for j in range(len(atu_multipolygons[i]['polys'])):
                plt.plot(atu_multipolygons[i]['polys'][j]['x'],atu_multipolygons[i]['polys'][j]['y'],'b-')
                #c=c+1
            ind=np.where(sxy['ID_atu_multipolygons']==atu_multipolygons[i]['ID_atu_multipolygons'])[0]
            c=c+ind.size
            plt.plot(sxy['x'][ind],sxy['y'][ind],'.')    

#%% Add TSA ID

sxy['ID_TSA']=0*sxy['x']
for i in range(sxy['x'].size):
    print(i)
    ix=np.where(zTSA['X'][0,:]==sxy['x'][i])[0]
    iy=np.where(zTSA['Y'][:,0]==sxy['y'][i])[0]
    id_tsa=zTSA['Data'][iy,ix]
    if id_tsa.size>0:
        sxy['ID_TSA'][i]=id_tsa
    #ind=np.where(lut_tsa['VALUE']==id_tsa)[0]
    #lut_tsa['Name'][ind]

#ind=np.where(lut_tsa['Name']=='Okanagan TSA')[0]

#%% Save to file
    
# Save multipolygon file
gu.opickle(Paths['Project'] + '\\Inputs\\Geospatial\\atu_multipolygons.pkl',atu_multipolygons)

# Save sparse sample file
gu.opickle(Paths['Project'] + '\\Inputs\\Geospatial\\sxy.pkl',sxy)

# Save geodataframe of polygons to shape file
gdf_atu_polygons=gdf_atu_polygons.set_geometry('geometry')
gdf_atu_polygons.to_file(filename=Paths['Project'] + '\\Inputs\\Geospatial\\atu_polygons.shp',driver="ESRI Shapefile")

# Save sparse grid as shapefile
flg=1
if flg==1:

    points=[]
    for k in range(sxy['x'].size):
        points.append(Point(sxy['x'][k],sxy['y'][k]))
    gdf_sxy=gpd.GeoDataFrame({'geometry':points,
                              'ID_atu_multipolygons':sxy['ID_atu_multipolygons'],
                              'ID_atu_polygons':sxy['ID_atu_polygons'],
                              'ID_TSA':sxy['ID_TSA']})
    gdf_sxy.crs=gdf_bm.crs   
    gdf_sxy.to_file(Paths['Project'] + '\\Inputs\\Geospatial\\sxy.shp')

print((time.time()-t0)/60)


#%% Extract inventory information for sparse grid sample over polygons

garc.collect()

atu_multipolygons=gu.ipickle(Paths['Project'] + '\\Inputs\\Geospatial\\atu_multipolygons.pkl')
sxy=gu.ipickle(Paths['Project'] + '\\Inputs\\Geospatial\\sxy.pkl')

for iLyr in range(len(InvLyrInfo)):
    
    # Define path
    path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    
    # Define layer
    lyr_nam=InvLyrInfo[iLyr]['Layer Name']
    print(lyr_nam)
    
    # Don't run this for planting - planting is done seperately below
    if lyr_nam=='RSLT_PLANTING_SVW':
        continue
    
    # Don't run this for FC silv - no valuable variables
    if lyr_nam=='RSLT_FOREST_COVER_SILV_SVW':
        continue
     
    # Initialize index to inventory
    IdxToInv=[None]*sxy['x'].size
        
    # Initialize inventory dictionary
    L=sxy['x'].size+1000000
    data={}
    data['IdxToSXY']=np.zeros(L,dtype=int)
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        if dtype=='<U20':
            data[fnam]=np.zeros(L,dtype=dtype)
        else:
            data[fnam]=-999*np.ones(L,dtype=dtype)
    cnt_inventory=0
    
    # Loop through features in layer
    t_start=time.time()

    # Scan through layer file to extract selected variables, and convert string
    # variables to numeric based on LUTs
    with fiona.open(path,layer=lyr_nam) as source:
        
        for feat in source:
            
            # Extract attributes and geometry
            prp=feat['properties']
            geom=feat['geometry']
            
            # Fill missing AT spatial with OP spatial            
            if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW'):
                
                # Populate missing ATU layer geometry with geometry from 
                # OPENING layer where possible.
                flg_geom_from_op=0
                if (geom==None):
                    ind=np.where(atu_mis['OPENING_ID']==prp['OPENING_ID'])[0]
                    if ind.size>0:
                        geom=at_geo_from_op[ind[0]]
                        flg_geom_from_op=1
                    else:
                        print('Checked for opening spatial, but no match')
            
            # Only continue if spatial info exists
            if (geom==None) | (geom==[]):
                continue
            
            # Extract multipolygon
            coords0=geom['coordinates']
        
            # loop through multipolygon
            for i in range(len(coords0)):
                
                # Extract multipolygon
                coords1=coords0[i]
                
                # loop through multipolygon
                for j in range(len(coords1)):
                    
                    # Extract polygon
                    coords2=np.asarray(coords1[j])
                    x_feat=coords2[:,0]
                    y_feat=coords2[:,1]
                    
                    # This should speed it up a bit
                    if np.max(x_feat)<np.min(sxy['x']):
                        continue
                    if np.min(x_feat)>np.max(sxy['x']):
                        continue
                    if np.max(y_feat)<np.min(sxy['y']): 
                        continue
                    if np.min(y_feat)>np.max(sxy['y']): 
                        continue
        
                    # Isolate cells within bounding box of the feature polygon
                    iBB=np.where( (sxy['x']>=np.min(x_feat)-1000) & (sxy['x']<=np.max(x_feat)+1000) & (sxy['y']>=np.min(y_feat)-1000) & (sxy['y']<=np.max(y_feat)+1000) )[0]
                    
                    # Only continue if overlaop
                    if iBB.size==0:
                        continue
        
                    # Isolate cells within the polygon
                    InPoly=gis.InPolygon(sxy['x'][iBB],sxy['y'][iBB],x_feat,y_feat)
                    
                    iInPoly=np.where(InPoly==1)[0]
        
                    # Index to cells within polygon
                    iKeep=iBB[iInPoly]
                    
                    # Only continue if overlap
                    if iKeep.size==0:     
                        continue
                                        
                    # Add attributes of overlapping grid cells to list
                    for k in range(iKeep.size):
                        
                        iKeepK=iKeep[k]
                        
                        # Update index to inventory
                        if IdxToInv[iKeepK]==None:
                            IdxToInv[iKeepK]={}
                            IdxToInv[iKeepK]['Index']=np.array([cnt_inventory])
                        else:
                            IdxToInv[iKeepK]['Index']=np.append(IdxToInv[iKeepK]['Index'],cnt_inventory)
                        
                        data['IdxToSXY'][cnt_inventory]=iKeepK
                        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
                            val=feat['properties'][fnam]
                            if val!=None:
                                if flag==0:
                                    # Numeric variable, leave as is
                                    data[fnam][cnt_inventory]=val
                                elif flag==1:
                                    # Convert string variable to numeric variable 
                                    # based on LUT
                                    data[fnam][cnt_inventory]=InvLyrInfo[iLyr]['LUT'][fnam][val]
                
                        cnt_inventory=cnt_inventory+1
    
    # Time
    t_ela=time.time()-t_start
    print(t_ela)
    
    #--------------------------------------------------------------------------
    # Truncate data
    #--------------------------------------------------------------------------
    
    # Truncate index to inventory
    IdxToInv=IdxToInv[0:sxy['x'].size+1]
            
    # Truncate variables at cnt_inventory
    data['IdxToSXY']=data['IdxToSXY'][0:cnt_inventory]
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        data[fnam]=data[fnam][0:cnt_inventory]
    
    # Replace time string with numeric variables
    # Note: Don't do this for planting - planting will be re-compiled below
    # to include planting info for projects that reported no geometries in the
    # planting layer. The Strings will be fixed then.
    if lyr_nam!='RSLT_PLANTING_SVW':
        data=invu.ExtractDateStringsFromRESULTS(lyr_nam,data)      
       
    #--------------------------------------------------------------------------    
    # Save to file
    #--------------------------------------------------------------------------
    
    gu.opickle(Paths['Project'] + '\\Inputs\\Geospatial\\' + InvLyrInfo[iLyr]['Layer Name'] + '.pkl',data)
    gu.opickle(Paths['Project'] + '\\Inputs\\Geospatial\\' + InvLyrInfo[iLyr]['Layer Name'] + '_IdxToInv.pkl',IdxToInv)
        

#%% Retrieve planting information 
# Some projects did not report spatial planting info. Without the spatial info
# in the planting layer, the initial import of the PL layer (above) will miss
# planting info in some AT polygons. Use the ACTIVITY TREATMENT UNIT ID as
# a crosswalk to retrieve all planting info for each activity.
# (10 min)

t0=time.time()

# Import coordinates   
sxy=gu.ipickle(Paths['Geospatial'] + '\\sxy.pkl')
 
# Get planting lyaer id
for iLyr in range(len(InvLyrInfo)):
    if InvLyrInfo[iLyr]['Layer Name']=='RSLT_PLANTING_SVW':
        break
    
# Import geodataframe, drop geometry and convert to dict
path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
lyr_nam=InvLyrInfo[iLyr]['Layer Name']
gdf_pl=gpd.read_file(path,layer=lyr_nam)        
d_pl=gu.DataFrameToDict(gdf_pl.drop(columns='geometry'))
        
# Open the planting sparse grid
#pl=gu.ipickle(Paths['TileGeospatial'] + '\\RSLT_PLANTING_SVW.pkl')
    
# Get keys for planting layer
key_pl=[]
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    key_pl.append(fnam)
            
# Open AT sparse grid
atu=gu.ipickle(Paths['Geospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')

pl_code=InvLyrInfo[0]['LUT']['SILV_BASE_CODE']['PL']
    
# Only proceed if planting occurs
ind_at=np.where(atu['SILV_BASE_CODE']==pl_code)[0]

# Initialize dictionary
L=1000000
pl={}
pl['IdxToSXY']=np.zeros(L,dtype=int)
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    if dtype=='<U20':
        pl[fnam]=np.zeros(L,dtype=dtype)
    else:
        pl[fnam]=-999*np.ones(L,dtype=dtype)

# Initialize index to inventory
IdxToInv=[None]*sxy['x'].size
    
# Populate planting layer
cnt_inventory=0
for i in range(ind_at.size):
    ind_pl=np.where(d_pl['ACTIVITY_TREATMENT_UNIT_ID']==atu['ACTIVITY_TREATMENT_UNIT_ID'][ind_at[i]])[0]
    for j in range(ind_pl.size):           
        
        idx=atu['IdxToSXY'][ind_at[i]]
        
        pl['IdxToSXY'][cnt_inventory]=idx
        
        # Update index to inventory
        if IdxToInv[idx]==None:
            IdxToInv[idx]={}
            IdxToInv[idx]['Index']=np.array([cnt_inventory])
        else:
            IdxToInv[idx]['Index']=np.append(IdxToInv[idx]['Index'],cnt_inventory)
        
        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
            val=d_pl[fnam][ind_pl[j]]
            if val==None:
                continue
            if flag==0:
                # Numeric variable, leave as is
                pl[fnam][cnt_inventory]=val
            elif flag==1:
                # Convert string variable to numeric variable 
                # based on LUT
                pl[fnam][cnt_inventory]=InvLyrInfo[iLyr]['LUT'][fnam][val]
        # Update counter
        cnt_inventory=cnt_inventory+1
    
# Truncate variables at cnt_inventory
pl['IdxToSXY']=pl['IdxToSXY'][0:cnt_inventory]        
for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
    pl[fnam]=pl[fnam][0:cnt_inventory]
    
# Convert date string to numeric
pl=invu.ExtractDateStringsFromRESULTS(lyr_nam,pl)
       
# Save    
gu.opickle(Paths['Geospatial'] + '\\RSLT_PLANTING_SVW.pkl',pl)
gu.opickle(Paths['Geospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl',IdxToInv)

t1=time.time()
print(t1-t0)

