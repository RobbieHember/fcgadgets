'''

FOREST CARBON BALANCE BY TILE - IMPORT INVENTORY

A single tile takes about 35 minutes to run. 

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
import gdal,osr
import pickle
from matplotlib import path
import matplotlib.pyplot as plt
import time
import gdal
from shapely.geometry import Polygon,Point
import rasterio
from rasterio import features
from rasterio.transform import from_origin

from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.utilities import utilities_tile as tu
from fcgadgets.cbrunner import cbrun_utilities

#%% Import paths

Paths=gu.ipickle(r'D:\Data\FCI_Projects\FCI_RollupByTile1\Paths.pkl')

#%% Import tile grid

# Import tile grid coordinates
gdf_tile=gpd.read_file(Paths['Project Root'] + '\\TilesBC.shp')

# Unique tiles (overlapping with land)
#ind=np.where(gdf_tile['WithLand']==1)[0]
#uTile=np.column_stack((gdf_tile.loc[ind,'i'].values,gdf_tile.loc[ind,'j'].values))

#%% Specify the tiles to run

TilesToRun_ij=np.zeros((1,2),dtype=int)
TilesToRun_ij[0,:]=[10,8]
#TilesToRun_ij[0,:]=[11,7]
#TilesToRun_ij[0,:]=[9,9]
#TilesToRun_ij[0,:]=[11,8]

#%% Import inventory layer information (names, variables, LUTs)

InvLyrInfo=invu.DefineInventoryLayersAndVariables()
for iLyr in range(len(InvLyrInfo)):
    lut=gu.ipickle(InvLyrInfo[iLyr]['Path'] + '\\LUTs_' + InvLyrInfo[iLyr]['Layer Name'] +'.pkl')
    for key in lut.keys():
        InvLyrInfo[iLyr]['LUT'][key]=lut[key]
    del lut

#%% Import look-up tables

lut_atu=gu.ipickle(Paths['Results'] + '\\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')

#%% Open crosswalk between missing AT geometries and opening geometries (if already exists)

atu_mis=gu.ipickle(Paths['Results'] + '\\atu_mis.pkl')
at_geo_from_op=gu.ipickle(Paths['Results'] + '\\at_geo_from_op.pkl')

#%% Compile inventory layers

for hTile in range(TilesToRun_ij.shape[0]):    
    
    #hTile=0
    
    # Define tile indices and tile name
    iTile=TilesToRun_ij[hTile,0]
    jTile=TilesToRun_ij[hTile,1]
    
    infoTile=tu.GetTileInfo(Paths,iTile,jTile)
    infoTile['X']=np.tile(np.reshape(infoTile['X'],(1,-1)),(infoTile['m'],1)) 
    infoTile['Y']=np.tile(np.reshape(infoTile['Y'],(-1,1)),(1,infoTile['n']))
    infoTile.keys()

    # Transform
    trf=from_origin(infoTile['xmin'],infoTile['ymax'],infoTile['cellsize'],infoTile['cellsize'])

    #------------------------------------------------------------------------------
    # Loop through layers
    #------------------------------------------------------------------------------
    
    for iLyr in range(7,8): #len(InvLyrInfo)):
        #iLyr=0
        
        # Loop through features in layer
        t_start=time.time()
    
        # Define path
        path=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    
        # Define layer
        lyr_nam=InvLyrInfo[iLyr]['Layer Name']
        print(lyr_nam)
        
        # Don't run this for planting - planting is done seperately below
        if lyr_nam=='RSLT_PLANTING_SVW':
            continue
        
        # Initialize index to inventory
        IdxToInv=[None]*infoTile['X'].size
        
        # Initialize dictionary
        L=9000000
        data={}
        data['IdxToGrid']=np.zeros(L,dtype=int)
        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
            if dtype=='<U20':
                data[fnam]=np.zeros(L,dtype=dtype)
            else:
                data[fnam]=-999*np.ones(L,dtype=dtype)       
        
        # Initialize counter
        cnt_inventory=0
        
        # Scan through layer file to extract selected variables, and convert string
        # variables to numeric based on LUTs
        with fiona.open(path,layer=lyr_nam) as source:    
            for feat in source:
                
                # Populate AT missing geometry with geometry from OPENING layer where
                # possible.
                if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW') & (feat['geometry']==None):
                    ind=np.where(atu_mis['OPENING_ID']==feat['properties']['OPENING_ID'])[0]
                    if ind.size>0:
                        feat['geometry']=at_geo_from_op[ind[0]]
                    else:
                        print('Checked for opening spatial, but no match')
            
                # If non-veg, continue
                if InvLyrInfo[iLyr]['Layer Name']=='VEG_COMP_LYR_R1_POLY':
                    if (feat['properties']['BCLCS_LEVEL_1']=='N'): 
                        continue
            
                # Don't conitnue if no spatial data
                if (feat['geometry']==None): 
                    continue        
        
                # Extract multipolygon
                flg_outside=0
                coords0=feat['geometry']['coordinates']                       
                for i in range(len(coords0)):
                    if flg_outside==1:
                        continue
                    coords1=coords0[i]                    
                    for j in range(len(coords1)):                        
                        if flg_outside==1:
                            continue
                        coords2=np.asarray(coords1[j])                    
                        x_feat=coords2[:,0]
                        y_feat=coords2[:,1]                    
                        # This should speed it up a bit
                        if np.max(x_feat)<infoTile['xmin']: 
                            flg_outside=1
                        if np.min(x_feat)>infoTile['xmax']: 
                            flg_outside=1
                        if np.max(y_feat)<infoTile['ymin']: 
                            flg_outside=1
                        if np.min(y_feat)>infoTile['ymax']: 
                            flg_outside=1 
                
                if flg_outside==1:
                    continue
                
                # Rasterize feature geometry
                df0=gpd.GeoDataFrame({'ID':[1,1],'geometry':[feat['geometry'],feat['geometry']]})                
                shapes=((geom,value) for geom,value in zip(df0.geometry,df0.ID))
                      
                z0=np.zeros(infoTile['X'].shape,dtype=float)
                burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=trf)
                #plt.imshow(burned)
                iKeep=np.where(burned.flatten()>0)[0]
                    
                # Only continue if overlap
                if iKeep.size==0: 
                    continue
                
                # Add attributes of overlapping grid cells to list
                for k in range(iKeep.size):
                    
                    iKeepK=iKeep[k]
                    
                    if IdxToInv[iKeepK]==None:
                        IdxToInv[iKeepK]={}
                        IdxToInv[iKeepK]['Index']=np.array([cnt_inventory])
                    else:
                        IdxToInv[iKeepK]['Index']=np.append(IdxToInv[iKeepK]['Index'],cnt_inventory)
                    
                    data['IdxToGrid'][cnt_inventory]=iKeepK
                    
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
                    
                    # Update counter
                    cnt_inventory=cnt_inventory+1
                
                #print(cnt_inventory)            
        
        # Truncate variables at cnt_inventory
        for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
            data[fnam]=data[fnam][0:cnt_inventory]
        data['IdxToGrid']=data['IdxToGrid'][0:cnt_inventory]
        
        # Time
        t_ela=time.time()-t_start
        print(t_ela)
    
        # Replace time string with numeric variables
        # Note: Don't do this for planting - planting will be re-compiled below
        # to include planting info for projects that reported no geometries in the
        # planting layer. The Strings will be fixed then.    
        data=invu.ExtractDateStringsFromRESULTS(lyr_nam,data)
        
        # Save        
        gu.opickle(infoTile['PathTileGeospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '.pkl',data)
        gu.opickle(infoTile['PathTileGeospatial'] + '\\' + InvLyrInfo[iLyr]['Layer Name'] + '_IdxToInv.pkl',IdxToInv)

    #--------------------------------------------------------------------------
    # Retrieve planting information 
    # Some projects did not report spatial planting info. Without the spatial info
    # in the planting layer, the initial import of the PL layer (above) will miss
    # planting info in some AT polygons. Use the ACTIVITY TREATMENT UNIT ID as
    # a crosswalk to retrieve all planting info for each activity.
    # (10 min)
    #--------------------------------------------------------------------------
    
    t0=time.time()
    
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
    atu=gu.ipickle(infoTile['PathTileGeospatial'] + '\\RSLT_ACTIVITY_TREATMENT_SVW.pkl')
    
    pl_code=lut_atu['SILV_BASE_CODE']['PL']
    
    # Only proceed if planting occurs
    ind_at=np.where(atu['SILV_BASE_CODE']==pl_code)[0]
    
    # Initialize dictionary
    L=1000000
    IdxToInv=[None]*infoTile['X'].size
    data={}
    data['IdxToGrid']=np.zeros(L,dtype=int)
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        if dtype=='<U20':
            data[fnam]=np.zeros(L,dtype=dtype)
        else:
            data[fnam]=-999*np.ones(L,dtype=dtype)
    
    # Populate planting layer
    cnt_inventory=0
    for i in range(ind_at.size):
        ind_pl=np.where(d_pl['ACTIVITY_TREATMENT_UNIT_ID']==atu['ACTIVITY_TREATMENT_UNIT_ID'][ind_at[i]])[0]
        
        ind_grd=atu['IdxToGrid'][ind_at[i]]
        
        for j in range(ind_pl.size):
            
            if IdxToInv[ind_grd]==None:
                IdxToInv[ind_grd]={}
                IdxToInv[ind_grd]['Index']=np.array([cnt_inventory])
            else:
                IdxToInv[ind_grd]['Index']=np.append(IdxToInv[ind_grd]['Index'],cnt_inventory)
            
            data['IdxToGrid'][cnt_inventory]=ind_grd
            
            for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
                val=d_pl[fnam][ind_pl[j]]
                if val==None:
                    continue
                if flag==0:
                    # Numeric variable, leave as is
                    data[fnam][cnt_inventory]=val
                elif flag==1:
                    # Convert string variable to numeric variable 
                    # based on LUT
                    data[fnam][cnt_inventory]=InvLyrInfo[iLyr]['LUT'][fnam][val]
            
            # Update counter
            cnt_inventory=cnt_inventory+1
    
    # Truncate variables at cnt_inventory
    for fnam,flag,dtype in InvLyrInfo[iLyr]['Field List']:
        data[fnam]=data[fnam][0:cnt_inventory]
    data['IdxToGrid']=data['IdxToGrid'][0:cnt_inventory]
    
    # Convert date string to numeric
    data=invu.ExtractDateStringsFromRESULTS(lyr_nam,data)
        
    t1=time.time()
    print(t1-t0)
       
    # Save    
    gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_PLANTING_SVW.pkl',data)
    gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_PLANTING_SVW_IdxToInv.pkl',IdxToInv)

