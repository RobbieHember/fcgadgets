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

from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.macgyver import utilities_tile as tu
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
#TilesToRun_ij[0,:]=[10,8]
#TilesToRun_ij[0,:]=[11,7]
TilesToRun_ij[0,:]=[9,9]
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

#%% Open crosswalk between missing AT geometries and opening geometries
# If this doesn't work, you need to run the script that creates the crosswalk

#atu_mis=gu.ipickle(Paths['Results'] + '\\atu_mis.pkl')
#at_geo_from_op=gu.ipickle(Paths['Results'] + '\\at_geo_from_op.pkl')
missing_geo_atu_list=gu.ipickle(Paths['Results'] + '\\missing_geo_atu_list.pkl')
missing_geo_op_geos=gu.ipickle(Paths['Results'] + '\\missing_geo_op_geos.pkl')
missing_geo_fc_geos=gu.ipickle(Paths['Results'] + '\\missing_geo_fc_geos.pkl')

#%% Open missing FC layer geometries that were retreived from VRI
# Forest cover spatial reporting didn’t get turned on until 2004ish when RESULTS got turned on.  
# Previous forest cover reporting was delivered to the government as paper or as a PDF then 
# digitized directly into the VRI (or predecessors).  

dMisFC=gu.ipickle(Paths['Results'] + '\\missing_geo_fc_list.pkl')

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
    
    for iLyr in range(len(InvLyrInfo)):
        #iLyr=0
        
        # Loop through features in layer
        t_start=time.time()
    
        # Define path
        path0=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    
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
        with fiona.open(path0,layer=lyr_nam) as source:    
            for feat in source:
                
                # Extract attributes and geometry
                prp=feat['properties']
                geom=feat['geometry']
            
                #if (geom!=None): 
                #    break
            
                # Populate missing ATU layer geometry with geometry from 
                # OPENING or FC layer where possible.
                flg_geom_from_op=0
                flg_geom_from_fc=0
                if (lyr_nam=='RSLT_ACTIVITY_TREATMENT_SVW') & (geom==None):
                        
                    # Check to see if the opening is listed in the AT missing dictionary
                    indMis=np.where( (missing_geo_atu_list['ACTIVITY_TREATMENT_UNIT_ID']==prp['ACTIVITY_TREATMENT_UNIT_ID']) )[0]
                
                    if indMis.size>0:
                    
                        idx2fc=missing_geo_atu_list['IdxToFC'][indMis[0]]
                    
                        if len(idx2fc)>0:
            
                            # Use forest cover geometries
                    
                            geom={}
                            geom['coordinates']=[]
                            for i in range(len(idx2fc)):
                                geo0=missing_geo_fc_geos[prp['OPENING_ID']][idx2fc[i]]
                                if type(geo0)==dict:
                                    geo1=geo0['coordinates']
                                    geom['coordinates'].append(geo1[0])
                                else:
                                    for j in range(len(geo0)):
                                        geo1=geo0[j]['coordinates']
                                        geom['coordinates'].append(geo1[0])
                        
                            flg_geom_from_fc=1
                        
                        #elif prp['OPENING_ID'] in missing_geo_op_geos==True:
                        
                        elif len(missing_geo_op_geos[prp['OPENING_ID']])>0:
                        
                            # Use opening geometry
                    
                            geom={}
                            geom['coordinates']=[]
                            geo0=missing_geo_op_geos[prp['OPENING_ID']]
                            if type(geo0)==dict:
                                geo1=geo0['coordinates']
                                geom['coordinates'].append(geo1[0])
                            else:
                                for j in range(len(geo0)):
                                    geo1=geo0[j]['coordinates']
                                    geom['coordinates'].append(geo1[0])
    
                            flg_geom_from_op=1
                        
                    else:
                    
                        # Could not use either FC or openign layer
                        #print('Missing spatial could not be recovered')
                        pass
            
                # Populate missing FC layer geometry with geometry from 
                # OPENING or VRI where possible.
                if (lyr_nam=='RSLT_FOREST_COVER_INV_SVW') & (geom==None) | (lyr_nam=='RSLT_FOREST_COVER_SILV_SVW') & (geom==None):
                    
                    iMis_fc=np.where( (dMisFC['Unique Openings with Missing FC Geom']==prp['OPENING_ID']) )[0]
                    iMis_fc=iMis_fc[0]
                    
                    if dMisFC['Geom from VRI'][iMis_fc]!=None:
                    
                        D=np.zeros(len(dMisFC['Geom from VRI'][iMis_fc]))
                        for iV in range(len(dMisFC['Geom from VRI'][iMis_fc])):
                            D[iV]=np.abs(prp['SILV_POLYGON_AREA']-dMisFC['Geom from VRI'][iMis_fc][iV]['Hectares'])
                        iMinD=np.where(D==np.min(D))[0]
                        iMinD=iMinD[0]
                
                        geom=dMisFC['Geom from VRI'][iMis_fc][iMinD]
                    
                        # QA: Look at polygon
                        #geom1=invu.GetPolygonsFromFionaFeature(geom)                
                        #gdf=gpd.GeoDataFrame(geom1,crs=gdf_bm.crs) 
                
                # If non-veg, continue
                if InvLyrInfo[iLyr]['Layer Name']=='VEG_COMP_LYR_R1_POLY':
                    if (prp['BCLCS_LEVEL_1']=='N'): 
                        continue
            
                # Don't conitnue if no spatial data
                if (geom==None): 
                    continue        
        
                # Extract multipolygon
                flg_outside=0
                coords0=geom['coordinates']                       
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
                
                # rasterio burn function requires type field
                if 'type' not in geom:
                    geom['type']='MultiPolygon'
                
                # Rasterize feature geometry
                df0=gpd.GeoDataFrame({'ID':[1,1],'geometry':[geom,geom]})                
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
                        val=prp[fnam]
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
    path0=InvLyrInfo[iLyr]['Path'] + '\\' + InvLyrInfo[iLyr]['File Name']
    lyr_nam=InvLyrInfo[iLyr]['Layer Name']
    gdf_pl=gpd.read_file(path0,layer=lyr_nam)        
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

    #--------------------------------------------------------------------------
    # Add FC Archive to FC Inventory dictionary
    #--------------------------------------------------------------------------
    
    # Import inputs
    #meta=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\Metadata.pkl')
    #meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
    fcinv=gu.ipickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_INV_SVW.pkl')
    fcsilv=gu.ipickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_SILV_SVW.pkl')
    
    # Save old versions as alterntive names in case you need to redo this/troubleshoot a problem
    flg=1
    if flg==1:
        gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_INV_SVW_before_adding_archive.pkl',fcinv)
        gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_SILV_SVW_before_adding_archive.pkl',fcsilv)
    
    # Load initial layers
    flg=0
    if flg==1:
        fcinv=gu.ipickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_INV_SVW_before_adding_archive.pkl')
        fcsilv=gu.ipickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_SILV_SVW_before_adding_archive.pkl')
    
    # Run script
    # *** requires meta - just do a quick and dirty mock up of it and required content to run AddArchive script ***
    meta={}
    meta['Paths']=Paths
    meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
    fcinv['IdxToSXY']=fcinv['IdxToGrid']
    fcsilv['IdxToSXY']=fcsilv['IdxToGrid']
    fcinv,fcsilv=invu.ForestCover_AddArchive(meta,fcinv,fcsilv)
    fcinv['IdxToGrid']=fcinv['IdxToSXY']
    fcsilv['IdxToGrid']=fcsilv['IdxToSXY']    
    
    # Save
    gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_INV_SVW.pkl',fcinv)
    gu.opickle(infoTile['PathTileGeospatial'] + '\\RSLT_FOREST_COVER_SILV_SVW.pkl',fcsilv)

