'''

TILE UTILITIES

'''

import geopandas as gpd
import numpy as np
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis

#%% 

def TileIdicesToName(iTile,jTile):
    
    if iTile<10:
        inam='0' + str(iTile)
    else:
        inam=str(iTile)
    
    if jTile<10:
        jnam='0' + str(jTile)
    else:
        jnam=str(jTile)
    
    name=inam + jnam
    
    return name

#%% GET TILE SPATIAL INFORMATION

def GetTileInfo(Paths,iTile,jTile):
    
    zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
    
    # Import geodatabase
    gdf_tile=gpd.read_file(Paths['Project Root'] + '\\TilesBC.shp')
    
    # Initialize dictionary
    info={}
    
    # Define string name
    if iTile<10: 
        iTileS='0' + str(iTile)
    else: 
        iTileS=str(iTile)    
    if jTile<10: 
        jTileS='0' + str(jTile)
    else: 
        jTileS=str(jTile)        
    info['Name']=iTileS + jTileS
    
    # Define tile path    
    info['PathTileGeospatial']=Paths['Project Root'] + '\\TileBC_' + info['Name'] + '\\Inputs\\Geospatial'    
    
    # Index to tile geodataframe
    ind_gdf_tile=(gdf_tile['i']==iTile) & (gdf_tile['j']==jTile)
    
    # Cell size
    info['cellsize']=gdf_tile['cellsize'][ind_gdf_tile].values
    
    # Tile width
    info['tilewidth']=gdf_tile['tilewidth'][ind_gdf_tile].values
    
    # Bounding box (minx,miny,maxx,maxy)
    bb=gdf_tile[ind_gdf_tile]['geometry'].bounds.values.flatten()
    info['Bounding Box']=bb

    # Get grid cells within ROI
    info['xv']=np.array([bb[0],bb[0],bb[2]-100,bb[2]-100,bb[0]])
    info['yv']=np.array([bb[1],bb[3]-100,bb[3]-100,bb[1],bb[1]])

    # Define min and max values
    info['xmin'],info['ymin'],info['xmax'],info['ymax']=gdf_tile[ind_gdf_tile].total_bounds
    
    info['X']=np.arange(info['xmin'],info['xmax'],info['cellsize'])
    info['Y']=np.arange(info['ymin'],info['ymax'],info['cellsize'])
    
    # Extent
    info['Extent']=[info['xmin'],info['xmax'],info['ymin'],info['ymax']]
    
    # Dimensions
    info['m']=info['X'].size
    info['n']=info['Y'].size
    
    # Limits
    info['xlim']=[info['xmin'],info['xmax']]
    info['ylim']=[info['ymin'],info['ymax']]
    
    #info=gu.Bunch(info)
    
    return info


#%% GET ROAD VECTORS FOR EACH TILE

def GetTileRoads():
    
    # Load roads
    gdf_road=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\roads.shp')
    
    # Import paths    
    Paths=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\Paths.pkl')

    # Import tile grid

    # Import tile grid coordinates
    gdf_tile=gpd.read_file(Paths['Project'] + '\\TilesBC.shp')

    for k in range(len(gdf_tile)):
        if gdf_tile.loc[k,'WithLand']==0:
            continue
        #if (gdf_tile.loc[k,'i']==9) & (gdf_tile.loc[k,'j']==9):
        #    break
        
        infoTile=tu.GetTileInfo(Paths,gdf_tile.loc[k,'i'],gdf_tile.loc[k,'j'])
        gdf_road_ij=gdf_road.cx[infoTile['xmin']:infoTile['xmax'],infoTile['ymin']:infoTile['ymax']]
        gdf_road_ij.to_file(infoTile['PathTileGeospatial'] + '\\Roads.shp',driver='ESRI Shapefile')
    
    return


#%% Get list of tiles intersecting major watersheds
    
def GetIntersectingTilesAndMajorWatersheds():

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.shp')

    wat=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='BC_MAJOR_WATERSHEDS')    
    u=wat['MAJOR_WATERSHED_SYSTEM'].unique()
    
    dID={}
    dName={}
    for iU in range(u.size):
        wat0=wat[wat['MAJOR_WATERSHED_SYSTEM']==u[iU]]
        wat0=wat0.reset_index()
        id=[]
        for iTile in range(len(tl)):
            for i in range(len(wat0)):
                if tl.loc[iTile,'geometry'].intersects(wat0.loc[i,'geometry']):
                    id.append(tl.loc[iTile,'ID'])
        dID[u[iU]]=np.unique(id)
        
        nam=[]
        for i in range(dID[u[iU]].size):
            iTile=tl.loc[dID[u[iU]][i],'i']
            jTile=tl.loc[dID[u[iU]][i],'j']
            
            # Define string name
            if iTile<10: 
                iTileS='0' + str(iTile)
            else: 
                iTileS=str(iTile)    
            if jTile<10: 
                jTileS='0' + str(jTile)
            else: 
                jTileS=str(jTile)     
            nam.append(iTileS + jTileS)
        dName[u[iU]]=np.unique(nam)
    
    flg=0
    if flg==1:
        list(dName.keys())
        dName['Chilcotin River']
    
    return dID,dName

#%% Get list of tiles intersecting TSAs
    
def GetIntersectingTilesAndTSAs():

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.shp')

    wat=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')
    
    u=wat['Name'].unique()
    dID={}
    dName={}
    for iU in range(u.size):
        wat0=wat[wat['Name']==u[iU]]
        wat0=wat0.reset_index()
        id=[]
        for iTile in range(len(tl)):
            for i in range(len(wat0)):
                if tl.loc[iTile,'geometry'].intersects(wat0.loc[i,'geometry']):
                    id.append(tl.loc[iTile,'ID'])
        dID[u[iU]]=np.unique(id)
    
        nam=[]
        for i in range(dID[u[iU]].size):
            iTile=tl.loc[dID[u[iU]][i],'i']
            jTile=tl.loc[dID[u[iU]][i],'j']
            
            # Define string name
            if iTile<10: 
                iTileS='0' + str(iTile)
            else: 
                iTileS=str(iTile)    
            if jTile<10: 
                jTileS='0' + str(jTile)
            else: 
                jTileS=str(jTile)     
            nam.append(iTileS + jTileS)
        dName[u[iU]]=np.unique(nam)        
    
    return dID,dName

#%% Get list of tiles intersecting TSAs
    
def GetIntersectingTilesAndDistricts():

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.shp')

    wat=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='ADM_NR_DISTRICTS_SP')
    
    u=wat['DISTRICT_NAME'].unique()
    dID={}
    dName={}
    for iU in range(u.size):
        wat0=wat[wat['DISTRICT_NAME']==u[iU]]
        wat0=wat0.reset_index()
        id=[]
        for iTile in range(len(tl)):
            for i in range(len(wat0)):
                if tl.loc[iTile,'geometry'].intersects(wat0.loc[i,'geometry']):
                    id.append(tl.loc[iTile,'ID'])
        dID[u[iU]]=np.unique(id)
    
        nam=[]
        for i in range(dID[u[iU]].size):
            iTile=tl.loc[dID[u[iU]][i],'i']
            jTile=tl.loc[dID[u[iU]][i],'j']
            
            # Define string name
            if iTile<10: 
                iTileS='0' + str(iTile)
            else: 
                iTileS=str(iTile)    
            if jTile<10: 
                jTileS='0' + str(jTile)
            else: 
                jTileS=str(jTile)     
            nam.append(iTileS + jTileS)
        dName[u[iU]]=np.unique(nam)
    
    return dID,dName