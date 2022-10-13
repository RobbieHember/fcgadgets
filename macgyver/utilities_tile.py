'''

TILE UTILITIES

'''

import geopandas as gpd
import numpy as np
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis

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

def GetTileInfo(metaTile,iTile,jTile):

    zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

    # Import geodatabase
    gdf_tile=gpd.read_file(metaTile['Paths']['Project Root'] + '\\TilesBC.geojson')

    # Initialize dictionary
    geos={}

    # Define string name
    if iTile<10:
        iTileS='0' + str(iTile)
    else:
        iTileS=str(iTile)
    if jTile<10:
        jTileS='0' + str(jTile)
    else:
        jTileS=str(jTile)
    geos['Name']=iTileS + jTileS

    # Define tile path
    geos['PathTileGeospatial']=metaTile['Paths']['Project Root'] + '\\TileBC_' + geos['Name'] + '\\Geospatial'

    # Index to tile geodataframe
    ind_gdf_tile=(gdf_tile['i']==iTile) & (gdf_tile['j']==jTile)

    # Tile ID
    geos['ID']=gdf_tile.loc[ind_gdf_tile,'ID'].values[0]

    # Cell size
    geos['Cellsize']=gdf_tile['Cellsize'][ind_gdf_tile].values

    # Tile width
    geos['MaxWidth']=gdf_tile['MaxWidth'][ind_gdf_tile].values

    # Bounding box (minx,miny,maxx,maxy)
    bb=gdf_tile[ind_gdf_tile]['geometry'].bounds.values.flatten()
    geos['Bounding Box']=bb

    # Get grid cells within ROI
    geos['xv']=np.array([bb[0],bb[0],bb[2]-100,bb[2]-100,bb[0]])
    geos['yv']=np.array([bb[1],bb[3]-100,bb[3]-100,bb[1],bb[1]])

    # Define min and max values
    geos['xmin'],geos['ymin'],geos['xmax'],geos['ymax']=gdf_tile[ind_gdf_tile].total_bounds

    geos['X']=np.arange(geos['xmin'],geos['xmax'],geos['Cellsize'])
    geos['Y']=np.arange(geos['ymin'],geos['ymax'],geos['Cellsize'])

    # Extent
    geos['Extent']=[geos['xmin'],geos['xmax'],geos['ymin'],geos['ymax']]

    # Dimensions
    geos['m']=geos['X'].size
    geos['n']=geos['Y'].size

    # Limits
    geos['xlim']=[geos['xmin'],geos['xmax']]
    geos['ylim']=[geos['ymin'],geos['ymax']]

    return geos

#%% GET ROAD VECTORS FOR EACH TILE

def GetTileRoads():

    # Load roads
    gdf_road=gpd.read_file(r'Z:\!Workgrp\Forest Carbon\Data\Basemaps\roads.shp')

    # Import paths
    Paths=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\Paths.pkl')

    # Import tile grid

    # Import tile grid coordinates
    gdf_tile=gpd.read_file(Paths['Project'] + '\\TilesBC.geojson')

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

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.geojson')

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

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.geojson')

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

    tl=gpd.read_file(r'D:\Data\FCI_Projects\FCI_RollupByTile1\TilesBC.geojson')

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