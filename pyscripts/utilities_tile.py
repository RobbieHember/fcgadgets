'''

TILE UTILITIES

'''

import geopandas as gpd
import numpy as np
from CustomFunctions.gis import *
from CustomFunctions import basic_functions as bf

'''============================================================================
GET TILE SPATIAL INFORMATION
============================================================================'''

def GetTileIJ_Spatial(Path,iTile,jTile):
    
    zTSA=OpenGdal(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Admin\tsa.tif')
    
    gdf_tile=gpd.read_file(Path)
    
    zTile={}
    
    # Bounding box (minx,miny,maxx,maxy)
    bb=gdf_tile[(gdf_tile['i']==iTile) & (gdf_tile['j']==jTile)]['geometry'].bounds.values.flatten()
    zTile['Bounding Box']=bb

    # Get grid cells within ROI
    zTile['xv']=np.array([bb[0],bb[0],bb[2]-100,bb[2]-100,bb[0]])
    zTile['yv']=np.array([bb[1],bb[3]-100,bb[3]-100,bb[1],bb[1]])
    
    # Get the bounding box
    zTile['ix']=np.where((zTSA.X[0,:]>=np.min(zTile['xv'])) & (zTSA.X[0,:]<=np.max(zTile['xv'])))[0]
    zTile['iy']=np.where((zTSA.Y[:,0]>=np.min(zTile['yv'])) & (zTSA.Y[:,0]<=np.max(zTile['yv'])))[0]
    ind=np.ix_(zTile['ix'],zTile['iy'])
    zTile['X']=zTSA.X[ind]
    zTile['Y']=zTSA.Y[ind]

    # Define min and max values
    zTile['xmin'],zTile['ymin'],zTile['xmax'],zTile['ymax']=gdf_tile[(gdf_tile['i']==iTile) & (gdf_tile['j']==jTile)].total_bounds
    
    # Extent
    zTile['Extent']=[zTile['xmin'],zTile['xmax'],zTile['ymin'],zTile['ymax']]
    
    zTile['m']=zTile['ix'].size
    zTile['n']=zTile['iy'].size
    
    zTile=bf.Bunch(zTile)
    
    return zTile