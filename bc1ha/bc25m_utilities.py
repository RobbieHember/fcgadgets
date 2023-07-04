
#%% Import modules

import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fiona
import time
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_query_gdb as qgdb

#%% Initialize project

def Init():

    meta={}

    # Set paths
    meta['Paths']={}
    meta['Paths']['bc25m']=r'C:\Users\rhember\Documents\Data\BC25m'
    meta['Paths']['bc25m Ref Grid']=meta['Paths']['bc25m'] + '\\LandCoverUse\\LandMask.tif'
    meta['Paths']['bc1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
    meta['Paths']['bc1ha Ref Grid']=meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandMask.tif'
    meta['Paths']['GDB']={}
    meta['Paths']['GDB']['GDB']=r'C:\Users\rhember\Documents\Data\Geodatabases'
    meta['Paths']['GDB']['LandCover']=r'C:\Users\rhember\Documents\Data\Geodatabases\LandCover\20230607\LandCover.gdb'
    meta['Paths']['GDB']['LandUse']=r'C:\Users\rhember\Documents\Data\Geodatabases\LandUse\20230501\LandUse.gdb'
    meta['Paths']['GDB']['Disturbance']=r'C:\Users\rhember\Documents\Data\Geodatabases\Disturbances\20230501\Disturbances.gdb'
    meta['Paths']['GDB']['Results']=r'C:\Users\rhember\Documents\Data\Geodatabases\Results\20230430\Results.gdb'
    meta['Paths']['GDB']['VRI']=r'C:\Users\rhember\Documents\Data\Geodatabases\VRI\20230401\VRI.gdb'
    meta['Paths']['Model']={}
    meta['Paths']['Model']['Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
    meta['Paths']['Model']['Taz Datasets']=r'C:\Users\rhember\Documents\Data\Taz Datasets'

    # Initiate geospatial info
    meta['Geos']={}

    # Import variable info
    meta['Geos']['Variable Info']=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\BC1ha Raster Variable List.xlsx')

    # Import coordinate reference system
    gdf_bm=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='TA_REGIONAL_DISTRICTS_SVW')
    meta['Geos']['crs']=gdf_bm.crs

    return meta

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 7.5 hours
#fiona.listlayers(r'C:\\Users\\rhember\\Documents\\Data\\Geodatabases\\LandUse\\20230501\\LandUse.gdb')

# *** Based on BC1ha database ***

#%% Look up tables

# *** Based on BC1ha database ***

