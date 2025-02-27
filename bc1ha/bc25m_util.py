
#%% Import modules

import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import time
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fiona
import pyproj
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb

#%% Initialize project
def Init_bc25m(meta):
	# Add paths
	meta['Paths']['bc25m']=r'C:\Data\BC25m'
	meta['Paths']['bc25m Ref Grid']=meta['Paths']['bc25m'] + '\\LandCover\\LandMask.tif'

	return meta

#%%
def DefineReferenceGrid(meta):
	# Create reference grid
	# Manually copy BC1ha LandMask to bc25m database folder, then resample below
	gis.ResampleRaster(meta['Paths']['bc25m Ref Grid'],4)
	return

#%% Build LUTs for variables with categorical data from source geodatabases
# Takes 7.5 hours
#fiona.listlayers(r'C:\\Users\\rhember\\Documents\\Data\\Geodatabases\\LandUse\\20230501\\LandUse.gdb')

# *** Based on BC1ha database ***

#%% Look up tables

# *** Based on BC1ha database ***

