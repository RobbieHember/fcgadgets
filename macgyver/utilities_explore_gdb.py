'''

EXPLORE GDB FILE

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
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Set path and layer

# Path to forest cover archive geodatabase
path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210401\Results.gdb'
#path=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20210401\VRI.gdb'
path=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'

# List layers
fiona.listlayers(path)

# Set layer
lyr='RSLT_ACTIVITY_TREATMENT_SVW'
#lyr='VEG_COMP_LYR_R1_POLY'

#%% Loop through gdb

with fiona.open(path,layer=lyr) as source:
    for feat in source:
        break

N=0
N_all=0
with fiona.open(path,layer='RSLT_FOREST_COVER_INV_SVW') as source:
    for feat in source:
        break
        
    N_all=N_all+1
        if feat['geometry']==None:
            N=N+1
        
        #prp=feat['properties']
        #if prp['OPENING_ID']==15241:
        #    if feat['geometry']!=None:
        #        print('found')


