'''

BC50KM CLIMATE

Build annual climate variables for a 50 km grid over BC.

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
import matplotlib.pyplot as plt
import time
from scipy.interpolate import griddata
import pyproj
import netCDF4 as nc
import gzip
from scipy.io import loadmat
from shapely.geometry import Polygon,Point
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_climate as clu

#%% Define grid

meta={}
meta['Paths']={}
meta['Paths']['Project']=r'D:\Data\FCI_Projects\Climate25km'

# Save
gu.opickle(meta['Paths']['Project'] + '\\Metadata.pkl',meta)

# Import raster grids
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

# Initialize geospatial info structure
geos={}

# Define regular grid sampling frequency
geos['rgsf']=500 # 50 km

# Extract subgrid
zTSA['Data']=zTSA['Data'][0::geos['rgsf'],0::geos['rgsf']]

# Populate geospatial structure with grid attributes
geos['X']=zTSA['X'][0::geos['rgsf'],0::geos['rgsf']]
geos['Y']=zTSA['Y'][0::geos['rgsf'],0::geos['rgsf']]
geos['m'],geos['n']=geos['X'].shape
geos['Data']=np.zeros((geos['m'],geos['n']),dtype=np.int8)

# Define additional inclusion criteria
geos['iMask']=np.where( (zTSA['Data']>0) & (zTSA['Data']<255) )

# Revise mask
geos['Data'][geos['iMask']]=1

# Generate sparse grid
geos['Sparse']={}
geos['Sparse']['X']=geos['X'][geos['iMask']]
geos['Sparse']['Y']=geos['Y'][geos['iMask']]
geos['Sparse']['ID_TSA']=zTSA['Data'][geos['iMask']]

# Save to pickle file
gu.opickle(meta['Paths']['Project'] + '\\geos.pkl',geos)

#%% Import monthly CMIP6 native grids and put them in project spatial reference system

def List_CMIP6_Availability(meta):

    # Path to CMIP6 data
    PathCMIP6=r'C:\Users\rhember\Documents\Data\CMIP6'
    
    # Model names
    ms=gu.ReadExcel(PathCMIP6 + '\\ModelSummary.xlsx')
    
    # Exclusion of models 
    namM=ms['Model']
    fsM=ms['File structure']
    fN=ms['F Number']
    
    # Future scenarios
    namS=['historical','ssp245','ssp585']
    
    # Variable list
    namV=['tasmin','tasmax','pr','rsds']
     
    data=[]
    
    for iM in range(namM.size):
        
        nM=namM[iM]
          
        for iR in range(0,100):
        
            flg_s=np.zeros(len(namS))
            
            for iS in range(len(namS)):
                
                nS=namS[iS]                
            
                flg_v=np.zeros(len(namV))
                
                for iV in range(len(namV)):
                
                    nV=namV[iV]
                    
                    if fsM[iM]=='Grouped':
                        if nS=='historical':
                            pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_185001-201412.nc'
                        else:
                            pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_201501-210012.nc'
                    else:
                        if nS=='historical':
                            pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_185101-185112.nc'
                        else:
                            pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_201501-201512.nc'
                    
                    if os.path.isfile(pin)==True:
                        flg_v[iV]=1
                        
                if np.sum(flg_v)==len(namV):
                    flg_s[iS]=1
            
            if np.sum(flg_s)==len(namS):
                d={}
                d['Model']=nM
                d['Run']=iR+1
                data.append(d)
                    
                

#%% Import CMIP6 data

d_cmip=clu.Import_CMIP6(meta,geos)

plt.close('all')
plt.plot(d_cmip['tasmin']['CESM2-FV2'][1])

# Calculate mean air temperature
data1={}
data1['tmean']

#%% Import CRU data

d_cru=clu.Import_CRU(meta,geos)

#%% Add bias-corrected future simulations to historical observations from CRU-TS

d_clm=clu.CombineMonthly(meta,d_cru,d_cmip)

#%% Run water balance model

d_clm=clu.RunWBM(d_clm)

