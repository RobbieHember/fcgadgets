'''

PREPARE CLIMATE ANOMALIES


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

geos={}
geos['Path']=r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies'

# Import raster grids
zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\FAIB_Standard.tif')

# Import BC boundary (for crs)
gdf_bc_boundary=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

# Define regular grid sampling frequency
#geos['rgsf']=40 # 4 km
geos['rgsf']=400 # 40 km

# Extract subgrid
X=zRef['X'][0::geos['rgsf'],0::geos['rgsf']].flatten()
Y=zRef['Y'][0::geos['rgsf'],0::geos['rgsf']].flatten()

# Elevation
zE=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif')
Elev=zE['Data'][0::geos['rgsf'],0::geos['rgsf']].flatten()

points=[]
for i in range(X.size):
    points.append(Point(X[i],Y[i]))
gdf_xy=gpd.GeoDataFrame({'geometry':points,'Elev':Elev})
gdf_xy.crs=gdf_bc_boundary.crs
gdf_xy.to_file(geos['Path'] + '\\geos.geojson',driver='GeoJSON')

# Convert to lat long
gdf_ll=gdf_xy.to_crs(4326)

fl=np.arange(0,len(gdf_ll))
Data=np.column_stack( (fl,fl,gdf_ll['geometry'].y.values,gdf_ll['geometry'].x.values,Elev) )
df=pd.DataFrame(columns=['ID1','ID2','Lat','Lon','Elev'],data=Data)
df.to_csv(geos['Path'] + '//coords.csv',index=False)

#%% Import climate data

gdf_xy=gpd.read_file(r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies\geos.geojson',driver='GeoJSON')

# Historical data
dH=gu.ReadCSV(r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies\From ClimateBC\coords_1901-2021MSY.csv')

id=np.unique(dH['ID1'])
tv=gu.tvec('m',1901,2021)
vL=['Tave','PPT','RH','Rad']
d={}
for v in vL:
    d[v]=np.zeros((tv.shape[0],id.size))
    for mo in range(1,13):
        if mo<10:
            nam='0' + str(mo)
        else:
            nam=str(mo)
        iT=np.where(tv[:,1]==mo)[0]
        d[v][iT,:]=np.reshape(np.reshape(dH[v + nam],(-1,1)), (int(tv.shape[0]/12),id.size) )
















# #%% Import monthly CMIP6 native grids and put them in project spatial reference system

# def List_CMIP6_Availability(meta):

#     # Path to CMIP6 data
#     PathCMIP6=r'C:\Users\rhember\Documents\Data\CMIP6'

#     # Model names
#     ms=gu.ReadExcel(PathCMIP6 + '\\ModelSummary.xlsx')

#     # Exclusion of models
#     namM=ms['Model']
#     fsM=ms['File structure']
#     fN=ms['F Number']

#     # Future scenarios
#     namS=['historical','ssp245','ssp585']

#     # Variable list
#     namV=['tasmin','tasmax','pr','rsds']

#     data=[]

#     for iM in range(namM.size):

#         nM=namM[iM]

#         for iR in range(0,100):

#             flg_s=np.zeros(len(namS))

#             for iS in range(len(namS)):

#                 nS=namS[iS]

#                 flg_v=np.zeros(len(namV))

#                 for iV in range(len(namV)):

#                     nV=namV[iV]

#                     if fsM[iM]=='Grouped':
#                         if nS=='historical':
#                             pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_185001-201412.nc'
#                         else:
#                             pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_201501-210012.nc'
#                     else:
#                         if nS=='historical':
#                             pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_185101-185112.nc'
#                         else:
#                             pin=PathCMIP6 + '\\' + nM + '\\' + nV + '_Amon_' + nM + '_' + nS + '_r' + str(iR+1) + 'i1p1f' + str(fN[iM]) + '_gn_201501-201512.nc'

#                     if os.path.isfile(pin)==True:
#                         flg_v[iV]=1

#                 if np.sum(flg_v)==len(namV):
#                     flg_s[iS]=1

#             if np.sum(flg_s)==len(namS):
#                 d={}
#                 d['Model']=nM
#                 d['Run']=iR+1
#                 data.append(d)



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

