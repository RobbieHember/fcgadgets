'''
PREPARE CLIMATE ANOMALIES
'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy.matlib as ml
import matplotlib.pyplot as plt
import time
from scipy.interpolate import griddata
import gc as garc
from scipy.io import loadmat
from shapely.geometry import Polygon,Point
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.gaia.gaia_utilities as clmu

#%% Define grid

geos={}
geos['Path']=r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies'

# Import raster grids
zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

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

# Export for ClimateBC
gdf_ll=gdf_xy.to_crs(4326) # Convert to lat long
fl=np.arange(0,len(gdf_ll))
Data=np.column_stack( (fl,fl,gdf_ll['geometry'].y.values,gdf_ll['geometry'].x.values,Elev) )
df=pd.DataFrame(columns=['ID1','ID2','Lat','Lon','Elev'],data=Data)
df.to_csv(geos['Path'] + '//coords.csv',index=False)

#%% Import climate data

gdf_xy=gpd.read_file(r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies\geos.geojson',driver='GeoJSON')

# Historical data
dH0=gu.ReadCSV(r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies\From ClimateBC\coords_1901-2021MSY.csv')

id=np.unique(dH0['ID1'])
tvH=gu.tvec('m',1901,2021)
vL=['Tave','PPT','RH','Rad']
dH={}
for v in vL:
    dH[v]=np.zeros((tvH.shape[0],id.size))
    for mo in range(1,13):
        if mo<10:
            nam='0' + str(mo)
        else:
            nam=str(mo)
        iT=np.where(tvH[:,1]==mo)[0]
        dH[v][iT,:]=np.reshape(np.reshape(dH0[v + nam],(-1,1)), (int(tvH.shape[0]/12),id.size) )

# Convert RH to VPD
es=clmu.GetEstar(dH['Tave'])
ea=(dH['RH']/100)*es
dH['vpd']=(es-ea)*10

dH['ta']=dH['Tave']
dH['rswd']=dH['Rad']
dH['prcp']=dH['PPT']

# Future data
dF0=gu.ReadCSV(r'C:\Users\rhember\Documents\Data\Climate\Gaia\Anomalies\From ClimateBC\coords_CanESM5_ssp370_2011-2100MP.csv')

idF=np.unique(dF0['ID1'])
tvF=gu.tvec('m',2011,2100)
vL=['Tave','Tmin','Tmax','PPT','Rad']
dF={}
for v in vL:
    dF[v]=np.zeros((tvF.shape[0],idF.size))
    for mo in range(1,13):
        if mo<10:
            nam='0' + str(mo)
        else:
            nam=str(mo)
        iT=np.where(tvF[:,1]==mo)[0]
        dF[v][iT,:]=np.reshape(np.reshape(dF0[v + nam],(-1,1)), (int(tvF.shape[0]/12),idF.size) )

# Predict future vpd from temperatures
es_min=clmu.GetEstar(dF['Tmin'])
es_max=clmu.GetEstar(dF['Tmax'])
ea=clmu.GetEstar(dF['Tmin']-2.5)
dF['vpd']=((es_min+es_max)/2-ea)*10

dF['RH']=0*dF['Tave']

dF['ta']=dF['Tave']
dF['rswd']=dF['Rad']
dF['prcp']=dF['PPT']

#%% Combine historical and future

tvC=gu.tvec('m',1801,2151)
dC={}
for v in ['ta','prcp','vpd','RH']:
    dC[v]=np.zeros((tvC.shape[0],dH[v].shape[1]))
    iT=np.where( (tvC[:,0]>=tvH[0,0]) & (tvC[:,0]<=tvH[-1,0]) )[0]
    dC[v][iT,:]=dH[v]

    iT1=np.where( (tvC[:,0]>=2021) & (tvC[:,0]<=tvF[-1,0]) )[0]
    iT2=np.where( (tvF[:,0]>=2021) )[0]
    dC[v][iT1,:]=dF[v][iT2,:]

    # Gap filling
    for mo in range(12):
        iTh=np.where( (tvC[:,0]<=1900) & (tvC[:,1]==mo+1) )[0]
        iTf=np.where( (tvC[:,0]>2100)  & (tvC[:,1]==mo+1) )[0]

        iTh2=np.where( (tvC[:,0]>=1901) & (tvC[:,0]<=1931) & (tvC[:,1]==mo+1) )[0]
        iTf2=np.where( (tvC[:,0]==2095) & (tvC[:,0]<=2100) & (tvC[:,1]==mo+1) )[0]
        for j in range(dC[v].shape[1]):
            dC[v][iTh,j]=np.mean(dC[v][iTh2,j])
            dC[v][iTf,j]=dC[v][iTf2,j]

# Radiation from normals
dC['rswd']=0*dC['ta'].copy()
for mo in range(12):
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
    dC['rswd'][mo::12,:]=z['Data'][0::geos['rgsf'],0::geos['rgsf']].flatten()


#%% Seasonal values

indS=np.where(np.mean(dC['ta'],axis=0)>-100)[0]

tv=np.arange(1801,2101,1)
dCs={}
for v in ['ta','prcp','rswd','vpd','RH']:
    dCs[v]=np.zeros((tv.size,dC[v].shape[1]))
    for iT in range(tv.size):
        ind=np.where( (tvC[:,0]==tv[iT]) & (tvC[:,1]>=5) & (tvC[:,1]<=9) )[0]
        dCs[v][iT]=np.mean(dC[v][ind,:],axis=0)

plt.close('all')
plt.plot(tv,np.mean(dCs['RH'][:,indS],axis=1),'-bo')

#%% Downscaling

#%% Calculate water balance

vi=dC.copy()

# Modify variables
con=clmu.HydroMetCon()

vi['rswn']=(1-con['Albedo']['Forest Coniferous'])*vi['rswd']
vi['LAI']=5.0
vi['Gs']=0.010
vi['Ga']=0.058

# Parameters
par={}
par['Method']='Combined'
par['Ws_max']=200.0 # mm
par['Tmin']=-3.0
par['Tmax']=3.0
par['Ei_FracMax']=0.15
par['Ei_ALMax']=5.0
par['ETp Method']='Penman-Monteith'
par['Daily_Interval']=5
par['Include Rainfall Fraction']='No'

# Calculate water balance iterating to get to a steady-state
vo=clmu.WBM(par,vi)











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

