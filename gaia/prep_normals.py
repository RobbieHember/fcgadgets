'''
PREPARE CLIMATE NORMALS
'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import scipy.io as spio
from scipy.io import loadmat
from shapely.geometry import Polygon,Point
import gc as garc
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.gaia.gaia_utilities as clmu
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import project config

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)

# Import land mask
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_9.tif')

plt.matshow(z['Data'],clim=[0,200])

#%% Prep input for ClimateBC

# ind=np.where( (zRef['Data']==1) )
# srs=gis.ImportSRSs()
# lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'][ind],zRef['Y'][ind])
# zE=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\elevation.tif')
# d={'ID1':np.arange(lon.size),'ID2':np.arange(lon.size),'Lat':lat,'Lon':lon,'Elev':zE['Data'][ind]}
# df=pd.DataFrame.from_dict(d)
# df.to_csv(r'C:\ClimateBC\InputFiles\bc.csv')

#%% Prepare monthly water balance normals

iMask=np.where(zRef['Data']>0)

# Initialize warm-season means
zWs=zRef.copy()
zWs['Data']=np.zeros(zRef['Data'].shape)
zETa=zRef.copy()
zETa['Data']=np.zeros(zRef['Data'].shape)
zETp=zRef.copy()
zETp['Data']=np.zeros(zRef['Data'].shape)

N_bin=15
bin=np.linspace(0,iMask[0].size,N_bin,dtype='int32')

for iBin in range(N_bin-1):
    print(iBin)
    indBin=np.arange(bin[iBin],bin[iBin+1],1).astype('int32')

    # Import climate data
    vi={}
    vi['ta']=np.zeros((12,indBin.size))
    vi['prcp']=np.zeros((12,indBin.size))
    vi['rswd']=np.zeros((12,indBin.size))
    vi['vpd']=np.zeros((12,indBin.size))
    for mo in range(12):
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1\BC1ha_tmean_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        vi['ta'][mo,:]=z['Data'][ iMask[0][indBin],iMask[1][indBin] ]

        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1\BC1ha_prcp_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        vi['prcp'][mo,:]=z['Data'][ iMask[0][indBin],iMask[1][indBin] ]

        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1\BC1ha_rswd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        vi['rswd'][mo,:]=z['Data'][ iMask[0][indBin],iMask[1][indBin] ]

        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Monthly\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1\BC1ha_vpd_mon_norm_1971to2000_si_hist_v1_' + str(mo+1) + '.tif')
        vi['vpd'][mo,:]=z['Data'][ iMask[0][indBin],iMask[1][indBin] ]

    N_repeat=2
    for v in vi.keys():
        vi[v]=np.repeat(vi[v],N_repeat,axis=0)

    # Modify variables
    con=clmu.HydroMetCon()
    vi['rswn']=(1-con['Albedo']['Forest Coniferous'])*vi['rswd']
    del vi['rswd']
    garc.collect()
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

    # Fix
    for v in vo.keys():
        vo[v]=vo[v][(N_repeat-1)*12:,:]

    # Populate
    zWs['Data'][ iMask[0][indBin],iMask[1][indBin] ]=np.mean(vo['Ws'][5:9,:],axis=0)
    zETa['Data'][ iMask[0][indBin],iMask[1][indBin] ]=np.mean(vo['ETa'][5:9,:],axis=0)
    zETp['Data'][ iMask[0][indBin],iMask[1][indBin] ]=np.mean(vo['ETp'][5:9,:],axis=0)

    garc.collect()

#%% Save

zWs['Data']=zWs['Data'].astype('int16')

zCWD=zWs.copy()
zCWD['Data']=zETp['Data']-zETa['Data']
zCWD['Data']=zCWD['Data']/np.mean(con['DIM'][5:9])*100
zCWD['Data']=zCWD['Data'].astype('int16')

zETp['Data']=zETp['Data']/np.mean(con['DIM'][5:9])*100
zETp['Data']=zETp['Data'].astype('int16')

gis.SaveGeoTiff(zWs,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
gis.SaveGeoTiff(zCWD,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_cwd_gs_norm_1971to2000.tif')
gis.SaveGeoTiff(zETp,r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_etp_gs_norm_1971to2000.tif')

#%%

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
plt.matshow(z['Data'])

#%%













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

