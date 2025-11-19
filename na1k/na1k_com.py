'''
NORTH AMERICA 1 KM COMMAND
'''
#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.na1k.na1k_util as u1k
import fcgadgets.gaia.gaia_util as gaia
gp=gu.SetGraphics('Manuscript')

#%% Initialize 1 km North American project structure
metaNA=u1k.Init()
zRef=gis.OpenGeoTiff(metaNA['Paths']['na1k Ref Grid'])

#%% Create reference grid from political boundaries
u1k.CreateReferenceGrid(metaNA)

#%% Process DEM (for input into ClimateNA)
u1k.ProcessDEM(metaNA)

#%% Import GSOD
u1k.ImportGSOD(metaNA)

#%% Update CRU
# Source: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/cruts.2503051245.v4.09/

#%% Climate normals (1971-2000)
u1k.ImportNormalsFromClimateNA(metaNA) # Import monthly climate normals from ClimateNA (temp and precip)
u1k.QA_ClimateNA(metaNA) # QA temperature from ClimateNA
u1k.CalcSolarRadiationNormals(metaNA) # Downward solar radiation
u1k.CalcSaturationVapourPressureNormal(metaNA) # Saturation vapour pressure normal
u1k.CalcActualVapourPressureNormalFromTemps(metaNA) # Actual vapour pressure (from ClimateNA temperatures)
u1k.CalcActualVapourPressureNormalFromCRU(metaNA) # Not in use
u1k.CalcActualVapourPressureNormalBiasCorrected(metaNA) # Atual vapour pressure (bias corrected)
u1k.CalcVapourPressureDeficitNormal(metaNA) # Vapour pressure deficit (using bias-corrected vapour pressure)

#%% Download daily data from Environment Canada
u1k.DownloadEnvironmentCanadaClimateData(metaNA)

#%% Gap-fill monthly data from Environment Canada
u1k.GapFillDailyECCD(metaNA)

#%% Look at subset of station data with long-term records
u1k.LongTermStationAssessment(metaNA)

#%% Import USHCN monthly data (crashing)
#u1k.ImportUSHCN(metaNA)

#%% Import AHCCD
d=u1k.ImportAHCCD(metaNA)

#%% Manually download monthly NCEP/NCAR reanalysis data
# https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface_gauss/

#%% Calculate monthly temperature anomalies
u1k.CalcTemperatureAnomaliesFromStationInterp(metaNA)
u1k.CalcTemperatureAnomaliesFromCRU(metaNA)

#%% Calculate monthly precipitation anomalies
u1k.CalcPrecipAnomaliesFromStationInterp(metaNA)

#%% Monthly solar radiation (NCEP/NCAR reanalysis)
u1k.CalcRadiationAnamlies(metaNA)

#%% Calculate monthly vapour pressure deficit anomalies
u1k.CalcVapourPressureDefictAnomaliesFromCRU(metaNA) # Long-term historical data
u1k.CalcVapourPressureDefictAnomaliesFromReanalysis(metaNA,2023) # Calculate monthly vapour pressure deficit anomalies (most recent year)

#%% Calculate basic variable seasonal summaries
u1k.CalcSummariesBasicVariables(metaNA)

#%% Calculate monthly surface water balance
u1k.CalcSurfaceWaterBalance(metaNA)
u1k.CalcSummariesWB(metaNA) # Calculate seasonal summaries

#%% Calculate time series averages by ecozone
u1k.CalcSummaryTimeSeries(metaNA)

#%% NDEP stats
u1k.Calc_NDEP_Stats()

#%% PCA of interpolated 50k air temperature grids
# iM=0
# tva=np.arange(1951,2025,1)
# m,n=130,128
# z=np.zeros( (tva.size,m,n) )
# for i,yr in enumerate(tva):
# 	z[i,:,:]=gu.ipickle(meta['Paths']['na1k'] + '\\Monthly\\Anomalies\\50k\\na1k_tmean_anom_' + str(yr) + '_' + str(iM+1) + '.tif')

# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler

# n_comp=70
# scaler=StandardScaler()
# X=z.reshape(z.shape[0],-1)
# X_scaled=scaler.fit_transform(X)
# pca=PCA(n_components=n_comp)
# X_pca=pca.fit_transform(X_scaled)
# ve=pca.explained_variance_
# pcv=pca.components_
# s=pca.transform(X_scaled)

# #
# i=60;
# plt.close('all');
# fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(24,10))
# ax[0].matshow(pcv[i,:].reshape(m,n))
# ax[1].plot(tva,gu.movingave(s[:,i],10,'center'),'-bo')
