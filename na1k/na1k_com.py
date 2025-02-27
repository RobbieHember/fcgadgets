'''
NORTH AMERICA 1 KM COMMAND
'''
#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_general as gu
import fcgadgets.na1k.na1k_util as u1k
import fcgadgets.gaia.gaia_util as gaia
gp=gu.SetGraphics('Manuscript')

#%% Initialize 1 km North American project structure
metaNA=u1k.Init()

#%% Create reference grid from political boundaries
u1k.CreateReferenceGrid(metaNA)

#%% Process DEM (for input into ClimateNA)
u1k.ProcessDEM(metaNA)

#%% Import GSOD
u1k.ImportGSOD(metaNA)

#%% Climate normals (1971-2000)
u1k.ImportNormalsFromClimateNA(metaNA) # Import monthly climate normals from ClimateNA (temp and precip)
u1k.QA_ClimateNA(metaNA) # QA temperature from ClimateNA
u1k.CalcSolarRadiationNormals(metaNA) # Downward solar radiation
u1k.CalcSaturationVapourPressureNormal(metaNA) # Saturation vapour pressure normal
u1k.CalcActualVapourPressureNormalFromTemps(metaNA) # Actual vapour pressure (from ClimateNA temperatures)
u1k.CalcActualVapourPressureNormalBiasCorrected(metaNA) # Atual vapour pressure (bias corrected)
u1k.CalcVapourPressureDeficitNormal(metaNA) # Vapour pressure deficit (using bias-corrected vapour pressure)

#%% Download daily data from Environment Canada
u1k.DownloadEnvironmentCanadaClimateData(metaNA)

#%% Gap-fill monthly data from Environment Canada
u1k.GapFillDailyECCD(metaNA)

#%% Import USHCN monthly data
u1k.ImportUSHCN(metaNA)

#%% Manually download monthly NCEP/NCAR reanalysis data
# https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface_gauss/

#%% Calculate monthly temperature anomalies
u1k.CalcTemperatureAnomalies(metaNA)

#%% Calculate monthly precipitation anomalies
u1k.CalcPrecipAnomalies(metaNA)

#%% Calculate monthly vapour pressure deficit anomalies
u1k.CalcVapourPressureDefictAnomaliesFromCRU(metaNA) # Long-term historical data
u1k.CalcVapourPressureDefictAnomaliesFromReanalysis(metaNA,2023) # Calculate monthly vapour pressure deficit anomalies (most recent year)

#%% Calculate monthly surface water balance
u1k.CalcSurfaceWaterBalance(metaNA)

#%% Calculate seasonal summarites
u1k.CalcSummariesBasicVariables(metaNA)
u1k.CalcSummariesWB(metaNA)

#%% Calculate time series averages by ecozone
u1k.CalcTimeSeriesByEZ(metaNA)

#%%
