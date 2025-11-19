#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import pandas as pd
import fiona
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.bc1ha.bc1ha_update as up1ha
import fcgadgets.bc1ha.bc1ha_climate as c1ha
import fcgadgets.bc5k.bc5k_util as u5k
import fcgadgets.na1k.na1k_util as u1k
import fcgadgets.gaia.gaia_util as gaia
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
meta['Graphics']['Print Figures']='On'
meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BC5K'
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Climate (5k)
u5k.CreateRefGrid_BC5k(meta) # Create reference grid for bc5k
u5k.Process_Ndep_ISIMIP(meta) # Prepare nitrogen deposition
u5k.ProcessCRU(meta) # Import CRU
u5k.ProcessCMIP6(meta) # Import CMIP6
u5k.BiasCorrectCMIP6_FromCRU(meta) # Bias-correct CMIP6 anomalies against CRU
u5k.Process20thCenturyReanalysis(meta) # Import 20thC reanlaysis
u5k.ProcessNA1K(meta) # Import NA1K data
u5k.CompileEnvironmentalData(meta) # Compile environment

#%% Graphics
u5k.Plot_HistoricalAssessment(meta)
u5k.Plot_EnvironmentalTimeSeries(meta)


