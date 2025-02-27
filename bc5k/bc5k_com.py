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
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Climate (5k)
u5k.CreateRefGrid_BC5k(meta) # Create reference grid for bc5k
u5k.Process_Ndep_ISIMIP(meta) # Prepare nitrogen deposition
u5k.ImportCRU(meta) # Import CRU
u5k.ImportCMIP6(meta) # Import CMIP6
u5k.BiasCorrectCMIP6_FromCRU(meta) # Bias-correct CMIP6 anomalies against CRU
u5k.Import20thCenturyReanalysis(meta) # Import 20thC reanlaysis

#%%
