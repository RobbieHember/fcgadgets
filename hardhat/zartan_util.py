#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import time
import warnings
import pandas as pd
import statsmodels.formula.api as smf
import lmfit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun_preprocess as prep
import fcgadgets.cbrunner.cbrun_postprocess as post
import fcgadgets.cbrunner.cbrun as cbr
import fcgadgets.macgyver.util_fcs_graphs as ufcs
import fcgadgets.macgyver.util_demo as udem
import fcgadgets.gaia.gaia_util as gaia
warnings.filterwarnings("ignore")
gp=gu.SetGraphics('Manuscript')



#%%

#%%