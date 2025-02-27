'''
REMOVED FATE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import util_general as gu

#%% Graphics parameters

gp=gu.SetGraphics('Manuscript')

#%% Import BaseCase HWP parameters
dDef0=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_RemovedFate.xlsx')

# Variables to keep
vrL=['RemovedMerchToChipperMill','RemovedMerchToPulpMill','RemovedMerchToPelletMill','RemovedMerchToLumberMill',
    'RemovedMerchToPlywoodMill','RemovedMerchToOSBMill','RemovedMerchToMDFMill','RemovedMerchToPolePostMill',
    'RemovedMerchToShakeShingleMill','RemovedMerchToLogExport','RemovedMerchToIPP','RemovedMerchToFirewood',
    'RemovedNonMerchToChipperMill','RemovedNonMerchToPulpMill','RemovedNonMerchToPelletMill','RemovedNonMerchToLumberMill','RemovedNonMerchToPlywoodMill',
    'RemovedNonMerchToOSBMill','RemovedNonMerchToMDFMill','RemovedNonMerchToPolePostMill','RemovedNonMerchToShakeShingleMill',
    'RemovedNonMerchToLogExport','RemovedNonMerchToIPP','RemovedNonMerchToFirewood',
    'RemovedDeadStemToChipperMill','RemovedDeadStemToPulpMill',
    'RemovedDeadStemToPelletMill','RemovedDeadStemToLumberMill','RemovedDeadStemToPlywoodMill','RemovedDeadStemToOSBMill',
    'RemovedDeadStemToMDFMill','RemovedDeadStemToPolePostMill','RemovedDeadStemToShakeShingleMill','RemovedDeadStemToLogExport',
    'RemovedDeadStemToIPP','RemovedDeadStemToFirewood']

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

# Region list
regL=['Coast',
      'Interior',
      'CT1',
      'CT2',
      'CT3',
      'GFS22',
      'Energy Production',
      'Burn Uneconomic',
      'Leave Uneconomic Standing',
      'Sawlogs and Pellets',
      'Energy Production Pellets']

# Scenario list
scnL=['BaseCase','S1','S2']

#%% Initialize with BaseCase parameters

d={}
d['Year']=tv
for scn in scnL:
    d[scn]={}
    for reg in regL:
        d[scn][reg]={}
        for i in range(dDef0['Name'].size):
            if np.isin(dDef0['Name'][i],vrL)==True:
                d[scn][reg][ dDef0['Name'][i] ]=np.ones(tv.size)*dDef0[reg][i]

#%% Scenario 1

# Shift away from pulp+pellets to LumberMills
# *** Probably not realistic! ***

sc='S1'

# Fractional reduction in direct transfers from forest ecosystem to pulp mill
# and pellet mill
fChange=0.5

# Phase in period
tStart=2021
tEnd=2030

for reg in regL:

    #--------------------------------------------------------------------------
    # Merch
    #--------------------------------------------------------------------------

    # Change in Pulp input
    nam='RemovedMerchToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Plot assumptions
    #plt.plot(tv,dS1[nam],'b-')
    # plt.plot(tv,1-yDelta,'b-')

    # Change in pellet input
    nam='RemovedMerchToPelletMill'
    y0_Pellet=d['BaseCase'][reg][nam][0].copy()
    y1_Pellet=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pellet,y1_Pellet,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pellet
    yD_Pellet=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Fractions of fibre sent to LumberMill, plywood, OSB
    fLumberMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2

    # Add to LumberMill input
    nam='RemovedMerchToLumberMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fLumberMill*(yD_Pulp+yD_Pellet)

    # Add to Plywood input
    nam='RemovedMerchToPlywoodMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)

    # Add to OSB input
    nam='RemovedMerchToOSBMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fOSBMill*(yD_Pulp+yD_Pellet)

    #--------------------------------------------------------------------------
    # Deads
    #--------------------------------------------------------------------------

    # Assumed change in Pulp input
    nam='RemovedDeadStemToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d[sc][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Plot assumptions
    #plt.plot(tv,dS1[nam],'b-')
    # plt.plot(tv,1-yDelta,'b-')

    # Assumed change in pellet input
    nam='RemovedDeadStemToPelletMill'
    y0_Pellet=d['BaseCase'][reg][nam][0].copy()
    y1_Pellet=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pellet,y1_Pellet,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pellet
    yD_Pellet=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Fractions of fibre sent to LumberMill, plywood, OSB
    fLumberMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2

    # Add to LumberMill input
    nam='RemovedDeadStemToLumberMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fLumberMill*(yD_Pulp+yD_Pellet)

    # Add to Plywood input
    nam='RemovedDeadStemToPlywoodMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)

    # Add to OSB input
    nam='RemovedDeadStemToOSBMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fOSBMill*(yD_Pulp+yD_Pellet)

    #--------------------------------------------------------------------------
    # Non-merch (same as BaseCase)
    #--------------------------------------------------------------------------

    # No change

#%% Scenario 2

# Transition from pulp to pellets as secondary

sc='S2'

# Fractional reduction in direct transfers from forest ecosystem to pulp mill
# and pellet mill
fChange=0.5

# Phase in period
tStart=2012
tEnd=2026

for reg in regL:

    #--------------------------------------------------------------------------
    # Merch
    #--------------------------------------------------------------------------

    # Change in Pulp input
    nam='RemovedMerchToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Change in pellet input
    nam='RemovedMerchToPelletMill'
    d[sc][reg][nam]=d[sc][reg][nam]+yD_Pulp

    #--------------------------------------------------------------------------
    # Deads
    #--------------------------------------------------------------------------

    # Assumed change in Pulp input
    nam='RemovedDeadStemToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d[sc][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Change in pellet input
    nam='RemovedDeadStemToPelletMill'
    d[sc][reg][nam]=d[sc][reg][nam]+yD_Pulp

    #--------------------------------------------------------------------------
    # Non-merch (same as BaseCase)
    #--------------------------------------------------------------------------

    # Assumed change in Pulp input
    nam='RemovedNonMerchToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d[sc][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Change in pellet input
    nam='RemovedNonMerchToPelletMill'
    d[sc][reg][nam]=d[sc][reg][nam]+yD_Pulp

#%% Save
gu.opickle(meta['Paths']['DB']['Harvest'] + '\\Variables_RemovedFate.pkl',d)

#%%