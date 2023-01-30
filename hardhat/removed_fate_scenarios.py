'''
REMOVED FATE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Graphics parameters

gp=gu.SetGraphics('Manuscript')

#%% Import BaseCase HWP parameters

dDef0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_RemovedFate.xlsx')

# Variables to keep
vrL=['RemovedMerchToChipperMill','RemovedMerchToPulpMill','RemovedMerchToPelletMill','RemovedMerchToSawMill',
    'RemovedMerchToPlywoodMill','RemovedMerchToOSBMill','RemovedMerchToMDFMill','RemovedMerchToPolePostMill',
    'RemovedMerchToShakeShingleMill','RemovedMerchToLogExport','RemovedMerchToIPP','RemovedMerchToFirewood',
    'RemovedNonMerchToChipperMill','RemovedNonMerchToPulpMill','RemovedNonMerchToPelletMill','RemovedNonMerchToSawMill','RemovedNonMerchToPlywoodMill',
    'RemovedNonMerchToOSBMill','RemovedNonMerchToMDFMill','RemovedNonMerchToPolePostMill','RemovedNonMerchToShakeShingleMill',
    'RemovedNonMerchToLogExport','RemovedNonMerchToIPP','RemovedNonMerchToFirewood',
    'RemovedSnagStemToChipperMill','RemovedSnagStemToPulpMill',
    'RemovedSnagStemToPelletMill','RemovedSnagStemToSawMill','RemovedSnagStemToPlywoodMill','RemovedSnagStemToOSBMill',
    'RemovedSnagStemToMDFMill','RemovedSnagStemToPolePostMill','RemovedSnagStemToShakeShingleMill','RemovedSnagStemToLogExport',
    'RemovedSnagStemToIPP','RemovedSnagStemToFirewood']

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

# Region list
regL=['Coast',
      'Interior',
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

# Shift away from pulp+pellets to sawmills
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

    # Fractions of fibre sent to sawmill, plywood, OSB
    fSawMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2

    # Add to Sawmill input
    nam='RemovedMerchToSawMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fSawMill*(yD_Pulp+yD_Pellet)

    # Add to Plywood input
    nam='RemovedMerchToPlywoodMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)

    # Add to OSB input
    nam='RemovedMerchToOSBMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fOSBMill*(yD_Pulp+yD_Pellet)

    #--------------------------------------------------------------------------
    # Snags
    #--------------------------------------------------------------------------

    # Assumed change in Pulp input
    nam='RemovedSnagStemToPulpMill'
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
    nam='RemovedSnagStemToPelletMill'
    y0_Pellet=d['BaseCase'][reg][nam][0].copy()
    y1_Pellet=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pellet,y1_Pellet,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pellet
    yD_Pellet=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Fractions of fibre sent to sawmill, plywood, OSB
    fSawMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2

    # Add to Sawmill input
    nam='RemovedSnagStemToSawMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fSawMill*(yD_Pulp+yD_Pellet)

    # Add to Plywood input
    nam='RemovedSnagStemToPlywoodMill'
    d[sc][reg][nam]=d[sc][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)

    # Add to OSB input
    nam='RemovedSnagStemToOSBMill'
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
    # Snags
    #--------------------------------------------------------------------------

    # Assumed change in Pulp input
    nam='RemovedSnagStemToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=fChange*d['BaseCase'][reg][nam][0]
    d[sc][reg][nam]=d[sc][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d[sc][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d[sc][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d[sc][reg][nam]

    # Change in pellet input
    nam='RemovedSnagStemToPelletMill'
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

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_RemovedFate.pkl',d)


#%% Plot

flg=0
if flg==1:
    reg='Interior'
    #reg='GFS22'
    sc='S2'

    plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,10)); ms=3; lw=1.25
    ax[0,0].plot(tv,d['BaseCase'][reg]['RemovedMerchToSawMill'],'b-',label='BaseCase',lw=lw)
    ax[0,0].plot(tv,d[sc][reg]['RemovedMerchToSawMill'],'g--',label='Scenario 2',lw=lw)
    ax[0,0].set(yticks=np.arange(0,1.2,0.2),ylabel='Merch to sawmill (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
    ax[0,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

    ax[0,1].plot(tv,d['BaseCase'][reg]['RemovedMerchToPulpMill'],'b-',label='BaseCase',lw=lw)
    ax[0,1].plot(tv,d[sc][reg]['RemovedMerchToPulpMill'],'g--',label='Scenario 2',lw=lw)
    ax[0,1].set(yticks=np.arange(0,1.2,0.2),ylabel='Merch to pulp mill (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)

    ax[1,0].plot(tv,d['BaseCase'][reg]['RemovedMerchToPelletMill'],'b-',label='BaseCase',lw=lw)
    ax[1,0].plot(tv,d[sc][reg]['RemovedMerchToPelletMill'],'g--',label='Scenario 2',lw=lw)
    ax[1,0].set(yticks=np.arange(0,1.2,0.2),ylabel='Merch to pellet mill (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)

    ax[1,1].plot(tv,d['BaseCase'][reg]['RemovedMerchToPlywoodMill']+d['BaseCase'][reg]['RemovedMerchToOSBMill'],'b-',label='BaseCase',lw=lw)
    ax[1,1].plot(tv,d[sc][reg]['RemovedMerchToPlywoodMill']+d[sc][reg]['RemovedMerchToOSBMill'],'g--',label='Scenario 2',lw=lw)
    ax[1,1].set(yticks=np.arange(0,1.2,0.2),ylabel='Merch to plywood mill (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)

    gu.axletters(ax,plt,0.028,0.9,LetterStyle='Caps',FontWeight='Bold') #

    fig.tight_layout()
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Removed Fate Scenarios\Removed Fate Scenarios ' + reg + 'Scenario ' + str(sc),'png',900)
