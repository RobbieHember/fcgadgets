'''
REMOVED FATE SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Graphics parameters

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

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

regL=['Coast','Interior','GFS22']
scnL=['BaseCase','S1']

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

sc='S1'

# Phase in period
tStart=2020
tEnd=2030
    
for reg in regL:
  
    #--------------------------------------------------------------------------
    # Merch
    #--------------------------------------------------------------------------
    
    # Assumed change in Pulp input
    nam='RemovedMerchToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=0.5*d['BaseCase'][reg][nam][0]
    d['S1'][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d['S1'][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d['S1'][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d['S1'][reg][nam]
    
    # Plot assumptions
    #plt.plot(tv,dS1[nam],'b-')
    # plt.plot(tv,1-yDelta,'b-')
    
    # Assumed change in pellet input
    nam='RemovedMerchToPelletMill'
    y0_Pellet=d['BaseCase'][reg][nam][0].copy()
    y1_Pellet=0.0*d['BaseCase'][reg][nam][0]
    d['S1'][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d['S1'][reg][nam][iT]=np.linspace(y0_Pellet,y1_Pellet,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d['S1'][reg][nam][iT]=y1_Pellet
    yD_Pellet=d['BaseCase'][reg][nam]-d['S1'][reg][nam]
    
    # Fractions of fibre sent to sawmill, plywood, OSB
    fSawMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2
    
    # Add to Sawmill input
    nam='RemovedMerchToSawMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fSawMill*(yD_Pulp+yD_Pellet)
    
    # Add to Plywood input
    nam='RemovedMerchToPlywoodMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)
    
    # Add to OSB input
    nam='RemovedMerchToOSBMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fOSBMill*(yD_Pulp+yD_Pellet)
    
    #--------------------------------------------------------------------------
    # Snags
    #--------------------------------------------------------------------------
    
    # Assumed change in Pulp input
    nam='RemovedSnagStemToPulpMill'
    y0_Pulp=d['BaseCase'][reg][nam][0].copy()
    y1_Pulp=0.25*d['BaseCase'][reg][nam][0]
    d['S1'][reg][nam]=d['S1'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d['S1'][reg][nam][iT]=np.linspace(y0_Pulp,y1_Pulp,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d['S1'][reg][nam][iT]=y1_Pulp
    yD_Pulp=d['BaseCase'][reg][nam]-d['S1'][reg][nam]
    
    # Plot assumptions
    #plt.plot(tv,dS1[nam],'b-')
    # plt.plot(tv,1-yDelta,'b-')
    
    # Assumed change in pellet input
    nam='RemovedSnagStemToPelletMill'
    y0_Pellet=d['BaseCase'][reg][nam][0].copy()
    y1_Pellet=0.0*d['BaseCase'][reg][nam][0]
    d['S1'][reg][nam]=d['BaseCase'][reg][nam].copy()
    iT=np.where( (tv>=tStart) & (tv<=tEnd) )[0]
    d['S1'][reg][nam][iT]=np.linspace(y0_Pellet,y1_Pellet,iT.size)
    iT=np.where( (tv>tEnd) )[0]
    d['S1'][reg][nam][iT]=y1_Pellet
    yD_Pellet=d['BaseCase'][reg][nam]-d['S1'][reg][nam]
    
    # Fractions of fibre sent to sawmill, plywood, OSB
    fSawMill=0.6
    fPlywoodMill=0.2
    fOSBMill=0.2
    
    # Add to Sawmill input
    nam='RemovedSnagStemToSawMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fSawMill*(yD_Pulp+yD_Pellet)
    
    # Add to Plywood input
    nam='RemovedSnagStemToPlywoodMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fPlywoodMill*(yD_Pulp+yD_Pellet)
    
    # Add to OSB input
    nam='RemovedSnagStemToOSBMill'
    d['S1'][reg][nam]=d['S1'][reg][nam]+fOSBMill*(yD_Pulp+yD_Pellet)
    
    #--------------------------------------------------------------------------
    # Non-merch (same as BaseCase)
    #--------------------------------------------------------------------------
    
    # No change


#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_RemovedFate.pkl',d)


#%% Plot

reg='Interior'
reg='GFS22'

plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,10)); ms=3; lw=1.25
ax[0,0].plot(tv,d['BaseCase'][reg]['RemovedMerchToSawMill'],'b-',label='BaseCase',lw=lw)
ax[0,0].plot(tv,d['S1'][reg]['RemovedMerchToSawMill'],'g--',label='Scenario 1',lw=lw)
ax[0,0].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Merch to sawmill (%)',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
ax[0,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

ax[0,1].plot(tv,d['BaseCase'][reg]['RemovedMerchToPulpMill'],'b-',label='BaseCase',lw=lw)
ax[0,1].plot(tv,d['S1'][reg]['RemovedMerchToPulpMill'],'g--',label='Scenario 1',lw=lw)
ax[0,1].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Merch to pulp mill (%)',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[0,1].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)

ax[1,0].plot(tv,d['BaseCase'][reg]['RemovedMerchToPelletMill'],'b-',label='BaseCase',lw=lw)
ax[1,0].plot(tv,d['S1'][reg]['RemovedMerchToPelletMill'],'g--',label='Scenario 1',lw=lw)
ax[1,0].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Merch to pellet mill (%)',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)

ax[1,1].plot(tv,d['BaseCase'][reg]['RemovedMerchToPlywoodMill']+d['BaseCase'][reg]['RemovedMerchToOSBMill'],'b-',label='BaseCase',lw=lw)
ax[1,1].plot(tv,d['S1'][reg]['RemovedMerchToPlywoodMill']+d['S1'][reg]['RemovedMerchToOSBMill'],'g--',label='Scenario 1',lw=lw)
ax[1,1].set(ylim=[0,1],yticks=np.arange(0,1.2,0.2),ylabel='Merch to plywood mill (%)',xlim=[1900,2100],xticks=np.arange(1500,2200,25),xlabel='Time, years')
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)

gu.axletters(ax,plt,0.028,0.9,LetterStyle='Caps',FontWeight='Bold') # 

fig.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Removed Fate Scenarios\Removed Fate Scenarios ' + reg,'png',900)
