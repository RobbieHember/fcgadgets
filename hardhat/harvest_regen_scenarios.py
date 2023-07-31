'''
HARVEST REGENERATION - SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from fcgadgets.macgyver import util_general as gu

# Graphics parameters
gp=gu.SetGraphics('Manuscript')

#%% Initilization

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

scnL=['BAU','Actual']

regL=['Coast',
      'Interior',
      'GFS22',
      'Energy Production',
      'Burn Uneconomic',
      'Leave Uneconomic Standing',
      'Sawlogs and Pellets',
      'Energy Production Pellets']

d={}
d['BAU']={}
d['Actual']={}

for scn in scnL:
    for reg in regL:
        d[scn][reg]={}

#%% Business as usual

scn='BAU'

for reg in regL:

    # Fraction planted

    d[scn][reg]['Frac Planted']=np.zeros(tv.size)

    f0=0.0; f0_yr=1940
    f1=0.0; f1_yr=1970

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Planted'][ind]=np.linspace(f0,f1,ind.size)
    d[scn][reg]['Frac Planted'][ind[-1]:]=f1

    # Planting Density

    d[scn][reg]['Planting Density']=np.zeros(tv.size)

    f0=1400; f0_yr=1850
    f1=1400; f1_yr=1970

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Planting Density'][ind]=np.linspace(f0,f1,ind.size)
    d[scn][reg]['Planting Density'][ind[-1]:]=f1

    # Fraction Class A

    d[scn][reg]['Frac Class A']=np.zeros(tv.size)

    f0=0.0; f0_yr=1980
    f1=0.0; f1_yr=2020

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Class A'][ind]=np.linspace(f0,f1,ind.size)

    f0=0.0; f0_yr=2020
    f1=0.0; f1_yr=2050

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Class A'][ind]=np.linspace(f0,f1,ind.size)

    d[scn][reg]['Frac Class A'][ind[-1]:]=f1

    # Average genetic gain

    d[scn][reg]['Ave Genetic Gain']=np.zeros(tv.size)

    f0=0.0; f0_yr=1980
    f1=0.0; f1_yr=2020

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Ave Genetic Gain'][ind]=np.linspace(f0,f1,ind.size)

    f0=0.0; f0_yr=2020
    f1=0.0; f1_yr=2050

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Ave Genetic Gain'][ind]=np.linspace(f0,f1,ind.size)

    d[scn][reg]['Ave Genetic Gain'][ind[-1]:]=f1

#%% Actual

scn='Actual'

for reg in regL:

    # Fraction planted

    d[scn][reg]['Frac Planted']=np.zeros(tv.size)

    f0=0.0; f0_yr=1940
    f1=0.8; f1_yr=2020

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Planted'][ind]=np.linspace(f0,f1,ind.size)

    f0=0.8; f0_yr=2020
    f1=1.0; f1_yr=2050

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Planted'][ind]=np.linspace(f0,f1,ind.size)

    d[scn][reg]['Frac Planted'][ind[-1]:]=f1

    # Planting Density

    d[scn][reg]['Planting Density']=np.zeros(tv.size)

    if reg=='Coast':
        f0=1200; f0_yr=1850
        f1=1200; f1_yr=1970
    else:
        f0=1500; f0_yr=1850
        f1=1500; f1_yr=1970
    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Planting Density'][ind]=np.linspace(f0,f1,ind.size)
    d[scn][reg]['Planting Density'][ind[-1]:]=f1

    # Fractin Class A

    d[scn][reg]['Frac Class A']=np.zeros(tv.size)

    f0=0.0; f0_yr=1980
    f1=0.75; f1_yr=2020

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Class A'][ind]=np.linspace(f0,f1,ind.size)

    f0=0.75; f0_yr=2020
    f1=1.00; f1_yr=2050

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Frac Class A'][ind]=np.linspace(f0,f1,ind.size)

    d[scn][reg]['Frac Class A'][ind[-1]:]=f1

    # Average genetic gain

    d[scn][reg]['Ave Genetic Gain']=np.zeros(tv.size)

    f0=0.0; f0_yr=1980
    f1=20.0; f1_yr=2020

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Ave Genetic Gain'][ind]=np.linspace(f0,f1,ind.size)

    f0=20.0; f0_yr=2020
    f1=25.0; f1_yr=2050

    iT1=np.where(tv==f0_yr)[0]
    iT2=np.where(tv==f1_yr)[0]
    ind=np.arange(iT1,iT2,1,dtype=int)
    d[scn][reg]['Ave Genetic Gain'][ind]=np.linspace(f0,f1,ind.size)

    d[scn][reg]['Ave Genetic Gain'][ind[-1]:]=f1

#%% Plot

iT1=np.where(tv==2020)[0]

plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,10)); aw=0.24; iw=0.08; ms=3; lw=1.25; cla=[0,0,0]; cl1=[0.6,0.9,0]; cl2=[0.27,0.47,0.79]; cl3=[0,0.5,0];
reg='Interior'
ax[0,0].plot([tv[iT1],tv[iT1]],[0,1],'k--',lw=lw,color=[0.7,0.7,0.7])
ax[0,0].plot(tv,d['Actual'][reg]['Frac Planted'],'k-',label='Actual',lw=lw,color=cl1)
ax[0,0].plot(tv,d['BAU'][reg]['Frac Planted'],'k--',label='BAU',lw=lw,color=cl2)
ax[0,0].set(ylabel='Proportion planted (%)',xlim=[1900,2100],xlabel='Time, years')
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

ax[0,1].plot([tv[iT1],tv[iT1]],[0,2000],'k--',lw=lw,color=[0.7,0.7,0.7])
ax[0,1].plot(tv,d['Actual']['Coast']['Planting Density'],'k:',label='Actual (Coast)',lw=lw,color=cl3)
ax[0,1].plot(tv,d['Actual']['Interior']['Planting Density'],'k-',label='Actual (Interior)',lw=lw,color=cl1)
#ax[0,1].plot(tv,d['BAU'][reg]['Planting Density'],'k--',label='BAU (Interior)',lw=lw,color=cl2)
ax[0,1].set(ylabel='Planting density (SPH)',xlim=[1900,2100],xlabel='Time, years')
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)
ax[0,1].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

ax[1,0].plot([tv[iT1],tv[iT1]],[0,1],'k--',lw=lw,color=[0.7,0.7,0.7])
ax[1,0].plot(tv,d['Actual'][reg]['Frac Class A'],'k-',label='Actual',lw=lw,color=cl1)
ax[1,0].plot(tv,d['BAU'][reg]['Frac Class A'],'k--',label='BAU',lw=lw,color=cl2)
ax[1,0].set(ylabel='Proportion Class A seed (%)',xlim=[1900,2100],xlabel='Time, years')
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)
#ax[1,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

ax[1,1].plot([tv[iT1],tv[iT1]],[0,25],'k--',lw=lw,color=[0.7,0.7,0.7])
ax[1,1].plot(tv,d['Actual'][reg]['Ave Genetic Gain'],'k-',label='Actual',lw=lw,color=cl1)
ax[1,1].plot(tv,d['BAU'][reg]['Ave Genetic Gain'],'k--',label='BAU',lw=lw,color=cl2)
ax[1,1].set(ylabel='Ave genetic gain',xlim=[1900,2100],xlabel='Time, years')
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)
#ax[1,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
plt.tight_layout()
gu.axletters(ax,plt,0.0175,0.9,FontColor=cla,Labels=['Proportion planted','Planting Density','Proportion Class A seed','Average genetic gain'],LabelSpacer=0.065)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Regeneration Scenarios\Harvest Regen Scenarios','png',900)

