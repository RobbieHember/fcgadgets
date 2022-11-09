'''
FATE OF FELLED MATERIALS - SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Graphics parameters

fs=6
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import BaseCase HWP parameters

dP0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_FelledFate.xlsx')

#%% Initialize all regions and all scenarios with the BaseCase values

regL=['Coast','Interior','GFS22','Energy Production','BurnUneconomic','LeaveUneconomicStanding','UseUnecomicForPellets','Energy Production Pellets']

dP={}
for reg in regL:
    dP[reg]={}
    for i in range(dP0['Name'].size):
        dP[reg][ dP0['Name'][i] ]=dP0[reg][i]
del dP0

# *** If you change the time vector, change it consistently for other variable scenarios - it needs to be consistent ***
tv=np.arange(1850,2101,1)

fl=np.ones(tv.size)

#%% Initialize database

vnL=['BiomassMerch_Removed','BiomassMerch_Burned','BiomassMerch_LeftOnSite','BiomassMerch_Piled', \
     'BiomassMerch_ToSnag','BiomassNonMerch_Removed','BiomassNonMerch_Burned','BiomassNonMerch_LeftOnSite', \
     'BiomassNonMerch_Piled','BiomassNonMerch_ToSnag','Snags_Removed','Snags_Burned','Snags_LeftOnSite', \
     'Snags_Piled','Snags_ToSnags']

d={}
d['Year']=tv
d['BaseCase']={}
d['S1']={}
d['S2']={}
for reg in regL:
    d['BaseCase'][reg]={}
    d['S1'][reg]={}
    d['S2'][reg]={}

#%% Populate BaseCase scenario

sc='BaseCase'

for rg in regL:
    d[sc][rg]['BiomassMerch_Removed']=fl*dP[rg]['BiomassMerch_Removed']
    d[sc][rg]['BiomassMerch_Piled']=fl*(1-dP[rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BiomassMerch_LeftOnSite']=fl*(1-dP[rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BiomassNonMerch_Removed']=fl*dP[rg]['BiomassNonMerch_Removed']
    d[sc][rg]['BiomassNonMerch_Piled']=fl*(1-dP[rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionNonMerchPiled']
    d[sc][rg]['BiomassNonMerch_LeftOnSite']=fl*(1-dP[rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionNonMerchLeftOnSite']

    d[sc][rg]['Snags_Removed']=fl*dP[rg]['Snags_Removed']
    d[sc][rg]['Snags_Piled']=fl*(1-dP[rg]['Snags_Removed'])*dP[rg]['SpecificFractionSnagsPiled']
    d[sc][rg]['Snags_LeftOnSite']=fl*(1-dP[rg]['Snags_Removed'])*dP[rg]['SpecificFractionSnagsLeftOnSite']

#%% Populate S1 scenario

sc='S1'

for rg in regL:

    # Initialize
    d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])

    iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
    iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

    d[sc][rg]['BiomassMerch_Removed'][iT1]=np.linspace(dP[rg]['BiomassMerch_Removed'],dP[rg]['S1_BiomassMerch_Removed'],iT1.size)
    d[sc][rg]['BiomassMerch_Removed'][iT2]=dP[rg]['S1_BiomassMerch_Removed']
    d[sc][rg]['BiomassMerch_Piled']=(1-d[sc][rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BiomassMerch_LeftOnSite']=(1-d[sc][rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BiomassNonMerch_Removed'][iT1]=np.linspace(dP[rg]['BiomassNonMerch_Removed'],dP[rg]['S1_BiomassNonMerch_Removed'],iT1.size)
    d[sc][rg]['BiomassNonMerch_Removed'][iT2]=dP[rg]['S1_BiomassNonMerch_Removed']
    d[sc][rg]['BiomassNonMerch_Piled']=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionNonMerchPiled']
    d[sc][rg]['BiomassNonMerch_LeftOnSite']=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionNonMerchLeftOnSite']

    d[sc][rg]['Snags_Removed'][iT1]=np.linspace(dP[rg]['Snags_Removed'],dP[rg]['S1_Snags_Removed'],iT1.size)
    d[sc][rg]['Snags_Removed'][iT2]=dP[rg]['S1_Snags_Removed']
    d[sc][rg]['Snags_Piled']=(1-d[sc][rg]['Snags_Removed'])*dP[rg]['SpecificFractionSnagsPiled']
    d[sc][rg]['Snags_LeftOnSite']=(1-d[sc][rg]['Snags_Removed'])*dP[rg]['SpecificFractionSnagsLeftOnSite']

#%% Populate S2 scenario
# Left on site fraction does not change, piled is converted to removals

sc='S2'
for rg in regL:

    d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])
    iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
    iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

    # Merch
    d[sc][rg]['BiomassMerch_LeftOnSite']=(1-d[sc][rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['BiomassMerch_Removed'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['BiomassMerch_Piled']=fPiled0
    d[sc][rg]['BiomassMerch_Piled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['BiomassMerch_Piled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['BiomassMerch_Piled']-d[sc][rg]['BiomassMerch_Piled'][0])
    d[sc][rg]['BiomassMerch_Removed']=d[sc][rg]['BiomassMerch_Removed']+fAddToRem

    # NonMerch
    d[sc][rg]['BiomassNonMerch_LeftOnSite']=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['BiomassNonMerch_Piled']=fPiled0
    d[sc][rg]['BiomassNonMerch_Piled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['BiomassNonMerch_Piled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['BiomassNonMerch_Piled']-d[sc][rg]['BiomassNonMerch_Piled'][0])
    d[sc][rg]['BiomassNonMerch_Removed']=d[sc][rg]['BiomassNonMerch_Removed']+fAddToRem

    # Snags
    d[sc][rg]['Snags_LeftOnSite']=(1-d[sc][rg]['Snags_Removed'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['Snags_Removed'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['Snags_Piled']=fPiled0
    d[sc][rg]['Snags_Piled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['Snags_Piled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['Snags_Piled']-d[sc][rg]['Snags_Piled'][0])
    d[sc][rg]['Snags_Removed']=d[sc][rg]['Snags_Removed']+fAddToRem

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_FelledFate.pkl',d)

#%% Plot

reg='Interior'
#reg='GFS22'

scn='S2'

plt.close('all'); fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(15,10)); aw=0.26; ms=3; lw=1.25; cla=[0,0,0]; cl1=[0.27,0.47,0.79]; cl2=[0.6,0.9,0];
ax[0,0].plot(tv,d['BaseCase'][reg]['BiomassMerch_Removed'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[0,0].plot(tv,d[scn][reg]['BiomassMerch_Removed'],'g--',label=scn,lw=lw,color=cl2)
ax[0,0].set(position=[0.065,0.69,aw,0.29],ylabel='Removal (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

ax[0,1].plot(tv,d['BaseCase'][reg]['BiomassMerch_Piled'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[0,1].plot(tv,d[scn][reg]['BiomassMerch_Piled'],'g--',label=scn,lw=lw,color=cl2)
ax[0,1].set(position=[0.395,0.69,aw,0.29],ylabel='Piled (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)

ax[0,2].plot(tv,d['BaseCase'][reg]['BiomassMerch_LeftOnSite'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[0,2].plot(tv,d[scn][reg]['BiomassMerch_LeftOnSite'],'g--',label=scn,lw=lw,color=cl2)
ax[0,2].set(position=[0.73,0.69,aw,0.29],ylabel='Left on site (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=1.5)

ax[1,0].plot(tv,d['BaseCase'][reg]['BiomassNonMerch_Removed'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[1,0].plot(tv,d[scn][reg]['BiomassNonMerch_Removed'],'g--',label=scn,lw=lw,color=cl2)
ax[1,0].set(position=[0.065,0.385,aw,0.29],ylabel='Removal (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)

ax[1,1].plot(tv,d['BaseCase'][reg]['BiomassNonMerch_Piled'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[1,1].plot(tv,d[scn][reg]['BiomassNonMerch_Piled'],'g--',label=scn,lw=lw,color=cl2)
ax[1,1].set(position=[0.395,0.385,aw,0.29],ylabel='Piled (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)

ax[1,2].plot(tv,d['BaseCase'][reg]['BiomassNonMerch_LeftOnSite'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[1,2].plot(tv,d[scn][reg]['BiomassNonMerch_LeftOnSite'],'g--',label=scn,lw=lw,color=cl2)
ax[1,2].set(position=[0.73,0.385,aw,0.29],ylabel='Left on site (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xticklabels=[''],xlabel='')
ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=1.5)

ax[2,0].plot(tv,d['BaseCase'][reg]['Snags_Removed'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[2,0].plot(tv,d[scn][reg]['Snags_Removed'],'g--',label=scn,lw=lw,color=cl2)
ax[2,0].set(position=[0.065,0.08,aw,0.29],ylabel='Removal (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xlabel='Time, years')
ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=1.5)

ax[2,1].plot(tv,d['BaseCase'][reg]['Snags_Piled'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[2,1].plot(tv,d[scn][reg]['Snags_Piled'],'g--',label=scn,lw=lw,color=cl2)
ax[2,1].set(position=[0.395,0.08,aw,0.29],ylabel='Piled (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xlabel='Time, years')
ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=1.5)

ax[2,2].plot(tv,d['BaseCase'][reg]['Snags_LeftOnSite'],'b-',label='BaseCase',lw=lw,color=cl1)
ax[2,2].plot(tv,d[scn][reg]['Snags_LeftOnSite'],'g--',label=scn,lw=lw,color=cl2)
ax[2,2].set(position=[0.73,0.08,aw,0.29],ylabel='Left on site (%)',xlim=[2000,2100],xticks=np.arange(1500,2200,20),xlabel='Time, years')
ax[2,2].yaxis.set_ticks_position('both'); ax[2,2].xaxis.set_ticks_position('both'); ax[2,2].tick_params(length=1.5)

gu.axletters(ax,plt,0.04,0.9,FontColor=cla,Labels=['Merchantable','Merchantable','Merchantable','Residuals','Residuals','Residuals','Snags','Snags','Snags'],LabelSpacer=0.075)

gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Felled Fate Scenarios\Felled Fate Scenarios ' + reg,'png',900)



