'''
FATE OF FELLED MATERIALS - SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from fcgadgets.macgyver import util_general as gu

#%% Graphics parameters

gp=gu.SetGraphics('Manuscript')

#%% Import BaseCase HWP parameters

dP0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_FelledFate.xlsx')

#%% Initialize all regions and all scenarios with the BaseCase values

regL=['Coast',
      'Interior',
      'GFS22',
      'Energy Production',
      'Burn Uneconomic',
      'Leave Uneconomic Standing',
      'Sawlogs and Pellets',
      'Energy Production Pellets']

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

vnL=['StemwoodMerchRemoved','StemwoodMerchBurned','StemwoodMerchLeftOnSite','StemwoodMerchPiled','StemwoodMerchToSnag', \
     'StemwoodNonMerchRemoved','StemwoodNonMerchBurned','StemwoodNonMerchLeftOnSite','StemwoodNonMerchPiled','StemwoodNonMerchToSnag', \
     'FoliageRemoved','FoliageBurned','FoliageLeftOnSite','FoliagePiled','FoliageToSnag', \
     'StemwoodNonMerchRemoved','BarkBurned','BarkLeftOnSite','BarkPiled','BarkToSnag', \
     'BranchRemoved','BranchBurned','BranchLeftOnSite','BranchPiled','BranchToSnag', \
     'SnagStemRemoved','SnagStemBurned','SnagStemLeftOnSite','SnagStemPiled','SnagStemToSnagStem', \
     'SnagBranchRemoved','SnagBranchBurned','SnagBranchLeftOnSite','SnagBranchPiled','SnagBranchToSnagBranch']

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
    d[sc][rg]['StemwoodMerchRemoved']=fl*dP[rg]['StemwoodMerchRemoved']
    d[sc][rg]['StemwoodMerchPiled']=fl*(1-dP[rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['StemwoodMerchLeftOnSite']=fl*(1-dP[rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['StemwoodNonMerchRemoved']=fl*dP[rg]['StemwoodNonMerchRemoved']
    d[sc][rg]['StemwoodNonMerchPiled']=fl*(1-dP[rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionNonMerchPiled']
    d[sc][rg]['StemwoodNonMerchLeftOnSite']=fl*(1-dP[rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionNonMerchLeftOnSite']

    d[sc][rg]['FoliageRemoved']=fl*dP[rg]['FoliageRemoved']
    d[sc][rg]['FoliagePiled']=fl*(1-dP[rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['FoliageLeftOnSite']=fl*(1-dP[rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BarkRemoved']=fl*dP[rg]['BarkRemoved']
    d[sc][rg]['BarkPiled']=fl*(1-dP[rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BarkLeftOnSite']=fl*(1-dP[rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BranchRemoved']=fl*dP[rg]['BranchRemoved']
    d[sc][rg]['BranchPiled']=fl*(1-dP[rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BranchLeftOnSite']=fl*(1-dP[rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['SnagStemRemoved']=fl*dP[rg]['SnagStemRemoved']
    d[sc][rg]['SnagStemPiled']=fl*(1-dP[rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionSnagStemPiled']
    d[sc][rg]['SnagStemLeftOnSite']=fl*(1-dP[rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionSnagStemLeftOnSite']

    d[sc][rg]['SnagBranchRemoved']=fl*dP[rg]['SnagBranchRemoved']
    d[sc][rg]['SnagBranchPiled']=fl*(1-dP[rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionSnagBranchPiled']
    d[sc][rg]['SnagBranchLeftOnSite']=fl*(1-dP[rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionSnagBranchLeftOnSite']

#%% Populate S1 scenario

sc='S1'

for rg in regL:

    # Initialize
    d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])

    iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
    iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

    d[sc][rg]['StemwoodMerchRemoved'][iT1]=np.linspace(dP[rg]['StemwoodMerchRemoved'],dP[rg]['S1_StemwoodMerchRemoved'],iT1.size)
    d[sc][rg]['StemwoodMerchRemoved'][iT2]=dP[rg]['S1_StemwoodMerchRemoved']
    d[sc][rg]['StemwoodMerchPiled']=(1-d[sc][rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['StemwoodMerchLeftOnSite']=(1-d[sc][rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['StemwoodNonMerchRemoved'][iT1]=np.linspace(dP[rg]['StemwoodNonMerchRemoved'],dP[rg]['S1_StemwoodNonMerchRemoved'],iT1.size)
    d[sc][rg]['StemwoodNonMerchRemoved'][iT2]=dP[rg]['S1_StemwoodNonMerchRemoved']
    d[sc][rg]['StemwoodNonMerchPiled']=(1-d[sc][rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionNonMerchPiled']
    d[sc][rg]['StemwoodNonMerchLeftOnSite']=(1-d[sc][rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionNonMerchLeftOnSite']

    d[sc][rg]['FoliageRemoved'][iT1]=np.linspace(dP[rg]['FoliageRemoved'],dP[rg]['S1_FoliageRemoved'],iT1.size)
    d[sc][rg]['FoliageRemoved'][iT2]=dP[rg]['S1_FoliageRemoved']
    d[sc][rg]['FoliagePiled']=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['FoliageLeftOnSite']=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BarkRemoved'][iT1]=np.linspace(dP[rg]['BarkRemoved'],dP[rg]['S1_BarkRemoved'],iT1.size)
    d[sc][rg]['BarkRemoved'][iT2]=dP[rg]['S1_BarkRemoved']
    d[sc][rg]['BarkPiled']=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BarkLeftOnSite']=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['BranchRemoved'][iT1]=np.linspace(dP[rg]['BranchRemoved'],dP[rg]['S1_BranchRemoved'],iT1.size)
    d[sc][rg]['BranchRemoved'][iT2]=dP[rg]['S1_BranchRemoved']
    d[sc][rg]['BranchPiled']=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    d[sc][rg]['BranchLeftOnSite']=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']

    d[sc][rg]['SnagStemRemoved'][iT1]=np.linspace(dP[rg]['SnagStemRemoved'],dP[rg]['S1_SnagStemRemoved'],iT1.size)
    d[sc][rg]['SnagStemRemoved'][iT2]=dP[rg]['S1_SnagStemRemoved']
    d[sc][rg]['SnagStemPiled']=(1-d[sc][rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionSnagStemPiled']
    d[sc][rg]['SnagStemLeftOnSite']=(1-d[sc][rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionSnagStemLeftOnSite']

    d[sc][rg]['SnagBranchRemoved'][iT1]=np.linspace(dP[rg]['SnagBranchRemoved'],dP[rg]['S1_SnagBranchRemoved'],iT1.size)
    d[sc][rg]['SnagBranchRemoved'][iT2]=dP[rg]['S1_SnagBranchRemoved']
    d[sc][rg]['SnagBranchPiled']=(1-d[sc][rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionSnagBranchPiled']
    d[sc][rg]['SnagBranchLeftOnSite']=(1-d[sc][rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionSnagBranchLeftOnSite']

#%% Populate S2 scenario
# Left on site fraction does not change, piled is converted to removals

sc='S2'
for rg in regL:

    d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])
    iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
    iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

    # Merch
    d[sc][rg]['StemwoodMerchLeftOnSite']=(1-d[sc][rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['StemwoodMerchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['StemwoodMerchPiled']=fPiled0
    d[sc][rg]['StemwoodMerchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['StemwoodMerchPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['StemwoodMerchPiled']-d[sc][rg]['StemwoodMerchPiled'][0])
    d[sc][rg]['StemwoodMerchRemoved']=d[sc][rg]['StemwoodMerchRemoved']+fAddToRem

    # NonMerch
    d[sc][rg]['StemwoodNonMerchLeftOnSite']=(1-d[sc][rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['StemwoodNonMerchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['StemwoodNonMerchPiled']=fPiled0
    d[sc][rg]['StemwoodNonMerchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['StemwoodNonMerchPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['StemwoodNonMerchPiled']-d[sc][rg]['StemwoodNonMerchPiled'][0])
    d[sc][rg]['StemwoodNonMerchRemoved']=d[sc][rg]['StemwoodNonMerchRemoved']+fAddToRem

    # Foliage
    d[sc][rg]['FoliageLeftOnSite']=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['FoliagePiled']=fPiled0
    d[sc][rg]['FoliagePiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['FoliagePiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['FoliagePiled']-d[sc][rg]['FoliagePiled'][0])
    d[sc][rg]['FoliageRemoved']=d[sc][rg]['FoliageRemoved']+fAddToRem

    # Bark
    d[sc][rg]['BarkLeftOnSite']=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['BarkPiled']=fPiled0
    d[sc][rg]['BarkPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['BarkPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['BarkPiled']-d[sc][rg]['BarkPiled'][0])
    d[sc][rg]['BarkRemoved']=d[sc][rg]['BarkRemoved']+fAddToRem

    # Branch
    d[sc][rg]['BranchLeftOnSite']=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['BranchPiled']=fPiled0
    d[sc][rg]['BranchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['BranchPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['BranchPiled']-d[sc][rg]['BranchPiled'][0])
    d[sc][rg]['BranchRemoved']=d[sc][rg]['BranchRemoved']+fAddToRem

    # SnagStem
    d[sc][rg]['SnagStemLeftOnSite']=(1-d[sc][rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['SnagStemRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['SnagStemPiled']=fPiled0
    d[sc][rg]['SnagStemPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['SnagStemPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['SnagStemPiled']-d[sc][rg]['SnagStemPiled'][0])
    d[sc][rg]['SnagStemRemoved']=d[sc][rg]['SnagStemRemoved']+fAddToRem

    # SnagBranch
    d[sc][rg]['SnagBranchLeftOnSite']=(1-d[sc][rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionMerchLeftOnSite']
    fPiled0=(1-d[sc][rg]['SnagBranchRemoved'])*dP[rg]['SpecificFractionMerchPiled']
    fPiled1=0.0
    d[sc][rg]['SnagBranchPiled']=fPiled0
    d[sc][rg]['SnagBranchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
    d[sc][rg]['SnagBranchPiled'][iT2]=fPiled1
    fAddToRem=-1*(d[sc][rg]['SnagBranchPiled']-d[sc][rg]['SnagBranchPiled'][0])
    d[sc][rg]['SnagBranchRemoved']=d[sc][rg]['SnagBranchRemoved']+fAddToRem

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_FelledFate.pkl',d)


