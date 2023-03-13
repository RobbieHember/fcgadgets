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


