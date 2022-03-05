'''
FATE OF FELLED MATERIALS - SCENARIOS

'''

#%% Import modules

import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

#%% Import default HWP parameters

dP0=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_Harvesting.xlsx')
dP={}
for i in range(dP0['Name'].size):
    dP[ dP0['Name'][i] ]=dP0['Value'][i]
del dP0

tv=np.arange(1850,2101,1)

fl=np.ones(tv.size)

rgL=['Coast','Interior']

#%% Initialize database

vnL=['BiomassMerch_Removed','BiomassMerch_Burned','BiomassMerch_LeftOnSite','BiomassMerch_Piled', \
     'BiomassMerch_ToSnag','BiomassNonMerch_Removed','BiomassNonMerch_Burned','BiomassNonMerch_LeftOnSite', \
     'BiomassNonMerch_Piled','BiomassNonMerch_ToSnag','Snags_Removed','Snags_Burned','Snags_LeftOnSite', \
     'Snags_Piled','Snags_ToSnags']

d={}
d['Year']=tv
d['Default']={'Coast':{},'Interior':{}}
d['S1']={'Coast':{},'Interior':{}}

#%% Populate default scenario

sc='Default'

for rg in rgL:
    d[sc][rg]['BiomassMerch_Removed']=fl*dP['BiomassMerch_Removed_' + rg]
    d[sc][rg]['BiomassMerch_Piled']=fl*(1-dP['BiomassMerch_Removed_' + rg])*dP['SpecificFractionMerchPiled_' + rg]
    d[sc][rg]['BiomassMerch_LeftOnSite']=fl*(1-dP['BiomassMerch_Removed_' + rg])*dP['SpecificFractionMerchLeftOnSite_' + rg]
    
    d[sc][rg]['BiomassNonMerch_Removed']=fl*dP['BiomassNonMerch_Removed_' + rg]
    d[sc][rg]['BiomassNonMerch_Piled']=fl*(1-dP['BiomassNonMerch_Removed_' + rg])*dP['SpecificFractionNonMerchPiled_' + rg]
    d[sc][rg]['BiomassNonMerch_LeftOnSite']=fl*(1-dP['BiomassNonMerch_Removed_' + rg])*dP['SpecificFractionNonMerchLeftOnSite_' + rg]
    
    d[sc][rg]['Snags_Removed']=fl*dP['Snags_Removed_' + rg]
    d[sc][rg]['Snags_Piled']=fl*(1-dP['Snags_Removed_' + rg])*dP['SpecificFractionSnagsPiled_' + rg]
    d[sc][rg]['Snags_LeftOnSite']=fl*(1-dP['Snags_Removed_' + rg])*dP['SpecificFractionSnagsLeftOnSite_' + rg]

#%% Populate S1 scenario

sc='S1'

for rg in rgL:
    
    # Initialize
    d[sc][rg]=copy.deepcopy(d['Default'][rg])
    
    iT1=np.where( (tv>=dP['S1_StartYear_' + rg]) & (tv<=dP['S1_TargetYear_' + rg]) )[0]
    iT2=np.where( (tv>dP['S1_TargetYear_' + rg]) )[0]
    
    d[sc][rg]['BiomassMerch_Removed'][iT1]=np.linspace(dP['BiomassMerch_Removed_' + rg],dP['S1_BiomassMerch_Removed_' + rg],iT1.size)
    d[sc][rg]['BiomassMerch_Removed'][iT2]=dP['S1_BiomassMerch_Removed_' + rg]
    d[sc][rg]['BiomassMerch_Piled']=(1-d[sc][rg]['BiomassMerch_Removed'])*dP['SpecificFractionMerchPiled_' + rg]
    d[sc][rg]['BiomassMerch_LeftOnSite']=(1-d[sc][rg]['BiomassMerch_Removed'])*dP['SpecificFractionMerchLeftOnSite_' + rg]
    
    d[sc][rg]['BiomassNonMerch_Removed'][iT1]=np.linspace(dP['BiomassNonMerch_Removed_' + rg],dP['S1_BiomassNonMerch_Removed_' + rg],iT1.size)
    d[sc][rg]['BiomassNonMerch_Removed'][iT2]=dP['S1_BiomassNonMerch_Removed_' + rg]
    d[sc][rg]['BiomassNonMerch_Piled']=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP['SpecificFractionNonMerchPiled_' + rg]
    d[sc][rg]['BiomassNonMerch_LeftOnSite']=(1-d[sc][rg]['BiomassNonMerch_Removed'])*dP['SpecificFractionNonMerchLeftOnSite_' + rg]
    
    d[sc][rg]['Snags_Removed'][iT1]=np.linspace(dP['Snags_Removed_' + rg],dP['S1_Snags_Removed_' + rg],iT1.size)
    d[sc][rg]['Snags_Removed'][iT2]=dP['S1_Snags_Removed_' + rg]
    d[sc][rg]['Snags_Piled']=(1-d[sc][rg]['Snags_Removed'])*dP['SpecificFractionSnagsPiled_' + rg]
    d[sc][rg]['Snags_LeftOnSite']=(1-d[sc][rg]['Snags_Removed'])*dP['SpecificFractionSnagsLeftOnSite_' + rg]

#%% Plot
    
plt.close('all')
plt.plot(tv,d['Default']['Interior']['BiomassMerch_Removed'],'b-')
plt.plot(tv,d['S1']['Interior']['BiomassMerch_Removed'],'g--')

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Variables_FelledFate.pkl',d)

