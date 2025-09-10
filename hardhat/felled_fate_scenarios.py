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
dP0=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\Parameters_FelledFate.xlsx')

#%% Initialize all regions and all scenarios with the BaseCase values
regL=['Coast',
      'Interior',
      'Full Utilization',
      'CT2',
      'CT3',
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

vnL=['GreenStemMerchRemoved','GreenStemMerchBurned','GreenStemMerchLeftOnSite','GreenStemMerchPiled','GreenStemMerchToSnag', \
     'GreenStemNonMerchRemoved','GreenStemNonMerchBurned','GreenStemNonMerchLeftOnSite','GreenStemNonMerchPiled','GreenStemNonMerchToSnag', \
     'FoliageRemoved','FoliageBurned','FoliageLeftOnSite','FoliagePiled','FoliageToSnag', \
     'BarkRemoved','BarkBurned','BarkLeftOnSite','BarkPiled','BarkToSnag', \
     'BranchRemoved','BranchBurned','BranchLeftOnSite','BranchPiled','BranchToSnag', \
     'DeadStemMerchRemoved','DeadStemMerchBurned','DeadStemMerchLeftOnSite','DeadStemMerchPiled','DeadStemMerchToDeadStemMerch', \
     'DeadStemNonMerchRemoved','DeadStemNonMerchBurned','DeadStemNonMerchLeftOnSite','DeadStemNonMerchPiled','DeadStemNonMerchToDeadStemNonMerch', \
     'DeadBranchRemoved','DeadBranchBurned','DeadBranchLeftOnSite','DeadBranchPiled','DeadBranchToDeadBranch']

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
	d[sc][rg]['GreenStemMerchRemoved']=fl*dP[rg]['GreenStemMerchRemoved']
	d[sc][rg]['GreenStemMerchPiled']=fl*dP[rg]['FractionGreenStemMerchPiled']
	d[sc][rg]['GreenStemMerchLeftOnSite']=fl*dP[rg]['FractionGreenStemMerchLeftOnSite']
	#d[sc][rg]['GreenStemMerchPiled']=fl*(1-dP[rg]['GreenStemMerchRemoved'])*dP[rg]['FractionGreenStemMerchPiled']
	#d[sc][rg]['GreenStemMerchLeftOnSite']=fl*(1-dP[rg]['GreenStemMerchRemoved'])*dP[rg]['FractionGreenStemMerchLeftOnSite']

	d[sc][rg]['GreenStemNonMerchRemoved']=fl*dP[rg]['GreenStemNonMerchRemoved']
	d[sc][rg]['GreenStemNonMerchPiled']=fl*dP[rg]['FractionGreenStemNonMerchPiled']
	d[sc][rg]['GreenStemNonMerchLeftOnSite']=fl*dP[rg]['FractionGreenStemNonMerchLeftOnSite']

	d[sc][rg]['FoliageRemoved']=fl*dP[rg]['FoliageRemoved']
	d[sc][rg]['FoliagePiled']=fl*dP[rg]['FractionFoliagePiled']
	d[sc][rg]['FoliageLeftOnSite']=fl*dP[rg]['FractionFoliageLeftOnSite']

	d[sc][rg]['BarkRemoved']=fl*dP[rg]['BarkRemoved']
	d[sc][rg]['BarkPiled']=fl*dP[rg]['FractionBarkPiled']
	d[sc][rg]['BarkLeftOnSite']=fl*dP[rg]['FractionBarkLeftOnSite']

	d[sc][rg]['BranchRemoved']=fl*dP[rg]['BranchRemoved']
	d[sc][rg]['BranchPiled']=fl*dP[rg]['FractionBranchPiled']
	d[sc][rg]['BranchLeftOnSite']=fl*dP[rg]['FractionBranchLeftOnSite']

	d[sc][rg]['DeadStemMerchRemoved']=fl*dP[rg]['DeadStemMerchRemoved']
	d[sc][rg]['DeadStemMerchPiled']=fl*dP[rg]['FractionDeadStemMerchPiled']
	d[sc][rg]['DeadStemMerchLeftOnSite']=fl*dP[rg]['FractionDeadStemMerchLeftOnSite']

	d[sc][rg]['DeadStemNonMerchRemoved']=fl*dP[rg]['DeadStemNonMerchRemoved']
	d[sc][rg]['DeadStemNonMerchPiled']=fl*dP[rg]['FractionDeadStemNonMerchPiled']
	d[sc][rg]['DeadStemNonMerchLeftOnSite']=fl*dP[rg]['FractionDeadStemNonMerchLeftOnSite']

	d[sc][rg]['DeadBranchRemoved']=fl*dP[rg]['DeadBranchRemoved']
	d[sc][rg]['DeadBranchPiled']=fl*dP[rg]['FractionDeadBranchPiled']
	d[sc][rg]['DeadBranchLeftOnSite']=fl*dP[rg]['FractionDeadBranchLeftOnSite']

	d[sc][rg]['DeadBarkRemoved']=fl*dP[rg]['DeadBarkRemoved']
	d[sc][rg]['DeadBarkPiled']=fl*dP[rg]['FractionDeadBarkPiled']
	d[sc][rg]['DeadBarkLeftOnSite']=fl*dP[rg]['FractionDeadBarkLeftOnSite']

#%% Populate S1 scenario

sc='S1'

for rg in regL:

	# Initialize
	d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])

	iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
	iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

	d[sc][rg]['GreenStemMerchRemoved'][iT1]=np.linspace(dP[rg]['GreenStemMerchRemoved'],dP[rg]['S1_GreenStemMerchRemoved'],iT1.size)
	d[sc][rg]['GreenStemMerchRemoved'][iT2]=dP[rg]['S1_GreenStemMerchRemoved']
	d[sc][rg]['GreenStemMerchPiled']=dP[rg]['FractionGreenStemMerchPiled']
	d[sc][rg]['GreenStemMerchLeftOnSite']=dP[rg]['FractionGreenStemMerchLeftOnSite']

	d[sc][rg]['GreenStemNonMerchRemoved'][iT1]=np.linspace(dP[rg]['GreenStemNonMerchRemoved'],dP[rg]['S1_GreenStemNonMerchRemoved'],iT1.size)
	d[sc][rg]['GreenStemNonMerchRemoved'][iT2]=dP[rg]['S1_GreenStemNonMerchRemoved']
	d[sc][rg]['GreenStemNonMerchPiled']=dP[rg]['FractionGreenStemNonMerchPiled']
	d[sc][rg]['GreenStemNonMerchLeftOnSite']=dP[rg]['FractionGreenStemNonMerchLeftOnSite']

	d[sc][rg]['FoliageRemoved'][iT1]=np.linspace(dP[rg]['FoliageRemoved'],dP[rg]['S1_FoliageRemoved'],iT1.size)
	d[sc][rg]['FoliageRemoved'][iT2]=dP[rg]['S1_FoliageRemoved']
	d[sc][rg]['FoliagePiled']=dP[rg]['FractionFoliagePiled']
	d[sc][rg]['FoliageLeftOnSite']=dP[rg]['FractionFoliageLeftOnSite']

	d[sc][rg]['BarkRemoved'][iT1]=np.linspace(dP[rg]['BarkRemoved'],dP[rg]['S1_BarkRemoved'],iT1.size)
	d[sc][rg]['BarkRemoved'][iT2]=dP[rg]['S1_BarkRemoved']
	d[sc][rg]['BarkPiled']=dP[rg]['FractionBarkPiled']
	d[sc][rg]['BarkLeftOnSite']=dP[rg]['FractionBarkLeftOnSite']

	d[sc][rg]['BranchRemoved'][iT1]=np.linspace(dP[rg]['BranchRemoved'],dP[rg]['S1_BranchRemoved'],iT1.size)
	d[sc][rg]['BranchRemoved'][iT2]=dP[rg]['S1_BranchRemoved']
	d[sc][rg]['BranchPiled']=dP[rg]['FractionBranchPiled']
	d[sc][rg]['BranchLeftOnSite']=dP[rg]['FractionBranchLeftOnSite']

	d[sc][rg]['DeadStemMerchRemoved'][iT1]=np.linspace(dP[rg]['DeadStemMerchRemoved'],dP[rg]['S1_DeadStemMerchRemoved'],iT1.size)
	d[sc][rg]['DeadStemMerchRemoved'][iT2]=dP[rg]['S1_DeadStemMerchRemoved']
	d[sc][rg]['DeadStemMerchPiled']=dP[rg]['FractionDeadStemMerchPiled']
	d[sc][rg]['DeadStemMerchLeftOnSite']=dP[rg]['FractionDeadStemMerchLeftOnSite']

	d[sc][rg]['DeadStemNonMerchRemoved'][iT1]=np.linspace(dP[rg]['DeadStemNonMerchRemoved'],dP[rg]['S1_DeadStemNonMerchRemoved'],iT1.size)
	d[sc][rg]['DeadStemNonMerchRemoved'][iT2]=dP[rg]['S1_DeadStemNonMerchRemoved']
	d[sc][rg]['DeadStemNonMerchPiled']=dP[rg]['FractionDeadStemNonMerchPiled']
	d[sc][rg]['DeadStemNonMerchLeftOnSite']=dP[rg]['FractionDeadStemNonMerchLeftOnSite']

	d[sc][rg]['DeadBranchRemoved'][iT1]=np.linspace(dP[rg]['DeadBranchRemoved'],dP[rg]['S1_DeadBranchRemoved'],iT1.size)
	d[sc][rg]['DeadBranchRemoved'][iT2]=dP[rg]['S1_DeadBranchRemoved']
	d[sc][rg]['DeadBranchPiled']=dP[rg]['FractionDeadBranchPiled']
	d[sc][rg]['DeadBranchLeftOnSite']=dP[rg]['FractionDeadBranchLeftOnSite']

	d[sc][rg]['DeadBarkRemoved'][iT1]=np.linspace(dP[rg]['DeadBarkRemoved'],dP[rg]['S1_DeadBarkRemoved'],iT1.size)
	d[sc][rg]['DeadBarkRemoved'][iT2]=dP[rg]['S1_DeadBarkRemoved']
	d[sc][rg]['DeadBarkPiled']=dP[rg]['FractionDeadBarkPiled']
	d[sc][rg]['DeadBarkLeftOnSite']=dP[rg]['FractionDeadBarkLeftOnSite']

#%% Populate S2 scenario
# Left on site fraction does not change, piled is converted to removals

sc='S2'
for rg in regL:

	d[sc][rg]=copy.deepcopy(d['BaseCase'][rg])
	iT1=np.where( (tv>=dP[rg]['S1_StartYear']) & (tv<=dP[rg]['S1_TargetYear']) )[0]
	iT2=np.where( (tv>dP[rg]['S1_TargetYear']) )[0]

	# Merch
	d[sc][rg]['GreenStemMerchLeftOnSite']=(1-d[sc][rg]['GreenStemMerchRemoved'])*dP[rg]['FractionGreenStemMerchLeftOnSite']
	fPiled0=(1-d[sc][rg]['GreenStemMerchRemoved'])*dP[rg]['FractionGreenStemMerchPiled']
	fPiled1=0.0
	d[sc][rg]['GreenStemMerchPiled']=fPiled0
	d[sc][rg]['GreenStemMerchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['GreenStemMerchPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['GreenStemMerchPiled']-d[sc][rg]['GreenStemMerchPiled'][0])
	d[sc][rg]['GreenStemMerchRemoved']=d[sc][rg]['GreenStemMerchRemoved']+fAddToRem

	# NonMerch
	d[sc][rg]['GreenStemNonMerchLeftOnSite']=(1-d[sc][rg]['GreenStemNonMerchRemoved'])*dP[rg]['FractionGreenStemNonMerchLeftOnSite']
	fPiled0=(1-d[sc][rg]['GreenStemNonMerchRemoved'])*dP[rg]['FractionGreenStemNonMerchPiled']
	fPiled1=0.0
	d[sc][rg]['GreenStemNonMerchPiled']=fPiled0
	d[sc][rg]['GreenStemNonMerchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['GreenStemNonMerchPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['GreenStemNonMerchPiled']-d[sc][rg]['GreenStemNonMerchPiled'][0])
	d[sc][rg]['GreenStemNonMerchRemoved']=d[sc][rg]['GreenStemNonMerchRemoved']+fAddToRem

	# Foliage
	d[sc][rg]['FoliageLeftOnSite']=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['FractionFoliageLeftOnSite']
	fPiled0=(1-d[sc][rg]['FoliageRemoved'])*dP[rg]['FractionFoliagePiled']
	fPiled1=0.0
	d[sc][rg]['FoliagePiled']=fPiled0
	d[sc][rg]['FoliagePiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['FoliagePiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['FoliagePiled']-d[sc][rg]['FoliagePiled'][0])
	d[sc][rg]['FoliageRemoved']=d[sc][rg]['FoliageRemoved']+fAddToRem

	# Bark
	d[sc][rg]['BarkLeftOnSite']=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['FractionBarkLeftOnSite']
	fPiled0=(1-d[sc][rg]['BarkRemoved'])*dP[rg]['FractionBarkPiled']
	fPiled1=0.0
	d[sc][rg]['BarkPiled']=fPiled0
	d[sc][rg]['BarkPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['BarkPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['BarkPiled']-d[sc][rg]['BarkPiled'][0])
	d[sc][rg]['BarkRemoved']=d[sc][rg]['BarkRemoved']+fAddToRem

	# Branch
	d[sc][rg]['BranchLeftOnSite']=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['FractionBranchLeftOnSite']
	fPiled0=(1-d[sc][rg]['BranchRemoved'])*dP[rg]['FractionBranchPiled']
	fPiled1=0.0
	d[sc][rg]['BranchPiled']=fPiled0
	d[sc][rg]['BranchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['BranchPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['BranchPiled']-d[sc][rg]['BranchPiled'][0])
	d[sc][rg]['BranchRemoved']=d[sc][rg]['BranchRemoved']+fAddToRem

	# DeadStemMerch
	d[sc][rg]['DeadStemMerchLeftOnSite']=(1-d[sc][rg]['DeadStemMerchRemoved'])*dP[rg]['FractionDeadStemMerchLeftOnSite']
	fPiled0=(1-d[sc][rg]['DeadStemMerchRemoved'])*dP[rg]['FractionDeadStemMerchPiled']
	fPiled1=0.0
	d[sc][rg]['DeadStemMerchPiled']=fPiled0
	d[sc][rg]['DeadStemMerchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['DeadStemMerchPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['DeadStemMerchPiled']-d[sc][rg]['DeadStemMerchPiled'][0])
	d[sc][rg]['DeadStemMerchRemoved']=d[sc][rg]['DeadStemMerchRemoved']+fAddToRem

	# DeadBranch
	d[sc][rg]['DeadBranchLeftOnSite']=(1-d[sc][rg]['DeadBranchRemoved'])*dP[rg]['FractionDeadBranchLeftOnSite']
	fPiled0=(1-d[sc][rg]['DeadBranchRemoved'])*dP[rg]['FractionDeadBranchPiled']
	fPiled1=0.0
	d[sc][rg]['DeadBranchPiled']=fPiled0
	d[sc][rg]['DeadBranchPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['DeadBranchPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['DeadBranchPiled']-d[sc][rg]['DeadBranchPiled'][0])
	d[sc][rg]['DeadBranchRemoved']=d[sc][rg]['DeadBranchRemoved']+fAddToRem

	# DeadBark
	d[sc][rg]['DeadBarkLeftOnSite']=(1-d[sc][rg]['DeadBarkRemoved'])*dP[rg]['FractionDeadBarkLeftOnSite']
	fPiled0=(1-d[sc][rg]['DeadBarkRemoved'])*dP[rg]['FractionDeadBarkPiled']
	fPiled1=0.0
	d[sc][rg]['DeadBarkPiled']=fPiled0
	d[sc][rg]['DeadBarkPiled'][iT1]=np.linspace(fPiled0[iT1[0]],fPiled1,iT1.size)
	d[sc][rg]['DeadBarkPiled'][iT2]=fPiled1
	fAddToRem=-1*(d[sc][rg]['DeadBarkPiled']-d[sc][rg]['DeadBarkPiled'][0])
	d[sc][rg]['DeadBarkRemoved']=d[sc][rg]['DeadBarkRemoved']+fAddToRem

#%% Save
gu.opickle(meta['Paths']['DB']['Harvest'] + '\\Variables_FelledFate.pkl',d)

#%%


