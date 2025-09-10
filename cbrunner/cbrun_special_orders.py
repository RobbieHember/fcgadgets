'''
SPECIAL ORDERS
'''

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_general as gu
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def CT_Johnstone2002(meta,pNam,vi,iScn,iT,iGY,NetGrowth):
	if meta[pNam]['Project']['Code Project']=='Demo_Harv_ThinDensePine':
		GN_Total=NetGrowth[:,iGY['StemMerch']]+NetGrowth[:,iGY['StemNonMerch']]
		if (vi['tv'][iT]<1952):
			# Alter the ratio of merch to non-merch
			NetGrowth[:,iGY['StemMerch']]=0.3*GN_Total
			NetGrowth[:,iGY['StemNonMerch']]=0.7*GN_Total
		if (vi['tv'][iT]>1952) & (vi['tv'][iT]<=1997):
			rBK=NetGrowth[:,iGY['Bark']]/NetGrowth[:,iGY['StemMerch']]
			rBR=NetGrowth[:,iGY['Branch']]/NetGrowth[:,iGY['StemMerch']]
			rF=NetGrowth[:,iGY['Foliage']]/NetGrowth[:,iGY['StemMerch']]
			if iScn==0:
				GG_StemMerch_SO=1.3
				GG_StemNonMerch_SO=0.29
				M_StemMerch_SO=0.5
				M_StemNonMerch_SO=0.88
				NetGrowth[:,iGY['StemMerch']]=GG_StemMerch_SO-M_StemMerch_SO
				NetGrowth[:,iGY['StemNonMerch']]=GG_StemNonMerch_SO-M_StemNonMerch_SO
				NetGrowth[:,iGY['Bark']]=rBK*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Branch']]=rBR*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Foliage']]=rF*NetGrowth[:,iGY['StemMerch']]
			elif iScn==1:
				GG_StemMerch_SO=1.38
				GG_StemNonMerch_SO=0.3
				M_StemMerch_SO=0.15
				M_StemNonMerch_SO=0.45
				NetGrowth[:,iGY['StemMerch']]=GG_StemMerch_SO-M_StemMerch_SO
				NetGrowth[:,iGY['StemNonMerch']]=GG_StemNonMerch_SO-M_StemNonMerch_SO
				NetGrowth[:,iGY['Bark']]=rBK*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Branch']]=rBR*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Foliage']]=rF*NetGrowth[:,iGY['StemMerch']]
	return NetGrowth