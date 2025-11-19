#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import time
import warnings
import pandas as pd
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun_preprocess as prep
import fcgadgets.cbrunner.cbrun_postprocess as post
import fcgadgets.cbrunner.cbrun as cbr
import fcgadgets.macgyver.util_fcs_graphs as ufcs
import fcgadgets.macgyver.util_demo as udem
import fcgadgets.gaia.gaia_util as gaia
warnings.filterwarnings("ignore")
gp=gu.SetGraphics('Manuscript')

#%% Configure project
meta=u1ha.Init()
pNam='Zartan'
meta['Paths'][pNam]={}
meta['Paths'][pNam]['Data']=r'D:\Modelling Projects\Zartan'
meta=cbu.ImportProjectConfig(meta,pNam)
meta['Graphics']['Print Figures']='On'
meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Demo\Demo_FNM'

#%% Prepare inputs
prep.Write_BatchTIPSY_Input_File(meta,pNam)
prep.PrepareInventoryFromSpreadsheet(meta,pNam)
prep.BuildEventChronologyFromSpreadsheet(meta,pNam)
prep.PrepGrowthCurvesForCBR(meta,pNam)

#%% Import growth curves

txtDat=np.loadtxt(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Output.out',skiprows=4)
df=pd.DataFrame(txtDat,columns=['A','V_G','V_N0','V_N125','V_N175','M'])
Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)
N_Age=Age.size
N_GC=int(df.shape[0]/N_Age)

d0={}
for k in df.columns:
	d0[k]=np.reshape(df[k].values,(N_Age,N_GC),order='F')

fl=np.ones(Age.size)

d1={}
d1['V_G']=np.array([0])
d1['V_N0']=np.array([0])
d1['A']=np.array([0])
d1['SI']=np.array([0])
d1['N']=np.array([0])
d1['RM']=np.array([0])
d1['GG']=np.array([0])
for i in range(d0['A'].shape[1]):
	dP=gu.ReadExcel(meta['Paths'][pNam]['Data'] + '\\Inputs\\GrowthCurvesTIPSY_Parameters.xlsx',sheet_name='Sheet1',skiprows=6)

	d1['A']=np.append(d1['A'],Age)
	d1['V_G']=np.append(d1['V_G'],d0['V_G'][:,i])
	d1['V_N0']=np.append(d1['V_N0'],d0['V_N0'][:,i])

	d1['SI']=np.append(d1['SI'],dP['i1'][i]*fl)
	d1['N']=np.append(d1['N'],dP['init_density'][i]*fl)
	if dP['regeneration_method'][i]=='P':
		RM=1
	else:
		RM=0
	d1['RM']=np.append(d1['RM'],RM*fl)
	d1['GG']=np.append(d1['GG'],np.nan_to_num(dP['gain1'][i])*fl)

for k in d1.keys():
	d1[k]=d1[k][1:]

#%% Remove zero growth in old age
# ikp=np.where( (d['A']>10) & (d['Gswm']>0) )[0]
# for k in d.keys():
# 	d[k]=d[k][ikp]

#%% Chapman Richards equation (Zeide 1993)
def funCR0(A,**b0):
	B1=b0['b1']
	B2=b0['b2']
	B3=b0['b3']
	yhat=B1*(1-np.exp(-B2*A))**B3
	yhat=np.nan_to_num(yhat)
	return yhat

def funCR0s(A,b1):
	B1=b1['b1']
	B2=b1['b2']
	B3=b1['b3']
	yhat=B1*(1-np.exp(-B2*A))**B3
	yhat=np.nan_to_num(yhat)
	return yhat

#%% Stats modelling
def FitYield(meta,d1):

	uSI=np.unique(d1['SI'])
	uN=np.unique(d1['N'])
	uRM=np.unique(d1['RM'])
	uGG=np.unique(d1['GG'])

	rs={}
	for iSI,si in enumerate(uSI):
		rs[si]={}
		for iN,n in enumerate(uN):
			rs[si][n]={}
			for iRM,rm in enumerate(uRM):
				rs[si][n][rm]={}
				for iGG,gg in enumerate(uGG):
					rs[si][n][rm][gg]={}

					# Filter
					ikp=np.where( (d1['SI']==si) & (d1['N']==n) & (d1['RM']==rm) & (d1['GG']==gg) )[0]
					d2={}
					for k in d1.keys():
						d2[k]=d1[k][ikp]
				
					# Initial parameters
					b0={}
					b0['b1']=1200
					b0['b2']=0.05
					b0['b3']=10
				
					# Fit model
					fitmodel=lmfit.Model(funCR0,independent_vars=['A'])
					params=fitmodel.make_params()
					for k,v in b0.items():
						params.add(k,value=v)
					rs0=fitmodel.fit(d2['V_G'],params,A=d2['A'])
				
					# Output dictionary
					rs[si][n][rm][gg]['R2']=rs0.rsquared
					rs[si][n][rm][gg]['BE']={}
					for k in rs0.params.keys():
						rs[si][n][rm][gg]['BE'][k]=np.round(rs0.params[k].value,decimals=4)
					rs[si][n][rm][gg]['TIPSY']=d2['V_G']

	return rs

rs=FitYield(meta,d1)

plt.close('all')
plt.plot(Age,funCR0s(Age,rs[25][1000][0][0]['BE']),'b-')
plt.plot(Age,funCR0s(Age,rs[30][1000][0][0]['BE']),'g-')
plt.plot(Age,funCR0s(Age,rs[30][1000][1][0]['BE']),'g--')
plt.plot(Age,funCR0s(Age,rs[30][1000][0][10]['BE']),'r:')

plt.plot(Age,rs[25][1000][0][0]['TIPSY'],'b--')
plt.plot(Age,rs[30][1000][0][0]['TIPSY'],'g--')

#%%