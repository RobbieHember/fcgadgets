'''
CO2 EXPERIMENTS SUMMARY
'''

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy.optimize import curve_fit
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu

#%% Import data
dC=gu.ReadExcel(r'C:\Data\CO2_Experiments\co2_experiments_R20200329.xlsx')

# Just keep summary list
ind=np.where(dC['Summary']=='Yes')[0]
for k in dC.keys():
	dC[k]=dC[k][ind]
n=dC['CO2_C'].size

#%% Functions
def LogFunc(x,b1,b2):
	return b1*(1+b2*np.log(x/280))

#%% Calculate standardized response
rs={}
rs['co2_ref']=280
rs['co2']=np.arange(rs['co2_ref'],1000,1)
rs['b_opt']=np.nan*np.ones((n,2))
rs['yhat']=np.nan*np.ones((rs['co2'].size,n))
rs['yhat_ex']=np.nan*np.ones((n,2))
for i in range(n):
	x=np.array([dC['CO2_C'][i],dC['CO2_T'][i]])
	y=np.array([1,dC['RR'][i]])
	try:
		p,pcov=curve_fit(LogFunc,x,y)
		rs['b_opt'][i,:]=p
		rs['yhat'][:,i]=LogFunc(rs['co2'],p[0],p[1])
		rs['yhat_ex'][i,:]=LogFunc(x,p[0],p[1])
	except:
		pass

iRef=np.where(rs['co2']==rs['co2_ref'])[0]
rs['yhat_z']=rs['yhat']/np.tile(rs['yhat'][iRef,:],(rs['co2'].size,1))

# Remove outliers
co2_today=400
iToday=np.where(rs['co2']==co2_today)[0]
rs['rr400']=rs['yhat_z'][iToday[0],:]
th_Outlier=2*np.nanstd(rs['rr400'])
dC['Outlier']=np.zeros(n)
iOutlier=np.where(rs['rr400']>th_Outlier)[0]
dC['Outlier'][iOutlier]=1
rs['yhat_z_or']=rs['yhat_z'].copy()
rs['rr400_or']=rs['rr400'].copy()
rs['yhat_z_or'][:,iOutlier]=np.nan
rs['rr400_or'][iOutlier]=np.nan
print(str(iOutlier.size) + ' experiments excluded because they were deemed to be outliers.')

#%% CO2 response function from field plots
rsG=gu.ipickle(meta['Paths']['Model']['Parameters'] + '\\Parameters_gromo_gg2.pkl')

y_FP=np.zeros( (rs['co2'].size,len(meta['Modules']['gromo']['Species CD'])) )
for i,s in enumerate(meta['Modules']['gromo']['Species CD']):
	A=75;Sd=0;Hi=0;Dab=0;Dad=0;Daf=0;Dap=0;Ta_z=0;Wa_z=0;Nd_z=0;Ca_z=0
	sc={};sc['PL']=0;sc['SW']=0;sc['BL']=0;sc['FDC']=0;sc['FDI']=0;sc['HWC']=0;sc['HWI']=0;sc['CW']=0;sc['SE']=0;sc['AT']=0;sc['SB']=0;sc['SS']=0;sc['AE']=0;sc['O']=0;
	sc[s]=1
	Tn_z=rsG['BySpecies'][s]['Mean']['Tn_z']
	Wn_z=rsG['BySpecies'][s]['Mean']['Wn_z']
	Nd_z=rsG['BySpecies'][s]['Mean']['Nd_z']
	
	x=np.array([280,380])
	yhat=np.array([0,0],dtype='float')
	Ca_z=(x[0]-rsG['Zscore Stats']['Ca']['mu'])/rsG['Zscore Stats']['Ca']['sig']
	yhat[0]=ugm.PredictGG2(A,sc['PL'],sc['SW'],sc['BL'],sc['FDC'],sc['FDI'],sc['HWC'],sc['HWI'],sc['CW'],sc['SE'],sc['AT'],sc['SB'],sc['SS'],sc['AE'],sc['O'],Sd,Hi,Dab,Dad,Daf,Dap,Tn_z,Wn_z,Ta_z,Wa_z,Nd_z,Ca_z,rsG['Param']['BE'])
	Ca_z=(x[1]-rsG['Zscore Stats']['Ca']['mu'])/rsG['Zscore Stats']['Ca']['sig']
	yhat[1]=ugm.PredictGG2(A,sc['PL'],sc['SW'],sc['BL'],sc['FDC'],sc['FDI'],sc['HWC'],sc['HWI'],sc['CW'],sc['SE'],sc['AT'],sc['SB'],sc['SS'],sc['AE'],sc['O'],Sd,Hi,Dab,Dad,Daf,Dap,Tn_z,Wn_z,Ta_z,Wa_z,Nd_z,Ca_z,rsG['Param']['BE'])
	rr=yhat[1]/yhat[0]
	y=np.array([1,rr])
	p,pcov=curve_fit(LogFunc,x,y)
	y_FP[:,i]=LogFunc(rs['co2'],p[0],p[1])

#%% Summary plot
def Plot_CO2_Response_Standard(rs):
	N_all=np.sum(~np.isnan(rs['yhat_z_or'][0,:]))
	mu_all=np.nanmean(rs['yhat_z_or'],axis=1)
	sd_all=np.nanstd(rs['yhat_z_or'],axis=1)
	sig2_all=2*sd_all/np.sqrt(N_all)
	med_all=np.nanmedian(rs['yhat_z_or'],axis=1)
	mu_Seedling=np.nanmean(rs['yhat_z_or'][:,dC['Source']=='Wullschleger et al. 1995'],axis=1)
	mu_FieldGrown=np.nanmean(rs['yhat_z_or'][:,dC['Source']=='Norby et al. 1999'],axis=1)
	mu_Plot=np.nanmean(rs['yhat_z_or'][:,dC['Experimental Setup']=='FACE'],axis=1)
	
	lw1=1
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
	ax.fill_between(rs['co2'],mu_all-sd_all,mu_all+sd_all,color=[1,0.87,0.83],alpha=0.5,linewidth=0)
	ax.fill_between(rs['co2'],mu_all-sig2_all,mu_all+sig2_all,color=[1,0.74,0.63],alpha=0.25,linewidth=0)
	ax.plot(rs['co2'],mu_all,'k-',color=[0.65,0,0],linewidth=lw1,label='All mean')
	#ax.plot(rs['co2'],med_all,'k--',color=[0.65,0,0],linewidth=lw1,label='All median')
	
	ax.plot(rs['co2'],mu_Seedling,'k:',color=[0.5,0,1],linewidth=lw1,label='Seedlings (Wullschleger et al. 1995)')
	ax.plot(rs['co2'],mu_FieldGrown,'k--',color=[0.5,0,1],linewidth=lw1,label='Field-grown (Norby et al. 1999)')
	ax.plot(rs['co2'],mu_Plot,'k-',color=[0.5,0,1],linewidth=lw1,label='FACE mean')

	# Add field plots
	mu=np.mean(y_FP,axis=1)
	sd=np.std(y_FP,axis=1)
	mx=np.max(y_FP,axis=1)
	mn=np.min(y_FP,axis=1)
	ax.fill_between(rs['co2'],mn,mx,color=[0.8,1,0.5],alpha=0.5,linewidth=0)
	#ax.fill_between(rs['co2'],mu-sd,mu+sd,color=[0.8,1,0.5],alpha=0.5,linewidth=0)
	ax.plot(rs['co2'],mu,'k-',color=[0.5,0.9,0],linewidth=lw1,label='Field plots')
	
	ax.plot([275,1050],[1,1],'k-',linewidth=0.5)
	ax.set(xlabel='Atmospheric CO$_2$ concentration (ppm)',ylabel='Relative response',ylim=[0.5,2],xlim=[275,550]);
	ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25,fontsize=6)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\ResponseCO2_Comparison','png',900)
	return
Plot_CO2_Response_Standard(rs)
#%%

