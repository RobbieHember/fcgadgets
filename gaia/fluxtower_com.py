#%% Import modules
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.gaia.gaia_util as uga
import fcgadgets.gaia.treefrog_utils as utf
import fcexplore.surface_climate.surfclim_util as usc
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

meta['Graphics']['Print Figures']='On'
meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Flux Towers'
meta['Paths']['Data']=r'C:\Data\Eddy Covariance'
con=uga.HydroMetCon() # Constants
meta=utf.CompileSiteSpecs(meta) # Compile site conditions

#%% EC data
dEC=usc.GetEC(meta)
dECCC=usc.Import10amTemperatureFromECCC()

#%%
import statsmodels.formula.api as smf

y=np.zeros(12)
for iM in range(12):
	d1=dEC['TS']['JP-H02']
	d2={}
	d2['QL']=d1['Heat flux latent'][iM::12]
	d2['D'],mu,sig=gu.zscore(d1['VPD'][iM::12])
	d2['Rn'],mu,sig=gu.zscore(d1['Radiation net'][iM::12])
	df=pd.DataFrame(data=d2);
	df=df.reset_index(drop=True)
	df.dropna()
	md=smf.ols('QL~Rn',df);
	mdf=md.fit(); mdf.summary()
	y[iM]=mdf.params['Rn']
plt.close('all'); plt.plot(y,'-bo')



#%% Temperature correction
utf.DiurnalTemperatureCorrection(meta)

#%% Bowen ratio
meta=utf.Calc_BowenRatio(meta,dEC)

#%% Energy balance closure correction
meta,dEC=utf.Calc_EnergyBalanceClosureCorrection(meta,dEC)

#%% Plot predicted aerodynamic conductance as function of H and U
utf.Plot_Ga_ModelSurface(meta)

#%% Calculate climatologies of aerodynamic conductance
utf.Calc_Ga(meta,dEC)

#%% Latent heat flux

#site='00'
#site='SB-Old'
#site='Mix-F77'
site='Mix-F98'
#site='JP-H75'
#site='JP-Old'

# Aerodynamic conductance
#ind=np.where(meta['Param']['Parameters_ConductanceAerodynamic_ByLC']['LC Class']==lc)[0]
#Ga=meta['Param']['Parameters_ConductanceAerodynamic_ByLC']['Value'][ind]
#Ga=0.2
#Ga=0.058
Ga=0.02

# Canopy conductance
Gs={}
Gs['DF-H49']=np.array([9,5.5,10,7,5.75,4.8,4.6,3.5,4.25,5.5,6,9])/1000
Gs['DF-H88']=np.array([9,5.5,7.5,5.9,4.5,4.5,3.4,2.25,3.65,5.85,7,9])/1000
Gs['DF-H00']=np.array([0.25,0.62,1.05,1.1,1.25,1.5,1.1,0.65,0.75,0.53,0.35,0.2])/1000
Gs['SB-Old']=np.array([0.5,0.6,0.6,1.0,2.0,2.6,2.5,2.7,2.25,2.15,0.75,0.2])/1000
Gs['Mix-F77']=np.array([10,12,5,2.75,5.25,7.5,8.5,12,15,8,12,12])/1000
Gs['Mix-F98']=np.array([10,6,3,2.15,2.5,5.3,6.1,3.25,7.5,5.25,7,12])/1000
Gs['JP-Old']=np.array([0.75,2.0,1.25,1.1,1.05,1.35,1.95,1.8,3.0,2.1,2.3,2.3])/1000
Gs['JP-H75']=np.array([0.75,2.0,1.25,1.1,1.85,2.55,2.1,2.25,2.9,2.1,2.3,2.3])/1000
flg=0
if flg==1:
	plt.close('all');
	plt.plot(Gs['DF-H49'],'-go');
	#plt.plot(Gs['88'],'-gd');
	plt.plot(Gs['DF-H00'],'-gs')
	plt.plot(Gs['SB-Old'],'-bs')
	plt.plot(Gs['JP-Old'],'-rv')
	plt.plot(Gs['JP-H75'],'-r*')

# Evapotranspiratiom (mm month-1)
par={'ETp Method':'Penman-Monteith'}
vi={'Month':np.arange(1,13,1),
	'Ga':Ga,
	'Gs':Gs[site],
	'rn':dEC['N'][site]['Radiation net'],
	'tmean':dEC['N'][site]['Ta'],
	'vpd':dEC['N'][site]['VPD']}
ET=uga.GetETp(vi,'Sparse grid',par['ETp Method'],'Month') # mm month-1
ET=ET/(con['DayLength']*30.5) # mm s-1
ET=ET # kg s-1
QL=ET*con['Lam'] # Latent heat flux (W m-2)

# Plot
x=dEC['N'][site]['Heat flux latent']
y=QL
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
ax.plot(x,y,'o',ms=5,mec='w',mfc=[0.27,0.49,0.77],mew=0.75)
ax.plot(x[1],y[1],'ks',ms=5,mec='k',mfc='w',mew=0.75)
rs,txt=gu.GetRegStats(x,y)
ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',lw=0.5)
ax.text(60,5,rs['txt'],fontsize=5,ha='right')
ax.text(92,92,'1:1',fontsize=5,ha='center')
ax.set(xlabel='Flux tower (W m-2)',ylabel='Model (W m-2)',xticks=np.arange(-30,300,10),yticks=np.arange(-30,300,10),xlim=[0,110],ylim=[0,110])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Name','png',900)

#%%
v='VPD'
#v='Ta'
#v='Tsoil'
#v='Radiation net'

plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(15.5,6)); ms=3; lw=0.5
ax[0].plot(dEC['N']['DF-H49'][v],Gs['DF-H49'],'-bo',ms=ms,lw=lw,mfc='w',mec='b',label='DF-H49')
ax[0].plot(dEC['N']['DF-H49'][v][3],Gs['DF-H49'][3],'ks',mfc='w')
ax[0].plot(dEC['N']['DF-H88'][v],Gs['DF-H88'],'-gs',ms=ms,lw=lw,mfc='w',mec='g',label='DF-H88')
ax[0].plot(dEC['N']['DF-H00'][v],Gs['DF-H00'],'-cv',ms=ms,lw=lw,mfc='w',mec='c',label='DF-H00')
ax[0].plot(dEC['N']['DF-H00'][v][3],Gs['DF-H00'][3],'ks',mfc='w')
ax[0].set(ylabel='Water conductance (m s$^{-1}$)',ylim=[0,0.011])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
ax[0].legend(loc='upper right',facecolor=[1,1,1],frameon=False,fontsize=6)
ax[1].plot(dEC['N']['JP-Old'][v],Gs['JP-Old'],'-bo',ms=ms,lw=lw,mfc='w',label='JP-Old')
ax[1].plot(dEC['N']['JP-Old'][v][3],Gs['JP-Old'][3],'ks',mfc='w')
ax[1].plot(dEC['N']['JP-H75'][v],Gs['JP-H75'],'-ms',ms=ms,lw=lw,mfc='w',label='JP-H75')
ax[1].set(ylabel='Water conductance (m s$^{-1}$)',ylim=[0,0.011])
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
ax[1].legend(loc='upper right',facecolor=[1,1,1],frameon=False,fontsize=6)
ax[2].plot(dEC['N']['SB-Old'][v],Gs['SB-Old'],'-bo',ms=ms,lw=lw,mfc='w',label='SB-Old')
#plt.plot(dEC['N']['SB-Old'][v][3],Gs['SB-Old'][3],'ks')
ax[2].plot(dEC['N']['Mix-F77'][v],Gs['Mix-F77'],'-gs',ms=ms,lw=lw,mfc='w',label='Mix-F77')
ax[2].plot(dEC['N']['Mix-F77'][v][3],Gs['Mix-F77'][3],'ks',mfc='w')
ax[2].plot(dEC['N']['Mix-F98'][v],Gs['Mix-F98'],'-rv',ms=ms,lw=lw,mfc='w',label='Mix-F98')
ax[2].set(ylabel='Water conductance (m s$^{-1}$)',ylim=[0,0.016])
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])
ax[2].legend(loc='upper right',facecolor=[1,1,1],frameon=False,fontsize=6)
plt.tight_layout()
gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Conductance vs VPD','png',900)


#%% Relationship between GC and slope of Ts-NDVI relationship
x=np.array([-25,25])
y=np.array([0.0005,0.015])
rs,txt=gu.GetRegStats(x,y) # 0.00028*x+0.008
meta['f(Gc)']=rs['B']

#%% Explore

v='albedo'
v='Heat flux ground'
v='Heat flux latent'
v='U'
iM=7
y=np.array([np.nanmean(dEC['TS']['00'][v][iM::12]),
			np.nanmean(dEC['TS']['88'][v][iM::12]),
			np.nanmean(dEC['TS']['49'][v][iM::12])])
plt.close('all')
plt.plot(y,'-bo')


#%% Summary of radiation balance
d={}
for v in dEC['N']['DF-H49'].keys():
	d[v]=np.zeros(10)
	d['Site']=np.array([])
	for i,s in enumerate(dEC['N'].keys()):
		d[v][i]=np.round(np.mean(dEC['N'][s][v]),decimals=2)
		d['Site']=np.append(d['Site'],s)
df=pd.DataFrame.from_dict(d)
df=df.set_index('Site')
df.T.to_excel(meta['Paths']['Data'] + '\\Summary.xlsx')

#%%
v='Radiation longwave down'
plt.close('all')
plt.plot(dEC['N']['00'][v],'-ro')
plt.plot(dEC['N']['88'][v],'-go')
plt.plot(dEC['N']['49'][v],'-bo')

#%%
v='Radiation longwave up'
plt.close('all')
plt.plot(dEC['N']['00'][v],'-ro')
plt.plot(dEC['N']['88'][v],'-go')
plt.plot(dEC['N']['49'][v],'-bo')

#%%
v='Radiation longwave net'
plt.close('all')
plt.plot(dEC['N']['00'][v],'-ro')
plt.plot(dEC['N']['88'][v],'-go')
plt.plot(dEC['N']['49'][v],'-bo')

#%% Plot monthly deltas (impact of forest loss)

plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(15.5,7)); ms=3; lw=0.5; x=np.arange(1,13,1); xlim=[0.5,12.5]
v='albedo';
ax[0,0].plot(x,dEC['N']['00'][v]-dEC['N']['49'][v],'-bo',ms=ms,lw=lw)
ax[0,0].set(ylabel='Albedo',xticks=x,xlim=xlim)
v='Radiation shortwave net';
rswn=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[0,1].plot(x,rswn,'-bo',ms=ms,lw=lw)
ax[0,1].set(ylabel='R shortwave net (W m-2)',xticks=x,xlim=xlim)
v='Radiation longwave net';
rlwn=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[0,2].plot(x,rlwn,'-bo',ms=ms,lw=lw)
ax[0,2].set(ylabel='R longwave net (W m-2)',xticks=x,xlim=xlim)
v='Radiation net';
rn=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[1,0].plot(x,rn,'-bo',ms=ms,lw=lw)
#ax[1,0].plot(x,rswn+rlwn,'--cs',ms=ms,lw=lw,mfc='w')
ax[1,0].set(ylabel='R net (W m-2)',xticks=x,xlim=xlim)
v='Heat flux sensible'#v='Ta'
qh=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[1,1].plot(x,qh,'-bo',ms=ms,lw=lw)
v='Heat flux ground'#v='Ta'
qh=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[1,1].plot(x,qh,'--cs',ms=ms,lw=lw,mfc='w')
ax[1,1].set(ylabel='QH (W m-2)',xticks=x,xlim=xlim)
#v='T_CNR1'; ax[1,1].plot(dEC['N']['00'][v]-dEC['N']['49'][v],'--cs',ms=ms,lw=lw,mfc='w')
#ax[1,1].set(ylabel='Ta (K)',xticks=x,xlim=xlim)
v='Heat flux latent';
ql=dEC['N']['00'][v]-dEC['N']['49'][v]
ax[1,2].plot(x,ql,'-bo',ms=ms,lw=lw)
#ax[1,0].plot(x,rn-ql,'--cs',ms=ms,lw=lw,mfc='w')
ax[1,0].plot(x,rn-ql-qh,'--cs',ms=ms,lw=lw,mfc='w')
ax[1,2].set(ylabel='QL (W m-2)',xticks=x,xlim=xlim)
for i in range(2):
	for j in range(3):
		ax[i,j].plot([0,13],[0,0],'k-',lw=1,color=[0.8,0.8,0.8])
		ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Albedo_ByBGC_ComparisonWithMODIS','png',900)
