'''
HARVEST SUMMARY BY TIMBER MARK - ANALYSIS
'''
#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fcgadgets.hardhat.harvest_by_tm_util as hut
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_util as u1ha
import statsmodels.formula.api as smf

#%% Import parameters
meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')
gradeL=list(meta['LUT']['Derived']['LogGrade'].keys())

#%% Import data
flg=0
if flg==1:
	# Import HBS by TM and year
	dTMY=hut.HarvestByTMY(meta)

	# Get area from forest tenure layer
	dTMY,dOP=hut.GetAreaFromFTEN(meta,dTMY)

	# Create TM DB from the TMY DB (sum all years wtihin TMs)
	dTM=hut.Harvest_ByTM_FromTMY(meta,dTMY)

	# Import Waste System data
	dTMY,dTM=hut.ImportWaste(meta,dTMY,dTM)

	# Add timber cruise data
	dTMY,dTM=hut.ImportTimberCruise(meta,dTMY,dTM)

	# Import RESULTS OP layer
	dTM=hut.ImportRESULTS(meta,dOP,dTM)

	# Create derived variables
	dTMY,dTM=hut.CreateDerivedVariables(meta,dTMY,dTM)

	# Add raster variables
	dTMY,dTM=hut.GetVariablesFromRasterDB(meta,dTMY,dTM)

	# Save
	gu.opickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTMAndYear.pkl',dTMY)
	gu.opickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTM.pkl',dTM)

	# Export to spreadsheet for review
	#df=pd.DataFrame(dTM)
	#df.to_excel(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTM.xlsx',index=False)
else:
	# Import saved data
	dTMY=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTMAndYear.pkl')
	dTM=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTM.pkl')

#%% QA 15% of TMs in HBS are missing area from FTEN layer (private land?)
ind1=np.where(dTMY['PLANNED_NET_BLOCK_AREA']==0)[0]
ind2=np.where(dTMY['V Logs m3']>0)[0]
ind1.size/ind2.size

#%% QA Relationship between net volumes from different sources
hut.NetVolumeComparison(meta,dTM)

#%% Plot time series of harvest volume
hut.PlotHarvestVolume_TS(meta,dTMY)

#%% Piling rate
hut.CalcPilingRateByBGCZone(meta,dTM,save='On',plot='On')

#%% Regression Waste Rate
hut.WasteRateRegression(meta,dTM,save='On',print_table='On',plot='On')

#%% Look at frequence of tenure type
u,N=gu.CountByCategories(dTM['Tenure Type'],'Percent')
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
ax.barh(np.arange(u.size),N)
ax.set(yticks=np.arange(u.size),yticklabels=u)

#%% Look at frequency of SSC
u,N=gu.CountByCategories(dTM['SSC'],'Percent')
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,14));
ax.barh(np.arange(u.size),N)
ax.set(yticks=np.arange(u.size),yticklabels=u)

#%% QA - Look at stumpage by grade
for gr in gradeL:
	dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/(453*dTM['V Logs Grade ' + gr + ' Abs m3']/1000)
	dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/dTM['V Logs Grade ' + gr + ' Abs m3']
ind=np.where(dTM['Region']=='Coast')[0]
d={}
for k in dTM.keys():
	try:
		ind=np.where( (dTM['Region']=='Coast') & (dTM[k]>0) & (dTM[k]<1000000) )[0]
		d[k]=np.nanmean(dTM[k][ind])
	except:
		pass

#%%
d=hut.StatsByDistrict(meta,dTM)

#%%
d=hut.StatsByBGCZone(meta,dTM)

#%%
d=hut.CalcHarvestVolumeByBGCZ(meta,dTMY,dTM)

#%%
hut.StatsBySILV_SYSTEM_CODE(meta,dTM)

#%%
zone='MS'
hut.StatsBySILV_SYSTEM_CODE_ForBGCZone(meta,dTM,zone)

#%% Relationship between Waste wood and net volume
list(dTM.keys())
ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['WR']>=0) & (dTM['WR']<2) )[0]
x=dTM['V Felled m3/ha'][ind]
y=dTM['Waste Total m3/ha'][ind]
#rs,txt=gu.GetRegStats(x,y,'No Intercept')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
#ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
#ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
#ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(ylabel='Waste wood (m3/ha)',xlabel='Net volume (m3/ha)',xlim=[0,1500],ylim=[0,300])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Relationship between waste wood from HBS and WS
list(dTM.keys())
ind=np.where( (dTM['V Logs Waste m3/ha']>0) & (dTM['V Logs Waste m3/ha']<400) & (dTM['Waste Total m3/ha']>=0) & (dTM['Waste Total m3/ha']<400) )[0]
x=dTM['V Logs Waste m3/ha'][ind]
y=dTM['Waste Total m3/ha'][ind]
rs,txt=gu.GetRegStats(x,y,'No Intercept')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
ax.text(300,300,'1:1',fontsize=10,ha='center',va='center')
ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
ax.text(385,60,txt,fontsize=8,ha='right',va='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(ylabel='Waste wood (m3/ha)',xlabel='Net volume (m3/ha)',xlim=[0,400],ylim=[0,400])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Relationship between WR and net volume
list(dTM.keys())
ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['WR']>=0) & (dTM['WR']<2) )[0]
x=dTM['V Logs Abs m3/ha'][ind]
y=dTM['WR'][ind]
#rs,txt=gu.GetRegStats(x,y,'No Intercept')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
#ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
#ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
#ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,1500],ylim=[0,2])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Relationship between WR and age
list(dTM.keys())
ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['WR']>=0) & (dTM['WR']<2) )[0]
x=dTM['age_vri02'][ind]
y=dTM['WR'][ind]
#rs,txt=gu.GetRegStats(x,y,'No Intercept')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
#ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
#ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
#ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,300],ylim=[0,2])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% WR vs. % dead (analysis of salvage)
ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']<2500) )[0]
x=dTM['Cruise_Pct Dead Net'][ind]
bw=20; bin=np.arange(0,120,bw)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
y=dTM['WR'][ind]
#y=dTM['Sawlog Ratio'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'-bo',ms=5,mfc=[0.27,0.49,0.74],mec=[0.27,0.49,0.74],color=[0.27,0.49,0.74])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Dead volume fraction (%)',ylabel='Sawlog ratio (%)')
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Yield vs. % dead (analysis of salvage)

list(dTM.keys())

ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']<1500) )[0]
x=dTM['Cruise_Pct Dead Net'][ind]
bw=20; bin=np.arange(0,120,bw)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
y=dTM['Cruise_V Gross (m3/ha)'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'go',ms=5,mfc='w',mec='g',label='Gross volume from cruise')
y=dTM['Cruise_V Net (m3/ha)'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'rs',ms=5,mfc='w',mec='r',label='Net volume from cruise')
y=dTM['V Logs Abs m3/ha'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'b^',ms=4,mfc='w',mec='b',label='Delivered volume from HBS')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Dead volume component (%)',ylabel='Net volume (m3/ha)')
ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Regression of sawlog ratio
def RegressionSR(meta,dTM):
	md=smf.ols("SR~Time+Age+PctDead+C(BGCZ)+C(SSC)+C(Ten)",data=df)
	rs=md.fit(maxiter=100)
	print(rs.summary())
	b=rs.params
	
	# Plot SR model behaviour
	Age=np.arange(1,300)
	AgeD=100; pdD=10; yrD=2020; bgcD='C(BGCZ)[T.SBS]'; sscD='C(SSC)[T.CLEAR]'; tenD='C(Ten)[T.Forest Licence]'
	plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,8));
	PctDead=0; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,0].plot(Age,y,'k-',label='% dead = 0')
	PctDead=100; Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,0].plot(Age,y,'k--',label='% dead = 100')
	ax[0,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[0,0].set(xlabel='Age, years',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	PctDead=pdD; Year=2010; bgc=bgcD; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,1].plot(Age,y,'k-',label='2010')
	PctDead=pdD; Year=2024; bgc=bgcD; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,1].plot(Age,y,'k--',label='2024')
	ax[0,1].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[0,1].set(xlabel='Age, years',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MH]'; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,2].plot(Age,y,'k-',label='MH')
	PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.CWH]'; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,2].plot(Age,y,'k--',label='CWH')
	PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.IDF]'; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,2].plot(Age,y,'k:',label='IDF')
	PctDead=pdD; Year=yrD; bgc='C(BGCZ)[T.MS]'; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[0,2].plot(Age,y,'k-.',label='MS')
	ax[0,2].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[0,2].set(xlabel='Age, years',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.CLEAR]'; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,0].plot(Age,y,'k-',label='CLEAR')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.IMCUT]'; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,0].plot(Age,y,'k--',label='IMCUT')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.RETEN]'; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,0].plot(Age,y,'k:',label='RETEN')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SEEDT]'; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,0].plot(Age,y,'k-.',label='SEEDT')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc='C(SSC)[T.SELEC]'; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,0].plot(Age,y,'k-',lw=1.5,label='SELEC')
	ax[1,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[1,0].set(xlabel='Age, years',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forest Licence]'
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,1].plot(Age,y,'k-',label='Forest Licence')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Forestry Licence to Cut]'
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,1].plot(Age,y,'k--',label='Forestry Licence to Cut')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.SB TSL S20 single mark]'
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,1].plot(Age,y,'k:',label='SB SSL S20 Single Mark')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Tree Farm Licence]'
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,1].plot(Age,y,'k-.',label='Tree Farm Licence')
	PctDead=pdD; Year=yrD; bgc=bgcD; ssc=sscD; ten='C(Ten)[T.Woodlot Licence]'
	y=np.maximum(0,b['Intercept']+b['Age']*Age+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,1].plot(Age,y,'k-',lw=1.5,label='Woodlot Licence')
	ax[1,1].set(xlabel='Age, years',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[1,1].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	PctDead=np.arange(0,100,1); Year=yrD; bgc=bgcD; ssc=sscD; ten=tenD
	y=np.maximum(0,b['Intercept']+b['Age']*60+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,2].plot(PctDead,y,'k-',label='Age = 60')
	y=np.maximum(0,b['Intercept']+b['Age']*200+b['Time']*Year+b['PctDead']*PctDead+b[ssc]+b[ten]+b[bgc])
	ax[1,2].plot(PctDead,y,'k--',label='Age = 200')
	ax[1,2].set(xlabel='Proportion of dead trees (%)',ylabel='Sawlog ratio (%)',ylim=[0,100])
	ax[1,2].legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	gu.axletters(ax,plt,0.03,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	return