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

#%% Recover TSA
#meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
zTM=gis.OpenGeoTiff(r'C:\Data\BC1ha\FTEN_CUT_BLOCK_POLY_SVW\TIMBER_MARK.tif')
zTSA=gis.OpenGeoTiff(r'C:\Data\BC1ha\FADM_TSA\TSA_NUMBER_DESCRIPTION.tif')
from scipy import stats

N_mis=0
u=np.unique(dTM['TM'])
dTM['ID TSA']=np.zeros(dTM['TM'].size,dtype='int16')
for iU in range(u.size):
	try:
		id=meta['LUT']['FTEN_CUT_BLOCK_POLY_SVW']['TIMBER_MARK'][u[iU]]
	except:
		print('missing')
		N_mis=N_mis+1
		continue
	ind1=np.where(zTM['Data']==id)
	md=stats.mode(zTSA['Data'][ind1])[0]
	ind2=np.where(dTM['TM']==u[iU])[0]
	dTM['ID TSA'][ind2]=md
gu.opickle(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTM_New.pkl',dTM)

#%% QA 15% of TMs in HBS are missing area from FTEN layer (private land?)
ind1=np.where(dTMY['PLANNED_NET_BLOCK_AREA']==0)[0]
ind2=np.where(dTMY['V Logs m3']>0)[0]
ind1.size/ind2.size

#%% QA Relationship between net volumes from different sources
hut.NetVolumeComparison(meta,dTM)

#%% Frequency of grades for each region
hut.PlotYieldByGrade(meta,dTM)

#%% Frequency of tenure type
u,N=gu.CountByCategories(dTM['Tenure Type'],'Percent')
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
ax.barh(np.arange(u.size),N)
ax.set(yticks=np.arange(u.size),yticklabels=u)

#%% Frequency of SSC
u,N=gu.CountByCategories(dTM['SSC'],'Percent')
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,14));
ax.barh(np.arange(u.size),N)
ax.set(yticks=np.arange(u.size),yticklabels=u)

#%% Plot time series of harvest volume
hut.PlotHarvestVolume_TS(meta,dTMY)

#%% Piling rate
hut.CalcPilingRateByBGCZone(meta,dTM,save='On',plot='On')

#%% Regression Waste Ratio
hut.WasteFractionRegression(meta,dTM,save='On',print_table='On',plot='On')

#%% Regression Sawlog Ratio
hut.SawlogFractionRegression(meta,dTM,save='On',print_table='On',plot='On')

#%% QA - Look at stumpage by grade
for gr in gradeL:
	#dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/(453*dTM['V Logs Grade ' + gr + ' Abs m3']/1000)
	dTM['Stumpage Grade ' + gr]=dTM['Stump Logs Grade ' + gr + ' $']/dTM['V Logs Grade ' + gr + ' m3']

ind=np.where( (dTM['Stumpage Grade 1']>0) & (dTM['V Logs Grade 1 Abs m3']>0) )[0]
plt.hist(dTM['Stumpage Grade 1'][ind])

ind=np.where( (dTM['Stumpage Grade J']>0) & (dTM['V Logs Grade J Abs m3']>0) )[0]
plt.hist(dTM['Stumpage Grade J'][ind])

d={}
for k in dTM.keys():
	try:
		ind=np.where( (dTM['Region']=='Coast') & (dTM[k]>0) & (dTM[k]<1000000) )[0]
		d[k]=np.nanmean(dTM[k][ind])
	except:
		pass

d={}
for k in dTM.keys():
	try:
		ind=np.where( (dTM['Region']=='Interior') & (dTM[k]>0) & (dTM[k]<1000000) )[0]
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
df=hut.StatsBySILV_SYSTEM_CODE(meta,dTM,'Interior')

#df.to_excel(r'C:\Data\Harvest\ByTM_Boundary.xlsx')
#df.to_excel(r'C:\Data\Harvest\ByTM_Interior.xlsx')

#%%
def StatsBySILV_SYSTEM_CODE(meta,dTM,reg):
	d={}
	for k in meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'].keys():
		#
		ind=np.where( (dTM['PLANNED_NET_BLOCK_AREA']>0) & \
			(dTM['ID TSA']==meta['LUT']['FADM_TSA']['TSA_NUMBER_DESCRIPTION']['Boundary TSA']) & \
			(dTM['Region']==reg) & \
			(dTM['Cruise_V Gross (m3/ha)']>0) & \
			(dTM['Sawlog Fraction']>=0) & (dTM['Sawlog Fraction']<=1) & \
			(dTM['V Logs m3/ha']>0) & (dTM['V Logs m3/ha']<3000) & (dTM['SILV_SYSTEM_CODE % 1']>90) & (dTM['SILV_SYSTEM_CODE ID 1']==meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'][k]) )[0]
		if ind.size>3:
			d[k]=[ind.size,
			  np.round(np.nanmean(dTM['Cruise_V Gross (m3/ha)'][ind]),decimals=0),
			  np.round(np.nanmean(dTM['Cruise_V Net (m3/ha)'][ind]),decimals=0)]
		else:
			d[k]=[0,0,0]
	df=pd.DataFrame(d).T
	df.columns=['Sample size (# TMs)','Gross Volume (m3/ha)','Net Volume (m3/ha)']
	df=df.astype(int)
	df=df.reset_index()
	df=df.rename(columns={'index':'SSC'})
	return df

aw=((158*238)+(5*286))/(158+5)
aw*0.45*0.5

#%%
zone='MS'
hut.StatsBySILV_SYSTEM_CODE_ForBGCZone(meta,dTM,zone)

#%% Gross vs. net volume

ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['Cruise_V Net (m3/ha)']<3000) & (dTM['Cruise_V Gross (m3/ha)']>0) & (dTM['Cruise_V Gross (m3/ha)']<3000) & (dTM['V Logs Abs m3/ha']<2500) )[0]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot([0,1200],[0,1200],'k-',lw=2,color=[0.8,0.8,0.8])
x=dTM['Cruise_V Gross (m3/ha)'][ind]
y=dTM['Cruise_V Net (m3/ha)'][ind]
ax.plot(x,y,'k.',ms=5,mfc='b',mec='w',mew=0.5)

rs,txt=gu.GetRegStats(x,y,'No Intercept')
ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])

ax.set(xlabel='Gross volume (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,1200],ylim=[0,1200])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)


#%% Waste Ratio vs. % dead (analysis of salvage)
ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<2500) )[0]
x=dTM['Cruise_Pct Dead Net'][ind]
bw=20; bin=np.arange(0,120,bw)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
y=dTM['Waste Ratio'][ind]
#y=dTM['Sawlog Ratio'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'-bo',ms=5,mfc=[0.27,0.49,0.74],mec=[0.27,0.49,0.74],color=[0.27,0.49,0.74])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Dead volume fraction (%)',ylabel='Waste ratio (%)')
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

#%% Trying to identify whether grade 7 and 8 are associated with dead stands

#ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']<1500) )[0]
#x=dTM['Cruise_Pct Dead Net'][ind]
d={}
for gr in gradeL:
	ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Grade ' + gr + ' m3/ha']>25) )[0]
	d[gr]=np.nanmean(dTM['Cruise_Pct Dead Net'][ind])


u,N=gu.CountByCategories(dTM['Tenure Type'],'Percent')
y=np.zeros(u.size)
for i in range(u.size):
	ind=np.where( (dTM['Tenure Type']==u[i]) & (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Grade 8 m3/ha']<10000) )[0]
	y[i]=np.nanmean(dTM['V Logs Grade 7 m3/ha'][ind])





# #%% Relationship between Waste wood and net volume
# list(dTM.keys())
# ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['Waste Ratio']>=0) & (dTM['Waste Ratio']<2) )[0]
# x=dTM['V Felled m3/ha'][ind]
# y=dTM['Waste Total m3/ha'][ind]
# #rs,txt=gu.GetRegStats(x,y,'No Intercept')

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
# ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
# #ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
# #ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
# #ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
# #ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
# ax.set(ylabel='Waste wood (m3/ha)',xlabel='Net volume (m3/ha)',xlim=[0,1500],ylim=[0,300])
# plt.tight_layout()
# #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

# #%% Relationship between waste wood from HBS and WS
# list(dTM.keys())
# ind=np.where( (dTM['V Logs Waste m3/ha']>0) & (dTM['V Logs Waste m3/ha']<400) & (dTM['Waste Total m3/ha']>=0) & (dTM['Waste Total m3/ha']<400) )[0]
# x=dTM['V Logs Waste m3/ha'][ind]
# y=dTM['Waste Total m3/ha'][ind]
# rs,txt=gu.GetRegStats(x,y,'No Intercept')

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
# ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
# ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
# ax.text(300,300,'1:1',fontsize=10,ha='center',va='center')
# ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
# ax.text(385,60,txt,fontsize=8,ha='right',va='center')
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
# ax.set(ylabel='Waste wood (m3/ha)',xlabel='Net volume (m3/ha)',xlim=[0,400],ylim=[0,400])
# plt.tight_layout()
# #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

# #%% Relationship between WR and net volume
# list(dTM.keys())
# ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['Waste Ratio']>=0) & (dTM['Waste Ratio']<2) )[0]
# x=dTM['V Logs Abs m3/ha'][ind]
# y=dTM['Waste Ratio'][ind]
# #rs,txt=gu.GetRegStats(x,y,'No Intercept')

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
# ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
# #ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
# #ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
# #ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
# #ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
# ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,1500],ylim=[0,2])
# plt.tight_layout()
# #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

# #%% Relationship between WR and age
# list(dTM.keys())
# ind=np.where( (dTM['V Logs Abs m3/ha']>0) & (dTM['V Logs Abs m3/ha']<1500) & (dTM['WR']>=0) & (dTM['WR']<2) )[0]
# x=dTM['age_vri02'][ind]
# y=dTM['WR'][ind]
# #rs,txt=gu.GetRegStats(x,y,'No Intercept')

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
# ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
# #ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
# #ax.text(1300,1300,'1:1',fontsize=10,ha='center',va='center')
# #ax.plot(rs['xhat Line'],rs['yhat Line'],'r--',lw=1.5,color=[0,0,0])
# #ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
# ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume (m3/ha)',xlim=[0,300],ylim=[0,2])
# plt.tight_layout()
# #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)
