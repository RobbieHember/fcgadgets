'''
HARVEST SUMMARY - COMMAND
'''
#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fcgadgets.hardhat.harvest_util as hut
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import statsmodels.formula.api as smf
gp=gu.SetGraphics('Manuscript')

#%% Import parameters
meta=u1ha.Init()

# List of grades
gradeL=list(meta['LUT']['Derived']['LogGrade'].keys())

#%% Processing of harvest DB by timber mark
flg=0
if flg==1:
	# Import HBS by TM and year
	dTMY=hut.CompileHarvest_ByTMY(meta)

	# Get area from forest tenure layer
	dTMY,dOP=hut.GetHarvestAreaFromFTEN(meta,dTMY)

	# Create TM DB from the TMY DB (sum all years wtihin TMs)
	dTM=hut.CalcHarvest_ByTM_FromTMY(meta,dTMY)

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
	#df=pd.DataFrame(dTMAndYear)
	#df.to_excel(meta['Paths']['DB']['Harvest'] + '\\HavestSummary_ByTMAndYear.xlsx',index=False)
else:
	# Import saved data
	dTMY=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HavestComp1_ByTMAndYear.pkl')
	dTM=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HavestComp1_ByTM.pkl')

#%% Harvest compilation 1 by time
hut.CalcHarvest_ByYear(meta)

#%% Plot time series of harvest volume
hut.Plot_HarvestVolumeTimeSeries(meta,dTMY)

#%% QA 15% of TMs in HBS are missing area from FTEN layer (private land?)
ind1=np.where(dTMY['PLANNED_NET_BLOCK_AREA']==0)[0]
ind2=np.where(dTMY['V Logs m3']>0)[0]
ind1.size/ind2.size

#%% QA Relationship between net volumes from different sources
hut.Plot_NetVolumeComparison(meta,dTM)

#%% Plot frequency of grades for each region
hut.Plot_YieldByGrade(meta,dTM)

#%% Frequency of tenure type
hut.Plot_TenureType_Frequency(meta,dTM)

#%% Frequency of SSC
hut.Plot_SilvicultureSystemCode_Frequency(meta,dTM)

#%% Relationship between total annual yield and low grades
dT=hut.StatsByTime(meta,dTMY)

x=(dT['Sum']['V Logs Grade 7 m3']+dT['Sum']['V Logs Grade 8 m3'])/1e6
ind=np.where( (x>0) )[0]
plt.close('all')
plt.plot(dT['Time'][ind],x[ind],'-ro')
plt.plot(dT['Time'][ind],dT['Sum']['V Felled m3'][ind]/1e6,'-bo')

#%% Relationship between grades 7-8 and total yield by District
hut.Plot_RelationshipBtwnGrades7and8andTotalYield(dTM)

#%%
plt.close('all')
plt.plot(dT['Mean']['Cruise_Pct Dead Net'],dT['Sum']['V Felled m3']/1e6,'-bo')

#%% Piling rate by BGC zone
hut.CalcPilingRateByBGCZone(meta,dTM,save='On',plot='On')

#%% Regression Waste Ratio
hut.WasteFractionRegression(meta,dTM,save='On',print_table='On',plot='On')

#%% Regression Sawlog Ratio
hut.SawlogFractionRegression(meta,dTM,save='On',print_table='On',plot='On')

#%%
d=hut.CalcStatsByDistrict(meta,dTM)

#%%
d=hut.CalcStatsByBGCZone(meta,dTM)

#%%
d=hut.CalcHarvestVolumeByBGCZ(meta,dTMY,dTM)

#%%
df=hut.CalcStatsBySILV_SYSTEM_CODE_And_Region(meta,dTM,'Interior')

#%%
zone='MS'
hut.CalcStatsBySILV_SYSTEM_CODE_ForBGCZone(meta,dTM,zone)

#%% Gross vs. net volume
hut.Plot_GrossVsNetVolume(dTM)

#%% Waste Ratio vs. % dead (analysis of salvage)
hut.Plot_WasteRatioVsPercentDead(dTM)

#%% Yield vs. % dead (analysis of salvage)
hut.Plot_YieldVsPercentDead(meta,dTM)

#%% Trying to identify whether grade 7 and 8 are associated with dead stands

#ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Abs m3/ha']<1500) )[0]
#x=dTM['Cruise_Pct Dead Net'][ind]
d={}
for gr in gradeL:
	ind=np.where( (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Grade ' + gr + ' m3/ha']>25) )[0]
	d[gr]=np.nanmean(dTM['Cruise_Pct Dead Net'][ind])


u,N=gu.CountByCategories(dTM['Tenure Type'],'Percent')
d={}
for i in range(u.size):
	ind=np.where( (dTM['Tenure Type']==u[i]) & (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Grade 7 m3/ha']<10000) & (dTM['V Logs Grade 8 m3/ha']<10000) )[0]
	d[u[i]]=np.nanmean(dTM['V Logs Grade 7 m3/ha'][ind]+dTM['V Logs Grade 8 m3/ha'][ind])

#%% Percent dead by tenure type
def Calc_PercentDeadByTenureType(dTM):
	d={}
	for i in range(u.size):
		ind0=np.where( (dTM['Tenure Type']==u[i]) & (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) )[0]
		ind1=np.where( (dTM['Tenure Type']==u[i]) & (dTM['Region']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['Cruise_Pct Dead Net']>50) )[0]
		try:
			d[u[i]]=ind1.size/ind0.size*100
		except:
			d[u[i]]=np.nan
	return

#%% QA - Look at stumpage by grade
d=hut.CalcStumpageByGrade(dTM)

#%% Trends in the proportion of young stands being harvested
hut.SummarizeAgeAtHarvest(dTM)

#%%