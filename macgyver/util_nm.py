'''
UTILITIES - NUTRIENT MANAGEMENT COMPLETED
'''

#%% Import Python modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import copy
import warnings
import matplotlib as mpl
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_inventory as uinv
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun as cbr
import fcgadgets.macgyver.util_fcs_graphs as ufcs
import fcgadgets.macgyver.util_fcs_qa as uqa
warnings.filterwarnings("ignore")
gp=gu.SetGraphics('Manuscript')

#%%
def ImportModelResults():
	pNam='BCFCS_NMC'
	meta=gu.ipickle(r'C:\Data\BCFCS\BCFCS_NMC\Inputs\Metadata.pkl')
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	mos=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNam)
	mos[pNam]['Delta']={}
	mos[pNam]['Delta']['Low Harvest']={'iB':0,'iP':1}
	mos[pNam]['Delta']['High Harvest']={'iB':2,'iP':3}
	mos=cbu.Import_MOS_ByScnComparisonAndStrata(meta,pNam,mos)
	mos=cbu.Import_MOS_ByScnAndStrata_Area(meta,pNam,mos)
	dmec=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\dmec.pkl')
	#thlb=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\thlb.pkl')
	#nPS='All'
	#nSS='All'
	#nYS='All'
	#iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
	#iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
	#iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
	meta[pNam]['Project']['Multi']=1e6

	# Add future projected implementation
	meta[pNam]['AIL CAP']=30000
	mosWF=ProjectFuture(meta,pNam,mos)
	return pNam,meta,tv,mos,mosWF,dmec

#%% Save summary of DMEC and results to spreadsheet for troubleshooting, investigation, QA
def QA_PrintSummarySparseToSpreadsheet(meta,pNam,iScn,dmec,E):
	#iScn=1

	# Get size
	n=0
	for i in range(len(dmec[iScn])):
		n=n+dmec[iScn][i]['Year'].size

	d={}
	d['Index']=np.zeros(n)
	d['Year']=np.zeros(n)
	d['DOY']=np.zeros(n)
	d['Event Type']=np.array(['' for _ in range(n)],dtype=object)
	d['Event Source']=np.array(['Inventory' for _ in range(n)],dtype=object)
	d['Mortality']=np.zeros(n)
	d['Year Focal']=np.zeros(n)
	#d['ASET Focal']=np.array(['' for _ in range(n)],dtype=object)
	#d['ASET']=np.array(['' for _ in range(n)],dtype=object)
	#d['Inciting NOSE']=np.array(['' for _ in range(n)],dtype=object)
	d['FSC Focal']=np.array(['' for _ in range(n)],dtype=object)
	d['FSC']=np.array(['' for _ in range(n)],dtype=object)
	d['OPENING']=np.array(['' for _ in range(n)],dtype=object)
	d['BGCZ']=np.array(['' for _ in range(n)],dtype=object)
	d['E_Atmosphere_SubstitutionExcluded Delta']=np.zeros(n)
	d['Pl Spc1 CD']=np.array(['' for _ in range(n)],dtype=object)
	d['Pl Spc1 %']=np.zeros(n)
	d['Pl SPH']=np.zeros(n)

	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\OPENING_ID.tif')
	zOPID=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
	zBGCZ=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

	c=0
	for i in range(len(dmec[iScn])):

		ind=np.where( (dmec[iScn][i]['ID Event Type']==meta['LUT']['Event']['Nutrient App Aerial']) )[0]
		if ind.size>0:
			year_focal=dmec[iScn][i]['Year'][ind[-1]]
			fsc_focal=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][ind[-1]])[0]
		else:
			year_focal=0
			fsc_focal='None'

		for j in range(dmec[iScn][i]['Year'].size):
			d['Index'][c]=i
			d['Year'][c]=dmec[iScn][i]['Year'][j]
			d['DOY'][c]=dmec[iScn][i]['DOY'][j]
			d['Mortality'][c]=dmec[iScn][i]['Mortality Factor'][j]
			d['Event Type'][c]=cbu.lut_n2s(meta['LUT']['Event'],dmec[iScn][i]['ID Event Type'][j])[0]
			if dmec[iScn][i]['Event Source'][j]==2:
				d['Event Source'][c]='Added'
			if dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][j]!=9999:
				d['FSC'][c]=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][j])[0]
			d['Year Focal'][c]=int(year_focal)
			d['FSC Focal'][c]=fsc_focal
			d['OPENING'][c]=zOPID[i]
			d['BGCZ'][c]=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],zBGCZ[i])[0]
			d['E_Atmosphere_SubstitutionExcluded Delta'][c]=int(E['Delta']['E_Atmosphere_SubstitutionExcluded'][i])
			d['Pl Spc1 CD'][c]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],dmec[iScn][i]['PL_SPECIES_CD1'][j])[0]
			d['Pl Spc1 %'][c]=dmec[iScn][i]['PL_SPECIES_PCT1'][j]
			d['Pl SPH'][c]=dmec[iScn][i]['Planted SPH'][j]
			c=c+1
	df=pd.DataFrame.from_dict(d)
	df.to_excel(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummarySparseSample.xlsx',index=False)
	return

#%%
def SummarizeServicePlan(meta,pNam,mos):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	v='E_Atmosphere_SubstitutionExcluded'
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) & (tv<=2050) )[0]
	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
	d={}
	d['Total']={}
	d['PerHa']={}
	for cNam in mos[pNam]['Delta'].keys():
		d['Total'][cNam]={}
		d['PerHa'][cNam]={}
		for iSS in range(meta[pNam]['Project']['Strata']['Spatial']['Unique CD'].size):
			cd=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
			y=np.sum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT,0,iSS,iYS])/1e6
			d['Total'][cNam][cd]=np.round(y,decimals=2)
			y=np.sum(mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,0,iSS,iYS])
			d['PerHa'][cNam][cd]=np.round(y,decimals=2)
	return d

#%%
def QA_Plot_CarbonTimeSeriesByStandSelect(meta,pNam,dmec,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	for iBat in range(meta[pNam]['Project']['N Batch']):
		indBat=cbu.IndexToBatch(meta[pNam],iBat)
		d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,iBat)
		d1=cbu.LoadSingleOutputFile(meta,pNam,1,0,iBat)	

		ind=np.where(np.isin(indBat,kwargs['Index'])==True)[0]
		if ind.size==0:
			continue
		for i in range(ind.size):
			plt.close('all'); fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(18,12));
			ax[0,0].plot(tv,d0['A'][:,ind[i]],'b-')
			ax[0,0].plot(tv,d1['A'][:,ind[i]],'g--')
			ax[0,0].set(ylabel='Age')
			ax[0,1].plot(tv,d0['C_Biomass'][:,ind[i]],'b-',label='Baseline')
			ax[0,1].plot(tv,d1['C_Biomass'][:,ind[i]],'g--',label='Actual')
			ax[0,1].set(ylabel='Biomass')
			ax[0,1].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)

			ax[1,0].plot(tv,d0['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'ob-',ms=3)
			ax[1,0].plot(tv,d1['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'sg--',ms=3)
			ax[1,0].set(ylabel='Wildfire emissions')
			ax[1,1].plot(tv,d0['V_ToMill_MerchTotal'][:,ind[i]],'ob-')
			ax[1,1].plot(tv,d1['V_ToMill_MerchTotal'][:,ind[i]],'sg--')
			ax[1,1].set(ylabel='Volume to mill, merch live')

			ax[2,0].plot(tv,d0['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'b-')
			ax[2,0].plot(tv,d1['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'g--')
			ax[2,0].set(ylabel='Open burning emissions')

			ax[2,1].plot(tv,d0['C_Litter'][:,ind[i]]+d0['C_Soil'][:,ind[i]]+d0['C_DeadWood'][:,ind[i]],'b-')
			ax[2,1].plot(tv,d1['C_Litter'][:,ind[i]]+d1['C_Soil'][:,ind[i]]+d0['C_DeadWood'][:,ind[i]],'g--')
			ax[2,1].set(ylabel='Litter+Soil+Dead Wood')

			ax[3,0].plot(tv,d0['E_Atmosphere_SubstitutionExcluded'][:,ind[i]],'b-')
			ax[3,0].plot(tv,d1['E_Atmosphere_SubstitutionExcluded'][:,ind[i]],'g--')
			ax[3,0].set(ylabel='Annual GHG balance')
			ax[3,1].plot(tv,d0['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,ind[i]],'b-')
			ax[3,1].plot(tv,d1['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,ind[i]],'g--')
			ax[3,1].set(ylabel='Cumulative GHG balance')

			ind1=np.where(dmec[1][indBat[ind[i]]]['ASET']<9999)[0]
			yr1=dmec[1][indBat[ind[i]]]['Year'][ind1[-1]]
			plt.tight_layout()
			type=cbu.lut_n2s(meta['LUT']['Derived']['ASET'],meta[pNam]['Project']['ASET'][indBat[ind[i]]])[0]
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA\\' + type + '_' + str(indBat[ind[i]]) + '_' + str(int(yr1)),'png',200)

	return

#%%
def PlotSummaryCombined(meta,mos,pNamC,pNamF,cNam):
	tv=np.arange(meta[pNamC]['Project']['Year Start Saving'],meta[pNamC]['Project']['Year End']+1,1)
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,9));

	# Annual
	pNam=pNamC
	v='E_Atmosphere_SubstitutionIncluded'
	yC=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi']
	yC[(tv<1971)]=0
	v='V_ToMill_MerchTotal'
	vC=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi']
	vC[(tv<1971)]=0
	pNam=pNamF
	v='E_Atmosphere_SubstitutionIncluded'
	yF=gu.movingave(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi'],50,'center')
	yF[(tv<2023)]=0
	v='V_ToMill_MerchTotal'
	vF=gu.movingave(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi'],50,'center')
	vF[(tv<2023)]=0
	ax[0,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,0].plot(tv,vC+vF,'k-',lw=1.25,color=[0.6,0.9,0],label='Completed + climate action planning')
	ax[0,0].plot(tv,vC,'g--',lw=1.25,color=[0.29,0.47,0.77],label='Completed')
	ax[0,0].set(xticks=np.arange(1950,2200,20),ylabel='Annual $\Delta$ volume yield (Mm$^3$ yr$^{-1}$)',xlabel='Time, years',xlim=[1975-0.5,2100+0.5],ylim=[-0.7,0.7])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,0].legend(frameon=False,loc='lower left',facecolor=[1,1,1],labelspacing=0.25)
	ax[0,1].plot(tv,np.zeros(tv.size),'k-',lw=1.25,color=[0.8,0.8,0.8])
	ax[0,1].plot(tv,yC+yF,'k-',lw=1.25,color=[0.6,0.9,0])
	ax[0,1].plot(tv,yC,'g--',lw=1,color=[0.29,0.47,0.77],label='Completed')
	ax[0,1].set(xticks=np.arange(1950,2200,20),ylabel='Annual $\Delta$E (MtCO$_2$e yr$^{-1}$)',xlabel='Time, years',xlim=[1975-0.5,2100+0.5],ylim=[-2.5,1])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	# Cumulative
	pNam=pNamC
	v='E_Atmosphere_SubstitutionIncluded'
	yC=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi']
	yC[(tv<1971)]=0
	v='V_ToMill_MerchTotal'
	vC=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi']
	vC[(tv<1971)]=0
	pNam=pNamF
	v='E_Atmosphere_SubstitutionIncluded'
	yF=gu.movingave(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi'],50,'center')
	yF[(tv<2023)]=0
	v='V_ToMill_MerchTotal'
	vF=gu.movingave(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNamC]['Project']['Multi'],50,'center')
	vF[(tv<2023)]=0
	ax[1,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,0].plot(tv,np.cumsum(vC+vF),'k-',lw=1.25,color=[0.6,0.9,0],label='Completed + Climate Action Plan')
	ax[1,0].plot(tv,np.cumsum(vC),'g--',lw=1.25,color=[0.29,0.47,0.77],label='Completed')
	ax[1,0].set(xticks=np.arange(1950,2200,20),ylabel='Cumulative $\Delta$ volume yield (Mm$^3$)',xlabel='Time, years',xlim=[1975-0.5,2100+0.5],ylim=[0,18])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(frameon=False,loc='upper center',facecolor=[1,1,1],labelspacing=0.25)
	ax[1,1].plot(tv,np.zeros(tv.size),'k-',lw=1.25,color=[0.8,0.8,0.8])
	ax[1,1].plot(tv,np.cumsum(yC+yF),'k-',lw=1.25,color=[0.6,0.9,0])
	ax[1,1].plot(tv,np.cumsum(yC),'g--',lw=1,color=[0.29,0.47,0.77],label='Completed')
	ax[1,1].set(xticks=np.arange(1950,2200,20),ylabel='Cumulative $\Delta$E (MtCO$_2$e)',xlabel='Time, years',xlim=[1975-0.5,2100+0.5],ylim=[-110,0])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\SummaryCombined','png',200)
	return

#%%
def Plot_AreaTreated_TimeSeriesByFSC(meta,pNam):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_SummaryByTimeAndFSC.pkl')

	# Treatment area differs from number of treatments due to repeat applications on the same area
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_Count.tif')['Data']
	print(np.where(z>0)[0].size)
	print(np.where(z>1)[0].size)
	cl=np.random.random((d['FSC'].size,3))

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,5.75));
	A_Cumulative=np.zeros(d['Year'].size)
	for iFSC in range(d['FSC'].size):
		plt.bar(d['Year'],d['Area'][:,iFSC]/1000,0.8,bottom=A_Cumulative,facecolor=cl[iFSC,:],label=d['FSC'][iFSC])
		A_Cumulative=A_Cumulative+d['Area'][:,iFSC]/1000
	ax.set(xticks=np.arange(1950,2025+1,5),ylabel='Treatment area (hectares x 1000)',xlabel='Time, years',xlim=[1975-0.5,d['Year'][-1]+0.5],ylim=[0,47])
	plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=3)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
# 	elif kwargs['Period']=='CF':
# 		tv=np.arange(meta[pNamC]['Project']['Year Start Saving'],meta[pNamC]['Project']['Year End']+1,1)
# 		iT=np.where(ail['Year']<=2022)[0]
# 		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(24,6))
# 		A_Cumulative=np.zeros(iT.size)
# 		for iFSC in range(ail['FSC Unique'].size):
# 			plt.bar(ail['Year'][iT],ail['A Unique'][iT,iFSC]/1000,0.8,bottom=A_Cumulative,facecolor=cl[iFSC,:],label=ail['FSC Unique'][iFSC])
# 			A_Cumulative=A_Cumulative+ail['A Unique'][iT,iFSC]/1000
# 		iPS=0; iSS=0;iYS=0
# 		plt.bar(tv,mos[pNamF]['Scenarios'][1]['Sum']['Area_Nutrient App Aerial']['Ensemble Mean'][:,iPS,iSS,iYS]*meta[pNamF]['Project']['AEF']/1000,0.8,fc=[0.7,0.7,0.7],label='CAP')
# 		ax.set(xticks=np.arange(1950,2225+1,10),ylabel='Implimentation (hectares x 1000)',xlabel='Time, years',xlim=[ail['Year'][0]-0.75,2101],ylim=[0,65])
# 		plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=7)
# 		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
# 		plt.tight_layout()
# 	else:
# 		pass
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TimeSeries_ByFSC','png',900)
	return

#%%
def Plot_AreaTreated_TimeSeriesByFSC_WithFuture(meta,pNam):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_SummaryByTimeAndFSC.pkl')
	cl=np.random.random((d['FSC'].size,3))

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,5.75));
	A_Cumulative=np.zeros(d['Year'].size)
	for iFSC in range(d['FSC'].size):
		plt.bar(d['Year'],d['Area'][:,iFSC]/1000,0.8,bottom=A_Cumulative,facecolor=cl[iFSC,:],label=d['FSC'][iFSC])
		A_Cumulative=A_Cumulative+d['Area'][:,iFSC]/1000

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where(tv>meta[pNam]['YearCurrent'])[0]
	#plt.bar(tv[iT],np.ones(iT.size)*meta[pNam]['AIL CAP']/1000,0.6,fc=[0.7,0.7,0.7],label='FIP Climate Action Plan')
	plt.plot(tv[iT],np.ones(iT.size)*meta[pNam]['AIL CAP']/1000,'-ko',ms=2,mfc='k',mec='k',label='CAP')

	ax.set(xticks=np.arange(1950,2025+1,5),ylabel='Treatment area (hectares x 1000)',xlabel='Time, years',xlim=[1975-0.5,2100+0.5],ylim=[0,47])
	plt.legend(frameon=False,loc='lower right',facecolor=[1,1,1],labelspacing=0.25,ncol=5)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TimeSeries_ByFSC_WithFuture','png',900)
	return

#%%
def Plot_AreaTreated_FrequencyByFSC(meta,t0,t1):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_SummaryByTimeAndFSC.pkl')
	iT=np.where( (d['Year']>=t0) & (d['Year']<=t1) )[0]

	A=np.sum(d['Area'][iT,:],axis=0)/1000
	ord=np.flip(np.argsort(A))
	A=A[ord]
	cd=d['FSC'][ord]
	ind=np.where(A>0)[0]
	A=A[ind]
	cd=cd[ind]
	bin=np.arange(1,A.size+1)
	A_max=np.max(A)

	if A_max<100:
		yl=[0,A_max+5]
		yo=0.95
	elif A_max>250:
		yl=[0,A_max+50]
		yo=4
	else:
		yl=[0,A_max+20]
		yo=1.35

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10.5));
	plt.bar(bin,A,0.85)
	for i in range(bin.size):
		y=A[i]
		ax.text(bin[i],y+yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Implementation level (Kha)',xticklabels=cd,xlim=[0.25,bin.size+0.75],ylim=yl) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyByFSC_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def Plot_EmissionsCumu_TimeSeriesFromCurrentYear(meta,mos,pNam,cNam,year):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	flg=0
	if flg==1:
		# Add Nutrient management
		pNamNMC='BCFCS_NMC'
		pth=r'D:\Data\FCI_Projects\BCFCS_NMC\Inputs\Metadata.pkl'
		metaNM=gu.ipickle(pth)
		yNM=np.nan*np.ones( (tv.size,metaNM[pNamNMC]['Project']['N Stand']) )
		for iBat in range(metaNM[pNamNMC]['Project']['N Batch']):
			indBat=cbu.IndexToBatch(metaNM[pNamNMC],iBat)
			indY=np.where(metaNM[pNamNMC]['Project']['Strata']['Year']['ID'][indBat]==2022)[0]
			d0=cbu.LoadSingleOutputFile(metaNM,pNamNMC,0,0,iBat)
			d1=cbu.LoadSingleOutputFile(metaNM,pNamNMC,1,0,iBat)
			y0=(d0['E_Atmosphere_SubstitutionExcluded'][:,indY]+d0['E_Atmosphere_SubstitutionIncluded'][:,indY])/2
			y1=(d1['E_Atmosphere_SubstitutionExcluded'][:,indY]+d1['E_Atmosphere_SubstitutionIncluded'][:,indY])/2
			yNM[:,indBat[indY]]=y1-y0
	
	vStat='Ensemble Mean'
	oper='Sum'
	iT2=np.where( (tv>=2020) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(year))[0][0]

		# Half substitutions
		#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,iPS,0,iYS])/2)/1e6

		# Without substitutions
		y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]/1e6

		ysum=ysum+np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)

	plt.plot(tv[iT2],np.cumsum(ysum[iT2]),'--',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	flg=0
	if flg==1:
		# Add Nutrient management
		#pNam=pNamNMC
		#iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='2022')[0][0]
		#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,0,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,0,0,iYS])/2)/1e6
		mu=np.nanmean(yNM,axis=1)
		iT3=np.where(tv<2022)[0]
		mu[iT3]=0
		A=33500
		yNM1=A*np.cumsum(mu)/1e6
		plt.plot(tv[iT2],yNM1[iT2],'-',ms=3,color=[0.5,0,1],mec=[0.5,0,1],mfc='w',lw=1,mew=0.75,label='Nutrient Management')
		plt.plot(tv[iT2],yNM1[iT2]+np.cumsum(ysum[iT2]),'--',ms=3,color='k',mec='k',mfc='k',lw=1.5,mew=0.75,label='Total')
	
	ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(2022,-1,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(2022,1,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2100,5),yticks=np.arange(-20,20,1),ylabel='Cumulative $\Delta$E (MtCO$_2$e)',xlabel='Time, years',ylim=[-10,2],xlim=[2020,2075])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumu_TimeSeries_FromYear' + str(year),'png',900)#%%
	return

#%%
def SummarizeBenefitForEachStand(meta,pNam,mos):
	# Sum
	vL=['E_Atmosphere_SubstitutionIncluded','E_Atmosphere_SubstitutionIncludedHalf','E_Atmosphere_SubstitutionExcluded']
	t0=1960
	t1=2050
	E_sum={}
	iScn=0;E_sum['Baseline']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Sum')
	iScn=1;E_sum['Actual']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Sum')
	E_sum['Delta']={}
	for k in E_sum['Baseline'].keys():
		E_sum['Delta'][k]=E_sum['Actual'][k]-E_sum['Baseline'][k]
	
	# Mean
	vL=['A','C_Biomass','C_Litter','C_Soil','C_Soil_OHorizon','C_DeadWood','C_G_Gross',
		'C_G_Net','C_M_Reg','C_M_Dist','C_LF','C_ToMillMerchGreen',
		'C_ToMillNonMerchGreen','C_ToMillMerchDead','C_ToMillNonMerchDead','C_ToPileBurnTot','E_Domestic_ForestSector_OpenBurning',
		'E_Atmosphere_SubstitutionIncluded','E_Atmosphere_SubstitutionIncludedHalf','E_Atmosphere_SubstitutionExcluded']
	E_mu={}
	iScn=0;E_mu['Baseline']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Mean')
	iScn=1;E_mu['Actual']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Mean')
	E_mu['Delta']={}
	for k in E_mu['Baseline'].keys():
		E_mu['Delta'][k]=E_mu['Actual'][k]-E_mu['Baseline'][k]
	return E_sum,E_mu

#%%
def Plot_TreatmentFrequency_BySpecies(meta,pNam,t0,t1):
	d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummarySpeciesComposition_' + str(t0) + 'to' + str(t1) + '.pkl')
	
	cd=np.array(list(d.keys()))
	N=np.array(list(d.values()))
	ord=np.flip(np.argsort(N))
	cd=cd[ord]
	N=N[ord]
	Np=N/np.sum(N)*100
	ind=np.where(Np>0.1)[0]
	cd=cd[ind]
	N=N[ind]
	Np=Np[ind]
	bin=np.arange(1,ind.size+1)
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10));
	plt.bar(bin,Np,0.85,facecolor=[0.25,0.5,1])
	# for i in range(bin.size):
	# 	y=N[i]
	# 	ax.text(bin[i],y+yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Frequency (%)',xticklabels=cd,xlim=[0.25,bin.size+0.75]) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyBySpecies_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def UpdateDatabase(meta,mos,pNam,cNam,YearCurrent):

	ac1='Silviculture investments'
	ac2='Non-obligation stand establishment'
	vStat='Ensemble Mean'
	oper='Sum'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	n=5000
	d={}
	d['Action Category 1']=np.array(['' for _ in range(n)],dtype=object)
	d['Action Category 2']=np.array(['' for _ in range(n)],dtype=object)
	d['Action Category 3']=np.array(['' for _ in range(n)],dtype=object)
	d['Funding Source Code']=np.array(['' for _ in range(n)],dtype=object)
	d['Scope Investments']=np.array(['' for _ in range(n)],dtype=object)
	d['Year Implemented']=np.array(['' for _ in range(n)],dtype=object)
	d['Area treated (ha/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2030 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2040 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2050 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2070 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2100 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta annual emissions 2120 wo subst (tCO2e/yr)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2030 wo subst (tCO2e)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2040 wo subst (tCO2e)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2050 wo subst (tCO2e)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2070 wo subst (tCO2e)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2100 wo subst (tCO2e)']=np.zeros(n,dtype='int32')
	d['Delta cumulative emissions 2120 wo subst (tCO2e)']=np.zeros(n,dtype='int32')

	d['Delta sum cost silviculture 2030 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost silviculture 2040 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost silviculture 2050 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost silviculture 2070 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost silviculture 2100 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost silviculture 2120 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2030 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2040 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2050 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2070 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2100 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum cost total 2120 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2030 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2040 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2050 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2070 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2100 (CAD)']=np.zeros(n,dtype='int32')
	d['Delta sum net revenue sector 2120 (CAD)']=np.zeros(n,dtype='int32')

	d['Mitigation value upfront undiscounted wo subst 2030 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront undiscounted wo subst 2040 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront undiscounted wo subst 2050 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront undiscounted wo subst 2070 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront undiscounted wo subst 2100 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront undiscounted wo subst 2120 (CAD/tCO2e)']=np.zeros(n,dtype='float')

	d['Mitigation value upfront discounted wo subst 2030 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront discounted wo subst 2040 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront discounted wo subst 2050 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront discounted wo subst 2070 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront discounted wo subst 2100 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value upfront discounted wo subst 2120 (CAD/tCO2e)']=np.zeros(n,dtype='float')

	d['Mitigation value net sector undiscounted wo subst 2030 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector undiscounted wo subst 2040 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector undiscounted wo subst 2050 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector undiscounted wo subst 2070 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector undiscounted wo subst 2100 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector undiscounted wo subst 2120 (CAD/tCO2e)']=np.zeros(n,dtype='float')

	d['Mitigation value net sector discounted wo subst 2030 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector discounted wo subst 2040 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector discounted wo subst 2050 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector discounted wo subst 2070 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector discounted wo subst 2100 (CAD/tCO2e)']=np.zeros(n,dtype='float')
	d['Mitigation value net sector discounted wo subst 2120 (CAD/tCO2e)']=np.zeros(n,dtype='float')

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysL=meta[pNam]['Project']['Strata']['Year']['Unique CD']
	ssL=meta[pNam]['Project']['Strata']['Spatial']['Unique CD']
	cnt=0
	for ps in psL:
		for ss in ssL:
			for ys in ysL:
				iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==ps)[0]
				iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==ys)[0]
				iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==ss)[0]

				d['Action Category 1'][cnt]=ac1
				d['Action Category 2'][cnt]=ac2
				d['Action Category 3'][cnt]=ps
				d['Funding Source Code'][cnt]=ss
				d['Year Implemented'][cnt]=ys
				d['Scope Investments'][cnt]='Completed'
				d['Area treated (ha/yr)'][cnt]=0

				iT=np.where( (tv==2030) )[0][0]
				d['Delta annual emissions 2030 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])
				iT=np.where( (tv==2040) )[0]
				d['Delta annual emissions 2040 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])
				iT=np.where( (tv==2050) )[0]
				d['Delta annual emissions 2050 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])
				iT=np.where( (tv==2070) )[0]
				d['Delta annual emissions 2070 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])
				iT=np.where( (tv==2100) )[0]
				d['Delta annual emissions 2100 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])
				iT=np.where( (tv==2120) )[0]
				d['Delta annual emissions 2120 wo subst (tCO2e/yr)'][cnt]=np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS])

				iT=np.where( (tv>=1960) & (tv<=2030) )[0]
				d['Delta cumulative emissions 2030 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				iT=np.where( (tv>=1960) & (tv<=2040) )[0]
				d['Delta cumulative emissions 2040 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				iT=np.where( (tv>=1960) & (tv<=2050) )[0]
				d['Delta cumulative emissions 2050 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				iT=np.where( (tv>=1960) & (tv<=2070) )[0]
				d['Delta cumulative emissions 2070 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				iT=np.where( (tv>=1960) & (tv<=2100) )[0]
				d['Delta cumulative emissions 2100 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				iT=np.where( (tv>=1960) & (tv<=2120) )[0]
				d['Delta cumulative emissions 2120 wo subst (tCO2e)'][cnt]=np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][iT,iPS,iSS,iYS]))
				cnt=cnt+1

	# Truncate
	for k in d.keys():
		d[k]=d[k][0:cnt]

	# Save to file
	df=pd.DataFrame.from_dict(d)
	df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Projects\Forest Carbon Summary\Data' + '\\' + str(YearCurrent) + '\\BC-FCS Database Summary By Year AC and FS.xlsx')

	return

#%%
# Raidial lines and CPT text
#yrs=[1990,2020,2050,2070,2100]
#ha=['right','right','left','left','left','left']
#xa=[1.8,1.35,1.2,1.06,1.035,1.07]

def Plot_MitigationValueUndiscounted(meta,mos,pNam,iPS,iSS,iYS):
	#tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,10)); fs2=7
	ax.plot([0,0],[-2000,2000],'k-',lw=0.5,color='k')
	ax.plot([-10000,10000],[0,0],'k-',lw=0.5,color='k')
	rings=np.arange(100,1200,100)
	for ring in rings:
		ax.plot(0,0,'o',ms=ring,mec=[0.8,0.8,0.8],mfc='none',mew=0.5)
	cpt=np.array([160,80,40,20,10,5,1],dtype=float)
	#x_lab=np.array([-3000,-2800,-2400,-1625,-910,-470,-300])
	#Slope=np.zeros(cpt.size)
	for iB in range(cpt.size):
		ax.plot([0,-10000],[0,-10000/cpt[iB]],'k-',color=[0.8,0.8,0.8])
	for iB in range(cpt.size):
		ax.plot([0,10000],[0,-10000/cpt[iB]],'k-',color=[0.8,0.8,0.8])
	#ax.text(xC,-y_new,int(cpt[iB]),fontsize=fs2)

	cNam='Low Harvest'
	dE=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dC=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.plot(-dRN,dE,'k-',color=[0,0,1],lw=1.5,label='High harvest (w/o substitutions)')

	cNam='Low Harvest'
	dE=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncludedHalf']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.plot(-dRN,dE,'k:',color=[0,1,1],lw=1.5,label='High harvest (with substitutions)')

	cNam='High Harvest'
	dE=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dC=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	#ax.plot(-1*dC,dE,'k-',color=[0.27,0.49,0.77],lw=1.5,label='Cost')
	ax.plot(-dRN,dE,'k--',color=[0,1,0],lw=1.5,label='Low harvest (w/o substitutions)')

	ax.set(xticks=np.arange(-12000,18000,1000),xlabel='Cumulative cost (Million CAD)',ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[-5000,2000],ylim=[-120,10])
	ax.legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.66,1) bbox_to_anchor=(0.3,0.9,0,0),
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
	ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	plt.tight_layout()

	flg=0
	if flg==1:
		def onclick(event):
			global ix, iy
			ix,iy = event.xdata, event.ydata
			global xt,yt
			xt.append(ix)
			yt.append(iy)
			return xt,yt

		xt=[]; yt=[]
		for i in range(7):
			cid=fig.canvas.mpl_connect('button_press_event', onclick)

	xt=np.array([-2037., -1978., -1824., -1417.,  -886.,  -484.,  -112.])
	yt=np.array([-13., -25., -45., -71., -89., -96., -99.])
	for i in range(len(xt)):
		ax.text(xt[i],yt[i],int(cpt[i]),fontsize=fs2,color=[0,0,0])
	for i in range(len(xt)):
		ax.text(-xt[i],yt[i],int(cpt[i]),fontsize=fs2,color=[0,0,0])

	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\MitigationValue_CompletedAndCAP_Undiscounted','png',900);

	return

#%%
def Plot_MitigationValueDiscounted(meta,mos,pNam,iPS,iSS,iYS):

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,10)); fs2=7
	ax.plot([0,0],[-2000,2000],'k-',lw=0.5,color='k')
	ax.plot([-10000,10000],[0,0],'k-',lw=0.5,color='k')
	rings=np.arange(100,1200,100)
	for ring in rings:
		ax.plot(0,0,'o',ms=ring,mec=[0.8,0.8,0.8],mfc='none',mew=0.5)
	cpt=np.array([160,80,40,20,10,5,1],dtype=float)
	#x_lab=np.array([-3000,-2800,-2400,-1625,-910,-470,-300])
	#Slope=np.zeros(cpt.size)
	for iB in range(cpt.size):
		ax.plot([0,-10000],[0,-10000/cpt[iB]],'k-',color=[0.8,0.8,0.8])
	#ax.text(xC,-y_new,int(cpt[iB]),fontsize=fs2)

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	r_disc=meta['Param']['BE']['Biophysical']['Discount Rate Emissions']
	t_disc=np.maximum(0,tv-meta[pNam]['Project']['Year Project'])
	t_disc=np.tile(t_disc,(v0['E_Atmosphere_SubstitutionExcluded'].shape[1],1)).T

	v0['E_Atmosphere_SubstitutionIncluded_Disc']=v0['E_Atmosphere_SubstitutionIncluded'].copy()/((1+r_disc)**t_disc)

	cNam='Low Harvest'
	E=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
	dE=np.cumsum()
	dC=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.plot(-dRN,dE,'k-',color=[0,0,1],lw=1.5,label='High harvest (w/o substitutions)')

	cNam='Low Harvest'
	dE=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncludedHalf']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.plot(-dRN,dE,'k:',color=[0,1,1],lw=1.5,label='High harvest (with substitutions)')

	cNam='High Harvest'
	dE=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dC=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	dRN=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][:,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	#ax.plot(-1*dC,dE,'k-',color=[0.27,0.49,0.77],lw=1.5,label='Cost')
	ax.plot(-dRN,dE,'k--',color=[0,1,0],lw=1.5,label='Low harvest (w/o substitutions)')

	ax.set(xticks=np.arange(-12000,18000,1000),xlabel='Cumulative cost (Million CAD)',ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[-8000,700],ylim=[-120,10])
	ax.legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.66,1) bbox_to_anchor=(0.3,0.9,0,0),
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
	ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	plt.tight_layout()

	flg=0
	if flg==1:
		def onclick(event):
			global ix, iy
			ix,iy = event.xdata, event.ydata
			global xt,yt
			xt.append(ix)
			yt.append(iy)
			return xt,yt

		xt=[]; yt=[]
		for i in range(7):
			cid=fig.canvas.mpl_connect('button_press_event', onclick)

	xt=np.array([-2037., -1978., -1824., -1417.,  -886.,  -484.,  -112.])
	yt=np.array([-13., -25., -45., -71., -89., -96., -99.])
	for i in range(len(xt)):
		ax.text(xt[i],yt[i],int(cpt[i]),fontsize=fs2,color=[0,0,0])

	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\MitigationValue_CompletedAndCAP_Discounted','png',900);

	return

#%%
def DefineStrata(meta,pNam,dmec,lsat,status_pt,status_t,status_s):

	#--------------------------------------------------------------------------
	# By project type (funding source code)
	# This uses the unique FSCs: ufcs.NMC_AIL_TS_ByFSC(meta)
	#--------------------------------------------------------------------------
	if status_pt=='On':
		uFSC=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\ail.pkl')['FSC Unique']
		meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=np.append('All',uFSC)
		meta[pNam]['Project']['Strata']['Project Type']['Unique ID']=np.append(0,np.zeros(uFSC.size))
		for i in range(uFSC.size):
			meta[pNam]['Project']['Strata']['Project Type']['Unique ID'][i+1]=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'][uFSC[i]]

	#--------------------------------------------------------------------------
	# By time (last year of implementation)
	#--------------------------------------------------------------------------
	if status_t=='On':

		# Define years to include
		t=np.arange(meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']+1,1)

		meta[pNam]['Project']['Strata']['Year']['Unique ID']=np.append(0,t)
		meta[pNam]['Project']['Strata']['Year']['Unique CD']=np.append('All',t.astype(str))
		iScn=1
		for iStand in range(meta[pNam]['Project']['N Stand']):
			ind=np.where( (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Nutrient App Aerial']) )[0]
			if ind.size>0:
				Year0=np.floor(dmec[iScn][iStand]['Year'][ind[-1]])
				iT=np.where(t==Year0)[0]
				if iT.size==0:
					continue
				meta[pNam]['Project']['Strata']['Year']['ID'][iStand]=t[iT]

	#--------------------------------------------------------------------------
	# By space (BGC zone)
	#--------------------------------------------------------------------------
	if status_s=='On':
		# Define zones to include
		u=np.unique(lsat['ID_BGCZ'])
		cd=np.array([])
		for iU in range(u.size):
			 cd=np.append(cd,u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU]))

		cd=np.array(['CWH','SBS','ICH','ESSF'])
		id=np.zeros(cd.size)
		for j in range(id.size):
			 id[j]=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][cd[j]]

		meta[pNam]['Project']['Strata']['Spatial']['Unique ID']=np.append(0,id)
		meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=np.append('All',cd)
		meta[pNam]['Project']['Strata']['Spatial']['ID']=lsat['ID_BGCZ']

	return meta

#%%
def TabulateSpeciesComposition(meta,pNam,t0,t1):
	zFE=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_YearLast.tif')['Data']
	ind1=np.where( (zFE>=t0) & (zFE<=t1) )
	sph=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\VRI_LIVE_STEMS_PER_HA.tif')['Data'][ind1]
	d=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'].copy()
	for k in d.keys():
		d[k]=0
	for iS in range(6):
		spc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_' + str(iS+1) + '.tif')['Data'][ind1]
		frac=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_' + str(iS+1) + '.tif')['Data'][ind1].astype('float')/100
		uS=np.unique(spc[spc>0])
		for iU in range(uS.size):
			ind2=np.where( (spc==uS[iU]) )
			cd=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],uS[iU])[0]
			d[cd]=d[cd]+np.sum(sph[ind2]*frac[ind2])
	# Save
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummarySpeciesComposition_' + str(t0) + 'to' + str(t1) + '.pkl',d)
	return

#%%
def Plot_AreaTreated_FrequencyByBGC(meta,pNam,t0,t1):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')['Data']
	d=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].copy()
	for k in d.keys():
		d[k]=0
	for i in range(10):
		try:
			zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(i+1) + '_Year.tif')['Data']
		except:
			continue

		ikp=np.where( (zYr>=t0) & (zYr<=t1) )
		if ikp[0].size==0:
			continue

		y=zBGC[ikp]
		u=np.unique(y)
		u=u[u>0]
		for iU in range(u.size):
			ind1=np.where( (y==u[iU]) )[0]
			cd=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			d[cd]=d[cd]+ind1.size
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummaryTreatmentAreaByBGCZone_' + str(t0) + 'to' + str(t1) + '.pkl',d)

	# Plot
	#d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummaryTreatmentAreaByBGCZone_' + str(t0) + 'to' + str(t1) + '.pkl')
	cd=np.array(list(d.keys()))
	A=np.array(list(d.values()))/1000
	ikp=np.where(A>0.1)[0]
	cd=cd[ikp]
	A=A[ikp]
	ord=np.flip(np.argsort(A))
	A=A[ord]
	cd=cd[ord]
	ind=np.where(A>0)[0]
	A=A[ind]
	cd=cd[ind]
	bin=np.arange(1,A.size+1)
	A_max=np.max(A)

	yo=1.05
	yl=[0,1.1*A_max]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10.5));
	plt.bar(bin,A,0.85)
	for i in range(bin.size):
		y=A[i]
		ax.text(bin[i],y*yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Implementation level (Kha)',xticklabels=cd,xlim=[0.25,bin.size+0.75],ylim=yl) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyByBGC_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def AreaTreated_FrequencyByLCC(meta,pNam,t0,t1):
	zNM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_YearLast.tif')['Data']
	zLCC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif')['Data']
	d=meta['LUT']['Derived']['lc_comp1'].copy()
	for k in d.keys():
		ind=np.where( (zNM>=t0) & (zNM<=t1) & (zLCC==meta['LUT']['Derived']['lc_comp1'][k]) )
		d[k]=ind[0].size/1e3

	# Plot
	cd=np.array(list(d.keys()))
	A=np.array(list(d.values()))
	#ikp=np.where(A>0.1)[0]
	#cd=cd[ikp]
	#A=A[ikp]
	ord=np.flip(np.argsort(A))
	A=A[ord]
	cd=cd[ord]
	ind=np.where(A>0)[0]
	A=A[ind]
	cd=cd[ind]
	bin=np.arange(1,A.size+1)
	A_max=np.max(A)

	yo=1.05
	yl=[0,1.1*A_max]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10.5));
	plt.bar(bin,A,0.85)
	for i in range(bin.size):
		y=A[i]
		ax.text(bin[i],y*yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Implementation level (Kha)',xticklabels=cd,xlim=[0.25,bin.size+0.75],ylim=yl) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyByLCC_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def Plot_EmissionsCumu_TimeSeriesFromCurrentYear(meta,mos,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	iT2=np.where( (tv>=2015) & (tv<=2130) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	cl=[0,0,0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Without substitutions
	y1=mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	y1[iT1]=0
	plt.plot(tv[iT2],np.cumsum(y1[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='High harvest')

	y2=mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	y2[iT1]=0
	plt.plot(tv[iT2],np.cumsum(y2[iT2]),'--',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='Low harvest')

	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	#plt.plot(tv[iT2],np.cumsum(y[iT2]),'--',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75)

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(2016,-2.4,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(2016,1.4,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2100,5),yticks=np.arange(-20,20,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlabel='Time, years',ylim=[-3,2],xlim=[2015,2075])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumu_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return

#%%
def Plot_EmissionsAnnualPerHectare_TimeSeriesFromCurrentYear(meta,mos,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	iT2=np.where( (tv>=2010) & (tv<=2100) )[0]

	cl=[0,0,0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])

	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	# Without substitutions
	y1=mos[pNam]['Delta']['Low Harvest']['ByStrata']['Mean']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]
	y1[iT1]=0
	y2=mos[pNam]['Delta']['High Harvest']['ByStrata']['Mean']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]
	y2[iT1]=0

	plt.plot(tv[iT2],y1[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='High harvest')
	plt.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='Low harvest')

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(2012,-13,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(2012,13,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2010,2200,10),yticks=np.arange(-2000,2000,2),ylabel='Annual GHG impact (tCO$_2$e ha$^{-1}$)',xlabel='Time, years',ylim=[-15,15],xlim=[2010,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnualPerHectare_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return


#%%
def Plot_EmissionsCumuPerHectare_TimeSeriesFromCurrentYear(meta,mos,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	iT2=np.where( (tv>=2010) & (tv<=2100) )[0]

	cl=[0,0,0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])

	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	# Half substitutions
	#y=((mos[pNam]['Delta'][cNam]['ByStrata']['Mean']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,iPS,0,iYS])/2)

	# Without substitutions
	y=mos[pNam]['Delta']['Low Harvest']['ByStrata']['Mean']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]
	y[iT1]=0
	plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='High harvest')

	y=mos[pNam]['Delta']['High Harvest']['ByStrata']['Mean']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]
	y[iT1]=0
	plt.plot(tv[iT2],np.cumsum(y[iT2]),'--',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='Low harvest')

	ax.legend(loc='upper left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(2097,-30,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',ha='right',va='center',color=[0.8,0.8,0.8])
	ax.text(2097,30,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',ha='right',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2010,2200,10),yticks=np.arange(-2000,2000,10),ylabel='Cumulative GHG impact (tCO$_2$e ha$^{-1}$)',xlabel='Time, years',ylim=[-100,50],xlim=[2010,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumuPerHectare_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return

#%%
def Plot_EmissionsAnnual_TimeSeriesFromCurrentYear(meta,mos,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	iT2=np.where( (tv>=2020) & (tv<=2130) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	cl=[0,0,0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Half substitutions
	#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,iPS,0,iYS])/2)/1e6

	# Without substitutions
	y1=mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	y1[iT1]=0
	plt.plot(tv[iT2],y1[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='High harvest')

	y2=mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	y2[iT1]=0
	plt.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label='Low harvest')

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(2022,-0.5,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(2022,0.5,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2100,5),yticks=np.arange(-20,20,0.5),ylabel='Annual GHG impact (MtCO$_2$e yr${^-1}$)',xlabel='Time, years',ylim=[-1,1],xlim=[2020,2075])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnual_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%

	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeriesCompleted(meta,pNam,mos,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=0.5,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Without substitutions
	y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],y[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1966,-2.2,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1966,2.2,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.5),xlabel='Time, years',ylim=[-2.5,2.5],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ax2.plot(tv[iT2],y[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,20),ylabel='Cumulative GHG impact (MtCO$_2$e)',ylim=[-100,100],xlim=[1965,2100])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_Completed','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeriesCurrentYear(meta,pNam,mos,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=0.5,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	y1=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	y1[iT1]=0
	ax.plot(tv[iT2],y1[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.text(1970,-0.4,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1970,0.4,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',ylim=[-0.5,0.5],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y2=np.cumsum(y1)
	ax2.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',ylim=[-4,4],xlim=[1965,2100])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_CurrentYear','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP(meta,pNam,mos,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],y[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Completed')

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1966,-2.2,'Net removals',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(1966,2.2,'Net emissions',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.5),xlabel='Time, years',ylim=[-2.5,2.5],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1)
	
	ax2=ax.twinx()
	y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ax2.plot(tv[iT2],y[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,20),ylabel='Cumulative GHG impact (MtCO$_2$e)',ylim=[-100,100],xlim=[1965,2100])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP','png',900)
	return

#%%
def Plot_YieldAnnAndCumu_TimeSeriesCompleted(meta,pNam,mos,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=0.5,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],y[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	ax.text(1966,-0.7,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1966,0.7,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.2),xlabel='Time, years',ylim=[-0.8,0.8],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ax2.plot(tv[iT2],y[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,10),ylabel='',ylim=[-40,40],xlim=[1965,2100])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\YieldAnnAndCumu_TimeSeries_Completed_' + cNam,'png',900)
	return

#%%
def Plot_YieldAnnAndCumu_TimeSeriesCompletedAndCAP(meta,pNam,mos,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=0.5,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],y[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1966,-0.7,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1966,0.7,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.2),xlabel='Time, years',ylim=[-0.8,0.8],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ax2.plot(tv[iT2],y[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,5),ylabel='',ylim=[-40,40],xlim=[1965,2100])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\YieldAnnAndCumu_TimeSeries_CompletedAndCAP_' + cNam,'png',900)
	return

#%%
def QA_Plot_CarbonTimeSeriesByStand(meta,pNam,dmec):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	for iBat in range(meta[pNam]['Project']['N Batch']):
		indBat=cbu.IndexToBatch(meta[pNam],iBat)
		d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,iBat)
		d1=cbu.LoadSingleOutputFile(meta,pNam,1,0,iBat)

		for i in range(0,indBat.size,25):
			plt.close('all'); fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(18,10));
			ax[0,0].plot(tv,d0['C_Biomass'][:,i],'b-',label='Baseline')
			ax[0,0].plot(tv,d1['C_Biomass'][:,i],'g--',label='Actual')
			ax[0,0].set(ylabel='Biomass')
			ax[0,0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
			ax[0,1].plot(tv,d0['C_DeadWood'][:,i],'b-')
			ax[0,1].plot(tv,d1['C_DeadWood'][:,i],'g--')
			ax[0,1].set(ylabel='Dead wood')
			ax[1,0].plot(tv,d0['E_Domestic_ForestSector_Wildfire'][:,i],'ob-',ms=3)
			ax[1,0].plot(tv,d1['E_Domestic_ForestSector_Wildfire'][:,i],'sg--',ms=3)
			ax[1,0].set(ylabel='Wildfire emissions')
			ax[1,1].plot(tv,d0['V_ToMill_MerchTotal'][:,i],'ob-')
			ax[1,1].plot(tv,d1['V_ToMill_MerchTotal'][:,i],'sg--')
			ax[1,1].set(ylabel='Volume to mill, merch live')
			ax[2,0].plot(tv,d0['E_Domestic_ForestSector_OpenBurning'][:,i],'b-')
			ax[2,0].plot(tv,d1['E_Domestic_ForestSector_OpenBurning'][:,i],'g--')
			ax[2,0].set(ylabel='Open burning\nemissions')
			#ax[2,1].plot(tv,d0['E_EnergySC_Bioenergy'][:,i],'b-')
			#ax[2,1].plot(tv,d1['E_EnergySC_Bioenergy'][:,i],'g--')
			#ax[2,1].set(ylabel='Bioenergy emissions')
			ax[2,1].plot(tv,d0['C_Litter'][:,i]+d0['C_Soil'][:,i],'b-')
			ax[2,1].plot(tv,d1['C_Litter'][:,i]+d1['C_Soil'][:,i],'g--')
			ax[2,1].set(ylabel='Litter+Soil')
			ax[3,0].plot(tv,d0['E_Atmosphere_SubstitutionExcluded'][:,i],'b-')
			ax[3,0].plot(tv,d1['E_Atmosphere_SubstitutionExcluded'][:,i],'g--')
			ax[3,0].set(ylabel='Annual GHG balance')
			ax[3,1].plot(tv,d0['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,i],'b-')
			ax[3,1].plot(tv,d1['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,i],'g--')
			ax[3,1].set(ylabel='Cumulative GHG balance')

			plt.tight_layout()
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA\\' + str(indBat[i]),'png',200)
	return

#%%
def ImportCombined():

	# Import Completed
	pNamC='BCFCS_NMC'
	pth=r'C:\Data\BCFCS\BCFCS_NMC\Inputs\Metadata.pkl'
	meta=gu.ipickle(pth)
	tv=np.arange(meta[pNamC]['Project']['Year Start Saving'],meta[pNamC]['Project']['Year End']+1,1)
	mos=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNamC)
	cNam='NM1'
	mos[pNamC]['Delta']={}
	mos[pNamC]['Delta'][cNam]={'iB':0,'iP':1}
	mos=cbu.Import_MOS_ByScnComparisonAndStrata(meta,pNamC,mos)
	mos=cbu.Import_MOS_ByScnAndStrata_Area(meta,pNamC,mos)
	#thlb=gu.ipickle(meta['Paths'][pNamC]['Data'] + '\\Inputs\\thlb.pkl')
	#lsat=gu.ipickle(meta['Paths'][pNamC]['Data'] + '\\Inputs\\lsat.pkl')
	#dmec=gu.ipickle(meta['Paths'][pNamC]['Data'] + '\\Inputs\\dmec.pkl')
	meta['Graphics']['Print Figures']='On'
	meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS_NM'
	nPS='All'
	nSS='All'
	nYS='All'
	iPS=np.where(meta[pNamC]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
	iSS=np.where(meta[pNamC]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
	iYS=np.where(meta[pNamC]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
	iT=np.where( (tv>=1850) & (tv<=2150) )[0]
	
	#Import future
	pNamF='BCFCS_NMF'
	pth=r'C:\Data\BCFCS\BCFCS_NMF\Inputs\Metadata.pkl'
	meta_tmp=gu.ipickle(pth)
	meta[pNamF]=meta_tmp[pNamF]
	meta['Paths'][pNamF]=meta_tmp['Paths'][pNamF]
	del meta_tmp
	mos_tmp=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNamF)
	mos[pNamF]=mos_tmp[pNamF]
	del mos_tmp
	mos[pNamF]['Delta']={}
	mos[pNamF]['Delta'][cNam]={'iB':0,'iP':1}
	mos=cbu.Import_MOS_ByScnComparisonAndStrata(meta,pNamF,mos)
	mos=cbu.Import_MOS_ByScnAndStrata_Area(meta,pNamF,mos)

	# Delete historical effects from future project
	flg=0
	if flg==1:
		iT2=np.where(tv<=2022)[0]
		for k1 in mos[pNamF]['Delta'][cNam]['ByStrata']['Sum'].keys():
			for k2 in mos[pNamF]['Delta'][cNam]['ByStrata']['Sum'][k1].keys():
				mos[pNamF]['Delta'][cNam]['ByStrata']['Sum'][k1][k2][iT2,:,:,:]=0
				mos[pNamF]['Delta'][cNam]['ByStrata']['Mean'][k1][k2][iT2,:,:,:]=0
		for k1 in mos[pNamF]['Scenarios'][1]['Sum'].keys():
			for k2 in mos[pNamF]['Scenarios'][1]['Sum'][k1].keys():
				mos[pNamF]['Scenarios'][0]['Sum'][k1][k2][iT2,:,:,:]=0
				mos[pNamF]['Scenarios'][0]['Mean'][k1][k2][iT2,:,:,:]=0
				mos[pNamF]['Scenarios'][1]['Sum'][k1][k2][iT2,:,:,:]=0
				mos[pNamF]['Scenarios'][1]['Mean'][k1][k2][iT2,:,:,:]=0

	# Combine completed and future
	pNamCF='BCFCS_NMCF'
	mos=cbu.CombineProjectMOSs(meta,mos,pNamCF,['BCFCS_NMC','BCFCS_NMF'])
	meta[pNamCF]=meta[pNamC]

	meta[pNamCF]['Project']['Multi']=1e6

	return pNamC,pNamF,meta,mos,iPS,iSS,iYS,iT,cNam

#%%
def AreaImpactedByDisturbance(meta,pNam,t0,t1):
	zNM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_YearLast.tif')['Data']
	zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')['Data']
	zWF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\WildfireComp1_YearLast.tif')['Data']
	
	t=np.arange(1960,meta[pNam]['YearCurrent']+1)
	d={'H':np.zeros(t.size),'WF':np.zeros(t.size)}
	for iT in range(t.size):
		print(t[iT])
		ind=np.where( (zNM>=t0) & (zNM<=t1) & (zH>zNM) & (zH==t[iT]) )
		d['H'][iT]=ind[0].size
		ind=np.where( (zNM>=t0) & (zNM<=t1) & (zWF==t[iT]) )
		d['WF'][iT]=ind[0].size
	
	plt.close('all'); plt.bar(t,d['H'],0.75,facecolor='b')
	plt.close('all'); plt.bar(t,d['WF'],0.75,facecolor='b')
	return
#%%
def AgeAtTimeOfApplication(meta,pNam,t0,t1):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zA=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2015\\PROJ_AGE_1.tif')['Data']
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')['Data']

	bw=10; bin=np.arange(bw,200+bw,bw)

	d=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].copy()
	for k in d.keys():
		d[k]=np.zeros(bin.size)

	for i in range(10):
		try:
			zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\FE-CA_' + str(i+1) + '_Year.tif')['Data']
		except:
			continue

		ikp=np.where( (zYr>=t0) & (zYr<=t1) )
		if ikp[0].size==0:
			continue

		A=zA[ikp]
		Z=zBGC[ikp]
		u=np.unique(Z)
		u=u[u>0]
		for iU in range(u.size):
			cd=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			for iA in range(bin.size):
				ind1=np.where( (Z==u[iU]) & (np.abs(A-bin[iA])<bw/2) )[0]
				d[cd][iA]=d[cd][iA]+ind1.size
	#gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummaryTreatmentAreaByBGCZone_' + str(t0) + 'to' + str(t1) + '.pkl',d)

	# Plot
	for zn in d.keys():
		if np.sum(d[zn])<10:
			continue
		cd=np.array(list(d.keys()))
		A=d[zn]/1000
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10.5))
		plt.bar(np.arange(bin.size),A,0.75,facecolor=[0.25,0.5,1])
		for i in range(bin.size):
			ax.text(i,A[i]+0.1,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
		ax.set(xticks=np.arange(bin.size),xticklabels=bin,ylabel='Frequency (ha)',xlim=[-0.75,bin.size-0.25],ylim=[0,np.max(A)+1])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyByAge_' + '_' + zn + '_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def ProjectFuture(meta,pNam,mos):
	# This function takes the current year per-hectare time series and applies it into the future at a projected AIL
	# Notes: Cumulative values are incorrect - do not use!

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]

	dCP=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_CostsAndPrices.xlsx','Data')

	mosWF=copy.deepcopy(mos)
	for cNam in mos[pNam]['Delta'].keys():
		for v in mos[pNam]['Delta'][cNam]['ByStrata']['Mean'].keys():
			#d[k]=np.zeros(tv.size)
			y0=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,:,:,iYS]
			tv0=tv[iT]
			for yr in range(meta[pNam]['YearCurrent']+1,tv[-1],1):
				iT1=np.where( (tv>=yr) )[0]
				iT2=np.where( (tv0>=meta[pNam]['YearCurrent']) )[0]
				iT2=iT2[0:iT1.size]
				#mosWF[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT1,:,:,0]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT1,:,:,0]+y0[iT2,:,:]
				mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,:,:,0]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,:,:,0]+meta[pNam]['AIL CAP']*y0[iT2,:,:]

		# Fix costs
		iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
		iTe=np.where(dCP['Year']==meta[pNam]['YearCurrent'])[0]
		y=meta[pNam]['AIL CAP']*(dCP['Cost Nutrient Purchase (CAD/ha)'][iTe]+dCP['Cost Nurtrient Application (CAD/ha)'][iTe]+dCP['Cost Nutrient Overhead (CAD/ha)'][iTe])
		mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Total']['Ensemble Mean'][iT,:,:,:]+y
		mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Nutrient Management']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Cost Nutrient Management']['Ensemble Mean'][iT,:,:,:]+y
		mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum']['Revenue Net']['Ensemble Mean'][iT,:,:,:]-y

		# The cumulative values are bogus
		#vSTZ=[]
		#for v in vSTZ:
		#	mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean']
	return mosWF

#%%
def ProjectServicePlan(meta,pNam,mos):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	YearProjection=meta[pNam]['YearCurrent']+1
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
	mosSPP=copy.deepcopy(mos)
	for cNam in mos[pNam]['Delta'].keys():
		for v in mos[pNam]['Delta'][cNam]['ByStrata']['Mean'].keys():
			mosSPP[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,:,:,:]=0
			y0=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,:,:,iYS]
			tv0=tv[iT]
			iT1=np.where( (tv>=YearProjection) )[0]
			iT2=np.where( (tv0>=meta[pNam]['YearCurrent']) )[0]
			iT2=iT2[0:iT1.size]
			mosSPP[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,:,:,0]=meta[pNam]['AIL CAP']*y0[iT2,:,:]
	return mosSPP

#%%
def Plot_PerHectareEmissionsByBGCZone(meta,pNam,d):
	cd=list(d['PerHa']['Low Harvest'].keys())
	bw=0.3
	bin=np.arange(1,len(cd)+1)
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,7));
	y=list(d['PerHa']['Low Harvest'].values())
	plt.bar(bin-bw/1.9,y,bw,facecolor=[0,0.5,0],label='Low harvest')
	y=list(d['PerHa']['High Harvest'].values())
	plt.bar(bin+bw/1.9,y,bw,facecolor=[0.6,0.9,0],label='High harvest')
	ax.legend(loc='lower left',labelspacing=0.5,facecolor=[1,1,1],frameon=False,ncols=1,fontsize=6);
	ax.set(xticks=bin,xticklabels=cd,ylabel='Cumulative GHG impact 2050 (tCO$_2$e ha${^-1}$)',xlim=[0.5,bin.size+0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\PerHectareEmissionsByBGCZone','png',200)
	return

#%%
def SummarizeYieldCumulative(meta,pNam,mos):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	cl1=[0.27,0.49,0.77]; cl2=[0.6,0.9,0]; lw=1; xl=[1950,2150]
	plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,10))
	
	v='V_ToMill_MerchTotal'
	ax[0,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,0].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw,label='Low harvest')
	ax[0,0].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw,label='High harvest')
	ax[0,0].set(ylabel='Removals, merch (Mm$^3$)',xlim=xl)
	ax[0,0].legend(loc='upper center',facecolor=[1,1,1],frameon=False,ncols=1,fontsize=8)

	v='V_ToMill_NonMerchTotal'
	ax[0,1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,1].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[0,1].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[0,1].set(ylabel='Removals, non-merch (Mm$^3$)',xlim=xl)
	
	#v='Cost Silviculture Total'
	v='Cost Nutrient Management'
	ax[0,2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,2].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k-',color=cl1,lw=lw)
	ax[0,2].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k--',color=cl2,lw=lw)
	ax[0,2].set(ylabel='Cost, nutrient app (Billion CAD)',xlim=xl)
	
	v='Cost Total'
	ax[1,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,0].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k-',color=cl1,lw=lw)
	ax[1,0].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k--',color=cl2,lw=lw)
	ax[1,0].set(ylabel='Cost, total (Billion CAD)',xlim=xl)
	
	v='Revenue Gross'
	ax[1,1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,1].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k-',color=cl1,lw=lw)
	ax[1,1].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k--',color=cl2,lw=lw)
	ax[1,1].set(ylabel='Revenue, gross (Billion CAD)',xlim=xl)
	
	v='Revenue Net'
	ax[1,2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,2].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k-',color=cl1,lw=lw)
	ax[1,2].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']/1e3),'k--',color=cl2,lw=lw)
	ax[1,2].set(ylabel='Revenue, net (Billion CAD)',xlim=xl)

	for i in range(2):
		for j in range(3):
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	return

#%%
def SummarizeFluxesCumulative(meta,pNam,mos):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	cl1=[0.27,0.49,0.77]; cl2=[0.6,0.9,0]; lw=1; xl=[1950,2150]
	plt.close('all'); fig,ax=plt.subplots(3,4,figsize=gu.cm2inch(22,10))
	
	v='E_Domestic_ForestSector_NPP'
	ax[0,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,0].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw,label='Low harvest')
	ax[0,0].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw,label='High harvest')
	ax[0,0].set(ylabel='-1 x NPP (MtCO$_2$e)',xlim=xl)
	ax[0,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False,ncols=1,fontsize=8)
	
	v='E_Domestic_ForestSector_RH'
	ax[0,1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,1].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[0,1].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[0,1].set(ylabel='RH (MtCO$_2$e)',xlim=xl)
	
	v='E_Domestic_ForestSector_NEE'
	ax[0,2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,2].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[0,2].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[0,2].set(ylabel='NEE (MtCO$_2$e)',xlim=xl)
	
	v='E_Domestic_ForestSector_Wildfire'
	ax[0,3].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[0,3].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[0,3].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[0,3].set(ylabel='Wildfire E (MtCO$_2$e)',xlim=xl)
	
	v='E_Domestic_ForestSector_OpenBurning'
	ax[1,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,0].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[1,0].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[1,0].set(ylabel='Pile burning (MtCO$_2$e)',xlim=xl)
	
	v='E_ForestOperations'
	ax[1,1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,1].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[1,1].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[1,1].set(ylabel='Operations (MtCO$_2$e)',xlim=xl)
	
	v='E_Bioenergy'
	ax[1,2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,2].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[1,2].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[1,2].set(ylabel='Bioenergy (MtCO$_2$e)',xlim=xl)
	
	v='E_Domestic_ForestSector_HWP'
	ax[1,3].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,3].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[1,3].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[1,3].set(ylabel='Product decay (MtCO$_2$e)',xlim=xl)
	
	v='E_Substitution_Material'
	ax[2,0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[2,0].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[2,0].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[2,0].set(ylabel='Material subst (MtCO$_2$e)',xlim=xl)
	
	v='E_Substitution_Energy'
	ax[2,1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[2,1].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[2,1].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[2,1].set(ylabel='Energy subst (MtCO$_2$e)',xlim=xl)
	
	v='E_Atmosphere_SubstitutionExcluded'
	ax[2,2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[2,2].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[2,2].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[2,2].set(ylabel='Net E w/o subst (MtCO$_2$e)',xlim=xl)
	
	v='E_Atmosphere_SubstitutionIncludedHalf'
	ax[2,3].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax[2,3].plot(tv,np.cumsum(mos[pNam]['Delta']['Low Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k-',color=cl1,lw=lw)
	ax[2,3].plot(tv,np.cumsum(mos[pNam]['Delta']['High Harvest']['ByStrata']['Sum'][v]['Ensemble Mean'][:,0,0,0]/meta[pNam]['Project']['Multi']),'k--',color=cl2,lw=lw)
	ax[2,3].set(ylabel='Net E with subst (MtCO$_2$e)',xlim=xl)
	
	for i in range(3):
		for j in range(4):
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.04,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	return

#%%
def NA_CalcNUE(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	iPS=0; iSS=0; iYS=0
	doseN=200 # NUE applied
	df=pd.DataFrame()
	Names=[]
	for sc in mos[pNam]['Delta'].keys():
		if 'cNam' in kwargs.keys():
			if np.isin(sc,kwargs['cNam'])==False:
				continue
		Names.append(sc)
		v='C_AboveGroundBiomass'
		Temp=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Biomass']['Ensemble Mean']-mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean']
		d={}
		v='C_AboveGroundBiomass'
		y=1000*np.mean(Temp)/doseN
		d[v[2:]]=np.round(y,decimals=2)
		v='C_Root'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v[2:]]=np.round(y,decimals=2)
		v='C_DeadWood'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v[2:]]=np.round(y,decimals=2)
		v='C_Litter'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v[2:]]=np.round(y,decimals=2)
		v='C_Soil'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v[2:]]=np.round(y,decimals=2)
		v='C_Forest'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v[2:]]=np.round(y,decimals=2)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)
	#df.index.name='Variable'
	df.columns=Names
	#df=df.sort_index(axis=0)
	if 'save' in kwargs.keys():
		df.to_excel(meta['Paths'][pNam]['Data'] + '\\Outputs\\TabularSummaryDelta_' + kwargs['table_name'] + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']) + '.xlsx')
	return df

#%%
def CumulativelativeGHGBenefit_WithCI(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]
	iPS=0; iSS=0; iYS=0

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,8)); lw=0.75; cl1=[0.96,0.92,1]; cl2=[0.9,0.82,0.96]
	# Without harvesting
	sc='Coast No Harvest'
	#sc='Interior No Harvest'
	v='E_Atmosphere_SubstitutionExcluded'
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[0,0].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='95% percentile range')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
		ax[0,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1,label='3 x S.E. range')
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
		ax[0,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1,label='50% percentile range')
	mu=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[0,0].set(yticks=np.arange(-40,50,10),ylabel='$\Delta$E\n(tCO$_2$e ha$^-$$^1$ YearProjection$^-$$^1$)',ylim=[-35,35],xlim=[ tv[iT[0]],tv[iT[-1]] ])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,0].legend(loc=kwargs['LegendPosition'],frameon=False,facecolor=None)
	ax[0,0].annotate('Application',xy=(2023,5),xytext=(2023,24),arrowprops=dict(arrowstyle="->"),ha='center')
	ax[0,0].text(2045,20,'Sink',ha='center',va='center',fontsize=8,style='italic',weight='bold',color=[0.5,0.5,0.5])
	ax[0,0].text(2045,-20,'Source',ha='center',va='center',fontsize=8,style='italic',weight='bold',color=[0.5,0.5,0.5])
	#ax[0].annotate('Manufacture, transport & N2O emissions',xy=(2020,-6),xytext=(2020,-16),arrowprops=dict(arrowstyle="->"),ha='center')
	
	v='E_Atmosphere_SubstitutionExcluded_Cumulative_from_tref'
	ax[0,1].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[0,1].yaxis.set_ticks_position('both');ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].set(yticks=np.arange(-260,160,20),ylabel='Cumuulative $\Delta$E\n(tCO$_2$e ha$^-$$^1$)',ylim=[-120,20],xlim=[ tv[iT[0]],tv[iT[-1]] ]);

	#------------------------------------------------------------------------------
	# With harvesting
	#------------------------------------------------------------------------------
	sc='Coast With Harvest'
	#sc='Interior With Harvest'
	v='E_Atmosphere_SubstitutionExcluded'
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[1,0].set(yticks=np.arange(-40,50,10),ylabel='$\Delta$E\n(tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',ylim=[-35,35],xlim=[2010,2080]);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].annotate('Post-application\nharvest',xy=(2047,7),xytext=(2047,22),arrowprops=dict(arrowstyle="->"),ha='center')

	ax[1,1].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	v='E_Atmosphere_SubstitutionExcluded_Cumulative_from_tref'
	lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].set(yticks=np.arange(-260,160,20),ylabel='Cumulative $\Delta$E\n(tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',ylim=[-120,20],xlim=[2010,2080]);
	
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\CumulativeGHGBenefit_WithCI','png',900)
	return