'''
UTILITIES - NON-OBLIGATION STAND ESTABLISHMENT (NOSE)
'''

#%% Import Python modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import warnings
import copy
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
	pNam='BCFCS_NOSEC'
	meta=gu.ipickle(r'C:\Data\BCFCS\BCFCS_NOSEC\Inputs\Metadata.pkl')
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	mos=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNam)
	cNam='NOSE1'
	mos[pNam]['Delta']={}
	mos[pNam]['Delta'][cNam]={'iB':0,'iP':1}
	mos=cbu.Import_MOS_ByScnComparisonAndStrata(meta,pNam,mos)
	mos=cbu.Import_MOS_ByScnAndStrata_Area(meta,pNam,mos)
	dmec=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\dmec.pkl')
	#thlb=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\thlb.pkl')
	nPS='All'
	nSS='All'
	nYS='All'
	iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
	meta[pNam]['Project']['Multi']=1e6
	return pNam,meta,tv,mos,dmec,cNam,iPS,iSS,iYS

#%%
def ImportModelResults():
	pNam='BCFCS_NOSEC'
	meta=gu.ipickle(r'C:\Data\BCFCS\BCFCS_NOSEC\Inputs\Metadata.pkl')
	meta[pNam]['Project']['Multi']=1e6
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	mos=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNam)
	cNam='NOSE1'
	mos[pNam]['Delta']={}
	mos[pNam]['Delta'][cNam]={'iB':0,'iP':1}
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
	return pNam,meta,tv,mos,dmec

#%% Save summary of DMEC and results to spreadsheet for troubleshooting, investigation, QA
def QA_PrintSummarySparseToSpreadsheet(meta,iScn,dmec,E):
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
	d['ASET Focal']=np.array(['' for _ in range(n)],dtype=object)
	d['ASET']=np.array(['' for _ in range(n)],dtype=object)
	d['Inciting NOSE']=np.array(['' for _ in range(n)],dtype=object)
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

		ind=np.where( (dmec[iScn][i]['ASET']!=9999) & (dmec[iScn][i]['ASET']!=0) )[0]
		if ind.size>0:
			year_focal=dmec[iScn][i]['Year'][ind[-1]]
			aset_focal=cbu.lut_n2s(meta['LUT']['Derived']['ASET'],dmec[iScn][i]['ASET'][ind[-1]])[0]
			fsc_focal=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][ind[-1]])[0]
		else:
			aset_focal='None'
			fsc_focal='None'

		FlagIncite=np.zeros(dmec[iScn][i]['Year'].size)
		ind=np.where(dmec[iScn][i]['Index to Event Inciting NOSE']!=9999)[0]
		FlagIncite[dmec[iScn][i]['Index to Event Inciting NOSE'][ind]]=1

		for j in range(dmec[iScn][i]['Year'].size):
			d['Index'][c]=i
			d['Year'][c]=dmec[iScn][i]['Year'][j]
			d['DOY'][c]=dmec[iScn][i]['DOY'][j]
			d['Mortality'][c]=dmec[iScn][i]['Mortality Factor'][j]
			d['Event Type'][c]=cbu.lut_n2s(meta['LUT']['Event'],dmec[iScn][i]['ID Event Type'][j])[0]
			if dmec[iScn][i]['Event Source'][j]==2:
				d['Event Source'][c]='Added'
			if (dmec[iScn][i]['ASET'][j]!=9999) & (dmec[iScn][i]['ASET'][j]!=0):
				d['ASET'][c]=cbu.lut_n2s(meta['LUT']['Derived']['ASET'],dmec[iScn][i]['ASET'][j])[0]
			if dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][j]!=9999:
				d['FSC'][c]=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],dmec[iScn][i]['SILV_FUND_SOURCE_CODE'][j])[0]
			d['Year Focal'][c]=int(year_focal)
			d['ASET Focal'][c]=aset_focal
			d['FSC Focal'][c]=fsc_focal
			if FlagIncite[j]==1:
				d['Inciting NOSE'][c]='Inciting'
			d['OPENING'][c]=zOPID[i]
			d['BGCZ'][c]=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],zBGCZ[i])[0]
			d['E_Atmosphere_SubstitutionExcluded Delta'][c]=int(E['Delta']['E_Atmosphere_SubstitutionExcluded'][i])
			if (dmec[iScn][i]['ASET'][j]!=9999) & (dmec[iScn][i]['ASET'][j]!=0):
				d['Pl Spc1 CD'][c]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],dmec[iScn][i]['PL_SPECIES_CD1'][j])[0]
				d['Pl Spc1 %'][c]=dmec[iScn][i]['PL_SPECIES_PCT1'][j]
				d['Pl SPH'][c]=dmec[iScn][i]['Planted SPH'][j]
			c=c+1

	#df=pd.DataFrame.from_dict(d)
	pthout=r'C:\Data\BCFCS\BCFCS_NOSEC\Inputs\SummarySparseSample.xlsx'
	gu.PrintDict(d,pthout,SheetName='Sheet1')

	return

#%%
def Tabulate_CurrentYear_ForServicePlan(meta,pNam,mos):

	cd=np.array(list(meta['LUT']['Derived']['ASET'].keys()))

	# No need for area, just take mean from MOS
	#dI=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ASET_SummaryByTime.pkl')
	#iT=np.where(dI['tv']==meta[pNam]['YearCurrent'])[0]
	#A_Treat=dI['A NOSE'][iT[0],:]
	
	v='E_Atmosphere_SubstitutionExcluded'
	#v='E_Atmosphere_SubstitutionIncluded'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=2015) & (tv<=2050) )[0]
	dE={}
	dEpHa={}
	for i in range(cd.size):
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==cd[i])[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		y=np.sum(mos[pNam]['Delta']['NOSE1']['ByStrata']['Sum'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
		dE[cd[i]]=np.round(y/meta[pNam]['Project']['Multi'],decimals=2)
		y=np.sum(mos[pNam]['Delta']['NOSE1']['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
		dEpHa[cd[i]]=y
	
	#d['NOSE']=d['Salvage and Planting Post Beetle']+d['Salvage and Planting Post Other']+d['Underplanting']
	#d
	return dE,dEpHa

#%%
def QA_Plot_CarbonTimeSeriesByStand(meta,pNam,dmec):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	for iBat in range(meta[pNam]['Project']['N Batch']):
		indBat=cbu.IndexToBatch(meta[pNam],iBat)
		d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,iBat)
		d1=cbu.LoadSingleOutputFile(meta,pNam,1,0,iBat)	
		for k in meta['LUT']['Derived']['ASET'].keys():
			ind=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID'][indBat]==meta['LUT']['Derived']['ASET'][k]) )[0]

			if ind.size==0:
				continue
			#if k!='Knockdown and Planting':
			#	continue

			cnt=0
			for i in range(0,ind.size,int(ind.size/10)):

				if cnt>15:
					# Limit the sample to a total of 25 per project type
					continue

				plt.close('all'); fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(18,10));
				ax[0,0].plot(tv,d0['C_Biomass'][:,ind[i]],'b-',label='Baseline')
				ax[0,0].plot(tv,d1['C_Biomass'][:,ind[i]],'g--',label='Actual')
				ax[0,0].set(ylabel='Biomass')
				ax[0,0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
				ax[0,1].plot(tv,d0['C_DeadWood'][:,ind[i]],'b-')
				ax[0,1].plot(tv,d1['C_DeadWood'][:,ind[i]],'g--')
				ax[0,1].set(ylabel='Dead wood')
				ax[1,0].plot(tv,d0['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'ob-',ms=3)
				ax[1,0].plot(tv,d1['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'sg--',ms=3)
				ax[1,0].set(ylabel='Wildfire emissions')
				ax[1,1].plot(tv,d0['V_ToMill_MerchTotal'][:,ind[i]],'ob-')
				ax[1,1].plot(tv,d1['V_ToMill_MerchTotal'][:,ind[i]],'sg--')
				ax[1,1].set(ylabel='Volume to mill, merch live')
				ax[2,0].plot(tv,d0['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'b-')
				ax[2,0].plot(tv,d1['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'g--')
				ax[2,0].set(ylabel='Open burning\nemissions')
				#ax[2,1].plot(tv,d0['E_EnergySC_Bioenergy'][:,ind[i]],'b-')
				#ax[2,1].plot(tv,d1['E_EnergySC_Bioenergy'][:,ind[i]],'g--')
				#ax[2,1].set(ylabel='Bioenergy emissions')
				ax[2,1].plot(tv,d0['C_Litter'][:,ind[i]]+d0['C_Soil'][:,ind[i]],'b-')
				ax[2,1].plot(tv,d1['C_Litter'][:,ind[i]]+d1['C_Soil'][:,ind[i]],'g--')
				ax[2,1].set(ylabel='Litter+Soil')
				ax[3,0].plot(tv,d0['E_Atmosphere_SubstitutionExcluded'][:,ind[i]],'b-')
				ax[3,0].plot(tv,d1['E_Atmosphere_SubstitutionExcluded'][:,ind[i]],'g--')
				ax[3,0].set(ylabel='Annual GHG balance')
				ax[3,1].plot(tv,d0['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,ind[i]],'b-')
				ax[3,1].plot(tv,d1['E_Atmosphere_SubstitutionExcluded_Cumulative'][:,ind[i]],'g--')
				ax[3,1].set(ylabel='Cumulative GHG balance')

				ind1=np.where(dmec[1][indBat[ind[i]]]['ASET']<9999)[0]
				yr1=dmec[1][indBat[ind[i]]]['Year'][ind1[-1]]
				plt.tight_layout()
				gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA\\' + k + '_' + str(indBat[ind[i]]) + '_' + str(int(yr1)),'png',200)
				cnt=cnt+1
	return

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
			#ax[0,1].plot(tv,d0['V_MerchDead'][:,ind[i]],'b-',label='Baseline')
			#ax[0,1].plot(tv,d1['V_MerchDead'][:,ind[i]],'g--',label='Actual')
			ax[0,1].set(ylabel='Biomass')
			ax[0,1].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)

			ax[1,0].plot(tv,d0['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'ob-',ms=3)
			ax[1,0].plot(tv,d1['E_Domestic_ForestSector_Wildfire'][:,ind[i]],'sg--',ms=3)
			ax[1,0].set(ylabel='Wildfire emissions')

			ax[1,1].plot(tv,d0['V_ToMill_MerchTotal'][:,ind[i]],'ob-')
			ax[1,1].plot(tv,d1['V_ToMill_MerchTotal'][:,ind[i]],'sg--')
			#ax[1,1].plot(tv,d0['V_ToMill_MerchDead'][:,ind[i]],'ob-')
			#ax[1,1].plot(tv,d1['V_ToMill_MerchDead'][:,ind[i]],'sg--')
			#ax[1,1].plot(tv,d0['C_ToMillTotal'][:,ind[i]],'ob-')
			#ax[1,1].plot(tv,d1['C_ToMillTotal'][:,ind[i]],'sg--')
			#ax[1,1].plot(tv,d0['C_ToMillMerchDead'][:,ind[i]],'ob-')
			#ax[1,1].plot(tv,d1['C_ToMillMerchDead'][:,ind[i]],'sg--')
			#ax[1,1].plot(tv,d0['C_Felled'][:,ind[i]],'ob-')
			#ax[1,1].plot(tv,d1['C_Felled'][:,ind[i]],'sg--')
			ax[1,1].set(ylabel='Volume to mill, merch')

			ax[2,0].plot(tv,d0['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'b-')
			ax[2,0].plot(tv,d1['E_Domestic_ForestSector_OpenBurning'][:,ind[i]],'g--')
			ax[2,0].set(ylabel='Open burning emissions')

			#ax[2,1].plot(tv,d0['C_Litter'][:,ind[i]]+d0['C_Soil'][:,ind[i]]+d0['C_DeadWood'][:,ind[i]],'b-')
			#ax[2,1].plot(tv,d1['C_Litter'][:,ind[i]]+d1['C_Soil'][:,ind[i]]+d1['C_DeadWood'][:,ind[i]],'g--')
			ax[2,1].plot(tv,d0['C_DeadWood'][:,ind[i]],'b-')
			ax[2,1].plot(tv,d1['C_DeadWood'][:,ind[i]],'g--')
			#ax[2,1].plot(tv,d0['C_DeadStemMerch'][:,ind[i]],'b-')
			#ax[2,1].plot(tv,d1['C_DeadStemMerch'][:,ind[i]],'g--')
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
def Plot_AreaTreated_TimeSeriesByASET(meta,pNam,type,**kwargs):
	#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\Planting_SummaryByTime.pkl')
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTime.pkl')

	# Summary for Service plan
	iT=np.where(d['tv']==2023)[0]
	cnt=0
	b={}
	for k in meta['LUT']['Derived']['ASET'].keys():
		b[k]=d['A NOSE'][iT[0],cnt]
		cnt=cnt+1

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[22,9]

	flg=0
	if flg==1:
		# Look at total amount of underplanting since 2007
		iT=np.where(d['tv']>=2007)[0]
		meta['LUT']['Derived']['ASET']['Underplanting']
		meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']
		print(np.sum(d['A All'][iT,5:7])/1e3)
		plt.bar(d['tv'],d['A NOSE'][iT,5]/1e3,0.85)
	
	iT=np.where( (d['tv']>=1960) & (d['tv']<=2050) )[0]
	#cl=np.array([[0.7,0.67,0.64],[0.9,0.87,0.84],[0.6,0.75,1],[0,0,0.5],[1,1,0],[0.6,1,0],[0.15,0.75,0],[0,0.4,0],[0.7,0.65,0.9],[1,0.75,0.25],[0.75,0,0],[0,1,1],[0.55,0.55,0.55]])	
	
	if (type=='NO') | (type=='NOSE'):
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]));
		A_Cumulative=np.zeros(d['tv'].size)
		cnt=0
		for k in meta['LUT']['Derived']['ASET'].keys():
			ind=np.where(meta['LUT']['Raw']['ASET']['Name']==k)[0][0]
			cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
			plt.bar(d['tv'],d['A NOSE'][:,cnt]/1e3,0.8,bottom=A_Cumulative,facecolor=cl,label=k)
			A_Cumulative=A_Cumulative+d['A NOSE'][:,cnt]/1e3; 
			cnt=cnt+1
		#plt.plot(d['tv'],d['A'/1e3,'ks',ms=2.5,mec='k',mfc='w',mew=0.5)
		ax.set(xticks=np.arange(1950,2225+1,5),ylabel='Implementation level (Kha yr$^{-1}$)',
			   xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[d['tv'][iT][0]-0.75,d['tv'][iT][-1]+0+.75],ylim=[0,140]) #
		plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=1,fontsize=6)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TimeSeriesByASET_NOSE','png',900)
	elif (type=='Licensees'):
		d['A Licensees']=d['A All']-d['A NOSE']
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]));
		A_Cumulative=np.zeros(d['tv'].size)
		cnt=0
		for k in meta['LUT']['Derived']['ASET'].keys():
			ind=np.where(meta['LUT']['Raw']['ASET']['Name']==k)[0][0]
			cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
			plt.bar(d['tv'],d['A Licensees'][:,cnt]/1e3,0.8,bottom=A_Cumulative,facecolor=cl,label=k)
			A_Cumulative=A_Cumulative+d['A Licensees'][:,cnt]/1e3;
			cnt=cnt+1
		#plt.plot(d['tv'],d['A'/1e3,'ks',ms=2.5,mec='k',mfc='w',mew=0.5)
		ax.set(xticks=np.arange(1950,2225+1,10),ylabel='Implementation level (Kha yr$^{-1}$)',
			   xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[d['tv'][iT][0]-0.75,d['tv'][iT][-1]+0+.75],ylim=[0,240]) #
		plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=1,fontsize=6)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TimeSeriesByASET_Licensees','png',900)
	else:
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]));
		A_Cumulative=np.zeros(d['tv'].size)
		cnt=0
		for k in meta['LUT']['Derived']['ASET'].keys():
			ind=np.where(meta['LUT']['Raw']['ASET']['Name']==k)[0][0]
			cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
			plt.bar(d['tv'],d['A All'][:,cnt]/1e3,0.8,bottom=A_Cumulative,facecolor=cl,label=k)
			A_Cumulative=A_Cumulative+d['A All'][:,cnt]/1e3; 
			cnt=cnt+1
		#plt.plot(d['tv'],d['A'/1e3,'ks',ms=2.5,mec='k',mfc='w',mew=0.5)
		ax.set(xticks=np.arange(1950,2225+1,10),ylabel='Implementation level (Kha yr$^{-1}$)',
			   xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[d['tv'][iT][0]-0.75,d['tv'][iT][-1]+0+.75],ylim=[0,240]) #
		plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=1,fontsize=6)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TimeSeriesByASET_NO_And_Licencees','png',900)

	return

#%%
def Plot_AIL_NOSE_ByProjectType_FromModel(meta,pNam,mos,tv,iScn):
	nSS='All'
	nYS='All'
	iScn=1
	iT=np.where( (tv>=1960) & (tv<=2050) )[0]
	#cl=np.random.random((len(meta['LUT']['Derived']['ASET'].keys()),3))
	cl=np.array([[0.7,0.67,0.64],[0.9,0.87,0.84],[0.6,0.75,1],[0,0,0.5],[1,1,0],[0.6,1,0],[0.15,0.75,0],[0,0.4,0],[0.7,0.65,0.9],[1,0.75,0.25],[0.75,0,0],[0,1,1],[0.55,0.55,0.55]])
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7.5));
	A_Cumulative=np.zeros(tv.size)
	cnt=0
	for nPS in meta['LUT']['Derived']['ASET'].keys():
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
		if nPS=='Fill Planting':
			A=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Fill Planting']['Ensemble Mean'][:,0,iSS,iYS]/1000#*meta[pNam]['Project']['AEF']
		else:
			A=mos[pNam]['Scenarios'][iScn]['Sum']['Area_Planting']['Ensemble Mean'][:,iPS,iSS,iYS]/1000#*meta[pNam]['Project']['AEF']
		if nPS!='Harvest and Planting':
			plt.bar(tv,A,0.8,bottom=A_Cumulative,facecolor=cl[cnt,:],label=nPS)
			A_Cumulative=A_Cumulative+A; 
		cnt=cnt+1
	ax.set(xticks=np.arange(1950,2225+1,10),ylabel='Implementation level (Kha yr$^{-1}$)',
		   xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[tv[iT][0]-0.75,tv[iT][-1]+0+.75],ylim=[0,250]) #
	plt.legend(frameon=False,loc='upper right',facecolor=[1,1,1],labelspacing=0.25,ncol=1,fontsize=6)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_AIL_ByPT_FromModel','png',900)
	return

#%%
def Plot_AreaTreated_FrequencyByFSC(meta,t0,t1):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_SummaryByTimeAndFSC.pkl')
	iT=np.where( (d['Year']>=t0) & (d['Year']<=t1) )[0]

	ind=np.where(np.isin(d['FSC'],meta['Param']['Raw']['FSC']['NO List Name'])==True)[0]
	d['FSC']=d['FSC'][ind]
	d['Area']=d['Area'][:,ind]

	ind=np.where(np.sum(d['Area'],axis=0)>100)[0]
	d['FSC']=d['FSC'][ind]
	d['Area']=d['Area'][:,ind]

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
def Plot_AreaTreated_FrequencyByASET(meta,pNam,t0,t1):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	d=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ASET_SummaryByTime.pkl')

	cd=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	bin=np.arange(1,cd.size+1)
	iT=np.where( (d['tv']>=t0) & (d['tv']<=t1) )[0]
	A_Treat=d['A NOSE'][iT[0],:]
	ord=np.argsort(A_Treat);A_Treat=A_Treat[ord];cd=cd[ord] # Put in order

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,7))
	ax.bar(bin,A_Treat/1e3,0.8,facecolor=[0.25,0.5,1])
	for i in range(bin.size):
		txt0 = ('{:,}'.format(int(A_Treat[i])))
		txt=str(txt0) + '\n(' + str(np.round(A_Treat[i]/np.sum(A_Treat)*100,decimals=1)) + '%)'
		ax.text(bin[i],A_Treat[i]/1e3+1,txt,ha='center',fontsize=7)
	ax.set(ylabel='Area treated\n(hectares x 1000 per year)',xticks=bin,xticklabels=cd,ylim=[0,np.max(A_Treat/1e3)+4],xlim=[0.5,bin.size+0.5])
	plt.xticks(rotation=30,ha='right')
	plt.tight_layout()
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_FrequencyByASET_' + str(t0) + 'to' + str(t1),'png',900)

	return

#%%
# def NOSE_PlotEmissionsAnn(meta,mos,pNam,cNam):
# 	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
# 	psL=['Salvage and Planting','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting'] #'Harvest and Planting NSR Backlog'
# 	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,9))
# 	ax.plot(tv,np.zeros(tv.size),'k-',lw=3,color=[0.8,0.8,0.8])
# 	tm=np.mean(tv[1:].reshape(-1,10),axis=1)
# 	ysum=np.zeros(tm.size)
# 	for nPS in psL: #meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:
# 		nSS='All';nYS='All'
# 		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
# 		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
# 		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
# 		y=((mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded']['Ensemble Mean'][:,iPS,iSS,iYS])/2)*meta[pNam]['Project']['AEF']/1e6
# 		y=np.mean(y[1:].reshape(-1,10),axis=1)
# 		ysum=ysum+y
# 		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
# 		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]	
# 		iT2=np.where( (tm>=1980) & (tm<=2110) )[0]
# 		plt.plot(tm[iT2],y[iT2],'-o',ms=4,color=cl,mec=cl,mfc=cl,lw=1.25,mew=1.25,label=nPS)
# 	plt.plot(tm[iT2],ysum[iT2],'-o',ms=3,color=[0.65,0.65,0.65],mec=[0.65,0.65,0.65],mfc='w',lw=1.25,mew=1.25,label='All Non-obligation Stand Establishment')
# 	flg=0
# 	if flg==1:
# 		yNM=mos[pNamNMC]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded']['Ensemble Mean'][:,0,0,0]*meta[pNamNMC]['Project']['AEF']/1e6
# 		yNM=np.mean(yNM[1:].reshape(-1,10),axis=1)
# 		plt.plot(tm[iT2],yNM[iT2],'-^',ms=3,color=[0.5,0,1],mec=[0.5,0,1],mfc='w',lw=0.75,mew=0.75,label='Nutrient Management')
# 		plt.plot(tm[iT2],yNM[iT2]+ysum[iT2],'-s',ms=3,color='k',mec='k',mfc='k',lw=0.75,mew=0.75,label='Total')
# 	ax.set(xticks=np.arange(1800,2200,10),ylabel='Annual $\Delta$ emissions (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[-1.8,0.6],xlim=[1979.5,2110.5])
# 	leg=ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
# 	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
# 	plt.tight_layout()
# 	if meta['Graphics']['Print Figures']=='On':
# 		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summary_DeltaE_Ann','png',	900)
# 	return

#%%
def NOSE_PlotEmissionsCumu(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	pSL=['Salvage and Planting','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting'] #'Harvest and Planting NSR Backlog'
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=3,color=[0.8,0.8,0.8])
	y=np.zeros((tv.size,len(pSL)))
	for i in range(len(pSL)): #meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:
		nPS=pSL[i]
		nSS='All';nYS='All'
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y[:,i]=((mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded_Cumulative_from_tref']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded_Cumulative_from_tref']['Ensemble Mean'][:,iPS,iSS,iYS])/2)*meta[pNam]['Project']['AEF']/1e6
		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
		plt.plot(tv,y[:,i],'-',ms=3,color=cl,mec=cl,mfc='w',lw=1.5,mew=0.75,label=nPS)
	plt.plot(tv,np.sum(y,axis=1),'-',color=[0.65,0.65,0.65],mec=[0.65,0.65,0.65],mfc='w',lw=2,mew=0.75,label='All Non-obligation Stand Establishment')
	ax.set(xticks=np.arange(1800,2200,5),ylabel='Cumulative $\Delta$ emissions (MtCO$_2$e)',xlabel='Time, years',yticks=np.arange(-200,200,10),xlim=[1979.5,2110.5],ylim=[-80,5])
	leg=ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summary_DeltaE_Cumulative','png',900)
	return

#%%
def Plot_EmissionsCumu_TimeSeries_FromCurrentYear(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	vStat='Ensemble Mean'
	oper='Sum'
	iT2=np.where( (tv>=2015) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=1,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

		# Half substitutions
		#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,iPS,0,iYS])/2)*meta[pNam]['Project']['AEF']/1e6

		# Without substitutions
		y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]/1e6

		ysum=ysum+np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)

	plt.plot(tv[iT2],np.cumsum(ysum[iT2]),'--',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')


	ax.legend(loc='lower left',frameon=True,facecolor='w',edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-1,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,1,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2150,5),yticks=np.arange(-20,20,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlabel='Time, years',ylim=[-4,1.5],xlim=[meta[pNam]['YearCurrent']-7,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumu_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return

#%%
def Plot_EmissionsCumuPerHectare_TimeSeries_FromCurrentYear(meta,mos,pNam,cNam):
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
	oper='Mean'
	iT2=np.where( (tv>=2010) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_mu=np.zeros((tv.size,len(psL)))
	cnt=0
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

		# Half substitutions
		#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,iPS,0,iYS])/2)*meta[pNam]['Project']['AEF']

		# Without substitutions
		y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]

		y_mu[:,cnt]=np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
		cnt=cnt+1

	plt.plot(tv[iT2],np.cumsum(np.mean(y_mu[iT2],axis=1)),'--',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	flg=0
	if flg==1:
		# Add Nutrient management
		#pNam=pNamNMC
		#iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='2022')[0][0]
		#y=((mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,0,0,iYS]+mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionIncluded'][vStat][:,0,0,iYS])/2)*meta[pNam]['Project']['AEF']/1e6
		mu=np.nanmean(yNM,axis=1)
		iT3=np.where(tv<2022)[0]
		mu[iT3]=0
		A=33500
		yNM1=A*np.cumsum(mu)/1e6
		plt.plot(tv[iT2],yNM1[iT2],'-',ms=3,color=[0.5,0,1],mec=[0.5,0,1],mfc='w',lw=1,mew=0.75,label='Nutrient Management')
		plt.plot(tv[iT2],yNM1[iT2]+np.cumsum(y_mu[iT2]),'--',ms=3,color='k',mec='k',mfc='k',lw=1.5,mew=0.75,label='Total')
	
	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-300,'Net removals',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,300,'Net emissions',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2010,2200,10),yticks=np.arange(-2000,2000,50),ylabel='Cumulative GHG impact (tCO$_2$e ha$^{-1}$)',xlabel='Time, years',ylim=[-350,350],xlim=[meta[pNam]['YearCurrent']-7,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumuPerHectare_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return

#%%
def Plot_EmissionsAnnual_TimeSeries_FromCurrentYear(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	vStat='Ensemble Mean'
	operSpace='Sum'
	iT2=np.where( (tv>=2015) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=1,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata'][operSpace]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
		plt.plot(tv[iT2],y[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
	plt.plot(tv[iT2],ysum[iT2],'-',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')
	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-0.25,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,0.25,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2150,5),yticks=np.arange(-20,20,0.1),ylabel='Annual GHG impact (MtCO$_2$e yr$^{-1}$)',xlabel='Time, years',ylim=[-0.3,0.3],xlim=[meta[pNam]['YearCurrent']-7,2075])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnn_TimeSeries_FromYear' + str(meta[pNam]['YearCurrent']),'png',900)#%%
	return

#%%
def Plot_EmissionsAnn_TimeSeries_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	vStat='Ensemble Mean'
	oper='Sum'
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
		plt.plot(tv[iT2],y[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
	plt.plot(tv[iT2],ysum[iT2],'-',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	ax.legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1968,-2.15,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1968,1.55,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.5),ylabel='Annual GHG impact (MtCO$_2$e yr$^{-1}$)',xlabel='Time, years',ylim=[-2.5,2],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnn_TimeSeries_Completed','png',900)
	return

#%%
def Plot_EmissionsCumu_TimeSeries_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	vStat='Ensemble Mean'
	oper='Sum'
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_Atmosphere_SubstitutionExcluded'][vStat][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
	plt.plot(tv[iT2],np.cumsum(ysum[iT2]),'-',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	ax.legend(loc='center left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1968,-60,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1968,20,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,10),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlabel='Time, years',ylim=[-70,30],xlim=[1965,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsCumu_TimeSeries_Completed','png',900)
	return

#%%
def QA_CalcBenefitForEachStand(meta,pNam,mos):
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
def TabulateSpeciesComposition(meta,pNam,t0,t1):
	d=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'].copy()
	for iP in range(6):
		print(iP)
		yr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
		indT=np.where( (yr>=t0) & (yr<=t1) )
		yr=yr[indT]
		fsc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data'][indT]
		sph=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SPH_Planted.tif')['Data'][indT]
		for iS in range(1,7):
			spc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_PL_SPECIES_CD' + (str(iS)) + '.tif')['Data'][indT]
			frac=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_PL_SPECIES_PCT' + (str(iS)) + '.tif')['Data'][indT].astype('float')/100
			uS=np.unique(spc[spc>0])
			for iU in range(uS.size):
				ind=np.where( (np.isin(fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) & (spc==uS[iU]) )
				cd=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],uS[iU])[0]
				d[cd]=d[cd]+np.sum(sph[ind]*frac[ind])
	# Save
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\SummarySpeciesComposition_' + str(t0) + 'to' + str(t1) + '.pkl',d)
	return

#%% Identify ecosystem restoration projects
def IdentifyEcosystemRestoration(meta):

	cs=np.array([])
	for iP in range(6):
		print(iP)
		yr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_Year.tif')['Data']
		fsc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
		ob1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_SILV_OBJECTIVE_CODE_1.tif')['Data']
		spc1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iP+1) + '_PL_SPECIES_CD1.tif')['Data']

		ind=np.where( (np.isin(fsc,meta['Param']['Raw']['FSC']['NO List ID'])==True) & (yr>2000) )
		cs0=np.column_stack(([fsc[ind],ob1[ind],spc1[ind]]))
		if cs.size==0:
			cs=cs0
		else:
			cs=np.vstack((cs,cs0))

	u=np.unique(cs,axis=0)
	n=u.shape[0]
	d={}
	d['FSC']=np.array(['' for _ in range(n)],dtype=object)
	d['Objective']=np.array(['' for _ in range(n)],dtype=object)
	d['Spc1']=np.array(['' for _ in range(n)],dtype=object)
	for i in range(u.shape[0]):
		d['FSC'][i]=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'],u[i,0])[0]
		try:
			d['Objective'][i]=cbu.lut_n2s(meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_OBJECTIVE_CODE_1'],u[i,1])[0]
		except:
			pass
		try:
			d['Spc1'][i]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],u[i,2])[0]
		except:
			pass
	df=pd.DataFrame.from_dict(d)
	df.to_excel(r'C:\Data\BCFCS\BCFCS_NOSEC\Inputs\ScanForRestoration.xlsx')

	return

#%%
def TabulateUnderplantingByBurnSeverityClass(meta):
	zBS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\BurnSevComp_SevClassLast.tif')['Data']
	zU=np.zeros(zBS.shape,dtype='int8')
	for iY in range(6):
		z0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(iY+1) + '_ASET.tif')['Data']
		ind=np.where(z0==meta['LUT']['Derived']['ASET']['Underplanting'])
		zU[ind]=1
	
	ind=np.where( (zBS>0) & (zBS<5) & (zU==1) )
	N_tot=ind[0].size
	
	d={}
	for k in meta['LUT']['Derived']['burnsev_comp1'].keys():
		ind=np.where( (zBS==meta['LUT']['Derived']['burnsev_comp1'][k]) & (zU==1) )
		d[k]=np.round(ind[0].size/N_tot*100,decimals=1)
	return d

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
def UpdateDatabase(meta,pNam,mos,mosWF,metaNM,pNamNM,mosNM,mosNM_WF):

	ac1='Silviculture investments'
	vStat='Ensemble Mean'
	operSpace='Sum'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Variable list
	vL=['E_Atmosphere_SubstitutionExcluded','E_Atmosphere_SubstitutionIncludedHalf','E_LULUCF','V_ToMill_MerchTotal','Cost Total','Revenue Gross','Revenue Net']

	# Initialize
	n=25000
	d={}
	d['Action Category 1']=np.array(['' for _ in range(n)],dtype=object)
	d['Action Category 2']=np.array(['' for _ in range(n)],dtype=object)
	d['Action Category 3']=np.array(['' for _ in range(n)],dtype=object)
	d['Funding Source Code']=np.array(['' for _ in range(n)],dtype=object)
	d['Comparison Name']=np.array(['' for _ in range(n)],dtype=object)
	d['Operation Spatial']=np.array(['' for _ in range(n)],dtype=object)
	d['Operation Temporal']=np.array(['' for _ in range(n)],dtype=object)
	d['Ensemble Statistic']=np.array(['' for _ in range(n)],dtype=object)
	d['Scope Investments']=np.array(['' for _ in range(n)],dtype=object)
	#d['Multiplier']=meta[pNam]['Project']['Multi']*np.ones(n,dtype='int32')
	d['Year Implemented']=np.array(['' for _ in range(n)],dtype=object)
	d['Time Start']=np.zeros(n,dtype='int32')
	d['Time End']=np.zeros(n,dtype='int32')
	d['Area Treated (ha/yr)']=np.zeros(n,dtype='int32')
	for v in vL:
		d[v]=np.zeros(n,dtype='float')

	# Non obligation stand establishment

	cNam='NOSE1'
	ac2='Non-obligation stand establishment'
	#psL=['All','Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	#  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	psL=meta[pNam]['Project']['Strata']['Project Type']['Unique CD']
	ysL=meta[pNam]['Project']['Strata']['Year']['Unique CD']
	ssL=meta[pNam]['Project']['Strata']['Spatial']['Unique CD']
	cnt=0
	for ps in psL:
		for ss in ssL:
			for ys in ysL:

				iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==ps)[0]
				iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==ys)[0]
				iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==ss)[0]

				# Completed

				# Add annual fluxes in specific reporting years
				yrL=np.array([2020,2030,2050,2100])
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=yr
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(mos[pNam]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']),decimals=1)
					cnt=cnt+1

				# Add cumulative fluxes between 1960 and future reporting years
				yrL=np.arange(2030,2140,10)
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=1960
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(np.sum(mos[pNam]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])),decimals=1)
					cnt=cnt+1

				# Completed + CAP

				# Add annual fluxes in specific reporting years
				yrL=np.array([2020,2030,2050,2100])
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed plus CAP'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=yr
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(mosWF[pNam]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']),decimals=1)
					cnt=cnt+1

				# Add cumulative fluxes between year of implementation (minus 5 to cover open burning) and future reporting years
				yrL=np.arange(2030,2140,10)
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed plus CAP'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=1960
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(np.sum(mosWF[pNam]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])),decimals=1)
					cnt=cnt+1

	# Add nutrient management
	ac2='Aerial Nutrient Application'
	cNam='Low Harvest'
	psL=metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']
	ssL=metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']
	ysL=metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']
	for ps in psL:
		for ss in ssL:
			for ys in ysL:
				iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']==ps)[0]
				iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']==ys)[0]
				iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']==ss)[0]

				# Completed

				# Add annual fluxes in specific reporting years
				yrL=np.array([2020,2030,2050,2100])
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=yr
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(mosNM[pNamNM]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']),decimals=1)
					cnt=cnt+1

				# Add cumulative fluxes between year of implementation (minus 5 to cover open burning) and future reporting years
				yrL=np.arange(2030,2140,10)
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=1960
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(np.sum(mosNM[pNamNM]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])),decimals=1)
					cnt=cnt+1

				# Completed + climate action planning

				# Add annual fluxes in specific reporting years
				yrL=np.array([2020,2030,2050,2100])
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed plus CAP'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=yr
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(mosNM_WF[pNamNM]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']),decimals=1)
					cnt=cnt+1

				# Add cumulative fluxes between year of implementation (minus 5 to cover open burning) and future reporting years
				yrL=np.arange(2030,2140,10)
				for yr in yrL:
					d['Action Category 1'][cnt]=ac1
					d['Action Category 2'][cnt]=ac2
					d['Action Category 3'][cnt]=ps
					d['Funding Source Code'][cnt]=ss
					d['Year Implemented'][cnt]=ys
					d['Comparison Name'][cnt]=cNam
					d['Scope Investments'][cnt]='Completed plus CAP'
					d['Area Treated (ha/yr)'][cnt]=0
					d['Operation Spatial'][cnt]=operSpace
					d['Operation Temporal'][cnt]='Sum'
					d['Ensemble Statistic'][cnt]='Ensemble Mean'
					d['Time Start'][cnt]=1960
					d['Time End'][cnt]=yr
					iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
					for v in vL:
						d[v][cnt]=np.round(np.nan_to_num(np.sum(mosNM_WF[pNamNM]['Delta'][cNam]['ByStrata'][operSpace][v][vStat][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])),decimals=1)
					cnt=cnt+1

	# Truncate
	for k in d.keys():
		d[k]=d[k][0:cnt]

	# Save to file
	#df=pd.DataFrame.from_dict(d)
	#df.to_excel()
	pthout=r'C:\Users\rhember\OneDrive - Government of BC\Projects\Forest Carbon Summary\Data' + '\\R' + str(meta[pNam]['YearCurrent']+1) + '\\BC-FOR Climate Mitigation Database Summary.xlsx'
	gu.PrintDict(d,pthout,SheetName='Summary')

	return

#%%
def DefineStrata(meta,pNam,dmec,lsat):

	# By project type
	cd=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	id=np.array(list(meta['LUT']['Derived']['ASET'].values()))
	meta[pNam]['Project']['Strata']['Project Type']['Unique ID']=np.append(0,id)
	meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=np.append('All',cd)
	#meta[pNam]['Project']['Strata']['Project Type']['ID']=np.zeros(meta[pNam]['Project']['N Stand'],dtype='int8') # This is populated by inventory processing

	# By time (last year of implementation)
	t=np.arange(meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']+1,1)
	meta[pNam]['Project']['Strata']['Year']['Unique CD']=np.append('All',t.astype(str))
	meta[pNam]['Project']['Strata']['Year']['Unique ID']=np.append(0,t)
	iScn=1
	for iStand in range(meta[pNam]['Project']['N Stand']):
		ind=np.where( (dmec[iScn][iStand]['Index to Event Inciting NOSE']!=9999) )[0]
		if ind.size>0:
			Year0=np.floor(dmec[iScn][iStand]['Year'][ind[-1]])
			iT=np.where(t==Year0)[0]
			if iT.size==0:
				continue
			meta[pNam]['Project']['Strata']['Year']['ID'][iStand]=t[iT]

	# By space
	flgSpace='ByFSC'
	if flgSpace=='ByFSC':
		# By Funding Source Code
		List=['FTM','FRP','FCE','FIP','FCM','LFP','FTL','FES','GA','FED']
		ind=np.where(np.isin(meta['Param']['Raw']['FSC']['NO List Name'],List)==True)[0]
		cd=meta['Param']['Raw']['FSC']['NO List Name'][ind]
		id=meta['Param']['Raw']['FSC']['NO List ID'][ind]
		meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=np.append('All',cd)
		meta[pNam]['Project']['Strata']['Spatial']['Unique ID']=np.append(0,id)
		iScn=1
		for iStand in range(meta[pNam]['Project']['N Stand']):
			ind=np.where( (dmec[iScn][iStand]['Index to Event Inciting NOSE']!=9999) )[0]
			if ind.size==0:
				continue
			if ind.size>1:
				ind=ind[-1]
			meta[pNam]['Project']['Strata']['Spatial']['ID'][iStand]=dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][ind]
	elif flgSpace=='ByBGC':
		# By space (BGC zone)
		u=np.unique(lsat['ID_BGCZ'])
		cd=np.array([])
		for iU in range(u.size):
			 cd=np.append(cd,u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU]))
		meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=np.append('All',cd)
		meta[pNam]['Project']['Strata']['Spatial']['Unique ID']=np.append(0,u)
		meta[pNam]['Project']['Strata']['Spatial']['ID']=lsat['ID_BGCZ']
	else:
		pass
	return meta

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_CurrentYear(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.text(1970,-0.25,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1970,0.25,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1965,2100],ylim=[-0.3,0.3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y2=np.cumsum(ysum)
	ax2.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-4,4])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_CurrentYear','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_ServicePlanProjectionYear(meta,mos,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent']+1,meta[pNam]['YearCurrent']+1],[-20,20],'k--',lw=2,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	cNam='NOSE1'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,0]/1e6
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.text(1970,-0.25,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1970,0.25,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1965,2100],ylim=[-0.3,0.3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	y2=np.cumsum(ysum)
	ax2.plot(tv[iT2],y2[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-4,4])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_ServicePlanProjectionYear','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_ServicePlanProjectionYear_WithNM(meta,mos,mosNM,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent']+1,meta[pNam]['YearCurrent']+1],[-20,20],'k--',lw=2,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	cNam='NOSE1'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,0]/1e6
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0

	# Add Nutrient management
	yNM=np.nan_to_num(mosNM['BCFCS_NMC']['Delta']['Low Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,0]/1e6)
	iTnm=np.where( (tv<meta[pNam]['YearCurrent']) )[0]
	yNM[iTnm]=0
	ysum=ysum+yNM

	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	#ax.plot(tv[iT2],yNM[iT2],'--',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.text(1970,-0.25,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1970,0.25,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1965,2100],ylim=[-0.4,0.4])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	#ax.plot(tv[iT2],np.cumsum(yNM)[iT2],'--',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75)
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-6,6])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	# Print number for Service Plan Projection
	iT=np.where( (tv==2050) )[0]
	print(ysumc[iT])

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_ServicePlanProjectionYear_WithNM','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-1.6,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.6,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-2.2,2.2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,10),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-70,70])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_Completed','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-3.2,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.7,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-3.5,2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,20),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-220,2/3.5*220])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-1.6,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.6,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-2.2,2.2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,10),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-70,70])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_Completed','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_Completed_WithNM(meta,mos,mosNM,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	cNam='NOSE1'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	# Add nutrient management
	ysum=ysum+mosNM['BCFCS_NMC']['Delta']['Low Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,0]/1e6
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-4.6,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.6,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-5,2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,50),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-320,2/5*320])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_Completed_WithNM','png',900)
	return

#%%
def Plot_EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP_WithNM(meta,mos,mosNM,pNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	cNam='NOSE1'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6

	# Add nutrient management
	y=mosNM['BCFCS_NMC']['Delta']['Low Harvest']['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6
	ysum=ysum+y

	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-4.6,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.6,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-5,2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,50),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-320,2/5*320])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsAnnAndCumu_TimeSeries_CompletedAndCAP_WithNM','png',900)
	return

#%%
def Plot_YieldAnnAndCumu_TimeSeries_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	v='V_ToMill_MerchTotal'

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-0.8,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,2.7,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-1,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,10),ylabel='',xlim=[1965,2100],ylim=[-50,150])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\YieldAnnAndCumu_TimeSeries_Completed','png',900)
	return

#%%
def Plot_YieldAnnAndCumu_TimeSeries_CompletedAndCAP(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]

	v='V_ToMill_MerchTotal'

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	ax.text(1967,-0.85,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,2.85,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-1,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionExcluded']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,10),ylabel='',xlim=[1965,2100],ylim=[-50,150])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\YieldAnnAndCumu_TimeSeries_CompletedAndCAP','png',900)
	return

#%%
def ProjectFuture(meta,pNam,mos):
	# This function takes the current year per-hectare time series and applies it into the future at a projected AIL
	# Notes: Cumulative values are incorrect - do not use!

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
	#iPS=np.where( meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=='All' )[0][0]
	#iSS=np.where( meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All' )[0][0]
	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]

	dCP=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_CostsAndPrices.xlsx','Data')

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	fracPT={'Salvage and Planting Post Beetle':0.1,
		   'Knockdown and Planting':0.05,
		   'Road Rehabilitation':0.0,
		   'Underplanting':0.7,
		   'Replanting':0.0,
		   'Fill Planting':0.15,
		   'Ecosystem Restoration':0}

	mosWF=copy.deepcopy(mos)
	for cNam in mos[pNam]['Delta'].keys():
		for v in mos[pNam]['Delta'][cNam]['ByStrata']['Mean'].keys():
			for nPS in psL:
				iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
				ail_pt=meta[pNam]['AIL CAP']*fracPT[nPS]
				y0=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,:,iYS]
				tv0=tv[iT]
				for yr in range(meta[pNam]['YearCurrent']+1,tv[-1],1):
					iT1=np.where( (tv>=yr) )[0]
					iT2=np.where( (tv0>=meta[pNam]['YearCurrent']) )[0]
					iT2=iT2[0:iT1.size]
					mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,iPS,:,0]=mosWF[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,iPS,:,0]+ail_pt*y0[iT2,:]

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
	# This function takes the current year per-hectare time series and applies it into the future at a projected AIL
	# Notes: Cumulative values are incorrect - do not use!

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	YearProjection=meta[pNam]['YearCurrent']+1
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	fracPT={'Salvage and Planting Post Beetle':0.1,
		   'Knockdown and Planting':0.05,
		   'Road Rehabilitation':0.0,
		   'Underplanting':0.7,
		   'Replanting':0.0,
		   'Fill Planting':0.15,
		   'Ecosystem Restoration':0}

	mosSPP=copy.deepcopy(mos)
	for cNam in mos[pNam]['Delta'].keys():
		for v in mos[pNam]['Delta'][cNam]['ByStrata']['Mean'].keys():
			mosSPP[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,:,:,:]=0
			for nPS in psL:
				iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
				ail_pt=meta[pNam]['AIL CAP']*fracPT[nPS]
				y0=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,:,iYS]
				tv0=tv[iT]
				iT1=np.where( (tv>=YearProjection) )[0]
				iT2=np.where( (tv0>=meta[pNam]['YearCurrent']) )[0]
				iT2=iT2[0:iT1.size]
				mosSPP[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT1,iPS,:,0]=ail_pt*y0[iT2,:]
	return mosSPP

#%%