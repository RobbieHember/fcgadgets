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
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun as cbr
import fcgadgets.macgyver.util_fcs_graphs as ufcs
import fcgadgets.macgyver.util_fcs_qa as uqa
warnings.filterwarnings("ignore")
gp=gu.SetGraphics('Manuscript')

#%%
def ImportModelResults():
	pNam='BCFCS_NOSE'
	meta=gu.ipickle(r'C:\Data\BCFCS\BCFCS_NOSE\Inputs\Metadata.pkl')
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
	pNam='BCFCS_NOSE'
	meta=gu.ipickle(r'D:\Modelling Projects\BCFCS_NOSE\Inputs\Metadata.pkl')#meta=gu.ipickle(r'C:\Data\BCFCS\BCFCS_NOSE\Inputs\Metadata.pkl')
	meta[pNam]['Project']['Multi']=1e6
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	print('Importing GHG emissions and economics by scenario and stratum')
	mos=cbu.Import_MOS_ByScnAndStrata_GHGEcon(meta,pNam)
	cNam='NOSE'
	mos[pNam]['Delta']={}
	mos[pNam]['Delta'][cNam]={'iB':0,'iP':1}
	print('Importing GHG emissions and economics for scenario comparisons')
	mos=cbu.Import_MOS_ByScnComparisonAndStrata(meta,pNam,mos)
	print('Importing DMEC')
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
def QA_Print_SummarySparseToSpreadsheet(meta,pNam,iScn,dmec,mos_ByStand,cNam):
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
	d['E_NSB Delta']=np.zeros(n)
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
			d['E_NSB Delta'][c]=int(mos_ByStand['Delta'][cNam]['1971-2050']['Sum']['E_NSB'][i])
			if (dmec[iScn][i]['ASET'][j]!=9999) & (dmec[iScn][i]['ASET'][j]!=0):
				d['Pl Spc1 CD'][c]=cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],dmec[iScn][i]['PL_SPECIES_CD1'][j])[0]
				d['Pl Spc1 %'][c]=dmec[iScn][i]['PL_SPECIES_PCT1'][j]
				d['Pl SPH'][c]=dmec[iScn][i]['Planted SPH'][j]
			c=c+1

	#df=pd.DataFrame.from_dict(d)
	pthout=meta['Paths'][pNam]['Data'] + '\\Inputs\\SummarySparseSample.xlsx'
	gu.PrintDict(d,pthout,SheetName='Sheet1')

	return

#%%
def QA_Plot_Carbon_TS_ByStand(meta,pNam,dmec,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	xlim=[1950,2050]
	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	if 'pt' in kwargs.keys():
		ptL=kwargs['pt']
	else:
		ptL=meta['LUT']['Derived']['ASET'].keys()

	for iBat in range(meta[pNam]['Project']['N Batch']):
		indBat=cbu.IndexToBatch(meta[pNam],iBat)
		d0=cbu.LoadSingleOutputFile(meta,pNam,0,0,iBat)
		d1=cbu.LoadSingleOutputFile(meta,pNam,1,0,iBat)
		for k in ptL:
			ind=np.where( (meta[pNam]['Project']['Strata']['Project Type']['ID'][indBat]==meta['LUT']['Derived']['ASET'][k]) )[0]

			if ind.size==0:
				continue
			#if k!='Knockdown and Planting':
			#	continue

			cnt=0
			max_per_pt=30
			ivl=1
			for i in range(0,ind.size,1):

# 				if cnt>max_per_pt:
# 					# Limit the sample to a total of x per project type
# 					continue

				plt.close('all'); fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(18,10));
				ax[0,0].plot(tv,d0['C_Biomass'][:,ind[i]],'b-',color=cl1,label='Baseline')
				ax[0,0].plot(tv,d1['C_Biomass'][:,ind[i]],'g--',color=cl2,label='Actual')
				ax[0,0].set(ylabel='Biomass',xlim=xlim)
				ax[0,0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
				ax[0,1].plot(tv,d0['C_DeadWood'][:,ind[i]],'b-',color=cl1)
				ax[0,1].plot(tv,d1['C_DeadWood'][:,ind[i]],'g--',color=cl2)
				ax[0,1].set(ylabel='Dead wood',xlim=xlim)
				ax[1,0].plot(tv,d0['E_Wildfire_ForestSector_Total'][:,ind[i]],'ob-',color=cl1,ms=3)
				ax[1,0].plot(tv,d1['E_Wildfire_ForestSector_Total'][:,ind[i]],'sg--',color=cl2,ms=3)
				ax[1,0].set(ylabel='Wildfire emissions',xlim=xlim)
				ax[1,1].plot(tv,d0['V_ToMill_MerchTotal'][:,ind[i]],'ob-',color=cl1)
				ax[1,1].plot(tv,d1['V_ToMill_MerchTotal'][:,ind[i]],'sg--',color=cl2)
				ax[1,1].set(ylabel='Volume to mill, merch live',xlim=xlim)
				ax[2,0].plot(tv,d0['E_OpenBurning_ForestSector_Total'][:,ind[i]],'b-',color=cl1)
				ax[2,0].plot(tv,d1['E_OpenBurning_ForestSector_Total'][:,ind[i]],'g--',color=cl2)
				ax[2,0].set(ylabel='Open burning\nemissions',xlim=xlim)
				#ax[2,1].plot(tv,d0['E_EnergySC_Bioenergy'][:,ind[i]],'b-')
				#ax[2,1].plot(tv,d1['E_EnergySC_Bioenergy'][:,ind[i]],'g--')
				#ax[2,1].set(ylabel='Bioenergy emissions')
				ax[2,1].plot(tv,d0['C_Litter'][:,ind[i]]+d0['C_Soil'][:,ind[i]],'b-',color=cl1)
				ax[2,1].plot(tv,d1['C_Litter'][:,ind[i]]+d1['C_Soil'][:,ind[i]],'g--',color=cl2)
				ax[2,1].set(ylabel='Litter+Soil',xlim=xlim)
				ax[3,0].plot(tv,d0['E_NSB'][:,ind[i]],'b-',color=cl1)
				ax[3,0].plot(tv,d1['E_NSB'][:,ind[i]],'g--',color=cl2)
				ax[3,0].set(ylabel='Annual GHG balance',xlim=xlim)
				ax[3,1].plot(tv,np.cumsum(d0['E_NSB'][:,ind[i]]),'b-',color=cl1)
				ax[3,1].plot(tv,np.cumsum(d1['E_NSB'][:,ind[i]]),'g--',color=cl2)
				ax[3,1].set(ylabel='Cumulative GHG balance',xlim=xlim)

				ind1=np.where(dmec[1][indBat[ind[i]]]['ASET']<9999)[0]
				yr1=dmec[1][indBat[ind[i]]]['Year'][ind1[-1]]
				plt.tight_layout()
				gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA\\' + k + '_' + str(indBat[ind[i]]) + '_' + str(int(yr1)),'png',200)
				cnt=cnt+1
	return

#%%
def QA_Plot_Carbon_TS_ByStandSelect(meta,pNam,dmec,**kwargs):
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

			ax[1,0].plot(tv,d0['E_Wildfire_ForestSector_Total'][:,ind[i]],'ob-',ms=3)
			ax[1,0].plot(tv,d1['E_Wildfire_ForestSector_Total'][:,ind[i]],'sg--',ms=3)
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

			ax[2,0].plot(tv,d0['E_OpenBurning_ForestSector_Total'][:,ind[i]],'b-')
			ax[2,0].plot(tv,d1['E_OpenBurning_ForestSector_Total'][:,ind[i]],'g--')
			ax[2,0].set(ylabel='Open burning emissions')

			#ax[2,1].plot(tv,d0['C_Litter'][:,ind[i]]+d0['C_Soil'][:,ind[i]]+d0['C_DeadWood'][:,ind[i]],'b-')
			#ax[2,1].plot(tv,d1['C_Litter'][:,ind[i]]+d1['C_Soil'][:,ind[i]]+d1['C_DeadWood'][:,ind[i]],'g--')
			ax[2,1].plot(tv,d0['C_DeadWood'][:,ind[i]],'b-')
			ax[2,1].plot(tv,d1['C_DeadWood'][:,ind[i]],'g--')
			#ax[2,1].plot(tv,d0['C_DeadStemMerch'][:,ind[i]],'b-')
			#ax[2,1].plot(tv,d1['C_DeadStemMerch'][:,ind[i]],'g--')
			ax[2,1].set(ylabel='Litter+Soil+Dead Wood')

			ax[3,0].plot(tv,d0['E_NSB'][:,ind[i]],'b-')
			ax[3,0].plot(tv,d1['E_NSB'][:,ind[i]],'g--')
			ax[3,0].set(ylabel='Annual GHG balance')
			ax[3,1].plot(tv,np.cumsum(d0['E_NSB'][:,ind[i]]),'b-')
			ax[3,1].plot(tv,np.cumsum(d1['E_NSB'][:,ind[i]]),'g--')
			ax[3,1].set(ylabel='Cumulative GHG balance')

			ind1=np.where(dmec[1][indBat[ind[i]]]['ASET']<9999)[0]
			yr1=dmec[1][indBat[ind[i]]]['Year'][ind1[-1]]
			plt.tight_layout()
			type=cbu.lut_n2s(meta['LUT']['Derived']['ASET'],meta[pNam]['Project']['ASET'][indBat[ind[i]]])[0]
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA\\' + type + '_' + str(indBat[ind[i]]) + '_' + str(int(yr1)),'png',200)

	return

#%%
def Plot_AreaTreated_TS_ByASET(meta,pNam,type,**kwargs):
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
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TS_ByASET_NOSE','png',900)
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
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TS_ByASET_Licensees','png',900)
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
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TS_ByASET_NO_And_Licencees','png',900)

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
def Plot_AreaTreated_Frequency_ByFSC(meta,t0,t1):
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
	plt.bar(bin,A,0.85,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	for i in range(bin.size):
		y=A[i]
		ax.text(bin[i],y+yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Implementation level (Kha)',xticklabels=cd,xlim=[0.25,bin.size+0.75],ylim=yl) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_Frequency_ByFSC_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def Plot_AreaTreated_Frequency_ByASET(meta,pNam,t0,t1):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	d=gu.ipickle(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\ASET_SummaryByTime.pkl')

	cd=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	bin=np.arange(1,cd.size+1)
	iT=np.where( (d['tv']>=t0) & (d['tv']<=t1) )[0]
	A_Treat=d['A NOSE'][iT[0],:]
	ord=np.argsort(A_Treat);A_Treat=A_Treat[ord];cd=cd[ord] # Put in order

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,7))
	ax.bar(bin,A_Treat/1e3,0.8,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	for i in range(bin.size):
		txt0 = ('{:,}'.format(int(A_Treat[i])))
		txt=str(txt0) + '\n(' + str(np.round(A_Treat[i]/np.sum(A_Treat)*100,decimals=1)) + '%)'
		ax.text(bin[i],A_Treat[i]/1e3+1,txt,ha='center',fontsize=7)
	ax.set(ylabel='Area treated\n(hectares x 1000 per year)',xticks=bin,xticklabels=cd,ylim=[0,np.max(A_Treat/1e3)+7],xlim=[0.5,bin.size+0.5])
	plt.xticks(rotation=30,ha='right')
	plt.tight_layout()
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_Frequency_ByASET_' + str(t0) + 'to' + str(t1),'png',900)

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
# 		y=((mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded']['Ensemble Mean'][:,iPS,iSS,iYS])/2)*meta[pNam]['Project']['AEF']/1e6
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
	pSL=['Salvage and Planting Post Beatle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting'] #'Harvest and Planting NSR Backlog'
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=3,color=[0.8,0.8,0.8])
	y=np.zeros((tv.size,len(pSL)))
	for i in range(len(pSL)): #meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:
		nPS=pSL[i]
		nSS='All';nYS='All'
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==nYS)[0][0]
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		y[:,i]=((mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB_Cumulative_from_tref']['Ensemble Mean'][:,iPS,iSS,iYS]+mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_Atmosphere_SubstitutionIncluded_Cumulative_from_tref']['Ensemble Mean'][:,iPS,iSS,iYS])/2)*meta[pNam]['Project']['AEF']/1e6
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
def Plot_Emissions_TS_Cumu_FromCurrentYear_ByASET(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=2015) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=1,color=[0,0,0])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting',
	  'Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='FIP')[0][0]

		#y=mos[pNam]['Delta'][cNam]['ByStrata'][oper]['E_NSB'][vStat][:,iPS,0,iYS]/1e6
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)

	plt.plot(tv[iT2],np.cumsum(ysum[iT2]),'--',ms=3,color=[0.25,0.25,0.25],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	ax.legend(loc='lower left',frameon=True,facecolor='w',edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-1.5,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,1.5,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2150,5),yticks=np.arange(-20,20,1),ylabel='Cumulative GHG impact (MtCO$_2$e)',
		xlabel='Time, years',ylim=[-4,2],xlim=[meta[pNam]['YearCurrent']-7,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_Cumu_FromYear' + str(meta[pNam]['YearCurrent']) + '_ByASET','png',900)#%%
	return

#%%
def Plot_EmissionsPerHectare_TS_Annual_FromCurrentYear_ByASET(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=meta[pNam]['YearCurrent']-15) & (tv<=2100) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_mu=np.zeros((tv.size,len(psL)))
	cnt=0
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]

		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		y_mu[:,cnt]=np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],y[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
		cnt=cnt+1

	plt.plot(tv[iT2],np.mean(y_mu[iT2],axis=1),'--',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ylim=[-40,200]
	ax.text(meta[pNam]['YearCurrent']-5,ylim[0]+15,'Net removals',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,ylim[0]-15,'Net emissions',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2010,2200,10),yticks=np.arange(-2000,2000,20),
		ylabel='Annual GHG impact (tCO$_2$e ha$^{-1}$)',xlabel='Time, years',
		ylim=ylim,xlim=[tv[iT2[0]],tv[iT2[-1]]])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsPerHectare_TS_Annual_FromYear' + str(meta[pNam]['YearCurrent']) + '_ByASET','png',900)#%%
	return

#%%
def Plot_EmissionsPerHectare_TS_Cumu_FromCurrentYear_ByASET(meta,mos,pNam,cNam,fsc):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=2010) & (tv<=2100) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_mu=np.zeros((tv.size,len(psL)))
	cnt=0
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]

		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		y_mu[:,cnt]=np.nan_to_num(y)

		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]

		plt.plot(tv[iT2],np.cumsum(y[iT2]),'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
		cnt=cnt+1

	plt.plot(tv[iT2],np.cumsum(np.mean(y_mu[iT2],axis=1)),'--',ms=3,color=[0,0,0],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')

	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-300,'Net removals',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,300,'Net emissions',fontsize=12,fontweight='bold',va='center',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2010,2200,10),yticks=np.arange(-2000,2000,50),ylabel='Cumulative GHG impact (tCO$_2$e ha$^{-1}$)',xlabel='Time, years',ylim=[-375,375],xlim=[meta[pNam]['YearCurrent']-7,2100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsPerHectare_TS_Cumu_FromYear' + str(meta[pNam]['YearCurrent']) + '_ByASET','png',900)#%%
	return

#%%
def Plot_Emissions_TS_Annual_FromCurrentYear_ByASET(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=2015) & (tv<=2100) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=1,color=[0,0,0])
	ax.plot([meta[pNam]['Project']['Year Project'],meta[pNam]['Project']['Year Project']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]

		#y=mos[pNam]['Delta'][cNam]['ByStrata'][operSpace]['E_NSB'][vStat][:,iPS,0,iYS]/1e6
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
		ind=np.where(meta['LUT']['Raw']['ASET']['Name']==nPS)[0][0]
		cl=[meta['LUT']['Raw']['ASET']['c1'][ind],meta['LUT']['Raw']['ASET']['c2'][ind],meta['LUT']['Raw']['ASET']['c3'][ind]]
		plt.plot(tv[iT2],y[iT2],'-',ms=3,color=cl,mec=cl,mfc='w',lw=2,mew=0.75,label=nPS)
	plt.plot(tv[iT2],ysum[iT2],'--',ms=3,color=[0.25,0.25,0.25],mec=cl,mfc='w',lw=2,mew=0.75,label='Total')
	ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(meta[pNam]['YearCurrent']-5,-0.25,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.text(meta[pNam]['YearCurrent']-5,0.25,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(2020,2150,5),yticks=np.arange(-20,20,0.1),ylabel='Annual GHG impact (MtCO$_2$e yr$^{-1}$)',
		xlabel='Time, years',ylim=[-0.6,0.6],xlim=[meta[pNam]['YearCurrent']-7,2075])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_Annual_FromYear' + str(meta[pNam]['YearCurrent']) + '_ByASET','png',900)#%%
	return

#%%
def Plot_Emissions_TS_Ann_Completed_ByASET(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	  'Underplanting','Replanting','Fill Planting',
	  'Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
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
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_Ann_Completed_ByASET','png',900)
	return

#%%
def Plot_Emissions_TS_Cumu_Completed_ByASET(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	vStat='Ensemble Mean'
	oper='Sum'
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,10))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([2030,2030],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-2000,2000],'k--',lw=0.5,color=[0,0,0])
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting',
	  'Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]

		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
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
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_Cumu_Completed_ByASET','png',900)
	return

#%%
def QA_CalcBenefitForEachStand(meta,pNam,mos):
	# Sum
	vL=['E_NEB','E_NSB','E_NAB']
	t0=1960
	t1=2050
	E_sum={}
	iScn=0;E_sum['Baseline']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Sum')
	iScn=1;E_sum['Actual']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Sum')
	E_sum['Delta']={}
	for k in E_sum['Baseline'].keys():
		E_sum['Delta'][k]=E_sum['Actual'][k]-E_sum['Baseline'][k]
	
	# Mean
	vL=['A','C_Biomass','C_Litter','C_Soil','C_DeadWood','C_G_Gross',
		'C_G_Net','C_M_Reg','C_M_Dist','C_LF','C_ToPileBurnTot','E_OpenBurning_ForestSector_Domestic',
		'E_NEB','E_NSB','E_NAB']
	#,'C_ToMillMerchGreen','C_ToMillNonMerchGreen','C_ToMillMerchDead','C_ToMillNonMerchDead'
	E_mu={}
	iScn=0;E_mu['Baseline']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Mean')
	iScn=1;E_mu['Actual']=cbu.Calc_MOS_Map(meta,pNam,iScn,t0=t0,t1=t1,Variables=vL,Operation='Mean')
	E_mu['Delta']={}
	for k in E_mu['Baseline'].keys():
		E_mu['Delta'][k]=E_mu['Actual'][k]-E_mu['Baseline'][k]
	return E_sum,E_mu

#%%
def Tabulate_SpeciesComposition(meta,pNam,t0,t1):
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
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Summary_SpeciesComposition_' + str(t0) + 'to' + str(t1) + '.pkl',d)
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
	df.to_excel(r'C:\Data\BCFCS\BCFCS_NOSE\Inputs\ScanForRestoration.xlsx')

	return

#%%
def Tabulate_UnderplantingByBurnSeverityClass(meta):
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
def Plot_AreaTreated_Frequency_BySpecies(meta,pNam,t0,t1):
	d=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Summary_SpeciesComposition_' + str(t0) + 'to' + str(t1) + '.pkl')
	
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
	plt.bar(bin,Np,0.85,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	# for i in range(bin.size):
	# 	y=N[i]
	# 	ax.text(bin[i],y+yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Frequency (%)',xticklabels=cd,xlim=[0.25,bin.size+0.75]) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_Frequency_BySpecies_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def Plot_AreaTreated_Frequency_ByBGC(meta,pNam,t0,t1):
	zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
	zBGC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')['Data']
	d=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].copy()
	for k in d.keys():
		d[k]=0
	for i in range(6):
		zYr=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_Year.tif')['Data']
		zFSC=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\PL_All_' + str(i+1) + '_SILV_FUND_SOURCE_CODE.tif')['Data']
		ikp=np.where( (zYr>=t0) & (zYr<=t1) & (np.isin(zFSC,meta['Param']['Raw']['FSC']['NO List ID'])==True) )
		if ikp[0].size==0:
			continue

		y=zBGC[ikp]
		u=np.unique(y)
		u=u[u>0]
		for iU in range(u.size):
			ind1=np.where( (y==u[iU]) )[0]
			cd=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			d[cd]=d[cd]+ind1.size
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\AreaTreated_Frequency_ByBGC_' + str(t0) + 'to' + str(t1) + '.pkl',d)

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
	plt.bar(bin,A,0.85,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	for i in range(bin.size):
		y=A[i]
		ax.text(bin[i],y*yo,str(np.round(A[i]/np.sum(A)*100,decimals=1)) + '%',fontsize=7,ha='center')
	ax.set(xticks=bin,ylabel='Implementation level (Kha)',xticklabels=cd,xlim=[0.25,bin.size+0.75],ylim=yl) # ,yticks=np.arange(0,110,10)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_Frequency_ByBGC_' + str(t0) + 'to' + str(t1),'png',900)
	return

#%%
def UpdateDatabase(meta,pNam,mos,mosWF,metaNM,pNamNM,mosNM,mosNM_WF):

	ac1='Silviculture investments'
	vStat='Ensemble Mean'
	operSpace='Sum'
	t_Start=1971

	# List of reporting years
	rep_yrL=[2030,2050,2070,2100]

	# Variable list
	vL=['E_NSB','V_ToMill_MerchTotal','Cost Total','Revenue Gross','Revenue Net']

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

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

	#--------------------------------------------------------------------------
	# Non obligation stand establishment
	#--------------------------------------------------------------------------

	cNam='NOSE'
	ac2='Non-obligation stand establishment'
	#psL=['All','Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	#  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'

	psL=meta[pNam]['Project']['Strata']['Project Type']['Unique CD']
	ysL=meta[pNam]['Project']['Strata']['Year']['Unique CD']
	ssL=['All']#meta[pNam]['Project']['Strata']['Spatial']['Unique CD']
	osL=['All','FIP'] #meta[pNam]['Project']['Strata']['Other']['Unique CD']

	cnt=0
	for ps in psL:
		for ss in ssL:
			for ys in ysL:
				for os in osL:
					print(ps + ' ' + ss + ' ' + ys + ' ' + os)

					iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==ps)[0]
					iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==ys)[0]
					iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==ss)[0]
					iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==os)[0]
	
					# Completed cumulative
					for yr in rep_yrL:
						d['Action Category 1'][cnt]=ac1
						d['Action Category 2'][cnt]=ac2
						d['Action Category 3'][cnt]=ps
						d['Funding Source Code'][cnt]=os
						d['Year Implemented'][cnt]=ys
						d['Comparison Name'][cnt]=cNam
						d['Scope Investments'][cnt]='Completed'
						d['Area Treated (ha/yr)'][cnt]=0
						d['Operation Spatial'][cnt]=operSpace
						d['Operation Temporal'][cnt]='Sum'
						d['Ensemble Statistic'][cnt]='Ensemble Mean'
						d['Time Start'][cnt]=t_Start
						d['Time End'][cnt]=yr
						for v in vL:
							y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)[operSpace][vStat]/meta[pNam]['Project']['Multi']
							if ys!='All':
								ind=np.where(tv<int(ys)-7)[0];
								y[ind]=0
							iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
							d[v][cnt]=np.round(np.nan_to_num(np.sum(y[iT])),decimals=1)
						cnt=cnt+1
	
					# Completed + CAP cumulative
					for yr in rep_yrL:
						d['Action Category 1'][cnt]=ac1
						d['Action Category 2'][cnt]=ac2
						d['Action Category 3'][cnt]=ps
						d['Funding Source Code'][cnt]=os
						d['Year Implemented'][cnt]=ys
						d['Comparison Name'][cnt]=cNam
						d['Scope Investments'][cnt]='Completed plus CAP'
						d['Area Treated (ha/yr)'][cnt]=0
						d['Operation Spatial'][cnt]=operSpace
						d['Operation Temporal'][cnt]='Sum'
						d['Ensemble Statistic'][cnt]='Ensemble Mean'
						d['Time Start'][cnt]=t_Start
						d['Time End'][cnt]=yr
						iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
						for v in vL:
							y=cbu.GetMosDeltaVar(meta,pNam,mosWF,cNam,v,iPS,iSS,iYS,iOS)[operSpace][vStat]/meta[pNam]['Project']['Multi']
							if ys!='All':
								ind=np.where(tv<int(ys)-7)[0];
								y[ind]=0
							iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
							d[v][cnt]=np.round(np.nan_to_num(np.sum(y[iT])),decimals=1)
						cnt=cnt+1

	#--------------------------------------------------------------------------
	# Add forest nutrient management
	#--------------------------------------------------------------------------

	ac2='Forest nutrient management'
	cNam='Moderate Harvest'
	psL=metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']
	ssL=['All']; #metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']
	ysL=metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']
	osL=['All','FIP'] #metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']

	for ps in psL:
		for ss in ssL:
			for ys in ysL:
				for os in osL:
					print(ps + ' ' + ss + ' ' + ys + ' ' + os)

					iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']==ps)[0]
					iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']==ys)[0]
					iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']==ss)[0]
					iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']==os)[0]

					# Completed cumulative
					for yr in rep_yrL:
						d['Action Category 1'][cnt]=ac1
						d['Action Category 2'][cnt]=ac2
						d['Action Category 3'][cnt]=ps
						d['Funding Source Code'][cnt]=os
						d['Year Implemented'][cnt]=ys
						d['Comparison Name'][cnt]=cNam
						d['Scope Investments'][cnt]='Completed'
						d['Area Treated (ha/yr)'][cnt]=0
						d['Operation Spatial'][cnt]=operSpace
						d['Operation Temporal'][cnt]='Sum'
						d['Ensemble Statistic'][cnt]='Ensemble Mean'
						d['Time Start'][cnt]=t_Start
						d['Time End'][cnt]=yr
						iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
						for v in vL:
							y=cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNam,v,iPS,iSS,iYS,iOS)[operSpace][vStat]/meta[pNam]['Project']['Multi']
							if ys!='All':
								ind=np.where(tv<int(ys))[0];
								y[ind]=0
							iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
							d[v][cnt]=np.round(np.nan_to_num(np.sum(y[iT])),decimals=1)
						cnt=cnt+1
	
					# Completed + climate action planning cumulative
					for yr in rep_yrL:
						d['Action Category 1'][cnt]=ac1
						d['Action Category 2'][cnt]=ac2
						d['Action Category 3'][cnt]=ps
						d['Funding Source Code'][cnt]=os
						d['Year Implemented'][cnt]=ys
						d['Comparison Name'][cnt]=cNam
						d['Scope Investments'][cnt]='Completed plus CAP'
						d['Area Treated (ha/yr)'][cnt]=0
						d['Operation Spatial'][cnt]=operSpace
						d['Operation Temporal'][cnt]='Sum'
						d['Ensemble Statistic'][cnt]='Ensemble Mean'
						d['Time Start'][cnt]=t_Start
						d['Time End'][cnt]=yr
						iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
						for v in vL:
							y=cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM_WF,cNam,v,iPS,iSS,iYS,iOS)[operSpace][vStat]/meta[pNam]['Project']['Multi']
							if ys!='All':
								ind=np.where(tv<int(ys))[0];
								y[ind]=0
							iT=np.where( (tv>=d['Time Start'][cnt]) & (tv<=d['Time End'][cnt]) )[0]
							d[v][cnt]=np.round(np.nan_to_num(np.sum(y[iT])),decimals=1)
						cnt=cnt+1

	# Truncate
	for k in d.keys():
		d[k]=d[k][0:cnt]

	# Save to file
	#df=pd.DataFrame.from_dict(d)
	#df.to_excel()
	pthout=r'C:\Users\rhember\Government of BC\External NRS Data Science and Modelling - General\Data\BC-FOR Climate Change Mitigation Database' + '\\R' + str(meta[pNam]['YearCurrent']+1) + '\\BC-FOR Climate Mitigation Database Summary.xlsx'
	gu.PrintDict(d,pthout,SheetName='Summary')

	return

#%%
def DefineStrata(meta,pNam,dmec,lsat):

	# By project type
	cd=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	id=np.array(list(meta['LUT']['Derived']['ASET'].values()))

	# Only keep a subset to avoid crashing in MOS calculation
	ind=np.where( (np.isin(cd,['Salvage and Planting Post Beetle',
							'Salvage and Planting Post Other',
							'Knockdown and Planting',
							'Underplanting',
							'Fill Planting',
							'Replanting'])==True) )[0]
	cd=cd[ind]
	id=id[ind]

	meta[pNam]['Project']['Strata']['Project Type']['Unique ID']=np.append(0,id)
	meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=np.append('All',cd)
	#meta[pNam]['Project']['Strata']['Project Type']['ID']=np.zeros(meta[pNam]['Project']['N Stand'],dtype='int8') # This is populated by inventory processing

	# By time (last year of implementation)
	#t=np.arange(meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']+1,1)
	t=np.arange(meta[pNam]['YearCurrent']-6,meta[pNam]['YearCurrent']+1,1)
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

	# By space (BGC zone)
	u=np.unique(lsat['ID_BGCZ'])
	cd=np.array([])
	for iU in range(u.size):
		 cd=np.append(cd,u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU]))

	cd=np.array(['IDF','SBPS','SBS','ICH','MS','ESSF','BWBS','CWH'])
	id=np.zeros(cd.size)
	for j in range(id.size):
		 id[j]=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][cd[j]]

	meta[pNam]['Project']['Strata']['Spatial']['Unique ID']=np.append(0,id)
	meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=np.append('All',cd)
	meta[pNam]['Project']['Strata']['Spatial']['ID']=lsat['ID_BGCZ']

	# By Other (FSC)
	cd=np.array(['FTM','FRP','FCE','FIP','FCM','LFP','FTL','FES','GA','FED'])
	id=np.zeros(cd.size)
	for i,code in enumerate(cd):
		id[i]=meta['LUT']['RSLT_ACTIVITY_TREATMENT_SVW']['SILV_FUND_SOURCE_CODE'][code]
	meta[pNam]['Project']['Strata']['Other']['Unique CD']=np.append('All',cd)
	meta[pNam]['Project']['Strata']['Other']['Unique ID']=np.append(0,id)
	iScn=1
	for iStand in range(meta[pNam]['Project']['N Stand']):
		ind=np.where( (dmec[iScn][iStand]['Index to Event Inciting NOSE']!=9999) )[0]
		if ind.size==0:
			continue
		if ind.size>1:
			ind=ind[-1]
		meta[pNam]['Project']['Strata']['Other']['ID'][iStand]=dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][ind]

	return meta

#%%
def Plot_Emissions_TS_AnnAndCumu_CurrentYear(meta,mos,pNam,cNam,fsc):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-10) )[0]
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		#iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]

		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1965,2100],ylim=[-0.5,0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(tv[iT2[0]]+2,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(tv[iT2[0]]+2,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y2=np.cumsum(ysum)
	ax2.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,2),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-10,10])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_CurrentYear','png',900)
	return

#%%
def Plot_Emissions_TS_AnnAndCumu_CurrentYear_NOSE_AND_FNM(meta,mos,pNam,cNam,pNamNM,metaNM,mosNM):
	v='E_NSB'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
	cl3=meta['Graphics']['Colours']['rgb']['Purple Light']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Non-obligation stand establishment
	fsc='FIP'
	cNam='NOSE'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_ann=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]

		y0=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		y_ann=y_ann+np.nan_to_num(y0)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	y_ann[iT1]=0

	# Forest nutrient management
	cNamNM='Moderate Harvest'
	iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
	iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']==str(metaNM[pNamNM]['YearCurrent']))[0][0]
	iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
	y0=cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNamNM,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/metaNM[pNamNM]['Project']['Multi']
	iT1=np.where( (tv<metaNM[pNamNM]['YearCurrent']) )[0]
	y0[iT1]=0
	y_ann_all=y_ann+y0

	ax.plot(tv,y_ann,'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1971,2100],ylim=[-0.6,0.6])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y_cumu=np.cumsum(y_ann)
	y_cumu_all=np.cumsum(y_ann_all)
	ax2.plot(tv,y_cumu,'-.',ms=3,color=cl3,lw=2,label='Just NOSE')
	ax2.plot(tv,y_cumu_all,'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,2),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1971,2100],ylim=[-10,10])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_CurrentYear_NOSE_AND_FNM','png',900)
	return

#%%
def Plot_Emissions_TS_AnnAndCumu_Completed_NOSE_AND_FNM(meta,mos,pNam,cNam,pNamNM,metaNM,mosNM):
	v='E_NSB'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
	cl3=meta['Graphics']['Colours']['rgb']['Purple Light']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Non-obligation stand establishment
	fsc='All'
	cNam='NOSE'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_ann=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]

		y0=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		y_ann=y_ann+np.nan_to_num(y0)
	#iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	#y_ann[iT1]=0

	# Forest nutrient management
	cNamNM='Moderate Harvest'
	iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
	iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
	iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
	y0=cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNamNM,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/metaNM[pNamNM]['Project']['Multi']
	#iT1=np.where( (tv<metaNM[pNamNM]['YearCurrent']) )[0]
	#y0[iT1]=0
	y_ann_all=y_ann+y0

	ax.plot(tv,y_ann,'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.5),xlabel='Time, years',xlim=[1971,2100],ylim=[-3,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y_cumu=np.cumsum(y_ann)
	y_cumu_all=np.cumsum(y_ann_all)
	ax2.plot(tv,y_cumu,'-.',ms=3,color=cl3,lw=2,label='Just NOSE')
	ax2.plot(tv,y_cumu_all,'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,50),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1971,2100],ylim=[-250,250])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_Completed_NOSE_AND_FNM','png',900)
	return

#%%
def Plot_Emissions_TS_AnnAndCumu_CompletedPlusCAP_NOSE_AND_FNM(meta,mos,pNam,cNam,pNamNM,metaNM,mosNM):
	v='E_NSB'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
	cl3=meta['Graphics']['Colours']['rgb']['Purple Light']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	# Non-obligation stand establishment
	fsc='All'
	cNam='NOSE'
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	y_ann=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]

		y0=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		y_ann=y_ann+np.nan_to_num(y0)
	#iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	#y_ann[iT1]=0

	# Forest nutrient management
	cNamNM='Moderate Harvest'
	iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
	iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
	iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
	y0=cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNamNM,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/metaNM[pNamNM]['Project']['Multi']
	#iT1=np.where( (tv<metaNM[pNamNM]['YearCurrent']) )[0]
	#y0[iT1]=0
	y_ann_all=y_ann+y0

	ax.plot(tv,y_ann,'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,2),xlabel='Time, years',xlim=[1971,2100],ylim=[-10,10])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y_cumu=np.cumsum(y_ann)
	y_cumu_all=np.cumsum(y_ann_all)
	ax2.plot(tv,y_cumu,'-.',ms=3,color=cl3,lw=2,label='Just NOSE')
	ax2.plot(tv,y_cumu_all,'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-2000,2000,100),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1971,2100],ylim=[-500,500])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_CompletedPlusCAP_NOSE_AND_FNM','png',900)
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
		y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.text(1967,-1.6,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,1.6,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1965,2100],ylim=[-2.2,2.2])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,0,0,iYS]/1e6)
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
def Plot_Emissions_TS_AnnAndCumu_CompletedAndCAP(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

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
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,1),xlabel='Time, years',
		xlim=[1965,2100],ylim=[-9,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1967,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,100),
		 ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-450,150])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_CompletedAndCAP','png',900)
	return

#%%
def Plot_Emissions_TS_AnnAndCumu_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

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
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]

		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,1),xlabel='Time, years',
		xlim=[1965,2100],ylim=[-9,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1967,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1967,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	#y=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,0,0,iYS]/1e6)
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,100),
		 ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1965,2100],ylim=[-450,150])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_Completed','png',900)
	return

#%%
def Plot_Yield_TS_AnnAndCumu_CurrentYear(meta,mos,pNam,cNam,fsc):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	iT2=np.where( (tv>=1971) & (tv<=2100) )[0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'V_ToMill_MerchTotal',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,10),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[tv[iT2[0]],2100],ylim=[-0.4,0.4])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y2=np.cumsum(ysum)
	ax2.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2170,10),yticks=np.arange(-2,2,0.5),ylabel='Cumulative',xlim=[tv[iT2[0]],2100],ylim=[-1.5,1.5])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Yield_TS_AnnAndCumu_CurrentYear','png',900)
	return

#%%
def Plot_Yield_TS_AnnAndCumu_Completed(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1971) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Salvage and Planting Post Other',
	  'Knockdown and Planting','Road Rehabilitation',
	  'Straight-to-planting Post Other',
	  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'

	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'V_ToMill_MerchTotal',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][:,iPS,0,iYS]/1e6
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	#ax.legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1971,2100],ylim=[-1,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,20),ylabel='',xlim=[1971,2100],ylim=[-66.66,200])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Yield_TS_AnnAndCumu_Completed','png',900)
	return

#%%
def Plot_Yield_TS_AnnAndCumu_CompletedAndCAP(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1965) & (tv<=2100) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Salvage and Planting Post Other',
	  'Knockdown and Planting','Road Rehabilitation',
	  'Straight-to-planting Post Other',
	  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'V_ToMill_MerchTotal',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')

	ax.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,0.5),xlabel='Time, years',xlim=[1971,2100],ylim=[-1,3])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual yield impact (Mm$^3$ yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(1973,y_min1+0.06*dy,'Losses',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(1973,y_max1-0.06*dy,'Gains',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	ysumc=np.cumsum(ysum)
	ax2.plot(tv[iT2],ysumc[iT2],'-.',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2200,10),yticks=np.arange(-2000,2000,20),ylabel='',xlim=[1971,2100],ylim=[-66.66,200])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative yield impact (Mm$^3$)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Yield_TS_AnnAndCumu_CompletedAndCAP','png',900)
	return

#%%
def ProjectFuture(meta,pNam,mos,TimeHorizon,fracPT):
	# This function takes the current year per-hectare time series and applies it into the future at a projected AIL
	# Notes: Cumulative values are incorrect - do not use!

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]

	nS=[]
	for k in meta[pNam]['Project']['Strata'].keys():
		nS.append(meta[pNam]['Project']['Strata'][k]['Unique ID'].size)

	iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]

	dCP=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_CostsAndPrices.xlsx',sheet_name='Sheet1',skiprows=0)

	# Proportion of activity levels
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'

	vL=['E_NSB','V_ToMill_MerchTotal']#,'Cost Total','Revenue Net','Cost Nutrient Management']
	#vL=['E_NSB']

	mosWF=copy.deepcopy(mos)

	for cNam in mos[pNam]['Delta'].keys():
		for v in vL:

			# Force completed to be zero
			if TimeHorizon=='SP Projection':
				mosWF[pNam]['Delta'][cNam]['Data']['Sum'][v]['Ensemble Mean']=0*np.ones(mosWF[pNam]['Delta'][cNam]['Data']['Sum'][v]['Ensemble Mean'].shape)

			for nPS in psL:

				y2=np.zeros((tv.size,nS[0],nS[1],nS[2],nS[3]))
				iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
				ail_pt=meta[pNam]['AIL CAP']*fracPT[nPS]

				# Get mean emissions for specific project type and current year
				y0=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,0,iYS,0)['Mean']['Ensemble Mean']

				# Define future years of projection
				if TimeHorizon=='Full':
					FutureYears=np.arange(meta[pNam]['YearCurrent']+1,tv[-1],1)
				elif TimeHorizon=='SP Projection':
					FutureYears=np.arange(meta[pNam]['YearCurrent']+1,meta[pNam]['YearCurrent']+2,1)

				y1=np.zeros((tv.size,FutureYears.size))
				for i,yr in enumerate(FutureYears):
					iT1=np.where( (tv>=yr) )[0]
					iT2=np.where( (tv>=meta[pNam]['YearCurrent']-10) )[0]
					iT2=iT2[0:iT1.size]
					y1[iT1-10,i]=ail_pt*y0[iT2]
				y1=np.sum(y1,axis=1)
				y2[:,iPS,0,0,0]=y1

				mosWF[pNam]['Delta'][cNam]['Data']['Sum'][v]['Ensemble Mean']=mosWF[pNam]['Delta'][cNam]['Data']['Sum'][v]['Ensemble Mean']+y2[ mos[pNam]['MOS Index'] ]

		# Fix costs
		flg=0
		if flg==1:
			iT=np.where( (tv>=meta[pNam]['YearCurrent']) )[0]
			iTe=np.where(dCP['Year']==meta[pNam]['YearCurrent'])[0]
			y=meta[pNam]['AIL CAP']*(dCP['Cost Nutrient Purchase (CAD/ha)'][iTe]+dCP['Cost Nurtrient Application (CAD/ha)'][iTe]+dCP['Cost Nutrient Overhead (CAD/ha)'][iTe])
			mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Cost Total']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Cost Total']['Ensemble Mean'][iT,:,:,:]+y
			mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Cost Nutrient Management']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Cost Nutrient Management']['Ensemble Mean'][iT,:,:,:]+y
			mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Revenue Net']['Ensemble Mean'][iT,:,:,:]=mosWF[pNam]['Delta'][cNam]['Data']['Sum']['Revenue Net']['Ensemble Mean'][iT,:,:,:]-y

	return mosWF


#%%
def Plot_Emissions_TS_AnnAndCumu_CompletedAndCAP_WithScenarios(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT2=np.where( (tv>=1990) & (tv<=2050) )[0]

	cl1=[0.27,0.49,0.77]
	cl2=[0.5,0.8,0]
	xl=[1990,2050]

	plt.close('all'); fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(11,11))
	ax[0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	#ax[0].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-20,20],'k--',lw=0.5,color=[0,0,0])
	#ax[0].plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	#ax[0].plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	#ax[0].plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation',
	  'Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']

	ysumB=np.zeros(tv.size)
	ysumP=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		yB=cbu.GetMosScnVar(meta,pNam,mos,iB,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		yP=cbu.GetMosScnVar(meta,pNam,mos,iP,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysumB=ysumB+np.nan_to_num(yB)
		ysumP=ysumP+np.nan_to_num(yP)
	#y=mos[pNam]['Delta'][cNam]['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][:,0,0,iYS]/1e6
	ax[0].plot(tv[iT2],ysumB[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Baseline')
	ax[0].plot(tv[iT2],ysumP[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Actual')
	ax[0].set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,5),xlabel='Time, years',xlim=xl,ylim=[-5,20])
	ax[0].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0].set_ylabel('GHG emissions\n(MtCO$_2$e yr$^{-1}$)',color='k',weight='normal',fontsize=8)

	ax[1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ax[1].plot(tv[iT2],ysum[iT2],'-',ms=3,color=[0.5,0,1],mec=cl2,mfc='w',lw=2,mew=0.75,label='Annual')
	ax[1].set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,2),xlabel='Time, years',xlim=xl,ylim=[-6,2])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1].set_ylabel('Total direct GHG\nimpact (MtCO$_2$e yr$^{-1}$)',color='k',weight='normal',fontsize=8)

	ax[2].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax[2].plot(tv[iT2],ysum[iT2],'-',ms=3,color=[0.5,0,1],mec=cl2,mfc='w',lw=2,mew=0.75,label='Total direct impact')
	iRL=np.where( (tv>=2008) & (tv<=2017) )
	RL=np.ones(tv.size)*np.mean(ysum[iRL])
	iT3=np.where( (tv>=2017) )
	ax[2].plot(tv[iT2],RL[iT2],'--',ms=3,color=[0.5,0.5,0.5],mec=cl2,mfc='w',lw=2,mew=0.75,label='Reference level (2008-2017 mean)')
	ax[2].fill_between(tv[iT3],ysum[iT3],RL[iT3],color=[1,1,0.5],alpha=0.3,lw=0,label='Incremental impact')
	ax[2].set(xticks=np.arange(1960,2200,10),yticks=np.arange(-200,200,2),xlabel='Time, years',xlim=xl,ylim=[-6,2])
	ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[2].set_ylabel('Incremental GHG\nimpact (MtCO$_2$e yr$^{-1}$)',color='k',weight='normal',fontsize=8)
	ax[2].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
	gu.axletters(ax,plt,0.02,0.85,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_CompletedAndCAP_WithScenarios','png',900)
	return

#%%
def QA_Plot_Validation_NetGrowth_ByBGC(meta,mos,pNam):
	dO=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Observation_NetGrowth.pkl')
	iB=0
	iP=1
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']+20) & (tv<=meta[pNam]['YearCurrent']+40) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
	iPS=0
	iOS=0
	v='C_G_Net'
	xlab=['Observations','Model\nBaseline','Model\nAction']
	lab=[]
	tym=1.05
	tfs=10

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Light']
	cl3=meta['Graphics']['Colours']['rgb']['Green Dark']
	plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(20,9));

	zone='IDF'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,0].bar(1,dO['Ctot Net Mean'][iObs],facecolor=cl1,label='Observations')
	ax[0,0].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,0].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,0].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,0].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,0].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='SBS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,1].bar(1,yO,facecolor=cl1,label='Observations')
	ax[0,1].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,1].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,1].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,1].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='SBPS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,2].bar(1,yO,facecolor=cl1,label='Observations')
	ax[0,2].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,2].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,2].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,2].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='ESSF'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,0].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,0].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,0].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,0].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,0].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,0].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='MS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,1].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,1].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,1].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,1].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,1].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='ICH'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,2].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,2].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,2].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,2].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,2].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.045,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.05)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Validation_NetGrowth_ByBGC','png',900)
	return

#%%
def QA_Plot_Validation_NetGrowth_ByProjectType(meta,mos,pNam):
	dO=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Observation_NetGrowth.pkl')
	iB=0
	iP=1
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']+20) & (tv<=meta[pNam]['YearCurrent']+40) )[0]
	iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
	v='C_G_Net'
	xlab=['Observations','Model\nBaseline','Model\nAction']
	lab=[]
	tym=1.05
	tfs=10

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Light']
	cl3=meta['Graphics']['Colours']['rgb']['Green Dark']
	plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(20,9));

	pt='Underplanting'; lab.append(pt)
	iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==pt)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,0].bar(1,dO['Ctot Net Mean'][iObs],facecolor=cl1,label='Observations')
	ax[0,0].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,0].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,0].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,0].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,0].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='SBS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,1].bar(1,yO,facecolor=cl1,label='Observations')
	ax[0,1].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,1].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,1].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,1].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='SBPS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[0,2].bar(1,yO,facecolor=cl1,label='Observations')
	ax[0,2].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[0,2].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[0,2].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,2].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[0,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='ESSF'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,0].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,0].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,0].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,0].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,0].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,0].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='MS'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,1].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,1].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,1].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,1].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,1].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	zone='ICH'; lab.append(zone)
	iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
	yB=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iB,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	yP=np.mean(cbu.GetMosScnVar(meta,pNam,mos,iP,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	iObs=np.where(dO['Zone']==zone)[0][0]
	yO=dO['Ctot Net Mean'][iObs]
	ax[1,2].bar(1,yO,facecolor=cl1,label='Observations')
	ax[1,2].bar(2,yB,0.8,facecolor=cl2,label='Model Baseline')
	ax[1,2].bar(3,yP,0.8,facecolor=cl3,label='Model Actual')
	ax[1,2].text(2,tym*yB,str(int((yB-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,2].text(3,tym*yP,str(int((yP-yO)/yO*100)) + '%',ha='center',fontsize=tfs)
	ax[1,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xticks=[1,2,3],xticklabels=xlab,ylim=[0,1.2*yO])
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.055,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.035)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Validation_NetGrowth_ByProjectType','png',900)
	return

#%%
def Plot_Emissions_Cumu2050_ByFSC(meta,mos,pNam,cNam):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=2015) & (tv<=2050) )[0]

	bw=0.7
	FSC=meta[pNam]['Project']['Strata']['Other']['Unique CD']
	bin=np.arange(1,FSC.size+1,1)
	y=np.zeros(bin.size)
	for i,fsc in enumerate(meta[pNam]['Project']['Strata']['Other']['Unique CD']):
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
		y[i]=np.sum(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'][iT])/meta[pNam]['Project']['Multi']
	ord=np.argsort(y)
	y=y[ord]
	FSC=FSC[ord]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,7));
	ax.bar(bin,y,bw,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	ax.plot([-1,bin.size+1],[0,0],'k-')
	ax.set(xticks=bin,xticklabels=FSC,ylabel='Cumulative GHG impact 2050\n(MtCO$_2$e)',xlim=[0.5,bin.size+0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_Cumu2050_ByFSC','png',200)
	return

#%%
def Plot_EmissionsPerHa_Cumu2050_ByFSC(meta,mos,pNam,cNam):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=2015) & (tv<=2050) )[0]

	bw=0.7
	FSC=meta[pNam]['Project']['Strata']['Other']['Unique CD']
	bin=np.arange(1,FSC.size+1,1)
	y=np.zeros(bin.size)
	for i,fsc in enumerate(meta[pNam]['Project']['Strata']['Other']['Unique CD']):
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
		y[i]=np.sum(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	ord=np.argsort(y)
	y=y[ord]
	FSC=FSC[ord]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,7));
	ax.bar(bin,y,bw,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	ax.plot([-1,bin.size+1],[0,0],'k-')
	ax.set(xticks=bin,xticklabels=FSC,ylabel='Cumulative GHG impact 2050\n(tCO$_2$e ha$^{-1}$)',xlim=[0.5,bin.size+0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsPerHa_Cumu2050_ByFSC','png',200)
	return

#%%
def Plot_EmissionsPerHectare_CurrentYear_ByBGC(meta,mos,pNam,cNam,fsc):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=2015) & (tv<=2050) )[0]

	bw=0.7
	Zones=meta[pNam]['Project']['Strata']['Spatial']['Unique CD']
	bin=np.arange(1,Zones.size+1,1)
	y=np.zeros(bin.size)
	for i,zone in enumerate(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']):
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']==zone)[0][0]
		iYS=np.where( meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']) )[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
		y[i]=np.sum(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT])
	ord=np.argsort(y)
	y=y[ord]
	Zones=Zones[ord]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,7));
	ax.bar(bin,y,bw,facecolor=meta['Graphics']['Colours']['rgb']['Blue Dark'])
	ax.plot([-1,bin.size+1],[0,0],'k-')
	#ax.legend(loc='lower left',labelspacing=0.5,facecolor=[1,1,1],frameon=False,ncols=1,fontsize=6);
	ax.set(xticks=bin,xticklabels=Zones,ylabel='Cumulative GHG impact 2050 (tCO$_2$e ha$^{-1}$)',xlim=[0.5,bin.size+0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\EmissionsPerHectare_CurrentYear_ByBGCZone','png',200)
	return

#%%
def Plot_Fluxes_TS_Annual_FromCurrentYear_ByASET(meta,mos,pNam,cNam,fsc):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	ptL=['Salvage and Planting Post Beetle','Straight-to-planting Post Other',
	  'Knockdown and Planting','Road Rehab','Underplanting','Fill Planting', 'Replanting']
	for pt in meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==pt)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
	
		xlim=[meta[pNam]['YearCurrent']-10,2150]
		vL=['E_NEE_ForestSector_Total','E_Wildfire_ForestSector_Total','E_OpenBurning_ForestSector_Total','E_ForestryOps_Total','E_HWP_ForestSector_Total','E_NSB']
		labL=['NEE','wildfire','open burning','operations','products','net emissions']
		cnt=0
		plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,10))
		for i in range(2):
			for j in range(3):
				ax[i,j].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
				ax[i,j].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
				#ax[i,j].plot([2030,2030],[-200,200],'k--',lw=0.5,color=[0,0,0])
				ax[i,j].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
				y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,vL[cnt],iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
				if (np.sum(y)==0) | (np.sum(np.isnan(y))>0):
					continue
				ax[i,j].plot(tv,y,'-',color=meta['Graphics']['Colours']['rgb']['Blue Dark'],lw=1.25,ms=3)
				ylim=[np.min(y),1.1*np.max(y)]
				if ylim[0]<0:
					ylim[0]=1.1*ylim[0]
				else:
					ylim[0]=0.9*ylim[0]
				ax[i,j].set(xticks=np.arange(2020,2150,20),ylabel='$\Delta$ ' + labL[cnt] + ' (tCO$_2$e ha$^{-1}$ yr$^{-1}$)',
				   xlabel='Time, years',ylim=ylim,xlim=xlim)
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
				cnt=cnt+1
		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		if (np.sum(y)==0) | (np.sum(np.isnan(y))>0):
			continue
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_TS_Annual_FromYear' + str(meta[pNam]['YearCurrent']) + '_' + pt,'png',900)#%%
	return

#%%
def Plot_Fluxes_TS_Cumulative_FromCurrentYear_ByASET(meta,mos,pNam,cNam,fsc):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	for pt in meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==pt)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
	
		xlim=[meta[pNam]['YearCurrent']-10,2150]
		vL=['E_NEE_ForestSector_Total','E_Wildfire_ForestSector_Total','E_OpenBurning_ForestSector_Total','E_ForestryOps_Total','E_HWP_ForestSector_Total','E_NSB']
		labL=['NEE','wildfire','open burning','operations','products','net emissions']
		cnt=0
		plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(22,10))
		for i in range(2):
			for j in range(3):
				ax[i,j].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
				ax[i,j].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
				#ax[i,j].plot([2030,2030],[-200,200],'k--',lw=0.5,color=[0,0,0])
				ax[i,j].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
				y=np.cumsum(np.nan_to_num(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,vL[cnt],iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'])/meta[pNam]['Project']['Multi'])
				if (np.sum(y)==0) | (np.sum(np.isnan(y))>0):
					continue
				ax[i,j].plot(tv,y,'-',color=meta['Graphics']['Colours']['rgb']['Blue Dark'],lw=1.25,ms=3)
				ylim=[np.min(y),1.1*np.max(y)]
				if ylim[0]<0:
					ylim[0]=1.1*ylim[0]
				else:
					ylim[0]=0.9*ylim[0]
				ax[i,j].set(xticks=np.arange(2020,2150,20),ylabel='$\Delta$ ' + labL[cnt] + ' (MtCO$_2$e)',
				   xlabel='Time, years',ylim=ylim,xlim=xlim)
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
				cnt=cnt+1
		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		if (np.sum(y)==0) | (np.sum(np.isnan(y))>0):
			continue
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_TS_Cumulative_FromYear' + str(meta[pNam]['YearCurrent']) + '_' + pt,'png',900)#%%
	return

#%%
def Plot_NetGrowth_TS_Annual_FromCurrentYear_ByASET(meta,mos,pNam,cNam,fsc):

	for nPS in meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:

		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']==fsc)[0][0]
		tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
		iT=np.where( (tv>=meta[pNam]['YearCurrent']-10) & (tv<=2100) )[0]
		
		cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
		cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
		lw=1.5
		
		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8))
		v='C_Biomass'
		ax[0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
		ax[0].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
		ax[0].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
		iScn=mos[pNam]['Delta'][cNam]['iB'];
		y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		if np.isnan(np.sum(y[iT]))==True:
			continue
		ax[0].plot(tv[iT],y[iT],'b-',color=cl1,lw=lw,label='Baseline scenario')

		iScn=mos[pNam]['Delta'][cNam]['iP'];
		y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[0].plot(tv[iT],y[iT],'g--',color=cl2,lw=lw,label='Action scenario')

		ax[0].set(xticks=np.arange(2020,2150,20),ylabel='Biomass (tC ha$^{-1}$)',
		   xlabel='Time, years',ylim=[0,150],xlim=[tv[iT[0]],tv[iT[-1]]])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax[0].legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
		
		v='C_G_Net'
		ax[1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
		ax[1].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
		ax[1].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
		iScn=mos[pNam]['Delta'][cNam]['iB']; y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[1].plot(tv[iT],y[iT],'b-',color=cl1,lw=lw)
		iScn=mos[pNam]['Delta'][cNam]['iP']; y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[1].plot(tv[iT],y[iT],'g--',color=cl2,lw=lw)
		ax[1].set(xticks=np.arange(2020,2150,20),ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',
		   xlabel='Time, years',ylim=[-10,3.5],xlim=[tv[iT[0]],tv[iT[-1]]])
		ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\NetGrowth_TS_Ann_CurrentYear_' + nPS,'png',900)
	return

#%%
def Plot_Age_TS_Annual_FromCurrentYear_ByASET(meta,mos,pNam,cNam):
	#nPS='Underplanting'
	#nPS='Salvage and Planting Post Beetle'
	for nPS in meta[pNam]['Project']['Strata']['Project Type']['Unique CD']:

		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='FIP')[0][0]
		tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
		iT=np.where( (tv>=meta[pNam]['YearCurrent']-50) & (tv<=2050) )[0]
		
		cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
		cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
		lw=1.5
		
		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8))
		v='A'
		ax[0].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
		ax[0].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
		ax[0].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
		iScn=mos[pNam]['Delta'][cNam]['iB'];
		y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		if np.isnan(np.sum(y[iT]))==True:
			continue
		ax[0].plot(tv[iT],y[iT],'b-',color=cl1,lw=lw)

		iScn=mos[pNam]['Delta'][cNam]['iP'];
		y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[0].plot(tv[iT],y[iT],'g--',color=cl2,lw=lw)

		ax[0].set(xticks=np.arange(1800,2150,10),ylabel='Age (years)',
		   xlabel='Time, years',ylim=[0,150],xlim=[tv[iT[0]],tv[iT[-1]]])
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		v='C_G_Net'
		ax[1].plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.8,0.8,0.8])
		ax[1].plot([meta[pNam]['YearCurrent'],meta[pNam]['YearCurrent']],[-200,200],'k--',lw=0.75,color=[0,0,0])
		ax[1].plot([2050,2050],[-200,200],'k--',lw=0.75,color=[0,0,0])
		iScn=mos[pNam]['Delta'][cNam]['iB']; y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[1].plot(tv[iT],y[iT],'b-',color=cl1,lw=lw)
		iScn=mos[pNam]['Delta'][cNam]['iP']; y=cbu.GetMosScnVar(meta,pNam,mos,iScn,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean']
		ax[1].plot(tv[iT],y[iT],'g--',color=cl2,lw=lw)
		ax[1].set(xticks=np.arange(1800,2150,10),ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',
		   xlabel='Time, years',ylim=[-10,3.5],xlim=[tv[iT[0]],tv[iT[-1]]])
		ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Age_TS_Ann_CurrentYear_' + nPS,'png',900)
	return

#%%
def AreaAffectedByEvents_ByASET(meta,pNam,mos):
	ivlT=1
	iScn=1
	
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['YearCurrent']-40) & (tv<=meta[pNam]['YearCurrent']+100) )[0]
	
	ptL=['Underplanting','Salvage and Planting Post Beetle','Knockdown and Planting','Fill Planting','Replanting']
	for pt in ptL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==pt)[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='FIP')[0][0]
		ufcs.Area_DisturbedAndManaged(meta,pNam,mos,ivlT,iScn,iT,iPS,iSS,iYS,iOS,
									ncols=2,
									LegendLoc='upper right',
									multi=1e3)
	return

#%%
def Plot_AreaTreated_TS_Ann_ByFSC(meta,pNam):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTimeAndFSC.pkl')
	
	cl1=meta['Graphics']['Colours']['rgb']['Blue Light']
	#cl2=meta['Graphics']['Colours']['rgb']['Green Neon']
	bw=0.8

	ptL=['Salvage and Planting Post Beetle','Salvage and Planting Post Other','Knockdown and Planting',
	  'Underplanting','Straight-to-planting Post Other','Road Rehabilitation','Direct Seeding',
	  'Fill Planting','Replanting','Ecosystem Restoration']

	for fsc in d['Data'].keys():
		cnt=0
		ysum=0
		plt.close('all'); fig,ax=plt.subplots(5,2,figsize=gu.cm2inch(16,12))
		for i in range(5):
			for j in range(2):
				y=d['Data'][fsc][ ptL[cnt] ]/1e3
				ysum=ysum+np.sum(y)
				ax[i,j].bar(d['tv'],y,bw,facecolor=cl1)
				ax[i,j].set(xticks=np.arange(1960,2100,1),xlabel='Time, years',ylabel='Treatment area\n(Kha/year)',
					xlim=[2017.25,2025.75])
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
				cnt=cnt+1
		if ysum==0:
			continue
		plt.tight_layout()
		gu.axletters(ax,plt,0.02,0.81,FontColor=meta['Graphics']['gp']['cla'],
					 LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],
					 FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],
					 Labels=ptL)
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaTreated_TS_Ann_' + fsc,'png',900)
	return

#%%
def Plot_Emissions_TS_AnnAndCumu_ProjectionYear(meta,mos,pNam,cNam):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT1=np.where( (tv<meta[pNam]['YearCurrent']-9) )[0]
	iT2=np.where( (tv>=1971) & (tv<=2100) )[0]
	yr_Proj=meta[pNam]['YearCurrent']+1

	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Green Neon']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9))
	ax.plot(tv,np.zeros(tv.size),'k-',lw=2,color=[0.85,0.85,0.85])
	ax.plot([yr_Proj,yr_Proj],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2030,2030],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2050,2050],[-20,20],'k--',lw=0.5,color=[0,0,0])
	ax.plot([2070,2070],[-20,20],'k--',lw=0.5,color=[0,0,0])

	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Road Rehabilitation','Underplanting','Replanting','Fill Planting','Ecosystem Restoration'] #'Harvest and Planting NSR Backlog'
	ysum=np.zeros(tv.size)
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=cbu.GetMosDeltaVar(meta,pNam,mos,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
		ysum=ysum+np.nan_to_num(y)
	ysum[iT1]=0
	ax.plot(tv[iT2],ysum[iT2],'-',ms=3,color=cl1,mec=cl1,mfc='w',lw=2,mew=0.75,label='Annual')
	ax.set(xticks=np.arange(1960,2100,5),yticks=np.arange(-20,20,0.1),xlabel='Time, years',xlim=[1971,2100],ylim=[-0.5,0.5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.set_ylabel('Annual GHG impact (MtCO$_2$e yr$^{-1}$)',color=cl1,weight='bold',fontsize=9)
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(tv[iT2[0]]+2,y_min1+0.06*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(tv[iT2[0]]+2,y_max1-0.06*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	
	ax2=ax.twinx()
	y2=np.cumsum(ysum)
	ax2.plot(tv[iT2],y2[iT2],'--',ms=3,color=cl2,mec=cl2,mfc='w',lw=2,mew=0.75,label='Cumulative')
	ax2.set(xticks=np.arange(1960,2100,5),yticks=np.arange(-2000,2000,2),ylabel='Cumulative GHG impact (MtCO$_2$e)',xlim=[1971,2100],ylim=[-10,10])
	ax2.yaxis.set_ticks_position('right'); ax2.xaxis.set_ticks_position('both'); ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2.set_ylabel('Cumulative GHG impact (MtCO$_2$e)', color=cl2,weight='bold',fontsize=9)

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_TS_AnnAndCumu_ProjectionYear','png',900)
	return

#%%
def Calculate_FractionOverstoryRemoval(meta,pNam,mos):
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\ASET_SummaryByTimeAndFSC.pkl')
	iT=np.where(d['tv']==meta[pNam]['YearCurrent'])[0]
	dS={}
	dS['N OverstoryRemoval']=d['Data']['FIP']['Salvage and Planting Post Beetle'][iT]+d['Data']['FIP']['Salvage and Planting Post Other'][iT]+d['Data']['FIP']['Knockdown and Planting'][iT]
	dS['N Total']=d['Data']['FIP']['Salvage and Planting Post Beetle'][iT]+d['Data']['FIP']['Salvage and Planting Post Other'][iT]+ \
		d['Data']['FIP']['Knockdown and Planting'][iT]+d['Data']['FIP']['Underplanting'][iT]+d['Data']['FIP']['Road Rehabilitation'][iT]+ \
		d['Data']['FIP']['Fill Planting'][iT]+d['Data']['FIP']['Replanting'][iT]
	dS['Frac OverstoryRemoval']=dS['N OverstoryRemoval']/dS['N Total']
	return dS

#%%
def SensitivityToOverstoryRemovalFrac(meta,pNam,mos,dSor):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	yr_Proj=meta[pNam]['YearCurrent']+1
	iT=np.where( (tv>yr_Proj-10) & (tv<=2050) )[0]
	
	dS={}
	dS['frac_OR']=np.arange(0,0.6,0.1)
	dS['y']=np.zeros(dS['frac_OR'].size)
	for i in range(dS['frac_OR'].size):
		fracPT={'Salvage and Planting Post Beetle':0.5*dS['frac_OR'][i],
			   'Knockdown and Planting':0.5*dS['frac_OR'][i],
			   'Road Rehabilitation':0.0,
			   'Underplanting':1-dS['frac_OR'][i],
			   'Replanting':0.0,
			   'Fill Planting':0.0,
			   'Ecosystem Restoration':0}
		mosPY=unose.ProjectFuture(meta,pNam,mos,'SP Projection',fracPT)
	
		psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Underplanting']
		ysum=np.zeros(tv.size)
		for nPS in psL:
			iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
			iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
			iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
			iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
			y=cbu.GetMosDeltaVar(meta,pNam,mosPY,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
			dS['y'][i]=dS['y'][i]+np.sum(np.nan_to_num(y[iT]))
	
	cl1=meta['Graphics']['Colours']['rgb']['Blue Dark']
	cl2=meta['Graphics']['Colours']['rgb']['Purple Dark']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
	ax.plot([0,1],[0,0],'k-',lw=1,color=[0,0,0])
	ax.plot(dSor['Frac OverstoryRemoval']*np.ones(2),[-5,5],'k--',lw=0.75,color=cl2)
	ax.text(1.04*dSor['Frac OverstoryRemoval'],-1.25,'Overstory removal fraction\nin 2024',color=cl2)
	ax.plot(dS['frac_OR'],dS['y'],'k-',lw=2.25,color=cl1)

	rs,txt=gu.GetRegStats(dS['frac_OR'],dS['y'])
	#ax.plot(rs['xhat Line'],rs['yhat Line'],'k-')
	ax.text(rs['xhat Line'][3]-0.1,rs['yhat Line'][3],txt,color=cl1,ha='left')

	ax.set(xlabel='Fraction overstory removal',ylabel='Emissions impact (MtCO$_2$e)',xlim=[0,0.5],ylim=[-2,1])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	y_min1,y_max1=plt.ylim()
	dy=y_max1-y_min1
	ax.text(0.02,y_min1+0.07*dy,'Net removals',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	ax.text(0.02,y_max1-0.07*dy,'Net emissions',fontsize=12,fontweight='bold',fontstyle='italic',va='center',color=[0.8,0.8,0.8])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Emissions_SensitivityToOverstoryRemovalFrac','png',900)
	return

#%%
def Calc_ServicePlan_ProgressReport(meta,pNam,mos,metaNM,pNamNM,mosNM):

	v='E_NSB'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	dE={}
	dEpHa={}

	# Non-obligation stand establishment
	cNam='NOSE'
	osL=['FIP']
	asetL=np.array(list(meta['LUT']['Derived']['ASET'].keys()))
	iT=np.where( (tv>=meta[pNam]['YearCurrent']-10) & (tv<=2050) )[0]
	for i,aset in enumerate(asetL):
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==aset)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']==str(meta[pNam]['YearCurrent']))[0][0]
		iOS=np.where(np.isin(meta[pNam]['Project']['Strata']['Other']['Unique CD'],osL)==True)[0][0]

		# Total
		y=np.nan_to_num(np.sum(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'][iT]))/meta[pNam]['Project']['Multi']
		dE[aset]=np.round(y,decimals=2)

		# Per-hectare
		y=np.nan_to_num(np.sum(cbu.GetMosDeltaVar(meta,pNam,mos,cNam,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT]))
		dEpHa[aset]=y

	# Forest nutrient management
	cNam='Moderate Harvest'
	pt='Aerial Nutrient Application'

	iT=np.where( (tv>=metaNM[pNamNM]['YearCurrent']) & (tv<=2050) )[0]

	iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
	iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']==str(metaNM[pNamNM]['YearCurrent']))[0][0]
	iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]

	y=np.sum(cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'][iT])/metaNM[pNamNM]['Project']['Multi']
	dE[pt]=np.round(y,decimals=2)
	y=np.nan_to_num(np.sum(cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM,cNam,v,iPS,iSS,iYS,iOS)['Mean']['Ensemble Mean'][iT]))
	dEpHa[pt]=y

	return dE,dEpHa

#%%
def Calc_ServicePlan_Projection(pNam,meta,mosPY,pNamNM,metaNM,mosNM_PY):
	v='E_NSB'
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	dE={}

	# Non-obligation stand establishment
	cNam='NOSE'
	iT=np.where( (tv>=meta[pNam]['YearCurrent']-10) & (tv<=2050) )[0]
	psL=['Salvage and Planting Post Beetle','Knockdown and Planting','Underplanting']
	for nPS in psL:
		iPS=np.where(meta[pNam]['Project']['Strata']['Project Type']['Unique CD']==nPS)[0][0]
		iSS=np.where(meta[pNam]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
		iYS=np.where(meta[pNam]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
		iOS=np.where(meta[pNam]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
		y=np.sum(cbu.GetMosDeltaVar(meta,pNam,mosPY,cNam,'E_NSB',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'][iT]/meta[pNam]['Project']['Multi'])
		dE[nPS]=np.round(y,decimals=2)

	# Forest Nutrient Management
	cNam='Moderate Harvest'
	iT=np.where( (tv>=metaNM[pNamNM]['YearCurrent']) & (tv<=2050) )[0]
	iPS=np.where(metaNM[pNamNM]['Project']['Strata']['Project Type']['Unique CD']=='All')[0][0]
	iSS=np.where(metaNM[pNamNM]['Project']['Strata']['Spatial']['Unique CD']=='All')[0][0]
	iYS=np.where(metaNM[pNamNM]['Project']['Strata']['Year']['Unique CD']=='All')[0][0]
	iOS=np.where(metaNM[pNamNM]['Project']['Strata']['Other']['Unique CD']=='All')[0][0]
	y=np.sum(cbu.GetMosDeltaVar(metaNM,pNamNM,mosNM_PY,cNam,v,iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean'][iT])/meta[pNam]['Project']['Multi']
	dE['Aerial Nutrient Application']=np.round(y,decimals=2)

	return dE

#%%