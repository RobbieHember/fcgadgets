#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import graphviz
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun as cbr
import fcgadgets.macgyver.util_fcs_graphs as ufcs
import fcgadgets.macgyver.util_demo as udem
import fcgadgets.gaia.gaia_util as gaia
warnings.filterwarnings("ignore")
gp=gu.SetGraphics('Manuscript')

#%%
def QA_ConservationOfMassTest_FromDemo(meta,pNam,mos,iScn):

	# This was designed with the Harvest demo in mind

	#iScn=3 # Coastal project
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

	# Ecosystem
	GN0=mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Net']['Ensemble Mean'][:,0,0,0]
	GNR0=mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Net_Reg']['Ensemble Mean'][:,0,0,0]
	DW0=mos[pNam]['Scenarios'][iScn]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0]
	LIT0=mos[pNam]['Scenarios'][iScn]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0]
	B0=mos[pNam]['Scenarios'][iScn]['Mean']['C_Biomass']['Ensemble Mean'][:,0,0,0]
	npp=mos[pNam]['Scenarios'][iScn]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]
	gg=mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Gross']['Ensemble Mean'][:,0,0,0]
	gn=mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Gross']['Ensemble Mean'][:,0,0,0]
	lf=mos[pNam]['Scenarios'][iScn]['Mean']['C_LF']['Ensemble Mean'][:,0,0,0]
	m=mos[pNam]['Scenarios'][iScn]['Mean']['C_M']['Ensemble Mean'][:,0,0,0]
	C0=mos[pNam]['Scenarios'][iScn]['Mean']['C_Forest']['Ensemble Mean'][:,0,0,0]
	npp=mos[pNam]['Scenarios'][iScn]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]
	rh=mos[pNam]['Scenarios'][iScn]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0]
	rem=mos[pNam]['Scenarios'][iScn]['Mean']['C_ToMillTotal']['Ensemble Mean'][:,0,0,0]
	fire=mos[pNam]['Scenarios'][iScn]['Mean']['C_ToFire']['Ensemble Mean'][:,0,0,0]
	
	dC0=npp-rh-rem-fire
	dB0=npp-lf-m
	
	C1=np.cumsum(npp-rh-rem-fire)
	C1a=C1+C0[0]-C1[0]
	dC1a=np.append(0,np.diff(C0))
	
	B1=np.cumsum(npp-lf-m)
	B1a=B1+B0[0]
	dB1a=np.append(0,np.diff(B0))

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,8)); lw=1
	ax[0,0].plot(tv,C0,'b-',lw=lw,label='Carbon stock')
	ax[0,0].plot(tv,C1a,'r--',lw=lw,label='Reconstructed from fluxes')
	ax[0,0].legend(loc='lower left',frameon=False,facecolor=None)
	ax[0,0].set(ylabel='Total ecosystem carbon\n(tC ha$^{-1}$)',xlabel='Time, years')
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[0,1].plot(tv,dC1a,'b-',lw=lw,label='Reconstructed from stock')
	ax[0,1].plot(tv,dC0,'r--',lw=lw,label='Carbon flux')
	ax[0,1].legend(loc='lower left',frameon=False,facecolor=None)
	ax[0,1].set(ylabel='Net carbon balance\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years')
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,0].plot(tv,C0-C1a,'b-',color=[0.5,0,1],lw=lw)
	ax[1,0].set(ylabel='$\Delta$ total ecosystem carbon\n(tC ha$^{-1}$)',xlabel='Time, years')
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,1].plot(tv,dC0-dC1a,'b-',color=[0.5,0,1],lw=lw)
	ax[1,1].set(ylabel='$\Delta$ net carbon balance\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years')
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_ConservationOfMass_Ecosystem','png',900)
	
	# Biomass
	flg=0
	if flg==1:
		plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,8)); lw=1
		ax[0,0].plot(tv,B0,'b-',lw=lw)
		ax[0,0].plot(tv,B1a,'r--',lw=lw)
		ax[0,1].plot(tv,dB0,'b-',lw=lw)
		ax[0,1].plot(tv,dB1a,'r--',lw=lw)
		ax[1,0].plot(tv,np.zeros(tv.size),'k-',lw=lw)
		ax[1,0].plot(tv,B0-B1a,'b-',color=[0.5,0,1],lw=lw)
		ax[1,1].plot(tv,np.zeros(tv.size),'k-',lw=lw)
		ax[1,1].plot(tv,dB0-dB1a,'b-',color=[0.5,0,1],lw=lw)
	
	# Harvested woood products
	C0=mos[pNam]['Scenarios'][iScn]['Mean']['C_HWP']['Ensemble Mean'][:,0,0,0]
	hr=mos[pNam]['Scenarios'][iScn]['Mean']['C_ToMillTotal']['Ensemble Mean'][:,0,0,0]
	bbp=mos[pNam]['Scenarios'][iScn]['Mean']['C_BBP']['Ensemble Mean'][:,0,0,0]
	rhp=mos[pNam]['Scenarios'][iScn]['Mean']['C_RHP']['Ensemble Mean'][:,0,0,0]
	lex=mos[pNam]['Scenarios'][iScn]['Mean']['C_ToLogExport']['Ensemble Mean'][:,0,0,0]
	#plt.plot(tv,lex,'.b-')
	
	dC0=hr-bbp-rhp
	
	C1=np.cumsum(hr-bbp-rhp)
	C1a=C1+C0[0]-C1[0]
	dC1a=np.append(0,np.diff(C0))
	
	# QA Conservation of mass check between pools nd fluxes
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,8)); lw=1
	ax[0,0].plot(tv,C0,'b-',lw=lw,label='Carbon stock')
	ax[0,0].plot(tv,C1a,'r--',lw=lw,label='Reconstructed from fluxes')
	ax[0,0].legend(loc='lower left',frameon=False,facecolor=None)
	ax[0,0].set(ylabel='Total HWP carbon\n(tC ha$^{-1}$)',xlabel='Time, years')
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[0,1].plot(tv,dC1a,'b-',lw=lw,label='Reconstructed from stock')
	ax[0,1].plot(tv,dC0,'r--',lw=lw,label='Carbon flux')
	ax[0,1].legend(loc='lower left',frameon=False,facecolor=None)
	ax[0,1].set(ylabel='Net carbon balance\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years')
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,0].plot(tv,C0-C1a,'b-',color=[0.5,0,1],lw=lw)
	ax[1,0].set(ylabel='$\Delta$ total HWP carbon\n(tC ha$^{-1}$)',xlabel='Time, years')
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].plot(tv,np.zeros(tv.size),'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,1].plot(tv,dC0-dC1a,'b-',color=[0.5,0,1],lw=lw)
	ax[1,1].set(ylabel='$\Delta$ net carbon balance\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years')
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_ConservationOfMass_HWP','png',900)

	return

#%%
def Record_In_ChangeTrackerDB(meta,pNam,mos,tbs,cNamL):
	pth=meta['Paths']['Model']['Code'] + '\\QA\\ChangeTrackingDB.xlsx'
	xfile=openpyxl.load_workbook(pth)
	sheet=xfile.get_sheet_by_name('Sheet1')
	tm=date.today()

	for cNam in cNamL:

		# Find first empty column
		cL=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
		for c in cL:
			col=c + '2'
			if sheet[col].value==None:
				break

		for r in range(1,1000):

			if sheet['A' + str(r)].value=='Date of Entry':
				sheet[c + str(r)].value=tm
			if sheet['A' + str(r)].value=='Year':
				sheet[c + str(r)].value=tm.year
			if sheet['A' + str(r)].value=='Project Code':
				sheet[c + str(r)].value=meta[pNam]['Project']['Code Project']
			if sheet['A' + str(r)].value=='Scenario Comparison Code':
				sheet[c + str(r)].value=cNam
			if sheet['A' + str(r)].value=='Focal Year':
				sheet[c + str(r)].value=meta[pNam]['Project']['Year Project']
			#------------------------------------------------------------------
			# Add variables for each scenario
			#------------------------------------------------------------------
			v='ByScenario'
			scnL=['Baseline','Actual']
			idxL=['iB','iP']
			for iScn,namScn in enumerate(scnL):
				if sheet['B' + str(r)].value==namScn:
					idx=mos[pNam]['Delta'][cNam][ idxL[iScn] ]
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t10':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t50':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2030':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2040':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2050':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][idx]['E_NEB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t10':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t50':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2030':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2040':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2050':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][idx]['E_NSB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t10':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t50':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2030':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2040':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2050':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][idx]['E_NAB']
					if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2100':
						sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][idx]['E_NAB']
	
					if sheet['A' + str(r)].value=='C_Biomass_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_Biomass']
					if sheet['A' + str(r)].value=='C_DeadWood_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_DeadWood']
					if sheet['A' + str(r)].value=='C_Litter_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_Litter']
					if sheet['A' + str(r)].value=='C_Soil_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_Soil']
					if sheet['A' + str(r)].value=='C_Ecosystem_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_Forest']
					if sheet['A' + str(r)].value=='C_HWP_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_HWP']
					if sheet['A' + str(r)].value=='C_HWP_PercentRemainingInHWP_in_t100':
						#y1=tbs[v]['Mean_in_t+1'][idx]['C_HWP']
						y2=tbs[v]['Mean_in_t+100'][idx]['C_HWP']
						y3=tbs[v]['Sum_t+0_to_t+10'][idx]['C_ToMillTotal']
						sheet[c + str(r)].value=np.round(np.nan_to_num(y2/y3*100),decimals=1)
					if sheet['A' + str(r)].value=='C_Geological_in_t100':
						sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][idx]['C_Geological']

					if sheet['A' + str(r)].value=='Cost_Sum_t0_to_t100':
						sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][idx]['Cost Total'],decimals=2)
					if sheet['A' + str(r)].value=='RevenueGross_Sum_t0_to_t100':
						sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][idx]['Revenue Gross'],decimals=2)
					if sheet['A' + str(r)].value=='RevenueNet_Sum_t0_to_t100':
						sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][idx]['Revenue Net'],decimals=2)

					if sheet['A' + str(r)].value=='Cost_PerYield_t-5_to_t10':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][idx]['Cost Total']/tbs[v]['Sum_t-5_to_t+10'][idx]['V_ToMill_MerchTotal']),decimals=2)
					if sheet['A' + str(r)].value=='RevenueGross_PerYield_t-5_to_t10':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][idx]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+10'][idx]['V_ToMill_MerchTotal']),decimals=2)
					if sheet['A' + str(r)].value=='RevenueNet_PerYield_t-5_to_t10':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][idx]['Revenue Net']/tbs[v]['Sum_t-5_to_t+10'][idx]['V_ToMill_MerchTotal']),decimals=2)
					if sheet['A' + str(r)].value=='Cost_PerYield_t-5_to_t100':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][idx]['Cost Total']/tbs[v]['Sum_t-5_to_t+100'][idx]['V_ToMill_MerchTotal']),decimals=2)
					if sheet['A' + str(r)].value=='RevenueGross_PerYield_t-5_to_t100':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][idx]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+100'][idx]['V_ToMill_MerchTotal']),decimals=2)
					if sheet['A' + str(r)].value=='RevenueNet_PerYield_t-5_to_t100':
						sheet[c + str(r)].value=np.round(np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][idx]['Revenue Net']/tbs[v]['Sum_t-5_to_t+100'][idx]['V_ToMill_MerchTotal']),decimals=2)

					L=['MitigationValue_NAB_UpfrontCost_Disc_t-5_to_t100','MitigationValue_NAB_UpfrontCost_Undisc_t-5_to_t100','MitigationValue_NAB_NetRev_Undisc_t-5_to_t100',
						'MitigationValue_NAB_NetRev_Disc_t-5_to_t100','RF_Total_Mean_t0_to_t100','RF_Biochem_Mean_t0_to_t100','RF_Biophys_Mean_t0_to_t100','dTemp_Mean_t0_to_t100']
					if np.isin(sheet['A' + str(r)].value,L)==True:
						sheet[c + str(r)].value='NA'

			#----------------------------------------------------------------------
			# Add delta
			#----------------------------------------------------------------------
			if sheet['B' + str(r)].value=='Delta':
				v='ByScenario'
				iB=mos[pNam]['Delta'][cNam]['iB']
				iP=mos[pNam]['Delta'][cNam]['iP']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t10':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_t+10'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t50':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_t+50'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_t100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_t+100'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2030':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_2030'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2040':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_2040'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2050':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_2050'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NEB_Sum_t0_to_2100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][iP]['E_NEB']-tbs[v]['Sum_t+0_to_2100'][iB]['E_NEB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t10':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_t+10'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t50':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_t+50'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_t100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_t+100'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2030':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_2030'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2040':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_2040'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2050':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_2050'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NSB_Sum_t0_to_2100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][iP]['E_NSB']-tbs[v]['Sum_t+0_to_2100'][iB]['E_NSB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t10':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+10'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_t+10'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t50':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+50'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_t+50'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_t100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_t+100'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_t+100'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2030':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2030'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_2030'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2040':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2040'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_2040'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2050':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2050'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_2050'][iB]['E_NAB']
				if sheet['A' + str(r)].value=='E_NAB_Sum_t0_to_2100':
					sheet[c + str(r)].value=tbs[v]['Sum_t+0_to_2100'][iP]['E_NAB']-tbs[v]['Sum_t+0_to_2100'][iB]['E_NAB']

				if sheet['A' + str(r)].value=='C_Biomass_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_Biomass']-tbs[v]['Mean_in_t+100'][iB]['C_Biomass']
				if sheet['A' + str(r)].value=='C_DeadWood_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_DeadWood']-tbs[v]['Mean_in_t+100'][iB]['C_DeadWood']
				if sheet['A' + str(r)].value=='C_Litter_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_Litter']-tbs[v]['Mean_in_t+100'][iB]['C_Litter']
				if sheet['A' + str(r)].value=='C_Soil_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_Soil']-tbs[v]['Mean_in_t+100'][iB]['C_Soil']
				if sheet['A' + str(r)].value=='C_Ecosystem_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_Forest']-tbs[v]['Mean_in_t+100'][iB]['C_Forest']
				if sheet['A' + str(r)].value=='C_HWP_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_HWP']-tbs[v]['Mean_in_t+100'][iB]['C_HWP']
				if sheet['A' + str(r)].value=='C_HWP_PercentRemainingInHWP_in_t100':
					#y1=tbs[v]['Mean_in_t+1'][iP]['C_HWP']-tbs[v]['Mean_in_t+1'][iB]['C_HWP']
					y2=tbs[v]['Mean_in_t+100'][iP]['C_HWP']-tbs[v]['Mean_in_t+100'][iB]['C_HWP']
					y3=tbs[v]['Sum_t+0_to_t+10'][iP]['C_ToMillTotal']-tbs[v]['Sum_t+0_to_t+10'][iB]['C_ToMillTotal']
					sheet[c + str(r)].value=np.round(y2/y3*100,decimals=1)
				if sheet['A' + str(r)].value=='C_Geological_in_t100':
					sheet[c + str(r)].value=tbs[v]['Mean_in_t+100'][iP]['C_Geological']-tbs[v]['Mean_in_t+100'][iB]['C_Geological']
				if sheet['A' + str(r)].value=='Cost_Sum_t0_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][iP]['Cost Total']-tbs[v]['Sum_t+0_to_t+100'][iB]['Cost Total'],decimals=2)
				if sheet['A' + str(r)].value=='RevenueGross_Sum_t0_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][iP]['Revenue Gross']-tbs[v]['Sum_t+0_to_t+100'][iB]['Revenue Gross'],decimals=2)
				if sheet['A' + str(r)].value=='RevenueNet_Sum_t0_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v]['Sum_t+0_to_t+100'][iP]['Revenue Net']-tbs[v]['Sum_t+0_to_t+100'][iB]['Revenue Net'],decimals=2)

				if sheet['A' + str(r)].value=='Cost_PerYield_t-5_to_t10':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iB]['Cost Total']/tbs[v]['Sum_t-5_to_t+10'][iB]['V_ToMill_MerchTotal'])
					yP=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iP]['Cost Total']/tbs[v]['Sum_t-5_to_t+10'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)
				if sheet['A' + str(r)].value=='RevenueGross_PerYield_t-5_to_t10':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iB]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+10'][iB]['V_ToMill_MerchTotal'])
					yP=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iP]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+10'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)
				if sheet['A' + str(r)].value=='RevenueNet_PerYield_t-5_to_t10':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iB]['Revenue Net']/tbs[v]['Sum_t-5_to_t+10'][iB]['V_ToMill_MerchTotal'])
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+10'][iP]['Revenue Net']/tbs[v]['Sum_t-5_to_t+10'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)
				if sheet['A' + str(r)].value=='Cost_PerYield_t-5_to_t100':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iB]['Cost Total']/tbs[v]['Sum_t-5_to_t+100'][iB]['V_ToMill_MerchTotal'])
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iP]['Cost Total']/tbs[v]['Sum_t-5_to_t+100'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)
				if sheet['A' + str(r)].value=='RevenueGross_PerYield_t-5_to_t100':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iB]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+100'][iB]['V_ToMill_MerchTotal'])
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iP]['Revenue Gross']/tbs[v]['Sum_t-5_to_t+100'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)
				if sheet['A' + str(r)].value=='RevenueNet_PerYield_t-5_to_t100':
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iB]['Revenue Net']/tbs[v]['Sum_t-5_to_t+100'][iB]['V_ToMill_MerchTotal'])
					yB=np.nan_to_num(tbs[v]['Sum_t-5_to_t+100'][iP]['Revenue Net']/tbs[v]['Sum_t-5_to_t+100'][iP]['V_ToMill_MerchTotal'])
					sheet[c + str(r)].value=np.round(yP-yB,decimals=2)

				v='Scenario Comparison'
				if sheet['A' + str(r)].value=='RF_Total_Mean_t0_to_t100':
					sheet[c + str(r)].value=tbs[v][cNam]['Mean_t+0_to_t+100']['RF_Biochem']+\
						tbs[v][cNam]['Mean_t+0_to_t+100']['RF_AlbedoSurfaceShortwave']
				if sheet['A' + str(r)].value=='RF_Biochem_Mean_t0_to_t100':
					sheet[c + str(r)].value=tbs[v][cNam]['Mean_t+0_to_t+100']['RF_Biochem']
				if sheet['A' + str(r)].value=='RF_Biophys_Mean_t0_to_t100':
					sheet[c + str(r)].value=tbs[v][cNam]['Mean_t+0_to_t+100']['RF_AlbedoSurfaceShortwave']
				if sheet['A' + str(r)].value=='dTemp_Mean_t0_to_t100':
					sheet[c + str(r)].value=tbs[v][cNam]['Mean_t+0_to_t+100']['Temperature']

				if sheet['A' + str(r)].value=='MitigationValue_NAB_UpfrontCost_Disc_t-5_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v][cNam]['Sum_t-5_to_t+10']['Cost Total Disc']/(-1*tbs[v][cNam]['Sum_t+0_to_t+100']['E_NAB_Disc']),decimals=1)
				if sheet['A' + str(r)].value=='MitigationValue_NAB_NetRev_Disc_t-5_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v][cNam]['Sum_t-5_to_t+100']['Revenue Net Disc']/(-1*tbs[v][cNam]['Sum_t+0_to_t+100']['E_NAB_Disc']),decimals=1)

				if sheet['A' + str(r)].value=='MitigationValue_NAB_UpfrontCost_Undisc_t-5_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v][cNam]['Sum_t-5_to_t+10']['Cost Total']/(-1*tbs[v][cNam]['Sum_t+0_to_t+100']['E_NAB']),decimals=1)
				if sheet['A' + str(r)].value=='MitigationValue_NAB_NetRev_Undisc_t-5_to_t100':
					sheet[c + str(r)].value=np.round(tbs[v][cNam]['Sum_t-5_to_t+100']['Revenue Net']/(-1*tbs[v][cNam]['Sum_t+0_to_t+100']['E_NAB']),decimals=1)

		xfile.save(pth)
	return

#%% QA Check that E matches E
def QA_CheckEmissions(meta,pNam,mos,scnC):
	# Used in the harvest demo
	scnC='Coast'
	E_co2=mos[pNam]['Delta'][scnC]['ByStrata']['Mean']['E_CO2']['Ensemble Mean'][:,iPS,iSS,iYS]
	E_ch4=mos[pNam]['Delta'][scnC]['ByStrata']['Mean']['E_CH4']['Ensemble Mean'][:,iPS,iSS,iYS]*28
	E_n2o=mos[pNam]['Delta'][scnC]['ByStrata']['Mean']['E_N2O']['Ensemble Mean'][:,iPS,iSS,iYS]*298
	E1=E_co2+E_ch4+E_n2o
	E2=mos[pNam]['Delta']['Coast']['ByStrata']['Mean']['E_NAB']['Ensemble Mean'][:,iPS,iSS,iYS] # tCO2e/ha/year
	
	plt.close('all')
	#plt.plot(E2,E1,'bo')
	plt.plot(E2,'bo')
	plt.plot(E1,'rs')
	return

#%%
def QA_BenchmarkCBM_ClearcutCoast(meta,pNam,mos):
	# Import field plot data
	meta,fp,soc=ufp.ImportFieldPlotData(meta,type='Stand',include_soil='False')
	dBGC=ufp.CalcStatsByBGC(meta,fp,soc)
	dAC=ufp.CalcStatsByAgeClass(meta,fp)

	# *** NOTES ***
	# To force tree biomass to match CBM25, turn phase shift on and add slight decline with age
	# *** NOTES ***
	StatusTuned=''
	#StatusTuned='Tuned'

	# Import CBM data
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250109\CBM Outputs (Scenario_ID_5 & 6) revised.xls',sheet_name='CBM Outputs (tC)',skiprows=0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250120\Scenario ID 5 & 6 CBM Outputs rah.xlsx','tC',0)
	dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\Comparison with CBM-CFS3\20250305\Scenario ID 5 & 6 CBM Outputs rah.xlsx',sheet_name='tC',skiprows=0)
	# Synchronize over time
	dCBM={}
	for k in dCBM0.keys():
		dCBM[k]=np.nan*np.ones(tv.size)
		ind1=np.where( (dCBM0['Years']>=tv[0]) & (dCBM0['Years']<=tv[-1]) )[0]
		ind2=np.where( (tv>=dCBM0['Years'][0]) & (tv<=dCBM0['Years'][-1]) )[0]
		dCBM[k][ind2]=dCBM0[k][ind1]
	dCBM['Year']=tv
	#dCBM['Fire S_ID-5']=dCBM['NetCO2emissions_removals_CO2e S_ID-5']+dCBM['SumofCOProduction_CO2e S_ID-5']+dCBM['SumofCH4Production_CO2e S_ID-5']+dCBM['N2O_CO2e S_ID-5']
	#dCBM['Fire S_ID-6']=dCBM['NetCO2emissions_removals_CO2e S_ID-6']+dCBM['SumofCOProduction_CO2e S_ID-6']+dCBM['SumofCH4Production_CO2e S_ID-6']+dCBM['N2O_CO2e S_ID-6']
	dCBM['Fire S_ID-5']=dCBM['NetCO2emissions_removals_CO2e S_ID-5']-3.667*dCBM['Harvest Removals S_ID 5']+3.667*dCBM['Net Ecosystem Productivity S_ID-5']+dCBM['SumofCOProduction_CO2e S_ID-5']+dCBM['SumofCH4Production_CO2e S_ID-5']+dCBM['N2O_CO2e S_ID-5']
	dCBM['Fire S_ID-6']=dCBM['NetCO2emissions_removals_CO2e S_ID-6']-3.667*dCBM['Harvest Removals S_ID 6']+3.667*dCBM['Net Ecosystem Productivity S_ID-6']+dCBM['SumofCOProduction_CO2e S_ID-6']+dCBM['SumofCH4Production_CO2e S_ID-6']+dCBM['N2O_CO2e S_ID-6']
	dCBM['HR S_ID-5']=dCBM['Harvest Removals S_ID 5']
	dCBM['HR S_ID-6']=dCBM['Harvest Removals S_ID 6']
	#dCBM['NEB S_ID-5']=dCBM['Fire S_ID-5']+dCBM['Harvest Removals S_ID 5']-3.667*dCBM['Net Ecosystem Productivity S_ID-5']
	#dCBM['NEB S_ID-6']=dCBM['Fire S_ID-6']+dCBM['Harvest Removals S_ID 6']-3.667*dCBM['Net Ecosystem Productivity S_ID-6']
	dCBM['NEB S_ID-5']=dCBM['NetCO2emissions_removals_CO2e S_ID-5']
	dCBM['NEB S_ID-6']=dCBM['NetCO2emissions_removals_CO2e S_ID-6']

	# Prep fcgadgets results
	iP=3
	iB=2
	A=mos[pNam]['Scenarios'][iP]['Mean']['A']['Ensemble Mean'][:,0,0,0]

	# Index to field plot sample
	#iFP=np.where( (fp['LS']=='FDC') & (fp['Ecozone BC L2']==meta['LUT']['GP']['Ecozone BC L2']['CWHdm']) )[0]
	iFP=np.where( (fp['LS']=='FDC') & (fp['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) )[0]

	#--------------------------------------------------------------------------
	# Carbon pools
	#--------------------------------------------------------------------------

	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15,11)); xticks=np.arange(1800,2220,30)
	ax[0,0].plot(1910,327,'ks',lw=1,ms=3,mec='k',mfc='w',label='Trofymow et al\n2008')
	ax[0,0].plot(2023+50,168,'ks',lw=1,ms=4,mec='k',mfc='w')
	ax[0,0].plot(2023+(2002-1988),24,'r^',ms=3,mfc='w',label='Ferster et al\n2015') #
	ax[0,0].plot(2023+(2002-1949),177,'r^',ms=3,mfc='w')
	bw=10; bin=np.arange(10,100+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot L t0'][iFP],bw,bin)
	ax[0,0].plot(bin+2023,mu,'go',ms=3,mfc='w',label='Field plots')
	ax[0,0].plot(tv,dCBM['Biomass S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#ax[0,0].plot(tv-108,dCBM['Biomass S_ID-6'],'k--',color=meta['Graphics']['gp']['cl1'],label='CBM-CFS3') # Compare planted and natural CBM curve
	ax[0,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Biomass']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=xticks,ylabel='Tree biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2180])
	ax[0,0].legend(loc='lower right',fontsize=4,facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[0,1].plot(1910,76,'ks',lw=1,ms=3,mec='k',mfc='w',label='Trofymow et al (2008)')
	#ax[0,1].plot(tv,dBGC['CN']['CWH']['Ctot D t0']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='CMI Network (Standing)')
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0]#+mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[0,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0],'k-',color=meta['Graphics']['gp']['cl1'],label='Forest Carbon Gadgets (standing)')
	ax[0,1].plot(tv,dCBM['Deadwood S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[0,1].set(xticks=xticks,ylabel='Dead wood (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,400])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv,dBGC['Soil']['CWH']['SOC_ORG_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB\n(Shaw et al. 2018)')
	ax[1,0].plot(1910,111,'ks',lw=1,ms=3,mec='k',mfc='w',label='Trofymow et al (2008)')
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0]#-mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Forest Carbon Gadgets')
	ax[1,0].plot(tv,dCBM['Litter S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#ax[1,0].plot(tv,dCBM['Aboveground DOM S_ID-6'],'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM-CFS3')
	#a=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,dCBM['Litter S_ID-6']-,'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM-CFS3')
	ax[1,0].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_LitterVF']['Ensemble Mean'][:,0,0,0],'m:',label='FCG LitterVF')
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_LitterF']['Ensemble Mean'][:,0,0,0],'r-.',label='FCG LitterF')
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_LitterS']['Ensemble Mean'][:,0,0,0],'c:',label='FCG LitterS')
	ax[1,0].set(xticks=xticks,ylabel='Litter (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,240])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)

	ax[1,1].plot(tv,dBGC['Soil']['CWH']['SOC_MIN_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB (Shaw et al. 2018)')
	ax[1,1].plot(1910,227,'ks',lw=1,ms=3,mec='k',mfc='w',label='Trofymow et al (2008)')
	ax[1,1].plot(tv,dCBM['Soil S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Soil']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=xticks,ylabel='Soil (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,350])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].legend(loc='lower right',fontsize=5,facecolor=[1,1,1],frameon=False)

	ax[2,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Piles']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	#ax[2,0].plot(tv,dCBM['Total Litter Project'],'k--',color=meta['Graphics']['gp']['cl2'],label='CBM-CFS3')
	ax[2,0].set(xticks=xticks,ylabel='Piles (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150])
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[2,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)

	ax[2,1].plot(tv,dCBM['Total Ecosystem S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[2,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Forest']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	#ax[2,1].plot(tv,dCBM['Biomass S_ID-6']+dCBM['Deadwood S_ID-6']+dCBM['Litter S_ID-6']+dCBM['Soil S_ID-6'],'k--',label='CBM-CFS3') # Check consistency

	ax[2,1].set(xticks=xticks,ylabel='Total ecosystem (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150])
	ax[2,1].legend(loc='lower left',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Pools_' + StatusTuned,'png',900)

	iT1=np.where( (tv==2022) )[0]
	iT2=np.where( (tv==2024) )[0]
	y=dCBM['Biomass S_ID-6']; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Biomass']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)

	y=dCBM['Deadwood S_ID-6']; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)

	y=dCBM['Litter S_ID-6']; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)

	y=dCBM['Soil S_ID-6']; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Soil']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)

	y=dCBM['Total Ecosystem S_ID-6']; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Forest']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)

	y=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Soil_OHorizon']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Soil']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_Soil_OHorizon']['Ensemble Mean'][:,0,0,0]; dr=(y[iT2]-y[iT1])/y[iT1]*100
	print(dr)


	#--------------------------------------------------------------------------
	# Net Growth
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
	ax.plot(54,10.2/4,'r^',ms=3,mfc='w',mec='r',label='Ferster et al. (2015)')
	ax.plot((2004-1988),10/4,'rs',ms=3,mfc='w',mec='r')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot Net'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=2*se,color='k',ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',label='Field plots')
	ax.plot(tv[1:]-2023,np.diff(dCBM['Biomass S_ID-6']),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_G_Net']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax.set(xticks=np.arange(0,2200,20),yticks=np.arange(-10,10,0.5),ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Age, years',xlim=[0,120],ylim=[0,7])
	ax.legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NetGrowth_' + StatusTuned,'png',900)

	# Fit CBM25
	flg=0
	if flg==1:
		def funB(A,**b0):
			# Chapman Richards (Zeide 1993)
			yhat=b0['b1']*(1-np.exp(-b0['b2']*A))**b0['b3']
			return yhat
		iT=np.where( (tv>=1915) & (tv<=2022) )[0]
		d={'A':A[iT],'B':dCBM['Biomass S_ID-6'][iT]}
		b0={'b1':400,'b2':0.05,'b3':2}
		plt.close('all');plt.plot(d['A'],d['B'],'bo')
		plt.plot(d['A'],b0['b1']*(1-np.exp(-b0['b2']*d['A']))**b0['b3'],'k-')

		fitmodel=lmfit.Model(funB,independent_vars=['A'])
		params=fitmodel.make_params()
		for k,v in b0.items():
			params.add(k,value=v)
		rs=fitmodel.fit(d['B'],params,A=d['A'])
		b=[rs.params['b1'].value,rs.params['b2'].value,rs.params['b3'].value]
		plt.close('all');plt.plot(d['A'],d['B'],'bo')
		plt.plot(d['A'],b[0]*(1-np.exp(-b[1]*d['A']))**b[2],'k-')

		def funG(A,**b0):
			B1=b0['b1']
			B2=b0['b2']
			B3=b0['b3']
			B4=b0['b4']
			B5=b0['b5']
			yhat=(B1*(1+((B2*(A/B3)**B4-1)/np.exp(A/B3))))+(B5*A)
			yhat=np.maximum(0,np.nan_to_num(yhat))
			return yhat
		iT=np.where( (tv>=1915) & (tv<=2022) )[0]
		d={'A':A[iT][1:],'G':np.diff(dCBM['Biomass S_ID-6'][iT])}
		b0={}
		b0['b1']=2
		b0['b2']=0.5
		b0['b3']=9
		b0['b4']=3.5
		b0['b5']=0.001
		plt.close('all');plt.plot(d['A'],d['G'],'bo')
		plt.plot(d['A'],(b0['b1']*(1+((b0['b2']*(d['A']/b0['b3'])**b0['b4']-1)/np.exp(d['A']/b0['b3']))))+(b0['b5']*d['A']),'k-') #

		fitmodel=lmfit.Model(funG,independent_vars=['A'])
		params=fitmodel.make_params()
		for k,v in b0.items():
			params.add(k,value=v)
		rs=fitmodel.fit(d['G'],params,A=d['A'])
		b0=[rs.params['b1'].value,rs.params['b2'].value,rs.params['b3'].value,rs.params['b4'].value,rs.params['b5'].value]
		plt.close('all');plt.plot(d['A'],d['G'],'bo')
		plt.plot(d['A'],(b0[0]*(1+((b0[1]*(d['A']/b0[2])**b0[3]-1)/np.exp(d['A']/b0[2]))))+(b0[4]*d['A']),'k-') # +(b0['b5']*d['A'])

	#--------------------------------------------------------------------------
	# Carbon fluxes
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,14))
	ax[0,0].plot(tv,dCBM['Net Ecosystem Productivity S_ID-6'],'.k-',color=meta['Graphics']['gp']['cl1'],label='CBM-CFS3')
	ax[0,0].plot(tv-1,mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0],'.k--',color=meta['Graphics']['gp']['cl2'],label='FCS')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='Net ecoystem producitvity (tC/ha/yr)',xlabel='Time, years',xlim=[1880,2150])
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Wildfire
	iT=np.where(tv==1915)[0][0]
	ax[0,1].bar(1,dCBM['Fire S_ID-6'][iT])
	ax[0,1].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_Wildfire']['Ensemble Mean'][iT,0,0,0])
	# Slashpile burn
	iT=np.where(tv==2023)[0][0]
	ax[1,0].bar(1,dCBM['Fire S_ID-6'][iT])
	ax[1,0].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_OpenBurning']['Ensemble Mean'][iT,0,0,0])
	# Harvest removals
	ax[1,1].bar(1,np.nansum(dCBM['HR S_ID-6']))
	ax[1,1].bar(2,np.sum(mos[pNam]['Scenarios'][iP]['Mean']['C_ToMillTotal']['Ensemble Mean'][:,0,0,0]))
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Fluxes_' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Delta
	#--------------------------------------------------------------------------
	iT=np.where(tv>=2000)[0]
	plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,9))
	d1=dCBM['Net Ecosystem Productivity S_ID-6']-dCBM['Net Ecosystem Productivity S_ID-5']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[0,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv[iT]-1,d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ Net ecoystem producitvity\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,0].legend(loc='lower right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['Net Ecosystem Productivity S_ID-6']-dCBM['Net Ecosystem Productivity S_ID-5']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,1].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[0,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv[iT]-1,np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].set(xticks=np.arange(1800,2200,20),ylabel='Cumulative $\Delta$ net ecoystem\nproducitvity(tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,1].legend(loc='lower right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-6']-dCBM['NEB S_ID-5']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[1,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(tv[iT],d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ emissions\n(tCO2e ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[1,0].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-6']-dCBM['NEB S_ID-5']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv[iT],np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=np.arange(1800,2200,20),ylabel='Cumulative $\Delta$ emissions\n(tCO2e ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[1,1].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	iT=np.where( (tv>=2023) & (tv<=2023+30) )[0]
	iT=np.where( (tv>=2023) & (tv<=2023+100) )[0]
	s1=np.sum(d1[iT])
	s2=np.sum(d2[iT])
	print(int(s1))
	print(int(s2))
	print(int(s2-s1))
	print(int((s2-s1)/s1*100))


	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Delta_' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Stemwood
	#--------------------------------------------------------------------------

	dGY=gu.ReadExcel(r'D:\Modelling Projects\Demo_Harv_Clearcut\Inputs\TIPSY output.xlsx')
	dGY['Merch Ratio']=np.nan_to_num(dGY['Net volume 175']/dGY['Volume total'])
	dGY['Stemwood Merch Biomass']=dGY['Merch Ratio']*dGY['Stemwood biomass']
	dGY['Stemwood Non-merch Biomass']=(1-dGY['Merch Ratio'])*dGY['Stemwood biomass']

	iGY7=np.where(dGY['Scenario']==7)[0]
	iGY8=np.where(dGY['Scenario']==8)[0]
	iT1=np.where( (tv>=1915) & (tv<=2020) )[0]
	iT2=np.where( (tv>=2023) )[0]

	iT1a=np.where( (tv>=1915) & (tv<=2020) & (A>=1) & (A<=100) )[0]
	iGY7a=np.where( (dGY['Scenario']==7) & (dGY['Age']>=1) & (dGY['Age']<=100) )[0]
	iT2a=np.where( (tv>=2023) & (A>=1) & (A<=100) )[0]
	iGY8a=np.where( (dGY['Scenario']==8) & (dGY['Age']>=1) & (dGY['Age']<=100) )[0]

	# Boudwyn et al Douglas-fir - Montane Cordillera
	#b_m=1.07518*dGY['Net volume 175']**0.84843 # Merch (including stumps and tops)
	b_m=0.55*dGY['Net volume 175']**1.002 # Merch (including stumps and tops)
	#nonmerchfactor=0.86758+55.22243*b_m**-1.05248 # total/b_m
	nonmerchfactor=0.9127+32.32*b_m**-0.9449 # total/b_m
	b_n=nonmerchfactor*b_m-b_m

	b_n_FP=fp['Csw L t0']-fp['Csw175 L t0']

	#plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,9.5))
	ax[0,0].plot(dGY['Age'][iGY7],dGY['Net volume 175'][iGY7],'b-',color=meta['Graphics']['gp']['cl2'],label='TIPSY')
	ax[0,0].set(ylabel='Net volume 17.5cm\nutilization (m$^3$ ha$^{-1}$)',xlabel='Age, years',xlim=[0,200])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

	ax[0,1].plot(dGY['Age'][iGY7],0.5*dGY['Stemwood biomass'][iGY7],'b--',color=meta['Graphics']['gp']['cl2'],label='Total stemwood (TIPSY)')
	ax[0,1].plot(dGY['Age'][iGY7],0.5*dGY['Stemwood Merch Biomass'][iGY7],'b-',color=meta['Graphics']['gp']['cl2'],label='Merch. stemwood (TIPSY)')
	#ax[0,1].plot(A[iT1],dCBM['Biomass S_ID-5'][iT1],'c-',label='Total biomass (CBM-CFS3)')
	ax[0,1].plot(A[iT1],dCBM['Softwood Merchantable S_ID-5'][iT1],'c-',color=meta['Graphics']['gp']['cl1'],label='Merch stemwood (CBM25)')
	ax[0,1].plot(dGY['Age'][iGY7],0.5*b_m[iGY7],'r--',label='Merch. stemwood\n(Boudwyn et al. 2007)')
	ax[0,1].plot(dGY['Age'][iGY7],0.5*b_n[iGY7],'r:',label='Non-merch. stemwood\n(Boudwyn et al. 2007)')
	ax[0,1].set(ylabel='Biomass (tC ha$^{-1}$)',xlabel='Age, years',xlim=[0,100],ylim=[0,250]); ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].legend(loc='upper left',fontsize=5,facecolor=[1,1,1],frameon=False)

	ax[1,0].plot([0,800],[0.53,0.53],'-',color=[0.8,0.8,0.8],lw=2,label='Constant wood density\n= 0.5 (ODT m$^3$)')
	ax[1,0].plot(dGY['Net volume 175'][iGY7a],dGY['Stemwood Merch Biomass'][iGY7a]/dGY['Net volume 175'][iGY7a],'-b',color=meta['Graphics']['gp']['cl2'],label='TIPSY')
	ax[1,0].plot(dGY['Net volume 175'][iGY7a],2*dCBM['Softwood Merchantable S_ID-5'][iT1a]/dGY['Net volume 175'][iGY7a],'-c',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(dGY['Net volume 175'][iGY7a],b_m[iGY7a]/dGY['Net volume 175'][iGY7a],'r--',label='Boudwyn et al. (2007)')
	ax[1,0].set(ylabel='Merch. wood density\n(ODT m$^3$)',xlabel='Net volume 17.5 (m$^3$ ha$^{-1}$)',ylim=[0,4],xlim=[0,200])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)

	#ax[1,1].plot(dGY['Age'][iGY7],dGY['Stemwood Non-merch Biomass'][iGY7]/dGY['Stemwood biomass'][iGY7],'b-')
	ax[1,1].plot([0,22000],[0,0],'k-')
	ax[1,1].plot(dGY['Age'][iGY7],0.5*dGY['Stemwood Non-merch Biomass'][iGY7],'b-',color=meta['Graphics']['gp']['cl2'],label='TIPSY')
	ax[1,1].plot(dGY['Age'][iGY7],0.5*b_n[iGY7],'r--',label='Boudwyn et al. (2007)')
	bw=10; bin=np.arange(10,220+bw,bw);N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],b_n_FP[iFP],bw,bin)
	ax[1,1].plot(bin,mu,'ko',ms=3,lw=0.75,mec='k',mfc='w',label='Field plots')
	ax[1,1].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[1,1].set(ylabel='Non-merch stemwood\nbiomass (tC ha$^{-1}$)',xlabel='Age, years',xlim=[0,200])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Stemwood_' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Non-stemwood biomass pools
	#--------------------------------------------------------------------------

	# With CBM25
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10)); xticks=np.arange(-2000,2220,20); xlim=[-120,125]
	ax[0,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax[0,0].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[0,0].plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	ax[0,0].plot(tv-2023,dCBM['Foilage S_ID 6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	x=2002-np.array([2000,1988,1949]); y=np.array([0,2.6,12.7])
	ax[0,0].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	ax[0,0].set(xticks=xticks,ylabel='Foliage biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,25])
	ax[0,0].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Branch (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP],bw,bin)
	ax[0,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[0,1].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Branches (field plots)')
	ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl3'],label='Bark (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	ax[0,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl3'],ls='',lw=0.5,capsize=2)
	ax[0,1].plot(bin,mu,'k^',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl3'],label='Bark (field plots)')
	ax[0,1].set(xticks=xticks,ylabel='Bark and branches (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,80])
	ax[0,1].legend(loc='upper center',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0]+mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP]+fp['Cbk L t0'][iFP],bw,bin)
	ax[1,0].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[1,0].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	ax[1,0].plot(tv-2023,dCBM['Other S_ID_6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	x=2002-np.array([2000,1988,1949]); y=np.array([0,2.6,22.1+17.6])
	ax[1,0].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	ax[1,0].set(xticks=xticks,ylabel='Bark+branches (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,80])
	ax[1,0].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Root']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cr L t0'][iFP],bw,bin)
	ax[1,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[1,1].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	ax[1,1].plot(tv-2023,dCBM['Roots S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	x=2002-np.array([2000,1988,1949]); y=np.array([0.1,2.8,32.6])
	ax[1,1].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	ax[1,1].set(xticks=xticks,ylabel='Roots, total (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,80])
	ax[1,1].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NonStemwoodBiomass_' + StatusTuned,'png',900)

	plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10)); xticks=np.arange(0,2220,10)
	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Foliage (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Foliage (field plots)')
	ax.plot(tv-2023,dCBM['Foilage S_ID 6'],'k-',color=meta['Graphics']['gp']['cl1'],label='Foliage (CBM25)')
	ax.set(xticks=xticks,ylabel='Foliage biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=[0,125],ylim=[0,20])
	ax.legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

	age=np.array([1,2,3,4,20,55])
	lai=np.array([0.24,0.76,2.17,2.53,10,7.3]) # m2 m-2

	sla=65 # cm2 g-1
	sla=65/10000 # m2 g-1
	sla=sla*1e6 # m2 MgDM-1
	sla=sla/0.5 # m2 MgC
	cf=lai/sla # MgC m-2
	cf=cf*10000 # MgC ha-1
	ax.plot(age,cf,'ko',mfc='w',mec='k')

	# Just FCG25
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,6)); xticks=np.arange(0,2220,10)
	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Foliage (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Foliage (field plots)')

	ax.plot(tv-2023,dCBM['Foilage S_ID 6'],'k-',color=meta['Graphics']['gp']['cl1'],label='Foliage (CBM25)')

	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl1'],label='Branch (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl1'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl1'],label='Branches (field plots)')

	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl3'],label='Bark (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl3'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'k^',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl3'],label='Bark (field plots)')

	ax.set(xticks=xticks,ylabel='Biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=[0,125],ylim=[0,50])
	ax.legend(loc='upper left',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NonStemwoodBiomass_' + StatusTuned,'png',900)

	flg=0
	if flg==1:
		sla=6 # m2/kg Landsberg and Waring (1997)
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,6)); xticks=np.arange(0,2220,10)
		ax.plot(tv-1915,mos[pNam]['Scenarios'][iScn]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0]*1000/10000*sla,'k-',color=meta['Graphics']['gp']['cl1'],label='Natural (FCG25)')
		ax.plot(tv-2023,mos[pNam]['Scenarios'][iScn]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0]*1000/10000*sla,'k--',color=meta['Graphics']['gp']['cl2'],label='Planted (FCG25)')
		ax.set(xticks=xticks,ylabel='Leaf area index (m$^2$ m$^{-2}$)',xlabel='Time, years',xlim=[0,100],ylim=[0,14])
		ax.legend(loc='lower right',fontsize=4,facecolor=[1,1,1],frameon=False)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

	#--------------------------------------------------------------------------
	# Check conservation of mass
	#--------------------------------------------------------------------------
	plt.close('all')
	plt.plot(dCBM['NEB S_ID-6']/3.667,'b-')
	plt.plot(-np.diff(dCBM['Total Ecosystem S_ID-6']),'r--')

	y1=dCBM['Total Ecosystem S_ID-6']
	y2=-1*np.cumsum(dCBM['NEB S_ID-6'])/3.667+dCBM['Total Ecosystem S_ID-6'][0]
	d=y1-y2
	plt.close('all')
	plt.plot(d-d[0],'b-')

	return

#%%
def QA_BenchmarkCBM_RefUnderplanting(meta,pNam,mos):
	zone='IDF'

	# Import field plot data
	meta,fp,soc=ufp.ImportFieldPlotData(meta,type='Stand',include_soil='False')
	dBGC=ufp.CalcStatsByBGC(meta,fp,soc)
	dAC=ufp.CalcStatsByAgeClass(meta,fp)

	# Import CBM data
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250109\CBM Outputs (Scenario_ID_5 & 6) revised.xls',sheet_name='CBM Outputs (tC)',skiprows=0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250120\Demo_RefUnderFire_Scenario ID 2 & 3 CBM Outputs.xlsx',sheet_name='tC',skiprows=0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250127\Scenario ID 2 & 3 CBM Outputs -v2.xlsx',sheet_name='tC',skiprows=0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\Comparison with CBM-CFS3\20250127\Scenario ID 2 & 3 CBM Outputs -v2.xlsx',sheet_name='tC',skiprows=0)
	dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\Comparison with CBM-CFS3\20250330\Scenario ID 2 & 3 CBM Outputs -v2.xlsx',sheet_name='tC',skiprows=0)

	dCBM={}
	for k in dCBM0.keys():
		dCBM[k]=np.nan*np.ones(tv.size)
		ind1=np.where( (dCBM0['Years']>=tv[0]) & (dCBM0['Years']<=tv[-1]) )[0]
		ind2=np.where( (tv>=dCBM0['Years'][0]) & (tv<=dCBM0['Years'][-1]) )[0]
		dCBM[k][ind2]=dCBM0[k][ind1]
	dCBM['Year']=tv
	dCBM['Fire S_ID-2']=dCBM['NetCO2emissions_removals_CO2e S_ID-2']+3.667*dCBM['Net Ecosystem Productivity S_ID-2']
	dCBM['Fire S_ID-3']=dCBM['NetCO2emissions_removals_CO2e S_ID-3']+3.667*dCBM['Net Ecosystem Productivity S_ID-3']
	dCBM['NEB S_ID-2']=dCBM['NetCO2emissions_removals_CO2e S_ID-2']
	dCBM['NEB S_ID-3']=dCBM['NetCO2emissions_removals_CO2e S_ID-3']

	# Prep fcgadgets results
	iB=1
	iP=3
	A=mos[pNam]['Scenarios'][iP]['Mean']['A']['Ensemble Mean'][:,0,0,0]

	iFP=np.where( (fp['Ctot Net']>-20) & (fp['Ecozone BC L2']==meta['LUT']['GP']['Ecozone BC L2']['IDFdk']) )[0]

	StatusTuned=''
	#StatusTuned='_Tuned'

	#--------------------------------------------------------------------------
	# Carbon pools
	#--------------------------------------------------------------------------
	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,14))
	ax[0,0].plot(tv,dCBM['Biomass S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Biomass']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')

	# Add field plots
	bw=10; bin=np.arange(10,110+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot L t0'][iFP],bw,bin)
	ax[0,0].errorbar(bin+1897,mu,yerr=2*se,color='k',ls='',lw=0.25,capsize=2)
	ax[0,0].plot(bin+1897,mu,'ko',ms=3,mfc='w',lw=0.25,label='Field plots')

	ax[0,0].errorbar(bin+2017,mu,yerr=2*se,color='k',ls='',lw=0.25,capsize=2)
	ax[0,0].plot(bin+2017,mu,'ko',ms=3,mfc='w',lw=0.25)

	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='Tree biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150])
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	y=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0]#+mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	ax[0,1].plot(tv,dCBM['Deadwood S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[0,1].set(xticks=np.arange(1800,2200,20),ylabel='Dead wood (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,200])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv,dBGC['Soil'][zone]['SOC_ORG_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB\n(Shaw et al. 2018)')
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0]#-mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].plot(tv,dCBM['Litter S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#ax[1,0].plot(tv,dCBM['Aboveground DOM S_ID-3'],'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#a=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,dCBM['Litter S_ID-3']-,'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].set(xticks=np.arange(1800,2200,20),ylabel='Litter (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,200])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)

	ax[1,1].plot(tv,dBGC['Soil'][zone]['SOC_MIN_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB (Shaw et al. 2018)')
	ax[1,1].plot(tv,dCBM['Soil S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Soil']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=np.arange(1800,2200,20),ylabel='Soil (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,200])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

	ax[2,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Piles']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[2,0].set(xticks=np.arange(1800,2200,20),ylabel='Piles (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150],ylim=[0,200])
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[2,1].plot(tv,dCBM['Total Ecosystem S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[2,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Forest']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[2,1].set(xticks=np.arange(1800,2200,20),ylabel='Total ecosystem (tC ha$^{-1}$)',xlabel='Time, years',xlim=[1880,2150])
	ax[2,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Pools' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Net Growth
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
	ax.plot(tv[1:]-1897,np.diff(dCBM['Biomass S_ID-3']),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax.plot(tv-1897,mos[pNam]['Scenarios'][iP]['Mean']['C_G_Net']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	#ax.plot(tv[1:]-1897,np.diff(mos[pNam]['Scenarios'][iP]['Mean']['C_StemNonMerch']['Ensemble Mean'][:,0,0,0]),'k:',color=meta['Graphics']['gp']['cl2'],label='FCG non-merch stemwood')
	#ax.plot(tv[1:]-1897,np.diff(mos[pNam]['Scenarios'][iP]['Mean']['C_StemMerch']['Ensemble Mean'][:,0,0,0]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG merch stemwood')
	# Add field plots
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot Net'][iFP],bw,bin)
	ikp=np.where(N>=20)[0]
	ax.errorbar(bin[ikp],mu[ikp],yerr=2*se[ikp],color='k',ls='',lw=0.75,capsize=2)
	ax.plot(bin[ikp],mu[ikp],'ko',ms=3,mec='k',mfc='w',mew=1,label='Field plots')
	ax.set(xticks=np.arange(0,2200,20),yticks=np.arange(-10,10,0.5),ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Age, years',xlim=[0,125],ylim=[0,3])
	ax.legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NetGrowth' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Carbon fluxes
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,14))
	ax[0,0].plot(tv,dCBM['Net Ecosystem Productivity S_ID-3'],'.k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv-1,mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0],'.k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='Net ecoystem producitvity (tC/ha/yr)',xlabel='Time, years',xlim=[1880,2150])
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Wildfire
	iT=np.where(tv==1915)[0][0]
	ax[0,1].bar(1,dCBM['Fire S_ID-3'][iT])
	ax[0,1].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_Wildfire']['Ensemble Mean'][iT,0,0,0])
	# Slashpile burn
	iT=np.where(tv==2023)[0][0]
	ax[1,0].bar(1,dCBM['Fire S_ID-3'][iT])
	ax[1,0].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_OpenBurning']['Ensemble Mean'][iT,0,0,0])
	# Harvest removals
	#ax[1,1].bar(1,np.nansum(dCBM['HR S_ID-3']))
	#ax[1,1].bar(2,np.sum(mos[pNam]['Scenarios'][iP]['Mean']['C_ToMillTotal']['Ensemble Mean'][:,0,0,0]))
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Fluxes','png',900)

	#--------------------------------------------------------------------------
	# Delta
	#--------------------------------------------------------------------------
	iT=np.where(tv>=2022)[0]
	plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,8.5))
	d1=dCBM['Net Ecosystem Productivity S_ID-3']-dCBM['Net Ecosystem Productivity S_ID-2']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[0,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv[iT]-1,d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ net ecoystem producitvity\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['Net Ecosystem Productivity S_ID-3']-dCBM['Net Ecosystem Productivity S_ID-2']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv[iT]-1,np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].set(xticks=np.arange(1800,2200,20),ylabel='Cumulative $\Delta$ net ecoystem\nproducitvity (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-3']-dCBM['NEB S_ID-2']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[1,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(tv[iT],d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ emissions (tCO$_2$e ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[1,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-3']-dCBM['NEB S_ID-2']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,1].plot([2050,2050],[0,-1000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv[iT],np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=np.arange(1800,2200,10),ylabel='Cumulative $\Delta$ emissions\n(tCO$_2$e ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2055],ylim=[-250,0])
	ax[1,1].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	iT=np.where( (tv>=2022) & (tv<=2022+30) )[0]
	#iT=np.where( (tv>=2022) & (tv<=2022+100) )[0]
	s1=np.sum(d1[iT])
	s2=np.sum(d2[iT])
	print(int(s1))
	print(int(s2))
	print(int(s2-s1))
	print(int((s2-s1)/s1*100))

	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Delta' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Non-stemwood biomass pools
	#--------------------------------------------------------------------------

	# Import age functions from EP964
	xhat=np.arange(0,26,1)
	pd=1500/1000
	rs=gu.ipickle(r'C:\Data\EP964\Processed\LinearAgeFunction.pkl')

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10)); xticks=np.arange(-2000,2220,20); xlim=[-20,125]
	ax[0,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax[0,0].errorbar(bin,mu,yerr=2*se,color='k',ls='',lw=0.5,capsize=2)
	ax[0,0].plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec='k',label='Field plots')
	#ax[0,0].plot(tv-2023,dCBM['Foilage S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#x=2002-np.array([2000,1988,1949]); y=np.array([0,2.6,12.7])
	#ax[0,0].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	yhat=rs['f'].Intercept+rs['f'].A*xhat+rs['f'].PD*pd
	ax[0,0].plot(xhat,yhat,'m-',label='EP964')
	ax[0,0].set(xticks=xticks,ylabel='Foliage (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,10])
	ax[0,0].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Bark (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	ax[0,1].errorbar(bin,mu,yerr=2*se,color='k',ls='',lw=0.5,capsize=2)
	ax[0,1].plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec='k',label='Branches (field plots)')
	yhat=rs['bk'].Intercept+rs['bk'].A*xhat+rs['bk'].PD*pd
	ax[0,1].plot(xhat,yhat,'m-',label='EP964')

	#ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl3'],label='Bark (FCG25)')
	#bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	#ax[0,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl3'],ls='',lw=0.5,capsize=2)
	#ax[0,1].plot(bin,mu,'k^',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl3'],label='Bark (field plots)')
	ax[0,1].set(xticks=xticks,ylabel='Bark (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,10])
	#ax[0,1].legend(loc='upper center',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP]+fp['Cbk L t0'][iFP],bw,bin)
	ax[1,0].errorbar(bin,mu,yerr=2*se,color='k',ls='',lw=0.5,capsize=2)
	ax[1,0].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec='k',label='Field plots')
	yhat=rs['br'].Intercept+rs['br'].A*xhat+rs['br'].PD*pd
	ax[1,0].plot(xhat,yhat,'m-',label='EP964')
	#ax[1,0].plot(tv-2023,dCBM['Other S_ID_6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#x=2002-np.array([2000,1988,1949]); y=np.array([0,2.6,22.1+17.6])
	#ax[1,0].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	ax[1,0].set(xticks=xticks,ylabel='Branches (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,20])
	#ax[1,0].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_StemMerch']['Ensemble Mean'][:,0,0,0]+mos[pNam]['Scenarios'][iP]['Mean']['C_StemNonMerch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Csw L t0'][iFP],bw,bin)
	ax[1,1].errorbar(bin,mu,yerr=2*se,color='k',ls='',lw=0.5,capsize=2)
	ax[1,1].plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec='k',label='Field plots')
	yhat=rs['sw'].Intercept+rs['sw'].A*xhat+rs['sw'].PD*pd
	ax[1,1].plot(xhat,yhat,'m-',label='EP964')
	#ax[1,1].plot(tv-2023,dCBM['Roots S_ID-6'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#x=2002-np.array([2000,1988,1949]); y=np.array([0.1,2.8,32.6])
	#ax[1,1].plot(x,y,'r^',mfc='w',mec='r',lw=0.5,label='Fertster et al\n(2002)')
	ax[1,1].set(xticks=xticks,ylabel='Stemwood (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,60])
	#ax[1,1].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NonStemwoodBiomass_' + StatusTuned,'png',900)

	# With CBM25
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10)); xticks=np.arange(-2000,2220,20); xlim=[0,125]
	ax[0,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].plot(tv-2023,dCBM['Foilage S_ID 3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax[0,0].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[0,0].plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	yhat=rs['f'].Intercept+rs['f'].A*xhat+rs['f'].PD*pd
	ax[0,0].plot(xhat,yhat,'m-.',label='EP964')
	ax[0,0].set(xticks=xticks,yticks=np.arange(0,30,5),ylabel='Foliage biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,15])
	ax[0,0].legend(loc='center left',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='Branch (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP],bw,bin)
	ax[0,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[0,1].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Branch (field plots)')
	yhat=rs['br'].Intercept+rs['br'].A*xhat+rs['br'].PD*pd
	ax[0,1].plot(xhat,yhat,'m--',color=0.5*np.array(meta['Graphics']['gp']['cl2']),label='Branch (EP964)')
	ax[0,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k-',color=meta['Graphics']['gp']['cl3'],label='Bark (FCG25)')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	ax[0,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl3'],ls='',lw=0.5,capsize=2)
	ax[0,1].plot(bin,mu,'k^',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl3'],label='Bark (field plots)')
	yhat=rs['bk'].Intercept+rs['bk'].A*xhat+rs['bk'].PD*pd
	ax[0,1].plot(xhat,yhat,'m--',color=0.5*np.array(meta['Graphics']['gp']['cl3']),label='Bark (EP964)')
	ax[0,1].set(xticks=xticks,yticks=np.arange(0,30,5),ylabel='Bark and branches (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,15])
	ax[0,1].legend(loc='center left',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0]+mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].plot(tv-2023,dCBM['Other S_ID_3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP]+fp['Cbk L t0'][iFP],bw,bin)
	ax[1,0].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[1,0].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	yhat=(rs['br'].Intercept+rs['br'].A*xhat+rs['br'].PD*pd)+(rs['bk'].Intercept+rs['bk'].A*xhat+rs['bk'].PD*pd)
	ax[1,0].plot(xhat,yhat,'m-.',color=meta['Graphics']['gp']['cl3'],label='EP964')
	ax[1,0].set(xticks=xticks,ylabel='Bark+branches (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,30])
	ax[1,0].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,1].plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Root']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].plot(tv-2023,dCBM['Roots S_ID-3'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cr L t0'][iFP],bw,bin)
	ax[1,1].errorbar(bin,mu,yerr=2*se,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax[1,1].plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'],label='Field plots')
	ax[1,1].set(xticks=xticks,ylabel='Roots, total (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,40])
	ax[1,1].legend(loc='lower right',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NonStemwoodBiomass_WithCBM25','png',900)




	#--------------------------------------------------------------------------
	# Stemwood
	#--------------------------------------------------------------------------
	# Import TIPSY output
	dGY=gu.ReadExcel(r'D:\Modelling Projects\Demo_Ref_Underplant\Inputs\TIPSY output.xlsx')
	dGY['Merch Ratio']=np.nan_to_num(dGY['Net volume 125']/dGY['Volume total'])
	dGY['Stemwood Merch Biomass']=dGY['Merch Ratio']*dGY['Stemwood biomass']
	dGY['Stemwood Non-merch Biomass']=(1-dGY['Merch Ratio'])*dGY['Stemwood biomass']

	iGY2=np.where(dGY['Scenario']==2)[0]
	iGY8=np.where(dGY['Scenario']==8)[0]
	iT1=np.where( (tv>=1897) & (tv<=2016) )[0]
	iT2=np.where( (tv>=2024) )[0]

	# Boudwyn et al Lodgepole pine - Montane Cordillera
	b_m=0.76*dGY['Net volume 125'][iGY2]**0.93 # Merch (including stumps and tops)
	nonmerchfactor=0.85829+18.79727*b_m**-0.86026 # total/b_m
	b_n=nonmerchfactor*b_m-b_m

	b_n_FP=fp['Csw L t0']-fp['Csw125 L t0']
	iFP2=np.where( (np.isnan(b_n_FP)==False) & (b_n_FP>0) & (b_n_FP<1000) & (fp['Csw L t0']>0) & (fp['Ecozone BC L2']==meta['LUT']['GP']['Ecozone BC L2']['IDFdk']) )[0]

	#plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,9.5)); xlim=[0,125]
	ax[0,0].plot(dGY['Age'][iGY2],dGY['Net volume 125'][iGY2],'b-',label='TIPSY')
	ax[0,0].set(ylabel='Net volume 12.5cm\n(m$^3$ ha$^{-1}$)',xlabel='Age, years',xlim=xlim); ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

	ax[0,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood biomass'][iGY2],'b-',label='Total stemwood (TIPSY)')
	ax[0,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood Merch Biomass'][iGY2],'b--',label='Merch. stemwood (TIPSY)')
	ax[0,1].plot(A[iT1],dCBM['Biomass S_ID-2'][iT1],'c-',label='Total biomass (CBM25)')
	ax[0,1].plot(A[iT1],dCBM['Softwood Merchantable S_ID-2'][iT1],'c--',label='Merch stemwood (CBM25)')
	ax[0,1].plot(dGY['Age'][iGY2],0.5*b_m,'r--',label='Merch. stemwood\n(Boudwyn et al. 2008)')
	ax[0,1].set(ylabel='Biomass (tC ha$^{-1}$)',xlabel='Age, years',xlim=xlim); ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].legend(loc='upper left',fontsize=5,facecolor=[1,1,1],frameon=False)

	ax[1,0].plot([0,800],[0.53,0.53],'-',color=[0.8,0.8,0.8],lw=2,label='Constant wood density =0.5 ODT m$^3$')
	iT1a=np.where( (tv>=1897) & (tv<=2016) & (A>=1) & (A<=100) )[0]
	iGY2a=np.where( (dGY['Scenario']==2) & (dGY['Age']>=1) & (dGY['Age']<=101) )[0]
	ax[1,0].plot(dGY['Net volume 125'][iGY2a],dGY['Stemwood Merch Biomass'][iGY2a]/dGY['Net volume 125'][iGY2a],'--b',label='TIPSY')
	ax[1,0].plot(dGY['Net volume 125'][iGY2a],(2*dCBM['Softwood Merchantable S_ID-2'][iT1a])/dGY['Net volume 125'][iGY2a],'--c',label='CBM25')
	ax[1,0].plot(dGY['Net volume 125'][iGY2],b_m/dGY['Net volume 125'][iGY2],'r--',label='Boudwyn et al. (2007)')
	ax[1,0].set(ylabel='Merch. wood density\n(ODT m$^3$)',xlabel='Net volume 12.5 TIPSY',ylim=[0,4],xlim=xlim);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='upper right',fontsize=6,facecolor=[1,1,1],frameon=False)

	ax[1,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood Non-merch Biomass'][iGY2],'b-',label='TIPSY merch. stemwood')
	ax[1,1].plot(dGY['Age'][iGY2],0.5*b_n,'r--',label='Boudwyn et al. (2008)')
	bw=10; bin=np.arange(10,220+bw,bw);N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP2],b_n_FP[iFP2],bw,bin)
	ax[1,1].plot(bin,mu,'ko',ms=3,mfc='w',mew=0.5,label='Field plots')
	ax[1,1].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
	ax[1,1].set(ylabel='Non-merch stemwood\nbiomass (tC ha$^{-1}$)',xlabel='Age, years',xlim=xlim);
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Stemwood','png',900)

	return

#%%
def QA_BenchmarkCBM_FNM(meta,pNam,mos):
	zone='CWH'

	# Import field plot data
	meta,fp,soc=ufp.ImportFieldPlotData(meta,type='Stand',include_soil='False')
	dBGC=ufp.CalcStatsByBGC(meta,fp,soc)
	dAC=ufp.CalcStatsByAgeClass(meta,fp)

	# Import CBM data
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250109\CBM Outputs (Scenario_ID_5 & 6) revised.xls',sheet_name='CBM Outputs (tC)',0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250120\Demo_RefUnderFire_Scenario ID 2 & 3 CBM Outputs.xlsx',sheet_name='tC',0)
	#dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\20250130\CBM
	dCBM0=gu.ReadExcel(r'C:\Users\rhember\Government of BC\FCCS Team - Forest Carbon\Projects\BC-FCS\Quality Assurance and Benchmarks\R2025\04 Model-Model Comparison\Comparison with CBM-CFS3\20250130\CBM Outputs SC-ID 1 & 2.xlsx',sheet_name='With APP',skiprows=0)

	dCBM={}
	for k in dCBM0.keys():
		dCBM[k]=np.nan*np.ones(tv.size)
		ind1=np.where( (dCBM0['Years']>=tv[0]) & (dCBM0['Years']<=tv[-1]) )[0]
		ind2=np.where( (tv>=dCBM0['Years'][0]) & (tv<=dCBM0['Years'][-1]) )[0]
		dCBM[k][ind2]=dCBM0[k][ind1]
	dCBM['Year']=tv
	dCBM['Fire S_ID-1']=dCBM['NetCO2emissions_removals_CO2e S_ID-1']+3.667*dCBM['Net Ecosystem Productivity S_ID-1']
	dCBM['Fire S_ID-2']=dCBM['NetCO2emissions_removals_CO2e S_ID-2']+3.667*dCBM['Net Ecosystem Productivity S_ID-2']
	dCBM['NEB S_ID-1']=dCBM['NetCO2emissions_removals_CO2e S_ID-1']
	dCBM['NEB S_ID-2']=dCBM['NetCO2emissions_removals_CO2e S_ID-2']

	# Prep fcgadgets results
	iB=0
	iP=1
	A=mos[pNam]['Scenarios'][iP]['Mean']['A']['Ensemble Mean'][:,0,0,0]

	iFP=np.where( (fp['Ctot Net']>-20) & (fp['Ecozone BC L2']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) )[0]

	StatusTuned=''
	#StatusTuned='Tuned'

	#--------------------------------------------------------------------------
	# Carbon pools
	#--------------------------------------------------------------------------
	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,14)); xlim=[1985,2100]
	ax[0,0].plot(tv,dCBM['Biomass S_ID-2'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#ax[0,0].plot(tv,dCBM['Biomass S_ID-1'],'k--',color='k',label='CBM25')
	ax[0,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Biomass']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')

	# Add field plots
	#bw=10; bin=np.arange(10,110+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot L t0'][iFP],bw,bin)
	#ax[0,0].errorbar(bin+1897,mu,yerr=2*se,color='k',ls='',lw=0.25,capsize=2)
	#ax[0,0].plot(bin+1897,mu,'ko',ms=3,mfc='w',lw=0.25,label='Field plots')

	#ax[0,0].errorbar(bin+2017,mu,yerr=2*se,color='k',ls='',lw=0.25,capsize=2)
	#ax[0,0].plot(bin+2017,mu,'ko',ms=3,mfc='w',lw=0.25)

	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='Tree biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim)
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	y=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWood']['Ensemble Mean'][:,0,0,0]#+mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	ax[0,1].plot(tv,dCBM['Deadwood S_ID-2'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[0,1].set(xticks=np.arange(1800,2200,20),ylabel='Dead wood (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,200])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,0].plot(tv,dBGC['Soil'][zone]['SOC_ORG_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB\n(Shaw et al. 2018)')
	y=mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0]#-mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Litter']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].plot(tv,dCBM['Litter S_ID-2'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#ax[1,0].plot(tv,dCBM['Aboveground DOM S_ID-2'],'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	#a=mos[pNam]['Scenarios'][iP]['Mean']['C_DeadWoodDown']['Ensemble Mean'][:,0,0,0]
	#ax[1,0].plot(tv,dCBM['Litter S_ID-2']-,'k-.',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(tv,y,'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].set(xticks=np.arange(1800,2200,20),ylabel='Litter (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,200])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)

	ax[1,1].plot(tv,dBGC['Soil'][zone]['SOC_MIN_C_THA']['Mean']*np.ones(tv.size),'k-',lw=2,color=[0.8,0.8,0.8],label='Upland soil DB (Shaw et al. 2018)')
	ax[1,1].plot(tv,dCBM['Soil S_ID-2'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Soil']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=np.arange(1800,2200,20),ylabel='Soil (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,300])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)

	ax[2,0].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Piles']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[2,0].set(xticks=np.arange(1800,2200,20),ylabel='Piles (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim,ylim=[0,200])
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[2,1].plot(tv,dCBM['Total Ecosystem S_ID-2'],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[2,1].plot(tv,mos[pNam]['Scenarios'][iP]['Mean']['C_Forest']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[2,1].set(xticks=np.arange(1800,2200,20),ylabel='Total ecosystem (tC ha$^{-1}$)',xlabel='Time, years',xlim=xlim)
	ax[2,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Pools' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Net Growth
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
	ax.plot(tv[1:]-1988,np.diff(dCBM['Biomass S_ID-2']),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax.plot(tv-1988,mos[pNam]['Scenarios'][iP]['Mean']['C_G_Net']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	# Add field plots
	#bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Ctot Net'][iFP],bw,bin)
	#ikp=np.where(N>=20)[0]
	#ax.errorbar(bin[ikp],mu[ikp],yerr=2*se[ikp],color='k',ls='',lw=0.75,capsize=2)
	#ax.plot(bin[ikp],mu[ikp],'ko',ms=3,mec='k',mfc='w',mew=1)
	ax.set(xticks=np.arange(0,2200,20),yticks=np.arange(-10,10,0.5),ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[0,125],ylim=[0,7])
	ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NetGrowth' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Carbon fluxes
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,14))
	ax[0,0].plot(tv,dCBM['Net Ecosystem Productivity S_ID-2'],'.k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv-1,mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0],'.k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='Net ecoystem producitvity (tC/ha/yr)',xlabel='Time, years',xlim=xlim)
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Wildfire
	iT=np.where(tv==1915)[0][0]
	ax[0,1].bar(1,dCBM['Fire S_ID-2'][iT])
	ax[0,1].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_Wildfire']['Ensemble Mean'][iT,0,0,0])
	# Slashpile burn
	iT=np.where(tv==2023)[0][0]
	ax[1,0].bar(1,dCBM['Fire S_ID-2'][iT])
	ax[1,0].bar(2,mos[pNam]['Scenarios'][iP]['Mean']['E_Domestic_ForestSector_OpenBurning']['Ensemble Mean'][iT,0,0,0])
	# Harvest removals
	#ax[1,1].bar(1,np.nansum(dCBM['HR S_ID-2']))
	#ax[1,1].bar(2,np.sum(mos[pNam]['Scenarios'][iP]['Mean']['C_ToMillTotal']['Ensemble Mean'][:,0,0,0]))
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Fluxes','png',900)

	#--------------------------------------------------------------------------
	# Delta
	#--------------------------------------------------------------------------
	iT=np.where(tv>=2000)[0]
	plt.close('all');
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,8.5))
	d1=dCBM['Net Ecosystem Productivity S_ID-2']-dCBM['Net Ecosystem Productivity S_ID-1']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[0,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,0].plot(tv[iT]-1,d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ net ecoystem producitvity\n(tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['Net Ecosystem Productivity S_ID-2']-dCBM['Net Ecosystem Productivity S_ID-1']
	d2=(mos[pNam]['Scenarios'][iP]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iP]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])-(mos[pNam]['Scenarios'][iB]['Mean']['C_NPP']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['C_RH']['Ensemble Mean'][:,0,0,0])
	ax[0,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[0,1].plot(tv[iT]-1,np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[0,1].set(xticks=np.arange(1800,2200,20),ylabel='Cumulative $\Delta$ net ecoystem\nproducitvity (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[0,1].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-2']-dCBM['NEB S_ID-1']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',color=[0.8,0.8,0.8],lw=2)
	ax[1,0].plot(tv[iT],d1[iT],'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,0].plot(tv[iT],d2[iT],'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,0].set(xticks=np.arange(1800,2200,20),ylabel='$\Delta$ emissions (tCO$_2$e ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xlim=[2000,2150])
	ax[1,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	d1=dCBM['NEB S_ID-2']-dCBM['NEB S_ID-1']
	d2=mos[pNam]['Scenarios'][iP]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]-mos[pNam]['Scenarios'][iB]['Mean']['E_NEB']['Ensemble Mean'][:,0,0,0]
	ax[1,1].plot([0,2250],[0,0],'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,1].plot([2050,2050],[1000,-1000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax[1,1].plot(tv[iT],np.cumsum(d1[iT]),'k-',color=meta['Graphics']['gp']['cl1'],label='CBM25')
	ax[1,1].plot(tv[iT],np.cumsum(d2[iT]),'k--',color=meta['Graphics']['gp']['cl2'],label='FCG25')
	ax[1,1].set(xticks=np.arange(1800,2200,10),ylabel='Cumulative $\Delta$ emissions\n(tCO$_2$e ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2055],ylim=[-100,20])
	ax[1,1].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	#iT=np.where( (tv>=2023) & (tv<=2023+30) )[0]
	iT=np.where( (tv>=2023) & (tv<=2023+100) )[0]
	s1=np.sum(d1[iT])
	s2=np.sum(d2[iT])
	print(int(s1))
	print(int(s2))
	print(int(s2-s1))
	print(int((s2-s1)/s1*100))

	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Delta' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Non-stemwood biomass pools
	#--------------------------------------------------------------------------
	#plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,6)); xticks=np.arange(0,2220,10)
	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Foliage']['Ensemble Mean'][:,0,0,0],'k-',color=meta['Graphics']['gp']['cl2'],label='Foliage')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cf L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=sig,color=meta['Graphics']['gp']['cl2'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ko',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl2'])

	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Branch']['Ensemble Mean'][:,0,0,0],'k--',color=meta['Graphics']['gp']['cl1'],label='Branch')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbr L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=sig,color=meta['Graphics']['gp']['cl1'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'ks',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl1'])

	ax.plot(tv-2023,mos[pNam]['Scenarios'][iP]['Mean']['C_Bark']['Ensemble Mean'][:,0,0,0],'k-.',color=meta['Graphics']['gp']['cl3'],label='Bark')
	bw=10; bin=np.arange(10,200+bw,bw); N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP],fp['Cbk L t0'][iFP],bw,bin)
	ax.errorbar(bin,mu,yerr=sig,color=meta['Graphics']['gp']['cl3'],ls='',lw=0.5,capsize=2)
	ax.plot(bin,mu,'k^',ms=3,lw=0.5,mfc='w',mec=meta['Graphics']['gp']['cl3'])

	ax.set(xticks=xticks,ylabel='Biomass (tC ha$^{-1}$)',xlabel='Time, years',xlim=[0,125],ylim=[0,25])
	ax.legend(loc='upper left',fontsize=6,facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_NonStemwoodBiomass_' + StatusTuned,'png',900)

	#--------------------------------------------------------------------------
	# Stemwood
	#--------------------------------------------------------------------------
	flg=0
	if flg==1:
		# Import TIPSY output
		dGY=gu.ReadExcel(r'D:\Modelling Projects\Demo_Ref_Underplant\Inputs\TIPSY output.xlsx')
		dGY['Merch Ratio']=np.nan_to_num(dGY['Net volume 125']/dGY['Volume total'])
		dGY['Stemwood Merch Biomass']=dGY['Merch Ratio']*dGY['Stemwood biomass']
		dGY['Stemwood Non-merch Biomass']=(1-dGY['Merch Ratio'])*dGY['Stemwood biomass']
	
		iGY2=np.where(dGY['Scenario']==2)[0]
		iGY8=np.where(dGY['Scenario']==8)[0]
		iT1=np.where( (tv>=1897) & (tv<=2016) )[0]
		iT2=np.where( (tv>=2024) )[0]
	
		# Boudwyn et al Lodgepole pine - Montane Cordillera
		b_m=0.76*dGY['Net volume 125'][iGY2]**0.93 # Merch (including stumps and tops)
		nonmerchfactor=0.85829+18.79727*b_m**-0.86026 # total/b_m
		b_n=nonmerchfactor*b_m-b_m
	
		b_n_FP=fp['Csw L t0']-fp['Csw125 L t0']
		iFP2=np.where( (np.isnan(b_n_FP)==False) & (b_n_FP>0) & (b_n_FP<1000) & (fp['Csw L t0']>0) & (fp['Ecozone BC L2']==meta['LUT']['GP']['Ecozone BC L2']['IDFdk']) )[0]
	
		#plt.close('all');
		fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,8)); xlim=[0,125]
		ax[0,0].plot(dGY['Age'][iGY2],dGY['Net volume 125'][iGY2],'b-',label='TIPSY')
		ax[0,0].set(ylabel='Net volume 12.5cm\n(m$^3$ ha$^{-1}$)',xlabel='Age, years',xlim=xlim); ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False)
	
		ax[0,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood biomass'][iGY2],'b-',label='Total stemwood (TIPSY)')
		ax[0,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood Merch Biomass'][iGY2],'b--',label='Merch. stemwood (TIPSY)')
		ax[0,1].plot(A[iT1],dCBM['Biomass S_ID-1'][iT1],'c-',label='Total biomass (CBM25)')
		ax[0,1].plot(A[iT1],dCBM['Softwood Merchantable S_ID-1'][iT1],'c--',label='Merch stemwood (CBM25)')
		ax[0,1].plot(dGY['Age'][iGY2],0.5*b_m,'r--',label='Merch. stemwood\n(Boudwyn et al. 2008)')
		ax[0,1].set(ylabel='Biomass (tC ha$^{-1}$)',xlabel='Age, years',xlim=xlim); ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax[0,1].legend(loc='upper left',fontsize=5,facecolor=[1,1,1],frameon=False)
	
		ax[1,0].plot([0,800],[0.53,0.53],'-',color=[0.8,0.8,0.8],lw=2,label='Constant wood density =0.5 ODT m$^3$')
		iT1a=np.where( (tv>=1897) & (tv<=2016) & (A>=1) & (A<=100) )[0]
		iGY2a=np.where( (dGY['Scenario']==2) & (dGY['Age']>=1) & (dGY['Age']<=101) )[0]
		ax[1,0].plot(dGY['Net volume 125'][iGY2a],dGY['Stemwood Merch Biomass'][iGY2a]/dGY['Net volume 125'][iGY2a],'--b',label='TIPSY')
		ax[1,0].plot(dGY['Net volume 125'][iGY2a],(2*dCBM['Softwood Merchantable S_ID-1'][iT1a])/dGY['Net volume 125'][iGY2a],'--c',label='CBM25')
		ax[1,0].plot(dGY['Net volume 125'][iGY2],b_m/dGY['Net volume 125'][iGY2],'r--',label='Boudwyn et al. (2007)')
		ax[1,0].set(ylabel='Merch. wood density\n(ODT m$^3$)',xlabel='Net volume 12.5 TIPSY',ylim=[0,4],xlim=xlim);
		ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax[1,0].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
	
		ax[1,1].plot(dGY['Age'][iGY2],0.5*dGY['Stemwood Non-merch Biomass'][iGY2],'b-',label='TIPSY merch. stemwood')
		ax[1,1].plot(dGY['Age'][iGY2],0.5*b_n,'r--',label='Boudwyn et al. (2008)')
		bw=10; bin=np.arange(10,220+bw,bw);N,mu,med,sig,se=gu.discres(fp['Age Mean t0'][iFP2],b_n_FP[iFP2],bw,bin)
		ax[1,1].plot(bin,mu,'ko',ms=3,mfc='w',mew=0.5,label='Field plots')
		ax[1,1].legend(loc='upper right',fontsize=5,facecolor=[1,1,1],frameon=False)
		ax[1,1].set(ylabel='Non-merch stemwood\nbiomass(tC ha$^{-1}$)',xlabel='Age, years',xlim=xlim);
		ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	
		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BenchmarkCBM_Stemwood','png',900)

	return

#%%
def QA_Benchmark_HWP_ResidenceTime(meta):
	iScn=3
	C_in=np.sum(mos[pNam]['Scenarios'][iScn]['Mean']['C_ToLumber']['Ensemble Mean'][:,0,0,0]+ \
		mos[pNam]['Scenarios'][iScn]['Mean']['C_ToOSB']['Ensemble Mean'][:,0,0,0]+ \
		mos[pNam]['Scenarios'][iScn]['Mean']['C_ToMDF']['Ensemble Mean'][:,0,0,0]+ \
		mos[pNam]['Scenarios'][iScn]['Mean']['C_ToPlywood']['Ensemble Mean'][:,0,0,0])
	C_hwp=mos[pNam]['Scenarios'][iScn]['Mean']['C_InUse']['Ensemble Mean'][:,0,0,0]+mos[pNam]['Scenarios'][iScn]['Mean']['C_WasteSystems']['Ensemble Mean'][:,0,0,0]
	C_use=mos[pNam]['Scenarios'][iScn]['Mean']['C_InUse']['Ensemble Mean'][:,0,0,0]
	iT1=np.where(tv==meta[pNam]['Project']['Year Project'])[0]
	iT2=np.where(tv==meta[pNam]['Project']['Year Project']+100)[0]
	r=C_use/C_use[iT1]*100

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,5.5))
	ax.plot(tv,r,'k-')
	#plt.plot(tv,C_hwp/C_hwp[iT1],'k-')
	ax.plot([tv[iT2],tv[iT2]],[0,100],'r--')
	ax.set(ylabel='Percent remaing (%)',xlabel='Time, years',xlim=[2000,2150],ylim=[0,100])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Benchmark_HWP_ResidenceTime100','png',900)
	return

#%%
def QA_AssessBiomassProportions(meta,mos):
	iScn=3
	iT=np.where(mos[pNam]['Scenarios'][iScn]['Mean']['A']['Ensemble Mean'][:,0,0,0]==50)[0]
	mu={}
	mu['Bark']=np.mean(mos[pNam]['Scenarios'][iScn]['Mean']['C_Bark']['Ensemble Mean'][iT,0,0,0])
	mu['Branch']=np.mean(mos[pNam]['Scenarios'][iScn]['Mean']['C_Branch']['Ensemble Mean'][iT,0,0,0])
	mu['Foliage']=np.mean(mos[pNam]['Scenarios'][iScn]['Mean']['C_Foliage']['Ensemble Mean'][iT,0,0,0])
	mu['Stemwood']=np.mean(mos[pNam]['Scenarios'][iScn]['Mean']['C_StemMerch']['Ensemble Mean'][iT,0,0,0]+mos[pNam]['Scenarios'][iScn]['Mean']['C_StemNonMerch']['Ensemble Mean'][iT,0,0,0])
	mu['Root']=np.mean(mos[pNam]['Scenarios'][iScn]['Mean']['C_Root']['Ensemble Mean'][iT,0,0,0])
	Ctot=mu['Bark']+mu['Branch']+mu['Foliage']+mu['Stemwood']
	muP={}
	muP['Bark']=int(mu['Bark']/Ctot*100)
	muP['Branch']=int(mu['Branch']/Ctot*100)
	muP['Foliage']=int(mu['Foliage']/Ctot*100)
	muP['Stemwood']=int(mu['Stemwood']/Ctot*100)
	return

#%%
def Plot_QA_FateOfFelled(meta,d,iS,iH,TimeHorizon,dpi):
	fnam='Figure ' + str(meta['Graphics']['Fig Count']) + ' Schematic of Fibre Flow'
	g=graphviz.Digraph('Fibre Flow',filename=fnam,
					 graph_attr={'rankdir':'LR',
								 'compound':'true',
								 'pad':'0.1',
								 'nodesep':'0.1',
								 'format':'png',
								 'dpi':dpi},
					 node_attr={'shape':'box',
								'fixedsize':'true',
								'width':'1.2',
								'height':'0.5',
								'fontname':meta['Graphics']['Flowchart']['Font Name'],
								'fontcolor':meta['Graphics']['Flowchart']['Font Color'],
								'fontsize':meta['Graphics']['Flowchart']['Font Size'],
								'style':'filled',
								'fillcolor':meta['Graphics']['Flowchart']['Node Background Color'],
								'penwidth':meta['Graphics']['Flowchart']['Penwidth'],
								'center':'true'},
					edge_attr={})
	
	g.edge_attr.update(arrowhead='normal',arrowsize='0.6',penwidth='0.5',fontsize='9',fontname=meta['Graphics']['Flowchart']['Font Name'])
	
	tfb='Biomass of\nfelled trees'
	r='Root\nmortality'
	m='Merchantable\nwood'
	nm='Non-merchant.\naboveground\nbiomass'
	lds='Left\ndispersed'
	pil='Piles'
	rpil='Remnant\npiles'
	lex='Log\nexports'
	fw='Firewood\ncollection'
	mill='Mill\ntransport'
	lum='Lumber'
	pp='Pulp and Paper'
	pan='Panels'
	pel='Pellets'
	bef='Bioenergy facility'
	#rh='Residual\nharvest'
	dom='Dead organic\nmatter'
	atm='Atmosphere'
	ob='Open\nburning'
	rhe='Decay'
	rhp='Product\ndecay'
	combp='Product\nburning'
	res='Residential'

	Felled=d['C_Felled'][iH,iS]
	FelledM=d['C_FelledMerch'][iH,iS]
	FelledR=d['C_FelledRoots'][iH,iS]
	FelledNM=Felled-FelledM-FelledR
	print(FelledNM)

	dDW=(d['C_DeadWood'][iH,iS])-(d['C_DeadWood'][iH-1,iS])
	ToDOM=d['C_ToDOM'][iH,iS]
	LdsToDom=ToDOM-FelledR
	DomToAtm=np.sum(d['C_RH_FelledRootsDispersed'][iH:iH+TimeHorizon,iS])

	ToPile=d['C_ToPile'][iH,iS]
	ToPileM=d['C_ToPileMerch'][iH,iS]
	ToPileNM=ToPile-ToPileM
	print(ToPileNM)

	ToMillSS=d['C_ToMillMerchDead'][iH,iS]
	ToMillM=d['C_ToMillMerchGreen'][iH,iS]
	ToMillNM=d['C_ToMillNonMerchGreen'][iH,iS]

	ToFW=d['C_ToBBP_FirewoodDom'][iH,iS]
	BurnTot=np.sum(d['C_ToPileBurnTot'][iH:iH+1,iS])
	RemnantPile=ToPile-BurnTot
	#BurnNM=d['C_ToPileBurnTot'][iH+1,iS]-d['C_ToPileBurnMerch'][iH+1,iS]
	ToLogEx=d['C_ToLogExport'][iH,iS]
	ToLum=d['C_ToLumber'][iH,iS]
	ToAtm=0

	#print(Felled-ToDOM-ToPile-ToMillM-ToMillNM-ToMillSS-ToFW) # Confirm that all felled fibre is accounted for

	tfb2r=str(np.round(FelledR/Felled*100,decimals=1)) + '%'
	tfb2m=str(np.round(FelledM/Felled*100,decimals=1)) + '%\l'
	tfb2nm=str(np.round(FelledNM/Felled*100,decimals=1)) + '%\l'

	r2dom=str(np.round(FelledR/Felled*100,decimals=1)) + '%'

	m2lds=str(np.round((FelledM-ToMillM-ToPileM)/Felled*100,decimals=1)) + '%\l'
	m2pb=str(np.round(ToPileM/Felled*100,decimals=1)) + '%\l'
	m2mill=str(np.round(ToMillM/Felled*100,decimals=1)) + '%\l'
	m2fw=str(np.round(ToFW/Felled*100,decimals=1)) + '%\l'
	m2lex=str(np.round(ToLogEx/Felled*100,decimals=1)) + '%\l'

	nm2lds=str(np.round((FelledNM-ToPileNM-ToMillNM)/Felled*100,decimals=1)) + '%\l'
	nm2pb=str(np.round(ToPileNM/Felled*100,decimals=1)) + '%\l'
	nm2mill=str(np.round(ToMillNM/Felled*100,decimals=1)) + '%\l'

	lds2dom=str(np.round(LdsToDom/Felled*100,decimals=1)) + '%\l'

	mill2lum=str(np.round(ToLum/Felled*100,decimals=1)) + '%\l'
	mill2pp=str(np.round(d['C_ToPaper'][iH,iS]/Felled*100,decimals=1)) + '%\l'
	mill2pan=str(np.round((d['C_ToPlywood'][iH,iS]+d['C_ToMDF'][iH,iS]+d['C_ToOSB'][iH,iS])/Felled*100,decimals=1)) + '%\l'
	mill2pel=str(np.round((d['C_ToBBP_PelletExport'][iH,iS]+d['C_ToBBP_PelletDomRNG'][iH,iS]+d['C_ToBBP_PelletDomGrid'][iH,iS])/Felled*100,decimals=1)) + '%\l'
	mill2bef=str(np.round((d['C_ToBBP_PowerGrid'][iH,iS]+d['C_ToBBP_PowerFacilityExport'][iH,iS]+d['C_ToBBP_PowerFacilityDom'][iH,iS])/Felled*100,decimals=1)) + '%\l'

	fw2atm=str(np.round(ToFW/Felled*100,decimals=1)) + '%\l'

	pil2ob=str(np.round(BurnTot/Felled*100,decimals=1)) + '%\l'
	pil2rpil=str(np.round(RemnantPile/Felled*100,decimals=1)) + '%\l'
	b=BurnTot/Felled*100
	ob2atm=str(np.round(b,decimals=1)) + '%\l'
	ToAtm=ToAtm+b

	f=(RemnantPile-d['C_Piles'][iH+TimeHorizon,iS])/RemnantPile
	rpil2rhe=f*(RemnantPile/Felled*100)
	rpil2rhe_str=str(np.round(rpil2rhe,decimals=1)) + '%\l'

	Frac_To_DOM=FelledR/Felled+LdsToDom/Felled
	dom2rhe=(DomToAtm/ToDOM)*Frac_To_DOM*100
	dom2rhe_str=str(np.round(dom2rhe,decimals=1)) + '%'
	rhe2atm=str(np.round(dom2rhe+rpil2rhe,decimals=1)) + '%'
	ToAtm=ToAtm+dom2rhe+rpil2rhe

	a=np.sum(d['E_Domestic_ForestSector_HWP'][iH:iH+TimeHorizon,iS])/3.667
	b=a/Felled*100
	ToAtm=ToAtm+b
	rhp2atm=str(np.round(b,decimals=1)) + '%'

	a=np.sum(d['E_Domestic_Bioenergy'][iH:iH+TimeHorizon,iS])/3.667
	b=a/Felled*100
	ToAtm=ToAtm+b
	combp2atm=str(np.round(b,decimals=1)) + '%'

	with g.subgraph(name='cluster_eco') as eco:
		eco.attr(label='Ecosystem',labeljust="l",style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		eco.node(tfb)
		eco.node(dom)
		eco.node(lds,shape='circle',width='0.9')
		eco.node(pil)
		eco.node(rpil)
		eco.node(ob,shape='circle',width='0.9')
		eco.node(rhe,shape='circle',width='0.9')
		eco.node(r,shape='circle',width='0.9')
		eco.node(nm)
		eco.node(m)

	g.node(atm,color='#578f91',fillcolor='#a2f9fc',
		fontcolor='#274445')

	with g.subgraph(name='cluster_mills') as mil:
		mil.attr(label='Harvested Wood Products',labeljust="r",style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		mil.node(lex,shape='circle',width='0.9')
		mil.node(fw,shape='circle',width='0.9')
		mil.node(mill,shape='circle',width='0.9')
		mil.node(lum) #,shape='house'
		mil.node(pp)
		mil.node(pan)
		mil.node(pel)
		mil.node(bef)
		mil.node(res)
		#mil.node(bef,shape='house')
		mil.node(combp,shape='circle',width='0.9')
		mil.node(rhp,shape='circle',width='0.9')

	g.edge(tfb,m,label=tfb2m)
	g.edge(tfb,nm,label=tfb2nm)
	g.edge(tfb,r,label=tfb2r)
	
	g.edge(r,dom,label=r2dom)
	g.edge(dom,rhe,label=dom2rhe_str)
	g.edge(rhe,atm,label=rhe2atm)

	g.edge(rpil,rhe,label=rpil2rhe_str)

	g.edge(m,lds,label=m2lds)
	g.edge(m,pil,label=m2pb)

	g.edge(nm,lds,label=nm2lds)
	g.edge(nm,pil,label=nm2pb)
	g.edge(nm,mill,label=nm2mill)
	
	g.edge(pil,ob,label=pil2ob)
	g.edge(pil,rpil,label=pil2rpil)
	g.edge(ob,atm,label=ob2atm)
	g.edge(lds,dom,label=lds2dom)

	g.edge(m,lex,label=m2lex)
	g.edge(m,mill,label=m2mill)
	g.edge(m,fw,label=m2fw)
	g.edge(fw,res,label=m2fw)
	g.edge(mill,lum,label=mill2lum)
	g.edge(mill,pp,label=mill2pp)
	g.edge(mill,pan,label=mill2pan)
	g.edge(mill,pel,label=mill2pel)
	g.edge(mill,bef,label=mill2bef)

	g.edge(bef,combp,label='')
	g.edge(pel,combp,label='')
	g.edge(pp,combp,label='')
	g.edge(res,combp,label='')

	g.edge(lum,rhp,label='')
	g.edge(pp,rhp,label='')
	g.edge(pan,rhp,label='')

	g.edge(rhp,atm,label=rhp2atm)

	g.edge(combp,atm,label=combp2atm)

	g.node(atm,label='Atmosphere\n' + str(np.round(ToAtm,decimals=1)) + '%')

	display(g)
	g.format='png'
	#g.dpi=72 # Controlled by "graph_attr" above
	if meta['Graphics']['Print Figures']=='On':
		g.render(filename=meta['Graphics']['Print Figure Path'] + '\\HarvestSchematic')
	return

#%%