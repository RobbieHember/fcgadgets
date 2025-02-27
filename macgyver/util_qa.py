#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
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

