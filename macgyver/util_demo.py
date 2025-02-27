'''
DEMO UTILITIES
'''
#%% Import Python modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
from scipy import stats
import openpyxl
from datetime import date
from matplotlib.patches import Rectangle
from fcgadgets.macgyver import util_general as gu
from fcgadgets.cbrunner import cbrun_util as cbu

#%%
def GetSingleEnsembleResults(meta,pNam):
	v0=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		v0.append(cbu.LoadSingleOutputFile(meta,iScn,0,0))
	return v0

# #%%
# def CalculateAggregateVariables(meta,v1):
# 	for iScn in range(meta[pNam]['Project']['N Scenario']):
# 		# Calculate carbon content of dead wood, organic and mineral soil horizons following Shaw et al. (2017)
# 		v1[iScn]['SoilOrgH']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterVF']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterS']]
# 		v1[iScn]['SoilMinH']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilVF']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilF']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SoilS']]
# 		v1[iScn]['DeadWood']=v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SnagStem']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['SnagBranch']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterM']]+ \
# 			v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['LitterF']]

# 		v1[iScn]['C_BiomassAG']=np.sum(v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['BiomassAboveground']],axis=2)
# 		v1[iScn]['C_BiomassBG']=np.sum(v1[iScn]['C_Eco_Pools'][:,:,meta['Core']['iEP']['BiomassBelowground']],axis=2)
# 		v1[iScn]['C_BiomassAG'][np.isnan(v1[iScn]['C_BiomassAG'])]=0

# 		v1[iScn]['C_Eco']=np.sum(v1[iScn]['C_Eco_Pools'],axis=2)
# 		v1[iScn]['C_Eco']=np.sum(v1[iScn]['C_Eco_Pools'],axis=2)
# 		v1[iScn]['C_Eco'][np.isnan(v1[iScn]['C_Eco'])]=0

# 		#v1[iScn]['Sum']['C_Forest']=v1[iScn]['Sum']['C_Biomass']+v1[iScn]['Sum']['C_DeadWood']+v1[iScn]['Sum']['C_Litter']+v1[iScn]['Sum']['C_Soil']

# 	return v1

#%%
def ExportTableByScenario(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	# Key word arguments
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']
		iSS=kwargs['iSS']
		iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'multip' in kwargs.keys():
		multip=kwargs['multip']
	else:
		multip=1.0

	if 'Variables' in kwargs.keys():
		vL=kwargs['variables']
	else:
		vL=mos[pNam]['Scenarios'][0]['Mean'].keys()
	d0=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for v in vL:

			if v not in mos[pNam]['Scenarios'][iScn]['Mean'].keys():
				continue

			if kwargs['OperTime']=='Mean':
				d[v]=np.round(multip*np.mean(mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]),decimals=2)
			else:
				d[v]=np.round(multip*np.sum(mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]),decimals=2)

		if iScn==0:
			df=pd.DataFrame().from_dict(d,orient='index');
		else:
			df0=pd.DataFrame().from_dict(d,orient='index');
			df=pd.concat([df,df0],axis=1);
		d0.append(d)

	df.columns=[np.arange(1,df.columns.size+1)];

	if kwargs['Save']=='On':
		fout=meta['Paths'][pNam]['Data'] + '\\Outputs\\ByScenario_' + kwargs['NameTable'] + '_' + str(kwargs['t0']) + '-' + str(kwargs['t1']) + '_ProjectStrat' + str(iPS) + '_SpatialStrat' + str(iSS) + '_TimeStrat' + str(iYS) + '.xlsx'
		df.to_excel(fout)
	#gu.PrintDF(df,fout,SheetName='Sheet1')

	return df,d0

#%% Custom scenario comparison tabular export
# Exmaple:
# udem.ExportTableDelta(meta,pNam,mos,'ScnComp1','Mean','Mean',[2020,2030])
def ExportTableDelta(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	iPS=0; iSS=0; iYS=0

	OperSpace=kwargs['OperSpace']

	if 'Variables' in kwargs.keys():
		vL=kwargs['Variables']
	else:
		vL=mos[pNam]['Delta'][ list(mos[pNam]['Delta'].keys())[0] ]['ByStrata'][OperSpace].keys()

	df=pd.DataFrame()
	Names=[]
	for sc in kwargs['cNam']:

		Names.append(sc)

		d={}
		for v in vL:
			if v not in mos[pNam]['Delta'][sc]['ByStrata'][OperSpace].keys():
				continue
			if (kwargs['Units']=='Actual') | (kwargs['Units']=='All'):
				if kwargs['OperTime']=='Mean':
					y=np.mean(mos[pNam]['Delta'][sc]['ByStrata'][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				else:
					y=np.sum(mos[pNam]['Delta'][sc]['ByStrata'][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				d[v]=np.round(y,decimals=2)
			if (kwargs['Units']=='Relative') | (kwargs['Units']=='All'):
				if kwargs['OperTime']=='Mean':
					y=np.mean(mos[pNam]['Delta'][sc]['ByStrata'][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.mean(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				else:
					y=np.sum(mos[pNam]['Delta'][sc]['ByStrata'][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.sum(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][OperSpace][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				d[v + ' (%)']=np.round(y,decimals=2)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)
	df.index.name='Variable'
	df.columns=Names
	#df=df.sort_index(axis=0)
	if kwargs['Save']=='On':
		df.to_excel(meta['Paths'][pNam]['Data'] + '\\Outputs\\Delta_' + kwargs['NameTable'] + '_' + kwargs['OperSpace'] + '_' + kwargs['OperTime'] + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']) + '.xlsx')
	return df

#%%
def ExportTableScenariosAndDelta(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	if 'Variables' in kwargs.keys():
		vL=kwargs['Variables']
	else:
		vL=mos[pNam]['Delta'][kwargs['cNam'][0]]['ByStrata']['Mean'].keys()

	if 'Decimals' in kwargs.keys():
		dec=kwargs['Decimals']
	else:
		dec=2

	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0;

	if 'Multiplier' in kwargs.keys():
		multi=kwargs['Multiplier']
	else:
		multi=1.0

	df=pd.DataFrame()
	for sc in mos[pNam]['Delta'].keys():

		if 'cNam' in kwargs.keys():
			if np.isin(sc,kwargs['cNam'])==False:
				continue

		d={}
		for v in vL:
			iScn=mos[pNam]['Delta'][sc]['iB']
			if v not in mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']].keys():
				continue
			if kwargs['OperTime']=='Mean':
				y=multi*np.mean(mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
			else:
				y=multi*np.sum(mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
			d[v]=np.round(y,decimals=dec)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)

		d={}
		for v in vL:
			iScn=mos[pNam]['Delta'][sc]['iP']
			if v not in mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']].keys():
				continue
			if kwargs['OperTime']=='Mean':
				y=multi*np.mean(mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
			else:
				y=multi*np.sum(mos[pNam]['Scenarios'][iScn][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
			d[v]=np.round(y,decimals=dec)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)

		d={}
		for v in vL:
			if v not in mos[pNam]['Delta'][sc]['ByStrata'][kwargs['OperSpace']].keys():
				continue
			if (kwargs['Units']=='Actual') | (kwargs['Units']=='All'):
				if kwargs['OperTime']=='Mean':
					y=multi*np.mean(mos[pNam]['Delta'][sc]['ByStrata'][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				else:
					y=multi*np.sum(mos[pNam]['Delta'][sc]['ByStrata'][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				d[v]=np.round(y,decimals=dec)
			if (kwargs['Units']=='Relative') | (kwargs['Units']=='All'):
				if kwargs['OperTime']=='Mean':
					y=multi*np.mean(mos[pNam]['Delta'][sc]['ByStrata'][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.mean(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				else:
					y=multi*np.sum(mos[pNam]['Delta'][sc]['ByStrata'][kwargs['OperSpace']][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.sum(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				d[v + ' (%)']=np.round(y,decimals=dec)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)
	
	df.index.name='Variable'
	df.columns=['Baseline','Action','Delta']
	#df=df.sort_index(axis=0)
	if kwargs['Save']=='On':
		df.to_excel(meta['Paths'][pNam]['Data'] + '\\Outputs\\TabularDelta_' + kwargs['NameTable'] + '_' + kwargs['OperTime'] + 'Time_' + kwargs['OperSpace'] + 'Space_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']) + '.xlsx')
	return df

#%% Plot mean fluxes and mean pools over a specified time horizon
def PlotSchematicBalance(meta,pNam,mos,**kwargs):

	# Key word arguments
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iSP']
		iSS=kwargs['iSS']
		iYS=kwargs['iYS']
	else:
		iPS=0
		iSS=0
		iYS=0
	iB=mos[pNam]['Delta'][ kwargs['cNam'] ]['iB']
	iP=mos[pNam]['Delta'][ kwargs['cNam'] ]['iP']
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	y_b={}
	y_p={}
	y_d={}
	for k in mos[pNam]['Scenarios'][0]['Mean'].keys():

		if (k=='C_Forest'):
			y_b[k]=mos[pNam]['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[0],iPS,iSS,iYS]
			y_p[k]=mos[pNam]['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[0],iPS,iSS,iYS]
			y_d[k]=y_p[k]-y_b[k]
		elif (k=='C_HWP'):
			# Reference to the year before the project, otherwise, HWP will go down following harvest
			y_b[k]=mos[pNam]['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT[0]-1,iPS,iSS,iYS]
			y_p[k]=mos[pNam]['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT[0]-1,iPS,iSS,iYS]
			y_d[k]=y_p[k]-y_b[k]
		else:
			y_b[k]=np.sum(mos[pNam]['Scenarios'][iB]['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS])
			y_p[k]=np.sum(mos[pNam]['Scenarios'][iP]['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS])
			y_d[k]=y_p[k]-y_b[k]

		# Round
		y_b[k]=y_b[k].astype(int)
		y_p[k]=y_p[k].astype(int)
		y_d[k]=y_d[k].astype(int)

	bx_ec='none'
	bx_fs=9
	bx_fc=[0.93,0.93,0.93]
	bx2_fc=[0.9,0.9,0.9]
	bx_lower_h=0.47
	bx_ForestSector_Domestic_w=0.48
	bx_esc_w=0.15
	bx_atmo_bottom=0.88
	arrow_head_w=0.007
	arrow_lw=0.05
	fs_flux=6.5
	decim=0
	cl_f1=[0.08,0.3,0.55] #cl_f1=[0.2,0.7,0]
	cl_f2=[0.5,0,0]

	def GetSign(y):
		if y>0:
			x='+'
		else:
			x=''
		return x

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(23,10))

	# Background
	#ax.add_patch(Rectangle([0,1],0,1,fc=[0.9,0.9,0.6],ec='none'))

	# Atmosphere
	ax.add_patch(Rectangle([0.01,bx_atmo_bottom],0.98,0.1,fc=[0.9,0.95,1],ec=bx_ec))
	ax.text(0.5,0.935,'Atmosphere',size=bx_fs,ha='center',fontweight='bold',color=[0.08,0.3,0.55])
	vr='E_NAB'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Change in storage (tCO$_2$e): ' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.5,0.9,txt,size=fs_flux+1,ha='center')

	# Non-forest
	ax.add_patch(Rectangle([-0.24,0.01],0.17,bx_lower_h,fc=pltc.to_rgb('#c2a9d6'),ec=bx_ec))
	ax.text(-0.155,0.04,'Non-Forest',size=bx_fs,ha='center',color=pltc.to_rgb('#4e3c5c'),fontweight='bold')

	# Non-forest land
	ax.add_patch(Rectangle([-0.23,0.28],0.15,0.17,fc=pltc.to_rgb('#f7edff'),ec=bx_ec))
	ax.text(-0.155,0.33,'Non-forest\nland',size=bx_fs,ha='center',color=pltc.to_rgb('#4e3c5c'))
	#ax.text(-0.13,0.26,'Change in storage (tCO$_2$e):',size=fs_flux,ha='center')
	vr='C_Forest'
	a1=np.round(y_b[vr]*3.667,decimals=decim); a2=np.round(y_p[vr]*3.667,decimals=decim); a3=np.round(y_d[vr]*3.667,decimals=decim)
	txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	#ax.text(-0.12,0.23,txt,size=fs_flux,ha='center')

	ax.arrow(-0.07,0.43,0.07,0.0,head_width=1.7*arrow_head_w,head_length=0.007,fc=cl_f1,ec=cl_f1,lw=arrow_lw)
	ax.text(-0.03,0.37,'Afforestation\n0,0 (0)',color=cl_f1,size=fs_flux,ha='center')

	ax.arrow(0.01,0.32,-0.07,0.0,head_width=1.7*arrow_head_w,head_length=0.007,fc=cl_f2,ec=cl_f2,lw=arrow_lw)
	ax.text(-0.03,0.26,'Deforestation\n0,0 (0)',color=cl_f2,size=fs_flux,ha='center')

	ax.arrow(0.01,0.16,-0.07,0.0,head_width=1.7*arrow_head_w,head_length=0.007,fc=cl_f2,ec=cl_f2,lw=arrow_lw)
	ax.text(-0.03,0.075,'Sediment\nloss\n0,0 (0)',color=cl_f2,size=fs_flux,ha='center')

	# Freshwater and ocean
	ax.add_patch(Rectangle([-0.23,0.08],0.15,0.17,fc=pltc.to_rgb('#f7edff'),ec=bx_ec))
	ax.text(-0.155,0.14,'Water\nsystems',size=bx_fs,ha='center',color=pltc.to_rgb('#4e3c5c'))
	#ax.text(-0.13,0.26,'Change in storage (tCO$_2$e):',size=fs_flux,ha='center')
	vr='C_Forest'
	a1=np.round(y_b[vr]*3.667,decimals=decim); a2=np.round(y_p[vr]*3.667,decimals=decim); a3=np.round(y_d[vr]*3.667,decimals=decim)
	txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	#ax.text(-0.12,0.23,txt,size=fs_flux,ha='center')

	# ForestSector_Domestic
	ax.add_patch(Rectangle([0.01,0.01],bx_ForestSector_Domestic_w,bx_lower_h,fc=[0.85,0.9,0.85],ec=bx_ec))
	ax.text(0.25,0.04,'Forest Sector (Biogenic Carbon)',size=bx_fs,ha='center',color=[0.2,0.5,0],fontweight='bold')

	# Forest land
	ax.add_patch(Rectangle([0.02,0.1],bx_ForestSector_Domestic_w*0.51,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
	ax.text(0.15,0.29,'Forest Land',size=bx_fs,ha='center',color=[0.2,0.5,0])
	ax.text(0.15,0.26,'Change in storage (tCO$_2$e):',size=fs_flux,ha='center')
	vr='C_Forest'
	a1=np.round(y_b[vr]*3.667,decimals=decim); a2=np.round(y_p[vr]*3.667,decimals=decim); a3=np.round(y_d[vr]*3.667,decimals=decim)
	txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.15,0.23,txt,size=fs_flux,ha='center')

	# Harvested wood products
	ax.add_patch(Rectangle([0.36,0.1],0.12,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
	ax.text(0.42,0.27,'Harvested\nWood\nProducts',size=bx_fs,ha='center',color=[0.2,0.5,0])
	vr='C_HWP'
	a1=np.round(y_b[vr]*3.667,decimals=decim); a2=np.round(y_p[vr]*3.667,decimals=decim); a3=np.round(y_d[vr]*3.667,decimals=decim)
	txt='Change in\nstorage (tCO$_2$e):\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.42,0.175,txt,size=fs_flux,ha='center')

	# Lithosphere
	ax.add_patch( Rectangle([bx_ForestSector_Domestic_w+0.02,0.01],bx_esc_w*3+0.03+0.01,bx_lower_h,
						 fc=pltc.to_rgb('#b09f80'),ec=bx_ec))
	ax.text(0.75,0.04,'External Sectors (Fossil Fuels & Limestone)',size=bx_fs,ha='center',
		 color=pltc.to_rgb('#4f473a'),fontweight='bold')

	# Energy - Stationary Combustion
	ax.add_patch(Rectangle([bx_ForestSector_Domestic_w+0.02+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=pltc.to_rgb('#faeed9'),ec=bx_ec))
	ax.text(0.585,0.28,'Energy\nStationary\nCombustion',size=bx_fs,ha='center',color=pltc.to_rgb('#4f473a'))
	a1=-1*np.round(y_b['E_Substitution_EnergySC']+y_b['E_Domestic_EnergySC_ForestOperations'],decimals=decim)
	a2=-1*np.round(y_p['E_Substitution_EnergySC']+y_p['E_Domestic_EnergySC_ForestOperations'],decimals=decim)
	a3=-1*np.round(y_d['E_Substitution_EnergySC']+y_d['E_Domestic_EnergySC_ForestOperations'],decimals=decim)
	txt='Change in\nstorage (tCO2e):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.585,0.20,txt,size=fs_flux,ha='center',color='k')

	# Energy - Transportation
	ax.add_patch(Rectangle([bx_ForestSector_Domestic_w+0.02+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=pltc.to_rgb('#faeed9'),ec=bx_ec))
	ax.text(0.745,0.29,'Energy\nTransportation',size=bx_fs,ha='center',color=pltc.to_rgb('#4f473a'))
	a1=-1*np.round(y_b['E_Substitution_EnergyT']+y_b['E_Domestic_EnergyT_ForestOperations'],decimals=decim)
	a2=-1*np.round(y_p['E_Substitution_EnergyT']+y_p['E_Domestic_EnergyT_ForestOperations'],decimals=decim)
	a3=-1*np.round(y_d['E_Substitution_EnergyT']+y_d['E_Domestic_EnergyT_ForestOperations'],decimals=decim)
	txt='Change in\n storage (tCO2e):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.745,0.21,txt,size=fs_flux,ha='center',color='k')

	# IPPU
	ax.add_patch(Rectangle([bx_ForestSector_Domestic_w+0.02+0.01+bx_esc_w+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=pltc.to_rgb('#faeed9'),ec=bx_ec))
	ax.text(0.905,0.27,'Industrial\nProcesses &\nProduct Use',size=bx_fs,ha='center',color=pltc.to_rgb('#4f473a'))
	a1=-1*np.round(y_b['E_Substitution_IPPU']+y_b['E_Domestic_IPPU_ForestOperations'],decimals=decim)
	a2=-1*np.round(y_p['E_Substitution_IPPU']+y_p['E_Domestic_IPPU_ForestOperations'],decimals=decim)
	a3=-1*np.round(y_d['E_Substitution_IPPU']+y_d['E_Domestic_IPPU_ForestOperations'],decimals=decim)
	txt='Change in\n storage (tCO2e):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.905,0.185,txt,size=fs_flux,ha='center',color='k')

	#--------------------------------------------------------------------------
	# Fluxes
	#--------------------------------------------------------------------------

	# NEE
	vr='E_Domestic_ForestSector_NEE'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Net ecosystem\nexchange\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.025,0.77,txt,ha='left',size=fs_flux,color=cl_f1)
	ax.arrow(0.02,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	# Wildfire
	vr='E_Domestic_ForestSector_Wildfire'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Wildfire\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
	ax.text(0.125,0.52,txt,ha='right',size=fs_flux,color=cl_f2)
	ax.arrow(0.13,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Open burning
	vr='E_Domestic_ForestSector_OpenBurning'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Open burning\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.155,0.52,txt,ha='left',size=fs_flux,color=cl_f2)
	ax.arrow(0.15,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Denitrification
	vr='E_Domestic_ForestSector_Denit'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Denitrification\nN$_2$O emissions\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.245,0.71,txt,ha='right',size=fs_flux,color=cl_f2)
	ax.arrow(0.25,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	if pNam=='Demo_FNM':
		# Volatilization (hardwired)
		a1=np.round(0,decimals=1); a2=np.round(0.3,decimals=1); a3=np.round(0.3,decimals=1)
		txt='Volatilization of CO$_2$\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
		ax.text(0.265,0.78,txt,ha='left',size=fs_flux,color=cl_f2,bbox=dict(boxstyle='square,pad=-0.07', fc='w', ec='none'))
		ax.arrow(0.26,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)
	
		# Volatilization
		vr='E_Domestic_ForestSector_Volat'
		a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
		txt='Volatilization /\ndeposition of NH$_3$\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')*'
		ax.text(0.30,0.63,txt,ha='left',size=fs_flux,color=cl_f1)
		ax.arrow(0.295,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	# HWP fluxes
	vr='E_Domestic_ForestSector_HWP'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Product decay\nand combustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.405,0.58,txt,ha='right',va='top',size=fs_flux,color=cl_f2)
	ax.arrow(0.41,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Removals
	vr='C_ToMillTotal'
	mult=3.667
	a1=int(mult*y_b[vr]); a2=int(mult*y_p[vr]); a3=int(mult*y_d[vr])
	txt='Removals\n(tCO$_2$e)\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.3125,0.30,txt,ha='center',size=fs_flux)
	ax.arrow(0.265,0.275,0.09,0,head_width=0.01,head_length=arrow_head_w,fc='k',ec='k',lw=arrow_lw)

	# Bioenergy combustion
	vr='E_Bioenergy'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Bioenergy\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.515-0.08,0.58,txt,ha='left',va='top',size=fs_flux,color=cl_f2)
	ax.arrow(0.51-0.08,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# ESC operational emissions
	#vr='E_ForestOperations'
	vr='E_Domestic_EnergySC_ForestOperations'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel\nuse\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.53,0.64,txt,ha='left',va='top',size=fs_flux,color=cl_f2)
	ax.arrow(0.525,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Substitution energy
	vr='E_Substitution_EnergySC'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Substitutions\nfrom stationary\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.645,0.82,txt,ha='right',va='top',size=fs_flux,color=cl_f1)
	ax.arrow(0.65,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	# Transportation
	vr='E_Domestic_EnergyT_ForestOperations'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel\nuse\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.685,0.64,txt,ha='left',va='top',size=fs_flux,color=cl_f2)
	ax.arrow(0.68,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Substitution transportation
	vr='E_Substitution_EnergyT'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Substitutions\nfrom\ntransportation\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.765,0.82,txt,ha='left',va='top',size=fs_flux,color=cl_f1)
	ax.arrow(0.76,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	# IPPU
	vr='E_Domestic_IPPU_ForestOperations'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel\nuse\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.865,0.64,txt,ha='right',va='top',size=fs_flux,color=cl_f2)
	ax.arrow(0.87,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc=cl_f2,ec=cl_f2,lw=arrow_lw)

	# Substitution IPPU
	vr='E_Substitution_IPPU'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Substitutions\nfrom IPPU\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.885,0.82,txt,ha='left',va='top',size=fs_flux,color=cl_f1)
	ax.arrow(0.88,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	if pNam=='Demo_FNM':
		# Volatilization (hardwired, CO2 sequestration during urea production)
		a1=np.round(0,decimals=1); a2=np.round(-0.3,decimals=1); a3=np.round(-0.3,decimals=1)
		txt='CO$_2$ fixation\nduring urea\nproduction\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
		ax.text(0.975,0.53,txt,ha='right',size=fs_flux,color=cl_f1)
		ax.arrow(0.98,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc=cl_f1,ec=cl_f1,lw=arrow_lw)

	#plt.tight_layout()
	ax.set(position=[0.0,0,1.04,1],visible='Off',xticks=[],yticks=[])
	ax.axis('off')

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GHG_Balance_Schematic_S' + str(iP) + '_minus_S' + str(iB) + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']),'png',900)
	return ax

#%% Plot Cashflow
def PlotCashflow(meta,pNam,mos,**kwargs):
	# Key word arguments
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iSP']
		iSS=kwargs['iSS']
		iYS=kwargs['iYS']
	else:
		iPS=0
		iSS=0
		iYS=0
	cNam=kwargs['cNam']
	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	cl=np.array([[0.17,0.34,0.69],[0.55,0.9,0],[0.5,0,1],[0,1,1]])

	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); lw=0.75
	# Cost
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label=lab1)
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label=lab2)
	ax[0,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cost (CAD x 1000)')
	ax[0,0].legend(loc='upper right',frameon=False,facecolor=None)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta cost
	ax[0,1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata']['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw)
	ax[0,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ Cost (CAD x 1000)')
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Net Revenue
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label='Baseline')
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label='Project')
	ax[1,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Net revenue (CAD x 1000)')
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta net revenue
	ax[1,1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata']['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw)
	ax[1,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ net revenue (CAD x 1000)')
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Cumulative net revenue
	ax[2,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Revenue Net_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label='Baseline')
	ax[2,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Revenue Net_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label='Project')
	ax[2,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cumulative net revenue\n (CAD x 1000)')
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta cumulative net revenue
	ax[2,1].plot(tv[iT],np.zeros(iT.size),'-k',lw=3,color=[0.85,0.85,0.85])
	ax[2,1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata']['Mean']['Revenue Net_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw,label='Undiscounted')
	ax[2,1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata']['Mean']['Revenue Net Disc_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--k',color=0.5*cl[2,:],mfc=0.5*cl[2,:],mec=0.5*cl[2,:],lw=lw,label='Discounted')
	ax[2,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ cumulative net\n revenue (CAD x 1000)')
	ax[2,1].legend(loc='lower left',frameon=False,facecolor=None)
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Cashflow_' + cNam,'png',900)
	return

#%%
def PlotPools(meta,mos,pNam,**kwargs):

	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	vs=['C_Biomass','C_DeadWood','C_Litter','C_Soil','C_InUse','C_WasteSystems','C_Geological','E_NAB_Cumulative']
	vs2=['Biomass','Dead Wood','Litter','Soil','In-use Products','Dump and Landfill','Geological','Atmosphere (50% Subs.)']

	operS=kwargs['OperSpace']

	cl=np.array([[0.27,0.44,0.79],[0.55,0.9,0],[0.5,0,1],[0,1,1]])
	symb=['-','--','-.',':','-']
	if 'ScenarioIndexList' in kwargs.keys():
		# Generate one figure for a custom set of scenarios
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		for i in range(4):
			for j in range(2):
				for iScn in range(len(kwargs['ScenarioIndexList'])):
					s=kwargs['ScenarioIndexList'][iScn]
					be=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					lo=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].plot([0,4000],[0,0],'k-',color=[0.85,0.85,0.85],lw=2)
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)

					if (i==0) & (j==0):
						if 'LegendLoc' in kwargs.keys():
							legloc=kwargs['LegendLoc']
						else:
							legloc='lower left'
						ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
						if vs[cnt]=='E_NAB_Cumulative':
							us='\n(tCO2e/ha)'
						else:
							us='\n(tC/ha)'
				ax[i,j].set(ylabel=vs2[cnt] + us,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				ax[i,j].yaxis.set_ticks_position('both');
				ax[i,j].xaxis.set_ticks_position('both');
				ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Pools_CustomScenarioList','png',900)
	elif 'cNam' in kwargs.keys():
		# Generate one figure per scenario comparison
		cNam=kwargs['cNam']
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		sL=[mos[pNam]['Delta'][cNam]['iB'],mos[pNam]['Delta'][cNam]['iP']]
		for i in range(4):
			for j in range(2):
				for iScn in range(len(sL)):
					be=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					if iScn==0:
						be0=be
					else:
						be1=be
					lo=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].plot([0,4000],[0,0],'k-',color=[0.85,0.85,0.85],lw=2)
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])

				if 'TextDelta' in kwargs.keys():
					if kwargs['TextDelta']['Units']=='Actual':
						be_d=be1-be0
						uni=''
					elif kwargs['TextDelta']['Units']=='Relative':
						be_d=(be1-be0)/be0*100
						uni='%'
					else:
						print('Text delta units unrecgonized')
					for yr in kwargs['TextDelta']['Year']:
						iT2=np.where(tv[iT]==yr)[0]
						if iT2.size==0:
							continue
						if (be_d[iT2]>0) & (be0[iT2]<0) & (be1[iT2]<0):
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						elif (be_d[iT2]>0):
							tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						else:
							tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						ax[i,j].text(tv[iT][iT2],mlti*be1[iT2],tx,fontsize=6,fontweight='bold',ha='center')

				if (i==0) & (j==0):
					if 'LegendLoc' in kwargs.keys():
						legloc=kwargs['LegendLoc']
					else:
						legloc='lower left'
					ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				if vs[cnt]=='E_NAB_Cumulative':
					us='\n(tCO2e/ha)'
				else:
					us='\n(tC/ha)'
				ax[i,j].set(ylabel=vs2[cnt] + us,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				cnt=cnt+1
		plt.tight_layout()
		gu.axletters(ax,plt,0.02,0.82,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Pools_' + cNam,'png',900)

	return

#%%
def PlotFluxes(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	vs=['E_Domestic_ForestSector_NPP','C_G_Net_Reg','E_Domestic_ForestSector_RH','E_Domestic_ForestSector_OpenBurning','E_Domestic_ForestSector_Wildfire','E_Domestic_ForestSector_HWP','E_Substitution_Total','E_NAB']
	vs2=['NPP','Net growth','RH','Open burning','Wildfire','HWP','Substitutions','Net emissions']

	operS=kwargs['OperSpace']

	cl=np.array([[0.27,0.44,0.79],[0.55,0.9,0],[0.5,0,1],[0,1,1]])
	symb=['-','--','-.',':','-']
	if 'ScenarioIndexList' in kwargs.keys():
		# Generate one figure for a custom set of scenarios
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		for i in range(4):
			for j in range(2):
				for iScn in range(len(kwargs['ScenarioIndexList'])):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Domestic_ForestSector_NPP':
						adj=-1.0
					s=kwargs['ScenarioIndexList'][iScn]
					be=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					lo=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					if (i==0) & (j==0):
						if 'LegendLoc' in kwargs.keys():
							legloc=kwargs['LegendLoc']
						else:
							legloc='lower left'
						ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				#if vs[cnt]=='E_Domestic_ForestSector_NPP':
				#	ylab='-1 x ' + vs2[cnt] + '\n(tCO2e/ha/yr)'
				#else:
				ylab=vs2[cnt] + '\n(tCO2e/ha/yr)'
				ax[i,j].set(ylabel=ylab + ' (tCO2e/ha/yr)',xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_CustomScenarioList','png',900)
	elif 'cNam' in kwargs.keys():
		# Generate one figure per scenario comparison
		k=kwargs['cNam']
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		sL=[mos[pNam]['Delta'][k]['iB'],mos[pNam]['Delta'][k]['iP']]
		for i in range(4):
			for j in range(2):
				ymx=0
				for iScn in range(len(sL)):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Domestic_ForestSector_NPP':
						adj=-1.0
					be=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					if iScn==0:
						be0=be
					else:
						be1=be
					ymx=np.maximum(ymx,np.max(be))
					lo=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
				if (i==0) & (j==0):
					if 'LegendLoc' in kwargs.keys():
						legloc=kwargs['LegendLoc']
					else:
						legloc='lower left'
					ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')

				if 'TextDelta' in kwargs.keys():
					if kwargs['TextDelta']['Units']=='Actual':
						be_d=be1-be0
						uni=''
					elif kwargs['TextDelta']['Units']=='Relative':
						be_d=(be1-be0)/be0*100
						uni='%'
					else:
						print('Text delta units unrecgonized')
					for yr in kwargs['TextDelta']['Year']:
						iT2=np.where(tv[iT]==yr)[0]
						if iT2.size==0:
							continue
						if (be_d[iT2]>0) & (be0[iT2]<0) & (be1[iT2]<0):
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						elif (be_d[iT2]>0):
							tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						else:
							tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						ax[i,j].text(tv[iT][iT2],mlti*be1[iT2],tx,fontsize=6,fontweight='bold',ha='center')

				ylab=vs2[cnt] + '\n(tCO2e/ha/yr)'
				if np.isin(vs[cnt],['E_Substitution_Total','E_NAB'])==True:
					plt.plot([0,5000],[0,0],'k-',lw=0.25)
					ax[i,j].set(ylabel=ylab,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				else:
					ax[i,j].set(ylabel=ylab,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])],ylim=[0,1.1*ymx])

				cnt=cnt+1
		gu.axletters(ax,plt,0.02,0.82,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_' + k,'png',900)

	return

#%%
def PlotAge(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	vs=['A','A','A','A','A','A','A','A']
	vs2=['Age','Age','Age','Age','','','','']

	operS=kwargs['OperSpace']

	cl=np.array([[0.27,0.44,0.79],[0.55,0.9,0],[0.5,0,1],[0,1,1]])
	symb=['-','--','-.',':','-']
	if 'ScenarioIndexList' in kwargs.keys():
		# Generate one figure for a custom set of scenarios
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		for i in range(4):
			for j in range(2):
				for iScn in range(len(kwargs['ScenarioIndexList'])):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Domestic_ForestSector_NPP':
						adj=-1.0
					s=kwargs['ScenarioIndexList'][iScn]
					be=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					lo=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=adj*mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					if (i==0) & (j==0):
						if 'LegendLoc' in kwargs.keys():
							legloc=kwargs['LegendLoc']
						else:
							legloc='lower left'
						ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				#if vs[cnt]=='E_Domestic_ForestSector_NPP':
				#	ylab='-1 x ' + vs2[cnt] + '\n(tCO2e/ha/yr)'
				#else:
				ylab=vs2[cnt] + '\n(tCO2e/ha/yr)'
				ax[i,j].set(ylabel=ylab + ' (tCO2e/ha/yr)',xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_CustomScenarioList','png',900)
	elif 'cNam' in kwargs.keys():
		# Generate one figure per scenario comparison
		k=kwargs['cNam']
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		sL=[mos[pNam]['Delta'][k]['iB'],mos[pNam]['Delta'][k]['iP']]
		for i in range(4):
			for j in range(2):
				ymx=0
				for iScn in range(len(sL)):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Domestic_ForestSector_NPP':
						adj=-1.0
					be=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					if iScn==0:
						be0=be
					else:
						be1=be
					ymx=np.maximum(ymx,np.max(be))
					lo=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
				if (i==0) & (j==0):
					if 'LegendLoc' in kwargs.keys():
						legloc=kwargs['LegendLoc']
					else:
						legloc='lower left'
					ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')

				if 'TextDelta' in kwargs.keys():
					if kwargs['TextDelta']['Units']=='Actual':
						be_d=be1-be0
						uni=''
					elif kwargs['TextDelta']['Units']=='Relative':
						be_d=(be1-be0)/be0*100
						uni='%'
					else:
						print('Text delta units unrecgonized')
					for yr in kwargs['TextDelta']['Year']:
						iT2=np.where(tv[iT]==yr)[0]
						if iT2.size==0:
							continue
						if (be_d[iT2]>0) & (be0[iT2]<0) & (be1[iT2]<0):
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						elif (be_d[iT2]>0):
							tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						else:
							tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						ax[i,j].text(tv[iT][iT2],mlti*be1[iT2],tx,fontsize=6,fontweight='bold',ha='center')

				ylab='Age (years)'
				ax[i,j].set(ylabel=ylab,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])],ylim=[0,1.1*ymx])

				cnt=cnt+1
		gu.axletters(ax,plt,0.02,0.82,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Age_' + k,'png',900)

	return

#%%
def PlotFluxesBiomass(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	vs=['C_G_Gross','C_M_Reg','C_G_Net_Reg','C_G_Net']
	vs2=['Gross growth','Regular mortality','Net growth (regular)','Net growth (total)']

	operS=kwargs['OperSpace']

	cl=np.array([[0.27,0.44,0.79],[0.55,0.9,0],[0.5,0,1],[0,1,1]])

	symb=['-','--','-.',':','-']
	if 'ScenarioIndexList' in kwargs.keys():
		# Generate one figure for a custom set of scenarios
		cnt=0
		fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		for i in range(2):
			for j in range(2):
				ymx=0
				for iScn in range(len(kwargs['ScenarioIndexList'])):
					s=kwargs['ScenarioIndexList'][iScn]
					be=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					ymx=np.maximum(ymx,np.max(be))
					lo=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=mos[pNam]['Scenarios'][s][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					if (i==0) & (j==0):
						if 'LegendLoc' in kwargs.keys():
							legloc=kwargs['LegendLoc']
						else:
							legloc='lower left'
						ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				#if vs[cnt]=='E_Domestic_ForestSector_NPP':
				#	ylab='-1 x ' + vs2[cnt] + '\n(tC/ha/yr)'
				#else:
				ylab=vs2[cnt] + '\n(tC/ha/yr)'
				ax[i,j].set(ylabel=ylab + ' (tC/ha/yr)',xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])],ylim=[0,1.1*ymx])
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_CustomScenarioList','png',900)
	elif 'cNam' in kwargs.keys():
		# Generate one figure per scenario comparison
		k=kwargs['cNam']
		cnt=0
		fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		sL=[mos[pNam]['Delta'][k]['iB'],mos[pNam]['Delta'][k]['iP']]
		for i in range(2):
			for j in range(2):
				ymx=0
				for iScn in range(len(sL)):
					be=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
					ymx=np.maximum(ymx,np.max(be))
					if iScn==0:
						be0=be
					else:
						be1=be
					lo=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P025'][iT,iPS,iSS,iYS]
					hi=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P975'][iT,iPS,iSS,iYS]
					lo2=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P250'][iT,iPS,iSS,iYS]
					hi2=mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble P750'][iT,iPS,iSS,iYS]
					ax[i,j].fill_between(tv[iT],lo,hi,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					ax[i,j].fill_between(tv[iT],lo2,hi2,color=cl[iScn,:],alpha=Alpha,linewidth=0)
					if 'ScenarioLabels' in kwargs.keys():
						lab=kwargs['ScenarioLabels'][iScn]
					else:
						lab='Scenario ' + str(iScn+1)
					ax[i,j].plot(tv[iT],be,symb[iScn],color=cl[iScn,:],lw=1,label=lab)
					ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])

				if 'TextDelta' in kwargs.keys():
					if kwargs['TextDelta']['Units']=='Actual':
						be_d=be1-be0
						uni=''
					elif kwargs['TextDelta']['Units']=='Relative':
						be_d=(be1-be0)/be0*100
						uni='%'
					else:
						print('Text delta units unrecgonized')
					for yr in kwargs['TextDelta']['Year']:
						iT2=np.where(tv[iT]==yr)[0]
						if iT2.size==0:
							continue
						if (be_d[iT2]>0) & (be0[iT2]<0) & (be1[iT2]<0):
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						elif (be_d[iT2]>0):
							tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						else:
							tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						ax[i,j].text(tv[iT][iT2],mlti*be1[iT2],tx,fontsize=6,fontweight='bold',ha='center')

				if (i==0) & (j==0):
					if 'LegendLoc' in kwargs.keys():
						legloc=kwargs['LegendLoc']
					else:
						legloc='lower left'
					ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				ylab=vs2[cnt] + '\n(tC/ha/yr)'
				ax[i,j].set(ylabel=ylab,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])],ylim=[0,1.1*ymx])
				cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\FluxesBiomass_' + k,'png',900)
	return

#%% Plot NEP, RH, and NPP
def PlotDeltaNEE(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	cNam=kwargs['cNam']

	operS=kwargs['OperSpace']

	cl=np.array([[0.75,0,0],[0,0.85,0],[0.29,0.49,0.77]])
	ymin=0.0;ymax=0.0

	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
	for j in range(0,2):
		ax[j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')
	vn='C_RH'
	lo=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=cl[0,:],label='RH')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	vn='C_NPP'
	lo=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'-.',color=cl[1,:],label='NPP')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	vn='E_Domestic_ForestSector_NEE'
	lo=-mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=-mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=-mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=-mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],-mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=cl[2,:],label='NEP')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	ax[0].legend(loc="lower right",frameon=False,facecolor=None,edgecolor='w')
	ax[0].set(ylabel='Annual $\Delta$ (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	ymin=0.0
	ymax=0.0

	vn='C_RH'
	lo=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'--',color=cl[0,:])
	ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]]);
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))

	vn='C_NPP'
	lo=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=np.cumsum(3.667*mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'-.',color=cl[1,:])
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))

	vn='E_Domestic_ForestSector_NEE'
	lo=-np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=-np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=-np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=-np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=-np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'-',color=cl[2,:])
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))
	ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.035,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\NEE_Balance_' + cNam,'png',900)

	return

#%%
def PlotDeltaGHGB(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1='Baseline'
		lab2='Action'

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	if 'LegendLoc' in kwargs.keys():
		ll=kwargs['LegendLoc']
	else:
		ll='upper right'

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	cNam=kwargs['cNam']

	operS=kwargs['OperSpace']

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']

	clD=np.array([0.6,0.1,1])

	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha1=0.06; Alpha2=0.08; Alpha3=0.11
	for i in range(0,2):
		for j in range(0,2):
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both')
			ax[i,j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB][operS]['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=(0,0.5,1),label=lab1)
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iP][operS]['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=(0,0.6,0),label=lab2)
	ax[0,0].legend(loc=ll,frameon=False,facecolor=None,edgecolor='w')
	ax[0,0].set(ylabel='Emissions\n(tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iB][operS]['E_NAB_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=(0,0.5,1))
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iP][operS]['E_NAB_Cumulative']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=(0,0.6,0))
	ax[0,1].set(ylabel='Cumulative emissions\n(tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	mu=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS]
	lo=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P750'][iT,iPS,iSS,iYS]
	sig=1*mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble SD'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=clD,alpha=Alpha1,linewidth=0,label='95 C.I.')
	ax[1,0].fill_between(tv[iT],lo2,hi2,color=clD,alpha=Alpha2,linewidth=0,label='50 C.I.')
	ax[1,0].fill_between(tv[iT],mu-sig,mu+sig,color=0.5*clD,alpha=Alpha3,linewidth=0,label='S.D.')
	ax[1,0].plot(tv[iT],mu,'-',color=clD,label='NAB')

	ax[1,0].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NSB']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=0.25*clD,label='NSB')
	ax[1,0].legend(loc=ll,frameon=0)
	ax[1,0].set(ylabel='$\Delta$ emissions\n(tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);
	if 'yLimPad' in kwargs.keys():
		ax[1,0].set(ylim=[kwargs['yLimPad']*np.min(mu),kwargs['yLimPad']*np.max(mu)])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	mu=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS]
	lo=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P750'][iT,iPS,iSS,iYS]
	sig=1*mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble SD'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=clD,alpha=Alpha1,linewidth=0)
	ax[1,1].fill_between(tv[iT],lo2,hi2,color=clD,alpha=Alpha2,linewidth=0)
	ax[1,1].fill_between(tv[iT],mu-sig,mu+sig,color=clD,alpha=Alpha3,linewidth=0)
	ax[1,1].plot(tv[iT],mu,'-',color=clD,label='NAB')
	ax[1,1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NSB_Cumulative_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=0.25*clD,label='NSB')
	ax[1,1].legend(loc=ll,frameon=0)
	ax[1,1].set(ylabel='Cumulative $\Delta$ emissions\n(tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	if 'yLimPad' in kwargs.keys():
		ax[1,1].set(ylim=[kwargs['yLimPad']*np.min(mu),kwargs['yLimPad']*np.max(mu)])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	plt.tight_layout()
	gu.axletters(ax,plt,0.035,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GHG_Balance_' + cNam,'png',900)

	return

#%%
def PlotGHGBenefit(meta,pNam,mos,tv,iT,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	operS=kwargs['OperSpace']

	cNam=kwargs['cNam']

	cl=np.array([[0,0.5,1],[0,0.6,0],[0,1,1],[0.5,0,1]])
	cnt=1
	symb=['-','--','-.',':','-']
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(16,6)); Alpha=0.09
	for i in range(0,2):
		ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both')
		ax[i].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

	lo=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[0].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1) )
	ax[0].legend(loc="upper right",frameon=False,facecolor=None,edgecolor='w')
	ax[0].set(ylabel='$\Delta$GHG (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);
	lo=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mos[pNam]['Delta'][cNam]['ByStrata'][operS]['E_NAB_Cumulative_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1))
	ax[1].set(ylabel='Cumulative $\Delta$GHG (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	gu.axletters(ax,plt,0.035,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['gp']['AxesLetterStyle'],FontWeight=meta['Graphics']['gp']['AxesLetterFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GHG_Benefit_' + cNam,'png',900)

	return

#%% Summary bioenergy description
def SummaryBioenergy(meta,mos,sc,th):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['Project']['Year Project']) & (tv<=meta[pNam]['Project']['Year Project']+th) )[0]

	d={}
	d['EI PelletExport (tCO2e/ODT)']=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI PelletExport Boiler (tCO2e/GJ)']=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['ByStrata']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI PelletExport Boiler+Ops (tCO2e/GJ)']=np.sum( (mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['ByStrata']['E_EnergyT_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['ByStrata']['E_IPPU_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['ByStrata']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,0,0]) )/np.sum(mos[pNam]['Delta'][sc]['ByStrata']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI Pellet Manufacture (tCO2e/ODT Pellets)']=np.sum( (mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_ForestOperations']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['ByStrata']['E_IPPU_ForestOperations']['Ensemble Mean'][iT,0,0]) )/np.sum(mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI OperationForestry (tCO2e/ODT)']= \
			np.sum( (mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_ForestOperations']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['ByStrata']['E_EnergyT_ForestOperations']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['ByStrata']['E_IPPU_ForestOperations']['Ensemble Mean'][iT,0,0]) ) / \
			np.sum( (mos[pNam]['Delta'][sc]['ByStrata']['ODT Lumber']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT LogExport']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT Plywood']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT OSB']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT MDF']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT Paper']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletExport']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletDomGrid']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletDomRNG']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT PowerFacilityDom']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['ByStrata']['ODT PowerGrid']['Ensemble Mean'][iT,0,0]) )

	d['Energy efficiency (GJ/ODT)']=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['GJ PelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['ByStrata']['ODT PelletExport']['Ensemble Mean'][iT,0,0])

	# Displacement factor

	SubTot=-np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_Substitution_Total']['Ensemble Mean'][iT,0,0])
	SubE=-np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_Substitution_Energy']['Ensemble Mean'][iT,0,0])
	SubM=-np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_Substitution_Material']['Ensemble Mean'][iT,0,0])
	Wood=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['C_ToLumber']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['ByStrata']['C_ToPlywood']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['ByStrata']['C_ToMDF']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['ByStrata']['C_ToOSB']['Ensemble Mean'][iT,0,0])
	Bioenergy=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_Bioenergy']['Ensemble Mean'][iT,0,0])
	BioenergyPlusOps=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_EnergyT_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_IPPU_ForestOperations']['Ensemble Mean'][iT,0,0])
	BioenergyOpsAndHWPDecay=np.sum(mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_EnergyT_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_IPPU_ForestOperations']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['ByStrata']['E_EnergySC_Bioenergy']['Ensemble Mean'][iT,0,0])

	d['Displacement Factor Total (tC/tC)']=SubTot/BioenergyOpsAndHWPDecay
	d['Displacement Factor Energy (tC/tC)']=SubE/BioenergyPlusOps
	d['Displacement Factor Materials (tC/tC)']=(SubM/3.667)/Wood

	# Save
	df=pd.DataFrame.from_dict(d,orient='index')
	df.to_excel(meta['Paths']['Project'] + '\\Outputs\\SummaryBioenergy.xlsx')

	return d

#%%
def NA_CompareSpecifications_ChangeInPools(meta,pNam,mos):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=2000) & (tv<=2100))[0]
	iPS=0; iSS=0; iYS=0
	cl=np.array([[0.6,0.95,0],[0,0.75,0],[0.95,0.7,0.7],[0.55,0.3,0.3],[0.25,0.1,0.1],[0.27,0.49,0.76]])
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,10))
	ax[0,0].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	sc='NGS'
	AGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[0,0].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[0,0].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[0,0].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[0,0].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[0,0].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[0,0].legend(loc="upper right",frameon=False,facecolor=None)
	ax[0,0].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,35])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT'
	AGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[0,1].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[0,1].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[0,1].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[0,1].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[0,1].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[0,1].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[0,1].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,35])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT+T'
	AGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,0].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[1,0].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[1,0].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[1,0].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[1,0].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[1,0].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[1,0].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,35])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT+T+D'
	AGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['ByStrata']['Mean']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,1].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[1,1].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[1,1].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[1,1].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[1,1].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[1,1].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[1,1].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (tC ha$^{-1}$)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,35])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.02,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=['NGS','NGT','NGT+T','NGT+T+D'],LabelSpacer=0.021)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\CompareSpecifications_ChangeInPools','png',900);
	return

#%%
def PlotVolumeMerchLive(meta,pNam,mos,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	if 'LegendLoc' in kwargs.keys():
		LegendLoc=kwargs['LegendLoc']
	else:
		LegendLoc='lower right'

	cNam=kwargs['cNam']

	operS=kwargs['OperSpace']

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	v='V_MerchLive'
	y1=mos[pNam]['Scenarios'][iB][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	y2=mos[pNam]['Scenarios'][iP][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ymx=1.2*np.max([np.max(y1),np.max(y2)])

	vb='V_WholeStemLive'
	y1b=mos[pNam]['Scenarios'][iB][operS][vb]['Ensemble Mean'][iT,iPS,iSS,iYS]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]))
	if 'FillDelta' in kwargs.keys():
		ax.fill_between(tv[iT],y1,y2,color=[0.8,1,0.3],alpha=0.5,linewidth=0)
	ax.plot(tv[iT],y1,'-',color=[0,0,0],lw=1.25,label=lab1)
	if 'IncludeWholeStem' in kwargs.keys():
		ax.plot(tv[iT],y1b,':',color=[0,0,0],lw=1.25,label='Whole stem live')
	ax.plot(tv[iT],y2,'--',color=[0.6,0.9,0],lw=1.25,label=lab2)
	ax.legend(loc=LegendLoc,frameon=False,facecolor=None)

	if 'TextDelta' in kwargs.keys():
		be_d=(y2-y1)/y1*100
		for yr in kwargs['TextDelta']['Year']:
			iT2=np.where(tv[iT]==yr)[0]
			if iT2.size==0:
				continue
			if be_d[iT2]>0:
				tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			elif np.isnan(be_d[iT2])==True:
				tx=''
			else:
				tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			ax.text(tv[iT][iT2],1.05*y2[iT2],tx,fontsize=6,fontweight='bold',ha='center',color=[0.3,0.45,0])

	ax.set(ylabel=r'Live net merch. stemwood volume (m$^3$ ha$^-$$^1$)',xlabel=r'Time, years',ylim=[0,ymx],xlim=[np.min(tv[iT]),np.max(tv[iT])]);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

	if 'HarvestYears' in kwargs.keys():
		for yr in kwargs['HarvestYears']:
			ind=np.where(tv==yr)[0]
			ax.annotate('Harvest',xy=(yr,y1[ind]),xytext=(yr,y1[ind]+200),arrowprops={'color':'black','arrowstyle':'->'},ha='center');

	flg=0
	if flg==1:
		ax.annotate('Nutrient\napplication',xy=(meta[pNam]['Project']['Year Project'],350),xytext=(meta[pNam]['Project']['Year Project'],650),
			 arrowprops={'color':'black','arrowstyle':'->'},ha='center');

	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\VolumeMerchLive_' + cNam,'png',900)
	return

#%%
def PlotVolumeMerchTotal(meta,pNam,mos,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	cNam=kwargs['cNam']

	operS=kwargs['OperSpace']

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	v='V_MerchTotal'
	y1=mos[pNam]['Scenarios'][iB][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	y2=mos[pNam]['Scenarios'][iP][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ymx=1.2*np.max([np.max(y1),np.max(y2)])

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]))
	if 'FillDelta' in kwargs.keys():
		ax.fill_between(tv[iT],y1,y2,color=[0.8,1,0.3],alpha=0.5,linewidth=0)
	ax.plot(tv[iT],y1,'-',color=[0,0,0],lw=1.25,label=lab1)
	ax.plot(tv[iT],y2,'--',color=[0.6,0.9,0],lw=1.25,label=lab2)
	ax.legend(loc='lower right',frameon=False,facecolor=None)

	if 'TextDelta' in kwargs.keys():
		be_d=(y2-y1)/y1*100
		for yr in kwargs['TextDelta']['Year']:
			iT2=np.where(tv[iT]==yr)[0]
			if iT2.size==0:
				continue
			if be_d[iT2]>0:
				tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			elif np.isnan(be_d[iT2])==True:
				tx=''
			else:
				tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			ax.text(tv[iT][iT2],1.05*y2[iT2],tx,fontsize=6,fontweight='bold',ha='center',color=[0.3,0.45,0])

	ax.set(ylabel=r'Net merch. stemwood volume (m$^3$ ha$^-$$^1$)',xlabel=r'Time, years',ylim=[0,ymx],xlim=[np.min(tv[iT]),np.max(tv[iT])]);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.annotate('Harvest',xy=(meta[pNam]['Project']['Year Project']-36,1000),xytext=(meta[pNam]['Project']['Year Project']-36,1200),
			 arrowprops={'color':'black','arrowstyle':'->'},ha='center');
	flg=0
	if flg==1:
		ax.annotate('Nutrient\napplication',xy=(meta[pNam]['Project']['Year Project'],350),xytext=(meta[pNam]['Project']['Year Project'],650),
			 arrowprops={'color':'black','arrowstyle':'->'},ha='center');
	#gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\VolumeMerchTotal_' + cNam,'png',900)
	return

#%%
def PlotNEP(meta,pNam,mos,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	if 'LegendLoc' in kwargs.keys():
		LegendLoc=kwargs['LegendLoc']
	else:
		LegendLoc='lower right'

	cNam=kwargs['cNam']

	operS=kwargs['OperSpace']

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	v='E_Domestic_ForestSector_NEE'
	y1=-1*mos[pNam]['Scenarios'][iB][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/3.667
	y2=-1*mos[pNam]['Scenarios'][iP][operS][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/3.667
	ymn=np.min([np.min(y1),np.min(y2)])-1
	ymx=np.max([np.max(y1),np.max(y2)])+1

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(FigSize[0],FigSize[1]))
	ax.plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=2)
	if 'FillDelta' in kwargs.keys():
		ax.fill_between(tv[iT],y1,y2,color=[0.8,1,0.3],alpha=0.5,linewidth=0)

	ax.plot(tv[iT],y1,'-',color=[0,0,0],lw=1.25,label=lab1)
	ax.plot(tv[iT],y2,'--',color=[0.6,0.9,0],lw=1.25,label=lab2)
	ax.legend(loc=LegendLoc,frameon=False,facecolor=None)

	if 'TextDelta' in kwargs.keys():
		be_d=(y2-y1)/y1*100
		for yr in kwargs['TextDelta']['Year']:
			iT2=np.where(tv[iT]==yr)[0]
			if iT2.size==0:
				continue
			if be_d[iT2]>0:
				tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			elif np.isnan(be_d[iT2])==True:
				tx=''
			else:
				tx='' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			ax.text(tv[iT][iT2],1.05*y2[iT2],tx,fontsize=6,fontweight='bold',ha='center',color=[0.3,0.45,0])

	ax.set(ylabel=r'Net ecosystem productivity (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel=r'Time, years',ylim=[ymn,ymx],xlim=[np.min(tv[iT]),np.max(tv[iT])])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.annotate('Harvest',xy=(meta[pNam]['Project']['Year Project']-36,1000),xytext=(meta[pNam]['Project']['Year Project']-36,1200),
			 arrowprops={'color':'black','arrowstyle':'->'},ha='center');
	flg=0
	if flg==1:
		ax.annotate('Nutrient\napplication',xy=(meta[pNam]['Project']['Year Project'],350),xytext=(meta[pNam]['Project']['Year Project'],650),
			 arrowprops={'color':'black','arrowstyle':'->'},ha='center');
	#gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\NEP_' + cNam,'png',900)
	return

#%%
def MitigationSummary(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	r_discE=0.01
	r_discC=0.03
	t_disc=tv-meta[pNam]['Project']['Year Project']
	multip=1000
	cNam=kwargs['cNam']
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	iPS=0
	iSS=0
	iYS=0

	vL=['E_NSB','E_NAB','Revenue Net','Revenue Gross','Cost Total']
	d={}
	for cNam in mos[pNam]['Delta'].keys():
		d[cNam]={}
		for v in vL:
			if np.isin(v,['E_NSB','E_NAB'])==True:
				r_disc=r_discE
			else:
				r_disc=r_discC
			d[cNam][v]={}
			d[cNam][v]['ann']=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
			d[cNam][v]['ann disc']=mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/((1+r_disc)**t_disc[iT])
			d[cNam][v]['cumu']=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/multip)
			d[cNam][v]['cumu disc']=np.cumsum(mos[pNam]['Delta'][cNam]['ByStrata']['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/multip/((1+r_disc)**t_disc[iT]))

	v='ann'
	v='cumu'
	plt.close('all');
	plt.plot(tv[iT],-d[cNam]['Cost Total'][v]/1000,'r-')
	plt.plot(tv[iT],d[cNam]['Revenue Gross'][v]/1000,'b-')
	plt.plot(tv[iT],d[cNam]['Revenue Net'][v]/1000,'g-')

	plt.close('all');
	plt.plot(-d[cNam]['Cost Total']['cumu'],-d[cNam]['E_NAB']['cumu'],'b-')
	plt.plot(d[cNam]['Revenue Net']['cumu'],-d[cNam]['E_NAB']['cumu'],'g-')
	plt.plot(-d[cNam]['Cost Total']['cumu disc'],-d[cNam]['E_NAB']['cumu disc'],'c--')
	plt.plot(d[cNam]['Revenue Net']['cumu disc'],-d[cNam]['E_NAB']['cumu disc'],'r--')

	dC=d[cNam]['Cost Total']['cumu disc'][-1]
	dNR=d[cNam]['Revenue Net']['cumu disc'][-1]
	dE=d[cNam]['E_NAB']['cumu disc'][-1]

	print(dC/-dE)
	print(-dNR/-dE)

	return

#%%
def PlotVolumeAll(meta,pNam,mos,**kwargs):

	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	operS=kwargs['OperSpace']
	cNam=kwargs['cNam']

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	cl1=[0.27,0.49,0.77]
	cl2=[0.45,0,0.95]
	cl3=[1,0.5,0]
	cl4=[1,1,0]
	ms=4; lw=1;

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	y0=np.max(mos[pNam]['Scenarios'][iB]['Mean']['V_WholeStemTotal']['Ensemble Mean'][iT,iPS,iSS,iYS])
	y1=np.max(mos[pNam]['Scenarios'][iP]['Mean']['V_WholeStemTotal']['Ensemble Mean'][iT,iPS,iSS,iYS])
	yl=np.maximum(y0,y1)+15

	plt.close('all'); fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(FigSize[0],FigSize[1]))
	#--------------------------------------------------------------------------
	# Live volume
	#--------------------------------------------------------------------------
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_WholeStemLive']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,mew=lw,label='Whole stem modelled')
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_MerchLive']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,mew=lw,label='Merchantable modelled')

	flg=1
	if flg==1:
		# Johnstone 2002
		ax[0,0].plot(1953,235,'ko',mec=cl1,mfc='w',ms=ms,label='Whole stem\n(Johnstone 2002)')
		ax[0,0].plot(1998,277,'ko',mec=cl1,mfc='w',ms=ms)
	
		ax[0,0].plot(1953,76,'ks',mec=cl2,mfc='w',ms=ms,label='Merchantable\n(Johnstone 2002)')
		ax[0,0].plot(1998,233,'ks',mec=cl2,mfc='w',ms=ms)

	# Reconstruct from carbon
	flg=0
	if flg==1:
		lsat=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Scenario0001\\Inventory_Bat0001.pkl')
		V_WholeStemLive_FromC=mos[pNam]['Scenarios'][iB]['Mean']['C_Stemwood']['Ensemble Mean'][iT,iPS,iSS,iYS]/(lsat['Wood Density'][0,0])/meta['Param']['BEV']['Biophysical']['Carbon Content Wood']
		ax[0].plot(tv[iT],V_WholeStemLive_FromC,'--k',color=cl4,label='Whole stem (from C model)')
		
		V_MerchLive_FromC=mos[pNam]['Scenarios'][iB]['Mean']['C_StemwoodMerch']['Ensemble Mean'][iT,iPS,iSS,iYS]/(lsat['Wood Density'][0,0])/meta['Param']['BEV']['Biophysical']['Carbon Content Wood']
		ax[0].plot(tv[iT],V_MerchLive_FromC,'--k',color=cl3,label='Merchantable (from C model)')
	
	ax[0,0].set(ylabel='Live volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	ax[0,0].legend(loc='upper center',facecolor=[1,1,1],frameon=False,fontsize=5);
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	# Scenario 2
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_WholeStemLive']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,label='Whole stem')
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_MerchLive']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,label='Merchantable')

	# Johnstone 2002
	ax[1,0].plot(1953,108,'ko',mec=cl1,mfc='w',ms=ms)
	ax[1,0].plot(1998,324,'ko',mec=cl1,mfc='w',ms=ms)

	ax[1,0].plot(1953,41,'ks',mec=cl2,mfc='w',ms=ms)
	ax[1,0].plot(1998,290,'ks',mec=cl2,mfc='w',ms=ms)

	ax[1,0].set(ylabel='Live volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	#ax[1,0].legend(loc='center left',facecolor=[1,1,1],frameon=False);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	#--------------------------------------------------------------------------
	# Dead volume
	#--------------------------------------------------------------------------
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_WholeStemDead']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,label='Whole stem')
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_MerchDead']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,label='Merchantable')
	ax[0,1].set(ylabel='Dead volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	#ax[0,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,1].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_WholeStemDead']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,label='Whole stem')
	ax[1,1].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_MerchDead']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,label='Merchantable')
	ax[1,1].set(ylabel='Dead volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	#ax[1,1].legend(loc='center left',facecolor=[1,1,1],frameon=False);
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	#--------------------------------------------------------------------------
	# Live+Dead volume
	#--------------------------------------------------------------------------
	ax[0,2].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_WholeStemTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,label='Whole stem')
	ax[0,2].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['V_MerchTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,label='Merchantable')
	ax[0,2].set(ylabel='Live + dead volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	#ax[0,2].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])

	#ax[0,2].plot(1953,235,'ko',mec=cl1,mfc='w',ms=ms)
	#ax[0,2].plot(1998,277,'ko',mec=cl1,mfc='w',ms=ms)

	#ax[0,2].plot(1953,76,'ks',mec=cl2,mfc='w',ms=ms)
	#ax[0,2].plot(1998,233,'ks',mec=cl2,mfc='w',ms=ms)

	# Scenario 2
	ax[1,2].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_WholeStemTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],'-k',color=cl1,lw=lw,label='Whole stem')
	ax[1,2].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['V_MerchTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],'--k',color=cl2,lw=lw,label='Merchantable')
	ax[1,2].set(ylabel='Live + dead volume (m$^3$ ha$^{-1}$)',xlabel='Time, years',ylim=[0,yl])
	#ax[1,2].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])

	#ax[1,2].plot(1953,108,'ko',mec=cl1,mfc='w',ms=ms)
	#ax[1,2].plot(1998,324,'ko',mec=cl1,mfc='w',ms=ms)

	#ax[1,2].plot(1953,41,'ks',mec=cl2,mfc='w',ms=ms)
	#ax[1,2].plot(1998,290,'ks',mec=cl2,mfc='w',ms=ms)

	gu.axletters(ax,plt,0.04,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_VolumeTimeSeriesSummary','png',900)
	return

#%%
def Export_Summary_Tables(meta,pNam,mos):
	tbs={}
	SaveStatus='On'
	
	v='ByScenario'
	tbs[v]={}
	
	nam='Mean_t+0_to_t+100'
	t0=meta[pNam]['Project']['Year Project']+1
	t1=meta[pNam]['Project']['Year Project']+100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Mean_t+1_to_t+100'
	t0=meta[pNam]['Project']['Year Project']+1
	t1=meta[pNam]['Project']['Year Project']+100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Mean_in_t+1'
	t0=meta[pNam]['Project']['Year Project']+1
	t1=meta[pNam]['Project']['Year Project']+1
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)

	nam='Mean_in_t+100'
	t0=meta[pNam]['Project']['Year Project']+100
	t1=meta[pNam]['Project']['Year Project']+100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)

	nam='Mean_t+0_to_2050'
	t0=meta[pNam]['Project']['Year Project']
	t1=2050
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)

	nam='Mean_in_2050'
	t0=2050
	t1=2050
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Mean',t0=t0,t1=t1,Save=SaveStatus)

	nam='Sum_t-5_to_t+10'
	t0=meta[pNam]['Project']['Year Project']-5
	t1=meta[pNam]['Project']['Year Project']+10
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)

	nam='Sum_t-5_to_t+100'
	t0=meta[pNam]['Project']['Year Project']-5
	t1=meta[pNam]['Project']['Year Project']+100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_t+10'
	t0=meta[pNam]['Project']['Year Project']
	t1=meta[pNam]['Project']['Year Project']+10
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_t+50'
	t0=meta[pNam]['Project']['Year Project']
	t1=meta[pNam]['Project']['Year Project']+50
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_t+100'
	t0=meta[pNam]['Project']['Year Project']
	t1=meta[pNam]['Project']['Year Project']+100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_2030'
	t0=meta[pNam]['Project']['Year Project']
	t1=2030
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_2040'
	t0=meta[pNam]['Project']['Year Project']
	t1=2040
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_2050'
	t0=meta[pNam]['Project']['Year Project']
	t1=2050
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	nam='Sum_t+0_to_2100'
	t0=meta[pNam]['Project']['Year Project']
	t1=2100
	df,tbs[v][nam]=ExportTableByScenario(meta,pNam,mos,NameTable=nam,OperTime='Sum',t0=t0,t1=t1,Save=SaveStatus)
	
	v='Scenario Comparison'
	tbs[v]={}
	for cNam in mos[pNam]['Delta'].keys():
		tbs[v][cNam]={}
		nam='Mean_t+0_to_t+100'
		t0=meta[pNam]['Project']['Year Project']+1
		t1=meta[pNam]['Project']['Year Project']+100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Mean',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Mean_t+1_to_t+100'
		t0=meta[pNam]['Project']['Year Project']+1
		t1=meta[pNam]['Project']['Year Project']+100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Mean',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Mean_in_t+1'
		t0=meta[pNam]['Project']['Year Project']+1
		t1=meta[pNam]['Project']['Year Project']+1
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Mean',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Mean_in_t+100'
		t0=meta[pNam]['Project']['Year Project']+100
		t1=meta[pNam]['Project']['Year Project']+100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Mean',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t-5_to_t+10'
		t0=meta[pNam]['Project']['Year Project']-5
		t1=meta[pNam]['Project']['Year Project']+10
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()

		nam='Sum_t-5_to_t+100'
		t0=meta[pNam]['Project']['Year Project']-5
		t1=meta[pNam]['Project']['Year Project']+100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_t+10'
		t0=meta[pNam]['Project']['Year Project']
		t1=meta[pNam]['Project']['Year Project']+10
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_t+50'
		t0=meta[pNam]['Project']['Year Project']
		t1=meta[pNam]['Project']['Year Project']+50
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_t+100'
		t0=meta[pNam]['Project']['Year Project']
		t1=meta[pNam]['Project']['Year Project']+100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_2030'
		t0=meta[pNam]['Project']['Year Project']
		t1=2030
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_2040'
		t0=meta[pNam]['Project']['Year Project']
		t1=2040
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_2050'
		t0=meta[pNam]['Project']['Year Project']
		t1=2050
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()
	
		nam='Sum_t+0_to_2100'
		t0=meta[pNam]['Project']['Year Project']
		t1=2100
		tbs[v][cNam][nam]=ExportTableDelta(meta,pNam,mos,cNam=[cNam],NameTable=nam,OperSpace='Mean',OperTime='Sum',t0=t0,t1=t1,Units='Actual',Save=SaveStatus)[cNam].to_dict()

	return tbs
