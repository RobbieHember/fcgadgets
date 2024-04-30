'''
DEMO UTILITIES
'''
#%% Import Python modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
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
def ExportSummariesByScenario(meta,pNam,mos,**kwargs):
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

	for iScn in range(meta[pNam]['Project']['N Scenario']):
		d={}
		for k in mos[pNam]['Scenarios'][iScn]['Mean'].keys():
			if kwargs['operTime']=='Mean':
				d[k]=np.round(multip*np.mean(mos[pNam]['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS]),decimals=2)
			else:
				d[k]=np.round(multip*np.sum(mos[pNam]['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS]),decimals=2)

		if iScn==0:
			df=pd.DataFrame().from_dict(d,orient='index');
		else:
			df0=pd.DataFrame().from_dict(d,orient='index');
			df=pd.concat([df,df0],axis=1);

	df.columns=[np.arange(1,df.columns.size+1)];

	fout=meta['Paths'][pNam]['Data'] + '\\Outputs\\TabularSummary_' + kwargs['table_name'] + '_' + str(kwargs['t0']) + '-' + str(kwargs['t1']) + '_ProjectStrat' + str(iPS) + '_SpatialStrat' + str(iSS) + '_TimeStrat' + str(iYS) + '.xlsx'
	df.to_excel(fout)

	return df

#%% Custom scenario comparison tabular export
# Exmaple:
# udem.ExportTableDelta(meta,pNam,mos,'ScnComp1','Mean','Mean',[2020,2030])
def ExportTableDelta(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	if 'variables' in kwargs.keys():
		vL=kwargs['variables']
	else:
		vL=mos[pNam]['Delta'][ list(mos[pNam]['Delta'].keys())[0] ]['Data'].keys()
	df=pd.DataFrame()
	Names=[]
	for sc in mos[pNam]['Delta'].keys():
		if 'cnam' in kwargs.keys():
			if np.isin(sc,kwargs['cnam'])==False:
				continue
		Names.append(sc)
		d={}
		for v in vL:
			if v not in mos[pNam]['Delta'][sc]['Data'].keys():
				continue
			if (kwargs['units']=='Actual') | (kwargs['units']=='All'):
				if kwargs['oper']=='mean':
					y=np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				else:
					y=np.sum(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])
				d[v]=np.round(y,decimals=2)
			if (kwargs['units']=='Relative') | (kwargs['units']=='All'):
				if kwargs['oper']=='mean':
					y=np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.mean(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				else:
					y=np.sum(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/np.sum(mos[pNam]['Scenarios'][ mos[pNam]['Delta'][sc]['iB'] ][v]['Ensemble Mean'][iT,iPS,iSS,iYS])*100
				d[v + ' (%)']=np.round(y,decimals=2)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)
	df.index.name='Variable'
	df.columns=Names
	#df=df.sort_index(axis=0)
	if kwargs['save']=='On':
		df.to_excel(meta['Paths'][pNam]['Data'] + '\\Outputs\\TabularDelta_' + kwargs['table_name'] + '_' + kwargs['oper'] + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']) + '.xlsx')
	return df

#%% Plot mean fluxes and mean pools over a specified time horizon
def PlotSchematicAtmoGHGBal(meta,pNam,mos,**kwargs):

	# Key word arguments
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iSP']
		iSS=kwargs['iSS']
		iYS=kwargs['iYS']
	else:
		iPS=0
		iSS=0
		iYS=0
	iB=mos[pNam]['Delta'][ kwargs['cnam'] ]['iB']
	iP=mos[pNam]['Delta'][ kwargs['cnam'] ]['iP']
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	y_b={}
	y_p={}
	y_d={}
	for k in mos[pNam]['Scenarios'][0]['Mean'].keys():
		if (k=='C_Forest') | (k=='C_HWP') | (k=='C_ToMill'):
			y_b[k]=mos[pNam]['Scenarios'][iB]['Mean']['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iB]['Mean']['Mean'][k]['Ensemble Mean'][iT[0],iPS,iSS,iYS]
			y_p[k]=mos[pNam]['Scenarios'][iP]['Mean']['Mean'][k]['Ensemble Mean'][iT[-1],iPS,iSS,iYS]-mos[pNam]['Scenarios'][iP]['Mean']['Mean'][k]['Ensemble Mean'][iT[0],iPS,iSS,iYS]
			y_d[k]=y_p[k]-y_b[k]
		else:
			y_b[k]=np.sum(mos[pNam]['Scenarios'][iB]['Mean']['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS])
			y_p[k]=np.sum(mos[pNam]['Scenarios'][iP]['Mean']['Mean'][k]['Ensemble Mean'][iT,iPS,iSS,iYS])
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
	bx_Dom_FS_w=0.48
	bx_esc_w=0.15
	bx_atmo_bottom=0.88
	arrow_head_w=0.007
	arrow_lw=0.05
	fs_flux=6.5
	decim=1

	def GetSign(y):
		if y>0:
			x='+'
		else:
			x=''
		return x

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,10))

	# Background
	#ax.add_patch(Rectangle([0,1],0,1,fc=[0.98,0.97,0.96],ec='none'))

	# Atmosphere
	ax.add_patch(Rectangle([0.01,bx_atmo_bottom],0.98,0.1,fc=[0.9,0.95,1],ec=bx_ec))
	ax.text(0.5,0.935,'Atmosphere',size=bx_fs,ha='center',fontweight='bold',color=[0.08,0.3,0.55])
	vr='E_AGHGB_WSub'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Change in storage (tCO$_2$e/ha): ' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.5,0.9,txt,size=fs_flux+1,ha='center')

	# Dom_FS
	ax.add_patch(Rectangle([0.01,0.01],bx_Dom_FS_w,bx_lower_h,fc=[0.85,0.9,0.85],ec=bx_ec))
	ax.text(0.25,0.04,'Domestic Forest Sector (LULUCF)',size=bx_fs,ha='center',color=[0,0.5,0])

	# Forest land
	ax.add_patch(Rectangle([0.02,0.1],bx_Dom_FS_w*0.53,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
	ax.text(0.15,0.29,'Forest Land',size=bx_fs,ha='center',color=[0,0.5,0])
	ax.text(0.15,0.26,'Change in storage (tC/ha):',size=fs_flux,ha='center')
	vr='C_Forest'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt=str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.15,0.23,txt,size=fs_flux,ha='center')

	# Harvested wood products
	ax.add_patch(Rectangle([0.36,0.1],0.12,bx_lower_h-0.11,fc=[0.9,0.95,0.9],ec=bx_ec))
	ax.text(0.42,0.27,'Harvested\nWood\nProducts',size=bx_fs,ha='center',color=[0,0.5,0])
	vr='C_HWP'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Change in\nstorage (tC/ha):\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.42,0.185,txt,size=fs_flux,ha='center')

	# Lithosphere
	ax.add_patch( Rectangle([bx_Dom_FS_w+0.02,0.01],bx_esc_w*3+0.03+0.01,bx_lower_h,fc=[0.94,0.88,0.84],ec=bx_ec))
	ax.text(0.75,0.04,'Geological deposits (Fossil Fuels & Limestone)',size=bx_fs,ha='center',color=[0.5,0,0])

	# Energy - Stationary Combustion
	ax.add_patch(Rectangle([bx_Dom_FS_w+0.02+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
	ax.text(0.585,0.29,'Stationary\nCombustion',size=bx_fs,ha='center',color=[0.5,0,0])
	a1=np.round(y_d['E_Sub_ESC']/3.667,decimals=decim);
	a2=np.round(-1*y_d['E_Dom_ESC_ForOps']/3.667,decimals=decim);
	a3=-1*np.round((y_d['E_Sub_ESC']+y_d['E_Dom_ESC_ForOps'])/3.667,decimals=decim)
	txt='Change in\nstorage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.585,0.21,txt,size=fs_flux,ha='center')

	# Energy - Transportation
	ax.add_patch(Rectangle([bx_Dom_FS_w+0.02+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
	ax.text(0.745,0.29,'Transportation',size=bx_fs,ha='center',color=[0.5,0,0])
	a1=np.round(y_d['E_Sub_ET']/3.667,decimals=decim);
	a2=np.round(-1*y_d['E_Dom_ET_ForOps']/3.667,decimals=decim);
	a3=-1*np.round((y_d['E_Sub_ET']+y_d['E_Dom_ET_ForOps'])/3.667,decimals=decim)
	txt='Change in\n storage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.745,0.21,txt,size=fs_flux,ha='center')

	# IPPU
	ax.add_patch(Rectangle([bx_Dom_FS_w+0.02+0.01+bx_esc_w+0.01+bx_esc_w+0.01,0.1],bx_esc_w,bx_lower_h-0.11,fc=[1,0.95,0.9],ec=bx_ec))
	ax.text(0.905,0.27,'Industrial\nProcesses &\nProduct Use',size=bx_fs,ha='center',color=[0.5,0,0])
	a1=np.round(y_d['E_Sub_IPPU']/3.667,decimals=decim);
	a2=np.round(-1*y_d['E_Dom_IPPU_ForOps']/3.667,decimals=decim);
	a3=-1*np.round((y_d['E_Sub_IPPU']+y_d['E_Dom_IPPU_ForOps'])/3.667,decimals=decim)
	txt='Change in\n storage (tC/ha):\n' + GetSign(a1) + str(a1) + ',' + GetSign(a2) + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.905,0.185,txt,size=fs_flux,ha='center')

	#--------------------------------------------------------------------------
	# Fluxes
	#--------------------------------------------------------------------------

	# NEE
	vr='E_Dom_FS_NEE'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Net ecosystem\nexchange\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.025,0.79,txt,ha='left',size=fs_flux)
	ax.arrow(0.02,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Wildfire
	vr='E_Dom_FS_Wildfire'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Wildfire\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
	ax.text(0.125,0.52,txt,ha='right',size=fs_flux)
	ax.arrow(0.13,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Open burning
	vr='E_Dom_FS_OpenBurning'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Open burning\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.155,0.52,txt,ha='left',size=fs_flux)
	ax.arrow(0.15,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Denitrification
	vr='E_Dom_FS_Denit'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Denitrification\n' + str(a1) + ',' + str(a2) + ' (' + str(a3) + ')'
	ax.text(0.245,0.71,txt,ha='right',size=fs_flux)
	ax.arrow(0.25,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Volatilization
	vr='E_Dom_FS_Other'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Volatilization\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.275,0.71,txt,ha='left',size=fs_flux)
	ax.arrow(0.27,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# HWP fluxes
	vr='E_Dom_FS_HWP'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Product decay\nand combustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.415,0.58,txt,ha='right',va='top',size=fs_flux)
	ax.arrow(0.42,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Removals
	vr='C_ToMill'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	#a1=0.2; a2=0.5; a3=0.3
	txt='Removals\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.315,0.31,txt,ha='center',size=fs_flux)
	ax.arrow(0.275,0.275,0.08,0,head_width=0.01,head_length=arrow_head_w,fc='k',ec='k',lw=arrow_lw)

	# Bioenergy combustion
	vr='E_ESC_Bioenergy'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Bioenergy\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.525-0.08,0.58,txt,ha='left',va='top',size=fs_flux)
	ax.arrow(0.52-0.08,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# ESC operational emissions
	vr='E_Dom_ESC_ForOps'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel for\nstationary\ncombustion\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.55,0.64,txt,ha='left',va='top',size=fs_flux)
	ax.arrow(0.545,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Displacement energy
	vr='E_Sub_E'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Displacement \neffects of\nbioenergy\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.655,0.82,txt,ha='right',va='top',size=fs_flux)
	ax.arrow(0.66,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Transportation
	vr='E_Dom_ET_ForOps'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel for\ntransportation\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.695,0.64,txt,ha='left',va='top',size=fs_flux)
	ax.arrow(0.69,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# Displacement building materials
	vr='E_Sub_M'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Displacement \neffects of\nsolid wood\nmaterials\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.805,0.82,txt,ha='left',va='top',size=fs_flux)
	ax.arrow(0.8,bx_atmo_bottom,0,-1*(bx_atmo_bottom-bx_lower_h-0.02),head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	# IPPU
	vr='E_Dom_IPPU_ForOps'
	a1=np.round(y_b[vr],decimals=decim); a2=np.round(y_p[vr],decimals=decim); a3=np.round(y_d[vr],decimals=decim)
	txt='Fossil fuel\ncombustion\nand\nurea\nsequestration\n' + str(a1) + ',' + str(a2) + ' (' + GetSign(a3) + str(a3) + ')'
	ax.text(0.915,0.64,txt,ha='left',va='top',size=fs_flux)
	ax.arrow(0.91,bx_lower_h+0.01,0,bx_atmo_bottom-bx_lower_h-0.02,head_width=arrow_head_w,head_length=0.01,fc='k',ec='k',lw=arrow_lw)

	plt.tight_layout()
	ax.set(position=[0,0,1,1],visible='Off',xticks=[],yticks=[])
	ax.axis('off')

	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AGHGB Schematic_S' + str(iP) + 'minusS' + str(iB) + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']),'png',900)
	return

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
	cnam=kwargs['cnam']
	iB=mos[pNam]['Delta'][cnam]['iB']
	iP=mos[pNam]['Delta'][cnam]['iP']
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	cl=np.array([[0.17,0.34,0.69],[0.55,0.9,0],[0.5,0,1],[0,1,1]])

	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(20,11)); lw=0.75
	# Cost
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label=lab1)
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label=lab2)
	ax[0,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cost (CAD x 1000)')
	ax[0,0].legend(loc='upper right',frameon=False,facecolor=None)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta cost
	ax[0,1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata']['Mean']['Cost Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw)
	ax[0,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ Cost (CAD x 1000)')
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Net Revenue
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label='Baseline')
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label='Project')
	ax[1,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Net revenue (CAD x 1000)')
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta net revenue
	ax[1,1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata']['Mean']['Revenue Net']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw)
	ax[1,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ net revenue (CAD x 1000)')
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Cumulative net revenue
	ax[2,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean']['Revenue Net_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-bo',color=cl[0,:],mfc=cl[0,:],mec=cl[0,:],lw=lw,ms=4,label='Baseline')
	ax[2,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean']['Revenue Net_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--r^',color=cl[1,:],mfc=cl[1,:],mec=cl[1,:],lw=lw,ms=2,label='Project')
	ax[2,0].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='Cumulative net revenue\n (CAD x 1000)')
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Delta cumulative net revenue
	ax[2,1].plot(tv[iT],np.zeros(iT.size),'-k',lw=3,color=[0.85,0.85,0.85])
	ax[2,1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata']['Mean']['Revenue Net_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'-k',color=cl[2,:],mfc=cl[2,:],mec=cl[2,:],lw=lw,label='Undiscounted')
	ax[2,1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata']['Mean']['Revenue Net Disc_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS]/1000,'--k',color=0.5*cl[2,:],mfc=0.5*cl[2,:],mec=0.5*cl[2,:],lw=lw,label='Discounted')
	ax[2,1].set(xlim=[tv[iT[0]], tv[iT[-1]]],ylabel='$\Delta$ cumulative net\n revenue (CAD x 1000)')
	ax[2,1].legend(loc='lower left',frameon=False,facecolor=None)
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Cashflow_' + cnam,'png',900)
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

	vs=['C_Biomass','C_DeadWood','C_Litter','C_Soil','C_InUse','C_DumpLandfill','C_Geological','E_AGHGB_WSub_cumu']
	vs2=['Biomass','Dead Wood','Litter','Soil','In-use Products','Dump and Landfill','Geological','Atmosphere']

	operS=kwargs['operSpace']

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
						if vs[cnt]=='E_AGHGB_WSub_cumu':
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
	elif 'cnam' in kwargs.keys():
		# Generate one figure per scenario comparison
		cnam=kwargs['cnam']
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
		sL=[mos[pNam]['Delta'][cnam]['iB'],mos[pNam]['Delta'][cnam]['iP']]
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
						if (be_d[iT2]>0) & (be0[iT2]<0) & (be1[iT2]<0):
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						elif (be_d[iT2]>0):
							tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						else:
							tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + uni
							mlti=1.15
						ax[i,j].text(tv[iT][iT2],mlti*be1[iT2],tx,fontsize=7,fontweight='bold',ha='center')


				if (i==0) & (j==0):
					if 'LegendLoc' in kwargs.keys():
						legloc=kwargs['LegendLoc']
					else:
						legloc='lower left'
					ax[i,j].legend(loc=legloc,frameon=False,facecolor=None,edgecolor='w')
				if vs[cnt]=='E_AGHGB_WSub_cumu':
					us='\n(tCO2e/ha)'
				else:
					us='\n(tC/ha)'
				ax[i,j].set(ylabel=vs2[cnt] + us,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				cnt=cnt+1
		plt.tight_layout()
		gu.axletters(ax,plt,0.02,0.82,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Pools_' + cnam,'png',900)

	return

#%%
def PlotFluxes(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	vs=['E_Dom_FS_NPP','C_G_Net_Reg','E_Dom_FS_RH','E_Dom_FS_OpenBurning','E_Dom_FS_Wildfire','E_Dom_FS_HWP','E_Sub','E_AGHGB_WSub']
	vs2=['NPP','Net growth','RH','Open burning','Wildfire','HWP','Substitutions','Net emissions']

	operS=kwargs['operSpace']

	cl=np.array([[0.27,0.44,0.79],[0.55,0.9,0],[0.5,0,1],[0,1,1]])
	symb=['-','--','-.',':','-']
	if 'ScenarioIndexList' in kwargs.keys():
		# Generate one figure for a custom set of scenarios
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(20,12)); Alpha=0.09
		for i in range(4):
			for j in range(2):
				for iScn in range(len(kwargs['ScenarioIndexList'])):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Dom_FS_NPP':
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
				#if vs[cnt]=='E_Dom_FS_NPP':
				#	ylab='-1 x ' + vs2[cnt] + '\n(tCO2e/ha/yr)'
				#else:
				ylab=vs2[cnt] + '\n(tCO2e/ha/yr)'
				ax[i,j].set(ylabel=ylab + ' (tCO2e/ha/yr)',xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			cnt=cnt+1
		gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_CustomScenarioList','png',900)
	elif 'cnam' in kwargs.keys():
		# Generate one figure per scenario comparison
		k=kwargs['cnam']
		cnt=0
		fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(20,12)); Alpha=0.09
		sL=[mos[pNam]['Delta'][k]['iB'],mos[pNam]['Delta'][k]['iP']]
		for i in range(4):
			for j in range(2):
				for iScn in range(len(sL)):
					adj=1.0
					if vs[cnt]=='C_G_Net_Reg':
						adj=3.667
					if vs[cnt]=='E_Dom_FS_NPP':
						adj=-1.0
					be=adj*mos[pNam]['Scenarios'][sL[iScn]][operS][vs[cnt]]['Ensemble Mean'][iT,iPS,iSS,iYS]
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
				#if vs[cnt]=='E_Dom_FS_NPP':
				#	ylab='-1 x ' + vs2[cnt] + '\n(tCO2e/ha/yr)'
				#else:
				ylab=vs2[cnt] + '\n(tCO2e/ha/yr)'
				ax[i,j].set(ylabel=ylab,xlabel='Time, years',xlim=[np.min(tv[iT]),np.max(tv[iT])])
				cnt=cnt+1
		gu.axletters(ax,plt,0.02,0.82,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Fluxes_' + k,'png',900)

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

	cnam=kwargs['cnam']

	operS=kwargs['operSpace']

	cl=np.array([[0.75,0,0],[0,0.85,0],[0.29,0.49,0.77]])
	ymin=0.0;ymax=0.0

	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(FigSize[0],FigSize[1])); Alpha=0.09
	for j in range(0,2):
		ax[j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')
	vn='C_RH'
	lo=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=cl[0,:],label='RH')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	vn='C_NPP'
	lo=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'-.',color=cl[1,:],label='NPP')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	vn='E_Dom_FS_NEE'
	lo=-mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=-mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=-mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=-mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS]
	#ax[0].fill_between(tv[iT],lo,hi,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],-mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=cl[2,:],label='NEP')
	ymin=np.minimum(ymin,np.min(lo))
	ymax=np.maximum(ymax,np.max(hi))

	ax[0].legend(loc="lower right",frameon=False,facecolor=None,edgecolor='w')
	ax[0].set(ylabel='Annual $\Delta$ (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	ymin=0.0
	ymax=0.0

	vn='C_RH'
	lo=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[0,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'--',color=cl[0,:])
	ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]]);
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))

	vn='C_NPP'
	lo=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=np.cumsum(3.667*mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[1,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'-.',color=cl[1,:])
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))

	vn='E_Dom_FS_NEE'
	lo=-np.cumsum(mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P025'][iT,iPS,iSS,iYS])
	hi=-np.cumsum(mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P975'][iT,iPS,iSS,iYS])
	lo2=-np.cumsum(mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P250'][iT,iPS,iSS,iYS])
	hi2=-np.cumsum(mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble P750'][iT,iPS,iSS,iYS])
	mu=-np.cumsum(mos[pNam]['Delta'][cnam]['ByStrata'][operS][vn]['Ensemble Mean'][iT,iPS,iSS,iYS])
	#ax[1].fill_between(tv[iT],lo,hi,color=[0.75,.5,1],alpha=Alpha,linewidth=0,label='95 C.I.')
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[2,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mu,'-',color=cl[2,:])
	ymin=np.minimum(ymin,np.min(mu))
	ymax=np.maximum(ymax,np.max(mu))
	ax[1].set(ylabel='Cumulative $\Delta$ (tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',ylim=[ymin,ymax],xlim=[tv[iT[0]],tv[iT[-1]]]);
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.035,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\NEE_Balance_' + cnam,'png',900)

	return

#%%
def PlotDeltaGHGB(meta,mos,pNam,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]

	cnam=kwargs['cnam']

	operS=kwargs['operSpace']

	iB=mos[pNam]['Delta'][cnam]['iB']
	iP=mos[pNam]['Delta'][cnam]['iP']

	clD=np.array([0.6,0.1,1])

	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,9)); Alpha1=0.06; Alpha2=0.08; Alpha3=0.11
	for i in range(0,2):
		for j in range(0,2):
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both')
			ax[i,j].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean'][operS]['E_AGHGB_WOSub']['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=(0,0.5,1),label='Baseline')
	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean'][operS]['E_AGHGB_WOSub']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=(0,0.6,0),label='Project (With Subs.)')
	ax[0,0].legend(loc="upper right",frameon=False,facecolor=None,edgecolor='w')
	ax[0,0].set(ylabel='Annual E (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Mean'][operS]['E_AGHGB_WOSub_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS],'-',color=(0,0.5,1),label='Baseline SR')
	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Mean'][operS]['E_AGHGB_WOSub_cumu']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=(0,0.6,0),label='Baseline NSR')
	ax[0,1].set(ylabel='Cumulative E (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	mu=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble Mean'][iT,iPS,iSS,iYS]
	lo=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble P750'][iT,iPS,iSS,iYS]
	sig=1*mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble SD'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=clD,alpha=Alpha1,linewidth=0,label='95 C.I.')
	ax[1,0].fill_between(tv[iT],lo2,hi2,color=clD,alpha=Alpha2,linewidth=0,label='50 C.I.')
	ax[1,0].fill_between(tv[iT],mu-sig,mu+sig,color=0.5*clD,alpha=Alpha3,linewidth=0,label='S.D.')
	ax[1,0].plot(tv[iT],mu,'-',color=clD,label='Best estimate (W/O Subs.)')
	ax[1,0].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=0.5*clD,label='Best estimate (With Subs.)')
	#ax[1,0].legend(loc="upper right",frameon=0)
	ax[1,0].legend(loc="lower left",frameon=0)
	ax[1,0].set(ylabel='$\Delta$E (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);
	if 'yLimPad' in kwargs.keys():
		ax[1,0].set(ylim=[-kwargs['yLimPad']*np.max(np.abs(mu)),kwargs['yLimPad']*np.max(np.abs(mu))])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	mu=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub_cumu_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS]
	lo=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub_cumu_from_tref']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub_cumu_from_tref']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub_cumu_from_tref']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub_cumu_from_tref']['Ensemble P750'][iT,iPS,iSS,iYS]
	sig=1*mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WOSub']['Ensemble SD'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=clD,alpha=Alpha1,linewidth=0,label='95 C.I.')
	ax[1,1].fill_between(tv[iT],lo2,hi2,color=clD,alpha=Alpha2,linewidth=0,label='50 C.I.')
	ax[1,1].fill_between(tv[iT],mu-sig,mu+sig,color=clD,alpha=Alpha3,linewidth=0,label='S.D.')
	ax[1,1].plot(tv[iT],mu,'-',color=clD,label='Best estimate (With Subs.)')
	ax[1,1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS],'--',color=0.5*clD,label='Best estimate (W/O Subs.)')
	ax[1,1].set(ylabel='Cumulative $\Delta$E (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	plt.tight_layout()
	if 'yLimPad' in kwargs.keys():
		ax[1,1].set(ylim=[-kwargs['yLimPad']*np.max(np.abs(mu)),kwargs['yLimPad']*np.max(np.abs(mu))])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GHG_Balance_' + cnam,'png',900)

	return

#%%
def PlotGHGBenefit(meta,pNam,mos,tv,iT,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	operS=kwargs['operSpace']

	cnam=kwargs['cnam']

	cl=np.array([[0,0.5,1],[0,0.6,0],[0,1,1],[0.5,0,1]])
	cnt=1
	symb=['-','--','-.',':','-']
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(16,6)); Alpha=0.09
	for i in range(0,2):
		ax[i].yaxis.set_ticks_position('both'); ax[i].xaxis.set_ticks_position('both')
		ax[i].plot(tv[iT],0*np.ones(tv[iT].shape),'-',lw=3,color=(0.8,0.8,0.8),label='')

	lo=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[0].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[0].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[0].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub']['Ensemble Mean'][iT,iPS,iSS,iYS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1) )
	ax[0].legend(loc="upper right",frameon=False,facecolor=None,edgecolor='w')
	ax[0].set(ylabel='$\Delta$GHG (tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years',ylim=[np.min(lo),np.max(hi)]);
	lo=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble P975'][iT,iPS,iSS,iYS]
	lo2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble P250'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1].fill_between(tv[iT],lo,hi,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[1].fill_between(tv[iT],lo2,hi2,color=cl[cnt,:],alpha=Alpha,linewidth=0)
	ax[1].plot(tv[iT],mos[pNam]['Delta'][cnam]['ByStrata'][operS]['E_AGHGB_WSub_cumu_from_tref']['Ensemble Mean'][iT,iPS,iSS,iYS],symb[cnt],color=cl[cnt,:],label='SC ' + str(cnt+1))
	ax[1].set(ylabel='Cumulative $\Delta$GHG (tCO$_2$e ha$^-$$^1$)',xlim=[tv[iT[0]],tv[iT[-1]]],xlabel='Time, years');
	gu.axletters(ax,plt,0.035,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['gp']['AxesLetterStyle'],FontWeight=meta['Graphics']['gp']['AxesLetterFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GHG_Benefit_' + cnam,'png',900)

	return

#%% Summary bioenergy description
def SummaryBioenergy(meta,mos,sc,th):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=meta[pNam]['Project']['Year Project']) & (tv<=meta[pNam]['Project']['Year Project']+th) )[0]

	d={}
	d['EI PelletExport (tCO2e/ODT)']=np.sum(mos[pNam]['Delta'][sc]['Data']['E_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['Data']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI PelletExport Boiler (tCO2e/GJ)']=np.sum(mos[pNam]['Delta'][sc]['Data']['E_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['Data']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI PelletExport Boiler+Ops (tCO2e/GJ)']=np.sum( (mos[pNam]['Delta'][sc]['Data']['E_ESC_BioenergyPelletExport']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['Data']['E_ESC_ForOps']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['Data']['E_ET_ForOps']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['Data']['E_IPPU_ForOps']['Ensemble Mean'][iT,0,0]+ \
													   mos[pNam]['Delta'][sc]['Data']['E_Dom_FS_NEE']['Ensemble Mean'][iT,0,0]) )/np.sum(mos[pNam]['Delta'][sc]['Data']['GJ PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI Pellet Manufacture (tCO2e/ODT Pellets)']=np.sum( (mos[pNam]['Delta'][sc]['Data']['E_ESC_ForOps']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['Data']['E_IPPU_ForOps']['Ensemble Mean'][iT,0,0]) )/np.sum(mos[pNam]['Delta'][sc]['Data']['ODT PelletExport']['Ensemble Mean'][iT,0,0])
	d['EI OperationForestry (tCO2e/ODT)']= \
			np.sum( (mos[pNam]['Delta'][sc]['Data']['E_ESC_ForOps']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['Data']['E_ET_ForOps']['Ensemble Mean'][iT,0,0]+mos[pNam]['Delta'][sc]['Data']['E_IPPU_ForOps']['Ensemble Mean'][iT,0,0]) ) / \
			np.sum( (mos[pNam]['Delta'][sc]['Data']['ODT Lumber']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT LogExport']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT Plywood']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT OSB']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT MDF']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT Paper']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT PelletExport']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT PelletDomGrid']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT PelletDomRNG']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT PowerFacilityDom']['Ensemble Mean'][iT,0,0]+ \
			mos[pNam]['Delta'][sc]['Data']['ODT PowerGrid']['Ensemble Mean'][iT,0,0]) )

	d['Energy efficiency (GJ/ODT)']=np.sum(mos[pNam]['Delta'][sc]['Data']['GJ PelletExport']['Ensemble Mean'][iT,0,0])/np.sum(mos[pNam]['Delta'][sc]['Data']['ODT PelletExport']['Ensemble Mean'][iT,0,0])

	# Displacement factor

	SubTot=-np.sum(mos[pNam]['Delta'][sc]['Data']['E_Sub']['Ensemble Mean'][iT,0,0])
	SubE=-np.sum(mos[pNam]['Delta'][sc]['Data']['E_Sub_E']['Ensemble Mean'][iT,0,0])
	SubM=-np.sum(mos[pNam]['Delta'][sc]['Data']['E_Sub_M']['Ensemble Mean'][iT,0,0])
	Wood=np.sum(mos[pNam]['Delta'][sc]['Data']['C_ToLumber']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['Data']['C_ToPlywood']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['Data']['C_ToMDF']['Ensemble Mean'][iT,0,0]+ \
		mos[pNam]['Delta'][sc]['Data']['C_ToOSB']['Ensemble Mean'][iT,0,0])
	Bioenergy=np.sum(mos[pNam]['Delta'][sc]['Data']['E_ESC_Bioenergy']['Ensemble Mean'][iT,0,0])
	BioenergyPlusOps=np.sum(mos[pNam]['Delta'][sc]['Data']['E_ESC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_ESC_ForOps']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_ET_ForOps']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_IPPU_ForOps']['Ensemble Mean'][iT,0,0])
	BioenergyOpsAndHWPDecay=np.sum(mos[pNam]['Delta'][sc]['Data']['E_ESC_Bioenergy']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_ESC_ForOps']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_ET_ForOps']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_IPPU_ForOps']['Ensemble Mean'][iT,0,0]+ \
							mos[pNam]['Delta'][sc]['Data']['E_ESC_Bioenergy']['Ensemble Mean'][iT,0,0])

	d['Displacement Factor Total (tC/tC)']=SubTot/BioenergyOpsAndHWPDecay
	d['Displacement Factor Energy (tC/tC)']=SubE/BioenergyPlusOps
	d['Displacement Factor Materials (tC/tC)']=(SubM/3.667)/Wood

	# Save
	df=pd.DataFrame.from_dict(d,orient='index')
	df.to_excel(meta['Paths']['Project'] + '\\Outputs\\SummaryBioenergy.xlsx')

	return d

#%%
def NA_CompareSpecifications_ChangeInPools(meta,pNam,mos,tv):
	iPS=0; iSS=0; iYS=0
	iT=np.where((tv>=2000) & (tv<=2100))[0]
	cl=np.array([[0.6,0.95,0],[0,0.75,0],[0.95,0.7,0.7],[0.55,0.3,0.3],[0.25,0.1,0.1],[0.27,0.49,0.76]])
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,9))
	ax[0,0].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	sc='NGS'
	AGB=mos[pNam]['Delta'][sc]['Data']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['Data']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['Data']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['Data']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['Data']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[0,0].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[0,0].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[0,0].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[0,0].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[0,0].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[0,0].legend(loc="upper right",frameon=False,facecolor=None)
	ax[0,0].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (MgC/ha)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,25])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT'
	AGB=mos[pNam]['Delta'][sc]['Data']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['Data']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['Data']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['Data']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['Data']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[0,1].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[0,1].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[0,1].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[0,1].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[0,1].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[0,1].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[0,1].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (MgC/ha)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,25])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT+T'
	AGB=mos[pNam]['Delta'][sc]['Data']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['Data']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['Data']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['Data']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['Data']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,0].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[1,0].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[1,0].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[1,0].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[1,0].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[1,0].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[1,0].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (MgC/ha)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,25])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	sc='NGT+T+D'
	AGB=mos[pNam]['Delta'][sc]['Data']['C_Biomass']['Ensemble Mean'][iT,iPS,iSS,iYS]-mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	BGB=mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean'][iT,iPS,iSS,iYS]
	DW=mos[pNam]['Delta'][sc]['Data']['C_DeadWood']['Ensemble Mean'][iT,iPS,iSS,iYS]
	LT=mos[pNam]['Delta'][sc]['Data']['C_Litter']['Ensemble Mean'][iT,iPS,iSS,iYS]
	SOC=mos[pNam]['Delta'][sc]['Data']['C_Soil']['Ensemble Mean'][iT,iPS,iSS,iYS]
	FOR=mos[pNam]['Delta'][sc]['Data']['C_Forest']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].plot([0,10000],[0,0],'k-',color=[0.8,0.8,0.8],lw=3)
	ax[1,1].plot(tv[iT],AGB,'-',color=cl[0,:],label='Biomass (aboveground)',linewidth=1.5)
	ax[1,1].plot(tv[iT],BGB,'-',color=cl[1,:],label='Biomass (belowground)',linewidth=1.5)
	ax[1,1].plot(tv[iT],DW,'-',color=cl[2,:],label='Dead wood',linewidth=1.5)
	ax[1,1].plot(tv[iT],LT,'-',color=cl[3,:],label='Litter',linewidth=1.5)
	ax[1,1].plot(tv[iT],SOC,'-',color=cl[4,:],label='Soil organic matter',linewidth=1.5)
	ax[1,1].plot(tv[iT],FOR,'-',color=cl[5,:],label='Forest ecosystem',linewidth=1.5)
	ax[1,1].set(xticks=np.arange(2000,2110,10),ylabel=r'Change in stocks (MgC/ha)',xlabel='Time, years',xlim=[2000,2100],ylim=[0,25])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	gu.axletters(ax,plt,0.02,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=['NGS','NGT','NGT+T','NGT+T+D'],LabelSpacer=0.021)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\CompareSpecifications_ChangeInPools','png',900);
	return

#%%
def NA_CumulativeGHGBenefit_WithCI(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,8)); lw=0.75; cl1=[0.96,0.92,1]; cl2=[0.9,0.82,0.96]
	# Without harvesting
	sc='Coast No Harvest'
	#sc='Interior No Harvest'
	v='E_AGHGB_WOSub'
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[0,0].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='95% percentile range')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['Data'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
		ax[0,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1,label='3 x S.E. range')
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
		ax[0,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1,label='50% percentile range')
	mu=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[0,0].set(yticks=np.arange(-40,50,10),ylabel='$\Delta$E\n(tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',ylim=[-35,35],xlim=[ tv[iT[0]],tv[iT[-1]] ])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,0].legend(loc=kwargs['LegendPosition'],frameon=False,facecolor=None)
	ax[0,0].annotate('Application',xy=(2020,5),xytext=(2020,24),arrowprops=dict(arrowstyle="->"),ha='center')
	ax[0,0].text(2045,20,'Sink',ha='center',va='center',fontsize=8,style='italic',weight='bold',color=[0.5,0.5,0.5])
	ax[0,0].text(2045,-20,'Source',ha='center',va='center',fontsize=8,style='italic',weight='bold',color=[0.5,0.5,0.5])
	#ax[0].annotate('Manufacture, transport & N2O emissions',xy=(2020,-6),xytext=(2020,-16),arrowprops=dict(arrowstyle="->"),ha='center')
	
	v='E_AGHGB_WOSub_cumu_from_tref'
	ax[0,1].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['Data'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[0,1].yaxis.set_ticks_position('both');ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].set(yticks=np.arange(-260,160,20),ylabel='Cumuulative $\Delta$E\n(tCO$_2$e ha$^-$$^1$)',ylim=[-120,20],xlim=[ tv[iT[0]],tv[iT[-1]] ]);

	#------------------------------------------------------------------------------
	# With harvesting
	#------------------------------------------------------------------------------
	sc='Coast With Harvest'
	#sc='Interior With Harvest'
	v='E_AGHGB_WOSub'
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['Data'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[1,0].set(yticks=np.arange(-40,50,10),ylabel='$\Delta$E\n(tCO$_2$e ha$^-$$^1$ yr$^-$$^1$)',ylim=[-35,35],xlim=[2010,2080]);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].annotate('Post-application\nharvest',xy=(2043,7),xytext=(2043,22),arrowprops=dict(arrowstyle="->"),ha='center')

	ax[1,1].plot(tv[iT],np.zeros(iT.size),'-',color=[0.8,0.8,0.8],lw=1.5)
	v='E_AGHGB_WOSub_cumu_from_tref'
	lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=cl1,linewidth=0,zorder=1,label='Mean $\pm$ 1.0 S.D.')
	if kwargs['Error']=='3SE':
		se=3*mos[pNam]['Delta'][sc]['Data'][v]['Ensemble SE'][iT,iPS,iSS,iYS]
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]-se
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]+se
	elif kwargs['Error']=='50%':
		lo=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P250'][iT,iPS,iSS,iYS]
		hi=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble P750'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=cl2,linewidth=0,zorder=1)
	mu=mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].plot(tv[iT],mu,'-',color=[0.25,0,0.5],label='Mean',linewidth=lw)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,1].set(yticks=np.arange(-260,160,20),ylabel='Cumulative $\Delta$E\n(tCO$_2$e ha$^-$$^1$)',xlabel='Time, years',ylim=[-120,20],xlim=[2010,2080]);
	
	gu.axletters(ax,plt,0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\CumulativeGHGBenefit_WithCI','png',900)
	return

#%%
def PlotVolume(meta,pNam,mos,**kwargs):
	if 'iPS' in kwargs.keys():
		iPS=kwargs['iPS']; iSS=kwargs['iSS']; iYS=kwargs['iYS']
	else:
		iPS=0; iSS=0; iYS=0

	if 'FigSize' in kwargs.keys():
		FigSize=kwargs['FigSize']
	else:
		FigSize=[20,12]

	cnam=kwargs['cnam']

	operS=kwargs['operSpace']

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where((tv>=kwargs['t0']) & (tv<=kwargs['t1']))[0]

	iB=mos[pNam]['Delta'][cnam]['iB']
	iP=mos[pNam]['Delta'][cnam]['iP']
	if 'ScenarioLabels' in kwargs.keys():
		lab1=kwargs['ScenarioLabels'][0]
		lab2=kwargs['ScenarioLabels'][1]
	else:
		lab1=meta[pNam]['Scenario'][iB]['Scenario_CD']
		lab2=meta[pNam]['Scenario'][iP]['Scenario_CD']

	#v='C_Biomass'
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
			if be_d[iT2]>0:
				tx='+' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			elif np.isnan(be_d[iT2])==True:
				tx=''
			else:
				tx='-' + str(np.round(be_d[iT2[0]],decimals=1)) + '%'
			ax.text(tv[iT][iT2],1.05*y2[iT2],tx,fontsize=7,fontweight='bold',ha='center',color=[0.3,0.45,0])

	ax.set(ylabel=r'Total net merch. stemwood volume (m$^3$ ha$^-$$^1$)',xlabel=r'Time, years',ylim=[0,ymx],xlim=[np.min(tv[iT]),np.max(tv[iT])]);
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
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Volume_' + cnam,'png',900)
	return

#%%
def NA_CalcNUE(meta,pNam,mos,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iT=np.where( (tv>=kwargs['t0']) & (tv<=kwargs['t1']) )[0]
	doseN=200 # NUE applied
	df=pd.DataFrame()
	Names=[]
	for sc in mos[pNam]['Delta'].keys():
		if 'cnam' in kwargs.keys():
			if np.isin(sc,kwargs['cnam'])==False:
				continue
		Names.append(sc)
		v='C_AboveGroundBiomass'
		Temp=mos[pNam]['Delta'][sc]['Data']['C_Biomass']['Ensemble Mean']-mos[pNam]['Delta'][sc]['Data']['C_Root']['Ensemble Mean']
		d={}
		v='C_AboveGroundBiomass'
		y=1000*np.mean(Temp)/doseN
		d[v]=np.round(y,decimals=2)
		v='C_Root'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v]=np.round(y,decimals=2)
		v='C_DeadWood'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v]=np.round(y,decimals=2)
		v='C_Litter'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v]=np.round(y,decimals=2)
		v='C_Soil'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v]=np.round(y,decimals=2)
		v='C_Forest'
		y=1000*np.mean(mos[pNam]['Delta'][sc]['Data'][v]['Ensemble Mean'][iT,iPS,iSS,iYS])/doseN
		d[v]=np.round(y,decimals=2)
		df0=pd.DataFrame().from_dict(d,orient='index')
		df=pd.concat([df,df0],axis=1)
	#df.index.name='Variable'
	df.columns=Names
	#df=df.sort_index(axis=0)
	if 'save' in kwargs.keys():
		df.to_excel(meta['Paths'][pNam]['Data'] + '\\Outputs\\TabularSummaryDelta_' + kwargs['table_name'] + '_' + str(kwargs['t0']) + 'to' + str(kwargs['t1']) + '.xlsx')
	return df