'''
INVENTORY UTILITIES
'''
#%% Import modules
import numpy as np
import time
import copy
import pyproj
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_general as gu
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.na1k.na1k_util as u1k
import fcgadgets.taz.aspatial_stat_models as asm

#%% Import variables
def Process1_ImportVariables(meta,pNam):

	#--------------------------------------------------------------------------
	# Define regular grid sampling frequency and name of mask
	#--------------------------------------------------------------------------
	rgsf=str(meta['Geos']['RGSF'])
	
	if meta[pNam]['Project']['Code Project']=='BCFCS_LUC':
		mask='BCFCS_LUC'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_EvalAtPlots':
		mask='BCFCS_EValAtPlots'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_EvalAtCN':
		mask='BCFCS_EvalAtCN'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_NOSEC':
		mask='BCFCS_NOSEC'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_NOSECs':
		mask='BCFCS_NOSECs'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_NMC':
		mask='BCFCS_NMC'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_Eval':
		mask='BCFCS_Eval' 
	elif meta[pNam]['Project']['Code Project']=='BCFCS_EvalCoast':
		mask='BCFCS_EvalCoast'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_EvalInterior':
		mask='BCFCS_EvalInterior'
	elif meta[pNam]['Project']['Code Project']=='TSA_DawsonCreek':
		mask='TSA_DawsonCreek'
	elif meta[pNam]['Project']['Code Project']=='Landscape_NicolaRiverWatershed':
		mask='Landscape_NicolaRiverWatershed'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_CWH':
		mask='BCFCS_CWH'
	elif meta[pNam]['Project']['Code Project']=='BCFCS_SBS':
		mask='BCFCS_SBS'
	else:
		# The default
		mask='Province'

	#--------------------------------------------------------------------------
	# Define land surface attributes
	#--------------------------------------------------------------------------
	print('Preparing land surface attributes')
	t0=time.time()

	lsat={}

	# Land cover / land use
	vL=['LandCover_Comp1_1800','LandCover_Comp1_2019',
		'LandUse_Comp1_2019','LandUse_Comp1_2049_Scn1','LandUse_Comp1_2049_Scn2',
		'LandUseChange_Comp1_1800to2019_Year','LandUseChange_Comp1_1800to2019_Type',
		'LandUseChange_Comp1_2020to2049_Scn1_Year','LandUseChange_Comp1_2020to2049_Scn1_Type',
		'LandUseChange_Comp1_2020to2049_Scn2_Year','LandUseChange_Comp1_2020to2049_Scn2_Type']
	for v in vL:
		lsat[v]=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_' + v + '.pkl')

	# BGC zone (no longer using VRI because there are errors and gaps)
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE_GapFilled.tif')
	lsat['ID_BGCZ']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	lsat['ID_BGCZ'][lsat['ID_BGCZ']==0]=7
	uBGC=np.unique(lsat['ID_BGCZ'])

	# Region code
	lsat['Region Code']=meta['LUT']['Region']['Interior']*np.ones(meta[pNam]['Project']['N Stand'],dtype='int8')
	cd=[meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF'],
		meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH'],
		meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['MH']]
	ind=np.where(np.isin(lsat['ID_BGCZ'],cd)==True)
	lsat['Region Code'][ind]=meta['LUT']['Region']['Coast']

	# Age
	if 'Custom Age Source' in meta[pNam]['Project']:
		z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom Age Source'])
		lsat['Age']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	else:
		z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')
		lsat['Age']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

	# Site index
	if 'Custom SI Source' in meta[pNam]['Project']:
		z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom SI Source'])
		lsat['SI']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	else:
		# From BGC Zone table
		lsat['SI']=18*np.ones(meta[pNam]['Project']['N Stand'])
		u=np.unique(lsat['ID_BGCZ'])
		for iU in range(u.size):
			cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			ind=np.where(lsat['ID_BGCZ']==u[iU])[0]
			lsat['SI'][ind]=meta['Param']['BE']['ByBGCZ'][cd]['SI Calibration']

	# Natural establishment parameters (by BGC)
	# *** This gets overridden by parameters for post-wildfire natural regeneration ***
	lsat['SPH Init Natural']=1500*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
	lsat['Regen Delay Natural']=0*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
	for iU in range(uBGC.size):
		cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uBGC[iU])[0]
		ind=np.where(lsat['ID_BGCZ']==u[iU])[0]
		lsat['SPH Init Natural'][ind]=meta['Param']['BE']['ByBGCZ'][cd]['Natural Initial Tree Density']
		lsat['Regen Delay Natural'][ind]=meta['Param']['BE']['ByBGCZ'][cd]['Natural Regeneration Delay']

	# Pile burn rate parameters (by BGC)
	lsat['Pile Burn Rate']=0*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
	for iU in range(uBGC.size):
		cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],uBGC[iU])[0]
		ind=np.where(lsat['ID_BGCZ']==uBGC[iU])[0]
		lsat['Pile Burn Rate'][ind]=100*meta['Param']['BE']['ByBGCZ'][cd]['Pile Burn Rate']

	# Species
	if 'Custom Species Source' in meta[pNam]['Project']:
		# Use custom species
		z=gis.OpenGeoTiff(meta[pNam]['Project']['Custom Species Source'])
		lsat['Spc1_ID']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
		lsat['Spc2_ID']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc3_ID']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc4_ID']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc5_ID']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc1_P']=100*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc2_P']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc3_P']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc4_P']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
		lsat['Spc5_P']=9999*np.ones(lsat['Spc1_ID'].size,dtype='int16')
	else:
		# Use species from VRI
		for s in range(5):
			ss=str(s+1)
			z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_CD_' + ss + '.tif')
			lsat['Spc' + ss + '_ID']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
			z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\SPECIES_PCT_' + ss + '.tif')
			lsat['Spc' + ss + '_P']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

		# Clean species
		ind=np.where( (lsat['Spc1_ID']==0) )[0]
		lsat['Spc1_ID'][ind]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']
		lsat['Spc1_P'][ind]=100
		for s in range(5):
			ss=str(s+1)
			ind=np.where( (lsat['Spc' + ss + '_ID']==0) )[0]
			lsat['Spc' + ss + '_ID'][ind]=9999
			lsat['Spc' + ss + '_P'][ind]=9999

	# Wood density (kg/m3)
	lsat['Wood Density']=meta['Param']['BE']['Biophysical']['Density Wood']*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
	u=np.unique(lsat['Spc1_ID'])
	for iU in range(u.size):
		ind0=np.where(lsat['Spc1_ID']==u[iU])[0]
		ind1=np.where(meta['Param']['Raw']['WoodDensity']['Species CD']==cbu.lut_n2s(meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1'],u[iU])[0])[0]
		lsat['Wood Density'][ind0]=meta['Param']['Raw']['WoodDensity']['Wood Density (kg/m3)'][ind1[0]]

	# Tree density class
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\TreeDensityClass_Current.tif')
	lsat['Tree Density Class']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	lsat['Tree Density Class'][lsat['Tree Density Class']==0]=2

# 	# Operational adjustment factors
# 	if 'Custom OAF1' in meta[pNam]['Project']:
# 		lsat['OAF1']=(100*meta[pNam]['Project']['Custom OAF1'])*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
# 	else:
# 		lsat['OAF1']=(100*meta['Modules']['GYM']['OAF1 Default'])*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')
# 		u=np.unique(lsat['Tree Density Class'])
# 		for iU in range(u.size):
# 			ind0=np.where(lsat['Tree Density Class']==u[iU])[0]
# 			ind1=np.where(meta['Param']['Raw']['TreeDensityClassStatistics']['Name']==cbu.lut_n2s(meta['LUT']['Derived']['tdc'],u[iU])[0])[0]
# 			lsat['OAF1'][ind0]=100*meta['Param']['Raw']['TreeDensityClassStatistics']['OAF1'][ind1[0]]
	lsat['OAF1']=85*np.ones(meta[pNam]['Project']['N Stand'],dtype='int16')

	# Harvest year
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_Comp2_Year.tif')
	z=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	lsat['Harvest Year Comp2']=z
	
	# Harvest first year
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearFirst.tif')
	z=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	lsat['Year Harvest First']=z
	
	# Adjustments to the stocking and fertility of harvested stands
	flg=1
	if flg==1:
		iH=np.where(lsat['Year Harvest First']>0)
		lsat['OAF1'][iH]=100
		lsat['SI'][iH]=lsat['SI'][iH]+4
		#iH=np.where(lsat['Year Harvest First']==0)
		#lsat['OAF1'][iH]=45
		#lsat['SI'][iH]=lsat['SI'][iH]-1

	# Harvest probability map
	z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
	lsat['Prob Harvest (%/yr) x 1000']=gis.UpdateGridCellsize(z,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

	# Define timber harvesting land base
	lsat=DefineTHLB(meta,pNam,lsat)

	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Import disturbance and management event chronology
	#--------------------------------------------------------------------------

	print('Preparing Disturbance/management event chronology')
	t0=time.time()

	# Initiate disturbance-management event history
	dmec0=[None]*meta[pNam]['Project']['N Stand']

	# Initialize a counter for the number of events recorded for each stand
	cnt_e=np.zeros(meta[pNam]['Project']['N Stand'],'int16')

	# Initialize each stand with dummy events that weill be truncated after
	N_Events_Init=50
	for iStand in range(meta[pNam]['Project']['N Stand']):
		dmec0[iStand]={}
		dmec0[iStand]['Year']=9999*np.ones(N_Events_Init,dtype='float')
		dmec0[iStand]['DOY']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['ID Event Type']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['Event Source']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['Mortality Factor']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['Growth Factor']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['SILV_FUND_SOURCE_CODE']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['ASET']=9999*np.ones(N_Events_Init,dtype='int16')
		dmec0[iStand]['Planted SPH']=9999*np.ones(N_Events_Init,dtype='int16')
		for iSpc in range(6):
			dmec0[iStand]['PL_SPECIES_CD' + str(iSpc+1)]=9999*np.ones(N_Events_Init,dtype='int16')
			dmec0[iStand]['PL_SPECIES_PCT' + str(iSpc+1)]=9999*np.ones(N_Events_Init,dtype='int16')
			dmec0[iStand]['PL_SPECIES_GW' + str(iSpc+1)]=9999*np.ones(N_Events_Init,dtype='int16')
		if meta[pNam]['Project']['Special Attribution Method']=='NOSE':
			dmec0[iStand]['Index to Event Inciting NOSE']=9999*np.ones(N_Events_Init,dtype='int8')

	if meta[pNam]['Project']['DMEC Method']=='From Events':

		#----------------------------------------------------------------------
		# Events are taken from historical record
		#----------------------------------------------------------------------

		# Add land use change (1800-2019)
		iAffected=np.where(lsat['LandUseChange_Comp1_1800to2019_Year']>0)[0]
		for iS in iAffected:
			if lsat['LandUseChange_Comp1_1800to2019_Type'][iS]==0:
				continue
			cd=cbu.lut_n2s(meta['LUT']['Derived']['lclu_chng_comp1'],lsat['LandUseChange_Comp1_1800to2019_Type'][iS])[0]
			idType=meta['LUT']['Event'][cd]
			dmec0[iS]['Year'][cnt_e[iS]]=np.round(lsat['LandUseChange_Comp1_1800to2019_Year'][iS],decimals=2)
			dmec0[iS]['ID Event Type'][cnt_e[iS]]=idType
			dmec0[iS]['Event Source'][cnt_e[iS]]=1
			dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
			cnt_e[iS]=cnt_e[iS]+1
			# Add Pile Burn
			dmec0[iS]['Year'][cnt_e[iS]]=lsat['LandUseChange_Comp1_1800to2019_Year'][iS]+1
			dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Pile Burn']
			dmec0[iS]['Event Source'][cnt_e[iS]]=1
			dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
			cnt_e[iS]=cnt_e[iS]+1
			# Add regen failure
			dmec0[iS]['Year'][cnt_e[iS]]=lsat['LandUseChange_Comp1_1800to2019_Year'][iS]+2
			dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Regen Failure']
			dmec0[iS]['Event Source'][cnt_e[iS]]=1
			dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
			cnt_e[iS]=cnt_e[iS]+1

		# Add land use change (2020-2049)
		# *** This can affect silviculture reporting if left on! ***
		flg=0
		if flg==1:
			iAffected=np.where(lsat['LandUseChange_Comp1_2020to2049_Scn1_Year']>0)[0]
			for iS in iAffected:
				if lsat['LandUseChange_Comp1_2020to2049_Scn1_Year'][iS]==0:
					continue
				cd=cbu.lut_n2s(meta['LUT']['Derived']['lclu_chng_comp1'],lsat['LandUseChange_Comp1_2020to2049_Scn1_Type'][iS])[0]
				idType=meta['LUT']['Event'][cd]
				dmec0[iS]['Year'][cnt_e[iS]]=lsat['LandUseChange_Comp1_2020to2049_Scn1_Year'][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=idType
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1
				# Add Pile Burn
				dmec0[iS]['Year'][cnt_e[iS]]=lsat['LandUseChange_Comp1_2020to2049_Scn1_Year'][iS]+1
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Pile Burn']
				dmec0[iS]['Event Source'][cnt_e[iS]]=2
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				cnt_e[iS]=cnt_e[iS]+1
				# Add regen failure
				dmec0[iS]['Year'][cnt_e[iS]]=lsat['LandUseChange_Comp1_2020to2049_Scn1_Year'][iS]+2
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Regen Failure']
				dmec0[iS]['Event Source'][cnt_e[iS]]=2
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1

		# Add wildfire observations
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PROT_HISTORICAL_FIRE_POLYS_SP_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PROT_HISTORICAL_FIRE_POLYS_SP_DOY.pkl')
		zS=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PROT_HISTORICAL_FIRE_POLYS_SP_SevClass.pkl')
		p=meta['Param']['BE']['BurnSev']['PoC']
		for iY in range(6):
			iAffected=np.where(zY[iY]>0)[0]

			doy=gu.Clamp(zDOY[iY],1,365)

			# Where severity observations are missing, use results from regression analysis
			rn=np.random.random(zY[iY].size)
			Severity=np.zeros(rn.size)
			Severity[(rn<=p[0])]=meta['Param']['BE']['BurnSev']['M'][0]
			Severity[(rn>p[0]) & (rn<=p[1])]=meta['Param']['BE']['BurnSev']['M'][1]
			Severity[(rn>p[1]) & (rn<=p[2])]=meta['Param']['BE']['BurnSev']['M'][2]
			Severity[(rn>p[3])]=meta['Param']['BE']['BurnSev']['M'][3]

			for iS in iAffected:
				if zS[iY][iS]==meta['LUT']['Derived']['burnsev_comp1']['Unburned']:
					dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+doy[iS]/367,decimals=2)
					dmec0[iS]['DOY'][cnt_e[iS]]=doy[iS]
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Wildfire']
					dmec0[iS]['Event Source'][cnt_e[iS]]=1
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100*meta['Param']['BE']['BurnSev']['M'][0]
				elif zS[iY][iS]==meta['LUT']['Derived']['burnsev_comp1']['Low']:
					dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+doy[iS]/367,decimals=2)
					dmec0[iS]['DOY'][cnt_e[iS]]=doy[iS]
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Wildfire']
					dmec0[iS]['Event Source'][cnt_e[iS]]=1
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100*meta['Param']['BE']['BurnSev']['M'][1]
				elif zS[iY][iS]==meta['LUT']['Derived']['burnsev_comp1']['Medium']:
					dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+doy[iS]/367,decimals=2)
					dmec0[iS]['DOY'][cnt_e[iS]]=doy[iS]
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Wildfire']
					dmec0[iS]['Event Source'][cnt_e[iS]]=1
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100*meta['Param']['BE']['BurnSev']['M'][2]
				elif zS[iY][iS]==meta['LUT']['Derived']['burnsev_comp1']['High']:
					dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+doy[iS]/367,decimals=2)
					dmec0[iS]['DOY'][cnt_e[iS]]=doy[iS]
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Wildfire']
					dmec0[iS]['Event Source'][cnt_e[iS]]=1
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100*meta['Param']['BE']['BurnSev']['M'][3]
				else:
					dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+doy[iS]/367,decimals=2)
					dmec0[iS]['DOY'][cnt_e[iS]]=doy[iS]
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Wildfire']
					dmec0[iS]['Event Source'][cnt_e[iS]]=1
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100*Severity[iS]
				cnt_e[iS]=cnt_e[iS]+1

		# Add insect outbreak observations from Pest Compilation 1
		# The compilation adds insects into one year and one type file to save time.
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_InsectComp1_Year.pkl')
		zT=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_InsectComp1_Type.pkl')
		for iY in range(10):
			iAffected=np.where( (zY[iY]>0) )[0]
			for iS in iAffected:
				ind=np.where( (meta['Param']['Raw']['InsectComp1']['ID']==zT[iY][iS]) )[0][0]
				dmec0[iS]['Year'][cnt_e[iS]]=zY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event'][meta['Param']['Raw']['InsectComp1']['Insect Name'][ind]]
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=meta['Param']['Raw']['InsectComp1']['Mortality (%)'][ind]
				cnt_e[iS]=cnt_e[iS]+1

		# Add disturbances from RESULTS ATU layer
		# 'R' can relate to unburned fireguards, firefighting infrastructure -> setting as road deactivation for now
		vL=['C','D','F','I','L','P','S','W','R','B']
		namL=['Drought','Disease Root','Flooding Lightning Slides','Flooding Lightning Slides','Harvest','Mountain Pine Beetle','Harvest','Wind','Road Deactivation','Wildfire']
		iY=0 # always just one list item
		for i in range(len(vL)):

			# Don't add harvest for any project other than NOSEC (it leads to double-counting of harvest)
			if (vL[i]=='L') | (vL[i]=='S'):
				if pNam!='BCFCS_NOSEC':
					continue

			zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_DisturbanceFromATU_' + vL[i] + '_Year.pkl')
			iAffected=np.where( (zY[iY]>0) )[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=zY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event'][namL[i]]
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1

		flg=0
		if flg==1:
			# Add harvest observations (compilation 2)
			zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_Harvest_Comp2_Year.pkl')
			iAffected=np.where(zY>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=zY[iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Harvest']
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1
				# Pile burning
				if np.random.random(1)<lsat['Pile Burn Rate'][iS].astype(float)/100:
					dmec0[iS]['Year'][cnt_e[iS]]=zY[iS]+1
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Pile Burn']
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
					cnt_e[iS]=cnt_e[iS]+1

		# Add harvest observations (CC database)
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_VEG_CONSOLIDATED_CUT_BLOCKS_SP_Year.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				doy=1
				dmec0[iS]['Year'][cnt_e[iS]]=zY[iY][iS]+doy/367
				dmec0[iS]['DOY'][cnt_e[iS]]=doy
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Harvest']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1
				# Pile burning
				doy=2
				if np.random.random(1)<lsat['Pile Burn Rate'][iS].astype(float)/100:
					dmec0[iS]['Year'][cnt_e[iS]]=zY[iY][iS]+doy/367
					dmec0[iS]['DOY'][cnt_e[iS]]=doy
					dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Pile Burn']
					dmec0[iS]['Event Source'][cnt_e[iS]]=2
					dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
					cnt_e[iS]=cnt_e[iS]+1

		# Add planting observations
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_All_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_All_DOY.pkl')
		zFSC=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_All_SILV_FUND_SOURCE_CODE.pkl')
		zSPH=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_All_SPH_Planted.pkl')
		zASET=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_All_ASET.pkl')
		cd={}; pct={}; gw={}
		for iSpc in range(6):
			cd[iSpc+1]=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_SPECIES_CD' + str(iSpc+1) + '.pkl')
			pct[iSpc+1]=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_SPECIES_PCT' + str(iSpc+1) + '.pkl')
			gw[iSpc+1]=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_PL_SPECIES_GW' + str(iSpc+1) + '.pkl')
			for iY in range(6):
				# Remove zeros
				iAffected=np.where(cd[iSpc+1][iY]==0)[0]
				cd[iSpc+1][iY][iAffected]=9999
		for iY in range(6):
			iAffected=np.where( (zY[iY]>0) )[0] #  & (zASET[iY]!=meta['LUT']['Derived']['ASET']['Back-to-back Planting'])
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Planting']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				dmec0[iS]['SILV_FUND_SOURCE_CODE'][cnt_e[iS]]=zFSC[iY][iS]
				dmec0[iS]['ASET'][cnt_e[iS]]=zASET[iY][iS]
				dmec0[iS]['Planted SPH'][cnt_e[iS]]=zSPH[iY][iS]
				for iSpc in range(6):
					dmec0[iS]['PL_SPECIES_CD' + str(iSpc+1)][cnt_e[iS]]=cd[iSpc+1][iY][iS]
					dmec0[iS]['PL_SPECIES_PCT' + str(iSpc+1)][cnt_e[iS]]=pct[iSpc+1][iY][iS]
					dmec0[iS]['PL_SPECIES_GW' + str(iSpc+1)][cnt_e[iS]]=gw[iSpc+1][iY][iS]
				cnt_e[iS]=cnt_e[iS]+1

		# Add aerial fertilization observations
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_FE-CA_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_FE-CA_DOY.pkl')
		zFSC=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_FE-CA_SILV_FUND_SOURCE_CODE.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Nutrient App Aerial']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['SILV_FUND_SOURCE_CODE'][cnt_e[iS]]=zFSC[iY][iS]
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				cnt_e[iS]=cnt_e[iS]+1

		# Add knockdown
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-KD_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-KD_DOY.pkl')
		zFSC=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-KD_SILV_FUND_SOURCE_CODE.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Knockdown']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['SILV_FUND_SOURCE_CODE'][cnt_e[iS]]=zFSC[iY][iS]
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
				cnt_e[iS]=cnt_e[iS]+1

		# Add ripping
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-Rip_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-Rip_DOY.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Mechanical Site Prep']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				cnt_e[iS]=cnt_e[iS]+1

		# Add pile burn
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-PBURN_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-PBURN_DOY.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Pile Burn']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				cnt_e[iS]=cnt_e[iS]+1

		# Add prescribed burning
		zY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-BU-BROAD_Year.pkl')
		zDOY=gu.ipickle(meta['Paths']['bc1ha'] + '\\Sparse\\RGSF' + rgsf + '_Mask' + mask + '_SP-BU-BROAD_DOY.pkl')
		for iY in range(3):
			iAffected=np.where(zY[iY]>0)[0]
			for iS in iAffected:
				dmec0[iS]['Year'][cnt_e[iS]]=np.round(zY[iY][iS]+zDOY[iY][iS]/367,decimals=2)
				dmec0[iS]['DOY'][cnt_e[iS]]=zDOY[iY][iS]
				dmec0[iS]['ID Event Type'][cnt_e[iS]]=meta['LUT']['Event']['Prescribed Burn']
				dmec0[iS]['Event Source'][cnt_e[iS]]=1
				dmec0[iS]['Mortality Factor'][cnt_e[iS]]=0
				cnt_e[iS]=cnt_e[iS]+1

	elif meta[pNam]['Project']['DMEC Method']=='From Inventory Age':

		#----------------------------------------------------------------------
		# Carbon stocks in the project year are being constrained by observations
		# through input of specified age and productivity.
		#----------------------------------------------------------------------
		dType=meta['LUT']['Event']['Harvest']*np.ones(lsat['Age'].size)
		dType[np.where(np.random.random(lsat['Age'].size)<0.1)]=meta['LUT']['Event']['Wildfire']
		for iS in range(meta[pNam]['Project']['N Stand']):
			dmec0[iS]['Year'][cnt_e[iS]]=meta[pNam]['Project']['Year Project']-lsat['Age'][iS]
			dmec0[iS]['ID Event Type'][cnt_e[iS]]=dType[iS]
			dmec0[iS]['Event Source'][cnt_e[iS]]=1
			dmec0[iS]['Mortality Factor'][cnt_e[iS]]=100
			cnt_e[iS]=cnt_e[iS]+1

	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Exclude duplicate events
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Exclude duplicate events']=='On':
		 print('Removing duplicate events from dmec0')
		 t0=time.time()
		 dmec0=Exclude_Duplicate_Events(meta,pNam,dmec0)
		 print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Put events in order
	# *** Must be done before several other processing steps ***
	#--------------------------------------------------------------------------
	print('Putting events in order')
	t0=time.time()
	dmec0=PutEventsInOrder1(meta,pNam,dmec0)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Ensure that a stand-replacing disturbance precedes fertilization so that age
	#--------------------------------------------------------------------------
	# Assume previous disturbance must be missing. Assume fertilization occurs at age
	# 35. Assume previous disturbance was harvest.
	# *** Must be in order of calendar date first ***
	# *** This will not work if the previous harvest or fire had a severity < 100. ***
	# Only applies to cbrunner when fertilization is simulated from TIPSY
	if meta[pNam]['Project']['Ensure aerial fert is preceded by disturbance']=='On':
		print('Ensure stand-replacing disturbance precedes fertilization')
		t0=time.time()
		dmec0=Ensure_Fert_Preceded_By_Disturbance(meta,pNam,dmec0,lsat)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Ensure every stand has a modern disturbance
	# So that there is at least one event
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Ensure every stand has a modern disturbance']=='On':
		name_dist='Wildfire'
		severity=100
		print('Ensure every stand has a disturbance in the modern era')
		t0=time.time()
		dmec0=Ensure_Every_Stand_Has_Modern_Disturbance(meta,pNam,dmec0,name_dist,severity)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Gap-fill stands with no event history based on age from inventory input
	# So that there is at least one event
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Gap-fill DMEC with inventory age']=='On':
		print('Gap-fill DMEC with inventory age')
		t0=time.time()
		dmec0=GapFill_DMEC_WithAge(meta,pNam,dmec0,lsat)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# IDW - Western Spruce Budworm - Fix severity
	# The dmec was populated with the numeric severity ID. Mortality only occurs
	# following repeated outrbreak years.
	#--------------------------------------------------------------------------
	# if meta[pNam]['Project']['Fix severity of western spruce budworm']=='On':
	#	 print('Adjust severity of IDW')
	#	 t0=time.time()
	#	 dmec0=IDW_Fix_Severity(meta,dmec0)
	#	 print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Reduce number of growth curves by adjusting site index
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Revise SI to reduce num of growth curves']=='On':
		print('Reducing the number of growth curves by lowering the precision of site index')
		t0=time.time()
		lsat=ReduceVariationInSiteIndex(meta,pNam,lsat)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Clean species composition - TIPSY will not recognize certain codes
	# Decomissioned - now cleaned upon import
	#--------------------------------------------------------------------------
	# print('Cleaning species composition')
	# t0=time.time()
	# #meta,dmec0,vri,fcinv=Clean_Species_Composition(meta,dmec0,vri,fcinv)
	# print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Put events in order
	# *** Must be done before NOSE types are defined ***
	#--------------------------------------------------------------------------
	print('Putting events in order')
	t0=time.time()
	dmec0=PutEventsInOrder1(meta,pNam,dmec0)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Define type of non-obligation stand establishment (NOSE)
	# Do this last - you can't do this and then rearange the order of the DEMEC
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Special Attribution Method']=='NOSE':
		print('Defining types of non-obligation stand establishment')
		t0=time.time()
		meta,dmec0,lsat=Define_NOSE_ProjectType(meta,pNam,dmec0,lsat)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Put events in order
	#--------------------------------------------------------------------------
	print('Putting events in order')
	t0=time.time()
	dmec0=PutEventsInOrder2(meta,pNam,dmec0)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Remove excess size from each stand
	#--------------------------------------------------------------------------
	for iStand in range(meta[pNam]['Project']['N Stand']):
		ind=np.where(dmec0[iStand]['Year']!=9999)
		for k in dmec0[iStand].keys():
			dmec0[iStand][k]=dmec0[iStand][k][ind]

	#--------------------------------------------------------------------------
	# Exclude unidentified activities or disturbances
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Exclude unidentified events']=='On':
		print('Removing unrecognized events from DMEC')
		t0=time.time()
		dmec0=Exclude_Unidentified_Events(meta,pNam,dmec0)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Reduce growth when a wildfire occurs within 10 years of planting
	# Introduced for NOSE project
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Special Attribution Method']=='NOSE':
		for iStand in range(meta[pNam]['Project']['N Stand']):
			iPL=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) )[0]
			iWF=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Wildfire']) )[0]
			if (iPL.size>0) | (iWF.size>0):
				if iPL.size>1:
					iPL=iPL[-1]
				if iWF.size>1:
					iWF=iWF[-1]
				dT=dmec0[iStand]['Year'][iWF]-dmec0[iStand]['Year'][iPL]
				if (iWF>iPL) & (dT<10):
					dmec0[iStand]['Growth Factor'][iWF]=-50#dmec0[iStand]['Mortality Factor'][iWF]

	#--------------------------------------------------------------------------
	# Reduce growth of event that incites fill planting
	#--------------------------------------------------------------------------

# 	for iStand in range(meta[pNam]['Project']['N Stand']):
# 		iPL=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) )[0]
# 		iWF=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Wildfire']) )[0]
# 		if (iPL.size>0) | (iWF.size>0):
# 			if iPL.size>1:
# 				iPL=iPL[-1]
# 			if iWF.size>1:
# 				iWF=iWF[-1]
# 			dT=dmec0[iStand]['Year'][iWF]-dmec0[iStand]['Year'][iPL]
# 			if (iWF>iPL) & (dT<10):
# 				dmec0[iStand]['Growth Factor'][iWF]=-50#dmec0[iStand]['Mortality Factor'][iWF]

	#--------------------------------------------------------------------------
	# Expand DMEC to each scenario
	#--------------------------------------------------------------------------
	print('Expanding DMEC for each scenario')
	t0=time.time()
	dmec=[]
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dmec.append(copy.deepcopy(dmec0))
	del dmec0
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Remove events from inventory if specified by special order for each scenario
	#--------------------------------------------------------------------------
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		if 'Exclude Wildfire from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Wildfire from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Wildfire']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Harvest from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Harvest from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Harvest']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Insects from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Insects from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Mountain Pine Beetle']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Douglas-fir Beetle']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Balsam Beetle']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Spruce Beetle']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Western Spruce Budworm']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Planting from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Planting from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Planting']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Nutrient App from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Nutrient App from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Nutrient App Aerial']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Knockdown from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Knockdown from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Knockdown']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Site Prep from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Site Prep from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Mechanical Site Prep']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999
		if 'Exclude Pile Burn from Inventory' in meta[pNam]['Scenario'][iScn].keys():
			if meta[pNam]['Scenario'][iScn]['Exclude Pile Burn from Inventory']=='On':
				for iStand in range(meta[pNam]['Project']['N Stand']):
					for iE in range(dmec[iScn][iStand]['ID Event Type'].size):
						if dmec[iScn][iStand]['ID Event Type'][iE]==meta['LUT']['Event']['Pile Burn']:
							dmec[iScn][iStand]['ID Event Type'][iE]=9999

	if ('Exclude Wildfire from Inventory' in meta[pNam]['Scenario'][iScn].keys()) | ('Exclude Harvest from Inventory' in meta[pNam]['Scenario'][iScn].keys()) | ('Exclude Insects from Inventory' in meta[pNam]['Scenario'][iScn].keys()):
		print('Removing unrecognized events from DMEC')
		t0=time.time()
		dmec=Exclude_Unidentified_Events2(meta,pNam,dmec)
		print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	return meta,lsat,dmec

#%% Define type of stand establishment (salvage, knockdown, underplanting, NSR backlog)

# This function creates an index to the inciting disturbance:
# It is critical that the above steps (ensuring there is a previous
# stand-replacing disturbance) work before this will work properly.

# If it is a salvage logging project, change the inciting disturbance (harvest)
# to be "Harvest Salvage" so that it removes a higher percentage of snags. I
# checked that a combo of "Harvest Salvage" + "Pile Burn" is consistent with
# the custom harvest (with Pile Burn) used in the salvage demo (March 2021).

# The mortality corrections can be overridden by the adjustment of species-specific
# mortality (Adjust species-specific mortality='Off' to avoid this)
def Define_NOSE_ProjectType(meta,pNam,dmec0,lsat):

	# Initialize project type
	meta[pNam]['Project']['ASET']=meta['LUT']['Derived']['ASET']['Unknown']*np.ones(meta[pNam]['Project']['N Stand'],dtype='int8')
	meta[pNam]['Project']['Strata']['Project Type']['ID']=meta['LUT']['Derived']['ASET']['Unknown']*np.ones(meta[pNam]['Project']['N Stand'],dtype='int8')

	for iStand in range(meta[pNam]['Project']['N Stand']):

		# Index to stand establishment events (exclude direct seeding)
		iNOSE=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) &
				 (np.isin(dmec0[iStand]['SILV_FUND_SOURCE_CODE'],meta['Param']['Raw']['FSC']['NO List ID'])==True) )[0]

		if iNOSE.size==0:
			continue

		# Only consider the last instance
		if iNOSE.size>1:
			iNOSE=iNOSE[-1]
# 			# Focus on the last non-Fill Planting ASET unless Fill Planting is the only ASET
# 			iNotFP=np.where(dmec0[iStand]['ASET'][iNOSE]!=meta['LUT']['Derived']['ASET']['Fill Planting'])[0]
# 			if iNotFP.size>0:
# 				iNOSE=iNOSE[iNotFP[-1]]
# 			else:
# 				iNOSE=iNOSE[-1]

		# Define open spaces if an event needs to be added
		iOpen=np.where(dmec0[iStand]['Year']==9999)[0]
		if iOpen.size>1:
			iOpen=iOpen[0]

		# Extract year of focal NOSE event
		YearNOSE=dmec0[iStand]['Year'][iNOSE]

		if dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Back-to-back Planting']:

			#--------------------------------------------------------------
			# 'Back-to-back Planting'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Back-to-back Planting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Back-to-back Planting']

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Ecosystem Restoration']:

			#--------------------------------------------------------------
			# 'Ecosystem Restoration'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Ecosystem Restoration']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Ecosystem Restoration']

			# Fix the mortality for a wildfire that was added through RESULTS DENUDATIONS (BUT HAD NO BURN SEVERITY)
			ind=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Wildfire']) & (dmec0[iStand]['Mortality Factor']==0) )[0]
			if ind.size>0:
				dmec0[iStand]['Mortality Factor'][ind]=100

			iIncite=np.where( (dmec0[iStand]['Mortality Factor']==100) & (dmec0[iStand]['Year']<=YearNOSE) )[0]
			if iIncite.size==0:
				# Change mortality to 100%
				iIncite=np.where( (np.isin(dmec0[iStand]['ID Event Type'],['Wildfire','Flooding Lightning Slides','Disease Root','Harvest'])==True) & (dmec0[iStand]['Year']<=YearNOSE) )[0]
				if iIncite.size>1:
					iIncite=iIncite[-1]
				dmec0[iStand]['Mortality Factor'][iIncite]=100
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
			else:
				if iIncite.size>1:
					iIncite=iIncite[-1]
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Fill Planting']:

			#--------------------------------------------------------------
			# 'Fill Planting'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Fill Planting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Fill Planting']

			# Add month if missing
			if dmec0[iStand]['DOY'][iNOSE]==9999:
				dmec0[iStand]['DOY'][iNOSE]=152
			dmec0[iStand]['Year'][iNOSE]=dmec0[iStand]['Year'][iNOSE]+dmec0[iStand]['DOY'][iNOSE]/367

			# Some Fill Planting missing a harvest
			List=[meta['LUT']['Event']['Harvest'],meta['LUT']['Event']['Knockdown']]
			ind=np.where( (np.isin(dmec0[iStand]['ID Event Type'],List)==True) & (dmec0[iStand]['Year']<YearNOSE) )[0]
			if ind.size==0:
				yr=np.floor(dmec0[iStand]['Year'][iNOSE])-10
				doy=1
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Harvest']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				#dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

				yr=np.floor(dmec0[iStand]['Year'][iNOSE])-10
				doy=10
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Pile Burn']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=0
				#dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

				yr=np.floor(dmec0[iStand]['Year'][iNOSE])-10
				doy=152
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Planting']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=0
				#dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

			# Define inciting event
			# Just prior to fill planting, impose disturbance with 50% or 100% mortality
			Mortality=100
			GrowthFactor=-10
			if dmec0[iStand]['Planted SPH'][iNOSE]>900:
				Mortality=100
				GrowthFactor=-50
			yr=np.floor(dmec0[iStand]['Year'][iNOSE])-1
			doy=335
			dmec0[iStand]['Year'][iOpen]=yr+doy/367
			dmec0[iStand]['DOY'][iOpen]=doy
			dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Drought']
			dmec0[iStand]['Event Source'][iOpen]=2
			dmec0[iStand]['Mortality Factor'][iOpen]=Mortality
			dmec0[iStand]['Growth Factor'][iOpen]=GrowthFactor
			dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
			iOpen=iOpen+1

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Harvest and Planting']:

			#--------------------------------------------------------------
			# 'Harvest and Planting'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Harvest and Planting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Harvest and Planting']

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']:

			#--------------------------------------------------------------
			# 'Harvest and Planting NSR Backlog'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']
			iIncite=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) & (dmec0[iStand]['Year']<YearNOSE) )[0]
			if iIncite.size>1:
				iIncite=iIncite[-1]
			if iIncite.size>0:
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
			else:
				yr_lag=YearNOSE-10
				doy=1
				dmec0[iStand]['Year'][iOpen]=yr_lag+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Harvest']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Knockdown and Planting']:

			#--------------------------------------------------------------
			# 'Knockdown and Planting'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Knockdown and Planting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Knockdown and Planting']

			iKD=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Knockdown']) )[0]
			if iKD.size>1:
				iKD=iKD[-1]
			YearKD=dmec0[iStand]['Year'][iKD]

			# Define inciting event:

			# Preferentially search for mortality = 100%, then mortality > 0%
			iIncite=np.where( (dmec0[iStand]['Year']<YearKD) & (dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Harvest']) &
					(dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Knockdown']) & (dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Pile Burn']) &
					(dmec0[iStand]['Mortality Factor']==100) )[0]
			if iIncite.size==0:
				iIncite=np.where( (dmec0[iStand]['Year']<YearKD) & (dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Harvest']) &
						(dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Knockdown']) & (dmec0[iStand]['ID Event Type']!=meta['LUT']['Event']['Pile Burn']) &
						(dmec0[iStand]['Mortality Factor']>0) )[0]

			if iIncite.size>0:
				# Inciting event found from inventory
				if iIncite.size>1:
					iIncite=iIncite[-1]
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
			else:
				# Inciting event missing, introduce an inciting event
				yr=np.floor(YearKD)-1
				doy=10
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Wildfire']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

			# Introduce Pile Burn to counteract likely under-reporting
			List=[meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF'],meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBPS'],
			 meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['BWBS'],meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS']]
			if np.isin(lsat['ID_BGCZ'][iStand],List)==True:
				iPB=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Pile Burn']) & (dmec0[iStand]['Year']>=YearNOSE) )[0]
				if iPB.size==0:
					yr=np.floor(YearNOSE)
					doy=dmec0[iStand]['DOY'][iNOSE]-5 # 5 days before the planting
					dmec0[iStand]['Year'][iOpen]=yr+doy/367
					dmec0[iStand]['DOY'][iOpen]=doy
					dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Pile Burn']
					dmec0[iStand]['Event Source'][iOpen]=2
					dmec0[iStand]['Mortality Factor'][iOpen]=100
					iOpen=iOpen+1

		elif (dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']):

			#--------------------------------------------------------------
			# 'Salvage and Planting Post Beetle' or 'Salvage and Planting Post Other'
			# For non-ob salvage, assume high mortality in the natural disturbance inciting salvage
			#--------------------------------------------------------------

			meta[pNam]['Project']['ASET'][iStand]=dmec0[iStand]['ASET'][iNOSE]
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=dmec0[iStand]['ASET'][iNOSE]

			# Adjust site index to reflect the fact that NO salvage is performed on marginal areas
			lsat['SI'][iStand]=lsat['SI'][iStand]-2

			# Define harvest event
			iH=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) & (dmec0[iStand]['Year']<=YearNOSE) )[0]
			if iH.size==0:
				# If missing, introduce harvest event
				# Small number of instances
				# I think this is happing because NTEMS was considered in the ASET classification
				# but it is not being pulled into the model. Can be fixed.
				yr=YearNOSE-1
				doy=1
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Harvest']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				iH=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) & (dmec0[iStand]['Year']==yr) )[0]
				iOpen=iOpen+1

			if iH.size>1:
				iH=iH[-1]

			# Define inciting event
			iIncite=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Wildfire']) & (dmec0[iStand]['Year']<dmec0[iStand]['Year'][iH]) | \
						  (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Mountain Pine Beetle']) & (dmec0[iStand]['Year']<dmec0[iStand]['Year'][iH]) |
						  (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Spruce Beetle']) & (dmec0[iStand]['Year']<dmec0[iStand]['Year'][iH]) |
						  (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Douglas-fir Beetle']) & (dmec0[iStand]['Year']<dmec0[iStand]['Year'][iH]) |
						  (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Balsam Beetle']) & (dmec0[iStand]['Year']<dmec0[iStand]['Year'][iH]) )[0]

			if iIncite.size>1:
				iIncite=iIncite[-1]

			if iIncite.size>0:

				# Sometimes a previous event is detected, but it occurred many years ago
				dY=dmec0[iStand]['Year'][iH]-dmec0[iStand]['Year'][iIncite]
				if dY<10:
					# Recent inciting event detected
					dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
					if dmec0[iStand]['Mortality Factor'][iIncite]<75:
						# Assume government would not intervene unless at least 75% dead
						dmec0[iStand]['Mortality Factor'][iIncite]=75
				else:
					# No recent inciting event detected, introduce a recent event
					yr=dmec0[iStand]['Year'][iH]-3
					doy=1
					dmec0[iStand]['Year'][iOpen]=yr+doy/367
					dmec0[iStand]['DOY'][iOpen]=doy
					dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Mountain Pine Beetle']
					dmec0[iStand]['Event Source'][iOpen]=2
					dmec0[iStand]['Mortality Factor'][iOpen]=75
					dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
					iOpen=iOpen+1

			else:
				# No inciting event detected, introduce inciting event
				yr=dmec0[iStand]['Year'][iH]-3
				doy=1
				dmec0[iStand]['Year'][iOpen]=yr+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Mountain Pine Beetle']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=75
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']:

			#--------------------------------------------------------------
			# Straight-to-planting Post Beetles
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']
			iIncite=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Mountain Pine Beetle']) & (dmec0[iStand]['Year']<YearNOSE) )[0]
			if iIncite.size==0:
				yr_lag=YearNOSE-1
				doy=1
				dmec0[iStand]['Year'][iOpen]=yr_lag+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Mountain Pine Beetle']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1
			else:
				if iIncite.size>1:
					iIncite=iIncite[-1]
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
				if dmec0[iStand]['Mortality Factor'][iIncite]<100:
					dmec0[iStand]['Mortality Factor'][iIncite]=100

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Underplanting']:

			#--------------------------------------------------------------
			# Underplanting
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Underplanting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Underplanting']
			iIncite=np.where( (dmec0[iStand]['ID Event Type']==meta['LUT']['Event']['Wildfire']) & (dmec0[iStand]['Year']<=YearNOSE) )[0]
			if iIncite.size==0:
				# No inciting event found, introduce one in the year before
				yr_lag=YearNOSE-1
				doy=182
				dmec0[iStand]['Year'][iOpen]=yr_lag+doy/367
				dmec0[iStand]['DOY'][iOpen]=doy
				dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Wildfire']
				dmec0[iStand]['Event Source'][iOpen]=2
				dmec0[iStand]['Mortality Factor'][iOpen]=100
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
				iOpen=iOpen+1
			else:
				if iIncite.size>1:
					iIncite=iIncite[-1]
				dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iIncite
				if dmec0[iStand]['Mortality Factor'][iIncite]<100:
					dmec0[iStand]['Mortality Factor'][iIncite]=100

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Straight-to-planting Post Other']:

			#--------------------------------------------------------------
			# 'Straight-to-planting Post Other'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Other']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Straight-to-planting Post Other']

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Replanting']:

			#--------------------------------------------------------------
			# 'Replanting'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Replanting']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Replanting']

			# Introduce an event that incites the replanting
			yr=np.floor(YearNOSE)-1
			doy=335
			dmec0[iStand]['Year'][iOpen]=yr+doy/367
			dmec0[iStand]['DOY'][iOpen]=doy
			dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Drought']
			dmec0[iStand]['Event Source'][iOpen]=2
			dmec0[iStand]['Mortality Factor'][iOpen]=100
			dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
			iOpen=iOpen+1
		
		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Road Rehabilitation']:

			#--------------------------------------------------------------
			# 'Road Rehabilitation'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Road Rehabilitation']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Road Rehabilitation']

			# Introduce a land use change event that incites the road planting
			yr_lag=YearNOSE-10
			doy=1
			dmec0[iStand]['Year'][iOpen]=yr_lag+doy/367
			dmec0[iStand]['DOY'][iOpen]=doy
			dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['FL-TR']
			dmec0[iStand]['Event Source'][iOpen]=2
			dmec0[iStand]['Mortality Factor'][iOpen]=100
			dmec0[iStand]['Index to Event Inciting NOSE'][iNOSE]=iOpen
			iOpen=iOpen+1
			
			dmec0[iStand]['Year'][iOpen]=YearNOSE-9
			dmec0[iStand]['ID Event Type'][iOpen]=meta['LUT']['Event']['Regen Failure']
			dmec0[iStand]['Mortality Factor'][iOpen]=100

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Direct Seeding']:

			#--------------------------------------------------------------
			# 'Direct Seeding'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Direct Seeding']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Direct Seeding']

		elif dmec0[iStand]['ASET'][iNOSE]==meta['LUT']['Derived']['ASET']['Unknown']:

			#--------------------------------------------------------------
			# 'Unknown'
			#--------------------------------------------------------------
			meta[pNam]['Project']['ASET'][iStand]=meta['LUT']['Derived']['ASET']['Unknown']
			meta[pNam]['Project']['Strata']['Project Type']['ID'][iStand]=meta['LUT']['Derived']['ASET']['Unknown']

	return meta,dmec0,lsat

#%%
def UpdateGC_BaselineNOSE(meta,pNam,iScn,iStand,iYr,lsat,dmec,gc,iIncitingNOSE,cnt_gc):

	# Growth curve update (at the time of inciting event)
	if meta[pNam]['Project']['ASET'][iStand]!=meta['LUT']['Derived']['ASET']['Fill Planting']:
		# Don't update the growth curve if it is an event that incited fill planting
		dmec[iScn][iStand]['ID_GC'][iIncitingNOSE:]=dmec[iScn][iStand]['ID_GC'][iIncitingNOSE]+1

	gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
	gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['N']

	if meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Underplanting']:

		#----------------------------------------------------------------------
		# 'Underplanting'
		#----------------------------------------------------------------------

		# Old default (up to 2023)
		#gc[iScn][iStand]['init_density'][cnt_gc]=200
		#gc[iScn][iStand]['regen_delay'][cnt_gc]=5

		# New parameterization (2024)

		# Regen delay (based on converstation with CH, 2024)
		gc[iScn][iStand]['regen_delay'][cnt_gc]=3

		# Initial density (based on stats from FC layer)
		cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],lsat['ID_BGCZ'][iStand])[0]
		iBGC=np.where(meta['Param']['Raw']['NaturalRegenPostWildfire']['Zone / burn severity rating']==cd)[0]
		gc[iScn][iStand]['init_density'][cnt_gc]=200
		if iBGC.size>0:
			if dmec[iScn][iStand]['Mortality Factor'][iIncitingNOSE]==100*meta['Param']['BE']['BurnSev']['M'][0]:
				# Unburned
				gc[iScn][iStand]['init_density'][cnt_gc]=meta['Param']['Raw']['NaturalRegenPostWildfire']['Unburned'][iBGC]
			elif dmec[iScn][iStand]['Mortality Factor'][iIncitingNOSE]==100*meta['Param']['BE']['BurnSev']['M'][1]:
				# Low
				gc[iScn][iStand]['init_density'][cnt_gc]=meta['Param']['Raw']['NaturalRegenPostWildfire']['Low'][iBGC]
			elif dmec[iScn][iStand]['Mortality Factor'][iIncitingNOSE]==100*meta['Param']['BE']['BurnSev']['M'][1]:
				# Medium
				gc[iScn][iStand]['init_density'][cnt_gc]=meta['Param']['Raw']['NaturalRegenPostWildfire']['Medium'][iBGC]
			elif dmec[iScn][iStand]['Mortality Factor'][iIncitingNOSE]==100*meta['Param']['BE']['BurnSev']['M'][1]:
				# High
				gc[iScn][iStand]['init_density'][cnt_gc]=meta['Param']['Raw']['NaturalRegenPostWildfire']['High'][iBGC]
		else:
			gc[iScn][iStand]['init_density'][cnt_gc]=200

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Ecosystem Restoration']:

		#----------------------------------------------------------------------
		# 'Ecosystem Restoration'
		#----------------------------------------------------------------------
		if dmec[iScn][iStand]['ID Event Type'][iIncitingNOSE]==meta['LUT']['Event']['Flooding Lightning Slides']:
			gc[iScn][iStand]['init_density'][cnt_gc]=200 # Can't have zero in TIPSY
			gc[iScn][iStand]['regen_delay'][cnt_gc]=20
		elif dmec[iScn][iStand]['ID Event Type'][iIncitingNOSE]==meta['LUT']['Event']['Wildfire']:
			gc[iScn][iStand]['init_density'][cnt_gc]=200
			gc[iScn][iStand]['regen_delay'][cnt_gc]=3
		elif dmec[iScn][iStand]['ID Event Type'][iIncitingNOSE]==meta['LUT']['Event']['Harvest']:
			gc[iScn][iStand]['init_density'][cnt_gc]=1450
			gc[iScn][iStand]['regen_delay'][cnt_gc]=0

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Fill Planting']:

		#----------------------------------------------------------------------
		# 'Fill Planting'
		#----------------------------------------------------------------------
		gc[iScn][iStand]['init_density'][cnt_gc]=500
		gc[iScn][iStand]['regen_delay'][cnt_gc]=0

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Straight-to-planting Post Beetles']:

		#----------------------------------------------------------------------
		# 'Straight-to-planting Post Beetles'
		#----------------------------------------------------------------------
		gc[iScn][iStand]['init_density'][cnt_gc]=500
		gc[iScn][iStand]['regen_delay'][cnt_gc]=5

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Salvage and Planting Post Beetle']:

		#----------------------------------------------------------------------
		# 'Salvage and Planting Post Beetle'
		#----------------------------------------------------------------------
		gc[iScn][iStand]['init_density'][cnt_gc]=1400
		gc[iScn][iStand]['regen_delay'][cnt_gc]=2

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Salvage and Planting Post Other']:
		gc[iScn][iStand]['init_density'][cnt_gc]=1400
		gc[iScn][iStand]['regen_delay'][cnt_gc]=2

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Knockdown and Planting']:
		gc[iScn][iStand]['init_density'][cnt_gc]=1400
		gc[iScn][iStand]['regen_delay'][cnt_gc]=0

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']:
		gc[iScn][iStand]['init_density'][cnt_gc]=500
		gc[iScn][iStand]['regen_delay'][cnt_gc]=5

	elif meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Replanting']:
		gc[iScn][iStand]['init_density'][cnt_gc]=500
		gc[iScn][iStand]['regen_delay'][cnt_gc]=5

	gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
	gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
	gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
	gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
	gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
	gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
	gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
	gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
	gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
	gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]
	return dmec,gc

#%% Process project inputs 2
def Process2_PrepareGrowthCurves(meta,pNam,lsat,dmec):

	#--------------------------------------------------------------------------
	# Indicate which events are affected in each scenario
	# *** Make sure scenarios have been defined as "baseline" or "actual" ***
	#--------------------------------------------------------------------------

	print('Applying exclusion rules to baseline scenario.')

	for iScn in range(meta[pNam]['Project']['N Scenario']):

		for iStand in range(meta[pNam]['Project']['N Stand']):

			# Initialize indicator for each scenario
			dmec[iScn][iStand]['Scenario Affected']=np.zeros(dmec[iScn][iStand]['Year'].size,dtype='int16')

			if meta[pNam]['Project']['Special Attribution Method']=='Off':

				# All events occur in all scenarios
				for iT in range(dmec[iScn][iStand]['Year'].size):
					dmec[iScn][iStand]['Scenario Affected'][iT]=1

			elif meta[pNam]['Project']['Special Attribution Method']=='TDAF':

				meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([meta['LUT']['Event']['Direct Seeding'],
				   meta['LUT']['Event']['Knockdown'],
				   meta['LUT']['Event']['Mechanical Site Prep'],
				   meta['LUT']['Event']['Harvest'],
				   meta['LUT']['Event']['Harvest Salvage'],
				   meta['LUT']['Event']['Thinning'],
				   meta['LUT']['Event']['Aerial BTK Spray'],
				   meta['LUT']['Event']['Planting'],
				   meta['LUT']['Event']['Nutrient App Aerial'],
				   meta['LUT']['Event']['Pile Burn'],
				   meta['LUT']['Event']['Prescribed Burn']])

				for iT in range(dmec[iScn][iStand]['Year'].size):
					if np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==False:
						# All events impact
						dmec[iScn][iStand]['Scenario Affected'][iT]=1
					elif np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True:
						if np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==False:
							# Events not in the list are added
							dmec[iScn][iStand]['Scenario Affected'][iT]=1
					else:
						pass

			elif meta[pNam]['Project']['Special Attribution Method']=='LUC':

				meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([meta['LUT']['Event']['FL-CL'],
				   meta['LUT']['Event']['FL-PA'],
				   meta['LUT']['Event']['FL-RC'],
				   meta['LUT']['Event']['FL-TR'],
				   meta['LUT']['Event']['FL-EM'],
				   meta['LUT']['Event']['Pile Burn'],
				   meta['LUT']['Event']['Regen Failure']])

				for iT in range(dmec[iScn][iStand]['Year'].size):
					if np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==False:
						# All events impact
						dmec[iScn][iStand]['Scenario Affected'][iT]=1
					elif np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True:
						if np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==False:
							# Events not in the list are added
							dmec[iScn][iStand]['Scenario Affected'][iT]=1
					else:
						pass

			elif meta[pNam]['Project']['Special Attribution Method']=='BAU':
				meta[pNam]['Project']['Activities To Exclude From Baseline']=np.array([
				   meta['LUT']['Event']['Direct Seeding'],
				   meta['LUT']['Event']['Knockdown'],
				   meta['LUT']['Event']['Mechanical Site Prep'],
				   meta['LUT']['Event']['Harvest'],
				   meta['LUT']['Event']['Thinning'],
				   meta['LUT']['Event']['Aerial BTK Spray'],
				   meta['LUT']['Event']['Planting'],
				   meta['LUT']['Event']['Nutrient App Aerial'],
				   meta['LUT']['Event']['Pile Burn'],
				   meta['LUT']['Event']['Prescribed Burn']])
				meta[pNam]['Project']['Activities To Exclude From Actual']=np.array([
				   meta['LUT']['Event']['Mechanical Site Prep'],
				   meta['LUT']['Event']['Thinning'],
				   meta['LUT']['Event']['Aerial BTK Spray'],
				   meta['LUT']['Event']['Nutrient App Aerial'],
				   meta['LUT']['Event']['Prescribed Burn']])

				for iT in range(dmec[iScn][iStand]['Year'].size):
					# Non-obligation status
					StatusNO=np.isin(dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][iT],meta['Param']['Raw']['FSC']['NO List ID'])

					if (np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Actual'])==False) & (StatusNO==False):
						# All but excluded events
						dmec[iScn][iStand]['Scenario Affected'][iT]=1

					elif (np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True) & (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],meta[pNam]['Project']['Activities To Exclude From Baseline'])==False):
						# Events not in the list are added
						dmec[iScn][iStand]['Scenario Affected'][iT]=1
					else:
						pass

			elif meta[pNam]['Project']['Special Attribution Method']=='NOSE':

				#--------------------------------------------------------------
				# Non-obligation stand establishment
				#--------------------------------------------------------------

				# List of activities that will be excluded from baseline
				if meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Harvest and Planting NSR Backlog']:
					excludeL=['Planting','Direct Seeding','Knockdown','Mechanical Site Prep']
				else:
					excludeL=['Planting','Direct Seeding','Harvest','Knockdown','Pile Burn','Mechanical Site Prep']

				meta[pNam]['Project']['Activities To Exclude From Baseline']=np.zeros(len(excludeL))
				for i in range(len(excludeL)):
					meta[pNam]['Project']['Activities To Exclude From Baseline'][i]=meta['LUT']['Event'][excludeL[i]]

				# Index to stand establishment events
				#StatusNO=np.isin(dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'],meta['Param']['BE']['FSC']['NO List ID'])
				iNOSE=np.where( (dmec[iScn][iStand]['Index to Event Inciting NOSE']!=9999) )[0]

				if iNOSE.size>0:
					# If multiple planting events, focus on the last instance
					iNOSE=iNOSE[-1]

				# Index to inciting event
				iIncitingEvent=dmec[iScn][iStand]['Index to Event Inciting NOSE'][iNOSE]

				for iT in range(dmec[iScn][iStand]['Year'].size):
					if np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True:

						# Baseline scenario, only some events are included
						if (dmec[iScn][iStand]['Year'][iT]<=dmec[iScn][iStand]['Year'][iIncitingEvent]):
							dmec[iScn][iStand]['Scenario Affected'][iT]=1

						# Post-treatment natural disturbances need to be included in baseline
						if dmec[iScn][iStand]['Year'][iT]>=dmec[iScn][iStand]['Year'][iNOSE]:
							List=[meta['LUT']['Event']['Wildfire'],meta['LUT']['Event']['Mountain Pine Beetle'],meta['LUT']['Event']['Spruce Beetle'],meta['LUT']['Event']['Douglas-fir Beetle'],meta['LUT']['Event']['Balsam Beetle']]
							if (np.isin(dmec[iScn][iStand]['ID Event Type'][iT],List)==True):
								dmec[iScn][iStand]['Scenario Affected'][iT]=1
					else:
						# Project scenario, everything occurs in all scenarios
						dmec[iScn][iStand]['Scenario Affected'][iT]=1

			elif meta[pNam]['Project']['Special Attribution Method']=='Nutrient Management':

				#--------------------------------------------------------------
				# Nutrient management
				#--------------------------------------------------------------

				for iT in range(dmec[iScn][iStand]['Year'].size):
					if (np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True):
						if meta[pNam]['Project']['Code Project']=='BCFCS_NMF':
							if dmec[iScn][iStand]['ID Event Type'][iT]==meta['LUT']['Event']['Nutrient App Aerial']:
								if dmec[iScn][iStand]['Year'][iT]>=meta[pNam]['Project']['Year Project']:
									dmec[iScn][iStand]['Scenario Affected'][iT]=1
							else:
								dmec[iScn][iStand]['Scenario Affected'][iT]=1
						else:
							dmec[iScn][iStand]['Scenario Affected'][iT]=1
					else:
						if dmec[iScn][iStand]['ID Event Type'][iT]!=meta['LUT']['Event']['Nutrient App Aerial']:
							dmec[iScn][iStand]['Scenario Affected'][iT]=1

	#--------------------------------------------------------------------------
	# Parameterize growth curves
	#--------------------------------------------------------------------------

	print('Preparing growth curves')
	t0=time.time()

	# Initialize list
	gc=[None]*meta[pNam]['Project']['N Scenario']

	for iScn in range(meta[pNam]['Project']['N Scenario']):

		gc[iScn]=[None]*meta[pNam]['Project']['N Stand']

		for iStand in range(meta[pNam]['Project']['N Stand']):

			# Index to growth curve
			cnt_gc=0

			# Initialize growth curve identifiers in DMEC
			dmec[iScn][iStand]['ID_GC']=1*np.ones(dmec[iScn][iStand]['Year'].size)

			# Initialize growth curve info
			gc[iScn][iStand]={}
			for key in meta['Modules']['GYM']['GC_Variable_List']:
				if np.isin(key,['ID_Stand','ID_Scn','ID_GC','s1','p1','s2','p2','s3','p3','s4','p4','s5','p5','s6','p6','init_density','regen_delay','regeneration_method','bec_zone','FIZ'])==True:
					if key=='ID_Stand':
						# For big projects, this needs to be 32 bits
						gc[iScn][iStand][key]=9999*np.ones(12,dtype='int32')
					else:
						gc[iScn][iStand][key]=9999*np.ones(12,dtype='int16')
				else:
					gc[iScn][iStand][key]=9999*np.ones(12)

			if lsat['Region Code'][iStand]==meta['LUT']['Region']['Coast']:
				fiz=meta['LUT']['TIPSY']['FIZ']['C']
			else:
				fiz=meta['LUT']['TIPSY']['FIZ']['I']

			#--------------------------------------------------------------------------
			# Add pre-contact (spinup) growth curve
			#--------------------------------------------------------------------------
			gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
			gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
			gc[iScn][iStand]['ID_GC'][cnt_gc]=1
			gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['N']
			gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
			gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
			gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]
			gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
			gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
			gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
			gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
			gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
			gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
			gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
			gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]
			gc[iScn][iStand]['init_density'][cnt_gc]=lsat['SPH Init Natural'][iStand]
			gc[iScn][iStand]['regen_delay'][cnt_gc]=lsat['Regen Delay Natural'][iStand]
			gc[iScn][iStand]['oaf1'][cnt_gc]=lsat['OAF1'][iStand].astype(float)/100
			gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
			gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
			gc[iScn][iStand]['FIZ'][cnt_gc]=fiz
			cnt_gc=cnt_gc+1

			#--------------------------------------------------------------------------
			# Add events from disturbance/management event history
			#--------------------------------------------------------------------------
			for iYr in range(dmec[iScn][iStand]['Year'].size):

				# Calculate planting density
				PlantingDensity=int(dmec[iScn][iStand]['Planted SPH'][iYr])
				if lsat['Region Code'][iStand]==meta['LUT']['Region']['Interior']:
					PlantingDensity=np.minimum(2400,np.maximum(1400,PlantingDensity))
				else:
					PlantingDensity=np.minimum(2400,np.maximum(900,PlantingDensity))
				PlantingDensity=int(PlantingDensity)

				# Create a flag that indicates whether there are back-to-back planting
				# Back to back planting I think occurs in some cases because they go back
				# and add a bit. You can tell by looking at the treatment area - if the
				# second planting treatment area is tiny compared to the first, you could
				# ignore it I guess. Certainly not ideal, but I don't see a work around.
				# We are going to ignore the second planting for now.
				Flag_PlantingBackToBack=0
				if iYr>0:
					if (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):
						Yr=dmec[iScn][iStand]['Year'][iYr]
						indPrevPL=np.where( (dmec[iScn][iStand]['Year']==Yr-1) & (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Planting']) )[0]
						if indPrevPL.size>0:
							Flag_PlantingBackToBack=1

				# Index to previous disturbance for fertilization
				#IndPrevDistForFert=int(dmec[iStand]['IndPrevDistForFert'][iYr])

				# Non-obligation status
				StatusNO=np.isin(dmec[iScn][iStand]['SILV_FUND_SOURCE_CODE'][iYr],meta['Param']['Raw']['FSC']['NO List ID'])

				if (meta[pNam]['Project']['Special Attribution Method']=='Off') | (meta[pNam]['Project']['Special Attribution Method']=='TDAF') | (meta[pNam]['Project']['Special Attribution Method']=='Nutrient Management'):

					#--------------------------------------------------------------
					# Off or TDAF or Nutrient Management
					#--------------------------------------------------------------
					
					if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):

						dmec[iScn][iStand]['ID_GC'][iYr:]=np.max(gc[iScn][iStand]['ID_GC'][gc[iScn][iStand]['ID_GC']<9999])+1

						gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
						gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
						gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
						gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
						gc[iScn][iStand]['init_density'][cnt_gc]=PlantingDensity
						gc[iScn][iStand]['regen_delay'][cnt_gc]=0
						gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]

						if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=9999) & (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]!=9999):

							# *** Adjust site index if it is energy production ***
							if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']):
								gc[iScn][iStand]['i1'][cnt_gc]=30
								gc[iScn][iStand]['init_density'][cnt_gc]=2000

							# Using planting info if it exists
							gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
							gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
							gc[iScn][iStand]['gain1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr]
							gc[iScn][iStand]['selage1'][cnt_gc]=10
							gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
							gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
							gc[iScn][iStand]['gain2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr]
							gc[iScn][iStand]['selage2'][cnt_gc]=10
							gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
							gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
							gc[iScn][iStand]['gain3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr]
							gc[iScn][iStand]['selage3'][cnt_gc]=10
							gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
							gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
							gc[iScn][iStand]['gain4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr]
							gc[iScn][iStand]['selage4'][cnt_gc]=10
							gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
							gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
							gc[iScn][iStand]['gain5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr]
							gc[iScn][iStand]['selage5'][cnt_gc]=10
						else:
							# Otherwise assume best-available inventory spc. comp.
							gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
							gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
							gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
							gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
							gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
							gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
							gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
							gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
							gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
							gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]
						gc[iScn][iStand]['oaf1'][cnt_gc]=lsat['OAF1'][iStand].astype(float)/100
						gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
						gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
						gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

						# Update counter
						cnt_gc=cnt_gc+1

				elif (meta[pNam]['Project']['Special Attribution Method']=='BAU'):
					
					#--------------------------------------------------------------
					# Harvesting
					# Exclude incremental improvements in silviculture, including:
					#  - genetic gains
					#--------------------------------------------------------------
					
					if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']):

						dmec[iScn][iStand]['ID_GC'][iYr:]=np.max(gc[iScn][iStand]['ID_GC'][gc[iScn][iStand]['ID_GC']<9999])+1

						gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
						gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
						gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
						gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
						gc[iScn][iStand]['init_density'][cnt_gc]=PlantingDensity
						gc[iScn][iStand]['regen_delay'][cnt_gc]=0
						gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]

						if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=9999) & (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]!=9999):
							# Using planting info if it exists
							gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
							gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
							gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
							gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
							gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
							gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
							gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
							gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
							gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
							gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
						else:
							# Otherwise assume best-available inventory spc. comp.
							gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
							gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
							gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
							gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
							gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
							gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
							gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
							gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
							gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
							gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]
						gc[iScn][iStand]['oaf1'][cnt_gc]=lsat['OAF1'][iStand].astype(float)/100
						gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
						gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
						gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

						# Update counter
						cnt_gc=cnt_gc+1

				elif (meta[pNam]['Project']['Special Attribution Method']=='NOSE'):

					#--------------------------------------------------------------
					# Non-obligation stand establishment (NOSE)
					#--------------------------------------------------------------

					# Index to event that incited NO stand establishment
					iIncitingNOSE=dmec[iScn][iStand]['Index to Event Inciting NOSE'][iYr]

					#----------------------------------------------------------------------
					# Planting (non-obligation)
					#----------------------------------------------------------------------

					if (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']) & (StatusNO==True) & (iIncitingNOSE!=9999):

						#------------------------------------------------------
						# Attributes common to baseline and project scenarios
						#------------------------------------------------------

						gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
						gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
						gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]
						gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
						gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
						gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
						gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

						# Whitebark pine restoration project
						if dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PA']:
							gc[iScn][iStand]['i1'][cnt_gc]=14
							gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

						# Aspen planted in drybelt
						if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']) & (lsat['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']):
							#print('Working Drybelt Aspen')
							gc[iScn][iStand]['i1'][cnt_gc]=16
							gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

						# Mt Meager slide
						if dmec[iScn][iStand]['ID Event Type'][iIncitingNOSE]==meta['LUT']['Event']['Flooding Lightning Slides']:
							#print('Working Slides')
							gc[iScn][iStand]['i1'][cnt_gc]=18
							gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

						#----------------------------------------------------------------------
						# Scneario-specific attributes
						#----------------------------------------------------------------------

						if np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True:

							# Baseline scenarios:
							dmec,gc=UpdateGC_BaselineNOSE(meta,pNam,iScn,iStand,iYr,lsat,dmec,gc,iIncitingNOSE,cnt_gc)

							# Update counter
							cnt_gc=cnt_gc+1

						else:

							# Project scenarios:

							# When there is a time gap between the inciting event and the NO planting,
							# the project DMEC must adopt the baseline curve prior to planting
							if meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Ecosystem Restoration']:
								dmec,gc=UpdateGC_BaselineNOSE(meta,pNam,iScn,iStand,iYr,lsat,dmec,gc,iIncitingNOSE,cnt_gc)

								# Update counter
								cnt_gc=cnt_gc+1

								gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
								gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
								gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]
								gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
								gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
								gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
								gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

								# Super awkward that the special orders (below) need to be repeated!

								# Whitebark pine restoration project
								if dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PA']:
									gc[iScn][iStand]['i1'][cnt_gc]=14
									gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

								# Aspen planted in drybelt
								if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']) & (lsat['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['IDF']):
									#print('Working Drybelt Aspen')
									gc[iScn][iStand]['i1'][cnt_gc]=16
									gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

								# Mt Meager slide
								if dmec[iScn][iStand]['ID Event Type'][iIncitingNOSE]==meta['LUT']['Event']['Flooding Lightning Slides']:
									#print('Working Slides')
									gc[iScn][iStand]['i1'][cnt_gc]=18
									gc[iScn][iStand]['oaf1'][cnt_gc]=0.5

							# Growth curve update: Project scenarios with planting at iYr
							dmec[iScn][iStand]['ID_GC'][iYr:]=np.max(gc[iScn][iStand]['ID_GC'][gc[iScn][iStand]['ID_GC']<9999])+1
							gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
							gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
							gc[iScn][iStand]['init_density'][cnt_gc]=PlantingDensity

							# Make sure a Fill Plant leads to a fully stocked stand
							if meta[pNam]['Project']['ASET'][iStand]==meta['LUT']['Derived']['ASET']['Fill Planting']:
								gc[iScn][iStand]['init_density'][cnt_gc]=1850

							gc[iScn][iStand]['regen_delay'][cnt_gc]=0
							if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=9999) & (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]!=9999):
								# Using planting layer information
								gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
								gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
								gc[iScn][iStand]['gain1'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr])
								gc[iScn][iStand]['selage1'][cnt_gc]=10
								gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
								gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
								gc[iScn][iStand]['gain2'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr])
								gc[iScn][iStand]['selage2'][cnt_gc]=10
								gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
								gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
								gc[iScn][iStand]['gain3'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr])
								gc[iScn][iStand]['selage3'][cnt_gc]=10
								gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
								gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
								gc[iScn][iStand]['gain4'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr])
								gc[iScn][iStand]['selage4'][cnt_gc]=10
								gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
								gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
								gc[iScn][iStand]['gain5'][cnt_gc]=int(dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr])
								gc[iScn][iStand]['selage5'][cnt_gc]=10
							else:
								# Planting info not available, using inventory
								gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
								gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
								gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
								gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
								gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
								gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
								gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
								gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
								gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
								gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]

							# Update counter
							cnt_gc=cnt_gc+1

					#----------------------------------------------------------------------
					# Planting (obligation)
					#----------------------------------------------------------------------
					
					if (dmec[iScn][iStand]['Scenario Affected'][iYr]==1) & (dmec[iScn][iStand]['ID Event Type'][iYr]==meta['LUT']['Event']['Planting']) & (Flag_PlantingBackToBack==0) & (StatusNO==False):

						dmec[iScn][iStand]['ID_GC'][iYr:]=np.max(gc[iScn][iStand]['ID_GC'][gc[iScn][iStand]['ID_GC']<9999])+1

						gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
						gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
						gc[iScn][iStand]['ID_GC'][cnt_gc]=dmec[iScn][iStand]['ID_GC'][iYr]
						gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
						gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]
						gc[iScn][iStand]['init_density'][cnt_gc]=PlantingDensity
						gc[iScn][iStand]['regen_delay'][cnt_gc]=0
						gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
						gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
						gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
						gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

						if (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=9999) & (dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]!=0) & (dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]!=9999):
							gc[iScn][iStand]['s1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD1'][iYr]
							gc[iScn][iStand]['p1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT1'][iYr]
							gc[iScn][iStand]['gain1'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW1'][iYr]
							gc[iScn][iStand]['selage1'][cnt_gc]=10
							gc[iScn][iStand]['s2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD2'][iYr]
							gc[iScn][iStand]['p2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT2'][iYr]
							gc[iScn][iStand]['gain2'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW2'][iYr]
							gc[iScn][iStand]['selage2'][cnt_gc]=10
							gc[iScn][iStand]['s3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD3'][iYr]
							gc[iScn][iStand]['p3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT3'][iYr]
							gc[iScn][iStand]['gain3'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW3'][iYr]
							gc[iScn][iStand]['selage3'][cnt_gc]=10
							gc[iScn][iStand]['s4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD4'][iYr]
							gc[iScn][iStand]['p4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT4'][iYr]
							gc[iScn][iStand]['gain4'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW4'][iYr]
							gc[iScn][iStand]['selage4'][cnt_gc]=10
							gc[iScn][iStand]['s5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_CD5'][iYr]
							gc[iScn][iStand]['p5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_PCT5'][iYr]
							gc[iScn][iStand]['gain5'][cnt_gc]=dmec[iScn][iStand]['PL_SPECIES_GW5'][iYr]
							gc[iScn][iStand]['selage5'][cnt_gc]=10
						else:
							gc[iScn][iStand]['s1'][cnt_gc]=lsat['Spc1_ID'][iStand]
							gc[iScn][iStand]['p1'][cnt_gc]=lsat['Spc1_P'][iStand]
							gc[iScn][iStand]['s2'][cnt_gc]=lsat['Spc2_ID'][iStand]
							gc[iScn][iStand]['p2'][cnt_gc]=lsat['Spc2_P'][iStand]
							gc[iScn][iStand]['s3'][cnt_gc]=lsat['Spc3_ID'][iStand]
							gc[iScn][iStand]['p3'][cnt_gc]=lsat['Spc3_P'][iStand]
							gc[iScn][iStand]['s4'][cnt_gc]=lsat['Spc4_ID'][iStand]
							gc[iScn][iStand]['p4'][cnt_gc]=lsat['Spc4_P'][iStand]
							gc[iScn][iStand]['s5'][cnt_gc]=lsat['Spc5_ID'][iStand]
							gc[iScn][iStand]['p5'][cnt_gc]=lsat['Spc5_P'][iStand]

						# Update counter
						cnt_gc=cnt_gc+1

			# Add a planted stand if missing:
			# Some stands may have no recorded historical management, but then they
			# may be harvested on-the-fly. Add a second growth curve for the on-the-fly
			# planting.
			if cnt_gc==1:
				gc[iScn][iStand]['ID_Stand'][cnt_gc]=iStand
				gc[iScn][iStand]['ID_Scn'][cnt_gc]=iScn
				gc[iScn][iStand]['ID_GC'][cnt_gc]=np.max(gc[iScn][iStand]['ID_GC'][gc[iScn][iStand]['ID_GC']<9999])+1
				gc[iScn][iStand]['regeneration_method'][cnt_gc]=meta['LUT']['TIPSY']['regeneration_method']['P']
				gc[iScn][iStand]['i1'][cnt_gc]=lsat['SI'][iStand]
				gc[iScn][iStand]['init_density'][cnt_gc]=int(1500)
				gc[iScn][iStand]['regen_delay'][cnt_gc]=0
				gc[iScn][iStand]['oaf1'][cnt_gc]=meta['Modules']['GYM']['OAF1 Default']
				gc[iScn][iStand]['oaf2'][cnt_gc]=meta['Modules']['GYM']['OAF2 Default']
				gc[iScn][iStand]['bec_zone'][cnt_gc]=lsat['ID_BGCZ'][iStand]
				gc[iScn][iStand]['FIZ'][cnt_gc]=fiz

				# The natural species composition may not be realistic. Use regional
				# default planting composition

				gain=15

				if (lsat['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) | (lsat['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CDF']):
					# Coastal
					gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['FDC']
					gc[iScn][iStand]['p1'][cnt_gc]=70
					gc[iScn][iStand]['gain1'][cnt_gc]=gain
					gc[iScn][iStand]['selage1'][cnt_gc]=10
					gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['CW']
					gc[iScn][iStand]['p2'][cnt_gc]=15
					gc[iScn][iStand]['gain2'][cnt_gc]=gain
					gc[iScn][iStand]['selage2'][cnt_gc]=10
					gc[iScn][iStand]['s3'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['BA']
					gc[iScn][iStand]['p3'][cnt_gc]=15
					gc[iScn][iStand]['gain3'][cnt_gc]=gain
					gc[iScn][iStand]['selage3'][cnt_gc]=10
				elif (lsat['ID_BGCZ'][iStand]==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['ICH']):
					# Interior wetbelt
					gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['CW']
					gc[iScn][iStand]['p1'][cnt_gc]=70
					gc[iScn][iStand]['gain1'][cnt_gc]=gain
					gc[iScn][iStand]['selage1'][cnt_gc]=10
					gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['HW']
					gc[iScn][iStand]['p2'][cnt_gc]=15
					gc[iScn][iStand]['gain2'][cnt_gc]=gain
					gc[iScn][iStand]['selage2'][cnt_gc]=10
				else:
					# Interior
					gc[iScn][iStand]['s1'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']
					gc[iScn][iStand]['p1'][cnt_gc]=60
					gc[iScn][iStand]['gain1'][cnt_gc]=gain
					gc[iScn][iStand]['selage1'][cnt_gc]=10
					gc[iScn][iStand]['s2'][cnt_gc]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SW']
					gc[iScn][iStand]['p2'][cnt_gc]=40
					gc[iScn][iStand]['gain2'][cnt_gc]=gain
					gc[iScn][iStand]['selage2'][cnt_gc]=10

				if gc[iScn][iStand]['s1'][cnt_gc]==0:
					print(iScn)
					print(iStand)
					print(iYr)

	# Get rid of rows with no info
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		for iStand in range(meta[pNam]['Project']['N Stand']):
			# Don't use ID_Stand, because there could be a legit ID_Stand=9999
			ind=np.where(gc[iScn][iStand]['ID_GC']!=9999)[0]
			for key in meta['Modules']['GYM']['GC_Variable_List']:
				gc[iScn][iStand][key]=gc[iScn][iStand][key][ind]

	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Adjust mortality factors that only affect specific tree species
	# *** MOVED TO CBRUN_ANNPROC - CONSIDER DELETING ***
	#--------------------------------------------------------------------------
	# if meta[pNam]['Project']['Adjust species-specific mortality']=='On':
	#	 print('Adjusting mortality based on species-specific pests')
	#	 t0=time.time()
	#	 dmec=AdjustSpeciesSpecificMortality(meta,pNam,dmec,gc,meta[pNam]['Project']['Actual Indices'][0])
	#	 print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Extract a set of unique growth curves
	# Decompose the full set of stands into a subset of unique stand types.
	# Exclude the first three columns, as they are all different.
	#--------------------------------------------------------------------------
	print('Extracting unique growth curves')
	t0=time.time()
	ugc=ExtractUniqueGrowthCurves(meta,pNam,gc)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Export to BatchTIPSY parameter spreadsheet
	#--------------------------------------------------------------------------
	print('Exporting BatchTIPSY parameters to spreadsheet')
	t0=time.time()
	cbu.Write_BatchTIPSY_Input_Spreadsheet(meta,pNam,ugc)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Populate BatchTIPSY.exe input variable (.dat) file
	#--------------------------------------------------------------------------
	print('Creating BatchTIPSY.exe input varialbe (.dat) file')
	t0=time.time()
	cbu.Write_BatchTIPSY_Input_File(meta,pNam)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	return gc,ugc,lsat,dmec

#%% Process project inputs 3
def Process3_PrepInputsByBatch(meta,pNam,lsat,dmec,gc,ugc):

	#--------------------------------------------------------------------------
	# Prepare land surface attributes
	#--------------------------------------------------------------------------
	print('Preparing inventory input files')
	t0=time.time()
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		# Loop through batches, saving inventory to file
		for iBat in range(meta[pNam]['Project']['N Batch']):
			
			# Index to batch
			indBat=cbu.IndexToBatch(meta[pNam],iBat)
			N_StandsInBatch=len(indBat)
			
			# Initilize
			lsat1={}
			
			# Populate with arrays
			for k in lsat.keys():
				# Exclude nested dictionaries
				if (k=='THLB'):
					continue
				lsat1[k]=np.zeros((1,N_StandsInBatch),dtype='int16')
				try:
					lsat1[k][0,:]=lsat[k][indBat]
				except:
					print(k)
					print(lsat[k].shape)

			# Timber harvesting landbase (1=yes, 0=no)
			sLCLU=meta[pNam]['Scenario'][iScn]['Land Cover/Land Use Scenario']
			lsat1['THLB']=lsat['THLB'][sLCLU][:,indBat]

			# Temperature will be updated automatically
			lsat1['MAT']=4*np.ones((1,N_StandsInBatch))

			# Probability of harvest (%/yr) from spatial map
			lsat1['Prob Harvest (%/yr) x 1000']=lsat['Prob Harvest (%/yr) x 1000'][indBat]

			# Wood density (kg/m3)
			lsat1['Wood Density']=lsat['Wood Density'][indBat]

			# Sawtooth species-region samples
			if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
				lsat1['Srs1_ID']=meta['LUT']['Spc'][meta[pNam]['Scenario'][iScn]['SRS1_CD']]*np.ones((1,N_StandsInBatch),dtype='int16')
			else:
				lsat1['Srs1_ID']=9999*np.ones((1,N_StandsInBatch),dtype='int16')

			# Save
			gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl',lsat1)

	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Prepare disturbance/management event chronology
	#--------------------------------------------------------------------------

	print('Preparing DMEC input files')
	t0=time.time()
	for iEns in range(meta[pNam]['Project']['N Ensemble']):
		for iScn in range(meta[pNam]['Project']['N Scenario']):
			for iBat in range(meta[pNam]['Project']['N Batch']):

				# Index to batch
				indBat=cbu.IndexToBatch(meta[pNam],iBat)

				# Initialize dictionary
				ec={}
				ec['ID Event Type']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
				ec['Mortality Factor']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
				ec['Growth Factor']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
				ec['ID Growth Curve']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
				ec['ASET']=np.zeros((meta[pNam]['Year'].size,indBat.size,meta['Core']['Max Events Per Year']),dtype='int16')
				tv=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

				#----------------------------------------------------------
				# Spinup with constant return interval
				#----------------------------------------------------------

				if (meta[pNam]['Project']['Spinup Status']=='On'):

					lsat0=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')

					for iS in range(indBat.size):

						# Index to stand
						iStandFull=indBat[iS]

						# Pre-industrial disturbance interval

						# Old: regional
						#if lsat0['Region Code'][0,iS]==meta['LUT']['Region']['Coast']:
						#	ivl_pi=300
						#else:
						#	ivl_pi=125

						if meta[pNam]['Project']['Return Interval Source']=='Custom':
							# From custom input
							ivl_pi=meta[pNam]['Project']['Custom Return Interval']
						elif meta[pNam]['Project']['Return Interval Source']=='BGC Zone':
							# BGC Zone values
							cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],lsat0['ID_BGCZ'][0,iS])[0]
							ivl_pi=meta['Param']['BE']['ByBGCZ'][cd]['Disturbance Return Interval']
						else:
							print('Spin-up return interval source incorrect.')
							
						# Timing of transition between pre-industrial and modern periods
						try:
							YearRef=dmec[iScn][iStandFull]['Year'][0]
						except:
							YearRef=np.random.randint(1775,high=1920,size=1,dtype=int)
						AgeRef=ivl_pi

						YrRegCyc=np.arange(YearRef-AgeRef-100*ivl_pi,YearRef-AgeRef+ivl_pi,ivl_pi)
						Year=YrRegCyc[np.where(YrRegCyc>=meta[pNam]['Year'][0])[0]]
						ID_Type=meta['LUT']['Event']['Wildfire']*np.ones(Year.size)
						MortF=100*np.ones(Year.size)
						GrowthF=0*np.ones(Year.size)
						ID_GrowthCurve=1*np.ones(Year.size)
						ASET=np.zeros(Year.size)
						ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve,ASET)

				#------------------------------------------------------------------
				# Add modern era events
				#------------------------------------------------------------------

				for iS in range(indBat.size):

					# Index to stand
					iStandFull=indBat[iS]

					ind=np.where(dmec[iScn][iStandFull]['Scenario Affected']==1)[0]

					if ind.size>0:
						ID_Type=dmec[iScn][iStandFull]['ID Event Type'][ind]

						# Remove historical wildfire if it is being simulated
						if meta[pNam]['Scenario'][iScn]['Wildfire Sim Obs Status']=='On':
							ind2=np.where(ID_Type!=meta['LUT']['Event']['Wildfire'])[0]
							ind=ind[ind2]
							ID_Type=ID_Type[ind2]

						# *** SPECIAL ORDER ***
						if (meta[pNam]['Project']['Code Project']=='BCFCS_Wildfire23') & (iScn==0):
							iToss=np.where( (dmec[iScn][iStandFull]['Year'][ind]==2023) & (ID_Type==meta['LUT']['Event']['Wildfire']) )[0]
							dmec[iScn][iStandFull]['Year'][ind[iToss]]=-1
						# *** SPECIAL ORDER ***

						Year=dmec[iScn][iStandFull]['Year'][ind]
						MortF=dmec[iScn][iStandFull]['Mortality Factor'][ind]
						GrowthF=dmec[iScn][iStandFull]['Growth Factor'][ind]
						ID_GrowthCurve=dmec[iScn][iStandFull]['ID_GC'][ind]
						ASET=dmec[iScn][iStandFull]['ASET'][ind]
						ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve,ASET)

				#------------------------------------------------------------------
				# Add future scheduled NOSE
				#------------------------------------------------------------------

				if 'NOSE Future' in meta[pNam]['Project']:

					c,ia,ib=np.intersect1d(indBat,meta[pNam]['Project']['NOSE Future']['Stand Index'],return_indices=True)

					if ia.size>0:

						for iP in range(ia.size):

							YearIncite=meta[pNam]['Project']['NOSE Future']['Year'][ib[iP]]
							iT=np.where(tv==YearIncite)[0]
							iAvailable=np.where(ec['ID Event Type'][iT,ia[iP],:]==0)[0]
							ec['ID Event Type'][iT,ia[iP],iAvailable]=meta['LUT']['Event']['GC Switch']
							ec['Mortality Factor'][iT,ia[iP],iAvailable]=0
							ec['Growth Factor'][iT,ia[iP],iAvailable]=-10
							ec['ID Growth Curve'][iT,ia[iP],iAvailable]=1

							if np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True:
								TimeBetweenFireAndPlant=2
								iT=np.where(tv==YearIncite+TimeBetweenFireAndPlant)[0]
								iAvailable=np.where(ec['ID Event Type'][iT,ia[iP],:]==0)[0]

								ID_GC=np.max(gc[iScn][indBat[ia[iP]]]['ID_GC'])

								ec['ID Event Type'][iT,ia[iP],iAvailable]=meta['LUT']['Event']['Planting']
								ec['Mortality Factor'][iT,ia[iP],iAvailable]=0
								ec['Growth Factor'][iT,ia[iP],iAvailable]=10
								ec['ID Growth Curve'][iT,ia[iP],iAvailable]=ID_GC

				#------------------------------------------------------------------
				# Compress by indexing into the elements with information
				#------------------------------------------------------------------

				ec['idx']=np.where(ec['ID Event Type']>0)
				ec['ID Event Type']=ec['ID Event Type'][ec['idx']]
				ec['Mortality Factor']=ec['Mortality Factor'][ec['idx']]
				ec['Growth Factor']=ec['Growth Factor'][ec['idx']]
				ec['ID Growth Curve']=ec['ID Growth Curve'][ec['idx']]
				ec['ASET']=ec['ASET'][ec['idx']]

				gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl',ec)

	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Prepare growth curves
	#--------------------------------------------------------------------------

	print('Preparing growth curve input files')
	t0=time.time()
	cbu.PrepGrowthCurvesUniqueForCBR(meta,pNam,ugc)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Save data
	#--------------------------------------------------------------------------

	print('Saving input files')
	t0=time.time()
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Metadata.pkl',meta)
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Metadata_backup.pkl',meta)
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\dmec.pkl',dmec)
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\lsat.pkl',lsat)
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\gc.pkl',gc)
	gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\ugc.pkl',ugc)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	#--------------------------------------------------------------------------
	# Delete all output files
	#--------------------------------------------------------------------------

	print('Deleting any output files')
	t0=time.time()
	cbu.DeleteAllOutputFiles(meta,pNam)
	print(str(np.round((time.time()-t0)/60,decimals=1)) + ' min')

	return meta,lsat,dmec

#%% Timber harvesting land base
def DefineTHLB(meta,pNam,lsat):
	
	# Initialize THLB flags for each land cover/land use change scenario (THLB=1,Non-THLB=0)
	thlb={}
	
	# Scenario 0: No REARs are removed from THLB historically or in the future
	thlb['S0']=np.ones((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
	
	# Scenario 1: Historical REARs are removed from THLB, but there is no 
	# no subequent establishment of REARs
	thlb['S1']=np.ones((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')

	# Scenario 2: Historical REARs are removed from THLB, with subsequent establishment
	# of REARs also removed from THLB
	thlb['S2']=np.ones((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
	
	# Index to stands that are uneconomic
	iUneconomic=np.where(lsat['SI']<=5)[0]
	if iUneconomic.size>0:
		# Remove uneconomic stands from THLB
		for k in thlb.keys():
			thlb[k][:,iUneconomic]=0

	# Idenify stands that have been harvested
	# has_been_harvested=np.zeros(meta[pNam]['Project']['N Stand'])
	# for iStand in range(len(dmec)):
	#	 ind=np.where( (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) | (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest Salvage']) )[0]
	#	 if ind.size>0:
	#		 has_been_harvested[iStand]=1

	# Index to stands that have not been harvested
	iNoHarv=np.where( (lsat['Harvest Year Comp2']==0) & (lsat['SI']>5) )[0]

	# Use the ratio of THLB to non-THLB as an indicator of what will be harvested
	# among remaining primary forest
	ratio_thlb=22/62 # ratio of THLB to total forest (SOF)

	# Ratio of uneconomic to total (needed to adjust probability)
	corr=iUneconomic.size/meta[pNam]['Project']['N Stand']

	# Probability of evading harvest
	if iNoHarv.size>0:
		p_evade=(1-ratio_thlb-corr)*(meta[pNam]['Project']['N Stand']/iNoHarv.size)
	else:
		p_evade=(1-ratio_thlb-corr)

	# Random prediction of whether it will evade harvesting
	iRem=np.where(np.random.random(iNoHarv.size)<p_evade)[0]
	for k in thlb.keys():
		thlb[k][:,iNoHarv[iRem]]=0

	# np.sum(thlb['Actual'][0,:])/meta[pNam]['Project']['N Stand']

	#--------------------------------------------------------------------------
	# Transitions
	#--------------------------------------------------------------------------

	# Initialize year of transition
	thlb_YearTransitionOut=np.zeros(meta[pNam]['Project']['N Stand'])

	# Historical transitions
	ind=np.where( (lsat['LandCover_Comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (lsat['LandUse_Comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) )
	try:
		thlb_YearTransitionOut[ind]=lsat['LandUseChange_Comp1_1800to2019_Year'][ind]
	except:
		print(thlb_YearTransitionOut.shape)
		print(ind[0].size)
		print(lsat['LandUseChange_Comp1_1800to2019_Year'].shape)
		print(lsat['LandCover_Comp1_2019'].shape)

	ind=np.where( (lsat['LandCover_Comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (lsat['LandUse_Comp1_2019']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
	try:
		thlb_YearTransitionOut[ind]=lsat['LandUseChange_Comp1_1800to2019_Year'][ind]
	except:
		print(thlb_YearTransitionOut.shape)
		print(ind[0].size)
		print(lsat['LandUseChange_Comp1_1800to2019_Year'].shape)
		print(lsat['LandCover_Comp1_2019'].shape)

	# Apply historical transitions to Scenario 1 and 2
	for j in range(thlb_YearTransitionOut.size):
		if thlb_YearTransitionOut[j]>0:
			it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
			thlb['S1'][it,j]=0

	# Add future transitions to historical transitions
	ind=np.where( (lsat['LandCover_Comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (lsat['LandUse_Comp1_2019']!=meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) & (lsat['LandUse_Comp1_2049_Scn2']==meta['LUT']['Derived']['lu_comp1']['Conservation Natural']) )
	thlb_YearTransitionOut[ind]=lsat['LandUseChange_Comp1_2020to2049_Scn2_Year'][ind]
	ind=np.where( (lsat['LandCover_Comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (lsat['LandUse_Comp1_2019']!=meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) & (lsat['LandUse_Comp1_2049_Scn2']==meta['LUT']['Derived']['lu_comp1']['Conservation Consistent']) )
	thlb_YearTransitionOut[ind]=lsat['LandUseChange_Comp1_2020to2049_Scn2_Year'][ind]
	
	# Apply historical and future transitions to Scenario 2
	for j in range(thlb_YearTransitionOut.size):
		if thlb_YearTransitionOut[j]>0:
			it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
			thlb['S2'][it,j]=0

	# #--------------------------------------------------------------------------
	# # Baselines
	# #--------------------------------------------------------------------------
	# # Adjust the baseline so that simulated harvesting between 1995 and 2022 only
	# # occurs in areas where the THLB was affected by value diversification
	# for year in range(1990,2023,1):
	#	 iT=np.where(meta[pNam]['Year']==year)[0]
	#	 iS=np.where( (thlb['Baseline'][iT,:]==1) & (thlb['Actual'][iT,:]==1) )[1]
	#	 thlb['Baseline'][iT,iS]=0
	#	 iS=np.where( (thlb['Scn1 Baseline'][iT,:]==1) & (thlb['Scn1 Actual'][iT,:]==1) )[1]
	#	 thlb[iScn]['Scn1 Baseline'][iT,iS]=0
	
	lsat['THLB']=thlb
	return lsat

# def DefineTHLB(meta,pNam,lsat,dmec,lsc):

#	 thlb=[{}]*meta[pNam]['Project']['N Scenario']

#	 for iScn in range(meta[pNam]['Project']['N Scenario']):

#		 #--------------------------------------------------------------------
#		 # Initialize THLB flags (THLB=1,Non-THLB=0)
#		 #--------------------------------------------------------------------
#		 # Initially assume everything is in the THLB
#		 thlb[iScn]['Actual']=np.ones((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')
#		 thlb[iScn]['Baseline']=thlb[iScn]['Actual'].copy()

#		 thlb[iScn]['Scn1 Actual']=thlb[iScn]['Actual'].copy()
#		 thlb[iScn]['Scn1 Baseline']=thlb[iScn]['Actual'].copy()

#		 # Index to stands that are uneconomic
#		 iUneconomic=np.where(lsat['SI']<=5)[0]
#		 if iUneconomic.size>0:
#			 # Remove uneconomic stands from THLB
#			 thlb[iScn]['Actual'][:,iUneconomic]=0
#			 thlb[iScn]['Baseline'][:,iUneconomic]=0
#			 thlb[iScn]['Scn1 Actual'][:,iUneconomic]=0
#			 thlb[iScn]['Scn1 Baseline'][:,iUneconomic]=0

#		 # Idenify stands that have been harvested
#		 has_been_harvested=np.zeros(meta[pNam]['Project']['N Stand'])
#		 for iStand in range(len(dmec)):
#			 ind=np.where( (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest']) | (dmec[iScn][iStand]['ID Event Type']==meta['LUT']['Event']['Harvest Salvage']) )[0]
#			 if ind.size>0:
#				 has_been_harvested[iStand]=1

#		 # Index to stands that have not been harvested
#		 iNoHarv=np.where( (has_been_harvested==0) & (lsat['SI']>5) )[0]

#		 # Use the ratio of THLB to non-THLB as an indicator of what will be harvested
#		 # among remaining primary forest
#		 ratio_thlb=22/60 # ratio of THLB to total forest (SOF)

#		 # Ratio of uneconomic to total (needed to adjust probability)
#		 corr=iUneconomic.size/meta[pNam]['Project']['N Stand']

#		 # Probability of evading harvest
#		 if iNoHarv.size>0:
#			 p_evade=(1-ratio_thlb-corr)*(meta[pNam]['Project']['N Stand']/iNoHarv.size)
#		 else:
#			 p_evade=(1-ratio_thlb-corr)

#		 # Random prediction of whether it will evade harvesting
#		 iRem=np.where(np.random.random(iNoHarv.size)<p_evade)[0]
#		 thlb[iScn]['Actual'][:,iNoHarv[iRem]]=0
#		 thlb[iScn]['Baseline'][:,iNoHarv[iRem]]=0

#		 thlb[iScn]['Scn1 Actual'][:,iNoHarv[iRem]]=0
#		 thlb[iScn]['Scn1 Baseline'][:,iNoHarv[iRem]]=0

#		 # np.sum(thlb['Actual'][0,:])/meta[pNam]['Project']['N Stand']

#		 #------------------------------------------------------------------------------
#		 # Actual (New based on REAR layer)
#		 #------------------------------------------------------------------------------

#		 # Initialize year of transition
#		 thlb_YearTransitionOut=np.zeros(meta[pNam]['Project']['N Stand'])

#		 ind=np.where(lsat['THLB Layer']==1)[0]
#		 thlb_YearTransitionOut[ind]=1990

#		 # Conservation from land surface classification
#		 name=meta[pNam]['Scenario'][iScn]['Land Surface Scenario']
#		 if name!='Off':
#			 idx=LSC_Scenario_Crosswalk(lsc,name)
#			 Use=np.reshape(lsc['Scenarios'][idx]['Use'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))
#			 ind=np.where( Use==meta['LUT']['LSC']['Use']['Conservation Consistent'] )
#			 if ind[0].size>0:
#				 for i in range(ind[0].size):
#					 thlb_YearTransitionOut[ind[1][i]]=lsc['tv'][ind[0][i]]

#		 # Apply transition to actual THLB
#		 for j in range(thlb_YearTransitionOut.size):
#			 if thlb_YearTransitionOut[j]>0:
#				 it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
#				 thlb[iScn]['Actual'][it,j]=0
#				 thlb[iScn]['Scn1 Actual'][it,j]=0

#		 #------------------------------------------------------------------------------
#		 # Scn1 Actual (with deferrals + random areas to achieve 30 by 30)
#		 #------------------------------------------------------------------------------

#		 thlb_YearTransitionOut=np.zeros(meta[pNam]['Project']['N Stand'])

#		 ind=np.where(lsat['THLB Layer']==2)[0]
#		 thlb_YearTransitionOut[ind]=2023

#		 # Apply transition to actual THLB
#		 for j in range(thlb_YearTransitionOut.size):
#			 if thlb_YearTransitionOut[j]>0:
#				 it=np.where( (meta[pNam]['Year']>=thlb_YearTransitionOut[j]) )[0]
#				 thlb[iScn]['Scn1 Actual'][it,j]=0

#		 #------------------------------------------------------------------------------
#		 # Baselines
#		 #------------------------------------------------------------------------------
#		 # Adjust the baseline so that simulated harvesting between 1995 and 2022 only
#		 # occurs in areas where the THLB was affected by value diversification
#		 for year in range(1990,2023,1):
#			 iT=np.where(meta[pNam]['Year']==year)[0]
#			 iS=np.where( (thlb[iScn]['Baseline'][iT,:]==1) & (thlb[iScn]['Actual'][iT,:]==1) )[1]
#			 thlb[iScn]['Baseline'][iT,iS]=0
#			 iS=np.where( (thlb[iScn]['Scn1 Baseline'][iT,:]==1) & (thlb[iScn]['Scn1 Actual'][iT,:]==1) )[1]
#			 thlb[iScn]['Scn1 Baseline'][iT,iS]=0

#	 return thlb

#%% Get index to a scenario in the land surface class list

# def LSC_Scenario_Crosswalk(lsc,name):
#	 if 'Scenarios' not in lsc:
#		 return
#	 for i in range(len(lsc['Scenarios'])):
#		 if lsc['Scenarios'][i]['Name']==name:
#			 return i

#%% Get unique growth curves
def ExtractUniqueGrowthCurves(meta,pNam,gc):

	ugc={}
	ugc['GC_Variable_List']=np.array(meta['Modules']['GYM']['GC_Variable_List'])[3:]

	# Calculate unique stand types
	ugc['Full']=np.zeros((int(10e6),len(meta['Modules']['GYM']['GC_Variable_List'])))

	cnt=0
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		for iStand in range(meta[pNam]['Project']['N Stand']):
			gc0=gc[iScn][iStand]
			for iGC in range(gc0['ID_GC'].size):
				for k in range(len(meta['Modules']['GYM']['GC_Variable_List'])):
					key=meta['Modules']['GYM']['GC_Variable_List'][k]
					ugc['Full'][cnt,k]=gc0[key][iGC]
				cnt=cnt+1
	ugc['Full']=ugc['Full'][0:cnt,:]

	# Unique growth curves
	# The 'Inverse' variable acts as the crosswalk between the full and unique gc arrays
	ugc['Unique'],ugc['Index'],ugc['Inverse']=np.unique(ugc['Full'][:,3:],return_index=True,return_inverse=True,axis=0)

	return ugc

#%% Add changes in land surface classfication to DMEC

def AddLandSurfaceChangesToDMEC(meta,pNam,dmec,lsc):

	for iScn in range(meta[pNam]['Project']['N Scenario']):

		# Name of LS scenario for each scenario
		sNam=meta[pNam]['Scenario'][iScn]['Land Surface Scenario']

		# Index to LSC scenario
		for i in range(len(lsc['Scenarios'])):
			if lsc['Scenarios'][i]['Name']==sNam:
				idx=i
				break

		if sNam!='Off':

			#Cover=np.reshape(lsc['Scenarios'][idx]['Cover'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))
			Use=np.reshape(lsc['Scenarios'][idx]['Use'].copy(),(lsc['tv'].size,meta[pNam]['Project']['N Stand']))

			#----------------------------------------------------------------------
			# Fuel breaks
			#----------------------------------------------------------------------

			nam='Fuel Break'
			indS=np.unique(np.where( Use==meta['LUT']['LSC']['Use'][nam] )[1])
			if indS.size>0:
				for i in range(indS.size):

					iS=indS[i]
					iT=np.where(Use[:,iS]==meta['LUT']['LSC']['Use'][nam])[0][0]

					# Add harvest
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT])
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

					# Add Pile Burn
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+1)
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Pile Burn'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

					# Add planting
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2)
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Planting'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],0)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)
					dmec[iScn][iS]['PL_SPECIES_CD1'][-1]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']
					dmec[iScn][iS]['PL_SPECIES_PCT1'][-1]=100

					# Add harvest
					RotationLength=14
					for iR in range(1,10):
						dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2+iR*RotationLength)
						dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
						dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
						dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
						if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
							dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
						for v in meta['Core']['StringsToFill']:
							dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

			#----------------------------------------------------------------------
			# Energy Production
			#----------------------------------------------------------------------

			nam='Energy Production'
			indS=np.unique(np.where( Use==meta['LUT']['LSC']['Use'][nam] )[1])
			if indS.size>0:
				for i in range(indS.size):

					iS=indS[i]
					iT=np.where(Use[:,iS]==meta['LUT']['LSC']['Use'][nam])[0][0]

					# Add harvest
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT])
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

					# Add Pile Burn
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+1)
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Pile Burn'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

					# Add planting
					dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2)
					dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Planting'])
					dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],0)
					dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
					if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
						dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
					for v in meta['Core']['StringsToFill']:
						dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)
					dmec[iScn][iS]['PL_SPECIES_CD1'][-1]=meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']
					dmec[iScn][iS]['PL_SPECIES_PCT1'][-1]=100

					# Add harvest
					RotationLength=14
					for iR in range(1,10):
						dmec[iScn][iS]['Year']=np.append(dmec[iScn][iS]['Year'],lsc['tv'][iT]+2+iR*RotationLength)
						dmec[iScn][iS]['ID Event Type']=np.append(dmec[iScn][iS]['ID Event Type'],meta['LUT']['Event']['Harvest'])
						dmec[iScn][iS]['Mortality Factor']=np.append(dmec[iScn][iS]['Mortality Factor'],100)
						dmec[iScn][iS]['Growth Factor']=np.append(dmec[iScn][iS]['Growth Factor'],9999)
						if 'Index to Event Inciting NOSE' in dmec[iScn][iS]:
							dmec[iScn][iS]['Index to Event Inciting NOSE']=np.append(dmec[iScn][iS]['Index to Event Inciting NOSE'],9999)
						for v in meta['Core']['StringsToFill']:
							dmec[iScn][iS][v]=np.append(dmec[iScn][iS][v],9999)

	return dmec

#%% Exclude unidentified events
def Exclude_Unidentified_Events(meta,pNam,dmec0):
	N_Removed=0
	for iStand in range(meta[pNam]['Project']['N Stand']):
		if dmec0[iStand]==None:
			continue
		ind=np.where(dmec0[iStand]['ID Event Type']!=9999)[0]
		N_Removed=N_Removed+ind.size
		for key in dmec0[iStand]:
			dmec0[iStand][key]=dmec0[iStand][key][ind]
	print(str(N_Removed) + ' unidentifiable events removed.')
	return dmec0

def Exclude_Unidentified_Events2(meta,pNam,dmec):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		N_Removed=0
		for iStand in range(meta[pNam]['Project']['N Stand']):
			if dmec[iScn][iStand]==None:
				continue
			ind=np.where(dmec[iScn][iStand]['ID Event Type']!=9999)[0]
			N_Removed=N_Removed+ind.size
			for key in dmec[iScn][iStand]:
				dmec[iScn][iStand][key]=dmec[iScn][iStand][key][ind]
		print(str(N_Removed) + ' unidentifiable events removed.')
	return dmec

#%% Exclude duplicate events from DMEC
def Exclude_Duplicate_Events(meta,pNam,dmec):
	for iStand in range(meta[pNam]['Project']['N Stand']):
		if dmec[iStand]==None:
			continue
		for key in meta['LUT']['Event'].keys():
			indType=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event'][key]) )[0]
			if indType.size==0:
				continue
			uYear=np.unique(np.floor(dmec[iStand]['Year'][indType]))
			for iYear in range(uYear.size):
				ind=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event'][key]) & (np.floor(dmec[iStand]['Year'])==uYear[iYear]) )[0]
				if ind.size>1:
					dmec[iStand]['ID Event Type'][ind[:-1]]=9999
					dmec[iStand]['Year'][ind[:-1]]=9999
					dmec[iStand]['Mortality Factor'][ind[:-1]]=9999
	return dmec

#%% Put DMEC events in order
def PutEventsInOrder1(meta,pNam,dmec):
	for iStand in range(meta[pNam]['Project']['N Stand']):
		d=dmec[iStand].copy()
		Ord=np.argsort(d['Year'])
		for key in d.keys():
			d[key]=d[key][Ord]
		dmec[iStand]=d.copy()
	return dmec

#%%
def PutEventsInOrder2(meta,pNam,dmec):

	def count_lists(l):
		return sum(1 + count_lists(i) for i in l if isinstance(i,list))

	if count_lists(dmec)==0:
		# DMEC has not been split up by scenario yet
		for iStand in range(meta[pNam]['Project']['N Stand']):
			d=dmec[iStand].copy()
			Ord=np.argsort(d['Year'])

			# Fix index to inciting events
			if 'Index to Event Inciting NOSE' in d.keys():
				indOld=np.where( (d['Index to Event Inciting NOSE']>=0) & (d['Index to Event Inciting NOSE']<9999) )[0]
				if indOld.size==0:
					# Some stands may not be represented (e.g. direct seeding)
					continue
				#indNew=np.zeros(indOld.size,dtype='int')
				#for j in range(indOld.size):
				#	indNew[j]=d['Index to Event Inciting NOSE'][indOld[j]]
				indNew=np.where(Ord==d['Index to Event Inciting NOSE'][indOld])[0]
				if indNew.size>0:
					d['Index to Event Inciting NOSE'][indOld]=indNew
				d['Index to Event Inciting NOSE'][indOld]=indNew

			# Update order
			for key in d.keys():
				d[key]=d[key][Ord]
			dmec[iStand]=d.copy()
	else:
		for iScn in range(meta[pNam]['Project']['N Scenario']):
			for iStand in range(meta[pNam]['Project']['N Stand']):
				d=dmec[iScn][iStand].copy()
				Ord=np.argsort(d['Year'])
				# Fix index to inciting events
				if 'Index to Event Inciting NOSE' in d.keys():
					#indOld=np.where(d['Index to Event Inciting NOSE']>=0)[0]
					indOld=np.where( (d['Index to Event Inciting NOSE']>=0) & (d['Index to Event Inciting NOSE']<9999) )[0]
					indNew=np.where(Ord==d['Index to Event Inciting NOSE'][indOld])[0]
					if indNew.size>0:
						d['Index to Event Inciting NOSE'][indOld]=indNew
				for key in d.keys():
					d[key]=d[key][Ord]
				dmec[iScn][iStand]=d.copy()
	return dmec

#%% Ensure disturbance prededes fertilization so age is specified
# So that age at fert is specified.
def Ensure_Fert_Preceded_By_Disturbance(meta,pNam,dmec,ba):

	ListOfTestedDist=[meta['LUT']['Event']['Wildfire'],
					  meta['LUT']['Event']['Harvest'],
					  meta['LUT']['Event']['Knockdown'],
					  meta['LUT']['Event']['Mountain Pine Beetle'],
					  meta['LUT']['Event']['Balsam Beetle'],
					  meta['LUT']['Event']['Douglas-fir Beetle'],
					  meta['LUT']['Event']['Spruce Beetle']]

	for iStand in range(meta[pNam]['Project']['N Stand']):

		if dmec[iStand]==None:
			continue

		iFert=np.where( (dmec[iStand]['ID Event Type']==meta['LUT']['Event']['Nutrient App Aerial']) )[0]

		if iFert.size==0:
			continue

		iFert=iFert[0]

		# Index to events prior to first fertilization with 100% mortality
		ind=np.where( (dmec[iStand]['Year']<dmec[iStand]['Year'][iFert]) & (dmec[iStand]['Mortality Factor']==100) & np.isin(dmec[iStand]['ID Event Type'],ListOfTestedDist) )[0]

		if ind.size>0:
			continue

		print('Adding a harvest before nutrient application')

		# Assume mean of 38 + random variation (planting is 2 years after harvest)
		r=38+np.random.randint(-6,high=6)
		Year=dmec[iStand]['Year'][iFert]-r

		# Add harvest
		for k in dmec[iStand].keys():
			dmec[iStand][k]=np.append(dmec[iStand][k],9999)
		dmec[iStand]['Year'][-1]=Year
		dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Harvest']
		dmec[iStand]['Event Source'][-1]=2
		dmec[iStand]['Mortality Factor'][-1]=100
		dmec[iStand]['Growth Factor'][-1]=9999

		# Add Pile Burn
		for k in dmec[iStand].keys():
			dmec[iStand][k]=np.append(dmec[iStand][k],9999)
		dmec[iStand]['Year'][-1]=Year+1
		dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Pile Burn']
		dmec[iStand]['Event Source'][-1]=2
		dmec[iStand]['Mortality Factor'][-1]=100
		dmec[iStand]['Growth Factor'][-1]=9999

		# Add planting
		for k in dmec[iStand].keys():
			dmec[iStand][k]=np.append(dmec[iStand][k],9999)
		dmec[iStand]['Year'][-1]=Year+2
		dmec[iStand]['ID Event Type'][-1]=meta['LUT']['Event']['Planting']
		dmec[iStand]['Event Source'][-1]=2
		dmec[iStand]['Mortality Factor'][-1]=0
		dmec[iStand]['Growth Factor'][-1]=9999

	return dmec

#%% Reduce number of growth curves by adjusting site index
def ReduceVariationInSiteIndex(meta,pNam,lsat):
	trig=0
	for i in range(1,55,2):
		ind=np.where(lsat['SI']==i)[0]
		if trig==0:
			lsat['SI'][ind]=lsat['SI'][ind]+1
			trig=1
		else:
			lsat['SI'][ind]=lsat['SI'][ind]-1
			trig=0
	return lsat

#%% Ensure every stand has a modern disturbance event
def Ensure_Every_Stand_Has_Modern_Disturbance(meta,pNam,dmec,name_dist,severity):

	for iStand in range(meta[pNam]['Project']['N Stand']):

		if dmec[iStand]==None:
			continue

		if dmec[iStand]['Year'].size==0:
			#print(iStand)
			#break
			r=np.random.randint(1700,2000)
			dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],r)
			dmec[iStand]['ID Event Type']=np.append(dmec[iStand]['ID Event Type'],meta['LUT']['Event'][name_dist])
			dmec[iStand]['Event Source']=np.append(dmec[iStand]['Event Source'],2)
			dmec[iStand]['Mortality Factor']=np.append(dmec[iStand]['Mortality Factor'],np.array(severity,dtype='int16'))
			dmec[iStand]['Growth Factor']=np.append(dmec[iStand]['Growth Factor'],np.array(0,dtype='int16'))
			#if 'FCI Funded' in dmec[iStand]:
			#	dmec[iStand]['FCI Funded']=np.append(dmec[iStand]['FCI Funded'],np.array(0,dtype='int16'))
			for v in meta['Core']['StringsToFill']:
				dmec[iStand][v]=np.append(dmec[iStand][v],9999)

	return dmec

#%% Generate DMEC from estimate of stand age
def GapFill_DMEC_WithAge(meta,pNam,dmec,lsat):
	for iStand in range(meta[pNam]['Project']['N Stand']):
		
		if dmec[iStand]==None:
			continue
		
		Age=np.minimum(300,lsat['Age'])
		
		# if dmec[iStand]['Mortality Factor'].size>0:
		#	 if np.max(dmec[iStand]['Mortality Factor']==100):
		#		 # If there is a stand-replacing disturbance on record, no need to proceed
		#		 continue
	
		# Allowing gap-filling of really young stands can be very problematic
		th=30
		if Age[iStand]>=th:
			r=np.random.random(1)
			if r<0.33:
				type=meta['LUT']['Event']['Wildfire']
			else:
				type=meta['LUT']['Event']['Mountain Pine Beetle']
			for k in dmec[iStand].keys():
				dmec[iStand][k]=np.append(dmec[iStand][k],9999)
			dmec[iStand]['Year'][-1]=meta[pNam]['Project']['Year Project']-Age[iStand]
			dmec[iStand]['ID Event Type'][-1]=type
			dmec[iStand]['Event Source'][-1]=2
			dmec[iStand]['Mortality Factor'][-1]=np.array(100,dtype='int16')
			dmec[iStand]['Growth Factor'][-1]=9999
	return dmec

#%%
def ImportEnvironmentalData(meta,pNam):
	# Notes:
	# Used by GROMO
	# The zscoring needs to be done on the fly unless willing to save both for
	# growth and mortality?

	# Specific time period for environmental inputs
	tv=np.arange(1851,2151,1)

	# Import project geospatial info
	geos=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\geos.pkl')

	# Get crosswalk between project stands and NA1k database
	#srs=gis.ImportSRSs()
	#lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],geos['Sparse']['X'],geos['Sparse']['Y'])
	#x,y=srs['Proj']['NACID'](lon,lat)

	# Import bc5k grid info
	z5k=gis.OpenGeoTiff(meta['Paths']['bc5k Ref Grid'])
	idxBC5k=gis.GetGridIndexToPoints(z5k,geos['Sparse']['X'],geos['Sparse']['Y'])

	# Import normal climate
	zTn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_tmean_ann_norm_1971to2000.tif')
	Tn=gis.UpdateGridCellsize(zTn,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]
	zWn=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Summaries\\Normals\\bc1ha_ws_mjjas_norm_1971to2000.tif')
	Wn=gis.UpdateGridCellsize(zWn,meta['Geos']['RGSF'])['Data'][ meta['Geos']['iMask'] ]

	# Import CRU climate anomalies
	dCRU=gu.ipickle(meta['Paths']['DB']['CRU'] + '\\cru_anomalies.pkl')
	tvCRU=np.arange(1851,2023)

	# Import CO2 (ppm), source: https://www.pik-potsdam.de/~mmalte/rcps/ (ppm)
	dCO2=gu.ReadExcel(meta['Paths']['DB']['CO2'])

	# Scenarios
	for iScn in range(meta[pNam]['Project']['N Scenario']):

		# Specifiy scenario strings
		cmS=meta[pNam]['Scenario'][iScn]['gromo Scenario']
		if cmS=='ssp245':
			ndS='26'
		else:
			ndS='85'
		cmM=meta[pNam]['Scenario'][iScn]['gromo GCM']

		# Import future climate
		dCM=gu.ipickle(meta['Paths']['DB']['CMIP6'] + '\\cmip6_anomalies_' + cmS + '.pkl')[cmM]
		tvCM=np.arange(1950,2101,1)

		# Import N deposition
		dND=gu.ipickle(meta['Paths']['DB']['NDEP'] + '\\ISIMIP\\ndep_' + ndS + '.pkl')

		# Compile environmental data (keep it in native precision)
		envd={}
		envd['Tn']=Tn
		envd['Wn']=Wn
		envd['Ta']=np.zeros( (tv.size,meta[pNam]['Project']['N Stand']),dtype='int16' )
		envd['Wa']=np.zeros( (tv.size,meta[pNam]['Project']['N Stand']),dtype='int16' )
		envd['ND']=np.zeros( (tv.size,meta[pNam]['Project']['N Stand']),dtype='int16' )
		envd['CO2']=np.zeros( (tv.size,meta[pNam]['Project']['N Stand']),dtype='int16' )
		for iT in range(tv.size):

			# Add historical climate data (ends in 2022)
			indT=np.where(tvCRU==tv[iT])[0]
			if indT.size>0:
				envd['Ta'][iT,:]=dCRU['tmean']['wyr'][indT[0],:,:][idxBC5k]#.astype('float')*meta['Climate']['SF']['tmean']
				envd['Wa'][iT,:]=dCRU['ws']['mjjas'][indT[0],:,:][idxBC5k]

			# Add N deposition
			indT=np.where(dND['ndep']['tv']==tv[iT])[0]
			if indT.size>0:
				envd['ND'][iT,:]=dND['ndep']['Data'][indT[0],:,:][idxBC5k]
			if (indT.size==0) & (tv[iT]<1900):
				envd['ND'][iT,:]=np.mean(dND['ndep']['Data'][0:10,:,:],axis=0)[idxBC5k]
			if (indT.size==0) & (tv[iT]>2050):
				envd['ND'][iT,:]=np.mean(dND['ndep']['Data'][-10:-1,:,:],axis=0)[idxBC5k]

			# Add future climate
			if (tv[iT]>=2023):
				indT=np.where( (tvCM==tv[iT]) )[0]
				if indT.size>0:
					envd['Ta'][iT,:]=dCM['tmean']['ann'][indT[0],:,:][idxBC5k]
					envd['Wa'][iT,:]=dCM['ws']['mjjas'][indT[0],:,:][idxBC5k]
					#envd['Ea'][iT,:]=dCM['etp']['mjjas'][indT[0],:,:][idxBC5k]
				if (tv[iT]>=2100):
					rn=np.random.randint(-15,high=-1)
					envd['Ta'][iT,:]=dCM['tmean']['ann'][rn,:,:][idxBC5k]
					envd['Wa'][iT,:]=dCM['ws']['mjjas'][rn,:,:][idxBC5k]
					#envd['Ea'][iT,:]=dCM['etp']['mjjas'][rn,:,:][idxBC5k]

			# Add CO2 (convert to integer for storage)
			indC=np.where(dCO2['Year']==tv[iT])[0]
			envd['CO2'][iT,:]=dCO2['CO2 RCP45'][indC]/meta['Climate']['SF']['ca']

		for k in envd.keys():
			envd[k]=envd[k].astype('int16')
	
		# Save to file
		for iBat in range(meta[pNam]['Project']['N Batch']):
			indBat=cbu.IndexToBatch(meta[pNam],iBat)
			envdBat={}
			for k in envd.keys():
				if (k=='Tn') | (k=='Wn') | (k=='En'):
					envdBat[k]=envd[k][indBat]
				else:
					envdBat[k]=envd[k][:,indBat]
			gu.opickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Environment_Bat' + cbu.FixFileNum(iBat) + '.pkl',envdBat)

	return