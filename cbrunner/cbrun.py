import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import gc as garc
import time
import datetime
import fcgadgets.macgyver.util_general as gu
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.cbrunner.cbrun_annproc as annproc

#%% Run simulation
def MeepMeep(meta,pNam):

	# Track time
	tStart=time.time()

	# Identify scenarios to run
	# If 'Scenario Override' is absent, it defaults to running all scenarios.
	# Specifying a list of scenarios in 'Scenario Override' facilitates running
	# the simulation in multiple instances of Python.
	if 'Scenario Override' in meta[pNam]['Project']:
		ScenariosToRun=meta[pNam]['Project']['Scenario Override']['Scenario List']
	else:
		ScenariosToRun=range(meta[pNam]['Project']['N Scenario'])

	# Loop through batches
	for iBat in range(meta[pNam]['Project']['N Batch']):
		
		# Initialize flag indicating acrive batch
		flag_WorkingOnBatch=0

		# Index to batch
		meta[pNam]['Project']['indBat']=cbu.IndexToBatch(meta[pNam],iBat)

		# Loop through scenarios
		for iScn in ScenariosToRun:

			# Track the current scenario index
			meta[pNam]['iScn']=iScn

			# Loop through ensembles
			for iEns in range(meta[pNam]['Project']['N Ensemble']):

				# Track time
				t0=time.time()

				# Path to temporary "working on batch..." file -> tells other instances
				pth_WorkingOnBatch=meta['Paths'][pNam]['Data'] + '\\Outputs\\WorkingOnBatch_' + cbu.FixFileNum(iBat) + '.pkl'

				# Only proceed if the file does not exist running multiple instances
				if meta[pNam]['Project']['Skip Completed Runs']=='On':
					if os.path.exists(pth_WorkingOnBatch)==False:
						gu.opickle(pth_WorkingOnBatch,[])
						flag_WorkingOnBatch=1
					else:
						if flag_WorkingOnBatch==0:
							print(pth_WorkingOnBatch)
							continue

				# Report progress
				if (meta[pNam]['Project']['Scenario Source']=='Spreadsheet') | ('Skip Run Update' in meta[pNam].keys()):
					pass
				else:
					print('Running -> Scenario:' + cbu.FixFileNum(iScn) + ', Ensemble:' + cbu.FixFileNum(iEns) + ', Batch:' + cbu.FixFileNum(iBat) + ', Runtime:' + str(np.round((t0-tStart)/60,decimals=1)) + ' min')

				# Initialize stands
				meta,vi,vo=InitializeStands(meta,pNam,iScn,iEns,iBat)
				t1=time.time()
				meta[pNam]['Project']['Run Time Summary']['Stand initialization']=meta[pNam]['Project']['Run Time Summary']['Stand initialization']+(t1-t0)

				# Set location-specific parameters
				meta,vi=PrepareParametersForBatch(meta,pNam,vi,iEns,iBat,iScn)
				t2=time.time()
				meta[pNam]['Project']['Run Time Summary']['Set location-specific parameters']=meta[pNam]['Project']['Run Time Summary']['Set location-specific parameters']+(t2-t1)

				# Indices to ecosystem pools
				iEP=meta['Core']['iEP']

				#--------------------------------------------------------------
				# Start at a later date and adopt spinup from a previous run
				#--------------------------------------------------------------
				if (meta[pNam]['Project']['Early Record Recycling']=='On') & (iEns!=0):

					# Start at an advanced date using previous data for spinup period
					Year_start=np.where(meta[pNam]['Year']==meta[pNam]['Project']['Year Start Saving']-meta['Core']['Recycle Early Record Buffer'])[0][0]

					# Get output variables from first run of this batch
					it=np.where(meta[pNam]['Year']==meta[pNam]['Project']['Year Start Saving']-meta['Core']['Recycle Early Record Buffer']-1)[0]
					for k in vo.keys():
						if (k=='C_M_DistByAgent') | (k=='C_M_DistByAgentPct'):
							continue
						if vo[k].size==0:
							continue
						vo[k][it,:]=vo_full[k]

					# Set future periods to zero
					it=np.where(meta[pNam]['Year']>=meta[pNam]['Project']['Year Start Saving']-meta['Core']['Recycle Early Record Buffer'])[0]
					for k in vo.keys():
						if (k=='C_M_DistByAgent') | (k=='C_M_DistByAgentPct'):
							# Nested dictionaries
							# *** These are already shortened time periods so set whole array to zero
							for k2 in vo[k].keys():
								vo[k][k2]=0*vo[k][k2]
						else:

							if vo[k].size==0:
								continue

							# Not a nested dictionary
							vo[k][it,:]=0*vo[k][it,:]

				else:

					# Start from beginning
					Year_start=1
					vo_full=[]

				# Track time
				t3=time.time()

				# Biomass dynamics from Sawtooth
				if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
					for iS in range(meta[pNam]['Project']['Batch Size'][iBat]):
						vo=annproc.BiomassFromSawtooth(iScn,iS,vi,vo,meta,iEP)
				t4=time.time()
				meta[pNam]['Project']['Run Time Summary']['Running biomass dynamics from sawtooth']=meta[pNam]['Project']['Run Time Summary']['Running biomass dynamics from sawtooth']+(t4-t3)

				# Force an actual scenario to match the baseline scenario prior
				# to the project year (e.g. future NM project)
				if 'Scenario Adjustment' in meta[pNam]['Project'].keys():
					if iScn==meta[pNam]['Project']['Scenario Adjustment']['Adjust Scenario']:
						Year_start=meta[pNam]['Project']['Year Project']
						iScnA=meta[pNam]['Project']['Scenario Adjustment']['With Scenario']
						fin=meta['Paths'][pNam]['Output Scenario'][iScnA] + '\\ScenarioAdjustment_Variables_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
						vo=gu.ipickle(fin)
						os.remove(fin)
						#d0=gu.ipickle(fin)
						#vo['A']=d0['A']
						#vo['C_Eco_ByPool']=d0['C_Eco_ByPool']
						#vo['C_Pro_ByPool']=d0['C_Pro_ByPool']
						#vi['GC']['Active']=d0['GCA']
						#os.remove(fin)
						#del d0
						fin=meta['Paths'][pNam]['Output Scenario'][iScnA] + '\\ScenarioAdjustment_GCActive_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
						vi['GC']['Active']=gu.ipickle(fin)
						os.remove(fin)
						vi['EC']=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScnA] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
						vi['EC']=cbu.EventChronologyDecompress(meta,pNam,vi['EC'],iScn,iEns,iBat)
						vi['EC']['Mortality Factor']=vi['EC']['Mortality Factor'].astype(float)/100 # Convert mortality percent to fraction

				# Loop through time intervals (start in second time step)
				for iT in range(Year_start,meta[pNam]['Project']['N Time']):
					
					t_Inn0=time.time()

					# If there are scenario adjustments, save active growth curves at the time of project year
					if 'Scenario Adjustment' in meta[pNam]['Project'].keys():
						if meta[pNam]['Year'][iT]==meta[pNam]['Project']['Year Project']:
							if iScn==meta[pNam]['Project']['Scenario Adjustment']['With Scenario']:
								#vo_out={}
								#vo_out['A']=vo['A']
								#vo_out['C_Eco_ByPool']=vo['C_Eco_ByPool']
								#vo_out['C_Pro_ByPool']=vo['C_Pro_ByPool']
								#vo_out['GCA']=vi['GC']['Active']
								#gu.opickle(fout,vo_out)
								fout=meta['Paths'][pNam]['Output Scenario'][iScn] + '\\ScenarioAdjustment_Variables_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
								gu.opickle(fout,vo)
								fout=meta['Paths'][pNam]['Output Scenario'][iScn] + '\\ScenarioAdjustment_GCActive_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
								gu.opickle(fout,vi['GC']['Active'])

					# Biomass dynamics
					if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':
						# Biomass from growth and yield model
						vo=annproc.TreeBiomassDynamicsFromGYModel(meta,pNam,iScn,iBat,iT,vi,vo,iEP)
						t_Inn1=time.time()
						meta[pNam]['Project']['Run Time Summary']['Biomass from GY model']=meta[pNam]['Project']['Run Time Summary']['Biomass from GY model']+(t_Inn1-t_Inn0)

					# Grassland dynamics
					if (meta[pNam]['Scenario'][iScn]['Grass Module Status']=='On') & (meta[pNam]['Year'][iT]>=meta[pNam]['Scenario'][iScn]['Grass Module Year Start']):
						vo=annproc.GrassBiomassDynamics(meta,pNam,iScn,iBat,iT,vi,vo,iEP)

					# Calculate annual dead organic matter dynamics
					vo=annproc.DeadWoodLitterAndSoilDynamics(meta,pNam,iT,iBat,vi,vo,iEP)
					t_Inn2=time.time()
					meta[pNam]['Project']['Run Time Summary']['DOM']=meta[pNam]['Project']['Run Time Summary']['DOM']+(t_Inn2-t_Inn1)

					# Calculate effects of disturbance and management
					vo,vi=annproc.DisturbanceAndManagementEvents(meta,pNam,iT,iScn,iEns,iBat,vi,vo,iEP)
					t_Inn3=time.time()
					meta[pNam]['Project']['Run Time Summary']['Events']=meta[pNam]['Project']['Run Time Summary']['Events']+(t_Inn3-t_Inn2)

					# Calculate products sector (no need to run this before a certain date)
					if meta[pNam]['Year'][iT]>=meta['Core']['HWP Year Start']:
						vo=annproc.ProductDynamics(meta,pNam,iT,iBat,vi,vo)
						t_Inn4=time.time()
						meta[pNam]['Project']['Run Time Summary']['HWP']=meta[pNam]['Project']['Run Time Summary']['HWP']+(t_Inn4-t_Inn3)

				t5=time.time()
				meta[pNam]['Project']['Run Time Summary']['Full annual loop']=meta[pNam]['Project']['Run Time Summary']['Full annual loop']+(t5-t4)
				
				# Calculate fossil fuel emissions from operations and substitution effects
				vo=annproc.GeologicalDynamics(meta,pNam,vi,vo)
				t6=time.time()
				meta[pNam]['Project']['Run Time Summary']['Run geological']=meta[pNam]['Project']['Run Time Summary']['Run geological']+(t6-t5)

				# Export simulation results to file
				vo_full=ExportSimulation(meta,pNam,vi,vo,iScn,iEns,iBat,iEP,vo_full)
				t7=time.time()
				meta[pNam]['Project']['Run Time Summary']['Export results to file']=meta[pNam]['Project']['Run Time Summary']['Export results to file']+(t7-t6)

				# Delete 'working on' file
				#if meta[pNam]['Project']['Skip Completed Runs']=='On':
					#os.remove(pthWO)

				# Delete variables
				del vi,vo
				garc.collect()
				t8=time.time()
				meta[pNam]['Project']['Run Time Summary']['Collect garbage']=meta[pNam]['Project']['Run Time Summary']['Collect garbage']+(t8-t7)

	tFinal=time.time()
	if (meta[pNam]['Project']['Scenario Source']=='Spreadsheet') | ('Skip Run Update' in meta[pNam].keys()):
		pass
	else:
		print('Run completed ' + str(np.round((tFinal-tStart)/60,decimals=1)) + ' min')

	return meta

#%% Phase shift correction of net growth
def PhaseShiftNG(meta,y):
	flg=1
	if flg==1:
		x=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)
		xCor=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)
		iCor=np.where(x<=55)[0]
		xCor[iCor]=0.5*x[iCor]
		for j in range(y.shape[1]):
			for k in range(y.shape[2]):
				y1=np.interp(x,xCor,y[:,j,k])
				y[:,j,k]=y1
	return y

#%% Initialize stands
def InitializeStands(meta,pNam,iScn,iEns,iBat):

	#--------------------------------------------------------------------------
	# Input variables
	#--------------------------------------------------------------------------

	vi={}

	vi['tv']=meta[pNam]['Year']
	tv=meta[pNam]['Year']
	meta[pNam]['Project']['N Time']=len(tv)

	# Import land surface attributes
	vi['lsat']=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl')

	# Get spatial indices to BGC zones
	tmp=gu.IndicesFromUniqueArrayValues(vi['lsat']['ID_BGCZ'].flatten())
	vi['lsat']['bgcz idx']={}
	for num in tmp.keys():
		cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],num)[0]
		vi['lsat']['bgcz idx'][cd]=tmp[num]

	# Create a biomass to volume conversion factor
	if meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
		rho=vi['lsat']['Wood Density'].astype('float')
	else:
		rho=vi['lsat']['Wood Density'].astype('float')/100
	vi['lsat']['Biomass to Volume CF']=(1/rho)*(1/meta['Param']['BE']['Biophysical']['Carbon Content Wood'])

	if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
		# Create relative index from spatial map of harvest annual probability
		Pa_H_Max=0.005
		vi['lsat']['Harvest Index']=(vi['lsat']['Prob Harvest (%/yr) x 1000'].astype('float')/1000/100)/Pa_H_Max
		
		# Was it harvested (yes=1, no=0)
		vi['lsat']['Harvest Flag']=np.minimum(1,vi['lsat']['Year Harvest First'])

	# Species-specific adjustments to insect mortality based on % species comp
	vi['lsat']=cbu.PrepInsectMortalityPercentTreeSpeciesAffected(meta,pNam,vi['lsat'],iBat)

	# Update number of stands for batch
	meta[pNam]['Project']['N Stand Batch']=vi['lsat']['ID_BGCZ'].shape[1]

	# Import event chronology
	vi['EC']=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')
	vi['EC']=cbu.EventChronologyDecompress(meta,pNam,vi['EC'],iScn,iEns,iBat)

	# Convert mortality percent to fraction
	vi['EC']['Mortality Factor']=vi['EC']['Mortality Factor'].astype(float)/100

	# Import growth curves
	if (meta[pNam]['Project']['Biomass Module']=='BatchTIPSY'):

		vi['GC']={}

		# Import growth curve 1
		vi['GC'][1]=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\GrowthCurve1_Bat' + cbu.FixFileNum(iBat) + '.pkl')
		vi['GC'][1]=PhaseShiftNG(meta,vi['GC'][1])

		# Set active growth curve to growth curve 1
		vi['GC']['Active']=vi['GC'][1].copy().astype(float)*meta['Modules']['GYM']['Scale Factor']

		# Initialize an indicator of the active growth curve ID
		vi['GC']['ID_GCA']=np.ones(meta[pNam]['Project']['Batch Size'][iBat])

		# Import growth curve 2
		vi['GC'][2]=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\GrowthCurve2_Bat' + cbu.FixFileNum(iBat) + '.pkl')
		vi['GC'][2]=PhaseShiftNG(meta,vi['GC'][2])

		# Import growth curve 3
		try:
			vi['GC'][3]=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\GrowthCurve3_Bat' + cbu.FixFileNum(iBat) + '.pkl')
			vi['GC'][3]=PhaseShiftNG(meta,vi['GC'][3])
		except:
			vi['GC'][3]=0

		# Import growth curve 4 (optional)
		try:
			vi['GC'][4]=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\GrowthCurve4_Bat' + cbu.FixFileNum(iBat) + '.pkl')
			vi['GC'][4]=PhaseShiftNG(meta,vi['GC'][4])
		except:
			vi['GC'][4]=0

		# Import growth curve 5 (optional)
		try:
			vi['GC'][5]=gu.ipickle(meta['Paths'][pNam]['Input Scenario'][iScn] + '\\GrowthCurve5_Bat' + cbu.FixFileNum(iBat) + '.pkl')
		except:
			vi['GC'][5]=0

	else:
		vi['GC']={}
		vi['GC'][1]=0
		vi['GC'][2]=0
		vi['GC'][3]=0
		vi['GC'][4]=0
		vi['GC'][5]=0
		vi['GC']['Active']=0
		vi['GC']['ID_GCA']=np.ones(meta[pNam]['Project']['Batch Size'][iBat])

	#--------------------------------------------------------------------------
	# Initialize output variables
	#--------------------------------------------------------------------------

	# Output variables dictionary
	vo={}

	# Get dimensions
	m=meta[pNam]['Project']['N Time']
	n=meta[pNam]['Project']['Batch Size'][iBat]
	o=meta['Core']['N Pools Eco']

	# Land cover / land use
	if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
		vo['LandCover']=np.tile(vi['lsat']['LandCover_Comp1_1800'],(m,1)).astype('int8')
		vo['LandUse']=np.tile(meta['LUT']['Derived']['lu_comp1']['No Designation'],(m,n)).astype('int8')
	else:
		vo['LandCover']=np.zeros((m,n),dtype='int8')
		vo['LandUse']=np.zeros((m,n),dtype='int8')

	# Stand age (i.e. time since stand-replacing disturbance)
	vo['A']=np.zeros((m,n))
	vo['A_Harvest']=np.nan*np.ones((m,n)) # This needs to be nans

	# Species composition
	#vo['Spc_ID']=np.zeros((m,n,6),dtype='int8')
	#vo['Spc_Percent']=np.zeros((m,n,6),dtype='int8')

	# Stemwood merch volume
	vo['V_MerchLive']=np.zeros((m,n))
	vo['V_MerchDead']=np.zeros((m,n))
	vo['V_MerchTotal']=np.zeros((m,n))

	# Stemwood merch volume sent to mill
	vo['V_ToMillMerchLive']=np.zeros((m,n))
	vo['V_ToMillMerchDead']=np.zeros((m,n))
	vo['V_ToMillMerchTotal']=np.zeros((m,n))
	vo['V_ToMillNonMerch']=np.zeros((m,n))

	# Piece size
	vo['LogSizeEnhancement']=np.zeros((m,n))

	# Carbon density of ecosystem (Mg C ha-1)
	# -> 3-D matrix: Time x Stand x Carbon pool
	vo['C_Eco_ByPool']=np.zeros((m,n,o))

	# Carbon density of products sector (Mg C ha-1)
	# -> 3-D matrix: Time x Stand x Carbon pool
	vo['C_Pro_ByPool']=np.zeros((m,n,meta['Core']['N Pools Pro']))

	# Aggregate pools (Mg C ha-1) (polulated upon export)
	vo['C_Biomass']=np.array([])
	vo['C_Stemwood']=np.array([])
	vo['C_Foliage']=np.array([])
	vo['C_Branch']=np.array([])
	vo['C_Bark']=np.array([])
	vo['C_Root']=np.array([])
	vo['C_Piles']=np.array([])
	vo['C_Litter']=np.array([])
	vo['C_DeadWood']=np.array([])
	vo['C_Soil']=np.array([])
	vo['C_Soil_OHorizon']=np.array([])
	vo['C_InUse']=np.array([])
	vo['C_DumpLandfill']=np.array([])
	vo['C_Buildings']=np.zeros((m,n))

	# Carbon flux densities (Mg C ha-1 yr-1)
	vo['C_G_Gross_ByPool']=np.zeros((m,n,o))
	vo['C_G_Net_Reg_ByPool']=np.zeros((m,n,o))
	vo['C_M_Reg_ByPool']=np.zeros((m,n,o))
	vo['C_LF_ByPool']=np.zeros((m,n,o))
	vo['C_RH_ByPool']=np.zeros((m,n,o))

	# Aggregate pools (Mg C ha-1) (polulated upon export)
	vo['C_G_Gross']=np.array([])
	vo['C_M_Reg']=np.array([])
	vo['C_M_Dist']=np.zeros((m,n))
	vo['C_G_Net_Reg']=np.array([])
	vo['C_G_Net']=np.array([])
	vo['C_LF']=np.array([])
	vo['C_RH']=np.array([])

	vo['C_M_DistByAgent']={}
	vo['C_M_DistByAgentPct']={}
	for k in meta['LUT']['Event'].keys():
		id=meta['LUT']['Event'][k]
		vo['C_M_DistByAgent'][id]=np.zeros((m,n))
		vo['C_M_DistByAgentPct'][id]=np.zeros((m,n))

	# Keep track of carbon transfers (for economics)
	vo['C_Felled']=np.zeros((m,n))
	vo['C_FelledMerch']=np.zeros((m,n))
	vo['C_FelledRoots']=np.zeros((m,n))
	vo['C_ToDOM']=np.zeros((m,n))
	vo['C_ToPile']=np.zeros((m,n))
	vo['C_ToPileMerch']=np.zeros((m,n))
	vo['C_ToMillMerch']=np.zeros((m,n))
	vo['C_ToMillNonMerch']=np.zeros((m,n))
	vo['C_ToMillSnagStem']=np.zeros((m,n))
	vo['C_ToSlashpileBurnTot']=np.zeros((m,n))
	vo['C_ToSlashpileBurnNonMerch']=np.zeros((m,n))
	vo['C_ToLumber']=np.zeros((m,n))
	vo['C_ToPlywood']=np.zeros((m,n))
	vo['C_ToOSB']=np.zeros((m,n))
	vo['C_ToMDF']=np.zeros((m,n))
	vo['C_ToPaper']=np.zeros((m,n))
	vo['C_ToPowerFacilityDom']=np.zeros((m,n))
	vo['C_ToPowerFacilityExport']=np.zeros((m,n))
	vo['C_ToPowerGrid']=np.zeros((m,n))
	vo['C_ToPelletExport']=np.zeros((m,n))
	vo['C_ToPelletDomRNG']=np.zeros((m,n))
	vo['C_ToPelletDomGrid']=np.zeros((m,n))
	vo['C_ToFirewoodDom']=np.zeros((m,n))
	vo['C_ToFirewoodExport']=np.zeros((m,n))
	vo['C_ToLogExport']=np.zeros((m,n))

	# Emissions from wildfire and open burning (will be deleted upon export)
	vo['C_E_OpenBurningAsCO2']=np.zeros((m,n))
	vo['C_E_OpenBurningAsCH4']=np.zeros((m,n))
	vo['C_E_OpenBurningAsCO']=np.zeros((m,n))
	vo['C_E_OpenBurningAsN2O']=np.zeros((m,n))

	vo['C_E_WildfireAsCO2']=np.zeros((m,n))
	vo['C_E_WildfireAsCH4']=np.zeros((m,n))
	vo['C_E_WildfireAsCO']=np.zeros((m,n))
	vo['C_E_WildfireAsN2O']=np.zeros((m,n))

	vo['C_FromCoal']=np.zeros((m,n))
	vo['C_FromOil']=np.zeros((m,n))
	vo['C_FromGas']=np.zeros((m,n))
	vo['C_FromLimestone']=np.zeros((m,n))
	vo['C_Coal']=np.zeros((m,n))
	vo['C_Oil']=np.zeros((m,n))
	vo['C_Gas']=np.zeros((m,n))
	vo['C_Limestone']=np.zeros((m,n))

	# Emissions, domestic, forest sector, ecosystem
	vo['E_Dom_FS_NPP']=np.zeros((m,n))
	vo['E_Dom_FS_RH']=np.zeros((m,n))
	vo['E_Dom_FS_NEE']=np.zeros((m,n))
	vo['E_Dom_FS_Wildfire']=np.zeros((m,n))
	vo['E_Dom_FS_OpenBurning']=np.zeros((m,n))
	vo['E_Dom_FS_Denit']=np.zeros((m,n))
	vo['E_Dom_FS_Other']=np.zeros((m,n))
	vo['E_Dom_FS_HWP']=np.zeros((m,n))

	# Emissions, international, forest sector, ecosystem
	vo['E_Int_FS_NPP']=np.zeros((m,n))
	vo['E_Int_FS_RH']=np.zeros((m,n))
	vo['E_Int_FS_NEE']=np.zeros((m,n))
	vo['E_Int_FS_Wildfire']=np.zeros((m,n))
	vo['E_Int_FS_OpenBurning']=np.zeros((m,n))
	vo['E_Int_FS_Denit']=np.zeros((m,n))
	vo['E_Int_FS_Other']=np.zeros((m,n))
	vo['E_Int_FS_HWP']=np.zeros((m,n))

	# Emissions, domestic, energy stationary combustion, bioenergy
	vo['E_ESC_Bioenergy']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyPowerFacility']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyPowerGrid']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyPelletGrid']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyPelletRNG']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyPelletHydrogen']=np.zeros((m,n))
	vo['E_Dom_ESC_BioenergyFirewood']=np.zeros((m,n))

	# Emissions, internationl, energy stationary combustion, bioenergy
	vo['E_Int_ESC_BioenergyPowerFacility']=np.zeros((m,n))
	vo['E_Int_ESC_BioenergyPelletGrid']=np.zeros((m,n))
	vo['E_Int_ESC_BioenergyPelletRNG']=np.zeros((m,n))
	vo['E_Int_ESC_BioenergyPelletHydrogen']=np.zeros((m,n))
	vo['E_Int_ESC_BioenergyFirewood']=np.zeros((m,n))

	# Emissions, domestic, energy stationary combustion, forestry operations
	vo['E_Dom_ESC_ForOpsBurnCoal']=np.zeros((m,n))
	vo['E_Dom_ESC_ForOpsBurnOil']=np.zeros((m,n))
	vo['E_Dom_ESC_ForOpsBurnGas']=np.zeros((m,n))

	# Emissions, domestic, energy transporation, forestry operations
	vo['E_Dom_ET_ForOpsBurnCoal']=np.zeros((m,n))
	vo['E_Dom_ET_ForOpsBurnOil']=np.zeros((m,n))
	vo['E_Dom_ET_ForOpsBurnGas']=np.zeros((m,n))

	# Emissions, domestic, industrial Produciton and Product Use, forestry operations
	vo['E_Dom_IPPU_ForOpsBurningCoal']=np.zeros((m,n))
	vo['E_Dom_IPPU_ForOpsBurningOil']=np.zeros((m,n))
	vo['E_Dom_IPPU_ForOpsBurningGas']=np.zeros((m,n))

	if meta[pNam]['Project']['Biomass Module']=='Sawtooth':

		meta['Core']['Sawtooth']={}
		meta['Core']['Sawtooth']['DBH Classes']=np.linspace(0,100,50)
		vo['DBH_Class']=np.zeros((m,n,meta['Core']['Sawtooth']['DBH Classes'].size))

		# Stand density
		vo['N']=np.zeros((m,n))

		# Change in stand density (stems ha-1 yr-1)
		vo['N_R']=np.zeros((m,n))
		vo['N_M_Tot']=np.zeros((m,n))
		vo['N_M_Reg']=np.zeros((m,n))

		# Mean of tree attributes
		vo['TreeMean_A']=np.zeros((m,n))
		vo['TreeMean_H']=np.zeros((m,n))
		vo['TreeMean_D']=np.zeros((m,n))
		vo['TreeMean_Csw']=np.zeros((m,n))
		vo['TreeMean_Csw_G']=np.zeros((m,n))

		# Needed to supply affected carbon to disturbance module
		#vo['C_M_Tot']=np.zeros((m,n,o))

		# Needed to track dead merch volume
		vo['V_Merch_M']=np.zeros((m,n))

	# Atmosphere
	vo['Atm_CO2_In']=np.zeros((m,n))
	vo['Atm_CO2_Out']=np.zeros((m,n))
	vo['Atm_CH4_In']=np.zeros((m,n))
	vo['Atm_CH4_Out']=np.zeros((m,n))
	vo['Atm_N2O_In']=np.zeros((m,n))
	vo['Atm_N2O_Out']=np.zeros((m,n))

	#--------------------------------------------------------------------------
	# Configure batch-specific setttings
	#--------------------------------------------------------------------------

	# Nutrient application response yearly counter
	meta['Modules']['NutrientApp']['ResponseCounter']=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])
	meta['Modules']['NutrientApp']['ResponseCounterContinuous']=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])

	# Initialize flag for fixing negative net growth. When TIPSY yields negative
	# net growth, the fluxes of gross growth and mortality need adjustment.
	# This flag helps achieve that.
	meta[pNam]['Project']['FlagNegNetGrowth']=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])
	meta[pNam]['Project']['G_Net_PriorToBreakup']=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],7))

	#--------------------------------------------------------------------------
	# Initialize a log book that will record various diagnostics, warnings,
	# and error flags
	#--------------------------------------------------------------------------

	meta['Logbook']=list()

	return meta,vi,vo

#%% Import parameters
def PrepareParametersForBatch(meta,pNam,vi,iEns,iBat,iScn):

	#--------------------------------------------------------------------------
	# Populate final parameters with initial best estimates
	#--------------------------------------------------------------------------
	meta['Param']['BEV']=copy.deepcopy(meta['Param']['BE'])

	#--------------------------------------------------------------------------
	# Add error variance to parameters
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Uncertainty Status Biomass Turnover']=='On':
		for k in meta['Param']['BE']['BiomassTurnover'].keys():
			meta['Param']['BEV']['BiomassTurnover'][k]=meta['Param']['By Ensemble'][iEns]['BiomassTurnover'][k]

	if meta[pNam]['Project']['Uncertainty Status Decomposition']=='On':
		for k in meta['Param']['BE']['Decomposition'].keys():
			meta['Param']['BEV']['Decomposition'][k]=meta['Param']['By Ensemble'][iEns]['Decomposition'][k]

	if meta[pNam]['Project']['Uncertainty Status Inter Pool Fluxes']=='On':
		for k in meta['Param']['BE']['InterPoolFluxes'].keys():
			meta['Param']['BEV']['InterPoolFluxes'][k]=meta['Param']['By Ensemble'][iEns]['InterPoolFluxes'][k]

	if meta[pNam]['Project']['Uncertainty Status Harvest Utilization']=='On':
		EventList=['Harvest','Harvest Salvage']
		VariableList=['BiomassMerch_Removed','BiomassNonMerch_Removed','Snags_Removed', \
					  'BiomassMerch_Piled','BiomassNonMerch_Piled','Snags_Piled', \
					  'BiomassMerch_LeftOnSite','BiomassNonMerch_LeftOnSite','Snags_LeftOnSite']
		for Event in EventList:
			ID_Type=meta['LUT']['Event'][Event]
			for Variable in VariableList:
				meta['Param']['BE']['Events'][ID_Type][Variable]=meta['Param']['By Ensemble'][iEns]['Event'][ID_Type][Variable]

	if meta[pNam]['Project']['Uncertainty Status Substitution']=='On':
		for k in meta['Param']['By Ensemble'][iEns]['Substitution'].keys():
			meta['Param']['BEV']['Substitution'][k]=meta['Param']['By Ensemble'][iEns]['Substitution'][k]

	if meta[pNam]['Project']['Uncertainty Status Nutrient Application']=='On':
		if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
			for k in meta['Param']['BE']['NutrientApp'].keys():
				meta['Param']['BEV']['NutrientApp'][k]=meta['Param']['By Ensemble'][iEns]['NutrientApp'][k]

	#--------------------------------------------------------------------------
	# Biomass allometry (stand level)
	#--------------------------------------------------------------------------
	u=np.unique(vi['lsat']['ID_BGCZ'].flatten())
	MaritimeZones=['CDF','CWH','ICH']
	for k in meta['Param']['BE']['BiomassAllometrySL'].keys():
		if k=='Region': continue
		meta['Param']['BEV']['BiomassAllometrySL'][k]=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])
		for iU in range(u.size):
			cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			ind=np.where(vi['lsat']['ID_BGCZ'].flatten()==u[iU])[0]
			if np.isin(cd,MaritimeZones)==True:
				meta['Param']['BEV']['BiomassAllometrySL'][k][ind]=meta['Param']['BE']['BiomassAllometrySL'][k][0]
			else:
				meta['Param']['BEV']['BiomassAllometrySL'][k][ind]=meta['Param']['BE']['BiomassAllometrySL'][k][1]

	#--------------------------------------------------------------------------
	# Fate of Felled Material
	#--------------------------------------------------------------------------
	# Initialize
	meta['Param']['BEV']['FelledFate']={}
	for k in meta['Param']['BE']['FelledFate']['BaseCase']['Coast'].keys():
		meta['Param']['BEV']['FelledFate'][k]=np.zeros( (meta['Param']['BE']['FelledFate']['Year'].size,meta[pNam]['Project']['Batch Size'][iBat]) )

	# Papulate batch-specific parameters
	if meta[pNam]['Project']['Scenario Source']=='Portfolio':

		# Isolate felled fate scenario names within this batch
		Scenario=meta[pNam]['Project']['Portfolio']['Felled Fate Scenario'][meta[pNam]['Project']['indBat']]

		# Isolate region
		Region=meta[pNam]['Project']['Portfolio']['Region Code'][meta[pNam]['Project']['indBat']]

		# Unique scenario and region (must be converted to string)
		SR=np.column_stack((Scenario,Region))
		SR=SR.astype(str)
		u=np.unique(SR,axis=0)

		for iU in range(u.shape[0]):

			# Index to each scenario
			scn=u[iU,0]
			reg=u[iU,1]
			ind=np.where( (Scenario==scn) & (Region==reg) )[0]

			for k in meta['Param']['BE']['FelledFate'][scn][reg].keys():
				x=meta['Param']['BE']['FelledFate'][scn][reg][k]
				for i in range(ind.size):
					meta['Param']['BEV']['FelledFate'][k][:,ind[i]]=x

	elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':

		Scenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Felled Fate Change Scenario']
		Region=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Region Code']
		for k in meta['Param']['BE']['FelledFate'][Scenario][Region].keys():
			x=meta['Param']['BE']['FelledFate'][Scenario][Region][k]
			x=np.reshape(x,(-1,1))
			x=np.tile(x,(1,meta[pNam]['Project']['Batch Size'][iBat]))
			meta['Param']['BEV']['FelledFate'][k]=x

	elif meta[pNam]['Project']['Scenario Source']=='Script':

		# Isolate felled fate scenario names within this batch
		HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Felled Fate Historical Regime']
		ChangeScenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Felled Fate Change Scenario']

		if HistoricalRegime=='Regional Defaults':

			for reg in meta['LUT']['Region'].keys():
				ind=np.where( vi['lsat']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
				for k in meta['Param']['BE']['FelledFate'][ChangeScenario][reg].keys():
					x=meta['Param']['BE']['FelledFate'][ChangeScenario][reg][k]
					for i in range(ind.size):
						meta['Param']['BEV']['FelledFate'][k][:,ind[i]]=x
		else:

			for k in meta['Param']['BE']['FelledFate'][ChangeScenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['FelledFate'][ChangeScenario][HistoricalRegime][k]
				for i in range(meta[pNam]['Project']['indBat'].size):
					meta['Param']['BEV']['FelledFate'][k][:,i]=x

		# # Override historical regime parameters for stands that have land use = energy production
		# if meta[pNam]['Scenario'][iScn]['Land Surface Scenario']!='Off':
		#	 iEnergy=np.where(vi['lsat']['LSC']['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
		#	 if iEnergy[0].size>0:
		#		 for k in meta['Param']['BEV']['FelledFate'].keys():
		#			 meta['Param']['BEV']['FelledFate'][k][iEnergy]=meta['Param']['BE']['FelledFate'][ChangeScenario]['Energy Production'][k][0]

	#--------------------------------------------------------------------------
	# Removed Fate
	#--------------------------------------------------------------------------
	# Initialize
	meta['Param']['BEV']['RemovedFate']={}
	for k in meta['Param']['BE']['RemovedFate']['BaseCase']['Coast'].keys():
		meta['Param']['BEV']['RemovedFate'][k]=np.zeros( (meta['Param']['BE']['RemovedFate']['Year'].size,meta[pNam]['Project']['Batch Size'][iBat]) )

	# Papulate batch-specific parameters
	if meta[pNam]['Project']['Scenario Source']=='Portfolio':

		# Isolate Removal Fate scenario names within this batch
		Scenario=meta[pNam]['Project']['Portfolio']['Removed Fate Change Scenario'][meta[pNam]['Project']['indBat']]

		# Isolate region
		Region=meta[pNam]['Project']['Portfolio']['Region Code'][meta[pNam]['Project']['indBat']]

		# Unique scenario and region (must be converted to string)
		SR=np.column_stack( (Scenario,Region) )
		SR=SR.astype(str)
		u=np.unique(SR,axis=0)

		for iU in range(u.shape[0]):

			# Index to each scenario
			scn=u[iU,0]
			reg=u[iU,1]
			ind=np.where( (Scenario==scn) & (Region==reg) )[0]

			for k in meta['Param']['BE']['RemovedFate'][scn][reg].keys():
				x=meta['Param']['BE']['RemovedFate'][scn][reg][k]
				for i in range(ind.size):
					meta['Param']['BEV']['RemovedFate'][k][:,ind[i]]=x

	elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':

		Scenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Removed Fate Change Scenario']
		Region=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Region Code']
		for k in meta['Param']['BE']['RemovedFate'][Scenario][Region].keys():
			x=meta['Param']['BE']['RemovedFate'][Scenario][Region][k]
			x=np.reshape(x,(-1,1))
			x=np.tile(x,(1,meta[pNam]['Project']['Batch Size'][iBat]))
			meta['Param']['BEV']['RemovedFate'][k]=x

	elif meta[pNam]['Project']['Scenario Source']=='Script':

		# Isolate Removal Fate scenario names within this batch
		HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Removed Fate Historical Regime']
		ChangeScenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Removed Fate Change Scenario']

		if HistoricalRegime=='Regional Defaults':

			for reg in meta['LUT']['Region'].keys():
				ind=np.where( vi['lsat']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
				for k in meta['Param']['BE']['RemovedFate'][ChangeScenario][reg].keys():
					x=meta['Param']['BE']['RemovedFate'][ChangeScenario][reg][k]
					for i in range(ind.size):
						meta['Param']['BEV']['RemovedFate'][k][:,ind[i]]=x

		else:

			for k in meta['Param']['BE']['RemovedFate'][ChangeScenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['RemovedFate'][ChangeScenario][HistoricalRegime][k]
				for i in range(meta[pNam]['Project']['indBat'].size):
					meta['Param']['BEV']['RemovedFate'][k][:,i]=x

		# # Override historical regime parameters for stands that have land use = energy production
		# if meta[pNam]['Scenario'][iScn]['Land Surface Scenario']!='Off':
		#	 iEnergy=np.where(vi['lsat']['LSC']['Use']==meta['LUT']['LSC']['Use']['Energy Production'])
		#	 if iEnergy[0].size>0:
		#		 for k in meta['Param']['BEV']['RemovedFate'].keys():
		#			 meta['Param']['BEV']['RemovedFate'][k][iEnergy]=meta['Param']['BE']['RemovedFate'][ChangeScenario]['Energy Production'][k][0]

	#--------------------------------------------------------------------------
	# Override mill transfers (for salvage logging)
	# Shift transfers to pulp in salvage logging projects
	#--------------------------------------------------------------------------
	if 'Salvage Mill Transfers' in meta[pNam]['Scenario'][iScn]:
		if meta[pNam]['Scenario'][iScn]['Salvage Mill Transfers']=='On':
			for k1 in meta['Param']['BE']['RemovedFate'].keys():
				if k1=='Year':
					continue
				for k2 in meta['Param']['BE']['RemovedFate'][k1].keys():
					y=meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedMerchToLumberMill'].copy()
					meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedMerchToLumberMill']=0.0*y
					meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedMerchToPulpMill']=meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedMerchToPulpMill']+1.0*y

					y=meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedSnagStemToLumberMill'].copy()
					meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedSnagStemToLumberMill']=0.0*y
					meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedSnagStemToPulpMill']=meta['Param']['BE']['RemovedFate'][k1][k2]['RemovedSnagStemToPulpMill']+1.0*y

	#--------------------------------------------------------------------------
	# Product Types
	#--------------------------------------------------------------------------
	# Initialize
	meta['Param']['BEV']['ProductTypes']={}
	for k in meta['Param']['BE']['ProductTypes']['BaseCase']['Coast'].keys():
		meta['Param']['BEV']['ProductTypes'][k]=np.zeros( (meta['Param']['BE']['ProductTypes']['Year'].size,meta[pNam]['Project']['Batch Size'][iBat]) )

	# Papulate batch-specific parameters
	if meta[pNam]['Project']['Scenario Source']=='Portfolio':

		# Isolate scenario names within this batch
		Scenario=meta[pNam]['Project']['Portfolio']['HWP End Use Scenario'][meta[pNam]['Project']['indBat']]

		# Isolate region
		Region=meta[pNam]['Project']['Portfolio']['Region Code'][meta[pNam]['Project']['indBat']]

		# Unique scenario and region (must be converted to string)
		SR=np.column_stack( (Scenario,Region) )
		SR=SR.astype(str)
		u=np.unique(SR,axis=0)
		for iU in range(u.shape[0]):
			# Index to each scenario
			scn=u[iU,0]
			reg=u[iU,1]
			ind=np.where( (Scenario==scn) & (Region==reg) )[0]
			for k in meta['Param']['BE']['ProductTypes'][scn][reg].keys():
				x=meta['Param']['BE']['ProductTypes'][scn][reg][k]
				x=np.reshape(x,(-1,1))
				x=np.tile(x,(1,ind.size))
				meta['Param']['BEV']['ProductTypes'][k][:,ind]=x

	elif meta[pNam]['Project']['Scenario Source']=='Spreadsheet':
		HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Product Type Historical Regime']
		ChangeScenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Product Type Change Scenario']
		if HistoricalRegime=='Regional Defaults':
			reg=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Region Code']
			for k in meta['Param']['BE']['ProductTypes'][ChangeScenario][reg].keys():
				x=meta['Param']['BE']['ProductTypes'][ChangeScenario][reg][k]
				for i in range(ind.size):
					meta['Param']['BEV']['ProductTypes'][k][:,ind[i]]=x
		else:
			for k in meta['Param']['BE']['ProductTypes'][ChangeScenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['ProductTypes'][ChangeScenario][HistoricalRegime][k]
				for i in range(meta[pNam]['Project']['indBat'].size):
					meta['Param']['BEV']['ProductTypes'][k][:,i]=x
	else:
		HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Product Type Historical Regime']
		ChangeScenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Product Type Change Scenario']
		if HistoricalRegime=='Regional Defaults':
			for reg in meta['LUT']['Region'].keys():
				ind=np.where( vi['lsat']['Region Code'][0,:]==meta['LUT']['Region'][reg] )[0]
				for k in meta['Param']['BE']['ProductTypes'][ChangeScenario][reg].keys():
					x=meta['Param']['BE']['ProductTypes'][ChangeScenario][reg][k]
					for i in range(ind.size):
						meta['Param']['BEV']['ProductTypes'][k][:,ind[i]]=x
		else:
			for k in meta['Param']['BE']['ProductTypes'][ChangeScenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['ProductTypes'][ChangeScenario][HistoricalRegime][k]
				for i in range(meta[pNam]['Project']['indBat'].size):
					meta['Param']['BEV']['ProductTypes'][k][:,i]=x

	#--------------------------------------------------------------------------
	# Product disposal
	#--------------------------------------------------------------------------
	# *** No differentiation at this time. ***

	#--------------------------------------------------------------------------
	# Growth enhancement factors
	# *** RETIRED AND MOVED TO BGC Zone-specific ***
	#--------------------------------------------------------------------------
# 	if meta[pNam]['Scenario'][iScn]['Gaia Status']=='On':
# 		meta['Param']['BEV']['GrowthModifier']={}
# 		ind=np.where(vi['lsat']['Region Code'][0,:]==meta['LUT']['Region']['Interior'])[0]
# 		for i in range(meta['Param']['Raw']['GrowthModifier']['Name'].size):
# 			k=meta['Param']['Raw']['GrowthModifier']['Name'][i]
# 			meta['Param']['BEV']['GrowthModifier'][k]=meta['Param']['Raw']['GrowthModifier']['Coast'][i]*np.ones(meta[pNam]['Project']['indBat'].size)
# 			meta['Param']['BEV']['GrowthModifier'][k][ind]=meta['Param']['Raw']['GrowthModifier']['Interior'][i]

	#--------------------------------------------------------------------------
	# BGC Zone-specific
	#--------------------------------------------------------------------------
	meta['Param']['BEV']['ByBGCZ']={}
	id=vi['lsat']['ID_BGCZ'][0,:]
	u=np.unique(id)
	meta['Param']['BEV']['ByBGCZ']['Disease Severity']=np.zeros(meta[pNam]['Project']['indBat'].size)
	meta['Param']['BEV']['ByBGCZ']['Frost Severity']=np.zeros(meta[pNam]['Project']['indBat'].size)
	meta['Param']['BEV']['ByBGCZ']['Wind Severity']=np.zeros(meta[pNam]['Project']['indBat'].size)
	for k in meta['Param']['BE']['ByBGCZ']['CWH'].keys():
		if k=='SpcLead': continue
		meta['Param']['BEV']['ByBGCZ'][k]=np.zeros(meta[pNam]['Project']['indBat'].size)
		for iU in range(u.size):
			cd=cbu.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[iU])[0]
			ind=np.where(id==u[iU])[0]
			meta['Param']['BEV']['ByBGCZ'][k][ind]=meta['Param']['BE']['ByBGCZ'][cd][k]

			# Severities
			meta['Param']['BEV']['ByBGCZ']['Disease Severity'][ind]=np.random.beta(meta['Param']['BE']['ByBGCZ'][cd]['Disease Sev Alpha'],meta['Param']['BE']['ByBGCZ'][cd]['Disease Sev Beta'],size=ind.size)
			meta['Param']['BEV']['ByBGCZ']['Frost Severity'][ind]=np.random.beta(meta['Param']['BE']['ByBGCZ'][cd]['Frost Sev Alpha'],meta['Param']['BE']['ByBGCZ'][cd]['Frost Sev Beta'],size=ind.size)
			meta['Param']['BEV']['ByBGCZ']['Wind Severity'][ind]=np.random.beta(meta['Param']['BE']['ByBGCZ'][cd]['Wind Sev Alpha'],meta['Param']['BE']['ByBGCZ'][cd]['Wind Sev Beta'],size=ind.size)

	return meta,vi

#%% Export simulation results
def ExportSimulation(meta,pNam,vi,vo,iScn,iEns,iBat,iEP,vo_full):

	# Save project metadata
	if (iBat==0) & (iScn==0):
		gu.opickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Metadata.pkl',meta)

	# Save a full copy of the first ensemble for re-use in big projects
	if (meta[pNam]['Project']['Early Record Recycling']=='On') & (iEns==0):

		it=np.where(meta[pNam]['Year']==meta[pNam]['Project']['Year Start Saving']-meta['Core']['Recycle Early Record Buffer']-1)[0]

		vo_full={}
		for k in vo.keys():

			if (k=='C_M_DistByAgent') | (k=='C_M_DistByAgentPct'):
				continue

			if vo[k].size==0:
				continue

			vo_full[k]=vo[k][it,:]

		#vo_full=copy.deepcopy(vo)

	# Extract parameters
	bB=meta['Param']['BEV']['Biophysical']

	# Isolate time period that will be saved to file
	it=np.where(meta[pNam]['Year']>=meta[pNam]['Project']['Year Start Saving'])[0]
	for k in vo:
		# Skip mortality summary
		if type(vo[k])==dict:
			continue
		# Skip variables that were initiated but not yet populated
		if vo[k].size==0:
			continue
		vo[k]=vo[k][it,:]

	# Mortality summaries
	for k in vo['C_M_DistByAgent'].keys():
		vo['C_M_DistByAgent'][k]=vo['C_M_DistByAgent'][k][it,:]/meta['Core']['Scale Factor C_M_DistByAgent']
		vo['C_M_DistByAgentPct'][k]=vo['C_M_DistByAgentPct'][k][it,:]/meta['Core']['Scale Factor C_M_DistByAgent']

		vo['C_M_DistByAgent'][k]=vo['C_M_DistByAgent'][k].astype('int16')
		vo['C_M_DistByAgentPct'][k]=vo['C_M_DistByAgentPct'][k].astype('int16')

	# Emissions from wildfire:

	# Carbon dioxide flux (tCO2/ha/yr)
	E_CO2=bB['Ratio_CO2_to_C']*vo['C_E_WildfireAsCO2']

	# Carbon monoxide flux (tCO/ha/yr)
	E_CO=bB['Ratio_CO_to_C']*vo['C_E_WildfireAsCO']

	# Methan flux *(tCH4/ha/yr)
	E_CH4=bB['Ratio_CH4_to_C']*vo['C_E_WildfireAsCH4']

	# Nitrous oxide flux (tN2O/ha/yr)
	E_N2O=bB['EF_N2O_fromCO2']*E_CO2

	# Convert fluxes to CO2e using global warming potential estimates
	CO2e_E_AsCO2=1*E_CO2
	CO2e_E_AsCH4=bB['GWP_CH4_AR5']*E_CH4
	CO2e_E_AsCO=bB['GWP_CO_AR5']*E_CO
	CO2e_E_AsN2O=bB['GWP_N2O_AR5']*E_N2O

	vo['E_Dom_FS_Wildfire']=CO2e_E_AsCO2+CO2e_E_AsCH4+CO2e_E_AsCO+CO2e_E_AsN2O

	# Add to atmosphere
	vo['Atm_CO2_In']=vo['Atm_CO2_In']+E_CO2
	vo['Atm_CH4_In']=vo['Atm_CH4_In']+E_CH4
	vo['Atm_N2O_In']=vo['Atm_N2O_In']+E_N2O

	# Emissions from open burning:

	# Carbon dioxide flux (tCO2/ha/yr)
	E_CO2=bB['Ratio_CO2_to_C']*vo['C_E_OpenBurningAsCO2']

	# Carbon monoxide flux (tCO/ha/yr)
	E_CO=bB['Ratio_CO_to_C']*vo['C_E_OpenBurningAsCO']

	# Methan flux *(tCH4/ha/yr)
	E_CH4=bB['Ratio_CH4_to_C']*vo['C_E_OpenBurningAsCH4']

	# Nitrous oxide flux (tN2O/ha/yr)
	E_N2O=bB['EF_N2O_fromCO2']*E_CO2

	# Convert fluxes to CO2e using global warming potential estimates
	CO2e_E_AsCO2=1*E_CO2
	CO2e_E_AsCH4=bB['GWP_CH4_AR5']*E_CH4
	CO2e_E_AsCO=bB['GWP_CO_AR5']*E_CO
	CO2e_E_AsN2O=bB['GWP_N2O_AR5']*E_N2O

	vo['E_Dom_FS_OpenBurning']=CO2e_E_AsCO2+CO2e_E_AsCH4+CO2e_E_AsCO+CO2e_E_AsN2O

	# Add to atmosphere
	vo['Atm_CO2_In']=vo['Atm_CO2_In']+E_CO2
	vo['Atm_CH4_In']=vo['Atm_CH4_In']+E_CH4
	vo['Atm_N2O_In']=vo['Atm_N2O_In']+E_N2O

	# Delete unnecessary fire emission variables
	del vo['C_E_WildfireAsCO2']
	del vo['C_E_WildfireAsCO']
	del vo['C_E_WildfireAsCH4']
	del vo['C_E_WildfireAsN2O']
	del vo['C_E_OpenBurningAsCO2']
	del vo['C_E_OpenBurningAsCO']
	del vo['C_E_OpenBurningAsCH4']
	del vo['C_E_OpenBurningAsN2O']

	# Calculate carbon content of dead wood, organic and mineral soil horizons following Shaw et al. (2017)
	vo['C_Soil_OHorizon']=vo['C_Eco_ByPool'][:,:,iEP['LitterVF']]+vo['C_Eco_ByPool'][:,:,iEP['LitterS']]

	# Aggregate pools
	vo['C_Biomass']=np.sum(vo['C_Eco_ByPool'][:,:,iEP['BiomassTotal']],axis=2)
	vo['C_Stemwood']=vo['C_Eco_ByPool'][:,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][:,:,iEP['StemNonMerch']]
	vo['C_Foliage']=vo['C_Eco_ByPool'][:,:,iEP['Foliage']]
	vo['C_Branch']=vo['C_Eco_ByPool'][:,:,iEP['Branch']]
	vo['C_Bark']=vo['C_Eco_ByPool'][:,:,iEP['Bark']]
	vo['C_Root']=vo['C_Eco_ByPool'][:,:,iEP['RootFine']]+vo['C_Eco_ByPool'][:,:,iEP['RootCoarse']]
	vo['C_Piles']=np.sum(vo['C_Eco_ByPool'][:,:,iEP['Piled']],axis=2)
	vo['C_Litter']=np.sum(vo['C_Eco_ByPool'][:,:,iEP['Litter']],axis=2)
	vo['C_DeadWood']=np.sum(vo['C_Eco_ByPool'][:,:,iEP['DeadWood']],axis=2)
	vo['C_Soil']=np.sum(vo['C_Eco_ByPool'][:,:,iEP['Soil']],axis=2)
	vo['C_InUse']=np.sum(vo['C_Pro_ByPool'][:,:,meta['Core']['iPP']['InUse']],axis=2)
	vo['C_Buildings']=np.sum(vo['C_Pro_ByPool'][:,:,meta['Core']['iPP']['Buildings'] ],axis=2)
	vo['C_DumpLandfill']=np.sum(vo['C_Pro_ByPool'][:,:,meta['Core']['iPP']['DumpLandfill']],axis=2)
	del vo['C_Eco_ByPool']
	del vo['C_Pro_ByPool']

	# Aggregate fluxes
	vo['C_G_Gross']=np.sum(vo['C_G_Gross_ByPool'],axis=2)
	vo['C_G_Net_Reg']=np.sum(vo['C_G_Net_Reg_ByPool'],axis=2)
	vo['C_M_Reg']=np.sum(vo['C_M_Reg_ByPool'],axis=2)
	vo['C_LF']=np.sum(vo['C_LF_ByPool'],axis=2)
	vo['C_RH']=np.sum(vo['C_RH_ByPool'],axis=2)
	del vo['C_G_Gross_ByPool']
	del vo['C_G_Net_Reg_ByPool']
	del vo['C_M_Reg_ByPool']
	del vo['C_LF_ByPool']
	del vo['C_RH_ByPool']

	# Apply scale factor and convert to integer
	for k in vo.keys():
		# Skip mortality summary by agent
		if (k=='C_M_DistByAgent') | (k=='C_M_DistByAgentPct'):
			continue
		if vo[k].size==0:
			continue
		if vo[k].dtype=='int8':
			continue

		# These variable needs a larger scale factor
		if (k=='E_Dom_FS_HWP') | (k=='E_ESC_Comb') | (k=='E_ET_Comb') | (k=='E_IPPU_Comb'):
			vo[k]=vo[k]/meta['Core']['Scale Factor Export Big']
		else:
			vo[k]=vo[k]/meta['Core']['Scale Factor Export Small']

		if np.max(vo[k])<32767:
			vo[k]=vo[k].astype('int16')
		else:
			vo[k]=vo[k].astype('int32')

	# Mortality summary by agent
	for k in vo['C_M_DistByAgent'].keys():
		idx=np.where(vo['C_M_DistByAgent'][k]>0)
		M=vo['C_M_DistByAgent'][k][idx].copy()
		vo['C_M_DistByAgent'][k]={}
		vo['C_M_DistByAgent'][k]['idx']=idx
		vo['C_M_DistByAgent'][k]['M']=M

	# Save data
	fout=meta['Paths'][pNam]['Output Scenario'][iScn] + '\\Data_Scn' + cbu.FixFileNum(iScn) + \
		'_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
	gu.opickle(fout,vo)

	# Save disturbance/management event chronology
	# If events are added on the fly, they will only be accessable if resaved.
	if (meta[pNam]['Scenario'][iScn]['Harvest Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Harvest Sim Future Status']=='On') | \
		(meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Status']=='On') | \
		(meta[pNam]['Scenario'][iScn]['IBM Sim Pre-obs Status']=='On') | (meta[pNam]['Scenario'][iScn]['IBM Sim Future Status']=='On') | \
		(meta[pNam]['Scenario'][iScn]['Disease Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Disease Sim Future Status']=='On') | \
		(meta[pNam]['Scenario'][iScn]['Wind Sim Historical Status']=='On') | (meta[pNam]['Scenario'][iScn]['Wind Sim Future Status']=='On'):
		vi['EC']['idx']=np.where(vi['EC']['ID Event Type']>0)
		vi['EC']['ID Event Type']=vi['EC']['ID Event Type'][vi['EC']['idx']]
		vi['EC']['Mortality Factor']=np.array(100*vi['EC']['Mortality Factor'][vi['EC']['idx']],dtype='int16')
		vi['EC']['Growth Factor']=vi['EC']['Growth Factor'][vi['EC']['idx']]
		vi['EC']['ID Growth Curve']=vi['EC']['ID Growth Curve'][vi['EC']['idx']]
		fout=meta['Paths'][pNam]['Input Scenario'][iScn] + '\\Modified_Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl'
		gu.opickle(fout,vi['EC'])

	# Save diagnostics (turned off for now)
	flg=0
	if flg==1:
		#if (iBat==0) & (iScn==0):

		cnams=['Value','Variable']
		df=pd.DataFrame(columns=cnams)

		# Model version
		df.loc[0]=[meta['Paths']['Model Code'],'Version']

		# Date
		now=datetime.datetime.now()
		df.loc[1]=[now,'Run date']

		# Simulation time
		#meta['t_Sim']=time.time()-meta['t_SimStart']
		#df.loc[2]=[meta['t_Sim']/60,'Simulation time (min)']

		# Check conservation of mass (stock change = NEBP)
		# *** only done for the first stand ***
		NPP=np.sum(vo['C_NPP'][:,0,0:7],axis=1)
		RH=np.sum(vo['C_RH_ByPool'][:,0,7:16],axis=1)
		E=vo['C_E_FireAsCO2'][:,0]+vo['C_E_FireAsCO'][:,0]+vo['C_E_FireAsCH4'][:,0]+vo['C_E_FireAsN2O'][:,0]
		R=vo['C_ToMillMerch'][:,0]+vo['C_ToMillNonMerch'][:,0]
		x=np.sum(vo['C_Eco_ByPool'][:,0,0:16],axis=1)
		y=np.sum(vo['C_Eco_ByPool'][0,0,0:16],axis=0)+np.cumsum(NPP-RH-R-E)
		D_abs=np.mean(np.abs(y-x))
		D_rel=np.mean(np.abs(y-x)/np.maximum(0.000001,x)*100)
		df.loc[3]=[D_abs,'Mean absolute difference between stock change and NECB (MgC/ha)']
		df.loc[4]=[D_rel,'Mean relative difference between stock change and NECB (%)']
		#plt.close(6)
		#plt.figure(6)
		#plt.plot(x)
		#plt.plot(y)

		# Year that slow soil carbon pool reaches dynamic equilibrium
		def runningMean(x, N):
			y=np.zeros((len(x),))
			for ctr in range(len(x)):
				y[ctr]=np.sum(x[ctr:(ctr+N)])
			return y/N

		#plt.plot(runningMean(vo['C_Eco_ByPool[:,0,15],150))
		#x=np.diff(vo['C_Eco_ByPool'][:,0,15])/vo['C_Eco_ByPool'][0:-1,0,15]*100
		#rmx=runningMean(x,150)
		#plt.figure(33)
		#plt.plot(rmx)
		#iDE=np.min(np.where(rmx<0.001))
		#yrDE=vi['tv'][iDE]
		#tDE=meta.N_t-iDE
		#df.loc[5]=[yrDE,'Equilibrium in SoilS reached at (year)']
		#df.loc[6]=[tDE,'Equlibrium in Soil S reached at (time before end of simulation)']

		# Save
		pthoutD=meta['Paths'][pNam]['Output Scenario'][iScn] + '\\Diagnostics.xlsx'
		writer=pd.ExcelWriter(pthoutD)
		df.to_excel(writer,'Sheet1')
		writer.save()

	return vo_full
