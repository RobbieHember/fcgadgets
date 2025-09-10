import numpy as np
from fcgadgets.cbrunner import cbrun_util as cbu

#%% Nutrient application effects
def NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,comp):

	# Exctract parameters
	if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
		bNA=meta['Param']['BEV']['NutrientApplication']
	else:
		bNA=meta['Param']['By Ensemble']['NutrientApplication']

	# Pull out index to applications
	iApp=meta['Modules']['NutrientApplication']['iApplication']

	if comp=='UpdateCounterForResponse':
		meta['Modules']['NutrientApplication']['ResponseCounter'][iApp]=meta['Modules']['NutrientApplication']['ResponseCounter'][iApp]+1

	elif comp=='UpdateCounterForHarvestRestriction':
		iHR=meta['Modules']['NutrientApplication']['iHarvestRestriction']
		meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction'][iHR]=meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction'][iHR]+1
		#print(np.max(meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction'][iHR]))

	elif (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner') & (comp=='AbovegroundNetGrowth'):

		#----------------------------------------------------------------------
		# Start response counter
		#----------------------------------------------------------------------

		meta['Modules']['NutrientApplication']['ResponseCounter'][iApp]=1
		meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction'][iApp]=1

		# Update log-size enhancement factor
		vo['LogSizeEnhancement'][iT:,iApp]=vo['LogSizeEnhancement'][iT:,iApp]+1

		# Get age vector for GC's
		A_gc=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

		# Extract active growth curves, apply scale factor
		GCA_SP=vi['GC']['Active'][:,iApp,:].copy()

		for iStand in range(iApp.size):
			AgeAtApplication=int(vo['A'][iT,iApp[iStand]])

			if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
				# Not from spreadsheet
				ResponseDuration=bNA['ResponseDuration']

				# Response ratio of stemwood (before wood density effect)
				rrS=bNA['r_Stemwood']

				# Wood density effect
				rr_wd=bNA['r_WoodDensity']
				f_wd=rr_wd-1
				f_wd=np.array([f_wd,f_wd,0,0,0])

				# Response relative to stemwood
				# Can't do roots here because they are not tracked in growth curve dict
				rXS=np.array([1,
						1,
						bNA['Ratio_Foliage_to_Stemwood'],
						bNA['Ratio_Branch_to_Stemwood'],
						bNA['Ratio_Bark_to_Stemwood']])

				# Biomass response ratios (after wood density effect accounted for)
				rr=1+(rrS-1)*rXS+f_wd

				# Append biomass response with volume response
				rr=np.append(rr,rrS)

			else:
				# From spreadsheet (stands act as ensembles!)
				ResponseDuration=bNA['ResponseDuration'][iStand]

				# Response ratio of stemwood (before wood density effect)
				rrS=bNA['r_Stemwood'][iStand]

				# Wood density effect
				rr_wd=bNA['r_WoodDensity'][iStand]
				f_wd=rr_wd-1
				f_wd=np.array([f_wd,f_wd,0,0,0])

				# Response relative to stemwood
				# Can't do roots here because they are not tracked in growth curve dict
				rXS=np.array([1,
						1,
						bNA['Ratio_Foliage_to_Stemwood'][iStand],
						bNA['Ratio_Branch_to_Stemwood'][iStand],
						bNA['Ratio_Bark_to_Stemwood'][iStand]])

				# Biomass response ratios (after wood density effect accounted for)
				rr=1+(rrS-1)*rXS+f_wd

				# *** SPECIAL ORDER FOR NUTRIENT MANAGEMENT DEMO - SPECIFICATION COMPARISON ***
				if (meta[pNam]['Project']['Code Subproject']=='Compare Specifications') & (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']=='Project NGS'):
					rr[1::]=1

				# Append biomass response with volume response
				rr=np.append(rr,rrS)

			# Index to the response period that will be alterred
			iResponse=np.where( (A_gc>=AgeAtApplication) & (A_gc<AgeAtApplication+ResponseDuration) )[0]

			# This will crash if an application occurs within 10 years of the maximum
			# TIPSY age curve (eg 200 years) -> adjust response period so that it
			# does not crash
			if iResponse.size==int(ResponseDuration):
				for iRR in range(rr.size):
					GCA_SP[iResponse,iStand,iRR]=rr[iRR]*GCA_SP[iResponse,iStand,iRR]
			else:
				em='Error: Nutrient application not implemented - stand age exceeds max age of growth curves.'
				print(em)
				#print(AgeAtApplication)

		# Repopulate in input variable dictionary
		vi['GC']['Active'][:,iApp,:]=GCA_SP

	elif (comp=='BelowgroundNetGrowth') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Adjust root net growth
		#----------------------------------------------------------------------

		iEP=meta['Core']['iEP']

		# Responses
		rrRC=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootC_to_Stemwood']
		rrRF=1+(bNA['r_Stemwood']-1)*bNA['Ratio_RootF_to_Stemwood']

		# *** SPECIAL ORDER STARTS ***
		# For nutrient management demo - comparison of specifications
		if (meta[pNam]['Project']['Code Subproject']=='Compare Specifications') & (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']=='Project NGS'):
			rrRC=1
			rrRF=1
		# *** SPECIAL ORDER ENDS ***

		# Calculate net growth of roots from change in pools
		#Gnet_RC=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]
		#Gnet_RF=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]

		# Remove the false stimulus that arises from estimating roots from AG
		# biomass
		vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootCoarse']]=0.84*vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootCoarse']]
		vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootFine']]=0.84*vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootFine']]
		#Gnet_RC[iApp]=0.65*(2-bNA['r_Stemwood'])*Gnet_RC[iApp]
		#Gnet_RF[iApp]=0.65*(2-bNA['r_Stemwood'])*Gnet_RF[iApp]

		# Stimulate root growth from N application
		#Gnet_RC[iApp]=rrRC*Gnet_RC[iApp]
		#Gnet_RF[iApp]=rrRF*Gnet_RF[iApp]

		# Recalculate current-year root biomass where stimulated by N application
		#vo['C_Eco_Pools'][iT,iApp,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT-1,iApp,iEP['RootCoarse']]+Gnet_RC[iApp]
		#vo['C_Eco_Pools'][iT,iApp,iEP['RootFine']]=vo['C_Eco_Pools'][iT-1,iApp,iEP['RootFine']]+Gnet_RF[iApp]

		# Populate net growth of root biomass
		vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootCoarse']]=rrRC*vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootCoarse']]
		vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootFine']]=rrRF*vo['C_G_Net_Reg_ByPool'][iT,iApp,iEP['RootFine']]
		#vo['C_G_Net_Reg'][iT,iApp,iEP['RootCoarse']]=Gnet_RC[iApp]
		#vo['C_G_Net_Reg'][iT,iApp,iEP['RootFine']]=Gnet_RF[iApp]

	elif (comp=='StopResponse') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Stop stimulation counter when it passes the response duration
		#----------------------------------------------------------------------

		iStop=np.where( meta['Modules']['NutrientApplication']['ResponseCounter']>bNA['ResponseDuration'] )[0]
		if iStop.size>0:
			meta['Modules']['NutrientApplication']['ResponseCounter'][iStop]=0

	elif (comp=='StopHarvestRestriction') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Stop harvest restrictions
		#----------------------------------------------------------------------
		if 'HarvestRestrictTSNA' in meta[pNam]['Scenario'][iScn].keys():
			iStop=np.where( meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction']>meta[pNam]['Scenario'][iScn]['HarvestRestrictTSNA'] )[0]
			if iStop.size>0:
				meta['Modules']['NutrientApplication']['ResponseCounterHarvestRestriction'][iStop]=0

	elif (comp=='Mortality') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Adjust mortality
		#----------------------------------------------------------------------

		if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
			vo['C_M_Reg_ByPool'][iT,iApp,0:7]=bNA['rPrime_TreeMortality']*vo['C_M_Reg_ByPool'][iT,iApp,0:7]
		else:
			vo['C_M_Reg_ByPool'][iT,iApp,0:7]=np.tile(bNA['rPrime_TreeMortality'][iApp],(7,1)).T*vo['C_M_Reg_ByPool'][iT,iApp,0:7]

	elif (comp=='Litterfall') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Adjust biomass turnover
		#----------------------------------------------------------------------

		iEP=meta['Core']['iEP']

		# *** SPECIAL ORDER FOR NUTRIENT MANAGEMENT DEMO - SPECIFICATION COMPARISON ***
		if (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']=='Project NGS') | (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']=='Project NGT'):
			pass
		else:
			vo['C_LF_ByPool'][iT,iApp,iEP['Foliage']]=bNA['rPrime_Litterfall']*vo['C_LF_ByPool'][iT,iApp,iEP['Foliage']]
			vo['C_LF_ByPool'][iT,iApp,iEP['Branch']]=bNA['rPrime_Litterfall']*vo['C_LF_ByPool'][iT,iApp,iEP['Branch']]
			vo['C_LF_ByPool'][iT,iApp,iEP['Bark']]=bNA['rPrime_Litterfall']*vo['C_LF_ByPool'][iT,iApp,iEP['Bark']]
			vo['C_LF_ByPool'][iT,iApp,iEP['RootCoarse']]=bNA['rPrime_TurnoverRootCoarse']*vo['C_LF_ByPool'][iT,iApp,iEP['RootCoarse']]
			vo['C_LF_ByPool'][iT,iApp,iEP['RootFine']]=bNA['rPrime_TurnoverRootFine']*vo['C_LF_ByPool'][iT,iApp,iEP['RootFine']]

	elif comp=='Emissions':

		#----------------------------------------------------------------------
		# Emissions from manufacture of ammonia and urea
		#----------------------------------------------------------------------

		# Urea dose (kgUrea/ha)
		DoseUrea=bNA['DoseUrea_Standard']

		# Nitrogen dose (kgN/ha)
		DoseN=bNA['Ratio_N_to_Urea']*DoseUrea

		# Emissions from production of ammonia (tCO2e/ha)
		E_ProdNH3_per_t=bNA['EmissionFromAmmoniaProduction_NRCAN']*bNA['TonneNH3PerTonneUrea']
		E_ProdNH3=E_ProdNH3_per_t*(DoseUrea/1000)

		# Emissions from production of Urea from ammonia (tCO2e/ha)
		MMBtu_per_app=(DoseUrea/1000)*bNA['UreaEnergyConsumption']
		therm_per_app=MMBtu_per_app/bNA['MMBtu_per_therm']
		E_ProdUrea=bNA['EmissionFromUreaProduction_per_therm']*therm_per_app

		vo['E_ForestryOps_EnergySC_Domestic_Gas'][iT,iApp]=vo['E_ForestryOps_EnergySC_Domestic_Gas'][iT,iApp] + \
			(E_ProdNH3+E_ProdUrea)

		#----------------------------------------------------------------------
		# Emissions from transportation (tCO2e/ha)
		#----------------------------------------------------------------------

		vo['E_ForestryOps_EnergyT_Domestic_Oil'][iT,iApp]=vo['E_ForestryOps_EnergyT_Domestic_Oil'][iT,iApp] + \
			(bNA['EmissionFromRailBargeTruck_Workbook']+bNA['EmissionFromHelicopter_SP10'])

		#----------------------------------------------------------------------
		# Volatilization
		#----------------------------------------------------------------------

		# CO2 emissions following application, Tier 1 approach, IPCC 2006, 11.4.1 (tCO2e/ha)
		# 0.2*(430/1000)*(1/0.27) = 0.32 tCO2e/ha
		#E_Vol=bNA['Ratio_C_to_Urea']*(DoseUrea/1000)*(1/meta['Param']['BEV']['Biophysical']['Ratio_C_to_CO2'])
		# Assume sequestration and volatilization of CO2 cancel out

		E_vol=0.3

		vo['E_ForestryOps_IPPU_Domestic_Gas'][iT,iApp]=vo['E_ForestryOps_IPPU_Domestic_Gas'][iT,iApp] - E_vol

		vo['E_Volat_ForestSector_Domestic'][iT,iApp]=vo['E_Volat_ForestSector_Domestic'][iT,iApp] + E_vol

		#----------------------------------------------------------------------
		# Denitrification: N2O emissions following application, Tier 1
		# approach, IPCC 2006, 11.4.1 (tCO2e/ha)
		#----------------------------------------------------------------------

		vo['E_Denit_ForestSector_Domestic'][iT,iApp]=vo['E_Denit_ForestSector_Domestic'][iT,iApp] + \
			(bNA['EmissionFactor_N2O_Jassaletal2008']*(DoseN/1000)*bNA['Ratio_N2OAsN_to_N2O']*meta['Param']['BEV']['Biophysical']['GWP_N2O'])

		#----------------------------------------------------------------------
		# Exterior area (volatilization/deposition effects)
		#----------------------------------------------------------------------

		flg=1

		if 'Nutrient Application Footprint Status' in meta[pNam]['Scenario'][meta[pNam]['iScn']]:
			if meta[pNam]['Scenario'][meta[pNam]['iScn']]['Nutrient Application Footprint Status']!='On':
				flg=0

		if (meta[pNam]['Project']['External Footprint Effect Status']=='On') & (flg==1):

			# Dose (kgN/ha)
			DoseN=bNA['DoseUrea_Standard']*bNA['Ratio_N_to_Urea']

			# Emissions (NH3-N ha-1)
			EmissionNH3_as_N=bNA['EA Fraction volatilized']*DoseN

			# Canopy uptake (kgN)
			CanopyUptakeN=EmissionNH3_as_N*bNA['EA Forest deposition fraction']*bNA['EA Leaf uptake fraction']

			# Root uptake (kgN)
			RootUptakeN=EmissionNH3_as_N*bNA['EA Forest deposition fraction']*bNA['EA Throughfall fraction']*bNA['EA RootUptakeFraction']

			# GHG benefit from canopy uptake (MgC/ha)
			GHG_Benefit_CanupyUptake_Coast=bNA['EA_NUEu_NGTT_Coast']*CanopyUptakeN/1000
			GHG_Benefit_CanupyUptake_Interior=CanopyUptakeN/1000*bNA['EA_NUEu_NGTT_Interior']

			# GHG benefit from root uptake (MgC/ha)
			GHG_Benefit_RootUptake_Coast=RootUptakeN/1000*bNA['EA_NUEa_NGTTD_Coast']
			GHG_Benefit_RootUptake_Interior=RootUptakeN/1000*bNA['EA_NUEa_NGTTD_Interior']

			# Total GHG benefit (MgC/ha)
			GHG_Benefit_Tot_Coast=GHG_Benefit_CanupyUptake_Coast+GHG_Benefit_RootUptake_Coast
			GHG_Benefit_Tot_Interior=GHG_Benefit_CanupyUptake_Interior+GHG_Benefit_RootUptake_Interior

			EA_GHG_Benefit=GHG_Benefit_Tot_Coast*bNA['EA Fraction of footprint coast']+GHG_Benefit_Tot_Interior*bNA['EA Fraction of footprint interior']

			# Annaul GHG benefit over response duration (tCO2e/ha/yr)
			#EA_GHG_Benefit=EA_GHG_Benefit/bNA['ResponseDuration']

			# Convert to CO2e (tCO2e/ha/yr)
			EA_GHG_Benefit=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*EA_GHG_Benefit

			# Subtract from ecosystem LULUCF emissions
			# *** it will crash in the last time step ***
			try:
				vo['E_Volat_ForestSector_Domestic'][iT+1,iApp]=vo['E_Volat_ForestSector_Domestic'][iT+1,iApp]-EA_GHG_Benefit
			except:
				pass

	elif (comp=='HeterotrophicRespiration') & (meta[pNam]['Project']['Nutrient Application Module']=='cbrunner'):

		#----------------------------------------------------------------------
		# Adjust rate of heterotrophic consumpion
		#----------------------------------------------------------------------

		# *** SPECIAL ORDER FOR NUTRIENT MANAGEMENT DEMO - SPECIFICATION COMPARISON ***
		if (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']!='Project NGS') & (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']!='Project NGT') & (meta[pNam]['Scenario'][ meta[pNam]['iScn'] ]['Scenario_CD']!='Project NGT+T'):
			rr=bNA['r_Decomp']
			meta[pNam]['Project']['HC']['LitterVF'][0,iApp]=rr*meta[pNam]['Project']['HC']['LitterVF'][0,iApp]
			meta[pNam]['Project']['HC']['LitterF'][0,iApp]=rr*meta[pNam]['Project']['HC']['LitterF'][0,iApp]
			meta[pNam]['Project']['HC']['LitterM'][0,iApp]=rr*meta[pNam]['Project']['HC']['LitterM'][0,iApp]
			meta[pNam]['Project']['HC']['LitterS'][0,iApp]=rr*meta[pNam]['Project']['HC']['LitterS'][0,iApp]
			meta[pNam]['Project']['HC']['SoilVF'][0,iApp]=rr*meta[pNam]['Project']['HC']['SoilVF'][0,iApp]
			meta[pNam]['Project']['HC']['SoilF'][0,iApp]=rr*meta[pNam]['Project']['HC']['SoilF'][0,iApp]
			meta[pNam]['Project']['HC']['SoilS'][0,iApp]=rr*meta[pNam]['Project']['HC']['SoilS'][0,iApp]

	return vi,vo,meta

#%% Simulate probability of harvesting on the fly
def ScheduleNutrientApplication(meta,pNam,vi,vo,iT,iScn,iEns,iBat):

	# Create a random number
	rn=np.random.random(meta[pNam]['Project']['Batch Size'][iBat])

	# Regional probabilities
	Po_Sat_Coast=meta['Param']['BEV']['NutrientApplication']['ProbOccSatAppFutureCoast']
	Po_Sat_Interior=meta['Param']['BEV']['NutrientApplication']['ProbOccSatAppFutureInterior']

	# Introduce time trend
	flg='Time trend'
	if flg=='Time trend':
		# Plot example
		flg_eg=0
		if flg_eg==1:
			tv=np.arange(1901,2100,1)
			Po_Sat_Coast=0.024
			Rate=-0.005*Po_Sat_Coast
			Po_Coast=np.maximum(0,Po_Sat_Coast+Rate*np.maximum(1,tv-2021))
			plt.close('all'); plt.plot(tv,Po_Coast,'r-',lw=1.5)

		Rate=-0.005*Po_Sat_Coast
		Po_Coast=np.maximum(0,Po_Sat_Coast+Rate*np.maximum(1,vi['tv']-2021))

		Rate=-0.005*Po_Sat_Interior
		Po_Interior=np.maximum(0,Po_Sat_Interior+Rate*np.maximum(1,vi['tv']-2021))
		#plt.close('all')
		#plt.plot(vi['tv'],Po_Interior,'r-',lw=1.5 )
	else:
		# Contstant
		Po_Coast=Po_Sat_Coast*np.ones(vi['tv'].size)
		Po_Interior=Po_Sat_Interior*np.ones(vi['tv'].size)

	# Time-dependent models to compensate for aging forests (paper)
	#Po_Coast=np.maximum(Po_Sat_Coast,Po_Sat_Coast+0.0005*np.maximum(1,vi['tv']-2021) )
	#Po_Interior=np.maximum(Po_Sat_Interior,Po_Sat_Interior+0.0005*np.maximum(1,vi['tv']-2021) )

	# Find eligible coastal stands to fertilize
	indS_Coast=np.where( (meta['Modules']['NutrientApplication']['ResponseCounter']==0) & \
			(vo['A'][iT,:]>=30) & \
			(vo['A'][iT,:]<=75) & \
			(rn<Po_Coast[iT]) & \
			(vo['V_MerchLive'][iT,:]>10) & \
			(np.isin(vi['lsat']['ID_BGCZ'][0,:],meta['Modules']['NutrientApplication']['BGC Zone Exclusion ID'])==False) & \
			(np.isin(vi['lsat']['ID_BGCZ'][0,:],meta['Modules']['NutrientApplication']['Coastal Zones ID'])==True) )[0]

	indS_Interior=np.where( (meta['Modules']['NutrientApplication']['ResponseCounter']==0) & \
			(vo['A'][iT,:]>=30) & \
			(vo['A'][iT,:]<=75) & \
			(rn<Po_Interior[iT]) & \
			(vo['V_MerchLive'][iT,:]>10) & \
			(np.isin(vi['lsat']['ID_BGCZ'][0,:],meta['Modules']['NutrientApplication']['BGC Zone Exclusion ID'])==False) & \
			(np.isin(vi['lsat']['ID_BGCZ'][0,:],meta['Modules']['NutrientApplication']['Coastal Zones ID'])==False) )[0]

	indS=np.append(indS_Coast,indS_Interior)

	if indS.size>0:
		for i in range(indS.size):
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
			if iAvailable.size>0:
				iE=iAvailable[0]
				vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Nutrient App Aerial']
				vi['EC']['Mortality Factor'][iT,indS[i],iE]=0
				vi['EC']['Growth Factor'][iT,indS[i],iE]=0
				#vi['EC']['ID Growth Curve'][iT,indS[i],iE]=np.max(vi['EC']['ID Growth Curve'][0:iT,indS[i],:])

	return vi
