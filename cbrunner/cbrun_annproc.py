#%% Import python modules
import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_general as gu
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.hardhat.nutrient_application as napp
import fcgadgets.hardhat.trenchfoot as tft
import fcgadgets.taz.aspatial_stat_models as asm
import fcgadgets.gaia.gaia_util as gaia
import warnings
import time
from scipy import stats

#%% Biomass dynamics
def TreeBiomassDynamicsFromGYModel(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

	# Update stand age
	vo['A'][iT,:]=vo['A'][iT-1,:]+1

	#--------------------------------------------------------------------------
	# Net growth of aboveground biomass
	#--------------------------------------------------------------------------

	# Index to growth curves at current age of stand
	iAge=np.minimum(vo['A'][iT,:],meta['Modules']['GYM']['BatchTIPSY Maximum Age'])-1

	# Convert to integer
	iAge=iAge.astype(int)

	# Extract net growth
	# Notes: I tried to make this faster using unravel but it was way slower
	NetGrowth=np.zeros( (meta[pNam]['Project']['Batch Size'][iBat],6) )
	for iS in range(meta[pNam]['Project']['Batch Size'][iBat]):
		NetGrowth[iS,:]=vi['GC']['Active'][iAge[iS],iS,:].copy()

	# Modify GY module estimates to represent variability in tree growth
	if meta[pNam]['Scenario'][iScn]['Gaia Status']=='On':
		NetGrowth=gaia.GYModelModifier(meta,pNam,NetGrowth,vo,iT)

	# Net growth of total stemwood
	Gnet_Stem=NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['StemMerch']]+NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['StemNonMerch']]

	# Net growth of foliage
	NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['Foliage']]=Gnet_Stem*(meta['Param']['BEV']['BiomassAllometrySL']['Gf1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gf2']-meta['Param']['BEV']['BiomassAllometrySL']['Gf1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gf3']*vo['A'][iT,:]))

	# Net growth of branches
	NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['Branch']]=Gnet_Stem*(meta['Param']['BEV']['BiomassAllometrySL']['Gbr1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbr2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbr1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbr3']*vo['A'][iT,:]))

	# Net growth of bark
	NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['Bark']]=Gnet_Stem*(meta['Param']['BEV']['BiomassAllometrySL']['Gbk1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbk2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbk1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbk3']*vo['A'][iT,:]))

	# Add net growth to output variable structure
	# Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions
	# of C_G_Net - don't change it.
	vo['C_G_Net_Reg'][iT,:,0:5]=NetGrowth[:,0:5]

	# Total net growth of root biomass (Li et al. 2003, Eq. 4)
	G_Net_Root_Total=0.22*np.sum(vo['C_G_Net_Reg'][iT,:,0:5],axis=1)

	# Fine root fraction should decline with total biomass, but estimating it based on
	# size causes a lot of problems. Think about creating new equation.
	vo['C_G_Net_Reg'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
	vo['C_G_Net_Reg'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total

	#--------------------------------------------------------------------------
	# Nutrient application effects to net growth
	#--------------------------------------------------------------------------

	# Index to stands that are stimulated by nutrient application
	meta['Modules']['NutrientApp']['iApplication']=np.where( meta['Modules']['NutrientApp']['ResponseCounter']>0 )[0]

	# Adjust N application response counter
	if meta['Modules']['NutrientApp']['iApplication'].size>0:
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'UpdateCounter')

	# Adjust root net growth
	if meta['Modules']['NutrientApp']['iApplication'].size>0:
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'BelowgroundNetGrowth')

	#--------------------------------------------------------------------------
	# Add net growth to biomass pools
	#--------------------------------------------------------------------------

	# Stemwood
	vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemMerch']]+vo['C_G_Net_Reg'][iT,:,iEP['StemMerch']])
	vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['StemNonMerch']]+vo['C_G_Net_Reg'][iT,:,iEP['StemNonMerch']])

	# Foliage
	vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Foliage']]+vo['C_G_Net_Reg'][iT,:,iEP['Foliage']])

	# Branches
	vo['C_Eco_Pools'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Branch']]+vo['C_G_Net_Reg'][iT,:,iEP['Branch']])

	# Bark
	vo['C_Eco_Pools'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['Bark']]+vo['C_G_Net_Reg'][iT,:,iEP['Bark']])

	# Coarse roots
	vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootCoarse']]+vo['C_G_Net_Reg'][iT,:,iEP['RootCoarse']])

	# Fine roots
	vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=np.maximum(0,vo['C_Eco_Pools'][iT-1,:,iEP['RootFine']]+vo['C_G_Net_Reg'][iT,:,iEP['RootFine']])

	#--------------------------------------------------------------------------
	# Add net growth to live merch volume and update other volume variables
	# Notes:
	# 1) Regualr mortality to be added to dead volume in DOM function ***
	# 2) Total volume updated at the start of disturfance event function.
	#--------------------------------------------------------------------------

	# Update live stemwood merchantable volume
	vo['V_MerchLive'][iT,:]=vo['V_MerchLive'][iT-1,:]+NetGrowth[:,meta['Modules']['GYM']['GC Input Indices']['StemMerchV']]

	# Update dead stemwood merchantable volume
	vo['V_MerchDead'][iT,:]=vo['V_MerchDead'][iT-1,:]

	#--------------------------------------------------------------------------
	# Biomass turnover
	#--------------------------------------------------------------------------

	# Biomass loss due to regular mortality
	bBT=meta['Param']['BEV']['BiomassTurnover']
	fM=bBT['Mreg0']+(bBT['Mreg1']-bBT['Mreg0'])*(1/(1+np.exp(-bBT['MregShape']*(vo['A'][iT,:]-bBT['MregInflect']))))
	vo['C_M_Reg'][iT,:,0:7]=np.tile(fM,(7,1)).T*vo['C_Eco_Pools'][iT,:,0:7]

	# Adjust mortality to account for N application response
	if meta['Modules']['NutrientApp']['iApplication'].size>0:
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Mortality')

	#--------------------------------------------------------------------------
	# If TIPSY indicates negative net growth (e.g., alder break-up), revise so
	# that net growth during the period of negative values equals the value in the
	# previous timestep.
	#--------------------------------------------------------------------------

	# Keep track of what regular mortality would be before it is affected by
	# correction for catastrophic mortality
	C_M_Reg=vo['C_M_Reg'][iT,:,:].copy()

	flg=0
	if flg==1:

		# Define a threshold level of negative growth so it isn't triggered by noise
		NegGrowthThreshold=-0.25
		#NegGrowthThreshold=-1e6

		# Find stands with negative net growth
		iNegNetG=np.where(vo['C_G_Net_Reg'][iT,:,0]<NegGrowthThreshold)[0]

		if iNegNetG.size>0:

			# If it is the first instance of negative net growth: 91) record net growth
			# of the preceeding timestep and (2) set flag = 1.
			iSwitchFlag=np.where(meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG]==0)[0]
			meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG[iSwitchFlag]]=1
			meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG[iSwitchFlag],0:7]=vo['C_G_Net_Reg'][iT-1,iNegNetG[iSwitchFlag],0:7]

			d=vo['C_G_Net_Reg'][iT,iNegNetG,0:7]-meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG,:]

			CToTransfer=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],7))
			CToTransfer[iNegNetG,:]=-1*d

			vo['C_G_Net_Reg'][iT,:,0:7]=vo['C_G_Net_Reg'][iT,:,0:7]+CToTransfer
			vo['C_M_Reg'][iT,:,0:7]=vo['C_M_Reg'][iT,:,0:7]+CToTransfer

			# # Logbook entry
			# for i in range(iNegNetG.size):
			#	 txt='Scenario:' + str(iScn) + ', Stand:' + str(iNegNetG[i]) + ', Time:' + str(iT) + ', Negative net growth received from GY model is being adjusted.'
			#	 meta['Logbook'].append(txt)

	#--------------------------------------------------------------------------
	# Litterfall
	#--------------------------------------------------------------------------

	# Setting turnover as a function of age will decouple NPP from net growth.
	Aref=meta['Param']['BEV']['BiomassTurnover']['BiomassTurnoverAgeRef']
	fA=-meta['Param']['BEV']['BiomassTurnover']['BiomassTurnoverAgeDependence']

	# Calculate foliage biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Foliage'])
	vo['C_LF'][iT,:,iEP['Foliage']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]

	# Calculate branch biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Branch'])
	vo['C_LF'][iT,:,iEP['Branch']]=tr*vo['C_Eco_Pools'][iT,:,iEP['Branch']]

	# Calculate bark biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Bark'])
	vo['C_LF'][iT,:,iEP['Bark']]=meta['Param']['BEV']['BiomassTurnover']['Bark']*vo['C_Eco_Pools'][iT,:,iEP['Bark']]

	# Calculate coarse root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootCoarse'])
	vo['C_LF'][iT,:,iEP['RootCoarse']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]

	# Calculate fine root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootFine'])
	vo['C_LF'][iT,:,iEP['RootFine']]=tr*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]

	# Adjust litterfall to account for N application response
	if meta['Modules']['NutrientApp']['iApplication'].size>0:
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Litterfall')

	#--------------------------------------------------------------------------
	# Update summary variables
	#--------------------------------------------------------------------------

	# Update gross growth
	vo['C_G_Gross'][iT,:,:]=vo['C_G_Net_Reg'][iT,:,:]+C_M_Reg

	# Update NPP
	# *** Compiled after the simulation is complete. ***
	#vo['C_NPP'][iT,:,:]=vo['C_G_Net_Reg'][iT,:,:]+C_M_Reg+vo['C_LF'][iT,:,:]

	return vo

#%% Dead organic matter and soil organic matter dynamics
def DeadWoodLitterAndSoilDynamics(meta,pNam,iT,iBat,vi,vo,iEP):

	# Extract parameters
	bIPF=meta['Param']['BEV']['InterPoolFluxes']
	bDec=meta['Param']['BEV']['Decomposition']

	#--------------------------------------------------------------------------
	# Flux of carbon between biomass pools and dead organic matter pools
	# Current DOM pools = previous DOM pools + biomass turnover
	#--------------------------------------------------------------------------

	# Transfer biomass turnover to very fast litter pool
	vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterVF']]+ \
		bIPF['FoliageLitToLitterVF']*vo['C_LF'][iT,:,iEP['Foliage']]+ \
		bIPF['RootFineLitToLitterVF']*vo['C_LF'][iT,:,iEP['RootFine']]+ \
		bIPF['FoliageMorToLitterVF']*vo['C_M_Reg'][iT,:,iEP['Foliage']]

	# Transfer biomass turnover to fast litter pool
	vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterF']]+ \
		bIPF['BranchLitToLitterF']*vo['C_LF'][iT,:,iEP['Branch']]+ \
		bIPF['BarkLitToLitterF']*vo['C_LF'][iT,:,iEP['Bark']]+ \
		bIPF['RootCoarseLitToLitterF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
		bIPF['BranchMorToLitterF']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
		bIPF['StemNonMerchLitToLitterF']*vo['C_LF'][iT,:,iEP['StemNonMerch']]

	# Transfer biomass turnover to medium litter pool
	vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterM']]+ \
		bIPF['BarkLitToLitterM']*vo['C_LF'][iT,:,iEP['Bark']]+ \
		bIPF['BarkMorToLitterM']*vo['C_M_Reg'][iT,:,iEP['Bark']]

	# Transfer biomass turnover to slow litter pool
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT-1,:,iEP['LitterS']]+ \
		bIPF['BarkLitToLitterS']*vo['C_LF'][iT,:,iEP['Bark']]+ \
		bIPF['BarkMorToLitterS']*vo['C_M_Reg'][iT,:,iEP['Bark']]

	# Transfer biomass turnover to snag stemwood pool
	vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagStem']]+ \
		bIPF['StemMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemMerch']]+ \
		bIPF['StemNonMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemNonMerch']]

	# Transfer mortality from live merch volume to dead merch volume
	vo['V_MerchDead'][iT,:]=vo['V_MerchDead'][iT,:]+vi['lsat']['Biomass to Volume CF']*bIPF['StemMerchMorToSnagStem']*vo['C_M_Reg'][iT,:,iEP['StemMerch']]

	# Transfer biomass turnover to snag branch pool
	vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['SnagBranch']]+ \
		bIPF['BranchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Branch']]+ \
		bIPF['BranchMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Branch']]+ \
		bIPF['BarkMorToSnagBranch']*vo['C_M_Reg'][iT,:,iEP['Bark']]+ \
		bIPF['StemNonMerchLitToSnagBranch']*vo['C_LF'][iT,:,iEP['Bark']]

	# Transfer biomass turnover to very fast soil pool
	vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilVF']]+ \
		bIPF['RootFineLitToSoilVF']*vo['C_LF'][iT,:,iEP['RootFine']] +\
		bIPF['RootCoarseMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]+ \
		bIPF['RootFineMorToSoilVF']*vo['C_M_Reg'][iT,:,iEP['RootFine']]

	# Transfer biomass turnover to fast soil pool
	vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilF']]+ \
		bIPF['RootCoarseLitToSoilF']*vo['C_LF'][iT,:,iEP['RootCoarse']]+ \
		bIPF['RootCoarseMorToSoilF']*vo['C_M_Reg'][iT,:,iEP['RootCoarse']]

	# Transfer biomass turnover to slow soil pool
	vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT-1,:,iEP['SoilS']]

	#--------------------------------------------------------------------------
	# Update piles
	# The above section updates DOM pools, but not piled pools
	#--------------------------------------------------------------------------

	vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledStemMerch']]
	vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledStemNonMerch']]
	vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledBranch']]
	vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledBark']]
	vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledSnagStem']]
	vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT-1,:,iEP['PiledSnagBranch']]

	#--------------------------------------------------------------------------
	# Decomposition
	# Notes:
	# 1) Vectorizing the parts of this that can be vectorized was attempted, but
	# runtime actually increased a bit.
	# 2) *** Currently no decomposition effect on dead volume ***
	#--------------------------------------------------------------------------

	# Prepare air temperature for respiration calculation
	Tref=10
	fT=(vi['lsat']['MAT']-Tref)/10

	# Respiration rate - Note that these terms do not equal the
	# atmosphere-bound efflux from heterotrophic respiration as a fraction is
	# emitted to the atmosphere and the remaining fraction is reorganized
	# within the ecosystem.
	meta[pNam]['Project']['R_LitterVF']=bDec['LitterVF_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']].flatten()*bDec['LitterVF_Q10']**fT
	meta[pNam]['Project']['R_LitterF']=bDec['LitterF_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterF']].flatten()*bDec['LitterF_Q10']**fT
	meta[pNam]['Project']['R_LitterM']=bDec['LitterM_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterM']].flatten()*bDec['LitterM_Q10']**fT
	meta[pNam]['Project']['R_LitterS']=bDec['LitterS_R10']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']].flatten()*bDec['LitterS_Q10']**fT
	meta[pNam]['Project']['R_SnagStem']=bDec['SnagStem_R10']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']].flatten()*bDec['SnagStem_Q10']**fT
	meta[pNam]['Project']['R_SnagBranch']=bDec['SnagBranch_R10']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']].flatten()*bDec['SnagBranch_Q10']**fT
	meta[pNam]['Project']['R_SoilVF']=bDec['SoilVF_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']].flatten()*bDec['SoilVF_Q10']**fT
	meta[pNam]['Project']['R_SoilF']=bDec['SoilF_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilF']].flatten()*bDec['SoilF_Q10']**fT
	meta[pNam]['Project']['R_SoilS']=bDec['SoilS_R10']*vo['C_Eco_Pools'][iT,:,iEP['SoilS']].flatten()*bDec['SoilS_Q10']**fT

	meta[pNam]['Project']['R_PiledStemMerch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']].flatten()*bDec['Piled_Q10']**fT
	meta[pNam]['Project']['R_PiledStemNonMerch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']].flatten()*bDec['Piled_Q10']**fT
	meta[pNam]['Project']['R_PiledBranch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']].flatten()*bDec['Piled_Q10']**fT
	meta[pNam]['Project']['R_PiledBark']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledBark']].flatten()*bDec['Piled_Q10']**fT
	meta[pNam]['Project']['R_PiledSnagStem']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']].flatten()*bDec['Piled_Q10']**fT
	meta[pNam]['Project']['R_PiledSnagBranch']=bDec['Piled_R10']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']].flatten()*bDec['Piled_Q10']**fT

	# Adjust decomposition to account for N application response
	if meta['Modules']['NutrientApp']['iApplication'].size>0:
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'HeterotrophicRespiration')

	# Remove respired carbon from source DOM pools
	vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]-meta[pNam]['Project']['R_LitterVF']
	vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]-meta[pNam]['Project']['R_LitterF']
	vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]-meta[pNam]['Project']['R_LitterM']
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-meta[pNam]['Project']['R_LitterS']
	vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-meta[pNam]['Project']['R_SnagStem']
	vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-meta[pNam]['Project']['R_SnagBranch']
	vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]-meta[pNam]['Project']['R_SoilVF']
	vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]-meta[pNam]['Project']['R_SoilF']
	vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]-meta[pNam]['Project']['R_SoilS']

	vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]-meta[pNam]['Project']['R_PiledStemMerch']
	vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]-meta[pNam]['Project']['R_PiledStemNonMerch']
	vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]-meta[pNam]['Project']['R_PiledBranch']
	vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]-meta[pNam]['Project']['R_PiledBark']
	vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]-meta[pNam]['Project']['R_PiledSnagStem']
	vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]-meta[pNam]['Project']['R_PiledSnagBranch']

	# Re-define decayed fast litter
	vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+bIPF['SnagBranchToLitterF']*meta[pNam]['Project']['R_SnagBranch']

	# Re-define decayed medium litter
	vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+bIPF['SnagStemToLitterM']*meta[pNam]['Project']['R_SnagStem']

	# Re-define decayed slow litter
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterVFToLitterS']*meta[pNam]['Project']['R_LitterVF']
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterFToLitterS']*meta[pNam]['Project']['R_LitterF']
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]+bIPF['LitterMToLitterS']*meta[pNam]['Project']['R_LitterM']

	# Re-define decayed slow soil
	vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilVFToSoilS']*meta[pNam]['Project']['R_SoilVF']
	vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+bIPF['SoilFToSoilS']*meta[pNam]['Project']['R_SoilF']

	# Heterotrophic respiration
	vo['C_RH'][iT,:,iEP['LitterVF']]=bIPF['LitterVFToCO2']*meta[pNam]['Project']['R_LitterVF']
	vo['C_RH'][iT,:,iEP['LitterF']]=bIPF['LitterFToCO2']*meta[pNam]['Project']['R_LitterF']
	vo['C_RH'][iT,:,iEP['LitterM']]=bIPF['LitterMToCO2']*meta[pNam]['Project']['R_LitterM']
	vo['C_RH'][iT,:,iEP['LitterS']]=bIPF['LitterSToCO2']*meta[pNam]['Project']['R_LitterS']
	vo['C_RH'][iT,:,iEP['SnagStem']]=bIPF['SnagStemToCO2']*meta[pNam]['Project']['R_SnagStem']
	vo['C_RH'][iT,:,iEP['SnagBranch']]=bIPF['SnagBranchToCO2']*meta[pNam]['Project']['R_SnagBranch']
	vo['C_RH'][iT,:,iEP['SoilVF']]=bIPF['SoilVFToCO2']*meta[pNam]['Project']['R_SoilVF']
	vo['C_RH'][iT,:,iEP['SoilF']]=bIPF['SoilFToCO2']*meta[pNam]['Project']['R_SoilF']
	vo['C_RH'][iT,:,iEP['SoilS']]=bIPF['SoilSToCO2']*meta[pNam]['Project']['R_SoilS']

	vo['C_RH'][iT,:,iEP['PiledStemMerch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledStemMerch']
	vo['C_RH'][iT,:,iEP['PiledStemNonMerch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledStemNonMerch']
	vo['C_RH'][iT,:,iEP['PiledBranch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledBranch']
	vo['C_RH'][iT,:,iEP['PiledBark']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledBark']
	vo['C_RH'][iT,:,iEP['PiledSnagStem']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledSnagStem']
	vo['C_RH'][iT,:,iEP['PiledSnagBranch']]=bIPF['PiledToCO2']*meta[pNam]['Project']['R_PiledSnagBranch']

	#--------------------------------------------------------------------------
	# Physical transfer
	#--------------------------------------------------------------------------

	PT_LitterS=bDec['LitterS_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
	PT_SnagStem=bDec['SnagStem_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
	PT_SnagBranch=bDec['SnagBranch_PhysTransRate']*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]

	# Remove carbon that is physically transferred
	vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-PT_LitterS
	vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-PT_SnagStem
	vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-PT_SnagBranch

	# Remove volume that is physically transferred from snags to litter
	vo['V_MerchDead'][iT,:]=vo['V_MerchDead'][iT,:]-vi['lsat']['Biomass to Volume CF']*PT_SnagStem

	# Add decomposed carbon to more decomposed DOM pools
	vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]+PT_LitterS
	vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+PT_SnagBranch
	vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+PT_SnagStem

	#--------------------------------------------------------------------------
	# litter decomposition
	# *** Added for fertilization study ***
	#--------------------------------------------------------------------------

	#vo['C_Eco_Pools'][iT,:,28]=meta[pNam]['Project']['R_LitterVF']+meta[pNam]['Project']['R_LitterF']+meta[pNam]['Project']['R_LitterM']+meta[pNam]['Project']['R_LitterS']

	return vo

#%% Disturbance and management events
def DisturbanceAndManagementEvents(meta,pNam,iT,iScn,iEns,iBat,vi,vo,iEP):

	# Update total (live+dead) stemwood merchantable volume
	vo['V_MerchTotal'][iT,:]=vo['V_MerchLive'][iT,:]+vo['V_MerchDead'][iT,:]

	# Predict wind (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Wind Status Historical']=='On') & (meta[pNam]['Scenario'][iScn]['Wind Status Future']=='On'):
		vi=asm.PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])
	if (meta[pNam]['Scenario'][iScn]['Wind Status Historical']=='On') & (meta[pNam]['Scenario'][iScn]['Wind Status Future']!='On'):
		if (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
			vi=asm.PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])
	if (meta[pNam]['Scenario'][iScn]['Wind Status Historical']!='On') & (meta[pNam]['Scenario'][iScn]['Wind Status Future']=='On'):
		if (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
			vi=asm.PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])
			
	# Predict disease (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Disease Status Historical']=='On') & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi=asm.PredictDisease_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])
	if (meta[pNam]['Scenario'][iScn]['Disease Status Future']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi=asm.PredictDisease_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])

	# Predict wildfire (on the fly)
	if meta[pNam]['Scenario'][iScn]['Wildfire Scn Pre-obs']!=-9999:
		if vi['tv'][iT]<1920:			
			vi=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns)
	if meta[pNam]['Scenario'][iScn]['Wildfire Scn Obs Period']!=-9999:
		if (vi['tv'][iT]>=1920) & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']): 
			vi=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns)
	if meta[pNam]['Scenario'][iScn]['Wildfire Scn Future']!=-9999:
		if vi['tv'][iT]>=meta[pNam]['Project']['Year Project']:
			vi=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns)

	# Predict beetles (on the fly)
	if meta[pNam]['Scenario'][iScn]['Beetle Scn Pre-obs']!=-9999:
		if vi['tv'][iT]<1950:
			vi=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns)
	if meta[pNam]['Scenario'][iScn]['Beetle Scn Obs Period']!=-9999:
		if (vi['tv'][iT]>=1950) & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']): 
			vi=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns)
	if meta[pNam]['Scenario'][iScn]['Beetle Scn Future']!=-9999:
		if vi['tv'][iT]>=meta[pNam]['Project']['Year Project']:
			vi=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns)

	# Predict harvesting (on the fly)
	if meta[pNam]['Scenario'][iScn]['Harvest Status Historical']=='On':
		if vi['tv'][iT]<meta[pNam]['Scenario'][iScn]['Harvest Year Transition']:
			Period='Historical'
			vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['V_MerchTotal'][iT,:],Period)
	if meta[pNam]['Scenario'][iScn]['Harvest Status Future']=='On':
		if vi['tv'][iT]>=meta[pNam]['Scenario'][iScn]['Harvest Year Transition']:
			Period='Future'
			vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['V_MerchTotal'][iT,:],Period)

	# Predict frost (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Frost Status Historical']=='On') & (vi['tv'][iT]<1950):
		vi=asm.PredictFrost_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])
	if (meta[pNam]['Scenario'][iScn]['Frost Status Future']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi=asm.PredictFrost_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:])

	# Predict future nutrient application (on the fly)
	if meta[pNam]['Scenario'][iScn]['Nutrient Application Status']=='On':
		if vi['tv'][iT]>=meta[pNam]['Project']['Year Project']:
			vi=napp.ScheduleNutrientApplication(meta,pNam,vi,vo,iT,iScn,iEns,iBat)
	
	# Predict future non-obligation stand establishment (on the fly)
	if meta[pNam]['Scenario'][iScn]['NOSE Status']=='On':
		if vi['tv'][iT]>=meta[pNam]['Project']['Year Project']:
			vi=tft.PredictNOSE_OnTheFly(meta,pNam,iScn,iBat,vi,iT)

	# Check to see how many events occur in this time step (don't do more than necessary)
	NumEventsInTimeStep=np.sum(np.sum(vi['EC']['ID Event Type'][iT,:,:]>0,axis=0)>0)

	# Initialize indicator of aerial nutrient application
	flag_nutrient_app=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])

	# Keep track of total tree biomass at the start of each annual time step
	C_Biomass_t0=np.sum(vo['C_Eco_Pools'][iT,:,iEP['BiomassTotal']],axis=0)

	# Loop through events in year
	for iE in range(NumEventsInTimeStep):

		# Event type IDs for the iE'th event of the year
		ID_Type=vi['EC']['ID Event Type'][iT,:,iE].copy()

		# Indexes to each unique event type
		idx_Type=gu.IndicesFromUniqueArrayValues(ID_Type)

		# Total affected biomass carbon
		MortalityFactor=vi['EC']['Mortality Factor'][iT,:,iE].copy()

		#----------------------------------------------------------------------
		# Get event-specific parameters
		#----------------------------------------------------------------------
		b={}
		for k1 in meta['Param']['BEV']['Events'][1].keys():
			b[k1]=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])
			for k2 in idx_Type.keys():
				if k2==0:
					continue
				b[k1][idx_Type[k2]]=meta['Param']['BEV']['Events'][k2][k1]
		
		#----------------------------------------------------------------------
		# Change in land cover / land use
		#----------------------------------------------------------------------
		for k in meta['Param']['Raw']['LUC']['LC Final'].keys():
			if meta['LUT']['Event'][k] in idx_Type:
				cd=meta['Param']['Raw']['LUC']['LC Final'][k]
				vo['LandCover'][iT:,idx_Type[meta['LUT']['Event'][k]]]=meta['LUT']['Derived']['lc_comp1'][cd]
				cd=meta['Param']['Raw']['LUC']['LU Final'][k]
				vo['LandUse'][iT:,idx_Type[meta['LUT']['Event'][k]]]=meta['LUT']['Derived']['lu_comp1'][cd]
			
		#----------------------------------------------------------------------
		# Record stands with aerial nutrient application
		#----------------------------------------------------------------------
		#iApp=np.where( (ID_Type==meta['LUT']['Event']['Nutrient App Aerial']) )[0]
		#flag_nutrient_app[iApp]=1
		if meta['LUT']['Event']['Nutrient App Aerial'] in idx_Type:
			flag_nutrient_app[idx_Type[meta['LUT']['Event']['Nutrient App Aerial']]]=1	  

		#----------------------------------------------------------------------
		# Adjust event-specific parameters to reflect time- and region-specific
		# fate of felled material
		#----------------------------------------------------------------------
		# Index to harvesting
		#iHarvest=np.where( (ID_Type==meta['LUT']['Event']['Harvest']) | (ID_Type==meta['LUT']['Event']['Harvest Salvage']) )[0]
		
		# Adjust fate of felled material parameters
		#if iHarvest.size>0:
		if meta['LUT']['Event']['Harvest'] in idx_Type:
			
			iHarvest=idx_Type[meta['LUT']['Event']['Harvest']]
			
			# Index to time-dependent fate of felled materials
			iT_P=np.where(meta['Param']['BE']['FelledFate']['Year']==meta[pNam]['Year'][iT])[0]

			# Simulations may exceed the timeframe of the felled fate parameters
			# If so, set to the last year
			if iT_P.size==0:
				iT_P=-1

			for k in meta['Param']['BEV']['FelledFate'].keys():
				b[k][iHarvest]=meta['Param']['BEV']['FelledFate'][k][iT_P,iHarvest]
			
			# Also update age at harvest
			vo['A_Harvest'][iT,iHarvest]=vo['A'][iT,iHarvest]

		#----------------------------------------------------------------------
		# Net-down insect mortality to reflect the proportion of host species
		#----------------------------------------------------------------------
		if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
			for nam in meta['Param']['Raw']['DisturbanceSpeciesAffected']['Insect Name']:
				id=meta['LUT']['Event'][nam]
				if id in idx_Type.keys():
					ndf=vi['lsat']['Insect Mortality Percent Tree Species Affected'][nam].astype('float')/100
					ind=idx_Type[id]
					MortalityFactor[ind]=MortalityFactor[ind]*ndf[ind]
		
		#----------------------------------------------------------------------
		# Define the amount of each pool that is affected by the event
		#----------------------------------------------------------------------
		# Affected biomass carbon
		if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
			Affected_StemwoodMerch=vo['C_M_Tot'][iT,:,iEP['StemMerch']]
			Affected_StemwoodNonMerch=vo['C_M_Tot'][iT,:,iEP['StemNonMerch']]
			Affected_Foliage=vo['C_M_Tot'][iT,:,iEP['Foliage']]
			Affected_Branch=vo['C_M_Tot'][iT,:,iEP['Branch']]
			Affected_Bark=vo['C_M_Tot'][iT,:,iEP['Bark']]
			Affected_RootCoarse=vo['C_M_Tot'][iT,:,iEP['RootCoarse']]
			Affected_RootFine=vo['C_M_Tot'][iT,:,iEP['RootFine']]
		else:
			Affected_StemwoodMerch=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]
			Affected_StemwoodNonMerch=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]
			Affected_Foliage=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['Foliage']]
			Affected_Branch=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['Branch']]
			Affected_Bark=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['Bark']]
			Affected_RootCoarse=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]
			Affected_RootFine=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['RootFine']]

		# Partition bark into merch and non-merch components
		#Affected_BarkMerch=0.85*Affected_Bark
		#Affected_BarkNonMerch=(1.00-0.85)*Affected_Bark

		# Sum up total affected non-merchantable biomass
		#Affected_TotNonMerch=Affected_StemwoodNonMerch+Affected_Branch+Affected_BarkNonMerch
		Affected_TotNonMerch=Affected_StemwoodNonMerch+Affected_Branch

		# Snags
		Affected_SnagStem=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]
		Affected_SnagBranch=MortalityFactor*vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]

		# All biomass
		Affected_BiomassTotal=Affected_StemwoodMerch+Affected_StemwoodNonMerch+Affected_Foliage+Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine

		# All biomass and dead wood
		if meta['LUT']['Event']['Harvest'] in idx_Type:
			Affected_BiomassAndSnags=Affected_BiomassTotal+Affected_SnagStem+Affected_SnagBranch
			vo['C_Felled'][iT,iHarvest]=Affected_BiomassAndSnags[iHarvest]
			vo['C_FelledMerch'][iT,iHarvest]=Affected_StemwoodMerch[iHarvest]+Affected_Bark[iHarvest]+Affected_SnagStem[iHarvest]
			vo['C_FelledRoots'][iT,iHarvest]=Affected_RootCoarse[iHarvest]+Affected_Bark[iHarvest]+Affected_RootFine[iHarvest]

		# Live merch. stemwood volume
		Affected_VolumeStemMerchLive=MortalityFactor*vo['V_MerchLive'][iT,:]

		# Dead merch. stemwood volume
		Affected_VolumeStemMerchDead=MortalityFactor*vo['V_MerchDead'][iT,:]

		#----------------------------------------------------------------------
		# Calculate mortality
		#----------------------------------------------------------------------
		# Total mortality
		vo['C_M_Dist'][iT,:]=vo['C_M_Dist'][iT,:]+Affected_BiomassTotal

		# Mortality by category (lumping stands together)
		for k in idx_Type.keys():
			if k==0:
				continue
			ind=idx_Type[k]
			vo['C_M_ByAgent'][k][iT,ind]=vo['C_M_ByAgent'][k][iT,ind]+Affected_BiomassTotal[ind]

		#----------------------------------------------------------------------
		# Remove affected amount from each pool
		#----------------------------------------------------------------------
		if meta[pNam]['Project']['Biomass Module']!='Sawtooth':

			# Remove carbon from affected biomass pools
			vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemMerch']]-Affected_StemwoodMerch
			vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['StemNonMerch']]-Affected_StemwoodNonMerch
			vo['C_Eco_Pools'][iT,:,iEP['Foliage']]=vo['C_Eco_Pools'][iT,:,iEP['Foliage']]-Affected_Foliage
			vo['C_Eco_Pools'][iT,:,iEP['Branch']]=vo['C_Eco_Pools'][iT,:,iEP['Branch']]-Affected_Branch
			vo['C_Eco_Pools'][iT,:,iEP['Bark']]=vo['C_Eco_Pools'][iT,:,iEP['Bark']]-Affected_Bark
			vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]-Affected_RootCoarse
			vo['C_Eco_Pools'][iT,:,iEP['RootFine']]=vo['C_Eco_Pools'][iT,:,iEP['RootFine']]-Affected_RootFine

			# Remove carbon from snag pools
			vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]-Affected_SnagStem
			vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]-Affected_SnagBranch

			# Remove stemwood merch volume
			vo['V_MerchLive'][iT,:]=np.maximum(0,vo['V_MerchLive'][iT,:]-Affected_VolumeStemMerchLive)
			vo['V_MerchDead'][iT,:]=np.maximum(0,vo['V_MerchDead'][iT,:]-Affected_VolumeStemMerchDead)

		#----------------------------------------------------------------------
		# We have not explicity tracked a partition of snags into merch and
		# non merch, but these categories may be removed in different degrees
		# Move a proportion of snag stemwood to non-merch biomass
		#----------------------------------------------------------------------
		#Transfer_SnagStemToNonMerch=meta['Param']['BEV']['InterPoolFluxes']['FractionSnagStemThatIsNonMerch']*Affected_SnagStem
		#Affected_SnagStem=Affected_SnagStem-Transfer_SnagStemToNonMerch
		#Affected_TotNonMerch=Affected_TotNonMerch+Transfer_SnagStemToNonMerch

		#----------------------------------------------------------------------
		# Carbon that is removed (ie sent to mills)
		#----------------------------------------------------------------------

		# Merch biomass to mill - of the total amount of biomass affected,
		vo['C_ToMillMerch'][iT,:]=vo['C_ToMillMerch'][iT,:]+b['StemwoodMerchRemoved']*Affected_StemwoodMerch
		vo['C_ToMillMerch'][iT,:]=vo['C_ToMillMerch'][iT,:]+b['BarkRemoved']*Affected_Bark

		# NonMerch biomass to mill - of the total amount of biomass affected,
		# what fraction of non-merch biomass was sent to the mill?
		# - NonMerch = NonMerchStem + Foliage + Branch + Bark
		vo['C_ToMillNonMerch'][iT,:]=vo['C_ToMillNonMerch'][iT,:]+b['StemwoodNonMerchRemoved']*Affected_TotNonMerch

		# Snag stemwood to mill
		vo['C_ToMillSnagStem'][iT,:]=vo['C_ToMillSnagStem'][iT,:]+b['SnagStemRemoved']*Affected_SnagStem

		# Snag branches to mill
		vo['C_ToMillNonMerch'][iT,:]=vo['C_ToMillNonMerch'][iT,:]+b['SnagBranchRemoved']*Affected_SnagBranch

		#----------------------------------------------------------------------
		# Volume that is removed (sent to mills)
		#----------------------------------------------------------------------

		if meta[pNam]['Project']['Biomass Module']!='Sawtooth':
			vo['V_ToMillMerchLive'][iT,:]=vo['V_ToMillMerchLive'][iT,:]+b['StemwoodMerchRemoved']*Affected_VolumeStemMerchLive
			vo['V_ToMillMerchDead'][iT,:]=vo['V_ToMillMerchDead'][iT,:]+b['SnagStemRemoved']*Affected_VolumeStemMerchDead
		else:
			vo['V_ToMillMerchLive'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_ToMillMerch'][iT,:]
			vo['V_ToMillMerchDead'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_ToMillSnagStem'][iT,:]

		# Total merch stemwood volume removed
		vo['V_ToMillMerchTotal'][iT,:]=vo['V_ToMillMerchLive'][iT,:]+vo['V_ToMillMerchDead'][iT,:]

		# Non-mech volume removed
		vo['V_ToMillNonMerch'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_ToMillNonMerch'][iT,:]

		#----------------------------------------------------------------------
		# Carbon that is left dispersed on site (after felling or wind storms)
		#----------------------------------------------------------------------

		# Foliage transferred directly to very fast litter
		c=b['FoliageLeftOnSite']*Affected_Foliage
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+c

		# Stem, branch and bark carbon transfered to medium and fast litter pools
		c=b['BranchLeftOnSite']*Affected_Branch
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+c

		c=b['BarkLeftOnSite']*Affected_Bark
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+c

		c=b['SnagBranchLeftOnSite']*Affected_SnagBranch
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+c

		c=b['StemwoodMerchLeftOnSite']*Affected_StemwoodMerch
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+c

		c=b['StemwoodNonMerchLeftOnSite']*Affected_StemwoodNonMerch
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+b['StemwoodNonMerchLeftOnSite']*Affected_StemwoodNonMerch
		vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+c

		c=b['SnagStemLeftOnSite']*Affected_SnagStem
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]+c

		# Roots transferred directly to DOM
		c=0.5*Affected_RootFine
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]+c

		c=0.5*Affected_RootFine
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]+c

		c=0.5*Affected_RootCoarse
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]+c

		c=0.5*Affected_RootCoarse
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]+c

		#----------------------------------------------------------------------
		# Carbon that is piled
		#----------------------------------------------------------------------

		PiledStemMerch=(1.0-meta['Param']['BE']['ProductDisposal']['PiledStemwoodToFirewoodDom'])*b['StemwoodMerchPiled']*Affected_StemwoodMerch
		PiledStemNonMerch=b['StemwoodNonMerchPiled']*Affected_StemwoodNonMerch
		PiledBranch=b['BranchPiled']*Affected_Branch
		PiledBark=b['BarkPiled']*Affected_Bark
		PiledFoliage=b['FoliagePiled']*Affected_Foliage
		PiledSnagStem=(1.0-meta['Param']['BE']['ProductDisposal']['PiledStemwoodToFirewoodDom'])*b['SnagStemPiled']*Affected_SnagStem
		PiledSnagBranch=b['SnagBranchPiled']*Affected_SnagBranch

		# Add to piles
		vo['C_ToPile'][iT,:]=PiledStemMerch+PiledStemNonMerch+PiledBranch+PiledBark+PiledFoliage+PiledSnagStem+PiledSnagBranch
		vo['C_ToPileMerch'][iT,:]=PiledStemMerch+PiledSnagStem
		vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]+PiledStemMerch
		vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]+PiledStemNonMerch
		vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]+PiledBranch
		vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]+PiledBark+PiledFoliage
		#vo['C_Eco_Pools'][iT,:,iEP['PiledFoliage']]=vo['C_Eco_Pools'][iT,:,iEP['PiledFoliage']]+b['FoliagePiled']*Affected_Foliage
		vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]+PiledSnagStem
		vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]+PiledSnagBranch

		# A small fraction of piled wood is collected for firewood
		vo['C_ToFirewoodDom'][iT,:]=vo['C_ToFirewoodDom'][iT,:]+meta['Param']['BE']['ProductDisposal']['PiledStemwoodToFirewoodDom']*(PiledStemMerch+PiledSnagStem)

		# Piles that are burned
		StemMerch=b['PiledStemMerchBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]
		StemNonMerch=b['PiledStemNonMerchBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]
		Branch=b['PiledBranchBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]
		Bark=b['PiledBarkBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]
		SnagStem=b['PiledSnagStemBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]
		SnagBranch=b['PiledSnagBranchBurned']*vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]
		Total_Burned=StemMerch+StemNonMerch+Branch+Bark+SnagStem+SnagBranch
		NonMerch_Burned=StemNonMerch+Branch+Bark+SnagBranch

		# Remove affected carbon
		vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemMerch']]-StemMerch
		vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledStemNonMerch']]-StemNonMerch
		vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBranch']]-Branch
		vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]=vo['C_Eco_Pools'][iT,:,iEP['PiledBark']]-Bark
		vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagStem']]-SnagStem
		vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['PiledSnagBranch']]-SnagBranch

		# Add to fire emissions
		vo['C_E_OpenBurningAsCO2'][iT,:]=vo['C_E_OpenBurningAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
		vo['C_E_OpenBurningAsCH4'][iT,:]=vo['C_E_OpenBurningAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
		vo['C_E_OpenBurningAsCO'][iT,:]=vo['C_E_OpenBurningAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned

		# If it is slashpile burning, record it
		vo['C_ToSlashpileBurnTot'][iT,:]=vo['C_ToSlashpileBurnTot'][iT,:]+Total_Burned
		vo['C_ToSlashpileBurnNonMerch'][iT,:]=vo['C_ToSlashpileBurnNonMerch'][iT,:]+NonMerch_Burned

		#----------------------------------------------------------------------
		# Carbon that is moved from biomass to snags
		#----------------------------------------------------------------------

		# Stem biomass to snag stem
		vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]+b['StemwoodMerchToSnagStem']*Affected_StemwoodMerch
		vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]+b['StemwoodNonMerchToSnagStem']*Affected_StemwoodNonMerch
		vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]=vo['C_Eco_Pools'][iT,:,iEP['SnagStem']]+b['BarkToSnagStem']*Affected_Bark

		# Branch biomass to snag branch
		vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]=vo['C_Eco_Pools'][iT,:,iEP['SnagBranch']]+b['BranchToSnagBranch']*Affected_Bark

		#----------------------------------------------------------------------
		# Biomass and DOM burned in wildfire
		#----------------------------------------------------------------------

		BurnedStemwoodMerch=b['StemwoodMerchBurned']*Affected_StemwoodMerch
		BurnedStemwoodNonMerch=b['StemwoodNonMerchBurned']*Affected_StemwoodNonMerch
		BurnedBark=b['BarkBurned']*Affected_Bark
		BurnedBranch=b['BranchBurned']*Affected_Branch
		BurnedFoliage=b['FoliageBurned']*Affected_Foliage
		BurnedSnagStem=b['SnagStemBurned']*Affected_SnagStem
		BurnedSnagBranch=b['SnagBranchBurned']*Affected_SnagBranch

		BurnedLitterVF=MortalityFactor*b['LitterVFBurned']*vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]
		BurnedLitterF=MortalityFactor*b['LitterFBurned']*vo['C_Eco_Pools'][iT,:,iEP['LitterF']]
		BurnedLitterM=MortalityFactor*b['LitterMBurned']*vo['C_Eco_Pools'][iT,:,iEP['LitterM']]
		BurnedLitterS=MortalityFactor*b['LitterSBurned']*vo['C_Eco_Pools'][iT,:,iEP['LitterS']]
		BurnedSoilVF=MortalityFactor*b['SoilVFBurned']*vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]
		BurnedSoilF=MortalityFactor*b['SoilFBurned']*vo['C_Eco_Pools'][iT,:,iEP['SoilF']]
		BurnedSoilS=MortalityFactor*b['SoilSBurned']*vo['C_Eco_Pools'][iT,:,iEP['SoilS']]

		vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterVF']]-BurnedLitterVF
		vo['C_Eco_Pools'][iT,:,iEP['LitterF']]=vo['C_Eco_Pools'][iT,:,iEP['LitterF']]-BurnedLitterF
		vo['C_Eco_Pools'][iT,:,iEP['LitterM']]=vo['C_Eco_Pools'][iT,:,iEP['LitterM']]-BurnedLitterM
		vo['C_Eco_Pools'][iT,:,iEP['LitterS']]=vo['C_Eco_Pools'][iT,:,iEP['LitterS']]-BurnedLitterS
		vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilVF']]-BurnedSoilVF
		vo['C_Eco_Pools'][iT,:,iEP['SoilF']]=vo['C_Eco_Pools'][iT,:,iEP['SoilF']]-BurnedSoilF
		vo['C_Eco_Pools'][iT,:,iEP['SoilS']]=vo['C_Eco_Pools'][iT,:,iEP['SoilS']]-BurnedSoilS

		BurnedBiomass=BurnedStemwoodMerch+BurnedStemwoodNonMerch+BurnedBark+BurnedBranch+BurnedFoliage
		BurnedDeadWood=BurnedSnagStem+BurnedSnagBranch
		BurnedLitter=BurnedLitterVF+BurnedLitterF+BurnedLitterM+BurnedLitterS
		BurnedSoil=BurnedSoilVF+BurnedSoilF+BurnedSoilS
		BurnedTotal=BurnedBiomass+BurnedDeadWood+BurnedLitter+BurnedSoil

		vo['C_E_WildfireAsCO2'][iT,:]=vo['C_E_WildfireAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*BurnedTotal
		vo['C_E_WildfireAsCH4'][iT,:]=vo['C_E_WildfireAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*BurnedTotal
		vo['C_E_WildfireAsCO'][iT,:]=vo['C_E_WildfireAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*BurnedTotal

		#----------------------------------------------------------------------
		# Update stand age
		#----------------------------------------------------------------------

		if (meta[pNam]['Project']['Biomass Module']!='Sawtooth') & (meta[pNam]['Project']['Partial Mortality Affects Age']=='On'):

			# List of exception event types (that will remain at the same age)
			Exceptions_to_Partial_Mortality=[  ]

			# Index to types where age will change
			ind=np.where( (np.isin(ID_Type,Exceptions_to_Partial_Mortality)==False) )[0]

			# Assume oldest trees were most affected, reduce age in prportion
			# with mortality rate
			# Not always realistic, but see how net growth is affected.
			vo['A'][iT,ind]=vo['A'][iT,ind]*(1-MortalityFactor)

		# Ensure stand-replacing events reset stand age to 0
		vo['A'][iT,(MortalityFactor==1)]=0

		# Ensure planting resets to age 0
		#vo['A'][iT,(ID_Type==meta['LUT']['Event']['Planting']) | (ID_Type==meta['LUT']['Event']['Direct Seeding'])]=0
		if meta['LUT']['Event']['Planting'] in idx_Type.keys():
			vo['A'][iT,idx_Type[meta['LUT']['Event']['Planting']]]=0

		#----------------------------------------------------------------------
		# Transition to new growth curve
		#----------------------------------------------------------------------
		# Only applies to BatchTIPSY
		if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':
			for iGC in range(meta['Modules']['GYM']['N Growth Curves']):

				# Don't alter growth curve for fertilization - index to non-fertilization events
				ind=np.where( (vi['EC']['ID Growth Curve'][iT,:,iE]==meta['Modules']['GYM']['ID GC Unique'][iGC]) & (ID_Type!=meta['LUT']['Event']['Nutrient App Aerial']) )[0]

				if ind.size>0:
					# Notes: This crashes often when the GY model has been accidently run with the wrong project info
					# - because the project path in BatchTipsy is wrong.
					vi['GC']['Active'][:,ind,:]=vi['GC'][ meta['Modules']['GYM']['ID GC Unique'][iGC] ][:,ind,:].astype(float)*meta['Modules']['GYM']['Scale Factor']
					vi['GC']['ID_GCA'][ind]=int(meta['Modules']['GYM']['ID GC Unique'][iGC])

		#----------------------------------------------------------------------
		# Impose regen failure
		#----------------------------------------------------------------------
		if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':
			#iFailure=np.where(ID_Type==meta['LUT']['Event']['Regen Failure'])[0]
			if meta['LUT']['Event']['Regen Failure'] in idx_Type.keys():
				iFailure=idx_Type[meta['LUT']['Event']['Regen Failure']]
				vi['GC']['Active'][:,iFailure,:]=0
				vi['GC'][1][:,iFailure,:]=0
			
			if meta[pNam]['Scenario'][iScn]['NOSE Status']=='On':
				#iPoorGrowth=np.where(ID_Type==meta['LUT']['Event']['Regen at 25% Growth'])[0]
				if meta['LUT']['Event']['Regen at 25% Growth'] in idx_Type.keys():
					iPoorGrowth=idx_Type[meta['LUT']['Event']['Regen at 25% Growth']]
					vi['GC']['Active'][:,iPoorGrowth,:]=0.25*vi['GC']['Active'][:,iPoorGrowth,:]
					vi['GC'][1][:,iPoorGrowth,:]=0.25*vi['GC'][1][:,iPoorGrowth,:]

		#------------------------------------------------------------------
		# Update net growth (in response to lethal events)
		#
		# Partial disturbances likely have a lasting impact on gross growth
		# by removing growing stock. Without ad hoc
		# adjustment, this will not be reflected in the growth curve from TIPSY.
		#
		# This will apply to instances where Mortality=100% only if the
		# user does not specify a change in growth curve.
		# With this functionality, the user has the options of:
		#  1) Changing the growth curve
		#  2) Assuming a gradual recovery to the previous curve, as defined
		#	 by the default recovery half life for that disturbance type.
		#------------------------------------------------------------------

		# *** This is not up to date - will crash if attempted - needs update. ***
		flg=0

		if (flg==1) & (meta[pNam]['Project']['Biomass Module']!='Sawtooth'):

			# Half life
			hl=meta['Param']['BEV']['bDist_GrowthRecovery_HL']

			# Only proceed if:
			#   1) the disturbance has a lasting growth impact/recovery
			#   2) the user is not specifying a change in growth curve
			if (hl>0) & (flg_gc_change==0):

				# Extract net growth for active growth curve
				NetGrowth=vi['GCA'][:,iS,:].copy().astype(float)*meta['Modules']['GYM']['Scale Factor']

				# Growth pre-event
				G_pre=NetGrowth.copy()

				# Growth post-event
				G_post=(1-Biomass_Affected_Frac)*NetGrowth.copy()

				# Difference in growth
				dG=G_pre-G_post

				# Age vector
				A=np.arange(1,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

				# Time since disturbance
				TSD=np.tile(A-vo['A'][iT,iS],(G_pre.shape[1],1)).T

				# Relative effect
				fTSD=1/(1+np.exp(-0.5*(TSD-hl)))
				fTSD[0:int(vo['A'][iT,iS]),:]=0

				# Growth recovery
				G_Recovery=fTSD*dG

				# New estimate of net growth
				NetGrowthNew=G_post+G_Recovery

				# Add back to dictionary
				NetGrowthNew=NetGrowthNew*meta['Modules']['GYM']['Scale Factor']
				NetGrowthNew=NetGrowthNew.astype(np.int16)

				vi['GCA'][:,iS,:]=NetGrowthNew

		#----------------------------------------------------------------------
		# Apply growth factors (in response to non-lethal events)
		# Growth factor should comes in DMEC as percent effect
		# 10 = 10% increase, -10 = 10% decrease
		#----------------------------------------------------------------------

		flg=1
		if (flg==1) & (meta[pNam]['Project']['Biomass Module']!='Sawtooth'):

			GrowthFactor=vi['EC']['Growth Factor'][iT,:,iE].copy()

			indAdj=np.where( (GrowthFactor!=9999) & (GrowthFactor!=0) )[0]

			if indAdj.size>0:

				Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

				# Convert growth factor to respoonse ratio
				GrowthResponseRatio=GrowthFactor.astype(float)/100+1.0

				NetGrowth=vi['GC']['Active'].copy()

				for iAdj in range(indAdj.size):

					if ID_Type[iAdj]==meta['LUT']['Event']['IDW']:

						# Only a temporary change in growth (see severity class table)
						# *** I have not checked to see if this is working properly ***
						ResponsePeriod=6

						iResponse=np.where( (Age>=vo['A'][iT,iAdj]) & (Age<=vo['A'][iT,iAdj]+ResponsePeriod) )[0]

						NetGrowth[iResponse,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[iResponse,iAdj,:]

					else:

						# A permanent change in growth
						NetGrowth[:,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[:,iAdj,:]

				vi['GC']['Active']=NetGrowth

	#--------------------------------------------------------------------------
	# Relative gravimetric mortality (%)
	#--------------------------------------------------------------------------

	for k in meta['LUT']['Event'].keys():
		id=meta['LUT']['Event'][k]
		vo['C_M_Pct_ByAgent'][id][iT,:]=vo['C_M_ByAgent'][id][iT,:]/C_Biomass_t0*100

	#--------------------------------------------------------------------------
	# Aerial nutrient application events
	#--------------------------------------------------------------------------

	# Generate index to stands that were fertilized
	meta['Modules']['NutrientApp']['iApplication']=np.where(flag_nutrient_app==1)[0]

	if meta['Modules']['NutrientApp']['iApplication'].size>0:

		# Adjust net growth of aboveground biomass
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'AbovegroundNetGrowth')

		# Adjust emissions
		vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Emissions')

	return vo,vi

#%%
def ProductDynamics(meta,pNam,iT,iBat,vi,vo):

	#--------------------------------------------------------------------------
	# Index to time-dependent parameters
	#--------------------------------------------------------------------------

	# This will be the same between end use and removed fate because the time vectors
	# are the same.
	# *** Don't change one time vector without changing the other. ***
	iT_P=np.where(meta['Param']['BE']['ProductTypes']['Year']==meta[pNam]['Year'][iT])[0]

	# Removed fate parameters
	bRF={}
	for k in meta['Param']['BEV']['RemovedFate'].keys():
		if iT_P.size!=0:
			bRF[k]=meta['Param']['BEV']['RemovedFate'][k][iT_P,:].flatten()
		else:
			# Simulations may extend beyond the 2100 limit of the time-dependent HWP variables
			bRF[k]=meta['Param']['BEV']['RemovedFate'][k][-1,:].flatten()

	# Product type parameters
	bPT={}
	for k in meta['Param']['BEV']['ProductTypes'].keys():
		if iT_P.size!=0:
			bPT[k]=meta['Param']['BEV']['ProductTypes'][k][iT_P,:].flatten()
		else:
			# Simulations may extend beyond the 2100 limit of the time-dependent HWP variables
			bPT[k]=meta['Param']['BEV']['ProductTypes'][k][-1,:].flatten()

	# Product disposal parameters
	bPD=meta['Param']['BEV']['ProductDisposal']

	#--------------------------------------------------------------------------
	# Carbon transferred from forest to mills or direct to end-uses
	#--------------------------------------------------------------------------

	C_ToPulpMill=bRF['RemovedMerchToPulpMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToPulpMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToPulpMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToPelletMill=bRF['RemovedMerchToPelletMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToPelletMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToPelletMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToLumberMill=bRF['RemovedMerchToLumberMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToLumberMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToLumberMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToPlywoodMill=bRF['RemovedMerchToPlywoodMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToPlywoodMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToPlywoodMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToOSBMill=bRF['RemovedMerchToOSBMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToOSBMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToOSBMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToMDFMill=bRF['RemovedMerchToMDFMill']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToMDFMill']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToMDFMill']*vo['C_ToMillSnagStem'][iT,:]

	C_ToLogExport=bRF['RemovedMerchToLogExport']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToLogExport']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToLogExport']*vo['C_ToMillSnagStem'][iT,:]

	C_ToPowerGrid=bRF['RemovedMerchToIPP']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToIPP']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToIPP']*vo['C_ToMillSnagStem'][iT,:]

	C_ToFirewoodCollection=bRF['RemovedMerchToFirewood']*vo['C_ToMillMerch'][iT,:] + \
				 bRF['RemovedNonMerchToFirewood']*vo['C_ToMillNonMerch'][iT,:] + \
				 bRF['RemovedSnagStemToFirewood']*vo['C_ToMillSnagStem'][iT,:]

	#C_ChipperMill=meta['Param']['BEV']['RemovedFate']['RemovedMerchToChipperMill']*vo['C_ToMillMerch'][iT,:] + \
	#			 meta['Param']['BEV']['RemovedFate']['RemovedNonMerchToChipperMill']*vo['C_ToMillNonMerch'][iT,:] + \
	#			 meta['Param']['BEV']['RemovedFate']['RemovedSnagStemToChipperMill']*vo['C_ToMillSnagStem'][iT,:]

	#C_PolePostMill=bRF['RemovedMerchToPolePostMill']*vo['C_ToMillMerch'][iT,:] + \
	#			 bRF['RemovedNonMerchToPolePostMill']*vo['C_ToMillNonMerch'][iT,:] + \
	#			 bRF['RemovedSnagStemToPolePostMill']*vo['C_ToMillSnagStem'][iT,:]

	#C_ShakeShingleMill=bRF['RemovedMerchToShakeShingleMill']*vo['C_ToMillMerch'][iT,:] + \
	#			 bRF['RemovedNonMerchToShakeShingleMill']*vo['C_ToMillNonMerch'][iT,:] + \
	#			 bRF['RemovedSnagStemToShakeShingleMill']*vo['C_ToMillSnagStem'][iT,:]

	#--------------------------------------------------------------------------
	# Carbon transferred from mill to mill
	#--------------------------------------------------------------------------

	C_ToPulpMill=C_ToPulpMill+C_ToLumberMill*bPT['LumberMillToPulpMill']
	C_ToMDFMill=C_ToMDFMill+C_ToLumberMill*bPT['LumberMillToMDFMill']
	C_ToPelletMill=C_ToPelletMill+C_ToLumberMill*bPT['LumberMillToPelletMill']

	#--------------------------------------------------------------------------
	# Carbon transferred from Lumber mill logs to log Export
	#--------------------------------------------------------------------------

	C_ToLogExport=C_ToLogExport+bPT['LumberMillToLogExport']*C_ToLumberMill

	#--------------------------------------------------------------------------
	# Carbon transferred from log Export to firewood
	#--------------------------------------------------------------------------

	C_ToFirewoodExport=bPT['LogExportToFirewood']*C_ToLogExport

	#--------------------------------------------------------------------------
	# Carbon transferred to paper
	#--------------------------------------------------------------------------

	C_ToPaper=bPT['PulpMillToPaper']*C_ToPulpMill

	#--------------------------------------------------------------------------
	# Carbon transferred to pulp effluent
	#--------------------------------------------------------------------------

	C_ToPulpEffluent=bPT['PulpMillToEffluent']*C_ToPulpMill

	#--------------------------------------------------------------------------
	# Carbon transferred to pellets
	#--------------------------------------------------------------------------

	C_ToPelletExport=bPT['PelletMillToPelletExport']*C_ToPelletMill
	C_ToPelletDomGrid=bPT['PelletMillToDomGrid']*C_ToPelletMill
	C_ToPelletDomRNG=bPT['PelletMillToDomRNG']*C_ToPelletMill

	#--------------------------------------------------------------------------
	# Carbon transferred to facility power
	#--------------------------------------------------------------------------

	C_ToPowerFacilityDom=bPT['LumberMillToPowerFacility']*C_ToLumberMill + \
		bPT['PulpMillToPowerFacility']*C_ToPulpMill + \
		bPT['PlywoodMillToPowerFacility']*C_ToPlywoodMill + \
		bPT['OSBMillToPowerFacility']*C_ToOSBMill
	C_ToPowerFacilityExport=bPT['LogExportToPowerFacility']*C_ToLogExport

	#--------------------------------------------------------------------------
	# Carbon transferred to grid by independent power producers
	#--------------------------------------------------------------------------

	C_ToPowerGrid=C_ToPowerGrid+bPT['LumberMillToIPP']*C_ToLumberMill + \
		bPT['PulpMillToIPP']*C_ToPulpMill + \
		bPT['PlywoodMillToIPP']*C_ToPlywoodMill + \
		bPT['OSBMillToIPP']*C_ToOSBMill

	#--------------------------------------------------------------------------
	# Log-size effect
	#--------------------------------------------------------------------------

	#	bLogSizeEffect=0.02
	#
	#	iLSE=np.where(vo['LogSizeEnhancement'][iT,:]>0)[0]
	#
	#	if iLSE.size>0:
	#
	#		C_Transfer=bLogSizeEffect*C_Paper[iLSE]
	#		C_Paper[iLSE]=C_Paper[iLSE]-C_Transfer
	#		C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
	#
	#		C_Transfer=bLogSizeEffect*C_Pellet[iLSE]
	#		C_Pellet[iLSE]=C_Pellet[iLSE]-C_Transfer
	#		C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
	#
	#		C_Transfer=bLogSizeEffect*C_PowerFacilityDom[iLSE]
	#		C_PowerFacilityDom[iLSE]=C_PowerFacilityDom[iLSE]-C_Transfer
	#		C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
	#
	#		C_Transfer=bLogSizeEffect*C_PowerFacilityFor[iLSE]
	#		C_PowerFacilityFor[iLSE]=C_PowerFacilityFor[iLSE]-C_Transfer
	#		C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer
	#
	#		C_Transfer=bLogSizeEffect*C_PowerGrid[iLSE]
	#		C_PowerGrid[iLSE]=C_PowerGrid[iLSE]-C_Transfer
	#		C_LumberMill[iLSE]=C_LumberMill[iLSE]+C_Transfer

	#--------------------------------------------------------------------------
	# Production of single-family homes
	#--------------------------------------------------------------------------

	C_LumberMillToSFH=C_ToLumberMill*bPT['LumberMillToSFH']
	C_PlywoodMillToSFH=C_ToPlywoodMill*bPT['PlywoodMillToSFH']
	C_OSBMillToSFH=C_ToOSBMill*bPT['OSBMillToSFH']
	C_MDFMillToSFH=C_ToMDFMill*bPT['MDFMillToSFH']
	C_LogExportToSFH=C_ToLogExport*bPT['LogExportToSFH']
	#C_ShakeShingleToSFH=C_ShakeShingleMill*bPD['ShakeShingleMillToSFH']

	#--------------------------------------------------------------------------
	# Production of multi-family homes
	#--------------------------------------------------------------------------

	C_LumberMillToMFH=C_ToLumberMill*bPT['LumberMillToMFH']
	C_PlywoodMillToMFH=C_ToPlywoodMill*bPT['PlywoodMillToMFH']
	C_OSBMillToMFH=C_ToOSBMill*bPT['OSBMillToMFH']
	C_MDFMillToMFH=C_ToMDFMill*bPT['MDFMillToMFH']
	C_LogExportToMFH=C_ToLogExport*bPT['LogExportToMFH']
	#C_ShakeShingleToMFH=C_ShakeShingleMill*bPT['ShakeShingleMillToMFH']

	#--------------------------------------------------------------------------
	# Production of commercial buildings
	#--------------------------------------------------------------------------

	C_LumberMillToCom=C_ToLumberMill*bPT['LumberMillToCom']
	C_PlywoodMillToCom=C_ToPlywoodMill*bPT['PlywoodMillToCom']
	C_OSBMillToCom=C_ToOSBMill*bPT['OSBMillToCom']
	C_MDFMillToCom=C_ToMDFMill*bPT['MDFMillToCom']
	C_LogExportToCom=C_ToLogExport*bPT['LogExportToCom']
	#C_ShakeShingleToCom=C_ShakeShingleMill*bPT['ShakeShingleMillToCom']

	#--------------------------------------------------------------------------
	# Production of furniture
	#--------------------------------------------------------------------------

	C_LumberMillToFurn=C_ToLumberMill*bPT['LumberMillToFurn']
	C_PlywoodMillToFurn=C_ToPlywoodMill*bPT['PlywoodMillToFurn']
	C_OSBMillToFurn=C_ToOSBMill*bPT['OSBMillToFurn']
	C_MDFMillToFurn=C_ToMDFMill*bPT['MDFMillToFurn']
	C_LogExportToFurn=C_ToLogExport*bPT['LogExportToFurn']

	#--------------------------------------------------------------------------
	# Production of shipping containers
	#--------------------------------------------------------------------------

	C_LumberMillToShip=C_ToLumberMill*bPT['LumberMillToShip']
	C_PlywoodMillToShip=C_ToPlywoodMill*bPT['PlywoodMillToShip']
	C_OSBMillToShip=C_ToOSBMill*bPT['OSBMillToShip']
	C_MDFMillToShip=C_ToMDFMill*bPT['MDFMillToShip']
	C_LogExportToShip=C_ToLogExport*bPT['LogExportToShip']

	#--------------------------------------------------------------------------
	# Production of repairs
	#--------------------------------------------------------------------------

	C_LumberMillToRepairs=C_ToLumberMill*bPT['LumberMillToRepairs']
	C_PlywoodMillToRepairs=C_ToPlywoodMill*bPT['PlywoodMillToRepairs']
	C_OSBMillToRepairs=C_ToOSBMill*bPT['OSBMillToRepairs']
	C_MDFMillToRepairs=C_ToMDFMill*bPT['MDFMillToRepairs']
	C_LogExportToRepairs=C_ToLogExport*bPT['LogExportToRepairs']

	#--------------------------------------------------------------------------
	# Production of other
	#--------------------------------------------------------------------------

	C_LumberMillToOther=C_ToLumberMill*bPT['LumberMillToOther']
	C_PlywoodMillToOther=C_ToPlywoodMill*bPT['PlywoodMillToOther']
	C_OSBMillToOther=C_ToOSBMill*bPT['OSBMillToOther']
	C_MDFMillToOther=C_ToMDFMill*bPT['MDFMillToOther']
	C_LogExportToOther=C_ToLogExport*bPT['LogExportToOther']
	#C_PolePostMillToOther=PolePostMill*bPT['PolePostMillToOther']

	#--------------------------------------------------------------------------
	# Track sales (for economic modelling)
	# *** Conservation of mass is not being retained in the below set of variables
	# - log Export retain the full amount of log Export, but some of it was
	# also added to foreign power facility and foreign firewood.
	# Economics only consider log Export, substitution effects only consider
	# foreign power facility ad foreign firewood. ***
	#--------------------------------------------------------------------------

	vo['C_ToLumber'][iT,:]=C_LumberMillToSFH+C_LumberMillToMFH+C_LumberMillToCom+C_LumberMillToFurn+C_LumberMillToShip+C_LumberMillToRepairs+C_LumberMillToOther

	vo['C_ToPlywood'][iT,:]=C_PlywoodMillToSFH+C_PlywoodMillToMFH+C_PlywoodMillToCom+C_PlywoodMillToFurn+C_PlywoodMillToShip+C_PlywoodMillToRepairs+C_PlywoodMillToOther

	vo['C_ToOSB'][iT,:]=C_OSBMillToSFH+C_OSBMillToMFH+C_OSBMillToCom+C_OSBMillToFurn+C_OSBMillToShip+C_OSBMillToRepairs+C_OSBMillToOther

	vo['C_ToMDF'][iT,:]=C_MDFMillToSFH+C_MDFMillToMFH+C_MDFMillToCom+C_MDFMillToFurn+C_MDFMillToShip+C_MDFMillToRepairs+C_MDFMillToOther

	vo['C_ToPaper'][iT,:]=C_ToPaper

	vo['C_ToPowerFacilityDom'][iT,:]=C_ToPowerFacilityDom

	vo['C_ToPowerFacilityExport'][iT,:]=C_ToPowerFacilityExport

	vo['C_ToPowerGrid'][iT,:]=C_ToPowerGrid

	vo['C_ToPelletExport'][iT,:]=C_ToPelletExport
	vo['C_ToPelletDomGrid'][iT,:]=C_ToPelletDomGrid
	vo['C_ToPelletDomRNG'][iT,:]=C_ToPelletDomRNG

	# *** There is firewood taken directly from forest ecosystems (see events function) ***
	vo['C_ToFirewoodDom'][iT,:]=vo['C_ToFirewoodDom'][iT,:]+C_ToFirewoodCollection

	vo['C_ToFirewoodExport'][iT,:]=C_ToFirewoodExport

	vo['C_ToLogExport'][iT,:]=C_ToLogExport

	#--------------------------------------------------------------------------
	# Carbon transferred from mills to in-use products
	#--------------------------------------------------------------------------

	# Transfer mill fibre to single-family homes
	iP=meta['Core']['iPP']['SFH']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToSFH + \
		C_PlywoodMillToSFH + \
		C_OSBMillToSFH + \
		C_MDFMillToSFH + \
		C_LogExportToSFH

	# Transfer mill fibre to multi-family homes
	iP=meta['Core']['iPP']['MFH']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToMFH + \
		C_PlywoodMillToMFH + \
		C_OSBMillToMFH + \
		C_MDFMillToMFH + \
		C_LogExportToMFH

	# Transfer mill fibre to commercial
	iP=meta['Core']['iPP']['Comm']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToCom + \
		C_PlywoodMillToCom + \
		C_OSBMillToCom + \
		C_MDFMillToCom + \
		C_LogExportToCom

	# Transfer mill fibre to furniture
	iP=meta['Core']['iPP']['Furn']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToFurn + \
		C_PlywoodMillToFurn + \
		C_OSBMillToFurn + \
		C_MDFMillToFurn + \
		C_LogExportToFurn

	# Transfer mill fibre to shipping
	iP=meta['Core']['iPP']['Ship']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToShip + \
		C_PlywoodMillToShip + \
		C_OSBMillToShip + \
		C_MDFMillToShip + \
		C_LogExportToShip

	# Transfer mill fibre to repairs
	iP=meta['Core']['iPP']['Repairs']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToRepairs + \
		C_PlywoodMillToRepairs + \
		C_OSBMillToRepairs + \
		C_MDFMillToRepairs + \
		C_LogExportToRepairs

	# Transfer mill fibre to other
	iP=meta['Core']['iPP']['Other']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP] + \
		C_LumberMillToOther + \
		C_PlywoodMillToOther + \
		C_OSBMillToOther + \
		C_MDFMillToOther + \
		C_LogExportToOther

	# Transfer pulp mill fibre to paper
	iP=meta['Core']['iPP']['Paper']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPaper

	# Transfer mill fibre to power facitity, domestic
	iP=meta['Core']['iPP']['PowerFacilityDom']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPowerFacilityDom

	# Transfer mill fibre to power facitity, foreign
	iP=meta['Core']['iPP']['PowerFacilityExport']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPowerFacilityExport

	# Transfer mill fibre to pellet Export
	iP=meta['Core']['iPP']['PelletExport']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletExport

	# Transfer mill fibre to pellet domestic grid
	iP=meta['Core']['iPP']['PelletDomGrid']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletDomGrid

	# Transfer mill fibre to pellet domestic RNG
	iP=meta['Core']['iPP']['PelletDomRNG']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPelletDomRNG

	# Transfer domestic firewood to domestic firewood pool
	iP=meta['Core']['iPP']['FirewoodDom']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToFirewoodCollection

	# Transfer foreign firewood to foreign firewood pool
	iP=meta['Core']['iPP']['FirewoodExport']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToFirewoodExport

	# Transfer pulp mill carbon to pulp-mill effluent
	iP=meta['Core']['iPP']['EffluentPulp']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT-1,:,iP]+C_ToPulpEffluent

	#--------------------------------------------------------------------------
	# Update dump and landfill reservoirs
	#--------------------------------------------------------------------------

	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpWood']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpWood']]
	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['DumpPaper']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['DumpPaper']]
	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodDegradable']]
	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillWoodNonDegradable']]
	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperDegradable']]
	vo['C_Pro_Pools'][iT,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]=vo['C_Pro_Pools'][iT-1,:,meta['Core']['iPP']['LandfillPaperNonDegradable']]

	#--------------------------------------------------------------------------
	# Single-family homes --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['SFH']
	C_retired=bPD['SFH_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['SFHToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['SFHToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['SFHToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Multi-family homes --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['MFH']
	C_retired=bPD['MFH_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['MFHToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['MFHToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['MFHToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Commercial building --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Comm']
	C_retired=bPD['Comm_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['CommToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['CommToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['CommToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Furniture --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Furn']
	C_retired=bPD['Furn_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['FurnToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['FurnToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['FurnToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Shipping --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Ship']
	C_retired=bPD['Ship_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['ShipToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['ShipToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['ShipToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Repairs --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Repairs']
	C_retired=bPD['Repairs_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['RepairsToDumpWood']*C_retired

	# Transfer carbon to landfill (degradble)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['RepairsToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['RepairsToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Other --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Other']
	C_retired=bPD['Other_tr']*vo['C_Pro_Pools'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['OtherToDumpWood']*C_retired

	# Transfer carbon to landfill (degradble)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['OtherToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['OtherToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Paper --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover (with adjustment for recycling)
	iP=meta['Core']['iPP']['Paper']
	C_retired=(1-bPD['PaperRecycleRate'])*bPD['Paper_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_retired

	# Transfer to dump
	iP=meta['Core']['iPP']['DumpPaper']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['PaperToDumpPaper']*C_retired

	# Transfer to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['PaperToLandfillPaper']*bPD['ToLandfillPaperDegradableFrac']*C_retired

	# Transfer to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] + bPD['PaperToLandfillPaper']*(1-bPD['ToLandfillPaperDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Emissions from combustion during domestic power generation
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['PowerFacilityDom']
	C_emitted=bPD['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyPowerFacilityDom'][iT,:]=vo['E_CO2e_ESC_BioenergyPowerFacilityDom'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion during foreign power generation
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['PowerFacilityExport']
	C_emitted=bPD['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyPowerFacilityExport'][iT,:]=vo['E_CO2e_ESC_BioenergyPowerFacilityExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet Export
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['PelletExport']
	C_emitted=bPD['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyPelletExport'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet domestic grid
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['PelletDomGrid']
	C_emitted=bPD['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyPelletDomGrid'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletDomGrid'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet domestic RNG
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['PelletDomRNG']
	C_emitted=bPD['Energy_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyPelletDomRNG'][iT,:]=vo['E_CO2e_ESC_BioenergyPelletDomRNG'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of domestic firewood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['FirewoodDom']
	C_emitted=bPD['Firewood_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyFirewoodDom'][iT,:]=vo['E_CO2e_ESC_BioenergyFirewoodDom'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of foreign firewood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['FirewoodExport']
	C_emitted=bPD['Firewood_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP] - C_emitted

	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*C_emitted
	E_CO2e_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_CO2e_ESC_Bioenergy'][iT,:]=vo['E_CO2e_ESC_Bioenergy'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)
	vo['E_CO2e_ESC_BioenergyFirewoodExport'][iT,:]=vo['E_CO2e_ESC_BioenergyFirewoodExport'][iT,:] + (E_CO2e_AsCO2+E_CO2e_AsCH4)

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from pulp effluent
	#--------------------------------------------------------------------------

	# Emissions from pulp effluent (CO2 from aerobic decomposition)
	iP=meta['Core']['iPP']['EffluentPulp']
	C_emitted=bPD['EffluentPulp_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Remove emitted carbon
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

	# Add emitted carbon to CO2 emission "pool"
	# Emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	#vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of dump wood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['DumpWood']
	C_emitted=bPD['DumpWood_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Removal
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*C_emitted

	# Add to emissions (CO2 emission from aerobic decomposition)
	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	#vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of dump paper
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['DumpPaper']
	C_emitted=bPD['DumpPaper_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Removal
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*C_emitted

	# Add to emissions (CO2 emission from aerobic decomposition)
	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:] + E_CO2e_AsCO2

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	#vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of landfill degradable wood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	C_emitted=bPD['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Removal
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-C_emitted

	# Add to emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['LandfillDegradableFracEmitCO2']*C_emitted

	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCO2

	# Adjustment for proportion of degradable landfills with gas collection systems,
	# efficiency of system, and methane oxided to CO2 from the landfill cover
	gcsp=bPD['LandfillMethaneEmit_GasColSysProp']
	c1=1-gcsp
	c2=1-bPD['LandfillMethaneEmit_GasColSysEffic']
	E_C_AsCH4=C_emitted*(c1-bPD['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-bPD['LandfillMethaneOxidizedToCO2']*c1)

	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*E_C_AsCH4

	E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_CH4

	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCH4

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of landfill degradable paper
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['LandfillPaperDegradable']
	c_emitted=bPD['LandfillWoodDegradable_tr']*vo['C_Pro_Pools'][iT,:,iP]

	# Removal
	vo['C_Pro_Pools'][iT,:,iP]=vo['C_Pro_Pools'][iT,:,iP]-c_emitted

	# Add to emissions
	E_CO2e_AsCO2=meta['Param']['BEV']['Biophysical']['Ratio_CO2_to_C']*bPD['LandfillDegradableFracEmitCO2']*C_emitted

	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCO2

	# Adjustment for proportion of degradable landfills with gas collection systems,
	# efficiency of system, and methane oxided to CO2 from the landfill cover
	E_C_AsCH4=C_emitted*(c1-bPD['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-bPD['LandfillMethaneOxidizedToCO2']*c1)
	E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_C_AsCH4

	E_CH4=meta['Param']['BEV']['Biophysical']['Ratio_CH4_to_C']*E_C_AsCH4

	E_CO2e_AsCH4=meta['Param']['BEV']['Biophysical']['GWP_CH4_AR5']*E_CH4

	vo['E_CO2e_LULUCF_HWP'][iT,:]=vo['E_CO2e_LULUCF_HWP'][iT,:]+E_CO2e_AsCH4

	vo['Atm_CO2_In'][iT,:]=vo['Atm_CO2_In'][iT,:]+E_CO2e_AsCO2
	vo['Atm_CH4_In'][iT,:]=vo['Atm_CH4_In'][iT,:]+E_CH4

	return vo

#%%
def GeologicalDynamics(meta,pNam,vi,vo):

	#==========================================================================
	# Operational emissions
	#==========================================================================

	# Extract parameters
	bB=meta['Param']['BEV']['Biophysical']
	bS=meta['Param']['BEV']['Substitution']

	# Electricity conversion efficiency ratio
	Ratio_EC=bB['Electrical Conversion Efficiency of Pellet Electricity Plant (>25MW)']/bB['Electrical Conversion Efficiency of Coal Electricity Plant']

	# Total removals (m3/ha)
	V_Tot=vo['V_ToMillMerchTotal']+vo['V_ToMillNonMerch']

	# Total removals (ODT/ha)
	ODT_Tot=V_Tot*bB['Density Wood']

	#--------------------------------------------------------------------------
	# Resource extraction
	#--------------------------------------------------------------------------

	# Construction and maintenance of roads
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Road Construction']*V_Tot

	# Cruising and reconnaissance
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Cruise And Recon']*V_Tot

	# Felling and processing logs
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Felling Process Logs']*V_Tot

	# Skidding trees to landing
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Skidding To Landing']*V_Tot

	# Piling and sorting
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Piling And Sorting Logs']*V_Tot

	# Loading logs
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Loading Logs At Landing']*V_Tot

	# Chipping (non-merch)
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Chipping']*vo['V_ToMillNonMerch']

	# Hauling (forest ecosystem to mill)
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+ODT_Tot/bB['Moisture Content Wood']*bB['Emission Intensity Transport Truck']*bB['Distance Forest To Mill (One Way)']

	# Site preparation
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Site Prep']*V_Tot

	# Sowing seeds
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Sowing']*V_Tot

	# Planting
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Planting']*V_Tot

	# Surveying
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Surveying']*V_Tot

	#--------------------------------------------------------------------------
	# Mill operations
	#--------------------------------------------------------------------------

	# Unloading
	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['EI Unloading At Mill']*V_Tot

	# Sawing and processing lumber
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing Lumber']*(vo['C_ToLumber']/bB['Carbon Content Wood'])

	# Sawing and processing plywood
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing Plywood']*(vo['C_ToPlywood']/bB['Carbon Content Wood'])

	# Sawing and processing OSB
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing OSB']*(vo['C_ToOSB']/bB['Carbon Content Wood'])

	# Sawing and processing MDF
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Sawing Processing MDF']*(vo['C_ToMDF']/bB['Carbon Content Wood'])

	# Sawing and processing pulp
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Processing Pulp']*(vo['C_ToPaper']/bB['Carbon Content Wood'])

	# Pellets
	C_Pellet=vo['C_ToPelletExport']+vo['C_ToPelletDomGrid']+vo['C_ToPelletDomRNG']

	# Size reduction of pellets
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Size Reduction']*C_Pellet/bB['Carbon Content Wood']

	# Drying pellets
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Drying']*C_Pellet/bB['Carbon Content Wood']

	# Pelletizing
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Pelletizing']*C_Pellet/bB['Carbon Content Wood']

	# Pellet sieving
	vo['E_CO2e_ESC_OperForBurnOil']=vo['E_CO2e_ESC_OperForBurnOil']+bB['EI Pellet Seiving']*C_Pellet/bB['Carbon Content Wood']

	#--------------------------------------------------------------------------
	# Transport Mill -> Distribution Hub (Vancouver)
	# *** Exclude raw log exports, as they originate at mills on the ocean ***
	#--------------------------------------------------------------------------

	Mass=(vo['C_ToLumber']+vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF']+vo['C_ToPaper']+vo['C_ToPelletExport']+vo['C_ToPelletDomRNG'])/bB['Carbon Content Wood']/bB['Moisture Content Lumber']

	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['Distance Mill To Distribution Hub']*bB['Emission Intensity Transport Rail']*Mass

	#--------------------------------------------------------------------------
	# Transport Hub -> Market
	#--------------------------------------------------------------------------

	# Solid wood products

	Mass=(vo['C_ToLumber']+vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF']+vo['C_ToPaper'])/bB['Carbon Content Wood']/bB['Moisture Content Lumber']

	E_HubToMarket_Water=bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction Solid Wood Product Water Dest 1']*bB['Distance Solid Wood Product Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

	E_HubToMarket_Rail=bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction Solid Wood Product Rail Dest 1']*bB['Distance Solid Wood Product Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass

	E_HubToMarket_Truck=bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
		bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
		bB['Fraction Solid Wood Product Truck Dest 1']*bB['Distance Solid Wood Product Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass

	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail+E_HubToMarket_Truck

	# Log exports

	Mass=(vo['C_ToLogExport'])/bB['Carbon Content Wood']/bB['Moisture Content Wood']

	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+bB['Distance LogExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

	# Pellets

	Mass=vo['C_ToPelletExport']/bB['Carbon Content Wood']

	E_HubToMarket_Water=bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

	E_HubToMarket_Rail=bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass

	# E_HubToMarket_Truck=bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
	#	 bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
	#	 bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass

	vo['E_CO2e_ET_OperForBurnOil']=vo['E_CO2e_ET_OperForBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail

	#==========================================================================
	# Substitution effects
	#==========================================================================

	#----------------------------------------------------------------------
	# Domestic facility power generation (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	Yield_PowerFacilityDom=vo['C_ToPowerFacilityDom']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerFacilityDom=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityDom

	GJ_PowerFacilityDom=GJ_PowerFacilityDom*Ratio_EC

	E_Sub_CoalForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityDom
	E_Sub_DieselForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityDom
	E_Sub_GasForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityDom
	E_Sub_OilForBioenergy_PowerFacilityDom=bS['PowerFacilityDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityDom

	#----------------------------------------------------------------------
	# Foreign facility power generation (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	# Yield (green tonnes/ha)
	Yield_PowerFacilityExport=vo['C_ToPowerFacilityExport']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerFacilityExport=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityExport

	GJ_PowerFacilityExport=GJ_PowerFacilityExport*Ratio_EC

	E_Sub_CoalForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityExport
	E_Sub_DieselForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityExport
	E_Sub_GasForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityExport
	E_Sub_OilForBioenergy_PowerFacilityExport=bS['PowerFacilityExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityExport

	#----------------------------------------------------------------------
	# Independent power producers (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	# Yield (green tonnes/ha)
	Yield_PowerGrid=vo['C_ToPowerGrid']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerGrid=bB['Energy Content Wood (0% moisture)']*Yield_PowerGrid

	GJ_PowerGrid=GJ_PowerGrid*Ratio_EC

	E_Sub_CoalForBioenergy_PowerGrid=bS['PowerGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerGrid
	E_Sub_DieselForBioenergy_PowerGrid=bS['PowerGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerGrid
	E_Sub_GasForBioenergy_PowerGrid=bS['PowerGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerGrid
	E_Sub_OilForBioenergy_PowerGrid=bS['PowerGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerGrid

	#----------------------------------------------------------------------
	# Pellet exports (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dried tonnes/ha)
	Yield_PelletExport=vo['C_ToPelletExport']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletExport=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletExport

	GJ_PelletExport=GJ_PelletExport*Ratio_EC

	E_Sub_CoalForBioenergy_PelletExport=bS['PelletExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletExport
	E_Sub_DieselForBioenergy_PelletExport=bS['PelletExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletExport
	E_Sub_GasForBioenergy_PelletExport=bS['PelletExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletExport
	E_Sub_OilForBioenergy_PelletExport=bS['PelletExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletExport

	#----------------------------------------------------------------------
	# Pellet domestic grid (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dired tonnes/ha)
	Yield_PelletDomGrid=vo['C_ToPelletDomGrid']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletDomGrid=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomGrid

	GJ_PelletDomGrid=GJ_PelletDomGrid*Ratio_EC

	E_Sub_CoalForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomGrid
	E_Sub_DieselForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomGrid
	E_Sub_GasForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomGrid
	E_Sub_OilForBioenergy_PelletDomGrid=bS['PelletDomGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomGrid

	#----------------------------------------------------------------------
	# Pellet domestic RNG (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dired tonnes/ha)
	Yield_PelletDomRNG=vo['C_ToPelletDomRNG']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletDomRNG=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomRNG

	E_Sub_CoalForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomRNG
	E_Sub_DieselForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomRNG
	E_Sub_GasForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomRNG
	E_Sub_OilForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomRNG

	#----------------------------------------------------------------------
	# Domestic firewood (MgC/ha) to (air dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (air dried tonnes/ha)
	Yield_FirewoodDom=vo['C_ToFirewoodDom']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_FirewoodDom=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodDom

	E_Sub_CoalForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodDom
	E_Sub_DieselForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodDom
	E_Sub_GasForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodDom
	E_Sub_OilForBioenergy_FirewoodDom=bS['FirewoodDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodDom

	#----------------------------------------------------------------------
	# Foreign firewood (MgC/ha) to (air dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (air dried tonnes/ha)
	Yield_FirewoodExport=vo['C_ToFirewoodExport']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_FirewoodExport=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodExport

	E_Sub_CoalForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodExport
	E_Sub_DieselForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodExport
	E_Sub_GasForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodExport
	E_Sub_OilForBioenergy_FirewoodExport=bS['FirewoodExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodExport

	#--------------------------------------------------------------------------
	# Substitution of fossil fuels for bioenergy
	# *** Save as positive and then change sign in post-processing ***
	#--------------------------------------------------------------------------

	vo['E_CO2e_SUB_CoalForBioenergy']= \
		(E_Sub_CoalForBioenergy_PowerFacilityDom + \
		 E_Sub_CoalForBioenergy_PowerFacilityExport + \
		 E_Sub_CoalForBioenergy_PowerGrid + \
		 E_Sub_CoalForBioenergy_PelletExport + \
		 E_Sub_CoalForBioenergy_PelletDomGrid + \
		 E_Sub_CoalForBioenergy_PelletDomRNG + \
		 E_Sub_CoalForBioenergy_FirewoodDom + \
		 E_Sub_CoalForBioenergy_FirewoodExport)

	vo['E_CO2e_SUB_OilForBioenergy']= \
		(E_Sub_OilForBioenergy_PowerFacilityDom + \
		 E_Sub_OilForBioenergy_PowerFacilityExport + \
		 E_Sub_OilForBioenergy_PowerGrid + \
		 E_Sub_OilForBioenergy_PelletExport + \
		 E_Sub_OilForBioenergy_FirewoodDom + \
		 E_Sub_OilForBioenergy_FirewoodExport + \
		 E_Sub_DieselForBioenergy_PowerFacilityDom + \
		 E_Sub_DieselForBioenergy_PowerFacilityExport + \
		 E_Sub_DieselForBioenergy_PowerGrid + \
		 E_Sub_DieselForBioenergy_PelletExport + \
		 E_Sub_DieselForBioenergy_PelletDomGrid + \
		 E_Sub_DieselForBioenergy_PelletDomRNG + \
		 E_Sub_DieselForBioenergy_FirewoodDom + \
		 E_Sub_DieselForBioenergy_FirewoodExport)

	vo['E_CO2e_SUB_GasForBioenergy']= \
		(E_Sub_GasForBioenergy_PowerFacilityDom + \
		 E_Sub_GasForBioenergy_PowerFacilityExport + \
		 E_Sub_GasForBioenergy_PowerGrid + \
		 E_Sub_GasForBioenergy_PelletExport + \
		 E_Sub_GasForBioenergy_PelletDomGrid + \
		 E_Sub_GasForBioenergy_PelletDomRNG + \
		 E_Sub_GasForBioenergy_FirewoodDom + \
		 E_Sub_GasForBioenergy_FirewoodExport)

	vo['E_CO2e_SUB_PowerFacilityDom']= \
		(E_Sub_CoalForBioenergy_PowerFacilityDom + \
		 E_Sub_OilForBioenergy_PowerFacilityDom + \
		 E_Sub_GasForBioenergy_PowerFacilityDom)

	vo['E_CO2e_SUB_PowerFacilityExport']= \
		(E_Sub_CoalForBioenergy_PowerFacilityExport + \
		 E_Sub_OilForBioenergy_PowerFacilityExport + \
		 E_Sub_GasForBioenergy_PowerFacilityExport)

	vo['E_CO2e_SUB_PowerGrid']= \
		(E_Sub_CoalForBioenergy_PowerGrid + \
		 E_Sub_OilForBioenergy_PowerGrid + \
		 E_Sub_GasForBioenergy_PowerGrid)

	vo['E_CO2e_SUB_PelletExport']= \
		(E_Sub_CoalForBioenergy_PelletExport + \
		 E_Sub_OilForBioenergy_PelletExport + \
		 E_Sub_GasForBioenergy_PelletExport)

	vo['E_CO2e_SUB_PelletDomGrid']= \
		(E_Sub_CoalForBioenergy_PelletDomGrid + \
		 E_Sub_OilForBioenergy_PelletDomGrid + \
		 E_Sub_GasForBioenergy_PelletDomGrid)

	vo['E_CO2e_SUB_PelletDomRNG']= \
		(E_Sub_CoalForBioenergy_PelletDomRNG + \
		 E_Sub_OilForBioenergy_PelletDomRNG + \
		 E_Sub_GasForBioenergy_PelletDomRNG)

	vo['E_CO2e_SUB_FirewoodDom']= \
		(E_Sub_CoalForBioenergy_FirewoodDom + \
		 E_Sub_OilForBioenergy_FirewoodDom + \
		 E_Sub_GasForBioenergy_FirewoodDom)

	vo['E_CO2e_SUB_FirewoodExport']= \
		(E_Sub_CoalForBioenergy_FirewoodExport + \
		 E_Sub_OilForBioenergy_FirewoodExport + \
		 E_Sub_GasForBioenergy_FirewoodExport)

	#----------------------------------------------------------------------
	# Substitution for structural wood produts (tCO2e/ha)
	#----------------------------------------------------------------------

	# Sawnwood (t DM)
	ODT_Sawnwood=(1/bB['Carbon Content Wood'])*vo['C_ToLumber']

	# Panel (t DM)
	ODT_Panel=(1/bB['Carbon Content Wood'])*(vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF'])

	# Residuals (tonnes)
	#Residuals=0

	# Substitution of concrete for structural wood
	fS_Concrete=bS['SawnwoodFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForSawnwood']*ODT_Sawnwood
	fP_Concrete=bS['PanelFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingConcrete']*bS['DisplacementRatio_ConcreteForResiduals']*Residuals
	vo['ODT Concrete']=fS_Concrete+fP_Concrete#-fR

	# Substitution of steel for structural wood
	fS_Steel=bS['SawnwoodFracDisplacingSteel']*bS['DisplacementRatio_SteelForSawnwood']*ODT_Sawnwood
	fP_Steel=bS['PanelFracDisplacingSteel']*bS['DisplacementRatio_SteelForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingSteel']*bS['DisplacementRatio_SteelForResiduals']*Residuals
	vo['ODT Steel']=fS_Steel+fP_Steel#-fR

	# Substitution of aluminum for structural wood
	fS_Aluminum=bS['SawnwoodFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForSawnwood']*ODT_Sawnwood
	fP_Aluminum=bS['PanelFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForResiduals']*Residuals
	vo['ODT Aluminum']=fS_Aluminum+fP_Aluminum#-fR

	# Substitution of plastics for structural wood
	fS_Plastic=bS['SawnwoodFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForSawnwood']*ODT_Sawnwood
	fP_Plastic=bS['PanelFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForResiduals']*Residuals
	vo['ODT Plastic']=fS_Plastic+fP_Plastic#-fR

	# Substitution of textiles for structural wood
	fS_Textile=bS['SawnwoodFracDisplacingTextile']*bS['DisplacementRatio_TextileForSawnwood']*ODT_Sawnwood
	fP_Textile=bS['PanelFracDisplacingTextile']*bS['DisplacementRatio_TextileForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingTextile']*bS['DisplacementRatio_TextileForResiduals']*Residuals
	vo['ODT Textile']=fS_Textile+fP_Textile#-fR

	# Emissions from productoin of structural materials
	vo['E_CO2e_SUB_Concrete']=bB['Emission Intensity Concrete']*vo['ODT Concrete']
	vo['E_CO2e_SUB_Steel']=bB['Emission Intensity Steel']*vo['ODT Steel']
	vo['E_CO2e_SUB_Aluminum']=bB['Emission Intensity Aluminum']*vo['ODT Aluminum']
	vo['E_CO2e_SUB_Plastic']=bB['Emission Intensity Plastic']*vo['ODT Plastic']
	vo['E_CO2e_SUB_Textile']=bB['Emission Intensity Textile']*vo['ODT Textile']

	# Emissions from sawnwood and Panel
	vo['E_CO2e_SUB_Sawnwood']=bB['Emission Intensity Concrete']*fS_Concrete+ \
		bB['Emission Intensity Steel']*fS_Steel + \
		bB['Emission Intensity Aluminum']*fS_Aluminum + \
		bB['Emission Intensity Plastic']*fS_Plastic + \
		bB['Emission Intensity Textile']*fS_Textile

	vo['E_CO2e_SUB_Panel']=bB['Emission Intensity Concrete']*fP_Concrete+ \
		bB['Emission Intensity Steel']*fP_Steel + \
		bB['Emission Intensity Aluminum']*fP_Aluminum + \
		bB['Emission Intensity Plastic']*fP_Plastic + \
		bB['Emission Intensity Textile']*fP_Textile

	# Emissions from structural matierials, tallied by feedstock (and calcination)
	vo['E_CO2e_SUB_CoalForWood']= \
		bS['FracConcreteEmissionsFromCoal']*vo['E_CO2e_SUB_Concrete']+ \
		bS['FracSteelEmissionsFromCoal']*vo['E_CO2e_SUB_Steel']+ \
		bS['FracAluminumEmissionsFromCoal']*vo['E_CO2e_SUB_Aluminum']+ \
		bS['FracPlasticEmissionsFromCoal']*vo['E_CO2e_SUB_Plastic']+ \
		bS['FracTextileEmissionsFromCoal']*vo['E_CO2e_SUB_Textile']

	vo['E_CO2e_SUB_OilForWood']= \
		bS['FracConcreteEmissionsFromOil']*vo['E_CO2e_SUB_Concrete']+ \
		bS['FracSteelEmissionsFromOil']*vo['E_CO2e_SUB_Steel']+ \
		bS['FracAluminumEmissionsFromOil']*vo['E_CO2e_SUB_Aluminum']+ \
		bS['FracPlasticEmissionsFromOil']*vo['E_CO2e_SUB_Plastic']+ \
		bS['FracTextileEmissionsFromOil']*vo['E_CO2e_SUB_Textile']

	vo['E_CO2e_SUB_GasForWood']= \
		bS['FracConcreteEmissionsFromGas']*vo['E_CO2e_SUB_Concrete']+ \
		bS['FracSteelEmissionsFromGas']*vo['E_CO2e_SUB_Steel']+ \
		bS['FracAluminumEmissionsFromGas']*vo['E_CO2e_SUB_Aluminum']+ \
		bS['FracPlasticEmissionsFromGas']*vo['E_CO2e_SUB_Plastic']+ \
		bS['FracTextileEmissionsFromGas']*vo['E_CO2e_SUB_Textile']

	vo['E_CO2e_SUB_ConcreteFromCalcination']=bS['FracConcreteEmissionsFromCalcination']*vo['E_CO2e_SUB_Concrete']
	vo['E_CO2e_SUB_ConcreteFromNonCalcination']=(1-bS['FracConcreteEmissionsFromCalcination'])*vo['E_CO2e_SUB_Concrete']

	#--------------------------------------------------------------------------
	# Back-calculate production of fossil fuel consumption from operational use
	# and substitution effects (GJ)
	# *** Moved to post-processing ***
	#--------------------------------------------------------------------------

	return vo

#%%
def GrassBiomassDynamics(meta,pNam,iScn,iBat,iT,vi,vo,iEP):
	# Aboveground biomass (foliage)
	vo['C_G_Gross'][iT,:,iEP['Foliage']]=vo['C_G_Gross'][iT,:,iEP['Foliage']]+1.5
	vo['C_LF'][iT,:,iEP['Foliage']]=vo['C_LF'][iT,:,iEP['Foliage']]+1.5

	# Belowground biomass
	if vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]>15:
		vo['C_G_Gross'][iT,:,iEP['RootCoarse']]=+1.5
	else:
		vo['C_G_Gross'][iT,:,iEP['RootCoarse']]=+1.75

	vo['C_LF'][iT,:,iEP['RootCoarse']]=+1.5

	# Net growth
	G_net=vo['C_G_Gross'][iT,:,iEP['RootCoarse']]-vo['C_LF'][iT,:,iEP['RootCoarse']]

	#vo['C_G_Net_Reg'][iT,:,iEP['RootCoarse']]=vo['C_G_Net_Reg'][iT,:,iEP['RootCoarse']]+G_net

	vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]=vo['C_Eco_Pools'][iT,:,iEP['RootCoarse']]+G_net

	return vo

