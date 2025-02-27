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

	# Extract incides for GY model data structure
	iGY=meta['Modules']['GYM']['GC Input Indices']

	# Extract net growth
	# Notes: I tried to make this faster using unravel but it was way slower
	NetGrowth=np.zeros( (meta[pNam]['Project']['Batch Size'][iBat],6) )
	for iS in range(meta[pNam]['Project']['Batch Size'][iBat]):
		NetGrowth[iS,:]=vi['GC']['Active'][iAge[iS],iS,:].copy()

	# Modify GY module estimates to represent variability in tree growth
	if meta[pNam]['Scenario'][iScn]['Gaia Status']=='On':
		NetGrowth=gaia.GYModelModifier(meta,pNam,NetGrowth,vo,iT)

	#--------------------------------------------------------------------------
	# *** < SPECIAL ORDER - CT Study (Johnston 2002) > ***
	if meta[pNam]['Project']['Code Project']=='Demo_Harv_ThinDensePine':
		GN_Total=NetGrowth[:,iGY['StemMerch']]+NetGrowth[:,iGY['StemNonMerch']]
		if (vi['tv'][iT]<1952):
			# Alter the ratio of merch to non-merch
			NetGrowth[:,iGY['StemMerch']]=0.3*GN_Total
			NetGrowth[:,iGY['StemNonMerch']]=0.7*GN_Total
		if (vi['tv'][iT]>1952) & (vi['tv'][iT]<=1997):
			rBK=NetGrowth[:,iGY['Bark']]/NetGrowth[:,iGY['StemMerch']]
			rBR=NetGrowth[:,iGY['Branch']]/NetGrowth[:,iGY['StemMerch']]
			rF=NetGrowth[:,iGY['Foliage']]/NetGrowth[:,iGY['StemMerch']]
			if iScn==0:
				GG_StemMerch_SO=1.3
				GG_StemNonMerch_SO=0.29
				M_StemMerch_SO=0.5
				M_StemNonMerch_SO=0.88
				NetGrowth[:,iGY['StemMerch']]=GG_StemMerch_SO-M_StemMerch_SO
				NetGrowth[:,iGY['StemNonMerch']]=GG_StemNonMerch_SO-M_StemNonMerch_SO
				NetGrowth[:,iGY['Bark']]=rBK*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Branch']]=rBR*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Foliage']]=rF*NetGrowth[:,iGY['StemMerch']]
			elif iScn==1:
				GG_StemMerch_SO=1.38
				GG_StemNonMerch_SO=0.3
				M_StemMerch_SO=0.15
				M_StemNonMerch_SO=0.45
				NetGrowth[:,iGY['StemMerch']]=GG_StemMerch_SO-M_StemMerch_SO
				NetGrowth[:,iGY['StemNonMerch']]=GG_StemNonMerch_SO-M_StemNonMerch_SO
				NetGrowth[:,iGY['Bark']]=rBK*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Branch']]=rBR*NetGrowth[:,iGY['StemMerch']]
				NetGrowth[:,iGY['Foliage']]=rF*NetGrowth[:,iGY['StemMerch']]
	# *** < END SPECIAL ORDER - CT Study (Johnstone 2002) > ***
	#--------------------------------------------------------------------------

	# Net growth of whole stemwood
	NetGrowth_WholeStem=NetGrowth[:,iGY['StemMerch']]+NetGrowth[:,iGY['StemNonMerch']]

	# *** SPECIAL ORDER ***
	# Used to force alignment with CBM25 (so that DOM pools can be compared)
	#NetGrowth_WholeStem=NetGrowth_WholeStem-0.01*vo['A'][iT,:]
	# *** SPECIAL ORDER ***

	flg=1
	if flg==1:
		# Net growth of foliage
		NetGrowth[:,iGY['Foliage']]=NetGrowth_WholeStem*(meta['Param']['BEV']['BiomassAllometrySL']['Gf1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gf2']-meta['Param']['BEV']['BiomassAllometrySL']['Gf1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gf3']*vo['A'][iT,:]))

		# Net growth of branches
		NetGrowth[:,iGY['Branch']]=NetGrowth_WholeStem*(meta['Param']['BEV']['BiomassAllometrySL']['Gbr1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbr2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbr1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbr3']*vo['A'][iT,:]))
	
		# Net growth of bark
		NetGrowth[:,iGY['Bark']]=NetGrowth_WholeStem*(meta['Param']['BEV']['BiomassAllometrySL']['Gbk1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbk2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbk1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbk3']*vo['A'][iT,:]))

	# Exploration of spcies-specific allometry
	flg=0
	if flg==1:
		# np.minimum(50,p[0]+(p[1]-p[0])*np.exp(-p[2]*dHat['Age']))
		p=[9.79194082e-02, 2.19770107e+02,4.92025804e-01]
		fun=np.minimum(3,p[0]+(p[1]-p[0])*np.exp(-p[2]*vo['A'][iT,:]))
		NetGrowth[:,iGY['Bark']]=np.maximum(0.1,NetGrowth_WholeStem)*fun
	
		p=[0.08881094, 7.76655694, 0.22940918]
		fun=np.minimum(3,p[0]+(p[1]-p[0])*np.exp(-p[2]*vo['A'][iT,:]))
		NetGrowth[:,iGY['Branch']]=np.maximum(0.1,NetGrowth_WholeStem)*fun
	
		p=[0.05464871,47.95376571,0.35506178]
		fun=np.minimum(3,p[0]+(p[1]-p[0])*np.exp(-p[2]*vo['A'][iT,:]))
		NetGrowth[:,iGY['Foliage']]=np.maximum(0.1,NetGrowth_WholeStem)*fun

	# *** SPECIAL ORDER ***
	# Used to force alignment with CBM25 (so that DOM pools can be compared)
	#NetGrowth[:,iGY['Foliage']]=0.5*NetGrowth[:,iGY['Foliage']]
	# *** SPECIAL ORDER ***

	# Add net growth to output variable structure
	# Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions
	# of C_G_Net - don't change it.
	vo['C_G_Net_Reg_ByPool'][iT,:,0:5]=NetGrowth[:,0:5]

	# Total net growth of root biomass (Li et al. 2003, Eq. 4)
	G_Net_Root_Total=0.22*np.sum(vo['C_G_Net_Reg_ByPool'][iT,:,0:5],axis=1)

	# Fine root fraction should decline with total biomass, but estimating it 	d on
	# size causes a lot of problems. Think about creating new equation.
	vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
	vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total

	#--------------------------------------------------------------------------
	# Nutrient application effects to net growth
	#--------------------------------------------------------------------------

	if meta[pNam]['Project']['Biomass Module']!='gromo':
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
	vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['StemMerch']])
	vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['StemNonMerch']])

	# Foliage
	vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Foliage']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Foliage']])

	# Branches
	vo['C_Eco_ByPool'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Branch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Branch']])

	# Bark
	vo['C_Eco_ByPool'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Bark']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Bark']])

	# Coarse roots
	vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['RootCoarse']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']])

	# Fine roots
	vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['RootFine']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootFine']])

	#--------------------------------------------------------------------------
	# Add net growth to live merch volume and update other volume variables
	# Notes:
	# 1) Regualr mortality to be added to dead volume in DOM function ***
	# 2) Total volume updated at the start of disturfance event function.
	#--------------------------------------------------------------------------

	# Update live whole stem volume
	V_Gnet=vi['lsat']['Biomass to Volume CF']*NetGrowth_WholeStem
	vo['V_WholeStemLive'][iT,:]=np.maximum(0,vo['V_WholeStemLive'][iT-1,:]+V_Gnet)

	# Update live stemwood merchantable volume (Use carbon curve instead of the volume directly from TIPSY (ensure consistency))
	V_Gnet=vi['lsat']['Biomass to Volume CF']*NetGrowth[:,iGY['StemMerch']]
	vo['V_MerchLive'][iT,:]=np.maximum(0,vo['V_MerchLive'][iT-1,:]+V_Gnet)

	#--------------------------------------------------------------------------
	# Biomass turnover
	#--------------------------------------------------------------------------

	# Biomass loss due to regular mortality
	#flg_Mortality='Age-dependent'
	flg_Mortality='Constant'

	if flg_Mortality=='Constant':
		# Constant mortality rate
		vo['C_M_Reg_ByPool'][iT,:,0:7]=meta['Param']['BEV']['BiomassTurnover']['Mreg_Kurzetal2009']*vo['C_Eco_ByPool'][iT,:,0:7]
	elif flg_Mortality=='Age-dependent':
		bBT=meta['Param']['BEV']['BiomassTurnover']
		flg2=0
		if flg2==1:
			b0={'b1':0.1,'b2':5,'b3':20,'b4':2.5};
			A=np.arange(0,200,1); y=(b0['b1']*(1+((b0['b2']*(A/b0['b3'])**b0['b4']-1)/np.exp(A/b0['b3']))))/100;
			plt.close('all'); plt.plot(A,y,'b-')
		fM=(bBT['Mreg1']*(1+((bBT['Mreg2']*(vo['A'][iT,:]/bBT['Mreg3'])**bBT['Mreg4']-1)/np.exp(vo['A'][iT,:]/bBT['Mreg3']))))/100
		vo['C_M_Reg_ByPool'][iT,:,0:7]=np.tile(fM,(7,1)).T*vo['C_Eco_ByPool'][iT,:,0:7]

	#--------------------------------------------------------------------------
	# *** SPECIAL ORDER - CT Study (Johnstone 2002) ***
	if meta[pNam]['Project']['Code Project']=='Demo_Harv_ThinDensePine':
		if (vi['tv'][iT]>1952) & (vi['tv'][iT]<=1997):
			if iScn==0:
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemMerch']]=M_StemMerch_SO
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemNonMerch']]=M_StemNonMerch_SO
			elif iScn==1:
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemMerch']]=M_StemMerch_SO
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemNonMerch']]=M_StemNonMerch_SO
	# *** SPECIAL ORDER - CT Study (Johnstone 2002) ***
	#--------------------------------------------------------------------------

	if meta[pNam]['Project']['Biomass Module']!='gromo':
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
	C_M_Reg_ByPool=vo['C_M_Reg_ByPool'][iT,:,:].copy()

	flg=0
	if flg==1:

		# Define a threshold level of negative growth so it isn't triggered by noise
		NegGrowthThreshold=-0.25
		#NegGrowthThreshold=-1e6

		# Find stands with negative net growth
		iNegNetG=np.where(vo['C_G_Net_Reg_ByPool'][iT,:,0]<NegGrowthThreshold)[0]

		if iNegNetG.size>0:

			# If it is the first instance of negative net growth: 91) record net growth
			# of the preceeding timestep and (2) set flag = 1.
			iSwitchFlag=np.where(meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG]==0)[0]
			meta[pNam]['Project']['FlagNegNetGrowth'][iNegNetG[iSwitchFlag]]=1
			meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG[iSwitchFlag],0:7]=vo['C_G_Net_Reg_ByPool'][iT-1,iNegNetG[iSwitchFlag],0:7]

			d=vo['C_G_Net_Reg_ByPool'][iT,iNegNetG,0:7]-meta[pNam]['Project']['G_Net_PriorToBreakup'][iNegNetG,:]

			CToTransfer=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],7))
			CToTransfer[iNegNetG,:]=-1*d

			vo['C_G_Net_Reg_ByPool'][iT,:,0:7]=vo['C_G_Net_Reg_ByPool'][iT,:,0:7]+CToTransfer
			vo['C_M_Reg_ByPool'][iT,:,0:7]=vo['C_M_Reg_ByPool'][iT,:,0:7]+CToTransfer

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
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Foliage_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Foliage']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]

	# Calculate branch biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Branch_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Branch']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Branch']]

	# Calculate bark biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Bark_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Bark']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Bark']]

	# Calculate coarse root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootCoarse_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]

	# Calculate fine root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootFine_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['RootFine']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]

	if meta[pNam]['Project']['Biomass Module']!='gromo':
		# Adjust litterfall to account for N application response
		if meta['Modules']['NutrientApp']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Litterfall')

	#--------------------------------------------------------------------------
	# Update summary variables
	#--------------------------------------------------------------------------

	# Update gross growth
	vo['C_G_Gross_ByPool'][iT,:,:]=vo['C_G_Net_Reg_ByPool'][iT,:,:]+C_M_Reg_ByPool

	# Update NPP
	# *** Compiled after the simulation is complete. ***
	#vo['C_NPP'][iT,:,:]=vo['C_G_Net_Reg_ByPool'][iT,:,:]+C_M_Reg_ByPool+vo['C_LF_ByPool'][iT,:,:]

	return vo

#%%
def TreeBiomassDynamicsFromGROMO(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

	# Update stand age
	vo['A'][iT,:]=vo['A'][iT-1,:]+1

	# Extract incides for GY model data structure
	iGY=meta['Modules']['GYM']['GC Input Indices']

	# Predictor variables
	A=vo['A'][iT,:]
	Ogp=0
	DAB=0; DAD=0; DAF=0; DAP=0
	#Hi=1*(vi['lsat']['Harvest Year Comp2']>0)
	Hi=0
	PL=1

	#--------------------------------------------------------------------------
	# Gross growth of stemwood
	#--------------------------------------------------------------------------

	bG=meta['Param']['Raw']['GROMO']['gg3']['Param']['BE']
	zs=meta['Param']['Raw']['GROMO']['gg3']['ZscoreStats']

	Tn=(vi['ENV']['Tn'][iT,:]-zs['Tn']['mu'])/zs['Tn']['sig']
	Wn=(vi['ENV']['Wn'][iT,:]-zs['Wn']['mu'])/zs['Wn']['sig']
	Ta=(vi['ENV']['Ta'][iT,:]-zs['Ta']['mu'])/zs['Ta']['sig']
	Wa=(vi['ENV']['Wa'][iT,:]-zs['Wa']['mu'])/zs['Wa']['sig']
	Nd=(vi['ENV']['ND'][iT,:]-zs['Nd']['mu'])/zs['Nd']['sig']
	Ca=(vi['ENV']['CO2'][iT]-zs['Ca']['mu'])/zs['Ca']['sig']

	B1=bG['b1']+bG['Ogp1']*Ogp+bG['Hi1']*Hi+bG['pl1']*PL+ \
		bG['dab1']*DAB+bG['dad1']*DAD+bG['daf1']*DAF+bG['dap1']*DAP+ \
		bG['Tn1']*Tn+bG['Wn1']*Wn+bG['Ta1']*Ta+bG['Wa1']*Wa+ \
		bG['TnTa1']*Tn*Ta+bG['WnWa1']*Wn*Wa+bG['Nd1']*Nd+bG['Ca1']*Ca
	B2=bG['b2']
	B3=bG['b3']
	B4=bG['b4']
	B5=bG['b5']+bG['Ogp5']*Ogp+bG['Hi5']*Hi+bG['pl5']*PL+bG['Tn5']*Tn+bG['Wn5']*Wn+bG['Ta5']*Ta+bG['Wa5']*Wa+bG['TnTa5']*Tn*Ta+bG['WnWa5']*Wn*Wa+bG['Nd5']*Nd+bG['Ca5']*Ca

	G_Gross_StemTot=(B1*(1+((B2*(A/B3)**B4-1)/np.exp(A/B3))))+(B5*A/100)
	#print(np.mean(G_Gross_StemTot))

	G_Gross_StemTot=np.maximum(0.1,G_Gross_StemTot)

	# Adjustment of stocking in harvested areas
	#G_Gross_StemTot=G_Gross_StemTot+0.60*Hi

	# Non-merch to merch ratio
	#Vtot_lag1=(vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']])/meta['Param']['BEV']['Biophysical']['Density Wood']/meta['Param']['BEV']['Biophysical']['Carbon Content Wood']
	fStemNonMerch=0.08#np.minimum(1,(1-0.086)*0.086+(1/np.exp(0.012*(Vtot_lag1-14.4))))

	# Update gross growth
	vo['C_G_Gross_ByPool'][iT,:,iEP['StemMerch']]=(1-fStemNonMerch)*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['StemNonMerch']]=fStemNonMerch*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['Bark']]=0.14*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['Branch']]=0.27*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['Foliage']]=0.17*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['RootCoarse']]=0.14*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['RootFine']]=0.17*G_Gross_StemTot

	#--------------------------------------------------------------------------
	# Mortality of stemwood
	#--------------------------------------------------------------------------

	bM=meta['Param']['Raw']['GROMO']['m2']['Param']['BE']
	zs=meta['Param']['Raw']['GROMO']['m2']['ZscoreStats']
	DT=1

	Csw_Tot_Lag1=vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']]
	A_z=(A-zs['A']['mu'])/zs['A']['sig']
	Tn_z=(vi['ENV']['Tn'][iT,:]-zs['Tn']['mu'])/zs['Tn']['sig']
	Wn_z=(vi['ENV']['Wn'][iT,:]-zs['Wn']['mu'])/zs['Wn']['sig']
	Ta_z=(vi['ENV']['Ta'][iT,:]-zs['Ta']['mu'])/zs['Ta']['sig']
	Wa_z=(vi['ENV']['Wa'][iT,:]-zs['Wa']['mu'])/zs['Wa']['sig']
	Nd_z=(vi['ENV']['ND'][iT,:]-zs['Nd']['mu'])/zs['Nd']['sig']

	M_Reg_StemTot=bM['Intercept']+bM['C(LS)[T.PL]']*PL+bM['OGP']*Ogp+bM['HI']*Hi+bM['DT']*DT+bM['B']*Csw_Tot_Lag1+bM['A_z']*A_z+ \
		bM['Tn_z']*Tn_z+bM['Wn_z']*Wn_z+bM['Ta_z']*Ta_z+bM['Wa_z']*Wa_z+bM['Nd_z']*Nd_z+ \
		bM['B:C(LS)[T.PL]']+bM['Tn_z:Ta_z']+bM['Wn_z:Wa_z']+bM['B:Tn_z']+bM['B:Wn_z']+bM['B:Ta_z']+bM['B:Wa_z']

	M_Reg_StemTot=np.maximum(0,M_Reg_StemTot)

	# Constrain mortality by the carbon pool at the start of the time step
	M_Reg_StemTot=np.minimum(M_Reg_StemTot,Csw_Tot_Lag1)

	#--------------------------------------------------------------------------
	# Net growth
	#--------------------------------------------------------------------------

	# Net growth by pool
	NetGrowthByPool=np.zeros( (meta[pNam]['Project']['Batch Size'][iBat],6) )

	# Net growth of total stemwood
	G_Net_StemTot=G_Gross_StemTot-M_Reg_StemTot

	# Net growth of merch stemwood and non-merch stemwood
	#Vtot_lag1=(vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']])/meta['Param']['BEV']['Biophysical']['Density Wood']/meta['Param']['BEV']['Biophysical']['Carbon Content Wood']
	fStemNonMerch=0.08 #np.minimum(1,(1-0.086)*0.086+(1/np.exp(0.012*(Vtot_lag1-14.4))))
	NetGrowthByPool[:,iEP['StemMerch']]=(1-fStemNonMerch)*G_Net_StemTot
	NetGrowthByPool[:,iEP['StemNonMerch']]=fStemNonMerch*G_Net_StemTot

	# Net growth of foliage
	NetGrowthByPool[:,iGY['Foliage']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gf1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gf2']-meta['Param']['BEV']['BiomassAllometrySL']['Gf1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gf3']*vo['A'][iT,:]))

	# Net growth of branches
	NetGrowthByPool[:,iGY['Branch']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gbr1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbr2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbr1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbr3']*vo['A'][iT,:]))

	# Net growth of bark
	NetGrowthByPool[:,iGY['Bark']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gbk1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbk2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbk1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbk3']*vo['A'][iT,:]))

	# Add net growth to output variable structure
	# Oddly, using meta['iEP']['BiomassAboveground'] will invert the dimensions
	# of C_G_Net - don't change it.
	vo['C_G_Net_Reg_ByPool'][iT,:,0:5]=NetGrowthByPool[:,0:5]

	# Total net growth of root biomass (Li et al. 2003, Eq. 4)
	G_Net_Root_Total=0.22*np.sum(vo['C_G_Net_Reg_ByPool'][iT,:,0:5],axis=1)

	# Fine root fraction should decline with total biomass, but estimating it based on
	# size causes a lot of problems. Think about creating new equation.
	vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']]=(1-0.072)*G_Net_Root_Total
	vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootFine']]=0.072*G_Net_Root_Total

	#--------------------------------------------------------------------------
	# Nutrient application effects to net growth
	#--------------------------------------------------------------------------

	if meta[pNam]['Project']['Biomass Module']!='gromo':
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
	vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['StemMerch']])
	vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['StemNonMerch']])

	# Foliage
	vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Foliage']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Foliage']])

	# Branches
	vo['C_Eco_ByPool'][iT,:,iEP['Branch']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Branch']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Branch']])

	# Bark
	vo['C_Eco_ByPool'][iT,:,iEP['Bark']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['Bark']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['Bark']])

	# Coarse roots
	vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['RootCoarse']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']])

	# Fine roots
	vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]=np.maximum(0,vo['C_Eco_ByPool'][iT-1,:,iEP['RootFine']]+vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootFine']])

	#--------------------------------------------------------------------------
	# Add net growth to live merch volume and update other volume variables
	# Notes:
	# 1) Regualr mortality to be added to dead volume in DOM function ***
	# 2) Total volume updated at the start of disturfance event function.
	#--------------------------------------------------------------------------

	# Update live stemwood merchantable volume
	vo['V_MerchLive'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]
	vo['V_WholeStemLive'][iT,:]=vi['lsat']['Biomass to Volume CF']*np.sum(vo['C_Eco_ByPool'][iT,:,iEP['StemTotal']],axis=0)

	#--------------------------------------------------------------------------
	# Biomass turnover
	#--------------------------------------------------------------------------

	# Biomass loss due to regular mortality
	vo['C_M_Reg_ByPool'][iT,:,0:7]=vo['C_G_Gross_ByPool'][iT,:,0:7]-vo['C_G_Net_Reg_ByPool'][iT,:,0:7]

	# Old:
	# 	bBT=meta['Param']['BEV']['BiomassTurnover']
	# 	flg=0
	# 	if flg==1:
	# 		b0={'b1':0.1,'b2':5,'b3':20,'b4':2.5};
	# 		A=np.arange(0,200,1); y=(b0['b1']*(1+((b0['b2']*(A/b0['b3'])**b0['b4']-1)/np.exp(A/b0['b3']))))/100; plt.close('all'); plt.plot(A,y,'b-')
	# 	fM=(bBT['Mreg1']*(1+((bBT['Mreg2']*(vo['A'][iT,:]/bBT['Mreg3'])**bBT['Mreg4']-1)/np.exp(vo['A'][iT,:]/bBT['Mreg3']))))/100
	# 	vo['C_M_Reg_ByPool'][iT,:,0:7]=np.tile(fM,(7,1)).T*vo['C_Eco_ByPool'][iT,:,0:7]

	# *** SPECIAL ORDER - CT Study (Johnstone 2002) ***
	if meta[pNam]['Project']['Code Project']=='Demo_Harv_ThinDensePine':
		if iScn==0:
			if (vi['tv'][iT]>1952) & (vi['tv'][iT]<=1997):
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemMerch']]=1.08
		elif iScn==1:
			if (vi['tv'][iT]>1952) & (vi['tv'][iT]<=1997):
				vo['C_M_Reg_ByPool'][iT,:,iEP['StemMerch']]=0.40

	if meta[pNam]['Project']['Biomass Module']!='gromo':
		# Adjust mortality to account for N application response
		if meta['Modules']['NutrientApp']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Mortality')

	#--------------------------------------------------------------------------
	# Litterfall
	#--------------------------------------------------------------------------

	# Setting turnover as a function of age will decouple NPP from net growth.
	Aref=meta['Param']['BEV']['BiomassTurnover']['BiomassTurnoverAgeRef']
	fA=-meta['Param']['BEV']['BiomassTurnover']['BiomassTurnoverAgeDependence']

	# Calculate merch stemwood biomass turnover due to litterfall
	# This is zero.
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['StemMerch_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['StemMerch']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]

	# Calculate non-merch stemwood biomass turnover due to litterfall
	# This is zero.
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['StemNonMerch_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['StemNonMerch']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]

	# Calculate foliage biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Foliage_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Foliage']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]

	# Calculate branch biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Branch_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Branch']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Branch']]

	# Calculate bark biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['Bark_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['Bark']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['Bark']]

	# Calculate coarse root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootCoarse_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]

	# Calculate fine root biomass turnover due to litterfall
	tr=np.maximum(0.001,fA*(vo['A'][iT,:]-Aref)+meta['Param']['BEV']['BiomassTurnover']['RootFine_RateLF'])
	vo['C_LF_ByPool'][iT,:,iEP['RootFine']]=tr*vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]

	if meta[pNam]['Project']['Biomass Module']!='gromo':
		# Adjust litterfall to account for N application response
		if meta['Modules']['NutrientApp']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Litterfall')

	return vo

#%% Dead organic matter and soil organic matter dynamics
def DeadWoodLitterAndSoilDynamics(meta,pNam,iT,iBat,vi,vo,iEP):

	# Extract parameters
	bBT=meta['Param']['BEV']['BiomassTurnover']
	bD=meta['Param']['BEV']['Decomposition']

	#--------------------------------------------------------------------------
	# Update pools of current time step
	#--------------------------------------------------------------------------

	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['DeadStemMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['DeadStemNonMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]=vo['C_Eco_ByPool'][iT-1,:,iEP['DeadBark']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['DeadBranch']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT-1,:,iEP['LitterVF']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT-1,:,iEP['LitterF']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT-1,:,iEP['LitterM']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT-1,:,iEP['LitterS']]
	vo['C_Eco_ByPool'][iT,:,iEP['SoilVF']]=vo['C_Eco_ByPool'][iT-1,:,iEP['SoilVF']]
	vo['C_Eco_ByPool'][iT,:,iEP['SoilF']]=vo['C_Eco_ByPool'][iT-1,:,iEP['SoilF']]
	vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]=vo['C_Eco_ByPool'][iT-1,:,iEP['SoilS']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledStemMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledStemNonMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledBranch']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledBark']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledDeadStem']]
	vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]=vo['C_Eco_ByPool'][iT-1,:,iEP['PiledDeadBranch']]

	#--------------------------------------------------------------------------
	# Transfer litterfall carbon to DOM pools
	#--------------------------------------------------------------------------

	# Transfer biomass turnover to very fast litter pool
	ToLitterVF=bBT['FoliageLitToLitterVF']*vo['C_LF_ByPool'][iT,:,iEP['Foliage']]+ \
		bBT['RootFineLitToLitterVF']*vo['C_LF_ByPool'][iT,:,iEP['RootFine']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]+ToLitterVF

	# Transfer biomass turnover to fast litter pool
	ToLitterF=bBT['BarkLitToLitterF']*vo['C_LF_ByPool'][iT,:,iEP['Bark']]+ \
		bBT['BranchLitToLitterF']*vo['C_LF_ByPool'][iT,:,iEP['Branch']]+ \
		bBT['RootCoarseLitToLitterF']*vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+ToLitterF

	# Transfer biomass turnover to medium litter pool
	ToLitterM=bBT['BarkLitToLitterM']*vo['C_LF_ByPool'][iT,:,iEP['Bark']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+ToLitterM

	# Transfer biomass turnover to slow litter pool
	ToLitterS=bBT['BarkLitToLitterS']*vo['C_LF_ByPool'][iT,:,iEP['Bark']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]+ToLitterS

	#--------------------------------------------------------------------------
	# Transfer mortality carbon to DOM pools
	#--------------------------------------------------------------------------

	ToDeadStemMerch=bBT['StemMerchMorToDeadStemMerch']*vo['C_M_Reg_ByPool'][iT,:,iEP['StemMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]+ToDeadStemMerch

	ToDeadStemNonMerch=bBT['StemNonMerchMorToDeadStemNonMerch']*vo['C_M_Reg_ByPool'][iT,:,iEP['StemNonMerch']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]+ToDeadStemNonMerch

	ToDeadBark=bBT['BarkMorToDeadBark']*vo['C_M_Reg_ByPool'][iT,:,iEP['Bark']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]+ToDeadBark

	ToDeadBranch=bBT['BranchMorToDeadBranch']*vo['C_M_Reg_ByPool'][iT,:,iEP['Branch']]
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]+ToDeadBranch

	ToLitterVF=bBT['FoliageMorToLitterVF']*vo['C_M_Reg_ByPool'][iT,:,iEP['Foliage']]
	vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]+ToLitterVF

	ToSoilVF=bBT['RootCoarseMorToSoilVF']*vo['C_M_Reg_ByPool'][iT,:,iEP['RootCoarse']]+ \
		bBT['RootFineLitToSoilVF']*vo['C_LF_ByPool'][iT,:,iEP['RootFine']] +\
		bBT['RootFineMorToSoilVF']*vo['C_M_Reg_ByPool'][iT,:,iEP['RootFine']]
	vo['C_Eco_ByPool'][iT,:,iEP['SoilVF']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilVF']]+ToSoilVF

	ToSoilF=bBT['RootCoarseLitToSoilF']*vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]+ \
		bBT['RootCoarseMorToSoilF']*vo['C_M_Reg_ByPool'][iT,:,iEP['RootCoarse']]
	vo['C_Eco_ByPool'][iT,:,iEP['SoilF']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilF']]+ToSoilF

	ToSoilS=np.zeros(ToSoilF.size)
	vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]+ToSoilS

	# QA check that litterfall = transfers to DOM
	flg=0
	if flg==1:
		iS=0
		a=np.sum(vo['C_LF_ByPool'][iT,iS,:]+vo['C_M_Reg_ByPool'][iT,iS,:])
		b=np.sum(ToLitterVF[iS]+ToLitterF[iS]+ToLitterM[iS]+ToLitterS[iS]+ToDeadStemMerch[iS]+ToDeadStemNonMerch[iS]+ToDeadBranch[iS]+ToSoilVF[iS]+ToSoilF[iS]+ToSoilS[iS])
		print(str(a) + ' ' + str(b) + ' '  + str(a-b))

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

	# Heterotrophic consumption
	# This is put in module "HC" dictionary so that it can be passed to the
	# nutrient application module for adjustment.
	for k in meta['Core']['Name Pools Dead']:
		meta[pNam]['Project']['HC'][k]=bD[k + '_R10']*vo['C_Eco_ByPool'][iT,:,iEP[k]].flatten()*bD[k + '_Q10']**fT

	HC_Total=0
	for k in meta[pNam]['Project']['HC'].keys():
		HC_Total=HC_Total+meta[pNam]['Project']['HC'][k][0,0]

	if meta[pNam]['Project']['Biomass Module']!='gromo':
		# Adjust decomposition to account for N application response
		if meta['Modules']['NutrientApp']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'HeterotrophicRespiration')

	# Remove heterotrophic consumption of carbon from source DOM pools
	for k in meta['Core']['Name Pools Dead']:
		vo['C_Eco_ByPool'][iT,:,iEP[k]]=vo['C_Eco_ByPool'][iT,:,iEP[k]]-meta[pNam]['Project']['HC'][k]

	# Heterotrophic respiration (emissions of CO2 from decomposition)
	rh={}
	RH_Total=0
	for k in meta['Core']['Name Pools Dead']:
		vo['C_RH_ByPool'][iT,:,iEP[k]]=bD[k + 'ToCO2']*meta[pNam]['Project']['HC'][k]

		rh[k]=bD[k + 'ToCO2']*meta[pNam]['Project']['HC'][k]
		RH_Total=RH_Total+rh[k][0,0]

	# Re-distribute the non-respired proportion of heterotrphic consumption to sink DOM pools

	DeadStemMerchToLitterM=bD['DeadStemMerchToLitterM']*meta[pNam]['Project']['HC']['DeadStemMerch']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+DeadStemMerchToLitterM

	DeadStemNonMerchToLitterM=bD['DeadStemNonMerchToLitterM']*meta[pNam]['Project']['HC']['DeadStemNonMerch']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+DeadStemNonMerchToLitterM

	DeadBarkToLitterF=bD['DeadBarkToLitterF']*meta[pNam]['Project']['HC']['DeadBark']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+DeadBarkToLitterF

	DeadBranchToLitterF=bD['DeadBranchToLitterF']*meta[pNam]['Project']['HC']['DeadBranch']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+DeadBranchToLitterF

	LitterVFToLitterS=bD['LitterVFToLitterS']*meta[pNam]['Project']['HC']['LitterVF']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]+LitterVFToLitterS

	LitterFToLitterS=bD['LitterFToLitterS']*meta[pNam]['Project']['HC']['LitterF']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]+LitterFToLitterS

	LitterMToLitterS=bD['LitterMToLitterS']*meta[pNam]['Project']['HC']['LitterM']
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]+LitterMToLitterS

	SoilVFToSoilS=bD['SoilVFToSoilS']*meta[pNam]['Project']['HC']['SoilVF']
	vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]+SoilVFToSoilS

	SoilFToSoilS=bD['SoilFToSoilS']*meta[pNam]['Project']['HC']['SoilF']
	vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]+SoilFToSoilS

	# QA - conservation of mass
	flg=0
	if flg==1:
		Red_Total=DeadBranchToLitterF+DeadStemMerchToLitterM+DeadStemNonMerchToLitterM+ \
			LitterVFToLitterS+LitterFToLitterS+LitterMToLitterS+ \
			SoilVFToSoilS+SoilFToSoilS

		print( str(HC_Total) + ' ' + str(RH_Total) + ' ' + str(Red_Total[0,0]) + ' ' + str(HC_Total-RH_Total-Red_Total[0,0]) )

	# Physical transfers
	PT_DeadStemMerch=bD['DeadStemMerch_PhysTransRate']*vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]
	PT_DeadStemNonMerch=bD['DeadStemNonMerch_PhysTransRate']*vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]
	PT_DeadBranch=bD['DeadBranch_PhysTransRate']*vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]
	PT_DeadBark=bD['DeadBark_PhysTransRate']*vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]
	PT_LitterS=bD['LitterS_PhysTransRate']*vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]

	# Remove carbon that is physically transferred
	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]-PT_DeadStemMerch
	vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]-PT_DeadStemNonMerch
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]-PT_DeadBark
	vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]-PT_DeadBranch
	vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterS']]-PT_LitterS

	# Add decomposed carbon to more decomposed DOM pools
	vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilS']]+PT_LitterS
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+PT_DeadBranch
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+PT_DeadBark
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+PT_DeadStemMerch
	vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+0.5*PT_DeadStemNonMerch
	vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+0.5*PT_DeadStemNonMerch

	#--------------------------------------------------------------------------
	# Update dead volume
	#--------------------------------------------------------------------------
	vo['V_MerchDead'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]
	vo['V_WholeStemDead'][iT,:]=vi['lsat']['Biomass to Volume CF']*(vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]+vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']])

	#--------------------------------------------------------------------------
	# Track Fate of Roots and Dispersed Slash to DOM
	#--------------------------------------------------------------------------
	if meta[pNam]['Project']['Track fate of felled to DOM Status']=='On':
		mnam='Fate of Roots and Dispersed Slash to DOM'

		# Update the pools
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,:]=meta['Modules'][mnam]['C_Eco_ByPool'][iT-1,:,:]

		# Respiration rate - Note that these terms do not equal the
		# atmosphere-bound efflux from heterotrophic respiration as a fraction is
		# emitted to the atmosphere and the remaining fraction is reorganized
		# within the ecosystem.
		R_DeadStemMerch=bD['DeadStemMerch_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']].flatten()*bD['DeadStemMerch_Q10']**fT
		R_DeadStemNonMerch=bD['DeadStemNonMerch_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']].flatten()*bD['DeadStemNonMerch_Q10']**fT
		R_DeadBranch=bD['DeadBranch_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']].flatten()*bD['DeadBranch_Q10']**fT
		R_LitterVF=bD['LitterVF_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterVF']].flatten()*bD['LitterVF_Q10']**fT
		R_LitterF=bD['LitterF_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']].flatten()*bD['LitterF_Q10']**fT
		R_LitterM=bD['LitterM_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']].flatten()*bD['LitterM_Q10']**fT
		R_LitterS=bD['LitterS_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']].flatten()*bD['LitterS_Q10']**fT
		R_SoilVF=bD['SoilVF_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilVF']].flatten()*bD['SoilVF_Q10']**fT
		R_SoilF=bD['SoilF_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilF']].flatten()*bD['SoilF_Q10']**fT
		R_SoilS=bD['SoilS_R10']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']].flatten()*bD['SoilS_Q10']**fT

		# Remove respired carbon from source DOM pools
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterVF']]-R_LitterVF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]-R_LitterF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]-R_LitterM
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]-R_LitterS
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]-R_DeadStemMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]-R_DeadStemNonMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]-R_DeadBranch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilVF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilVF']]-R_SoilVF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilF']]-R_SoilF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]-R_SoilS

		# Re-define decayed fast litter
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]+bD['DeadBranchToLitterF']*R_DeadBranch

		# Re-define decayed medium litter
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]+bD['DeadStemMerchToLitterM']*R_DeadStemMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]+bD['DeadStemNonMerchToLitterM']*R_DeadStemNonMerch

		# Re-define decayed slow litter
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]+bD['LitterVFToLitterS']*R_LitterVF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]+bD['LitterFToLitterS']*R_LitterF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]+bD['LitterMToLitterS']*R_LitterM

		# Re-define decayed slow soil
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]+bD['SoilVFToSoilS']*R_SoilVF
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]+bD['SoilFToSoilS']*R_SoilF

		# Heterotrophic respiration
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['LitterVF']]=bD['LitterVFToCO2']*R_LitterVF
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['LitterF']]=bD['LitterFToCO2']*R_LitterF
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['LitterM']]=bD['LitterMToCO2']*R_LitterM
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['LitterS']]=bD['LitterSToCO2']*R_LitterS
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['DeadStemMerch']]=bD['DeadStemMerchToCO2']*R_DeadStemMerch
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['DeadStemNonMerch']]=bD['DeadStemNonMerchToCO2']*R_DeadStemNonMerch
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['DeadBranch']]=bD['DeadBranchToCO2']*R_DeadBranch
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['SoilVF']]=bD['SoilVFToCO2']*R_SoilVF
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['SoilF']]=bD['SoilFToCO2']*R_SoilF
		meta['Modules'][mnam]['C_RH_ByPool'][iT,:,iEP['SoilS']]=bD['SoilSToCO2']*R_SoilS

		# Physical transfer
		PT_LitterS=bD['LitterS_PhysTransRate']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]
		PT_DeadStemMerch=bD['DeadStemMerch_PhysTransRate']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]
		PT_DeadStemNonMerch=bD['DeadStemNonMerch_PhysTransRate']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]
		PT_DeadBranch=bD['DeadBranch_PhysTransRate']*meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]

		# Remove carbon that is physically transferred
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterS']]-PT_LitterS
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]-PT_DeadStemMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]-PT_DeadStemNonMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]-PT_DeadBranch

		# Add decomposed carbon to more decomposed DOM pools
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['SoilS']]+PT_LitterS
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]+PT_DeadBranch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]+PT_DeadStemMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterF']]+0.5*PT_DeadStemNonMerch
		meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,:,iEP['LitterM']]+0.5*PT_DeadStemNonMerch

	return vo

#%% Disturbance and management events
def DisturbanceAndManagementEvents(meta,pNam,iT,iScn,iEns,iBat,vi,vo,iEP):

	# Track initial carbon pools for QA purposes
	flg=0
	if flg==1:
		Biomass0=np.sum(vo['C_Eco_ByPool'][iT,0,iEP['BiomassTotal']])
		Litter0=np.sum(vo['C_Eco_ByPool'][iT,0,iEP['Litter']])
		Soil0=np.sum(vo['C_Eco_ByPool'][iT,0,iEP['Soil']])
		DeadWood0=np.sum(vo['C_Eco_ByPool'][iT,0,iEP['DeadWood']])
		Piles0=np.sum(vo['C_Eco_ByPool'][iT,0,iEP['Piled']])

	# Update total (live+dead) stemwood merchantable volume
	vo['V_MerchTotal'][iT,:]=vo['V_MerchLive'][iT,:]+vo['V_MerchDead'][iT,:]
	vo['V_WholeStemTotal'][iT,:]=vo['V_WholeStemLive'][iT,:]+vo['V_WholeStemDead'][iT,:]

	# Keep tabs on whether a catastrophic event has been added
	flag_sim_event=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])

	# Predict wildfire (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Status']=='On') & (vi['tv'][iT]<1920):
		vi,flag_sim_event=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['Wildfire Sim Obs Status']=='On') & (vi['tv'][iT]>=1920) & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)

	# Predict mountain pine beetle (on the fly)
	if (meta[pNam]['Scenario'][iScn]['IBM Sim Pre-obs Status']=='On') & (vi['tv'][iT]<1950):
		vi,flag_sim_event=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['IBM Sim Obs Status']=='On') & (vi['tv'][iT]>=1950) & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['IBM Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event)

	# Predict disease (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Disease Sim Historical Status']=='On') & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictDisease_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['A'][iT,:],flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['Disease Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictDisease_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['A'][iT,:],flag_sim_event)

	# Predict wind (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Wind Sim Historical Status']=='On') & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:],flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['Wind Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:],flag_sim_event)

	# Predict frost (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Frost Sim Historical Status']=='On') & (vi['tv'][iT]<meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictFrost_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:],flag_sim_event)
	if (meta[pNam]['Scenario'][iScn]['Frost Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi,flag_sim_event=asm.PredictFrost_OnTheFly(meta,pNam,vi,iT,iEns,vo['A'][iT,:],flag_sim_event)

	# Predict harvesting (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Harvest Sim Historical Status']=='On') & (vi['tv'][iT]<meta[pNam]['Scenario'][iScn]['Harvest Sim Year Transition']):
		Period='Historical'
		vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['V_MerchTotal'][iT,:],Period)
	if (meta[pNam]['Scenario'][iScn]['Harvest Sim Future Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Scenario'][iScn]['Harvest Sim Year Transition']):
		Period='Future'
		vi=asm.PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,vo['V_MerchTotal'][iT,:],Period)

	# Predict future nutrient application (on the fly)
	if (meta[pNam]['Scenario'][iScn]['Nutrient App Sim Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi=napp.ScheduleNutrientApplication(meta,pNam,vi,vo,iT,iScn,iEns,iBat)
	
	# Predict future non-obligation stand establishment (on the fly)
	if (meta[pNam]['Scenario'][iScn]['NOSE Sim Status']=='On') & (vi['tv'][iT]>=meta[pNam]['Project']['Year Project']):
		vi=tft.PredictNOSE_OnTheFly(meta,pNam,iScn,iBat,vi,iT)

	# Check to see how many events occur in this time step (don't do more than necessary)
	NumEventsInTimeStep=np.sum(np.sum(vi['EC']['ID Event Type'][iT,:,:]>0,axis=0)>0)

	# Initialize indicator of aerial nutrient application
	flag_nutrient_app=np.zeros(meta[pNam]['Project']['Batch Size'][iBat])

	# Keep track of total tree biomass at the start of each annual time step
	C_Biomass_t0=np.sum(vo['C_Eco_ByPool'][iT,:,iEP['BiomassTotal']],axis=0)

	# Loop through events in year
	for iE in range(NumEventsInTimeStep):

		# Event type IDs for the iE'th event of the year
		ID_Type=vi['EC']['ID Event Type'][iT,:,iE].copy()

		# Indexes to each unique event type
		idx_Type=gu.IndicesFromUniqueArrayValues(ID_Type)

		# Total affected biomass carbon
		LiveImpactFactor=vi['EC']['Mortality Factor'][iT,:,iE].copy()

		# Define a seperate impact factor for dead carbon pools
		DeadImpactFactor=LiveImpactFactor.copy()

		# Age correction factor
		# By default, assume oldest trees were most affected, reduce age in prportion
		# with mortality rate. See exceptions below.
		AgeCorrection=1-LiveImpactFactor

		# Growth factor
		GrowthFactor=vi['EC']['Growth Factor'][iT,:,iE].copy()

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
		if meta['LUT']['Event']['Nutrient App Aerial'] in idx_Type:
			flag_nutrient_app[idx_Type[meta['LUT']['Event']['Nutrient App Aerial']]]=1

		#----------------------------------------------------------------------
		# Adjustments to harvest events
		#----------------------------------------------------------------------
		if meta['LUT']['Event']['Harvest'] in idx_Type:
			
			iHarvest=idx_Type[meta['LUT']['Event']['Harvest']]

			# Adjust event-specific parameters to reflect time- and region-specific
			# fate of felled material
			
			# Index to time-dependent fate of felled materials
			iT_P=np.where(meta['Param']['BE']['FelledFate']['Year']==meta[pNam]['Year'][iT])[0]

			# Simulations may exceed the timeframe of the felled fate parameters
			# If so, set to the last year
			if iT_P.size==0:
				iT_P=-1
			for k in meta['Param']['BEV']['FelledFate'].keys():
				b[k][iHarvest]=meta['Param']['BEV']['FelledFate'][k][iT_P,iHarvest]

			# Update age at harvest
			vo['A_Harvest'][iT,iHarvest]=vo['A'][iT,iHarvest]

			# Track occurrence of harvest for wildfire modelling
			mnam='Disturbance Effects on Wildfire Occurrence'
			meta['Modules'][mnam]['Harvest Flag'][iHarvest]=1

			# Track occurrence of salvage harvesting for wildfire modelling
			meta['Modules'][mnam]['Salvage Harvest Flag'][iHarvest]=meta['Modules'][mnam]['Mountain Pine Beetle Flag'][iHarvest]

			#--------------------------------------------------------------------------
			# Update albedo surface shortwave RF
			#--------------------------------------------------------------------------
			pSet='AlbedoSurfaceShortwaveRF_HarvestResponseByBGCZone'

			#vi['tv']=np.arange(1,2100)
			#iT=2022
			dt=np.tile(vi['tv']-vi['tv'][iT],(iHarvest.size,1)).T
			b0=np.tile(meta['Param']['BEV'][pSet]['Intercept'][iHarvest],(vi['tv'].size,1))
			b1=np.tile(meta['Param']['BEV'][pSet]['Slope'][iHarvest],(vi['tv'].size,1))
			b2=np.tile(meta['Param']['BEV'][pSet]['Initial'][iHarvest],(vi['tv'].size,1))
			yhat0=b0+b1*dt
			Response=b2*np.ones(dt.shape)
			Response[vi['tv']>=vi['tv'][iT],:]=np.minimum(b2[vi['tv']>=vi['tv'][iT],:],yhat0[vi['tv']>=vi['tv'][iT],:])
			Response=Response-np.tile(Response[iT-1],(vi['tv'].size,1))
			vo['RF_AlbedoSurfaceShortwave'][iT:,iHarvest]=vo['RF_AlbedoSurfaceShortwave'][iT:,iHarvest]+Response[iT:,:]

# 			yhat0=beta[zone][0]+beta[zone][1]*(tv-t_harv)
# 			yhat=beta[zone][2]*np.ones(tv.size)
# 			yhat[tv>=t_harv]=np.minimum(beta[zone][2],yhat0[tv>=t_harv])
# 			yhat=yhat-yhat[0]

		#----------------------------------------------------------------------
		# Adjustments to pre-commercial thinning events
		#----------------------------------------------------------------------
		if meta['LUT']['Event']['PCT'] in idx_Type:

			iPCT=idx_Type[meta['LUT']['Event']['PCT']]

			# Age correction
			AgeCorrection[iPCT]=1.0

		#----------------------------------------------------------------------
		# Adjustments to Mountain Pine Beetle events
		#----------------------------------------------------------------------
		if meta['LUT']['Event']['Mountain Pine Beetle'] in idx_Type:

			iIBM=idx_Type[meta['LUT']['Event']['Mountain Pine Beetle']]

			# Age is less affected than biomass (see field plot summary)
			AgeCorrection[iIBM]=1-0.6*LiveImpactFactor[iIBM]
			#AgeCorrection[iIBM]=1-LiveImpactFactor[iIBM]
			#AgeCorrection[iIBM]=1.0

			# MPB detection in young stands does not have impact
			LiveImpactFactor[ iIBM[ vo['A'][iT,iIBM]<20 ] ]=0
			AgeCorrection[ iIBM[ vo['A'][iT,iIBM]<20 ] ]=1.0

			# Track occurrence of IBM for wildfire modelling
			meta['Modules']['Disturbance Effects on Wildfire Occurrence']['Mountain Pine Beetle Flag'][iIBM]=1

		#----------------------------------------------------------------------
		# Adjustments to wildfire events
		#----------------------------------------------------------------------
		if meta['LUT']['Event']['Wildfire'] in idx_Type:
			iWF=idx_Type[meta['LUT']['Event']['Wildfire']]

			# Track occurrence of IBM for wildfire modelling
			meta['Modules']['Disturbance Effects on Wildfire Occurrence']['Wildfire Flag'][iWF]=1

		#----------------------------------------------------------------------
		# Adjustments to wind events
		#----------------------------------------------------------------------
		if meta['LUT']['Event']['Wind'] in idx_Type:
			iWind=idx_Type[meta['LUT']['Event']['Wind']]

			# Reduce mortality from wind (the 'dead factor' will affect dead trees,
			# but the effect on live trees is reduced)
			#LiveImpactFactor[iWind]=LiveImpactFactor[iWind]
			LiveImpactFactor[iWind]=0.0
			AgeCorrection[iWind]=1.0#*LiveImpactFactor[iWind]

		#----------------------------------------------------------------------
		# Net-down insect mortality to reflect the proportion of host species
		# This was perhaps also very slow and inefficent!
		#----------------------------------------------------------------------
		if meta[pNam]['Project']['Special Attribution Method']!='NOSE':
			# *** This is causing problems in the NOSE project - exclude from that project ***
			for nam in meta['Param']['Raw']['DisturbanceSpeciesAffected']['Insect Name']:
				id=meta['LUT']['Event'][nam]
				if id in idx_Type.keys():
					# Net down factor
					ndf=vi['lsat']['Insect Mortality Percent Tree Species Affected'][nam].astype('float')/100
					ind=idx_Type[id]
					LiveImpactFactor[ind]=LiveImpactFactor[ind]*ndf[ind]

		#----------------------------------------------------------------------
		# Define the amount of each pool that is affected by the event
		#----------------------------------------------------------------------

		# Affected biomass carbon
		if meta[pNam]['Project']['Biomass Module']=='Sawtooth':
			Affected_StemMerch=vo['C_M_Tot'][iT,:,iEP['StemMerch']]
			Affected_StemNonMerch=vo['C_M_Tot'][iT,:,iEP['StemNonMerch']]
			Affected_Foliage=vo['C_M_Tot'][iT,:,iEP['Foliage']]
			Affected_Branch=vo['C_M_Tot'][iT,:,iEP['Branch']]
			Affected_Bark=vo['C_M_Tot'][iT,:,iEP['Bark']]
			Affected_RootCoarse=vo['C_M_Tot'][iT,:,iEP['RootCoarse']]
			Affected_RootFine=vo['C_M_Tot'][iT,:,iEP['RootFine']]
		else:
			Affected_StemMerch=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]
			Affected_StemNonMerch=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]
			Affected_Bark=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['Bark']]
			Affected_Branch=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['Branch']]
			Affected_Foliage=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]
			Affected_RootCoarse=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]
			Affected_RootFine=LiveImpactFactor*vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]

		# Detect mortality factors > 1.0
		#ind=np.where(LiveImpactFactor>1)[0]
		#if ind.size>0:
		#	print('Warning, mortality factor > 1.0 detected.')

		# Affected all biomass
		Affected_BiomassTotal=Affected_StemMerch+Affected_StemNonMerch+Affected_Foliage+ \
			Affected_Branch+Affected_Bark+Affected_RootCoarse+Affected_RootFine

		# Partition bark into merch and non-merch components (done for removals)
		#MerchFraction=np.nan_to_num(np.minimum(1,np.maximum(0,Affected_StemMerch/(Affected_StemMerch+Affected_StemNonMerch))))
		#Affected_BarkMerch=Affected_Bark*MerchFraction
		#Affected_BarkNonMerch=Affected_Bark*(1-MerchFraction)

		# Affected dead wood
		Affected_DeadStemMerch=DeadImpactFactor*(1-b['DeadStemMerchToDeadStemMerch'])*vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]
		Affected_DeadStemNonMerch=DeadImpactFactor*(1-b['DeadStemNonMerchToDeadStemNonMerch'])*vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]
		Affected_DeadBark=DeadImpactFactor*(1-b['DeadBarkToDeadBark'])*vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]
		Affected_DeadBranch=DeadImpactFactor*(1-b['DeadBranchToDeadBranch'])*vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]
		Affected_DeadWood=Affected_DeadStemMerch+Affected_DeadStemNonMerch+Affected_DeadBark+Affected_DeadBranch

		# All biomass and dead wood
		if meta['LUT']['Event']['Harvest'] in idx_Type:
			Affected_All=Affected_BiomassTotal+Affected_DeadStemMerch+Affected_DeadStemNonMerch+ \
				Affected_DeadBranch+Affected_DeadBark
			vo['C_Felled'][iT,iHarvest]=Affected_All[iHarvest]
			vo['C_FelledMerch'][iT,iHarvest]=Affected_StemMerch[iHarvest]+Affected_Bark[iHarvest]+Affected_DeadStemMerch[iHarvest]
			vo['C_FelledRoots'][iT,iHarvest]=Affected_RootCoarse[iHarvest]+Affected_Bark[iHarvest]+Affected_RootFine[iHarvest]

		# Live merch. stemwood volume
		Affected_VolumeStemMerchLive=LiveImpactFactor*vo['V_MerchLive'][iT,:]
		Affected_VolumeStemWholeLive=LiveImpactFactor*vo['V_WholeStemLive'][iT,:]

		# Dead merch. stemwood volume
		Affected_VolumeStemMerchDead=DeadImpactFactor*vo['V_MerchDead'][iT,:]
		Affected_VolumeStemWholeDead=DeadImpactFactor*vo['V_WholeStemDead'][iT,:]

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
			vo['C_M_DistByAgent'][k][iT,ind]=vo['C_M_DistByAgent'][k][iT,ind]+Affected_BiomassTotal[ind]

		#----------------------------------------------------------------------
		# Remove affected amount from each pool
		#----------------------------------------------------------------------

		if meta[pNam]['Project']['Biomass Module']!='Sawtooth':

			# Remove carbon from affected biomass pools
			vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['StemMerch']]-Affected_StemMerch
			vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['StemNonMerch']]-Affected_StemNonMerch
			vo['C_Eco_ByPool'][iT,:,iEP['Bark']]=vo['C_Eco_ByPool'][iT,:,iEP['Bark']]-Affected_Bark
			vo['C_Eco_ByPool'][iT,:,iEP['Branch']]=vo['C_Eco_ByPool'][iT,:,iEP['Branch']]-Affected_Branch
			vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]=vo['C_Eco_ByPool'][iT,:,iEP['Foliage']]-Affected_Foliage
			vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]=vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]-Affected_RootCoarse
			vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]=vo['C_Eco_ByPool'][iT,:,iEP['RootFine']]-Affected_RootFine

			# Remove carbon from dead wood pools
			vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]-Affected_DeadStemMerch
			vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]-Affected_DeadStemNonMerch
			vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]-Affected_DeadBark
			vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]-Affected_DeadBranch

			# Remove stemwood volume
			vo['V_MerchLive'][iT,:]=np.maximum(0,vo['V_MerchLive'][iT,:]-Affected_VolumeStemMerchLive)
			vo['V_WholeStemLive'][iT,:]=np.maximum(0,vo['V_WholeStemLive'][iT,:]-Affected_VolumeStemWholeLive)
			vo['V_MerchDead'][iT,:]=np.maximum(0,vo['V_MerchDead'][iT,:]-Affected_VolumeStemMerchDead)
			vo['V_WholeStemDead'][iT,:]=np.maximum(0,vo['V_WholeStemDead'][iT,:]-Affected_VolumeStemWholeDead)

		#----------------------------------------------------------------------
		# Carbon that is removed (ie sent to mills)
		#----------------------------------------------------------------------

		# Green material to mill - of the total amount of biomass affected,
		vo['C_ToMillMerchGreen'][iT,:]=vo['C_ToMillMerchGreen'][iT,:]+b['GreenStemMerchRemoved']*Affected_StemMerch+ \
			b['BarkRemoved']*(0.9*Affected_Bark)
		vo['C_ToMillNonMerchGreen'][iT,:]=vo['C_ToMillNonMerchGreen'][iT,:]+b['GreenStemNonMerchRemoved']*Affected_StemNonMerch+ \
			b['BarkRemoved']*(0.1*Affected_Bark)

		# Dead wood to mill
		vo['C_ToMillMerchDead'][iT,:]=vo['C_ToMillMerchDead'][iT,:]+b['DeadStemMerchRemoved']*Affected_DeadStemMerch+ \
			+b['DeadBarkRemoved']*(0.9*Affected_DeadBark)
		vo['C_ToMillNonMerchDead'][iT,:]=vo['C_ToMillNonMerchDead'][iT,:]+b['DeadStemNonMerchRemoved']*Affected_DeadStemNonMerch+ \
			b['DeadBarkRemoved']*(0.1*Affected_DeadBark)

		#----------------------------------------------------------------------
		# Volume that is removed (sent to mills)
		#----------------------------------------------------------------------

		if meta[pNam]['Project']['Biomass Module']!='Sawtooth':
			vo['V_ToMill_MerchGreen'][iT,:]=vo['V_ToMill_MerchGreen'][iT,:]+b['GreenStemMerchRemoved']*Affected_VolumeStemMerchLive
			vo['V_ToMill_MerchDead'][iT,:]=vo['V_ToMill_MerchDead'][iT,:]+b['DeadStemMerchRemoved']*Affected_VolumeStemMerchDead
		else:
			vo['V_ToMill_MerchGreen'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_ToMillMerchGreen'][iT,:]
			vo['V_ToMill_MerchDead'][iT,:]=vi['lsat']['Biomass to Volume CF']*vo['C_ToMillMerchDead'][iT,:]

		# Total merch stemwood volume removed
		vo['V_ToMill_MerchTotal'][iT,:]=vo['V_ToMill_MerchGreen'][iT,:]+vo['V_ToMill_MerchDead'][iT,:]

		# Non-mech volume removed
		vo['V_ToMill_NonMerchTotal'][iT,:]=vi['lsat']['Biomass to Volume CF']*(vo['C_ToMillNonMerchGreen'][iT,:]+vo['C_ToMillNonMerchDead'][iT,:])

		#----------------------------------------------------------------------
		# Carbon that is moved from biomass to dead wood
		#----------------------------------------------------------------------

		StemMerch_M=b['GreenStemMerchToDeadStemMerch']*Affected_StemMerch
		StemNonMerch_M=b['GreenStemNonMerchToDeadStemNonMerch']*Affected_StemNonMerch
		Bark_M=b['BarkToDeadBark']*Affected_Bark
		Branch_M=b['BranchToDeadBranch']*Affected_Branch
		C_ToDeadWood=StemMerch_M+StemNonMerch_M+Bark_M+Branch_M

		vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemMerch']]+StemMerch_M
		vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadStemNonMerch']]+StemNonMerch_M
		vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBark']]+Bark_M
		vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['DeadBranch']]+Branch_M

		#vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+StemMerch_M+StemNonMerch_M+Bark_M+Branch_M

		#----------------------------------------------------------------------
		# Carbon that is left dispersed on site (after felling or wind storms)
		# and root carbon that is transferred to DOM pools
		#----------------------------------------------------------------------

		# Biomasss to DOM

		c=b['GreenStemMerchLeftOnSite']*Affected_StemMerch
		C_BiomassToDOM=c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+c

		c=b['GreenStemNonMerchLeftOnSite']*Affected_StemNonMerch
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+c

		c=b['BarkLeftOnSite']*Affected_Bark
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+c

		c=b['BranchLeftOnSite']*Affected_Branch
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+c

		c=b['FoliageLeftOnSite']*Affected_Foliage
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]+c

		c=0.5*Affected_RootFine
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterVF']]+c

		c=0.5*Affected_RootFine
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['SoilVF']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilVF']]+c

		c=0.5*Affected_RootCoarse
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+c

		c=0.5*Affected_RootCoarse
		C_BiomassToDOM=C_BiomassToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['SoilF']]=vo['C_Eco_ByPool'][iT,:,iEP['SoilF']]+c

		# Dead Wood to DOM

		c=b['DeadStemMerchLeftOnSite']*Affected_DeadStemMerch
		C_DeadWoodToDOM=c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+c

		c=b['DeadStemNonMerchLeftOnSite']*Affected_DeadStemNonMerch
		C_DeadWoodToDOM=C_DeadWoodToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterM']]+c

		c=b['DeadBranchLeftOnSite']*Affected_DeadBranch
		C_DeadWoodToDOM=C_DeadWoodToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+c

		c=b['DeadBarkLeftOnSite']*Affected_DeadBark
		C_DeadWoodToDOM=C_DeadWoodToDOM+c
		vo['C_ToDOM'][iT,:]=vo['C_ToDOM'][iT,:]+c
		vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,:,iEP['LitterF']]+c

		# Track fate of felled to DOM
		if meta[pNam]['Project']['Track fate of felled to DOM Status']=='On':
			if meta['LUT']['Event']['Harvest'] in idx_Type:
				mnam='Fate of Roots and Dispersed Slash to DOM'
				c=b['GreenStemMerchLeftOnSite']*Affected_StemMerch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]+c

				c=b['GreenStemNonMerchLeftOnSite']*Affected_StemNonMerch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]+c
		
				c=b['BarkLeftOnSite']*Affected_Bark
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]+c

				c=b['BranchLeftOnSite']*Affected_Branch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]+c

				c=b['FoliageLeftOnSite']*Affected_Foliage
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterVF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterVF']]+c

				c=b['DeadStemMerchLeftOnSite']*Affected_DeadStemMerch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]+c

				c=b['DeadStemNonMerchLeftOnSite']*Affected_DeadStemNonMerch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterM']]+c

				c=b['DeadBranchLeftOnSite']*Affected_DeadBranch
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]+c

				c=b['DeadBarkLeftOnSite']*Affected_DeadBark
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]+c
		
				c=0.5*Affected_RootFine
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterVF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterVF']]+c

				c=0.5*Affected_RootFine
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['SoilVF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['SoilVF']]+c
		
				c=0.5*Affected_RootCoarse
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['LitterF']]+c

				c=0.5*Affected_RootCoarse
				meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['SoilF']]=meta['Modules'][mnam]['C_Eco_ByPool'][iT,iHarvest,iEP['SoilF']]+c

		#----------------------------------------------------------------------
		# Carbon that is piled
		#----------------------------------------------------------------------

		#PiledStemMerch=(1.0-meta['Param']['BE']['ProductDisposal']['PiledStemwoodToFirewoodDom'])*b['GreenStemMerchPiled']*Affected_StemMerch
		PiledStemMerch=b['GreenStemMerchPiled']*Affected_StemMerch
		PiledStemNonMerch=b['GreenStemNonMerchPiled']*Affected_StemNonMerch
		PiledBranch=b['BranchPiled']*Affected_Branch
		PiledBark=b['BarkPiled']*Affected_Bark
		PiledFoliage=b['FoliagePiled']*Affected_Foliage
		PiledDeadStemMerch=b['DeadStemMerchPiled']*Affected_DeadStemMerch
		PiledDeadStemNonMerch=b['DeadStemNonMerchPiled']*Affected_DeadStemNonMerch
		PiledDeadBark=b['DeadBarkPiled']*Affected_DeadBark
		PiledDeadBranch=b['DeadBranchPiled']*Affected_DeadBranch

		C_BiomassToPile=PiledStemMerch+PiledStemNonMerch+PiledBark+PiledBranch+PiledFoliage
		C_DeadWoodToPile=PiledDeadStemMerch+PiledDeadStemNonMerch+PiledDeadBark+PiledDeadBranch
		vo['C_ToPile'][iT,:]=vo['C_ToPile'][iT,:]+C_BiomassToPile+C_DeadWoodToPile
		vo['C_ToPileMerch'][iT,:]=vo['C_ToPileMerch'][iT,:]+PiledStemMerch

		# Add piled carbon to pile pools
		vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]+PiledStemMerch
		vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]+PiledStemNonMerch
		vo['C_Eco_ByPool'][iT,:,iEP['PiledFoliage']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledFoliage']]+PiledFoliage
		vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]+PiledBranch
		vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]+PiledBark
		vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]+PiledDeadStemMerch+PiledDeadStemNonMerch+PiledDeadBark
		vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]+PiledDeadBranch

		# A small fraction of piled wood is collected for firewood
		vo['C_ToBBP_FirewoodDom'][iT,:]=vo['C_ToBBP_FirewoodDom'][iT,:]+meta['Param']['BE']['ProductDisposal']['PiledStemwoodToFirewoodDom']*(PiledStemMerch+PiledDeadStemMerch+PiledDeadStemNonMerch)

		# Piles that are burned
		StemMerch_Burned=b['PiledStemMerchBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]
		StemNonMerch_Burned=b['PiledStemNonMerchBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]
		Branch_Burned=b['PiledBranchBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]
		Bark_Burned=b['PiledBarkBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]
		Foliage_Burned=b['PiledFoliageBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledFoliage']]
		DeadStem_Burned=b['PiledDeadStemBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]
		DeadBranch_Burned=b['PiledDeadBranchBurned']*vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]

		Total_Burned=StemMerch_Burned+StemNonMerch_Burned+Branch_Burned+Bark_Burned+ \
			Foliage_Burned+DeadStem_Burned+DeadBranch_Burned

		NonMerch_Burned=StemNonMerch_Burned+Bark_Burned+Branch_Burned+Foliage_Burned+ \
			DeadStem_Burned+DeadBranch_Burned

		DeadBurned=DeadStem_Burned+DeadBranch_Burned

		# Remove burned carbon from piles
		vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledStemMerch']]-StemMerch_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledStemNonMerch']]-StemNonMerch_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledBranch']]-Branch_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledBark']]-Bark_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledFoliage']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledFoliage']]-Foliage_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadStem']]-DeadStem_Burned
		vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]=vo['C_Eco_ByPool'][iT,:,iEP['PiledDeadBranch']]-DeadBranch_Burned

		# Add burned pile carbon to fire emissions
		vo['C_E_OpenBurningAsCO2'][iT,:]=vo['C_E_OpenBurningAsCO2'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*Total_Burned
		vo['C_E_OpenBurningAsCH4'][iT,:]=vo['C_E_OpenBurningAsCH4'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*Total_Burned
		vo['C_E_OpenBurningAsCO'][iT,:]=vo['C_E_OpenBurningAsCO'][iT,:]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*Total_Burned

		# If it is pile burning, record it
		vo['C_ToPileBurnTot'][iT,:]=vo['C_ToPileBurnTot'][iT,:]+Total_Burned
		vo['C_ToPileBurnMerch'][iT,:]=vo['C_ToPileBurnMerch'][iT,:]+Total_Burned-NonMerch_Burned

		vo['C_ToFire'][iT,:]=vo['C_ToFire'][iT,:]+Total_Burned

		#----------------------------------------------------------------------
		# Biomass and DOM burned in wildfire
		#----------------------------------------------------------------------

		if meta['LUT']['Event']['Wildfire'] in idx_Type:

			# Biomass burned
			BurnedStemwoodMerch=b['GreenStemMerchBurned'][iWF]*Affected_StemMerch[iWF]
			BurnedStemwoodNonMerch=b['GreenStemNonMerchBurned'][iWF]*Affected_StemNonMerch[iWF]
			BurnedBark=b['BarkBurned'][iWF]*Affected_Bark[iWF]
			BurnedBranch=b['BranchBurned'][iWF]*Affected_Branch[iWF]
			BurnedFoliage=b['FoliageBurned'][iWF]*Affected_Foliage[iWF]
	
			# Dead wood burned
			# "Affected" dead wood has already accounted for the fraction that is unaffected
			# so there is no fraction applied
			BurnedDeadStemMerch=Affected_DeadStemMerch[iWF]
			BurnedDeadStemNonMerch=Affected_DeadStemNonMerch[iWF]
			BurnedDeadBark=Affected_DeadBark[iWF]
			BurnedDeadBranch=Affected_DeadBranch[iWF]
	
			# Litter and soil burned
			BurnedLitterVF=DeadImpactFactor[iWF]*b['LitterVFBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['LitterVF']]
			BurnedLitterF=DeadImpactFactor[iWF]*b['LitterFBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['LitterF']]
			BurnedLitterM=DeadImpactFactor[iWF]*b['LitterMBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['LitterM']]
			BurnedLitterS=DeadImpactFactor[iWF]*b['LitterSBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['LitterS']]
			BurnedSoilVF=DeadImpactFactor[iWF]*b['SoilVFBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['SoilVF']]
			BurnedSoilF=DeadImpactFactor[iWF]*b['SoilFBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['SoilF']]
			BurnedSoilS=DeadImpactFactor[iWF]*b['SoilSBurned'][iWF]*vo['C_Eco_ByPool'][iT,iWF,iEP['SoilS']]
	
			# Remove burned carbon from DOM pools - burned biomass and dead wood has
			# already been removed from biomass and dead wood
			vo['C_Eco_ByPool'][iT,iWF,iEP['LitterVF']]=vo['C_Eco_ByPool'][iT,iWF,iEP['LitterVF']]-BurnedLitterVF
			vo['C_Eco_ByPool'][iT,iWF,iEP['LitterF']]=vo['C_Eco_ByPool'][iT,iWF,iEP['LitterF']]-BurnedLitterF
			vo['C_Eco_ByPool'][iT,iWF,iEP['LitterM']]=vo['C_Eco_ByPool'][iT,iWF,iEP['LitterM']]-BurnedLitterM
			vo['C_Eco_ByPool'][iT,iWF,iEP['LitterS']]=vo['C_Eco_ByPool'][iT,iWF,iEP['LitterS']]-BurnedLitterS
			vo['C_Eco_ByPool'][iT,iWF,iEP['SoilVF']]=vo['C_Eco_ByPool'][iT,iWF,iEP['SoilVF']]-BurnedSoilVF
			vo['C_Eco_ByPool'][iT,iWF,iEP['SoilF']]=vo['C_Eco_ByPool'][iT,iWF,iEP['SoilF']]-BurnedSoilF
			vo['C_Eco_ByPool'][iT,iWF,iEP['SoilS']]=vo['C_Eco_ByPool'][iT,iWF,iEP['SoilS']]-BurnedSoilS
	
			BurnedBiomass=BurnedStemwoodMerch+BurnedStemwoodNonMerch+BurnedBark+BurnedBranch+BurnedFoliage
			BurnedDeadWood=BurnedDeadStemMerch+BurnedDeadStemNonMerch+BurnedDeadBark+BurnedDeadBranch
			BurnedLitter=BurnedLitterVF+BurnedLitterF+BurnedLitterM+BurnedLitterS
			BurnedSoil=BurnedSoilVF+BurnedSoilF+BurnedSoilS
			BurnedTotal=BurnedBiomass+BurnedDeadWood+BurnedLitter+BurnedSoil
	
			# Used for Conservation of Mass Test
			if (vi['tv'][iT]==2023) & (iE>=0) & (iScn==3):
				print(BurnedTotal[0])
			vo['C_ToFire'][iT,iWF]=vo['C_ToFire'][iT,iWF]+BurnedTotal
	
			vo['C_E_WildfireAsCO2'][iT,iWF]=vo['C_E_WildfireAsCO2'][iT,iWF]+meta['Param']['BEV']['Biophysical']['CombFrac_CO2']*BurnedTotal
			vo['C_E_WildfireAsCH4'][iT,iWF]=vo['C_E_WildfireAsCH4'][iT,iWF]+meta['Param']['BEV']['Biophysical']['CombFrac_CH4']*BurnedTotal
			vo['C_E_WildfireAsCO'][iT,iWF]=vo['C_E_WildfireAsCO'][iT,iWF]+meta['Param']['BEV']['Biophysical']['CombFrac_CO']*BurnedTotal

		#----------------------------------------------------------------------
		# Update stand age
		#----------------------------------------------------------------------

		if meta[pNam]['Project']['Biomass Module']!='Sawtooth':
			if (meta[pNam]['Project']['Partial Mortality Affects Age']=='On'):
				vo['A'][iT,:]=AgeCorrection*vo['A'][iT,:]
			else:
				# Ensure stand-replacing events reset stand age to 0
				vo['A'][iT,(LiveImpactFactor==1)]=0

		# Ensure planting resets to age 0 (except in fill-planting)
		if meta['LUT']['Event']['Planting'] in idx_Type.keys():
			iPL=idx_Type[meta['LUT']['Event']['Planting']]

			# Don't reset age of fill-planting events - the age is already reset by the inciting disturbance
			iPL=iPL[ vi['EC']['ASET'][iT,iPL,iE]!=meta['LUT']['Derived']['ASET']['Fill Planting'] ]

			vo['A'][iT,iPL]=0

		#----------------------------------------------------------------------
		# Transition to new growth curve
		#----------------------------------------------------------------------

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
		# Impose effects on growth through growth curves
		#----------------------------------------------------------------------
		if meta[pNam]['Project']['Biomass Module']=='BatchTIPSY':

			if meta['LUT']['Event']['Regen Failure'] in idx_Type.keys():
				iFailure=idx_Type[meta['LUT']['Event']['Regen Failure']]
				vi['GC']['Active'][:,iFailure,:]=0
				vi['GC'][1][:,iFailure,:]=0

			if meta[pNam]['Scenario'][iScn]['NOSE Sim Status']=='On':
				if meta['LUT']['Event']['Regen at 25% Growth'] in idx_Type.keys():
					iPoorGrowth=idx_Type[meta['LUT']['Event']['Regen at 25% Growth']]
					vi['GC']['Active'][:,iPoorGrowth,:]=0.25*vi['GC']['Active'][:,iPoorGrowth,:]
					vi['GC'][1][:,iPoorGrowth,:]=0.25*vi['GC'][1][:,iPoorGrowth,:]

		#----------------------------------------------------------------------
		# Apply prescribed growth factors (in response to non-lethal events)
		# Growth factor should comes in DMEC as percent effect
		# +10 = 10% increase, -10 = 10% decrease
		#----------------------------------------------------------------------

		if (meta[pNam]['Project']['Biomass Module']!='Sawtooth') & (meta[pNam]['Project']['Biomass Module']!='gromo'):

			indAdj=np.where( (GrowthFactor!=9999) & (GrowthFactor!=0) )[0]

			if indAdj.size>0:

				Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

				# Convert growth factor to respoonse ratio
				GrowthResponseRatio=GrowthFactor.astype(float)/100+1.0

				NetGrowth=vi['GC']['Active'].copy()

				for iAdj in indAdj:
					if ID_Type[iAdj]==meta['LUT']['Event']['Western Spruce Budworm']:
						# Only a temporary change in growth (see severity class table)
						ResponsePeriod=6
						iResponse=np.where( (Age>=vo['A'][iT,iAdj]) & (Age<=vo['A'][iT,iAdj]+ResponsePeriod) )[0]
						NetGrowth[iResponse,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[iResponse,iAdj,:]
					else:
						# A permanent change in growth
						NetGrowth[:,iAdj,:]=GrowthResponseRatio[iAdj]*NetGrowth[:,iAdj,:]

				vi['GC']['Active']=NetGrowth

		#----------------------------------------------------------------------
		# Apply growth factors (in response to non-lethal events)
		# +10 = 10% increase, -10 = 10% decrease
		#----------------------------------------------------------------------

		if (meta[pNam]['Project']['Biomass Module']!='Sawtooth') & (meta[pNam]['Project']['Biomass Module']!='gromo'):

			Age=np.arange(0,meta['Modules']['GYM']['BatchTIPSY Maximum Age']+1,1)

			# Wildfire
			if meta['LUT']['Event']['Wildfire'] in idx_Type.keys():
				indA=idx_Type[meta['LUT']['Event']['Wildfire']]
				if indA.size>0:
					NetGrowth=vi['GC']['Active'].copy()
					ResponsePeriod=25
					ResponseRatio=0.25
					Shape=3
					for iA in indA:
						iResponse=np.where( (Age>=vo['A'][iT,iA]) & (Age<=vo['A'][iT,iA]+ResponsePeriod) )[0]
						n=iResponse.size
						x=np.arange(1,n+1)
						fD=ResponseRatio+(1-ResponseRatio)*(1-np.exp(-Shape*(x/ResponsePeriod)))
						fD=np.tile(np.reshape(fD,(-1,1)),(1,6))
						NetGrowth[iResponse,iA,:]=fD*NetGrowth[iResponse,iA,:]

				flg=0
				if flg==1:
					ResponsePeriod=25
					n=5
					GRR=0.25
					x=np.arange(0,n)
					plt.close('all'); plt.plot(x,GRR+(1-GRR)*(1-np.exp(-3*(x/ResponsePeriod))),'.b-')

				vi['GC']['Active']=NetGrowth

			# Mountain pine beetle
			if meta['LUT']['Event']['Mountain Pine Beetle'] in idx_Type.keys():
				indA=idx_Type[meta['LUT']['Event']['Mountain Pine Beetle']]
				if indA.size>0:
					NetGrowth=vi['GC']['Active'].copy()
					ResponsePeriod=20
					ResponseRatio=0.3
					Shape=3
					for iA in indA:
						iResponse=np.where( (Age>=vo['A'][iT,iA]) & (Age<=vo['A'][iT,iA]+ResponsePeriod) )[0]
						n=iResponse.size
						x=np.arange(1,n+1)
						fD=ResponseRatio+(1-ResponseRatio)*(1-np.exp(-Shape*(x/ResponsePeriod)))
						fD=np.tile(np.reshape(fD,(-1,1)),(1,6))
						NetGrowth[iResponse,iA,:]=fD*NetGrowth[iResponse,iA,:]
	
				flg=0
				if flg==1:
					ResponsePeriod=25
					n=5
					GRR=0.1
					x=np.arange(0,n)
					plt.close('all'); plt.plot(x,GRR+(1-GRR)*(1-np.exp(-3*(x/ResponsePeriod))),'.b-')
	
				vi['GC']['Active']=NetGrowth

			# Disease
			if meta['LUT']['Event']['Disease Root'] in idx_Type.keys():
				indA=idx_Type[meta['LUT']['Event']['Disease Root']]
				if indA.size>0:
					NetGrowth=vi['GC']['Active'].copy()
					ResponsePeriod=25
					ResponseRatio=0.75
					Shape=3
					for iA in indA:
						iResponse=np.where( (Age>=vo['A'][iT,iA]) & (Age<=vo['A'][iT,iA]+ResponsePeriod) )[0]
						n=iResponse.size
						x=np.arange(1,n+1)
						fD=ResponseRatio+(1-ResponseRatio)*(1-np.exp(-Shape*(x/ResponsePeriod)))
						fD=np.tile(np.reshape(fD,(-1,1)),(1,6))
						NetGrowth[iResponse,iA,:]=fD*NetGrowth[iResponse,iA,:]
				vi['GC']['Active']=NetGrowth

			# Wind
			if meta['LUT']['Event']['Wind'] in idx_Type.keys():
				indA=idx_Type[meta['LUT']['Event']['Wind']]
				if indA.size>0:
					NetGrowth=vi['GC']['Active'].copy()
					ResponsePeriod=5
					ResponseRatio=0.75
					Shape=3
					for iA in indA:
						iResponse=np.where( (Age>=vo['A'][iT,iA]) & (Age<=vo['A'][iT,iA]+ResponsePeriod) )[0]
						n=iResponse.size
						x=np.arange(1,n+1)
						fD=ResponseRatio+(1-ResponseRatio)*(1-np.exp(-Shape*(x/ResponsePeriod)))
						fD=np.tile(np.reshape(fD,(-1,1)),(1,6))
						NetGrowth[iResponse,iA,:]=fD*NetGrowth[iResponse,iA,:]
				vi['GC']['Active']=NetGrowth

	#--------------------------------------------------------------------------
	# Relative gravimetric mortality (%)
	#--------------------------------------------------------------------------

	for k in meta['LUT']['Event'].keys():
		id=meta['LUT']['Event'][k]
		vo['C_M_DistByAgentPct'][id][iT,:]=vo['C_M_DistByAgent'][id][iT,:]/C_Biomass_t0*100

	#--------------------------------------------------------------------------
	# Aerial nutrient application events
	#--------------------------------------------------------------------------

	if meta[pNam]['Project']['Biomass Module']!='gromo':

		# Generate index to stands that were fertilized
		meta['Modules']['NutrientApp']['iApplication']=np.where(flag_nutrient_app==1)[0]
	
		if meta['Modules']['NutrientApp']['iApplication'].size>0:
	
			# Adjust net growth of aboveground biomass
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'AbovegroundNetGrowth')
	
			# Adjust emissions
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,'Emissions')

	#--------------------------------------------------------------------------
	# Update time since wildfire counter
	#--------------------------------------------------------------------------
	mnam='Disturbance Effects on Wildfire Occurrence'
	# Update counter
	meta['Modules'][mnam]['Wildfire Counter']=meta['Modules'][mnam]['Wildfire Flag']*(meta['Modules'][mnam]['Wildfire Counter']+1)

	# Stop counter once it reaches recovery time
	iStop=meta['Modules'][mnam]['Wildfire Counter']>=meta['Param']['BE']['WildfireDisturbanceEffects']['Wildfire Response Recovery Time']
	meta['Modules'][mnam]['Wildfire Flag'][iStop]=0
	meta['Modules'][mnam]['Wildfire Counter'][iStop]=0

	return vo,vi

#%%
def ProductDynamics(meta,pNam,iT,iBat,vi,vo):

	#--------------------------------------------------------------------------
	# Some projects may want to switch Removed Fate Regime half way through a simulation
	# E.g. CT regime for a commercial thinning event in 1952, followed by regional
	# default in 1998.
	#--------------------------------------------------------------------------
	if 'Switch Felled Fate Regime' in meta[pNam]['Scenario'][meta[pNam]['iScn']].keys():
		if vi['tv'][iT]==meta[pNam]['Scenario'][meta[pNam]['iScn']]['Switch Felled Fate Regime Year']:
			Scenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Felled Fate Change Scenario']
			HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Switch Felled Fate Regime']
			if HistoricalRegime=='Regional Defaults':
				HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Region Code']
			for k in meta['Param']['BE']['FelledFate'][Scenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['FelledFate'][Scenario][HistoricalRegime][k]
				x=np.reshape(x,(-1,1))
				x=np.tile(x,(1,meta[pNam]['Project']['Batch Size'][iBat]))
				meta['Param']['BEV']['FelledFate'][k]=x

	if 'Switch Removed Fate Regime' in meta[pNam]['Scenario'][meta[pNam]['iScn']].keys():
		if vi['tv'][iT]==meta[pNam]['Scenario'][meta[pNam]['iScn']]['Switch Removed Fate Regime Year']:
			Scenario=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Removed Fate Change Scenario']
			HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Switch Removed Fate Regime']
			if HistoricalRegime=='Regional Defaults':
				HistoricalRegime=meta[pNam]['Scenario'][meta[pNam]['iScn']]['Region Code']
			for k in meta['Param']['BE']['RemovedFate'][Scenario][HistoricalRegime].keys():
				x=meta['Param']['BE']['RemovedFate'][Scenario][HistoricalRegime][k]
				x=np.reshape(x,(-1,1))
				x=np.tile(x,(1,meta[pNam]['Project']['Batch Size'][iBat]))
				meta['Param']['BEV']['RemovedFate'][k]=x

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

	# Biophysical parameters
	bB=meta['Param']['BEV']['Biophysical']

	#--------------------------------------------------------------------------
	# Update product pools
	#--------------------------------------------------------------------------

	vo['C_Pro_ByPool'][iT,:,:]=vo['C_Pro_ByPool'][iT-1,:,:]

	#--------------------------------------------------------------------------
	# Carbon transferred from forest to mills or direct to end-uses
	#--------------------------------------------------------------------------

	fx={}

	fx['C_ToPulpMill']=bRF['RemovedMerchToPulpMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToPulpMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToPulpMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToPulpMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToPelletMill']=bRF['RemovedMerchToPelletMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToPelletMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToPelletMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToPelletMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToLumberMill']=bRF['RemovedMerchToLumberMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToLumberMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToLumberMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToLumberMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToPlywoodMill']=bRF['RemovedMerchToPlywoodMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToPlywoodMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToPlywoodMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToPlywoodMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToOSBMill']=bRF['RemovedMerchToOSBMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToOSBMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToOSBMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToOSBMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToMDFMill']=bRF['RemovedMerchToMDFMill']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToMDFMill']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToMDFMill']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToMDFMill']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToLogExport']=bRF['RemovedMerchToLogExport']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToLogExport']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToLogExport']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToLogExport']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToBBP_PowerGrid']=bRF['RemovedMerchToIPP']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToIPP']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToIPP']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToIPP']*vo['C_ToMillNonMerchDead'][iT,:]

	fx['C_ToBBP_FirewoodCollection']=bRF['RemovedMerchToFirewood']*vo['C_ToMillMerchGreen'][iT,:] + \
				 bRF['RemovedNonMerchToFirewood']*vo['C_ToMillNonMerchGreen'][iT,:] + \
				 bRF['RemovedDeadStemToFirewood']*vo['C_ToMillMerchDead'][iT,:] + \
				 bRF['RemovedDeadStemToFirewood']*vo['C_ToMillNonMerchDead'][iT,:]

	#--------------------------------------------------------------------------
	# Carbon transferred from mill to mill
	#--------------------------------------------------------------------------

	C_LumberToPulp=bPT['LumberMillToPulpMill']*fx['C_ToLumberMill']
	C_LumberToMDF=bPT['LumberMillToMDFMill']*fx['C_ToLumberMill']
	C_LumberToPellet=bPT['LumberMillToPelletMill']*fx['C_ToLumberMill']
	C_LumberToExport=bPT['LumberMillToLogExport']*fx['C_ToLumberMill']

	#fx['C_ToLumberMill']=fx['C_ToLumberMill']-C_LumberToPulp-C_LumberToMDF-C_LumberToPellet-C_LumberToExport

	fx['C_ToPulpMill']=fx['C_ToPulpMill']+C_LumberToPulp
	fx['C_ToMDFMill']=fx['C_ToMDFMill']+C_LumberToMDF
	fx['C_ToPelletMill']=fx['C_ToPelletMill']+C_LumberToPellet
	fx['C_ToLogExport']=fx['C_ToLogExport']+C_LumberToExport

	#--------------------------------------------------------------------------
	# Carbon transferred from log Export to firewood
	#--------------------------------------------------------------------------

	C_ToBBP_FirewoodExport=bPT['LogExportToFirewood']*fx['C_ToLogExport']

	#--------------------------------------------------------------------------
	# Carbon transferred from pulp to paper and effluent
	#--------------------------------------------------------------------------

	C_ToPaper=bPT['PulpMillToPaper']*fx['C_ToPulpMill']
	C_ToPulpEffluent=bPT['PulpMillToEffluent']*fx['C_ToPulpMill']

	#--------------------------------------------------------------------------
	# Carbon transferred to pellets
	#--------------------------------------------------------------------------

	C_ToBBP_PelletExport=bPT['PelletMillToPelletExport']*fx['C_ToPelletMill']
	C_ToBBP_PelletDomGrid=bPT['PelletMillToDomGrid']*fx['C_ToPelletMill']
	C_ToBBP_PelletDomRNG=bPT['PelletMillToDomRNG']*fx['C_ToPelletMill']

	#--------------------------------------------------------------------------
	# Carbon transferred to facility power
	#--------------------------------------------------------------------------

	C_ToBBP_PowerFacilityDom=bPT['LumberMillToPowerFacility']*fx['C_ToLumberMill'] + \
		bPT['PulpMillToPowerFacility']*fx['C_ToPulpMill'] + \
		bPT['PlywoodMillToPowerFacility']*fx['C_ToPlywoodMill'] + \
		bPT['OSBMillToPowerFacility']*fx['C_ToOSBMill']

	C_ToBBP_PowerFacilityExport=bPT['LogExportToPowerFacility']*fx['C_ToLogExport']

	#--------------------------------------------------------------------------
	# Carbon transferred to grid by independent power producers
	#--------------------------------------------------------------------------

	fx['C_ToBBP_PowerGrid']=fx['C_ToBBP_PowerGrid']+bPT['LumberMillToIPP']*fx['C_ToLumberMill'] + \
		bPT['PulpMillToIPP']*fx['C_ToPulpMill'] + \
		bPT['PlywoodMillToIPP']*fx['C_ToPlywoodMill'] + \
		bPT['OSBMillToIPP']*fx['C_ToOSBMill']

	#--------------------------------------------------------------------------
	# Production of single-family homes
	#--------------------------------------------------------------------------

	C_LumberMillToSFH=fx['C_ToLumberMill']*bPT['LumberMillToSFH']
	C_PlywoodMillToSFH=fx['C_ToPlywoodMill']*bPT['PlywoodMillToSFH']
	C_OSBMillToSFH=fx['C_ToOSBMill']*bPT['OSBMillToSFH']
	C_MDFMillToSFH=fx['C_ToMDFMill']*bPT['MDFMillToSFH']
	C_LogExportToSFH=fx['C_ToLogExport']*bPT['LogExportToSFH']

	#--------------------------------------------------------------------------
	# Production of multi-family homes
	#--------------------------------------------------------------------------

	C_LumberMillToMFH=fx['C_ToLumberMill']*bPT['LumberMillToMFH']
	C_PlywoodMillToMFH=fx['C_ToPlywoodMill']*bPT['PlywoodMillToMFH']
	C_OSBMillToMFH=fx['C_ToOSBMill']*bPT['OSBMillToMFH']
	C_MDFMillToMFH=fx['C_ToMDFMill']*bPT['MDFMillToMFH']
	C_LogExportToMFH=fx['C_ToLogExport']*bPT['LogExportToMFH']

	#--------------------------------------------------------------------------
	# Production of commercial buildings
	#--------------------------------------------------------------------------

	C_LumberMillToCom=fx['C_ToLumberMill']*bPT['LumberMillToCom']
	C_PlywoodMillToCom=fx['C_ToPlywoodMill']*bPT['PlywoodMillToCom']
	C_OSBMillToCom=fx['C_ToOSBMill']*bPT['OSBMillToCom']
	C_MDFMillToCom=fx['C_ToMDFMill']*bPT['MDFMillToCom']
	C_LogExportToCom=fx['C_ToLogExport']*bPT['LogExportToCom']

	#--------------------------------------------------------------------------
	# Production of furniture
	#--------------------------------------------------------------------------

	C_LumberMillToFurn=fx['C_ToLumberMill']*bPT['LumberMillToFurn']
	C_PlywoodMillToFurn=fx['C_ToPlywoodMill']*bPT['PlywoodMillToFurn']
	C_OSBMillToFurn=fx['C_ToOSBMill']*bPT['OSBMillToFurn']
	C_MDFMillToFurn=fx['C_ToMDFMill']*bPT['MDFMillToFurn']
	C_LogExportToFurn=fx['C_ToLogExport']*bPT['LogExportToFurn']

	#--------------------------------------------------------------------------
	# Production of shipping containers
	#--------------------------------------------------------------------------

	C_LumberMillToShip=fx['C_ToLumberMill']*bPT['LumberMillToShip']
	C_PlywoodMillToShip=fx['C_ToPlywoodMill']*bPT['PlywoodMillToShip']
	C_OSBMillToShip=fx['C_ToOSBMill']*bPT['OSBMillToShip']
	C_MDFMillToShip=fx['C_ToMDFMill']*bPT['MDFMillToShip']
	C_LogExportToShip=fx['C_ToLogExport']*bPT['LogExportToShip']

	#--------------------------------------------------------------------------
	# Production of repairs
	#--------------------------------------------------------------------------

	C_LumberMillToRepairs=fx['C_ToLumberMill']*bPT['LumberMillToRepairs']
	C_PlywoodMillToRepairs=fx['C_ToPlywoodMill']*bPT['PlywoodMillToRepairs']
	C_OSBMillToRepairs=fx['C_ToOSBMill']*bPT['OSBMillToRepairs']
	C_MDFMillToRepairs=fx['C_ToMDFMill']*bPT['MDFMillToRepairs']
	C_LogExportToRepairs=fx['C_ToLogExport']*bPT['LogExportToRepairs']

	#--------------------------------------------------------------------------
	# Production of other
	#--------------------------------------------------------------------------

	C_LumberMillToOther=fx['C_ToLumberMill']*bPT['LumberMillToOther']
	C_PlywoodMillToOther=fx['C_ToPlywoodMill']*bPT['PlywoodMillToOther']
	C_OSBMillToOther=fx['C_ToOSBMill']*bPT['OSBMillToOther']
	C_MDFMillToOther=fx['C_ToMDFMill']*bPT['MDFMillToOther']
	C_LogExportToOther=fx['C_ToLogExport']*bPT['LogExportToOther']

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
	vo['C_ToBBP_PowerGrid'][iT,:]=fx['C_ToBBP_PowerGrid']
	vo['C_ToBBP_PowerFacilityDom'][iT,:]=C_ToBBP_PowerFacilityDom
	vo['C_ToBBP_PowerFacilityExport'][iT,:]=C_ToBBP_PowerFacilityExport
	vo['C_ToBBP_PelletExport'][iT,:]=C_ToBBP_PelletExport
	vo['C_ToBBP_PelletDomGrid'][iT,:]=C_ToBBP_PelletDomGrid
	vo['C_ToBBP_PelletDomRNG'][iT,:]=C_ToBBP_PelletDomRNG
	vo['C_ToBBP_FirewoodDom'][iT,:]=vo['C_ToBBP_FirewoodDom'][iT,:]+fx['C_ToBBP_FirewoodCollection'] # *** There is firewood taken directly from forest ecosystems (see events function) ***
	vo['C_ToBBP_FirewoodExport'][iT,:]=C_ToBBP_FirewoodExport
	vo['C_ToLogExport'][iT,:]=fx['C_ToLogExport']

	#--------------------------------------------------------------------------
	# Carbon transferred from mills to in-use products
	#--------------------------------------------------------------------------

	# Transfer mill fibre to single-family homes
	iP=meta['Core']['iPP']['SFH']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToSFH + \
		C_PlywoodMillToSFH + \
		C_OSBMillToSFH + \
		C_MDFMillToSFH + \
		C_LogExportToSFH

	# Transfer mill fibre to multi-family homes
	iP=meta['Core']['iPP']['MFH']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToMFH + \
		C_PlywoodMillToMFH + \
		C_OSBMillToMFH + \
		C_MDFMillToMFH + \
		C_LogExportToMFH

	# Transfer mill fibre to commercial
	iP=meta['Core']['iPP']['Comm']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToCom + \
		C_PlywoodMillToCom + \
		C_OSBMillToCom + \
		C_MDFMillToCom + \
		C_LogExportToCom

	# Transfer mill fibre to furniture
	iP=meta['Core']['iPP']['Furn']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToFurn + \
		C_PlywoodMillToFurn + \
		C_OSBMillToFurn + \
		C_MDFMillToFurn + \
		C_LogExportToFurn

	# Transfer mill fibre to shipping
	iP=meta['Core']['iPP']['Ship']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToShip + \
		C_PlywoodMillToShip + \
		C_OSBMillToShip + \
		C_MDFMillToShip + \
		C_LogExportToShip

	# Transfer mill fibre to repairs
	iP=meta['Core']['iPP']['Repairs']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToRepairs + \
		C_PlywoodMillToRepairs + \
		C_OSBMillToRepairs + \
		C_MDFMillToRepairs + \
		C_LogExportToRepairs

	# Transfer mill fibre to other
	iP=meta['Core']['iPP']['Other']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + \
		C_LumberMillToOther + \
		C_PlywoodMillToOther + \
		C_OSBMillToOther + \
		C_MDFMillToOther + \
		C_LogExportToOther

	# Transfer pulp mill fibre to paper
	iP=meta['Core']['iPP']['Paper']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]+C_ToPaper

	# Transfer pulp mill carbon to pulp-mill effluent
	iP=meta['Core']['iPP']['EffluentPulp']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]+C_ToPulpEffluent

	#--------------------------------------------------------------------------
	# Single-family homes --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['SFH']
	C_retired=bPD['SFH_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['SFHToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['SFHToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['SFHToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Multi-family homes --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['MFH']
	C_retired=bPD['MFH_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['MFHToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['MFHToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['MFHToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Commercial building --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Comm']
	C_retired=bPD['Comm_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['CommToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['CommToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['CommToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Furniture --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Furn']
	C_retired=bPD['Furn_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['FurnToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['FurnToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['FurnToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Shipping --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Ship']
	C_retired=bPD['Ship_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['ShipToDumpWood']*C_retired

	# Transfer carbon to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['ShipToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['ShipToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Repairs --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Repairs']
	C_retired=bPD['Repairs_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['RepairsToDumpWood']*C_retired

	# Transfer carbon to landfill (degradble)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['RepairsToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['RepairsToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Other --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['Other']
	C_retired=bPD['Other_tr']*vo['C_Pro_ByPool'][iT-1,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer carbon to dump wood
	iP=meta['Core']['iPP']['DumpWood']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['OtherToDumpWood']*C_retired

	# Transfer carbon to landfill (degradble)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['OtherToLandfillWood']*bPD['ToLandfillWoodDegradableFrac']*C_retired

	# Transfer carbon to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['OtherToLandfillWood']*(1-bPD['ToLandfillWoodDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Paper --> dump and landfill
	#--------------------------------------------------------------------------

	# Turnover (with adjustment for recycling)
	iP=meta['Core']['iPP']['Paper']
	C_retired=(1-bPD['PaperRecycleRate'])*bPD['Paper_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Remove carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] - C_retired

	# Transfer to dump
	iP=meta['Core']['iPP']['DumpPaper']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['PaperToDumpPaper']*C_retired

	# Transfer to landfill (degradable)
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['PaperToLandfillPaper']*bPD['ToLandfillPaperDegradableFrac']*C_retired

	# Transfer to landfill (non-degradable)
	iP=meta['Core']['iPP']['LandfillWoodNonDegradable']
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP] + bPD['PaperToLandfillPaper']*(1-bPD['ToLandfillPaperDegradableFrac'])*C_retired

	#--------------------------------------------------------------------------
	# Emissions from combustion during grid power
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PowerGrid'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PowerGrid'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Domestic_EnergySC_Bioenergy_PowerGrid'][iT,:]=vo['E_Domestic_EnergySC_Bioenergy_PowerGrid'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PowerGrid'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion during domestic facility power generation
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PowerFacilityDom'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PowerFacilityDom'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Domestic_EnergySC_Bioenergy_PowerFacility'][iT,:]=vo['E_Domestic_EnergySC_Bioenergy_PowerFacility'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PowerFacilityDom'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion during international power generation
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PowerFacilityExport'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PowerFacilityExport'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Internat_EnergySC_Bioenergy_PowerFacility'][iT,:]=vo['E_Internat_EnergySC_Bioenergy_PowerFacility'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PowerFacilityExport'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet Export
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PelletExport'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PelletExport'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Internat_EnergySC_Bioenergy_PelletGrid'][iT,:]=vo['E_Internat_EnergySC_Bioenergy_PelletGrid'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PelletExport'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet domestic grid
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PelletDomGrid'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PelletDomGrid'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Domestic_EnergySC_Bioenergy_PelletGrid'][iT,:]=vo['E_Domestic_EnergySC_Bioenergy_PelletGrid'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PelletDomGrid'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of pellet domestic RNG
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_PelletDomRNG'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_PelletDomRNG'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Domestic_EnergySC_Bioenergy_PelletRNG'][iT,:]=vo['E_Domestic_EnergySC_Bioenergy_PelletRNG'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_PelletDomRNG'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of domestic firewood
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_FirewoodDom'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_FirewoodDom'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Domestic_EnergySC_Bioenergy_Firewood'][iT,:]=vo['E_Domestic_EnergySC_Bioenergy_Firewood'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_FirewoodDom'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from combustion of firewood that was exported
	#--------------------------------------------------------------------------

	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*vo['C_ToBBP_FirewoodExport'][iT,:]
	E_CH4=bB['Ratio_CH4_to_C']*(1-bPD['EnergyCombustionFracEmitCO2'])*vo['C_ToBBP_FirewoodExport'][iT,:]
	E_AsCH4=meta['Param']['BE']['Biophysical']['GWP_CH4_AR5']*E_CH4
	vo['E_Internat_EnergySC_Bioenergy_Firewood'][iT,:]=vo['E_Internat_EnergySC_Bioenergy_Firewood'][iT,:] + (E_AsCO2+E_AsCH4)
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+vo['C_ToBBP_FirewoodExport'][iT,:]

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Emissions from pulp effluent
	#--------------------------------------------------------------------------

	# Emissions from pulp effluent (CO2 from aerobic decomposition)
	iP=meta['Core']['iPP']['EffluentPulp']
	C_emitted=bPD['EffluentPulp_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Remove emitted carbon
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]-C_emitted

	# Add emitted carbon to CO2 emission "pool"
	# Emissions
	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['EnergyCombustionFracEmitCO2']*C_emitted
	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:] + E_AsCO2
	vo['C_BBP'][iT,:]=vo['C_BBP'][iT,:]+C_emitted

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of dump wood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['DumpWood']
	C_emitted=bPD['DumpWood_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Removal
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]-C_emitted

	E_AsCO2=bB['Ratio_CO2_to_C']*C_emitted

	# Add to emissions (CO2 emission from aerobic decomposition)
	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:] + E_AsCO2
	vo['C_RHP'][iT,:]=vo['C_RHP'][iT,:]+C_emitted

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of dump paper
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['DumpPaper']
	C_emitted=bPD['DumpPaper_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Removal
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]-C_emitted

	E_AsCO2=bB['Ratio_CO2_to_C']*C_emitted

	# Add to emissions (CO2 emission from aerobic decomposition)
	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:] + E_AsCO2
	vo['C_RHP'][iT,:]=vo['C_RHP'][iT,:]+C_emitted

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of landfill degradable wood
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['LandfillWoodDegradable']
	C_emitted=bPD['LandfillWoodDegradable_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Removal
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]-C_emitted

	# Add to emissions
	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['LandfillDegradableFracEmitCO2']*C_emitted

	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:]+E_AsCO2
	vo['C_RHP'][iT,:]=vo['C_RHP'][iT,:]+C_emitted

	# Adjustment for proportion of degradable landfills with gas collection systems,
	# efficiency of system, and methane oxided to CO2 from the landfill cover
	gcsp=bPD['LandfillMethaneEmit_GasColSysProp']
	c1=1-gcsp
	c2=1-bPD['LandfillMethaneEmit_GasColSysEffic']
	E_C_AsCH4=C_emitted*(c1-bPD['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-bPD['LandfillMethaneOxidizedToCO2']*c1)

	E_CH4=bB['Ratio_CH4_to_C']*E_C_AsCH4

	E_AsCH4=bB['GWP_CH4_AR5']*E_CH4

	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:]+E_AsCH4

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

	#--------------------------------------------------------------------------
	# Decomposition of landfill degradable paper
	#--------------------------------------------------------------------------

	# Turnover
	iP=meta['Core']['iPP']['LandfillPaperDegradable']
	C_emitted=bPD['LandfillWoodDegradable_tr']*vo['C_Pro_ByPool'][iT,:,iP]

	# Removal
	vo['C_Pro_ByPool'][iT,:,iP]=vo['C_Pro_ByPool'][iT,:,iP]-C_emitted

	# Add to emissions
	E_AsCO2=bB['Ratio_CO2_to_C']*bPD['LandfillDegradableFracEmitCO2']*C_emitted

	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:]+E_AsCO2
	vo['C_RHP'][iT,:]=vo['C_RHP'][iT,:]+C_emitted

	# Adjustment for proportion of degradable landfills with gas collection systems,
	# efficiency of system, and methane oxided to CO2 from the landfill cover
	E_C_AsCH4=C_emitted*(c1-bPD['LandfillMethaneOxidizedToCO2']*c1)+C_emitted*gcsp*(c2-bPD['LandfillMethaneOxidizedToCO2']*c1)
	E_AsCH4=bB['GWP_CH4_AR5']*E_C_AsCH4

	E_CH4=bB['Ratio_CH4_to_C']*E_C_AsCH4

	E_AsCH4=bB['GWP_CH4_AR5']*E_CH4

	vo['E_Domestic_ForestSector_HWP'][iT,:]=vo['E_Domestic_ForestSector_HWP'][iT,:]+E_AsCH4

	# Track if radiative forcing status is on
	if meta[pNam]['Project']['Radiative Forcing Status']=='On':
		vo['E_CO2'][iT,:]=vo['E_CO2'][iT,:]+E_AsCO2
		vo['E_CH4'][iT,:]=vo['E_CH4'][iT,:]+E_CH4

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
	V_Tot=vo['V_ToMill_MerchTotal']+vo['V_ToMill_NonMerchTotal']

	# Total removals (ODT/ha)
	ODT_Tot=V_Tot*bB['Density Wood']

	#--------------------------------------------------------------------------
	# Resource extraction
	#--------------------------------------------------------------------------

	# Construction and maintenance of roads
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Road Construction']*V_Tot

	# Cruising and reconnaissance
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Cruise And Recon']*V_Tot

	# Felling and processing logs
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Felling Process Logs']*V_Tot

	# Skidding trees to landing
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Skidding To Landing']*V_Tot

	# Piling and sorting
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Piling And Sorting Logs']*V_Tot

	# Loading logs
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Loading Logs At Landing']*V_Tot

	# Chipping (non-merch)
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Chipping']*vo['V_ToMill_NonMerchTotal']

	# Hauling (forest ecosystem to mill)
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+ODT_Tot/bB['Moisture Content Wood']*bB['Emission Intensity Transport Truck']*bB['Distance Forest To Mill (One Way)']

	# Site preparation
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Site Prep']*V_Tot

	# Sowing seeds
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Sowing']*V_Tot

	# Planting
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Planting']*V_Tot

	# Surveying
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Surveying']*V_Tot

	#--------------------------------------------------------------------------
	# Mill operations
	#--------------------------------------------------------------------------

	# Unloading
	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['EI Unloading At Mill']*V_Tot

	# Sawing and processing lumber
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Sawing Processing Lumber']*(vo['C_ToLumber']/bB['Carbon Content Wood'])

	# Sawing and processing plywood
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Sawing Processing Plywood']*(vo['C_ToPlywood']/bB['Carbon Content Wood'])

	# Sawing and processing OSB
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Sawing Processing OSB']*(vo['C_ToOSB']/bB['Carbon Content Wood'])

	# Sawing and processing MDF
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Sawing Processing MDF']*(vo['C_ToMDF']/bB['Carbon Content Wood'])

	# Sawing and processing pulp
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Processing Pulp']*(vo['C_ToPaper']/bB['Carbon Content Wood'])

	# Pellets
	C_Pellet=vo['C_ToBBP_PelletExport']+vo['C_ToBBP_PelletDomGrid']+vo['C_ToBBP_PelletDomRNG']

	# Size reduction of pellets
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Pellet Size Reduction']*C_Pellet/bB['Carbon Content Wood']

	# Drying pellets
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Pellet Drying']*C_Pellet/bB['Carbon Content Wood']

	# Pelletizing
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Pellet Pelletizing']*C_Pellet/bB['Carbon Content Wood']

	# Pellet sieving
	vo['E_Domestic_EnergySC_ForestOperationsBurnOil']=vo['E_Domestic_EnergySC_ForestOperationsBurnOil']+bB['EI Pellet Seiving']*C_Pellet/bB['Carbon Content Wood']

	#--------------------------------------------------------------------------
	# Transport Mill -> Distribution Hub (Vancouver)
	# *** Exclude raw log exports, as they originate at mills on the ocean ***
	#--------------------------------------------------------------------------

	Mass=(vo['C_ToLumber']+vo['C_ToPlywood']+vo['C_ToOSB']+vo['C_ToMDF']+vo['C_ToPaper']+vo['C_ToBBP_PelletExport']+vo['C_ToBBP_PelletDomRNG'])/bB['Carbon Content Wood']/bB['Moisture Content Lumber']

	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['Distance Mill To Distribution Hub']*bB['Emission Intensity Transport Rail']*Mass

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

	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail+E_HubToMarket_Truck

	# Log exports

	Mass=(vo['C_ToLogExport'])/bB['Carbon Content Wood']/bB['Moisture Content Wood']

	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+bB['Distance LogExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

	# Pellets

	Mass=vo['C_ToBBP_PelletExport']/bB['Carbon Content Wood']

	E_HubToMarket_Water=bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass+ \
		bB['Fraction PelletExport Water Dest 1']*bB['Distance PelletExport Water Dest 1']*bB['Emission Intensity Transport Water (Bulk)']*Mass

	E_HubToMarket_Rail=bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass+ \
		bB['Fraction PelletExport Rail Dest 1']*bB['Distance PelletExport Rail Dest 1']*bB['Emission Intensity Transport Rail']*Mass

	# E_HubToMarket_Truck=bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
	#	 bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass+ \
	#	 bB['Fraction PelletExport Truck Dest 1']*bB['Distance PelletExport Truck Dest 1']*bB['Emission Intensity Transport Truck']*Mass

	vo['E_Domestic_EnergyT_ForestOperationsBurnOil']=vo['E_Domestic_EnergyT_ForestOperationsBurnOil']+E_HubToMarket_Water+E_HubToMarket_Rail

	#==========================================================================
	# Substitution effects
	#==========================================================================

	#----------------------------------------------------------------------
	# Domestic facility power generation (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	Yield_PowerFacilityDom=vo['C_ToBBP_PowerFacilityDom']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerFacilityDom=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityDom

	GJ_PowerFacilityDom=GJ_PowerFacilityDom*Ratio_EC

	E_Domestic_Substitution_CoalForBioenergy_PowerFacility=bS['PowerFacilityDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityDom
	E_Domestic_Substitution_DieselForBioenergy_PowerFacility=bS['PowerFacilityDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityDom
	E_Domestic_Substitution_GasForBioenergy_PowerFacility=bS['PowerFacilityDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityDom
	E_Domestic_Substitution_OilForBioenergy_PowerFacility=bS['PowerFacilityDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityDom

	#----------------------------------------------------------------------
	# International facility power generation (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	# Yield (green tonnes/ha)
	Yield_PowerFacilityExport=vo['C_ToBBP_PowerFacilityExport']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerFacilityExport=bB['Energy Content Wood (0% moisture)']*Yield_PowerFacilityExport

	GJ_PowerFacilityExport=GJ_PowerFacilityExport*Ratio_EC

	E_Internat_Substitution_CoalForBioenergy_PowerFacility=bS['PowerFacilityExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerFacilityExport
	E_Internat_Substitution_DieselForBioenergy_PowerFacility=bS['PowerFacilityExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerFacilityExport
	E_Internat_Substitution_GasForBioenergy_PowerFacility=bS['PowerFacilityExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerFacilityExport
	E_Internat_Substitution_OilForBioenergy_PowerFacility=bS['PowerFacilityExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerFacilityExport

	#----------------------------------------------------------------------
	# Independent power producers (MgC/ha) to (green tonne/ha)
	#----------------------------------------------------------------------

	# Yield (green tonnes/ha)
	Yield_PowerGrid=vo['C_ToBBP_PowerGrid']/bB['Density Wood']/bB['Moisture Content Wood']

	# Yield to energy (GJ/ha)
	GJ_PowerGrid=bB['Energy Content Wood (0% moisture)']*Yield_PowerGrid

	GJ_PowerGrid=GJ_PowerGrid*Ratio_EC

	E_Domestic_Substitution_CoalForBioenergy_ElectricityGrid=bS['PowerGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PowerGrid
	E_Domestic_Substitution_DieselForBioenergy_ElectricityGrid=bS['PowerGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PowerGrid
	E_Domestic_Substitution_GasForBioenergy_ElectricityGrid=bS['PowerGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PowerGrid
	E_Domestic_Substitution_OilForBioenergy_ElectricityGrid=bS['PowerGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PowerGrid

	#----------------------------------------------------------------------
	# Pellet exports (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dried tonnes/ha)
	Yield_PelletExport=vo['C_ToBBP_PelletExport']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletExport=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletExport

	GJ_PelletExport=GJ_PelletExport*Ratio_EC

	E_Internat_Substitution_CoalForBioenergy_Pellet=bS['PelletExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletExport
	E_Internat_Substitution_DieselForBioenergy_Pellet=bS['PelletExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletExport
	E_Internat_Substitution_GasForBioenergy_Pellet=bS['PelletExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletExport
	E_Internat_Substitution_OilForBioenergy_Pellet=bS['PelletExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletExport

	#----------------------------------------------------------------------
	# Pellet domestic grid (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dired tonnes/ha)
	Yield_PelletDomGrid=vo['C_ToBBP_PelletDomGrid']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletDomGrid=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomGrid

	GJ_PelletDomGrid=GJ_PelletDomGrid*Ratio_EC

	E_Domestic_Substitution_CoalForBioenergy_PelletElectricityGrid=bS['PelletDomGridFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomGrid
	E_Domestic_Substitution_DieselForBioenergy_PelletElectricityGrid=bS['PelletDomGridFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomGrid
	E_Domestic_Substitution_GasForBioenergy_PelletElectricityGrid=bS['PelletDomGridFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomGrid
	E_Domestic_Substitution_OilForBioenergy_PelletElectricityGrid=bS['PelletDomGridFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomGrid

	#----------------------------------------------------------------------
	# Pellet domestic RNG (MgC/ha) to (kiln dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (kiln dired tonnes/ha)
	Yield_PelletDomRNG=vo['C_ToBBP_PelletDomRNG']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_PelletDomRNG=bB['Energy Content Wood (Kiln-dried)']*Yield_PelletDomRNG

	E_Domestic_Substitution_CoalForBioenergy_PelletRNG=bS['PelletDomRNGFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_PelletDomRNG
	E_Domestic_Substitution_DieselForBioenergy_PelletRNG=bS['PelletDomRNGFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_PelletDomRNG
	E_Domestic_Substitution_GasForBioenergy_PelletRNG=bS['PelletDomRNGFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_PelletDomRNG
	E_Substitution_OilForBioenergy_PelletDomRNG=bS['PelletDomRNGFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_PelletDomRNG

	#----------------------------------------------------------------------
	# Domestic firewood (MgC/ha) to (air dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (air dried tonnes/ha)
	Yield_FirewoodDom=vo['C_ToBBP_FirewoodDom']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_FirewoodDom=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodDom

	E_Domestic_Substitution_CoalForBioenergy_Firewood=bS['FirewoodDomFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodDom
	E_Domestic_Substitution_DieselForBioenergy_Firewood=bS['FirewoodDomFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodDom
	E_Domestic_Substitution_GasForBioenergy_Firewood=bS['FirewoodDomFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodDom
	E_Domestic_Substitution_OilForBioenergy_Firewood=bS['FirewoodDomFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodDom

	#----------------------------------------------------------------------
	# International firewood (MgC/ha) to (air dried tonne/ha)
	#----------------------------------------------------------------------

	# Yield (air dried tonnes/ha)
	Yield_FirewoodExport=vo['C_ToBBP_FirewoodExport']/bB['Density Wood']

	# Yield to energy (GJ/ha)
	GJ_FirewoodExport=bB['Energy Content Wood (0% moisture)']*Yield_FirewoodExport

	E_Internat_Substitution_CoalForBioenergy_Firewood=bS['FirewoodExportFracDisplacingCoal']*bB['Emission Intensity Coal']/1000*GJ_FirewoodExport
	E_Internat_Substitution_DieselForBioenergy_Firewood=bS['FirewoodExportFracDisplacingDiesel']*bB['Emission Intensity Diesel']/1000*GJ_FirewoodExport
	E_Internat_Substitution_GasForBioenergy_Firewood=bS['FirewoodExportFracDisplacingNaturalGas']*bB['Emission Intensity Natural Gas']/1000*GJ_FirewoodExport
	E_Internat_Substitution_OilForBioenergy_Firewood=bS['FirewoodExportFracDisplacingOil']*bB['Emission Intensity Oil']/1000*GJ_FirewoodExport

	#--------------------------------------------------------------------------
	# Substitution of fossil fuels for bioenergy
	# *** Save as positive and then change sign in post-processing ***
	#--------------------------------------------------------------------------

	vo['E_Domestic_Substitution_CoalForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_CoalForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_CoalForBioenergy_PelletElectricityGrid + \
		 E_Domestic_Substitution_CoalForBioenergy_PelletRNG + \
		 E_Domestic_Substitution_CoalForBioenergy_Firewood)

	vo['E_Internat_Substitution_CoalForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_CoalForBioenergy_PowerFacility + \
		 E_Internat_Substitution_CoalForBioenergy_Pellet + \
		 E_Internat_Substitution_CoalForBioenergy_Firewood)

	vo['E_Domestic_Substitution_OilForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_OilForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_OilForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_OilForBioenergy_Firewood + \
		 E_Domestic_Substitution_DieselForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_DieselForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_DieselForBioenergy_PelletElectricityGrid + \
		 E_Domestic_Substitution_DieselForBioenergy_PelletRNG + \
		 E_Domestic_Substitution_DieselForBioenergy_Firewood)

	vo['E_Internat_Substitution_OilForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_OilForBioenergy_PowerFacility + \
		 E_Internat_Substitution_OilForBioenergy_Pellet + \
		 E_Internat_Substitution_OilForBioenergy_Firewood + \
		 E_Internat_Substitution_DieselForBioenergy_PowerFacility + \
		 E_Internat_Substitution_DieselForBioenergy_Pellet + \
		 E_Internat_Substitution_DieselForBioenergy_Firewood)

	vo['E_Domestic_Substitution_GasForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_GasForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_GasForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_GasForBioenergy_PelletElectricityGrid + \
		 E_Domestic_Substitution_GasForBioenergy_PelletRNG + \
		 E_Domestic_Substitution_GasForBioenergy_Firewood)

	vo['E_Internat_Substitution_GasForBioenergy']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_GasForBioenergy_PowerFacility + \
		 E_Internat_Substitution_GasForBioenergy_Pellet + \
		 E_Internat_Substitution_GasForBioenergy_Firewood)

	vo['E_Domestic_Substitution_PowerFacility']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_OilForBioenergy_PowerFacility + \
		 E_Domestic_Substitution_GasForBioenergy_PowerFacility)

	vo['E_Internat_Substitution_PowerFacility']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_CoalForBioenergy_PowerFacility + \
		 E_Internat_Substitution_OilForBioenergy_PowerFacility + \
		 E_Internat_Substitution_GasForBioenergy_PowerFacility)

	vo['E_Domestic_Substitution_ElectricityGrid']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_OilForBioenergy_ElectricityGrid + \
		 E_Domestic_Substitution_GasForBioenergy_ElectricityGrid)

	vo['E_Internat_Substitution_Pellet']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_CoalForBioenergy_Pellet + \
		 E_Internat_Substitution_OilForBioenergy_Pellet + \
		 E_Internat_Substitution_GasForBioenergy_Pellet)

	vo['E_Domestic_Substitution_PelletElectricityGrid']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_PelletElectricityGrid + \
		 E_Domestic_Substitution_OilForBioenergy_PelletElectricityGrid + \
		 E_Domestic_Substitution_GasForBioenergy_PelletElectricityGrid)

	vo['E_Domestic_Substitution_PelletRNG']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_PelletRNG + \
		 E_Substitution_OilForBioenergy_PelletDomRNG + \
		 E_Domestic_Substitution_GasForBioenergy_PelletRNG)

	vo['E_Domestic_Substitution_Firewood']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Domestic_Substitution_CoalForBioenergy_Firewood + \
		 E_Domestic_Substitution_OilForBioenergy_Firewood + \
		 E_Domestic_Substitution_GasForBioenergy_Firewood)

	vo['E_Internat_Substitution_Firewood']= (1-bS['Economic Contraction Fraction'])* \
		 (E_Internat_Substitution_CoalForBioenergy_Firewood + \
		 E_Internat_Substitution_OilForBioenergy_Firewood + \
		 E_Internat_Substitution_GasForBioenergy_Firewood)

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
	vo['ODT_Concrete']=fS_Concrete+fP_Concrete#-fR

	# Substitution of steel for structural wood
	fS_Steel=bS['SawnwoodFracDisplacingSteel']*bS['DisplacementRatio_SteelForSawnwood']*ODT_Sawnwood
	fP_Steel=bS['PanelFracDisplacingSteel']*bS['DisplacementRatio_SteelForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingSteel']*bS['DisplacementRatio_SteelForResiduals']*Residuals
	vo['ODT_Steel']=fS_Steel+fP_Steel#-fR

	# Substitution of aluminum for structural wood
	fS_Aluminum=bS['SawnwoodFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForSawnwood']*ODT_Sawnwood
	fP_Aluminum=bS['PanelFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingAluminum']*bS['DisplacementRatio_AluminumForResiduals']*Residuals
	vo['ODT_Aluminum']=fS_Aluminum+fP_Aluminum#-fR

	# Substitution of plastics for structural wood
	fS_Plastic=bS['SawnwoodFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForSawnwood']*ODT_Sawnwood
	fP_Plastic=bS['PanelFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingPlastic']*bS['DisplacementRatio_PlasticForResiduals']*Residuals
	vo['ODT_Plastic']=fS_Plastic+fP_Plastic#-fR

	# Substitution of textiles for structural wood
	fS_Textile=bS['SawnwoodFracDisplacingTextile']*bS['DisplacementRatio_TextileForSawnwood']*ODT_Sawnwood
	fP_Textile=bS['PanelFracDisplacingTextile']*bS['DisplacementRatio_TextileForPanel']*ODT_Panel
	#fR=0#bS['ResidualsFracDisplacingTextile']*bS['DisplacementRatio_TextileForResiduals']*Residuals
	vo['ODT_Textile']=fS_Textile+fP_Textile#-fR

	# Emissions from productoin of structural materials
	vo['E_Internat_Substitution_Concrete']=(1-bS['Economic Contraction Fraction'])*bB['Emission Intensity Concrete']*vo['ODT_Concrete']
	vo['E_Internat_Substitution_ConcreteFromCalcination']=(1-bS['Economic Contraction Fraction'])*bS['FracConcreteEmissionsFromCalcination']*vo['E_Internat_Substitution_Concrete']
	vo['E_Internat_Substitution_ConcreteFromNonCalcination']=(1-bS['Economic Contraction Fraction'])*(1-bS['FracConcreteEmissionsFromCalcination'])*vo['E_Internat_Substitution_Concrete']
	vo['E_Internat_Substitution_Steel']=(1-bS['Economic Contraction Fraction'])*bB['Emission Intensity Steel']*vo['ODT_Steel']
	vo['E_Internat_Substitution_Aluminum']=(1-bS['Economic Contraction Fraction'])*bB['Emission Intensity Aluminum']*vo['ODT_Aluminum']
	vo['E_Internat_Substitution_Plastic']=(1-bS['Economic Contraction Fraction'])*bB['Emission Intensity Plastic']*vo['ODT_Plastic']
	vo['E_Internat_Substitution_Textile']=(1-bS['Economic Contraction Fraction'])*bB['Emission Intensity Textile']*vo['ODT_Textile']

	# Emissions from sawnwood and Panel
	vo['E_Internat_Substitution_Sawnwood']=(1-bS['Economic Contraction Fraction'])*(bB['Emission Intensity Concrete']*fS_Concrete+ \
		bB['Emission Intensity Steel']*fS_Steel + \
		bB['Emission Intensity Aluminum']*fS_Aluminum + \
		bB['Emission Intensity Plastic']*fS_Plastic + \
		bB['Emission Intensity Textile']*fS_Textile)

	vo['E_Internat_Substitution_Panel']=(1-bS['Economic Contraction Fraction'])*(bB['Emission Intensity Concrete']*fP_Concrete+ \
		bB['Emission Intensity Steel']*fP_Steel + \
		bB['Emission Intensity Aluminum']*fP_Aluminum + \
		bB['Emission Intensity Plastic']*fP_Plastic + \
		bB['Emission Intensity Textile']*fP_Textile)

	# Emissions from structural matierials, tallied by feedstock (and calcination)
	vo['E_Internat_Substitution_CoalForSolidWood']= \
		bS['FracConcreteEmissionsFromCoal']*vo['E_Internat_Substitution_Concrete']+ \
		bS['FracSteelEmissionsFromCoal']*vo['E_Internat_Substitution_Steel']+ \
		bS['FracAluminumEmissionsFromCoal']*vo['E_Internat_Substitution_Aluminum']+ \
		bS['FracPlasticEmissionsFromCoal']*vo['E_Internat_Substitution_Plastic']+ \
		bS['FracTextileEmissionsFromCoal']*vo['E_Internat_Substitution_Textile']

	vo['E_Internat_Substitution_OilForSolidWood']= \
		bS['FracConcreteEmissionsFromOil']*vo['E_Internat_Substitution_Concrete']+ \
		bS['FracSteelEmissionsFromOil']*vo['E_Internat_Substitution_Steel']+ \
		bS['FracAluminumEmissionsFromOil']*vo['E_Internat_Substitution_Aluminum']+ \
		bS['FracPlasticEmissionsFromOil']*vo['E_Internat_Substitution_Plastic']+ \
		bS['FracTextileEmissionsFromOil']*vo['E_Internat_Substitution_Textile']

	vo['E_Internat_Substitution_GasForSolidWood']= \
		bS['FracConcreteEmissionsFromGas']*vo['E_Internat_Substitution_Concrete']+ \
		bS['FracSteelEmissionsFromGas']*vo['E_Internat_Substitution_Steel']+ \
		bS['FracAluminumEmissionsFromGas']*vo['E_Internat_Substitution_Aluminum']+ \
		bS['FracPlasticEmissionsFromGas']*vo['E_Internat_Substitution_Plastic']+ \
		bS['FracTextileEmissionsFromGas']*vo['E_Internat_Substitution_Textile']

	#--------------------------------------------------------------------------
	# Back-calculate production of fossil fuel consumption from operational use
	# and substitution effects (GJ)
	# *** Moved to post-processing ***
	#--------------------------------------------------------------------------

	return vo

#%%
def GrassBiomassDynamics(meta,pNam,iScn,iBat,iT,vi,vo,iEP):
	# Aboveground biomass (foliage)
	vo['C_G_Gross_ByPool'][iT,:,iEP['Foliage']]=vo['C_G_Gross_ByPool'][iT,:,iEP['Foliage']]+1.5
	vo['C_LF_ByPool'][iT,:,iEP['Foliage']]=vo['C_LF_ByPool'][iT,:,iEP['Foliage']]+1.5

	# Belowground biomass
	if vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]>15:
		vo['C_G_Gross_ByPool'][iT,:,iEP['RootCoarse']]=+1.5
	else:
		vo['C_G_Gross_ByPool'][iT,:,iEP['RootCoarse']]=+1.75

	vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]=+1.5

	# Net growth
	G_net=vo['C_G_Gross_ByPool'][iT,:,iEP['RootCoarse']]-vo['C_LF_ByPool'][iT,:,iEP['RootCoarse']]

	#vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']]=vo['C_G_Net_Reg_ByPool'][iT,:,iEP['RootCoarse']]+G_net

	vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]=vo['C_Eco_ByPool'][iT,:,iEP['RootCoarse']]+G_net

	return vo

