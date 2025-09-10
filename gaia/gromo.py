

#%%
def TreeBiomassDynamicsFromGROMO1(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

	# Update stand age
	vo['A'][iT,:]=vo['A'][iT-1,:]+1

	# Extract incides for GY model data structure
	iGY=meta['Modules']['GYM']['GC Input Indices']

	# Predictor variables
	A=vo['A'][iT,:]
	DAB=0; DAD=0; DAF=0; DAP=0;
	IH=0;
	IPN=1;
	#IH=1*(vi['lsat']['Harvest Year Comp2']>0)
	FDC=1

	#--------------------------------------------------------------------------
	# Gross growth of stemwood
	#--------------------------------------------------------------------------

	flg_GG='Chen'
	if flg_GG=='Chapman Richards':
		# Chapman Richards
		bG=meta['Param']['Raw']['GROMO']['gg1']['Chapman Richards']['Param']['BE']
		B1=bG['fdc1']+bG['Sd1']*IPN+bG['Hi1']*IH+bG['SdHi1']*IPN*IH+bG['dab1']*DAB+bG['dad1']*DAD+bG['daf1']*DAF+bG['dap1']*DAP
		B2=bG['fdc2']
		B3=bG['b3']
		G_Gross_StemTot=B1*B2*B3*np.exp(-B2*A)*(1-np.exp(-B2*A))**(B3-1)
	elif flg_GG=='Chen':
		# Chen
		bG=meta['Param']['Raw']['GROMO']['gg1']['Chen']['Param']['BE']
		B1=bG['fdc1']+bG['Sd1']*IPN+bG['Hi1']*IH+bG['SdHi1']*IPN*IH+bG['dab1']*DAB+bG['dad1']*DAD+bG['daf1']*DAF+bG['dap1']*DAP
		B2=bG['fdc2']
		B3=bG['b3']
		B4=bG['b4']
		G_Gross_StemTot=B1*(1+((B2*(A/B3)**B4-1)/np.exp(A/B3)))

	G_Gross_StemTot=np.maximum(0.1,G_Gross_StemTot)

	# Non-merch to merch ratio
	#Vtot_lag1=(vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']])/meta['Param']['BEV']['Biophysical']['Density Wood']/meta['Param']['BEV']['Biophysical']['Carbon Content Wood']
	fStemNonMerch=0.08#np.minimum(1,(1-0.086)*0.086+(1/np.exp(0.012*(Vtot_lag1-14.4))))

	# Update gross growth
	vo['C_G_Gross_ByPool'][iT,:,iEP['StemMerch']]=(1-fStemNonMerch)*G_Gross_StemTot
	vo['C_G_Gross_ByPool'][iT,:,iEP['StemNonMerch']]=fStemNonMerch*G_Gross_StemTot

	#--------------------------------------------------------------------------
	# Mortality of stemwood
	#--------------------------------------------------------------------------

	bM=meta['Param']['Raw']['GROMO']['m1']['Param']['BE']
	zs=meta['Param']['Raw']['GROMO']['m1']['Zscore Stats']

	B=vo['C_Eco_ByPool'][iT-1,:,iEP['StemMerch']]+vo['C_Eco_ByPool'][iT-1,:,iEP['StemNonMerch']]
	B_z=(B-zs['B']['mu'])/zs['B']['sig']
	A_z=(A-zs['A']['mu'])/zs['A']['sig']
	A_z_sq=A_z**2

	BTR_Reg_StemTot=bM['Intercept']+bM['C(LS)[T.FDC]']*FDC+bM['IPN']*IPN+bM['IH']*IH+bM['B_z']*B_z+bM['A_z']*A_z+bM['A_z_sq']*A_z_sq

	M_Reg_StemTot=BTR_Reg_StemTot/100*B

	M_Reg_StemTot=np.maximum(0,M_Reg_StemTot)

	# Constrain mortality by the carbon pool at the start of the time step
	M_Reg_StemTot=np.minimum(M_Reg_StemTot,B)

	vo['C_M_Reg_ByPool'][iT,:,0:7]=np.tile(BTR_Reg_StemTot/100,(7,1)).T*vo['C_Eco_ByPool'][iT-1,:,0:7]

	#--------------------------------------------------------------------------
	# Net growth
	#--------------------------------------------------------------------------

	# Net growth by pool
	NetGrowthByPool=np.zeros( (meta[pNam]['Project']['Batch Size'][iBat],6) )

	# Net growth of total stemwood
	G_Net_StemTot=G_Gross_StemTot-M_Reg_StemTot

	# Net growth of merch stemwood and non-merch stemwood
	fStemNonMerch=0.08 #np.minimum(1,(1-0.086)*0.086+(1/np.exp(0.012*(Vtot_lag1-14.4))))
	NetGrowthByPool[:,iEP['StemMerch']]=(1-fStemNonMerch)*G_Net_StemTot
	NetGrowthByPool[:,iEP['StemNonMerch']]=fStemNonMerch*G_Net_StemTot

	# Net growth of foliage
	#NetGrowthByPool[:,iGY['Foliage']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gf1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gf2']-meta['Param']['BEV']['BiomassAllometrySL']['Gf1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gf3']*vo['A'][iT,:]))

	# Net growth of branches
	#NetGrowthByPool[:,iGY['Branch']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gbr1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbr2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbr1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbr3']*vo['A'][iT,:]))

	# Net growth of bark
	#NetGrowthByPool[:,iGY['Bark']]=G_Net_StemTot*(meta['Param']['BEV']['BiomassAllometrySL']['Gbk1']+(meta['Param']['BEV']['BiomassAllometrySL']['Gbk2']-meta['Param']['BEV']['BiomassAllometrySL']['Gbk1'])*np.exp(-meta['Param']['BEV']['BiomassAllometrySL']['Gbk3']*vo['A'][iT,:]))

	NetGrowthByPool[:,iGY['Foliage']]=0.08*G_Net_StemTot
	NetGrowthByPool[:,iGY['Branch']]=0.14*G_Net_StemTot
	NetGrowthByPool[:,iGY['Bark']]=0.05*G_Net_StemTot

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
		meta['Modules']['NutrientApplication']['iApplication']=np.where( meta['Modules']['NutrientApplication']['ResponseCounter']>0 )[0]
	
		# Adjust N application response counter
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'UpdateCounter')
	
		# Adjust root net growth
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'BelowgroundNetGrowth')

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
	#vo['C_M_Reg_ByPool'][iT,:,0:7]=vo['C_G_Gross_ByPool'][iT,:,0:7]-vo['C_G_Net_Reg_ByPool'][iT,:,0:7]

	if meta[pNam]['Project']['Biomass Module']!='gromo1':
		# Adjust mortality to account for N application response
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'Mortality')

	#--------------------------------------------------------------------------
	# Update gross growth
	#--------------------------------------------------------------------------
	vo['C_G_Gross_ByPool'][iT,:,0:7]=vo['C_G_Net_Reg_ByPool'][iT,:,0:7]+vo['C_M_Reg_ByPool'][iT,:,0:7]
	#vo['C_G_Gross_ByPool'][iT,:,iEP['Bark']]=0.14*G_Gross_StemTot
	#vo['C_G_Gross_ByPool'][iT,:,iEP['Branch']]=0.27*G_Gross_StemTot
	#vo['C_G_Gross_ByPool'][iT,:,iEP['Foliage']]=0.17*G_Gross_StemTot
	#vo['C_G_Gross_ByPool'][iT,:,iEP['RootCoarse']]=0.14*G_Gross_StemTot
	#vo['C_G_Gross_ByPool'][iT,:,iEP['RootFine']]=0.17*G_Gross_StemTot

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
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'Litterfall')

	return vo

#%%
def TreeBiomassDynamicsFromGROMO2(meta,pNam,iScn,iBat,iT,vi,vo,iEP):

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
		meta['Modules']['NutrientApplication']['iApplication']=np.where( meta['Modules']['NutrientApplication']['ResponseCounter']>0 )[0]
	
		# Adjust N application response counter
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'UpdateCounter')
	
		# Adjust root net growth
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'BelowgroundNetGrowth')

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
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'Mortality')

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
		if meta['Modules']['NutrientApplication']['iApplication'].size>0:
			vi,vo,meta=napp.NutrientApplicationResponse(meta,pNam,vi,vo,iT,iScn,'Litterfall')

	return vo