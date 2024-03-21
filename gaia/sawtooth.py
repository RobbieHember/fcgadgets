#%%
def BiomassFromSawtooth(meta,pNam,iScn,iS,vi,vo,iEP):

	t0=time.time()

	# Turn warnings off to avoid divide by zero warning
	warnings.filterwarnings('ignore')

	# *************************************************************************
	# What if we only track the current and previous time step - memory useage would
	# be low enough to vectorize stands
	# Requires constant population of stand level summary variables at the end
	# of each time step.
	#**************************************************************************

	# Define demensions
	N_time=meta[pNam]['Project']['N Time']

	# Maximum number of trees per stand
	N_tree=meta['Param']['BE']['Sawtooth']['Core']['Max SPH'].astype(int)

	# Multiplier to convert kg C to Mg C
	cm=1000

	#--------------------------------------------------------------------------
	# Initialize tree-level variables
	#--------------------------------------------------------------------------

	# List of variables
	vL=['ID_SRS','ID_Decid','A','H','D','N_R','Csw','Csw_Larger','Csw_G',
		'N_M_Reg','N_M_Fir','N_M_Ins','N_M_Pat','N_M_Har','N_M_Win',
		'Csw_M_Reg','Csw_M_Fir','Csw_M_Ins','Csw_M_Pat','Csw_M_Har','Csw_M_Win']

	# Populate
	tl={}
	for v in vL:
		if v[0:2]=='ID':
			tl[v]=np.zeros((N_time,N_tree),dtype=np.int)
		else:
			tl[v]=np.zeros((N_time,N_tree),dtype=np.float)

	#--------------------------------------------------------------------------
	# Populate species ID
	#--------------------------------------------------------------------------

	# Generate random number vector
	rp=np.random.permutation(N_tree).astype(int)

	# Species 1
	n1=np.ceil(vi['lsat']['SRS1_PCT'][0,iS]/100*N_tree).astype(int)
	tl['ID_SRS'][:,rp[0:n1]]=vi['lsat']['SRS1_ID'][0,iS]

	# Species 2
	if vi['lsat']['SRS2_PCT'][0,iS]>0:
		n2=np.ceil(vi['lsat']['SRS2_Pct'][0,iS]/100*N_tree).astype(int)
		tl['ID_SRS'][:,rp[n1:n2]]=vi['lsat']['SRS2_ID'][0,iS]

	# Species 3
	if vi['lsat']['SRS3_PCT'][0,iS]>0:
		tl['ID_SRS'][:,rp[n2:]]=vi['lsat']['SRS3_ID'][0,iS]

	#--------------------------------------------------------------------------
	# Initialize parameter vectors
	#--------------------------------------------------------------------------

	bA={}
	for k in meta['Param']['BE']['Sawtooth']['Allom'].keys():
		bA[k]=np.zeros(N_tree)

	bR={}
	for k in meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ].keys():
		bR[k]=np.zeros(N_tree)

	bM={}
	for k in meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ].keys():
		bM[k]=np.zeros(N_tree)

	bG={}
	for k in meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ].keys():
		bG[k]=np.zeros(N_tree)

	#--------------------------------------------------------------------------
	# Populate parameter vectors based on species to each tree based on fractions
	# from inventory.
	#--------------------------------------------------------------------------

	uS=np.unique(tl['ID_SRS'])
	for iU in range(uS.size):

		ind0=np.where(tl['ID_SRS'][0,:]==uS[iU])[0]

		ind1=np.where(meta['Param']['BE']['Sawtooth']['Allom']['SRS_ID']==uS[iU])[0]
		for k in bA.keys():
			if (k=='SRS_ID') | (k=='SRS_CD'):
				continue
			bA[k][ind0]=meta['Param']['BE']['Sawtooth']['Allom'][k][ind1]

		ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ]['SRS_ID']==uS[iU])[0]
		for k in bR.keys():
			if (k=='SRS_ID') | (k=='SRS_CD'):
				continue
			bR[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq R'][ meta[pNam]['Scenario'][iScn]['Eq R CD'] ][k][ind1]

		ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ]['SRS_ID']==uS[iU])[0]
		for k in bM.keys():
			if (k=='SRS_ID') | (k=='SRS_CD'):
				continue
			bM[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq M'][ meta[pNam]['Scenario'][iScn]['Eq M CD'] ][k][ind1]

		ind1=np.where(meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ]['SRS_ID']==uS[iU])[0]
		for k in bG.keys():
			if (k=='SRS_ID') | (k=='SRS_CD'):
				continue
			bG[k][ind0]=meta['Param']['BE']['Sawtooth']['Eq G'][ meta[pNam]['Scenario'][iScn]['Eq G CD'] ][k][ind1]

	#--------------------------------------------------------------------------
	# Germination
	#--------------------------------------------------------------------------

	# In the first time step, assume 1...seed_n seedlings germinate with initial
	# biomass drawn from a normal distribution, with a very small mean
	# biomass and a coefficient of variation of 10 percent. Force values to be positive.
	# Once a few live trees have been Invscribed, recruitment rate should
	# rise quickly because stand biomass (the primary driver of recruitment)
	# is low. This should lead to repopulation of the stand assuming mortality
	# does not exceed recruitment, and growth is viable.

	# Planted:
	# seed_n=1100;
	# seed_mu=0.06; % (kg C tree-1)

	# Natural:
	seed_n=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant SPH'].astype(int)
	seed_mass_mu=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass Mean'] # (kg C tree-1)
	seed_mass_sd=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass SD'] # (kg C tree-1)
	seed_min=meta['Param']['BE']['Sawtooth']['Core']['Natural Germinant Mass Min'] # (kg C tree-1)
	Csw_seed=np.random.normal(seed_mass_mu,seed_mass_sd,(1,seed_n))
	Csw_seed=np.maximum(seed_min,Csw_seed)

	# Populate
	tl['Csw'][0,0:seed_n]=Csw_seed

	# *** This needs updating ***
	tl['H'][0,0:seed_n]=bA['Cag2H1'][0:seed_n]*((1-np.exp(-bA['Cag2H2'][0:seed_n]*tl['Csw'][0,0:seed_n])))**(1/(1-bA['Cag2H3'][0:seed_n]))

	#--------------------------------------------------------------------------
	# Reshape allometric parameters
	#--------------------------------------------------------------------------

	# Foliage biomass
	Csw2Cf1=np.tile( np.reshape(bA['Csw2Cf1'],(1,-1)),(N_time,1) )
	Csw2Cf2=np.tile( np.reshape(bA['Csw2Cf2'],(1,-1)),(N_time,1) )

	# Branch biomass
	Csw2Cbk1=np.tile( np.reshape(bA['Csw2Cbk1'],(1,-1)),(N_time,1) )
	Csw2Cbk2=np.tile( np.reshape(bA['Csw2Cbk2'],(1,-1)),(N_time,1) )

	# Branch biomass
	Csw2Cbr1=np.tile( np.reshape(bA['Csw2Cbr1'],(1,-1)),(N_time,1) )
	Csw2Cbr2=np.tile( np.reshape(bA['Csw2Cbr2'],(1,-1)),(N_time,1) )

	#--------------------------------------------------------------------------
	# Function for calculating total biomass from stemwood biomass
	#--------------------------------------------------------------------------

	def ExpandBiomass(Csw_tot,Csw_merch):

		# Merchantable stemwood
		SWm=Csw_merch*meta['Param']['BE']['Sawtooth']['Core']['Merch Stemwood Fraction']

		# Non-merchantable stemwood
		SWnm=Csw_tot-SWm

		# Foliage biomass
		F=Csw_tot*Csw2Cf1*Csw_tot**Csw2Cf2

		# Bark biomass
		Bk=Csw_tot*Csw2Cbk1*Csw_tot**Csw2Cbk2

		# Branch biomass
		Br=Csw_tot*Csw2Cbr1*Csw_tot**Csw2Cbr2

		# Total aboveground biomass
		AG=Csw_tot+F+Bk+Br

		# Total root biomass (Li et al. 2003)

		# Conifer
		R=0.222*AG

		# Deciduous (equation from Li et al. XXXX)
		ind=np.where(tl['ID_Decid']==1)[0]
		if ind.size!=0:
			R[ind]=1.576*AG[ind]**0.615

		# Fine root biomass (equation from Li et al. XXXX)
		Rf=R*(0.072+0.354*np.exp(-0.06*(2*R)))

		# Coarse root biomass
		Rc=R-Rf

		# Total biomass
		Tot=AG+R

		return SWm,SWnm,F,Bk,Br,Rc,Rf,Tot

	#--------------------------------------------------------------------------
	# Loop through time intervals (start in second time step)
	#--------------------------------------------------------------------------

	for iT in range(1,N_time):

		#----------------------------------------------------------------------
		# Apply first disturbance cycle to subsequent disturbance cycles during
		# the spinup period to save time.
		# *** Not working properly - just needs maintenance ***
		#----------------------------------------------------------------------

#		if (vi['tv'][iT]>meta[pNam]['Project']['SpinupSpanFastTrack'][iScn]['Start']+1) & (vi['tv'][iT]<=meta[pNam]['Project']['SpinupSpanFastTrack'][iScn]['End']+1):
#
#			iT0=iT-meta[pNam]['Project']['Spinup Disturbance Return Inverval']
#
#			tl['A'][iT,:]=tl['A'][iT0,:]
#			tl['H'][iT,:]=tl['H'][iT0,:]
#			tl['D'][iT,:]=tl['D'][iT0,:]
#			tl['N_R'][iT,:]=tl['N_R'][iT0,:]
#			tl['N_M_Reg'][iT,:]=tl['N_M_Reg'][iT0,:]
#			tl['Csw'][iT,:]=tl['Csw'][iT0,:]
#			tl['Csw_Larger'][iT,:]=tl['Csw_Larger'][iT0,:]
#			tl['Csw_G'][iT,:]=tl['Csw_G'][iT0,:]
#			tl['Csw_M_Reg'][iT,:]=tl['Csw_M_Reg'][iT0,:]
#			vo['A'][iT,iS]=vo['A'][iT0,iS]
#
#			continue

		#----------------------------------------------------------------------

		# Index to live trees
		iLive=np.array(np.where(~np.isnan(tl['Csw'][iT-1,:]))).flatten()
		nLive=iLive.size

		# Index to dead trees
		#iDead=np.array(np.where(np.isnan(tl['Csw'][iT-1,:])))
		iDead=np.where(np.isnan(tl['Csw'][iT-1,:])==True)[0]
		nDead=iDead.size

		# Create random number vectors that will be used for calculation of
		# annual probabilities of recruitment and mortality
		rLive=np.random.uniform(0,1,nLive)
		rDead=np.random.uniform(0,1,nDead)

		# Update stand age (i.e., time since stand-replacing disturbance)
		vo['A'][iT,iS]=vo['A'][iT-1,iS]+1

		# Update tree age
		tl['A'][iT,iLive]=tl['A'][iT-1,iLive]+1

		tmp=tl['D'][iT-1,iLive]
		tmp=tmp[tmp>=0]
		try:
			kde=stats.gaussian_kde(tmp)
		except:
			print(tmp.shape)
		try:
			vo['DBH_Class'][iT,iS,:]=kde(meta['Core']['Sawtooth']['DBH Classes'])
		except:
			print(vo['DBH_Class'][iT,iS,:].shape)
			print(meta['Core']['Sawtooth']['DBH Classes'].shape)

		#data=np.random.random(100)
		#kde=stats.gaussian_kde(data)
		#a=kde(np.linspace(-2,2,100))

		#----------------------------------------------------------------------
		# Calculate predictor variables used in the equations of recruitment,
		# mortality and growth
		#----------------------------------------------------------------------

		# Tree age
		A=tl['A'][iT-1,:]

		# Aboveground biomass of individual trees from t-1 (kg C tree-1)
		Csw=tl['Csw'][iT-1,:]

		# Tree height (m)
		H=tl['H'][iT-1,:]

		# Stand age
		SA=np.nanmean(tl['A'][iT-1,:])

		# Stand-level biomass from t-1 (Mg C ha-1)
		SCsw=np.nansum(tl['Csw'][iT-1,:])/cm

		# Stand density from t-1 (stems ha-1)
		SN=nLive

		# Biomass of larger trees (Mg C ha-1)
		ListToSort=np.zeros((N_tree,2))
		ListToSort[:,0]=np.arange(0,N_tree,1)
		ListToSort[:,1]=Csw/cm
		SortedList=ListToSort[ListToSort[:,1].argsort()]
		SortedList=np.flip(SortedList,0)
		tmp=np.reshape(np.nancumsum(SortedList[:,1])-SortedList[:,1]/1000,(N_tree,1))
		SortedList2=np.append(SortedList,tmp,axis=1)

		Csw_Larger=np.zeros(N_tree)
		Csw_Larger[SortedList2[:,0].astype(int)]=SortedList2[:,2]

		tl['Csw_Larger'][iT,:]=Csw_Larger

		#----------------------------------------------------------------------
		# Probability of recruitment
		#----------------------------------------------------------------------

		# Only do this if there are living trees
		if iDead.size!=0:

			if meta[pNam]['Scenario'][iScn]['Eq R CD']=='Def1':
				SCsw_z=(SCsw-bR['SB_mu'])/bR['SB_sig']
				lgit=bR['Int']+bR['SB']*SCsw_z
				Pr=(np.exp(lgit)/(1+np.exp(lgit)))

			# Establish trees based on annual probability of recruitment.
			# Initial values of biomass and height are set low, arbitrary until
			# they can be set according to actual observations
			# at end of first year.

			iRec=np.where(Pr[iDead]>=rDead)[0]

			tl['A'][iT:,iDead[iRec]]=1
			tl['H'][iT:,iDead[iRec]]=0.1
			tl['Csw'][iT:,iDead[iRec]]=0.05
			tl['N_R'][iT,iDead[iRec]]=1

		#----------------------------------------------------------------------
		# Growth of stemwood biomass (kg C tree-1 yr-1)
		#----------------------------------------------------------------------

		if meta[pNam]['Scenario'][iScn]['Eq G CD']=='Def1':

			# Standardization
			LnCsw_z=(np.log(Csw)-bG['LnB_mu'])/bG['LnB_sig']
			Csw_z=(Csw-bG['B_mu'])/bG['B_sig']
			SA_z=(A-bG['SA_mu'])/bG['SA_sig']
			SCswLT_z=(Csw_Larger-bG['SBLT_mu'])/bG['SBLT_sig']
			SCsw_z=(SCsw-bG['SB_mu'])/bG['SB_sig']

			# Add all effects to intercept
			yhat=bG['Int'] + bG['LnB']*LnCsw_z + bG['B']*Csw_z + bG['SA']*SA_z + bG['SBLT']*SCswLT_z + bG['SB']*SCsw_z

			# Back-transform and apply log correction
			yhat=np.exp(yhat)

			# Cap unrealistic growth
			#yhat[yhat>G_max]=G_max

			# Populate tree level structure with growth predictions
			tl['Csw_G'][iT,iLive]=yhat[iLive]

		#----------------------------------------------------------------------
		# Update state variables
		#----------------------------------------------------------------------

		# Update stemwood biomass (kg C tree-1)
		tl['Csw'][iT,iLive]=tl['Csw'][iT-1,iLive]+tl['Csw_G'][iT,iLive]

		# Update tree height (m)
		# *** This needs updating ***
		tl['H'][iT,:]=bA['Cag2H1']*((1-np.exp(-bA['Cag2H2']*tl['Csw'][iT,:]))**(1/(1-bA['Cag2H3'])))

		# Update diameter (cm)
		tl['D'][iT,:]=(tl['Csw'][iT,:]/0.5/meta['Param']['BE']['Sawtooth']['Core']['D_to_Bsw_b0'])**(1/meta['Param']['BE']['Sawtooth']['Core']['D_to_Bsw_b1'])

		#----------------------------------------------------------------------
		# Probability of tree mortality (regular)
		#----------------------------------------------------------------------

		if meta[pNam]['Scenario'][iScn]['Eq M CD']=='Def1':

			Csw_z=(Csw-bM['B_mu'])/bM['B_sig']
			Csw2_z=(Csw**2-bM['B2_mu'])/bM['B2_sig']
			SA_z=(A-bM['SA_mu'])/bM['SA_sig']
			SCswLT_z=(Csw_Larger-bM['SBLT_mu'])/bM['SBLT_sig']
			SCsw_z=(SCsw-bM['SB_mu'])/bM['SB_sig']

			lgit=bM['Int'] + bM['B']*Csw_z + bM['B2']*Csw2_z + bM['SA']*SA_z + bM['SBLT']*SCswLT_z + bM['SB']*SCsw_z

			Pm_Sim_Reg=1*(np.exp(lgit)/(1+np.exp(lgit)))
			Pm_Sim_Ins=np.zeros(N_tree)
			Pm_Sim_Pat=np.zeros(N_tree)

		#----------------------------------------------------------------------
		# Remove biomass of trees that died
		#----------------------------------------------------------------------

		# Remove biomass of trees that died directly from competition and
		# environmental conditions
		iKill=np.where(Pm_Sim_Reg[iLive]>=rLive)[0]
		Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
		tl['Csw_M_Reg'][iT,iLive[iKill]]=Csw_M
		tl['N_M_Reg'][iT,iLive[iKill]]=1
		tl['A'][iT:,iLive[iKill]]=np.nan
		tl['H'][iT:,iLive[iKill]]=np.nan
		tl['D'][iT+1:,iLive[iKill]]=np.nan
		tl['Csw'][iT:,iLive[iKill]]=np.nan
		rLive[iKill]=10; # Update rLive to avoid double-counting mortality

		# Remove biomass of trees that die due to insect attack
		iKill=np.where(Pm_Sim_Ins[iLive]>=rLive)[0]
		Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
		tl['Csw_M_Ins'][iT,iLive[iKill]]=Csw_M
		tl['N_M_Ins'][iT,iLive[iKill]]=1
		tl['A'][iT:,iLive[iKill]]=np.nan
		tl['H'][iT:,iLive[iKill]]=np.nan
		tl['D'][iT+1:,iLive[iKill]]=np.nan
		tl['Csw'][iT:,iLive[iKill]]=np.nan
		rLive[iKill]=10

		# Remove biomass of trees that die due to pathogen infection
		iKill=np.where(Pm_Sim_Pat[iLive]>=rLive)[0]
		Csw_M=tl['Csw'][iT,iLive[iKill]].copy()
		tl['Csw_M_Pat'][iT,iLive[iKill]]=Csw_M
		tl['N_M_Pat'][iT,iLive[iKill]]=1
		tl['A'][iT:,iLive[iKill]]=np.nan
		tl['H'][iT:,iLive[iKill]]=np.nan
		tl['D'][iT+1:,iLive[iKill]]=np.nan
		tl['Csw'][iT:,iLive[iKill]]=np.nan
		rLive[iKill]=10

		#----------------------------------------------------------------------
		# Mortality from inventory data sources
		#----------------------------------------------------------------------

		# Loop through events in year
		for iE in range(meta['Core']['Max Events Per Year']):

			# Event type IDs for the iE'th event of the year
			ID_Type=vi['EC']['ID Event Type'][iT,iS,iE].copy()

			# Total affected biomass carbon
			MF=vi['EC']['Mortality Factor'][iT,iS,iE].copy()

			if MF>0:

				iLive=np.where(np.isnan(tl['Csw'][iT,:])==False)[0]

				nKill=int(MF*iLive.size)

				if ID_Type==meta['LUT']['Event']['Wildfire']:
					# Wildfire-------------------------------------------------

					iKill=np.arange(0,nKill)
					tl['Csw_M_Fir'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
					tl['N_M_Fir'][iT,iLive[iKill]]=1
					vo['A'][iT,iS]=0

				elif (ID_Type==meta['LUT']['Event']['Harvest']) | (ID_Type==meta['LUT']['Event']['Harvest Salvage']):
					# Harvest--------------------------------------------------

					iKill=np.arange(0,nKill)
					tl['Csw_M_Har'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
					tl['N_M_Har'][iT,iLive[iKill]]=1
					vo['A'][iT,iS]=0

				elif (ID_Type==meta['LUT']['Event']['Sawtooth Commercial Thinning']):
					# Commercial thinning--------------------------------------

					idx_sorted=np.flip(np.argsort(tl['Csw'][iT,iLive]))
					iKill=idx_sorted[0:nKill]
					tl['Csw_M_Har'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
					tl['N_M_Har'][iT,iLive[iKill]]=1

				elif (ID_Type==meta['LUT']['Event']['Beetles']) | (ID_Type==meta['LUT']['Event']['Mountain Pine Beetle']) | (ID_Type==meta['LUT']['Event']['Defoliators']):
					# Insects--------------------------------------------------

					iKill=np.arange(0,nKill)
					tl['Csw_M_Ins'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
					tl['N_M_Ins'][iT,iLive[iKill]]=1

				elif (ID_Type==meta['LUT']['Event']['Sawtooth IDW']):
					# Sawtooth IDW---------------------------------------------

					Csw_interm=np.percentile(tl['Csw'][iT,iLive],33)
					idx_sorted=np.argsort( np.abs(tl['Csw'][iT,iLive]-Csw_interm) )
					iKill=idx_sorted[0:nKill]
					tl['Csw_M_Ins'][iT,iLive[iKill]]=tl['Csw'][iT,iLive[iKill]].copy()
					tl['N_M_Ins'][iT,iLive[iKill]]=1

				else:
					pass

				# Update tree variables
				tl['A'][iT:,iLive[iKill]]=0
				tl['H'][iT:,iLive[iKill]]=np.nan
				tl['D'][iT+1:,iLive[iKill]]=np.nan
				tl['Csw'][iT:,iLive[iKill]]=np.nan

	#--------------------------------------------------------------------------
	# Calculate biomass of other tissues from stemwood biomass
	#--------------------------------------------------------------------------

	# Calculate merch stemwood biomass
	iMerch=np.where(tl['D']>12.5)

	Csw_merch=np.zeros(tl['Csw'].shape); Csw_merch[iMerch]=tl['Csw'][iMerch]
	Csw_merch_M_Reg=np.zeros(tl['Csw'].shape); Csw_merch_M_Reg[iMerch]=tl['Csw_M_Reg'][iMerch]
	Csw_merch_M_Fir=np.zeros(tl['Csw'].shape); Csw_merch_M_Fir[iMerch]=tl['Csw_M_Fir'][iMerch]
	Csw_merch_M_Ins=np.zeros(tl['Csw'].shape); Csw_merch_M_Ins[iMerch]=tl['Csw_M_Ins'][iMerch]
	Csw_merch_M_Pat=np.zeros(tl['Csw'].shape); Csw_merch_M_Pat[iMerch]=tl['Csw_M_Pat'][iMerch]
	Csw_merch_M_Har=np.zeros(tl['Csw'].shape); Csw_merch_M_Har[iMerch]=tl['Csw_M_Har'][iMerch]

	# Standing biomass
	tl['Cswm'],tl['Cswnm'],tl['Cf'],tl['Cbk'],tl['Cbr'],tl['Crc'],tl['Crf'],tl['Ctot']=ExpandBiomass(tl['Csw'],Csw_merch)

	# Biomass loss due to regular mortality
	tl['Cswm_M_Reg'],tl['Cswnm_M_Reg'],tl['Cf_M_Reg'],tl['Cbk_M_Reg'],tl['Cbr_M_Reg'],tl['Crc_M_Reg'],tl['Crf_M_Reg'],tl['Ctot_M_Reg']=ExpandBiomass(tl['Csw_M_Reg'],Csw_merch_M_Reg)

	# Biomass loss due to wildfire
	tl['Cswm_M_Fir'],tl['Cswnm_M_Fir'],tl['Cf_M_Fir'],tl['Cbk_M_Fir'],tl['Cbr_M_Fir'],tl['Crc_M_Fir'],tl['Crf_M_Fir'],tl['Ctot_M_Fir']=ExpandBiomass(tl['Csw_M_Fir'],Csw_merch_M_Fir)

	# Biomass loss due to insects
	tl['Cswm_M_Ins'],tl['Cswnm_M_Ins'],tl['Cf_M_Ins'],tl['Cbk_M_Ins'],tl['Cbr_M_Ins'],tl['Crc_M_Ins'],tl['Crf_M_Ins'],tl['Ctot_M_Ins']=ExpandBiomass(tl['Csw_M_Ins'],Csw_merch_M_Ins)

	# Biomass loss due to pathogens
	tl['Cswm_M_Pat'],tl['Cswnm_M_Pat'],tl['Cf_M_Pat'],tl['Cbk_M_Pat'],tl['Cbr_M_Pat'],tl['Crc_M_Pat'],tl['Crf_M_Pat'],tl['Ctot_M_Pat']=ExpandBiomass(tl['Csw_M_Pat'],Csw_merch_M_Pat)

	# Biomass loss due to harvest
	tl['Cswm_M_Har'],tl['Cswnm_M_Har'],tl['Cf_M_Har'],tl['Cbk_M_Har'],tl['Cbr_M_Har'],tl['Crc_M_Har'],tl['Crf_M_Har'],tl['Ctot_M_Har']=ExpandBiomass(tl['Csw_M_Har'],Csw_merch_M_Har)

	# Biomass loss due to wind
	#tl['Cswm_M_Win'],tl['Cswnm_M_Win'],tl['Cf_M_Win'],tl['Cbk_M_Win'],tl['Cbr_M_Win'],tl['Crc_M_Win'],tl['Crf_M_Win'],tl['Ctot_M_Win']=ExpandBiomass(tl['Csw_M_Win'])

	#--------------------------------------------------------------------------
	# Biomass turnover (litterfall)
	#--------------------------------------------------------------------------

	# Foliage turnover
	tl['Cf_LF']=0.100*tl['Cf']

	# Branch turnover
	tl['Cbr_LF']=0.035*tl['Cbr']

	# Bark turnover
	tl['Cbk_LF']=0.035*tl['Cbk']

	# Fine root turnover
	tl['Crf_LF']=0.641*tl['Crf']

	# Coarse root turnover
	tl['Crc_LF']=0.02*tl['Crc']

	#--------------------------------------------------------------------------
	# Populate stand variables (individual-tree statistics)
	#--------------------------------------------------------------------------

	# Mean tree age
	vo['TreeMean_A'][:,iS]=np.nanmean(tl['A'],axis=1)

	# Mean tree diamter at breast height
	vo['TreeMean_D'][:,iS]=np.nanmean(tl['D'],axis=1)

	# Mean tree height (m)
	vo['TreeMean_H'][:,iS]=np.nanmean(tl['H'],axis=1)

	# Mean tree stemwood biomass (kg C)
	vo['TreeMean_Csw'][:,iS]=np.maximum(0,np.nanmean(tl['Csw'],axis=1))

	# Mean stemwood biomass growth (kg C yr-1)
	vo['TreeMean_Csw_G'][:,iS]=np.maximum(0,np.nanmean(tl['Csw_G'],axis=1))

	#--------------------------------------------------------------------------
	# Populate stand-level demographics
	#--------------------------------------------------------------------------

	# Stand density (stems ha-1)
	vo['N'][:,iS]=np.sum(~np.isnan(tl['Csw']),axis=1)

	# Stand density at beginning of time step
	N_time0=np.append(vo['N'][0,iS],vo['N'][0:-1,iS])

	# Demographic recruitment rate (% yr-1)
	vo['N_R'][:,iS]=np.minimum(100,np.sum(tl['N_R'],axis=1)/N_time0*100)

	# Demographic mortality rate (% yr-1)

	# Need to update N0 or else stand-replacing disturbances will have lower
	# relative mortality rate
	# *** The other solution is to move the prescribed disturbances above the
	# background mortality ***
	#N0=N0-Sum_N_M_Sim_Reg-Sum_N_M_Sim_Ins-Sum_N_M_Sim_Pat;

	N_M_Tot=tl['N_M_Reg']+tl['N_M_Fir']+tl['N_M_Ins']+tl['N_M_Pat']+tl['N_M_Har']+tl['N_M_Win']

	vo['N_M_Tot'][:,iS]=np.maximum(0,np.minimum(100,np.sum(N_M_Tot,axis=1)/N_time0*100))

	#	vo['N_M_Reg'][:,iS]=np.maximum(0,np.sum(tl['N_M_Reg'],axis=1)/N_time0*100)
	#	vo['N_M_Fir'][:,iS]=np.maximum(0,np.sum(tl['N_M_Fir'],axis=1)/N_time0*100)
	#	vo['N_M_Ins'][:,iS]=np.maximum(0,np.sum(tl['N_M_Ins'],axis=1)/N_time0*100)
	#	vo['N_M_Pat'][:,iS]=np.maximum(0,np.sum(tl['N_M_Pat'],axis=1)/N_time0*100)
	#	vo['N_M_Win'][:,iS]=np.maximum(0,np.sum(tl['N_M_Win'],axis=1)/N_time0*100)
	#	vo['N_M_Har'][:,iS]=np.maximum(0,np.sum(tl['N_M_Har'],axis=1)/N_time0*100)

	#--------------------------------------------------------------------------
	# Populate stand-level biomass
	#--------------------------------------------------------------------------

	# Gross growth of biomass (Mg C ha-1 yr-1)
	vo['C_G_Gross'][1:,iS,0]=np.nansum(np.maximum(0,np.diff(tl['Cswm'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,1]=np.nansum(np.maximum(0,np.diff(tl['Cswnm'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,2]=np.nansum(np.maximum(0,np.diff(tl['Cf'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,3]=np.nansum(np.maximum(0,np.diff(tl['Cbr'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,4]=np.nansum(np.maximum(0,np.diff(tl['Cbk'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,5]=np.nansum(np.maximum(0,np.diff(tl['Crc'],axis=0)),axis=1)/cm
	vo['C_G_Gross'][1:,iS,6]=np.nansum(np.maximum(0,np.diff(tl['Crf'],axis=0)),axis=1)/cm

	# Biomass turnover (Mg C ha-1 yr-1)
	vo['C_LF'][:,iS,2]=np.nansum(tl['Cf_LF'],axis=1)/cm
	vo['C_LF'][:,iS,3]=np.nansum(tl['Cbr_LF'],axis=1)/cm
	vo['C_LF'][:,iS,4]=np.nansum(tl['Cbk_LF'],axis=1)/cm
	vo['C_LF'][:,iS,5]=np.nansum(tl['Crc_LF'],axis=1)/cm
	vo['C_LF'][:,iS,6]=np.nansum(tl['Crf_LF'],axis=1)/cm

	# Biomass loss due to regular mortality (Mg C ha-1 yr-1)
	vo['C_M_Reg'][:,iS,0]=np.nansum(tl['Cswm_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,1]=np.nansum(tl['Cswnm_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,2]=np.nansum(tl['Cf_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,3]=np.nansum(tl['Cbr_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,4]=np.nansum(tl['Cbk_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,5]=np.nansum(tl['Crc_M_Reg'],axis=1)/cm
	vo['C_M_Reg'][:,iS,6]=np.nansum(tl['Crf_M_Reg'],axis=1)/cm

	# Biomass loss due to all mortality (Mg C ha-1 yr-1)
	vo['C_M_Tot'][:,iS,0]=np.nansum(tl['Cswm_M_Reg']+tl['Cswm_M_Fir']+tl['Cswm_M_Ins']+tl['Cswm_M_Pat']+tl['Cswm_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,1]=np.nansum(tl['Cswnm_M_Reg']+tl['Cswnm_M_Fir']+tl['Cswnm_M_Ins']+tl['Cswnm_M_Pat']+tl['Cswnm_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,2]=np.nansum(tl['Cf_M_Reg']+tl['Cf_M_Fir']+tl['Cswnm_M_Ins']+tl['Cswnm_M_Pat']+tl['Cswnm_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,3]=np.nansum(tl['Cbr_M_Reg']+tl['Cbr_M_Fir']+tl['Cbr_M_Ins']+tl['Cbr_M_Pat']+tl['Cbr_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,4]=np.nansum(tl['Cbk_M_Reg']+tl['Cbk_M_Fir']+tl['Cbk_M_Ins']+tl['Cbk_M_Pat']+tl['Cbk_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,5]=np.nansum(tl['Crc_M_Reg']+tl['Crc_M_Fir']+tl['Crc_M_Ins']+tl['Crc_M_Pat']+tl['Crc_M_Har'],axis=1)/cm
	vo['C_M_Tot'][:,iS,6]=np.nansum(tl['Crf_M_Reg']+tl['Crf_M_Fir']+tl['Crf_M_Ins']+tl['Crf_M_Pat']+tl['Crf_M_Har'],axis=1)/cm

	#	vo['C_M_Fir'][:,iS]=np.nansum(tl['Ctot_M_Fir'],axis=1)/cm
	#	vo['C_M_Ins'][:,iS]=np.nansum(tl['Ctot_M_Ins'],axis=1)/cm
	#	vo['C_M_Pat'][:,iS]=np.nansum(tl['Ctot_M_Pat'],axis=1)/cm
	#	vo['C_M_Win'][:,iS]=np.nansum(tl['Ctot_M_Win'],axis=1)/cm
	#	vo['C_M_Har'][:,iS]=np.nansum(tl['Ctot_M_Har'],axis=1)/cm

	# Net growth by tissue (Mg C ha-1 yr-1)
	vo['C_G_Net'][:,iS,0:7]=vo['C_G_Gross'][:,iS,0:7]-vo['C_M_Reg'][:,iS,0:7]

	# Ecosystem carbon pools by tissue (Mg C ha-1 yr-1)
	vo['C_Eco_Pools'][:,iS,0]=np.maximum(0,np.nansum(tl['Cswm'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,1]=np.maximum(0,np.nansum(tl['Cswnm'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,2]=np.maximum(0,np.nansum(tl['Cf'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,3]=np.maximum(0,np.nansum(tl['Cbr'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,4]=np.maximum(0,np.nansum(tl['Cbk'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,5]=np.maximum(0,np.nansum(tl['Crc'],axis=1))/cm
	vo['C_Eco_Pools'][:,iS,6]=np.maximum(0,np.nansum(tl['Crf'],axis=1))/cm

	#--------------------------------------------------------------------------
	# Populate stand-level merch stemwood volume
	#--------------------------------------------------------------------------

	# Merchantable stemwood volume of live trees
	vo['V_MerchLive'][:,iS]=np.maximum(0,np.nansum(tl['Cswm'],axis=1))/cm/meta['Param']['BEV']['Biophysical']['Ratio_Wood_C_to_DM']/0.45

	#t1=time.time()
	#print((t1-t0))

	return vo