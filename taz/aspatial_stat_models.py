'''
Aspatial Statistical Models of Disturbance Events
'''

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import beta
from fcgadgets.macgyver import util_general as gu
from fcgadgets.cbrunner import cbrun_util as cbu
import warnings
warnings.filterwarnings("ignore")

#%%
def PredictWildfire_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event):

	if (meta[pNam]['Year'][iT]>=meta[pNam]['Project']['Year Project']) & (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Scn ID']==-9999):
		# No future wildfire
		return vi

	# Occurrence
	sH='H' + str(meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Scn ID'])
	sF='F' + str(meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Scn ID'])

	if (meta[pNam]['Scenario'][iScn]['Wildfire Sim Pre-obs Scn ID']==-9999) & (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Scn ID']!=-9999):
		# The project is trying to include the future but it crashes because no history was set. 
		# Just hardwire the historical scenario - it will not actually be used.
		sH='H2'

	if (meta[pNam]['Scenario'][iScn]['Wildfire Sim Future Scn ID']==-9999):
		# If future is turned off, still need something so that pre-obs doesn't crash
		sF='F2'

	iT_wf=np.where(meta['Modules']['Taz']['wfss']['tv_scn']==meta[pNam]['Year'][iT])[0]
	if iT_wf.size==0:
		iT_wf=np.where(meta['Modules']['Taz']['wfss']['tv_scn']==meta['Modules']['Taz']['wfss']['tv_scn'][-1])[0]

	# Calculate annual probability of wildfire occurrence from BGC zone- and time-dependent Pareto functions of area burned
	Po=np.zeros(meta[pNam]['Project']['indBat'].size)
	for zone in vi['lsat']['bgcz idx'].keys():
		Po_Det=meta['Modules']['Taz']['wfss']['ByZone'][zone]['Po Det Scenarios'][sH][sF][iT_wf]
		beta=meta['Modules']['Taz']['wfss']['ByZone'][zone]['Beta_Pareto_Cal'].copy()
		b0=meta['Modules']['Taz']['wfss']['ByZone'][zone]['Pareto_scale_to_match_Po_mu'][0]
		b1=meta['Modules']['Taz']['wfss']['ByZone'][zone]['Pareto_scale_to_match_Po_mu'][1]
		Scale=b1*Po_Det+b0
		beta[1]=-np.abs(Scale)
		beta[2]=np.abs(Scale)
		N_t=1
		Po[vi['lsat']['bgcz idx'][zone]]=stats.pareto.rvs(beta[0],loc=beta[1],scale=beta[2],size=N_t)

	# Adjustments based on disturbance history
	mnam='Disturbance Effects on Wildfire Occurrence'
	Po1=Po.copy()

	# Effects of wildfire on probabilty of wildfire occurrence
	flg=0
	if flg==1:
		# Testing
		#a=0.07;b=1.0;t0=0;t1=84;t=np.arange(0,t1,1);plt.close('all'); plt.plot(t,a+(b-a)*(t/t1),'b.')
		adj=meta['Param']['BE']['WildfireDisturbanceEffects']['Wildfire Response']+(1.0-meta['Param']['BE']['WildfireDisturbanceEffects']['Wildfire Response'])*(meta['Modules'][mnam]['Wildfire Counter']/meta['Param']['BE']['WildfireDisturbanceEffects']['Wildfire Response Recovery Time'])
		ind=np.where(meta['Modules'][mnam]['Wildfire Flag']==1)[0]
		Po1[ind]=Po[ind]*adj[ind]

	# Effects of harvest on probabilty of wildfire occurrence
	flg=0
	if flg==1:
		adj=1+meta['Modules'][mnam]['Harvest Flag']*(meta['Param']['BE']['WildfireDisturbanceEffects']['Harvest Response']-1)
		Po1=np.maximum(Po1,Po*adj)

	# Effects of MPB on probabilty of wildfire occurrence
	flg=0
	if flg==1:
		ind=np.where(meta['Modules'][mnam]['Mountain Pine Beetle Flag']==1)[0]
		Po1[ind]=Po[ind]*meta['Param']['BE']['WildfireDisturbanceEffects']['Mountain Pine Beetle Response']

	# Effects of salvage harvest on probabilty of wildfire occurrence
	flg=0
	if flg==1:
		ind=np.where(meta['Modules'][mnam]['Salvage Harvest Flag']==1)[0]
		Po1[ind]=Po[ind]*meta['Param']['BE']['WildfireDisturbanceEffects']['Salvage Harvest Response']

	# Calculate occurrence based on random number generator
	rn=np.random.random(meta[pNam]['Project']['indBat'].size)
	iOccur=np.where(rn<Po1)[0]
	if iOccur.size>0:

		# Update simulated event tracking flag
		flag_sim_event[iOccur]=1

		# Severity
		rn=np.random.random(meta[pNam]['Project']['indBat'].size)
		p=meta['Param']['BE']['BurnSev']['PoC']
		Severity=np.zeros(rn.size)
		Severity[(rn<=p[0])]=meta['Param']['BE']['BurnSev']['M'][0]
		Severity[(rn>p[0]) & (rn<=p[1])]=meta['Param']['BE']['BurnSev']['M'][1]
		Severity[(rn>p[1]) & (rn<=p[2])]=meta['Param']['BE']['BurnSev']['M'][2]
		Severity[(rn>p[3])]=meta['Param']['BE']['BurnSev']['M'][3]

		for i in range(iOccur.size):
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,iOccur[i],:]==0)[0]
			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,iOccur[i],iE]=meta['LUT']['Event']['Wildfire']
				vi['EC']['Mortality Factor'][iT,iOccur[i],iE]=Severity[i]
				vi['EC']['ID Growth Curve'][iT,iOccur[i],iE]=1

	return vi,flag_sim_event

#%%
def PredictIBM_OnTheFly(meta,pNam,vi,iT,iScn,iEns,flag_sim_event):
	Po=meta['Modules']['Taz']['IBM Pareto']['Po'][iT,iEns]
	rn=np.random.random(meta[pNam]['Project']['indBat'].size)
	iOccur=np.where( (rn<Po) & (flag_sim_event==0) & (vi['lsat']['Insect Mortality Percent Tree Species Affected']['Mountain Pine Beetle']>0) )[0]
	#iOccur=np.where( (rn<Po) & (flag_sim_event==0) )[0]
	if iOccur.size>0:
		# Update simulated event tracking flag
		flag_sim_event[iOccur]=1
		for i in range(iOccur.size):
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,iOccur[i],:]==0)[0]
			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,iOccur[i],iE]=meta['LUT']['Event']['Mountain Pine Beetle']
				vi['EC']['Mortality Factor'][iT,iOccur[i],iE]=meta['Modules']['Taz']['IBM Pareto']['Mortality Fraction']
				vi['EC']['ID Growth Curve'][iT,iOccur[i],iE]=1
	return vi,flag_sim_event

#%%
def PredictHarvesting_OnTheFly(meta,pNam,vi,iT,iScn,iEns,V_Merch,Period):

	# This can crash in the last 2 time steps
	if vi['tv'][iT]>=vi['tv'][-2]:
		return vi

	# Indicator of THLB (THLB=1, Non-THLB=0)
	flag_thlb=vi['lsat']['THLB'][iT,:]

	# *** SPECIAL ORDER ***
	# # Indicator of energy production (don't harvest on the fly, it is pre-scheduled)
	# if meta[pNam]['Project']['Land Surface Class Dependent']!='No':
	#	 iT_lsc=np.where(vi['lsat']['LSC']['tv']==meta[pNam]['Year'][iT])[0]
	#	 if iT_lsc.size>0:
	#		 flag_ep=1*( (vi['lsat']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Fuel Break']) | (vi['lsat']['LSC']['Use'][iT_lsc,:]==meta['LUT']['LSC']['Use']['Energy Production']) )
	#	 else:
	#		 flag_ep=np.ones(flag_thlb.size,dtype=int)
	# else:
	#	 flag_ep=np.ones(flag_thlb.size,dtype=int)

	if Period=='Historical':

		# Inflection point
		Po_H_Inf=500

		# Shape parameter
		Po_H_Shape=meta['Param']['BEV']['ByBGCZ']['Harvest Po Shape']

		# Saturation value
		#bH=[0.0007,5.05,0.32,1975]
		#bH=[0.00085,5.15,0.32,1975]
		# 		if meta[pNam]['Year'][iT]<1933:
		# 			bH=[0.002,5.0,0.32,1977]
		# 			f1=bH[0]*np.maximum(0,(meta[pNam]['Year'][iT]-1800)/100)**bH[1]
		# 			f2=(1/(1+np.exp(bH[2]*(meta[pNam]['Year'][iT]-bH[3]))))
		# 		else:
		bH=[0.0007,5.15,0.32,1975]
		f1=bH[0]*np.maximum(0,(meta[pNam]['Year'][iT]-1800)/100)**bH[1]
		f2=(1/(1+np.exp(bH[2]*(meta[pNam]['Year'][iT]-bH[3]))))
		Po_H_Sat=f1*f2

		# Plot
		flg=0
		if flg==1:
			t=np.arange(1700,2001,1)
			f1=bH[0]*np.maximum(0,(t-1800)/100)**bH[1]
			f2=(1/(1+np.exp(bH[2]*(t-bH[3]))))
			Po_H_Sat=f1*f2
			plt.plot(t,Po_H_Sat,'g--',lw=1.5)

	else:

		# Future period

		if 'Po Harvest Inf' in meta[pNam]['Scenario'][iScn]:
			# Default has been overriden with project-specific value
			Po_H_Inf=meta[pNam]['Scenario'][iScn]['Po Harvest Inf']
		else:
			# Use BGC zone-specific defaults
			Po_H_Inf=meta['Param']['BEV']['ByBGCZ']['Harvest Po Inf']

		if 'Po Harvest Shape' in meta[pNam]['Scenario'][iScn]:
			# Default has been overriden with project-specific value
			Po_H_Shape=meta[pNam]['Scenario'][iScn]['Po Harvest Shape']
		else:
			# Use BGC zone-specific defaults
			Po_H_Shape=meta['Param']['BEV']['ByBGCZ']['Harvest Po Shape']

		if 'Po Harvest Sat Pct' in meta[pNam]['Scenario'][iScn]:
			# Default has been overriden with project-specific value
			Po_H_Sat=(meta[pNam]['Scenario'][iScn]['Po Harvest Sat Pct']/100)*vi['lsat']['Harvest Index']
		else:
			# Use BGC zone-specific defaults
			Po_H_Sat=(meta['Param']['BEV']['ByBGCZ']['Harvest Po Sat Pct']/100)*vi['lsat']['Harvest Index']

		# QA - Plot function:
		flg=0
		if flg==1:
			beta=[0.03,-0.02,400]
			V_Merch=np.arange(1,1200)
			Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
	
			plt.close('all')
			fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
			ax.plot(V_Merch,Po*100,'k-',linewidth=0.75,label='Harvest on-the-fly model 1')
			ax.set(position=[0.1,0.12,0.87,0.86],xlim=[0,800],xticks=np.arange(0,1300,100),xlabel='Merchantable volume (m$^3$ ha$^-$$^1$)', \
				   ylim=[0,5],ylabel='Annual probability of harvest (%)')
			ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
			ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
			gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\taz_ann_prob_harvest','png',500)

# 		if pNam=='BCFCS_NMF':
# 			# Plot example
# 			flg_eg=0
# 			if flg_eg==1:
# 				tv=np.arange(1901,2100,1)
# 				Po_Sat=5
# 				Rate=0.01*Po_Sat
# 				Po=np.maximum(0,Po_Sat+Rate*np.maximum(1,tv-2021))
# 				plt.close('all'); plt.plot(tv,Po,'b-',lw=1.5)
# 				Po_Sat=5
# 				Rate=0.02*Po_Sat
# 				Po=np.maximum(0,Po_Sat+Rate*np.maximum(1,tv-2021))
# 				plt.plot(tv,Po,'c--',lw=1.5)
# 	
# 			if np.isin(iScn,meta[pNam]['Project']['Actual Indices'])==True:
# 				Po_Sat=5
# 				Rate=0.01*Po_Sat
# 				Po_H_Sat=0.01*np.maximum(0,Po_Sat+Rate*np.maximum(1,vi['tv'][iT]-2021))
# 			else:
# 				Po_Sat=5
# 				Rate=0.02*Po_Sat
# 				Po_H_Sat=0.01*np.maximum(0,Po_Sat+Rate*np.maximum(1,vi['tv'][iT]-2021))

	# Annual probability of occurrence
	Po=Po_H_Sat*(1/(1+np.exp(-Po_H_Shape*(V_Merch-Po_H_Inf))))

	# Don't allow harvesting below a specific volume
	#Po[V_Merch<150]=0

	# Random number
	rn=np.random.random(V_Merch.size)

	# Don't allow historical harvest simulations to occur on
	# the footprint of recorded historical harvesting
	if Period=='Historical':
		flag_HistHarv=-1*(np.minimum(1,vi['lsat']['Year Harvest First'])-1)
	else:
		flag_HistHarv=np.ones(rn.size)

	# Don't allow harvesting of stands that were recently fertilized
	flag_NA=np.ones(rn.size)
	#if Period=='Future':
	#	flag_NA=(meta['Modules']['NutrientApp']['ResponseCounterContinuous']==0) | (meta['Modules']['NutrientApp']['ResponseCounterContinuous']>meta['Param']['BEV']['NutrientApp']['HarvestRestrictTSNA']).astype(int)
	if 'HarvestRestrictTSNA' in meta[pNam]['Scenario'][iScn]:
		if Period=='Future':
			flag_NA=(meta['Modules']['NutrientApp']['ResponseCounterContinuous']==0) | (meta['Modules']['NutrientApp']['ResponseCounterContinuous']>meta[pNam]['Scenario'][iScn]['HarvestRestrictTSNA']).astype(int)

	# Occurrence
	#Oc=flag_ep*flag_thlb*flag_HistHarv*np.floor(np.minimum(1,Po/rn))
	#Oc=flag_ep*flag_thlb*np.floor(np.minimum(1,Po/rn))
	Oc=flag_thlb*flag_NA*np.floor(np.minimum(1,Po/rn))

	# Index to occurrence
	indS=np.where(Oc==1)[0]

	if indS.size>0:
		for i in range(indS.size):
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,indS[i],:]==0)[0]
			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Harvest']
				vi['EC']['Mortality Factor'][iT,indS[i],iE]=1.0
				vi['EC']['ID Growth Curve'][iT,indS[i],iE]=1
				iE=iAvailable[0]+1 # changing this to zero will cause the harvest to be overwritten
				vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Pile Burn']
				vi['EC']['Mortality Factor'][iT,indS[i],iE]=0
				vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2
				iE=iAvailable[0]+2
				vi['EC']['ID Event Type'][iT,indS[i],iE]=meta['LUT']['Event']['Planting']
				vi['EC']['Mortality Factor'][iT,indS[i],iE]=0
				vi['EC']['ID Growth Curve'][iT,indS[i],iE]=2

	return vi

#%%
def PredictDisease_OnTheFly(meta,pNam,vi,iT,iScn,iEns,Age,flag_sim_event):

	# Occurrence
	flg=0
	if flg==1:
		beta=[-0.035,100]
		Age=np.arange(1,500)
		Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
		fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
		ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of disease')
		ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

	if 'Disease Occ Sat Pct' in meta[pNam]['Scenario'][iScn].keys():
		# Use custom values specified for each scenario
		a1=meta[pNam]['Scenario'][iScn]['Disease Occ Sat Pct']
		a2=meta[pNam]['Scenario'][iScn]['Disease Occ Shape']
		a3=meta[pNam]['Scenario'][iScn]['Disease Occ Inf']
		Po=(a1/100)*(1/(1+np.exp(-a2*(Age-a3))))
	else:
		# Use BGC zone values
		Po=(meta['Param']['BEV']['ByBGCZ']['Disease Occ Sat Pct']/100)*(1/(1+np.exp(-meta['Param']['BEV']['ByBGCZ']['Disease Occ Shape']*(Age-meta['Param']['BEV']['ByBGCZ']['Disease Occ Inf']))))

	# Control disturbances before harvest
# 	if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
# 		yr=meta[pNam]['Year'][iT]
# 		yrH=vi['lsat']['Year Harvest First']
# 		dY=yr-yrH
# 		ind=np.where( (dY>-300) & (dY<=0) )[0]
# 		Po[ind]=0

	# Don't simulate where there has already been a major disturbance event
	#Po[flag_sim_event==1]=0

	rn=np.random.random(Age.size)
	iOccur=np.where(rn<Po)[0]

	if iOccur.size>0:

		# Update simulated event tracking flag
		flag_sim_event[iOccur]=1

		# Populate
		for iA in iOccur:
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]

			# Severity
			rn=np.random.random(1)[0]
			P=[meta['Param']['BEV']['ByBGCZ']['Disease Sev Prob Low'][iA],meta['Param']['BEV']['ByBGCZ']['Disease Sev Prob Medium'][iA],meta['Param']['BEV']['ByBGCZ']['Disease Sev Prob High'][iA]]
			if (rn<P[0]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Disease Sev Mort Low'][iA]
			elif (rn>=P[0]) & (rn<P[0]+P[1]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Disease Sev Mort Medium'][iA]
			elif (rn>=P[0]+P[1]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Disease Sev Mort High'][iA]
			else:
				Severity=0.0

			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Disease Root']
				vi['EC']['Mortality Factor'][iT,iA,iE]=Severity
				vi['EC']['ID Growth Curve'][iT,iA,iE]=1
			else:
				print('No space left in event chronology for on-the-fly event!')

	return vi,flag_sim_event

#%%
def PredictWind_OnTheFly(meta,pNam,vi,iT,iEns,Age,flag_sim_event):

	# Occurrence
	flg=0
	if flg==1:
		beta=[-0.035,150,0.005]
		Age=np.arange(1,500)
		Po=beta[2]*(1/(1+np.exp(beta[0]*(Age-beta[1]))))
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
		ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of breakup')
		ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

	Po=(meta['Param']['BEV']['ByBGCZ']['Wind Occ Sat Pct']/100)*(1/(1+np.exp(-meta['Param']['BEV']['ByBGCZ']['Wind Occ Shape']*(Age-meta['Param']['BEV']['ByBGCZ']['Wind Occ Inf']))))

	# Control occurrence before harvest
# 	if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
# 		yr=meta[pNam]['Year'][iT]
# 		yrH=vi['lsat']['Year Harvest First']
# 		dY=yr-yrH
# 		ind=np.where( (dY>-300) & (dY<=0) )[0]
# 		Po[ind]=0

	# Don't simulate where there has already been a major disturbance event
	#Po[flag_sim_event==1]=0

	rn=np.random.random(Age.size)
	iOccur=np.where(rn<Po)[0]

	if iOccur.size>0:
		# Populate
		for iA in iOccur:
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]

			# Severity
			rn=np.random.random(1)[0]
			P=[meta['Param']['BEV']['ByBGCZ']['Wind Sev Prob Low'][iA],meta['Param']['BEV']['ByBGCZ']['Wind Sev Prob Medium'][iA],meta['Param']['BEV']['ByBGCZ']['Wind Sev Prob High'][iA]]
			if (rn<P[0]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Wind Sev Mort Low'][iA]
			elif (rn>=P[0]) & (rn<P[0]+P[1]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Wind Sev Mort Medium'][iA]
			elif (rn>=P[0]+P[1]):
				Severity=meta['Param']['BEV']['ByBGCZ']['Wind Sev Mort High'][iA]
			else:
				Severity=0.0

			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Wind']
				vi['EC']['Mortality Factor'][iT,iA,iE]=Severity
				vi['EC']['ID Growth Curve'][iT,iA,iE]=1
			else:
				print('No space left in event chronology for on-the-fly event!')

	return vi,flag_sim_event

#%%
def PredictFrost_OnTheFly(meta,pNam,vi,iT,iEns,Age,flag_sim_event):

	# Occurrence
	flg=0
	if flg==1:
		beta=[0.5,5]
		Age=np.arange(1,500)
		Po=1/(1+np.exp(beta[0]*(Age-beta[1])))
		fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
		ax.plot(Age,Po,'k-',linewidth=0.75,label='Default model')
		ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of disease')
		ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

	Po=(meta['Param']['BEV']['ByBGCZ']['Frost Po Sat Pct']/100)*(1/(1+np.exp(meta['Param']['BEV']['ByBGCZ']['Frost Po Shape']*(Age-meta['Param']['BEV']['ByBGCZ']['Frost Po Inf']))))

	# Don't simulate where there has already been a major disturbance event
	Po[flag_sim_event==1]=0

	rn=np.random.random(Age.size)
	iOccur=np.where(rn<Po)[0]

	if iOccur.size>0:

		# Control disturbances before harvest
		if meta[pNam]['Project']['Scenario Source']!='Spreadsheet':
			yr=meta[pNam]['Year'][iT]
			yrH=vi['lsat']['Year Harvest First']
			dY=yr-yrH
			ind=np.where( (dY>-300) & (dY<=0) )[0]
			Po[ind]=0

		# Severity
		# Severity (Uniform distribution)
		#Severity=np.random.random(1)[0] # Uniform distribution

		# Severity (Beta distribution)
		Severity=meta['Param']['BEV']['ByBGCZ']['Frost Severity']
		#Severity=np.random.beta(meta['Param']['BEV']['ByBGCZ']['Frost Sev Alpha'],meta['Param']['BEV']['ByBGCZ']['Frost Sev Beta'],size=1)[0]

		# Populate
		for iA in iOccur:
			iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]
			if iAvailable.size>0:
				iE=iAvailable[0]+0
				vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Frost Snow Ice Hail']
				vi['EC']['Mortality Factor'][iT,iA,iE]=Severity[iA]
				vi['EC']['ID Growth Curve'][iT,iA,iE]=1
			else:
				print('No space left in event chronology for on-the-fly event!')

	return vi,flag_sim_event

#%% Modelling severity of wind and disease
# Import modules
flg=0
if flg==1:
	x=np.linspace(0,1,100)
	plt.close('all'); fig,ax=plt.subplots(1)
	a=2;b=8; ax.plot(x,beta.pdf(x,a,b),'b-',lw=1,label='beta pdf')
	a=4.7;b=3.3; ax.plot(x,beta.pdf(x,a,b),'g-',lw=1,label='beta pdf')
	a=8;b=2; ax.plot(x,beta.pdf(x,a,b),'r-',lw=1,label='beta pdf')

# 	a=4;b=4; yhat=beta.pdf(x,a,b);
# 	ind=np.where(np.abs(yhat-np.mode(yhat))==np.min(np.abs(yhat-np.mode(yhat))))[0]
# 	print(x[ind])

	r=np.random.beta(a,b,size=1000)
	plt.hist(r)

