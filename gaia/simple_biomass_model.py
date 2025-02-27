'''
SIMPLE BIOMASS MODEL
'''

#%% Import modules
import numpy as np
#from pymcmcstat.MCMC import MCMC
import fcgadgets.macgyver.util_general as gu
gp=gu.SetGraphics('Manuscript')
global par,con

#%%
def SimpleBiomassModel(reg,sts,con,dH):

	pr=sts[reg]['Param']

	d={}
	d['Year']=np.arange(con['Year_Start'],2023,1)
	for v in ['A','B','GG','M_Reg','M_ND','M_H','M_Tot','GN_Reg','GN_Tot','A_Harv','V_Harv']:
		d[v]=np.zeros((d['Year'].size,con['N Stand']))

	# Probability of harvest
	d['p_Harv']=np.zeros(d['Year'].size)
	ind1=np.where( (d['Year']>=dH['Year'][0]) & (d['Year']<=dH['Year'][-1]) )[0]
	ind2=np.where( (dH['Year']>=d['Year'][0]) & (dH['Year']<=d['Year'][-1]) )[0]
	d['p_Harv'][ind1]=np.nan_to_num( (dH['Data']['Area']['CC'][reg][ind2]*1000)/(sts[reg]['Area']*1e6) )
	#plt.plot(d['p_Harv'],'ob-')

	# Probability of MPB
	dMPB=sts[reg]['MPB Area']
	d['p_MPB']=np.zeros(d['Year'].size)
	if reg=='Subhumid continental':
		ind1=np.where( (d['Year']>=dMPB['Year'][0]) & (d['Year']<=dMPB['Year'][-1]) )[0]
		ind2=np.where( (dMPB['Year']>=d['Year'][0]) & (dMPB['Year']<=d['Year'][-1]) )[0]
		d['p_MPB'][ind1]=np.nan_to_num( (dMPB['Area'][ind2])/(sts[reg]['Area']*1e6) )

	# Spatial variability in growth (a multiplier ranging from 0 to 2)
	GF_Spatial=np.maximum(0.02,np.random.randint(200,size=con['N Stand']).astype('float')/100)

	for i,yr in enumerate(d['Year']):
		if i==0:
			continue

		# Update age
		d['A'][i,:]=d['A'][i-1,:]+1

		# Growth
		d['GG'][i,:]=(pr['G_B1']*(1+((pr['G_B2']*(d['A'][i,:]/pr['G_B3'])**pr['G_B4']-1)/np.exp(d['A'][i,:]/pr['G_B3']))))+(pr['G_B5']*d['A'][i,:])
		d['GG'][i,:]=pr['G_adj']*d['GG'][i,:]

		# Add spatial variability
		d['GG'][i,:]=GF_Spatial*d['GG'][i,:]

		# Adjust growth based on modifier
		if con['Growth Modifier Status']=='On':
			# Age dependent
			#m1=(pr['G_m1_young']+pr['G_m1_old'])-gu.Clamp(pr['G_m0']+(pr['G_Age_young']-d['A'][i,:])/(pr['G_Age_young']-pr['G_Age_old']),pr['G_m1_old'],pr['G_m1_young'])
			G_m1=1.0
			G_t1=2020
			if pr['G_m0']<=1:
				GF=gu.Clamp(pr['G_m0']+(G_m1-pr['G_m0'])/(G_t1-pr['G_t0'])*(yr-pr['G_t0']),pr['G_m0'],G_m1)
			else:
				GF=gu.Clamp(pr['G_m0']-(pr['G_m0']-G_m1)/(G_t1-pr['G_t0'])*(yr-pr['G_t0']),G_m1,pr['G_m0'])
			d['GG'][i,:]=GF*d['GG'][i,:]

		# Ensure growth > 0
		d['GG'][i,:]=np.maximum(0.05,d['GG'][i,:])

		# Add growth to biomass
		d['B'][i,:]=d['B'][i-1,:]+d['GG'][i,:]

		# Regular mortality
		BTR=pr['M_Intercept']+pr['M_BGC']+ \
			pr['M_B']*d['B'][i,:]+ \
			pr['M_B:BGC']*d['B'][i,:]+ \
			pr['M_A']*d['A'][i,:]+ \
			pr['M_OGP']*con['OGP']+ \
			pr['M_HI']*con['HI']+ \
			pr['M_DT']*con['DT']
		d['M_Reg'][i,:]=np.maximum(0,BTR/100*d['B'][i,:])

		# Modify regular mortality
		if con['Mortality Modifier Status']=='On':
			M_t1=2010
			MF=np.maximum(0,1-np.maximum(0,(yr-pr['M_Reg_t0'])/(M_t1-pr['M_Reg_t0'])))
			MF=pr['M_Reg_m1']+(1-pr['M_Reg_m1'])*MF
			d['M_Reg'][i,:]=np.maximum(0,MF*d['M_Reg'][i,:])

		# Remove mortality from biomass
		d['B'][i,:]=d['B'][i,:]-d['M_Reg'][i,:]

		# Regular net growth
		d['GN_Reg'][i,:]=d['GG'][i,:]-d['M_Reg'][i,:]

		# Natural disturbance mortality
		P_ND=pr['P_NatDist']
		rn_NatDist=np.random.random(con['N Stand'])
		if con['Mortality Modifier Status']=='On':
			M_t1=2010
			MF=np.maximum(0,1-np.maximum(0,(yr-pr['M_ND_t0'])/(M_t1-pr['M_ND_t0'])))
			MF=pr['M_ND_m1']+(1-pr['M_ND_m1'])*MF
			P_ND=MF*P_ND
		ind=np.where(rn_NatDist<P_ND)[0]
		d['M_ND'][i, ind ]=d['B'][i, ind ]
		d['A'][i, ind ]=0
		d['B'][i, ind ]=0

		# Mountain pine beetle
		if reg=='Subhumid continental':
			rn_MPB=np.random.random(con['N Stand'])
			ind=np.where(rn_MPB<d['p_MPB'][i])[0]
			d['M_ND'][i, ind ]=0.5*d['B'][i, ind ]
			d['A'][i, ind ]=0.5*d['A'][i, ind ]
			d['B'][i, ind ]=0.5*d['B'][i, ind ]

		# Harvest
		if con['Harvest Status']=='On':
			if (d['Year'][i]>=1901):
				# Confine harvesting to the most productive, stocked stands
				indH1=np.where( (GF_Spatial>1.0) & (d['B'][i,:]>45) & (rn_NatDist>pr['P_NatDist']) )[0]
				#indH1=np.where( (GF_Spatial>0) & (d['B'][i,:]>0) )[0]
				if indH1.size>0:
					# Adjust probability of harvest to match observations
					# Notes: This factor varies with N_Stand - not sure why!
					AEF_Harv=pr['Multi_Harv']*con['N Stand']/indH1.size
					p_Harv0=AEF_Harv*d['p_Harv'][i]
					rn_Harv=np.random.random(indH1.size)
					indH2=np.where( (rn_Harv<p_Harv0) )[0]
					B_Fel=d['B'][i, indH1[indH2] ]
					#d['M_H'][i, indH1[indH2] ]=d['M_H'][i, indH1[indH2] ]+B_Fel
					d['M_H'][i, indH1[indH2] ]=B_Fel
					d['A'][i, indH1[indH2] ]=0
					d['B'][i, indH1[indH2] ]=0
					d['A_Harv'][i, indH1[indH2] ]=1
					d['V_Harv'][i, indH1[indH2] ]=pr['FelledBiomassToSawlogRatio']*B_Fel/con['CarbonContent']/con['DensityWood']

		d['M_Tot'][i,:]=d['M_Reg'][i,:]+d['M_ND'][i,:]+d['M_H'][i,:]
		d['GN_Tot'][i,:]=d['GG'][i,:]-d['M_Tot'][i,:]

	return d
