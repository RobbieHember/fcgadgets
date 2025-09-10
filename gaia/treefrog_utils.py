#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.gaia.gaia_util as uga
import fcexplore.surface_climate.surfclim_util as usc
gp=gu.SetGraphics('Manuscript')

#%% Site conditions
def CompileSiteSpecs(meta):

	sL=['DF-H49','DF-H88','DF-H00','SB-Old','Mix-F77','Mix-F98','PB-Old','PB-H75','PB-H02']

	Site={}
	for s in sL:
		Site[s]={}

	# Canopy height (m)
	v='Height Canopy'
	s='DF-H00'; Site[s][v]=0.5
	s='DF-H88'; Site[s][v]=16
	s='DF-H49'; Site[s][v]=40
	s='SB-Old'; Site[s][v]=40
	s='Mix-F77'; Site[s][v]=40
	s='Mix-F98'; Site[s][v]=40
	s='PB-Old'; Site[s][v]=40
	s='PB-H75'; Site[s][v]=40
	s='PB-H02'; Site[s][v]=40

	meta['Site']=Site
	return meta


#%% Bowen ratio
def Calc_BowenRatio(meta,dEC):
	meta['Bowen Ratio']={}
	for site in dEC['TS'].keys():
		#site='DF-H49'
		meta['Bowen Ratio'][site]=np.zeros(12)
		for iM in range(12):
			meta['Bowen Ratio'][site][iM]=dEC['N'][site]['Heat flux sensible'][iM]/dEC['N'][site]['Heat flux latent'][iM]
	return meta

#%% Energy balance closure correction
def Calc_EnergyBalanceClosureCorrection(meta,dEC):
	meta['EBC Parameters']={}
	for site in dEC['TS'].keys():
		#site='DF-H49'

		x=dEC['N'][site]['Radiation net']
		y=dEC['N'][site]['Heat flux latent']+dEC['N'][site]['Heat flux sensible']
		rs,txt=gu.GetRegStats(x,y)
		meta['EBC Parameters'][site]=rs['B']

		y1=(y-rs['B'][0])/rs['B'][1]
		rs1,txt1=gu.GetRegStats(x,y1)

		dEC['N'][site]['Heat flux latent C']=np.zeros(12)
		dEC['N'][site]['Heat flux sensible C']=np.zeros(12)
		dEC['TS'][site]['Heat flux latent C']=dEC['TS'][site]['Heat flux latent'].copy()
		dEC['TS'][site]['Heat flux sensible C']=dEC['TS'][site]['Heat flux sensible'].copy()
		for iM in range(12):
			# Normals
			QL=dEC['N'][site]['Heat flux latent'][iM]
			QH=dEC['N'][site]['Heat flux sensible'][iM]
			Qtot=QL+QH
			Qtot_C=(Qtot-rs['B'][0])/rs['B'][1]
			deltaQtot=Qtot_C-Qtot
			fQL=1-meta['Bowen Ratio'][site][iM]/2
			fQH=meta['Bowen Ratio'][site][iM]/2
			dEC['N'][site]['Heat flux latent C'][iM]=QL+fQL*deltaQtot
			dEC['N'][site]['Heat flux sensible C'][iM]=QH+fQH*deltaQtot

			# Time series
			QL=dEC['TS'][site]['Heat flux latent'][iM::12]
			QH=dEC['TS'][site]['Heat flux sensible'][iM::12]
			Qtot=QL+QH
			Qtot_C=(Qtot-rs['B'][0])/rs['B'][1]
			deltaQtot=Qtot_C-Qtot
			fQL=1-meta['Bowen Ratio'][site][iM]/2
			fQH=meta['Bowen Ratio'][site][iM]/2
			dEC['TS'][site]['Heat flux latent C'][iM::12]=QL+fQL*deltaQtot
			dEC['TS'][site]['Heat flux sensible C'][iM::12]=QH+fQH*deltaQtot

		flg=0
		if flg==1:
			plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
			ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
			ax.plot(x,y,'o',ms=5,mec='w',mfc=[0.27,0.49,0.77],mew=0.75)
			ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',lw=0.5)
			ax.text(170,10,rs['txt'],fontsize=5,ha='right')
			ax.plot(rs1['xhat Line'],rs1['yhat Line'],'g--',lw=0.5)
			ax.text(160,160,'1:1',fontsize=5,ha='center')
			ax.set(xlabel='Net radiation (W m-2)',ylabel='Latent + sensible heat flux (W m-2)',xticks=np.arange(-40,300,20),yticks=np.arange(-40,300,20),xlim=[-20,180],ylim=[-20,180])
			ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
			plt.tight_layout()
			#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Name','png',900)
	return meta,dEC

#%%
def Plot_Ga_ModelSurface(meta):
	bwH=1; binH=np.arange(0,50+bwH,bwH)
	bwU=0.1; binU=np.arange(0,2.0+bwU,bwU)
	y=np.zeros((binH.size,binU.size))
	for i,H in enumerate(binH):
		for j,U in enumerate(binU):
			y[i,j]=uga.AerodynamicConductance(H,U)
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
	c=plt.contour(binU,binH,y)
	ax.set(xlabel='Wind speed (m s$^{-1}$)',ylabel='Vegetation height (m)')
	ax.clabel(c)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AerodynamicConductance','png',900)
	return
#%%
def Calc_Ga(meta,dEC):
	dGa={}
	for s in dEC['N'].keys():
		dEC['N'][s]['Ga']=uga.AerodynamicConductance(meta['Site'][s]['Height Canopy'],dEC['N'][s]['U'])
		dGa[s]=np.mean(dEC['N'][s]['Ga'])
	pd.DataFrame.from_dict(dGa).to_excel(meta['Paths']['Data'] + '\\Ga_Predicted.xlsx')
	return

#%%