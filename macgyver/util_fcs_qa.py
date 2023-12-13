
#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.patches import Rectangle
import gc as garc
import copy
import statsmodels.formula.api as smf
import warnings
import time
import copy
from scipy import stats
import matplotlib.colors
import matplotlib.ticker as ticker
#from matplotlib import animation
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcexplore.field_plots.Processing.psp_util as ugp

#%% Age class distribution (by BGC Zone)
def QA_FullComparisonAgeDistByBGC_CN(meta,gpt):
	x=np.arange(0,501,1); yt=np.arange(0,2,0.1)
	ord=np.flip(np.argsort(meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)']))
	lab=np.array(['' for _ in range(ord.size)],dtype=object)
	plt.close('all'); fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(22,11)); cnt=0
	for j in range(3):
		for i in range(3):
			zone=meta['Param']['BE']['BGC Zone Averages']['Name'][ord[cnt]]
			#ind=np.where( (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) & (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
			#kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
			#p=kde(x); y1=p/np.sum(p)*100
			#ax[i,j].plot(x,y1,'k-',lw=0.75,color=[0.5,0,1],label='VRI ground plots')
	
			ind=np.where( (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) & (gpt['PTF CN']==1) & (gpt['Age Mean t0']>=0) )[0]
			kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
			p=kde(x); y2=p/np.sum(p)*100
			ax[i,j].fill_between(x,y2,alpha=0.15)
			ax[i,j].plot(x,y2,'b-',lw=0.75,color=[0.27,0.44,0.79],label='CMI + NFI network\ntree cores')
			#ym=np.maximum(np.max(y1),np.max(y2)); indYT=np.where(yt>ym)[0]
			ym=np.max(y2); indYT=np.where(yt>ym)[0]
			mu=np.mean(gpt['Age Mean t0'][ind])
			ax[i,j].plot([mu,mu],[0,20],'b--',color=[0.27,0.44,0.79],label='Mean')
			ax[i,j].set(ylabel='Frequency (%)',xlabel='Stand age, years',xlim=[0,500],ylim=[0,yt[indYT[1]]])
			if (i==0) & (j==0):
				ax[i,j].legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			plt.tight_layout()
			lab[cnt]=zone
			cnt=cnt+1
	gu.axletters(ax,plt,0.025,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Labels=lab,LabelSpacer=0.035)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareAgeDistByBGC_CN_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullCompareBiomassDynamicsAve_CN(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):

		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()

		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]

		dO_mu={}
		dO_se={}
		for v in dObs['data']['CNV'].keys():
			dO_mu[v]=np.nansum(Area*dObs['data']['CNV'][v]['mu'])/np.nansum(Area)
			dO_se[v]=np.nansum(Area*dObs['data']['CNV'][v]['se'])/np.nansum(Area)
		dO_mu['Ctot G Tot']=dO_mu['Ctot G Surv']+dO_mu['Ctot G Recr']
		dO_se['Ctot G Tot']=dO_se['Ctot G Surv']+dO_se['Ctot G Recr']

		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		dM_mu={}
		dM_se={}
		for v in dMod[0]['SBS'].keys():
			mu=np.zeros(lab.size)
			se=np.zeros(lab.size)
			for i in range(dObs['code'].size):
				tmp=dMod[iScn][dObs['code'][i]]
				mu[i]=np.mean(tmp[v][iT])
				se[i]=2*np.std(tmp[v][iT])/np.sqrt(iT.size)
			dM_mu[v]=np.sum(Area*mu)/np.sum(Area)
			dM_se[v]=np.sum(Area*se)/np.sum(Area)

		lab=['Gross\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #,'Harvest\nmortality'
		cl=meta['Graphics']['GP Comp']; barw=0.38
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(9,7))
		ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
		ax.bar(1-barw/2-0.01,dO_mu['Ctot G Tot'],barw,facecolor=cl['bl'],label='Ground plot observations')
		ax.bar(1+barw/2+0.01,dM_mu['C_G_Gross_Tot'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
		ax.errorbar(1-barw/2-0.01,dO_mu['Ctot G Tot'],yerr=dO_se['Ctot G Tot'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(1+barw/2+0.01,dM_mu['C_G_Gross_Tot'],yerr=dM_se['C_G_Gross_Tot'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		
		ax.bar(2-barw/2-0.01,-dO_mu['Ctot Mort Nat'],barw,facecolor=cl['bl'])
		ax.bar(2+barw/2+0.01,-dM_mu['C_M_Nat'],barw,facecolor=cl['gl'])
		ax.errorbar(2-barw/2-0.01,-dO_mu['Ctot Mort Nat'],yerr=dO_se['Ctot Mort Nat'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(2+barw/2+0.01,-dM_mu['C_M_Nat'],yerr=dM_se['C_M_Nat'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		
		ax.bar(3-barw/2-0.01,-dO_mu['Ctot Mort Harv'],barw,facecolor=cl['bl'])
		ax.bar(3+barw/2+0.01,-dM_mu['C_M_Harv'],barw,facecolor=cl['gl'])
		ax.errorbar(3-barw/2-0.01,-dO_mu['Ctot Mort Harv'],yerr=dO_se['Ctot Mort Harv'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(3+barw/2+0.01,-dM_mu['C_M_Harv'],yerr=dM_se['C_M_Harv'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		
		ax.bar(4-barw/2-0.01,dO_mu['Ctot Net'],barw,facecolor=cl['bl'])
		ax.bar(4+barw/2+0.01,dM_mu['C_G_Net_Tot'],barw,facecolor=cl['gl'])
		ax.errorbar(4-barw/2-0.01,dO_mu['Ctot Net'],yerr=dO_se['Ctot Net'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(4+barw/2+0.01,dM_mu['C_G_Net_Tot'],yerr=dM_se['C_G_Net_Tot'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		
		ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-2,3,0.5),ylabel='Carbon balance of trees (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-1.5,2])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareBiomassDynamicsAve_CN_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullCompareBiomassByBGC_CNV(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()

		d={}
		d['obs mu']=dObs['data']['CNV']['Ctot L t0']['mu']
		d['obs se']=dObs['data']['CNV']['Ctot L t0']['se']
		# Add modelled data
		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		d['mod mu']=np.zeros(d['obs mu'].size)
		d['mod se']=np.zeros(d['obs se'].size)
		for i in range(dObs['code'].size):
			d['mod mu'][i]=np.mean(dMod[iScn][dObs['code'][i]]['C_Biomass_Tot'][iT])
			d['mod se'][i]=2*np.std(dMod[iScn][dObs['code'][i]]['C_Biomass_Tot'][iT])/np.sqrt(iT.size)
		
		# Put in order
		ord=np.argsort(d['obs mu'])
		lab=np.flip(lab[ord])
		u=u[ord]
		for v in d:
			d[v]=np.flip(d[v][ord])
		
		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
		
		# Area weighting
		for v in d:
			d[v]=np.append(d[v],np.nansum(d[v]*Area)/np.nansum(Area))
		lab=np.append(lab,'Weighted\naverage')
		u=np.append(u,0.0)
		
		# Percent difference
		yp=d['mod mu']
		yo=d['obs mu']
		Dp=(yp-yo)/yo*100
		
		cl=meta['Graphics']['GP Comp']; barw=0.32
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,6))
		ax.bar(np.arange(u.size)-barw/2-0.01,d['obs mu'],barw,facecolor=cl['bl'],label='Ground plots')
		ax.bar(np.arange(u.size)+barw/2+0.01,d['mod mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
		ax.errorbar(np.arange(u.size)-barw/2-0.01,d['obs mu'],yerr=d['obs se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(np.arange(u.size)+barw/2+0.01,d['mod mu'],yerr=d['mod se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		for i in range(u.size):
			if Dp[i]>=0:
				a='+'
			else:
				a=''
			ax.text(i+barw/2+0.01,yp[i]+d['mod se'][i]+6,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
		#for i in range(u.size):
		#ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
		ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
		plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareBiomassByBGC_CNV_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullCompareGrowthGrossByBGC_CN(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()
		d={}
		d['obs mu']=dObs['data']['CN']['Ctot G Surv']['mu']+dObs['data']['CN']['Ctot G Recr']['mu']
		d['obs se']=dObs['data']['CN']['Ctot G Surv']['se']+dObs['data']['CN']['Ctot G Recr']['se']
		
		# Add modelled data
		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		d['mod mu']=np.zeros(d['obs mu'].size)
		d['mod se']=np.zeros(d['obs se'].size)
		for i in range(dObs['code'].size):
			d['mod mu'][i]=np.mean(dMod[iScn][dObs['code'][i]]['C_G_Gross_Tot'][iT])
			d['mod se'][i]=2*np.std(dMod[iScn][dObs['code'][i]]['C_G_Gross_Tot'][iT])/np.sqrt(iT.size)
		
		# Put in order
		ord=np.argsort(np.nan_to_num(d['obs mu']))
		lab=np.flip(lab[ord])
		u=u[ord]
		for v in d:
			d[v]=np.flip(d[v][ord])
		
		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
		
		# Area weighting
		for v in d:
			d[v]=np.append(d[v],np.nansum(d[v]*Area)/np.nansum(Area))
		lab=np.append(lab,'Weighted\naverage')
		u=np.append(u,0.0)
		
		# Percent difference
		yp=d['mod mu']
		yo=d['obs mu']
		Dp=(yp-yo)/yo*100
		
		cl=meta['Graphics']['GP Comp']; barw=0.32
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,6))
		ax.bar(np.arange(u.size)-barw/2-0.01,d['obs mu'],barw,facecolor=cl['bl'],label='Ground plots')
		ax.bar(np.arange(u.size)+barw/2+0.01,d['mod mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
		ax.errorbar(np.arange(u.size)-barw/2-0.01,d['obs mu'],yerr=d['obs se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(np.arange(u.size)+barw/2+0.01,d['mod mu'],yerr=d['mod se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		for i in range(u.size):
			if (np.abs(Dp[i])>10000) | (np.isnan(yo[i])==True):
				continue
			if Dp[i]>=0:
				a='+'
			else:
				a=''
			ax.text(i+barw/2+0.01,yp[i]+d['mod se'][i]+0.25,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
		#for i in range(u.size):
		#ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
		ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Gross growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,5])
		plt.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareGrowthGrossByBGC_CNV_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullCompareGrowthNetByBGC_CN(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()
		d={}
		d['obs mu']=dObs['data']['CN']['Ctot Net']['mu']
		d['obs se']=dObs['data']['CN']['Ctot Net']['se']
		
		# Add modelled data
		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		d['mod mu']=np.zeros(d['obs mu'].size)
		d['mod se']=np.zeros(d['obs se'].size)
		for i in range(dObs['code'].size):
			d['mod mu'][i]=np.mean(dMod[iScn][dObs['code'][i]]['C_G_Net_Tot'][iT])
			d['mod se'][i]=2*np.std(dMod[iScn][dObs['code'][i]]['C_G_Net_Tot'][iT])/np.sqrt(iT.size)
		
		# Put in order
		tmp=d['obs mu'].copy()
		ind=np.where(np.isnan(tmp)==True)[0]
		tmp[ind]=-1000
		ord=np.argsort(tmp)
		lab=np.flip(lab[ord])
		u=u[ord]
		for v in d.keys():
			d[v]=np.flip(d[v][ord])
		
		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
		
		# Area weighting
		for v in d:
			d[v]=np.append(d[v],np.nansum(d[v]*Area)/np.nansum(Area))
		lab=np.append(lab,'Weighted\naverage')
		u=np.append(u,0.0)
		
		# Percent difference
		yp=d['mod mu']
		yo=d['obs mu']
		Dp=(yp-yo)/yo*100
		
		cl=meta['Graphics']['GP Comp']; barw=0.32
		
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,6))
		ax.plot([-1,100],[0,0],'k-',lw=1)
		ax.bar(np.arange(u.size)-barw/2-0.01,d['obs mu'],barw,facecolor=cl['bl'],label='Ground plots')
		ax.bar(np.arange(u.size)+barw/2+0.01,d['mod mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
		ax.errorbar(np.arange(u.size)-barw/2-0.01,d['obs mu'],yerr=d['obs se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(np.arange(u.size)+barw/2+0.01,d['mod mu'],yerr=d['mod se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		for i in range(u.size):
			if (np.abs(Dp[i])>10000) | (np.isnan(yo[i])==True):
				continue
			if Dp[i]>=0:
				a='+'
			else:
				a=''
			ax.text(i+barw/2+0.01,yp[i]+d['mod se'][i]+0.25,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
		#for i in range(u.size):
		#ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
		ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[-3,3])
		plt.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareGrowthNetByBGC_CN_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullCompareMortalityByBGC_CN(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()
		d={}
		d['obs mu']=dObs['data']['CN']['Ctot Mort Nat']['mu']+dObs['data']['CN']['Ctot Mort Harv']['mu']
		d['obs se']=dObs['data']['CN']['Ctot Mort Nat']['se']+dObs['data']['CN']['Ctot Mort Harv']['se']
		
		# Add modelled data
		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		d['mod mu']=np.zeros(d['obs mu'].size)
		d['mod se']=np.zeros(d['obs se'].size)
		for i in range(dObs['code'].size):
			d['mod mu'][i]=np.mean(dMod[iScn][dObs['code'][i]]['C_M_Tot'][iT])
			d['mod se'][i]=2*np.std(dMod[iScn][dObs['code'][i]]['C_M_Tot'][iT])/np.sqrt(iT.size)
		
		# Put in order
		ord=np.argsort(np.nan_to_num(d['obs mu']))
		lab=np.flip(lab[ord])
		u=u[ord]
		for v in d:
			d[v]=np.flip(d[v][ord])
		
		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
		
		# Area weighting
		for v in d:
			d[v]=np.append(d[v],np.nansum(d[v]*Area)/np.nansum(Area))
		lab=np.append(lab,'Weighted\naverage')
		u=np.append(u,0.0)
		
		# Percent difference
		yp=d['mod mu']
		yo=d['obs mu']
		Dp=(yp-yo)/yo*100
		
		cl=meta['Graphics']['GP Comp']; barw=0.32
		plt.close('all');fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,6))
		ax.bar(np.arange(u.size)-barw/2-0.01,d['obs mu'],barw,facecolor=cl['bl'],label='Ground plots')
		ax.bar(np.arange(u.size)+barw/2+0.01,d['mod mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
		ax.errorbar(np.arange(u.size)-barw/2-0.01,d['obs mu'],yerr=d['obs se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(np.arange(u.size)+barw/2+0.01,d['mod mu'],yerr=d['mod se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		for i in range(u.size):
			if (np.abs(Dp[i])>10000) | (np.isnan(yo[i])==True):
				continue
			if Dp[i]>=0:
				a='+'
			else:
				a=''
			ax.text(i+barw/2+0.01,yp[i]+d['mod se'][i]+0.25,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
		#for i in range(u.size):
		#ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
		ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Mortality (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,6])
		plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.2)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareMortByBGC_CN_' + str(iScn+1),'png',900)
	
	return

#%%
def QA_FullCompareSOCByBGC(meta,pNam,tv,dObs0,dMod0):
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		dObs=copy.deepcopy(dObs0)
		dMod=copy.deepcopy(dMod0)
		lab=dObs['code'].copy()
		u=dObs['id'].copy()

		# Add modelled data
		iT=np.where( (tv>=2000) & (tv<=2018) )[0]
		dM={}
		dM['Tot_mu']=np.zeros(lab.size)
		dM['Org_mu']=np.zeros(lab.size)
		dM['Min_mu']=np.zeros(lab.size)
		dM['Tot_se']=np.zeros(lab.size)
		dM['Org_se']=np.zeros(lab.size)
		dM['Min_se']=np.zeros(lab.size)
		for i in range(lab.size):
			dM['Tot_mu'][i]=np.mean(dMod[iScn][lab[i]]['C_Soil_Tot'][iT])
			dM['Org_mu'][i]=np.mean(dMod[iScn][lab[i]]['C_Soil_OHorizon'][iT])
			dM['Min_mu'][i]=dM['Tot_mu'][i]-dM['Org_mu'][i]
			dM['Tot_se'][i]=2*np.std(dMod[iScn][lab[i]]['C_Soil_Tot'][iT])/np.sqrt(iT.size)
			dM['Org_se'][i]=2*np.std(dMod[iScn][lab[i]]['C_Soil_OHorizon'][iT])/np.sqrt(iT.size)
			dM['Min_se'][i]=dM['Tot_se'][i]-dM['Org_se'][i]

		# Put in order
		ord=np.argsort(dObs['data']['Soil']['TOT_C_THA']['mu'])
		lab=np.flip(lab[ord])
		for v in dObs['data']['Soil'].keys():
			for k in dObs['data']['Soil'][v].keys():
				dObs['data']['Soil'][v][k]=np.flip(dObs['data']['Soil'][v][k][ord])
		for v in dM.keys():
			dM[v]=np.flip(dM[v][ord])

		Area=np.zeros(lab.size)
		for i in range(lab.size):
			ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
			Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
		
		# Area weighting
		for v in dObs['data']['Soil'].keys():
			dObs['data']['Soil'][v]['mu']=np.append(dObs['data']['Soil'][v]['mu'],np.sum(dObs['data']['Soil'][v]['mu']*Area)/np.sum(Area))
			dObs['data']['Soil'][v]['se']=np.append(dObs['data']['Soil'][v]['se'],np.sum(dObs['data']['Soil'][v]['se']*Area)/np.sum(Area))
		for v in dM.keys():
			dM[v]=np.append(dM[v],np.sum(dM[v]*Area)/np.sum(Area))
		lab=np.append(lab,'Weighted\naverage')
		u=np.append(u,0.0)
		
		# Percent difference
		yp=dM['Tot_mu']
		yo=dObs['data']['Soil']['TOT_C_THA']['mu']
		Dp=(yp-yo)/yo*100
		
		cl=meta['Graphics']['GP Comp']; barw=0.32
		plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,6))
		ax.bar(np.arange(u.size)-barw/2-0.01,dObs['data']['Soil']['MIN_C_THA']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations, mineral horizon (Shaw et al. 2018)')
		ax.bar(np.arange(u.size)-barw/2-0.01,dObs['data']['Soil']['ORG_C_THA']['mu'],barw,facecolor=cl['bd'],bottom=dObs['data']['Soil']['MIN_C_THA']['mu'],label='Ground plot observations, organic horizon (Shaw et al. 2018)')
		ax.bar(np.arange(u.size)+barw/2+0.01,dM['Min_mu'],barw,facecolor=cl['gl'],label='Prediction, mineral horizon (FCS)')
		ax.bar(np.arange(u.size)+barw/2+0.01,dM['Org_mu'],barw,facecolor=cl['gd'],bottom=dM['Min_mu'],label='Prediction, organic horizon (FCS)')
		ax.errorbar(np.arange(u.size)-barw/2-0.01,dObs['data']['Soil']['TOT_C_THA']['mu'],yerr=dObs['data']['Soil']['TOT_C_THA']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax.errorbar(np.arange(u.size)+barw/2+0.01,dM['Tot_mu'],yerr=dM['Tot_se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		for i in range(u.size):
			if (np.abs(Dp[i])>10000) | (np.isnan(yo[i])==True):
				continue
			if Dp[i]>=0:
				a='+'
			else:
				a=''
			ax.text(i+barw/2+0.01,yp[i]+dM['Tot_se'][i]+15,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=5)
		ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,450],xticks=np.arange(u.size),
		 xticklabels=lab,ylabel='Soil organic carbon (tC ha$^{-1}$ yr$^{-1}$)')
		plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25,fontsize=6)
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareSOCByBGC_' + str(iScn+1),'png',900)
	
	return

#%%
def QA_FullComparison_AgeResponsesBiomassAndNetGrowth_ByReg_CNY(meta,pNam,dObs,dMod):
	lw=0.5; ms=3; cl=meta['Graphics']['GP Comp']
	ptf='PTF CN'
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(22,7))
		reg='Coast'
		ax1[0].plot(dObs['bin'],dObs['data'][ptf][reg]['Ctot L t0']['mu'],'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],label='Observed biomass',zorder=1)
		ax1[0].plot(dObs['bin'],dMod[iScn][reg]['C_Biomass_Tot']['mu'],'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],label='Predicted biomass',zorder=1)
		ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dObs['bw']),yticks=np.arange(0,500,50),xlim=[0,250+dObs['bw']],ylim=[0,250])
		ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax2=ax1[0].twinx()
		ax2.bar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot Net']['mu'],0.45*dObs['bw'],ec='none',fc=cl['bl'],zorder=-1,label='Observed net growth')
		ax2.bar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Net_Tot']['mu'] ,0.45*dObs['bw'],ec='none',fc=cl['gl'],zorder=-1,label='Predicted net growth')
		ax2.plot([0,500],[0,0],'-k',lw=lw)
		ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
		ax1[0].set_zorder(ax2.get_zorder()+1)
		ax1[0].patch.set_visible(False)
		ax2.tick_params(length=meta['Graphics']['gp']['tickl'])
		reg='Interior'
		ax1[1].plot(dObs['bin'],dObs['data'][ptf][reg]['Ctot L t0']['mu'],'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],zorder=1)
		ax1[1].plot(dObs['bin'],dMod[iScn][reg]['C_Biomass_Tot']['mu'],'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],zorder=1)
		ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dObs['bw']),yticks=np.arange(0,500,50),xlim=[0,250+dObs['bw']],ylim=[0,250])
		ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax3=ax1[1].twinx()
		ax3.bar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot Net']['mu'],0.45*dObs['bw'],ec='none',fc=cl['bl'],label='Observed net growth',zorder=-1)
		ax3.bar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Net_Tot']['mu'],0.45*dObs['bw'],ec='none',fc=cl['gl'],label='Predicted net growth',zorder=-1)
		ax3.plot([0,500],[0,0],'-k',lw=lw)
		ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
		ax1[1].set_zorder(ax3.get_zorder()+1)
		ax1[1].patch.set_visible(False)
		ax3.tick_params(length=meta['Graphics']['gp']['tickl'])
		ax1[0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
		ax3.legend(loc='upper center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
		gu.axletters(ax1,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullComparison_AgeResponsesBiomassAndNetGrowth_ByReg_CNY_' + str(iScn+1),'png',900)
	return

#%%
def QA_FullComparison_AgeResponsesGrossGrowthAndMortality_ByReg_CNY(meta,pNam,dObs,dMod):
	cl=meta['Graphics']['GP Comp']
	ptf='PTF CN'
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(22,7.25))
		reg='Coast'
		ax1[0].bar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot G Surv']['mu']+dObs['data'][ptf][reg]['Ctot G Recr']['mu'],0.4*dObs['bw'],ec='none',fc=cl['bl'],label='Observed gross growth')
		ax1[0].errorbar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot G Surv']['mu']+dObs['data'][ptf][reg]['Ctot G Recr']['mu'],yerr=dObs['data'][ptf][reg]['Ctot G Surv']['se']+dObs['data'][ptf][reg]['Ctot G Recr']['se'],color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[0].bar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Gross_Tot']['mu'],0.4*dObs['bw'],ec='none',fc=cl['gl'],label='Predicted gross growth')
		ax1[0].errorbar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Gross_Tot']['mu'],yerr=dMod[iScn][reg]['C_G_Gross_Tot']['se'],color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		
		ax1[0].bar(dObs['bin']-(0.45*dObs['bw']/2),-dObs['data'][ptf][reg]['Ctot Mort']['mu'],0.4*dObs['bw'],ec='none',fc=cl['bd'],label='Observed mortality')
		ax1[0].errorbar(dObs['bin']-(0.45*dObs['bw']/2),-dObs['data'][ptf][reg]['Ctot Mort']['mu'],yerr=dObs['data'][ptf][reg]['Ctot Mort']['se'],color=0.5*cl['bd'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[0].bar(dObs['bin']+(0.45*dObs['bw']/2),-dMod[iScn][reg]['C_M_Tot']['mu'],0.4*dObs['bw'],ec='none',fc=cl['gd'],label='Predicted mortality')
		ax1[0].errorbar(dObs['bin']+(0.45*dObs['bw']/2),-dMod[iScn][reg]['C_M_Tot']['mu'],yerr=dMod[iScn][reg]['C_M_Tot']['se'],color=0.5*cl['gd'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[0].plot([0,300],[0,0],'k-',lw=0.5)
		ax1[0].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dObs['bw']),yticks=np.arange(-10,500,1),xlim=[0,250+dObs['bw']],ylim=[-3,7])
		ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax1[0].legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25,fontsize=6)
		
		reg='Interior'
		ax1[1].bar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot G Surv']['mu']+dObs['data'][ptf][reg]['Ctot G Recr']['mu'],0.4*dObs['bw'],ec='none',fc=cl['bl'])
		ax1[1].errorbar(dObs['bin']-(0.45*dObs['bw']/2),dObs['data'][ptf][reg]['Ctot G Surv']['mu']+dObs['data'][ptf][reg]['Ctot G Recr']['mu'],yerr=dObs['data'][ptf][reg]['Ctot G Surv']['se']+dObs['data'][ptf][reg]['Ctot G Recr']['se'],color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[1].bar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Gross_Tot']['mu'],0.4*dObs['bw'],ec='none',fc=cl['gl'])
		ax1[1].errorbar(dObs['bin']+(0.45*dObs['bw']/2),dMod[iScn][reg]['C_G_Gross_Tot']['mu'],yerr=dMod[iScn][reg]['C_G_Gross_Tot']['se'],color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

		ax1[1].bar(dObs['bin']-(0.45*dObs['bw']/2),-dObs['data'][ptf][reg]['Ctot Mort']['mu'],0.4*dObs['bw'],ec='none',fc=cl['bd'])
		ax1[1].errorbar(dObs['bin']-(0.45*dObs['bw']/2),-dObs['data'][ptf][reg]['Ctot Mort']['mu'],yerr=dObs['data'][ptf][reg]['Ctot Mort']['se'],color=0.5*cl['bd'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[1].bar(dObs['bin']+(0.45*dObs['bw']/2),-dMod[iScn][reg]['C_M_Tot']['mu'],0.4*dObs['bw'],ec='none',fc=cl['gd'])
		ax1[1].errorbar(dObs['bin']+(0.45*dObs['bw']/2),-dMod[iScn][reg]['C_M_Tot']['mu'],yerr=dMod[iScn][reg]['C_M_Tot']['se'],color=0.5*cl['gd'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
		ax1[1].plot([0,300],[0,0],'k-',lw=0.5)
		ax1[1].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dObs['bw']),yticks=np.arange(-10,500,1),xlim=[0,250+dObs['bw']],ylim=[-3,7])
		ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
		gu.axletters(ax1,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullComparison_AgeResponsesGrossGrowthAndMortality_ByReg_CNY_' + str(iScn+1),'png',900)
	return

#%%

def QA_AgeResponseBiomassNetGrowth_CN(meta,pNam,iScn,gpt):
    #E=[None]*meta[pNam]['Project']['N Scenario']
    #for iScn in range(meta[pNam]['Project']['N Scenario']):
    #E[iScn]={'Coast':{},'Interior':{}}
    bw=25; bin=np.arange(bw,250+bw,bw)
    #ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Mod C_M_Harv']==0) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) & (gpt[iScn]['Mod C_M_Harv']==0) )[0]
    
    xO=gpt[iScn]['Age Med t0'][ind]    
    xM=gpt[iScn]['Mod A t0'][ind]
    yO=gpt[iScn]['Ctot L t0'][ind]
    yM=gpt[iScn]['Mod C_Biomass_Tot t0'][ind]
    N,mu1,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu2,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu1+mu2)==False)[0]
    E['Coast']['Biomass'],txt=gu.GetRegStats(mu1[ikp],mu2[ikp])
    #E[0]=np.sqrt(np.sum((mu1[ikp]-mu2[ikp])**2))
    #E[0]=np.nanmedian((mu2[ikp]-mu1[ikp])/mu1[ikp]*100)        
    
    yO=gpt[iScn]['Ctot Net'][ind]
    yM=gpt[iScn]['Mod C_G_Net'][ind]
    N,mu3,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu4,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu3+mu4)==False)[0]
    E['Coast']['Net Growth'],txt=gu.GetRegStats(mu3[ikp],mu4[ikp])
    #E[1]=np.sqrt(np.sum((mu3[ikp]-mu4[ikp])**2))
    #E[1]=np.nanmedian((mu4[ikp]-mu3[ikp])/mu3[ikp]*100)
    
    ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Mod C_M_Harv']==0) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt[iScn]['Age Med t0'][ind]
    xM=gpt[iScn]['Mod A t0'][ind]
    yO=gpt[iScn]['Ctot L t0'][ind]
    yM=gpt[iScn]['Mod C_Biomass_Tot t0'][ind]
    N,mu5,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu6,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu5+mu6)==False)[0]
    E['Interior']['Biomass'],txt=gu.GetRegStats(mu5[ikp],mu6[ikp])
    #E[2]=np.sqrt(np.sum((mu5[ikp]-mu6[ikp])**2))
    #E[2]=np.nanmedian((mu5[ikp]-mu4[ikp])/mu4[ikp]*100)
    
    yO=gpt[iScn]['Ctot Net'][ind]
    yM=gpt[iScn]['Mod C_G_Net'][ind]
    N,mu7,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu8,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu7+mu8)==False)[0]
    E['Interior']['Net Growth'],txt=gu.GetRegStats(mu7[ikp],mu8[ikp])
    #E[3]=np.sqrt(np.sum((mu7[ikp]-mu8[ikp])**2))
    #E[3]=np.nanmedian((mu8[ikp]-mu7[ikp])/mu7[ikp]*100)
    
    if 'Just get errors' not in meta[pNam].keys():            
        lw=0.5; ms=3; cl=meta['Graphics']['GP Comp']
        plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,6))
        # Coast    
        ax1[0].plot(bin,mu1,'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],label='Observed biomass',zorder=1)    
        ax1[0].plot(bin,mu2,'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],label='Predicted biomass',zorder=1)
        ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,50),xlim=[0,250+bw],ylim=[0,400])
        ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])        
        ax2=ax1[0].twinx()
        ax2.bar(bin-(0.45*bw/2),mu3,0.4*bw,ec='none',fc=cl['bl'],zorder=-1)
        ax2.bar(bin+(0.45*bw/2),mu4,0.45*bw,ec='none',fc=cl['gl'],zorder=-1)
        ax2.plot([0,500],[0,0],'-k',lw=lw)
        ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
        ax1[0].set_zorder(ax2.get_zorder()+1)
        ax1[0].patch.set_visible(False)
        ax2.tick_params(length=meta['Graphics']['gp']['tickl'])                
        # Interior    
        ax1[1].plot(bin,mu5,'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],zorder=1)    
        ax1[1].plot(bin,mu6,'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],zorder=1)
        ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,20),xlim=[0,250+bw],ylim=[0,200])
        ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])        
        ax3=ax1[1].twinx()    
        ax3.bar(bin-(0.45*bw/2),mu7,0.4*bw,ec='none',fc=cl['bl'],label='Observed net growth',zorder=-1)    
        ax3.bar(bin+(0.45*bw/2),mu8,0.45*bw,ec='none',fc=cl['gl'],label='Predicted net growth',zorder=-1)    
        ax3.plot([0,500],[0,0],'-k',lw=lw)
        ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
        ax1[1].set_zorder(ax3.get_zorder()+1)
        ax1[1].patch.set_visible(False)
        ax3.tick_params(length=meta['Graphics']['gp']['tickl'])
        ax1[0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
        ax3.legend(loc='upper center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
        gu.axletters(ax1,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
        plt.tight_layout()
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseBiomassNetGrowth_CN_' + str(iScn+1),'png',900)
    
    return E

#%% Age
def EvalAtPlots_AgeByBGCZ_CNV(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CNV']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Age VRI t0']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Age VRI t0']['mu']
        y=d['Mod A t0']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Age VRI t0']['N']>=30) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(215,30,txt,fontsize=10,color='k',ha='right')
        ax.text(230,230,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed age (years)',ylabel='Predicted age (years)',xlim=[0,250],ylim=[0,250])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Scatterplot_CNV_' + str(iScn+1),'png',900)
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Age VRI t0']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod A t0']['mu']
        yo=d['Age VRI t0']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Age VRI t0']['mu'],barw,facecolor=cl['bl'],label='Ground plots')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod A t0']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Age VRI t0']['mu'],yerr=d['Age VRI t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod A t0']['mu'],yerr=d['Mod A t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod A t0']['se'][i]+6,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Age (years)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Barchart_CNV_' + str(iScn+1),'png',900)

    return

#%%
def EvalAtPlots_BiomassByBGC_CNV(meta,pNam,gpt):

    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CNV']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot L t0']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # # Scatterplot
        # flg=0
        # if flg==1:
        #     x=d['Ctot L t0']['mu']
        #     y=d['Mod C_Biomass_Tot t0']['mu']
        #     ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot L t0']['N']>=30) )[0]
        #     rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        #     fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        #     ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #     #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        #     for i in range(ikp.size):
        #         ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        #     ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        #     ax.text(200,20,txt,fontsize=10,color='k',ha='right')
        #     ax.text(190,190,'1:1',fontsize=8,ha='center')
        #     ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (tC ha$^{-1}$)',ylabel='Predicted biomass (tC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
        #     #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        #     ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        #     if meta['Graphics']['Print Figures']=='On':
        #         gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_CNV_' + str(iScn+1),'png',900)
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot L t0']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_Biomass_Tot t0']['mu']
        yo=d['Ctot L t0']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl['bl'],label='Ground plots')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],yerr=d['Mod C_Biomass_Tot t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_Biomass_Tot t0']['se'][i]+6,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_CNV_' + str(iScn+1),'png',900)

        # CN only
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
        
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot L t0']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot L t0']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_Biomass_Tot t0']['mu']
        yo=d['Ctot L t0']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl['bl'],label='Ground plots')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],yerr=d['Mod C_Biomass_Tot t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_Biomass_Tot t0']['se'][i]+6,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_CN_' + str(iScn+1),'png',900)
    return

#%% Biomass (YSM)
def EvalAtPlots_BiomassByBGC_YSM(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF YSM']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot L t0']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot L t0']['mu']
        y=d['Mod C_Biomass_Tot t0']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot L t0']['N']>=30) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(200,20,txt,fontsize=10,color='k',ha='right')
        ax.text(190,190,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (tC ha$^{-1}$)',ylabel='Predicted biomass (tC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_YSM_' + str(iScn+1),'png',900)
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot L t0']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_Biomass_Tot t0']['mu']
        yo=d['Ctot L t0']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl['bl'],label='Ground plots')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],yerr=d['Mod C_Biomass_Tot t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_Biomass_Tot t0']['se'][i]+6,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,150])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_YSM_' + str(iScn+1),'png',900)

    return

#%% Biomass (TIPSY stemwood)
def EvalAtPlots_StemwoodFromTIPSYByBGC_CNV(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CNV']==1) &
                             (gpt[iScn]['Csw L t0']>=0) & (gpt[iScn]['Csw L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Csw L t0']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Csw L t0']['mu']
        y=d['GY Csw t0']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot L t0']['N']>=30) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
        print(rs['Mean Dif'])
        print(rs['Dif (%)'])
    
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(200,20,txt,fontsize=10,color='k',ha='right')
        ax.text(190,190,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (tC ha$^{-1}$)',ylabel='Predicted biomass (tC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_BiomassSW_ByBGCZone_Scatterplot_' + str(iScn+1),'png',900)
    
        # # Plot bar chart
    
        # # Put in order
        # ord=np.argsort(d['Ctot L t0']['mu'])
        # lab=np.flip(lab[ord])
        # for v in d:
        #     for k in d[v].keys():
        #         d[v][k]=np.flip(d[v][k][ord])
    
        # Area=np.zeros(lab.size)
        # for i in range(lab.size):
        #     ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
        #     Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # # Area weighting
        # for v in d:
        #     d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
        #     d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        # lab=np.append(lab,'Weighted\naverage')
        # u=np.append(u,0.0)
    
        # cl=meta['Graphics']['GP Comp']
        # barw=0.32
        # plt.close('all')
        # fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        # ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl['bl'],label='Ground plots')
        # ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        # ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        # ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],yerr=d['Mod C_Biomass_Tot t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        # #for i in range(u.size):
        # #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        # ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
        # plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        # ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        # gu.PrintFig(meta['Paths']['Figures'] + '\\QA_Biomass_ByBGCZone_BarChart_' + str(iScn+1),'png',900)

    return

#%% Biomass components
def EvalAtPlots_BiomassComponents_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.array([1])
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
        for i in range(u.size):
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        D=np.zeros(5)
        D[0]=(d['Mod C_Stemwood_Tot t0']['mu'][0]-d['Csw L t0']['mu'][0])/d['Csw L t0']['mu'][0]*100
        D[1]=(d['Mod C_Bark_Tot t0']['mu'][0]-d['Cbk L t0']['mu'][0])/d['Cbk L t0']['mu'][0]*100
        D[2]=(d['Mod C_Branch_Tot t0']['mu'][0]-d['Cbr L t0']['mu'][0])/d['Cbr L t0']['mu'][0]*100
        D[3]=(d['Mod C_Foliage_Tot t0']['mu'][0]-d['Cf L t0']['mu'][0])/d['Cf L t0']['mu'][0]*100
        D[4]=(d['Mod C_Root_Tot t0']['mu'][0]-d['Cr L t0']['mu'][0])/d['Cr L t0']['mu'][0]*100
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6))
        ax.bar(np.arange(5),D,facecolor=[0.8,0.8,0.8])
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(D.size),xticklabels=['Stemwood','Bark','Branch','Foliage','Roots'],ylabel='Mean difference (%)',xlim=[-0.5,D.size-0.5],ylim=[0,100])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassComponents_CN_' + str(iScn+1),'png',900)

    return

#%% Biomass dynammics summary (total CO2e)
def EvalAtPlots_BiomassDynamicsAve_TotCO2e_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000) &
                             (gpt[iScn]['Ctot G Tot']>0) & (gpt[iScn]['Ctot G Tot']<30))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot G Tot']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Area weighting
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Area of forest
        #ind=np.where(zLCC1['Data']==lut_1ha['lcc1']['Forest Land'])
        #A_Tot=ind[0].size/1e6*3.67
        A_Tot=62000000/1e6*3.67
    
        cl=meta['Graphics']['GP Comp']
        lab=['Gross\ngrowth','Mortality','Net\ngrowth'] #,'Harvest\nmortality'
    
        barw=0.38
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6))
        ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
        ax.bar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1]*A_Tot,barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1]*A_Tot,barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1]*A_Tot,yerr=d['Ctot G Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1]*A_Tot,yerr=d['Mod C_G_Gross_Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.bar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1]*A_Tot,barw,facecolor=cl['bl'])
        ax.bar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1]*A_Tot,barw,facecolor=cl['gl'])
        ax.errorbar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1]*A_Tot,yerr=d['Ctot Mort+Lost']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1]*A_Tot,yerr=d['Mod C_M_Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.bar(3-barw/2-0.01,d['Ctot Net']['mu'][-1]*A_Tot,barw,facecolor=cl['bl'])
        ax.bar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1]*A_Tot,barw,facecolor=cl['gl'])
        ax.errorbar(3-barw/2-0.01,d['Ctot Net']['mu'][-1]*A_Tot,yerr=d['Ctot Net']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1]*A_Tot,yerr=d['Mod C_G_Net']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-500,600,100),ylabel='Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)',xlim=[0.5,3.5],ylim=[-500,500])
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        ax.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsTotal_CN_' + str(iScn+1),'png',900)

    return

#%% Biomass dynammics summary (CN average)
def EvalAtPlots_BiomassDynamicsAve_CN(meta,pNam,gpt):

    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) & (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ctot Net']>=-1000) & (gpt[iScn]['Ctot Net']<1000) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot Net']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Area weighting
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
            
        # z=u1ha.Import_Raster(meta,[],['lcc1_c','bgcz'])
        # dA={}
        # for k in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
        #     ind=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][k]) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
        #     dA[k]=ind[0].size 
        
        # Area=np.zeros(lab.size)
        # for k in meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys():
        #     ind1=np.where(lab==k)[0]
        #     Area[ind1]=dA[k]
            
        wa={}
        for v in d:
            wa[v]={}
            wa[v]['mu']=np.nansum(d[v]['mu']*Area)/np.nansum(Area)
            wa[v]['se']=np.nansum(d[v]['se']*Area)/np.nansum(Area)
    
        #wa['Ctot G Tot']['mu']=wa['Ctot G Surv']['mu']+wa['Ctot G Recr']['mu']
    
        lab=['Gross\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #,'Harvest\nmortality'
        cl=meta['Graphics']['GP Comp']
        barw=0.38
        plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(9,7))
        ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
        ax.bar(1-barw/2-0.01,wa['Ctot G Tot']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(1+barw/2+0.01,wa['Mod C_G_Gross_Tot']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(1-barw/2-0.01,wa['Ctot G Tot']['mu'],yerr=wa['Ctot G Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(1+barw/2+0.01,wa['Mod C_G_Gross_Tot']['mu'],yerr=wa['Mod C_G_Gross_Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.bar(2-barw/2-0.01,-wa['Ctot Mort Nat']['mu'],barw,facecolor=cl['bl'])
        ax.bar(2+barw/2+0.01,-wa['Mod C_M_Nat']['mu'],barw,facecolor=cl['gl'])
        ax.errorbar(2-barw/2-0.01,-wa['Ctot Mort Nat']['mu'],yerr=wa['Ctot Mort Nat']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(2+barw/2+0.01,-wa['Mod C_M_Nat']['mu'],yerr=wa['Mod C_M_Nat']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        
        ax.bar(3-barw/2-0.01,-wa['Ctot Mort Harv']['mu'],barw,facecolor=cl['bl'])
        ax.bar(3+barw/2+0.01,-wa['Mod C_M_Harv']['mu'],barw,facecolor=cl['gl'])
        ax.errorbar(3-barw/2-0.01,-wa['Ctot Mort Harv']['mu'],yerr=wa['Ctot Mort Harv']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(3+barw/2+0.01,-wa['Mod C_M_Harv']['mu'],yerr=wa['Mod C_M_Harv']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.bar(4-barw/2-0.01,wa['Ctot Net']['mu'],barw,facecolor=cl['bl'])
        ax.bar(4+barw/2+0.01,wa['Mod C_G_Net']['mu'],barw,facecolor=cl['gl'])
        ax.errorbar(4-barw/2-0.01,wa['Ctot Net']['mu'],yerr=wa['Ctot Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(4+barw/2+0.01,wa['Mod C_G_Net']['mu'],yerr=wa['Mod C_G_Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-2,3,0.5),ylabel='Carbon balance of trees (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-1.5,2])
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        ax.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage__CN_' + str(iScn+1),'png',900)

    return

#%% Biomass dynammics summary (average YSM)
def EvalAtPlots_BiomassDynamicsAve_YSM(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF YSM']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000) &
                             (gpt[iScn]['Ctot G Tot']>0) & (gpt[iScn]['Ctot G Tot']<30))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot G Tot']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Area weighting
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        cl=meta['Graphics']['GP Comp']['cl']
        cle=[0.05,0.2,0.45]
        lab=['Gross\ngrowth','Natural\nmortality','Net\ngrowth'] #,'Harvest\nmortality'
    
        barw=0.38
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6))
        ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
        ax.bar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],yerr=d['Ctot G Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],yerr=d['Mod C_G_Gross_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.bar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],barw,facecolor=cl['bl'])
        ax.bar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],barw,facecolor=cl['gl'])
        ax.errorbar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],yerr=d['Ctot Mort+Lost']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],yerr=d['Mod C_M_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.bar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],barw,facecolor=cl['bl'])
        ax.bar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],barw,facecolor=cl['gl'])
        ax.errorbar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],yerr=d['Ctot Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],yerr=d['Mod C_G_Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    
        ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-4,6,1),ylabel='Carbon balance of trees (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,3.5],ylim=[-3,4])
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        ax.legend(frameon=False,loc='lower left',facecolor=[1,1,1],labelspacing=0.25)
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage_YSM_' + str(iScn+1),'png',900)

    return

#%% Gross growth (CN)
def EvalAtPlots_GrossGrowthByBGC_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000) &
                             (gpt[iScn]['Ctot G Tot']>0) & (gpt[iScn]['Ctot G Tot']<30))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot G Tot']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot G Tot']['mu']
        y=d['Mod C_G_Gross_Tot']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot G Tot']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=3,color=[0.8,0.8,0.8])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=6)
        ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(4.1,0.45,txt,fontsize=10,color='k',ha='right')
        ax.text(4.05,4.05,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed gross growth (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted gross growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0,4.5],ylim=[0,4.5])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_CN_' + str(iScn+1),'png',900)
        plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot G Tot']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_G_Gross_Tot']['mu']
        yo=d['Ctot G Tot']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],yerr=d['Ctot G Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],yerr=d['Mod C_G_Gross_Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_G_Gross_Tot']['se'][i]+0.07,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=5)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Gross growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,4.5])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_CN_' + str(iScn+1),'png',900)

    return

#%% Gross growth (YSM)
def EvalAtPlots_GrossGrowthByBGC_YSM(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF YSM']==1) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000) &
                             (gpt[iScn]['Ctot G Tot']>0) & (gpt[iScn]['Ctot G Tot']<30))[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot G Tot']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot G Tot']['mu']
        y=d['Mod C_G_Gross_Tot']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot G Tot']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=3,color=[0.8,0.8,0.8])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
        ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed gross growth (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted gross growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0,5.5],ylim=[0,5.5])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_YSM_' + str(iScn+1),'png',900)
        plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot G Tot']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_G_Gross_Tot']['mu']
        yo=d['Ctot G Tot']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],yerr=d['Ctot G Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],yerr=d['Mod C_G_Gross_Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_G_Gross_Tot']['se'][i]+0.07,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Gross growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,5])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_YSM_' + str(iScn+1),'png',900)

    return

#%% Mortality
def EvalAtPlots_MortalityByBGC_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot Mort+Lost']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot Mort+Lost']['mu']
        y=d['Mod C_M_Tot']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot Mort+Lost']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=3,color=[0.8,0.8,0.8])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=6)
        ax.plot(rs['xhat Line'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(2.8,0.45,txt,fontsize=10,color='k',ha='right')
        ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed mortality (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted mortality (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0,3],ylim=[0,3])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_CN_' + str(iScn+1),'png',900)
        plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot Mort+Lost']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_M_Tot']['mu']
        yo=d['Ctot Mort+Lost']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],yerr=d['Ctot Mort+Lost']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],yerr=d['Mod C_M_Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_M_Tot']['se'][i]+0.07,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Mortality (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,4])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Barchart_CN_' + str(iScn+1),'png',900)

    return

#%% Mortality (YSM)
def EvalAtPlots_MortalityByBGC_YSM(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF YSM']==1) &
                             (gpt[iScn]['Mod C_M_Tot']>=0) &
                             (gpt[iScn]['Ctot L t0']>=0) & (gpt[iScn]['Ctot L t0']<10000) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot Mort+Lost']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot Mort+Lost']['mu']
        y=d['Mod C_M_Tot']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot Mort+Lost']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
        ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed mortality (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted mortality (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0,3],ylim=[0,3])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM_' + str(iScn+1),'png',900)
        plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot Mort+Lost']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_M_Tot']['mu']
        yo=d['Ctot Mort+Lost']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],yerr=d['Ctot Mort+Lost']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],yerr=d['Mod C_M_Tot']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['Mod C_M_Tot']['se'][i]+0.07,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Mortality (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,6])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM_' + str(iScn+1),'png',900)

    return

#%% Net growth (CN)

def EvalAtPlots_GrowthNetByBGC_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot Net']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot Net']['mu']
        y=d['Mod C_G_Net']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot Mort+Lost']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([-200,500],[-200,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(1.8,-2.75,txt,fontsize=10,color='k',ha='right')
        ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-3,2],ylim=[-3,2])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_CN_' + str(iScn+1),'png',900)
        plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot Net']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_G_Net']['mu']
        yo=d['Ctot Net']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.plot([-1,20],[0,0],'k-',lw=0.5)
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],yerr=d['Ctot Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],yerr=d['Mod C_G_Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+';
            else:
                a='';
            if yp[i]>=0:
                a2=1
            else:
                a2=-1
            ax.text(i+barw/2+0.01,yp[i]+a2*(d['Mod C_G_Net']['se'][i]+0.25),a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[-3,2])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_CN_' + str(iScn+1),'png',900)

    return

#%% Net growth (YSM)

def EvalAtPlots_GrowthNetByBGC_YSM(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF YSM']==1) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Ctot Net']['N']>=3)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Ctot Net']['mu']
        y=d['Mod C_G_Net']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot Mort+Lost']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([-200,500],[-200,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(1.8,-2.75,txt,fontsize=10,color='k',ha='right')
        ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-3,5],ylim=[-3,5])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_YSM_' + str(iScn+1),'png',900)
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['Ctot Net']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['Mod C_G_Net']['mu']
        yo=d['Ctot Net']['mu']
        Dp=(yp-yo)/yo*100
    
        cl=meta['Graphics']['GP Comp']
        barw=0.32
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.plot([-1,20],[0,0],'k-',lw=0.5)
        ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],yerr=d['Ctot Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],yerr=d['Mod C_G_Net']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+';
            else:
                a='';
            if yp[i]>=0:
                a2=1
            else:
                a2=-1
            ax.text(i+barw/2+0.01,yp[i]+a2*(d['Mod C_G_Net']['se'][i]+0.25),a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    
        #for i in range(u.size):
        #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[-3,2])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_YSM_' + str(iScn+1),'png',900)
    return

#%% Net growth (TIPSY stemwood)

def EvalAtPlots_GrowthNetStemwoodTIPSY_CN(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Unique BGC zones
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) &
                             (gpt[iScn]['Csw Mort']<100) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Csw Net']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Csw Net']['mu']
        y=d['GY Csw Net']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Csw Net']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
        print(rs['Mean Dif'])
        #print(rs['Dif (%)'])
    
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([-200,500],[-200,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        ax.text(1.8,-0.75,txt,fontsize=10,color='k',ha='right')
        ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (tC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-1,2],ylim=[-1,2])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN_' + str(iScn+1),'png',900)
    
        # With no disturbance mortality
        d={}
        for v in gpt[iScn]['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
    
        for i in range(u.size):
            lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
            for v in gpt[iScn]['vaL']:
                ind=np.where( (gpt[iScn]['Ecozone BC L1']==u[i]) &
                             (gpt[iScn]['PTF CN']==1) &
                             (gpt[iScn]['Csw Mort']<0.2) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gpt[iScn][v][ind])
                d[v]['sd'][i]=np.nanstd(gpt[iScn][v][ind])
                d[v]['se'][i]=np.nanstd(gpt[iScn][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['Csw Net']['N']>=10)[0]
        for v in gpt[iScn]['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # Scatterplot
        x=d['Csw Net']['mu']
        y=d['GY Csw Net']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Csw Net']['N']>=10) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])
        print(rs['Mean Dif'])
        #print(rs['Dif (%)'])
    
        ax.plot(rs['xhat Lnie'],rs['yhat Line'],'r-',lw=1,label='Best fit')
    
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN_' + str(iScn+1),'png',900)

    return

#%%
def EvalAtPlots_SOCByBGC_ShawComp(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        u=np.unique(gpt[iScn]['Ecozone BC L1'][gpt[iScn]['Ecozone BC L1']>0])
        lab=np.array(['' for _ in range(u.size)],dtype=object)
    
        d={}
        for v in gpt[iScn]['soils']['vaL']:
            d[v]={}
            d[v]['N']=np.zeros(u.size)
            d[v]['mu']=np.zeros(u.size)
            d[v]['sd']=np.zeros(u.size)
            d[v]['se']=np.zeros(u.size)
            for i in range(u.size):
                lab[i]=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[i])[0]
                ind=np.where(gpt[iScn]['soils']['bgcz']==u[i])[0]
                if ind.size>0:
                    ind=np.where( (gpt[iScn]['soils']['bgcz']==u[i]) & (gpt[iScn]['soils']['TOT_C_THA']>0) )[0]
                    d[v]['N'][i]=ind.size
                    d[v]['mu'][i]=np.nanmean(gpt[iScn]['soils'][v][ind])
                    d[v]['sd'][i]=np.nanstd(gpt[iScn]['soils'][v][ind])
                    d[v]['se'][i]=np.nanstd(gpt[iScn]['soils'][v][ind])/np.sqrt(ind.size)
    
        # Remove classes with inadequate data
        ind=np.where(d['TOT_C_THA']['N']>=10)[0]
        for v in gpt[iScn]['soils']['vaL']:
            for k in d[v].keys():
                d[v][k]=d[v][k][ind]
        u=u[ind]
        lab=lab[ind]
    
        # # Scatterplot
        # x=d['TOT_C_THA']['mu']
        # y=d['C_Soil_Tot']['mu']
        # ikp=np.where( (np.isnan(x+y)==False) & (d['TOT_C_THA']['N']>=10) )[0]
        # rs,txt=gu.GetRegStats(x[ikp],y[ikp])
    
        # plt.close('all')
        # fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        # ax.plot([0,1000],[0,1000],'-k',lw=2,color=[0.8,0.8,0.8])
        # #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        # ax.plot(rs['xhat Lnie'],rs['yhat Line'],'k-',lw=1,label='Best fit')
        # for i in range(ikp.size):
        #     ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        # ax.text(330,30,txt,fontsize=10,color='k',ha='right')
        # ax.text(300,300,'1:1',fontsize=8,ha='center')
        # ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed SOC (tC ha$^{-1}$)',ylabel='Predicted SOC (tC ha$^{-1}$)',xlim=[0,350],ylim=[0,350])
        # #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        # ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        # if meta['Graphics']['Print Figures']=='On':
        #     gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_Scatterplot_' + str(iScn+1),'png',900)
        # plt.close('all')
    
        # Plot bar chart
    
        # Put in order
        ord=np.argsort(d['TOT_C_THA']['mu'])
        lab=np.flip(lab[ord])
        for v in d:
            for k in d[v].keys():
                d[v][k]=np.flip(d[v][k][ord])
    
        Area=np.zeros(lab.size)
        for i in range(lab.size):
            ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
            Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
    
        # Area weighting
        for v in d:
            d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
            d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
        lab=np.append(lab,'Weighted\naverage')
        u=np.append(u,0.0)
    
        # Percent difference
        yp=d['C_Soil_Tot']['mu']+d['C_Soil_OHorizon']['mu']
        yo=d['TOT_C_THA']['mu']
        Dp=(yp-yo)/yo*100
    
        barw=0.32
        cl=meta['Graphics']['GP Comp']
    
        plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
        ax.bar(np.arange(u.size)-barw/2-0.01,d['MIN_C_THA']['mu'],barw,facecolor=cl['bl'],label='Ground plot observations, mineral horizon (Shaw et al. 2018)')
        ax.bar(np.arange(u.size)-barw/2-0.01,d['ORG_C_THA']['mu'],barw,facecolor=cl['bd'],bottom=d['MIN_C_THA']['mu'],label='Ground plot observations, organic horizon (Shaw et al. 2018)')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['C_Soil_Tot']['mu'],barw,facecolor=cl['gl'],label='Prediction, mineral horizon (FCS)')
        ax.bar(np.arange(u.size)+barw/2+0.01,d['C_Soil_OHorizon']['mu'],barw,facecolor=cl['gd'],bottom=d['C_Soil_Tot']['mu'],label='Prediction, organic horizon (FCS)')
        ax.errorbar(np.arange(u.size)-barw/2-0.01,d['MIN_C_THA']['mu']+d['ORG_C_THA']['mu'],yerr=d['TOT_C_THA']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(u.size)+barw/2+0.01,d['C_Soil_Tot']['mu']+d['C_Soil_OHorizon']['mu'],yerr=d['C_Soil_Tot']['se']+d['C_Soil_OHorizon']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        for i in range(u.size):
            if Dp[i]>=0:
                a='+'
            else:
                a=''
            ax.text(i+barw/2+0.01,yp[i]+d['C_Soil_Tot']['se'][i]+d['C_Soil_OHorizon']['se'][i]+13,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=5)
        #    ax.text(i,8,str(d['TOT_C_THA']['N'][i].astype(int)),color=meta['Graphics']['gp']['cla'],ha='center',fontsize=7)
        #    #ax.text(i,30,str(d['Model SS'][i].astype(int)),color='c',ha='center',fontsize=8)
        ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,450],xticks=np.arange(u.size),
               xticklabels=lab,ylabel='Soil organic carbon (tC ha$^{-1}$ yr$^{-1}$)')
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25,fontsize=6)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_' + str(iScn+1),'png',900)

    return

#%%

def EvalAtPlots_AgeResponsesBiomassAndNetGrowth_ByReg_CN(meta,pNam,iScn,gpt,E):
    #E=[None]*meta[pNam]['Project']['N Scenario']
    #for iScn in range(meta[pNam]['Project']['N Scenario']):
    #E[iScn]={'Coast':{},'Interior':{}}
    bw=25; bin=np.arange(bw,250+bw,bw)        
    #ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Mod C_M_Harv']==0) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) & (gpt[iScn]['Mod C_M_Harv']==0) )[0]
    
    xO=gpt[iScn]['Age Med t0'][ind]    
    xM=gpt[iScn]['Mod A t0'][ind]
    yO=gpt[iScn]['Ctot L t0'][ind]
    yM=gpt[iScn]['Mod C_Biomass_Tot t0'][ind]
    N,mu1,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu2,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu1+mu2)==False)[0]
    E['Coast']['Biomass'],txt=gu.GetRegStats(mu1[ikp],mu2[ikp])
    #E[0]=np.sqrt(np.sum((mu1[ikp]-mu2[ikp])**2))
    #E[0]=np.nanmedian((mu2[ikp]-mu1[ikp])/mu1[ikp]*100)        
    
    yO=gpt[iScn]['Ctot Net'][ind]
    yM=gpt[iScn]['Mod C_G_Net'][ind]
    N,mu3,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu4,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu3+mu4)==False)[0]
    E['Coast']['Net Growth'],txt=gu.GetRegStats(mu3[ikp],mu4[ikp])
    #E[1]=np.sqrt(np.sum((mu3[ikp]-mu4[ikp])**2))
    #E[1]=np.nanmedian((mu4[ikp]-mu3[ikp])/mu3[ikp]*100)
    
    ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Mod C_M_Harv']==0) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt[iScn]['Age Med t0'][ind]
    xM=gpt[iScn]['Mod A t0'][ind]
    yO=gpt[iScn]['Ctot L t0'][ind]
    yM=gpt[iScn]['Mod C_Biomass_Tot t0'][ind]
    N,mu5,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu6,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu5+mu6)==False)[0]
    E['Interior']['Biomass'],txt=gu.GetRegStats(mu5[ikp],mu6[ikp])
    #E[2]=np.sqrt(np.sum((mu5[ikp]-mu6[ikp])**2))
    #E[2]=np.nanmedian((mu5[ikp]-mu4[ikp])/mu4[ikp]*100)
    
    yO=gpt[iScn]['Ctot Net'][ind]
    yM=gpt[iScn]['Mod C_G_Net'][ind]
    N,mu7,med,sig,se=gu.discres(xO,yO,bw,bin)
    N,mu8,med,sig,se=gu.discres(xM,yM,bw,bin)
    ikp=np.where(np.isnan(mu7+mu8)==False)[0]
    E['Interior']['Net Growth'],txt=gu.GetRegStats(mu7[ikp],mu8[ikp])
    #E[3]=np.sqrt(np.sum((mu7[ikp]-mu8[ikp])**2))
    #E[3]=np.nanmedian((mu8[ikp]-mu7[ikp])/mu7[ikp]*100)
    
    if 'Just get errors' not in meta[pNam].keys():            
        lw=0.5; ms=3; cl=meta['Graphics']['GP Comp']
        plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,6))
        # Coast    
        ax1[0].plot(bin,mu1,'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],label='Observed biomass',zorder=1)    
        ax1[0].plot(bin,mu2,'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],label='Predicted biomass',zorder=1)
        ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,50),xlim=[0,250+bw],ylim=[0,400])
        ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])        
        ax2=ax1[0].twinx()
        ax2.bar(bin-(0.45*bw/2),mu3,0.4*bw,ec='none',fc=cl['bl'],zorder=-1)    
        ax2.bar(bin+(0.45*bw/2),mu4,0.45*bw,ec='none',fc=cl['gl'],zorder=-1)
        ax2.plot([0,500],[0,0],'-k',lw=lw)
        ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
        ax1[0].set_zorder(ax2.get_zorder()+1)
        ax1[0].patch.set_visible(False)
        ax2.tick_params(length=meta['Graphics']['gp']['tickl'])                
        # Interior    
        ax1[1].plot(bin,mu5,'-ko',ms=ms,lw=lw,mew=lw,color=cl['bd'],mfc='w',mec=cl['bd'],zorder=1)    
        ax1[1].plot(bin,mu6,'--ks',ms=ms,lw=lw,mew=lw,color=cl['gd'],mfc='w',mec=cl['gd'],zorder=1)
        ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,20),xlim=[0,250+bw],ylim=[0,200])
        ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])        
        ax3=ax1[1].twinx()    
        ax3.bar(bin-(0.45*bw/2),mu7,0.4*bw,ec='none',fc=cl['bl'],label='Observed net growth',zorder=-1)    
        ax3.bar(bin+(0.45*bw/2),mu8,0.45*bw,ec='none',fc=cl['gl'],label='Predicted net growth',zorder=-1)    
        ax3.plot([0,500],[0,0],'-k',lw=lw)
        ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
        ax1[1].set_zorder(ax3.get_zorder()+1)
        ax1[1].patch.set_visible(False)
        ax3.tick_params(length=meta['Graphics']['gp']['tickl'])
        ax1[0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
        ax3.legend(loc='upper center',frameon=False,facecolor=None,edgecolor='w',fontsize=6)
        gu.axletters(ax1,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
        plt.tight_layout()
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseBiomassNetGrowth_CN_' + str(iScn+1),'png',900)
    
    return E

#%%
def EvalAtPlots_AgeResponsesGrossGrowthAndMortality_ByReg(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        bw=25; bin=np.arange(bw,250+bw,bw)
        cl=meta['Graphics']['GP Comp']
        plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,6.25))
        # Coast
        ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
        xO=gpt[iScn]['Age Med t0'][ind]
        xM=gpt[iScn]['Mod A t0'][ind]
        yO=gpt[iScn]['Ctot G Tot'][ind]
        yM=gpt[iScn]['Mod C_G_Gross_Tot'][ind]
        N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
        ax1[0].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=cl['bl'],label='Observed gross growth')
        ax1[0].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
        ax1[0].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=cl['gl'],label='Predicted gross growth')
        ax1[0].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        
        xO=gpt[iScn]['Age Med t0'][ind]
        xM=gpt[iScn]['Mod A t0'][ind]
        yO=gpt[iScn]['Ctot Mort+Lost'][ind]
        yM=gpt[iScn]['Mod C_M_Tot'][ind]
        N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
        ax1[0].bar(bin-(0.45*bw/2),-mu,0.4*bw,ec='none',fc=cl['bl'],label='Observed mortality')
        ax1[0].errorbar(bin-(0.45*bw/2),-mu,yerr=se,color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
        ax1[0].bar(bin+(0.45*bw/2),-mu,0.4*bw,ec='none',fc=cl['gl'],label='Predicted mortality')
        ax1[0].errorbar(bin+(0.45*bw/2),-mu,yerr=se,color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax1[0].plot([0,300],[0,0],'k-',lw=0.5)
        ax1[0].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-10,500,1),xlim=[0,250+bw],ylim=[-7,8])
        ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
        ax1[0].legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25,fontsize=6)
        # Interior
        ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
        xO=gpt[iScn]['Age Med t0'][ind]
        xM=gpt[iScn]['Mod A t0'][ind]
        yO=gpt[iScn]['Ctot G Tot'][ind]
        yM=gpt[iScn]['Mod C_G_Gross_Tot'][ind]
        N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
        ax1[1].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=cl['bl'])
        ax1[1].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
        ax1[1].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=cl['gl'])
        ax1[1].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        
        xO=gpt[iScn]['Age Med t0'][ind]
        xM=gpt[iScn]['Mod A t0'][ind]
        yO=gpt[iScn]['Ctot Mort+Lost'][ind]
        yM=gpt[iScn]['Mod C_M_Tot'][ind]
        N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
        ax1[1].bar(bin-(0.45*bw/2),-mu,0.4*bw,ec='w',fc=cl['bl'])
        ax1[1].errorbar(bin-(0.45*bw/2),-mu,yerr=se,color=0.5*cl['bl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
        ax1[1].bar(bin+(0.45*bw/2),-mu,0.4*bw,ec='w',fc=cl['gl'])
        ax1[1].errorbar(bin+(0.45*bw/2),-mu,yerr=se,color=0.5*cl['gl'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax1[1].plot([0,300],[0,0],'k-',lw=0.5)    
        ax1[1].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-10,500,1),xlim=[0,250+bw],ylim=[-6,4])
        ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
        gu.axletters(ax1,plt,0.04,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
        plt.tight_layout()
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseGrossGrowthAndMortality_CN_' + str(iScn+1),'png',900)
    return

# def EvalAgeResponsesMortality_ByReg(meta,pNam,gpt):
#     for iScn in range(meta[pNam]['Project']['N Scenario']):
#         bw=25; bin=np.arange(bw,250+bw,bw)
#         clO=np.array([0.8,0.85,1]); clM=np.array([0.85,1,0.8]);
#         plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,5.5))
#         # Coast
#         ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
#         xO=gpt[iScn]['Age Med t0'][ind]
#         xM=gpt[iScn]['Mod A t0'][ind]
#         yO=gpt[iScn]['Ctot Mort+Lost'][ind]
#         yM=gpt[iScn]['Mod C_M_Tot'][ind]
#         N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
#         ax1[0].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
#         ax1[0].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
#         N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
#         ax1[0].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
#         ax1[0].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
#         ax1[0].set(ylabel='Mortality (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,8])
#         ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
#         # Interior
#         ind=np.where( (gpt[iScn]['PTF CN']==1) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt[iScn]['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
#         xO=gpt[iScn]['Age Med t0'][ind]
#         xM=gpt[iScn]['Mod A t0'][ind]
#         yO=gpt[iScn]['Ctot Mort+Lost'][ind]
#         yM=gpt[iScn]['Mod C_M_Tot'][ind]
#         N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
#         ax1[1].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
#         ax1[1].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
#         N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
#         ax1[1].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
#         ax1[1].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
#         ax1[1].set(ylabel='Mortality (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,5])
#         ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
#         gu.axletters(ax1,plt,0.04,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
#         plt.tight_layout()
#         if meta['Graphics']['Print Figures']=='On':
#             gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseMortality_CN_' + str(iScn+1),'png',900)
#     return

#%%

def QA_Mortality_Regression(meta,pNam,gpt):
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        # Modelled harvest mortality is too high because salvage harvest is being counted as mortality
        # even though something else killed most of the trees.
        gpt2=copy.deepcopy(gpt)
        for k in gpt2.keys():
            try:
                gpt2[k]=gpt2[k][(gpt[iScn]['PTF CNY']==1)]
            except:
                pass
        del gpt2['vaL'],gpt2['soils']
    
        df=pd.DataFrame.from_dict(gpt2)
        df.columns=df.columns.str.replace(' ','_')
        #x=df[['debt_ratio','industry']]
        #y=df['cash_flow']
    
        #frm='Ctot_Mort ~ Age_t0 + C(IBB_Occurrence) + C(IBD_Occurrence) + C(IBM_Occurrence) + C(IBS_Occurrence) + C(IDW_Occurrence) + C(Wildfire_Occurrence) + C(Harv_Occurrence)'
        frm1='Ctot_Mort ~ Ctot_L_t0 + C(Occ_IBM) + C(Occ_Wildfire) + C(Occ_Harv)' #  + C(Occ_IBB) + C(Occ_IBS) + C(Occ_IDW)
        md1=smf.ols(formula=frm1,data=df)
        mr1=md1.fit()
        print(mr1.summary())
    
        frm2='Mod_C_M_Tot ~ Mod_C_Biomass_Tot_t0 + C(Occ_IBM) + C(Occ_Wildfire) + C(Occ_Harv)' #  + C(Occ_IBB)  + C(Occ_IBS) + C(Occ_IDW)
        md2=smf.ols(formula=frm2,data=df)
        mr2=md2.fit()
        print(mr2.summary())
    
        b1=mr1.params.values
        b2=mr2.params.values
        se1=mr1.bse.values
        se2=mr2.bse.values
        lab=['IBM','Wildfire','Harvest'] #'Int', ,'Initial\nbiomass''IBB',,'IBS','IDW'
    
        # Correct harvest to account for salvage
        b2[-2]=0.65*b2[-2]
    
        barw=0.32
        cl=meta['Graphics']['GP Comp']
    
        plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,9))
        ax.plot([-1,b1.size+1],[0,0],'k-',lw=0.5)
        ax.bar(np.arange(b1.size)-barw/2-0.01,b1,barw,facecolor=cl['bl'],label='Observations')
        ax.bar(np.arange(b1.size)+barw/2+0.01,b2,barw,facecolor=cl['gl'],label='Predictions (FCS)')
        ax.errorbar(np.arange(b1.size)-barw/2-0.01,b1,yerr=se1,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.errorbar(np.arange(b1.size)+barw/2+0.01,b2,yerr=se2,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
        ax.set(xticks=np.arange(1,b1.size-1,1),xticklabels=lab,ylabel='Regression coefficient (tC ha$^{-1}$ yr$^{-1}$)',
               xlim=[-0.5+1,b1.size-1-0.5],ylim=[-4,10.5])
        plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        plt.tight_layout()
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_Regression_' + str(iScn+1),'png',900)
    return

#%%

def Prepare_ObsVsModComparison(meta,pNam):

    gptS=[None]*meta[pNam]['Project']['N Scenario']
    for iScn in range(meta[pNam]['Project']['N Scenario']):
        
        # Import ground plot data
        metaGP,gpt=ugp.ImportGroundPlotData(meta,type='Stand')
        gpt['Ctot G Tot']=gpt['Ctot G Surv']+gpt['Ctot G Recr']
        gpt['Csw Net']=gpt['Csw G Surv']+gpt['Csw G Recr']-gpt['Csw Mort']
    
        # Import BC1ha data
        #z=u1ha.Import_Raster(meta,[],['refg','lcc1_c','harv_yr_con1'])
    
        # Import modelled data
        tvSaved=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
    
        vmL=['A','C_Biomass_Tot','C_Stemwood_Tot','C_Foliage_Tot','C_Branch_Tot',
             'C_Bark_Tot','C_Root_Tot','C_DeadWood_Tot','C_G_Gross_Tot','C_G_Net_Tot',
             'C_M_Reg_Tot','C_M_Dist','C_Soil_Tot','C_Soil_OHorizon','C_M_Harvest_Salvage','C_M_Harv']
        dM={}
        for v in vmL:
            dM[v]=np.zeros( (tvSaved.size,meta[pNam]['Project']['N Stand']) )        
        for iEns in range(meta[pNam]['Project']['N Ensemble']):
            for iBat in range(meta[pNam]['Project']['N Batch']):
                indBat=cbu.IndexToBatch(meta[pNam],iBat)
                d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
                for v in vmL:
                    if v=='C_M_Harvest_Salvage':
                        id=meta['LUT']['Event']['Harvest Salvage']
                        idx=d1['C_M_ByAgent'][id]['idx']
                        d1[v]=np.zeros( (tvSaved.size,indBat.size),dtype='float32')
                        d1[v][idx[0],idx[1]]=meta['Core']['Scale Factor C_M_ByAgent']*d1['C_M_ByAgent'][id]['M'].astype('float32')
                        dM[v][:,indBat]=dM[v][:,indBat]+d1[v].copy()
                    elif v=='C_M_Harv':                    
                        d1[v]=np.zeros( (tvSaved.size,indBat.size),dtype='float32')
                        id=meta['LUT']['Event']['Harvest Salvage']
                        idx=d1['C_M_ByAgent'][id]['idx']                    
                        d1[v][idx[0],idx[1]]=meta['Core']['Scale Factor C_M_ByAgent']*d1['C_M_ByAgent'][id]['M'].astype('float32')
                        dM[v][:,indBat]=dM[v][:,indBat]+d1[v].copy()
                        d1[v]=np.zeros( (tvSaved.size,indBat.size),dtype='float32')
                        id=meta['LUT']['Event']['Harvest']
                        idx=d1['C_M_ByAgent'][id]['idx']                    
                        d1[v][idx[0],idx[1]]=meta['Core']['Scale Factor C_M_ByAgent']*d1['C_M_ByAgent'][id]['M'].astype('float32')
                        dM[v][:,indBat]=dM[v][:,indBat]+d1[v].copy()
                    else:
                        dM[v][:,indBat]=dM[v][:,indBat]+d1[v]
                        
        for v in vmL:
            dM[v]=dM[v]/meta[pNam]['Project']['N Ensemble']
    
        # Calculate total mortality
        dM['C_M_Tot']=dM['C_M_Reg_Tot']+dM['C_M_Dist']-dM['C_M_Harvest_Salvage']
        dM['C_M_Nat']=dM['C_M_Reg_Tot']+dM['C_M_Dist']-dM['C_M_Harv']
        dM['C_G_Net']=dM['C_G_Gross_Tot']-dM['C_M_Tot']
    
        # Add simulations to gpt structure
        for k in dM.keys():
            gpt['Mod ' + k + ' t0']=np.nan*np.ones(gpt['Year t0'].size)
            gpt['Mod ' + k + ' t1']=np.nan*np.ones(gpt['Year t0'].size)
            gpt['Mod ' + k]=np.nan*np.ones(gpt['Year t0'].size)
    
        for i in range(gpt['X'].size):
            d=np.sqrt((meta['Geos']['Sparse']['X']-gpt['X'][i])**2+(meta['Geos']['Sparse']['Y']-gpt['Y'][i])**2)
            iS=np.where(d==np.min(d))[0]
            if iS.size==0:
                continue
            if np.isnan(gpt['Year t0'][i])==True:
                continue
            if np.isnan(gpt['Year t1'][i])==True:
                iT=np.where( (tvSaved==gpt['Year t0'][i]) )[0]
                gpt['Mod A t0'][i]=dM['A'][iT,iS]
                gpt['Mod C_Biomass_Tot t0'][i]=dM['C_Biomass_Tot'][iT,iS]
                gpt['Mod C_Stemwood_Tot t0'][i]=dM['C_Stemwood_Tot'][iT,iS]
                gpt['Mod C_Foliage_Tot t0'][i]=dM['C_Foliage_Tot'][iT,iS]
                gpt['Mod C_Branch_Tot t0'][i]=dM['C_Branch_Tot'][iT,iS]
                gpt['Mod C_Bark_Tot t0'][i]=dM['C_Bark_Tot'][iT,iS]
                gpt['Mod C_Root_Tot t0'][i]=dM['C_Root_Tot'][iT,iS]
                gpt['Mod C_DeadWood_Tot t0'][i]=dM['C_DeadWood_Tot'][iT,iS]
            else:
                iT=np.where( (tvSaved>=gpt['Year t0'][i]) & (tvSaved<=gpt['Year t1'][i]) )[0]
                gpt['Mod A t0'][i]=dM['A'][iT[0],iS]
                gpt['Mod C_Biomass_Tot t0'][i]=dM['C_Biomass_Tot'][iT[0],iS]
                gpt['Mod C_Stemwood_Tot t0'][i]=dM['C_Stemwood_Tot'][iT[0],iS]
                gpt['Mod C_Foliage_Tot t0'][i]=dM['C_Foliage_Tot'][iT[0],iS]
                gpt['Mod C_Branch_Tot t0'][i]=dM['C_Branch_Tot'][iT[0],iS]
                gpt['Mod C_Bark_Tot t0'][i]=dM['C_Bark_Tot'][iT[0],iS]
                gpt['Mod C_Root_Tot t0'][i]=dM['C_Root_Tot'][iT[0],iS]
                gpt['Mod C_Biomass_Tot t1'][i]=dM['C_Biomass_Tot'][iT[-1],iS]
                gpt['Mod C_DeadWood_Tot t0'][i]=dM['C_DeadWood_Tot'][iT[0],iS]
                gpt['Mod C_DeadWood_Tot t1'][i]=dM['C_DeadWood_Tot'][iT[-1],iS]
                gpt['Mod C_M_Tot'][i]=np.mean(dM['C_M_Tot'][iT,iS])
                gpt['Mod C_M_Harv'][i]=np.mean(dM['C_M_Harv'][iT,iS])
                gpt['Mod C_M_Nat'][i]=np.mean(dM['C_M_Nat'][iT,iS])
                gpt['Mod C_G_Gross_Tot'][i]=np.mean(dM['C_G_Gross_Tot'][iT,iS])
                gpt['Mod C_G_Net'][i]=np.mean(dM['C_G_Net'][iT,iS])
    
        # Import TIPSY stemwood biomass
        iGC=0
        gc=cbu.Import_BatchTIPSY_Output(meta,pNam,iScn,iGC)
    
        gpt['GY Csw t0']=np.nan*np.ones( gpt['Year t0'].size )
        gpt['GY Csw Net']=np.nan*np.ones( gpt['Year t0'].size )
        for iStand in range(meta['Geos']['Sparse']['X'].size):
            iS=np.where( (np.abs(gpt['X']-meta['Geos']['Sparse']['X'][iStand])<=100) & (np.abs(gpt['Y']-meta['Geos']['Sparse']['Y'][iStand])<=100) )[0]
            if iS.size==0:
                continue
            for j in range(iS.size):
                if (np.isnan(gpt['Age VRI t0'][iS[j]])==True) | (np.isnan(gpt['Delta t'][iS[j]])==True):
                    continue
                iA=np.minimum(200.,gpt['Age VRI t0'][iS[j]]).astype(int)
                gpt['GY Csw t0'][iS[j]]=gc['Csw'][iA,iStand]
                if (np.isnan(gpt['Year t1'][iS[j]])==False) & (iA<198):
                    iA2=np.arange(iA,np.minimum(200,iA+gpt['Delta t'][iS[j]]+1),1).astype('int')
                    gpt['GY Csw Net'][iS[j]]=np.mean(gc['Gsw_Net'][iA2,iStand])
    
        # List of variables for analysis
        gpt['vaL']=['Age VRI t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0','Cdw t0',
            'Ctot G Tot','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Mort Nat','Ctot Mort Harv','Ctot Net',
            'Mod A t0','Mod C_Biomass_Tot t0','Mod C_Biomass_Tot t1','Mod C_M_Tot','Mod C_G_Gross_Tot','Mod C_G_Net',
            'Mod C_M_Nat','Mod C_M_Harv',
            'Mod C_Stemwood_Tot t0','Mod C_Foliage_Tot t0','Mod C_Branch_Tot t0','Mod C_Bark_Tot t0','Mod C_Root_Tot t0',
            'Csw Net','GY Csw t0','GY Csw Net']
    
        # Soil organic carbon
    
        # Import soils (see soil_Shawetal2018_01_process.py
        soils=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')
    
        # Add simulations to soils data
        iT=np.where(tvSaved==2020)[0]
        soils['C_Soil_Tot']=np.nan*np.ones(soils['x'].size)
        soils['C_Soil_OHorizon']=np.nan*np.ones(soils['x'].size)
        cnt=0
        for i in range(meta['Geos']['Sparse']['X'].size):
            iS=np.where( (np.abs(soils['x']-meta['Geos']['Sparse']['X'][i])<=2000) & (np.abs(soils['y']-meta['Geos']['Sparse']['Y'][i])<=2000) )[0]
            if iS.size==0:
                continue
            soils['C_Soil_Tot'][iS]=dM['C_Soil_Tot'][iT,i]
            soils['C_Soil_OHorizon'][iS]=dM['C_Soil_OHorizon'][iT,i]
            cnt=cnt+1
    
        gpt['soils']=soils
    
        # Ground plots
        gpt['soils']['vaL']=['TOT_C_THA','MIN_C_THA','ORG_C_THA','C_Soil_Tot','C_Soil_OHorizon']
    
        # Add simulations to soils data
        iT=np.where(tvSaved==2020)[0]
        gpt['soils']['C_Soil_Tot']=np.nan*np.ones(gpt['soils']['x'].size)
        gpt['soils']['C_Soil_OHorizon']=np.nan*np.ones(gpt['soils']['x'].size)
        cnt=0
        for i in range(meta['Geos']['Sparse']['X'].size):
            iS=np.where( (np.abs(gpt['soils']['x']-meta['Geos']['Sparse']['X'][i])<=2000) & (np.abs(gpt['soils']['y']-meta['Geos']['Sparse']['Y'][i])<=2000) )[0]
            if iS.size==0:
                continue
            gpt['soils']['C_Soil_Tot'][iS]=dM['C_Soil_Tot'][iT,i]
            gpt['soils']['C_Soil_OHorizon'][iS]=dM['C_Soil_OHorizon'][iT,i]
            cnt=cnt+1
        
        # Add to scenarios
        gptS[iScn]=gpt
    
    # Save    
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\EvaluationAtGroundPlots.pkl',gptS)

    return

#%%
def QA_ProfileIBM(meta,pNam,gpt,iScn):

    # Import AOS data
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    ind=gis.GetGridIndexToPoints(zRef,gpt[iScn]['X'],gpt[iScn]['Y'])
    tvMPB=np.arange(1995,2022,1)
    ios=np.zeros( (tvMPB.size,ind[0].size) ,dtype='int8')
    for iT in range(tvMPB.size):
        z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(tvMPB[iT]) + '.tif')
        ios[iT,:]=z['Data'][ind].copy()
    
    # Calculate average outbreak profiles
    
    ios_mx=np.max(ios,axis=0)
    ind=np.where(ios_mx>=2)[0]
    sts={}
    sts['tso']=np.zeros(ind.size)
    for v in gpt[iScn].keys():
        if (v=='vaL') | (v=='soils'):
            continue
        sts[v]=np.zeros(ind.size)
        for iStand in range(ind.size):
            indOutbreak=np.where(ios[:,ind[iStand]]>0)[0]
            sts['tso'][iStand]=gpt[iScn]['Year t0'][ind[iStand]]-tvMPB[indOutbreak[0]]
            sts[v][iStand]=gpt[iScn][v][ind[iStand]]
    
    bw=4; bin=np.arange(-10,20,bw)
    prfl={}
    for k in sts.keys():
        prfl[k]={}
        prfl[k]['N'],prfl[k]['mu'],prfl[k]['med'],prfl[k]['sig'],prfl[k]['se']=gu.discres(sts['tso'],sts[k],bw,bin)
    
    # Plot
    
    it0=np.where(bin<0)[0]
    it1=np.where(bin>=10)[0]
    it2=np.where(bin>0)[0]
    
    # Plot observed profile
    
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,7)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
    ax[0].plot([0,0],[0,230],'k-',lw=2,color=[0.75,0.75,0.75])
    ax[0].plot(bin,prfl['Ctot L t0']['mu']+prfl['Cdw t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')
    
    mu0=np.mean(prfl['Ctot L t0']['mu'][it0])
    mu1=np.mean(prfl['Ctot L t0']['mu'][it1])
    ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
    ax[0].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
    ax[0].text(np.mean(bin[it1]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=6,style='normal',weight='bold',color='k')
    
    mu0=np.mean(prfl['Cdw t0']['mu'][it0])
    ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
    it1b=np.where(prfl['Cdw t0']['mu']==np.max(prfl['Cdw t0']['mu']) )[0]
    mu1=prfl['Cdw t0']['mu'][it1b[0]]
    ax[0].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
    ax[0].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=6,style='normal',weight='bold',color='k')
    
    ax[0].plot(bin,prfl['Ctot L t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
    ax[0].plot(bin,prfl['Cdw t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')
    
    ax[0].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
    ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
    
    # Plot modelled profile
    #fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,7.5)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
    ax[1].plot([0,0],[0,250],'k-',lw=2,color=[0.75,0.75,0.75])
    ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu']+prfl['Mod C_DeadWood_Tot t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')
    
    mu0=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it0])
    mu1=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it1])
    ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
    ax[1].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
    ax[1].text(np.mean(bin[it1]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=6,style='normal',weight='bold',color='k')
    
    mu0=np.mean(prfl['Mod C_DeadWood_Tot t0']['mu'][it0])
    ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
    it1b=np.where(prfl['Mod C_DeadWood_Tot t0']['mu']==np.max(prfl['Mod C_DeadWood_Tot t0']['mu']) )[0]
    mu1=prfl['Mod C_DeadWood_Tot t0']['mu'][it1b[0]]
    ax[1].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
    ax[1].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=6,style='normal',weight='bold',color='k')
    
    ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
    ax[1].plot(bin,prfl['Mod C_DeadWood_Tot t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')
    
    ax[1].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
    ax[1].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
    ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
    
    gu.axletters(ax,plt,0.045,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
    plt.tight_layout()
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Profile_IBM_' + str(iScn+1),'png',900)
    return

#%% Fire profiles
def QA_ProfileWildfire(meta,pNam,gpt,iScn):
    zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
    iGrd=gis.GetGridIndexToPoints(zRef,gpt[iScn]['X'],gpt[iScn]['Y'])
    tvFire=np.arange(1990,2022,1)
    yrFire=np.zeros( (tvFire.size,iGrd[0].size) ,dtype='int16')
    for iT in range(tvFire.size):
        z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tvFire[iT]) + '.tif')['Data']
        yrFire[iT,:]=z[iGrd]
    
    yrFire_mx=np.max(yrFire,axis=0)
    ind=np.where(yrFire_mx>0)[0]
    sts={}
    sts['tsf']=np.zeros(ind.size)
    for v in gpt[iScn].keys():
        if (v=='vaL') | (v=='soils'):
            continue
        sts[v]=np.zeros(ind.size)
        for iStand in range(ind.size):
            indO=np.where(yrFire[:,ind[iStand]]>0)[0]
            sts['tsf'][iStand]=gpt[iScn]['Year t0'][ind[iStand]]-tvFire[indO[0]]
            sts[v][iStand]=gpt[iScn][v][ind[iStand]]
    
    bw=4; bin=np.arange(-10,20,bw)
    prfl={}
    for k in sts.keys():
        prfl[k]={}
        prfl[k]['N'],prfl[k]['mu'],prfl[k]['med'],prfl[k]['sig'],prfl[k]['se']=gu.discres(sts['tsf'],sts[k],bw,bin)
    
    # Plot
    
    it0=np.where(bin<0)[0]
    it1=np.where(bin>0)[0]
    it2=np.where(bin>0)[0]
    
    # Plot observed profile
    
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,7)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
    ax[0].plot([0,0],[0,130],'k-',lw=2,color=[0.75,0.75,0.75])
    #ax[0].plot(bin,prfl['Ctot L t0']['mu']+prfl['Cdw t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')
    
    mu0=np.mean(prfl['Ctot L t0']['mu'][it0])
    mu1=np.mean(prfl['Ctot L t0']['mu'][it1])
    ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
    ax[0].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
    ax[0].text(np.mean(bin[it1]),mu1+12,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[1,:])
    
    mu0=np.mean(prfl['Cdw t0']['mu'][it0])
    ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
    it1b=np.where(prfl['Cdw t0']['mu']==np.max(prfl['Cdw t0']['mu']) )[0]
    mu1=prfl['Cdw t0']['mu'][it1b[0]]
    #ax[0].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
    #ax[0].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[2,:])
    
    ax[0].plot(bin,prfl['Ctot L t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
    ax[0].plot(bin,prfl['Cdw t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')
    
    ax[0].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
    ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=gp['tickl'])
    
    # Plot modelled profile
    #fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,7.5)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
    ax[1].plot([0,0],[0,250],'k-',lw=2,color=[0.75,0.75,0.75])
    ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu']+prfl['Mod C_DeadWood_Tot t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')
    
    mu0=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it0])
    mu1=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it1])
    ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
    ax[1].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
    ax[1].text(np.mean(bin[it1]),mu1+12,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[1,:])
    
    mu0=np.mean(prfl['Mod C_DeadWood_Tot t0']['mu'][it0])
    ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
    it1b=np.where(prfl['Mod C_DeadWood_Tot t0']['mu']==np.max(prfl['Mod C_DeadWood_Tot t0']['mu']) )[0]
    mu1=prfl['Mod C_DeadWood_Tot t0']['mu'][it1b[0]]
    #ax[1].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
    #ax[1].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[2,:])
    
    ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
    ax[1].plot(bin,prfl['Mod C_DeadWood_Tot t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')
    
    ax[1].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
    ax[1].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
    ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=gp['tickl'])
    
    gu.axletters(ax,plt,0.04,0.92,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')
    plt.tight_layout()
    gu.PrintFig(meta['Paths']['Figures'] + '\\QA_Wildfire_Profiles_' + str(iScn+1),'png',900)
    return