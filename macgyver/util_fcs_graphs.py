
#%% Import modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.patches import Rectangle
import gc as garc
import warnings
import time
import copy
import docx
import matplotlib.patches as patches
from IPython.display import display,HTML,Image
from IPython.display import Markdown as md
import graphviz
from diagrams import Diagram,Cluster,Edge
from diagrams import aws
from diagrams.aws import analytics as aws_analytics
from diagrams.aws import network as aws_network
from diagrams.aws import compute as aws_compute
from diagrams.aws import database as aws_database
from diagrams.onprem import vcs
import matplotlib.colors
import matplotlib as mpl
import matplotlib.ticker as ticker
#from matplotlib import animation
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_utils as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcgadgets.macgyver.util_fcs_qa as uqa
import fcgadgets.hardhat.harvest_util as hut
import fcexplore.field_plots.processing.fp_util as ufp

#%% Update counters
def UpdateTabCount(meta):
	meta['Graphics']['Tab Count']=meta['Graphics']['Tab Count']+1
	return meta
def UpdateFigCount(meta):
	meta['Graphics']['Fig Count']=meta['Graphics']['Fig Count']+1
	return meta
#%% Add caption
def FigureCaption(meta,caption):
	txt='<b>Figure ' + str(meta['Graphics']['Fig Count']) + '</b>. ' + caption
	meta=UpdateFigCount(meta)
	display(md(txt))
	return meta
def TableCaption(meta,caption):
	txt='<b>Table ' + str(meta['Graphics']['Tab Count']) + '</b>. ' + caption
	meta=UpdateTabCount(meta)
	display(md(txt))
	return meta
#%% Table
def Table(x):
	try:
		df=pd.read_excel(x)
	except:
		df=pd.DataFrame(x)
	df=df.fillna("")
	df=df.round(decimals=2)
	df=df.style.set_table_styles([dict(selector='th',props=[('text-align','left')])])
	df.set_properties(**{'text-align':'left'}).hide()
	df.set_properties(**{'text-align':'left'},**{'width':'300px'}).hide()
	display(df)
	return

#%%
def RunOutputGraphics(meta,pNam,lsat,mos,dR,iPS,iSS,iYS,tv,iT,gpt,dH_ByBGCZ):
	
	# Comparisons
	for c in mos[pNam]['Delta'].keys():
		if dR['ScenarioComp_GHGBalanceAnnual']=='On':
			ScenarioComp_GHGBalanceAnnual(meta,pNam,mos,tv,c,iT,iPS,iSS,iYS)
		if dR['ScenarioComp_GHGBalanceAnnualPlusCumu']=='On':
			ScenarioComp_GHGBalanceAnnualPlusCumu(meta,pNam,mos,c,iT,iPS,iSS,iYS)
		if dR['ScenarioComp_ForcingBarChart']=='On':
			ScenarioComp_ForcingBarChart(meta,pNam,mos,tv,c,iPS,iSS,iYS)
	
	E=[None]*meta[pNam]['Project']['N Scenario']
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		
		# Areas
		if dR['Area_DisturbedAndManaged']=='On':
			ivlT=1
			Area_DisturbedAndManaged(meta,pNam,mos,ivlT,iScn,iT,iPS,iSS,iYS)
		if dR['Area_HarvestTimeSeries']=='On':
			Area_HarvestTimeSeries(meta,pNam,mos,tv,iScn,iPS,iSS,iYS)
		if dR['Area_THLB_TimeSeries']=='On':
			#iScn=meta[pNam]['Project']['Actual Indices'][0]
			Area_THLB_TimeSeries(meta,pNam,tv,lsat,iScn)
			
		# GHG / Carbon
		if dR['GHGBalance']=='On':
			GHGBalance(meta,mos,pNam,tv,iScn,iT,iPS,iSS,iYS)
		if dR['GHGBalanceSimple']=='On':
			GHGBalanceSimple(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS)
		if dR['CarbonFluxTimeSeries']=='On':
			CarbonFluxTimeSeries(meta,pNam,mos,iT,iScn,iPS,iSS,iYS)
		if dR['CarbonPoolBarChart']=='On':
			CarbonPoolBarChart(meta,pNam,mos,tv,iScn,iPS,iSS,iYS)
		if dR['CarbonFluxesBarChart']=='On':
			CarbonFluxesBarChart(meta,pNam,mos,tv,iScn,iPS,iSS,iYS)
		if dR['CarbonFluxMortality']=='On':
			ivlT=1
			CarbonFluxMortality(meta,pNam,mos,tv,iScn,iT,ivlT,iPS,iSS,iYS)
		if dR['ComparisonWithPIR']=='On':
			#iScn=meta[pNam]['Project']['Actual Indices'][0]
			ComparisonWithPIR(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS)
		
		if dR['AgeClassDist']=='On':
			AgeClassDist(meta,pNam,iScn,iPS,iSS,iYS)
		if dR['HarvestVolumeTimeSeries']=='On':
			flg_ann=1; flg_dead=1; flg_hbs=1
			HarvestVolumeTimeSeries(meta,mos,pNam,tv,iScn,iT,iPS,iSS,iYS,flg_ann,flg_dead,flg_hbs,dH_ByBGCZ)
		if dR['HarvestVolumePerHectareTimeSeries']=='On':
			iT2=np.where( (tv>=1960) & (tv<=2025) )[0]
			#iScn=meta[pNam]['Project']['Actual Indices'][0]
			HarvestVolumePerHectareTimeSeries(meta,pNam,mos,tv,iScn,iT2,iPS,iSS,iYS,dH_ByBGCZ)
			#iScn=1
			#PlotVolumePerHectare(meta,mos,tv,iScn,iT2,iPS,iSS,iYS,AEF)
		if dR['MortalitySpectrum']=='On':
			iEns=0
			MortalitySpectrum(meta,pNam,iEns,iScn,iPS,iSS,iYS)
		if dR['NetGrowthTimeSeries']=='On':
			NetGrowthTimeSeries(meta,pNam,mos,tv,iT,iScn,iPS,iSS,iYS)
		if dR['WildfireProbHistTimeSeries']=='On':
			hw=WildfireProbHistTimeSeries(meta,pNam,iScn)
		if dR['AllVariableTimeSeries']=='On':
			PlotAllVariableTimeSeries(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS)
		#if dR['Summary of Approaches']=='On':
		#	SummarizeExistingInitiatives(meta)
		if dR['SaveTabData']=='On':
			t_start=1990
			t_end=2021
			sum_mult=meta[pNam]['Project']['AEF']/1e6
			df=udem.ExportSummariesByScenario(meta,pNam,mos,t_start,t_end,sum_mult=sum_mult)
		
		# QA
		E[iScn]={'Province':{},'Coast':{},'Interior':{}}
		if dR['EvalAtPlots_AgeResponsesBiomassAndNetGrowth_ByReg_CN']=='On':
			E[iScn]=uqa.EvalAtPlots_AgeResponsesBiomassAndNetGrowth_ByReg_CN(meta,pNam,iScn,gpt,E[iScn])
		if dR['EvalAtPlots_AgeResponsesGrossGrowthAndMortality_ByReg']=='On':
			uqa.EvalAtPlots_AgeResponsesGrossGrowthAndMortality_ByReg(meta,pNam,gpt)
		if dR['EvalAtPlots_BiomassDynamicsAve_CN']=='On':
			uqa.EvalAtPlots_BiomassDynamicsAve_CN(meta,pNam,gpt)
		if dR['EvalAtPlots_AgeByBGCZ_CNV']=='On':
			uqa.EvalAtPlots_AgeByBGCZ_CNV(meta,pNam,gpt)
		if dR['EvalAtPlots_BiomassByBGC_CNV']=='On':
			uqa.EvalAtPlots_BiomassByBGC_CNV(meta,pNam,gpt)
		if dR['EvalAtPlots_GrossGrowthByBGC_CN']=='On':
			uqa.EvalAtPlots_GrossGrowthByBGC_CN(meta,pNam,gpt)
		if dR['EvalAtPlots_MortalityByBGC_CN']=='On':
			uqa.EvalAtPlots_MortalityByBGC_CN(meta,pNam,gpt)
		if dR['EvalAtPlots_SOCByBGC_ShawComp']=='On':
			uqa.EvalAtPlots_SOCByBGC_ShawComp(meta,pNam,gpt)
		if dR['QA Profile IBM']=='On':
			uqa.QA_ProfileIBM(meta,pNam,gpt,iScn)
		if dR['QA Profile Wildfire']=='On':
			uqa.QA_ProfileWildfire(meta,pNam,gpt,iScn)
		
	return E

#%%
def Area_HarvestTimeSeries(meta,pNam,mos,tv,iScn,iPS,iSS,iYS):
	#dA=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestArea_TS.pkl')
	dA=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestTimeSeries.pkl')

	iT2=np.where( (tv>=1901) & (tv<=2100) )[0]
	ivlT=1
	A=cbu.SummarizeAreaAffected(meta,pNam,mos,tv,iScn,iPS,iSS,iYS,ivlT)
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,8));
	ax.plot(dA['Year'],dA['Province']['Area']['CC'],'k-',color=meta['Graphics']['gp']['cl1'],ms=2,label='Consolidated Cutblocks database')
	ax.plot(tv[iT2],A['Management']['Harvest']['Data'][iT2]/1e3,'k--',color=meta['Graphics']['gp']['cl2'],ms=2,label='Model')
	ax.set(xticks=np.arange(1800,2200,10),yticks=np.arange(0,400,50),ylabel='Harvest area (1000 hectares)',xlabel='Time, years',
		xlim=[1901,2100],ylim=[0,300])
	ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_All_All_All_Harvest_Area_TS','png',900)
	return 

#%% THLB
def Area_THLB_TimeSeries(meta,pNam,tv,lsat,iScn):

	xlim=[tv[0],tv[-1]]
	xt=np.arange(1800,2200,20)
	#sLCLU=meta[pNam]['Scenario'][iScn]['Land Cover/Land Use Scenario']
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
	ax.plot(meta[pNam]['Year'],meta[pNam]['Project']['AEF']*np.sum(lsat['THLB']['S0'],axis=1)/1e6,'k-',color=meta['Graphics']['gp']['cl1'],label='S0')
	ax.plot(meta[pNam]['Year'],meta[pNam]['Project']['AEF']*np.sum(lsat['THLB']['S1'],axis=1)/1e6,'k--',color=meta['Graphics']['gp']['cl2'],label='S1')
	ax.plot(meta[pNam]['Year'],meta[pNam]['Project']['AEF']*np.sum(lsat['THLB']['S2'],axis=1)/1e6,'k-.',color=meta['Graphics']['gp']['cl3'],label='S2')
	ax.set(xticks=xt,yticks=np.arange(0,26,2),ylabel='Timber harvesting landbase (million ha)',xlabel='Time, years',xlim=xlim,ylim=[0,28])
	ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn1_All_All_All_Area_THLB_TimeSeries','png',900)
	return
#%% Carbon pool summary
def CarbonPoolBarChart(meta,pNam,mos,tv,iScn,iPS,iSS,iYS):
	iT2=[]
	iT2.append(np.where( (tv>=1800) & (tv<=1820) )[0])
	iT2.append(np.where( (tv>=1990) & (tv<=2020) )[0])
	iT2.append(np.where( (tv>=2100) & (tv<=2120) )[0])

	bw=0.25; bw2=0.3;
	cl=np.array([[0.65,0.75,0.75],[0.7,0.8,0.9],[0.8,0.9,1]]).T
	lab=['1820','1990','2100']

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,6));

	# Biomass
	nam='C_Biomass'; x0=1; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
		ax[0].bar(x0,y[i],0.01,color=cl[i,:],label=lab[i])
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	# Dead wood
	nam='C_DeadWood'; x0=2; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	# LiT2ter
	nam='C_Litter'; x0=3; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	# Soil
	nam='C_Soil'; x0=4; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	# In-use products
	nam='C_InUse'; x0=5; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	# Dump and landfill
	nam='C_WasteSystems'; x0=6; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[0].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

	ax[0].set(position=[0.07,0.125,0.7,0.86],xticks=np.arange(1,7,1),xticklabels=['Biomass','Dead \nwood','Litter','Soil','In-use\nproducts','Dumps & \nlandfills'],xlabel='',xlim=[0.5,6.5], \
	  yticks=np.arange(0,120,2),ylabel='Carbon stock (GtC)',ylim=[0,11])
	ax[0].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])

	# Total
	nam='C_HWP'; x0=1; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
	for i in range(3):
		y[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum']['C_HWP']['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)+np.mean((mos[pNam]['Scenarios'][iScn]['Sum']['C_Forest']['Ensemble Mean'][iT2[i],iPS,iSS,iYS])/1e9)
		yl[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum']['C_HWP']['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)+np.mean((mos[pNam]['Scenarios'][iScn]['Sum']['C_Forest']['Ensemble P025'][iT2[i],iPS,iSS,iYS])/1e9)
		yh[i]=np.nanmean((mos[pNam]['Scenarios'][iScn]['Sum']['C_HWP']['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)+np.mean((mos[pNam]['Scenarios'][iScn]['Sum']['C_Forest']['Ensemble P975'][iT2[i],iPS,iSS,iYS])/1e9)
	ax[1].bar(x,y,bw,color=cl)
	if meta[pNam]['Project']['N Ensemble']>1:
		ax[1].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)
	ax[1].set(position=[0.84,0.125,0.14,0.86],xticks=np.arange(1,2,1),xticklabels=['Total'],yticks=np.arange(0,120,2),ylabel='Carbon stock (GtC)',xlabel='',ylim=[0,30],xlim=[0.5,1.5])
	ax[1].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])

	#ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_CarbonPoolBarChart','png',900)

	return

#%% flux summary
def CarbonFluxesBarChart(meta,pNam,mos,tv,iScn,iPS,iSS,iYS):
	#list(mos[pNam]['Scenarios'][0]['Sum'].keys())

	nam=np.array(['E_Domestic_ForestSector_NEE','E_Wildfire_ForestSector_Domestic','E_OpenBurning_ForestSector_Domestic','E_Domestic_ForestSector_HWP','E_Bioenergy','E_ForestOperations', \
				  'E_Substitution_Energy','E_Substitution_Material','E_Substitution_Total','E_NAB','E_NSB'])

	xlab=['Net\necosystem\nexchange','Wildfire\nemissions','Open\nburning\nemissions','Wood\nproduct\nemissions','Bioenergy\nemissions','Operational\nemissions', \
		  'Energy\nsubst.','Material\nsubst.','Total\nsubst.','Net GHG\nbalance\n(with subs.)','Net GHG\nbalance\n(w/o subs.)']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,8.5)); bw=0.2; bw2=0.25;
	ax.plot([0,nam.size+1],[0,0],'-',color=meta['Graphics']['gp']['cla'])

	iT2=np.where( (tv>=1865) & (tv<=2021) )[0]
	y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
	ax.bar(np.arange(1,y.size+1,1)-bw2,y,bw,color=[0.57,0.79,1],label='1865-2021')
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(np.arange(1,y.size+1,1)-bw2,y,yerr=[np.maximum(0,y-yl),np.maximum(0,yh-y)],color=[0.37,0.59,0.8],ls='',capsize=2)

	iT2=np.where( (tv>=1990) & (tv<=2021) )[0]
	y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
	ax.bar(np.arange(1,y.size+1,1),y,bw,color=[0.8,1,0.6],label='1990-2021')
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(np.arange(1,y.size+1,1),y,yerr=[np.maximum(0,y-yl),np.maximum(0,yh-y)],color=[0.6,0.8,0.4],ls='',capsize=2)

	iT2=np.where( (tv>=2022) & (tv<=2100) )[0]
	y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		y[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yl[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
		yh[i]=np.mean((mos[pNam]['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'])
	ax.bar(np.arange(1,y.size+1,1)+bw2,y,bw,color=[0.8,0.4,1],label='2023-2100')
	if meta[pNam]['Project']['N Ensemble']>1:
		try:
			ax.errorbar(np.arange(1,y.size+1,1)+bw2,y,yerr=[y-yl,yh-y],color=[0.6,0.2,0.8],ls='',capsize=2)
		except:
			pass

	ax.set(xlim=[0.5,nam.size+.5],xticks=np.arange(1,nam.size+1,1),xticklabels=xlab, \
		   ylabel='Mean GHG flux (MtCO$_2$e yr$^-$$^1$)') # ylim=[-70,70],yticks=np.arange(-100,100,10)
	ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	
	#ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_CarbonFluxesBarChart','png',900)
	return

#%% flux summary
def ScenarioComp_ForcingBarChart(meta,pNam,mos,tv,cNam,iPS,iSS,iYS):

	nam=np.array(['E_Domestic_ForestSector_NEE','E_Wildfire_ForestSector_Domestic','E_OpenBurning_ForestSector_Domestic','E_Domestic_ForestSector_HWP','E_Bioenergy','E_ForestOperations', \
				  'E_Substitution_Energy','E_Substitution_Material','E_Substitution_Total','E_NAB','E_NSB'])

	xlab=['Net\necosystem\nexchange','Wildfire\nemissions','Open\nburning\nemissions','Wood\nproduct\nemissions','Bioenergy\nemissions','Operational\nemissions', \
		  'Energy\nsubst.','Material\nsubst.','Total\nsubst.','Net GHG\nbalance\n(with subs.)','Net GHG\nbalance\n(w/o subs.)']

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,8.5)); bw=0.2; bw2=0.25;
	ax.plot([0,nam.size+1],[0,0],'k-',color=meta['Graphics']['gp']['cla'])

	iT2=np.where( (tv>=1865) & (tv<1990) )[0]
	ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		ymu[i]=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		ylo=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		yhi=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		cs=np.array([ylo,yhi])
		yh[i]=np.max(cs)
		yl[i]=np.min(cs)
	ax.bar(np.arange(1,ymu.size+1,1)-bw2,ymu,bw,color=[0.57,0.79,1],label='1865-1990')
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(np.arange(1,ymu.size+1,1)-bw2,ymu,yerr=[ymu-yl,yh-ymu],color=[0.37,0.59,0.8],ls='',capsize=2)

	iT2=np.where( (tv>=1990) & (tv<=2021) )[0]
	ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		ymu[i]=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		ylo=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		yhi=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		cs=np.array([ylo,yhi])
		yh[i]=np.max(cs)
		yl[i]=np.min(cs)
	ax.bar(np.arange(1,ymu.size+1,1),ymu,bw,color=[0.8,1,0.6],label='1991-2021')
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(np.arange(1,ymu.size+1,1),ymu,yerr=[np.maximum(0,ymu-yl),np.maximum(0,yh-ymu)],color=[0.6,0.8,0.4],ls='',capsize=2)

	iT2=np.where( (tv>=2022) & (tv<=2100) )[0]
	ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
	for i in range(nam.size):
		ymu[i]=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		ylo=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		yhi=np.mean(mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
		cs=np.array([ylo,yhi])
		yh[i]=np.max(cs)
		yl[i]=np.min(cs)
	ax.bar(np.arange(1,ymu.size+1,1)+bw2,ymu,bw,color=[0.8,0.4,1],label='2023-2100')
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(np.arange(1,ymu.size+1,1)+bw2,ymu,yerr=[ymu-yl,yh-ymu],color=[0.6,0.2,0.8],ls='',capsize=2)

	ax.set(xticks=np.arange(1,nam.size+1,1),xticklabels=xlab,ylabel='Forcing (MtCO$_2$e yr$^-$$^1$)',xlim=[0.5,nam.size+.5]) # ,ylim=[-100,100] yticks=np.arange(-100,200,10),
	ax.legend(loc='upper left',ncol=3,facecolor=[1,1,1],frameon=False);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()

	#ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_ScnComp_' + cNam + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_ForcingBarChart','png',900)

#%% Harvest removals
def HarvestVolumeTimeSeries(meta,mos,pNam,tv,iScn,iT,iPS,iSS,iYS,iOS,flg_ann,flg_dead,flg_hbs,dH_ByBGCZ):

	# Import FAIB public summary
	#H_FLNR=gu.ReadExcel(r'C:\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')
	dH=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestTimeSeries.pkl')

	H_HBS=gu.ipickle(r'C:\Data\Harvest\HBS\HBS_AnnualSummary.pkl')

	#matype='historical'
	matype='center'

	# Total harvest
	#y=(mos[pNam]['Scenarios'][iScn]['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,iPS,iSS,iYS]+
	#	mos[pNam]['Scenarios'][iScn]['Sum']['V_ToMill_NonMerchTotal']['Ensemble Mean'][:,iPS,iSS,iYS]+
	#	mos[pNam]['Scenarios'][iScn]['Sum']['V_ToMill_MerchDead']['Ensemble Mean'][:,iPS,iSS,iYS])/meta[pNam]['Project']['Multi']
	y_Live=cbu.GetMosScnVar(meta,pNam,mos,iScn,'V_ToMill_MerchTotal',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
	y_Dead=cbu.GetMosScnVar(meta,pNam,mos,iScn,'V_ToMill_MerchDead',iPS,iSS,iYS,iOS)['Sum']['Ensemble Mean']/meta[pNam]['Project']['Multi']
	y=y_Live+y_Dead
	y_ma=gu.movingave(y[iT],10,matype)

	#y_Dead=(mos[pNam]['Scenarios'][iScn]['Sum']['V_ToMill_MerchDead']['Ensemble Mean'][:,iPS,iSS,iYS])/meta[pNam]['Project']['Multi']
	y_Dead_ma=gu.movingave(y_Dead[iT],10,matype)

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));

	if flg_ann==1:
		ax.plot(tv[iT],y[iT],'-',color=[0.6,0.9,0],lw=0.25,label='Tot. merch. vol. (FCS)')
	ax.plot(tv[iT],y_ma,'-',color=[0.3,0.75,0],lw=1,label='Tot. merch vol. mov. ave. (FCS)')

	if flg_dead==1:
		ax.plot(tv[iT],y_Dead_ma,'--',color=[0,0.7,0],lw=1,label='Dead merch. vol. (FCS)')

	if flg_hbs==1:
		if pNam=='BCFCS_CWH':
			ax.plot(dH_ByBGCZ['Time']['CWH'],dH_ByBGCZ['Mm3/yr']['CWH'],color=[0,1,1],lw=1.25,label='Tot. merch. vol. (HBS)')
		elif pNam=='BCFCS_SBS':
			ax.plot(dH_ByBGCZ['Time']['SBS'],dH_ByBGCZ['Mm3/yr']['SBS'],color=[0,1,1],lw=1.25,label='Tot. merch. vol. (HBS)')
		else:
			ax.plot(dH['Year'],dH['Province']['Volume']['FLNR'],'k-',color=meta['Graphics']['gp']['cl1'],lw=1.25,label='Tot. merch. vol. (FLNR)')
			ax.plot(dH['Year'],dH['Province']['Volume']['HBS'],color=[0,1,1],lw=1.25,label='Tot. merch. vol. (HBS)')

	#ax.plot(tv[iT],y,'-',color=meta['Graphics']['gp']['cl2'],lw=0.5,label='Total')
	ax.set(position=[0.085,0.125,0.88,0.84],xlim=[tv[iT[0]],2100],xticks=np.arange(1800,2120,20),ylabel='Volume removed (Million m$^3$ yr$^-$$^1$)',xlabel='Time, years')
	ylim=ax.get_ylim()
	ax.set(ylim=[0,ylim[1]+0.03*ylim[1]])
	ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	#ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_Harvest_Volume_TS','png',900)

	return

#%% Volume per hectare of harvest
def HarvestVolumePerHectareTimeSeries(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS,dH_ByBGCZ):
	# Observations are relatively stable at 350 m3/ha over 1990-2018 (https://cfs.nrcan.gc.ca/statsprofile/forest/bc)

	dH=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestTimeSeries.pkl')

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));

	# Modelled
	ivlT=1;
	A=cbu.SummarizeAreaAffected(meta,pNam,mos,tv,iScn,iPS,iSS,iYS,ivlT)
	A_Harvest=A['Management']['Harvest']['Data']
	V_Harvest=mos[pNam]['Scenarios'][iScn]['Sum']['V_ToMill_MerchTotal']['Ensemble Mean'][:,iPS,iSS,iYS]
	rat=np.maximum(0,V_Harvest/A_Harvest)
	rat_ma=gu.movingave(rat,5,'historical')
	ax.plot(tv[iT],rat_ma[iT],'-go',lw=0.5,mew=0.5,mfc=[1,1,1],mec=[0.5,0.85,0],color=[0.5,0.85,0],ms=2,label='BC-FCS')

	if pNam=='BCFCS_CWH':
		ax.plot(dH_ByBGCZ['Time']['CWH'],dH_ByBGCZ['m3/ha']['CWH']*np.ones(dH_ByBGCZ['Time']['CWH'].size),color=[0,1,1],lw=1.25,label='Tot. merch. vol. (HBS)')
	else:
		y1=(dH['Province']['Volume']['FLNR']*1e6)/(dH['Province']['Area']['CC']*1e3)
		y2=(dH['Province']['Volume']['HBS']*1e6)/(dH['Province']['Area']['CC']*1e3)
		plt.plot(dH['Year'],y1,'-bs',color=[0.27,0.49,0.79],mfc=[0.27,0.49,0.79],mec=[0.27,0.49,0.79],ms=2,lw=0.5,mew=0.5,label='Harvest observations DB')
		plt.plot(dH['Year'],y2,'-b^',color=[0,1,1],mfc=[0.27,0.49,0.79],mec=[0.27,0.49,0.79],ms=2,lw=0.5,mew=0.5,label='Harvest observations DB (HBS)')

	ax.set(xticks=np.arange(tv[iT[0]],2250,5),xlabel='Time, years',yticks=np.arange(0,2000,100),ylabel='Harvest volume (m$^3$ ha$^-$$^1$)',
		xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[0,1000])
	ax.legend(loc='lower left',frameon=False,facecolor='w',edgecolor='w');
	ax.grid(False)
	plt.tight_layout()
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_Harvest_VolumePerHectare_TS','png',900)
	return

#%% Harvest rates (Actual vs. Forest Retention Scenario)
def PlotHarvestMultipleScenarios(meta,pNam,mos,tv,iT,iPS,iSS,iYS,iScn1,iScn2,flg_hbs):

	# Import FAIB public summary
	H_FLNR=gu.ReadExcel(r'C:\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));
	#ax.add_patch(Rectangle([2005,0],10,200,fc=[0.8,0.8,0.8],ec="none"))
	#ax.add_patch(Rectangle([2015,0],100,200,fc=[0.9,0.9,0.9],ec="none"))

	if flg_hbs==1:
		ax.plot(H_FLNR['Year'],H_FLNR['Total_harvest_millions_m3'],'k-',lw=1.25,label='Tot. merch. vol. (FLNR)')
		#ax.plot(H_HBS['Year'],H_HBS['V All m3']/1e6,'-cd')

	y=gu.movingave(mos[pNam]['Scenarios'][iScn1]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],10,'historical')/meta[pNam]['Project']['Multi']
	ax.plot(tv[iT],y,'-',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')

	y=gu.movingave(mos[pNam]['Scenarios'][iScn2]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][iT,iPS,iSS,iYS],10,'historical')/meta[pNam]['Project']['Multi']
	ax.plot(tv[iT],y,'--',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')

	ax.set(position=[0.09,0.12,0.88,0.84],xlim=[1860,2100],xticks=np.arange(1800,2120,20),ylim=[0,140], \
		   ylabel='Harvest volume (million m$^3$ yr$^-$$^1$)',xlabel='Time, years')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.legend(loc='upper left',frameon=False,facecolor='w',edgecolor='w');
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_HarvestVolumeMerchCompareScenarios_Scn' + str(iScn1+1) + '_vs_' + str(iScn2+1),'png',900)
	return

#%% Plot net sector GHG balance
def GHGBalance(meta,mos,pNam,tv,iScn,iT,iPS,iSS,iYS):

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,11)); wd=1;
	ax.add_patch(Rectangle([1920,-1000],103,2000,fc=[0.94,0.94,0.94],ec="none"))
	ax.plot(tv[iT],np.zeros(iT.size),'k-',color=meta['Graphics']['gp']['cla'])

	y_hwp=mos[pNam]['Scenarios'][iScn]['Sum']['E_Domestic_ForestSector_HWP']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
	ax.bar(tv[iT],y_hwp,wd,label='Product emissions',facecolor=[0.75,0.85,1])

	y_op=mos[pNam]['Scenarios'][iScn]['Sum']['E_ForestOperations']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
	ax.bar(tv[iT],y_op,wd,bottom=y_hwp,label='Fossil fuel emissions',facecolor=[0,0,0.7])

	y_ob=mos[pNam]['Scenarios'][iScn]['Sum']['E_OpenBurning_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
	ax.bar(tv[iT],y_ob,wd,bottom=y_hwp+y_op,label='Open burning emissions',facecolor=[1,0.1,0.1])

	y_wf=mos[pNam]['Scenarios'][iScn]['Sum']['E_Wildfire_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
	ax.bar(tv[iT],y_wf,wd,bottom=y_hwp+y_op+y_ob,label='Wildfire emissions',facecolor=[1,0.7,0.7])

	# Negative harvest removals
	#ax.bar(tv[iT],-1*(mos[pNam]['Scenarios'][iScn]['Sum']['C_ToMill']['Ensemble Mean'][iT,iPS,iSS,iYS])*3.667/meta[pNam]['Project']['Multi'],wd,facecolor=[0.65,.85,0.2],label='Removals (ecosystem to HWP)')
	ax.bar(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_Substitution_Total']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],wd,facecolor=[0.65,.85,0.2],label='Maximum substitution effects')

	ax.plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'k--',color=[0.5,0.5,0.5],mfc='w',mew=0.5,ms=2,label='Net atmospheric balance')
	ax.plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_NSB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'k-',color=[0,0,0],mfc='w',mew=0.75,ms=2,label='Net sector balance')
	ax.plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'k-.',color=[0.4,0.6,0.15],mfc='w',ms=2,mew=0.5,label='Net ecosystem exchange')

	if tv[iT[0]]>=1800:
		ax.set(xticks=np.arange(tv[iT[0]],2200,20),yticks=np.arange(-500,500,25),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[-350,250],xlim=[tv[iT[0]],tv[iT[-1]]])
	else:
		ax.set(position=[0.1,0.11,0.89,0.83],xticks=np.arange(tv[iT[0]],2200,50), \
			yticks=np.arange(-200,300,25),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[-180,220],xlim=[tv[iT[0]],tv[iT[-1]]])
	ax.legend(loc='lower left',fontsize=meta['Graphics']['gp']['fs_s'],frameon=False,facecolor='w',edgecolor='w');
	#plt.grid()
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.text(1865,150,'Emissions',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7])
	ax.text(1865,-150,'Removals',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7])
	ax.text(1970,100,'Modern\nperiod',fontsize=9,style='normal',weight='bold',color=[0.6,0.6,0.6],ha='center')
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_GHGBalance','png',900)
	return

#%% Plot net sector GHG balance (simple)
def GHGBalanceSimple(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS):
	v1='E_NAB'
	v2='E_NSB'
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(13,7.5));
	ax.plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'],lw=0.5)
	ylo1=(mos[pNam]['Scenarios'][iScn]['Sum'][v1]['Ensemble P250'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum'][v2]['Ensemble P250'][iT,iPS,iSS,iYS])/2/meta[pNam]['Project']['Multi']
	yhi1=(mos[pNam]['Scenarios'][iScn]['Sum'][v1]['Ensemble P750'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum'][v2]['Ensemble P750'][iT,iPS,iSS,iYS])/2/meta[pNam]['Project']['Multi']
	ylo2=(mos[pNam]['Scenarios'][iScn]['Sum'][v1]['Ensemble P025'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum'][v2]['Ensemble P025'][iT,iPS,iSS,iYS])/2/meta[pNam]['Project']['Multi']
	yhi2=(mos[pNam]['Scenarios'][iScn]['Sum'][v1]['Ensemble P975'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum'][v2]['Ensemble P975'][iT,iPS,iSS,iYS])/2/meta[pNam]['Project']['Multi']
	mu=(mos[pNam]['Scenarios'][iScn]['Sum'][v1]['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum'][v2]['Ensemble Mean'][iT,iPS,iSS,iYS])/2/meta[pNam]['Project']['Multi']

	cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
	ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
	ax.fill_between(tv[iT],ylo1,yhi1,color=meta['Graphics']['gp']['cl1'],alpha=0.3,lw=0,label='50% C.I.')
	ax.fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl1'],alpha=0.1,lw=0,label='95% C.I.')

	ax.plot(tv[iT],mu,'k-',color=meta['Graphics']['gp']['cl1'],mew=0.5,ms=2,label='Best estimate')
	leg=ax.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
	ax.text(tv[iT[0]]+7,165,'Source',fontsize=11,style='italic',weight='bold',color=[0.75,0.75,0.75],va='center')
	ax.text(tv[iT[0]]+7,-165,'Sink',fontsize=11,style='italic',weight='bold',color=[0.75,0.75,0.75],va='center')
	if tv[iT[0]]>=1800:
		ax.set(position=[0.1,0.11,0.89,0.83],xticks=np.arange(tv[iT[0]],2200,20), \
			   yticks=np.arange(-500,500,50),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]]) # ,ylim=[-300,200]
	else:
		ax.set(position=[0.1,0.11,0.89,0.83],xticks=np.arange(tv[iT[0]],2200,50), \
			   yticks=np.arange(-500,500,50),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]]) # ,ylim=[-300,200]
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_GHGBalance_AverageSubs','png',900);
	return

#%% Plot scenario comparisons
def ScenarioComp_GHGBalanceAnnualPlusCumu(meta,pNam,mos,cNam,iT,iPS,iSS,iYS):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	v1=['E_NAB','E_NSB']
	v2=['E_NAB_Cumulative_from_tref','E_NSB_Cumulative_from_tref']
	#v2='E_NAB_Cumulative'
	for iV in range(len(v1)):
		plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(18.8,10))
		ax[0,0].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ylo1=mos[pNam]['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi1=mos[pNam]['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ylo2=mos[pNam]['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		if np.isnan(ymn)==True:
			ymn=-1; ymx=1
		ax[0,0].fill_between(tv[iT],ylo1,yhi1,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0)
		ax[0,0].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0)
		ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'-',color=meta['Graphics']['gp']['cl1'],label='Baseline')
		ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'--',color=meta['Graphics']['gp']['cl2'],label='Actual')

		leg=ax[0,0].legend(loc='upper left',bbox_to_anchor=(0.1,0.45,0.5,0.5),frameon=False,facecolor=None,edgecolor='w');
		for text in leg.get_texts():
			plt.setp(text,color=meta['Graphics']['gp']['cla']);
		#ax[0,0].text(tv[iT[0]]+7,80,'Source',fontsize=7,style='italic',weight='bold',color=clt,va='center')
		#ax[0,0].text(tv[iT[0]]+7,-80,'Sink',fontsize=7,style='italic',weight='bold',color=clt,va='center')
		ax[0,0].set(xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50), \
		  ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],ylabel='GHG emissions\n(MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
		ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

		ax[0,1].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ylo1=mos[pNam]['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi1=mos[pNam]['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ylo2=mos[pNam]['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		if np.isnan(ymx)==True:
			ymx=1; ymn=-1
		ax[0,1].fill_between(tv[iT],ylo1,yhi1,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
		ax[0,1].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0)
		ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'-',color=meta['Graphics']['gp']['cl1'],label='Baseline')
		ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'--',color=meta['Graphics']['gp']['cl2'],label='Actual')
		for text in leg.get_texts():
			plt.setp(text,color=meta['Graphics']['gp']['cla'])
		ax[0,1].set(xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50), \
		  ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],ylabel='Cumulative GHG emissions\n(MtCO$_2$e)',xlabel='Time, years')
		ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

		ylo=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ymu=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo,ymu,yhi))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		if np.isnan(ymx)==True:
			ymx=1; ymn=-1
		ax[1,0].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ax[1,0].fill_between(tv[iT],ylo,yhi,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0,label='95% C.I.')
		ylo2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v1[iV]]['Ensemble P250'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v1[iV]]['Ensemble P750'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ax[1,0].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0,label='50% C.I.')
		ax[1,0].plot(tv[iT],ymu,'-',color=meta['Graphics']['gp']['cl3'],label='Mean')
		ax[1,0].set(xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50),
		  ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)], \
		  ylabel='$\Delta$ GHG emissions\n(MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
		ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
		#ax[1,0].text(tv[iT[0]]+7,0.6*np.abs(np.minimum(ymn,ymx)),'Net emission',fontsize=7,style='italic',weight='bold',color=clt,va='center')
		#ax[1,0].text(tv[iT[0]]+7,-0.6*np.abs(np.minimum(ymn,ymx)),'Net removal',fontsize=7,style='italic',weight='bold',color=clt,va='center')
		leg=ax[1,0].legend(loc='best',frameon=False,facecolor=None,edgecolor='w')
		for text in leg.get_texts():
			plt.setp(text,color=meta['Graphics']['gp']['cla']);

		ax[1,1].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ymu=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ylab='Cumulative $\Delta$GHG emissions\n(MtCO$_2$e)'
		#ylab='Cumulative $\Delta$GHG emissions\n(GtCO$_2$e)'
		ymu=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ylo=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo,ymu,yhi))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		if np.isnan(ymx)==True:
			ymx=1; ymn=-1
		ylo2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble P250'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v2[iV]]['Ensemble P750'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ax[1,1].fill_between(tv[iT],ylo,yhi,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
		ax[1,1].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0,label='50% C.I.')
		ax[1,1].plot(tv[iT],ymu,'-',color=meta['Graphics']['gp']['cl3'],label='Actual minus baseline')
		ax[1,1].set(xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50),ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],ylabel=ylab,xlabel='Time, years')
		ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

		gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
		plt.tight_layout()
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_ScnComp_' + cNam + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_ScenarioComp_GHGBalanceAnnualPlusCumu_' + v1[iV],'png',900)
	return

#%% GHG benefit (with Subs)
def ScenarioComp_GHGBalanceAnnual(meta,pNam,mos,tv,cNam,iT,iPS,iSS,iYS):

	vL=['E_NAB','E_NSB']
	lab=['WithSub','WithoutSub']

	iB=mos[pNam]['Delta'][cNam]['iB']
	iP=mos[pNam]['Delta'][cNam]['iP']
	for iv in range(len(vL)):
		v=vL[iv]
		plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.8,6));
		ax[0].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ylo1=mos[pNam]['Scenarios'][iB]['Sum'][v]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi1=mos[pNam]['Scenarios'][iB]['Sum'][v]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ylo2=mos[pNam]['Scenarios'][iP]['Sum'][v]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Scenarios'][iP]['Sum'][v]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		ax[0].fill_between(tv[iT],ylo1,yhi1,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0)
		ax[0].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0)
		ax[0].plot(tv[iT],mos[pNam]['Scenarios'][iB]['Sum'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'-',color=meta['Graphics']['gp']['cl1'],label='Baseline')
		ax[0].plot(tv[iT],mos[pNam]['Scenarios'][iP]['Sum'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'--',color=meta['Graphics']['gp']['cl2'],label='Actual')
		leg=ax[0].legend(loc='upper left',bbox_to_anchor=(0.08,0.48,0.5,0.5),frameon=False,facecolor=None,edgecolor='w')
		#for text in leg.get_texts():
		#	plt.setp(text, color=cla)
		ax[0].text(tv[iT[0]]+7,0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Source',fontsize=7,style='italic',weight='bold',color=meta['Graphics']['gp']['clt'],va='center')
		ax[0].text(tv[iT[0]]+7,-0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Sink',fontsize=7,style='italic',weight='bold',color=meta['Graphics']['gp']['clt'],va='center')
		ax[0].set(xticks=np.arange(tv[iT[0]],2200,50),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]],tv[iT[-1]]])
		#,ylim=[-160,220]yticks=np.arange(-200,300,40),
		ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
		
		ylo=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ymu=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		cs=np.column_stack((ylo,ymu,yhi))
		ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
		if np.isnan(ymx)==True:
			ymn=-1; ymx=1
		ax[1].plot(tv[iT],np.zeros(iT.size),color=meta['Graphics']['gp']['cla'])
		ax[1].fill_between(tv[iT],ylo,yhi,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0,label='95% C.I.')
		ylo2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble P250'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		yhi2=mos[pNam]['Delta'][cNam]['ByStrata']['Sum'][v]['Ensemble P750'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']
		ax[1].fill_between(tv[iT],ylo2,yhi2,color=meta['Graphics']['gp']['cl3'],alpha=meta['Graphics']['gp']['Alpha2'],lw=0,label='50% C.I.')
		ax[1].plot(tv[iT],ymu,'-',color=meta['Graphics']['gp']['cl3'],label='Mean')
		ax[1].set(xticks=np.arange(tv[iT[0]],2200,50),ylabel='$\Delta$GHG (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],xlim=[tv[iT[0]],tv[iT[-1]]])
		# ,yticks=np.arange(-200,300,20)
		ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
		ax[1].text(tv[iT[0]]+7,0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Net emission',fontsize=7,style='italic',weight='bold',color=meta['Graphics']['gp']['clt'],va='center');
		ax[1].text(tv[iT[0]]+7,-0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Net removal',fontsize=7,style='italic',weight='bold',color=meta['Graphics']['gp']['clt'],va='center');
		leg=ax[1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
		plt.tight_layout()
		gu.axletters(ax,plt,0.03,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight']);
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_ScnComp_' + cNam + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_ScenarioComp_GHGBalanceAnnual_' + lab[iv],'png',900);
	return

#%%
def ComparisonWithPIR(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS):
	dPIR=gu.ReadExcel(r'C:\Data\PIR\PIR Reformatted.xlsx',sheet_name='Sheet1',skiprows=0)
	iT=np.where( (tv>=1900) & (tv<=2023) )[0]

	plt.close('all')
	fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,8.5)); cl=np.array([[0.27,0.47,0.79],[0,1,1]]);
	ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',color=meta['Graphics']['gp']['cla'])

	ax[0,0].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'ob-',color=cl[0,:],mec=cl[0,:],mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='FCS')
	ax[0,0].plot(dPIR['Year'],dPIR['Forest Growth Minus Decay']/1e3,'rs-',mfc=[1,1,1],ms=2,mew=0.5,lw=meta['Graphics']['gp']['lw1'],label='PIR')
	ax[0,0].set(xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-300,300,50),
		   ylabel='Growth - Decay (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[-250,100])
	ax[0,0].legend(loc='best',frameon=False,facecolor='w',edgecolor='w')
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[0,1].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_Wildfire_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'ob-',color=cl[0,:],mec=cl[0,:],mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='Net ecosystem exchange')
	ax[0,1].plot(dPIR['Year'],dPIR['Wildfires']/1e3,'rs-',mfc=[1,1,1],mew=0.5,ms=2,lw=meta['Graphics']['gp']['lw1'],label='PIR')
	ax[0,1].set(xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(0,300,25),
		   ylabel='Wildfire (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[0,220])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	#ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_OpenBurning_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'ob-',mfc='w',mew=0.5,ms=2,label='Net ecosystem exchange')
	#ax[1,0].plot(dPIR['Year'],dPIR['Slash Pile Burning']/1e3,'rs-',mfc=[1,1,1],label='PIR')
	ax[1,0].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_Domestic_ForestSector_HWP']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'ob-',color=cl[0,:],mec=cl[0,:],mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='HWP decay (FCS)')
	ax[1,0].plot(tv[iT],(mos[pNam]['Scenarios'][iScn]['Sum']['E_Domestic_ForestSector_HWP']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Sum']['E_Bioenergy']['Ensemble Mean'][iT,iPS,iSS,iYS])/meta[pNam]['Project']['Multi'],'oc-',mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='HWP decay + bioenergy (FCS)')
	#ax[1,0].plot(tv[iT],gu.movingave(mos[pNam]['Scenarios'][iScn]['Sum']['C_ToMill']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi']*3.667,10,'center'),'og-',mfc='w',mew=0.5,ms=2,label='Harvest removals (FCAST)')
	ax[1,0].plot(dPIR['Year'],dPIR['Decomposition of Harvested Wood Products']/1e3,'rs-',mfc=[1,1,1],mew=0.5,lw=meta['Graphics']['gp']['lw1'],ms=2,label='HWP decay (PIR)')
	ax[1,0].set(xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-120,300,20),
		   ylabel='HWP emissions (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[0,100])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1,0].legend(loc='center left',frameon=False,facecolor='w',edgecolor='w')

	#ax[1,0].plot(Harvest['Year'],Harvest['Total_harvest_millions_m3']*0.46*0.5*3.667,'k-',lw=0.75,label='FLNRORD estimate')
	#ax.plot(H_HBS['Year'],H_HBS['V All m3']/1e6,'-cd')

	ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',color=meta['Graphics']['gp']['cla'])
	ax[1,1].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'ob-',color=cl[0,:],mec=cl[0,:],mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='With subs. (FCS)')
	ax[1,1].plot(tv[iT],mos[pNam]['Scenarios'][iScn]['Sum']['E_NSB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'],'sb-',color='c',mec='c',mfc='w',lw=meta['Graphics']['gp']['lw1'],mew=0.5,ms=2,label='W/O subs. (FCS)')
	ax[1,1].plot(dPIR['Year'],dPIR['Forest Management']/1e3,'rs-',ms=2,mfc='w',mew=0.5,lw=meta['Graphics']['gp']['lw1'],label='Forest Management (PIR)')
	ax[1,1].set(xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-350,300,50),
		   ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[-250,250])
	ax[1,1].legend(loc='upper center',frameon=False,facecolor='w',edgecolor='w')
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	#plt.grid()
	gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_ComparisonWithPIR','png',900)
	return

#%%
def Area_DisturbedAndManaged(meta,pNam,mos,ivlT,iScn,iT,iPS,iSS,iYS,iOS,**kwargs):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	A=cbu.SummarizeAreaAffected(meta,pNam,mos,tv,iScn,iPS,iSS,iYS,iOS,ivlT)
	yr_start=tv[iT[0]]
	if tv[iT[0]]>=1800:
		xtivl=20
	else:
		xtivl=50

	if iT.size<50:
		xtivl=5

	if 'multi' in kwargs.keys():
		multi=kwargs['multi']
		ylab='Area affected (ha x ' + str(int(multi)) + ')'
	else:
		multi=1e6
		ylab='Area affected (Mha)'

	if 'ncols' in kwargs.keys():
		ncols=kwargs['ncols']
	else:
		ncols=1

	if 'LegendLoc' in kwargs.keys():
		LegendLoc=kwargs['LegendLoc']
	else:
		LegendLoc='upper left'


	plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(22,10));

	pl_d=[None]*len(A['Natural'])
	nams_d=[None]*len(A['Natural'])
	cnt=0
	bottom=0
	for k in A['Natural'].keys():
		pl_d[cnt]=ax[0].bar(tv[iT],A['Natural'][k]['Data'][iT]/multi,ivlT,color=A['Natural'][k]['Color'],bottom=bottom)
		nams_d[cnt]=k
		bottom=bottom+A['Natural'][k]['Data'][iT]/multi
		cnt=cnt+1
	ax[0].legend(pl_d,nams_d,loc=LegendLoc,labelspacing=0.12,facecolor=[1,1,1],
			  frameon=False,fontsize=6,ncols=ncols)
	ax[0].set(yscale='linear',xticks=np.arange(np.min(tv),np.max(tv)+1,xtivl),ylabel=ylab,
		   xlim=[yr_start,np.max(tv[iT])],ylim=[0,1.05*np.max(bottom)])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	#ax[0].text((2020+1920)/2,1050,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])

	#ax[1].fill_between([1920,2020],[0,0],[1100,1100],color=[0.9,0.9,0.9],lw=0)
	pl_m=[None]*len(A['Management'])
	nams_m=[None]*len(A['Management']);
	cnt=0
	bottom=0
	for k in A['Management'].keys():
		pl_m[cnt]=ax[1].bar(tv[iT],A['Management'][k]['Data'][iT]/multi,ivlT,color=A['Management'][k]['Color'],bottom=bottom)
		nams_m[cnt]=k
		bottom=bottom+A['Management'][k]['Data'][iT]/multi
		cnt=cnt+1
	#ax[1].plot([2020,2020],[0,5500],'k--')
	ax[1].legend(pl_m,nams_m,loc=LegendLoc,labelspacing=0.12,facecolor=[1,1,1],frameon=False,ncols=ncols,fontsize=6)
	ax[1].set(yscale='linear',xticks=np.arange(np.min(tv),np.max(tv)+1,xtivl),xlabel='Time, years',ylabel=ylab,
		   xlim=[yr_start,np.max(tv[iT])],ylim=[0,1.05*np.max(bottom)])
	ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.014,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
	nam_os=meta[pNam]['Project']['Strata']['Other']['Unique CD'][iOS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\AreaAffected_ByEventType_' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_' + nam_os,'png',900)
	return

#%% Mortality summary
def CarbonFluxMortality(meta,pNam,mos,tv,iScn,iT,ivlT,iPS,iSS,iYS):

	y=[None]*7; c=-1
	c=c+1; y[c]={}; y[c]['Name']='Competition'; y[c]['Color']=[1,0.9,0.55]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Reg']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Wind'; y[c]['Color']=[0.27,0.49,0.74]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Wind']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Disease'; y[c]['Color']=[0.9,0.85,1]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Disease Root']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Disease Foliage']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Disease Stem']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Wildfire'; y[c]['Color']=[0.75,0,0]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Wildfire']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Beetles'; y[c]['Color']=[0.5,0.8,0]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Mountain Pine Beetle']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Douglas-fir Beetle']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Balsam Beetle']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Spruce Beetle']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Defoliators'; y[c]['Color']=[0.75,1,0.1]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Western Spruce Budworm']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Hemlock Looper']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#c=c+1; y[c]={}; y[c]['Name']='Harvest'; y[c]['Color']=[0.75,0.9,1]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Harvest']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Harvest Salvage']['Ensemble Mean'][iT,iPS,iSS,iYS]
	c=c+1; y[c]={}; y[c]['Name']='Harvest'; y[c]['Color']=[0.75,0.9,1]; y[c]['Data']=mos[pNam]['Scenarios'][iScn]['Mean']['C_M_Harvest']['Ensemble Mean'][iT,iPS,iSS,iYS]

	# Convert to x-year intervals
	y[0]['tv']=gu.BlockMean(tv[iT],ivlT)
	for i in range(len(y)):
		y[i]['Data']=gu.BlockMean(y[i]['Data'],ivlT)

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(16.5,8)); yr_start=tv[iT[0]];
	pl_d=[None]*len(y); nams_d=[None]*len(y);
	for i in range(len(y)):
		bottom=0;
		if i!=0:
			for j in range(i):
				bottom=bottom+y[j]['Data']
		pl_d[i]=ax.bar(y[0]['tv'],y[i]['Data'],ivlT,color=y[i]['Color'],bottom=bottom)
		nams_d[i]=y[i]['Name']
	ax.legend(pl_d,nams_d,loc='upper left',labelspacing=0.12,facecolor=[1,1,1],frameon=False) # ,bbox_to_anchor=(0.03,0.96)
	if yr_start<1800:
		xivl=50
	else:
		xivl=20
	ax.set(ylabel='Mortality (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Time, years',xticks=np.arange(np.min(y[0]['tv']),np.max(y[0]['tv'])+1,xivl),
		xlim=[yr_start,np.max(y[0]['tv'])],ylim=[0,1.05*np.max(bottom)]); # ,ylim=[0,3]yticks=np.arange(0,5,0.2),
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_CarbonFluxMortality','png',900)

#%% Plot time series
def AllVariableTimeSeries(meta,pNam,mos,tv,iScn,iT,iPS,iSS,iYS):
	
	pth=meta['Graphics']['Print Figure Path'] + '\\TS'
	if os.path.exists(pth)==False:
		os.mkdir(pth)
		
	for k in mos[pNam]['Scenarios'][iScn]['Sum']:
		plt.close('all')
		fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,12))
		mxn=np.zeros(2)
		for iScn in range(meta[pNam]['Project']['N Scenario']):
			#mu=mos[iScn]['v1'][k]['Mean']
			#mu=mos[iScn]['v1']['Mean'][k]['Ensemble Mean'];
			mu=mos[pNam]['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][:,iPS,iSS,iYS]
			mxn[0]=np.minimum(mxn[0],np.min(mu[iT])); 
			mxn[1]=np.maximum(mxn[1],np.max(mu[iT]));
			ax.plot(tv[iT],mu[iT],'-',lw=0.75)
		mu=np.mean(mxn)
		mxn[0]=mxn[0]-0.1*mu
		mxn[1]=mxn[1]+0.1*mu
		ax.set(xlim=[np.min(tv[iT]),np.max(tv[iT])],ylabel=k,xlabel='Time, years',aspect='auto'); # ylim=mxn,
		fig.patch.set_facecolor('w')
		plt.tight_layout()
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		if meta['Graphics']['Print Figures']=='On':
			gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\TS\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_ts_' + k + '_mean','png',150)
	return

#%% Plot other fluxes
def CarbonFluxTimeSeries(meta,pNam,mos,iT,iScn,iPS,iSS,iYS):
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	xt=np.arange(0,3000,50)
	plt.close('all'); fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(18,12))
	#ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=(mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Gross']['Ensemble P025'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_LF']['Ensemble P025'][iT,iPS,iSS,iYS])*3.667
	hi=(mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Gross']['Ensemble P975'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_LF']['Ensemble P975'][iT,iPS,iSS,iYS])*3.667
	mu=(mos[pNam]['Scenarios'][iScn]['Mean']['C_G_Gross']['Ensemble Mean'][iT,iPS,iSS,iYS]+mos[pNam]['Scenarios'][iScn]['Mean']['C_LF']['Ensemble Mean'][iT,iPS,iSS,iYS])*3.667
	ax[0,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,0].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[0,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[0,0].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[0,0].set(xticks=xt,xlabel='Time, years',ylabel='NPP\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	#ax[0,0].legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	#ax[0,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iScn]['Mean']['C_RH']['Ensemble P025'][iT,iPS,iSS,iYS]*3.667
	hi=mos[pNam]['Scenarios'][iScn]['Mean']['C_RH']['Ensemble P975'][iT,iPS,iSS,iYS]*3.667
	mu=mos[pNam]['Scenarios'][iScn]['Mean']['C_RH']['Ensemble Mean'][iT,iPS,iSS,iYS]*3.667
	ax[0,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,1].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[0,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[0,1].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[0,1].set(xticks=xt,xlabel='Time, years',ylabel='RH\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	
	ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,0].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_NEE']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[0,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[0,0].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[1,0].set(xticks=xt,xlabel='Time, years',ylabel='NEE\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	#ax[1,0].legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iScn]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iScn]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iScn]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,1].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_Wildfire_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[0,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[0,1].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[1,1].set(xticks=xt,xlabel='Time, years',ylabel='Wildfire emissions\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[2,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iScn]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iScn]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iScn]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[2,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[2,0].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_OpenBurning_ForestSector_Domestic']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[1,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[1,0].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[2,0].set(xticks=xt,xlabel='Time, years',ylabel='Open burning\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])

	ax[2,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iScn]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[2,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl2'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[2,1].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl2'],label='Actual scenario')
	#lo=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble P025'][iT,iPS,iSS,iYS]
	#hi=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble P975'][iT,iPS,iSS,iYS]
	#mu=mos[pNam]['Scenarios'][iB]['Mean']['E_Domestic_ForestSector_HWP']['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax[1,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax[1,1].plot(tv[iT],mu,'--',color=meta['Graphics']['gp']['cl1'],label='Baseline scenario')
	ax[2,1].set(xticks=xt,xlabel='Time, years',ylabel='Product emissions\n(tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	gu.axletters(ax,plt,0.025,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout()
	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_CarbonFluxTimeSeries','png',900)

#%% Plot pools
def PlotCarbonPoolTS(meta,pNam,mos,tv,iT,iB,iP,iPS,iSS,iYS):

	ms=2; Alpha=0.16; lw=1; cl0=[0.27,0.44,0.79]; cl1=[0.4,0.8,0]; cl2=[0,0.9,0.9]
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,10))

	v='C_Biomass'
	#ax[0,0].plot(tv[it],np.zeros(it.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,0].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl1'],label='Actual scenario')
	lo=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,0].fill_between(tv[iT],lo,hi,color=cl0,alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,0].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario')
	ax[0,0].set(position=[0.07,0.57,0.42,0.42],xticks=np.arange(0,2220,25), \
	  xlabel='Time, years',ylabel='Biomass (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[0,0].legend(loc='lower right',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both')
	ymin,ymax=ax[0,0].get_ylim(); ax[0,0].set_ylim(0.0,ymax)

	v='C_DeadWood'
	#ax[0,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,1].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl1'],label='Actual scenario')
	lo=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[0,1].fill_between(tv[iT],lo,hi,color=cl0,alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[0,1].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario')
	ax[0,1].set(position=[0.57,0.57,0.42,0.42],yscale='linear',xticks=np.arange(0,2220,25), \
	  xlabel='Time, years',ylabel='Dead wood (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both')
	ymin,ymax=ax[0,1].get_ylim(); ax[0,1].set_ylim(0.0,ymax)

	v='C_Litter'
	#ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,0].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl1'],label='Actual scenario')
	lo=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,0].fill_between(tv[iT],lo,hi,color=cl0,alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,0].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario')
	ax[1,0].set(position=[0.07,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
	  xlabel='Time, years',ylabel='Litter (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both')
	ymin,ymax=ax[1,0].get_ylim(); ax[1,0].set_ylim(0.0,ymax)

	v='C_Soil'
	#ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
	lo=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,1].plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl1'],label='Actual scenario')
	lo=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax[1,1].fill_between(tv[iT],lo,hi,color=cl0,alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax[1,1].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario')
	ax[1,1].set(position=[0.57,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
	  xlabel='Time, years',ylabel='Soil (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both')
	ymin,ymax=ax[1,1].get_ylim(); ax[1,1].set_ylim(0.0,ymax)

	gu.axletters(ax,plt,0.03,0.91)

	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_' + nam_ps + '_' + nam_ss + '_CarbonPools_Scns' + str(iB) + 'and' + str(iP),'png',900)

#%% Plot age
def PlotAgeTS(meta,pNam,mos,tv,iT,iB,iP,iPS,iSS,iYS):
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,10)); cl=np.array([[0.27,0.49,0.74],[0.6,0.95,0]])
	ax.plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)	
	for iScn in range(meta[pNam]['Project']['N Scenario']):
		mu=mos[pNam]['Scenarios'][iScn]['Mean']['A']['Ensemble Mean'][iT,iPS,iSS,iYS]	
		ax.plot(tv[iT],mu,'-',color=cl[iScn,:],label='Scenario ' + str(iScn+1))
	ax.set(xticks=np.arange(0,2220,25),xlabel='Time, years',ylabel='NEE (tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	#ax[0].legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	return

#%% Plot net growth
def NetGrowthTimeSeries(meta,pNam,mos,tv,iT,iScn,iPS,iSS,iYS):

	#ms=2; Alpha=0.16; lw=1; cl0=[0.27,0.44,0.79]; cl1=[0.4,0.8,0]; cl2=[0,0.9,0.9]
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,10))
	#ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)

	v='C_G_Net'
	lo=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	ax.fill_between(tv[iT],lo,hi,color=meta['Graphics']['gp']['cl1'],alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax.plot(tv[iT],mu,'-',color=meta['Graphics']['gp']['cl1'],label='Net growth (w/o instects)')

	v='C_M_Reg'
	lo1=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi1=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu1=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]
	#ax.fill_between(tv[iT],lo1,hi1,color='r',alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax.plot(tv[iT],mu1,'-',color='r',label='Mortality'+0.25)

	v='C_M_Dist'
	lo2=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS,iYS]
	hi2=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS,iYS]
	mu2=mos[pNam]['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS,iYS]

	# Mortality from beetles
	#Mb=C_M_ByAgent[iScn]['Beetles'][iT]+C_M_ByAgent[iScn]['IBM'][iT]+C_M_ByAgent[iScn]['IBS'][iT]+C_M_ByAgent[iScn]['IBB'][iT]

	lo3=lo1+lo2
	hi3=hi1+hi2
	mu3=gu.movingave(mu1+mu2,10,'historical')
	#ax.fill_between(tv[iT],lo3,hi3,color='y',alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	#ax.plot(tv[iT],mu3,'-',color='y',label='Mortality'+0.25)

	mth='center'
	lo4=gu.movingave(lo,7,mth)
	hi4=gu.movingave(hi,7,mth)
	mu4=gu.movingave(mu,7,mth)
	ax.fill_between(tv[iT],lo4,hi4,color='y',alpha=meta['Graphics']['gp']['Alpha1'],lw=0)
	ax.plot(tv[iT],mu4,'-',color='y',label='Net growth (with insects)')

	#------------

	#d=gu.ImportMat(r'G:\My Drive\Data\psp_sl_ts_bc.mat','qt')
	#ax.plot(d['tv'].flatten(),np.nanmean(d['NEBP'],axis=1),'-bo')

	ax.set(position=[0.07,0.07,0.84,0.84],
	  xticks=np.arange(0,2220,25), \
	  xlabel='Time, years', \
	  ylabel='Net growth (tC ha$^-$$^1$yr$^-$$^1$)',
	  xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
	ax.legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

	#gu.axletters(ax,plt,0.03,0.91,FontColor=cla,LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
	nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
	nam_ys=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iYS]
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_NetGrowthTimeSeries','png',900)
	return

#%% Mortality spectrum (frequency distribution)
def MortalitySpectrum(meta,pNam,iEns,iScn,iPS,iSS,iYS):
	# Import simualtions
	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	
	vInc=['Wildfire','Mountain Pine Beetle','Balsam Beetle','Douglas-fir Beetle','Spruce Beetle','Western Spruce Budworm','Hemlock Looper',
		  'Disease Root','Disease Foliage','Disease Stem','Wind','Frost Snow Ice Hail','Flooding Lightning Slides','Drought',
		  'Harvest','FL-CL','FL-PA','FL-RC','FL-EM','FL-TR']
	
	M={}
	M['Reg']=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))
	for v in vInc:
		M[v]=np.zeros((tv.size,meta[pNam]['Project']['N Stand']))
		
	for iBat in range(meta[pNam]['Project']['N Batch']):
		indBat=cbu.IndexToBatch(meta[pNam],iBat)
		d1=cbu.LoadSingleOutputFile(meta,pNam,iScn,iEns,iBat)
		for v in vInc:
			id=meta['LUT']['Event'][v]
			M[v][:,indBat]=d1['C_M_DistByAgentPct'][id].astype('float32')*meta['Core']['Scale Factor C_M_DistByAgent']
		M['Reg'][1::,indBat]=np.nan_to_num(d1['C_M_Reg'][1::,:]/d1['C_Biomass'][0::-1,:]*100)

	# Frequency distribution
	bw=1; bin=np.arange(0,100+bw,bw)
	Frq={}
	for k in M.keys():
		y=M[k].copy().flatten()
		Frq[k]=np.zeros(bin.size)
		for iBin in range(bin.size):
			ind=np.where(np.abs(y-bin[iBin])<=bw/2)[0]
			Frq[k][iBin]=ind.size/y.size*100

	# Plot
	lab=list(M.keys()); lab[0]='Competition'
	dat={}
	cnt=0
	Bottom=np.zeros(bin.size)
	for k in M.keys():
		dat[k]={}
		dat[k]['Data']=Frq[k]
		dat[k]['Bottom']=Bottom.copy()
		Bottom=Bottom+Frq[k].copy()
		cnt=cnt+1

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9));
	br=[None]*len(dat)
	cnt=0
	for k in dat.keys():
		if k=='Reg':
			cl=[0.85,0.85,0.85]
		else:
			ind=np.where(meta['Param']['Raw']['Events']['Name']==k)[0][0]
			cl=[meta['Param']['Raw']['Events']['clr'][ind],meta['Param']['Raw']['Events']['clg'][ind],meta['Param']['Raw']['Events']['clb'][ind]]
		br[cnt]=plt.bar(bin,dat[k]['Data'],bottom=dat[k]['Bottom'],width=bw,color=cl,label=lab[cnt])
		cnt=cnt+1
	ax.legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w',ncol=2)
	ax.set(ylabel=r'Frequency, with stacked overlap (%)',xlabel='Stand mortality rate (%)',xlim=[-0.5,100.5],xticks=np.arange(0,110,10),yscale='log')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_MortalitySpectrum','png',900)

	return M

#%% Plot map of investigation sites
def PlotSitesOfInterest(meta):

	# Load basemap
	gdf_bm=gpd.read_file(r'C:\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

	gdf_sxy=gpd.read_file(meta['Paths']['Geospatial'] + '\\geos.geojson')

	plt.close('all')
	fig,ax=plt.subplots(figsize=gu.cm2inch(7.8,6.6))
	#mngr=plt.get_current_fig_manager()
	#mngr.window.setGeometry(700,20,620,600)
	gdf_bm.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
	#tsa_boundaries.plot(ax=ax,facecolor='none',edgecolor=[0,0,0],linewidth=0.25)
	gdf_sxy.plot(ax=ax,markersize=1,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
	ax.grid(color='k',linestyle='-',linewidth=0.25)

	#iP=2000
	#for iD in range(len(nddat[iP])):
	#	x,y=nddat[iP][iD]['Geometry'].exterior.xy
	#	plt.plot(x,y,'r-')

	ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
	#plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)
	plt.close('all')

	return

#%% Comparison of anthropogenic component with CFS
def CompareAnthropogenicComponentWithCFS(meta,pNam,mos,iT,iPS,iSS,iYS):

	tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)
	dfCFS=gu.ReadExcel(r'C:\Data\Emission Reduction Projections\Summary of Reporting Initiatives.xlsx')

	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,5.5)); wd=0.8; lw=0.5; cla=[0.8,0.8,0.8]
	ax.bar(1,np.mean(dfCFS['BC FLFL+HWP Anthropogenic (NIR22)']),label='National GHG Inventory (CFS)',facecolor=[0.5,0.7,1])
	#ax.bar(2,np.mean(dfCFS['NEW']),label='National GHG Inventory (CFS old)',facecolor=[0.2,0.4,1])
	ax.text(1,np.mean(dfCFS['BC FLFL+HWP Anthropogenic (NIR22)'])-48,'National\nGHG inventory\n(BC CFS)',fontsize=7,style='italic',weight='bold',color=[0.5,0.7,1],ha='center')
	#ax.text(2,np.mean(dfCFS['NEW'])-45,'National\nGHG inventory\n(BC CFS new)',fontsize=9,style='italic',weight='bold',color=[0.2,0.4,1],ha='center')
	iT=np.where( (tv>=1990) & (tv<=2020) )[0]
	y1=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NSB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.bar(2,y1,wd,label='BC FCS W/O\nsubstitution\neffects',facecolor=[1,0.7,0.5])
	y_lo=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NSB']['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	y_hi=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NSB']['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(2,y1,yerr=y_hi-y1,color=[0.8,0.5,0.2],ls='',lw=1.5,capsize=2)
	ax.text(2,y1+25,'BC FCS\nW/O\nsubstitution',fontsize=7,style='italic',weight='bold',color=[1,0.7,0.5],ha='center')
	y2=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NAB']['Ensemble Mean'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	y_lo=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NAB']['Ensemble P025'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	y_hi=np.mean(mos[pNam]['Delta']['TDAF']['ByStrata']['Sum']['E_NAB']['Ensemble P975'][iT,iPS,iSS,iYS]/meta[pNam]['Project']['Multi'])
	ax.bar(3,y2,wd,label='BC FCS with\nsubstitution\neffects',facecolor=[1,0.4,0.2])
	if meta[pNam]['Project']['N Ensemble']>1:
		ax.errorbar(3,y2,yerr=y2-y_lo,color=[0.4,0.2,0],ls='',lw=1.5,capsize=2)
	ax.text(3,y2+35,'BC FCS\nwith\nsubstitution',fontsize=7,style='italic',weight='bold',color=[1,0.4,0.2],ha='center')
	ax.plot([0,5],[0,0],'k-',color=cla,lw=0.25)
	ax.text(0.5,125,'Emissions',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7],va='center')
	ax.text(0.5,-125,'Removals',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7],va='center')
	ax.set(position=[0.16,.05,0.8,0.88],xticks=[0],yticks=np.arange(-175,300,25),ylabel='Flux (MtCO$_2$e yr$^-$$^1$)',xticklabels={''},xlabel='',ylim=[-140,160],xlim=[0.35,3.65])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	#ax.set_title('British Columbia, 1990-2020',fontsize=11,weight='bold',color=cla)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_DirectAnthroFlux_BC','png',900)
	return

#%% Plot Summary of existing initiatives
def SummarizeExistingInitiatives(meta,pNam):
	d=gu.ReadExcel(r'C:\Data\Emission Reduction Projections\Summary of Reporting Initiatives.xlsx')
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8));
	plt.plot(d['Year'],np.zeros(d['Year'].size),'k-',lw=0.5)
	plt.plot(d['Year'],d['BC Total Emissions FLFL+HWP All Gases (NIR2022)'],'-ko',lw=1,color=[0.7,0.7,0.7],mec=[0.7,0.7,0.7],mfc=[0.7,0.7,0.7],label='Total emissions NIR 22')
	#plt.plot(d['Year'],d['Forest Management (PIR)'],'cs')
	plt.plot(d['Year'],d['BC PTT BAU Scenario (CFS)'],'-rs',label='PTT BAU (CFS)',lw=1)
	plt.plot(d['Year'],d['BC PTT Reference Level Scenario (CFS)'],'-b^',label='PTT RL (CFS)',lw=1)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
	ax.set(position=[0.1,0.1,0.82,0.82],xlim=[d['Year'][0],d['Year'][-1]],ylabel='Flux (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Summary of Existing Initiatives','png',900)

	return

#%% Look at historical wildfire from aspatial model
def WildfireProbHistTimeSeries(meta,pNam,iScn):

	hw={}
	hw['tv']=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)
	hw['Po']=np.zeros((meta[pNam]['Project']['N Time'],meta[pNam]['Project']['N Stand']),dtype='int8')

	for iEns in range(meta[pNam]['Project']['N Ensemble']):

		if meta[pNam]['Project']['Frozen Ensembles Status']=='Off':
			wf_sim=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Inputs\\Ensembles\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
		else:
			wf_sim=gu.ipickle(meta[pNam]['Project']['Frozen Ensembles Path'] + '\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')

		if 'idx' in wf_sim:
			idx=wf_sim['idx']
			hw['Po'][idx[0],idx[1]]=hw['Po'][idx[0],idx[1]]+wf_sim['Occurrence']

	hw['Po']=hw['Po']/meta[pNam]['Project']['N Ensemble']

	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7.5));
	ax.plot(hw['tv'],100*np.mean(hw['Po'],axis=1))
	ax.set(position=[0.07,0.1,0.92,0.82],ylabel='Prob wildfire (% yr$^-$$^1$)',xlabel='Time, years',xlim=[hw['tv'][0],hw['tv'][-1]])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_WildfireProbHistTimeSeries','png',900)

	return hw

#%% Age class distribution
def AgeClassDist(meta,pNam,iScn,iPS,iSS,iYS):

	# Import modelled data
	acd=gu.ipickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\MOS_ByStrata_AgeClassDist_Scn' + str(iScn+1) + '.pkl')

	# Import VRI
	zA=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI 2023\PROJ_AGE_1.tif')
	zA=zA['Data'].flatten()[0::5]
	acd['Data VRI']=np.zeros(acd['binA'].size)
	for iA in range(acd['binA'].size):
		ind=np.where( np.abs(zA-acd['binA'][iA])<=acd['bwA']/2 )[0]
		acd['Data VRI'][iA]=ind.size
	acd['Data VRI'][0]=0

	# Import age class distributions fitted against field plots
	dFP=gu.ReadExcel(meta['Paths']['DB']['Field Plot Observations'] + '\\FP_AgeDistributionsByBGCZone.xlsx',sheet_name='Sheet1',skiprows=0)

	lw=1.25
	tL=[1850,2020]
	cl=['b','g','r','c']

	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8));

	# Add redictions
	cnt=0
	for t in tL:
		iT=np.where(acd['binT']==t)[0]
		if iT.size>0:
			plt.plot(acd['binA'],gu.movingave(acd['Data'][iT[0],:,iPS,iSS,iYS]/np.sum(acd['Data'][iT[0],:,iPS,iSS,iYS])*100,10,'Centre'),'g-',color=cl[cnt],lw=lw,label='FCS ' + str(t))	
			cnt=cnt+1

	# Add observations
	if pNam=='BCFC_CWH':
		plt.plot(dFP['Age'],dFP['CWH'],'b-',color=[0.24,0.47,0.79],lw=2,label='Field plots (2000-2020)')
	else:
		plt.plot(acd['binA'],gu.movingave(acd['Data VRI']/np.sum(acd['Data VRI'])*100,10,'Centre'),'k--',color=[0.24,0.47,0.79],lw=2,label='Vegetation resource inventory')

	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax.legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
	ax.set(position=[0.1,0.1,0.82,0.82],ylabel='Frequency (%)',xlabel='Stand age, years since major disturbance',xlim=[0,400],ylim=[0,3])
	if meta['Graphics']['Print Figures']=='On':
		nam_ps=meta[pNam]['Project']['Strata']['Project Type']['Unique CD'][iPS]
		nam_ss=meta[pNam]['Project']['Strata']['Spatial']['Unique CD'][iSS]
		nam_ys=meta[pNam]['Project']['Strata']['Year']['Unique CD'][iYS]
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_Scn' + str(iScn+1) + '_' + nam_ps + '_' + nam_ss + '_' + nam_ys + '_AgeClassDist','png',900)

	return

#%%
def Plot_NSR_Area(meta):
	# Import data
	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\HarvestArea_TS.pkl')

	# Import NOSE
	nose=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NOSE_Sparse\Inputs\AIL.pkl')
	A_regen_wo_nose=np.zeros(d['tv'].size)
	iT=np.where( (d['tv']>=nose['Year'][0]) & (d['tv']<=nose['Year'][-1]) )[0]
	A_regen_wo_nose[iT]=nose['Area All']

	plt.close('all')
	#plt.plot(d['tv'],,'og-')
	plt.plot(d['tv'],d['Area Harvested Max']-d['Area Planted RESULTS'],'-ko')

	NSR=d['Area Harvested Max']-d['Area Planted RESULTS']
	NSR_wo_nose=d['Area Harvested Max']-d['Area Planted RESULTS']+A_regen_wo_nose

	plt.close('all'); cl=np.array([[0.17,0.35,0.7],[0.3,0.6,0.8],[0.5,0.9,1],[0,0.75,0],[0.6,1,0]])
	fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8));
	ax[0].plot(d['tv'],np.zeros(d['tv'].size),'-k',lw=0.5)
	ax[0].fill_between(d['tv'],NSR/1e3,NSR_wo_nose/1e3,color=[0.9,0.9,0.9],alpha=1,lw=0,label='NOSE programs')
	ax[0].plot(d['tv'],d['Area Harvested Max']/1e3,'-k',lw=0.5,ms=2,label='Area harvested')
	ax[0].plot(d['tv'],NSR_wo_nose/1e3,'-g',lw=0.5,ms=1.25,label='Change in NSR area without NOSE')
	ax[0].plot(d['tv'],NSR/1e3,'-b',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],lw=0.5,ms=1.25,label='Change in NSR area with NOSE')
	ax[0].set(xticks=np.arange(1800,2120,10),yticks=np.arange(-150,350,50),ylabel='Change in area (ha x 1000 yr$^-$$^1$)',xlabel='Time, years',ylim=[-125,300],xlim=[1949.5,2021.5])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
	ax[0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

	ax[1].fill_between(d['tv'],np.cumsum(NSR/1e6),np.cumsum(NSR_wo_nose/1e6),color=[0.9,0.9,0.9],alpha=1,lw=0)
	ax[1].plot(d['tv'],np.cumsum(d['Area Harvested Max'])/1e6,'-k',lw=0.5,ms=2,label='Area harvested')
	ax[1].plot(d['tv'],np.cumsum(NSR_wo_nose)/1e6,'-g',lw=0.5,ms=2,label='NSR without NOSE')
	ax[1].plot(d['tv'],np.cumsum((d['Area Harvested Max']-d['Area Planted RESULTS'])/1e6),'-b',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],ms=2,label='NSR with NOSE')
	ax[1].set(xticks=np.arange(1800,2120,10),yticks=np.arange(0,10000,1),ylabel='Cumulative area (Mha)',xlabel='Time, years',ylim=[0,12],xlim=[1949.5,2021.5])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
	ax[1].legend(loc='center left',facecolor=[1,1,1],frameon=False)
	gu.axletters(ax,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	plt.tight_layout();
	if meta['Graphics']['Print Figures']=='On':
		#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area\NSR_Harvest_Summary','png',900)
		pass

	return

#%% Plot Time series and by BGC zone
def Plot_Wildfire_Po_By_BGCZone(meta):

	wfss=gu.ipickle(meta['Paths']['DB']['Taz'] + '\\Wildfire\\Wildfire_Stats_Scenarios_By_BGCZ.pkl')
	plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(22,9))
	ax[0].bar(wfss['tv_obs'],wfss['A_Burned_Tot']/1e6,facecolor=meta['Graphics']['Colours']['rgb']['Blue Light'])
	ax[0].set(xticks=np.arange(1920,2040,10),ylabel='Area affected (Million ha yr$^-$$^1$)',xlabel='Time, years',xlim=[1919,2025])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
	ax[1].bar(np.arange(1,wfss['Po_obs_mu'].size+1),wfss['Po_obs_mu']*100,facecolor=meta['Graphics']['Colours']['rgb']['Blue Light'])
	ax[1].set(xticks=np.arange(1,wfss['Po_obs_mu'].size+1),
		xticklabels=wfss['bgcz_cdS'],xlabel='BGC zone',ylabel='Probability of occurrence (% yr$^-$$^1$)',xlim=[0.5,wfss['Po_obs_mu'].size+.5],ylim=[0,0.65])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
	plt.tight_layout()
	gu.axletters(ax,plt,0.015,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Wildfire_Po_By_BGCZone','png',900)

	d={}
	d['Zone']=wfss['bgcz_cdS']
	d['Po']=wfss['Po_obs_mu']*100
	pd.DataFrame.from_dict(d)

	return

#%%
def Plot_Wildfire_Po_By_BGCZoneAndNDT(meta):
	#d=gu.ReadExcel(meta['Paths']['bc1ha'] + '\\LUTs\\LUT_bgcz_ndt_combo.xlsx')
	d=meta['LUT']['Raw']['bgc-ndt']

	ord=np.argsort(d['Po (%/yr)'])
	for k in d.keys():
		d[k]=d[k][ord]

	ind=np.where( (d['Po (%/yr)']>0) & (d['Area']>=1000) )[0]
	for k in d.keys():
		d[k]=d[k][ind]

	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,6.75))
	plt.bar(np.arange(d['ID'].size),d['Po (%/yr)'],facecolor=[0.8,0.8,0.8])
	 #ax.set_xticklabels(rotation = 45)
	for tick in ax.get_xticklabels():
		tick.set_rotation(45)
	for i in range(d['ID'].size):
		ax.text(i,d['Po (%/yr)'][i]+0.03,int(d['Area'][i]/1000),rotation=90,ha='center',fontsize=meta['Graphics']['gp']['fs_m'])
	ax.set(xlim=[-0.75,d['ID'].size],xticks=np.arange(d['ID'].size),xticklabels=d['BGC-NDT'],xlabel='BGC zone / NDT combination',ylim=[0,0.65],
		   ylabel='Annual probability of wildfire (%/yr)')
	ax.tick_params(axis='both',labelsize=meta['Graphics']['gp']['fs_s'])
	plt.tight_layout()
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Probability of Wildfire Historical','png',900)
		pass
	return

#%%
def Plot_FelledFate_Scenarios(meta,scn,reg):
	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\Variables_FelledFate.pkl')
	for k1 in d.keys():
		if k1=='Year':
			continue
		for k2 in d[k1].keys():
			for k3 in d[k1][k2].keys():
				d[k1][k2][k3]=d[k1][k2][k3]*100;
	#reg='Interior'; #reg='GFS22'
	#scn='S2'
	aw=0.26;
	plt.close('all');
	fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(22,10));
	ax[0,0].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,0].plot(d['Year'],d[scn][reg]['GreenStemMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,0].set(position=[0.065,0.69,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5);
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False);

	ax[0,1].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,1].plot(d['Year'],d[scn][reg]['GreenStemMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,1].set(position=[0.395,0.69,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5);

	ax[0,2].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,2].plot(d['Year'],d[scn][reg]['GreenStemMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,2].set(position=[0.73,0.69,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=1.5);

	ax[1,0].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,0].plot(d['Year'],d[scn][reg]['GreenStemNonMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,0].set(position=[0.065,0.385,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5);

	ax[1,1].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,1].plot(d['Year'],d[scn][reg]['GreenStemNonMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,1].set(position=[0.395,0.385,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5);

	ax[1,2].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,2].plot(d['Year'],d[scn][reg]['GreenStemNonMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,2].set(position=[0.73,0.385,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=1.5);

	ax[2,0].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,0].plot(d['Year'],d[scn][reg]['DeadStemMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,0].set(position=[0.065,0.08,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=1.5);

	ax[2,1].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,1].plot(d['Year'],d[scn][reg]['DeadStemMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,1].set(position=[0.395,0.08,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=1.5);

	ax[2,2].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,2].plot(d['Year'],d[scn][reg]['DeadStemMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,2].set(position=[0.73,0.08,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,2].yaxis.set_ticks_position('both'); ax[2,2].xaxis.set_ticks_position('both'); ax[2,2].tick_params(length=1.5);
	lab=['Merchantable','Merchantable','Merchantable','Residuals','Residuals','Residuals','Snags','Snags','Snags'];

	gu.axletters(ax,plt,0.035,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.03)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Felled Fate Scenarios\Felled Fate Scenarios ' + reg,'png',900)
	return

#%%
def Plot_FelledFate_Scenarios(meta,scn,reg):
	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\Variables_FelledFate.pkl')
	for k1 in d.keys():
		if k1=='Year':
			continue
		for k2 in d[k1].keys():
			for k3 in d[k1][k2].keys():
				d[k1][k2][k3]=d[k1][k2][k3]*100;
	#reg='Interior'; #reg='GFS22'
	#scn='S2'
	aw=0.26;
	plt.close('all');
	fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(22,10));
	ax[0,0].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,0].plot(d['Year'],d[scn][reg]['GreenStemMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,0].set(position=[0.065,0.69,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5);
	ax[0,0].legend(loc='lower right',facecolor=[1,1,1],frameon=False);

	ax[0,1].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,1].plot(d['Year'],d[scn][reg]['GreenStemMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,1].set(position=[0.395,0.69,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5);

	ax[0,2].plot(d['Year'],d['BaseCase'][reg]['GreenStemMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[0,2].plot(d['Year'],d[scn][reg]['GreenStemMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[0,2].set(position=[0.73,0.69,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=1.5);

	ax[1,0].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,0].plot(d['Year'],d[scn][reg]['GreenStemNonMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,0].set(position=[0.065,0.385,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5);

	ax[1,1].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,1].plot(d['Year'],d[scn][reg]['GreenStemNonMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,1].set(position=[0.395,0.385,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5);

	ax[1,2].plot(d['Year'],d['BaseCase'][reg]['GreenStemNonMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[1,2].plot(d['Year'],d[scn][reg]['GreenStemNonMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[1,2].set(position=[0.73,0.385,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xticklabels='',xlabel='',xlim=[2000,2100]);
	ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=1.5);

	ax[2,0].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchRemoved'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,0].plot(d['Year'],d[scn][reg]['DeadStemMerchRemoved'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,0].set(position=[0.065,0.08,aw,0.29],ylabel='Removal (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=1.5);

	ax[2,1].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchPiled'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,1].plot(d['Year'],d[scn][reg]['DeadStemMerchPiled'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,1].set(position=[0.395,0.08,aw,0.29],ylabel='Piled (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=1.5);

	ax[2,2].plot(d['Year'],d['BaseCase'][reg]['DeadStemMerchLeftOnSite'],'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1']);
	ax[2,2].plot(d['Year'],d[scn][reg]['DeadStemMerchLeftOnSite'],'g--',label=scn,lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2']);
	ax[2,2].set(position=[0.73,0.08,aw,0.29],ylabel='Left on site (%)',xticks=np.arange(1500,2200,20),xlabel='Time, years',xlim=[2000,2100]);
	ax[2,2].yaxis.set_ticks_position('both'); ax[2,2].xaxis.set_ticks_position('both'); ax[2,2].tick_params(length=1.5);
	lab=['Merchantable','Merchantable','Merchantable','Residuals','Residuals','Residuals','Snags','Snags','Snags'];

	gu.axletters(ax,plt,0.035,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.03)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Felled Fate Scenarios\Felled Fate Scenarios ' + reg,'png',900)
	return


#%% 
def Plot_RemovedFate_Scenarios(meta,scn,reg):

	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\Variables_RemovedFate.pkl')
	#reg='Interior'  #reg='GFS22'
	#sc='S2'

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,10));
	ax[0,0].plot(d['Year'],d['BaseCase'][reg]['RemovedMerchToLumberMill'],'k-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[0,0].plot(d['Year'],d[scn][reg]['RemovedMerchToLumberMill'],'k--',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[0,0].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
	ax[0,0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

	ax[0,1].plot(d['Year'],d['BaseCase'][reg]['RemovedMerchToPulpMill'],'k-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[0,1].plot(d['Year'],d[scn][reg]['RemovedMerchToPulpMill'],'k--',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[0,1].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)

	ax[1,0].plot(d['Year'],d['BaseCase'][reg]['RemovedMerchToPelletMill'],'k-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[1,0].plot(d['Year'],d[scn][reg]['RemovedMerchToPelletMill'],'k--',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[1,0].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)

	ax[1,1].plot(d['Year'],d['BaseCase'][reg]['RemovedMerchToPlywoodMill']+d['BaseCase'][reg]['RemovedMerchToOSBMill'],'k-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[1,1].plot(d['Year'],d[scn][reg]['RemovedMerchToPlywoodMill']+d[scn][reg]['RemovedMerchToOSBMill'],'k--',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[1,1].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction (%)',xticks=np.arange(1500,2200,25),xlabel='Time, years',ylim=[0,1],xlim=[1900,2100])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)
	lab=['Merchantable fibre to lumber mill','Merchantable fibre to pulp mill','Merchantable fibre to pellet mill','Merchantable fibre to plywood mill']
	#gu.axletters(ax,plt,0.028,0.9,LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight']); #
	gu.axletters(ax,plt,0.035,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.025)
	fig.tight_layout();
	if meta['Graphics']['Print Figures']=='On':
		#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Removed Fate Scenarios\Removed Fate Scenarios ' + reg + 'Scenario ' + str(sc),'png',900)
		pass
	return

#%%
def Plot_HWP_EndUse_Scenarios(meta,scn,reg):

	d=gu.ipickle(meta['Paths']['DB']['Harvest'] + '\\Variables_ProductTypes.pkl')
	#reg='GFS22'

	# Pool groups
	BuildingL=['SFH','MFH','Com']
	NonBuildingL=['PulpMill','MDFMill','PelletMill','PowerFacility','IPP','LogExport','Furn','Ship','Repairs','Other']

	nm='LumberMillTo'

	yBD=np.zeros(d['Year'].size)
	yBS1=np.zeros(d['Year'].size)
	yBS2=np.zeros(d['Year'].size)
	for ip in BuildingL:
		yBD=yBD+d['BaseCase'][reg][nm + ip]
		yBS1=yBS1+d['S1'][reg][nm + ip]
		yBS2=yBS2+d['S2'][reg][nm + ip]

	yNBD=np.zeros(d['Year'].size)
	yNBS1=np.zeros(d['Year'].size)
	yNBS2=np.zeros(d['Year'].size)
	for ip in NonBuildingL:
		yNBD=yNBD+d['BaseCase'][reg][nm + ip]
		yNBS1=yNBS1+d['S1'][reg][nm + ip]
		yNBS2=yNBS2+d['S2'][reg][nm + ip]

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,6));
	ax[0].plot(d['Year'],yBD,'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[0].plot(d['Year'],yBS1,'g--',label='Scenario 1',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[0].plot(d['Year'],yBS2,'r-.',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl3'])
	ax[0].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction',xticks=np.arange(1500,2200,25),xlabel='Time, years',xlim=[1900,2100],ylim=[0,1])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
	ax[0].legend(loc='upper right',facecolor=[1,1,1],frameon=False)
	ax[1].plot(d['Year'],yNBD,'b-',label='BaseCase',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl1'])
	ax[1].plot(d['Year'],yNBS1,'g--',label='Scenario 1',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl2'])
	ax[1].plot(d['Year'],yNBS2,'r-.',label='Scenario 2',lw=meta['Graphics']['gp']['lw3'],color=meta['Graphics']['gp']['cl3'])
	ax[1].set(yticks=np.arange(0,1.2,0.2),ylabel='Transfer fraction',xticks=np.arange(1500,2200,25),xlabel='Time, years',xlim=[1900,2100],ylim=[0,1])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
	lab=['Fibre to building materials','Fibre to non-building materials']
	gu.axletters(ax,plt,0.028,0.92,LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.025) #
	fig.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\End Use Scenarios\End Use Scenarios ' + reg,'png',900)
		pass
	return

#%% Plot deteriministic component of wildfire scenarios
def Plot_WildfireScenarios(meta,zone):
	wfss=gu.ipickle(meta['Paths']['DB']['Taz'] + '\\Wildfire\Wildfire_Stats_Scenarios_By_BGCZ.pkl')

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,7)); lw=1
	rc=patches.Rectangle((1920,0),103,20,facecolor=[0.92,0.92,0.92])
	ax.add_patch(rc)
	# Observations
	ax.plot(wfss['tv_obs'],wfss['ByZone'][zone]['Area Wildfire']/wfss['ByZone'][zone]['Area Zone']*100,'k.',ms=2,lw=0.25,label='Observations')
	# Scenario 1
	ax.plot(wfss['tv_scn'],wfss['ByZone'][zone]['Po Det Scenarios']['H1']['F1'],'b--',color=[1,0.5,0],lw=lw,label='Scenario H1F1')
	# Scenario 2
	ax.plot(wfss['tv_scn'],wfss['ByZone'][zone]['Po Det Scenarios']['H2']['F2'],'b-.',color=[1,0,0],lw=lw,label='Scenario H2F2')
	# Scenario 3
	ax.plot(wfss['tv_scn'],wfss['ByZone'][zone]['Po Det Scenarios']['H3']['F3'],'b:',color=[0.45,0,0],lw=lw,label='Scenario H3F3')
	# Scenario 4
	ax.plot(wfss['tv_scn'],wfss['ByZone'][zone]['Po Det Scenarios']['H2']['F0'],'b-',color=[0.47,0.69,0.97],lw=2,label='Scenario H2F0')
	ax.plot(wfss['tv_scn'],wfss['ByZone'][zone]['Po Det Scenarios']['H2']['F2'],'b-.',color=[1,0,0],lw=lw)
	
	ax.annotate('Observation\nperiod',(1920+50,2.25),ha='center',fontweight='bold',color=[0.5,0.5,0.5])
	ax.legend(loc='upper left',frameon=False)
	ax.set(xticks=np.arange(wfss['tv_scn'][0],wfss['tv_scn'][-1]+100,100),ylabel='Probability of occurrence (% yr$^-$$^1$)',xlabel='Time, calendar year',
		   yscale='linear',yticks=np.arange(0,6,1),ylim=[0,5],xlim=[500-0.5,wfss['tv_scn'][-1]+1+0.5])
	#ax.get_yaxis().set_minor_formatter(ScalarFormatter())
	#ax.set_ylim([0,5])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)	
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Wildfire_Scenarios_ts_' + zone,'png',900)
	return


#%% Bar chart of total area burned
def Plot_AreaBurned(meta):
	lNam='PROT_HISTORICAL_FIRE_POLYS_SP'
	vNam='FIRE_YEAR'
	tv,A=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,6)

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6))
	ax.bar(tv,A/1e6,0.75,facecolor=[0.7,0,0])
	ax.set(xticks=np.arange(1920,2040,10),ylabel='Affected area (Million hectares/year)',xlabel='Time, years',yticks=np.arange(0,2.2,0.2),xlim=[1919.25,2023.75])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\Area Burned','png',900)
		#gu.PrintFig(PathFigures + '\\Wildfire_AreaAffectedal','png',900)
	return tv,A

#%%
def Plot_AIL_Silv(meta,cd):
	ail=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + cd + '_SummaryByTimeAndFSC.pkl')
	cl=np.random.random((ail['FSC Unique'].size,3))
	ylim=[0,1.1*np.max(ail['A Total'])/1000]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,5.75));
	A_Cumulative=np.zeros(ail['Year'].size)
	for iFSC in range(ail['FSC Unique'].size):
		plt.bar(ail['Year'],ail['A Unique'][:,iFSC]/1000,0.8,bottom=A_Cumulative,facecolor=cl[iFSC,:],label=ail['FSC Unique'][iFSC])
		A_Cumulative=A_Cumulative+ail['A Unique'][:,iFSC]/1000
	ax.set(xticks=np.arange(1950,2025+1,5),ylabel='Treatment area (hectares x 1000)',xlabel='Time, years',xlim=[1975-0.5,ail['Year'][-1]+0.5],ylim=ylim)
	plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=3)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()

	#if meta['Graphics']['Print Figures']=='On':
	#	gu.PrintFig(meta['Paths'][pNamC]['Figures'] + '\\AIL_' + kwargs['Period'],'png',900)
	return

#%%
def PlotBeetleComp1_TS(meta):
	#tv=np.arange(1951,2025,1)
	#name=meta['Param']['Raw']['InsectComp1']['Insect Name']
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_SummaryTS.pkl')
	tv=d['tv']
	lab=['Mountain Pine Beetle','Balsam Beetle','Douglas-fir Beetle','Spruce Beetle']
	bw=0.8;
	#cl=np.array([[0.85,0.85,1],[0.5,0.5,1],[0.25,0.25,0.9],[0,0,0.3]])
	c=np.array(meta['Graphics']['Colours']['rgb']['Blue Light'])
	cl=np.array([c,0.75*c,0.5*c,0.25*c])

	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,7.75))
	s='Mountain Pine Beetle'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[0,0].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:],label='Light')
	ax[0,0].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0],label='Moderate')
	ax[0,0].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1],label='Severe')
	ax[0,0].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2],label='Very severe')
	ax[0,0].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)
	ax[0,0].legend(loc='lower left',frameon=False,facecolor=None,edgecolor='w')
	s='Balsam Beetle'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[0,1].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:])
	ax[0,1].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0])
	ax[0,1].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1])
	ax[0,1].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2])
	ax[0,1].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)
	s='Douglas-fir Beetle'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[1,0].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:])
	ax[1,0].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0])
	ax[1,0].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1])
	ax[1,0].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2])
	ax[1,0].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)
	s='Spruce Beetle'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[1,1].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:])
	ax[1,1].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0])
	ax[1,1].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1])
	ax[1,1].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2])
	ax[1,1].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.025)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics'][ 'Print Figure Path'] + '\\AreaBeetleTimeSeries','png',900)
	return

#%%
def PlotDefoliatorComp1_TS(meta):
	#tv=np.arange(1951,2025,1)
	#name=meta['Param']['Raw']['InsectComp1']['Insect Name']
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\InsectComp1_SummaryTS.pkl')
	tv=d['tv']
	lab=['Western Spruce Budworm','Hemlock Looper']
	bw=0.8

	#cl=np.array([[0.85,0.85,1],[0.5,0.5,1],[0.25,0.25,0.9],[0,0,0.3]])
	c=np.array(meta['Graphics']['Colours']['rgb']['Blue Light'])
	cl=np.array([c,0.75*c,0.5*c,0.25*c])

	plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,4.7))
	s='Western Spruce Budworm'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[0].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:],label='Light')
	ax[0].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0],label='Moderate')
	ax[0].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1],label='Severe')
	ax[0].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2],label='Very severe')
	ax[0].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
	ax[0].legend(loc='center left',frameon=False,facecolor=None,edgecolor='w')
	s='Hemlock Looper'
	y=np.column_stack( (d['data'][s]['L'],d['data'][s]['M'],d['data'][s]['S'],d['data'][s]['V']) )/1000
	ax[1].bar(tv,y[:,0],bw,fc=cl[0,:],ec=cl[0,:])
	ax[1].bar(tv,y[:,1],bw,fc=cl[1,:],ec=cl[1,:],bottom=y[:,0])
	ax[1].bar(tv,y[:,2],bw,fc=cl[2,:],ec=cl[2,:],bottom=y[:,0]+y[:,1])
	ax[1].bar(tv,y[:,3],bw,fc=cl[3,:],ec=cl[3,:],bottom=y[:,0]+y[:,1]+y[:,2])
	ax[1].set(xticks=np.arange(1950,2025,5),ylabel='Area affected\n(ha x 1000 yr$^-$$^1$)',xlabel='Time, calendar year',xlim=[1950,2025])
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
	plt.tight_layout()
	gu.axletters(ax,plt,0.025,0.87,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'],Labels=lab,LabelSpacer=0.025)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics'][ 'Print Figure Path'] + '\\AreaDefoliatorTimeSeries','png',900)
	return

#%%
def AddParagraphs(txt,c):
	for i in range(txt['Counter'],txt['Counter']+c):
		if txt['style'][i]=='Title':
			txt0='''<p style="color:black;text-align:center;font-size:50px;padding-top:5px">''' + txt['body'][i] + '</p>'
		elif txt['style'][i]=='Subtitle':
			txt0='''<p style="color:black;text-align:center;font-size:16px;padding-top:5px">''' + txt['body'][i] + '</p>'
		elif txt['style'][i]=='Authorship':
			txt0='''<p style="color:black;text-align:center;font-size:12px;padding-top:5px">''' + txt['body'][i] + '</p>'
		elif txt['style'][i]=='Heading 1':
			txt0='<h1>' + txt['body'][i] + '</h1>'
		elif txt['style'][i]=='Heading 2':
			txt0='<h2>' + txt['body'][i] + '</h2>'
		elif txt['style'][i]=='Heading 3':
			txt0='<h3>' + txt['body'][i] + '</h3>'
		elif txt['style'][i]=='Heading 4':
			txt0='<h4>' + txt['body'][i] + '</h4>'
		elif txt['style'][i]=='List Paragraph':
			txt0='''<li style="margin-left:20px;margin-bottom:1px;margin-top:1px;">''' + txt['body'][i] + ''''</li>'''
		elif txt['style'][i]=='Header':
			txt0='<b><i>' + txt['body'][i] + '</i></b>'
		elif txt['style'][i]=='Box':
			txt0='''<p style="color:#094985;font-weight:normal;text-align:left;background-color:#e0f0ff;padding-top:5px;padding-bottom:5px;padding-left:5px;padding-right:5px">''' + txt['body'][i] + '</p>'
		else:
			txt0='''<p style="color:black;padding-bottom:0px">''' + txt['body'][i] + '</p>'
		display(HTML(txt0))
	txt['Counter']=txt['Counter']+c
	return txt

#%%
def ImportText(meta,pNam):
	# Import body text	
	doc=docx.Document(meta['Paths'][pNam]['Text'])
	txt={'body':[],'style':[],'hyperlinks':[],'Counter':0}
	for para in doc.paragraphs:
		txt['body'].append(para.text)
		txt['style'].append(para.style.name)
		#if para.hyperlinks==[]:
		#	txt['hyperlinks'].append(para.hyperlinks[0].address)
		#else:
		#	txt['hyperlinks'].append([])
	return txt

#%%
def Plot_FateOfRemovals(meta):
	fnam='Figure ' + str(meta['Graphics']['Fig Count']) + ' Schematic Fate of Harvest Removals'
	g=graphviz.Digraph('Fibre Flow',filename=fnam,
					 graph_attr={'rankdir':'TB',
								 'compound':'true',
								 'pad':'0.1',
								 'nodesep':'0.1'},
					 node_attr={'shape':'box',
								'fixedsize':'true',
								'width':'1.2',
								'height':'0.5',
								'fontname':meta['Graphics']['Flowchart']['Font Name'],
								'fontcolor':meta['Graphics']['Flowchart']['Font Color'],
								'fontsize':meta['Graphics']['Flowchart']['Font Size'],
								'style':'filled',
								'fillcolor':meta['Graphics']['Flowchart']['Node Background Color'],
								'penwidth':meta['Graphics']['Flowchart']['Penwidth'],
								'center':'true'},
					edge_attr={})
	
	g.edge_attr.update(arrowhead='normal',arrowsize='0.6',penwidth='0.5',fontsize='9',fontname=meta['Graphics']['Flowchart']['Font Name'])
	
	tot='Removals'
	lumb='Lumber mill'
	pulp='Pulp mill'
	plyw='Plywood mill'
	pane='OSB/MDF mill'
	pelm='Pellet mill'
	loge='Log exports'
	
	tot2lumb='X%\l'
	tot2pulp='X%\l'
	tot2plyw='X%\l'
	tot2pane='X%\l'
	tot2pelm='X%\l'
	tot2loge='X%\l'
	
	soli='Solid products'
	pape='Paper'
	facp='Facility power'
	grid='Grid power'
	pell='Pellet exports'
	
	lumb2soli='X%\l'
	lumb2facp='X%\l'
	lumb2grid='X%\l'
	lumb2pulp='X%\l'
	lumb2pane='X%\l'
	lumb2pelm='X%\l'
	pulp2pape='X%\l'
	pulp2facp='X%\l'
	pulp2grid='X%\l'
	plyw2soli='X%\l'
	pane2soli='X%\l'
	pelm2pell='X%\l'
	
	#eco.node(tot)
	with g.subgraph(name='cluster_mill') as mil:
		mil.attr(label='Mills',labeljust="l",style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		mil.node(lumb,shape='house')
		mil.node(pulp,shape='house')
		mil.node(plyw,shape='house')
		mil.node(pane,shape='house')
		mil.node(pelm,shape='house')
	with g.subgraph(name='cluster_prod') as pro:
		pro.attr(label='Products',labeljust="r",style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		pro.node(soli)
		pro.node(pape)
		pro.node(facp)
		pro.node(grid)
		pro.node(pell)
		pro.node(loge)
	g.edge(tot,lumb,label=tot2lumb)
	g.edge(tot,pulp,label=tot2pulp)
	g.edge(tot,plyw,label=tot2plyw)
	g.edge(tot,pane,label=tot2pane)
	g.edge(tot,pelm,label=tot2pelm)
	g.edge(tot,loge,label=tot2loge)
	g.edge(lumb,soli,label=lumb2soli)
	g.edge(lumb,facp,label=lumb2facp)
	g.edge(lumb,pulp,label=lumb2pulp)
	g.edge(lumb,pane,label=lumb2pane)
	g.edge(lumb,pelm,label=lumb2pelm)
	g.edge(pulp,pape,label=pulp2pape)
	g.edge(pulp,facp,label=pulp2facp)
	g.edge(pulp,grid,label=pulp2grid)
	g.edge(plyw,soli,label=plyw2soli)
	g.edge(pane,soli,label=pane2soli)
	
	g.edge(pelm,pell,label=pelm2pell)
	
	display(g)
	return

#%%
def Plot_Tier2Categories(meta,capt):
	fig=str(meta['Graphics']['Fig Count'])	
	g=graphviz.Graph('Taxonomy',filename='Figure ' + fig + ' ' + capt,
					 graph_attr={'rankdir':'TB',
								 'compound':'true',
								 'pad':'0.1',
								 'nodesep':'0.1'},
					 node_attr={'shape':'box',
										   'fixedsize':'true',
										   'width':'1.5',
										   'height':'0.7',
										   'fontname':meta['Graphics']['Flowchart']['Font Name'],
										   'fontsize':meta['Graphics']['Flowchart']['Font Size'],
										   'fontcolor':meta['Graphics']['Flowchart']['Font Color'],
										   'style':'filled',
										   'fillcolor':meta['Graphics']['Flowchart']['Node Background Color'],
										   'penwidth':meta['Graphics']['Flowchart']['Penwidth'],
										   'center':'true'})
	
	act=meta['Tables']['BCFCS_ForcingCategories']['Category'][0]
	ant=meta['Tables']['BCFCS_ForcingCategories']['Category'][1]
	iant=meta['Tables']['BCFCS_ForcingCategories']['Category'][2]
	nat=meta['Tables']['BCFCS_ForcingCategories']['Category'][3]
	lcc=meta['Tables']['BCFCS_ForcingCategories']['Category'][4]
	bau=meta['Tables']['BCFCS_ForcingCategories']['Category'][5]
	rea=meta['Tables']['BCFCS_ForcingCategories']['Category'][6]
	ih=meta['Tables']['BCFCS_ForcingCategories']['Category'][7]
	sil=meta['Tables']['BCFCS_ForcingCategories']['Category'][8]
	dis=meta['Tables']['BCFCS_ForcingCategories']['Category'][9]
	pro=meta['Tables']['BCFCS_ForcingCategories']['Category'][10]
	
	with g.subgraph(name='cluster_0') as c:
		c.attr(style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		c.node_attr.update(style='filled',color='white')
		c.node(nat)
		c.node(ant)
		c.node(iant)
		c.attr(label='Level 2',labeljust="l")
		
	with g.subgraph(name='cluster_1') as c:
		c.attr(style='filled',color=meta['Graphics']['Flowchart']['Cluster Background Color'],fontname=meta['Graphics']['Flowchart']['Font Name'])
		c.node_attr.update(style='filled', color='white')
		c.node(lcc)
		c.node(bau)
		c.node(rea)
		c.node(ih)
		c.node(sil)
		c.node(dis)
		c.node(pro)
		c.attr(label='Level 3',labeljust="l")
	
	g.edge(act,ant)
	g.edge(act,iant)
	g.edge(act,nat)
	
	g.edge(ant,lcc)
	g.edge(ant,bau)
	g.edge(ant,rea)
	g.edge(ant,ih)
	g.edge(ant,sil)
	g.edge(ant,dis)
	g.edge(ant,pro)
	display(g)
	return

#%%
def Plot_FlowChartHarvestByTimberMarkDB(meta,capt):
	fig=str(meta['Graphics']['Fig Count'])	
	g=graphviz.Graph('Harvest DB',filename='Figure ' + fig + ' ' + capt,
					 graph_attr={'rankdir':'TB',
								 'compound':'true',
								 'pad':'0.05',
								 'nodesep':'0.05'},
					 node_attr={'shape':'box',
										   'fixedsize':'true',
										   'width':'1.25',
										   'height':'0.7',
										   'fontname':meta['Graphics']['Flowchart']['Font Name'],
										   'fontsize':meta['Graphics']['Flowchart']['Font Size'],
										   'fontcolor':meta['Graphics']['Flowchart']['Font Color'],
										   'style':'filled',
										   'fillcolor':meta['Graphics']['Flowchart']['Node Background Color'],
										   'penwidth':meta['Graphics']['Flowchart']['Penwidth'],
										   'center':'true'})

	g.node('db1','Consolidated\nCutblocks',shape='cylinder',style='filled')
	g.node('db2','Forest\nTenures',shape='cylinder',style='filled')
	g.node('db3','RESULTS',shape='cylinder',style='filled')
	g.node('db4','Harvest Billing\nSystem',shape='cylinder',style='filled')
	g.node('db5','Timber Cruise\nCompilation',shape='cylinder',style='filled')
	g.node('db6','Waste & Residues\nSystem',shape='cylinder',style='filled')
	g.node('db7','BC Mill\nSurveyReport',shape='cylinder',style='filled')

	g.node('p1','Compilation',shape='ellipse')

	g.node('ByOP','Harvest by\nOpening\nDB',shape='cylinder',style='filled')
	g.node('ByTM','Harvest by\nTimber Mark\nDB',shape='cylinder',style='filled')

	g.node('REG Table','Harvest Stats\nby Region',shape='cylinder',style='filled')
	g.node('BGC Table','Harvest Stats\nby BGC zone',shape='cylinder',style='filled')
	g.node('SSC Table','Harvest Stats\nby SSC',shape='cylinder',style='filled')
	g.node('TS Table','Harvest Stats\nby Time',shape='cylinder',style='filled')

	for n in ['db1','db2','db3','db4','db5','db6','db7']:
		g.edge(n,'p1')
	g.edge('p1','ByOP')
	g.edge('p1','ByTM')
	g.edge('ByTM','REG Table')
	g.edge('ByTM','BGC Table')
	g.edge('ByTM','SSC Table')
	g.edge('ByTM','TS Table')
	display(g)
	return

#%%
def Tabulate_Prices(meta,year):
	dCP=gu.ReadExcel(meta['Paths']['Model']['Code'] + '\\Parameters\\Parameters_CostsAndPrices.xlsx',sheet_name='Data')
	iT=np.where(dCP['Year']==year)[0]
	d={}
	for k in dCP.keys():
		if k[0:5]=='Price':
			d[k[5:]]=dCP[k][iT]
	df=pd.DataFrame.from_dict(d)
	df=df.rename(index={0:str(dCP['Year'][iT[0]])}).T
	return d,df

#%% Plot map of GHG balance
# *** MOVED TO ROI_MAP ***
# def PlotMap(meta,pNam,iScn,mu_mod,v):

#%%
def PlotMapOfSparseXYSample(meta,pNam):
	if meta['Geos']['Sparse']['X'].size<50000:
		gdf_bcb=gpd.read_file(meta['Paths']['GDB']['LandUse'],layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
		gdf_sxy=gpd.read_file(meta['Paths'][pNam]['Data'] + '\\geos.geojson')
		plt.close('all'); fig,ax=plt.subplots(figsize=gu.cm2inch(16,16)); ms=2
		gdf_bcb.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
		gdf_sxy.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
		#ax.grid(color='none',linestyle='-',linewidth=0.25)
		ax.set(position=[0,0,1,1],xticks=[],yticks=[],xlim=meta['Geos']['Grid']['xlim'],ylim=meta['Geos']['Grid']['ylim'])
		ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.grid(meta['Graphics']['Map']['Map Grid Vis']); ax.axis(meta['Graphics']['Map']['Map Axis Vis'])
		plt.tight_layout()
		#gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Map_GroundPlots','png',900)
		#plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)
		#plt.close('all')
		#garc.collect()
	return

#%%
def Plot_AIL_PileBurn(meta,pNam):
	vNam='SP-BU-LAND'
	vNam='SP-PILE'
	vNam='SP-PBURN'
	#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PileBurn_SummaryByTimeAndFSC.pkl')
	#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PBURN_SummaryByTimeAndFSC.pkl')
	#d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\SP-PILE_SummaryByTimeAndFSC.pkl')
	d=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')

	iT=np.where( (d['Year']>=1960) & (d['Year']<=2050) )[0]
	cl=np.random.random((d['FSC Unique'].size,3))

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9));
	A_Cumulative=np.zeros(d['Year'].size)
	for i in range(d['FSC Unique'].size):
		plt.bar(d['Year'],d['A Unique'][:,i]/1e3,0.8,bottom=A_Cumulative,facecolor=cl[i,:],label=d['FSC Unique'][i])
		A_Cumulative=A_Cumulative+d['A Unique'][:,i]/1e3
	#plt.plot(d['tv'],d['A'/1e3,'ks',ms=2.5,mec='k',mfc='w',mew=0.5)
	ax.set(xticks=np.arange(1950,2225+1,5),ylabel='Implementation level (Kha yr$^{-1}$)',
		   xlabel='Time, years',yticks=np.arange(0,300,5),xlim=[d['Year'][iT][0]-0.75,d['Year'][iT][-1]+0+.75],ylim=[0,35]) ##
	plt.legend(frameon=False,loc='upper right',facecolor=[1,1,1],labelspacing=0.25,ncol=6)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_AIL_' + vNam + '_ByFSC','png',900)

	# Mix of all types
	vNam='SP-BU-LAND'
	d1=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')
	vNam='SP-PILE'
	d2=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')
	vNam='SP-PBURN'
	d3=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndFSC.pkl')

	iT=np.where( (d1['Year']>=1960) & (d1['Year']<=2050) )[0]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9));
	plt.bar(d1['Year'],d1['A Total']/1e3,0.8,facecolor=[0.6,0.6,0.6],label='Landings (SP-BU-LAND)')
	plt.bar(d1['Year'],d2['A Total']/1e3,0.8,bottom=d1['A Total']/1e3,facecolor=[0.2,0.2,0.2],label='Pile burn (SP-BU-PBURN)')
	plt.bar(d1['Year'],d3['A Total']/1e3,0.8,bottom=(d1['A Total']+d2['A Total'])/1e3,facecolor=[0.85,0.85,0.85],label='Piling (SP-BU-PILE, SP-BU-RPILE)')
	ax.set(xticks=np.arange(1950,2225+1,5),ylabel='Implementation level (Kha yr$^{-1}$)',
		   xlabel='Time, years',yticks=np.arange(0,300,5),xlim=[d['Year'][iT][0]-0.75,d['Year'][iT][-1]+0+.75],ylim=[0,45]) ##
	plt.legend(frameon=False,loc='upper right',facecolor=[1,1,1],labelspacing=0.25,ncol=1)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_AIL_SP-BU-Summary','png',900)

	# By BGC zone - mix of all types
	d={}
	vNam='SP-BU-LAND'
	d[vNam]=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')
	vNam='SP-PILE'
	d[vNam]=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')
	vNam='SP-PBURN'
	d[vNam]=gu.ipickle(meta['Paths']['bc1ha'] + '\\RSLT_ACTIVITY_TREATMENT_SVW\\' + vNam + '_SummaryByTimeAndBGCZone.pkl')

	u=d1['BGC']
	for k in d.keys():
		cnt=0
		d[k]['Summary']=np.zeros(u.size)
		for zone in u:
			ind=np.where(d[k]['BGC']==zone)[0]
			d[k]['Summary'][cnt]=np.sum(d[k]['A'][:,ind])/1e3
			cnt=cnt+1

	zH=[None]*6
	for iY in range(6):
		zH[iY]=gis.OpenGeoTiff(r'C:\Data\BC1ha\VEG_CONSOLIDATED_CUT_BLOCKS_SP\HARVEST_YEAR_' + str(iY+1) + '_Year.tif')['Data']
	zZ=gis.OpenGeoTiff(r'C:\Data\BC1ha\BEC_BIOGEOCLIMATIC_POLY\ZONE_GapFilled.tif')
	A_Harvest=np.zeros(u.size)
	cnt=0
	for z in u:
		for iY in range(6):
			ind=np.where( (zZ['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][z]) & (zH[iY]>=1960) )
			A_Harvest[cnt]=A_Harvest[cnt]+ind[0].size
		cnt=cnt+1

	for k in d.keys():
		d[k]['Summary %']=d[k]['Summary']/A_Harvest*100
	
	v='Summary %'

	ord=np.argsort(d['SP-BU-LAND'][v]+d['SP-PBURN'][v]+d['SP-PILE'][v])
	u=u[ord]
	for k in d.keys():
		d[k][v]=d[k][v][ord]

	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,9));
	plt.bar(u,d['SP-BU-LAND'][v],0.8,facecolor=[0.6,0.6,0.6],label='Landings (SP-BU-LAND)')
	plt.bar(u,d['SP-PBURN'][v],0.8,bottom=d['SP-BU-LAND'][v],facecolor=[0.2,0.2,0.2],label='Pile burn (SP-BU-PBURN)')
	plt.bar(u,d['SP-PILE'][v],0.8,bottom=(d['SP-BU-LAND'][v]+d['SP-PBURN'][v]),facecolor=[0.85,0.85,0.85],label='Piling (SP-BU-PILE, SP-BU-RPILE)')
	ax.set(xticks=np.arange(1,u.size+1,1),ylabel='Implementation level (Kha yr$^{-1}$)',
		   xlabel='Time, years',yticks=np.arange(0,300,20),xlim=[0,u.size],ylim=[0,285]) ##
	plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=1)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + pNam + '_AIL_SP-BU_SummaryByBGCZone','png',900)

	return