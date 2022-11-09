
#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.patches import Rectangle
import gc as garc
import warnings
import time
#import matplotlib.colors
#from matplotlib import animation
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Plotting parameters

flag_printfigs='On'

gp=gu.SetGraphics('Presentation Dark')

#%% THLB
def Plot_THLB(meta,tv,AEF,thlb,iScn):

    #xlim=[1980,2040]; figs=[10,10]; xt=np.arange(1800,2200,10)
    xlim=[tv[0],tv[-1]]; figs=[14,8]; xt=np.arange(1800,2200,20)

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(figs[0],figs[1]));
    ax.plot(meta['Year'],AEF*np.sum(thlb[iScn]['Baseline WithDef'],axis=1)[0]*np.ones(meta['Year'].size)/1e6,'k-',color=gp['cl1'],label='Baseline')
    ax.plot(meta['Year'],AEF*np.sum(thlb[iScn]['Actual WithDef'],axis=1)/1e6,'k--',color=gp['cl2'],label='Actual')
    ax.set(xticks=xt,yticks=np.arange(0,26,2),ylabel='Timber harvesting landbase (million ha)',xlabel='Time, years',xlim=xlim,ylim=[0,24])
    ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    if flag_printfigs=='On':
        gu.PrintFig(meta['Paths']['Figures'] + '\\THLB','png',900)

    return

#%% Carbon pool summary

def PlotCarbonPoolBarChart(meta,tv,mos,AEF,iScn,iPS,iSS):

    iT2=[]
    iT2.append(np.where( (tv>=1800) & (tv<=1820) )[0])
    iT2.append(np.where( (tv>=1990) & (tv<=2020) )[0])
    iT2.append(np.where( (tv>=2100) & (tv<=2120) )[0])

    bw=0.25; bw2=0.3;
    cl=np.array([[0.65,0.75,0.75],[0.7,0.8,0.9],[0.8,0.9,1]]).T
    lab=['1820','1990','2100']

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,6));

    # Biomass
    nam='C_Biomass_Tot'; x0=1; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
        ax[0].bar(x0,y[i],0.01,color=cl[i,:],label=lab[i])
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    # Dead wood
    nam='C_DeadWood_Tot'; x0=2; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    # LiT2ter
    nam='C_Litter_Tot'; x0=3; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    # Soil
    nam='C_Soil_Tot'; x0=4; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    # In-use products
    nam='C_InUse_Tot'; x0=5; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    # Dump and landfill
    nam='C_DumpLandfill_Tot'; x0=6; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam]['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[0].bar(x,y,bw,color=cl)
    ax[0].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)

    ax[0].set(position=[0.07,0.125,0.7,0.86],xticks=np.arange(1,7,1),xticklabels=['Biomass','Dead \nwood','Litter','Soil','In-use\nproducts','Dumps & \nlandfills'],xlabel='',xlim=[0.5,6.5], \
      yticks=np.arange(0,120,2),ylabel='Carbon stock (GtC)',ylim=[0,11])
    ax[0].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=gp['tickl'])

    # Total
    nam='C_HWP_Tot'; x0=1; x=np.array([x0-bw2,x0,x0+bw2]); y=np.zeros(3); yl=np.zeros(3); yh=np.zeros(3)
    for i in range(3):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum']['C_HWP_Tot']['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)+np.mean((mos['Scenarios'][iScn]['Sum']['C_Forest_Tot']['Ensemble Mean'][iT2[i],iPS,iSS])/1e9*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum']['C_HWP_Tot']['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)+np.mean((mos['Scenarios'][iScn]['Sum']['C_Forest_Tot']['Ensemble P025'][iT2[i],iPS,iSS])/1e9*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum']['C_HWP_Tot']['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)+np.mean((mos['Scenarios'][iScn]['Sum']['C_Forest_Tot']['Ensemble P975'][iT2[i],iPS,iSS])/1e9*AEF)
    ax[1].bar(x,y,bw,color=cl)
    ax[1].errorbar(x,y,yerr=[y-yl,yh-y],color=0.8*cl[i,:],ls='',capsize=2)
    ax[1].set(position=[0.84,0.125,0.14,0.86],xticks=np.arange(1,2,1),xticklabels=['Total'],yticks=np.arange(0,120,2),ylabel='Carbon stock (GtC)',xlabel='',ylim=[0,26],xlim=[0.5,1.5])
    ax[1].legend(loc='upper left',facecolor=[1,1,1],frameon=False);
    ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=gp['tickl'])

    #ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_PoolSummaryBarChart_Scn' + str(iScn+1),'png',900)

    return

#%% flux summary

def PlotFluxesBarChart(meta,tv,mos,AEF,iScn,iPS,iSS):
    #list(mos['Scenarios'][0]['Sum'].keys())

    #iScn=0

    mos['Scenarios'][iScn]['Sum']['OPS']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations'].copy()
    mos['Scenarios'][iScn]['Sum']['OPS']['Ensemble Mean']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean']+mos['Scenarios'][iScn]['Sum']['E_CO2e_ET_Operations']['Ensemble Mean']+mos['Scenarios'][iScn]['Sum']['E_CO2e_IPPU_Operations']['Ensemble Mean']

    mos['Scenarios'][iScn]['Sum']['OPS']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations'].copy()
    mos['Scenarios'][iScn]['Sum']['OPS']['Ensemble P250']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations']['Ensemble P250']+mos['Scenarios'][iScn]['Sum']['E_CO2e_ET_Operations']['Ensemble P250']+mos['Scenarios'][iScn]['Sum']['E_CO2e_IPPU_Operations']['Ensemble P250']

    mos['Scenarios'][iScn]['Sum']['OPS']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations'].copy()
    mos['Scenarios'][iScn]['Sum']['OPS']['Ensemble P750']=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations']['Ensemble P750']+mos['Scenarios'][iScn]['Sum']['E_CO2e_ET_Operations']['Ensemble P750']+mos['Scenarios'][iScn]['Sum']['E_CO2e_IPPU_Operations']['Ensemble P750']

    nam=np.array(['E_CO2e_LULUCF_NEE','E_CO2e_LULUCF_Wildfire','E_CO2e_LULUCF_OpenBurning','E_CO2e_LULUCF_HWP','E_CO2e_ESC_Bioenergy','OPS', \
                  'E_CO2e_SUB_E','E_CO2e_SUB_M','E_CO2e_SUB_Tot','E_CO2e_AGHGB_WSub','E_CO2e_AGHGB_WOSub'])
    xlab=['Net\necosystem\nexchange','Wildfire\nemissions','Open\nburning\nemissions','Wood\nproduct\nemissions','Bioenergy\nemissions','Operational\nemissions', \
          'Energy\nsubst.','Material\nsubst.','Total\nsubst.','Net GHG\nbalance\n(with subs.)','Net GHG\nbalance\n(w/o subs.)']

    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,8.5)); bw=0.2; bw2=0.25;
    ax.plot([0,nam.size+1],[0,0],'-',color=gp['cla'])

    iT2=np.where( (tv>=1865) & (tv<=2021) )[0]
    y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS])/1e6*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS])/1e6*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS])/1e6*AEF)
    ax.bar(np.arange(1,y.size+1,1)-bw2,y,bw,color=[0.57,0.79,1],label='1865-2021')
    ax.errorbar(np.arange(1,y.size+1,1)-bw2,y,yerr=[y-yl,yh-y],color=[0.37,0.59,0.8],ls='',capsize=2)

    iT2=np.where( (tv>=1990) & (tv<=2021) )[0]
    y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS])/1e6*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS])/1e6*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS])/1e6*AEF)
    ax.bar(np.arange(1,y.size+1,1),y,bw,color=[0.8,1,0.6],label='1990-2021')
    ax.errorbar(np.arange(1,y.size+1,1),y,yerr=[y-yl,yh-y],color=[0.6,0.8,0.4],ls='',capsize=2)

    iT2=np.where( (tv>=2022) & (tv<=2100) )[0]
    y=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        y[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS])/1e6*AEF)
        yl[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P250'][iT2,iPS,iSS])/1e6*AEF)
        yh[i]=np.mean((mos['Scenarios'][iScn]['Sum'][nam[i]]['Ensemble P750'][iT2,iPS,iSS])/1e6*AEF)
    ax.bar(np.arange(1,y.size+1,1)+bw2,y,bw,color=[0.8,0.4,1],label='2023-2100')
    ax.errorbar(np.arange(1,y.size+1,1)+bw2,y,yerr=[y-yl,yh-y],color=[0.6,0.2,0.8],ls='',capsize=2)

    ax.set(position=[0.06,0.14,0.92,0.82],xlim=[0.5,nam.size+.5],xticks=np.arange(1,nam.size+1,1),xticklabels=xlab, \
           ylim=[-70,70],yticks=np.arange(-100,100,10),ylabel='Mean GHG flux (MtCO$_2$e yr$^-$$^1$)')
    ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

    #ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_FluxSummaryBarChart_Scn' + str(iScn+1),'png',900)

#%% flux summary

def PlotForcingBarChart(meta,tv,mos,AEF,cmp,iPS,iSS):

    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']={}
    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble Mean']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean'].copy()
    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble P025']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean'].copy()
    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble P975']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean'].copy()

    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble Mean']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ET_Operations']['Ensemble Mean']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_IPPU_Operations']['Ensemble Mean']
    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble P025']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble P025']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ET_Operations']['Ensemble P025']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_IPPU_Operations']['Ensemble P025']
    mos['Delta'][cmp]['ByStrata']['Sum']['OPS']['Ensemble P975']=mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ESC_Operations']['Ensemble P975']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_ET_Operations']['Ensemble P975']+mos['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_IPPU_Operations']['Ensemble P975']

    nam=np.array(['E_CO2e_LULUCF_NEE','E_CO2e_LULUCF_Wildfire','E_CO2e_LULUCF_OpenBurning','E_CO2e_LULUCF_HWP','E_CO2e_ESC_Bioenergy','OPS', \
                  'E_CO2e_SUB_E','E_CO2e_SUB_M','E_CO2e_SUB_Tot','E_CO2e_AGHGB_WSub','E_CO2e_AGHGB_WOSub'])

    xlab=['Net\necosystem\nexchange','Wildfire\nemissions','Open\nburning\nemissions','Wood\nproduct\nemissions','Bioenergy\nemissions','Operational\nemissions', \
          'Energy\nsubst.','Material\nsubst.','Total\nsubst.','Net GHG\nbalance\n(with subs.)','Net GHG\nbalance\n(w/o subs.)']

    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,8.5)); bw=0.2; bw2=0.25;
    ax.plot([0,nam.size+1],[0,0],'k-',color=gp['cla'])

    iT2=np.where( (tv>=1865) & (tv<1990) )[0]
    ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        ymu[i]=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS]/1e6*AEF)
        ylo=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS]/1e6*AEF)
        yhi=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS]/1e6*AEF)
        cs=np.array([ylo,yhi])
        yh[i]=np.max(cs)
        yl[i]=np.min(cs)
    ax.bar(np.arange(1,ymu.size+1,1)-bw2,ymu,bw,color=[0.57,0.79,1],label='1865-1990')
    ax.errorbar(np.arange(1,ymu.size+1,1)-bw2,ymu,yerr=[ymu-yl,yh-ymu],color=[0.37,0.59,0.8],ls='',capsize=2)

    iT2=np.where( (tv>=1990) & (tv<=2021) )[0]
    ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        ymu[i]=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS]/1e6*AEF)
        ylo=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS]/1e6*AEF)
        yhi=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS]/1e6*AEF)
        cs=np.array([ylo,yhi])
        yh[i]=np.max(cs)
        yl[i]=np.min(cs)
    ax.bar(np.arange(1,ymu.size+1,1),ymu,bw,color=[0.8,1,0.6],label='1991-2021')
    ax.errorbar(np.arange(1,ymu.size+1,1),ymu,yerr=[ymu-yl,yh-ymu],color=[0.6,0.8,0.4],ls='',capsize=2)

    iT2=np.where( (tv>=2022) & (tv<=2100) )[0]
    ymu=np.zeros(nam.size); yl=np.zeros(nam.size); yh=np.zeros(nam.size)
    for i in range(nam.size):
        ymu[i]=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble Mean'][iT2,iPS,iSS]/1e6*AEF)
        ylo=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P025'][iT2,iPS,iSS]/1e6*AEF)
        yhi=np.mean(mos['Delta'][cmp]['ByStrata']['Sum'][nam[i]]['Ensemble P975'][iT2,iPS,iSS]/1e6*AEF)
        cs=np.array([ylo,yhi])
        yh[i]=np.max(cs)
        yl[i]=np.min(cs)
    ax.bar(np.arange(1,ymu.size+1,1)+bw2,ymu,bw,color=[0.8,0.4,1],label='2023-2100')
    ax.errorbar(np.arange(1,ymu.size+1,1)+bw2,ymu,yerr=[ymu-yl,yh-ymu],color=[0.6,0.2,0.8],ls='',capsize=2)

    ax.set(position=[0.07,0.14,0.92,0.82],xlim=[0.5,nam.size+.5],xticks=np.arange(1,nam.size+1,1),xticklabels=xlab, \
           yticks=np.arange(-100,200,10),ylabel='Forcing (MtCO$_2$e yr$^-$$^1$)',ylim=[-100,100])
    ax.legend(loc='upper left',ncol=3,facecolor=[1,1,1],frameon=False);
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

    #ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_ForcingBarChart_' + cmp,'png',900)

#%% Harvest removals
def PlotHarvestVolume(meta,mos,tv,AEF,iScn,iT,iPS,iSS,flg_ann,flg_dead,flg_hbs):

    # Import FAIB public summary
    H_FLNR=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')

    H_HBS=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS_AnnualSummary.pkl')

    matype='historical'
    #matype='center'

    # Total harvest
    y_Tot=(mos['Scenarios'][iScn]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][:,iPS,iSS]+mos['Scenarios'][iScn]['Sum']['V_ToMillNonMerch']['Ensemble Mean'][:,iPS,iSS])/1e6*AEF
    y_Tot_ma=gu.movingave(y_Tot[iT],10,matype)

    y_Dead=(mos['Scenarios'][iScn]['Sum']['V_ToMillMerchDead']['Ensemble Mean'][:,iPS,iSS])/1e6*AEF
    y_Dead_ma=gu.movingave(y_Dead[iT],10,matype)

    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));

    if flg_ann==1:
        ax.plot(tv[iT],y_Tot[iT],'-',color=[0.6,0.9,0],lw=0.25,label='Tot. merch. vol. (FCS)')
    ax.plot(tv[iT],y_Tot_ma,'-',color=[0.3,0.75,0],lw=1,label='Tot. merch vol. mov. ave. (FCS)')

    if flg_dead==1:
        ax.plot(tv[iT],y_Dead_ma,'--',color=[0,0.7,0],lw=1,label='Dead merch. vol. (FCS)')

    if flg_hbs==1:
        ax.plot(H_FLNR['Year'],H_FLNR['Total_harvest_millions_m3'],'k-',color=gp['cl1'],lw=1.25,label='Tot. merch. vol. (FLNR)')
        ax.plot(H_HBS['Year'],H_HBS['V All (Mm3/yr)'],color=[0,1,1],lw=1.25,label='Tot. merch. vol. (HBS)')

    #ax.plot(tv[iT],y,'-',color=gp['cl2'],lw=0.5,label='Total')
    ax.set(position=[0.085,0.125,0.88,0.84],xlim=[tv[iT[0]],2100],xticks=np.arange(1800,2120,20),ylabel='Volume removed (Million m$^3$ yr$^-$$^1$)',xlabel='Time, years')
    ylim=ax.get_ylim()
    ax.set(ylim=[0,ylim[1]+0.03*ylim[1]])
    ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False);
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    #ax.text((2020+1960)/2,12,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_HarvestVolume_Scn' + str(iScn+1),'png',900)

#%% Harvest rates (Actual vs. Forest Retention Scenario)
def PlotHarvestMultipleScenarios(meta,mos,tv,AEF,iT,iPS,iSS,iScn1,iScn2,flg_hbs):

    # Import FAIB public summary
    H_FLNR=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')

    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));
    #ax.add_patch(Rectangle([2005,0],10,200,fc=[0.8,0.8,0.8],ec="none"))
    #ax.add_patch(Rectangle([2015,0],100,200,fc=[0.9,0.9,0.9],ec="none"))

    if flg_hbs==1:
        ax.plot(H_FLNR['Year'],H_FLNR['Total_harvest_millions_m3'],'k-',lw=1.25,label='Tot. merch. vol. (FLNR)')
        #ax.plot(H_HBS['Year'],H_HBS['V All m3']/1e6,'-cd')

    y=AEF*gu.movingave(mos['Scenarios'][iScn1]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][iT,iPS,iSS],10,'historical')/1e6
    ax.plot(tv[iT],y,'-',color=gp['cl1'],label='Baseline scenario')

    y=AEF*gu.movingave(mos['Scenarios'][iScn2]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][iT,iPS,iSS],10,'historical')/1e6
    ax.plot(tv[iT],y,'--',color=gp['cl2'],label='Actual scenario')

    ax.set(position=[0.09,0.12,0.88,0.84],xlim=[1860,2100],xticks=np.arange(1800,2120,20),ylim=[0,140], \
           ylabel='Harvest volume (million m$^3$ yr$^-$$^1$)',xlabel='Time, years')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    ax.legend(loc='upper left',frameon=False,facecolor='w',edgecolor='w');
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_HarvestVolumeMerchCompareScenarios_Scn' + str(iScn1+1) + '_vs_' + str(iScn2+1),'png',900)

#%% Volume per hectare of harvest
def PlotVolumePerHectare(meta,mos,tv,iScn,iT,iPS,iSS,AEF):

    # Observations are relatively stable at 350 m3/ha over 1990-2018 (https://cfs.nrcan.gc.ca/statsprofile/forest/bc)

    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));

    # Modelled
    ivlT=1;
    A=cbu.SummarizeAreaAffected(meta,mos,tv,iScn,iPS,iSS,AEF,ivlT)
    A_Harvest=A['Management'][0]['Data']
    V_Harvest=AEF*mos['Scenarios'][iScn]['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][:,iPS,iSS]
    rat=np.maximum(0,V_Harvest/A_Harvest)
    rat_ma=gu.movingave(rat,5,'historical')
    ax.plot(tv[iT],rat_ma[iT],'-go',lw=1,mfc=[0.6,0.96,0],mec=[0.6,0.96,0],color=[0.6,0.96,0],ms=3,label='FCS prediction')

    # Observations
    d={}
    d['Year']=np.arange(1850,2022,1)
    d['V']=np.zeros(d['Year'].size)
    d['A']=np.zeros(d['Year'].size)

    H_FLNR=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.xlsx')
    ind=np.where( (d['Year']>=H_FLNR['Year'][0]) & (d['Year']<=H_FLNR['Year'][-1]) )[0]
    d['V'][ind]=H_FLNR['Total_harvest_millions_m3']*1e6

    H_HBS=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS_AnnualSummary.pkl')
    ind=np.where( (d['Year']>=H_HBS['Year'][0]) & (d['Year']<=H_HBS['Year'][-1]) )[0]
    d['V'][ind]=np.maximum(d['V'][ind],H_HBS['V All (Mm3/yr)']*1e6)

    A_CCB=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromConCutblocksDB.pkl')
    ind=np.where( (d['Year']>=A_CCB['Year'][0]) & (d['Year']<=A_CCB['Year'][-1]) )[0]
    d['A'][ind]=A_CCB['Area Harvested']

    A_RES=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromRESULTS.pkl')
    ind=np.where( (d['Year']>=A_RES['Year'][0]) & (d['Year']<=A_RES['Year'][-1]) )[0]
    d['A'][ind]=np.maximum(d['A'][ind],A_RES['Area Harvested'])

    iT3=np.where( (d['Year']>=tv[iT[0]]) & (d['Year']<=tv[iT[-1]]) )[0]

    plt.plot(d['Year'][iT3],d['V'][iT3]/d['A'][iT3],'-bs',color=[0.27,0.49,0.79],mfc=[0.27,0.49,0.79],mec=[0.27,0.49,0.79],ms=3,lw=1,label='HBS + Harvest Area Analysis')

    ax.set(position=[0.1,0.12,0.88,0.84],xticks=np.arange(tv[iT[0]],2250,5),xlabel='Time, years', \
           yticks=np.arange(0,2000,100),ylabel='Harvest volume (m$^3$ ha$^-$$^1$)',
           xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[200,1220])
    ax.legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w');
    ax.grid(True)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_HarvestVolumePerHectare_Scn' + str(iScn+1),'png',900)

#%% Plot net sector GHG balance
def PlotGHGBalance(meta,mos,tv,AEF,iScn,iT,iPS,iSS):
    #iScn=0

    v='E_CO2e_AGHGB_WSub'
    v2='E_CO2e_AGHGB_WOSub'

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,8)); wd=1;
    ax.add_patch(Rectangle([1920,-200],100,550,fc=[0.94,0.94,0.94],ec="none"))
    ax.plot(tv[iT],np.zeros(iT.size),'k-',color=gp['cla'])

    y_hwp=mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
    ax.bar(tv[iT],y_hwp,wd,label='Product emissions',facecolor=[0.75,0.85,1])

    y_op=mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Operations']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
    ax.bar(tv[iT],y_op,wd,bottom=y_hwp,label='Fossil fuel emissions',facecolor=[0,0,0.7])

    y_ob=mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
    ax.bar(tv[iT],y_ob,wd,bottom=y_hwp+y_op,label='Open burning emissions',facecolor=[1,0.1,0.1])

    y_wf=mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
    ax.bar(tv[iT],y_wf,wd,bottom=y_hwp+y_op+y_ob,label='Wildfire emissions',facecolor=[1,0.7,0.7])

    # Negative harvest removals
    #ax.bar(tv[iT],-1*(mos['Scenarios'][iScn]['Sum']['C_ToMill']['Ensemble Mean'][iT,iPS,iSS])*3.667/1e6*AEF,wd,facecolor=[0.65,.85,0.2],label='Removals (ecosystem to HWP)')
    ax.bar(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_SUB_Tot']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,wd,facecolor=[0.65,.85,0.2],label='Substitution effects')

    ax.plot(tv[iT],mos['Scenarios'][iScn]['Sum'][v]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'k-',mfc='w',mew=0.5,ms=2,label='GHG balance (with subs.)')
    ax.plot(tv[iT],mos['Scenarios'][iScn]['Sum'][v2]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'k--',color=[0.45,0.45,0.45],mfc='w',mew=0.5,ms=2,label='GHG balance (w/o subs.)')

    ax.plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_NEE']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'k-.',color=[0.4,0.6,0.15],mfc='w',ms=2,mew=0.5,label='Net ecosystem exchange')

    ax.set(position=[0.075,0.085,0.9,0.86],xticks=np.arange(tv[iT[0]],2250,20),ylim=[-160,160],yticks=np.arange(-200,200,25),
           ylabel='GHG fluxes (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax.legend(loc='lower right',fontsize=gp['fs3'],frameon=False,facecolor='w',edgecolor='w');
    #plt.grid()
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    ax.text(1865,92.5,'Emissions',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7])
    ax.text(1865,-92.5,'Removals',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7])
    ax.text(1970,100,'Modern period',fontsize=9,style='normal',weight='bold',color=[0.6,0.6,0.6],ha='center')
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_GHGBalance_Scn' + str(iScn+1),'png',900)

#%% Plot scenario comparisons
def PlotGHGBalanceAndBenefitWithCumulative(meta,mos,tv,AEF,cmp,iT,iPS,iSS):

    iB=mos['Delta'][cmp]['iB']
    iP=mos['Delta'][cmp]['iP']

    v1=['E_CO2e_AGHGB_WSub','E_CO2e_AGHGB_WOSub']
    v2=['E_CO2e_AGHGB_WSub_cumu_from_tref','E_CO2e_AGHGB_WOSub_cumu_from_tref']
    #v2='E_CO2e_AGHGB_WSub_cumu'
    for iV in range(len(v1)):
        plt.close('all')
        fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(18.8,10))

        ax[0,0].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])
        ylo1=mos['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi1=mos['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        ylo2=mos['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi2=mos['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[0,0].fill_between(tv[iT],ylo1,yhi1,color=gp['cl1'],alpha=gp['Alpha2'],lw=0)
        ax[0,0].fill_between(tv[iT],ylo2,yhi2,color=gp['cl2'],alpha=gp['Alpha2'],lw=0)
        ax[0,0].plot(tv[iT],mos['Scenarios'][iB]['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'-',color=gp['cl1'],label='Baseline')
        ax[0,0].plot(tv[iT],mos['Scenarios'][iP]['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'--',color=gp['cl2'],label='Actual')

        leg=ax[0,0].legend(loc='upper left',bbox_to_anchor=(0.1,0.45,0.5,0.5),frameon=False,facecolor=None,edgecolor='w');
        for text in leg.get_texts():
            plt.setp(text,color=gp['cla']);
        #ax[0,0].text(tv[iT[0]]+7,80,'Source',fontsize=7,style='italic',weight='bold',color=clt,va='center')
        #ax[0,0].text(tv[iT[0]]+7,-80,'Sink',fontsize=7,style='italic',weight='bold',color=clt,va='center')
        ax[0,0].set(position=[0.07,0.57,0.4,0.4],xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50), \
          ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
        ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

        ax[0,1].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])
        ylo1=mos['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS]/1e9*AEF
        yhi1=mos['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS]/1e9*AEF
        ylo2=mos['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS]/1e9*AEF
        yhi2=mos['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS]/1e9*AEF
        cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[0,1].fill_between(tv[iT],ylo1,yhi1,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
        ax[0,1].fill_between(tv[iT],ylo2,yhi2,color=gp['cl2'],alpha=gp['Alpha2'],lw=0)
        ax[0,1].plot(tv[iT],mos['Scenarios'][iB]['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e9*AEF,'-',color=gp['cl1'],label='Baseline')
        ax[0,1].plot(tv[iT],mos['Scenarios'][iP]['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e9*AEF,'--',color=gp['cl2'],label='Actual')
        for text in leg.get_texts():
            plt.setp(text, color=gp['cla'])
        ax[0,1].set(position=[0.57,0.57,0.4,0.4],xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50), \
          ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],ylabel='Cumulative GHG balance (GtCO$_2$e)',xlabel='Time, years')
        ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

        ylo=mos['Delta'][cmp]['ByStrata']['Sum'][v1[iV]]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi=mos['Delta'][cmp]['ByStrata']['Sum'][v1[iV]]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        ymu=mos['Delta'][cmp]['ByStrata']['Sum'][v1[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
        cs=np.column_stack((ylo,ymu,yhi))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[1,0].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])
        ax[1,0].fill_between(tv[iT],ylo,yhi,color=gp['cl3'],alpha=gp['Alpha2'],lw=0,label='95% C.I.')

        ylo2=mos['Delta'][cmp]['ByStrata']['Sum'][v1[iV]]['Ensemble P250'][iT,iPS,iSS]/1e6*AEF
        yhi2=mos['Delta'][cmp]['ByStrata']['Sum'][v1[iV]]['Ensemble P750'][iT,iPS,iSS]/1e6*AEF
        ax[1,0].fill_between(tv[iT],ylo2,yhi2,color=gp['cl3'],alpha=gp['Alpha2'],lw=0,label='50% C.I.')

        ax[1,0].plot(tv[iT],ymu,'-',color=gp['cl3'],label='Mean')
        ax[1,0].set(position=[0.07,0.07,0.4,0.42],xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50),
          ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)], \
          ylabel='$\Delta$GHG (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
        ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
        #ax[1,0].text(tv[iT[0]]+7,0.6*np.abs(np.minimum(ymn,ymx)),'Net emission',fontsize=7,style='italic',weight='bold',color=clt,va='center')
        #ax[1,0].text(tv[iT[0]]+7,-0.6*np.abs(np.minimum(ymn,ymx)),'Net removal',fontsize=7,style='italic',weight='bold',color=clt,va='center')
        leg=ax[1,0].legend(loc='upper left',bbox_to_anchor=(0.65,0.47,0.5,0.5),frameon=False,facecolor=None,edgecolor='w')
        for text in leg.get_texts():
            plt.setp(text,color=gp['cla']);
        ylo=mos['Delta'][cmp]['ByStrata']['Sum'][v2[iV]]['Ensemble P025'][iT,iPS,iSS]/1e9*AEF
        yhi=mos['Delta'][cmp]['ByStrata']['Sum'][v2[iV]]['Ensemble P975'][iT,iPS,iSS]/1e9*AEF
        ymu=mos['Delta'][cmp]['ByStrata']['Sum'][v2[iV]]['Ensemble Mean'][iT,iPS,iSS]/1e9*AEF
        cs=np.column_stack((ylo,ymu,yhi))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[1,1].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])
        ax[1,1].fill_between(tv[iT],ylo,yhi,color=gp['cl3'],alpha=gp['Alpha1'],lw=0)

        ylo2=mos['Delta'][cmp]['ByStrata']['Sum'][v2[iV]]['Ensemble P250'][iT,iPS,iSS]/1e9*AEF
        yhi2=mos['Delta'][cmp]['ByStrata']['Sum'][v2[iV]]['Ensemble P750'][iT,iPS,iSS]/1e9*AEF
        ax[1,1].fill_between(tv[iT],ylo2,yhi2,color=gp['cl3'],alpha=gp['Alpha2'],lw=0,label='50% C.I.')

        ax[1,1].plot(tv[iT],ymu,'-',color=gp['cl3'],label='Actual minus baseline')
        ax[1,1].set(position=[0.57,0.07,0.4,0.42],xlim=[tv[iT[0]],tv[iT[-1]]],xticks=np.arange(tv[iT[0]],2200,50),ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)], \
          ylabel='Cumulative $\Delta$GHG (GtCO$_2$e)',xlabel='Time, years')
        ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])

        gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')

        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_DeltaGHG_And_GHGBalance_' + cmp + '_' + v1[iV],'png',900)


#%% GHG benefit (with Subs)
def PlotGHGBenefit(meta,mos,tv,AEF,cmp,iT,iPS,iSS):

    vL=['E_CO2e_AGHGB_WSub','E_CO2e_AGHGB_WOSub']
    lab=['WithSub','WithoutSub']

    iB=mos['Delta'][cmp]['iB']
    iP=mos['Delta'][cmp]['iP']

    for iv in range(len(vL)):
        v=vL[iv]
        plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.8,6));
        ax[0].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])
        ylo1=mos['Scenarios'][iB]['Sum'][v]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi1=mos['Scenarios'][iB]['Sum'][v]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        ylo2=mos['Scenarios'][iP]['Sum'][v]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi2=mos['Scenarios'][iP]['Sum'][v]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        cs=np.column_stack((ylo1,ylo2,yhi1,yhi2))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[0].fill_between(tv[iT],ylo1,yhi1,color=gp['cl1'],alpha=gp['Alpha2'],lw=0)
        ax[0].fill_between(tv[iT],ylo2,yhi2,color=gp['cl2'],alpha=gp['Alpha2'],lw=0)
        ax[0].plot(tv[iT],mos['Scenarios'][iB]['Sum'][v]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'-',color=gp['cl1'],label='Baseline')
        ax[0].plot(tv[iT],mos['Scenarios'][iP]['Sum'][v]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'--',color=gp['cl2'],label='Actual')
        leg=ax[0].legend(loc='upper left',bbox_to_anchor=(0.08,0.48,0.5,0.5),frameon=False,facecolor=None,edgecolor='w')
        #for text in leg.get_texts():
        #    plt.setp(text, color=cla)
        ax[0].text(tv[iT[0]]+7,0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Source',fontsize=7,style='italic',weight='bold',color=gp['clt'],va='center')
        ax[0].text(tv[iT[0]]+7,-0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Sink',fontsize=7,style='italic',weight='bold',color=gp['clt'],va='center')
        ax[0].set(position=[0.07,0.11,0.43,0.83],xticks=np.arange(tv[iT[0]],2200,50), \
          yticks=np.arange(-200,300,40),ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[-160,220],xlim=[tv[iT[0]],tv[iT[-1]]])
        ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=gp['tickl'])

        ylo=mos['Delta'][cmp]['ByStrata']['Sum'][v]['Ensemble P025'][iT,iPS,iSS]/1e6*AEF
        yhi=mos['Delta'][cmp]['ByStrata']['Sum'][v]['Ensemble P975'][iT,iPS,iSS]/1e6*AEF
        ymu=mos['Delta'][cmp]['ByStrata']['Sum'][v]['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF
        cs=np.column_stack((ylo,ymu,yhi))
        ymx=np.max(np.max(cs,axis=1)); ymn=np.min(np.min(cs,axis=1))
        ax[1].plot(tv[iT],np.zeros(iT.size),color=gp['cla'])

        ax[1].fill_between(tv[iT],ylo,yhi,color=gp['cl3'],alpha=gp['Alpha1'],lw=0,label='95% C.I.')
        ylo2=mos['Delta'][cmp]['ByStrata']['Sum'][v]['Ensemble P250'][iT,iPS,iSS]/1e6*AEF
        yhi2=mos['Delta'][cmp]['ByStrata']['Sum'][v]['Ensemble P750'][iT,iPS,iSS]/1e6*AEF
        ax[1].fill_between(tv[iT],ylo2,yhi2,color=gp['cl3'],alpha=gp['Alpha2'],lw=0,label='50% C.I.')
        ax[1].plot(tv[iT],ymu,'-',color=gp['cl3'],label='Mean')
        ax[1].set(position=[0.57,0.11,0.43,0.83],xticks=np.arange(tv[iT[0]],2200,50),yticks=np.arange(-200,300,20),
          ylabel='$\Delta$GHG (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',ylim=[ymn-0.01*np.abs(ymn),ymx+0.01*np.abs(ymx)],xlim=[tv[iT[0]],tv[iT[-1]]])
        ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=gp['tickl'])
        ax[1].text(tv[iT[0]]+7,0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Net emission',fontsize=7,style='italic',weight='bold',color=gp['clt'],va='center');
        ax[1].text(tv[iT[0]]+7,-0.65*np.minimum(np.abs(ymx),np.abs(ymn)),'Net removal',fontsize=7,style='italic',weight='bold',color=gp['clt'],va='center');
        leg=ax[1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')

        gu.axletters(ax,plt,0.03,0.92,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold');
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS];
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS];
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_DeltaGHGOnly_' + lab[iv] + '_' + cmp,'png',900);

    return

#%%
def PlotComparisonWithPIR(meta,mos,tv,AEF,iScn,iT,iPS,iSS):

    dPIR=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\PIR\PIR Reformatted.xlsx')

    iT=np.where( (tv>=1900) & (tv<=2022) )[0]

    plt.close('all')
    fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,8));
    ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-')

    ax[0,0].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_NEE']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'ob-',mfc='w',mew=0.5,ms=2,label='FCAST')
    ax[0,0].plot(dPIR['Year'],dPIR['Forest Growth Minus Decay']/1e3,'rs-',mfc=[1,1,1],label='PIR')
    ax[0,0].set(position=[0.08,0.59,0.41,0.4],xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-120,300,20),
           ylabel='Growth - Decay (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[-100,60])
    ax[0,0].legend(loc='lower left',frameon=False,facecolor='w',edgecolor='w')
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

    ax[0,1].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'ob-',mfc='w',mew=0.5,ms=2,label='Net ecosystem exchange')
    ax[0,1].plot(dPIR['Year'],dPIR['Wildfires']/1e3,'rs-',mfc=[1,1,1],label='PIR')
    ax[0,1].set(position=[0.57,0.59,0.41,0.4],xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(0,300,25),
           ylabel='Wildfire (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[0,220])
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

    #ax[1,0].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'ob-',mfc='w',mew=0.5,ms=2,label='Net ecosystem exchange')
    #ax[1,0].plot(dPIR['Year'],dPIR['Slash Pile Burning']/1e3,'rs-',mfc=[1,1,1],label='PIR')
    ax[1,0].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'ob-',mfc='w',mew=0.5,ms=2,label='HWP decay (FCS)')
    ax[1,0].plot(tv[iT],(mos['Scenarios'][iScn]['Sum']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][iT,iPS,iSS]+mos['Scenarios'][iScn]['Sum']['E_CO2e_ESC_Bioenergy']['Ensemble Mean'][iT,iPS,iSS])/1e6*AEF,'oc-',mfc='w',mew=0.5,ms=2,label='HWP decay + bioenergy (FCS)')
    #ax[1,0].plot(tv[iT],gu.movingave(mos['Scenarios'][iScn]['Sum']['C_ToMill']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF*3.667,10,'center'),'og-',mfc='w',mew=0.5,ms=2,label='Harvest removals (FCAST)')
    ax[1,0].plot(dPIR['Year'],dPIR['Decomposition of Harvested Wood Products']/1e3,'rs-',mfc=[1,1,1],label='HWP decay (PIR)')
    ax[1,0].set(position=[0.08,0.09,0.41,0.4],xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-120,300,20),
           ylabel='HWP (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[0,100])
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
    ax[1,0].legend(loc='center left',frameon=False,facecolor='w',edgecolor='w')

    #ax[1,0].plot(Harvest['Year'],Harvest['Total_harvest_millions_m3']*0.46*0.5*3.667,'k-',lw=0.75,label='FLNRORD estimate')
    #ax.plot(H_HBS['Year'],H_HBS['V All m3']/1e6,'-cd')

    ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.5)
    ax[1,1].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'ob-',mfc='w',mew=0.5,ms=2,label='With subs. (FCS)')
    ax[1,1].plot(tv[iT],mos['Scenarios'][iScn]['Sum']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF,'sb-',color='c',mec='c',mfc='w',mew=0.5,ms=2,label='W/O subs. (FCS)')
    ax[1,1].plot(dPIR['Year'],dPIR['Forest Management']/1e3,'rs-',lw=0.25,mfc='w',label='Forest Management (PIR)')
    ax[1,1].set(position=[0.57,0.09,0.41,0.4],xticks=np.arange(tv[iT[0]],2250,20),yticks=np.arange(-150,300,50),
           ylabel='GHG balance (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5],ylim=[-150,250])
    ax[1,1].legend(loc='upper center',frameon=False,facecolor='w',edgecolor='w')
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
    #plt.grid()
    gu.axletters(ax,plt,0.03,0.885,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\ComparisonWithPIR','png',900)

#%%
def PlotAreaDisturbed(meta,mos,tv,AEF,ivlT,iScn,iT,iPS,iSS):

    A=cbu.SummarizeAreaAffected(meta,mos,tv,iScn,iPS,iSS,AEF,ivlT)

    yr_start=tv[iT[0]]

    plt.close('all');
    fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(15,10));
    #ax[0].plot([2020,2020],[0,3400],'k--')
    #ax[0].fill_between([1920,2020],[0,0],[11000,11000],color=[0.9,0.9,0.9],lw=0)
    pl_d=[None]*len(A['Nat Dist']); nams_d=[None]*len(A['Nat Dist']);
    for i in range(len(A['Nat Dist'])):
        bottom=0;
        if i!=0:
            for j in range(i):
                bottom=bottom+A['Nat Dist'][j]['Data']/1000
        pl_d[i]=ax[0].bar(A['tv'],A['Nat Dist'][i]['Data']/1000,ivlT,color=A['Nat Dist'][i]['Color'],bottom=bottom)
        nams_d[i]=A['Nat Dist'][i]['Name']
    ax[0].legend(pl_d,nams_d,loc='upper left',bbox_to_anchor=(0.05,0.98),labelspacing=0.12,facecolor=[1,1,1],frameon=False);
    ax[0].set(position=[0.06,0.56,0.92,0.43],yscale='linear', \
      xticks=np.arange(np.min(A['tv']),np.max(A['tv'])+1,20),ylabel='Area affected (ha x 1000)',xlim=[yr_start,np.max(A['tv'])]);
    #ax[0].text((2020+1920)/2,1050,'Observation Period',ha='center',fontsize=6,style='normal',weight='normal',color=[0,0,0])

    #ax[1].fill_between([1920,2020],[0,0],[1100,1100],color=[0.9,0.9,0.9],lw=0)
    pl_m=[None]*len(A['Management']); nams_m=[None]*len(A['Management']);
    for i in range(len(A['Management'])):
        bottom=0;
        if i!=0:
            for j in range(i):
                bottom=bottom+A['Management'][j]['Data']/1000
        pl_m[i]=ax[1].bar(A['tv'],A['Management'][i]['Data']/1000,ivlT,color=A['Management'][i]['Color'],bottom=bottom)
        nams_m[i]=A['Management'][i]['Name']
    #ax[1].plot([2020,2020],[0,5500],'k--')
    ax[1].legend(pl_m,nams_m,loc='upper left',bbox_to_anchor=(0.05,0.98),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
    ax[1].set(position=[0.06,0.08,0.92,0.43],yscale='linear',
      xticks=np.arange(np.min(A['tv']),np.max(A['tv'])+1,20),xlabel='Time, years',ylabel='Area affected (ha x 1000)',xlim=[yr_start,np.max(A['tv'])]);
    ax[1].tick_params(length=gp['tickl'])
    gu.axletters(ax,plt,0.01,0.91)

    nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
    nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_AreaDisturbedAndManaged_Scn' + str(iScn+1),'png',900)

#%% Mortality summary
def MortalitySummary(meta,mos,C_M_ByAgent,tv,iScn,iT,ivlT,iPS,iSS):

    y=[None]*7; c=-1
    c=c+1; y[c]={}; y[c]['Name']='Competition'; y[c]['Color']=[0.9,0.85,1]; y[c]['Data']=mos['Scenarios'][iScn]['Mean']['C_M_Reg_Tot']['Ensemble Mean'][iT,iPS,iSS]
    c=c+1; y[c]={}; y[c]['Name']='Weather & disease'; y[c]['Color']=[1,0.95,0.8]; y[c]['Data']=C_M_ByAgent[iScn]['Mechanical'][iT]
    c=c+1; y[c]={}; y[c]['Name']='Wildfire'; y[c]['Color']=[0.75,0,0]; y[c]['Data']=C_M_ByAgent[iScn]['Wildfire'][iT]
    c=c+1; y[c]={}; y[c]['Name']='Beetles'; y[c]['Color']=[0.5,0.8,0]; y[c]['Data']=C_M_ByAgent[iScn]['IBM'][iT]+C_M_ByAgent[iScn]['Beetles'][iT]+C_M_ByAgent[iScn]['IBB'][iT]
    c=c+1; y[c]={}; y[c]['Name']='Weevils'; y[c]['Color']=[0.5,0.75,0.75]; y[c]['Data']=C_M_ByAgent[iScn]['Weevils'][iT]
    c=c+1; y[c]={}; y[c]['Name']='Defoliators'; y[c]['Color']=[0.65,1,0]; y[c]['Data']=C_M_ByAgent[iScn]['Defoliators'][iT]
    c=c+1; y[c]={}; y[c]['Name']='Harvest'; y[c]['Color']=[0.75,0.9,1]; y[c]['Data']=C_M_ByAgent[iScn]['Harvest'][iT]+C_M_ByAgent[iScn]['Harvest Salvage'][iT]

    # Convert to x-year intervals
    y[0]['tv']=gu.BlockMean(tv[iT],ivlT)
    for i in range(len(y)):
        y[i]['Data']=gu.BlockMean(y[i]['Data'],ivlT)

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(16.5,8)); yr_start=tv[iT[0]];
    pl_d=[None]*len(y); nams_d=[None]*len(y);
    for i in range(len(y)):
        bottom=0;
        if i!=0:
            for j in range(i):
                bottom=bottom+y[j]['Data']
        pl_d[i]=ax.bar(y[0]['tv'],y[i]['Data'],ivlT,color=y[i]['Color'],bottom=bottom)
        nams_d[i]=y[i]['Name']
    ax.legend(pl_d,nams_d,loc='upper left',bbox_to_anchor=(0.03,0.96),labelspacing=0.12,facecolor=[1,1,1],frameon=False)
    ax.set(position=[0.07,0.08,0.92,0.9],ylim=[0,2.8],yticks=np.arange(0,3,0.2),ylabel='Mortality (MgC/ha/yr)', \
           xlim=[yr_start,np.max(y[0]['tv'])],xticks=np.arange(np.min(y[0]['tv']),np.max(y[0]['tv'])+1,20),xlabel='Time, years');
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    if flag_printfigs=='On':
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_MortalityFlux_Scn' + str(iScn+1),'png',900)

#%% Plot time series

def PlotAllVariableTimeSeries(meta,mos,tv,iScn,iT,iPS,iSS):
    for k in mos['Scenarios'][iScn]['Sum']:
        plt.close('all')
        fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7))
        mxn=np.zeros(2)
        for iScn in range(meta['Project']['N Scenario']):
            #mu=mos[iScn]['v1'][k]['Mean']
            #mu=mos[iScn]['v1']['Mean'][k]['Ensemble Mean'];
            mu=mos['Scenarios'][iScn]['Mean'][k]['Ensemble Mean'][:,iPS,iSS]
            mxn[0]=np.minimum(mxn[0],np.min(mu[iT])); mxn[1]=np.maximum(mxn[1],np.max(mu[iT]));
            ax.plot(tv[iT],mu[iT],'-',lw=0.75)
        mu=np.mean(mxn)
        mxn[0]=mxn[0]-0.1*mu
        mxn[1]=mxn[1]+0.1*mu
        ax.set(position=[0.13,0.13,0.8,0.82],xlim=[np.min(tv[iT]),np.max(tv[iT])],ylim=mxn,
                   ylabel=k,xlabel='Time, years',aspect='auto');
        fig.patch.set_facecolor('w')
        nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
        nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
        gu.PrintFig(meta['Paths']['Figures'] + '\\TS\\' + nam_ps + '_' + nam_ss + '_ts_' + k + '_mean_' + str(iScn+1),'png',150)

#%% Plot other fluxes

def PlotCarbonFluxTS(meta,mos,tv,iT,iB,iP,iPS,iSS):

    plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,10))

    ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble Mean'][iT,iPS,iSS]
    ax[0,0].fill_between(tv[iT],lo,hi,color=gp['cl2'],alpha=gp['Alpha1'],lw=0)
    ax[0,0].plot(tv[iT],mu,'-',color=gp['cl2'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_NEE']['Ensemble Mean'][iT,iPS,iSS]
    ax[0,0].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[0,0].plot(tv[iT],mu,'--',color=gp['cl1'],label='Baseline scenario'+0.25)
    ax[0,0].set(position=[0.07,0.57,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='NEE (tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[0,0].legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

    ax[0,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][iT,iPS,iSS]
    ax[0,1].fill_between(tv[iT],lo,hi,color=gp['cl2'],alpha=gp['Alpha1'],lw=0)
    ax[0,1].plot(tv[iT],mu,'-',color=gp['cl2'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_Wildfire']['Ensemble Mean'][iT,iPS,iSS]
    ax[0,1].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[0,1].plot(tv[iT],mu,'--',color=gp['cl1'],label='Baseline scenario'+0.25)
    ax[0,1].set(position=[0.57,0.57,0.42,0.42],yscale='linear', \
      xticks=np.arange(0,2220,25),xlabel='Time, years',ylabel='Wildfire emissions (tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

    ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][iT,iPS,iSS]
    ax[1,0].fill_between(tv[iT],lo,hi,color=gp['cl2'],alpha=gp['Alpha1'],lw=0)
    ax[1,0].plot(tv[iT],mu,'-',color=gp['cl2'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_OpenBurning']['Ensemble Mean'][iT,iPS,iSS]
    ax[1,0].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[1,0].plot(tv[iT],mu,'--',color=gp['cl1'],label='Baseline scenario'+0.25)
    ax[1,0].set(position=[0.07,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Open burning emissions (tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])

    ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][iT,iPS,iSS]
    ax[1,1].fill_between(tv[iT],lo,hi,color=gp['cl2'],alpha=gp['Alpha1'],lw=0)
    ax[1,1].plot(tv[iT],mu,'-',color=gp['cl2'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean']['E_CO2e_LULUCF_HWP']['Ensemble Mean'][iT,iPS,iSS]
    ax[1,1].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[1,1].plot(tv[iT],mu,'--',color=gp['cl1'],label='Baseline scenario'+0.25)
    ax[1,1].set(position=[0.57,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Product emissions (tCO$_2$e ha$^-$$^1$yr$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
    gu.axletters(ax,plt,0.03,0.91,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')

    nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
    nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_CarbonFluxes_Scns' + str(iB) + 'and' + str(iP),'png',900)

#%% Plot pools
def PlotCarbonPoolTS(meta,mos,tv,iT,iB,iP,iPS,iSS):

    ms=2; Alpha=0.16; lw=1; cl0=[0.27,0.44,0.79]; cl1=[0.4,0.8,0]; cl2=[0,0.9,0.9]
    plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(20,10))

    v='C_Biomass_Tot'
    #ax[0,0].plot(tv[it],np.zeros(it.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[0,0].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[0,0].plot(tv[iT],mu,'-',color=gp['cl1'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[0,0].fill_between(tv[iT],lo,hi,color=cl0,alpha=gp['Alpha1'],lw=0)
    ax[0,0].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario'+0.25)
    ax[0,0].set(position=[0.07,0.57,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Biomass (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[0,0].legend(loc='lower right',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both')
    ymin,ymax=ax[0,0].get_ylim(); ax[0,0].set_ylim(0.0,ymax)

    v='C_DeadWood_Tot'
    #ax[0,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[0,1].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[0,1].plot(tv[iT],mu,'-',color=gp['cl1'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[0,1].fill_between(tv[iT],lo,hi,color=cl0,alpha=gp['Alpha1'],lw=0)
    ax[0,1].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario'+0.25)
    ax[0,1].set(position=[0.57,0.57,0.42,0.42],yscale='linear',xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Dead wood (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both')
    ymin,ymax=ax[0,1].get_ylim(); ax[0,1].set_ylim(0.0,ymax)

    v='C_Litter_Tot'
    #ax[1,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[1,0].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[1,0].plot(tv[iT],mu,'-',color=gp['cl1'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[1,0].fill_between(tv[iT],lo,hi,color=cl0,alpha=gp['Alpha1'],lw=0)
    ax[1,0].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario'+0.25)
    ax[1,0].set(position=[0.07,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Litter (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both')
    ymin,ymax=ax[1,0].get_ylim(); ax[1,0].set_ylim(0.0,ymax)

    v='C_Soil_Tot'
    #ax[1,1].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)
    lo=mos['Scenarios'][iP]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iP]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iP]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[1,1].fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax[1,1].plot(tv[iT],mu,'-',color=gp['cl1'],label='Actual scenario'+0.25)
    lo=mos['Scenarios'][iB]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iB]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iB]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax[1,1].fill_between(tv[iT],lo,hi,color=cl0,alpha=gp['Alpha1'],lw=0)
    ax[1,1].plot(tv[iT],mu,'--',color=cl0,label='Baseline scenario'+0.25)
    ax[1,1].set(position=[0.57,0.07,0.42,0.42],xticks=np.arange(0,2220,25), \
      xlabel='Time, years',ylabel='Soil (MgC ha$^-$$^1$)',xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both')
    ymin,ymax=ax[1,1].get_ylim(); ax[1,1].set_ylim(0.0,ymax)

    gu.axletters(ax,plt,0.03,0.91)

    nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
    nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_CarbonPools_Scns' + str(iB) + 'and' + str(iP),'png',900)

#%% Plot net growth
def PlotNetGrowthTS(meta,mos,tv,C_M_ByAgent,iT,iScn,iPS,iSS):

    #ms=2; Alpha=0.16; lw=1; cl0=[0.27,0.44,0.79]; cl1=[0.4,0.8,0]; cl2=[0,0.9,0.9]
    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,10))
    #ax[0,0].plot(tv[iT],np.zeros(iT.size),'k-',lw=0.75)

    v='C_G_Net_Tot'
    lo=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu=mos['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    ax.fill_between(tv[iT],lo,hi,color=gp['cl1'],alpha=gp['Alpha1'],lw=0)
    ax.plot(tv[iT],mu,'-',color=gp['cl1'],label='Net growth (w/o instects)')

    v='C_M_Reg_Tot'
    lo1=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi1=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu1=mos['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]
    #ax.fill_between(tv[iT],lo1,hi1,color='r',alpha=gp['Alpha1'],lw=0)
    #ax.plot(tv[iT],mu1,'-',color='r',label='Mortality'+0.25)

    v='C_M_Dist'
    lo2=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P025'][iT,iPS,iSS]
    hi2=mos['Scenarios'][iScn]['Mean'][v]['Ensemble P975'][iT,iPS,iSS]
    mu2=mos['Scenarios'][iScn]['Mean'][v]['Ensemble Mean'][iT,iPS,iSS]

    # Mortality from beetles
    Mb=C_M_ByAgent[iScn]['Beetles'][iT]+C_M_ByAgent[iScn]['IBM'][iT]+C_M_ByAgent[iScn]['IBS'][iT]+C_M_ByAgent[iScn]['IBB'][iT]

    lo3=lo1+lo2
    hi3=hi1+hi2
    mu3=gu.movingave(mu1+mu2,10,'historical')
    #ax.fill_between(tv[iT],lo3,hi3,color='y',alpha=gp['Alpha1'],lw=0)
    #ax.plot(tv[iT],mu3,'-',color='y',label='Mortality'+0.25)

    mth='center'
    lo4=gu.movingave(lo-Mb,7,mth)
    hi4=gu.movingave(hi-Mb,7,mth)
    mu4=gu.movingave(mu-Mb,7,mth)
    ax.fill_between(tv[iT],lo4,hi4,color='y',alpha=gp['Alpha1'],lw=0)
    ax.plot(tv[iT],mu4,'-',color='y',label='Net growth (with insects)')

    #------------

    d=gu.ImportMat(r'G:\My Drive\Data\psp_sl_ts_bc.mat','qt')
    ax.plot(d['tv'].flatten(),np.nanmean(d['NEBP'],axis=1),'-bo')

    ax.set(position=[0.07,0.07,0.84,0.84],
      xticks=np.arange(0,2220,25), \
      xlabel='Time, years', \
      ylabel='Net growth (tC ha$^-$$^1$yr$^-$$^1$)',
      xlim=[tv[iT[0]]-0.5,tv[iT[-1]]+0.5])
    ax.legend(loc='lower left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.06,0.92)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

    #gu.axletters(ax,plt,0.03,0.91,FontColor=cla,LetterStyle='Caps',FontWeight='Bold')
    nam_ps=meta['Project']['Strata']['Project']['Unique CD'][iPS]
    nam_ss=meta['Project']['Strata']['Spatial']['Unique CD'][iSS]
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + nam_ps + '_' + nam_ss + '_GrowthNet_' + str(iScn+1),'png',900)


#%% Mortality frequency distribution

def GetMortalityFrequencyDistribution(meta,iScn):

    iEns=0

    M=[None]*meta['Project']['N Scenario']
    #for iScn in range(meta['Project']['N Scenario']):

    tv=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)
    M[iScn]={}
    M[iScn]['Ma']={}
    M[iScn]['Mr']={}
    for k in meta['LUT']['Dist'].keys():
        M[iScn]['Ma'][k]=np.zeros((tv.size,meta['Project']['N Stand']))
        M[iScn]['Mr'][k]=np.zeros((tv.size,meta['Project']['N Stand']))
    M[iScn]['Ma']['Reg']=np.zeros((tv.size,meta['Project']['N Stand']))
    M[iScn]['Mr']['Reg']=np.zeros((tv.size,meta['Project']['N Stand']))

    for iBat in range(meta['Project']['N Batch']):

        indBat=cbu.IndexToBatch(meta,iBat)

        d1=cbu.LoadSingleOutputFile(meta,iScn,iEns,iBat)
        dmec=gu.ipickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl')

        # Uncompress event chronology if it has been compressed
        if 'idx' in dmec:
            idx=dmec['idx']
            tmp=dmec.copy()
            for v in ['ID_Type','MortalityFactor','GrowthFactor','ID_GrowthCurve']:
                dmec[v]=np.zeros((meta['N Time'],meta['N Stand'],meta['Core']['Max Events Per Year']),dtype='int16')
                dmec[v][idx[0],idx[1],idx[2]]=tmp[v]
            del tmp

        for iStandInBat in range(len(indBat)):
            iStand=indBat[iStandInBat]
            for iYr in range(tv.size):
                it=np.where(meta['Year']==tv[iYr])[0]
                if it.size==0:
                    continue
                for iE in range(meta['Core']['Max Events Per Year']):
                    id=dmec['ID_Type'][iYr,iStandInBat,iE]
                    if id==0:
                        continue
                    nam=cbu.lut_n2s(meta['LUT']['Dist'],id)[0]
                    M[iScn]['Ma'][nam][iYr,iStand]=d1['C_M_Dist'][iYr,iStandInBat]
                    M[iScn]['Mr'][nam][iYr,iStand]=d1['C_M_Dist'][iYr,iStandInBat]/np.maximum(0.0001,d1['C_Biomass'][iYr-1,iStandInBat])
                    M[iScn]['Ma']['Reg'][iYr,iStand]=d1['C_M_Reg'][iYr,iStandInBat]
                    M[iScn]['Mr']['Reg'][iYr,iStand]=d1['C_M_Reg'][iYr,iStandInBat]/np.maximum(0.0001,d1['C_Biomass'][iYr-1,iStandInBat])
        del d1,dmec
        garc.collect()

    return M

##%%
#
#iScn=2
##M=cbu.GetMortalityFrequencyDistribution(meta)
#M=GetMortalityFrequencyDistribution(meta,iScn)
#
## Frequency distribution
#bw=1; bin=np.arange(0,100+bw,bw)
#M[iScn]['Percent Freq']={}
#for k in M[iScn]['Mr'].keys():
#    y=M[iScn]['Mr'][k].copy().flatten()*100
#    M[iScn]['Percent Freq'][k]=np.zeros(bin.size)
#    for iBin in range(bin.size):
#        ind=np.where(np.abs(y-bin[iBin])<=bw/2)[0]
#        M[iScn]['Percent Freq'][k][iBin]=ind.size/y.size*100
#
## Organize and categorize so easy to plot
#cat=[['Reg'],['Mechanical'],['Harvest'],['Wildfire'],['IBM'],['Beetles','IBB','IBD','IBS'],['IDW']]
#lab=['Competition','Mechanical damage','Harvesting','Wildfire','Mountain pine beetle','Beetles (other)','W. spruce budworm']
#cl=[[0.27,0.49,0.78],[1,1,0],[0,0.9,0.9],[0.5,0,0],[0,0.85,0],[0,0.5,0],[1,0.5,0]]
#dat=[None]*len(cat)
#for i in range(len(cat)):
#    dat[i]={}
#    dat[i]['Data']=np.zeros(bin.size)
#    for j in range(len(cat[i])):
#        dat[i]['Data']=dat[i]['Data']+M[iScn]['Percent Freq'][cat[i][j]]
#    if i==0:
#        dat[i]['Bottom']=np.zeros(bin.size)
#    else:
#        dat[i]['Bottom']=dat[i-1]['Bottom']+dat[i-1]['Data']
#
#plt.close('all')
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9)); br=[None]*len(dat)
#for i in range(len(dat)):
#    br[i]=plt.bar(bin,dat[i]['Data'],bottom=dat[i]['Bottom'],width=bw,color=cl[i])
#plt.legend((br[0][0],br[1][0],br[2][0],br[3][0],br[4][0],br[5][0],br[6][0]),
#           lab,frameon=False)
#ax.set(position=[0.08,0.12,0.86,0.82],ylabel=r'Frequency',xlabel='Severity (% mortality)',xlim=[-0.5,100.5],yscale='log')
##gu.PrintFig(meta['Paths']['Figures'] + '\\StandMortalitySpectrum','png',900)
#
#
#plt.close('all')
#plt.plot(tv,np.mean(M[1]['Ma']['Mechanical'],axis=1),'-')

#%% Plot map of investigation sites

def PlotSitesOfInterest(meta):

    # Load basemap
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

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
    #    x,y=nddat[iP][iD]['Geometry'].exterior.xy
    #    plt.plot(x,y,'r-')

    ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
    #plt.savefig(PathProject + '\\SparseGrid_Map.png',format='png',dpi=900)
    plt.close('all')

    return

#%% QA - soil organic carbon stocks by BGC zone

def QA_SOC_ByBGCZone(meta,mu_mod):

    # Import soils (see soil_Shawetal2018_01_process.py
    soils=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

    # Lookup table for BGC zone (based on raster map that was used to populate values in the soils DB)
    lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

    # Import best-available inventory
    ba=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\ba.pkl')

    # Calculate mean by BGC zone
    u=np.unique(soils['becz'])
    lab=np.array(['' for _ in range(u.size)],dtype=object)
    soc={}
    soc['SS']=np.zeros(u.size)
    soc['mu']=np.zeros(u.size)
    soc['se']=np.zeros(u.size)
    soc['min_mu']=np.zeros(u.size)
    soc['org_mu']=np.zeros(u.size)
    for i in range(u.size):
        ind=np.where( (soils['becz']==u[i]) & (soils['TOT_C_THA']>0) )[0]
        soc['SS'][i]=ind.size
        soc['mu'][i]=np.nanmean(soils['TOT_C_THA'][ind])
        soc['se'][i]=np.nanstd(soils['TOT_C_THA'][ind])/np.sqrt(ind.size)
        soc['min_mu'][i]=np.nanmean(soils['MIN_C_THA'][ind])
        soc['org_mu'][i]=np.nanmean(soils['ORG_C_THA'][ind])
        ind=np.where(lutBGC['VALUE']==u[i])[0]
        if ind.size>0:
            lab[i]=lutBGC['ZONE'][ind][0]

    soc['Model Soil C']=np.zeros(u.size)
    soc['Model SS']=np.zeros(u.size)
    soc['Model SI']=np.zeros(u.size)
    soc['Model SI SPL']=np.zeros(u.size)
    for i in range(u.size):
        id=meta['LUT']['VRI']['BEC_ZONE_CODE'][lab[i]][0]
        ind_mod=np.where( (ba['BEC_ZONE_CODE']==id) & (ba['SI']>5) )[0]
        soc['Model SS'][i]=ind_mod.size
        soc['Model Soil C'][i]=np.nanmean(mu_mod['C_Soil_Tot'][ind_mod])
        soc['Model SI'][i]=np.nanmean(ba['SI'][ind_mod])
        ind_mod=np.where( (ba['BEC_ZONE_CODE']==id) & (ba['SI SPL']>5) )[0]
        soc['Model SI SPL'][i]=np.nanmean(ba['SI SPL'][ind_mod])

    # Put in order
    ord=np.argsort(soc['mu'])
    for k in soc:
        soc[k]=np.flip(soc[k][ord])
    lab=np.flip(lab[ord])
    soc['Zone']=lab

    # Plot bar chart
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size),soc['min_mu'],facecolor=[0.45,0.3,0.3],label='Mineral horizons')
    #ax.bar(np.arange(u.size),soc['org_mu'],facecolor=[0.15,0.05,0.05],bottom=soc['min_mu'],label='Organic horizon')
    ax.errorbar(np.arange(u.size),soc['min_mu'],yerr=soc['se'],color='k',fmt='none',capsize=2)
    ax.plot(np.arange(u.size),soc['Model Soil C'],'co',lw=0.75,ms=6,mfc='w')
    for i in range(u.size):
        ax.text(i,10,str(soc['SS'][i].astype(int)),color='w',ha='center',fontsize=8)
        ax.text(i,30,str(soc['Model SS'][i].astype(int)),color='c',ha='center',fontsize=8)
    ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,375],xticks=np.arange(u.size),
           xticklabels=lab,ylabel='Soil organic carbon (MgC ha$^{-1}$ yr$^{-1}$)')
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    gu.PrintFig(meta['Paths']['Figures'] + '\\QA_SOC_ByBGCZone_BarChart','png',900)

    # Scatterplot
    x=soc['min_mu']
    y=soc['Model Soil C']
    rs,txt=gu.GetRegStats(x,y)

    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
    ax.plot([0,1000],[0,1000],'-k',lw=3,color=[0.8,0.8,0.8])
    ax.plot(x,y,'ko',mfc='k',mec='w',lw=0.5,ms=6)
    ax.plot(rs['xhat'],rs['yhat'],'r-',lw=1,label='Best fit')
    ax.text(350,30,txt,fontsize=10,color='r',ha='right')
    ax.text(300,300,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed SOC (MgC ha$^{-1}$)',ylabel='Predicted SOC (MgC ha$^{-1}$)',xlim=[0,375],ylim=[0,375])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    gu.PrintFig(meta['Paths']['Figures'] + '\\QA_SOC_ByBGCZone_Scatterplot','png',900)

    return soc

#%% Comparison of anthropogenic component with CFS

def CompareAnthropogenicComponentWithCFS(meta,mos,iT,iPS,iSS,AEF):

    tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)
    dfCFS=pd.read_excel(r'C:\Users\rhember\Documents\Data\CFS gaff.xlsx')
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,5.5)); wd=0.8; lw=0.5; cla=[0.8,0.8,0.8]
    ax.bar(1,np.mean(dfCFS['OLD']),label='National GHG Inventory (CFS)',facecolor=[0.5,0.7,1])
    #ax.bar(2,np.mean(dfCFS['NEW']),label='National GHG Inventory (CFS old)',facecolor=[0.2,0.4,1])
    ax.text(1,np.mean(dfCFS['OLD'])-48,'National\nGHG inventory\n(BC CFS)',fontsize=7,style='italic',weight='bold',color=[0.5,0.7,1],ha='center')
    #ax.text(2,np.mean(dfCFS['NEW'])-45,'National\nGHG inventory\n(BC CFS new)',fontsize=9,style='italic',weight='bold',color=[0.2,0.4,1],ha='center')
    iT=np.where( (tv>=1990) & (tv<=2020) )[0]
    y1=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WOSub']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF)
    ax.bar(2,y1,wd,label='BC FCS W/O\nsubstitution\neffects',facecolor=[1,0.7,0.5])
    y_lo=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WOSub']['Ensemble P025'][iT,iPS,iSS]/1e6*AEF)
    y_hi=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WOSub']['Ensemble P975'][iT,iPS,iSS]/1e6*AEF)
    ax.errorbar(2,y1,yerr=y_hi-y1,color=[0.8,0.5,0.2],ls='',lw=1.5,capsize=2)
    ax.text(2,y1+25,'BC FCS\nW/O\nsubstitution',fontsize=7,style='italic',weight='bold',color=[1,0.7,0.5],ha='center')
    y2=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][iT,iPS,iSS]/1e6*AEF)
    y_lo=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble P025'][iT,iPS,iSS]/1e6*AEF)
    y_hi=np.mean(mos['Delta']['TDAF']['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble P975'][iT,iPS,iSS]/1e6*AEF)
    ax.bar(3,y2,wd,label='BC FCS with\nsubstitution\neffects',facecolor=[1,0.4,0.2])
    ax.errorbar(3,y2,yerr=y2-y_lo,color=[0.4,0.2,0],ls='',lw=1.5,capsize=2)
    ax.text(3,y2+35,'BC FCS\nwith\nsubstitution',fontsize=7,style='italic',weight='bold',color=[1,0.4,0.2],ha='center')
    ax.plot([0,5],[0,0],'k-',color=cla,lw=0.25)
    ax.text(0.5,125,'Emissions',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7],va='center')
    ax.text(0.5,-125,'Removals',fontsize=9,style='italic',weight='bold',color=[0.7,0.7,0.7],va='center')
    ax.set(position=[0.16,.05,0.8,0.88],xticks=[0],yticks=np.arange(-175,300,25),ylabel='Flux (MtCO$_2$e yr$^-$$^1$)',xticklabels={''},xlabel='',ylim=[-140,160],xlim=[0.35,3.65])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    #ax.set_title('British Columbia, 1990-2020',fontsize=11,weight='bold',color=cla)
    gu.PrintFig(meta['Paths']['Figures'] + '\\DirectAnthroFlux_BC','png',900)

    #%%

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,5.5)); wd=0.8;

    ax.bar(1,0,label='National GHG Inventory (CFS old)',facecolor=[0.5,0.7,1])
    ax.text(1,0-7,'National\nGHG inventories\n(equals zero)',fontsize=7,style='italic',weight='bold',color=[0.4,0.6,1],ha='center')

    ax.bar(2,5.9,label='National GHG Inventory (CFS old)',facecolor=[1,0.4,0.2])
    ax.errorbar(2,5.9,yerr=4.1,color=[0.7,0.2,0],ls='',lw=1.5,capsize=2)
    ax.text(2,5.9+6,'Bookkeeping &\ndynamic global\nvegetation\n models',fontsize=7,style='italic',weight='bold',color=[1,0.4,0.2],ha='center')

    ax.set(position=[0.16,.05,0.8,0.88],xticks=[0],yticks=np.arange(-180,300,5),
           ylabel='Flux (GtCO$_2$e yr$^-$$^1$)',xticklabels={''},xlabel='',xlim=[0.35,2.65],ylim=[-15,25])
    #ax.set_title('Global direct anthropogenic GHG emissions\n 2010-2019 (Nabuurs et al. 2022)',fontsize=11,weight='bold',color=[0.7,0.7,0.7])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    ax.plot([0,5],[0,0],'k-',color=cla,lw=0.25)
    gu.PrintFig(meta['Paths']['Figures'] + '\\DirectAnthroFlux_Global','png',900)
    return

#%% Plot Summary of existing initiatives

def SummarizeExistingInitiatives(meta):

    d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Emission Reduction Projections\Summary of Reporting Initiatives.xlsx')

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8));
    plt.plot(d['Year'],np.zeros(d['Year'].size),'k-',lw=0.5)
    plt.plot(d['Year'],d['BC Total Emissions FLFL+HWP All Gases (NIR2022)'],'-ko',lw=1,color=[0.7,0.7,0.7],mec=[0.7,0.7,0.7],mfc=[0.7,0.7,0.7],label='Total emissions NIR 22')
    #plt.plot(d['Year'],d['Forest Management (PIR)'],'cs')
    plt.plot(d['Year'],d['PTT BAU Scenario (CFS)'],'-rs',label='PTT BAU (CFS)',lw=1)
    plt.plot(d['Year'],d['PTT RL Scenario (CFS)'],'-b^',label='PTT RL (CFS)',lw=1)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
    ax.legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
    ax.set(position=[0.1,0.1,0.82,0.82],xlim=[d['Year'][0],d['Year'][-1]],ylabel='Flux (MtCO$_2$e yr$^-$$^1$)',xlabel='Time, years')
    gu.PrintFig(meta['Paths']['Figures'] + '\\Summary of Existing Initiatives','png',900)

    return

#%% Look at historical wildfire from aspatial model

def LookAtWildfireRecord(meta,iScn):

    hw={}
    hw['tv']=np.arange(meta['Project']['Year Start'],meta['Project']['Year End']+1,1)
    hw['Po']=np.zeros((meta['Project']['N Time'],meta['Project']['N Stand']),dtype='int8')

    for iEns in range(meta['Project']['N Ensemble']):

        if 'Use Frozen Ensembles' not in meta['Project']:
            wf_sim=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\Ensembles\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')
        else:
            wf_sim=gu.ipickle(meta['Project']['Use Frozen Ensembles'] + '\\wf_sim_Scn' + cbu.FixFileNum(iScn) + '_Ens' + cbu.FixFileNum(iEns) + '.pkl')

        if 'idx' in wf_sim:
            idx=wf_sim['idx']
            hw['Po'][idx[0],idx[1]]=hw['Po'][idx[0],idx[1]]+wf_sim['Occurrence']

    hw['Po']=hw['Po']/meta['Project']['N Ensemble']

    return hw

