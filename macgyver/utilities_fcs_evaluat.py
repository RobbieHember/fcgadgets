
#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.patches import Rectangle
import gc as garc
import warnings
import time
import copy
import matplotlib.colors
import matplotlib.ticker as ticker
#from matplotlib import animation
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcexplore.psp.Processing import psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Age

def Plot_EvalAge_CNV(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CNV']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Age VRI t0']['N']>=10)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(215,30,txt,fontsize=10,color='k',ha='right')
    ax.text(230,230,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed age (years)',ylabel='Predicted age (years)',xlim=[0,250],ylim=[0,250])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Scatterplot_CNV','png',meta['Graphics']['gp']['save fig dpi'])

    # Plot bar chart

    # Put in order
    ord=np.argsort(d['Age VRI t0']['mu'])
    lab=np.flip(lab[ord])
    for v in d:
        for k in d[v].keys():
            d[v][k]=np.flip(d[v][k][ord])

    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Age VRI t0']['mu'],barw,facecolor=cl[0,:],label='Ground plots')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod A t0']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Barchart_CNV','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%%

def Plot_EvalBiomass_CNV(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CNV']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot L t0']['N']>=10)[0]
    for v in gplt['vaL']:
        for k in d[v].keys():
            d[v][k]=d[v][k][ind]
    u=u[ind]
    lab=lab[ind]

    # Scatterplot
    flg=0
    if flg==1:
        x=d['Ctot L t0']['mu']
        y=d['Mod C_Biomass_Tot t0']['mu']
        ikp=np.where( (np.isnan(x+y)==False) & (d['Ctot L t0']['N']>=30) )[0]
        rs,txt=gu.GetRegStats(x[ikp],y[ikp])

        fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
        ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
        #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
        for i in range(ikp.size):
            ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
        ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
        ax.text(200,20,txt,fontsize=10,color='k',ha='right')
        ax.text(190,190,'1:1',fontsize=8,ha='center')
        ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (MgC ha$^{-1}$)',ylabel='Predicted biomass (MgC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
        #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
        if meta['Graphics']['Print Figures']=='On':
            gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_CNV','png',meta['Graphics']['gp']['save fig dpi'])

    # Plot bar chart

    # Put in order
    ord=np.argsort(d['Ctot L t0']['mu'])
    lab=np.flip(lab[ord])
    for v in d:
        for k in d[v].keys():
            d[v][k]=np.flip(d[v][k][ord])

    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl[0,:],label='Ground plots')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_CNV','png',meta['Graphics']['gp']['save fig dpi'])
        pass

    return

#%% Biomass (YSM)

def Plot_EvalBiomass_YSM(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF YSM']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot L t0']['N']>=10)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(200,20,txt,fontsize=10,color='k',ha='right')
    ax.text(190,190,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (MgC ha$^{-1}$)',ylabel='Predicted biomass (MgC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    # Plot bar chart

    # Put in order
    ord=np.argsort(d['Ctot L t0']['mu'])
    lab=np.flip(lab[ord])
    for v in d:
        for k in d[v].keys():
            d[v][k]=np.flip(d[v][k][ord])

    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl[0,:],label='Ground plots')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,150])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Biomass (TIPSY stemwood)

def Plot_EvalStemwoodFromTIPSY_CNV(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CNV']==1) &
                         (gplt['Csw L t0']>=0) & (gplt['Csw L t0']<10000))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Csw L t0']['N']>=10)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(200,20,txt,fontsize=10,color='k',ha='right')
    ax.text(190,190,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (MgC ha$^{-1}$)',ylabel='Predicted biomass (MgC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_BiomassSW_ByBGCZone_Scatterplot','png',meta['Graphics']['gp']['save fig dpi'])

    # # Plot bar chart

    # # Put in order
    # ord=np.argsort(d['Ctot L t0']['mu'])
    # lab=np.flip(lab[ord])
    # for v in d:
    #     for k in d[v].keys():
    #         d[v][k]=np.flip(d[v][k][ord])

    # Area=np.zeros(lab.size)
    # for i in range(lab.size):
    #     ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
    #     Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

    # # Area weighting
    # for v in d:
    #     d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
    #     d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
    # lab=np.append(lab,'Weighted\naverage')
    # u=np.append(u,0.0)

    # cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    # barw=0.32
    # plt.close('all')
    # fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    # ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],barw,facecolor=cl[0,:],label='Ground plots')
    # ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
    # ax.errorbar(np.arange(u.size)-barw/2-0.01,d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    # ax.errorbar(np.arange(u.size)+barw/2+0.01,d['Mod C_Biomass_Tot t0']['mu'],yerr=d['Mod C_Biomass_Tot t0']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    # #for i in range(u.size):
    # #    ax.text(i,8,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=8)
    # ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,250])
    # plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    # ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    # gu.PrintFig(meta['Paths']['Figures'] + '\\QA_Biomass_ByBGCZone_BarChart','png',900)

    return

#%% Biomass components

def Plot_Eval_BiomassComponents_CN(meta,gplt):
    u=np.array([1])
    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)
    for i in range(u.size):
        for v in gplt['vaL']:
            ind=np.where( (gplt['PTF CN']==1) & (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassComponents_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Biomass dynammics summary (total CO2e)

def Plot_Eval_BiomassDynamicsSummaryTotCO2e_CN(meta,gplt):

    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)
    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) &
                         (gplt['Ctot G Tot']>0) & (gplt['Ctot G Tot']<30))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot G Tot']['N']>=3)[0]
    for v in gplt['vaL']:
        for k in d[v].keys():
            d[v][k]=d[v][k][ind]
    u=u[ind]
    lab=lab[ind]

    # Area weighting
    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]
    for v in d:
        d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
        d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
    lab=np.append(lab,'Weighted\naverage')
    u=np.append(u,0.0)

    # Area of forest
    #ind=np.where(zLCC1['Data']==lut_1ha['lcc1']['Forest Land'])
    #A_Tot=ind[0].size/1e6*3.67
    A_Tot=62000000/1e6*3.67

    cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
    lab=['Gross\ngrowth','Mortality','Net\ngrowth'] #,'Harvest\nmortality'

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]]); barw=0.38
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6))
    ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
    ax.bar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1]*A_Tot,barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1]*A_Tot,barw,facecolor=cl[1,:],label='Predictions (FCS)')
    ax.errorbar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1]*A_Tot,yerr=d['Ctot G Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1]*A_Tot,yerr=d['Mod C_G_Gross_Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1]*A_Tot,barw,facecolor=cl[0,:])
    ax.bar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1]*A_Tot,barw,facecolor=cl[1,:])
    ax.errorbar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1]*A_Tot,yerr=d['Ctot Mort+Lost']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1]*A_Tot,yerr=d['Mod C_M_Tot']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(3-barw/2-0.01,d['Ctot Net']['mu'][-1]*A_Tot,barw,facecolor=cl[0,:])
    ax.bar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1]*A_Tot,barw,facecolor=cl[1,:])
    ax.errorbar(3-barw/2-0.01,d['Ctot Net']['mu'][-1]*A_Tot,yerr=d['Ctot Net']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1]*A_Tot,yerr=d['Mod C_G_Net']['se'][-1]*A_Tot,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-500,600,100),ylabel='Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)',xlim=[0.5,3.5],ylim=[-500,500])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    ax.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsTotal_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Biomass dynammics summary (CN average)

def Plot_Eval_BiomassDynamicsSummaryAverage_CN(meta,gplt):

    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)
    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) &
                         (gplt['Ctot G Tot']>0) & (gplt['Ctot G Tot']<30))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot G Tot']['N']>=3)[0]
    for v in gplt['vaL']:
        for k in d[v].keys():
            d[v][k]=d[v][k][ind]
    u=u[ind]
    lab=lab[ind]

    # Area weighting
    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]
    for v in d:
        d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
        d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
    lab=np.append(lab,'Weighted\naverage')
    u=np.append(u,0.0)

    cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
    cle=[0.05,0.2,0.45]
    lab=['Gross\ngrowth','Mortality','Net\ngrowth'] #,'Harvest\nmortality'

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]]); barw=0.38
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6))
    ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
    ax.bar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],barw,facecolor=cl[1,:],label='Predictions (FCS)')
    ax.errorbar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],yerr=d['Ctot G Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],yerr=d['Mod C_G_Gross_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],barw,facecolor=cl[0,:])
    ax.bar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],barw,facecolor=cl[1,:])
    ax.errorbar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],yerr=d['Ctot Mort+Lost']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],yerr=d['Mod C_M_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],barw,facecolor=cl[0,:])
    ax.bar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],barw,facecolor=cl[1,:])
    ax.errorbar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],yerr=d['Ctot Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],yerr=d['Mod C_G_Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-2,3,0.5),ylabel='Carbon balance of trees (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,3.5],ylim=[-2,2])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    ax.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage__CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Biomass dynammics summary (average YSM)

def Plot_EvalBiomassDynamicsSummary_Average_YSM(meta,gplt):
    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)
    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF YSM']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) &
                         (gplt['Ctot G Tot']>0) & (gplt['Ctot G Tot']<30))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot G Tot']['N']>=3)[0]
    for v in gplt['vaL']:
        for k in d[v].keys():
            d[v][k]=d[v][k][ind]
    u=u[ind]
    lab=lab[ind]

    # Area weighting
    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]
    for v in d:
        d[v]['mu']=np.append(d[v]['mu'],np.sum(d[v]['mu']*Area)/np.sum(Area))
        d[v]['se']=np.append(d[v]['se'],np.sum(d[v]['se']*Area)/np.sum(Area))
    lab=np.append(lab,'Weighted\naverage')
    u=np.append(u,0.0)

    cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
    cle=[0.05,0.2,0.45]
    lab=['Gross\ngrowth','Natural\nmortality','Net\ngrowth'] #,'Harvest\nmortality'

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]]); barw=0.38
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,6))
    ax.plot([0,5],[0,0],'k-',color=meta['Graphics']['gp']['cla'],lw=meta['Graphics']['gp']['lw1'])
    ax.bar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],barw,facecolor=cl[1,:],label='Predictions (FCS)')
    ax.errorbar(1-barw/2-0.01,d['Ctot G Tot']['mu'][-1],yerr=d['Ctot G Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(1+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'][-1],yerr=d['Mod C_G_Gross_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],barw,facecolor=cl[0,:])
    ax.bar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],barw,facecolor=cl[1,:])
    ax.errorbar(2-barw/2-0.01,-d['Ctot Mort+Lost']['mu'][-1],yerr=d['Ctot Mort+Lost']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(2+barw/2+0.01,-d['Mod C_M_Tot']['mu'][-1],yerr=d['Mod C_M_Tot']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.bar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],barw,facecolor=cl[0,:])
    ax.bar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],barw,facecolor=cl[1,:])
    ax.errorbar(3-barw/2-0.01,d['Ctot Net']['mu'][-1],yerr=d['Ctot Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(3+barw/2+0.01,d['Mod C_G_Net']['mu'][-1],yerr=d['Mod C_G_Net']['se'][-1],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

    ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,yticks=np.arange(-4,6,1),ylabel='Carbon balance of trees (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,3.5],ylim=[-3,4])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    ax.legend(frameon=False,loc='lower left',facecolor=[1,1,1],labelspacing=0.25)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Gross growth (CN)

def Plot_EvalGrossGrowth_CN(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) &
                         (gplt['Ctot G Tot']>0) & (gplt['Ctot G Tot']<30))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot G Tot']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
    #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
    for i in range(ikp.size):
        ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
    ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed gross growth (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted gross growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0,3.6],ylim=[0,3.6])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_CN','png',meta['Graphics']['gp']['save fig dpi'])
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
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Gross growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,3.5])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Gross growth (YSM)

def Plot_EvalGrossGrowth_YSM(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF YSM']==1) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) &
                         (gplt['Ctot G Tot']>0) & (gplt['Ctot G Tot']<30))[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot G Tot']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
    #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
    for i in range(ikp.size):
        ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
    ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed gross growth (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted gross growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0,5.5],ylim=[0,5.5])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_YSM','png',meta['Graphics']['gp']['save fig dpi'])
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
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot G Tot']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Gross_Tot']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Gross growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,5])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Mortality

def Plot_EvalMortality_CN(meta,gplt):

    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Mod C_M_Tot']>=0) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot Mort+Lost']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
    ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed mortality (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted mortality (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0,3],ylim=[0,3])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_CN','png',meta['Graphics']['gp']['save fig dpi'])
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
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Mortality (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,6])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Barchart_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Mortality (YSM)

def Plot_EvalMortality_YSM(meta,gplt):
    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF YSM']==1) &
                         (gplt['Mod C_M_Tot']>=0) &
                         (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<10000) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot Mort+Lost']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(2.8,0.3,txt,fontsize=10,color='k',ha='right')
    ax.text(2.7,2.7,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed mortality (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted mortality (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[0,3],ylim=[0,3])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM','png',meta['Graphics']['gp']['save fig dpi'])
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
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Mort+Lost']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_M_Tot']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Mortality (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,6])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Net growth (CN)

def Plot_EvalGrowthNet_CN(meta,gplt):
    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot Net']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(1.8,-2.75,txt,fontsize=10,color='k',ha='right')
    ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-3,2],ylim=[-3,2])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_CN','png',meta['Graphics']['gp']['save fig dpi'])
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
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.plot([-1,20],[0,0],'k-',lw=0.5)
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Net growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[-3,2])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Net growth (YSM)

def Plot_EvalGrowthNet_YSM(meta,gplt):
    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF YSM']==1) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Ctot Net']['N']>=3)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(1.8,-2.75,txt,fontsize=10,color='k',ha='right')
    ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-3,5],ylim=[-3,5])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_YSM','png',meta['Graphics']['gp']['save fig dpi'])

    # Plot bar chart

    # Put in order
    ord=np.argsort(d['Ctot Net']['mu'])
    lab=np.flip(lab[ord])
    for v in d:
        for k in d[v].keys():
            d[v][k]=np.flip(d[v][k][ord])

    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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

    cl=np.array([[0.85,0.85,0.95],[0.75,0.85,0.75]])
    barw=0.32
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.plot([-1,20],[0,0],'k-',lw=0.5)
    ax.bar(np.arange(u.size)-barw/2-0.01,d['Ctot Net']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['Mod C_G_Net']['mu'],barw,facecolor=cl[1,:],label='Predictions (FCS)')
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
    ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),xticklabels=lab,ylabel='Net growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[-3,2])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_YSM','png',meta['Graphics']['gp']['save fig dpi'])
    return

#%% Net growth (TIPSY stemwood)

def Plot_EvalGrowthNetStemwoodTIPSY_CN(meta,gplt):
    # Unique BGC zones
    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Csw Mort']<100) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Csw Net']['N']>=10)[0]
    for v in gplt['vaL']:
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
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    ax.text(1.8,-0.75,txt,fontsize=10,color='k',ha='right')
    ax.text(1.8,1.8,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed net growth (MgC ha$^{-1}$ yr$^{-1}$)',ylabel='Predicted net growth (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-1,2],ylim=[-1,2])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN','png',meta['Graphics']['gp']['save fig dpi'])

    # With no disturbance mortality
    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['Ground Plots'],'Ecozone BC L1',u[i])
        for v in gplt['vaL']:
            ind=np.where( (gplt['Ecozone BC L1']==u[i]) &
                         (gplt['PTF CN']==1) &
                         (gplt['Csw Mort']<0.2) )[0]
            d[v]['N'][i]=ind.size
            d[v]['mu'][i]=np.nanmean(gplt[v][ind])
            d[v]['sd'][i]=np.nanstd(gplt[v][ind])
            d[v]['se'][i]=np.nanstd(gplt[v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['Csw Net']['N']>=10)[0]
    for v in gplt['vaL']:
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

    ax.plot(rs['xhat'],rs['yhat'],'r-',lw=1,label='Best fit')

    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%% Soil organic carbon

def Plot_EvalSOC(meta,gplt):

    u=np.unique(gplt['Ecozone BC L1'][gplt['Ecozone BC L1']>0])
    lab=np.array(['' for _ in range(u.size)],dtype=object)

    d={}
    for v in gplt['soils']['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)
        for i in range(u.size):
            lab[i]=u1ha.lut_n2s(meta['u1ha']['bgcz'],u[i])[0]
            ind=np.where(gplt['soils']['becz']==u[i])[0]
            if ind.size>0:
                ind=np.where( (gplt['soils']['becz']==u[i]) & (gplt['soils']['TOT_C_THA']>0) )[0]
                d[v]['N'][i]=ind.size
                d[v]['mu'][i]=np.nanmean(gplt['soils'][v][ind])
                d[v]['sd'][i]=np.nanstd(gplt['soils'][v][ind])
                d[v]['se'][i]=np.nanstd(gplt['soils'][v][ind])/np.sqrt(ind.size)

    # Remove classes with inadequate data
    ind=np.where(d['TOT_C_THA']['N']>=10)[0]
    for v in gplt['soils']['vaL']:
        for k in d[v].keys():
            d[v][k]=d[v][k][ind]
    u=u[ind]
    lab=lab[ind]

    # Scatterplot
    x=d['TOT_C_THA']['mu']
    y=d['C_Soil_Tot']['mu']
    ikp=np.where( (np.isnan(x+y)==False) & (d['TOT_C_THA']['N']>=10) )[0]
    rs,txt=gu.GetRegStats(x[ikp],y[ikp])

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
    ax.plot([0,1000],[0,1000],'-k',lw=2,color=[0.8,0.8,0.8])
    #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
    ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    for i in range(ikp.size):
        ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
    ax.text(330,30,txt,fontsize=10,color='k',ha='right')
    ax.text(300,300,'1:1',fontsize=8,ha='center')
    ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed SOC (MgC ha$^{-1}$)',ylabel='Predicted SOC (MgC ha$^{-1}$)',xlim=[0,350],ylim=[0,350])
    #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_Scatterplot','png',meta['Graphics']['gp']['save fig dpi'])
    plt.close('all')

    # Plot bar chart

    # Put in order
    ord=np.argsort(d['TOT_C_THA']['mu'])
    lab=np.flip(lab[ord])
    for v in d:
        for k in d[v].keys():
            d[v][k]=np.flip(d[v][k][ord])

    Area=np.zeros(lab.size)
    for i in range(lab.size):
        ind1=np.where(meta['par']['By BGC']['Name']==lab[i])[0]
        Area[i]=meta['par']['By BGC']['Area Treed (Mha)'][ind1]

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
    cl=np.array([[0.9,0.9,1.0],[0.8,0.8,0.9],[0.9,1,0.9],[0.8,0.9,0.8]])

    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
    ax.bar(np.arange(u.size)-barw/2-0.01,d['MIN_C_THA']['mu'],barw,facecolor=cl[0,:],label='Ground plot observations, mineral horizon (Shaw et al. 2018)')
    ax.bar(np.arange(u.size)-barw/2-0.01,d['ORG_C_THA']['mu'],barw,facecolor=cl[1,:],bottom=d['MIN_C_THA']['mu'],label='Ground plot observations, organic horizon (Shaw et al. 2018)')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['C_Soil_Tot']['mu'],barw,facecolor=cl[2,:],label='Prediction, mineral horizon (FCS)')
    ax.bar(np.arange(u.size)+barw/2+0.01,d['C_Soil_OHorizon']['mu'],barw,facecolor=cl[3,:],bottom=d['C_Soil_Tot']['mu'],label='Prediction, organic horizon (FCS)')
    ax.errorbar(np.arange(u.size)-barw/2-0.01,d['MIN_C_THA']['mu']+d['ORG_C_THA']['mu'],yerr=d['TOT_C_THA']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(np.arange(u.size)+barw/2+0.01,d['C_Soil_Tot']['mu']+d['C_Soil_OHorizon']['mu'],yerr=d['C_Soil_Tot']['se']+d['C_Soil_OHorizon']['se'],color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    for i in range(u.size):
        if Dp[i]>=0:
            a='+'
        else:
            a=''
        ax.text(i+barw/2+0.01,yp[i]+d['C_Soil_Tot']['se'][i]+d['C_Soil_OHorizon']['se'][i]+13,a + str(Dp[i].astype(int)) + '%',color=meta['Graphics']['gp']['cla'],ha='center',fontsize=6)
    #    ax.text(i,8,str(d['TOT_C_THA']['N'][i].astype(int)),color=meta['Graphics']['gp']['cla'],ha='center',fontsize=7)
    #    #ax.text(i,30,str(d['Model SS'][i].astype(int)),color='c',ha='center',fontsize=8)
    ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,450],xticks=np.arange(u.size),
           xticklabels=lab,ylabel='Soil organic carbon (MgC ha$^{-1}$ yr$^{-1}$)')
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_Scatterplot','png',meta['Graphics']['gp']['save fig dpi'])

    return

#%%

def Plot_EvalAgeResponsesCoast(meta,gplt):

    # Coast
    ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==meta['Ground Plots']['LUT']['Ecozone BC L1']['CWH']) | (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==meta['Ground Plots']['LUT']['Ecozone BC L1']['MH']) )[0]

    # All
    #ind=np.where( (gplt['PTF CN']==1) | (gplt['PTF YSM']==1) )[0]

    x=gplt['Age VRI t0'][ind]
    xM=gplt['Mod A t0'][ind]
    bw=25; bin=np.arange(bw,300,bw)
    xhat=np.arange(1,301,1)

    lw=1; ms=4; mec='b'; mfc='w'; cl='b';
    plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10))
    y=gplt['Ctot L t0'][ind]
    N,mu,med,sig,se=gu.discres(x,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[0,0].plot(xhat,yhat,'b-',lw=0.5)
    ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Observations')

    y=gplt['Mod C_Biomass_Tot t0'][ind]
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[0,0].plot(xhat,yhat,'r-',lw=0.5)
    ax[0,0].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r',label='Model (FCS)')
    ax[0,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
    ax[0,0].legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')

    y=gplt['Ctot G Tot'][ind]
    N,mu,med,sig,se=gu.discres(x,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[0,1].plot(xhat,yhat,'b-',lw=0.5)
    ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Observations')

    y=gplt['Mod C_G_Gross_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[0,1].plot(xhat,yhat,'r-',lw=0.5)
    ax[0,1].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r',label='Model (FCS)')

    ax[0,1].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])


    y=gplt['Ctot Mort+Lost'][ind]
    N,mu,med,sig,se=gu.discres(x,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[1,0].plot(xhat,yhat,'b-',lw=0.5)
    ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')

    y=gplt['Ctot Mort+Lost'][ind]-gplt['Ctot Mort+Lost Fire'][ind]-gplt['Ctot Mort+Lost Insect'][ind]
    N,mu,med,sig,se=gu.discres(x,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    #ax[1,0].plot(xhat,yhat,'g--',lw=0.5)
    #ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')

    y=gplt['Mod C_M_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[1,0].plot(xhat,yhat,'r-',lw=0.5)
    ax[1,0].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r')

    ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

    ax[1,1].plot(xhat,0*xhat,'k-',lw=1.5,color=[0.8,0.8,0.8])
    y=gplt['Ctot Net'][ind]
    N,mu,med,sig,se=gu.discres(x,y,bw,bin)
    yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
    ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5)
    ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')

    y=gplt['Mod C_G_Net'][ind]
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
    yhat=np.interp(xhat,bin,mu)
    ax[1,1].plot(xhat,yhat,'r-',lw=0.5)
    ax[1,1].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r')

    flg=0
    if flg==1:
        y=gplt['GY Csw Net'][ind]
        N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
        yhat=np.interp(xhat,bin,mu)
        ax[1,1].plot(xhat,yhat,'c-',lw=0.5)
        ax[1,1].plot(bin,mu,'c^',ms=ms,lw=lw,color='c',mfc=mfc,mec='c')

    ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
    ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
    gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponses_Coast','png',meta['Graphics']['gp']['save fig dpi'])

    return

def Plot_EvalAgeResponsesInterior(meta,gplt):

    # Interior
    ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']!=meta['Ground Plots']['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=meta['Ground Plots']['LUT']['Ecozone BC L1']['MH']) & (gplt['Ecozone BC L1']!=meta['Ground Plots']['LUT']['Ecozone BC L1']['ICH']) )[0];

    x=gplt['Age VRI t0'][ind];
    xM=gplt['Mod A t0'][ind];
    bw=25; bin=np.arange(bw,300,bw);
    xhat=np.arange(1,301,1);

    lw=1; ms=4; mec='b'; mfc='w'; cl='b';
    plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10));
    y=gplt['Ctot L t0'][ind];
    N,mu,med,sig,se=gu.discres(x,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[0,0].plot(xhat,yhat,'b-',lw=0.5);
    ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Observations');

    y=gplt['Mod C_Biomass_Tot t0'][ind];
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[0,0].plot(xhat,yhat,'r-',lw=0.5);
    ax[0,0].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r',label='Model (FCS)');
    ax[0,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300]);
    ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl']);
    ax[0,0].legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w');

    y=gplt['Ctot G Tot'][ind];
    N,mu,med,sig,se=gu.discres(x,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[0,1].plot(xhat,yhat,'b-',lw=0.5);
    ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Observations');

    y=gplt['Mod C_G_Gross_Tot'][ind];
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[0,1].plot(xhat,yhat,'r-',lw=0.5);
    ax[0,1].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r',label='Model (FCS)');

    ax[0,1].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300]);
    ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl']);

    y=gplt['Ctot Mort+Lost'][ind];
    N,mu,med,sig,se=gu.discres(x,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[1,0].plot(xhat,yhat,'b-',lw=0.5);
    ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total');

    y=gplt['Ctot Mort+Lost'][ind]-gplt['Ctot Mort+Lost Fire'][ind]-gplt['Ctot Mort+Lost Insect'][ind];
    N,mu,med,sig,se=gu.discres(x,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    #ax[1,0].plot(xhat,yhat,'g--',lw=0.5)
    #ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')

    y=gplt['Mod C_M_Tot'][ind];
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[1,0].plot(xhat,yhat,'r-',lw=0.5);
    ax[1,0].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r');

    ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300]);
    ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl']);

    ax[1,1].plot(xhat,0*xhat,'k-',lw=1.5,color=[0.8,0.8,0.8]);
    y=gplt['Ctot Net'][ind];
    N,mu,med,sig,se=gu.discres(x,y,bw,bin);
    yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False]);
    ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5);
    ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total');

    y=gplt['Mod C_G_Net'][ind];
    N,mu,med,sig,se=gu.discres(xM,y,bw,bin);
    yhat=np.interp(xhat,bin,mu);
    ax[1,1].plot(xhat,yhat,'r-',lw=0.5);
    ax[1,1].plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r');

    flg=0
    if flg==1:
        y=gplt['GY Csw Net'][ind]
        N,mu,med,sig,se=gu.discres(xM,y,bw,bin)
        yhat=np.interp(xhat,bin,mu)
        ax[1,1].plot(xhat,yhat,'c-',lw=0.5)
        ax[1,1].plot(bin,mu,'c^',ms=ms,lw=lw,color='c',mfc=mfc,mec='c')

    ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300]);
    ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl']);
    ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w');
    gu.axletters(ax,plt,0.03,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold');
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponses_Interior','png',meta['Graphics']['gp']['save fig dpi']);

    return