
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
import matplotlib.colors
import matplotlib.ticker as ticker
#from matplotlib import animation
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_inventory as uinv
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
import fcexplore.psp.Processing.psp_util as ugp

#%% Age

def EvalAgeByBGCZ_CNV(meta,pNam,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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

    # # Scatterplot
    # x=d['Age VRI t0']['mu']
    # y=d['Mod A t0']['mu']
    # ikp=np.where( (np.isnan(x+y)==False) & (d['Age VRI t0']['N']>=30) )[0]
    # rs,txt=gu.GetRegStats(x[ikp],y[ikp])

    # fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
    # ax.plot([0,500],[0,500],'-k',lw=2,color=[0.75,0.75,0.75])
    # #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
    # for i in range(ikp.size):
    #     ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
    # ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    # ax.text(215,30,txt,fontsize=10,color='k',ha='right')
    # ax.text(230,230,'1:1',fontsize=8,ha='center')
    # ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed age (years)',ylabel='Predicted age (years)',xlim=[0,250],ylim=[0,250])
    # #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    # ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    # if meta['Graphics']['Print Figures']=='On':
    #     gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Scatterplot_CNV','png',900)

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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Age_ByBGCZone_Barchart_CNV','png',900)

    return

#%%

def EvalBiomassByBGC_CNV(meta,pNam,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
    #     ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    #     ax.text(200,20,txt,fontsize=10,color='k',ha='right')
    #     ax.text(190,190,'1:1',fontsize=8,ha='center')
    #     ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed biomass (MgC ha$^{-1}$)',ylabel='Predicted biomass (MgC ha$^{-1}$)',xlim=[0,220],ylim=[0,220])
    #     #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    #     ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    #     if meta['Graphics']['Print Figures']=='On':
    #         gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_CNV','png',900)

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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_CNV','png',900)
        pass

    return

#%% Biomass (YSM)

def EvalBiomassByBGC_YSM(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_Scatterplot_YSM','png',900)

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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Biomass_ByBGCZone_BarChart_YSM','png',900)

    return

#%% Biomass (TIPSY stemwood)

def EvalStemwoodFromTIPSYByBGC_CNV(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_BiomassSW_ByBGCZone_Scatterplot','png',900)

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

def Eval_BiomassComponents_CN(meta,gplt):
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassComponents_CN','png',900)

    return

#%% Biomass dynammics summary (total CO2e)

def Eval_BiomassDynamicsAve_TotCO2e_CN(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsTotal_CN','png',900)

    return

#%% Biomass dynammics summary (CN average)

def Eval_BiomassDynamicsAve_CN(meta,pNam,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
        Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage__CN','png',900)

    return

#%% Biomass dynammics summary (average YSM)

def Plot_BiomassDynamicsAve_YSM(meta,gplt):
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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        ind1=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[i])[0]
        Area[i]=meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind1]
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_BiomassDynamicsAverage_YSM','png',900)

    return

#%% Gross growth (CN)

def EvalGrossGrowthByBGC_CN(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_CN','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_CN','png',900)

    return

#%% Gross growth (YSM)

def EvalGrossGrowthByBGC_YSM(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Scatterplot_YSM','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthGross_ByBGCZone_Barchart_YSM','png',900)

    return

#%% Mortality

def EvalMortalityByBGC_CN(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_CN','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Barchart_CN','png',900)

    return

#%% Mortality (YSM)

def EvalMortalityByBGC_YSM(meta,gplt):

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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_ByBGCZone_Scatterplot_YSM','png',900)

    return

#%% Net growth (CN)

def EvalGrowthNetByBGC_CN(meta,gplt):
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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_CN','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_CN','png',900)

    return

#%% Net growth (YSM)

def EvalGrowthNetByBGC_YSM(meta,gplt):
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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Scatterplot_YSM','png',900)

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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_GrowthNet_ByBGCZone_Barchart_YSM','png',900)
    return

#%% Net growth (TIPSY stemwood)

def EvalGrowthNetStemwoodTIPSY_CN(meta,gplt):
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
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN','png',900)

    # With no disturbance mortality
    d={}
    for v in gplt['vaL']:
        d[v]={}
        d[v]['N']=np.zeros(u.size)
        d[v]['mu']=np.zeros(u.size)
        d[v]['sd']=np.zeros(u.size)
        d[v]['se']=np.zeros(u.size)

    for i in range(u.size):
        lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_TIPSY_GrowthNetSW_ByBGCZone_Scatterplot_CN','png',900)

    return

#%% Soil organic carbon

def EvalSOCByBGC_ShawComp(meta,pNam,gplt):

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
            lab[i]=u1ha.lut_n2s(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'],u[i])[0]
            ind=np.where(gplt['soils']['bgcz']==u[i])[0]
            if ind.size>0:
                ind=np.where( (gplt['soils']['bgcz']==u[i]) & (gplt['soils']['TOT_C_THA']>0) )[0]
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

    # # Scatterplot
    # x=d['TOT_C_THA']['mu']
    # y=d['C_Soil_Tot']['mu']
    # ikp=np.where( (np.isnan(x+y)==False) & (d['TOT_C_THA']['N']>=10) )[0]
    # rs,txt=gu.GetRegStats(x[ikp],y[ikp])

    # plt.close('all')
    # fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
    # ax.plot([0,1000],[0,1000],'-k',lw=2,color=[0.8,0.8,0.8])
    # #ax.plot(x[ikp],y[ikp],'ko',mfc=[0.29,0.49,0.78],mec=[0.29,0.49,0.78],lw=0.5,ms=5)
    # ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
    # for i in range(ikp.size):
    #     ax.text(x[ikp[i]],y[ikp[i]],lab[ikp[i]],color='k',ha='center',fontsize=8)
    # ax.text(330,30,txt,fontsize=10,color='k',ha='right')
    # ax.text(300,300,'1:1',fontsize=8,ha='center')
    # ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed SOC (MgC ha$^{-1}$)',ylabel='Predicted SOC (MgC ha$^{-1}$)',xlim=[0,350],ylim=[0,350])
    # #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    # ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    # if meta['Graphics']['Print Figures']=='On':
    #     gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_Scatterplot','png',900)
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
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_SOC_ByBGCZone_Scatterplot','png',900)

    return

#%%

def EvalAgeResponsesBiomassAndNetGrowth_ByReg(meta,pNam,gpt):

    bw=25; bin=np.arange(bw,250+bw,bw)
    lw=0.5; ms=3; clO=np.array([0.8,0.85,1]); clM=np.array([0.85,1,0.8]);
    plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,6))

    # Coast
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt['Age Med t0'][ind]
    #xO=gpt['Age VRI t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot L t0'][ind]
    yM=gpt['Mod C_Biomass_Tot t0'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[0].plot(bin,mu,'-ko',ms=ms,lw=lw,mew=lw,color=0.7*clO,mfc='w',mec=0.7*clO,label='Observed biomass',zorder=1)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[0].plot(bin,mu,'--ks',ms=ms,lw=lw,mew=lw,color=0.7*clM,mfc='w',mec=0.7*clM,label='Predicted biomass',zorder=1)
    ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,50),xlim=[0,250+bw],ylim=[0,400])
    ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])

    ax2=ax1[0].twinx()
    yO=gpt['Ctot Net'][ind]
    yM=gpt['Mod C_G_Net'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax2.bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO,zorder=-1)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax2.bar(bin+(0.45*bw/2),mu,0.45*bw,ec='none',fc=clM,zorder=-1)
    ax2.plot([0,500],[0,0],'-k',lw=lw)
    ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
    ax1[0].set_zorder(ax2.get_zorder()+1)
    ax1[0].patch.set_visible(False)
    ax2.tick_params(length=meta['Graphics']['gp']['tickl'])

    # Interior
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    #xO=gpt['Age VRI t0'][ind]
    xO=gpt['Age Med t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot L t0'][ind]
    yM=gpt['Mod C_Biomass_Tot t0'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[1].plot(bin,mu,'-ko',ms=ms,lw=lw,mew=lw,color=0.7*clO,mfc='w',mec=0.7*clO,zorder=1)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[1].plot(bin,mu,'--ks',ms=ms,lw=lw,mew=lw,color=0.7*clM,mfc='w',mec=0.7*clM,zorder=1)
    ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,20),xlim=[0,250+bw],ylim=[0,200])
    ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
    ax3=ax1[1].twinx()
    yO=gpt['Ctot Net'][ind]
    yM=gpt['Mod C_G_Net'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax3.bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO,label='Observed net growth',zorder=-1)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax3.bar(bin+(0.45*bw/2),mu,0.45*bw,ec='none',fc=clM,label='Predicted net growth',zorder=-1)
    ax3.plot([0,500],[0,0],'-k',lw=lw)
    ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
    ax1[1].set_zorder(ax3.get_zorder()+1)
    ax1[1].patch.set_visible(False)
    ax3.tick_params(length=meta['Graphics']['gp']['tickl'])
    ax1[0].legend(loc='lower right',frameon=False,facecolor=None,edgecolor='w')
    ax3.legend(loc='upper center',frameon=False,facecolor=None,edgecolor='w')
    gu.axletters(ax1,plt,0.04,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
    plt.tight_layout()
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseBiomassNetGrowth','png',900)
    return

def EvalAgeResponsesGrossGrowth_ByReg(meta,pNam,gpt):
    bw=25; bin=np.arange(bw,250+bw,bw)
    clO=np.array([0.8,0.85,1]); clM=np.array([0.85,1,0.8]);
    plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,5.5))
    # Coast
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt['Age Med t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot G Tot'][ind]
    yM=gpt['Mod C_G_Gross_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[0].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
    ax1[0].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[0].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
    ax1[0].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax1[0].set(ylabel='Gross growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,7])
    ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
    # Interior
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt['Age Med t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot G Tot'][ind]
    yM=gpt['Mod C_G_Gross_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[1].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
    ax1[1].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[1].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
    ax1[1].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax1[1].set(ylabel='Gross growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,4])
    ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
    gu.axletters(ax1,plt,0.04,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
    plt.tight_layout()
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseGrossGrowth','png',900)
    return

def EvalAgeResponsesMortality_ByReg(meta,pNam,gpt):
    bw=25; bin=np.arange(bw,250+bw,bw)
    clO=np.array([0.8,0.85,1]); clM=np.array([0.85,1,0.8]);
    plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,5.5))
    # Coast
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt['Age Med t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot Mort+Lost'][ind]
    yM=gpt['Mod C_M_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[0].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
    ax1[0].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[0].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
    ax1[0].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax1[0].set(ylabel='Mortality (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,8])
    ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
    # Interior
    ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
    xO=gpt['Age Med t0'][ind]
    xM=gpt['Mod A t0'][ind]
    yO=gpt['Ctot Mort+Lost'][ind]
    yM=gpt['Mod C_M_Tot'][ind]
    N,mu,med,sig,se=gu.discres(xO,yO,bw,bin)
    ax1[1].bar(bin-(0.45*bw/2),mu,0.4*bw,ec='none',fc=clO)
    ax1[1].errorbar(bin-(0.45*bw/2),mu,yerr=se,color=0.5*clO,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    N,mu,med,sig,se=gu.discres(xM,yM,bw,bin)
    ax1[1].bar(bin+(0.45*bw/2),mu,0.4*bw,ec='none',fc=clM)
    ax1[1].errorbar(bin+(0.45*bw/2),mu,yerr=se,color=0.5*clM,fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax1[1].set(ylabel='Mortality (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,500,1),xlim=[0,250+bw],ylim=[0,5])
    ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
    gu.axletters(ax1,plt,0.04,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
    plt.tight_layout()
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_AgeResponseMortality','png',900)
    return

#%%

def QA_Mortality_Regression(meta,gpt):
    # Modelled harvest mortality is too high because salvage harvest is being counted as mortality
    # even though something else killed most of the trees.
    gpt2=copy.deepcopy(gpt)
    for k in gpt2.keys():
        try:
            gpt2[k]=gpt2[k][(gpt['PTF CNY']==1)]
        except:
            pass
    del gpt2['vaL'],gpt2['soils']

    df=pd.DataFrame.from_dict(gpt2)
    df.columns=df.columns.str.replace(' ','_')
    #x=df[['debt_ratio','industry']]
    #y=df['cash_flow']

    #frm='Ctot_Mort ~ Age_t0 + C(IBB_Occurrence) + C(IBD_Occurrence) + C(IBM_Occurrence) + C(IBS_Occurrence) + C(IDW_Occurrence) + C(Wildfire_Occurrence) + C(Harv_Occurrence)'
    frm1='Ctot_Mort ~ Ctot_L_t0 + C(Occ_IBB) + C(Occ_IBM) + C(Occ_IBS) + C(Occ_IDW) + C(Occ_Wildfire) + C(Occ_Harv)'
    md1=smf.ols(formula=frm1,data=df)
    mr1=md1.fit()
    print(mr1.summary())

    frm2='Mod_C_M_Tot ~ Mod_C_Biomass_Tot_t0 + C(Occ_IBB) + C(Occ_IBM) + C(Occ_IBS) + C(Occ_IDW) + C(Occ_Wildfire) + C(Occ_Harv)'
    md2=smf.ols(formula=frm2,data=df)
    mr2=md2.fit()
    print(mr2.summary())

    b1=mr1.params.values
    b2=mr2.params.values
    se1=mr1.bse.values
    se2=mr2.bse.values
    lab=['IBB','IBM','IBS','IDW','Wildfire','Harvest'] #'Int', ,'Initial\nbiomass'

    # Correct harvest to account for salvage
    b2[-2]=0.65*b2[-2]

    barw=0.32
    cl=np.array([[0.84,0.87,1.0],[0.8,0.9,0.8]])

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,9))
    ax.plot([-1,b1.size+1],[0,0],'k-',lw=0.5)
    ax.bar(np.arange(b1.size)-barw/2-0.01,b1,barw,facecolor=cl[0,:],label='Observations')
    ax.bar(np.arange(b1.size)+barw/2+0.01,b2,barw,facecolor=cl[1,:],label='Predictions (FCS)')
    ax.errorbar(np.arange(b1.size)-barw/2-0.01,b1,yerr=se1,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.errorbar(np.arange(b1.size)+barw/2+0.01,b2,yerr=se2,color=meta['Graphics']['gp']['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
    ax.set(xticks=np.arange(1,b1.size-1,1),xticklabels=lab,ylabel='Regression coefficient (tC ha$^{-1}$ yr$^{-1}$)',
           xlim=[-0.5+1,b1.size-1-0.5],ylim=[-4,10.5])
    plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
    plt.tight_layout()
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_Mortality_Regression','png',900)
    return

#%%

def Prepare_ObsVsModComparison(meta,pNam):

    # Import ground plot data
    metaGP,gplt=ugp.ImportPlotData(meta,type='Stand')
    gplt['Ctot G Tot']=gplt['Ctot G Surv']+gplt['Ctot G Recr']
    gplt['Csw Net']=gplt['Csw G Surv']+gplt['Csw G Recr']-gplt['Csw Mort']

    # Import BC1ha data
    #z=u1ha.Import_Raster(meta,[],['refg','lcc1_c','harv_yr_con1'])

    # Import modelled data
    tvSaved=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

    vmL=['A','C_Biomass_Tot','C_Stemwood_Tot','C_Foliage_Tot','C_Branch_Tot',
         'C_Bark_Tot','C_Root_Tot','C_DeadWood_Tot','C_G_Gross_Tot','C_G_Net_Tot',
         'C_M_Reg_Tot','C_M_Dist','C_Soil_Tot','C_Soil_OHorizon','C_M_Harvest_Salvage']
    dM={}
    for v in vmL:
        dM[v]=np.zeros( (tvSaved.size,meta[pNam]['Project']['N Stand']) )

    tvFull=np.arange(meta[pNam]['Project']['Year Start'],meta[pNam]['Project']['Year End']+1,1)

    iT2=np.where(tvFull>=tvSaved[0])[0]
    for iEns in range(meta[pNam]['Project']['N Ensemble']):
        for iBat in range(meta[pNam]['Project']['N Batch']):
            indBat=cbu.IndexToBatch(meta[pNam],iBat)
            d1=cbu.LoadSingleOutputFile(meta,pNam,0,iEns,iBat)
            for v in vmL:
                if v=='C_M_Harvest_Salvage':
                    idx=d1['C_M_ByAgent']['Harvest Salvage']['idx']
                    d1[v]=np.zeros( (tvFull.size,indBat.size),dtype='float32')
                    d1[v][idx[0],idx[1]]=meta['Core']['Scale Factor C_M_ByAgent']*d1['C_M_ByAgent']['Harvest Salvage']['M'].astype('float32')
                    dM[v][:,indBat]=dM[v][:,indBat]+d1[v][iT2,:].copy()
                else:
                    dM[v][:,indBat]=dM[v][:,indBat]+d1[v]
    for v in vmL:
        dM[v]=dM[v]/meta[pNam]['Project']['N Ensemble']

    # Calculate total mortality
    dM['C_M_Tot']=dM['C_M_Reg_Tot']+dM['C_M_Dist']-dM['C_M_Harvest_Salvage']
    dM['C_G_Net']=dM['C_G_Gross_Tot']-dM['C_M_Tot']

    # Add simulations to gplt structure
    for k in dM.keys():
        gplt['Mod ' + k + ' t0']=np.nan*np.ones(gplt['Year t0'].size)
        gplt['Mod ' + k + ' t1']=np.nan*np.ones(gplt['Year t0'].size)
        gplt['Mod ' + k]=np.nan*np.ones(gplt['Year t0'].size)

    for i in range(meta['Geos']['Sparse']['X'].size):
        iS=np.where( (np.abs(gplt['X']-meta['Geos']['Sparse']['X'][i])<=500) & (np.abs(gplt['Y']-meta['Geos']['Sparse']['Y'][i])<=500) )[0]
        if iS.size==0:
            continue
        for j in range(iS.size):
            if np.isnan(gplt['Year t1'][iS[j]])==True:
                iT=np.where( (tvSaved==gplt['Year t0'][iS[j]]) )[0]
                gplt['Mod A t0'][iS[j]]=dM['A'][iT,i]
                gplt['Mod C_Biomass_Tot t0'][iS[j]]=dM['C_Biomass_Tot'][iT,i]
                gplt['Mod C_Stemwood_Tot t0'][iS[j]]=dM['C_Stemwood_Tot'][iT,i]
                gplt['Mod C_Foliage_Tot t0'][iS[j]]=dM['C_Foliage_Tot'][iT,i]
                gplt['Mod C_Branch_Tot t0'][iS[j]]=dM['C_Branch_Tot'][iT,i]
                gplt['Mod C_Bark_Tot t0'][iS[j]]=dM['C_Bark_Tot'][iT,i]
                gplt['Mod C_Root_Tot t0'][iS[j]]=dM['C_Root_Tot'][iT,i]
                gplt['Mod C_DeadWood_Tot t0'][iS[j]]=dM['C_DeadWood_Tot'][iT,i]
            else:
                iT=np.where( (tvSaved>=gplt['Year t0'][iS[j]]) & (tvSaved<=gplt['Year t1'][iS[j]]) )[0]
                gplt['Mod A t0'][iS[j]]=dM['A'][iT[0],i]
                gplt['Mod C_Biomass_Tot t0'][iS[j]]=dM['C_Biomass_Tot'][iT[0],i]
                gplt['Mod C_Stemwood_Tot t0'][iS[j]]=dM['C_Stemwood_Tot'][iT[0],i]
                gplt['Mod C_Foliage_Tot t0'][iS[j]]=dM['C_Foliage_Tot'][iT[0],i]
                gplt['Mod C_Branch_Tot t0'][iS[j]]=dM['C_Branch_Tot'][iT[0],i]
                gplt['Mod C_Bark_Tot t0'][iS[j]]=dM['C_Bark_Tot'][iT[0],i]
                gplt['Mod C_Root_Tot t0'][iS[j]]=dM['C_Root_Tot'][iT[0],i]
                gplt['Mod C_Biomass_Tot t1'][iS[j]]=dM['C_Biomass_Tot'][iT[-1],i]
                gplt['Mod C_DeadWood_Tot t0'][iS[j]]=dM['C_DeadWood_Tot'][iT[0],i]
                gplt['Mod C_DeadWood_Tot t1'][iS[j]]=dM['C_DeadWood_Tot'][iT[-1],i]
                gplt['Mod C_M_Tot'][iS[j]]=np.mean(dM['C_M_Tot'][iT,i])
                gplt['Mod C_G_Gross_Tot'][iS[j]]=np.mean(dM['C_G_Gross_Tot'][iT,i])
                gplt['Mod C_G_Net'][iS[j]]=np.mean(dM['C_G_Net'][iT,i])

    # Import TIPSY stemwood biomass
    iScn=0
    iGC=0
    gc=cbu.Import_BatchTIPSY_Output(meta,pNam,iScn,iGC)

    gplt['GY Csw t0']=np.nan*np.ones( gplt['Year t0'].size )
    gplt['GY Csw Net']=np.nan*np.ones( gplt['Year t0'].size )
    for iStand in range(meta['Geos']['Sparse']['X'].size):
        iS=np.where( (np.abs(gplt['X']-meta['Geos']['Sparse']['X'][iStand])<=100) & (np.abs(gplt['Y']-meta['Geos']['Sparse']['Y'][iStand])<=100) )[0]
        if iS.size==0:
            continue
        for j in range(iS.size):
            if (np.isnan(gplt['Age VRI t0'][iS[j]])==True) | (np.isnan(gplt['Delta t'][iS[j]])==True):
                continue
            iA=np.minimum(200.,gplt['Age VRI t0'][iS[j]]).astype(int)
            gplt['GY Csw t0'][iS[j]]=gc['Csw'][iA,iStand]
            if (np.isnan(gplt['Year t1'][iS[j]])==False) & (iA<198):
                iA2=np.arange(iA,np.minimum(200,iA+gplt['Delta t'][iS[j]]+1),1).astype('int')
                gplt['GY Csw Net'][iS[j]]=np.mean(gc['Gsw_Net'][iA2,iStand])

    # List of variables for analysis
    gplt['vaL']=['Age VRI t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0','Cdw t0','Ctot G Tot','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Net',
        'Mod A t0','Mod C_Biomass_Tot t0','Mod C_Biomass_Tot t1','Mod C_M_Tot','Mod C_G_Gross_Tot','Mod C_G_Net',
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

    gplt['soils']=soils

    # Ground plots
    gplt['soils']['vaL']=['TOT_C_THA','MIN_C_THA','ORG_C_THA','C_Soil_Tot','C_Soil_OHorizon']

    # Add simulations to soils data
    iT=np.where(tvSaved==2020)[0]
    gplt['soils']['C_Soil_Tot']=np.nan*np.ones(gplt['soils']['x'].size)
    gplt['soils']['C_Soil_OHorizon']=np.nan*np.ones(gplt['soils']['x'].size)
    cnt=0
    for i in range(meta['Geos']['Sparse']['X'].size):
        iS=np.where( (np.abs(gplt['soils']['x']-meta['Geos']['Sparse']['X'][i])<=2000) & (np.abs(gplt['soils']['y']-meta['Geos']['Sparse']['Y'][i])<=2000) )[0]
        if iS.size==0:
            continue
        gplt['soils']['C_Soil_Tot'][iS]=dM['C_Soil_Tot'][iT,i]
        gplt['soils']['C_Soil_OHorizon'][iS]=dM['C_Soil_OHorizon'][iT,i]
        cnt=cnt+1

    # Save
    gu.opickle(meta['Paths'][pNam]['Data'] + '\\Outputs\\EvaluationAtGroundPlots.pkl',gplt)

    return