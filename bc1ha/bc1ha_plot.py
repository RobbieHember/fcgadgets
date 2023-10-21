
#%% Import modules

import os
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fiona
import time
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcexplore.field_plots.Processing.psp_util as ugp
import fcgadgets.cbrunner.cbrun_util as cbu

#%%
def PlotVectorBaseMaps(meta,roi,ax):
    if meta['Graphics']['Map']['Show Bound Within']=='On':
        roi['gdf']['bound within'].plot(ax=ax,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
        roi['gdf']['bc_bound'].plot(ax=ax,edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    if meta['Graphics']['Map']['Show Lakes']=='On':
        roi['gdf']['lakes'].plot(ax=ax,facecolor=[0.82,0.88,1],label='Lakes',linewidth=0.25)
    if meta['Graphics']['Map']['Show Rivers']=='On':
        roi['gdf']['rivers'].plot(ax=ax,color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    if meta['Graphics']['Map']['Show Roads']=='On':
        roi['gdf']['road'].plot(ax=ax,edgecolor=[0.4,0,0],linewidth=0.25,label='Road',alpha=1,zorder=1)
    if meta['Graphics']['Map']['Show TPF and Cities']=='On':
        roi['gdf']['tpf'].plot(ax=ax,marker='s',edgecolor=[0.4,0,1],linewidth=0.25,label='TPF',alpha=1,zorder=1,markersize=15)
        if meta['Graphics']['Map']['Show Symbol Labels']=='On':
            for x,y,label in zip(roi['gdf']['tpf'].geometry.x,roi['gdf']['tpf'].geometry.y,roi['gdf']['tpf']['COMPANY_NAME']):
                ax.annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=6)
        roi['gdf']['cities'].plot(ax=ax,marker='o',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
        if meta['Graphics']['Map']['Show Symbol Labels']=='On':
            for x,y,label in zip(roi['gdf']['cities'].geometry.x,roi['gdf']['cities'].geometry.y,roi['gdf']['cities']['City Name']):
                ax.annotate(label,xy=(x,y),xytext=(4,3),textcoords="offset points",color=[0,0,0],fontsize=6)
    return roi

#%%
def Plot_LandCoverComp1(meta,roi,nam):
    z0=roi['grd'][nam]['Data']    
    lab=list(meta['LUT']['Derived']['lc_comp1'].keys())

    z1=len(lab)*np.ones(z0.shape,dtype='int8')
    for k in meta['LUT']['Derived']['lc_comp1'].keys():
        ind=np.where(z0==meta['LUT']['Derived']['lc_comp1'][k])
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
        else:
            z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]
    ind=np.where(roi['grd']['Data']==0)
    if ind[0].size>0:
        z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]+1
    else:
        z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.3,0.6,0,1),(0.65,1,0,1),(1,1,0.5,1),(0.8,0.7,1,1),(0.75,0.5,0.5,1),(0.75,0.75,0.75,1),(0.8,0,0,1),(0.95,0.95,0.95,1),(0.75,0.85,0.95,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    #ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
    
    return fig,ax

#%%
def Plot_LandUseComp1(meta,roi,nam):
    
    z0=roi['grd'][nam]['Data']
    lab=list(meta['LUT']['Derived']['lu_comp1'].keys())
    z1=len(lab)*np.ones(z0.shape,dtype='int8')
    for k in meta['LUT']['Derived']['lu_comp1'].keys():
        ind=np.where(z0==meta['LUT']['Derived']['lu_comp1'][k])
        z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]
    ind=np.where(roi['grd']['Data']==0)
    if ind[0].size>0:
        z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
    else:
        z1[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
    for k in meta['LUT']['Derived']['lu_comp1'].keys():
        z1[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
    
    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden
    cm=np.vstack( ((0.65,0,0,1),(1,0.95,0.5,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.25,0,0,1),(0.2,0.5,0,1),
                   (0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,0.75,0.2,1),(1,0.4,0.15,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
    return fig,ax

#%%
def Plot_LandUseComp1Panels(meta,roi):
    namL=['lu_comp1_2029s1','lu_comp1_2029s2','lu_comp1_2029s3','lu_comp1_2029s4']
    cnt=0
    plt.close('all'); fig,ax=plt.subplots(5,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    for nam in namL:
        z0=roi['grd'][nam]['Data']
        lab=list(meta['LUT']['Derived']['lu_comp1'].keys())
        z1=len(lab)*np.ones(z0.shape,dtype='int8')
        for k in meta['LUT']['Derived']['lu_comp1'].keys():
            ind=np.where(z0==meta['LUT']['Derived']['lu_comp1'][k])
            z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]
        ind=np.where(roi['grd']['Data']==0)
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['lu_comp1'][k]+1
        else:
            z1[-1,-1]=meta['LUT']['Derived']['lu_comp1'][k]+1
        for k in meta['LUT']['Derived']['lu_comp1'].keys():
            z1[0,meta['LUT']['Derived']['lu_comp1'][k]]=meta['LUT']['Derived']['lu_comp1'][k]
        
        N_vis=len(lab)
        N_hidden=1
        N_tot=N_vis+N_hidden
        cm=np.vstack( ((0.65,0,0,1),(0.9,0.8,0.4,1),(0.5,0.8,0.1,1),(0,0.2,0.5,1),(0.25,0,0,1),(0.2,0.5,0,1),(0.35,0.87,0.94,1),(0.85,0.95,0.7,1),(0.8,0.6,1,1),(1,1,0,1),(1,0,0,1),(0.4,0.7,0.79,1),(0.6,0.6,0.6,1),(1,1,1,1)) )
        cm=matplotlib.colors.ListedColormap(cm)    
        im=ax[cnt].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
        roi=PlotVectorBaseMaps(meta,roi,ax[0])
        ax[cnt].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
        ax[cnt].yaxis.set_ticks_position('both'); ax[cnt].xaxis.set_ticks_position('both'); ax[cnt].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[cnt].axis(meta['Graphics']['Map']['Map Axis Vis'])
        cnt=cnt+1
    gu.axletters(ax,plt,-0.025,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Skip=[4])

    ax[0].set(position=[-0.1,0.5,0.5,0.5])
    ax[1].set(position=[0.14,0.5,0.5,0.5])
    ax[2].set(position=[-0.1,0,0.5,0.5])
    ax[3].set(position=[0.14,0,0.5,0.5])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[4],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[4].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[0.5,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[4].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_LandUseComp1Panels','png',900)
    return fig,ax

#%%
def Plot_LandCoverLandUseChange(meta,roi,nam):
    #nam='luc1_early_type'
    z0=roi['grd'][nam]['Data']            
    lab=list(meta['LUT']['Derived']['lclu_chng_comp1'].keys())

    z1=(len(lab)+1)*np.ones(z0.shape,dtype='int8')
    for k in meta['LUT']['Derived']['lc_comp1'].keys():
        ind=np.where(z0==meta['LUT']['Derived']['lc_comp1'][k])
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]
        else:
            z1[0,meta['LUT']['Derived']['lc_comp1'][k]]=meta['LUT']['Derived']['lc_comp1'][k]
    ind=np.where(roi['grd']['Data']==0)
    if ind[0].size>0:
        z1[ind]=meta['LUT']['Derived']['lc_comp1'][k]+1
    else:
        z1[-1,-1]=meta['LUT']['Derived']['lc_comp1'][k]+1
    #z0[0,0:7]=np.arange(1,8)
    
    N_vis=len(lab)
    N_hidden=2
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.95,0.9,0.7,1),(0.75,0.9,0.55,1),(0.85,0.4,0.4,1),(0.4,0.75,1,1),(0.65,0.1,0,1),(0.4,0.6,0.8,1),(0.87,0.87,0.87,1),(0.15,0.15,0.15,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
    
    return fig,ax

#%%
def Plot_Geomorphons(meta,roi):
    
    lab=['Flat','Peak','Ridge','Shoulder','Spur','Slope','Hollow','Footslope','Valley','Pit']
    z1=roi['grd']['geomorph']['Data']
    ind=np.where(z1==0); z1[ind]=len(lab)+1
    ind=np.where(roi['grd']['Data']==0)
    if ind[0].size>0:
        z1[ind]=len(lab)+1
    
    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden
    cm=np.vstack( ((0.9,0.9,0.9,1),(0.25,0,0,1),
                   (0.75,0,0,1),(1,0.5,0.5,1),
                   (1,1,0.5,1),(0.65,1,0,1),
                   (0.65,0.87,0.24,1),(0.85,0.95,0.7,1),
                   (0.8,0.9,1,1),(0,0.25,0.5,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_Geomorphons','png',900)
    return fig,ax

#%%
def Plot_UplandWetlandForest(meta,roi,nam):

    z0=roi['grd'][nam]['Data']
    lab=list(meta['LUT']['Derived']['upwetf'].keys())
    z1=len(lab)*np.ones(z0.shape,dtype='int8')
    for k in meta['LUT']['Derived']['upwetf'].keys():
        ind=np.where(z0==meta['LUT']['Derived']['upwetf'][k])
        if ind[0].size>0:
            z1[ind]=meta['LUT']['Derived']['upwetf'][k]
        else:
            z1[0,meta['LUT']['Derived']['upwetf'][k]]=meta['LUT']['Derived']['upwetf'][k]
    ind=np.where(roi['grd']['Data']==0)
    if ind[0].size>0:
        z1[ind]=meta['LUT']['Derived']['upwetf'][k]+1
    else:
        z1[-1,-1]=meta['LUT']['Derived']['upwetf'][k]+1

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.3,0.6,0,1),(0.75,0.95,1,1),(0.2,0.65,1,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_' + nam,'png',900)
    plt.show()
    return fig,ax

#%%
def Plot_BurnSeverity(meta,roi):
    z0=roi['grd']['bsr_sc']['Data']    
    lab=list(meta['LUT']['Derived']['burnsev_comp1'].keys())
    z1=z0
    
    ind=np.where(roi['grd']['Data']==0)
    z1[ind]=6
    z1[0,0:6]=np.arange(1,7)

    N_vis=4
    N_hidden=2
    N_tot=N_vis+N_hidden
    
    lab.append('')
    
    cm=np.vstack( ((0.85,1,0.55,1),(1,0.9,0.65,1),(1,0.5,0.4,1),(0.5,0,0,1),(0.92,0.92,0.92,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    #ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_burnsev','png',900)
    
    return fig,ax

#%%
def Plot_HarvestYear(meta,roi,nam):
    z0=roi['grd']['harv_yr_' + nam]['Data']
    bw=10; bin=np.arange(1855,2025+bw,bw);
    
    N_vis=bin.size
    N_hidden=3
    N_tot=N_vis+N_hidden
    
    z1=N_vis*np.ones(z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==1) & (z0==0)]=i+2
    z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=i+1    
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3
    z1[0,0:N_tot]=np.arange(0,N_tot,1)

    lab=bin.astype(str)

    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.25,0.25,0.25,1),(0.75,0.75,0.75,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm,clim=(0,N_tot)) #
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis+1,1),ticks=np.arange(0.5,N_vis,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
        
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_harv_yr_' + nam,'png',900)
    return fig,ax

#%%

def Plot_HarvestRetentionComp1(meta,roi):

    lab=list(meta['LUT']['Derived']['harvret1'].keys())
    z1=roi['grd']['harvret1']['Data']
    z1[roi['grd']['Data']==0]=9
    z1[0,0:9]=np.arange(1,10)

    N_vis=8
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.6,1,0,1),(1,0.5,0,1),(0.5,0.8,1,1),(0,0.5,0,1),(0.75,0,1,1),(0.4,0.4,0.4,1),(0.6,0.6,0.6,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_harvest_retention1','png',900)

    return fig,ax

#%%
def Plot_REARs(meta,roi):
    if 'rears' in roi['grd'].keys():
        z1=roi['grd']['rears']['Data']
    else:
        z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\REARs_CompPlusProp.tif')['Data']
        z1=gis.UpdateGridCellsize(z1,meta['Graphics']['Map']['RGSF'])['Data']
    z1[(roi['grd']['Data']==1) & (z1==0) & (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])]=3
    z1[(roi['grd']['Data']==1) & (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=4
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=5

    lab=['Protected','Protected (proposed)','Unprotected forest','Non-forest land']

    N_vis=4
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.7,0.6,1,1),(0.25,0.85,0.9,1),(0,0.4,0,1),(0.83,0.86,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_REARs','png',900)

    return fig,ax

#%%
def Plot_Harvest_Po(meta,roi):
    z=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\HarvestProbability.tif')
    z=gis.UpdateGridCellsize(z,meta['Graphics']['Map']['RGSF'])
    sf=1000
    # Apply scale factor to convert to (%/yr)
    #z['Data']=z['Data'].astype('float')/1000
    bw=0.1; bin=np.arange(0,0.5+bw,bw)
    z1=np.ones( z['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z['Data']-bin[i]*sf)<=bw*sf/2)
        if ind[0].size==0:
            z1[ind[0][i],ind[1][i]]=i
        else:
            z1[ind]=i
    ind=np.where( z['Data']>=bin[i]*sf ); z1[ind]=i
    ind=np.where( roi['grd']['Data']!=1 ); z1[ind]=i+1

    lab=["%.2f" % x for x in bin]
    lab=np.append(lab,np.array(['Water']))

    N_vis=bin.size
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=plt.cm.get_cmap('viridis',bin.size)
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(18,18*0.5*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
        roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    ax[0].set(position=[0,0,0.5,1],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[0.39,0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],0.02,N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    # Model
    V_Merch=np.arange(1,1200)
    beta=[0.005,-0.04,400]
    Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    ax[2].plot(V_Merch,Po*100,'k-',linewidth=1,label='Harvest on-the-fly model 1')
    beta=[0.005,-0.04,500]
    Po=beta[0]*(1/(1+np.exp(beta[1]*(V_Merch-beta[2]))))
    ax[2].plot(V_Merch,Po*100,'g--',linewidth=1,label='Harvest on-the-fly model 2')
    ax[2].set(position=[0.55,0.12,0.4,0.86],xticks=np.arange(0,1300,100),xlabel='Merchantable volume (m$^3$ ha$^-$$^1$)',ylabel='Annual probability of harvest (%)',xlim=[0,800],ylim=[0,0.7])
    ax[2].legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')
    ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])

    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\taz_ann_prob_harvest','png',900)
    plt.show()
    return

#%%
def Plot_BGC_Zone(meta,roi):

    z0=roi['grd']['bgcz']['Data']
    lab0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys())
    cl0=np.column_stack([meta['LUT']['Raw']['bgc_zone']['R'],meta['LUT']['Raw']['bgc_zone']['G'],meta['LUT']['Raw']['bgc_zone']['B']])
    id0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].values())
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)

    N_vis=int(np.max(z1))
    N_hidden=2
    N_tot=N_vis+N_hidden
    ind=np.where(z1==0); z1[ind]=N_vis+1
    lab1=lab1[1:]
    cl1=cl1[1:,:]

    z1[0,0]=N_vis+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=N_vis+2

    lab1=np.append(lab1,[''])

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_vis);
    for i in range(N_vis):
        cm.colors[i,0:3]=cl1[i,:]
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab1)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_bgcz','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_Ownership(meta,roi):

    z0=roi['grd']['own']['Data']
    lab0=list(meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION'].keys())
    #cl0=np.column_stack([meta['LUT']['Raw']['bgc_zone']['R'],meta['LUT']['Raw']['bgc_zone']['G'],meta['LUT']['Raw']['bgc_zone']['B']])
    #cl0=np.random.random((len(lab0),3))
    cl0=np.array([[0.95,0.9,1],[0.9,0.85,1],[0.85,0.8,1],[0.2,0.5,0],[0.75,0.7,1],[0.7,0.65,1],[0.3,0.6,0],[0.6,0.55,1],[0.55,0.5,1],
                  [1,1,0.8],[0.95,0.95,0.7],[0.95,0.95,0.6],[0.95,0.95,0.5],[0.95,0.95,0.4],[0.9,0.9,0.3],[0.9,0.9,0.2],[0.85,0.85,0.1],                  
                  [0,0.75,0],[0,0.5,0],[0,0.25,0],[0.85,0.65,0.75],[0.7,0.5,0.6],[0.6,0.8,0.8],[0.7,0.9,0.8],[0.8,1,0.9],[0.9,1,0.65],[0.9,1,0.6],
                  [0.85,0.95,1],[0.8,0.9,1],[0.75,0.85,1],[0.7,0.8,1],[0.65,0.75,1],
                  [1,0.9,0.8],[1,0.8,0.7],[1,0.7,0.6],[0.9,0.9,0.9]])
    id0=list(meta['LUT']['F_OWN']['OWNERSHIP_DESCRIPTION'].values())
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)

    N_vis=int(np.max(z1))
    N_hidden=2
    N_tot=N_vis+N_hidden
    ind=np.where(z1==0); z1[ind]=N_vis+1
    lab1=lab1[1:]
    cl1=cl1[1:,:]

    z1[0,0]=N_vis+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=N_vis+2

    lab1=np.append(lab1,[''])

    cm=plt.cm.get_cmap('viridis',N_vis);
    for i in range(N_vis):
        cm.colors[i,0:3]=cl1[i,:]
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.01,0.03,0.98])
    cb.ax.set(yticklabels=lab1)
    cb.ax.tick_params(labelsize=4,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[0.66,0.45,0.01,0.54]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_ownership','png',900)
    return fig,ax

#%%
def Plot_SalvageLogging(meta,roi):

    lab=['Salvage (>50% dead)','Salvage (>10% dead)','Harvest (<10% dead)','No harvesting','Non-forest land','Outside']
    lab=lab[0:-1]
    z1=roi['grd']['harv_salv']['Data']
    z1[roi['grd']['Data']==0]=6

    # Number of colours and number of colours excluded from colorbar
    N_vis=4
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.9,0,0,1),(1,0.85,0,1),(0.85,1,0.65,1),(0.6,0.6,0.6,1), (0.93,0.93,0.93,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)

    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_harv_salv','png',900)

    return fig,ax

#%%
def Plot_MAP(meta,roi):

    z0=roi['grd']['prcp_ann_n']['Data']

    bw=200; bin=np.arange(0,2400+bw,bw)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

    for i in range(N_vis):
        z1[0,i]=i+1

    lab=['']*(N_tot-1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_prcp_ann_n','png',900)

    return fig,ax

#%%
def Plot_MAT(meta,roi):

    z0=roi['grd']['tmean_ann_n']['Data']

    bw=10; bin=np.arange(-30,120+bw,bw)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i+1
    ind=np.where(z0>=bin[i]); z1[ind]=i+1
    z1[1,1]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+3

    for i in range(N_vis):
        z1[0,i]=i+1

    lab=['']*(N_tot-1)
    lab[0:N_vis]=np.array(bin/10).astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_tmean_ann_n','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_Fire2023(meta,roi):

    z1=4*np.ones(roi['grd']['Data'].shape,dtype='int8')

    ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==1) ); z1[ind]=1
    ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==1) & (roi['grd']['fire_2023']['Data']>0) ); z1[ind]=2
    ind=np.where( (roi['grd']['lc_comp1_2019']['Data']!=1) ); z1[ind]=3
    ind=np.where( (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water']) ); z1[ind]=4

    lab=['Forest land','Forest land (affected)','Non-forest land']

    N_vis=len(lab)
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.vstack( ((0.8,0.8,0.8,1),(0.85,0,0,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    gdf['cities'][gdf['cities']['Territory']=='BC'].plot(ax=ax[0],marker='s',edgecolor=[0,0,0],facecolor=[1,1,0],lw=0.25,markersize=9,alpha=1,zorder=2)
    for x,y,label in zip(gdf['cities'].geometry.x,gdf['cities'].geometry.y,gdf['cities'].Name):
        ax[0].annotate(label,xy=(x,y),xytext=(3,2),textcoords="offset points",color=[0,0,0],fontsize=4)

    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    if meta['Graphics']['Print Figures']=='On':
        gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_fire2023','png',900)

    return fig,ax

#%% Plot major timber producing facilities
def PlotMills(ax,gdf,mtypeL,labels):
    lw=0.75
    for mtype in mtypeL:
        #mtype='PLP'
        ind=np.where( (gdf['PRODUCT_CODE']==mtype) )[0]
        if mtype=='LBR':
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_BOARD_FT']/453
            ms=500*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='g',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif (mtype=='PLY') | (mtype=='VNR') | (mtype=='OSB') | (mtype=='PNL'):
            # One PNL (panel) mill WestPine MDF in Quesnell - it is MDF
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_SQ_FT']/885
            ms=300*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='c',facecolor='c',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLT':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=1.0*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='r',facecolor='r',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='CHP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_BDUS']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='y',facecolor='y',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='PLP':
            y=gdf.iloc[ind]['EST_AN_CAP_000_TONNES']*2
            ms=0.75*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='k',lw=lw,markersize=ms,alpha=0.35,zorder=2)
        elif mtype=='LVL':
            # Laminated veneer lumber
            y=gdf.iloc[ind]['EST_AN_CAP_MLN_CUBIC_FT']#/885
            ms=200*y
            gdf.iloc[ind].plot(ax=ax,marker='o',edgecolor='k',facecolor='g',lw=lw,markersize=ms,alpha=0.35,zorder=2)

        if labels=='On':
            for x,y,label in zip(gdf.iloc[ind].geometry.x,gdf.iloc[ind].geometry.y,gdf.iloc[ind].COMPANY_NAME):
                ax.annotate(label,xy=(x,y),xytext=(5,4),textcoords="offset points")
    return ax

#%% Plot FECA year

def Plot_FECA_Year(meta,roi):

    z0=roi['grd']['feca_yr']['Data']
    bw=5; bin=np.arange(1975,2025+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(z0==0)]=i+1
    z1[(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2

    N_vis=bin.size+3
    N_hidden=3
    N_tot=N_vis+N_hidden

    lab=bin.astype(str)

    #cm=plt.cm.get_cmap('viridis',N_vis)
    cm=plt.cm.get_cmap('plasma',N_vis)
    cm=np.vstack( (cm.colors,(0,0,0,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_feca_yr','png',900)
    return fig,ax

#%%
def Plot_Spc1(meta,roi):

    # Compress categories
    z0=roi['grd']['spc1_ntems']['Data']
    lab0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].keys()))
    id0=np.array(list(meta['LUT']['Derived']['spc1_ntems'].values()))
    cl0=np.column_stack([meta['LUT']['Raw']['spc1_ntems']['R'],meta['LUT']['Raw']['spc1_ntems']['G'],meta['LUT']['Raw']['spc1_ntems']['B']])
    z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
    z1=z1+1

    n=int(np.max(z1))
    z1[0,0]=n+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=n+2

    lab1=np.append(lab1,[''])

    N_vis=n
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_vis);
    for i in range(N_vis):
        cm.colors[i,0:3]=cl1[i,:]
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab1)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_spc1_ntems','png',900)
    plt.show()
    return fig,ax

#%%
def Plot_PROJ_AGE_1(meta,roi):

    z0=roi['grd']['PROJ_AGE_1']['Data']

    bw=20
    bin=np.arange(0,220,bw)
    id=np.arange(1,bin.size+1)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(roi['grd']['PROJ_AGE_1']['Data']>bin[i]); z1[ind]=id[-1]
    z1[(roi['grd']['PROJ_AGE_1']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_PROJ_AGE_1','png',900)

    return fig,ax

#%%
def Plot_SI(meta,roi):

    #z0=roi['grd']['SITE_INDEX']['Data']
    #z0=roi['grd']['si_spl_fd']['Data']
    #z0=roi['grd']['si_spl_hw']['Data']
    z0=roi['grd']['si_spl_ntems']['Data']

    bw=2
    bin=np.arange(10,40+bw,bw)
    id=np.arange(1,bin.size+1)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(z0>bin[i]); z1[ind]=id[-1]
    z1[(z0==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
    #z1[1,1]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

    lab=['']*(N_vis+1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_SITE_INDEX','png',900)
    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_si_spl_hw','png',900)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_si_spl_ntems','png',900)

    return fig,ax

#%%
def Plot_Age_NTEMS(meta,roi):

    z0=roi['grd']['age_ntem']['Data']

    bw=20
    bin=np.arange(0,220,bw)
    id=np.arange(1,bin.size+1)

    N_vis=bin.size
    N_hidden=2
    N_tot=N_vis+N_hidden

    z1=N_vis*np.ones(z0.shape)
    for i in range(N_vis):
        ind=np.where(np.abs(z0-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=id[i]
        else:
            z1[0,i]=id[i]
    ind=np.where(z0>bin[i]); z1[ind]=id[-1]
    z1[(z0==0) | (roi['grd']['lc_comp1_2019']['Data']!=meta['LUT']['Derived']['lc_comp1']['Forest'])]=id[-1]+1
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=id[-1]+2

    lab=['']*(N_tot-1)
    lab[0:N_vis]=bin.astype(str)

    cm=plt.cm.get_cmap('viridis',N_vis)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])

    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_Age_NTEM','png',900)

    return fig,ax

#%%

def Plot_PlantedMask(meta,roi):

    lab=[]
    z1=2*np.ones(roi['grd']['plam']['Data'].shape)
    z1[(roi['grd']['Data']==1) & (roi['grd']['plam']['Data']>0)]=1; lab.append('Planted')
    z1[(roi['grd']['Data']==1) & (roi['grd']['plam']['Data']==0)]=2; lab.append('Not planted')
    z1[(roi['grd']['Data']!=1)]=3;

    # Number of colours and number of colours excluded from colorbar
    N_vis=2
    N_hidden=1
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0,0.6,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_planted','png',900)

    return fig,ax

#%%
def Plot_Elev(meta,roi):

    z=roi['grd']['elev']['Data']
    if (roi['Type']=='ByTSA') | (roi['Type']=='ByRegDis'):
        z[roi['grd']['Data']==0]=0

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['elev']['Data'],extent=roi['grd']['Extent'],cmap='Greys')
    roi=PlotVectorBaseMaps(meta,roi,ax[0])

    # Plot fuel treatments
    #fre['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0.5,0],linewidth=1.5,label='Opening',alpha=1)

    # Plot lumber mills, pulp mills and chipper mills
    #mtypeL=['LBR','PLP','PLT']
    #mtypeL=list(gdf['tpf']['gdf']['PRODUCT_CODE'].unique())
    #ax=u1ha.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

    roi['gdf']['cities'].plot(ax=ax[0],marker='s',edgecolor=[0.75,0.3,0],facecolor=[1,0.6,0],lw=1,markersize=20,alpha=1,zorder=2)
    for x,y,label in zip(roi['gdf']['cities'].geometry.x,roi['gdf']['cities'].geometry.y,roi['gdf']['cities'].Name):
        ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0.9,0.45,0],fontsize=4)

    #roi['gdf']['popp'].plot(ax=ax[0],marker='s',edgecolor=None,facecolor=[0.9,0.45,0],lw=0.5,markersize=3,alpha=1,zorder=2)
    #for x,y,label in zip(roi['gdf']['popp'].geometry.x,roi['gdf']['popp'].geometry.y,roi['gdf']['popp'].NAME):
    #    ax[0].annotate(label,xy=(x,y),xytext=(2,1.5),textcoords="offset points",color=[0.9,0.45,0],fontsize=4)

    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos=[meta['Graphics']['Map']['Legend X'],0.5,meta['Graphics']['Map']['Legend Width'],0.4]
    ax[1].set(position=pos);

    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_elev','png',900)
    return fig,ax

#%%
def Plot_LC20_CEC(meta,roi):

    dCEC=gu.ReadExcel(meta['Paths']['Model']['Parameters'] + '\\LUT_lc_cec_Compressed.xlsx')

    N_vis=dCEC['Name'].size
    N_hidden=1
    N_tot=N_vis+N_hidden

    z1=roi['grd']['lc_cec_2020']['Data']
    z1[(roi['grd']['Data']==0)]=np.max(z1)
    z1[0,0:N_vis]=np.arange(1,N_vis+1,1)
    
    # Labels
    lab=list(meta['LUT']['Derived']['lc_cec_c'].keys())
    #lab.append('')
    #lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_vis);
    for i in range(N_vis):
        cm.colors[i,0:3]=[ dCEC['R'][i]/255,dCEC['G'][i]/255,dCEC['B'][i]/255 ]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_lc_cec_2020','png',900)

    return fig,ax

#%%
def Plot_TreeDensityClass(meta,roi):

    lab=['Treed sparse','Treed open','Treed dense','Shrubland','Grassland','Non-treed']
    z1=roi['grd']['tdc_wsg']['Data']
    #z1[roi['grd']['Data']==0]=np.max(z1)
    z1[0,0:6]=np.arange(1,7,1)

    N_vis=6
    N_hidden=1
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.75,0.95,0.6,1),(0,0.7,0,1),(0,0.25,0,1),(0.8,0.6,0.6,1),(1,1,0.5,1),(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_tdc_wsg','png',900)

    return fig,ax

#%%
def Plot_SI(meta,roi):

    z=roi['grd']['SITE_INDEX']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['SITE_INDEX']['Data'],extent=roi['grd']['Extent'],cmap='magma',clim=[5,22])

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].grid(False)

    # # Add relief shading
    # z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    # ls=LightSource(azdeg=90,altdeg=45)
    # ve=0.1
    # hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    # ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

#%%
def Plot_PFI(meta,roi):

    # Grid
    bw=20; bin=np.arange(0,140+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['pfi_c']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['pfi_c']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['pfi_c']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('cividis',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

#%%

def Plot_PFI(meta,roi):

    zB=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')
    zB['Data']=0.5*0.5*zB['Data']
    zB['Data']=zB['Data'].astype('int16')

    # Grid
    bw=10; bin=np.arange(0,160+bw,bw);
    z1=(bin.size)*np.ones( zB['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( zB['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    ind=np.where(zB['Data']<0); z1[ind]=0
    ind=np.where(zB['Data']>bin[i]); z1[ind]=i
    z1[0,0]=i+1
    z1[1,1]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('cividis',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);
    #ind=np.where(roi['gdf']['cc']['gdf']['HARVEST_YEAR']==1995)
    roi['gdf']['cc']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.25,0.75,1],linewidth=1,label='Opening')
    roi['gdf']['fcres']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.65,1,0.1],linewidth=1,linestyle='--',label='Opening')
    #for idx, row in roi['gdf']['cc']['gdf'].iterrows():
        #print(row.keys())
        #break
    #    plt.text(row['geometry'].centroid.x,row['geometry'].centroid.y,str(row['HARVEST_YEAR']),ha='center',fontsize=5,color=[0.25,0.75,1])
    return fig,ax

#%%

def Plot_SoilOrganicCarbon_GSOC(meta,roi):

    # Grid
    bw=10; bin=np.arange(0,100+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['gsoc']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['gsoc']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['gsoc']['Data']>=bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['gsoc']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gsoc','png',900)

    return fig,ax

#%%

def Plot_CrownCover(meta,roi):

    # Grid
    bw=10; bin=np.arange(0,100+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['crownc']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['crownc']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['crownc']['Data']>bin[i])
    z1[ind]=0
    z1[(roi['grd']['Data']==1) & ( roi['grd']['crownc']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_crownc','png',900)

    return fig,ax

#%%
def Plot_SPH_Live(meta,roi):

    # Grid
    bw=400; bin=np.arange(0,2400,bw);
    z1=(bin.size)*np.ones( roi['grd']['sphlive']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['sphlive']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(roi['grd']['sphlive']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['sphlive']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['ax2 pos long']);

    # cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    # cb.ax.set(yticklabels=lab)
    # cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    # for i in range(0,N_vis):
    #     ax[1].plot([0,100],[i/(N_vis-N_hidden),i/(N_vis-N_hidden)],'k-',linewidth=0.5)
    # ax[1].set(position=meta['Graphics']['ax2 pos long']);
    # pos2=copy.copy(meta['Graphics']['pos2'])
    # pos2[1]=0.6
    # pos2[3]=0.24
    # ax[1].set(position=pos2)

    return fig,ax

#%%

def Plot_LUC_Year(meta,roi):

    z0=roi['grd']['luc1_yr']['Data']

    bw=10; bin=np.arange(1850,2020+bw,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(np.abs( z0-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
    #z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_luc1_yr','png',900)
    return fig,ax

#%%
def Plot_WildfireYear(meta,roi):

    # Grid
    bw=10; bin=np.arange(1910,2020+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['fire_yr']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['fire_yr']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
    #z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_fire_yr','png',900)

    return fig,ax

#%%

def Plot_GFC_LossYear(meta,roi):

    z0=roi['grd']['gfcly']['Data']
    #z0=roi['grd']['gfcly_filt']['Data']

    bw=1; bin=np.arange(2001,2022,bw);
    z1=(bin.size)*np.ones( z0.shape)
    for i in range(bin.size):
        ind=np.where(z0==bin[i])
        z1[ind]=i
    z1[(z0==0)]=i+1
    #z1[(roi['grd']['Data']==0)]=i+2
    z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
    L=i+2

    lab=np.array(bin).astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.85,0.85,0.85,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)

    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gfcly','png',900)
    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_gfcly_filt','png',900)

    return fig,ax

#%%
def Plot_SoilWaterContent(meta,roi):

    bw=20; bin=np.arange(0,200+bw,bw);
    z1=bin.size*np.ones( roi['grd']['ws_gs_n']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['ws_gs_n']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    z1[(roi['grd']['ws_gs_n']['Data']>200)]=i
    z1[(roi['grd']['lcc1']==0) | (roi['grd']['lcc1']==meta['LUT']['Derived']['lc_comp1']['Water']) | (roi['grd']['Data']!=1)]=i+1

    lab=bin.astype(str)
    lab=np.append(lab,'Hidden')

    N_vis=bin.size
    N_hidden=1
    N_tot=N_vis+N_hidden

    cm=np.column_stack((roi['grd']['ws_gs_n']['cm']['cl1'],roi['grd']['ws_gs_n']['cm']['cl2'],roi['grd']['ws_gs_n']['cm']['cl3'],np.ones(roi['grd']['ws_gs_n']['cm']['bin'].size)))
    cm=np.vstack( (cm,(1,1,1,1)) ) # (0.83137,0.81569,0.78431,1)
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1);
    zmx=np.max(z1);
    cb_ivl=(zmx-zmn)/N_tot;
    cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_SoilWaterContent_Normal','png',300)

    return fig,ax

#%% Empty forest mask
# *** For use with custom vector layers ***
def Plot_ForestMask(meta,roi):
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2    
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)    
    ms=2
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    # *** Add Vector layers ***
    
    # points=[]
    # for k in range(x.size):
    #     points.append(Point(x[k],y[k]))
    # gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    # gdf1.crs=roi['crs']
    # gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    # gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
        
    ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
    #gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_GroundPlots','png',900)
    return fig,ax

#%%

def Plot_GroundPlots(meta,roi):
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2
    
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    ms=2
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
       
    ind=np.where(z1==2)
    A=ind[0].size/1e3
    
    ivl=32
    ind=np.where(z1[0::ivl,0::ivl]==2)
    x=roi['grd']['X'][0::ivl,0::ivl][ind]
    y=roi['grd']['Y'][0::ivl,0::ivl][ind]
    rho1=x.size/A
    print(A)
    print(rho1)
    
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
        
    ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_GroundPlots','png',900)
    return

#%%
def Plot_GroundPlotsPanels(meta,roi):
    import fcexplore.field_plots.Processing.psp_util as ugp
    meta2,gpt=ugp.ImportGroundPlotData(meta,type='Stand')
    
    plt.close('all'); fig,ax=plt.subplots(1,4,figsize=gu.cm2inch(17.5,15))
    
    # # BGC Zone
    
    # z0=roi['grd']['bgcz']['Data']
    # lab0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys())
    # cl0=np.column_stack([meta['LUT']['Raw']['bgc_zone']['R'],meta['LUT']['Raw']['bgc_zone']['G'],meta['LUT']['Raw']['bgc_zone']['B']])
    # id0=list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].values())
    # z1,lab1,cl1=gis.CompressCats(z0,id0,lab0,cl0)
    
    # N_vis=int(np.max(z1))
    # N_hidden=2
    # N_tot=N_vis+N_hidden
    # ind=np.where(z1==0); z1[ind]=N_vis+1
    # lab1=lab1[1:]
    # cl1=cl1[1:,:]
    
    # z1[0,0]=N_vis+1
    # z1[(roi['grd']['Data']==0) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=N_vis+2
    
    # lab1=np.append(lab1,[''])
    
    # cm=plt.cm.get_cmap('viridis',N_vis);
    # for i in range(N_vis):
    #     cm.colors[i,0:3]=cl1[i,:]
    # cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    # cm=matplotlib.colors.ListedColormap(cm)
    
    # im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    # if meta['Graphics']['Map']['Show Bound Within']=='On':
    #     roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    # if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
    #     roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    # if meta['Graphics']['Map']['Show Lakes']=='On':
    #     roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    # if meta['Graphics']['Map']['Show Rivers']=='On':
    #     roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    # if meta['Graphics']['Map']['Show Roads']=='On':
    #     roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    
    # ax[0].set(position=[0,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    # ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    # zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    # cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    # cb=plt.colorbar(im,cax=ax[4],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    # ax[4].set(position=[0.71,0.6,0.05,0.14])
    # cb.ax.set(yticklabels=lab1)
    # cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    # cb.outline.set_edgecolor('w')
    # for i in range(cb_bnd.size):
    #     ax[4].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
    # ax[4].set(position=[0.4,0.68,0.02,0.3])
      
    # VRI plots
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    ms=2
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=[0,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    ind=np.where( (gpt['PTF CNV']==1) & (gpt['PTF CN']==0) )[0]
    x=gpt['X'][ind]; y=gpt['Y'][ind]
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[0],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
    
    # CMI plots
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    ms=2
    im=ax[1].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[1].set(position=[0.5,0.5,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[1].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    ind=np.where( (gpt['PTF CN']==1) & (np.isnan(gpt['Year t1'])==True) )[0]
    x=gpt['X'][ind]; y=gpt['Y'][ind]
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[1],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Single measurement')
    
    ind=np.where( (gpt['PTF CN']==1) & (np.isnan(gpt['Year t1'])==False) )[0]
    x=gpt['X'][ind]; y=gpt['Y'][ind]
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[1],marker='s',markersize=0.7,facecolor=[0.25,0.5,0.75],edgecolor=None,linewidth=0.75,alpha=1,label='With remeasurement')
    ax[1].legend(loc='lower left',facecolor=[1,1,1],frameon=False)
          
    # YSM plots
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    ms=2
    im=ax[2].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[2].set(position=[0,0,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[2].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    ind=np.where( (gpt['PTF YSM']==1) & (np.isnan(gpt['Year t1'])==True) )[0]
    x=gpt['X'][ind]; y=gpt['Y'][ind]
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[2],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
    
    ind=np.where( (gpt['PTF YSM']==1) & (np.isnan(gpt['Year t1'])==False) )[0]
    x=gpt['X'][ind]; y=gpt['Y'][ind]
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[2],marker='s',markersize=0.7,facecolor=[0.25,0.5,0.75],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
    
    # SOC    
    z1=roi['grd']['Data'].copy()
    ind=np.where(roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Forest'])
    z1[ind]=2
    cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    ms=2
    im=ax[3].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[3].set(position=[0.5,0,0.5,0.5],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[3].yaxis.set_ticks_position('both'); ax[3].xaxis.set_ticks_position('both'); ax[3].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[3].axis(meta['Graphics']['Map']['Map Axis Vis'])
    
    soils=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')
    x=soils['x']; y=soils['y']
    points=[]
    for k in range(x.size):
        points.append(Point(x[k],y[k]))
    gdf1=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
    gdf1.crs=roi['crs']
    gdf1=gpd.overlay(gdf1,roi['gdf']['bound within'],how='intersection')
    gdf1.plot(ax=ax[3],markersize=0.7,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')
    
    gu.axletters(ax,plt,0.0,0.9,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')
    
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_GroundPlotsFull2','png',900)
    return

#%%

def Plot_RangeTenure(meta,roi):

    lab=['Forest with grazing tenure','Forest with haycutting tenure','Forest with no range tenure','Non-forest land','Non land']
    lab=lab[0:-1]
    z1=roi['grd']['rangecon']['Data']
    z1[roi['grd']['Data']==0]=5

    # Number of colours and number of colours excluded from colorbar
    N_vis=4
    N_hidden=1
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.85,1,0.5,1),(0,0.5,0,1),(0.6,0.6,0.6,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_range_consol','png',900)

    return fig,ax

#%%

def Plot_biomass_glob(meta,roi):

    # Grid
    bw=20; bin=np.arange(0,200+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['biomass_glob']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['biomass_glob']['Data']-bin[i])<=bw/2)
        if ind[0].size>0:
            z1[ind]=i
        else:
            z1[0,i]=i
    ind=np.where(roi['grd']['biomass_glob']['Data']>bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['biomass_glob']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_vis=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi=PlotVectorBaseMaps(meta,roi,ax[0])
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_vis):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_biomass_glob','png',900)

    return fig,ax

#%%
def Plot_DistanceFromRoad(meta,roi):

    #z=roi['grd']['d2fac']['Data']
    z=roi['grd']['d2road']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    N_vis=10
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z,extent=roi['grd']['Extent'],cmap='viridis')

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    # Add relief shading
    #z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    #ls=LightSource(azdeg=90,altdeg=45)
    #ve=0.1
    #hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    #ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_vis-(N_hidden-1),1),ticks=np.arange(0.5,N_vis+1.5,1))
    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)

    return fig,ax

#%%
def Plot_LU_Change_FromCEC(meta,roi,type):

    #type='def'
    if type=='def':
        lab=['Forest to Settlement','Forest to Cropland','Forest to Barren Ground','Land']# Deforestation
    elif type=='aff':
        lab=['Settlement to Forest','Cropland to Forest','Barren Ground to Forest','Land'] # Afforestation
    else:
        pass

    z0=roi['grd']['luc_' + type + '_cec']['Data']
    z1=z0
    z1[(z0==0)]=4
    z1[(roi['grd']['Data']!=1) | (roi['grd']['lc_comp1_2019']['Data']==meta['LUT']['Derived']['lc_comp1']['Water'])]=5

    # Number of colours and number of colours excluded from colorbar
    N_vis=3
    N_hidden=2
    N_tot=N_vis+N_hidden

    # Colormap
    cm=np.vstack( ((0.5,0.0,0.0,1),(1.0,0,0,1),(1.0,0.5,0,1),(0.83,0.83,0.83,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,extent=roi['grd']['Extent'],cmap=cm)
    if meta['Graphics']['Map']['Show Bound Within']=='On':
        roi['gdf']['bound within'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    if meta['Graphics']['Map']['Show Bound Land Mask']=='On':
        roi['gdf']['bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    if meta['Graphics']['Map']['Show Lakes']=='On':
        roi['gdf']['lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],label='Lakes',linewidth=2.25)
    if meta['Graphics']['Map']['Show Rivers']=='On':
        roi['gdf']['rivers'].plot(ax=ax[0],color=[0.6,0.8,1],label='Rivers',linewidth=0.25)
    if meta['Graphics']['Map']['Show Roads']=='On':
        roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.25,label='Road',alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='')
    ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

    zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
    cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
    cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
    ax[1].set(position=[0.71,0.6,0.05,0.14])
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=meta['Graphics']['Map']['Legend Font Size'],length=0)
    cb.outline.set_edgecolor('w')
    for i in range(cb_bnd.size):
        ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)

    pos2=[meta['Graphics']['Map']['Legend X'],0.99-N_vis*meta['Graphics']['Map']['Legend Text Space'],meta['Graphics']['Map']['Legend Width'],N_vis*meta['Graphics']['Map']['Legend Text Space']]
    ax[1].set(position=pos2)
    gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\' + roi['Name'] + '_lu_change_' + type + '_cec','png',900)

    return fig,ax

#%% Forest health outbreak maps

def Plot_OutbreakMap(meta,roi):
    pest_cd='IBS'
    tv=np.arange(1999,2022,1)
    zG=np.zeros(zRef['Data'].shape,dtype='int8')
    for iT in range(tv.size):
    
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_' + pest_cd + '_SeverityClass_' + str(tv[iT]) + '.tif')
    
        zG[(z['Data']>1)]=1
    
        # Grid and labels
        z1=6*np.ones(zRef['Data'].shape)
        z1[(z['Data']==1)]=0
        z1[(z['Data']==2)]=1
        z1[(z['Data']==3)]=2
        z1[(z['Data']==4)]=3
        z1[(z['Data']==5)]=4
        z1[(z['Data']==0) & (zRef['Data']==1)]=6
        z1[(z['Data']==0) & (zG==1)]=5
        z1[(zRef['Data']==0)]=7
    
        ind=np.where(z1==7)
        for i in range(7):
            z1[ind[0][i],ind[1][i]]=i
        #z1[ind[0][1],ind[1][1]]=4
    
        lab=['Trace','Light','Moderate','Severe','Very severe','Previously affected','Unaffected','N/A']
    
        # Number of colours and number of colours excluded from colorbar
        N_vis=7
        N_hidden=1
        N_tot=N_vis+N_hidden
    
        # Colormap
        cm=np.vstack( ((1,1,0,1),(1,0.5,0,1),(1,0,0,1),(0.7,0,0,1),(0.35,0,0,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
        cm=matplotlib.colors.ListedColormap(cm)
    
        plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
        im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm)
        gdf['bc_bound']['gdf'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
        ax[0].set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
        ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both')
        ax[0].grid(meta['Graphics']['ax1 gridvis'])
        ax[0].axis(meta['Graphics']['ax1 vis'])
    
        cb_ivl=(N_tot-1)/N_tot
        cbivls=np.arange( cb_ivl , N_tot , cb_ivl)
        cbivls_low=np.arange( 0 , N_tot , cb_ivl)
    
        cb=plt.colorbar(im,cax=ax[1],cmap=cm,
                        boundaries=np.arange(0,cbivls_low[N_vis],cb_ivl),
                        ticks=np.arange(cb_ivl/2,N_tot-1,cb_ivl) )
        ax[1].set(position=[0.71,0.65,0.045,0.2])
    
        cb.ax.set(yticklabels=lab)
        cb.ax.tick_params(labelsize=11,length=0)
        cb.outline.set_edgecolor('w')
        for i in range(cbivls.size):
            ax[1].plot([0,100],[cbivls[i],cbivls[i]],'w-',linewidth=2)
    
        gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Beetles\Maps\Map_' + pest_cd + '_Severity_' + str(tv[iT]),'png',300)
    return