
#%% Import modules

import os
import numpy as np
import gc
import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import copy
import fiona
import time
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.query_vector_db as qv
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1hau

#%% Path management

meta={}
meta['Paths']={}
meta['Paths']['BC1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
meta['Paths']['Forest Inventory Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930\Disturbances.gdb'
#meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation'

#%% Plotting parameters

meta['Graphics']={}
meta['Graphics']['figsize1']=[900,700]
meta['Graphics']['pos1']=[0.04,0.02,0.74,0.95]
meta['Graphics']['pos2']=[0.79,0.6,0.03,0.35]
meta['Graphics']['pos2_long']=[0.79,0.5,0.03,0.45]

params_graphic=cbu.Import_GraphicsParameters('bc1ha_1')
plt.rcParams.update(params_graphic)

#%% Import base maps

bm,tsa,road,district=bc1hau.Import_BaseMaps()

# Load dataset with CRS
meta['crs']=bm['gdf_bm'].crs

#%% Define region of interest

flg_roi='ByTSA'
#flg_roi='ByLatLon'

if flg_roi=='ByTSA':
    
    t0=time.time()
    roi={}
    roi['Type']='ByTSA'
    
    # Pick the TSAs to include
    #roi['TSA List']=['Soo TSA']
    #roi['TSA List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
    roi['TSA List']=['Williams Lake TSA']
    #roi['TSA List']=['100 Mile House TSA']
    #roi['TSA List']=list(tsa['key']['Name'])
    
elif flg_roi=='ByLatLon':
    
    roi={}
    roi['Type']='ByLatLon'
    
    # Hanceville fire
    flg=0
    if flg==1:
        roi['Centre']=[-122.92,51.92]
        roi['Radius']=40*1000
        meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'
    
    # Elephant Hill fire
    flg=1
    if flg==1:
        roi['Centre']=[-121.15,51.15]
        roi['Radius']=45*1000 # metres
        meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Elephant Hill ROI'
    
    # Yung Lake, Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.07,51.308]
        roi['Radius']=10*1000 # metres

# Prepare region of interest
roi=bc1hau.DefineROI(roi,tsa,bm,road)
t1=time.time()
print((t1-t0)/60)

#%% Import rasters over ROI

# Always import these
lc2=bc1hau.Import_Raster_Over_ROI(meta,'lc2',roi)
btm=bc1hau.Import_Raster_Over_ROI(meta,'btm',roi)

flg_cut_yr=1
if flg_cut_yr==1:
    cut_yr=bc1hau.Import_Raster_Over_ROI(meta,'cut_yr',roi)

flg_bsr=1
if flg_bsr==1:
    bsr=bc1hau.Import_Raster_Over_ROI(meta,'bsr',roi)

flg_bgcz=1
if flg_bgcz==1:
    bgcz=bc1hau.Import_Raster_Over_ROI(meta,'bgcz',roi)

flg_wf=1
if flg_bgcz==1:
    wf=bc1hau.Import_Raster_Over_ROI(meta,'wf',roi)
    
flg_age=1
if flg_age==1:
    age1=bc1hau.Import_Raster_Over_ROI(meta,'age1',roi)
    
flg_sph=1
if flg_sph==1:
    sphlive=bc1hau.Import_Raster_Over_ROI(meta,'sphlive',roi)  
    sphdead=bc1hau.Import_Raster_Over_ROI(meta,'sphdead',roi)  

# Adjust roi mask to exlcude areas not burned
flg=1
if flg==1:
    ind=np.where(wf['grd']['Data']==0)
    roi['Mask']['Data'][ind]=0

#%% Import required vector geodatabases

# Wildfire permimiter
flg=0
if flg==1:
    wfp=qv.GetWildfirePerimiter(meta,2017,2022)
    wfp=bc1hau.ClipGDF_ByROI(wfp,roi)

# Non-ob reforestation polygons
flg=0
if flg==1:
    atup={}
    atup['gdf']=gpd.read_file(r'D:\Data\FCI_Projects\SummaryReforestationNonOb\Geospatial\atu_polygons.geojson')
    atup['gdf clip']=bc1hau.ClipGDF_ByROI(atup['gdf'],roi)
    atup['gdf overlay']=gpd.overlay(atup['gdf clip'],wfp,how='intersection')

# Import planting
flg=0
if flg==1:
    pl_from_op={}
    pl_from_op['Year Start']=2018
    pl_from_op['Year End']=2021
    pl_from_op['gdf']=qr.GetOpeningsWithPlanting(pl_from_op['Year Start'],pl_from_op['Year End'])
    pl_from_op['gdf']=bc1hau.ClipGDF_ByROI(pl_from_op['gdf'],roi)
    
    # Planting with spatial from AT layer
    pls={}
    pls['gdf']=bc1hau.GetPlantingWithinROI(2018,2021,roi)

flg=0
if flg==1:
    #t0=time.time()
    # This takes 17 min!!
    #vri={}
    #vri['gdf']=bc1hau.Get_Vectors_For_ROI(roi,'vri',1900,2200)
    #t1=time.time()
    #print(t1-t0)
    #vri['gdf2']=bc1hau.ClipGDF_ByROI(vri['gdf'],roi)
    
    vri['gdf'].plot()
    vri['gdf2'].plot()

# Import surveyed areas
flg=0
if flg==1:
    su={}
    su['gdf']=bc1hau.GetSurveyWithinROI(2018,2021,roi)

# Import openings within ROI
flg=0
if flg==1:
    t0=time.time()
    op={}
    op['gdf']=bc1hau.GetOpeningsWithinROI(roi)
    t1=time.time()
    print(t1-t0)
    
    op['gdf']=bc1hau.ClipGDF_ByROI(op['gdf'],roi)

#%% Plot ROI mask

plt.close('all')
fig,ax=bc1hau.Plot_ROI_Mask(meta,roi,lc2,bm)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].loc[atup['gdf overlay']['Year']>=2018].plot(ax=ax[0],facecolor=[0,0.8,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta['Paths']['Figures'] + '\\Planted areas','png',900)

#%% Plot BTM

def Plot_ROI_BTM(btm):
    
    plt.close('all')

    # Grid
    bin=np.unique(btm['Data1'])
    
    N_bin=bin.size
    N_hidden=1
    N_color=N_bin+N_hidden
    
    z1=np.ones(btm['Data1'].shape)
    for i in range(bin.size):
        z1[(btm['Data1']==bin[i])]=i
    
    z1[(roi['Mask']['Data']!=1)]=N_bin

    # Labels
    lab=[]
    for i in range(N_bin):
        try:
            lab.append(btm['lab1'][i])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_bin);
    for i in range(N_bin):
        cm.colors[i,0:3]=btm['cl1'][i,:]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1,clim=(0,N_color),extent=lc2['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)
    #cb=plt.colorbar(im,cax=ax[1])

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/N_bin,i/N_bin],'k-',linewidth=0.5)
    pos2=meta['Graphics']['pos2']
    pos2[1]=0.6
    pos2[3]=0.24
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BTM(btm)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta['Paths']['Figures'] + '\\BTM','png',300)


#%% Plot Burn Severity within the TSA mask

def Plot_ROI_BSR(bsr):
    
    z1=np.zeros(bsr['grd']['Data'].shape)
    ic=np.where(bsr['key']['Code']=='High')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=1
    ic=np.where(bsr['key']['Code']=='Medium')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=2
    ic=np.where(bsr['key']['Code']=='Low')[0]; z1[np.where(bsr['grd']['Data']==bsr['key']['ID'][ic])]=3   
    ind=np.where(z1==0);
    z1[ind]=3
    ind=np.where(roi['Mask']['Data']==0)
    z1[ind]=4
    
    lab=['High','Medium','Low','Unburned','']

    # Number of colours and number of colours excluded from colorbar
    N_color=5
    N_hidden=0

    # Colormap
    cm=np.vstack( ( (0.5,0,0,1),(1,0.25,0.25,1),(1,0.75,0.75,1),(0.96,0.96,0.96,1),(1,1,1,1) ) )
    cm=matplotlib.colors.ListedColormap(cm)

    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1,clim=(0,N_color),extent=lc2['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden-1),i/(N_color-N_hidden-1)],'k-',linewidth=0.5)
    pos2=meta['Graphics']['pos2']
    pos2[1]=0.8
    pos2[3]=0.14
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BSR(bsr)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta['Paths']['Figures'] + '\\Planted areas and BSR','png',300)


#%% PLOT BGC Zones within the TSA mask

def Plot_BGCZone_WithinROI(meta,bgcz):

    # Grid
    bin=np.unique(bgcz['grd']['Data'])
    
    N_bin=bin.size
    N_hidden=1
    N_color=N_bin+N_hidden
    
    z1=np.ones(bgcz['grd']['Data'].shape)
    for i in range(bin.size):
        z1[(bgcz['grd']['Data']==bin[i])]=i
    
    z1[(roi['Mask']['Data']!=1)]=N_bin

    # Labels
    lab=[]
    for i in range(N_bin):
        try:
            lab.append(bgcz['key'].ZONE[bgcz['key'].VALUE==bin[i]].values[0])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_bin)
    #cm=plt.cm.get_cmap('gray',N_bin)
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,N_color),extent=bgcz['grd']['Extent'],cmap=cm)
    #tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    #ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+0.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/N_bin,i/N_bin],'k-',linewidth=0.5)
    pos2=copy.copy(meta['Graphics']['pos2'])
    pos2[1]=0.6
    pos2[3]=0.24
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_BGCZone_WithinROI(meta,bgcz)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor=[1,1,1],edgecolor=[0.5,0,1],linewidth=1.25,label='Planting',alpha=0.25)

#%% PLOT age from VRI

def Plot_Age1_WithinROI(meta,age1):
    
    # Grid
    bw=20; bin=np.arange(0,220,bw); 
    z1=(bin.size)*np.ones( age1['grd']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( age1['grd']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( age1['grd']['Data']==0)]=i+1
    z1[(roi['Mask']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    N_color=bin.size+3
    N_hidden=3
    
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=age1['grd']['Extent'],cmap=cm)
    tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['pos2_long']);
    
    return fig,ax

plt.close('all')
fig,ax=Plot_Age1_WithinROI(meta,age1)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)

#%% Plot SPH from VRI

def Plot_sph_WithinROI(meta,sphlive):
    
    # Grid
    bw=500; bin=np.arange(0,2500,bw); 
    z1=(bin.size)*np.ones( sphlive['grd']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( sphlive['grd']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(sphlive['grd']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( sphlive['grd']['Data']==0)]=i+1
    z1[(roi['Mask']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    N_color=bin.size+3
    N_hidden=3
    
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=sphlive['grd']['Extent'],cmap=cm)
    #tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['pos2_long']);
    pos2=copy.copy(meta['Graphics']['pos2'])
    pos2[1]=0.6
    pos2[3]=0.24
    ax[1].set(position=pos2)
    
    return fig,ax

plt.close('all')
fig,ax=Plot_sph_WithinROI(meta,sphlive)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)
#gu.PrintFig(meta['Paths']['Figures'] + '\\sphlive','png',300)


#%% Plot SPH dead from VRI

def Plot_sph_WithinROI(meta,sphdead):
    
    # Grid
    bw=200; bin=np.arange(0,3200,200); 
    z1=(bin.size)*np.ones( sphdead['grd']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( sphdead['grd']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(sphlive['grd']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( sphdead['grd']['Data']==0)]=i+1
    z1[(roi['Mask']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    N_color=bin.size+3
    N_hidden=3
    
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=sphdead['grd']['Extent'],cmap=cm)
    tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['pos2_long']);
    
    return fig,ax

plt.close('all')
fig,ax=Plot_sph_WithinROI(meta,sphdead)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)
#gu.PrintFig(meta['Paths']['Figures'] + '\\sphdead','png',300)


#%% PLOT harvested year within the TSA mask

def Plot_ROI_HarvestYear(cut_yr):
    
    # Grid
    bw=5; bin=np.arange(1960,2025,bw); 
    z1=(bin.size)*np.ones( cut_yr['grd']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( cut_yr['grd']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( cut_yr['grd']['Data']==0)]=i+1
    z1[(roi['Mask']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    #cm=plt.cm.get_cmap('viridis',i)
    cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    N_color=bin.size+3
    N_hidden=3
    
    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,meta['Graphics']['figsize1'][0],meta['Graphics']['figsize1'][1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=cut_yr['grd']['Extent'],cmap=cm)
    tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta['Graphics']['pos1'],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta['Graphics']['pos2']);
    
    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_HarvestYear(cut_yr)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)

#%% Reforestation Monitoring Site Selection Maps

def Plot_ROI_Climate():

    zT_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
    zT=tsa['grd'].copy()
    zT['Data']=zT_tmp['Data']
    del zT_tmp
    zT=gis.ClipRaster(zT,roi['xlim'],roi['ylim'])  
    zT['Data']=zT['Data'].astype('float')/10
    gc.collect()    
    
    zW_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
    zW=tsa['grd'].copy()
    zW['Data']=zW_tmp['Data']
    del zW_tmp
    zW=gis.ClipRaster(zW,roi['xlim'],roi['ylim'])  
    zW['Data']=zW['Data'].astype('float')
    gc.collect()
    #plt.matshow(zW['Data']);plt.colorbar()

    clm={}
    clm['Tmin']=[-10,-10,-10,-5,-5,-5]
    clm['W']=[150,65+(150-65)/2,65,120,25+(120-27)/2,27]
    clm['Tmin buffer']=0.25
    clm['W buffer']=5

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
    ax.plot(zT['Data'][0::10,0::10].flatten(),zW['Data'][0::10,0::10].flatten(),'.',markerfacecolor=[0.85,0.85,0.85],markeredgecolor='None')
    iBurn=np.where( (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='High')[0] ]) | (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='Medium')[0] ]) )
    ax.plot(zT['Data'][iBurn].flatten()[0::10],zW['Data'][iBurn].flatten()[0::10],'.',markerfacecolor=[0.8,0.7,0.4],markeredgecolor='None')
    for i in range(len(clm['Tmin'])):
        ax.plot(clm['Tmin'][i],clm['W'][i],'s',markersize=15,markeredgecolor='k',mfc='None')
    ax.set(position=[0.14,0.14,0.8,0.8],xlim=[-14,2],ylim=[0,200],xlabel='Minimum monthly temperature (\circC)',ylabel='Soil water content (mm)')
    
    # Grid
    z1=np.zeros(zT['Data'].shape)
    for i in range(len(clm['Tmin'])):
        ind=np.where( (np.abs(zT['Data']-clm['Tmin'][i])<clm['Tmin buffer']) & (np.abs(zW['Data']-clm['W'][i])<clm['W buffer']) )
        z1[ind]=i
    ind=np.where(z1==0); z1[ind]=i+1
    L=8

    # Labels
    lab=['Cold/Wet','Cold/Mesic','Cold/Dry','Warm/Wet','Warm/Mesic','Warm/Dry','','']

    # Colormap
    #cm=plt.cm.get_cmap('viridis',6)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( ((0,0,0.5,1),(0.1,0.5,1,1),(0.5,0.75,1,1),(1,0.75,0.25,1),(1,0,0,1),(0.5,0,0,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,950,750)
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L),extent=zT['Extent'],cmap=cm)
    bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)
    #ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L,1),ticks=np.arange(0.5,L+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,L):
        ax[1].plot([0,100],[i/(L-1),i/(L-1)],'k-',linewidth=0.5)
    ax[1].set(position=[0.8,0.04,0.025,0.5])

    return

#%% Plot climate

def Plot_ROI_Climate():

    zT_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
    zT=tsa['grd'].copy()
    zT['Data']=zT_tmp['Data']
    del zT_tmp
    zT=gis.ClipRaster(zT,roi['xlim'],roi['ylim'])  
    zT['Data']=zT['Data'].astype('float')/10
    gc.collect()    
    
    zW_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
    zW=tsa['grd'].copy()
    zW['Data']=zW_tmp['Data']
    del zW_tmp
    zW=gis.ClipRaster(zW,roi['xlim'],roi['ylim'])  
    zW['Data']=zW['Data'].astype('float')
    gc.collect()
    #plt.matshow(zW['Data']);plt.colorbar()

    plt.close('all')
    plt.plot(zT['Data'][0::10,0::10].flatten(),zW['Data'][0::10,0::10].flatten(),'.',markerfacecolor=[0.85,0.85,0.85],markeredgecolor='None')
    iBurn=np.where( (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='High')[0] ]) | (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='Medium')[0] ]) )
    plt.plot(zT['Data'][iBurn].flatten()[0::10],zW['Data'][iBurn].flatten()[0::10],'.',markerfacecolor=[0.8,0.7,0.4],markeredgecolor='None')
    plt.plot(-10,150,'s',markersize=15,markeredgecolor='k',mfc='None')
    plt.plot(-10,65+(150-65)/2,'s',markersize=15,markeredgecolor='k',mfc='None')
    plt.plot(-10,65,'s',markersize=15,markeredgecolor='k',mfc='None')
    
    plt.plot(-5,120,'s',markersize=15,markeredgecolor='k',mfc='None')
    plt.plot(-5,25+(120-25)/2,'s',markersize=15,markeredgecolor='k',mfc='None')
    plt.plot(-5,25,'s',markersize=15,markeredgecolor='k',mfc='None')
    
    # Grid
    z1=np.zeros(zT['Data'].shape)
    # Cold and wet
    
    ind=np.where( (np.abs(zT['Data']--9.5)<0.25) & (np.abs(zW['Data']-100)<5) )
    z1[ind]=1
    # Cold and mesic
    ind=np.where( (np.abs(zT['Data']--9)<0.25) & (np.abs(zW['Data']-75)<5) )
    z1[ind]=2
    # Cold and dry
    ind=np.where( (np.abs(zT['Data']--8)<0.25) & (np.abs(zW['Data']-39)<5) )
    z1[ind]=3
    # Warm and wet
    ind=np.where( (np.abs(zT['Data']--5.5)<0.25) & (np.abs(zW['Data']-100)<5) )
    z1[ind]=4
    # Warm and mesic
    ind=np.where( (np.abs(zT['Data']--5.5)<0.25) & (np.abs(zW['Data']-75)<5) )
    z1[ind]=5
    # Warm and dry
    ind=np.where( (np.abs(zT['Data']--5.5)<0.25) & (np.abs(zW['Data']-39)<5) )
    z1[ind]=6
    L=8
    plt.matshow(zT['Data'])
    
    # Grid
    bin=np.arange(8,17,0.5); bw=0.5;
    z1=(bin.size)*np.ones(zC['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs(zC['Data']-bin[i])<=bw/2)
        z1[ind]=i
    #z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']!=4)]=i+1
    z1[(roi['Mask']['Data']!=1)]=i+2
    L=i+2

    # Labels
    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,950,750)
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zC['Extent'],cmap=cm)
    bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)
    #ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L,1),ticks=np.arange(0.5,L+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,L):
        ax[1].plot([0,100],[i/(L-1),i/(L-1)],'k-',linewidth=0.5)
    ax[1].set(position=[0.8,0.04,0.025,0.5])

    return

##%%
#
## Add burn seveirty rating 
#gdf_bsr_roi.plot(ax=ax[0],facecolor=[1,0,0],edgecolor=[1,0,0],linewidth=1,label='BSR',alpha=0.25)
#
#gdf_at_roi.plot(ax=ax[0],facecolor=[0,0.5,1],edgecolor=[0,0.5,1],linewidth=0.5,label='AT',hatch='////',alpha=0.25)
#
## Plot planted areas
#gdf_fci=gpd.read_file(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\atu_polygons.geojson')
#gdf_fci2=gdf_fci[gdf_fci['SILV_BASE_CODE']=='PL']
#gdf_fci2.plot(ax=ax[0],facecolor=[0.5,1,0],edgecolor=[0,0.5,0],linewidth=0.5,label='FCI',hatch='///',alpha=0.25)
#
## Surveys (zanzibar)
#dS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Surveys\ZANZIBAR FCI EXPORT.xlsx')
#import re
#dS['x']=np.zeros(dS['Latitude'].size)
#dS['y']=np.zeros(dS['Latitude'].size)
#for i in range(dS['Latitude'].size):
#    try:
#        lat=re.split('[°\'"]+',dS['Latitude'][i])
#        lon=re.split('[°\'"]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    try:
#        lat=re.split('[.]+',dS['Latitude'][i])
#        lon=re.split('[.]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    try:
#        lat=re.split('[\'"]+',dS['Latitude'][i])
#        lon=re.split('[\'"]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    
#points=[]
#for k in range(dS['x'].size):
#    points.append(Point(dS['x'][k],dS['y'][k])) 
#gdf_su=gpd.GeoDataFrame({'geometry':points,'Opening ID':dS['Opening ID']} )
#gdf_su=gdf_su.set_geometry('geometry')
#gdf_su.crs={'init':'epsg:4326'}
#gdf_su=gdf_su.to_crs(bm['gdf_bm'].crs)
#
#gdf_su.plot(ax=ax[0],marker='o',markersize=14,edgecolor='k',facecolor='w',linewidth=1)
#
## Surveys (All)
#dS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Surveys\FCI Post-wildfire Survey Summary.xlsx','2018')
#import re
#dS['x']=np.zeros(dS['Latitude'].size)
#dS['y']=np.zeros(dS['Latitude'].size)
#for i in range(dS['Latitude'].size):
#    try:
#        lat=re.split('[°\'"]+',dS['Latitude'][i])
#        lon=re.split('[°\'"]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    try:
#        lat=re.split('[.]+',dS['Latitude'][i])
#        lon=re.split('[.]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    try:
#        lat=re.split('[\'"]+',dS['Latitude'][i])
#        lon=re.split('[\'"]+',dS['Longitude'][i])
#        dS['y'][i]=np.array(lat[0],dtype=float)+np.array(lat[1],dtype=float)/60+np.array(lat[2],dtype=float)/(60*60)
#        dS['x'][i]=-np.array(lon[0],dtype=float)-np.array(lon[1],dtype=float)/60-np.array(lon[2],dtype=float)/(60*60)
#    except:
#        pass
#    
#points=[]
#for k in range(dS['x'].size):
#    points.append(Point(dS['x'][k],dS['y'][k])) 
#gdf_su=gpd.GeoDataFrame({'geometry':points,'Opening ID':dS['Opening ID']} )
#gdf_su=gdf_su.set_geometry('geometry')
#gdf_su.crs={'init':'epsg:4326'}
#gdf_su=gdf_su.to_crs(bm['gdf_bm'].crs)
#
#gdf_su.plot(ax=ax[0],marker='s',markersize=14,edgecolor='k',facecolor='w',linewidth=1)











#%% Plot AOS year within the TSA mask

dSC=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_PEST_SEVERITY_CODE.xlsx')
sco=['L','M','S','V']
tv=np.arange(1990,2021,1)
A=np.zeros((tv.size,len(sco)))
for iT in range(tv.size):
    print(tv[iT])
    z=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_IDL_SeverityClass_' + str(tv[iT]) + '.tif')
    for iSC in range(len(sco)):
        ind=np.where(dSC['PEST_SEVERITY_CODE']==sco[iSC])[0]
        ind=np.where(z['Data']==dSC['ID'][ind])[0]
        A[iT,iSC]=ind.size

plt.close('all')
plt.plot(tv,A[:,1],'-s')        
plt.plot(tv,A[:,2],'-o')
plt.plot(tv,A[:,3],'-d')


zLC2_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\lc2.tif')
zLC2=tsa['grd'].copy()
lc2['grd']['Data']=zLC2_tmp['Data']
del zLC2_tmp
zLC2=gis.ClipRaster(zLC2,xlim,ylim)
gc.collect()

# Grid
bw=5; bin=np.arange(1960,2025,bw); 
z1=(bin.size)*np.ones(zH['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs(zH['Data']-bin[i])<=bw/2)
    z1[ind]=i
z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']==4) & (zH['Data']==0)]=i+1
z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']!=4)]=i+2
z1[(roi['Mask']['Data']!=1)]=i+3
L=i+3

# Labels
lab=bin.astype(str)

# Colormap
#cm=plt.cm.get_cmap('viridis',i)
cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(0.8,0.8,0.8,1),(0.9,0.9,0.9,1),(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,950,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
ax[0].grid(False)
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L-1,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-2),i/(L-2)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.025,0.35])
gu.PrintFig(r'G:\My Drive\Figures\CutYear','png',500)






#%% Map of N deposition from 2017 aerial fertilization

pth=r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\Inputs\Nitrogen Deposition\ndep_2017.tif'
z0=gis.OpenGeoTiff(pth)

Mask=np.zeros(z0['Data'].shape)
id=np.where(zBTM.lab1=='Old Forest (FO)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Range Lands (RANG)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Recently Burned (BURN)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Recently Logged (LOG)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Selectively Logged (LOGS)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Shurbs (SHRB)')[0]; Mask[(zBTM['Data1']==id)]=1
id=np.where(zBTM.lab1=='Young Forest (FY)')[0]; Mask[(zBTM['Data1']==id)]=1

N0=np.sum(z0['Data'].flatten())/25000
N1=np.sum(z0['Data'].flatten()*Mask.flatten())/25000
N1/N0

N1=np.sum(z0['Data'].flatten()*Mask.flatten())
3.667*35*0.8*N1/1000/2500

#Mask=np.ones(z0['Data'].shape)

# Grid
L=8
z1=L*np.ones(tsa['grd']['Data'].shape,dtype='int8')
ind=np.where((Mask==1) & (z0['Data']>0.001) & (z0['Data']<=0.01)); z1[ind]=1
ind=np.where((Mask==1) & (z0['Data']>0.01) & (z0['Data']<=0.05)); z1[ind]=2
ind=np.where((Mask==1) & (z0['Data']>0.05) & (z0['Data']<=0.1)); z1[ind]=3
ind=np.where((Mask==1) & (z0['Data']>0.1) & (z0['Data']<=0.2)); z1[ind]=4
ind=np.where((Mask==1) & (z0['Data']>0.2) & (z0['Data']<=0.5)); z1[ind]=5
ind=np.where((Mask==1) & (z0['Data']>0.5) & (z0['Data']<=1)); z1[ind]=6
ind=np.where((Mask==1) & (z0['Data']>1) ); z1[ind]=7
z1[(Mask==0)]=8
z1[(z0['Data']==0)]=8
#id=np.where(zBTM.lab1=='Unknown')[0]; z1[(zBTM['Data1']==id)]=8
#id=np.where(zBTM.lab1=='Outside B.C. (NONE)')[0]; z1[(zBTM['Data1']==id)]=8

pth=r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\Inputs\Nitrogen Deposition\ndep_2017_cat.tif'
zo=z0.copy()
zo['Data']=z1.astype('int16')
gis.SaveGeoTiff(zo,pth)

# Labels
lab=['0.005','0.025','0.075','0.150','0.350','0.750','1.000']

# Colormap
cm=plt.cm.get_cmap('plasma',7)
cm=np.vstack( (cm.colors,(1,1,1,1)) ) # ,(1,1,1,1)
cm=matplotlib.colors.ListedColormap(cm)

from rasterio.plot import show
show(z1)
                              
# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,600,600)
im=ax[0].matshow(z1[0::50,0::50])
im=ax[0].matshow(z1[0::50,0::50],clim=(0,L),extent=tsa['grd']['Extent'],cmap=cm)
#tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.5)
ax[0].set(position=[0,0,1,1],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto',xticks=xlim,yticks=ylim)
ax[0].grid(False)
ax[0].tick_params(left=False, right=False)

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L-1,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=13,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-2),i/(L-2)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.54,0.03,0.3])
gu.PrintFig(r'G:\My Drive\Figures\Fertilization\ndep_from_aerial_fert_2017','png',500)

