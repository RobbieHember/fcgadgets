
#%% Import modules

import os
import numpy as np
import gc
import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import fiona
import time
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.query_vector_db as qv
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1ha

#%% Set figure properties

params_graphic=cbu.Import_GraphicsParameters('bc1ha_1')
plt.rcParams.update(params_graphic)

figsize1=[900,700]
pos1=[0.04,0.02,0.74,0.95]
pos2=[0.79,0.6,0.03,0.35]


#%% Import base maps

bm,tsa,road,district=bc1ha.Import_BaseMaps()

#%% Define region of interest

# By TSA
t0=time.time()
roi={}
roi['Type']='ByTSA'
# Pick the TSAs to include
#roi['TSA List']=['Soo TSA']
#roi['TSA List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
roi['TSA List']=['Williams Lake TSA']
#roi['TSA List']=list(tsa['key']['Name'])
roi=bc1ha.DefineROI(roi,tsa,bm,road)
t1=time.time()
print((t1-t0)/60)

# By Lat and Long
roi={}
roi['Type']='ByLatLon'
#roi['Centre']=[-123.1308, 51.8157]
#roi['Radius']=2500
roi['Centre']=[-123.107309,51.964553]
roi['Radius']=500*1000
roi=bc1ha.DefineROI(roi,tsa,bm,road)

#%% Import rasters over ROI

lc2=bc1ha.Import_Raster_Over_ROI('lc2',roi)

btm=bc1ha.Import_Raster_Over_ROI('btm',roi)

cut_yr=bc1ha.Import_Raster_Over_ROI('cut_yr',roi)

bsr=bc1ha.Import_Raster_Over_ROI('bsr',roi)

#%% Import geodatabases

#t0=time.time()

# This takes 17 min!!
#vri={}
#vri['gdf']=bc1ha.Get_Vectors_For_ROI(roi,'vri',1900,2200)
#t1=time.time()
#print(t1-t0)
#vri['gdf2']=bc1ha.ClipGDF_ByROI(vri['gdf'],roi)

vri['gdf'].plot()
vri['gdf2'].plot()

#%% Import planting

pl_from_op={}
pl_from_op['Year Start']=2018
pl_from_op['Year End']=2021
pl_from_op['gdf']=qr.GetOpeningsWithPlanting(pl_from_op['Year Start'],pl_from_op['Year End'])
pl_from_op['gdf']=bc1ha.ClipGDF_ByROI(pl_from_op['gdf'],roi)

# Planting with spatial from AT layer
pls={}
pls['gdf']=bc1ha.GetPlantingWithinROI(2018,2021,roi)


#%% Import surveyed areas

su={}
su['gdf']=bc1ha.GetSurveyWithinROI(2018,2021,roi)

#%% Import openings within ROI

t0=time.time()
op={}
op['gdf']=bc1ha.GetOpeningsWithinROI(roi)
t1=time.time()
print(t1-t0)

op['gdf']=bc1ha.ClipGDF_ByROI(op['gdf'],roi)


#%% Plot

fig,ax=Plot_BSR_WithinROI(bsr)

op['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.25,0.25,0.25],linewidth=0.75,label='Openings')
su['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.75,0.5,1],linewidth=2,label='Surveys')

pls['gdf'].plot(ax=ax[0],linestyle='--',facecolor='None',edgecolor=[0,0.6,0],linewidth=1.5,label='Planting')
pl_from_op['gdf'].plot(ax=ax[0],linestyle='--',facecolor='None',edgecolor=[0,0.6,0],linewidth=1.5,label='Planting from Opening')

#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\ReforestationPairedPlots\Hanceville1_BSR','png',300)


fig,ax=Plot_ROI_HarvestYear(cut_yr)
op['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.25,0.25,0.25],linewidth=0.75,label='Openings')
su['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.75,0.5,1],linewidth=2,label='Surveys')
pls['gdf'].plot(ax=ax[0],linestyle='--',facecolor='None',edgecolor=[0,0.6,0],linewidth=1.5,label='Planting')
pl_from_op['gdf'].plot(ax=ax[0],linestyle='--',facecolor='None',edgecolor=[0,0.6,0],linewidth=1.5,label='Planting from Opening')
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\ReforestationPairedPlots\Hanceville1_HarvestYear','png',300)


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
    mngr.window.setGeometry(100,100,figsize1[0],figsize1[1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=cut_yr['grd']['Extent'],cmap=cm)
    tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=pos1,xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=pos2);
    
    return fig,ax

#%% Plot Burn Severity within the TSA mask

def Plot_BSR_WithinROI(bsr):

    #z=np.zeros(zLC2['Data'].shape,dtype='int8')
    #ind=np.where( (bs['FIRE_YEAR']==yr) & (bs['BURN_SEVERITY_RATING']==lut['LUT BS']['BURN_SEVERITY_RATING']['High']) )[0]
    #ind=np.unravel_index(bs['IdxToGrid'][ind],z.shape,'C'); z[ind]=1
    #ind=np.where( (bs['FIRE_YEAR']==yr) & (bs['BURN_SEVERITY_RATING']==lut['LUT BS']['BURN_SEVERITY_RATING']['Medium']) )[0]
    #ind=np.unravel_index(bs['IdxToGrid'][ind],z.shape,'C'); z[ind]=2
    #ind=np.where( (bs['FIRE_YEAR']==yr) & (bs['BURN_SEVERITY_RATING']==lut['LUT BS']['BURN_SEVERITY_RATING']['Low']) )[0]
    #ind=np.unravel_index(bs['IdxToGrid'][ind],z.shape,'C'); z[ind]=3
    #zBS=z.copy(); del z;
    ## Grid and labels
    #lab=[]
    #z1=3*np.ones(zLC2['Data'].shape)
    #z1[(zBS==1)]=0; lab.append('High' )
    #z1[(zBS==2)]=1; lab.append('Medium')
    #z1[(zBS==3)]=2; lab.append('Low')
    
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

    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,figsize1[0],figsize1[1])
    im=ax[0].matshow(z1,clim=(0,N_color),extent=lc2['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=pos1,xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=pos2)

    return fig,ax

#%% PLOT ROI mask

def Plot_ROI_Mask():

    # Grid and labels
    lab=[]
    z1=np.ones(lc2['grd']['Data'].shape)
    z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']==4)]=0; lab.append('Treed')
    z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']!=4)]=1; lab.append('Non-treed')
    z1[(roi['Mask']['Data']!=1) & (lc2['grd']['Data']!=1)]=2; lab.append('hidden')
    z1[(roi['Mask']['Data']!=1) & (lc2['grd']['Data']==1)]=3; lab.append('hidden')

    # Number of colours and number of colours excluded from colorbar
    N_color=4
    N_hidden=2

    # Colormap
    cm=np.vstack( ((0.7,0.7,0.7,1),(0.8,0.8,0.8,1),(0.93,0.93,0.93,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,figsize1[0],figsize1[1])
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,N_color),extent=lc2['grd']['Extent'],cmap=cm)
    bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=pos1,xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-1,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=pos2)

    return fig,ax


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







#def GetBSR():
#
#    pth=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'
#    lyr='VEG_BURN_SEVERITY_SP'
#    
#    gl=[None]*int(1e5)
#    cnt=0
#    with fiona.open(pth,layer=lyr) as source:
#        for feat in source:
#            prp=feat['properties']
#            if (feat['geometry']==None):
#                continue
#            if (prp['FIRE_YEAR']!=2017):
#                continue
#            if (prp['BURN_SEVERITY_RATING']=='Unburned') | (prp['BURN_SEVERITY_RATING']=='Low'):
#                continue
#            gl[cnt]=feat
#            cnt=cnt+1
#
#    # Truncate
#    gl=gl[0:cnt-1]
#
#    # Convert to geodataframe
#    gdf_bsr=gpd.GeoDataFrame.from_features(gl,crs=bm['gdf_bm'].crs)
#
#    # Isolate features that intersect ROI
#    gdf_bsr_roi=gpd.overlay(gdf_bsr,tsa['gdf'][np.isin(tsa['gdf'].Name,tsaList)],how='intersection')
#
#    # Convert to geographic coordinate system
#    #gdf_at=gdf_at.to_crs({'init':'epsg:4326'})
#
#    # Save
#    #gdf_at.to_file(filename=meta['Paths']['Project'] + '\\Geospatial\\fcinv.geojson',driver='GeoJSON')
#
#    return


#%% PLOT BGC Zones within the TSA mask

def Plot_BGCZone_WithinROI():

    # BGC classification
    zBGC_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif')
    zBGC=tsa['grd'].copy()
    zBGC['Data']=zBGC_tmp['Data']
    del zBGC_tmp
    zBGC['key']=pd.read_excel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif.vat.xlsx')
    zBGC=gis.ClipRaster(zBGC,xlim,ylim)
    gc.collect()

    # Grid
    u=np.unique(zBGC['Data'])
    z1=np.ones(zBGC['Data'].shape)
    for i in range(u.size):
        z1[(zBGC['Data']==u[i])]=i
    L=i+1
    z1[(roi['Mask']['Data']!=1)]=L

    # Labels
    lab=[]
    for i in range(len(u)):
        try:
            lab.append(zBGC['key'].ZONE[zBGC['key'].VALUE==u[i]].values[0])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',L)
    #cm.colors=np.concatenate((cm.colors,np.array([1,1,1,1]).reshape(1,4)),axis=0)

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2)
    mngr=plt.get_current_fig_manager()
    mngr.window.setGeometry(100,100,750,750)
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zBGC['Extent'],cmap=cm)
    tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
    #ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L+1,1),ticks=np.arange(0.5,L+0.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    for i in range(0,L):
        ax[1].plot([0,100],[i/L,i/L],'k-',linewidth=0.5)
    ax[1].set(position=[0.06,0.04,0.025,0.4])

    return fig,ax




## Add invenotry data
#
## Inventory data
#sgrd=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\Coordinates.pkl')
#atu=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\RSLT_ACTIVITY_TREATMENT_SVW.pkl')
#pest=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\PEST_INFESTATION_POLY.pkl')
#vri=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupByTile1\TileBC_0909\Inputs\Geospatial\VEG_COMP_LYR_R1_POLY.pkl')
#lut_at=gu.ipickle(r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20200430\LUTs_RSLT_ACTIVITY_TREATMENT_SVW.pkl')
#lut_pest=gu.ipickle(r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20200430\LUTs_PEST_INFESTATION_POLY.pkl')
#xg=sgrd['X'].flatten()
#yg=sgrd['Y'].flatten()
#
#
#ind=np.where( (atu['SILV_BASE_CODE']==lut_at['SILV_BASE_CODE']['PL']) & (atu['Year']>=2017) )[0]
#
#ax[0].plot(xg[atu['IdxToGrid'][ind]],yg[atu['IdxToGrid'][ind]],'o',markeredgecolor=[0,1,0],markerfacecolor=[0.7,1,0.7],markersize=3,linewidth=0.25)




#%% PLOT stand age within the TSA mask

zAge_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\age1.tif')
zAge=tsa['grd'].copy()
zAge['Data']=zAge_tmp['Data']
zAge=gis.ClipRaster(zAge,xlim,ylim)
del zAge_tmp
gc.collect()

# Grid
bin=np.arange(0,210,10); bw=10;
z1=(bin.size)*np.ones(zAge['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs(zAge['Data']-bin[i])<=bw/2)
    z1[ind]=i
z1[(roi['Mask']['Data']==1) & (lc2['grd']['Data']!=4)]=i+1
z1[(roi['Mask']['Data']!=1)]=i+2
L=i+2

ind=np.where( (roi['Mask']['Data']==1) & (lc2['grd']['Data']==4) )
np.mean(zAge['Data'][ind])

# Labels
lab=bin.astype(str)

# Colormap
#cm=plt.cm.get_cmap('viridis',i)
cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
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

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L,1),ticks=np.arange(0.5,L+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/(L-1),i/(L-1)],'k-',linewidth=0.5)
ax[1].set(position=[0.8,0.04,0.025,0.5])





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




#%% PLOT wildfire within the TSA mask

zF1=tsa['grd'].copy()
zF1=gis.ClipRaster(zF1,xlim,ylim)
zF1['Data']=0*zF1['Data']
for i in range(1,5):
    zF1_tmp=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_year' + str(i) + '.tif')    
    zF1_tmp=gis.ClipRaster(zF1_tmp,xlim,ylim)
    ind=np.where(zF1_tmp['Data']==2017)
    zF1['Data'][ind]=zF1['Data'][ind]+zF1_tmp['Data'][ind]
    del zF1_tmp    
    gc.collect()

# Grid
z1=np.ones(zF1['Data'].shape)
z1[(zF1['Data']==2017)]=2
L=3
z1[(roi['Mask']['Data']!=1)]=L

# Labels
lab=[]
for i in range(len(u)):
    try:
        lab.append(zF1['key'].ZONE[zF1['key'].VALUE==u[i]].values[0])
    except:
        lab.append('')

# Colormap
cm=plt.cm.get_cmap('viridis',L)
#cm.colors=np.concatenate((cm.colors,np.array([1,1,1,1]).reshape(1,4)),axis=0)
cm.colors[0,:]=[0.8,0.8,0.8,1]
cm.colors[1,:]=[0.5,0,0,1]
cm.colors[2,:]=[1,1,1,1]

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2)
mngr=plt.get_current_fig_manager()
mngr.window.setGeometry(100,100,750,750)
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zF1['Extent'],cmap=cm)
tsa['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto')
#ax[0].set(position=[0.04,0.02,0.92,0.96],xlim=[tsa['grd'].minx,tsa['grd'].maxx],ylim=[tsa['grd'].miny,tsa['grd'].maxy],aspect='auto')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,L+1,1),ticks=np.arange(0.5,L+0.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=11,length=0)
for i in range(0,L):
    ax[1].plot([0,100],[i/L,i/L],'k-',linewidth=0.5)
ax[1].set(position=[0.06,0.04,0.025,0.4])





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

