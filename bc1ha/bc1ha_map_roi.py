
#%% Import modules

import os
import numpy as np
import gc
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.cbook import get_sample_data
from matplotlib.colors import LightSource
import geopandas as gpd
import pandas as pd
import copy
import fiona
import time
from shapely.geometry import Polygon,Point,LineString,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1hau
import fcgadgets.macgyver.utilities_query_gdb as qgdb

#%% Path management

meta_bc1ha={}
meta_bc1ha['Paths']={}
meta_bc1ha['Paths']['BC1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
meta_bc1ha['Paths']['Forest Inventory Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210930\Disturbances.gdb'
#meta_bc1ha['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation'
#meta_bc1ha['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS'
#meta_bc1ha['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\LICS Hanceville'

#%% Plotting parameters

meta_bc1ha['Graphics']={}
meta_bc1ha['Graphics']['sidespace']=0.45
meta_bc1ha['Graphics']['figwidth']=20
meta_bc1ha['Graphics']['ax1 pos']=[0,0,1-meta_bc1ha['Graphics']['sidespace'],1]
meta_bc1ha['Graphics']['ax1 vis']='off'
meta_bc1ha['Graphics']['ax1 gridvis']=False
meta_bc1ha['Graphics']['ax2 pos']=[1-meta_bc1ha['Graphics']['sidespace']+meta_bc1ha['Graphics']['ax1 pos'][0]+0.005,0.5,0.02,0.3]
meta_bc1ha['Graphics']['ax2 pos long']=[1-meta_bc1ha['Graphics']['sidespace']+meta_bc1ha['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

params_graphic=cbu.Import_GraphicsParameters('bc1ha_1')
plt.rcParams.update(params_graphic)

#%% Import base maps (full province)

gdf=bc1hau.Import_GDBs_ProvinceWide()

# Load dataset with CRS
meta_bc1ha['crs']=gdf['bc_land']['gdf'].crs

#%% Define region of interest

roi={}

# Note: roads can slow down by hour, set as [] to skip

#flg_roi='ByTSA'
flg_roi='ByLatLon'

t0=time.time()
if flg_roi=='ByTSA':

    roi['Type']='ByTSA'

    # Pick the TSAs to include
    # Search: gdf['tsa']['key']
    #roi['TSA List']=['Soo TSA']
    #roi['TSA List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
    roi['TSA List']=['Revelstoke TSA']
    #roi['TSA List']=['Kootenay Lake TSA']
    #roi['TSA List']=['100 Mile House TSA']
    #roi['TSA List']=['Merritt TSA','Kamloops TSA','100 Mile House TSA','Okanagan TSA','Williams Lake TSA','Lillooet TSA','Boundary TSA'] # ,'Arrow TSA','Revelstoke TSA'
    #roi['TSA List']=list(gdf['tsa']['key']['Name'])

elif flg_roi=='ByLatLon':

    roi['Type']='ByLatLon'

    # Hanceville fire
    flg=0
    if flg==1:
        roi['Centre']=[-122.92,51.92]
        roi['Radius']=40*1000
        meta_bc1ha['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Hanceville ROI'

    # Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.15,51.15]
        roi['Radius']=45*1000 # metres
        meta_bc1ha['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Elephant Hill ROI'

    # Yung Lake, Elephant Hill fire
    flg=0
    if flg==1:
        roi['Centre']=[-121.07,51.308]
        roi['Radius']=10*1000 # metres

    # Yahk
    flg=0
    if flg==1:
        roi['Centre']=[-116.086296,49.106646]
        roi['Radius']=100*1000 # metres

    # IWB
    flg=1
    if flg==1:
        roi['Centre']=[-117.5,50.8]
        roi['Radius']=100*1000 # metres

    # LICS
    flg=0
    if flg==1:
        lics=gu.ipickle(r'D:\Data\FCI_Projects\LICS Site List\LICS Sites.pkl')
        id=5; roi['Centre']=[lics[id]['Lon'],lics[id]['Lat']]; roi['Name']=lics[id]['Loc']
        roi['Radius']=15*1000 # metres

# Prepare region of interest
# Roads is empty because it takes a long time
roi=bc1hau.DefineROI(roi,gdf)
t1=time.time()
print((t1-t0)/60)

#%% Import rasters over ROI

#vList=['lc2','btm','elev','bgcz','cut_yr','bsr','wf','age1','sphlive','sphdead']
#vList=['temp_norm','ws_norm','lc2','btm','elev','bgcz','cut_yr']
#vList=['lc2','btm','elev','becz','age1','cut_yr']
vList=['lc2','btm','elev','soc','age1','si','temp_norm','ws_norm','cut_yr']
#vList=['cut_yr']
#vList=['si']
roi=bc1hau.Import_Raster_Over_ROI(meta_bc1ha,roi,vList)

#%% Import required vector geodatabases

vList=['op','ogsr']
roi=bc1hau.Import_GDB_Over_ROI(meta_bc1ha,roi,vList)

#%% Plot ROI mask

def Plot_ROI_Mask():
    plt.close('all')
    fig,ax=bc1hau.Plot_ROI_Mask(meta_bc1ha,roi,gdf)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=1.25,label='Road',alpha=1)
    roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1)
    roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=45)
    #wf['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,0],linewidth=1,label='Opening',alpha=1)
    roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=2,label='Opening',alpha=1)
    #fcres['gdf'].plot(ax=ax[0],facecolor=[0.75,0.9,0.75],edgecolor=[0,1,0],linewidth=1,label='Reserves',alpha=1,linestyle='--')
    #atu['gdf'].plot(ax=ax[0],facecolor=[0,1,1],edgecolor=[0,0,1],linewidth=1,label='Planting',alpha=0.25)
    #cc['gdf'].plot(ax=ax[0],facecolor=[1,0,0],edgecolor=[1,0,0],linewidth=1,label='Cut',alpha=0.25,linestyle='--')
    #gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\Planted areas','png',900)

#%% Plot elevation

def Plot_Elev(meta_bc1ha,roi):

    z=roi['grd']['elev']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['elev']['Data'],extent=roi['grd']['elev']['Extent'],cmap='terrain')

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')
    ax[0].grid(False)

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    #cb.ax.set(yticklabels=lab)
    #cb.ax.tick_params(labelsize=6,length=0)
    #for i in range(0,N_color):
    #    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_Elev(meta_bc1ha,roi)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=3,label='Road',alpha=1)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.75,label='Road',alpha=1)
roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=45)
#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,0],linewidth=0.5,label='Opening',alpha=1)
roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=0.5,label='Opening',alpha=1)

gp=gpd.read_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson')
ind=np.where( (sl['Age_t0']>125) & (sl['Csw_L_t0']<50) )[0]
gp.iloc[ind].plot(ax=ax[0],facecolor='None',marker='s',edgecolor=[1,0,0],linewidth=1,markersize=12,label='Opening',alpha=1)

ind=np.where( (sl['Age_t0']>125) & (sl['Csw_L_t0']>250) )[0]
gp.iloc[ind].plot(ax=ax[0],facecolor='None',marker='s',edgecolor=[0,1,0],linewidth=1,markersize=12,label='Opening',alpha=1)

a=gpd.read_file(r'C:\Users\rhember\Documents\Data\GroundPlots\DellaSala et al 2022 IWB\data\v10\outputs.gdb')
a=a.to_crs(roi['crs'])
a.plot(ax=ax[0],facecolor='None',marker='^',edgecolor=[1,1,0],linewidth=1.25,markersize=14,label='Opening',alpha=1)


#ind=np.where(roi['gdf']['op']['gdf']['OPENING_ID']==1760606)[0]
#roi['gdf']['op']['gdf'].iloc[ind].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,1],linewidth=2.5,label='Opening',alpha=1)

#wf['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,0],linewidth=1,label='Opening',alpha=1)
#fcres['gdf'].plot(ax=ax[0],facecolor=[0.75,0.9,0.75],edgecolor=[0,1,0],linewidth=1,label='Reserves',alpha=1,linestyle='--')
#atu['gdf'].plot(ax=ax[0],facecolor=[0,1,1],edgecolor=[0,0,1],linewidth=1,label='Planting',alpha=0.25)
#cc['gdf'].plot(ax=ax[0],facecolor=[0,0,1],edgecolor=[0,0,1],linewidth=1,label='Cut',alpha=0.25,linestyle='--')
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\Elevation','png',300)

#%% Plot SI

def Plot_SI(meta_bc1ha,roi):

    z=roi['grd']['si']['Data']
    if roi['Type']=='ByTSA':
        z[roi['grd']['Data']==0]=0

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(roi['grd']['si']['Data'],extent=roi['grd']['si']['Extent'],cmap='magma',clim=[5,22])

    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')
    ax[0].grid(False)

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['Cellsize'],roi['grd']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1])#,boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    #cb.ax.set(yticklabels=lab)
    #cb.ax.tick_params(labelsize=6,length=0)
    #for i in range(0,N_color):
    #    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_SI(meta_bc1ha,roi)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=1.25,label='Road',alpha=1)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=0.5,label='Road',alpha=1)
roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=45)
roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,1,0],linewidth=1,label='Opening',alpha=1)
roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,1,1],linewidth=1,label='Opening',alpha=1)

#%% Plot BTM

def Plot_ROI_BTM(roi):

    # Grid
    bin=np.unique(roi['grd']['btm']['Compressed']['Data'])

    N_bin=bin.size
    N_hidden=1
    N_color=N_bin+N_hidden

    z1=np.ones(roi['grd']['btm']['Compressed']['Data'].shape,dtype=int)
    for i in range(bin.size):
        z1[(roi['grd']['btm']['Compressed']['Data']==bin[i])]=i

    z1[(roi['grd']['Data']!=1)]=N_bin

    # Labels
    lab=[]
    for i in range(N_bin):
        try:
            lab.append(roi['grd']['btm']['Compressed']['lab'][i])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_bin);
    for i in range(N_bin):
        cm.colors[i,0:3]=roi['grd']['btm']['Compressed']['cl1'][i,:]
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_color),extent=roi['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    #roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf']['rivers'].plot(ax=ax[0],color=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf']['bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
#    try:
#        roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
#    except:
#        pass
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')
    ax[0].grid(False)
    #cb=plt.colorbar(im,cax=ax[1])

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['elev']['Cellsize'],roi['grd']['elev']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.4,clim=(np.min(hs),np.max(hs)))

    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-1,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'w-',linewidth=1.5)
    pos2=meta_bc1ha['Graphics']['ax2 pos']
    pos2[1]=0.5
    pos2[3]=0.4
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BTM(roi)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='y',linewidth=2.25,label='Road',alpha=1)
roi['gdf']['road'].plot(ax=ax[0],edgecolor='k',linewidth=1.5,label='Road',alpha=1)
roi['gdf']['tpf'].plot(ax=ax[0],marker='^',edgecolor='c',facecolor=[0.5,1,1],markersize=75)
roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,1,0],linewidth=0.5,label='Opening',alpha=1)
#op['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.8,0.8,1],linewidth=2,label='Wildfire',alpha=1)
#fcinv['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0.4,0.4,1],linewidth=1,label='Wildfire',alpha=1)
#fcres['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0.5,0.25],linewidth=1,label='Wildfire',alpha=1)
#atu['gdf'].plot(ax=ax[0],facecolor='None',ls='--',edgecolor=[1,0,0],linewidth=2,label='Planting',alpha=1)
#wf['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[1,0.5,0.5],linewidth=2,label='Planting',alpha=1)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\BTM_' + roi['Name'],'png',300)

#%% PLOT BGC Zones within the TSA mask

def Plot_BGCZone(roi):

    # Grid
    bin=np.unique(roi['grd']['bgcz']['Data'])

    N_bin=bin.size
    N_hidden=1
    N_color=N_bin+N_hidden

    z1=np.ones(grd['bgcz']['Data'].shape)
    for i in range(bin.size):
        z1[(grd['bgcz']['Data']==bin[i])]=i

    z1[(roi['Mask']['Data']!=1)]=N_bin

    # Labels
    lab=[]
    for i in range(N_bin):
        try:
            lab.append(grd['bgcz']['key'].ZONE[grd['bgcz']['key'].VALUE==bin[i]].values[0])
        except:
            lab.append('')

    # Colormap
    cm=plt.cm.get_cmap('viridis',N_bin)
    #cm=plt.cm.get_cmap('gray',N_bin)
    cm=np.vstack( (cm.colors,(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,N_color),extent=grd['bgcz']['Extent'],cmap=cm)
    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    try:
        roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    except:
        pass
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')

    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+0.5,1))
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=11,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'k-',linewidth=0.5)
    pos2=copy.copy(meta_bc1ha['Graphics']['ax2 pos'])
    pos2[1]=0.6
    pos2[3]=0.24
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_BGCZone(roi)
#wf['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[1,1,1],edgecolor=[0.5,0,1],linewidth=1.25,label='Planting',alpha=0.25)

#%% PLOT age from VRI

def Plot_Age1_WithinROI(roi):

    # Grid
    bw=20; bin=np.arange(0,220,bw);
    z1=(bin.size)*np.ones( roi['grd']['age1']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['age1']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    #ind=np.where(roi['grd']['age1']['Data']>bin[i]); z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['age1']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap
    cm=plt.cm.get_cmap('viridis',i)
    cm=np.vstack( (cm.colors,(0.9,0.9,0.9,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_color=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].grid(False)
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_Age1_WithinROI(roi)
roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],linewidth=0.5,facecolor='none')
#roi['gdf']['ogsr']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[1,0.5,0],linewidth=0.75,label='Opening',alpha=1)
#roi['gdf']['op']['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=0.75,label='Opening',alpha=1)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)

#%% PLOT soc

def Plot_SOC_WithinROI(roi):

    # Grid
    bw=20; bin=np.arange(0,320,bw);
    z1=(bin.size)*np.ones( roi['grd']['soc']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['soc']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['soc']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_SOC_WithinROI(roi)

#%% Plot Live SPH from VRI

def Plot_sph_WithinROI(meta_bc1ha,grd):

    # Grid
    bw=500; bin=np.arange(0,2500,bw);
    z1=(bin.size)*np.ones( grd['sphlive']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( grd['sphlive']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(grd['sphlive']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( grd['sphlive']['Data']==0)]=i+1
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=grd['sphlive']['Extent'],cmap=cm)
    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos long']);
    pos2=copy.copy(meta_bc1ha['Graphics']['pos2'])
    pos2[1]=0.6
    pos2[3]=0.24
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_sph_WithinROI(meta_bc1ha,grd)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\sphlive','png',300)

#%% Plot SPH dead from VRI

def Plot_sph_WithinROI(meta_bc1ha,grd):

    # Grid
    bw=200; bin=np.arange(0,3200,200);
    z1=(bin.size)*np.ones( grd['sphdead']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( grd['sphdead']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    ind=np.where(grd['sphdead']['Data']>bin[i])
    z1[ind]=i
    z1[(roi['Mask']['Data']==1) & ( grd['sphdead']['Data']==0)]=i+1
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=grd['sphdead']['Extent'],cmap=cm)
    roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    ax[0].set(position=meta_bc1ha['Graphics']['pos1'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['pos2_long']);

    return fig,ax

plt.close('all')
fig,ax=Plot_sph_WithinROI(meta_bc1ha,grd)
wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\sphdead','png',300)


#%% PLOT harvested year within the TSA mask

def Plot_HarvestYear(roi):

    # Grid
    bw=5; bin=np.arange(1960,2025,bw);
    z1=(bin.size)*np.ones( roi['grd']['cut_yr']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['cut_yr']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    #z1[(roi['grd']['Data']==1) & ( roi['grd']['cut_yr']['Data']==0)]=i+1
    #z1[(roi['grd']['Data']!=1)]=i+2
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
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf']['rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    try:
        roi['gdf']['roads'].plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    except:
        pass
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')
    ax[0].grid(False)

    #cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'k-',linewidth=0.5)
    ax[1].set(position=meta_bc1ha['Graphics']['ax2 pos']);

    return fig,ax

plt.close('all')
fig,ax=Plot_HarvestYear(roi)
#wf['gdf'].plot(ax=ax[0],facecolor='None',edgecolor=[0.85,0,0],lw=1,label='Wildfire',alpha=1)
#op['gdf'].plot(ax=ax[0],facecolor='None',ls='-',edgecolor=[0,0,0],lw=0.5,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.25,label='Planting',alpha=1)


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

    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1,clim=(0,N_color),extent=lc2['grd']['Extent'],cmap=cm)
    #bm['gdf_bc_bound'].plot(ax=ax[0],edgecolor=[0,0,0],facecolor='none',linewidth=0.25)
    roi['gdf_lakes'].plot(ax=ax[0],facecolor=[0.82,0.88,1],edgecolor=[0.7*0.82,0.7*0.88,0.7*1],linewidth=0.25,label='Water')
    roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    roi['gdf_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[0,0,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],aspect='auto')
    ax[0].grid(False)

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color,1),ticks=np.arange(0.5,N_color+1.5,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    for i in range(0,N_color):
        ax[1].plot([0,100],[i/(N_color-N_hidden-1),i/(N_color-N_hidden-1)],'k-',linewidth=0.5)
    pos2=meta_bc1ha['Graphics']['ax2 pos']
    pos2[1]=0.8
    pos2[3]=0.14
    ax[1].set(position=pos2)

    return fig,ax

plt.close('all')
fig,ax=Plot_ROI_BSR(bsr)
#wfp.plot(ax=ax[0],facecolor='None',edgecolor=[0,0,0],linewidth=1.5,label='Wildfire',alpha=1)
#atup['gdf overlay'].plot(ax=ax[0],facecolor=[0,0,0],edgecolor=[0,0.5,0],linewidth=1.25,label='Planting',alpha=0.25)
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\Planted areas and BSR','png',300)


#%% Plot soil water content

def Plot_SoilWaterContent(roi):

    # Grid
    bw=20; bin=np.arange(0,200+bw,bw);
    z1=(bin.size)*np.ones( roi['grd']['ws_norm']['Data'].shape)
    for i in range(bin.size):
        ind=np.where(np.abs( roi['grd']['ws_norm']['Data']-bin[i])<=bw/2)
        z1[ind]=i
    z1[(roi['grd']['Data']==1) & ( roi['grd']['ws_norm']['Data']==0)]=i+1
    z1[(roi['grd']['Data']!=1)]=i+2
    L=i+2

    lab=bin.astype(str)

    # Colormap

    cm=np.column_stack((roi['grd']['ws_norm']['cm']['cl1'],roi['grd']['ws_norm']['cm']['cl2'],roi['grd']['ws_norm']['cm']['cl3'],np.ones(roi['grd']['ws_norm']['cm']['bin'].size)))
    cm=np.vstack( (cm,(0.83137,0.81569,0.78431,1),(0.9,0.95,1,1),(1,1,1,1)) )
    cm=matplotlib.colors.ListedColormap(cm)

    N_color=bin.size+3
    N_hidden=3

    # Plot
    plt.close('all')
    fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta_bc1ha['Graphics']['figwidth'],(1-meta_bc1ha['Graphics']['sidespace'])*meta_bc1ha['Graphics']['figwidth']*roi['grd']['yxrat']))
    im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=roi['grd']['Extent'],cmap=cm)
    #roi['gdf']['tsa'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
    #roi['gdf_rivers'].plot(ax=ax[0],linecolor=[0,0,0.7],label='Water',linewidth=0.25)
    #try:
    #    roi['gdf_roads'].plot(ax=ax[0],facecolor='none',edgecolor=[1,1,0],label='Roads',linewidth=0.75,alpha=1,zorder=1)
    #except:
    #    pass
    ax[0].set(position=meta_bc1ha['Graphics']['ax1 pos'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'],xticklabels='',yticklabels='',aspect='auto')
    ax[0].grid(False)

    # Add relief shading
    z=roi['grd']['elev']['Data']; dx,dy=roi['grd']['elev']['Cellsize'],roi['grd']['elev']['Cellsize']
    ls=LightSource(azdeg=90,altdeg=45)
    ve=0.1
    hs=ls.hillshade(z,vert_exag=ve,dx=dx,dy=dy)
    ax[0].matshow(hs,extent=roi['grd']['elev']['Extent'],cmap='Greys',alpha=0.45,clim=(np.min(hs),np.max(hs)))

    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
    cb.ax.set(yticklabels=lab)
    cb.ax.tick_params(labelsize=6,length=0)
    cb.outline.set_edgecolor('w')
    for i in range(0,N_color):
        ax[1].plot([0,100],[i,i],'w-',linewidth=0.5)
    pos=[meta_bc1ha['Graphics']['ax2 pos'][0],0.5,0.03,0.4]
    ax[1].set(position=pos);

    #cb.ax.tick_params(labelsize=6,length=0)
    #for i in range(0,N_color):
    #    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'k-',linewidth=0.5)

    return fig,ax

plt.close('all')
fig,ax=Plot_SoilWaterContent(roi)
#gu.PrintFig(meta_bc1ha['Paths']['Figures'] + '\\SoilWaterContentNorm','png',300)


