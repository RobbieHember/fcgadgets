
#%% Import modules

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import PatchCollection
import geopandas as gpd
import pandas as pd
import copy
import fiona
import time
from scipy.interpolate import griddata
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1hau

#%% Path management

meta={}
meta['Paths']={}
meta['Paths']['BC1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
meta['Paths']['Forest Inventory Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
#meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation'

#%% Plotting parameters

meta['Graphics']={}
meta['Graphics']['figwidth']=16

#meta['Graphics']['sidespace']=0.25
meta['Graphics']['sidespace']=0

meta['Graphics']['ax1 pos']=[0,0,1-meta['Graphics']['sidespace']-0.01,1]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

gp=gu.SetGraphics('Manuscript')

#%% Import base maps

gdf=bc1hau.Import_GDBs_ProvinceWide()

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
#zBTM=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif')

#%% Forest health outbreak maps

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

#%% Plot mills

u=gdf['tpf']['gdf']['PRODUCT_CODE'].unique()
gdf['tpf']['gdf'].keys()

plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,14*zBTM['yxrat']))
#im=ax[0].matshow(z1,clim=(0,L+1),extent=zBTM['Extent'],cmap=cm)
gdf['bc_bound']['gdf'].plot(ax=ax[0],color=None,linewidth=0.25,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.05,0.05,0.9,0.9],xlim=zBTM['xlim'],ylim=zBTM['ylim'],aspect='auto')
ax[0].grid(False)
ax[0].axis('off')
ax[1].set(position=[0.06,0.15,0.05,0.65]);

mtypeL=['LBR','PLP','CHP','PLT']
ax=bc1hau.PlotMills(ax,gdf['tpf']['gdf'],mtypeL,labels='On')

#%% Find a ROI (example Revelstoke Complex)

flg=0
if flg==1:
    zBTM=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif')
    zMask=np.zeros(zBTM['Data'].shape)
    ind=np.where( (zBTM['Data']!=9) & (zBTM['Data']!=15) & (zBTM['Data']!=0) )
    zMask[ind]=1
    #plt.close('all')
    #plt.matshow(zMask)

    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,14*zBTM['yxrat']))
    im=ax.matshow(zMask,extent=zBTM['Extent'])

    ch_aoi=gpd.read_file(r'D:\Data\FCI_Projects\CaribouRecovery\Geospatial\Received 2022-09-22\REVY_HERDS.geojson')
    ch_aoi.plot(ax=ax)

#%% Plot site productivity

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\SPL\Site_Prod_Pl.tif')

plt.close('all')
plt.matshow(z['Data'],clim=[5,24])

#%% Plot probability of harvest

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\HarvestProbability.tif')

# Apply scale factor to convert to (%/yr)
z['Data']=z['Data'].astype('float')/100

bw=0.05; bin=np.arange(0,0.55,bw)
z1=np.ones( z['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs( z['Data']-bin[i])<=bw/2)
    if ind[0].size==0:
        z1[ind[0][i],ind[1][i]]=i
    else:
        z1[ind]=i
ind=np.where( z['Data']>=bin[i] ); z1[ind]=i
ind=np.where( zRef['Data']!=1 ); z1[ind]=i+1

lab=["%.2f" % x for x in bin]
lab=np.append(lab,np.array(['Water','Water']))

# Number of colours and number of colours excluded from colorbar
N_vis=bin.size+1
N_hidden=1
N_tot=N_vis+N_hidden

# Colormap
cm=plt.cm.get_cmap('viridis',bin.size)
cm=np.vstack( (cm.colors,(1,1,1,1)) )
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
ax[1].set(position=[0.71,0.6,0.045,0.32])

cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=11,length=0)
cb.outline.set_edgecolor('w')
for i in range(cbivls.size):
    ax[1].plot([0,100],[cbivls[i],cbivls[i]],'w-',linewidth=2)

gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest_Probability','png',300)

#%% Plot soil organic carbon

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soc_tot_forest_Shawetal2018.tif')

bw=20; bin=np.arange(0,320,bw)
z1=np.ones( z['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs( z['Data']-bin[i])<=bw/2)
    z1[ind]=i
ind=np.where( z['Data']>=bin[i] ); z1[ind]=i
ind=np.where( zMask!=1 ); z1[ind]=i+1

L=i+1

lab=bin.astype(str)

# Colormap
cm=plt.cm.get_cmap('viridis',i)
#cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

N_color=bin.size+1
N_hidden=1

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,14*zRef['yxrat']))
im=ax[0].matshow(z1,clim=(0,L+1),extent=z['Extent'],cmap=cm)
gdf['bc_bound']['gdf'].plot(ax=ax[0],color=None,linewidth=0.25,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.05,0.05,0.9,0.9],xlim=z['xlim'],ylim=z['ylim'],aspect='auto')
ax[0].grid(False)
ax[0].axis('off')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=10,length=0)
cb.outline.set_edgecolor('w')
for i in range(0,N_color):
    ax[1].plot([0,100],[i,i],'w-',linewidth=0.5)
ax[1].set(position=[0.06,0.15,0.05,0.65]);
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\SOC_Shawetal2018','png',300)

#%% Plot GlobBiomass

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
z['Data']=0.5*z['Data'] # Convert to tC/ha

bw=10; bin=np.arange(0,160,bw)
z1=np.ones( z['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs( z['Data']-bin[i])<=bw/2)
    z1[ind]=i
ind=np.where( z['Data']>=bin[i] ); z1[ind]=i
ind=np.where( zMask!=1 ); z1[ind]=i+1

L=i+1

lab=bin.astype(str)

# Colormap
cm=plt.cm.get_cmap('viridis',i)
#cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

N_color=bin.size+1
N_hidden=1

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,14*zRef['yxrat']))
im=ax[0].matshow(z1,clim=(0,L+1),extent=z['Extent'],cmap=cm)
gdf['bc_bound']['gdf'].plot(ax=ax[0],color=None,linewidth=0.25,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.05,0.05,0.9,0.9],xlim=z['xlim'],ylim=z['ylim'],aspect='auto')
ax[0].grid(False)
ax[0].axis('off')

cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color-N_hidden,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=10,length=0)
cb.outline.set_edgecolor('w')
for i in range(0,N_color):
    ax[1].plot([0,100],[i,i],'w-',linewidth=0.5)
ax[1].set(position=[0.06,0.15,0.05,0.65]);
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\GlobBiomass','png',300)

#%% Plot nitrogen deposition

#nam='1985'
nam='2016'

pth=r'C:\Users\rhember\Documents\Data\Nitrogen Deposition\Ackerman\ndep_tot_' + nam + '_clip_projected.tif'
z=gis.OpenGeoTiff(pth)
z['X']=z['X'][:,0:-1]

zND=tsa['grd'].copy()
ivl=10
zND['X']=zND['X'][0::ivl,0::ivl]
zND['Y']=zND['Y'][0::ivl,0::ivl]
zND['Data']=griddata( (z['X'].flatten(),z['Y'].flatten()) ,z['Data'].flatten(),(zND['X'],zND['Y']),method='cubic')

# plot

bw=0.2; bin=np.round(np.arange(0,4,bw),decimals=2);
z1=np.ones( zND['Data'].shape)
for i in range(bin.size):
    ind=np.where(np.abs( zND['Data']-bin[i])<=bw/2)
    z1[ind]=i
ind=np.where( zND['Data']>=bin[i] ); z1[ind]=i
ind=np.where( zMask[0::ivl,0::ivl]!=1 ); z1[ind]=i+1
#z1[(roi['Mask']['Data']==1) & ( age1['grd']['Data']==0)]=i+1
#z1[(roi['Mask']['Data']!=1)]=i+2
L=i+1

lab=bin.astype(str)

# Colormap
cm=plt.cm.get_cmap('viridis',i)
#cm=plt.cm.get_cmap('plasma',i)
cm=np.vstack( (cm.colors,(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

N_color=bin.size+1
N_hidden=1

# Plot
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(8,8*zRef['yxrat']))
im=ax[0].matshow(z1[0::1,0::1],clim=(0,L+1),extent=zND['Extent'],cmap=cm)
bm['gdf_bc_bound'].plot(ax=ax[0],color=None,edgecolor=[0,0,0],facecolor='none')
ax[0].set(position=[0.05,0.05,0.9,0.9],xlim=zND['xlim'],ylim=zND['ylim'],aspect='auto')
ax[0].grid(False)
ax[0].axis('off')
cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_color-(N_hidden-1),1),ticks=np.arange(0.5,N_color+1.5,1))
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=6,length=0)
cb.outline.set_edgecolor('w')
for i in range(0,N_color):
    ax[1].plot([0,100],[i/(N_color-N_hidden),i/(N_color-N_hidden)],'w-',linewidth=0.5)
ax[1].set(position=[0.84,0.35,0.04,0.6]);

gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\N Deposition\ndep_fromAckerman_' + nam,'png',300)

