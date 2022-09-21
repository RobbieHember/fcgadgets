
#%% Import modules

import os
import numpy as np
import gc
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import copy
import fiona
import time
from scipy.interpolate import griddata
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
#import fcgadgets.macgyver.utilities_query_gdbs as qv
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
meta['Graphics']['sidespace']=0.25
meta['Graphics']['figwidth']=16
meta['Graphics']['ax1 pos']=[0.04,0.02,1-meta['Graphics']['sidespace']-0.01,0.95]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

params_graphic=cbu.Import_GraphicsParameters('bc1ha_1')
plt.rcParams.update(params_graphic)

#%% Import base maps

gdf=bc1hau.Import_GDBs_ProvinceWide()

# Load dataset with CRS
#meta_bc1ha['crs']=gdf['bc_land']['gdf'].crs

#%% Get mask

zBTM=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif')
zMask=np.zeros(zBTM['Data'].shape)
ind=np.where( (zBTM['Data']!=9) & (zBTM['Data']!=15) & (zBTM['Data']!=0) )
zMask[ind]=1
#plt.close('all')
#plt.matshow(zMask)


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
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,14*zBTM['yxrat']))
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
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14,14*zBTM['yxrat']))
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
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(8,8*zBTM['yxrat']))
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

